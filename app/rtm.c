///
/// @copyright Copyright 2017- Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// \b simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details->
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author Issam Said.
/// @file src/main.c
/// @brief Main program of the wave simulator.
///
/// Contains the main program that runs the Reverse Time Migration (RTM) using
/// simwave.
///
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <cuda_runtime.h>
#include <simwave/parser.h>
#include <simwave/simwave.h>
#include <simwave/wave_tb.h>
#include <simwave/wtime.h>
#include <simwave/macros.h>
void wave_tb_data_close_open(tb_data_t * data,
                             const int nb_thread_groups,
const int shotid) {
  char tmp[512];
  sprintf(tmp, "mkdir -p %s", OUTDIR);
  system(tmp);
  sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SNAP_BASE, shotid);

  for (int i = 0; i < nb_thread_groups; i++) {
    fclose(data->gfd[i]);
    CHK(data->gfd[i] == NULL, "failed to open snapshot file");
  }

  for (int i = 0; i < nb_thread_groups; i++) {
    data->gfd[i] = fopen(tmp,"wb+");
    CHK(data->gfd[i] == NULL, "failed to open snapshot file");
  }
}

void cmemcpy_omp(char *dst,
                 char *src,
                 const size_t size) {
  size_t my_start, my_size;
  int tid,nth;
  #pragma omp parallel default(shared) private(tid,nth,my_start,my_size)
  {
    tid = omp_get_thread_num();
    nth = omp_get_num_threads();

    my_start = (tid*size)/nth;
    my_size = ((tid+1)*size)/nth - my_start;

    memcpy(dst + my_start,src + my_start, my_size);
  }
}

void memcpy_omp(float* u, float *v,sismap_t*s) {
  const int nnx = s->dimx + 2* s->sx;
  const int nny = s->dimy + 2* s->sy;
  const int nnxy = nnx*nny;
  float * restrict ux;
  float * restrict vx;
  #pragma omp parallel for collapse(2) private(ux,vx)
  for (int zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
    for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
      int zmax = zmin+BLOCKZ;
      if (zmax > s->dimz) zmax = s->dimz;
      int ymax = ymin+BLOCKY;
      if (ymax > s->dimy) ymax = s->dimy;
      for(int z = zmin; z < zmax ; z++) {
        for(int y = ymin; y < ymax; y++) {
          ux = &(  u[1ULL* (z+s->sz) * nnxy            + (y+s->sy) * nnx + s->sx]);
          vx = &(  v[1ULL* (z+s->sz) * nnxy            + (y+s->sy) * nnx + s->sx]);
//          #pragma simd
          for(int x = 0; x < s->dimx; x++) {
            ux[x] = vx[x];
          }
        }
      }
    }
  }
}

///
/// Reverse Time Migration on CPU.
///
///
void run_rtm_cpu (sismap_t *s, float* vel,  float *source, float *pml_tab) {
  /// contains the fields pressure value at time step t.
  float* u0;
  /// contains the fields pressure value at time step t+1.
  float* u1;
  /// forward wave-field for RTM.
  float* fwd;
  float** fwd_all;
  /// seismic traces for a given shot.
  float *sismos;
  /// image and illumination of each shot.
  float *ilm_shot, *img_shot;
  /// PML tmp tab.
  float *pml_tmp;
  /// wave-field arrays.
  CREATE_BUFFER_ONLY(u0, s->size);
  array_openmp_init(u0,s);
  CREATE_BUFFER_ONLY(u1, s->size);
  array_openmp_init(u1,s);
  CREATE_BUFFER_ONLY(fwd, s->size);
  array_openmp_init(fwd,s);

  CREATE_BUFFER(sismos, s->rcv_len*(s->time_steps+1));
  CREATE_BUFFER_ONLY(img_shot, s->size_img);
  array_openmp_inner_init(img_shot,s);

  CREATE_BUFFER_ONLY(ilm_shot, s->size_img);
  array_openmp_inner_init(ilm_shot,s);
  CREATE_BUFFER(pml_tmp, s->size_eff);
  shot_t *shot;

  double t0,t1,t2,t_snap,t_prop,t_sismos,t_image;
  wtime_init();

  int nb_fwd_snap = (s->time_steps+1 + s->nb_snap - 1) / s->nb_snap;
  printf("s->time_steps %d, nb_snap %d\n",s->time_steps, s->nb_snap);
  printf("Nb of snapshots in FWD: %d\n", nb_fwd_snap);

  if (s->mode == 2) {
    fwd_all = calloc(nb_fwd_snap,sizeof(float*));
    for (int i=0;i<nb_fwd_snap;i++) {
      CREATE_BUFFER_ONLY(fwd_all[i],s->size);
      array_openmp_init(fwd_all[i],s);
    }
    printf("MODE MEM\n");
  } else {
    printf("MODE I/O\n");
  }

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    MSG("Start of Shot (%d)",sidx);

    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, true, s->modeling);
    /// load the seismic traces for the shot.
    wave_read_sismos(s, shot, sismos);

    wave_min_max("sismos", sismos, s->rcv_len*(s->time_steps+1));

    /// reset buffers for the shot (forward).
    s->snap_idx = 0;
    NULIFY_BUFFER(u0, s->size);
    NULIFY_BUFFER(u1, s->size);
    NULIFY_BUFFER(pml_tmp, s->size_eff);
    /// forward modeling.

    t1 = wtime();
    t_snap   = 0.0;
    t_prop   = 0.0;
    t_sismos = 0.0;
    t_image  = 0.0;

    t0 = wtime();
/*
    if (s->mode == 2) {
        if ((0)%s->nb_snap == 0) {
          memcpy_omp(fwd_all[s->snap_idx],u0,s);
          if (s->verbose) MSG("... saving snapshot %d t=%u", s->snap_idx, 0);
          s->snap_idx ++;
        }
    } else {
      wave_save_snapshot(s, shot, u0, 0);
    }
    */
    t_snap += wtime() - t0;

    wave_update_source(s, shot, u0, source[0]);
    for(int t = 0; t < s->time_steps; ++t) {
      t0 = wtime();
      wave_update_fields_block_bis(s, u0, u1, vel, pml_tmp, pml_tab);
      t_prop += wtime() - t0;

      t0 = wtime();
      if (s->mode == 2) {
        if ((t+1)%s->nb_snap == 0) {
          memcpy_omp(fwd_all[s->snap_idx],u1,s);
          if (s->verbose) MSG("... saving snapshot %d t=%u", s->snap_idx, t+1);
          s->snap_idx ++;
        }
      } else {
        wave_save_snapshot(s, shot, u1, t+1);
      }
      t_snap += wtime() - t0;

      t0 = wtime();
      wave_update_source(s, shot, u1, source[t+1]);
      t_sismos += wtime() - t0;

      WAVE_SWAP_POINTERS(u0, u1);
    }
    t2 = wtime();

    MSG("forward timer (%d)",sidx);
    MSG("Total:        %f (s)",t2-t1);
    MSG("SNAP:         %f (s)",t_snap);
    MSG("PROP:         %f (s)",t_prop);
    MSG("SISMOS:       %f (s)",t_sismos);
    MSG("Speed:        %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t2-t1));
    MSG("PropSpeed:    %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t_prop));


    /// reset buffers for the shot (backward).
    s->snap_idx = s->snap_idx-1;
    NULIFY_BUFFER(u0,       s->size);
    NULIFY_BUFFER(u1,       s->size);
    NULIFY_BUFFER(pml_tmp,  s->size_eff);
    NULIFY_BUFFER(ilm_shot, s->size_img);
    NULIFY_BUFFER(img_shot, s->size_img);

    t1 = wtime();
    t_snap   = 0.0;
    t_prop   = 0.0;
    t_sismos = 0.0;
    t_image  = 0.0;

    /// backward modeling and imaging.
    t0 = wtime();
    wave_inject_sismos(s, u0, s->time_steps, sismos);
    t_sismos += wtime() - t0;

    for(int t = s->time_steps-1; t>=0 ; --t) {
      t0 = wtime();
      if (s->mode == 2) {
        if ((t+1)%s->nb_snap == 0) {
          memcpy_omp(fwd,fwd_all[s->snap_idx],s);
          if (s->verbose) MSG("... reading snapshot %d t=%u", s->snap_idx, t+1);
          s->snap_idx --;
        }
    }  else {
        wave_read_snapshot(s, shot, fwd, t+1);
      }
      t_snap += wtime() - t0;

      t0 = wtime();
      wave_image_condition_block(s, u0, fwd, img_shot, ilm_shot, t+1);
      t_image += wtime() - t0;

      t0 = wtime();
      wave_update_fields_block_bis(s, u0, u1, vel, pml_tmp, pml_tab);
      t_prop += wtime() - t0;

      t0 = wtime();
      wave_inject_sismos(s, u1, t, sismos);
      t_sismos += wtime() - t0;

      WAVE_SWAP_POINTERS(u0, u1);
    }


    t2 = wtime();

    MSG("backward timer (%d)",sidx);
    MSG("Total:        %f (s)",t2-t1);
    MSG("SNAP:         %f (s)",t_snap);
    MSG("PROP:         %f (s)",t_prop);
    MSG("SISMOS:       %f (s)",t_sismos);
    MSG("IMAGE COND:   %f (s)",t_image);
    MSG("Speed:        %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t2-t1));
    MSG("PropSpeed:    %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t_prop));

    /// save the img/ilm of the shot:
    wave_save_img(s, shot, img_shot, ilm_shot);
    /// release/close the resources related to the current shot.
    shot_release(shot);

    MSG("End of Shot (%d)",sidx);
  }

  if (s->mode == 2) {
  // TODO DEALLOCATE FWD_ALL
    DELETE_BUFFER(fwd_all);
  }
  /// free the simulation buffers.
  DELETE_BUFFER(u0);
  DELETE_BUFFER(u1);
  DELETE_BUFFER(fwd);
  DELETE_BUFFER(img_shot);
  DELETE_BUFFER(ilm_shot);
  DELETE_BUFFER(sismos);
  DELETE_BUFFER(pml_tmp);
}

///
/// Reverse Time Migration on CPU.
///
///
void run_rtm_tb_cpu (sismap_t *s, float* vel,  float *source, float *pml_tab,
                     parser* p) {
  /// contains the fields pressure value at time step t.
  float* u0;
  /// contains the fields pressure value at time step t+1.
  float* u1;
  /// forward wave-field for RTM.
  float* fwd,*fwd_io;
  /// seismic traces for a given shot.
  float *sismos;
  /// image and illumination of each shot.
  float *ilm_shot, *img_shot;
  /// PML tmp tab.
  float *pml_tmp;
  /// wave-field arrays.
printf("allcoation\n");
CREATE_BUFFER(u0,  s->size);
  CREATE_BUFFER(u1,  s->size);
  CREATE_BUFFER(img_shot, s->size_img);
  CREATE_BUFFER(ilm_shot, s->size_img);
  CREATE_BUFFER(pml_tmp, s->size_eff);
  shot_t *shot;
printf("all\n");
  tb_t * ctx         = (tb_t*)       malloc(sizeof(tb_t));
  tb_data_t * data   = (tb_data_t*)  malloc(sizeof(tb_data_t));
  tb_timer_t * timer = (tb_timer_t*) malloc(sizeof(tb_timer_t));
  wtime_init();

  // setup tb_data
  printf("tb_init\n");
  wave_tb_init(ctx,s,p);
//    source = realloc(source,sizeof(float)*(s->time_steps+1));

  printf("tb_info\n");
  wave_tb_info(ctx);
  printf("tb_timer\n");
  wave_tb_timer_init(timer,ctx->thread_group_size,ctx->num_thread_groups);

  printf("Nb of snapshots in FWD: %d\n",((ctx->t_len + ctx->fwd_steps-1) / ctx->fwd_steps));

  if (ctx->mode == 1 || ctx->mode == 2) {
    CREATE_BUFFER(fwd, 1ULL*s->size * ((ctx->t_len + ctx->fwd_steps-1) / ctx->fwd_steps));
  } else {
    printf("nnx %d, nnz %d, diam_width %d ntg %d\n",ctx->nnx,ctx->nnz, ctx->diam_width, ctx->num_thread_groups);
    CREATE_BUFFER(fwd, 1ULL*ctx->nnx * ctx->nnz * ctx->diam_width* ctx->num_thread_groups);
  }

  printf("rcv_len %d, time_steps %d\n",s->rcv_len,s->time_steps);
  CREATE_BUFFER(sismos, s->rcv_len*(s->time_steps+1));

  MSG("loop over the shots between %d and %d",s->first,s->last);

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    MSG("Processing shot %d (BEGIN)",sidx);
    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, true, s->modeling);
    /// load the seismic traces for the shot.
    wave_read_sismos(s, shot, sismos);

    wave_min_max("sismos", sismos, s->rcv_len*(s->time_steps+1));

    /// reset buffers for the shot (forward).
    NULIFY_BUFFER(u0, s->size);
    NULIFY_BUFFER(u1, s->size);
    if (ctx->mode == 1 || ctx->mode == 2) {
      NULIFY_BUFFER(fwd,s->size * ((ctx->t_len + ctx->fwd_steps-1) / ctx->fwd_steps));
      fwd_io = NULL;
    } else {
      NULIFY_BUFFER(fwd,1ULL*ctx->nnx * ctx->nnz * ctx->diam_width* ctx->num_thread_groups);
    }
    NULIFY_BUFFER(pml_tmp, s->size_eff);

    wave_tb_data_init(data,ctx,s,ctx->num_thread_groups,shot->id,
                      1ULL*ctx->nnx * ctx->nnz * ctx->diam_width);

    data->fwd    = fwd;
    data->img    = img_shot;
    data->ilm    = ilm_shot;

    /////////
    // FWD //
    /////////
    MSG("FWD STEP %d",sidx);
    wave_tb_timer_clear(timer);
    data->flag_fwd = 1;
    data->flag_bwd = 0;

    /// forward modeling.
    wave_tb_data_set_src(data,s,shot->srcidx,source);

    wave_tb_forward(ctx,data,timer,u0,u1,vel);

    wave_tb_data_unset_src(data);

    wave_tb_timer_info(timer,ctx->nb_stencils_total_fwd,ctx->nb_stencils_main);

//wave_tb_data_close_open(data,ctx->num_thread_groups,shot->id);

    /////////
    // BWD //
    /////////
    MSG("BWD STEP %d",sidx);
    wave_tb_timer_clear(timer);
    data->flag_fwd = 0;
    data->flag_bwd = 1;

    /// reset buffers for the shot (backward).
    NULIFY_BUFFER(u0,       s->size);
    NULIFY_BUFFER(u1,       s->size);
    NULIFY_BUFFER(pml_tmp,  s->size_eff);
    NULIFY_BUFFER(ilm_shot, s->size_img);
    NULIFY_BUFFER(img_shot, s->size_img);

    wave_tb_data_set_rcv(data,s,sismos);

    wave_inject_sismos(s, u0, s->time_steps, sismos);
    wave_update_fields(s, u0, u1, vel, pml_tmp, pml_tab);
    wave_inject_sismos(s, u1, s->time_steps-1 , sismos);
    wave_tb_backward(ctx,data, timer, u0,u1,vel);

    wave_tb_data_unset_rcv(data);

    wave_tb_timer_info(timer,ctx->nb_stencils_total_bwd,ctx->nb_stencils_main);

    /// save the img/ilm of the shot:
    wave_save_img(s, shot, img_shot, ilm_shot);
    /// release/close the resources related to the current shot.
    shot_release(shot);

    wave_tb_data_free(data,ctx->num_thread_groups);

    MSG("Processing shot %d (END)",sidx);
  }

  wave_tb_free(ctx);
  wave_tb_timer_free(timer);

  free(ctx);   ctx = NULL;
  free(data);  data = NULL;
  free(timer); timer = NULL;
  /// free the simulation buffers.
  DELETE_BUFFER(fwd);
  DELETE_BUFFER(u0);
  DELETE_BUFFER(u1);
  DELETE_BUFFER(img_shot);
  DELETE_BUFFER(ilm_shot);
  DELETE_BUFFER(sismos);
  DELETE_BUFFER(pml_tmp);
}

/*
///
/// Reverse Time Migration on GPU.
///
///
void run_rtm_gpu(sismap_t *s, float* vel, float *source, float *pml_tab) {
  /// u0 on the GPU.
  float* d_u0;
  /// u1 on the GPU.
  float* d_u1;
  /// vel on the GPU.
  float* d_vel;
  /// sismos on the GPU.
  float* d_sismos, *sismos;
  /// the fwd wave-field on the GPU.
  float *d_fwd, *fwd;
  /// the final image on the GPU.
  float *d_img_shot, *d_ilm_shot;
  /// PML GPU arrays.
  float *d_pml_tab, *d_pml_tmp;
  gpu_wave_set(s->device);
  GPU_CHK(cudaMalloc((void**)&d_u0,  s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_u1,  s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_fwd, s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_img_shot, s->size_img*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_ilm_shot, s->size_img*sizeof(float)));
  GPU_CHK(cudaMemset(d_u0,  0, s->size*sizeof(float)));
  GPU_CHK(cudaMemset(d_u1,  0, s->size*sizeof(float)));
  GPU_CHK(cudaMemset(d_fwd, 0, s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_vel, s->size_eff*sizeof(float)));
  GPU_CHK(cudaMemcpy(d_vel, vel,
                     s->dimx*s->dimy*s->dimz*sizeof(float),
                     cudaMemcpyHostToDevice));
  GPU_CHK(cudaMalloc((void**)&d_sismos,
                     s->rcv_len*s->time_steps*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_pml_tab, (s->dimx+2)*(s->dimy+2)*
                                         (s->dimz+2)*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_pml_tmp,
                             s->dimx*s->dimy*s->dimz*sizeof(float)));
  GPU_CHK(cudaMemcpy(d_pml_tab, pml_tab,
                    (s->dimx+2)*(s->dimy+2)*(s->dimz+2)*sizeof(float),
                    cudaMemcpyHostToDevice));
  CREATE_BUFFER(sismos, s->rcv_len*s->time_steps);
  CREATE_BUFFER(fwd, s->size);
  float *tmp;
  CREATE_BUFFER(tmp,  s->size_img);
  gpu_wave_init(s);

  if (s->verbose) gpu_wave_info(s);

  shot_t *shot;

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, false, s->modeling);
    /// load the seismic traces for the shot.
    gpu_wave_read_sismos(s, shot, d_sismos, sismos);
    /// reset the buffers for the shot.
    s->snap_idx = 0;
    GPU_CHK(cudaMemset(d_u0,     0,     s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_u1,     0,     s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_pml_tmp, 0, s->size_eff*sizeof(float)));
    /// forward modeling.
    for(int t = 0; t < s->time_steps; ++t) {
      gpu_wave_update_source(s, shot, d_u0, source[t]);
      gpu_wave_update_fields(s, d_u0, d_u1, d_vel, d_pml_tmp, d_pml_tab);
      gpu_wave_save_snapshot(s, shot, d_u1, fwd, t);
      GPU_WAVE_SWAP_POINTERS(d_u0, d_u1);
    }
    /// reset buffers for the shot (backward).
    GPU_CHK(cudaDeviceSynchronize());
    s->snap_idx = s->snap_idx-1;
    GPU_CHK(cudaMemset(d_u0, 0, s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_u1, 0, s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_pml_tmp, 0, s->size_eff*sizeof(float)));
    GPU_CHK(cudaMemset(d_img_shot, 0, s->size_img*sizeof(float)));
    GPU_CHK(cudaMemset(d_ilm_shot, 0, s->size_img*sizeof(float)));
    /// backward modeling and imaging.
    for(int t = s->time_steps-1; t>=0 ; --t) {
      gpu_wave_read_snapshot(s, shot, d_fwd, fwd, t);
      gpu_wave_inject_sismos(s, d_u0, t, d_sismos);
      gpu_wave_update_fields(s, d_u0, d_u1, d_vel, d_pml_tmp, d_pml_tab);
      gpu_wave_image_condition(s, d_u1, d_fwd, d_img_shot, d_ilm_shot, t);
      WAVE_SWAP_POINTERS(d_u0, d_u1);
    }
    /// save the img and ilm of the shot:
    gpu_wave_save_img(s, shot, d_img_shot, d_ilm_shot, tmp);
    /// release/close the resources related to the current shot.
    shot_release(shot);
  }
  DELETE_BUFFER(fwd);
  DELETE_BUFFER(tmp);
  DELETE_BUFFER(sismos);
  GPU_CHK(cudaFree(d_u0));
  GPU_CHK(cudaFree(d_u1));
  GPU_CHK(cudaFree(d_vel));
  GPU_CHK(cudaFree(d_fwd));
  GPU_CHK(cudaFree(d_img_shot));
  GPU_CHK(cudaFree(d_ilm_shot));
  GPU_CHK(cudaFree(d_sismos));
  GPU_CHK(cudaFree(d_pml_tmp));
  GPU_CHK(cudaFree(d_pml_tab));
  gpu_wave_release(s);
  gpu_wave_unset();
}
*/

/// @brief The main function of the first part of \b simwave
/// @param argc the number of user's options
/// @param argv contains the user options
/// @return 0 on success
///
///
/// User options parser:
/// - \b p is an options @ref parser
/// - it parses the user's command line or from file options
/// - check if saving snapshots is enabled
/// - get the snapshot frequency
/// - check if the GPU results have to be compared to those of the CPU
/// - check if verbose mode is enabled
/// - check if the execution on the CPU is disabled
///
/// GPU environment:
/// - create an GPU resources descriptor (@ref gpu_engine_t)
/// - print informations about the GPU environment
/// - allocate resources before the GPU runs
/// - deallocate resources after the GPU runs
///
/// CPU wave descriptor (@ref wave_t):
/// - print informations about the wave
/// - simulation on the CPU if enabled
///
/// GPU wave descriptor (@ref gpu_wave_t):
/// - simulation on the GPU (default behavior)
/// - check the GPU results if asked by the user
int main(int argc, char* argv[]) {
  /// structure to maintain the user choices.
  sismap_t *s = (sismap_t*)malloc(sizeof(sismap_t));
  /// create a parser.
  parser *p = parser_create("Reverse Time Migration using simwave");
  /// parse command line arguments.
  PARSER_BOOTSTRAP(p);
  parser_parse(p, argc, argv);
  s->verbose    = parser_get_bool(p, "verbose");
  s->cpu        = parser_get_bool(p, "cpu");
  s->time_steps = parser_get_int(p, "iter");
  s->cfl        = parser_get_float(p, "cfl");
  s->fmax       = parser_get_float(p, "fmax");
  s->vel_file   = parser_get_string(p, "in");
  s->vel_dimx   = parser_get_int(p, "n1");
  s->vel_dimy   = parser_get_int(p, "n2");
  s->vel_dimz   = parser_get_int(p, "n3");
  s->dx         = parser_get_int(p, "dx");
  s->dy         = parser_get_int(p, "dy");
  s->dz         = parser_get_int(p, "dz");
  s->dcdp       = parser_get_int(p, "dcdp");
  s->dline      = parser_get_int(p, "dline");
  s->drcv       = parser_get_int(p, "drcv");
  s->dshot      = parser_get_int(p, "dshot");
  s->ddepth     = parser_get_int(p, "ddepth");
  s->device     = parser_get_int(p, "device");
  s->first      = parser_get_int(p, "first");
  s->last       = parser_get_int(p, "last");
  s->src_depth  = parser_get_int(p, "src_depth");
  s->rcv_depth  = parser_get_int(p, "rcv_depth");
  s->modeling   = false;
  s->nb_snap    = parser_get_int(p, "nbsnap");
  s->mode       = parser_get_int(p, "mode");

  /// contains the velocity values of the traversed mediums.
  float* vel;
  /// contains the terms of the source.
  float* source;
  /// contains the PML coefficients.
  float* pml_tab;
  /// get velocity min max from file and setup numerics.
  wave_init_numerics(s);
  /// initialize the velocity and the compute sizes.
  wave_init_dimensions(s);
  wave_init_damp(s);
  /// initialize the geometry.
  wave_init_acquisition(s);
  /// initialize the simulation buffers.
  if (s->cpu) {
    CREATE_BUFFER(vel, s->size_eff);
  } else {
    CREATE_BUFFER_ONLY(vel, s->size_eff);
    array_openmp_inner_init(vel,s);
  }
  CREATE_BUFFER(source, s->time_steps+1);
  CREATE_BUFFER(pml_tab, (s->dimx+2)*(s->dimy+2)*(s->dimz+2));
  /// load/generate the velocity model.
  velocity_load_model(s, vel);
  /// compute PML parameters.
  pml_compute_coefs(s, pml_tab);
  /// generate the Ricker source.
  source_ricker_wavelet(s, source);
  source[s->time_steps] = 0.0f; // an extra time step for girih.
  /// print info if needed.
  if (s->verbose) wave_print(s);
  /// run RTM on CPU or GPU.
  if(s->cpu) {
    run_rtm_tb_cpu(s, vel, source, pml_tab,p);
  } else {
    run_rtm_cpu(s, vel, source, pml_tab);
//    run_rtm_gpu(s, vel, source, pml_tab);
  }
  /// free the simulation buffers.
  DELETE_BUFFER(vel);
  DELETE_BUFFER(source);
  DELETE_BUFFER(pml_tab);
  /// release simwave.
  wave_release(s);
  /// release the simulation structure.
  free(s);
  /// delete the parser.
  parser_delete(p);
  return EXIT_SUCCESS;
}
