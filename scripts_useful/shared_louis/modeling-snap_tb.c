///
/// @copyright Copyright 2017. All rights reserved.
/// This file is part of \b stencil.
///
/// \b stencil is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// \b stencil is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details->
///
/// You should have received a copy of the GNU General Public License
/// along with \b stencil.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author .
/// @file src/main.c
/// @brief Main program of the wave simulator.
///
/// Contains the main program that runs the Seismic Modeling using stencil.
///
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
//#include <cuda_runtime.h>
#include <simwave/parser.h>
#include <simwave/simwave.h>
#include <simwave/wave_tb.h>
#include <simwave/wtime.h>
#include <simwave/macros.h>
/// Modeling on CPU.
///
///


void run_modeling_cpu(sismap_t *s, float* vel,  float *source, float *pml_tab)  {
  /// contains the fields pressure value at time step t.
  float* u0;
  /// contains the fields pressure value at time step t+1.
  float* u1;
  /// seismic traces for a given shot.
  float *sismos;
  /// PML temporary tab.
  float *pml_tmp;
  MSG("... !SB MODE! ...");
  CREATE_BUFFER_ONLY(u0, s->size);
  CREATE_BUFFER_ONLY(u1, s->size);
  array_openmp_init(u0,s);
  array_openmp_init(u1,s);

  CREATE_BUFFER(sismos, s->rcv_len*(s->time_steps+1));
  CREATE_BUFFER(pml_tmp, s->size_eff);
  shot_t *shot;

  double t0,t1,t2,t_prop,t_sismos;
  wtime_init();

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, true, s->modeling);
    /// reset some buffers for the shot.
    NULIFY_BUFFER(u0, s->size);
    NULIFY_BUFFER(u1, s->size);
    NULIFY_BUFFER(pml_tmp, s->size_eff);
    NULIFY_BUFFER(sismos, s->rcv_len*(s->time_steps+1));
    t1 = wtime();
    /// forward modeling.
#if 0
    for(int t = 0; t <= s->time_steps-1; ++t) {
      wave_update_source(s, shot, u0, source[t]);
      wave_update_fields_block_bis(s, u0, u1, vel, pml_tmp, pml_tab);
      #ifdef __DEBUG
      wave_save_fwd_dbg(s, shot, u1, t%s->nb_snap==0);
      #endif // __DEBUG
      wave_extract_sismos(s, u1, t, sismos);
      WAVE_SWAP_POINTERS(u0, u1);
    }
#else
    t1 = wtime();
    t_prop   = 0.0;
    t_sismos = 0.0;

    t0 = wtime();
    wave_extract_sismos(s, u0, 0, sismos);
    t_sismos += wtime() - t0;
    #ifdef __DEBUG
    wave_save_fwd_dbg(s, shot, u0, 0%s->nb_snap==0);
//    wave_save_snapshot(s, shot, u0, 0%s->nb_snap==0);
    #endif // __DEBUG
    for(int t = 0; t <= s->time_steps-1; ++t) {
/*	//CUSTOM SNAPSHOTTING
        int snap_freq = 1000;
        if ((t+1) % snap_freq == 0) {
            MSG("SNAPSHOT CUSTOM");
            FILE *snap_fd;
            char *snap_fd_name = (char*)malloc(20*sizeof(char));

            sprintf(snap_fd_name, "snapshot_%u", t+1);
            snap_fd = fopen(snap_fd_name, "wb+");

            CHK(snap_fd == NULL, "failed to open custom snapshot file");
            CHK(fwrite(u1, sizeof(float), s->size, snap_fd) != s->size,
                "failed to write custom snapshot");
            CHK(fclose(snap_fd) != 0, "failed to close custom snapshot file");

            if (s->verbose) MSG("... saving snapshot n°%u (size %d)", t+1, s->size);

            free(snap_fd_name);
        }
*/
	//@MODIF activer modeling avec MODELING=1
	//char* value = getenv("MODELING");
	//int intValue;
	//if (value != NULL && sscanf(value, "%d", &intValue) == 1 && intValue == 1){
        	wave_update_source(s, shot, u0, source[t]);
	//}

        t0 = wtime();
        wave_update_fields_block_bis(s, u0, u1, vel, pml_tmp, pml_tab);
        t_prop += wtime() - t0;

        #ifdef __DEBUG
        wave_save_fwd_dbg(s, shot, u1, (t+1)%s->nb_snap==0);
        //      wave_save_snapshot(s, shot, u1, (t+1)%s->nb_snap==0);
        #endif // __DEBUG

        t0 = wtime();
        wave_extract_sismos(s, u1, t+1, sismos);
        t_sismos += wtime() - t0;

        WAVE_SWAP_POINTERS(u0, u1);
    }
#endif
    t2 = wtime();
    MSG("forward timer");
    MSG("Total:        %f (s)",t2-t1);
    MSG("PROP:         %f (s)",t_prop);
    MSG("SISMOS:       %f (s)",t_sismos);
    MSG("Speed:        %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t2-t1));
    MSG("PropSpeed:    %f Mstencils/s",1.0*s->time_steps*s->size_eff/1e6/(t_prop) );
/// save the seismic traces for the shot.
    wave_save_sismos(s, shot, sismos);
    /// release/close the resources related to the current shot.
    shot_release(shot);
  }
  /// free the simulation buffers.
  DELETE_BUFFER(u0);
  DELETE_BUFFER(u1);
  DELETE_BUFFER(sismos);
  DELETE_BUFFER(pml_tmp);
}
/*
/// Modeling on GPU.
///
///
void run_modeling_gpu(sismap_t *s, float* vel, float *source, float *pml_tab) {
  /// seismic traces for a given shot.
  float *sismos;
  /// An image for @ref u0 on the GPU.
  float *d_u0;
  /// An image for @ref u1 on the GPU.
  float *d_u1;
  /// An image for @ref vel on the GPU.
  float *d_vel;
  /// An image for @ref sismos on the GPU.
  float *d_sismos;
  /// PML GPU arrays.
  float *d_pml_tab, *d_pml_tmp;
  gpu_wave_set(s->device);

  GPU_CHK(cudaMalloc((void**)&d_u0, s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_u1, s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_sismos,
                     s->rcv_len*s->time_steps*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_vel, s->size_eff*sizeof(float)));
  GPU_CHK(cudaMemcpy(d_vel, vel,
                     s->dimx*s->dimy*s->dimz*sizeof(float),
                     cudaMemcpyHostToDevice));
 	GPU_CHK(cudaMalloc((void**)&d_pml_tab, (s->dimx+2)*(s->dimy+2)*
                                         (s->dimz+2)*sizeof(float)));
 	GPU_CHK(cudaMalloc((void**)&d_pml_tmp,
                             s->dimx*s->dimy*s->dimz*sizeof(float)));
	GPU_CHK(cudaMemcpy(d_pml_tab, pml_tab,
                    (s->dimx+2)*(s->dimy+2)*(s->dimz+2)*sizeof(float),
                    cudaMemcpyHostToDevice));
  CREATE_BUFFER(sismos, s->rcv_len*s->time_steps);

  #ifdef __DEBUG
  float *tmp;
  CREATE_BUFFER(tmp, s->size);
  #endif // __DEBUG

  gpu_wave_init(s);

  if (s->verbose) gpu_wave_info(s);

  shot_t *shot;

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, false, s->modeling);
    /// reset the buffers for the shot.
    GPU_CHK(cudaMemset(d_u0, 0, s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_u1, 0, s->size*sizeof(float)));
    GPU_CHK(cudaMemset(d_sismos, 0, s->rcv_len*s->time_steps*sizeof(float)));
    GPU_CHK(cudaMemset(d_pml_tmp, 0, s->size_eff*sizeof(float)));
    /// forward modeling.
    for(int t = 0; t <= s->time_steps-1; ++t) {
      gpu_wave_update_source(s, shot,  d_u0, source[t]);
      gpu_wave_update_fields(s, d_u0, d_u1, d_vel, d_pml_tmp, d_pml_tab);
      #ifdef __DEBUG
      gpu_wave_save_fwd_dbg(s, shot, d_u1, tmp, t%s->nb_snap==0);
      #endif // __DEBUG
      gpu_wave_extract_sismos(s, d_u1, t, d_sismos);
      GPU_WAVE_SWAP_POINTERS(d_u0, d_u1);
    }
    /// save the seismic traces for the shot.
    gpu_wave_save_sismos(s, shot, d_sismos, sismos);
    /// release/close the resources related to the current shot.
    shot_release(shot);
  }
  /// release buffers.
  GPU_CHK(cudaFree(d_u0));
  GPU_CHK(cudaFree(d_u1));
  GPU_CHK(cudaFree(d_vel));
  GPU_CHK(cudaFree(d_sismos));
  GPU_CHK(cudaFree(d_pml_tmp));
  GPU_CHK(cudaFree(d_pml_tab));
  #ifdef __DEBUG
  DELETE_BUFFER(tmp);
  #endif // __DEBUG
  DELETE_BUFFER(sismos);
  gpu_wave_release(s);
  gpu_wave_unset();
}
*/
/// Modeling on CPU.
///
///
void run_modeling_tb_cpu(sismap_t *s, float* vel,  float *source, parser *p) {
  /// contains the fields pressure value at time step t.
  float* u0;
  /// contains the fields pressure value at time step t+1.
  float* u1;
  /// seismic traces for a given shot.
  float *sismos;
  /// PML temporary tab.
  float *pml_tmp;
  CREATE_BUFFER(u0, s->size);
  CREATE_BUFFER(u1, s->size);
  CREATE_BUFFER(pml_tmp, s->size_eff);
  shot_t *shot;

  wtime_init();

  tb_t * ctx         = (tb_t*)       malloc(sizeof(tb_t));
  tb_data_t * data   = (tb_data_t*)  malloc(sizeof(tb_data_t));
  tb_timer_t * timer = (tb_timer_t*) malloc(sizeof(tb_timer_t));

  wave_tb_init(ctx,s,p);
  source = realloc(source,sizeof(float)*(s->time_steps+1));

  /*
  //@MODIF: NUMA First touch policy, do it if FIRST_TOUCH==1
  char* value = getenv("FIRST_TOUCH");
  int intValue;
  if (value != NULL && sscanf(value, "%d", &intValue) == 1 && intValue == 1){
    int zmin, zmax;
    int ymin, ymax;
    int xmin, xmax;
#pragma omp parallel for collapse(3) private(zmin,zmax,ymin,ymax,xmin,xmax)
    for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
      for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
        for(int xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
          zmax = zmin+BLOCKZ;
          ymax = ymin+BLOCKY;
          xmax = xmin+BLOCKX;
          for(int z = zmin; z < zmax ; z++) {
            for(int y = ymin; y < ymax; y++) {
              for(int x = xmin; x < xmax; x++) {
                u0[x] = 0.;
	        u1[x] = 0.;
	      }
            }
          }
        }
      }
    }
  }
  */

  wave_tb_info(ctx);
  wave_tb_timer_init(timer,ctx->thread_group_size,ctx->num_thread_groups);
  CREATE_BUFFER(sismos, s->rcv_len*(s->time_steps+1));

  MSG("loop over the shots");
  printf("rcv_len %d, time_steps %d\n",s->rcv_len,s->time_steps);

  /// loop over the shots.
  for (int sidx = s->first; sidx <= s->last; sidx++) {
    MSG("Processing shot %d",sidx);
    /// retrieve the shot descriptor.
    shot = s->shots[sidx];
    /// initialize the current shot.
    shot_init(shot, true, s->modeling);
    /// reset some buffers for the shot.
    NULIFY_BUFFER(u0, s->size);
    NULIFY_BUFFER(u1, s->size);
    NULIFY_BUFFER(pml_tmp, s->size_eff);
    NULIFY_BUFFER(sismos, s->rcv_len*(s->time_steps+1)); // add one time step

    // setup tb_data
    wave_tb_data_init(data,ctx,s,ctx->num_thread_groups,shot->id,
                      1ULL*ctx->nnx * ctx->nny * ctx->diam_width);
    data->flag_fwd = 1;
    data->flag_bwd = 0;
    data->fwd = NULL;
    data->ilm = NULL;
    data->img = NULL;

    //  Pavel's modification
//    data->fwd = 1.0;

    /// forward modeling.
    wave_tb_timer_clear(timer);

    wave_tb_data_set_src(data,s,shot->srcidx,source);
    wave_tb_data_set_rcv(data,s,sismos);
    wave_tb_data_info(data);

    wave_tb_forward(ctx,data,timer,u0,u1,vel);

    //@MODIF : trying to snapshot after all computations
    wave_tb_save_lastshot(s, shot, u0, u1, 104);
    //wave_save_snapshot(s, shot, u1, 104);


    wave_tb_data_unset_src(data);
    wave_tb_data_unset_rcv(data);

    wave_tb_timer_info(timer,ctx->nb_stencils_total_fwd,ctx->nb_stencils_main);

    /// save the seismic traces for the shot.
    wave_save_sismos(s, shot, sismos);

    /// release/close the resources related to the current shot.
    shot_release(shot);
    wave_tb_data_free(data,ctx->num_thread_groups);
  }

  wave_tb_free(ctx);
  wave_tb_timer_free(timer);
  /// free the simulation buffers.

  free(ctx);   ctx = NULL;
  free(data);  data = NULL;
  free(timer); timer = NULL;

  DELETE_BUFFER(u0);
  DELETE_BUFFER(u1);
  DELETE_BUFFER(sismos);
  DELETE_BUFFER(pml_tmp);
}

/// @brief The main function of the first part of \b stencil
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
  parser *p = parser_create("Seismic Modelling using stencil");
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
  s->ddepth     = parser_get_int(p, "ddepth");
  s->drcv       = parser_get_int(p, "drcv");
  s->dshot      = parser_get_int(p, "dshot");
  s->device     = parser_get_int(p, "device");
  s->first      = parser_get_int(p, "first");
  s->last       = parser_get_int(p, "last");
  s->src_depth  = parser_get_int(p, "src_depth");
  s->rcv_depth  = parser_get_int(p, "rcv_depth");
  s->modeling   = true;
  /// Nyquist sampling.
  s->nb_snap    = parser_get_int(p, "nbsnap");
    s->mode       = parser_get_int(p, "mode");

  /// contains the velocity values of the traversed mediums.
  float* vel;
  /// contains the terms of the source.
  float* source;
  /// contains the PML coefficients.
  float* pml_tab;
  /// get velocity min max from file and setup numerics.
  MSG("... cpu=: %d\n", s->cpu);
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
  /// generate the ricker source.
  source_ricker_wavelet(s, source);
  source[s->time_steps] = 0.0f; // an extra time step for girih.
  /// print info if needed.
  if (s->verbose) wave_print(s);
  /// run RTM on CPU or GPU.
  if(s->cpu) {
//    run_modeling_cpu(s, vel, source, pml_tab);
    run_modeling_tb_cpu(s, vel, source, p);
  } else {
    run_modeling_cpu(s, vel, source, pml_tab);

//    run_modeling_gpu(s, vel, source, pml_tab);
  }
  /// free the simulation buffers.
  DELETE_BUFFER(vel);
  DELETE_BUFFER(source);
  DELETE_BUFFER(pml_tab);
  /// release stencil.
  wave_release(s);
  /// release the sismap structure.
  free(s);
  /// delete the parser.
  parser_delete(p);

  MSG("END\n");
  return EXIT_SUCCESS;
}
