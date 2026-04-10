/// Contains the main program that runs the Reverse Time Migration (RTM) using ///
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/sysinfo.h>
#include <errno.h>
#include <stencil/parser.h>
#include <stencil/stencil.h>
#include <stencil/wave_tb.h>
#include <stencil/wave.h>
#include <stencil/wtime.h>
#include <time.h>
#include <stencil/macros.h>


//#include <cuda_runtime.h>

void wave_tb_data_close_open(tb_data_t * data,const int nb_thread_groups,const int shotid) {
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
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;

	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const int sx = s->sx;
	const int sy = s->sy;
	const int sz = s->sz;
	const int nnx = dimx + 2 * sx;
	const int nny = dimy + 2 * sy;
	const int nnz = dimz + 2 * sz;
	const long int nnxy=(long int)nnx * nny; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)nny * nnz;
	float * restrict ux;
	float * restrict vx;
	/////////////
#pragma omp parallel for collapse(3) schedule(dynamic) private(xmin, xmax, zmin, zmax, ymin, ymax, ux,vx)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);
                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
                        ux = &(  u[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx = &(  v[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        #pragma omp simd
                        for (int z = zmin; z < zmax; z++) {
                        	ux[z] = vx[z];
                        }
                    }
                }
            }
        }
    }
    /////////////
}

///
/// Reverse Time Migration on CPU.
///
///
void run_rtm_cpu(sismap_t *s,float* vel,float *source,float *pml_tab) {
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
//  printf("debug point 3\n");

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

//    printf("debug point 1\n");
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
////      wave_extract_sismos(s, u1, t, sismos);
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
    MSG("Speed:        %f GStencils/s",1.0*s->time_steps*s->size_eff/1e9/(t2-t1));
    MSG("PropSpeed:    %f GStencils/s",1.0*s->time_steps*s->size_eff/1e9/(t_prop));


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

//      t0 = wtime();
//      wave_update_fields_block_bis(s, u0, u1, vel, pml_tmp, pml_tab);
//      t_prop += wtime() - t0;

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
    MSG("Speed:        %f GStencils/s",1.0*s->time_steps*s->size_eff/1e9/(t2-t1));
    MSG("PropSpeed:    %f GStencils/s",1.0*s->time_steps*s->size_eff/1e9/(t_prop));

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
/// 1st order rtm
void run_rtm_1st_cpu(sismap_t *s, float* vel,float *inv_rho,float *source, float *pml_tab) {
    /// contains the fields pressure value at time step t.
    float* u0;
    /// contains the fields (particle velocity across x direction) value at time step t.
    float* vx;
    /// contains the fields (particle velocity across y direction) value at time step t.
    float* vy;
    /// contains the fields (particle velocity across z direction) value at time step t.
    float* vz;
    /// forward wave-field for RTM.
    float* fwd;
    float** fwd_all;
    /// seismic traces for a given shot.
    float *sismos;
    /// image and illumination of each shot.
    float *ilm_shot, *img_shot;
    /// PML tmp tab.
    /// wave-field arrays.
    CREATE_BUFFER_ONLY(u0, s->size);
    array_openmp_init(u0,s);
    CREATE_BUFFER_ONLY(vx, s->size);
    array_openmp_init(vx,s);
    CREATE_BUFFER_ONLY(vy, s->size);
    array_openmp_init(vy,s);
    CREATE_BUFFER_ONLY(vz,s->size);
    array_openmp_init(vz,s);
    CREATE_BUFFER_ONLY(fwd, s->size);
    array_openmp_init(fwd,s);

    CREATE_BUFFER(sismos, s->rcv_len*(s->time_steps+1));
    CREATE_BUFFER_ONLY(img_shot, s->size_img);
    array_openmp_inner_init(img_shot,s);

    CREATE_BUFFER_ONLY(ilm_shot, s->size_img);
    array_openmp_inner_init(ilm_shot,s);
    shot_t *shot;

    double t0,t1,t2,t_snap,t_prop,t_sismos,t_image;
    wtime_init();

    int nb_fwd_snap = (s->time_steps+1 + s->nb_snap - 1) / s->nb_snap;
    printf("s->time_steps %d, nb_snap %d\n",s->time_steps, s->nb_snap);
    printf("Nb of snapshots in FWD: %d\n", nb_fwd_snap);
//  printf("debug point 3\n");

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
        NULIFY_BUFFER(u0,s->size);
        NULIFY_BUFFER(vx,s->size);
        NULIFY_BUFFER(vy,s->size);
        NULIFY_BUFFER(vz,s->size);

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
            wave_update_fields_block_1st(s,u0,vx,vy,vz,vel,inv_rho);
            t_prop += wtime() - t0;

            t0 = wtime();
            if (s->mode == 2) {
                if ((t+1)%s->nb_snap == 0) {
                    memcpy_omp(fwd_all[s->snap_idx],u0,s);
                    if (s->verbose) MSG("... saving snapshot %d t=%u", s->snap_idx, t+1);
                    s->snap_idx ++;
                }
            } else {
                wave_save_snapshot(s,shot,u0,t+1);
            }
            t_snap += wtime() - t0;

            t0 = wtime();
            wave_update_source(s, shot, u0, source[t]);
            t_sismos += wtime() - t0;

//            WAVE_SWAP_POINTERS(u0, u1);
        }
        t2 = wtime();

        MSG("forward timer (%d)",sidx);
        MSG("Total:        %f (s)",t2-t1);
        MSG("SNAP:         %f (s)",t_snap);
        MSG("PROP:         %f (s)",t_prop);
        MSG("SISMOS:       %f (s)",t_sismos);
        MSG("Speed:        %f GStencils/s",2.0*s->time_steps*s->size_eff/1e9/(t2-t1));
        MSG("PropSpeed:    %f GStencils/s",2.0*s->time_steps*s->size_eff/1e9/(t_prop));


        /// reset buffers for the shot (backward).
        s->snap_idx = s->snap_idx-1;
        NULIFY_BUFFER(u0,s->size);
        NULIFY_BUFFER(vx,s->size);
        NULIFY_BUFFER(vy,s->size);
        NULIFY_BUFFER(vz,s->size);
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
            t_snap += wtime()-t0;

            t0 = wtime();
            if ((t+1)%s->nb_snap == 0) {
            	wave_image_condition_block_xyz(s,u0,fwd,img_shot,ilm_shot,t+1);
            }
            t_image += wtime() - t0;

            t0 = wtime();
            wave_update_fields_block_1st(s,u0,vx,vy,vz,vel,inv_rho);
            t_prop += wtime()-t0;

            t0 = wtime();
            wave_inject_sismos(s,u0,t,sismos);
            t_sismos += wtime() - t0;

//            WAVE_SWAP_POINTERS(u0, u1);
        }


        t2 = wtime();

        MSG("backward timer (%d)",sidx);
        MSG("Total:        %f (s)",t2-t1);
        MSG("SNAP:         %f (s)",t_snap);
        MSG("PROP:         %f (s)",t_prop);
        MSG("SISMOS:       %f (s)",t_sismos);
        MSG("IMAGE COND:   %f (s)",t_image);
        MSG("Speed:        %f GStencils/s",2.0*s->time_steps*s->size_eff/1e9/(t2-t1));
        MSG("PropSpeed:    %f GStencils/s",2.0*s->time_steps*s->size_eff/1e9/(t_prop));

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
    DELETE_BUFFER(vx);
    DELETE_BUFFER(vy);
    DELETE_BUFFER(vz);
    DELETE_BUFFER(fwd);
    DELETE_BUFFER(img_shot);
    DELETE_BUFFER(ilm_shot);
    DELETE_BUFFER(sismos);
}

/// Reverse Time Migration on CPU.///

void run_rtm_1st_tb_cpu(sismap_t *s,float *vel,float *inv_rho, float *source, float *pml_tab,parser *p) {
    /// contains the fields pressure value at time step t.
    float* u0;
    /// contains the fields (particle velocity across x direction) value at time step t.
    float* vx;
    /// contains the fields (particle velocity across y direction) value at time step t.
    float* vy;
    /// contains the fields (particle velocity across z direction) value at time step t.
    float* vz;

    /// forward wave-field for RTM.
    float *fwd, *fwd_io;
    /// seismic traces for a given shot.
    float *sismos;
    /// image and illumination of each shot.
    float *ilm_shot, *img_shot;
    /// PML tmp tab.
    float *pml_tmp;
    /// wave-field arrays.
    printf("allocation\n");
    CREATE_BUFFER(u0, s->size);
    CREATE_BUFFER(vx, s->size);
    CREATE_BUFFER(vy, s->size);
    CREATE_BUFFER(vz, s->size);
    CREATE_BUFFER(img_shot, s->size_img);
    CREATE_BUFFER(ilm_shot, s->size_img);
    CREATE_BUFFER(pml_tmp, s->size_eff);
    shot_t *shot;
    wtime_init();
    printf("all\n");

    tb_t *ctx = (tb_t *) malloc(sizeof(tb_t));
    tb_data_t *data = (tb_data_t *) malloc(sizeof(tb_data_t));
    tb_timer_t *timer = (tb_timer_t *) malloc(sizeof(tb_timer_t));

    // setup tb_data
    printf("tb_init\n");
    wave_tb_init(ctx,s,p);
	Parameters *P = (Parameters*) calloc(1, sizeof(Parameters));
	wave_tb_init_p(ctx,s,P);
//    source = realloc(source,sizeof(float)*(s->time_steps+1));

    printf("tb_info\n");
    wave_tb_info(ctx);

    printf("tb_timer\n");
    wave_tb_timer_init(timer, ctx->thread_group_size, ctx->num_thread_groups);

//    printf("Nb of snapshots in FWD: %d\n", ((ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps));

//		int max_ifwd = (2 * P->nt + ctx->fwd_steps - 1) / ctx->fwd_steps;
//	int max_ifwd = (2 * ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps  +1;
	int max_ifwd = (2 * ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps  +1+1;
	printf("Nb of snapshots in FWD: %d\n", max_ifwd );
	size_t fwd_size = s->size * max_ifwd;

    if (ctx->mode == 1 || ctx->mode == 2) {
		ctx->fwd_size=fwd_size;
		CREATE_BUFFER(fwd,fwd_size);
		MSG("Nb of snapshots max_ifwd=%d,fwd_size=%llu",max_ifwd,fwd_size);
////		exit(1);
		///////////////////////
//        CREATE_BUFFER(fwd, 1ULL * s->size * ((ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps));
//        MSG("s->size=%d",s->size);
//        MSG("( ( (ctx->t_len + ctx->fwd_steps - 1)/ctx->fwd_steps)=%d",((ctx->t_len + ctx->fwd_steps-1)/ctx->fwd_steps) );
//		MSG("fwd_size=%llu",1ULL * s->size * ((ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps) );
//        exit(1);
    } else {
        printf("nnx %d, nnz %d, diam_width %d ntg %d\n", ctx->nnx, ctx->nnz, ctx->diam_width, ctx->num_thread_groups);
        CREATE_BUFFER(fwd, 1ULL * ctx->nnx * ctx->nnz * ctx->diam_width * ctx->num_thread_groups);
    }

    printf("rcv_len %d, time_steps %d\n", s->rcv_len, s->time_steps);
    CREATE_BUFFER(sismos, s->rcv_len * (s->time_steps + 1));
//    printf("s->rcv_len*(s->time_steps+1)=%d",s->rcv_len*(s->time_steps+1));
//    exit(1);

    MSG("loop over the shots between %d and %d", s->first, s->last);

    /// loop over the shots.
    for (int sidx = s->first; sidx <= s->last; sidx++) {
        MSG("Processing shot %d (BEGIN)",sidx);
        /// retrieve the shot descriptor.
        shot = s->shots[sidx];
        /// initialize the current shot.
        shot_init(shot, true, s->modeling);
        /// load the seismic traces for the shot.
        wave_read_sismos(s, shot, sismos);

        wave_min_max("sismos", sismos, s->rcv_len * (s->time_steps + 1));

        /// reset buffers for the shot (forward).
        NULIFY_BUFFER(u0, s->size);
        NULIFY_BUFFER(vx, s->size);
        NULIFY_BUFFER(vy, s->size);
        NULIFY_BUFFER(vz, s->size);
        if (ctx->mode == 1 || ctx->mode == 2) {
//			NULIFY_BUFFER(fwd, s->size * ((ctx->t_len + ctx->fwd_steps - 1) / ctx->fwd_steps)); //orig
			NULIFY_BUFFER(fwd, fwd_size);
            fwd_io = NULL;
        } else {
            NULIFY_BUFFER(fwd, 1ULL * ctx->nnx * ctx->nnz * ctx->diam_width * ctx->num_thread_groups);
        }
        NULIFY_BUFFER(pml_tmp, s->size_eff);

        wave_tb_data_init(data, ctx, s, ctx->num_thread_groups, shot->id,
                          1ULL * ctx->nnx * ctx->nnz * ctx->diam_width);

        data->fwd = fwd;
        data->img = img_shot;
        data->ilm = ilm_shot;

        /////////
        // FWD //
        /////////
        MSG("FWD STEP %d", sidx);
        wave_tb_timer_clear(timer);
        data->flag_fwd = 1;
        data->flag_bwd = 0;
        data->rec_sismos=s->rec_sismos;

        Parameters *P = (Parameters*) calloc(1, sizeof(Parameters));

        /// forward modeling.
        wave_tb_data_set_src(data, s, shot->srcidx, source);
        wave_tb_data_set_rcv(data,s,sismos);
        wave_tb_forward_1st(ctx,data,P,timer,u0,vx,vy,vz,vel,inv_rho);
//		exit(1);

        wave_tb_data_unset_src(data);

        MSG("before wave_tb_timer_info");
        MSG("hallo!");
        MSG("ctx->nb_stencils_total_fwd=%llu ",ctx->nb_stencils_total_fwd);
        MSG("ctx->nb_stencils_main=%llu",ctx->nb_stencils_main);
        wave_tb_timer_info(timer, ctx->nb_stencils_total_fwd, ctx->nb_stencils_main);

        wave_tb_data_close_open(data,ctx->num_thread_groups,shot->id);

        /////////
        // BWD //
        /////////
        MSG("BWD STEP %d", sidx);
        wave_tb_timer_clear(timer);
        data->flag_fwd = 0;
        data->flag_bwd = 1;

        /// reset buffers for the shot (backward).
        NULIFY_BUFFER(u0, s->size);
        NULIFY_BUFFER(vx, s->size);
        NULIFY_BUFFER(vy, s->size);
        NULIFY_BUFFER(vz, s->size);
        NULIFY_BUFFER(pml_tmp, s->size_eff);
        NULIFY_BUFFER(ilm_shot, s->size_img);
        NULIFY_BUFFER(img_shot, s->size_img);

        wave_tb_data_set_rcv(data, s, sismos);

		MSG("ctx->nb_stencils_main=%f",ctx->nb_stencils_main);
        wave_inject_sismos(s, u0, P->nt, sismos);
//        wave_update_fields(s, u0, u1, vel, pml_tmp, pml_tab);
//        wave_update_fields_1st(s,u0,vx,vy,vz, vel, pml_tmp, pml_tab);
		wave_update_fields_block_1st(s,u0,vx,vy,vz,vel,inv_rho);
        wave_tb_backward_1st(ctx,data,P,timer, u0,vx,vy,vz,vel,inv_rho);

        wave_tb_data_unset_rcv(data);

        wave_tb_timer_info(timer, ctx->nb_stencils_total_bwd, ctx->nb_stencils_main);

        /// save the img/ilm of the shot:
        wave_save_img(s, shot, img_shot, ilm_shot);
        /// release/close the resources related to the current shot.
        shot_release(shot);

        wave_tb_data_free(data, ctx->num_thread_groups);
        free(P);

        MSG("Processing shot %d (END)", sidx);
    }

    wave_tb_free(ctx);
    wave_tb_timer_free(timer);

    free(ctx);
    ctx = NULL;
    free(data);
    data = NULL;
    free(timer);
    timer = NULL;

    /// free the simulation buffers.
    DELETE_BUFFER(fwd);
    DELETE_BUFFER(u0);
    DELETE_BUFFER(vx);
    DELETE_BUFFER(vy);
    DELETE_BUFFER(vz);
    DELETE_BUFFER(img_shot);
    DELETE_BUFFER(ilm_shot);
    DELETE_BUFFER(sismos);
    DELETE_BUFFER(pml_tmp);
}

int main(int argc, char *argv[]) {
//	MSG("Hi");
	printf("Hi\n");
	MSG("Before running RTM");
//	time_t rawtime;
//	struct tm *timeinfo;
//	char buffer[80];
//	// Start of program
//	time(&rawtime);
//	timeinfo = localtime(&rawtime);
//	strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
//	printf("Program started at: %s\n",buffer);

    /// structure to maintain the user choices.
    sismap_t *s = (sismap_t *) malloc(sizeof(sismap_t));
    /// create a parser.
    parser *p = parser_create("Reverse Time Migration using stencil");
    /// parse command line arguments.
    PARSER_BOOTSTRAP(p);
    parser_parse(p, argc, argv);
    s->verbose = parser_get_bool(p, "verbose");
    s->cpu = parser_get_bool(p, "cpu");
    s->time_steps = parser_get_int(p, "iter");
    s->dt = parser_get_float(p, "dt");
    s->cfl = parser_get_float(p, "cfl");
    s->fmax = parser_get_float(p, "fmax");
    s->vel_file = parser_get_string(p, "in");
    s->vel_dimx = parser_get_int(p, "n1");
    s->vel_dimy = parser_get_int(p, "n2");
    s->vel_dimz = parser_get_int(p, "n3");
    s->dx = parser_get_float(p, "dx");
    s->dy = parser_get_float(p, "dy");
    s->dz = parser_get_float(p, "dz");
    s->dcdp = parser_get_int(p, "dcdp");
    s->dline = parser_get_int(p,"dline");
    s->drcv = parser_get_int(p,"drcv"); // code works only with drcv=1. due to sismos recording specifique.
    s->dshot = parser_get_int(p, "dshot");
    s->ddepth = parser_get_int(p, "ddepth");
    s->device = parser_get_int(p, "device");
    s->first = parser_get_int(p, "first");
    s->last = parser_get_int(p, "last");
    s->src_depth = parser_get_int(p, "src_depth");
    s->rcv_depth = parser_get_int(p, "rcv_depth");
    s->modeling = false;
    s->nb_snap = parser_get_int(p, "nbsnap");
    s->mode = parser_get_int(p, "mode");
    s->order = parser_get_int(p, "order");
    s->rec_sismos= parser_get_int(p, "rec_sismos");

    // Read cache blocking parameters for SB method //
	s->blockx=parser_get_int(p,"cbx");
	s->blocky=parser_get_int(p,"cby");
	s->blockz=parser_get_int(p,"cbz");

    int ncpus=get_nprocs();
    printf("ncpus : %d\n",ncpus);
    printf("# THREADS : %d\n",omp_get_max_threads());

    /// contains the velocity values of the traversed mediums.
    float* vel;
    /// contains the density values of the traversed mediums.
	float* rho;
	/// contains the inverse density values of the traversed mediums.
	float* inv_rho;
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

    printf("cpu=%d, time_steps=%d, size_eff=%ld, dimx=%d, dimy=%d, dimz=%d\n",
               s->cpu, s->time_steps, s->size_eff, s->dimx, s->dimy, s->dimz);

    /// initialize the simulation buffers.
    if (s->cpu) {
        CREATE_BUFFER(vel, s->size_eff);
        CREATE_BUFFER(rho,s->size_eff);
        CREATE_BUFFER(inv_rho,s->size_eff);
    } else {
		MSG("... !SB MODE! ...");
		MSG("BLOCKX=%d, BLOCKY=%d, BLOCKZ=%d\n",s->blockx,s->blocky,s->blockz);
		///////////////////////////////////////
        CREATE_BUFFER_ONLY(vel, s->size_eff);
        array_openmp_inner_init(vel, s);
        CREATE_BUFFER_ONLY(rho, s->size_eff);
		array_openmp_inner_init(rho, s);
		CREATE_BUFFER_ONLY(inv_rho,s->size_eff);
		array_openmp_inner_init(inv_rho, s);
    }
    printf("Before CREATE_BUFFER(source \n");
    CREATE_BUFFER(source, s->time_steps + 1);
    CREATE_BUFFER(pml_tab, (s->dimx + 2) * (s->dimy + 2) * (s->dimz + 2));

    /// load/generate the velocity model.
    velocity_const_model2(s,vel);
//    velocity_2layer_model(s,vel);
//    velocity_load_salt3d(s,vel);

    /// load/generate the density model.
	density_const_model(s,rho,inv_rho);

	/// dump velocity and density grids.
	dump_vel(s,vel,rho);

	/// generate coefficient matrix.
	fill_coef_matrix(s,vel,rho);
	/// dump coefficient matrix.
	dump_coef(s,vel);

    /// generate the Ricker source.
    if (s->order==1) {
        MSG("source 1st order");
//        source_ricker_wavelet(s, source);
        source_ricker_wavelet_1st(s, source);
//        source_ricker_wavelet_2nd(s, source);
    } else {
        MSG("source 2nd order");
        source_ricker_wavelet_2nd(s, source);
    }
    source[s->time_steps] = 0.0f; // an extra time step for girih.
    /// print info if needed.
    if (s->verbose) wave_print(s);

    /// run RTM on CPU or GPU.
    if (s->cpu) {
        if (s->order==1) {
            MSG("run 1st order TB RTM");
            run_rtm_1st_tb_cpu(s,vel,inv_rho,source, pml_tab, p);
        }
    } else {
        if (s->order==1) {
            MSG("run 1st order SB RTM");
            run_rtm_1st_cpu(s,vel,inv_rho,source, pml_tab);
        }
//    run_rtm_gpu(s, vel, source, pml_tab);
    }
    /// free the simulation buffers.
    DELETE_BUFFER(vel);
    DELETE_BUFFER(rho);
    DELETE_BUFFER(inv_rho);
    DELETE_BUFFER(source);
    DELETE_BUFFER(pml_tab);
    /// release stencil.
    wave_release(s);
    /// release the simulation structure.
    free(s);
    /// delete the parser.
    parser_delete(p);
    MSG("RTM HALAS");
    return EXIT_SUCCESS;
}
