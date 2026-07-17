#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/sismap.h>
#include <stencil/shot.h>
#include <stencil/gpu_wave.h>
#include <stencil/modeling_gpu.h>
/// Modeling on GPU.
extern "C" void run_modeling_1st_sb_gpu(sismap_t *s,float* vel,float* inv_rho,float *source, float *pml_tab) {
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

  /** CUDA profiling events. cudaEventElapsedTime() returns milliseconds.*/
  cudaEvent_t event_start;
  cudaEvent_t event_stop;

  GPU_CHK(cudaEventCreate(&event_start));
  GPU_CHK(cudaEventCreate(&event_stop));

  gpu_wave_set(s->device);

  GPU_CHK(cudaMalloc((void**)&d_u0, s->size*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_u1, s->size*sizeof(float)));
  
  GPU_CHK(cudaMalloc((void**)&d_sismos,s->rcv_len*s->time_steps*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_vel, s->size_eff*sizeof(float)));
  GPU_CHK(cudaMemcpy(d_vel, vel,s->dimx*s->dimy*s->dimz*sizeof(float),cudaMemcpyHostToDevice));
  GPU_CHK(cudaMalloc((void**)&d_pml_tab, (s->dimx+2)*(s->dimy+2)*(s->dimz+2)*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&d_pml_tmp,s->dimx*s->dimy*s->dimz*sizeof(float)));
  GPU_CHK(cudaMemcpy(d_pml_tab, pml_tab,(s->dimx+2)*(s->dimy+2)*(s->dimz+2)*sizeof(float),cudaMemcpyHostToDevice));
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
