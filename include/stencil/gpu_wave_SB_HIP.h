
#ifndef STENCIL_GPU_WAVE_SB_HIP_H
#define STENCIL_GPU_WAVE_SB_HIP_H

#include <stencil/sismap.h>
#include <stencil/shot.h>

#ifdef __cplusplus
extern "C" {
#endif

/* HIP execution test */
int hip_test(void);

/* Device management */
void gpu_wave_set(int device);
void gpu_wave_unset(void);
void gpu_wave_init(sismap_t *s);
void gpu_wave_release(sismap_t *s);
void gpu_wave_info(sismap_t *s);

/* Wave propagation */
void gpu_wave_update_fields(sismap_t *s, float *d_u0, float *d_u1,
                            float *d_roc2, float *phi, float *eta);
void gpu_wave_update_source(sismap_t *s, shot_t *shot, float *d_u0,
                            float sterm);

/* Seismic traces */
void gpu_wave_extract_sismos(sismap_t *s, float *d_u0, unsigned int t,
                             float *d_sismos);
void gpu_wave_inject_sismos(sismap_t *s, float *d_u0, unsigned int t,
                            float *d_sismos);

/* Imaging */
void gpu_wave_image_condition(sismap_t *s, float *d_u1, float *d_fwd,
                              float *d_img, float *d_ilm, unsigned int t);
void gpu_wave_image_gather(sismap_t *s, float *d_img, float *d_ilm,
                           float *gimg);

/* Snapshots */
void gpu_wave_save_snapshot(sismap_t *s, shot_t *shot, float *d_u1,
                            float *u1, unsigned int t);
void gpu_wave_read_snapshot(sismap_t *s, shot_t *shot, float *d_fwd,
                            float *fwd, unsigned int t);

/* Numerical validation */
unsigned int gpu_wave_check_fields(sismap_t *s, float *d_u0, float *u1,
                                   float *u0, float e);

/* Seismogram I/O */
void gpu_wave_save_sismos(sismap_t *s, shot_t *shot, float *d_sismos,
                          float *sismos);
void gpu_wave_read_sismos(sismap_t *s, shot_t *shot, float *d_sismos,
                          float *sismos);

/* Debug output */
void gpu_wave_save_fwd_dbg(sismap_t *s, shot_t *shot, float *d_u1,
                           float *u1, unsigned int t);
void gpu_wave_save_bwd_dbg(sismap_t *s, shot_t *shot, float *d_u1,
                           float *u1, unsigned int t);

/* Image output */
void gpu_wave_save_img(sismap_t *s, shot_t *shot, float *d_img,
                       float *d_ilm, float *tmp);

#ifdef __cplusplus
}
#endif

#endif /* STENCIL_GPU_WAVE_SB_HIP_H */
