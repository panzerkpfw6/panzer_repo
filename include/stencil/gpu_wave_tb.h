#ifndef STENCIL_GPU_WAVE_TB_H
#define STENCIL_GPU_WAVE_TB_H

#include <stencil/sismap.h>
#include <stencil/shot.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct gpu_tb_ctx
{
    int device_id;

    int nx;
    int ny;
    int nz;

    int nt;
    int nrcv;

    size_t nxyz;

    size_t field_bytes;
    size_t source_bytes;
    size_t coef_bytes;
    size_t receiver_index_bytes;
    size_t sismos_bytes;
    size_t allocated_bytes;

    int coef_count;

    float *d_vel;
    float *d_rho;
    float *d_inv_rho;
    float *d_roc2;

    float *d_source;

    float *d_p1;
    float *d_p2;
    float *d_p3;

    float *d_v1;
    float *d_v2;
    float *d_v3;

    float *d_coefx;
    float *d_coefy;
    float *d_coefz;

    unsigned int *d_rcv;
    float *d_sismos;

} gpu_tb_ctx_t;

int gpu_wave_tb_init(
    gpu_tb_ctx_t *ctx,
    const sismap_t *s,
    int coef_count);

int gpu_wave_tb_allocate(
    gpu_tb_ctx_t *ctx);

int gpu_wave_tb_copy_static_data(
    gpu_tb_ctx_t *ctx,
    const float *vel,
    const float *inv_rho,
    const float *source,
    const unsigned int *rcv);

int gpu_wave_tb_zero(
    gpu_tb_ctx_t *ctx);

void gpu_wave_tb_info(
    const gpu_tb_ctx_t *ctx,
    const sismap_t *s);

void gpu_wave_tb_release(
    gpu_tb_ctx_t *ctx);


#ifdef __cplusplus
}
#endif

#endif