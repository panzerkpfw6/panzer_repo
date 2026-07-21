#ifndef STENCIL_MODELING_TB_GPU_H
#define STENCIL_MODELING_TB_GPU_H

#include <stencil/sismap.h>

#ifdef __cplusplus
extern "C" {
#endif
/**
 * Run first-order MWD-TB acoustic modeling on a CUDA GPU.
 *
 * @param s        Seismic-modeling configuration.
 * @param vel      Host velocity/roc2 model.
 * @param inv_rho  Host inverse-density model.
 * @param source   Host source wavelet.
 * @param p        Command-line/parser configuration containing TB parameters.
 */
void run_modeling_1st_tb_gpu(
    sismap_t *s,
    float *vel,
    float *inv_rho,
    float *source,
    parser *p);

#ifdef __cplusplus
}
#endif

#endif