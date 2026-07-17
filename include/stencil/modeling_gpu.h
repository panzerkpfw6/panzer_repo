#ifndef STENCIL_MODELING_GPU_H
#define STENCIL_MODELING_GPU_H

#include <stencil/sismap.h>

#ifdef __cplusplus
extern "C" {
#endif

void run_modeling_1st_sb_gpu(
    sismap_t *s,
    float *vel,
    float* inv_rho,
    float *source,
    float *pml_tab
);

#ifdef __cplusplus
}
#endif

#endif