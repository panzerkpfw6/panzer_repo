#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>

#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/parser.h>
#include <stencil/sismap.h>
#include <stencil/shot.h>
#include <stencil/gpu_wave.h>
#include <stencil/gpu_wave_tb.h>
#include <stencil/modeling_tb_gpu.h>

/**
 * First-order MWD-TB modeling on GPU.
 *
 * Stage 2:
 *   - select and validate the CUDA device;
 *   - initialize the GPU-TB context;
 *   - allocate GPU buffers;
 *   - initialize wavefields and receiver output;
 *   - report allocation information;
 *   - release GPU resources.
 *
 * Static model transfer and numerical propagation are not yet complete.
 */
 
extern "C" void run_modeling_1st_tb_gpu(
    sismap_t *s,
    float *vel,
    float *inv_rho,
    float *source,
    parser *p)
{
    (void)p;

    gpu_tb_ctx_t ctx;
    int status = 0;

    memset(&ctx, 0, sizeof(ctx));

    MSG("... entering 1st-order TB GPU driver");

    /* Select and validate the CUDA device.*/
    gpu_wave_set(s->device);

    /* Initialize GPU context.*/
    status = gpu_wave_tb_init(&ctx,s);
    if (status != 0) {
        MSG("... GPU TB context initialization failed");
        goto cleanup;
    }

    // Allocate GPU memory.
    status = gpu_wave_tb_allocate(&ctx);
    if (status != 0) {
        MSG("... GPU TB memory allocation failed");
        goto cleanup;
    }
    
    gpu_wave_tb_zero(&ctx);

    // Copy all static model data to the GPU.
    status = gpu_wave_tb_copy_static_data(
        &ctx,
        vel,
        inv_rho,
        source,
        s->rcv);

    if (status != 0) {
        MSG("... GPU TB static-data transfer failed");
        goto cleanup;
    }

    // Zero only dynamic wavefields and receiver output.
    status = gpu_wave_tb_zero(&ctx);

    if (status != 0) {
        MSG("... GPU TB wavefield initialization failed");
        goto cleanup;
    }

    // Print GPU information.
    if (s->verbose) {
        gpu_wave_tb_info(&ctx, s);
    }


    gpu_wave_tb_release(&ctx);

    // Release CUDA context.
    gpu_wave_unset();

    MSG("... leaving 1st-order TB GPU driver");
}
