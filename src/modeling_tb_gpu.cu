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
#include <stencil/wave_tb.h>
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
    /*
     * Host-side output buffer.
     *
     * The actual wavefields live on the GPU in gctx.
     */
    float *sismos = NULL;

    /*
     * Host TB descriptors.
     *
     * We keep these because they contain:
     *   - diamond geometry
     *   - thread/wavefront parameters
     *   - source metadata
     *   - receiver metadata
     *   - timing information
     */
    tb_t *ctx = NULL;
    tb_data_t *data = NULL;
    tb_timer_t *timer = NULL;

    shot_t *shot = NULL;

    (void)p;

    gpu_tb_ctx_t gctx;
    int status = 0;
    int gpu_selected = 0;

    memset(&gctx, 0, sizeof(gctx));

    MSG("... entering 1st-order TB GPU driver");

    /*
     * ------------------------------------------------------------
     * 1. Initialize normal CPU-side TB metadata.
     * ------------------------------------------------------------
    */
    wtime_init();

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
    
    
    // Copy all static model data to the GPU.
    status = gpu_wave_tb_copy_static_data(
        &ctx,
        s,vel,
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

    cleanup:
    gpu_wave_tb_release(&ctx);

    // Release CUDA context.
    gpu_wave_unset();

    MSG("... leaving 1st-order TB GPU driver");
}
