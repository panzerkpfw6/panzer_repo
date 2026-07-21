#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>
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
 * Stage 1 purpose:
 *
 *   1. Validate CPU-to-CUDA function linkage.
 *   2. Validate GPU runtime selection.
 *   3. Validate that the new CUDA files are included in libstencil_cuda.
 *   4. Establish the initialization/release lifecycle.
 *
 * No wavefield allocation or numerical propagation is performed yet.
 */
extern "C" void run_modeling_1st_tb_gpu(
    sismap_t *s,
    float *vel,
    float *inv_rho,
    float *source,
    parser *p)
{
    CHK(
        s == NULL,
        "run_modeling_1st_tb_gpu received a NULL sismap pointer");

    CHK(
        vel == NULL,
        "run_modeling_1st_tb_gpu received a NULL velocity pointer");

    CHK(
        inv_rho == NULL,
        "run_modeling_1st_tb_gpu received a NULL inverse-density pointer");

    CHK(
        source == NULL,
        "run_modeling_1st_tb_gpu received a NULL source pointer");

    CHK(
        p == NULL,
        "run_modeling_1st_tb_gpu received a NULL parser pointer");

    MSG("... entering 1st-order TB GPU driver");

    /*
     * Select the GPU requested by the existing sismap configuration.
     *
     * gpu_wave_set() already checks:
     *
     *   - whether CUDA devices exist;
     *   - whether s->device is a valid device index.
     */
    gpu_wave_set(s->device);

    /*
     * Stage 1 initialization.
     */
    gpu_wave_tb_init(s);

    if (s->verbose) {
        gpu_wave_tb_info(s);

        MSG("... GPU TB model dimensions: %u x %u x %u",
            s->dimx,
            s->dimy,
            s->dimz);

        MSG("... GPU TB time steps: %u",
            s->time_steps);

        MSG("... GPU TB number of shots: first=%d last=%d",
            s->first,
            s->last);

        MSG("... GPU TB numerical propagation is not implemented yet");
    }

    /*
     * Stage 1 cleanup.
     */
    gpu_wave_tb_release();

    /*
     * Keep this call last.
     *
     * The existing implementation of gpu_wave_unset() calls
     * cudaDeviceReset(), so no CUDA work should follow it.
     */
    gpu_wave_unset();

    MSG("... leaving 1st-order TB GPU driver");
}