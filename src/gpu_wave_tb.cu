#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>
#include <cuda_runtime.h>

#include <stencil/macros.h>
#include <stencil/sismap.h>
#include <stencil/gpu_wave_tb.h>

/*
 * Stage 1 state.
 *
 * This is deliberately minimal. Later this file will own:
 *
 *   - GPU-TB context initialization;
 *   - CUDA propagation kernels;
 *   - diamond geometry;
 *   - wavefront scheduling;
 *   - source injection;
 *   - receiver extraction.
 */
static int gpu_tb_initialized = 0;

extern "C" void gpu_wave_tb_init(sismap_t *s)
{
    CHK(
        s == NULL,
        "gpu_wave_tb_init received a NULL sismap pointer");

    gpu_tb_initialized = 1;

    if (s->verbose) {
        MSG("... GPU TB numerical subsystem initialized");
    }
}

extern "C" void gpu_wave_tb_info(const sismap_t *s)
{
    CHK(
        s == NULL,
        "gpu_wave_tb_info received a NULL sismap pointer");

    cudaDeviceProp props;

    GPU_CHK(cudaGetDeviceProperties(
        &props,
        s->device));

    MSG("... GPU TB implementation: Stage 1 skeleton");
    MSG("... GPU TB device: %s", props.name);
    MSG("... GPU TB device index: %d", s->device);
    MSG("... GPU TB compute capability: %d.%d",
        props.major,
        props.minor);
}

extern "C" void gpu_wave_tb_release(void)
{
    if (!gpu_tb_initialized) {
        return;
    }

    /*
     * No GPU-TB allocations exist in Stage 1.
     *
     * Later stages will release the GPU-TB context and its device
     * buffers here.
     */
    gpu_tb_initialized = 0;
}