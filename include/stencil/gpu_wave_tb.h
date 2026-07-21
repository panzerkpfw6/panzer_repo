#ifndef STENCIL_GPU_WAVE_TB_H
#define STENCIL_GPU_WAVE_TB_H

#include <stencil/sismap.h>
#include <stencil/shot.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * Initialize the GPU-TB numerical subsystem.
     *
     * Stage 1 only validates that the new CUDA translation unit is compiled
     * and linked correctly.
     */
    void gpu_wave_tb_init(sismap_t *s);

    /**
     * Print information about the GPU-TB implementation.
     */
    void gpu_wave_tb_info(const sismap_t *s);

    /**
     * Release resources owned by the GPU-TB numerical subsystem.
     *
     * Stage 1 does not yet allocate resources, but the function is provided
     * now so the final driver lifecycle is already defined.
     */
    void gpu_wave_tb_release(void);

#ifdef __cplusplus
}
#endif

#endif /* STENCIL_GPU_WAVE_TB_H */