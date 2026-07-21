#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>

#include <stencil/macros.h>
#include <stencil/sismap.h>
#include <stencil/gpu_wave_tb.h>


/*
 * Free a CUDA pointer and reset it to NULL.
 *
 * cudaFree(NULL) is legal, but checking the pointer makes cleanup
 * easier to inspect and avoids unnecessary runtime calls.
 */
#define GPU_TB_FREE(pointer)            \
    do                                  \
    {                                   \
        if ((pointer) != NULL)          \
        {                               \
            GPU_CHK(cudaFree(pointer)); \
            (pointer) = NULL;           \
        }                               \
    } while (0)


extern "C" int gpu_wave_tb_init(
    gpu_tb_ctx_t *ctx,
    const sismap_t *s,
    int coef_count)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_init received a NULL context pointer");

    CHK(
        s == NULL,
        "gpu_wave_tb_init received a NULL sismap pointer");

    CHK(
        coef_count <= 0,
        "gpu_wave_tb_init received an invalid coefficient count");

    /*
     * Make every device pointer NULL before any allocation occurs.
     * This guarantees that gpu_wave_tb_release() is safe after a
     * partially successful allocation.
     */
    memset(ctx, 0, sizeof(*ctx));

    ctx->device_id = s->device;

    /*
     * Check these names against sismap.h.
     *
     * Based on your current code and output, the computational dimensions
     * are expected to be stored in dimx, dimy and dimz.
     */
    ctx->nx = s->dimx;
    ctx->ny = s->dimy;
    ctx->nz = s->dimz;

    ctx->nt = s->iter;
    ctx->nrcv = s->rcv_len;

    CHK(
        ctx->nx <= 0 ||
        ctx->ny <= 0 ||
        ctx->nz <= 0,
        "gpu_wave_tb_init received invalid model dimensions");

    CHK(
        ctx->nt <= 0,
        "gpu_wave_tb_init received an invalid number of time steps");

    CHK(
        ctx->nrcv < 0,
        "gpu_wave_tb_init received an invalid number of receivers");

    ctx->nxyz =
        (size_t)ctx->nx *
        (size_t)ctx->ny *
        (size_t)ctx->nz;

    ctx->field_bytes =
        ctx->nxyz * sizeof(float);

    ctx->coef_count = coef_count;

    ctx->coef_bytes =
        (size_t)ctx->coef_count * sizeof(float);

    ctx->receiver_index_bytes =
        (size_t)ctx->nrcv *
        sizeof(unsigned int);

    /*
     * Store one float for every receiver and every time step.
     *
     * If your working SB implementation records fewer time samples,
     * replace ctx->nt here with the exact SB output-sample count.
     */
    ctx->sismos_bytes =
        (size_t)ctx->nrcv *
        (size_t)ctx->nt *
        sizeof(float);

    GPU_CHK(cudaSetDevice(ctx->device_id));

    if (s->verbose)
    {
        MSG("... GPU TB numerical subsystem initialized");
        MSG(
            "... GPU TB context dimensions: %d x %d x %d",
            ctx->nx,
            ctx->ny,
            ctx->nz);
        MSG(
            "... GPU TB grid points: %zu",
            ctx->nxyz);
    }

    return 0;
}


extern "C" int gpu_wave_tb_allocate(
    gpu_tb_ctx_t *ctx)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_allocate received a NULL context pointer");

    CHK(
        ctx->field_bytes == 0,
        "gpu_wave_tb_allocate received an uninitialized context");

    MSG("... allocating GPU TB arrays");

    /*
     * Static model arrays.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_vel),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_rho),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_inv_rho),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_roc2),
        ctx->field_bytes));

    /*
     * First-order acoustic wavefields.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_p1),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_p2),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_p3),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_v1),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_v2),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_v3),
        ctx->field_bytes));

    /*
     * Finite-difference coefficients.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_coefx),
        ctx->coef_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_coefy),
        ctx->coef_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_coefz),
        ctx->coef_bytes));

    /*
     * Receiver geometry and seismogram output.
     */
    if (ctx->nrcv > 0)
    {
        GPU_CHK(cudaMalloc(
            reinterpret_cast<void **>(&ctx->d_rcv),
            ctx->receiver_index_bytes));

        GPU_CHK(cudaMalloc(
            reinterpret_cast<void **>(&ctx->d_sismos),
            ctx->sismos_bytes));
    }

    /*
     * Four model arrays:
     *   vel, rho, inv_rho, roc2
     *
     * Six wavefield arrays:
     *   p1, p2, p3, v1, v2, v3
     */
    ctx->allocated_bytes =
        (size_t)10 * ctx->field_bytes +
        (size_t)3 * ctx->coef_bytes +
        ctx->receiver_index_bytes +
        ctx->sismos_bytes;

    MSG("... GPU TB allocation successful");

    return 0;
}


extern "C" int gpu_wave_tb_copy_static_data(
    gpu_tb_ctx_t *ctx,
    const float *vel,
    const float *rho,
    const float *inv_rho,
    const float *roc2,
    const float *coefx,
    const float *coefy,
    const float *coefz,
    const unsigned int *rcv)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_copy_static_data received a NULL context pointer");

    CHK(
        vel == NULL,
        "gpu_wave_tb_copy_static_data received a NULL velocity pointer");

    CHK(
        rho == NULL,
        "gpu_wave_tb_copy_static_data received a NULL density pointer");

    CHK(
        inv_rho == NULL,
        "gpu_wave_tb_copy_static_data received a NULL inverse-density pointer");

    CHK(
        roc2 == NULL,
        "gpu_wave_tb_copy_static_data received a NULL roc2 pointer");

    CHK(
        coefx == NULL ||
        coefy == NULL ||
        coefz == NULL,
        "gpu_wave_tb_copy_static_data received a NULL coefficient pointer");

    CHK(
        ctx->d_vel == NULL,
        "GPU TB arrays have not been allocated");

    /*
     * Copy static model arrays.
     */
    GPU_CHK(cudaMemcpy(
        ctx->d_vel,
        vel,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_rho,
        rho,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_inv_rho,
        inv_rho,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_roc2,
        roc2,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    /*
     * Copy finite-difference coefficients.
     */
    GPU_CHK(cudaMemcpy(
        ctx->d_coefx,
        coefx,
        ctx->coef_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_coefy,
        coefy,
        ctx->coef_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_coefz,
        coefz,
        ctx->coef_bytes,
        cudaMemcpyHostToDevice));

    /*
     * Copy receiver indices only when receivers exist.
     */
    if (ctx->nrcv > 0)
    {
        CHK(
            rcv == NULL,
            "gpu_wave_tb_copy_static_data received a NULL receiver pointer");

        GPU_CHK(cudaMemcpy(
            ctx->d_rcv,
            rcv,
            ctx->receiver_index_bytes,
            cudaMemcpyHostToDevice));
    }

    MSG("... GPU TB static arrays copied to the device");

    return 0;
}


extern "C" int gpu_wave_tb_zero(
    gpu_tb_ctx_t *ctx)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_zero received a NULL context pointer");

    CHK(
        ctx->d_p1 == NULL ||
        ctx->d_p2 == NULL ||
        ctx->d_p3 == NULL ||
        ctx->d_v1 == NULL ||
        ctx->d_v2 == NULL ||
        ctx->d_v3 == NULL,
        "gpu_wave_tb_zero received unallocated wavefield pointers");

    /*
     * Initialize all first-order wavefields to zero.
     */
    GPU_CHK(cudaMemset(
        ctx->d_p1,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_p2,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_p3,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_v1,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_v2,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_v3,
        0,
        ctx->field_bytes));

    if (ctx->d_sismos != NULL)
    {
        GPU_CHK(cudaMemset(
            ctx->d_sismos,
            0,
            ctx->sismos_bytes));
    }

    /*
     * Detect asynchronous CUDA errors before returning.
     */
    GPU_CHK(cudaDeviceSynchronize());

    MSG("... GPU TB wavefields initialized to zero");

    return 0;
}


extern "C" void gpu_wave_tb_info(
    const gpu_tb_ctx_t *ctx,
    const sismap_t *s)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_info received a NULL context pointer");

    CHK(
        s == NULL,
        "gpu_wave_tb_info received a NULL sismap pointer");

    cudaDeviceProp props;

    size_t free_bytes = 0;
    size_t total_bytes = 0;

    GPU_CHK(cudaGetDeviceProperties(
        &props,
        ctx->device_id));

    GPU_CHK(cudaMemGetInfo(
        &free_bytes,
        &total_bytes));

    MSG("... GPU TB implementation: Stage 2 memory setup");
    MSG("... GPU TB device: %s", props.name);
    MSG("... GPU TB device index: %d", ctx->device_id);

    MSG(
        "... GPU TB compute capability: %d.%d",
        props.major,
        props.minor);

    MSG(
        "... GPU TB model dimensions: %d x %d x %d",
        ctx->nx,
        ctx->ny,
        ctx->nz);

    MSG(
        "... GPU TB grid points: %zu",
        ctx->nxyz);

    MSG(
        "... GPU TB time steps: %d",
        ctx->nt);

    MSG(
        "... GPU TB receivers: %d",
        ctx->nrcv);

    MSG(
        "... GPU TB one-field size: %.3f MiB",
        (double)ctx->field_bytes /
        (1024.0 * 1024.0));

    MSG(
        "... GPU TB receiver output size: %.3f MiB",
        (double)ctx->sismos_bytes /
        (1024.0 * 1024.0));

    MSG(
        "... GPU TB tracked allocation: %.3f MiB",
        (double)ctx->allocated_bytes /
        (1024.0 * 1024.0));

    MSG(
        "... GPU memory free/total: %.3f / %.3f MiB",
        (double)free_bytes /
        (1024.0 * 1024.0),
        (double)total_bytes /
        (1024.0 * 1024.0));

    if (s->verbose)
    {
        MSG(
            "... GPU TB number of shots: first=%u last=%u",
            s->first,
            s->last);
    }
}


extern "C" void gpu_wave_tb_release(
    gpu_tb_ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }

    /*
     * Static model arrays.
     */
    GPU_TB_FREE(ctx->d_vel);
    GPU_TB_FREE(ctx->d_rho);
    GPU_TB_FREE(ctx->d_inv_rho);
    GPU_TB_FREE(ctx->d_roc2);

    /*
     * Wavefield arrays.
     */
    GPU_TB_FREE(ctx->d_p1);
    GPU_TB_FREE(ctx->d_p2);
    GPU_TB_FREE(ctx->d_p3);

    GPU_TB_FREE(ctx->d_v1);
    GPU_TB_FREE(ctx->d_v2);
    GPU_TB_FREE(ctx->d_v3);

    /*
     * Coefficient arrays.
     */
    GPU_TB_FREE(ctx->d_coefx);
    GPU_TB_FREE(ctx->d_coefy);
    GPU_TB_FREE(ctx->d_coefz);

    /*
     * Receiver arrays.
     */
    GPU_TB_FREE(ctx->d_rcv);
    GPU_TB_FREE(ctx->d_sismos);

    ctx->allocated_bytes = 0;

    MSG("... GPU TB memory released");
}


#undef GPU_TB_FREE