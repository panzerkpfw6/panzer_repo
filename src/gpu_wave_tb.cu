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
    const sismap_t *s)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_init received a NULL context pointer");

    CHK(
        s == NULL,
        "gpu_wave_tb_init received a NULL sismap pointer");

    memset(ctx, 0, sizeof(*ctx));

    ctx->device_id = s->device;

<<<<<<< HEAD
    ctx->nx = (int)s->dimx;
    ctx->ny = (int)s->dimy;
    ctx->nz = (int)s->dimz;

    ctx->nt = s->time_steps;
    ctx->nrcv = (int)s->rcv_len;

    /*
     * The active first-order stencil uses four coefficients
     * on either side of the staggered derivative.
     */
    ctx->stencil_radius = 4;
=======
    ctx->nx = s->dimx;
    ctx->ny = s->dimy;
    ctx->nz = s->dimz;

    ctx->nt = s->time_steps;
    ctx->nrcv = s->rcv_len;
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1

    CHK(
        ctx->nx <= 0 ||
        ctx->ny <= 0 ||
        ctx->nz <= 0,
        "gpu_wave_tb_init received invalid dimensions");

    CHK(
        ctx->nt <= 0,
        "gpu_wave_tb_init received an invalid time-step count");

    /*
     * If nrcv is unsigned in gpu_tb_ctx_t, this check is unnecessary.
     * It is harmless only when nrcv is signed.
     */
    CHK(
        ctx->nrcv < 0,
        "gpu_wave_tb_init received an invalid receiver count");

    ctx->nxyz =
        (size_t)ctx->nx *
        (size_t)ctx->ny *
        (size_t)ctx->nz;

    /*
     * Keep this aligned with the host arrays passed to
     * run_modeling_1st_tb_cpu().
     *
     * The CPU driver allocates u0/vx/vy/vz using s->size.
     * Therefore, use s->size rather than nx*ny*nz when they differ.
     */
    ctx->field_bytes =
        (size_t)s->size * sizeof(float);

    ctx->source_bytes =
        (size_t)ctx->nt * sizeof(float);
<<<<<<< HEAD
=======

    ctx->coef_count = coef_count;

    ctx->coef_bytes =
        (size_t)ctx->coef_count * sizeof(float);
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1

    ctx->receiver_index_bytes =
        (size_t)ctx->nrcv *
        sizeof(unsigned int);

<<<<<<< HEAD
=======
    /*
     * Store one float for every receiver and every time step.
     *
     * If the final solver records only decimated output samples,
     * replace ctx->nt with the actual number of recorded samples.
     */
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
    ctx->sismos_bytes =
        (size_t)ctx->nrcv *
        (size_t)(ctx->nt + 1) *
        sizeof(float);

    GPU_CHK(cudaSetDevice(ctx->device_id));

    if (s->verbose)
    {
        MSG("... GPU TB numerical subsystem initialized");

        MSG(
            "... GPU TB physical dimensions: %d x %d x %d",
            ctx->nx,
            ctx->ny,
            ctx->nz);

        MSG(
            "... GPU TB allocated field elements: %zu",
            (size_t)s->size);
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
        "gpu_wave_tb_allocate received an invalid field size");

    CHK(
        ctx->source_bytes == 0,
        "gpu_wave_tb_allocate received an invalid source size");

    MSG("... allocating GPU TB arrays");

    /*
     * Four dynamic first-order fields.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_u0),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_vx),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_vy),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_vz),
        ctx->field_bytes));

    /*
     * Static physical parameters.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_roc2),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_inv_rho),
        ctx->field_bytes));

    /*
     * Damping arrays.
     *
     * We initially allocate full-field arrays because that is the
     * simplest safe representation. This can be reduced later if
     * ctx->dampx/dampy/dampz are one-dimensional arrays.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_dampx),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_dampy),
        ctx->field_bytes));

    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_dampz),
        ctx->field_bytes));

    /*
<<<<<<< HEAD
     * Source and receivers.
=======
     * Source time function.
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_source),
        ctx->source_bytes));

    /*
     * First-order acoustic wavefields.
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
     */
    GPU_CHK(cudaMalloc(
        reinterpret_cast<void **>(&ctx->d_source),
        ctx->source_bytes));

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
<<<<<<< HEAD
     * Nine full fields:
     *
     *   u0, vx, vy, vz
     *   roc2, inv_rho
     *   dampx, dampy, dampz
     */
    ctx->allocated_bytes =
        (size_t)9 * ctx->field_bytes +
        ctx->source_bytes +
=======
     * Four full model arrays:
     *   vel, rho, inv_rho, roc2
     *
     * Six full wavefield arrays:
     *   p1, p2, p3, v1, v2, v3
     *
     * Plus:
     *   source, coefficients, receiver indices and seismograms.
     */
    ctx->allocated_bytes =
        (size_t)10 * ctx->field_bytes +
        ctx->source_bytes +
        (size_t)3 * ctx->coef_bytes +
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
        ctx->receiver_index_bytes +
        ctx->sismos_bytes;

    MSG("... GPU TB allocation successful");

    return 0;
}


extern "C" int gpu_wave_tb_copy_static_data(
    gpu_tb_ctx_t *ctx,
<<<<<<< HEAD
    const sismap_t *s,
    const float *roc2,
=======
    const float *vel,
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
    const float *inv_rho,
    const float *source,
    const unsigned int *rcv)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_copy_static_data received a NULL context pointer");

    CHK(
<<<<<<< HEAD
        s == NULL,
        "gpu_wave_tb_copy_static_data received a NULL sismap pointer");

    CHK(
        roc2 == NULL,
        "gpu_wave_tb_copy_static_data received a NULL roc2 pointer");

    CHK(
        inv_rho == NULL,
        "gpu_wave_tb_copy_static_data received a NULL inv_rho pointer");

    CHK(
        source == NULL,
        "gpu_wave_tb_copy_static_data received a NULL source pointer");

    CHK(
        s->dampx == NULL ||
        s->dampy == NULL ||
        s->dampz == NULL,
        "gpu_wave_tb_copy_static_data received NULL damping arrays");

    if (ctx->nrcv > 0)
    {
        CHK(
            rcv == NULL,
            "gpu_wave_tb_copy_static_data received a NULL receiver array");
    }

    MSG("... copying GPU TB static numerical arrays");

    GPU_CHK(cudaMemcpy(
        ctx->d_roc2,
        roc2,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
=======
        vel == NULL,
        "gpu_wave_tb_copy_static_data received a NULL velocity pointer");

    CHK(
        inv_rho == NULL,
        "gpu_wave_tb_copy_static_data received a NULL inverse-density pointer");

    CHK(
        source == NULL,
        "gpu_wave_tb_copy_static_data received a NULL source pointer");

    CHK(
        ctx->d_vel == NULL,
        "gpu_wave_tb_copy_static_data received an unallocated d_vel");

    CHK(
        ctx->d_inv_rho == NULL,
        "gpu_wave_tb_copy_static_data received an unallocated d_inv_rho");

    CHK(
        ctx->d_source == NULL,
        "gpu_wave_tb_copy_static_data received an unallocated d_source");

    CHK(
        ctx->field_bytes == 0,
        "gpu_wave_tb_copy_static_data received an invalid field size");

    CHK(
        ctx->source_bytes == 0,
        "gpu_wave_tb_copy_static_data received an invalid source size");

    if (ctx->nrcv > 0)
    {
        CHK(
            rcv == NULL,
            "gpu_wave_tb_copy_static_data received a NULL receiver array");

        CHK(
            ctx->d_rcv == NULL,
            "gpu_wave_tb_copy_static_data received an unallocated d_rcv");
    }

    MSG(
        "... copying GPU TB static arrays: field=%zu bytes, source=%zu bytes",
        ctx->field_bytes,
        ctx->source_bytes);

    GPU_CHK(cudaMemcpy(
        ctx->d_vel,
        vel,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
        ctx->d_inv_rho,
        inv_rho,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    /*
     * This assumes s->dampx, dampy and dampz each contain s->size
     * elements. Verify their actual allocation before running kernels.
     */
    GPU_CHK(cudaMemcpy(
<<<<<<< HEAD
        ctx->d_dampx,
        s->dampx,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_dampy,
        s->dampy,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_dampz,
        s->dampz,
        ctx->field_bytes,
        cudaMemcpyHostToDevice));

    GPU_CHK(cudaMemcpy(
        ctx->d_source,
        source,
        ctx->source_bytes,
        cudaMemcpyHostToDevice));

=======
        ctx->d_source,
        source,
        ctx->source_bytes,
        cudaMemcpyHostToDevice));

>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
    if (ctx->nrcv > 0)
    {
        GPU_CHK(cudaMemcpy(
            ctx->d_rcv,
            rcv,
            ctx->receiver_index_bytes,
            cudaMemcpyHostToDevice));
    }

<<<<<<< HEAD
    GPU_CHK(cudaDeviceSynchronize());

    MSG("... GPU TB static numerical arrays copied");
=======
    /*
     * cudaMemcpy is synchronous for ordinary pageable host memory,
     * but this also exposes any previous CUDA runtime error here.
     */
    GPU_CHK(cudaDeviceSynchronize());

    MSG("... GPU TB static arrays copied to the device");
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1

    return 0;
}


extern "C" int gpu_wave_tb_zero(
    gpu_tb_ctx_t *ctx)
{
    CHK(
        ctx == NULL,
        "gpu_wave_tb_zero received a NULL context pointer");

    CHK(
        ctx->d_u0 == NULL ||
        ctx->d_vx == NULL ||
        ctx->d_vy == NULL ||
        ctx->d_vz == NULL,
        "gpu_wave_tb_zero received unallocated dynamic fields");

<<<<<<< HEAD
=======
    /*
     * Initialize all dynamic first-order wavefields to zero.
     *
     * Do not clear d_vel, d_inv_rho, d_source or d_rcv here because
     * they already contain copied static data.
     */
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
    GPU_CHK(cudaMemset(
        ctx->d_u0,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_vx,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_vy,
        0,
        ctx->field_bytes));

    GPU_CHK(cudaMemset(
        ctx->d_vz,
        0,
        ctx->field_bytes));

    if (ctx->d_sismos != NULL &&
        ctx->sismos_bytes > 0)
    {
        GPU_CHK(cudaMemset(
            ctx->d_sismos,
            0,
            ctx->sismos_bytes));
    }

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

    MSG("... GPU TB implementation: Stage 3A numerical data setup");
    MSG("... GPU TB device: %s", props.name);
    MSG("... GPU TB device index: %d", ctx->device_id);

    MSG(
        "... GPU TB compute capability: %d.%d",
        props.major,
        props.minor);

    MSG(
        "... GPU TB dimensions: %d x %d x %d",
        ctx->nx,
        ctx->ny,
        ctx->nz);

    MSG(
        "... GPU TB field allocation: %.3f MiB",
        (double)ctx->field_bytes /
        (1024.0 * 1024.0));

    MSG(
<<<<<<< HEAD
        "... GPU TB source allocation: %.6f MiB",
=======
        "... GPU TB source size: %.6f MiB",
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
        (double)ctx->source_bytes /
        (1024.0 * 1024.0));

    MSG(
<<<<<<< HEAD
        "... GPU TB receiver output: %.3f MiB",
=======
        "... GPU TB receiver output size: %.3f MiB",
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1
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

    MSG(
        "... GPU TB number of shots: first=%d last=%d",
        s->first,
        s->last);
}


extern "C" void gpu_wave_tb_release(
    gpu_tb_ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }

    GPU_TB_FREE(ctx->d_u0);
    GPU_TB_FREE(ctx->d_vx);
    GPU_TB_FREE(ctx->d_vy);
    GPU_TB_FREE(ctx->d_vz);

    GPU_TB_FREE(ctx->d_roc2);
    GPU_TB_FREE(ctx->d_inv_rho);

<<<<<<< HEAD
    GPU_TB_FREE(ctx->d_dampx);
    GPU_TB_FREE(ctx->d_dampy);
    GPU_TB_FREE(ctx->d_dampz);
=======
    /*
     * Source time function.
     */
    GPU_TB_FREE(ctx->d_source);

    /*
     * Wavefield arrays.
     */
    GPU_TB_FREE(ctx->d_p1);
    GPU_TB_FREE(ctx->d_p2);
    GPU_TB_FREE(ctx->d_p3);
>>>>>>> d61c5e3146b8862ad57e002cfcde2472247c2fb1

    GPU_TB_FREE(ctx->d_source);

    GPU_TB_FREE(ctx->d_rcv);
    GPU_TB_FREE(ctx->d_sismos);

    ctx->allocated_bytes = 0;

    MSG("... GPU TB memory released");
}


#undef GPU_TB_FREE