#include <hip/hip_runtime.h>
#include <cstdio>

__global__ void hip_test_kernel(float *data, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        data[i] *= 2.0f;
    }
}

extern "C" int hip_test()
{
    constexpr int n = 1024;
    constexpr size_t bytes = n * sizeof(float);

    float *d_data = nullptr;

    hipError_t err = hipMalloc(
        reinterpret_cast<void **>(&d_data), bytes
    );

    if (err != hipSuccess) {
        std::fprintf(stderr, "hipMalloc failed: %s\n",
                     hipGetErrorString(err));
        return 1;
    }

    err = hipMemset(d_data, 0, bytes);
    if (err != hipSuccess) {
        std::fprintf(stderr, "hipMemset failed: %s\n",
                     hipGetErrorString(err));
        hipFree(d_data);
        return 1;
    }

    hipLaunchKernelGGL(
        hip_test_kernel,
        dim3((n + 255) / 256),
        dim3(256),
        0,
        0,
        d_data,
        n
    );

    err = hipGetLastError();
    if (err == hipSuccess) {
        err = hipDeviceSynchronize();
    }

    hipError_t free_err = hipFree(d_data);
    if (free_err != hipSuccess) {
        std::fprintf(stderr, "hipFree failed: %s\n",
                    hipGetErrorString(free_err));
        return 1;
    }

    if (err != hipSuccess) {
        std::fprintf(stderr, "HIP kernel failed: %s\n",
                     hipGetErrorString(err));
        return 1;
    }

    std::printf("HIP kernel executed successfully\n");
    return 0;
}


/// @brief Computes the pressure on each element of the domain
/// using an OpenCL kernel on the GPU
/// @param gpu_w is a pointer to the wave descriptor to be updated
void gpu_wave_update_fields(sismap_t *s, float*,
                            float*, float*, float*, float* );

/// @brief Adds an impulse to the grid point situated at the
/// source location using an OpenCL kernel on the GPUn
/// @param gpu_w is a pointer to the wave descriptor to be updated
/// @param sterm is the source term to add at the current time step.
void gpu_wave_update_source(sismap_t *s, shot_t *shot, float *, float);

void gpu_wave_image_condition(sismap_t* s, float *d_u1,
                              float *d_fwd,
                              float *d_img, float *d_ilm, unsigned int t);
void gpu_wave_image_gather(sismap_t* s, float *d_img_shot,
                           float *d_ilm_shot, float *d_img);

/// @brief Switches the @ref gpu_wave_t::d_u0 and @ref gpu_wave_t::d_u1 pointers
/// @param gpu_w is a pointer to the wave descriptor to be updated
#define GPU_WAVE_SWAP_POINTERS(u0, u1) \
  float* tmp = u0;                     \
  u0 = u1;                             \
  u1 = tmp

/// @brief Retrieves data from the GPU and saves a snapshot of the
/// wave fields. This can be improved by using asynchronous communications
/// @param gpu_w is a pointer to the wave descriptor to be save
void gpu_wave_save_snapshot(sismap_t *s, shot_t *shot,
                            float *, float *, unsigned int t);

void gpu_wave_read_snapshot(sismap_t* s, shot_t *shot,
                            float *d_fwd, float *fwd, unsigned int t);

/// @brief Retrieves the sismos from the GPU and saves them on disk.
void gpu_wave_save_sismos(sismap_t* s, shot_t *shot,
                          float *d_sismos, float *sismos);

/// @brief read the sismos and send to the GPU.
void gpu_wave_read_sismos(sismap_t* s, shot_t *shot,
                          float *d_sismos, float *sismos);

void gpu_wave_extract_sismos(sismap_t* s,
                             float *d_u1, unsigned int t, float *d_sismos);

void gpu_wave_inject_sismos(sismap_t* s,
                            float *d_u1, unsigned int t, float *d_sismos);

/// @brief Creates a GPU resources descriptor
/// @param p is the parser that may contain a list of additional options
void gpu_wave_set(int);

void gpu_wave_unset();

void gpu_wave_init(sismap_t *s);

/// @brief Releases all the CUDA resources properly
void gpu_wave_release(sismap_t *s);

/// @brief Prints informations about the used OpenCL resources
void gpu_wave_info(sismap_t *s);

void gpu_wave_save_fwd_dbg(sismap_t* s,
                           shot_t *shot,
                           float *d_u1, float* u1, unsigned int t);
void gpu_wave_save_bwd_dbg(sismap_t* s,
                           shot_t *shot,
                           float *d_u1, float* u1, unsigned int t);
void gpu_wave_save_img(sismap_t* s,
                       shot_t *shot,
                       float *d_img, float *d_ilm, float *tmp);