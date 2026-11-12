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