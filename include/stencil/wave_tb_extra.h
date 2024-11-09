// wave_tb.h

#ifndef WAVE_TB_H  // Include guard to prevent multiple inclusions
#define WAVE_TB_H

#include <stdint.h>
#include <sched.h>
#include <omp.h>
#include "stencil/macros.h"
#include "stencil/sismap.h"
#include "stencil/parser.h"
#include "stencil/wtime.h"

// Typedefs for float types
typedef float hFloat;
typedef float real_t;

// Profiling structure
typedef struct {
    double compute, communicate, send_recv, wait, total, others, ts_main, ts_others;
} Profile;

// Halo information structure
typedef struct {
    int shape[3], recv_b[3], recv_e[3], send_b[3], send_e[3], is_contiguous, size;
} Halo;

// MPI topology structure
typedef struct {
    int right, left, up, down, front, back;
    int shape[3], is_periodic[3], rank_coords[3];
} mpi_topology;

// CLU context structure
typedef struct {
    int nnz, nny, nnx;
    uint64_t ln_domain;
} CLU_ctx;

// Function signature for CLU function
#define CLU_SIG (const CLU_ctx clu_ctx, const int zb, const int ze, const int j, const int k, \
        const real_t *coef, hFloat *u, const hFloat *v, const hFloat *roc2)

typedef void (*clu_func_t)CLU_SIG;
// Function signatures for kernel functions
#define KERNEL_SIG (const int shape[3], const int zb, const int yb, const int xb, const int ze, const int ye, const int xe, \
    const real_t *coef, hFloat *p11, hFloat *p12, hFloat *p13, const hFloat *p21, const hFloat *p22, const hFloat *p23, const hFloat *roc2, \
    int mod, stencil_ctx stencil_ctx)

// Stencil context structure
typedef struct {
    int bs_z;
    int bs_y;
    int thread_group_size;
    int th_z, th_y, th_x, th_c;
    cpu_set_t **bind_masks;
    int setsize;
    int use_manual_cpu_bind;
    clu_func_t clu_func;
    int num_wf;
    double *t_wf_comm, *t_wait, *t_wf_prologue, *t_wf_main, *t_wf_epilogue;
    real_t idz, idy, idx, idzyx_sum;
} stencil_ctx;
#define KERNEL_MWD_SIG (const int shape[3], const int zb, const int yb_r, const int xb, const int ze, const int ye_r, const int xe, \
    const real_t *coef, hFloat *u1, hFloat *u2, hFloat *u3, hFloat *v1, hFloat *v2, hFloat *v3, const hFloat *roc2, \
    float *dampx, float *dampy, float *dampz, int t_dim, int b_inc, int e_inc, int NHALO, int tb, int te, stencil_ctx stencil_ctx, int mtid)

typedef void (*spt_blk_func_t)KERNEL_SIG;
typedef void (*mwd_func_t)KERNEL_MWD_SIG;

// Stencil structure
struct Stencil {
    const char *name;
    int r;
    int time_order;
    int nd;
    enum Stencil_Shapes shape;
    enum Stencil_Coefficients coeff;
    enum Stencil_Type type;
    spt_blk_func_t spt_blk_func;
    spt_blk_func_t stat_sched_func;
    mwd_func_t mwd_func;
};

// Stencil information structure (similar to Stencil)
struct StencilInfo {
    const char *name;
    int r;
    int time_order;
    int nd;
    enum Stencil_Shapes shape;
    enum Stencil_Coefficients coeff;
    enum Stencil_Type type;
};


// Parameters structure for stencil
typedef struct {
    int alignment, verbose, stencil_shape[3];
    uint64_t n_stencils, ln_domain, ln_stencils;
    int target_ts, target_kernel;
    int mpi_rank, mpi_size;
    int n_tests, nt;
    int verify;
    int source_pt[3];
    int debug;
    int num_threads;
    int use_omp_stat_sched;
    int lstencil_shape[3], ldomain_shape[3], gb[3], ge[3], lsource_pt[3], has_source;
    int notuning;
    int call_combined_function;
    hFloat *U1, *U2, *U3, *U4, *U5, *source;
    hFloat *rU1, *rU2, *rU3, *rU4, *rU5;
    float *dampx, *dampy, *dampz;
    real_t *coef;
    real_t *src_exc_coef;
    int th_block;
    int th_stride;
    stencil_ctx stencil_ctx;
    int source_point_enabled;
    int cache_size;
    int halo_concat;
    int t_dim, larger_t_dim, is_last, mwd_type;
    Halo hu[3], hv[3];
    int wavefront;
    uint64_t idiamond_pro_epi_logue_updates;
    uint64_t wf_blk_size, wf_larger_blk_size;
    Halo h[3];
    mpi_topology t;
    Profile prof;
    struct Stencil stencil;
    real_t g_coef[11];
    int array_padding;
    int in_auto_tuning;
    int orig_thread_group_size;
} Parameters;

#endif  // WAVE_TB_H
