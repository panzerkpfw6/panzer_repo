#define _GNU_SOURCE
#include <sched.h>
#include <sys/sysinfo.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include "stencil/macros.h"
#include "stencil/sismap.h"
#include "stencil/parser.h"
#include "stencil/wave_tb.h"
#include "stencil/wtime.h"

////////////////////////
// Paste from stencil-dev. A bit different notattion for TB.
typedef float hFloat ;
typedef float real_t;
#define MPI_real_t MPI_FLOAT
enum Stencil_Shapes{
	STAR,
	TTI,
	BOz
};

enum Stencil_Coefficients{
	CONSTANT_COEFFICIENT,
	VARIABLE_COEFFICIENT,
	VARIABLE_COEFFICIENT_AXSYM,
	VARIABLE_COEFFICIENT_NOSYM,
	SOLAR_COEFFICIENT
};

enum Stencil_Type{
	REGULAR,
	SOLAR
};

// Fields in solar kernel grid cell
enum Solar_Fields{
	ALL_FIELDS,
	H_FIELD,
	E_FIELD,
};


#ifndef _OPENMP
#define  omp_get_num_threads() (1)
#define omp_set_nested(a) {}
#define omp_get_thread_num() (0)
#define omp_get_max_threads() (1)
#endif

// Profiling
typedef struct{
    double compute, communicate, send_recv, wait, total, others, ts_main, ts_others;
}Profile;

// Halo information
typedef struct{
    int shape[3], recv_b[3], recv_e[3], send_b[3], send_e[3], is_contiguous, size;
}Halo;

// MPI topology information
typedef struct{
    int right, // z+
    left,  // z-
    up,    // y+
    down,  // y-
    front, // x+
    back;  // x-
    int shape[3],is_periodic[3],rank_coords[3];
    //  MPI_Request wait_req[8];
    //  MPI_Status wait_stat16[16],wait_stat8[8],wait_stat4[4];
}mpi_topology;

typedef struct{
    int nnz, nny, nnx;
    uint64_t ln_domain;
}CLU_ctx;

#define CLU_SIG (const CLU_ctx clu_ctx, const int zb, const int ze, const int j, const int k, \
		const real_t *  coef, hFloat *  u, \
		const hFloat *  v, const hFloat *  roc2)
typedef void (*clu_func_t)CLU_SIG;

// contezt information passed to the stencil kernel
typedef struct{
    int bs_z; //depricated
    int bs_y; // for spatial blocking in Y at the standard methods
    int thread_group_size;
    int th_z, th_y, th_x, th_c; // number of threads per dimension in z, y, and x, and per component

    // cpu binding masks
    cpu_set_t **bind_masks;
    int setsize;
    int use_manual_cpu_bind;

    // for separate stride-1 functions
    clu_func_t clu_func;

    int num_wf; // number of wavefront updats per iteration
    // wavefront profiling
    double *t_wf_comm, *t_wait, *t_wf_prologue, *t_wf_main, *t_wf_epilogue, *wf_num_resolved_diamonds, *t_group_wait;
    real_t idz, idy, idx, idzyx_sum; //@Rqched 1/dz 1/dy 1/dx

}stencil_ctx;

////// Kernels and time steppers data structures
#define KERNEL_SIG ( const int shape[3], const int zb, const int yb,  const int xb, const int ze, const int ye, const int xe,\
		const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13, const hFloat *  p21, const hFloat *  p22, const hFloat *  p23, const hFloat *  roc2, int mod, stencil_ctx stencil_ctx);
#define KERNEL_MWD_SIG ( const int shape[3], const int zb, const int yb_r, const int xb, const int ze, const int ye_r, const int xe, \
		const real_t *  coef, hFloat *  u1,hFloat *  u2,hFloat *  u3, \
		hFloat *  v1, hFloat *  v2,hFloat *  v3, const hFloat *  roc2,\
        float *  dampx, float *  dampy,float *  dampz,\
        int t_dim, int b_inc, int e_inc, int NHALO, int tb, int te, stencil_ctx stencil_ctx, int mtid)
typedef void (*spt_blk_func_t)KERNEL_SIG;
typedef void (*mwd_func_t)KERNEL_MWD_SIG;
//////

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
struct StencilInfo {
    const char *name;
    int r;
    int time_order;
    int nd;
    enum Stencil_Shapes shape;
    enum Stencil_Coefficients coeff;
    enum Stencil_Type type;
};

// contezt information
typedef struct{
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
    int lstencil_shape[3], ldomain_shape[3], gb[3], ge[3], lsource_pt[3], has_source; //MPI ranks' global indices, and local source locations

    int notuning; //@KADIR returns early from autuning functions
    int call_combined_function; //@KADIR calls Kadir's function that does everything

    //  int stencil_radius, is_constant_coefficient;
    //  enum Stencil_Types stencil_type;

    hFloat *  U1, *  U2, *  U3,*  U4, *  U5, *  source;
    hFloat *  rU1, *  rU2, *  rU3,*  rU4, *  rU5;

    // damping ABCs
    float * dampx;
    float * dampy;
    float * dampz;
    real_t *  coef;

    real_t *  src_exc_coef; //@KADIR: coef used in source ezcitation. length is number of time steps (nt)
    // Use source instead of src_exc_conf after verification results

    // parameters for internal thread affinity
    int th_block;
    int th_stride;

    // Holds the value of cache blocking across Y axis
    stencil_ctx stencil_ctx;

    // to enable/disable source point update
    int source_point_enabled;

    int cache_size; // Last level cache usable size in KB for blocking

    // to enable concatenating halo information before communication
    int halo_concat;

    // Specific data for the diamond method
    int t_dim, larger_t_dim, is_last, mwd_type;
    Halo hu[3], hv[3];

    int wavefront;
    uint64_t idiamond_pro_epi_logue_updates;
    uint64_t wf_blk_size, wf_larger_blk_size;

    Halo h[3]; // Halo information for z,Y, and x directions
    mpi_topology t;
    Profile prof;

    struct Stencil stencil;

    // list of coefficients to be used in stencil operators
    real_t g_coef[11];

    int array_padding;

    int in_auto_tuning;
    int orig_thread_group_size; // to distingquish whether thread group size is set by the user

}Parameters;
////////////////////////