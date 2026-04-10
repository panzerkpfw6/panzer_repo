///
/// @copyright Copyright 2024- Pavel Plotnitskii. All rights reserved.
/// This file is part of the \b stencil project.
///
/// \b stencil is free software: you may redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// The stencil project is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with the \b stencil project. If not, see <http://www.gnu.org/licenses/>.
///
/// @author Pavel Plotnitskii
#ifndef __STENCIL_WAVE_TB_H_
#define __STENCIL_WAVE_TB_H_

#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdbool.h>
#include <stencil/config.h>
#include <stencil/sismap.h>
#include <stencil/shot.h>
#include <stencil/parser.h>
#include <stencil/mlbs.h>
#include <stencil/bwriter.h>
////#include <stencil/wave_tb_extra.h>
#include <sys/sysinfo.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include "stencil/macros.h"
////#include "stencil/wave_tb.h"
#include "stencil/wtime.h"

///////////////Definitions/////////////////////////////////

// Define the types for MPI real numbers (assuming MPI_FLOAT is defined in your environment)
#define MPI_real_t MPI_FLOAT

typedef struct tb_s tb_t;
typedef struct tb_data_s tb_data_t;
typedef struct tb_timer_s tb_timer_t;

// Enumeration types for stencil shapes, coefficients, and types
enum Stencil_Shapes {
    STAR,
    TTI,
    BOz
};

enum Stencil_Coefficients {
    CONSTANT_COEFFICIENT,
    VARIABLE_COEFFICIENT,
    VARIABLE_COEFFICIENT_AXSYM,
    VARIABLE_COEFFICIENT_NOSYM,
    SOLAR_COEFFICIENT
};

enum Stencil_Type {
    REGULAR,
    SOLAR
};

// Fields in solar kernel grid cell
enum Solar_Fields {
    ALL_FIELDS,
    H_FIELD,
    E_FIELD
};

// Profiling structure
typedef struct {
    double compute, communicate, send_recv, wait, total, others, ts_main, ts_others;
} Profile;

// Halo information structure
typedef struct {
    int shape[3], recv_b[3], recv_e[3], send_b[3], send_e[3], is_contiguous, size;
} Halo;

// MPI topology information structure
typedef struct {
    int right, left, up, down, front, back;
    int shape[3], is_periodic[3], rank_coords[3];
} mpi_topology;

// CLU context structure
typedef struct {
    int nnz, nny, nnx;
    uint64_t ln_domain;
} CLU_ctx;

// Function pointer type for CLU kernel function
#define CLU_SIG (const CLU_ctx clu_ctx, const int zb, const int ze, const int j, const int k, \
                const real_t *coef, hFloat *u, const hFloat *v, const hFloat *roc2)
typedef void (*clu_func_t)CLU_SIG;
double get_wall_time();
double get_cpu_time();

///////////////////////////////////////////////////////////////

// Context information passed to the stencil kernel
typedef struct {
    int bs_z; // deprecated
    int bs_y; // for spatial blocking in Y at the standard methods
    int thread_group_size;
    int th_z, th_y, th_x, th_c; // number of threads per dimension
    int fwd_steps;// snapshotting step
    float dz, dy, dx; // spacing per dimension
    int nz, ny, nx; // grid size without halo area
    unsigned long fwd_size;

    // cpu binding masks
    cpu_set_t **bind_masks;
    int setsize;
    int use_manual_cpu_bind;

    // For separate stride-1 functions
    clu_func_t clu_func;

    int num_wf; // number of wavefront updates per iteration
    float dt;   // time step

    // wavefront profiling
    double *t_wf_comm, *t_wait, *t_wf_prologue, *t_wf_main, *t_wf_epilogue, *wf_num_resolved_diamonds, *t_group_wait;
    real_t idz, idy, idx, idzyx_sum; // inverse spacing

} stencil_ctx;

////// Kernels and time steppers data structures
#define KERNEL_SIG ( const int shape[3], const int zb, const int yb,  const int xb, const int ze, const int ye, const int xe,\
		const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13, const hFloat *  p21, const hFloat *  p22, const hFloat *  p23, const hFloat *  roc2, int mod, stencil_ctx stencil_ctx);
#define KERNEL_MWD_SIG ( const int shape[3], const int zb, const int yb_r, \
                            const int xb, const int ze, const int ye_r, const int xe, \
                            const real_t *coef, hFloat *u1, hFloat *u2, hFloat *u3,   \
                            hFloat *v1, hFloat *v2, hFloat *v3, \
							const hFloat *roc2, const hFloat *inv_rho, \
                            float *dampx, float *dampy, float *dampz,      \
                            int t_dim, int b_inc, int e_inc, int NHALO,    \
                            int tb, int te,int t0,int ifwd,stencil_ctx stencil_ctx, int mtid,tb_data_t * data)
typedef void (*spt_blk_func_t)KERNEL_SIG;
typedef void (*mwd_func_t)KERNEL_MWD_SIG;


// Stencil structure definition
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

// Stencil info structure (for meta-information)
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

    hFloat *  U1, *  U2, *  U3,*  U4, *  source;
    const float * U5;
    const float * U6;
//    hFloat *  rU1, *  rU2, *  rU3,*  rU4, *  rU5;

    // damping ABCs
    float * dampx;
    float * dampy;
    float * dampz;
    real_t * coef;

    real_t * src_exc_coef; //@KADIR: coef used in source ezcitation. length is number of time steps (nt)
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
    int t_dim,larger_t_dim,is_last,mwd_type,t_len;
    Halo hu[3],hv[3];

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
    tb_data_t * data;
}Parameters;

// Global parameter variable
extern Parameters *gp;  // Global parameter structure within a node
extern real_t *recv_rec;  // Array for receiver recording
extern size_t *irecv_rec; // Index array into recv_rec
extern size_t isrc_exc;  // Number of source excitations performed so far


int get_ntg(Parameters p);  // Function declaration
// Function prototypes
void femwd_iso_ref_1st( const int shape[3], const int zb, const int yb_r0, const int xb,
                        const int ze, const int ye_r0, const int xe,
                        const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13,
                        hFloat *  p21, hFloat *  p22, hFloat *  p23,
						const hFloat *  roc2,const hFloat * inv_rho,
                        float * dampx,float * dampy,float * dampz,
                        int t_dim, int b_inc, int e_inc,int NHALO,
                        int tb, int te,int t0,int ifwd,
						stencil_ctx stencil_ctx, int mtid,tb_data_t * data);
void femwd_iso_ref_1st_grok( const int shape[3], const int zb, const int yb_r0, const int xb,
                        const int ze, const int ye_r0, const int xe,
                        const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13,
                        hFloat *  p21, hFloat *  p22, hFloat *  p23,
						const hFloat *  roc2,const hFloat * inv_rho,
                        float * dampx,float * dampy,float * dampz,
                        int t_dim, int b_inc, int e_inc,int NHALO,
                        int tb, int te,int t0,int ifwd,
						stencil_ctx stencil_ctx, int mtid,tb_data_t * data);
void intra_diamond_mwd_comp_std(Parameters *p, int yb_r, int ye_r, int b_inc, int e_inc, int tb, int te, int tid,int t_coord);
void intra_diamond_mwd_comp(Parameters *p, int yb_r, int ye_r, int b_inc, int e_inc, int tb, int te, int tid,int t0,int ifwd);
void dynamic_intra_diamond_ts_combined(Parameters *p);
void dynamic_intra_diamond_ts_combined_backward(Parameters *p);
void reset_timers(Profile * p);
void reset_wf_timers(Parameters * p);
void cpu_bind_init(Parameters *p);

///////////////////////////////////////////////////////////////

typedef struct tb_s tb_t;
typedef struct tb_data_s tb_data_t;
typedef struct tb_timer_s tb_timer_t;

struct tb_s {
  // threads
  int num_threads;
  int num_thread_groups;
  int thread_group_size;
  int th_x;
  int th_y;
  int th_z; // number of threads per dimension in x, y, and z, and per component
  char * affinity_file;

  // thread affinity
  cpu_set_t **bind_masks;
  size_t setsize;

  // for tiling blocking
  int time_steps;
  int t_dim;  //
  int num_wf; // number of wavefront updates per iteration
  int diam_width;

  // scheduling
  int t_len;
  int y_len_l;
  int y_len_r;

  // stride
  int nnx;
  int nny;
  int nnz;

  // grid size
  int stencilx;
  int stencily;
  int stencilz;

  // coef (copy pointer)
  float *coefx;
  float *coefy;
  float *coefz;

  // halo
  int r;

  volatile int *t_pos;
  volatile int *avail_list;

  // imaging step
  int fwd_steps;
  unsigned long fwd_size;

  // estimation
  unsigned long long nb_stencils_main;
  unsigned long long nb_stencils_total_fwd;
  unsigned long long nb_stencils_total_bwd;

  // damping PML
  float * dampx;
  float * dampy;
  float * dampz;

  // function pointer
  int mode;
  void (*kernel_spatial_blocking_1st)(const int nnx, const int nny, const int nnz,
                                  const int xb,  int yb_r, const int* zb,
                                  const int xe,  int ye_r, const int* ze,
                                  const float * restrict coefx,
                                  const float * restrict coefy,
                                  const float * restrict coefz,
                                  const float * restrict dampx,
                                  const float * restrict dampy,
                                  const float * restrict dampz,
                                  float * restrict u_r,
                                  float * restrict vx_r,
                                  float * restrict vy_r,
                                  float * restrict vz_r,
                                  const float * restrict roc2,
                                  const int t_dim,
                                  int b_inc,int e_inc,
                                  const int stencilr,
                                  const int tb,const int te,
                                  const int thread_group_size,const int groupid,
                                  const int setsize,cpu_set_t ** bind_masks,
                                  const tb_data_t* data,
                                  const int t0,const int ifwd,
                                  tb_timer_t* timer);

  void (*kernel_tiling_blocking_1st)(const int nnx, const int nny, const int nnz,
                                 const int xb, const int yb_r0, const int zb,
                                 const int xe, const int ye_r0, const int ze,
                                 const float * restrict coefx,
                                 const float * restrict coefy,
                                 const float * restrict coefz,
                                 const float * restrict dampx,
                                 const float * restrict dampy,
                                 const float * restrict dampz,
                                 float * restrict u_r,
                                 float * restrict vx_r,
                                 float * restrict vy_r,
                                 float * restrict vz_r,
                                 const float * restrict roc2,
                                 const int t_dim,
                                 int b_inc, int e_inc,
                                 const int stencilr,
                                 const int tb, const int te,
                                 const int num_wf,
                                 const int thread_group_size, const int threadx,
                                 const int thready, const int threadz,
                                 const int groupid,
                                 const int setsize, cpu_set_t ** bind_masks,
                                 tb_data_t* data,
                                 const int t0, const int ifwd,
                                 tb_timer_t* timer);
    void (*kernel_spatial_blocking)(const int nnx, const int nny, const int nnz,
                                    const int xb,  int yb_r, const int* zb,
                                    const int xe,  int ye_r, const int* ze,
                                    const float * restrict coefx,
                                    const float * restrict coefy,
                                    const float * restrict coefz,
                                    const float * restrict dampx,
                                    const float * restrict dampy,
                                    const float * restrict dampz,
                                    float * restrict u_r,
                                    float * restrict v_r,
                                    const float * restrict roc2,
                                    const int t_dim,
                                    int b_inc,int e_inc,
                                    const int stencilr,
                                    const int tb,const int te,
                                    const int thread_group_size,const int groupid,
                                    const int setsize,cpu_set_t ** bind_masks,
                                    const tb_data_t* data,
                                    const int t0,const int ifwd,
                                    tb_timer_t* timer);

    void (*kernel_tiling_blocking)(const int nnx, const int nny, const int nnz,
                                   const int xb, const int yb_r0, const int zb,
                                   const int xe, const int ye_r0, const int ze,
                                   const float * restrict coefx,
                                   const float * restrict coefy,
                                   const float * restrict coefz,
                                   const float * restrict dampx,
                                   const float * restrict dampy,
                                   const float * restrict dampz,
                                   float * restrict u, float * restrict v,
                                   const float * restrict roc2,
                                   const int t_dim,
                                   int b_inc, int e_inc,
                                   const int stencilr,
                                   const int tb, const int te,
                                   const int num_wf,
                                   const int thread_group_size, const int threadx,
                                   const int thready, const int threadz,
                                   const int groupid,
                                   const int setsize, cpu_set_t ** bind_masks,
                                   tb_data_t* data,
                                   const int t0, const int ifwd,
                                   tb_timer_t* timer);
};

struct tb_data_s {
  float * roc2;
  const float* source;
  unsigned int src_idx;
  int src_depth;
  int src_x;
  int src_y;
  int src_z;
  int order;

  float dx,dy,dz;
  float dt;

  float* sismos;
  int rcv_depth;
  unsigned int rcv_len;
  const unsigned int *rcv;
  int *ix, *iy;

  float * wave;

  // image related data
  float * fwd;
  int flag_fwd; // 1 if save fwd snap
  int flag_bwd; // 1 if load fwd snap and build image
  int flag_img;
  int rec_sismos; // 1 if record sismos, 0 if not

  float * img;
  float * ilm;

  int fwd_steps;


  // mode 3
  int *gfwd_y0;
  FILE ** gfd;
  size_t groupsize;
  char gfd_name[512];


  mlbs_t* mlbs0;
  mlbs_t* mlbs1;

  bwriter_t* bwriter0;
  bwriter_t* bwriter1;
};

struct tb_timer_s {
  int num_thread_groups;
  int thread_group_size;
  int num_threads;

  double compute;
  double total;
  double others;
  double ts_main;        // mainloop (diamonds)
  double ts_others;      // prologue and epilogue
  double *t_wait;        // waiting time for openmp barrier inside TB
  double *t_wf_prologue; // wavefront spatial blocking (prologue)
  double *t_wf_mainloop; // wavefront tiling blocking  (mainloop)
  double *t_wf_epilogue; // wavefront spatial blocking (epilogue)
  double *t_wf_snapio; // wavefront tiling blocking  (mainloop)
  double *t_group_wait;  // time used for waiting a new diamond
  double *wf_num_resolved_diamonds;

};

void wave_tb_timer_init(tb_timer_t * timer,
                        const int thread_group_size,
                        const int num_thread_groups);

void wave_tb_timer_clear(tb_timer_t * timer);

void wave_tb_timer_free(tb_timer_t * timer);

void wave_tb_timer_info(tb_timer_t * timer,
                        const int64_t nb_stencils_total,
                        const int64_t nb_stencils_main);

void wave_tb_info(tb_t * ctx);

void wave_tb_init(tb_t* ctx,
                  sismap_t* s,
                  parser *p);

void wave_tb_init_p(tb_t *ctx,
                  sismap_t *s,
                  Parameters *p);

void wave_tb_free(tb_t* ctx);

void wave_tb_save_lastshot(sismap_t* s,
                           shot_t *shot,
                           float* u0,
                           float * u1);

void wave_tb_save_lastshot_1st(sismap_t* s,
                           shot_t *shot,
                           float* u0);

void wave_tb_data_init(tb_data_t * data,
                       tb_t *tb,
                       sismap_t *s,
                       const int nb_thread_groups,
                       const int shotid,
                       size_t groupsize);

void wave_tb_data_set_src(tb_data_t * data,
                          sismap_t *s,
                          const unsigned int src_idx,
                          float * source);

void wave_tb_data_unset_src(tb_data_t * data);

void wave_tb_data_set_rcv(tb_data_t * data,
                          sismap_t *s,
                          float* sismos);

void wave_tb_data_unset_rcv(tb_data_t * data);

void wave_tb_data_set_wave(tb_data_t * data,
                           sismap_t *s);
void wave_tb_data_unset_wave(tb_data_t * data);

void wave_tb_data_dump_wave(tb_data_t *data,
                            sismap_t* s);

void wave_tb_data_info(tb_data_t* data);

void wave_tb_data_free(tb_data_t * data,
                       const int nb_thread_groups);

void wave_tb_forward(tb_t* ctx,
                     tb_data_t* data,
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict v0,
                     const float * restrict roc2);

void wave_tb_backward(tb_t* ctx,
                      tb_data_t* data,
                      tb_timer_t* timer,
                      float * restrict u0,
                      float * restrict v0,
                      const float * restrict roc2);

void wave_tb_forward_1st(tb_t* ctx,
                     tb_data_t* data,
					 Parameters *p,
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict vx,
                     float * restrict vy,
                     float * restrict vz,
                     const float * restrict roc2,
					 const float *restrict inv_rho);

void wave_tb_backward_1st(tb_t* ctx,
                     tb_data_t* data,
					 Parameters *p,
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict vx,
                     float * restrict vy,
                     float * restrict vz,
                     const float * restrict roc2,
					 const float *restrict inv_rho);

#endif
