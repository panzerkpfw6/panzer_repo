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

// time measure
#include <sys/time.h>
#include <time.h>

double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
	return (double)clock()/CLOCKS_PER_SEC;
}

////////////////////////
///// Paste from stencil-dev. A bit different notattion for TB.
const Myint   FDM_O1_8_2_LSTENCIL = 4 ;
const Myfloat FDM_O1_8_2_A1      = 1225./1024. ;
const Myfloat FDM_O1_8_2_A2       = -245./3072. ;
const Myfloat FDM_O1_8_2_A3       = 49./5120. ;
const Myfloat FDM_O1_8_2_A4       = -5./7168. ;
const Myint   NB_OP_O2_8          = 12 ;

const Myfloat FDM_O2_8_2_A0       = -205./72. ;
const Myfloat FDM_O2_8_2_A1       = 8./5. ;
const Myfloat FDM_O2_8_2_A2       = -1/5. ;
const Myfloat FDM_O2_8_2_A3       = 8./315. ;
const Myfloat FDM_O2_8_2_A4       = -1/560. ;


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
    float dz, dy, dx; // spacing per dimension in z, y, and x

    // cpu binding masks
    cpu_set_t **bind_masks;
    int setsize;
    int use_manual_cpu_bind;

    // for separate stride-1 functions
    clu_func_t clu_func;

    int num_wf; // number of wavefront updats per iteration
    float dt;
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

    hFloat *  U1, *  U2, *  U3,*  U4, *  source;
    const float * U5;
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


Parameters *gp; //@KADIR global parameter within a node
real_t* recv_rec; //@KADIR array for receiver recording
size_t* irecv_rec; //@KADIR indez into the recv_rec
size_t isrc_exc; //@KADIR number of source ezcitations performed so far

//from diamond_ts.c
#define ST_BUSY (0)
#define ST_NOT_BUSY (1)
typedef struct{
    volatile  int *t_pos;
    int *state;
} Diam_Sched_State;
volatile int *avail_list;
volatile uint64_t head, tail;
int diam_width;
Diam_Sched_State st;
int y_len_l, y_len_r;
int t_len;
int mpi_size;
#define F2H(f) (f)
#define H2F(h) (h)
#define T_POS_L(y) (st.t_pos[(((y)+(y_len_l))%(y_len_l))])
#define T_POS_R(y) (st.t_pos[(((y)+(y_len_r))%(y_len_r))])
Parameters *gp; //@KADIR global parameter within a node
static inline void update_state(int y_coord, Parameters *p){
    //@KADIR EXECUTED IN DIAMOND
    int sh;
    st.t_pos[y_coord]++; // advance the current tile in time
    if(p->is_last != 1) {
        sh = ((st.t_pos[y_coord]%2 == 0) ? 1 : -1);// define the dependency direction
        // add the current tile to the ready queue if its dependency is satisfied
        if( (T_POS_L(y_coord+sh) >= st.t_pos[y_coord]) & (st.t_pos[y_coord] < t_len) )
        {
            avail_list[head%y_len_r] = y_coord;
            head++;
        }
        // add the dependent tile to the ready queue if its other dependency is satisfied
        if( (T_POS_L(y_coord-sh) == st.t_pos[y_coord]) & (T_POS_L(y_coord-sh) < t_len) )
        {
            avail_list[head%y_len_r] = (y_coord - sh + y_len_l)%y_len_l; // add the dependent neighbor to the list if the dependency is satisfied
            head++;
        }

    } else { // last process (and single process case)

        if(st.t_pos[y_coord]%2 == 0){ // right row case
            // add the current diamond to the ready queue if dependencies are satisfied
            if(st.t_pos[y_coord] < t_len){
                // if left-half diamond, no dependencies. Add to the list
                if(y_coord == y_len_l-1){
                    avail_list[head%y_len_r] = y_coord;
                    head++;
                } else if(T_POS_R(y_coord+1) >= st.t_pos[y_coord]) {
                    //the reset have the same circular dependence (ezcept the right-half diamond) if:
                    // 1) the current tile did not reach the end of the temporal dimension
                    // 2) the right neighbor is at least at the same time step
                    avail_list[head%y_len_r] = y_coord;
                    head++;
                }
            } // check validity in range of temporal dimension

            // add the dependent diamond to the ready queue if other dependencies are satisfied:
            if (T_POS_R(y_coord-1) < t_len){
                // add the right-half diamond automatically when the left most diamond is updated
                if(y_coord == 0){ // no dependencies. Add to the list
                    st.t_pos[y_len_r-1]++; // advance the right-half diamond in time
                    avail_list[head%y_len_r] = y_len_r-1;
                    head++;
                }
                else if(T_POS_R(y_coord-1) == st.t_pos[y_coord]) {
                    // 1) the neighbor did not reach the end of the temporal dimension
                    // 2) the left neighbor is at the same time step
                    // 3) is not the right-half diamond
                    avail_list[head%y_len_r] = (y_coord - 1 + y_len_r)%y_len_r; // add the dependent neighbor to the list if the dependency is satisfied
                    head++;
                }
            } // check validity in temporal dimension
        } //end right row case

        else if(st.t_pos[y_coord]%2 == 1){ // left row
            // add the current diamond to the ready queue if dependencies are satisfied:
            if( (T_POS_R(y_coord-1) >= st.t_pos[y_coord]) && (st.t_pos[y_coord] < t_len)  && (y_coord != y_len_r-1) ) {
                // 1) the left neighbor is at least at the same time step
                // 2) the current diamond did not reach the end of the temporal dimension
                // 3) is not the right-half diamond
                avail_list[head%y_len_r] = y_coord;
                head++;
            }

            // add the dependent diamond to the ready queue if other dependencies are satisfied:
            if( (T_POS_R(y_coord+1) == st.t_pos[y_coord]) && (T_POS_R(y_coord+1) < t_len) && (y_coord != y_len_l-1) ) {
                // 1) the right neighbor is at the same time step
                // 2) the neighbor did not reach the end of the temporal dimension
                // 3) is not the right most diamond in space
                avail_list[head%y_len_r] = (y_coord + 1 + y_len_r)%y_len_r; // add the dependent neighbor to the list if the dependency is satisfied
                head++;
            }
        } // end left row case

    } //end is_last process case
}

// Function definitions below
void femwd_iso_ref_2nd( const int shape[3], const int zb, const int yb_r0, const int xb, const int ze, const int ye_r0, const int xe,
                        const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13,
                        hFloat *  p21, hFloat *  p22, hFloat *  p23,const hFloat *  roc2,
                        float * dampx,float * dampy,float * dampz,
                        int t_dim, int b_inc, int e_inc, int NHALO, int tb, int te, stencil_ctx stencil_ctx, int mtid)
{
#pragma omp parallel shared(shape, stencil_ctx, roc2, coef, mtid, tb, te, t_dim, NHALO,recv_rec,irecv_rec) \
firstprivate(b_inc, e_inc) \
num_threads(stencil_ctx.thread_group_size)
    {
        int lstencil=NHALO;// @pavel  allocate variable lstencil
        //    KA_SHOW(FUNC_NAME);
        //    printf("%p\n",  TEMPLATE(femwd, FUNC_NAME));
        int tgs, nwf, th_nwf, tid, gtid, xi, yb, ye, ib, ie, kt, t, ix, iy, iz, q, r, err;
        double t_start;

        const int nny =shape[1];
        const int nnz =shape[0];
        const unsigned long nnzy = 1UL * nnz * nny;
        const unsigned long nnyz = nnzy;

        uint64_t  ln_domain = ((uint64_t) 1)* shape[0]*shape[1]*shape[2];

        tgs = stencil_ctx.thread_group_size;
        nwf = stencil_ctx.num_wf;

        tid = 0;
        gtid = 0;
        ////////////////////////////////
        int nx,ny,nz;
        nz=shape[0]-2*lstencil;
        ny=shape[1]-2*lstencil;
        nx=shape[2]-2*lstencil;

        // load damping parameters
        int NDAMP;
        NDAMP=20; //
        float alpha = 0.2;
        float tabdamp[NDAMP];
        float dampx[nx+2*lstencil];
        float dampy[ny+2*lstencil];
        float dampz[nz+2*lstencil];

        for (int i = 1; i <= NDAMP; i++) {
            tabdamp[NDAMP-i] = exp( -alpha * (1.0*i/NDAMP) * (1.0*i/NDAMP) );
        }

        for (int i = 0; i < 2*lstencil + nx; i++) { //pavel modification (loop)
            dampx[i] = 1.0;
        }
        for (int i = 0; i < 2*lstencil + ny; i++) { //pavel modification (loop)
            dampy[i] = 1.0;
        }
        for (int i = 0; i < 2*lstencil + nz; i++) { //pavel modification (loop)
            dampz[i] = 1.0;
        }
//
        for (int i = 0; i < NDAMP; i++) {
            dampx[i] = tabdamp[i];
            dampy[i] = tabdamp[i];
//        dampz[lstencil+i] = tabdamp[i];

            // original version
            dampx[nx+lstencil-1-i] = tabdamp[i];
            dampy[ny+lstencil-1-i] = tabdamp[i];
            dampz[nz+lstencil-1-i] = tabdamp[i];
            // pavel modification
//        dampx[nx+2*lstencil-1-i] = tabdamp[i];
//        dampy[ny+2*lstencil-1-i] = tabdamp[i];
//        dampz[nz+2*lstencil-1-i] = tabdamp[i];
        }
        // end of loading damping parameters


#if defined(_OPENMP)
        tid = omp_get_thread_num();
		gtid = tid + mtid * tgs;
#endif

        if(stencil_ctx.use_manual_cpu_bind == 1){
            err = sched_setaffinity(0, stencil_ctx.setsize, stencil_ctx.bind_masks[mtid*tgs+tid]);
            if(err==-1) printf("WARNING: Could not set CPU Affinity\n");
        }

        hFloat *  u1 = p11;
        hFloat *  u2 = p12;
        hFloat *  u3 = p13;
        hFloat *  v1 = p21;
        hFloat *  v2 = p22;
        hFloat *  v3 = p23;


        int th_z = stencil_ctx.th_z;
        int th_y = stencil_ctx.th_y;
        int th_x = stencil_ctx.th_x;

        // tid = tid_x*(th_z*th_y) + tid_y*th_z + tid_z
        int tid_z = tid%th_z;
        int tid_y = tid/th_z;
        int tid_x = tid/(th_z*th_y);

        int yb_r = yb_r0;
        int ye_r = ye_r0;

        if(stencil_ctx.th_y>1 ){
            if(b_inc !=0 && e_inc!=0){ // split only at full diamonds
                if (tid_y%2 == 0){ // left thread
                    ye_r = (yb_r + ye_r)/2;
                    e_inc = 0;
                } else{
                    yb_r = (yb_r + ye_r)/2;
                    b_inc = 0;
                }
            }else{// use the y-threads along x-axis make sure to use sufficient number of frontlines
                th_x *= th_y;
                tid_x = tid/th_z;
                if (nwf < th_x) nwf = th_x;
            }
        }

        int nbz = (ze-zb)/th_z;
        q = (int)((ze-zb)/th_z);
        r = (ze-zb)%th_z;
        if(tid_z < r) {
            ib = zb + tid_z * (q+1);
            ie = ib + (q+1);
        }else {
            ib = zb + r * (q+1) + (tid_z - r) * q;
            ie = ib + q;
        }

        th_nwf = nwf/th_x;

        int printed = 0; //@KADIR

        int end=0;

        for(xi=xb; xi<xe; xi+=nwf) { // wavefront loop (x direction)

            if(xe-xi <= nwf){
                nwf = xe-xi;
                end =1;
            }

            yb = yb_r;
            ye = ye_r;

            kt = xi;
            int kte=kt+nwf;


            for(t=tb; t< te; t++){ // Diamond blocking in time
                int mod = (t)%2;
                if (mod) {
                    u1 = p11 ;
                    v1 = p21 ;
                } else {
                    u1 = p21 ;
                    v1 = p11 ;
                }


//#pragma omp barrier
                const Myfloat inv_dz2 = 1. / (stencil_ctx.dz*stencil_ctx.dz) ;
                const Myfloat inv_dx2 = 1. / (stencil_ctx.dx*stencil_ctx.dx) ;
                const Myfloat inv_dy2 = 1. / (stencil_ctx.dy*stencil_ctx.dy);
                const Myfloat coef = stencil_ctx.dt;
                for(ix=kt; ix<kte; ix++){
                    if( ((ix)/th_nwf)%th_x == tid_x ) {
                        for(iy=yb; iy<ye; iy++) {
                            float *  pr0_v = &(u1[ix*nnyz+iy*nnz]);
                            float *  vx0_v = &(v1[ix*nnyz+iy*nnz]);
                            const float *  coef0_v = &(roc2[ix*nnyz+iy*nnz]);

#pragma ivdep
                            for(iz=ib; iz<ie; iz++) {
                                Myfloat d2_pr_x = (  FDM_O2_8_2_A0 *  pr0_v[ 0*nnyz + iz]
                                                     + FDM_O2_8_2_A1 * (pr0_v[-1*nnyz + iz] + pr0_v[ 1*nnyz + iz])
                                                     + FDM_O2_8_2_A2 * (pr0_v[-2*nnyz + iz] + pr0_v[ 2*nnyz + iz])
                                                     + FDM_O2_8_2_A3 * (pr0_v[-3*nnyz + iz] + pr0_v[ 3*nnyz + iz])
                                                     + FDM_O2_8_2_A4 * (pr0_v[-4*nnyz + iz] + pr0_v[ 4*nnyz + iz])) * inv_dx2 ;


                                Myfloat d2_pr_y = (  FDM_O2_8_2_A0 *  pr0_v[ 0*nnz + iz]
                                                     + FDM_O2_8_2_A1 * (pr0_v[-1*nnz + iz] + pr0_v[ 1*nnz + iz])
                                                     + FDM_O2_8_2_A2 * (pr0_v[-2*nnz + iz] + pr0_v[ 2*nnz + iz])
                                                     + FDM_O2_8_2_A3 * (pr0_v[-3*nnz + iz] + pr0_v[ 3*nnz + iz])
                                                     + FDM_O2_8_2_A4 * (pr0_v[-4*nnz + iz] + pr0_v[ 4*nnz + iz])) * inv_dy2 ;

                                Myfloat d2_pr_z = (  FDM_O2_8_2_A0 *  pr0_v[ 0 + iz]
                                                     + FDM_O2_8_2_A1 * (pr0_v[-1 + iz] + pr0_v[ 1 + iz])
                                                     + FDM_O2_8_2_A2 * (pr0_v[-2 + iz] + pr0_v[ 2 + iz])
                                                     + FDM_O2_8_2_A3 * (pr0_v[-3 + iz] + pr0_v[ 3 + iz])
                                                     + FDM_O2_8_2_A4 * (pr0_v[-4 + iz] + pr0_v[ 4 + iz])) * inv_dz2 ;

                                vx0_v[iz] = coef0_v[iz] * stencil_ctx.dt * (d2_pr_x + d2_pr_y + d2_pr_z) - vx0_v[iz] + (Myfloat)(2.0) * pr0_v[iz];
                                // ABCs
                                vx0_v[iz] = dampx[ix+lstencil] * vx0_v[iz] + (1 -dampx[ix+lstencil]) * pr0_v[iz];
                                vx0_v[iz] = dampy[iy+lstencil] * vx0_v[iz] + (1 -dampy[iy+lstencil]) * pr0_v[iz];
                                vx0_v[iz] = dampz[iz+lstencil] * vx0_v[iz] + (1 -dampz[iz+lstencil]) * pr0_v[iz];
                            }

                        }


                        if( (gp->source_point_enabled==1)
                            && (gp->lsource_pt[2] >= ib ) //@KADIR
                            && (gp->lsource_pt[2] <  ie ) //@KADIR
                            && (gp->lsource_pt[1] >= yb ) //@KADIR
                            && (gp->lsource_pt[1] <  ye ) //@KADIR
                            && (gp->lsource_pt[0] == ix ) //@KADIR
                                )
                        {
                            V1_xyz(gp->lsource_pt[0],gp->lsource_pt[1],gp->lsource_pt[2]) += gp->src_exc_coef[isrc_exc];
                            isrc_exc++;
                        }

                    }
                }
                // Update block size in Y
                if(t< t_dim){ // lower half of the diamond
                    yb += -b_inc;
                    ye += e_inc;
                }else{ // upper half of the diamond
                    yb += b_inc;
                    ye += -e_inc;
                }
                kte=max(kte-NHALO,xb);
                if (end==1) kte =xe;
                kt=max(kt-NHALO,xb);
                t_start = get_wall_time();
#pragma omp barrier
                stencil_ctx.t_wait[gtid] += get_wall_time() - t_start;
            } // diamond blocking in time (time loop)
        } // wavefront loop
    } // parallel region
}

void intra_diamond_mwd_comp_std(Parameters *p, int yb_r, int ye_r, int b_inc, int e_inc, int tb, int te, int tid){
    //@KADIR1 EXECUTED IN DIAMOND
    int t, x, xb[32], xe[32];
    int xb0,xe0;
    int yb, ye;
    int time_len = te-tb;
    double t1, t2, t3;
//    int lstencil; //pavel edition

    // wavefront prologue
    t1 = get_wall_time();
#if 0
    yb = yb_r;
	ye = ye_r;

	for(t=tb; t< te-1; t++){
		xb[t] = p->stencil.r;
		xe[t] = p->stencil.r*(time_len-(t-tb));
	}

	stat_sched_iso_ref_v2(p->ldomain_shape, p->stencil.r, yb, xb,
                			  p->lstencil_shape[0]+p->stencil.r, ye, xe, p->coef, p->U1,p->U1, p->U1, p->U2, p->U3,p->U4, p->U5, p->t_dim, b_inc, e_inc, p->stencil.r, tb, te-1, p->stencil_ctx, tid);
#endif
    t2 = get_wall_time();
    // main wavefront loop
    yb = yb_r;
    ye = ye_r;
    //xb0 = (te-tb)*p->stencil.r;
    xb0 = p->stencil.r;
    xe0 = p->ldomain_shape[2]-p->stencil.r;
//    lstencil=(p->ldomain_shape[2]-p->lstencil_shape[2])/2; //pavel edition
    p->stencil.mwd_func(p->ldomain_shape, p->stencil.r, yb, xb0,
                        p->lstencil_shape[0]+p->stencil.r, ye, xe0, p->coef, p->U1,p->U1, p->U1, p->U2, p->U3,p->U4, p->U5,
                        p->dampx,p->dampy,p->dampz,
                        p->t_dim, b_inc, e_inc, p->stencil.r, tb, te, p->stencil_ctx, tid);
    t3 = get_wall_time();

    // wavefront epilogue
#if 0
    yb = yb_r;
	ye = ye_r;

	if(tb < p->t_dim) { // lower half of the diamond
    yb -= b_inc;
    ye += e_inc;
  } else { // upper half of the diamond
    yb += b_inc;
    ye -= e_inc;
  }


	for(t=tb+1; t< te; t++){
		xe[t] = p->ldomain_shape[2]-p->stencil.r;
		xb[t] = p->ldomain_shape[2]-p->stencil.r - (t-tb)*p->stencil.r;
	}

	stat_sched_iso_ref_v2(p->ldomain_shape, p->stencil.r, yb, xb,
                			  p->lstencil_shape[0]+p->stencil.r, ye, xe, p->coef, p->U1,p->U1, p->U1, p->U2, p->U3,p->U4, p->U5, p->t_dim, b_inc, e_inc, p->stencil.r, tb+1, te, p->stencil_ctx, tid);
#endif
    p->stencil_ctx.t_wf_prologue[tid] += t2-t1;
    p->stencil_ctx.t_wf_main[tid]     += t3-t2;
    p->stencil_ctx.t_wf_epilogue[tid] += get_wall_time() - t3;
}

void dynamic_intra_diamond_ts_combined(Parameters *p) {
    //@5
    MSG("inside dynamic_intra_diamond_ts_combined");
    int t_dim = p->t_dim;
    diam_width = (t_dim+1) * 2 *p->stencil.r;
    if(p->stencil.type == REGULAR){
        t_len = 2*( (p->nt-2)/((t_dim+1)*2) ) - 1;
    } else if(p->stencil.type == SOLAR){
        t_len = 2*( (p->nt)/((t_dim+1)*2) ) - 1;
    }
    int num_thread_groups = get_ntg(*p);

    y_len_l = p->lstencil_shape[1] / (diam_width);
    y_len_r = y_len_l;
    if(p->is_last == 1) y_len_r++;

    int i, y, t;
    double t1,t2,t3,t4;
    int yb,ye;
    double db_t;

    // allocate scheduling variables
    st.t_pos = (int*) malloc(y_len_r*sizeof(int));
    st.state = (int*) malloc(y_len_r*sizeof(int));
    avail_list = (int*) malloc(y_len_r*sizeof(int));
    head=y_len_r;
    tail=0;
    // initialixe scheduling variables
    for(i=0; i<y_len_r; i++){
        st.t_pos[i] = 0;
        st.state[i] = ST_NOT_BUSY;
    }
#if defined(_OPENMP)
    omp_set_nested(1);
#endif
    MSG("689");
    gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] = F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]) + gp->src_exc_coef[isrc_exc]);//@KADIR
    isrc_exc++;

    // Prologue
    t1 = get_wall_time();
    if(p->in_auto_tuning == 0){
        //dynamic_intra_diamond_prologue(p);
        //@4.1
        if(p->stencil.type == REGULAR){
            //dynamic_intra_diamond_prologue_std(p);
            //@3.1
            // compute all the trapexoids
            int i, yb, ye;
            int ntg = get_ntg(*p);
#pragma omp parallel num_threads(ntg)
            {
                int b_inc = p->stencil.r;
                int e_inc = p->stencil.r;
                int tid = 0;
#if defined(_OPENMP)
                tid = omp_get_thread_num();
#endif
                MSG("712");
#pragma omp for schedule(dynamic) private(i,yb,ye)
                for(i=0; i<y_len_l; i++){
                    yb = p->stencil.r + i*diam_width;
                    ye = yb + diam_width;
                    //intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, p->t_dim, p->t_dim*2+1, tid);
                    //(Parameters *p, int yb, int ye, int b_inc, int e_inc, int tb, int te, int tid)
                    //@2
                    {
                        int tb = p->t_dim;
                        int te = p->t_dim*2+1;
                        if(p->stencil.type == REGULAR){
                            intra_diamond_mwd_comp_std(p, yb, ye, b_inc, e_inc, tb, te, tid);
                        }
                    }
                }
            }
        }
    }
    t2 = get_wall_time();

    // main loop
    //dynamic_intra_diamond_main_loop(p);
    {
        //@4.3.2
        int not_complete, th_y_coord, i;
        uint64_t il;
        int num_thread_groups = get_ntg(*p);
        uint64_t diam_size = y_len_l*(t_len-1)/2 + y_len_r*((t_len-1)/2 +1);
        int tid;
        double t1;

        int idz=0;

        if(p->in_auto_tuning == 0) {
            for(i=0; i<y_len_r; i++){
                avail_list[i] = i;
            }
        } else { // diversify the startup for shorter autotuning
            for(i=0; i<y_len_r; i++){
                if(i%2==0){
                    avail_list[i] = idz++;
                }
            }
            for(i=0; i<y_len_r; i++){
                if(i%2==1){
                    avail_list[i] = idz++;
                }
            }
        }

#pragma omp parallel num_threads(num_thread_groups) shared(head, tail) private(tid)
        {
            tid = 0;
#if defined(_OPENMP)
            tid = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic) private(il, th_y_coord, not_complete)//shared(head,tail)
            for (il=0; il<diam_size; il++){
                not_complete = 1;
                th_y_coord = -1;
                while(not_complete) {
                    t1 = get_wall_time();
                    while(head-tail<1); // spin-wait for available tasks
                    p->stencil_ctx.t_group_wait[tid] += (get_wall_time() - t1);
#pragma omp critical// (consumer)
                    {
#pragma omp flush (head, tail)
                        if(head-tail>0){ // make sure there is still available work
                            th_y_coord = avail_list[tail%y_len_r]; //acquire task
                            tail++;
                        }
                    }
                    if(th_y_coord>=0){
                        //intra_diamond_resolve(p, th_y_coord, tid);
                        //(Parameters *p, int y_coord, int tid)
                        {
                            int y_coord = th_y_coord;
                            //@4.3.1
                            int t_coord = st.t_pos[y_coord];
                            double t1, t2;
                            //intra_diamond_comp_using_location(p, y_coord, tid, t_coord);
                            //(Parameters *p, int y_coord, int tid, int t_coord)
                            {
                                //@3.3
                                int yb, ye, b_inc, e_inc;
                                if(p->stencil.type == REGULAR){
                                    //intra_diamond_get_info_std(p, y_coord, tid, t_coord, &yb, &ye, &b_inc, &e_inc);
                                    //(Parameters *p, int y_coord, int tid, int t_coord, int *yb, int *ye, int *b_inc, int *e_inc)
                                    {
                                        double diam_size;
                                        if( (p->is_last == 1) && (y_coord == y_len_l-1) && (t_coord%2 == 0) ){ // right most process & left-half diamond
                                            // left half computations
                                            yb = p->stencil.r + p->lstencil_shape[1] - p->stencil.r;
                                            ye = yb + p->stencil.r;
                                            b_inc = p->stencil.r;
                                            e_inc = 0;
                                            diam_size = 0.5;
                                        } else if( (p->is_last == 1) && (y_coord == y_len_r-1) && (t_coord%2 == 0) ){ // right most process & right-half diamond
                                            // right half computations
                                            b_inc = 0;
                                            e_inc = p->stencil.r;
                                            if(p->t.shape[1] > 1)
                                                yb = p->stencil.r + p->lstencil_shape[1] + 2*p->stencil.r;
                                            else // serial code case
                                                yb = p->stencil.r;
                                            ye = yb + p->stencil.r;
                                            diam_size = 0.5;
                                        } else{ // full diamond computation
                                            if(t_coord%2 == 0)// row shifted to the right
                                                yb = p->stencil.r + diam_width - p->stencil.r + y_coord*diam_width;
                                            else // row shifted to the left
                                                yb = p->stencil.r + diam_width/2 - p->stencil.r+ y_coord*diam_width;
                                            ye = yb + 2*p->stencil.r;
                                            b_inc = p->stencil.r;
                                            e_inc = p->stencil.r;
                                            diam_size = 1.0;
                                        }
                                        p->stencil_ctx.wf_num_resolved_diamonds[tid] += diam_size;
                                    }

                                    //intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, 0, p->t_dim*2+1, tid);
                                    //(Parameters *p, int yb, int ye, int b_inc, int e_inc, int tb, int te, int tid)
                                    {
                                        int tb = 0;
                                        int te = p->t_dim*2+1;
                                        //@2
                                        if(p->stencil.type == REGULAR){
                                            intra_diamond_mwd_comp_std(p, yb, ye, b_inc, e_inc, tb, te, tid);
                                        }
                                    }

                                }
                            }
                            p->stencil_ctx.t_wf_comm[tid] += t2-t1;
                        }
#pragma omp critical// (producer)
                        {
#pragma omp flush (head)
                            update_state(th_y_coord, p);
                        }
                        not_complete = 0;
                    }
                }
            }
        }
    }
    t3 = get_wall_time();
    // Epilogue
    if(p->in_auto_tuning == 0){
        //dynamic_intra_diamond_epilogue(p);
        //@4
        //dynamic_intra_diamond_epilogue_std(p);
        //@3
        int yb, ye, i;
        int ntg = get_ntg(*p);
#pragma omp parallel num_threads(ntg)
        {
            int b_inc = p->stencil.r;
            int e_inc = p->stencil.r;
            int yb_r = p->stencil.r + diam_width/2 - p->stencil.r;
            int ye_r = yb_r + 2*p->stencil.r;
            int tid = 0;
#if defined(_OPENMP)
            tid = omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic) private(i,yb,ye)
            for(i=0; i<y_len_l; i++){
                yb = yb_r + i*diam_width;
                ye = ye_r + i*diam_width;
                //intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, 0, p->t_dim+1, tid);
                {
                    int tb = 0;
                    int te = p->t_dim+1;
                    //@2
                    intra_diamond_mwd_comp_std(p, yb, ye, b_inc, e_inc, tb, te, tid);
                }
            }
        }
    }
    t4 = get_wall_time();

    p->prof.ts_main += (t3-t2);
    p->prof.ts_others += (t2-t1) + (t4-t3);
    // cleanup the state variables
    free((void *) st.t_pos);
    free(st.state);
    free((void *) avail_list);
}

void reset_timers(Profile * p){
    p->compute = 0.;
    p->communicate = 0.;
    p->send_recv = 0.;
    p->wait = 0.;
    p->total = 0.;
    p->others = 0.;
    p->ts_main = 0.;
    p->ts_others = 0.;
}
void reset_wf_timers(Parameters * p){
    int i;
    int num_thread_groups = get_ntg(*p);

    // reset if the wavefront profiling is allocated
    if( (p->wavefront != 0) && (p->target_ts == 2) ) {

        for(i=0; i<p->num_threads; i++) p->stencil_ctx.t_wait[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_main[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_comm[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_prologue[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_epilogue[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.wf_num_resolved_diamonds[i] = 0.0;
        for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_group_wait[i] = 0.0;
    }
}
void cpu_bind_init(Parameters *p){
    if(p->stencil_ctx.use_manual_cpu_bind==0)
        return;
    //printf("%s %d: using sched_affinity()\n", __FILE__, __LINE__);

    // Source for finding number of CPUs: https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
    cpu_set_t *mask;
    int ncpus;
    for ( ncpus = sizeof(cpu_set_t)/8; ncpus < 16*1024; ncpus <<= 1 ) {
        mask = CPU_ALLOC( ncpus );
        if ( !mask ) break;
        const size_t size = CPU_ALLOC_SIZE( ncpus );
        CPU_ZERO_S( size, mask );
        const int err = sched_getaffinity( 0, size, mask );
        if ( !err ) break;
        CPU_FREE( mask );
        mask = NULL;
        //if ( errno != EINVAL )  break; REALLY FIXME KADIR
    }
    if ( !mask )
        printf("Warning: Failed to obtain process affinity mask. Thread affinitixation is disabled.\n");

    p->stencil_ctx.setsize = CPU_ALLOC_SIZE(ncpus);
    p->stencil_ctx.bind_masks = (cpu_set_t**) malloc(p->num_threads*sizeof(cpu_set_t*));

    int i, ib, idz=0;
    ib=0;
#if __MIC__
    ib = 1;
#endif
    for(i=ib; i<p->num_threads*p->th_stride/p->th_block+ib;i++){
        if((i-ib)%p->th_stride < p->th_block){
            p->stencil_ctx.bind_masks[idz] = CPU_ALLOC( ncpus );
            CPU_ZERO_S(p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[idz]);
            CPU_SET_S(i,p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[idz]);
            idz++;
        }
    }

    int *phys_cpu = (int*) malloc(p->num_threads*sizeof(int));

    omp_set_nested(1);
    // Set the affinity to reduce the cost of first run
    int num_thread_groups = get_ntg(*p);
#pragma omp parallel num_threads(num_thread_groups) //PROC_BIND(spread)
    {
        int mtid = omp_get_thread_num();
#pragma omp parallel shared(mtid)  num_threads(p->stencil_ctx.thread_group_size) //PROC_BIND(master)
        {
            int tid = omp_get_thread_num();
            int gtid = tid + mtid * p->stencil_ctx.thread_group_size;
            int err = sched_setaffinity(0, p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[gtid]);
            if(err==-1) printf("WARNING: Could not set CPU Affinity of thread:%d error:%d\n", gtid, err);
            phys_cpu[gtid] = sched_getcpu();
        }
    }
#ifdef _GNU_SOURCE
    //printf("__________________________GNU SOURCE IS DEFINED\n");
#else
    printf("__________________________GNU SOURCE IS NOT NOT NOT DEFINED\n");
#endif
    printf("Threads binding (tid->OS tid):");
    for(i=0;i<p->num_threads;i++){
        printf(" %d->%d", i, phys_cpu[i]);
    }
    printf("\n");

    free(phys_cpu);
}

////////////////////////