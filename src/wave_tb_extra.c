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


////////////////////////