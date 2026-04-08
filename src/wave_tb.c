#define _GNU_SOURCE
#include <sched.h>
#include <sys/sysinfo.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include "stencil/macros.h"
#include "stencil/sismap.h"
#include "stencil/parser.h"
#include "stencil/wave_tb.h"
//#include "stencil/wave_tb_extra.h"
#include "stencil/wtime.h"
// time measure
#include <sys/time.h>
#include <time.h>

volatile int wave_tb_head;
volatile int wave_tb_tail;

#define NDAMP 20

#if 1
#define FUNC_BODY() {                                                       \
ux[i] = 2.0f * vx[i] - ux[i]                                                \
      + rx[i] * (coef0 * vx[i] + coefx[1] * (vx[i+1     ] + vx[i-1     ])   \
                               + coefy[1] * (vx[i+nnx   ] + vx[i-nnx   ])   \
                               + coefz[1] * (vx[i+nnxy  ] + vx[i-nnxy  ])   \
                               + coefx[2] * (vx[i+2     ] + vx[i-2     ])   \
                               + coefy[2] * (vx[i+2*nnx ] + vx[i-2*nnx ])   \
                               + coefz[2] * (vx[i+2*nnxy] + vx[i-2*nnxy])   \
                               + coefx[3] * (vx[i+3     ] + vx[i-3     ])   \
                               + coefy[3] * (vx[i+3*nnx ] + vx[i-3*nnx ])   \
                               + coefz[3] * (vx[i+3*nnxy] + vx[i-3*nnxy])   \
                               + coefx[4] * (vx[i+4     ] + vx[i-4     ])   \
                               + coefy[4] * (vx[i+4*nnx ] + vx[i-4*nnx ])   \
                               + coefz[4] * (vx[i+4*nnxy] + vx[i-4*nnxy])); \
ux[i] = dampx[i] * ux[i] + (1 - dampx[i]) * vx[i];                          \
ux[i] = dampy[j] * ux[i] + (1 - dampy[j]) * vx[i];                          \
ux[i] = dampz[k] * ux[i] + (1 - dampz[k]) * vx[i];                          \
}
#else
#define FUNC_BODY() {                                                       \
ux[i] = 2.0f * vx[i] - ux[i]                                                \
      + rx[i] * (coef0 * vx[i] + coefx[1] * (vx[i+1     ] + vx[i-1     ])   \
                               + coefy[1] * (vx[i+nnx   ] + vx[i-nnx   ])   \
                               + coefz[1] * (vx[i+nnxy  ] + vx[i-nnxy  ])   \
                               + coefx[2] * (vx[i+2     ] + vx[i-2     ])   \
                               + coefy[2] * (vx[i+2*nnx ] + vx[i-2*nnx ])   \
                               + coefz[2] * (vx[i+2*nnxy] + vx[i-2*nnxy])   \
                               + coefx[3] * (vx[i+3     ] + vx[i-3     ])   \
                               + coefy[3] * (vx[i+3*nnx ] + vx[i-3*nnx ])   \
                               + coefz[3] * (vx[i+3*nnxy] + vx[i-3*nnxy])   \
                               + coefx[4] * (vx[i+4     ] + vx[i-4     ])   \
                               + coefy[4] * (vx[i+4*nnx ] + vx[i-4*nnx ])   \
                               + coefz[4] * (vx[i+4*nnxy] + vx[i-4*nnxy])); \
}
#endif

#define FUNC_BODY_1st_ord_Psweep()  {                            \
ux[i] = ux[i] + rx[i] * (coefx[0]/data->dx * (vx_x[i] - vx_x[i-1])   \
                       + coefy[0]/data->dy * (vy_x[i] - vy_x[i-nnx])   \
                       + coefz[0]/data->dz * (vz_x[i] - vz_x[i-nnxy])   \
                       + coefx[1]/data->dx * (vx_x[i+1]-vx_x[i-2])   \
                       + coefy[1]/data->dy * (vy_x[i+1*nnx]-vy_x[i-2*nnx])   \
                       + coefz[1]/data->dz * (vz_x[i+1*nnxy]-vz_x[i-2*nnxy])   \
                       + coefx[2]/data->dx * (vx_x[i+2]-  vx_x[i-3])   \
                       + coefy[2]/data->dy * (vy_x[i+2*nnx] - vy_x[i-3*nnx])   \
                       + coefz[2]/data->dz * (vz_x[i+2*nnxy] - vz_x[i-3*nnxy])   \
                       + coefx[3]/data->dx * (vx_x[i+3]-  vx_x[i-4])   \
                       + coefy[3]/data->dy * (vy_x[i+3*nnx] -vy_x[i-4*nnx])   \
                       + coefz[3]/data->dz * (vz_x[i+3*nnxy] -vz_x[i-4*nnxy])); \
ux[i] = dampx[i] * ux[i];                          \
ux[i] = dampy[j] * ux[i];                          \
ux[i] = dampz[k] * ux[i];                           \
}

#define FUNC_BODY_1st_ord_Vsweep() {     \
vx_x[i] = vx_x[i] \
      + data->dt/data->dx * (     coefx[0] * (ux[i+1] - ux[i])    \
                    + coefy[1] * (ux[i+2] - ux[i-1])   \
                    + coefz[2] * (ux[i+3] - ux[i-2])   \
                    + coefz[3] * (ux[i+4] - ux[i-3]));  \
vy_x[i] = vy_x[i] \
      + data->dt/data->dy * (     coefx[0] * (ux[i+1*nnx] - ux[i])   \
                   +  coefy[1] * (ux[i+2*nnx] - ux[i-1*nnx])   \
                   +  coefz[2] * (ux[i+3*nnx] - ux[i-2*nnx])   \
                   +  coefz[3] * (ux[i+4*nnx] - ux[i-3*nnx]));  \
vz_x[i] = vz_x[i] \
      + data->dt/data->dz * (    coefx[0] * (ux[i+1*nnxy] - ux[i])   \
                   + coefy[1] * (ux[i+2*nnxy] - ux[i-1*nnxy])   \
                   + coefz[2] * (ux[i+3*nnxy] - ux[i-2*nnxy])   \
                   + coefz[3] * (ux[i+4*nnxy] - ux[i-3*nnxy])); \
}

///////////////////////////////////////////////////////////////////////
//
//  wave_tb_extra.h code
//

double get_wall_time() {
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return time.tv_sec + time.tv_nsec / 1000000000.0;
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

#ifndef _OPENMP
#define  omp_get_num_threads() (1)
#define omp_set_nested(a) {}
#define omp_get_thread_num() (0)
#define omp_get_max_threads() (1)
#endif


//////
//struct Stencil {
//    const char *name;
//    int r;
//    int time_order;
//    int nd;
//    enum Stencil_Shapes shape;
//    enum Stencil_Coefficients coeff;
//    enum Stencil_Type type;
//    spt_blk_func_t spt_blk_func;
//    spt_blk_func_t stat_sched_func;
//    mwd_func_t mwd_func;
//};
//struct StencilInfo {
//    const char *name;
//    int r;
//    int time_order;
//    int nd;
//    enum Stencil_Shapes shape;
//    enum Stencil_Coefficients coeff;
//    enum Stencil_Type type;
//};

//// contezt information
//typedef struct{
//    int alignment, verbose, stencil_shape[3];
//    uint64_t n_stencils, ln_domain, ln_stencils;
//    int target_ts, target_kernel;
//    int mpi_rank, mpi_size;
//    int n_tests, nt;
//    int verify;
//    int source_pt[3];
//    int debug;
//    int num_threads;
//    int use_omp_stat_sched;
//    int lstencil_shape[3], ldomain_shape[3], gb[3], ge[3], lsource_pt[3], has_source; //MPI ranks' global indices, and local source locations
//
//    int notuning; //@KADIR returns early from autuning functions
//    int call_combined_function; //@KADIR calls Kadir's function that does everything
//
//    //  int stencil_radius, is_constant_coefficient;
//    //  enum Stencil_Types stencil_type;
//
//    hFloat *  U1, *  U2, *  U3,*  U4, *  source;
//    const float * U5;
//    const float * U6;
////////    hFloat *  rU1, *  rU2, *  rU3,*  rU4, *  rU5;
//
//    // damping ABCs
//    float * dampx;
//    float * dampy;
//    float * dampz;
//    real_t *  coef;
//
//    real_t *  src_exc_coef; //@KADIR: coef used in source ezcitation. length is number of time steps (nt)
//    // Use source instead of src_exc_conf after verification results
//
//    // parameters for internal thread affinity
//    int th_block;
//    int th_stride;
//
//    // Holds the value of cache blocking across Y axis
//    stencil_ctx stencil_ctx;
//
//    // to enable/disable source point update
//    int source_point_enabled;
//
//    int cache_size; // Last level cache usable size in KB for blocking
//
//    // to enable concatenating halo information before communication
//    int halo_concat;
//
//    // Specific data for the diamond method
//    int t_dim, larger_t_dim, is_last, mwd_type;
//    Halo hu[3], hv[3];
//
//    int wavefront;
//    uint64_t idiamond_pro_epi_logue_updates;
//    uint64_t wf_blk_size, wf_larger_blk_size;
//
//    Halo h[3]; // Halo information for z,Y, and x directions
//    mpi_topology t;
//    Profile prof;
//
//    struct Stencil stencil;
//
//    // list of coefficients to be used in stencil operators
//    real_t g_coef[11];
//
//    int array_padding;
//
//    int in_auto_tuning;
//    int orig_thread_group_size; // to distingquish whether thread group size is set by the user
//    tb_data_t * data;
//}Parameters;
//Parameters *gp; //@KADIR global parameter within a node


real_t* recv_rec; //@KADIR array for receiver recording
size_t* irecv_rec; //@KADIR indez into the recv_rec
size_t isrc_exc; // number of source ezcitations performed so far
size_t isrc_exc2; // number of source ezcitations performed so far, 2nd counter for rtm

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

int get_ntg(Parameters p){
    return (int) ceil(1.0*p.num_threads/p.stencil_ctx.thread_group_size);
}

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
                    avail_list[head%y_len_r]=y_coord;
                    head++;
                } else if(T_POS_R(y_coord+1)>=st.t_pos[y_coord]) {
                    //the reset have the same circular dependence (ezcept the right-half diamond) if:
                    // 1) the current tile did not reach the end of the temporal dimension
                    // 2) the right neighbor is at least at the same time step
                    avail_list[head%y_len_r]=y_coord;
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
static inline int get_t_coord(int y_coord) {
    // Ensure st.t_pos is allocated
    if (!st.t_pos) {
        fprintf(stderr, "st.t_pos is NULL\n");
        exit(1);
    }
    // Validate y_coord
    if (y_coord < 0 || y_coord >= y_len_r) {
        fprintf(stderr, "Invalid y_coord=%d, y_len_r=%d\n", y_coord, y_len_r);
        exit(1);
    }
    // Read t_pos[y_coord] thread-safely
    int t_coord;
    #pragma omp critical
    {
        t_coord = st.t_pos[y_coord];
    }
    return t_coord;
}
// Function definitions below
void femwd_iso_ref_1st( const int shape[3], const int zb, const int yb_r0,
                        const int xb, const int ze, const int ye_r0, const int xe,
                    const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13,
                    hFloat *  p21, hFloat *  p22, hFloat *  p23,
					const hFloat * roc2, const hFloat * inv_rho,
                    float * dampx,float * dampy,float * dampz,
                    int t_dim, int b_inc, int e_inc,int NHALO,
                    int tb, int te,int t0,int ifwd,
					stencil_ctx stencil_ctx,int mtid,tb_data_t * data)
{
#pragma omp parallel shared(shape, stencil_ctx, roc2, coef, mtid, tb, te, t_dim, NHALO,recv_rec,irecv_rec) \
firstprivate(b_inc, e_inc) \
num_threads(stencil_ctx.thread_group_size)
    {
//    	MSG("xb=%d,xe=%d",xb,xe);
        int lstencil=NHALO;// @pavel  allocate variable lstencil
        int tgs, nwf, th_nwf, tid, gtid, xi, yb, ye, ib, ie, kt, t,  q, r, err;
        double t_start;

        const int nnx =shape[2];
        const int nny =shape[1];
        const int nnz =shape[0];

        const unsigned long nnzy = 1UL * nnz * nny;
        const unsigned long nnyz = nnzy;
        const int64_t nnxyz=1ULL*nnx * nny * nnz;
        const int64_t nnxy=1ULL*nnx * nny;
        const int64_t nnyz_grid=1ULL*nnx * nny;

        // index notation for velocity array
        const int nnz_v=stencil_ctx.nz;
        const unsigned long nnyz_v=1UL*stencil_ctx.nz*stencil_ctx.ny;

//        MSG("nnx=%d,nny=%d,nnz=%d",nnx,nny,nnz);
//        MSG("stencil. nnx=%d,nny=%d,nnz=%d",stencil_ctx.nx,stencil_ctx.ny,stencil_ctx.nz);
//        MSG("nnxy=%d",nnxy);
//        exit(1);

        tgs = stencil_ctx.thread_group_size;
        nwf = stencil_ctx.num_wf;

        tid = 0;
        gtid = 0;
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
        int iz_=data->rcv_depth; //@pavel
        int end=0;

        // Precompute coefficients with dt for velocity updates
		const Myfloat dt_inv_dx = stencil_ctx.dt / (stencil_ctx.dx);
		const Myfloat dt_inv_dy = stencil_ctx.dt / (stencil_ctx.dy);
		const Myfloat dt_inv_dz = stencil_ctx.dt / (stencil_ctx.dz);

        const Myfloat inv_dx = 1./ (stencil_ctx.dx);
        const Myfloat inv_dy = 1. / (stencil_ctx.dy);
        const Myfloat inv_dz = 1. / (stencil_ctx.dz);

        // Load wavefield and imaging condition /////////////////////
        // RTM variables

        float *restrict ux;
        float *restrict vx;
		float *restrict wx;
		float *restrict imgx;
		float *restrict ilmx;
		float * __restrict v3_v;
		int kte;
		//////////////////////////////////////Backward//////////////////////////////////////
        if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // load fwd wavefield and compute IC
//        	MSG("bwd phase, ifwd=%d",ifwd);
        	yb = yb_r;
        	ye = ye_r;
        	kt = xb;
//        	kte=kt+nwf;
        	kte=xe;
        	hFloat *  v3=p13;
        	v3=	p13;
        	for(int t=tb; t< te; t++){ // Diamond blocking in time
				hFloat* output_buffer = NULL;  //@KADIR
				int mod = (t)%2;
//                MSG("t=%d",t);
//				MSG("bwd, ifwd=%d",ifwd);
				if(mod==0){// compute p from v
//					u1=	p21 ;
					for(int ix=kt; ix<kte; ix++){
						if( ((ix)/th_nwf)%th_x == tid_x ) {
							for(int iy=yb; iy<ye; iy++) {
								unsigned long int index=1ULL*(ix-NHALO)*nnyz_v+(iy-NHALO)*nnz_v;
//								MSG("IMG rec: ix=%d, iy=%d",ix,iy);
								vx=&(v3[1ULL*ix*nnyz+iy*nnz]);
								wx   = &(data->fwd[1ULL * ifwd*nnxyz + 1ULL*ix*nnyz + iy*nnz]);
								imgx = &(data->img[ index ]);
								ilmx = &(data->ilm[ index ]);

//								coef0_v = &(roc2[1ULL*ix*nnyz+iy*nnz]); //// original
//								coef0_v = &(roc2[1ULL*(ix-NHALO)*nnyz_v+(iy-NHALO)*nnz_v]);
//								unsigned long int index=1ULL*ix*nnyz+iy*nnz;
//								imgx = &(data->img[ index ]);
//								ilmx = &(data->ilm[ index ]);

#pragma ivdep
								for(int iz=ib; iz<ie; iz++) {
									imgx[iz] += vx[iz]*wx[iz];
									ilmx[iz] += wx[iz]*wx[iz];
//									imgx[iz] += 1;
//									ilmx[iz] += 2;
								}
							}
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
			} // diamond blocking in time (time loop)
        }
        //////////////////////////////////////
        // calculate wavefield update
        for(xi=xb; xi<xe; xi+=nwf) { // wavefront loop (x direction)
            if(xe-xi <= nwf){
                nwf = xe-xi;
                end =1;
            }
            yb = yb_r;
            ye = ye_r;
//            MSG("xb=%d,xe=%d,yb=%d,ye=%d",xb,xe,yb,ye);
            kt = xi;
            kte=kt+nwf;

            float * __restrict v1_v;
            float * __restrict v2_v;
            float * __restrict v3_v;
            float * __restrict u1_v;
            float * __restrict u2_v;
            float * __restrict u3_v;
            const float * __restrict coef0_v;
            const float * __restrict inv_rho_v;

            //
            int t_real=0;
            int tb_real=(tb)/2+1;
            int t0_real=(t0)/2+1;
            for(int t=tb; t< te; t++){ // Diamond blocking in time
                t_real=(t)/2+1;
                hFloat* output_buffer = NULL;  //@KADIR
                int mod = (t)%2;
//                MSG("t=%d",t);
                if(mod){ // compute v from p
                    u1 = p11 ; //p
                    u2 = p12 ;
                    u3 = p13 ;
                    v1 = p21 ; //vx
                    v2 = p22 ; //vy
                    v3 = p23 ; //vz
//#pragma omp barrier
                    const Myfloat coef=stencil_ctx.dt;
                    for(int ix=kt; ix<kte; ix++){    // X
                        if( ((ix)/th_nwf)%th_x == tid_x ) {
                            for(int iy=yb; iy<ye; iy++) {
								v1_v = &(v1[ix*nnyz+iy*nnz]);
								v2_v = &(v2[ix*nnyz+iy*nnz]);
								v3_v = &(v3[ix*nnyz+iy*nnz]);
								u1_v = &(u1[ix*nnyz+iy*nnz]);
								u2_v = &(u2[ix*nnyz+iy*nnz]);
								u3_v = &(u3[ix*nnyz+iy*nnz]);
								inv_rho_v = &(inv_rho[(ix-NHALO)*nnyz_v+(iy-NHALO)*nnz_v]);
#pragma ivdep
								for(int iz=ib; iz<ie; iz++) {
									const Myfloat xum4 = u1_v[-3*nnyz + iz];
									const Myfloat xum3 = u1_v[-2*nnyz + iz];
									const Myfloat xum2 = u1_v[-1*nnyz + iz];
									const Myfloat xum1 = u1_v[ 0*nnyz + iz];
									const Myfloat xu0  = u1_v[ 1*nnyz + iz];
									const Myfloat xup1 = u1_v[ 2*nnyz + iz];
									const Myfloat xup2 = u1_v[ 3*nnyz + iz];
									const Myfloat xup3 = u1_v[ 4*nnyz + iz];

									Myfloat d_pr_x  = ( ( FDM_O1_8_2_A1 * (xu0  - xum1)
														  + FDM_O1_8_2_A2 * (xup1 - xum2)
														  + FDM_O1_8_2_A3 * (xup2 - xum3)
														  + FDM_O1_8_2_A4 * (xup3 - xum4)) ) ;

									v1_v[iz] += inv_rho_v[iz]*dt_inv_dx* d_pr_x;

									const Myfloat yum4 = u1_v[-3*nnz + iz];
									const Myfloat yum3 = u1_v[-2*nnz + iz];
									const Myfloat yum2 = u1_v[-1*nnz + iz];
									const Myfloat yum1 = u1_v[ 0*nnz + iz];
									const Myfloat yu0  = u1_v[ 1*nnz + iz];
									const Myfloat yup1 = u1_v[ 2*nnz + iz];
									const Myfloat yup2 = u1_v[ 3*nnz + iz];
									const Myfloat yup3 = u1_v[ 4*nnz + iz];

									Myfloat d_pr_y  = ( ( FDM_O1_8_2_A1 * (yu0  - yum1)
														  + FDM_O1_8_2_A2 * (yup1 - yum2)
														  + FDM_O1_8_2_A3 * (yup2 - yum3)
														  + FDM_O1_8_2_A4 * (yup3 - yum4)) ) ;

									v2_v[iz] += inv_rho_v[iz]*dt_inv_dy* d_pr_y;

									const Myfloat zum4 = u1_v[-3 + iz];
									const Myfloat zum3 = u1_v[-2 + iz];
									const Myfloat zum2 = u1_v[-1 + iz];
									const Myfloat zum1 = u1_v[ 0 + iz];
									const Myfloat zu0  = u1_v[ 1 + iz];
									const Myfloat zup1 = u1_v[ 2 + iz];
									const Myfloat zup2 = u1_v[ 3 + iz];
									const Myfloat zup3 = u1_v[ 4 + iz];

									Myfloat d_pr_z  = ( ( FDM_O1_8_2_A1 * (zu0  - zum1)
														  + FDM_O1_8_2_A2 * (zup1 - zum2)
														  + FDM_O1_8_2_A3 * (zup2 - zum3)
														  + FDM_O1_8_2_A4 * (zup3 - zum4)) ) ;

									v3_v[iz] += inv_rho_v[iz]*dt_inv_dz * d_pr_z;
								}
							}
                        }
                    }
                } else{// compute p from v
                    u1=	p21 ; //vx
                    u2=	p22 ; //vy
                    u3=	p23 ; //vz
                    v1=	p11 ; //p
                    v2=	p12 ;
                    v3=	p13 ;
                    for(int ix=kt; ix<kte; ix++){
                        if( ((ix)/th_nwf)%th_x == tid_x ) {
                            for(int iy=yb; iy<ye; iy++) {
//                                printf("iy=%d\n",iy);
                                v1_v = &(v1[ix*nnyz+iy*nnz]);
                                v2_v = &(v2[ix*nnyz+iy*nnz]);
                                v3_v = &(v3[ix*nnyz+iy*nnz]);
                                u1_v = &(u1[ix*nnyz+iy*nnz]);
                                u2_v = &(u2[ix*nnyz+iy*nnz]);
                                u3_v = &(u3[ix*nnyz+iy*nnz]);
//                                coef0_v = &(roc2[ix*nnyz+iy*nnz]); original
                                coef0_v = &(roc2[(ix-NHALO)*nnyz_v+(iy-NHALO)*nnz_v]);
#pragma ivdep
                                for(int iz=ib; iz<ie; iz++) {
                                    const Myfloat xum4 = u1_v[-4*nnyz + iz];
                                    const Myfloat xum3 = u1_v[-3*nnyz + iz];
                                    const Myfloat xum2 = u1_v[-2*nnyz + iz];
                                    const Myfloat xum1 = u1_v[-1*nnyz + iz];
                                    const Myfloat xu0  = u1_v[ 0*nnyz + iz];
                                    const Myfloat xup1 = u1_v[ 1*nnyz + iz];
                                    const Myfloat xup2 = u1_v[ 2*nnyz + iz];
                                    const Myfloat xup3 = u1_v[ 3*nnyz + iz];

                                    Myfloat d_vx_x  = ( ( FDM_O1_8_2_A1 * (xu0  - xum1)
                                                          + FDM_O1_8_2_A2 * (xup1 - xum2)
                                                          + FDM_O1_8_2_A3 * (xup2 - xum3)
                                                          + FDM_O1_8_2_A4 * (xup3 - xum4)) * inv_dx) ;

                                    const Myfloat yum4 = u2_v[-4*nnz + iz];
                                    const Myfloat yum3 = u2_v[-3*nnz + iz];
                                    const Myfloat yum2 = u2_v[-2*nnz + iz];
                                    const Myfloat yum1 = u2_v[-1*nnz + iz];
                                    const Myfloat yu0  = u2_v[ 0*nnz + iz];
                                    const Myfloat yup1 = u2_v[ 1*nnz + iz];
                                    const Myfloat yup2 = u2_v[ 2*nnz + iz];
                                    const Myfloat yup3 = u2_v[ 3*nnz + iz];

                                    Myfloat d_vy_y  = ( ( FDM_O1_8_2_A1 * (yu0  - yum1)
                                                          + FDM_O1_8_2_A2 * (yup1 - yum2)
                                                          + FDM_O1_8_2_A3 * (yup2 - yum3)
                                                          + FDM_O1_8_2_A4 * (yup3 - yum4)) * inv_dy) ;

                                    const Myfloat zum4 = u3_v[-4 + iz];
                                    const Myfloat zum3 = u3_v[-3 + iz];
                                    const Myfloat zum2 = u3_v[-2 + iz];
                                    const Myfloat zum1 = u3_v[-1 + iz];
                                    const Myfloat zu0  = u3_v[ 0 + iz];
                                    const Myfloat zup1 = u3_v[ 1 + iz];
                                    const Myfloat zup2 = u3_v[ 2 + iz];
                                    const Myfloat zup3 = u3_v[ 3 + iz];

                                    Myfloat d_vz_z  = ( ( FDM_O1_8_2_A1 * (zu0  - zum1)
                                                          + FDM_O1_8_2_A2 * (zup1 - zum2)
                                                          + FDM_O1_8_2_A3 * (zup2 - zum3)
                                                          + FDM_O1_8_2_A4 * (zup3 - zum4)) * inv_dz);
                                    v3_v[iz] += coef0_v[iz] * (d_vx_x + d_vy_y + d_vz_z);
                                    v3_v[iz]*=dampx[ix] * dampy[iy] * dampz[iz];
                                }
                                if (data->flag_bwd == 1) {
									///////  add sismos
									//////////////////////////////////////////
									double time_term = t0_real - (t_real - tb_real);
									int64_t term1 = (int64_t)data->rcv_len * (int64_t)time_term;
									int64_t term2 = (int64_t)(ix - 4) * (nny - 2 * NHALO);
									int64_t sismos_ind = term1 + term2 + (iy - 4);
									v3_v[iz_]+=data->sismos[sismos_ind];
								}
////                                MSG("data->flag_fwd=%d \n",data->flag_fwd);
                                if (data->flag_fwd == 1 && data->src_depth!=-1 && data->rec_sismos==1 ) {
//                                	MSG("recording sismos\n");
									///////  save sismos
									////////////////////////////////////////
//////////									data->sismos[data->rcv_len*(t0_real+(t_real-tb_real))+(ix-4)*(nny-2*NHALO)+(iy-4)]=(v3_v[iz_]);
									////////////////////////////////////////
//									MSG("ix=%d,iy=%d,iz=%d \n",ix,iy,iz_);
	//                                MSG("t0_real=%d,t_real=%d,tb_real=%d",t0_real,t_real,tb_real);
									double time_term = t0_real + (t_real - tb_real);
									int64_t term1 = (int64_t)data->rcv_len * (int64_t)time_term;
									int64_t term2 = (int64_t)(ix - 4) * (nny - 2 * NHALO);
									int64_t sismos_ind=term1 + term2 + (iy - 4);
	//								MSG("s_ind=%lld,term1=%lld,time_term=%f,v3_v[iz_]=%f \n",sismos_ind,term1,time_term,v3_v[iz_]);
									data->sismos[sismos_ind]=v3_v[iz_];
                                }
                            }
                            ///
                            if (data->flag_fwd == 1) {
								///////  add source
								if( (gp->source_point_enabled==1)
									&& (gp->lsource_pt[2] >= ib ) //@KADIR
									&& (gp->lsource_pt[2] <  ie ) //@KADIR
									&& (gp->lsource_pt[1] >= yb ) //@KADIR
									&& (gp->lsource_pt[1] <  ye ) //@KADIR
									&& (gp->lsource_pt[0] == ix ) )
								{
/////									ux[data->src_x] += data->source[t0+(t-tb)];// original
////									gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] = F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]) + gp->src_exc_coef[isrc_exc]);//@KADIR
//									gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+(gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] +=gp->src_exc_coef[isrc_exc];
//									if(0)  printf("DIA\tts:%d idzU:-- valU:%.4f src_exc_coef:%.4f coef:%g %g %g %g %g\ti(%d-%d) %d/%d\n", isrc_exc, H2F(gp->U1[((1ULL)*((gp->lsource_pt[2])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[0]))]),  gp->src_exc_coef[isrc_exc], coef[0], coef[1], coef[2], coef[3], coef[4],
//												  ib, ie, omp_get_thread_num(), omp_get_num_threads());

//									int src_index=(t0+(t-tb))/2;
//									MSG( "t0=%d,t=%d,tb=%d,src_index=%d,gp->src_exc_coef=%f",t0,t,tb,src_index,gp->src_exc_coef[ src_index ] );
//									MSG("isrc_exc=%d",isrc_exc);
									gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+(gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] +=gp->src_exc_coef[isrc_exc];
									isrc_exc++;
								}
                            }
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
        //////////////////////////////////////Forward//////////////////////////////////////
        if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // load fwd wavefield and compute IC
//        	MSG("fwd phase, ifwd=%d",ifwd);
			yb = yb_r;
			ye = ye_r;
			kt = xb;
//			kte=kt+nwf;
			kte=xe;
			hFloat *v3=p13;
			v3=p13;
//			MSG("fwd, ifwd=%d",ifwd);
			for(int t=tb; t< te; t++){ // Diamond blocking in time
				hFloat* output_buffer = NULL;  //@KADIR
				int mod = (t)%2;
//                MSG("t=%d",t);
				if(mod==0){// compute p from v
					u1=	p21 ;
////////////////////////
					for(int ix=kt; ix<kte; ix++){
						if( ((ix)/th_nwf)%th_x == tid_x ) {
							for(int iy=yb; iy<ye; iy++) {
								size_t index = 1ULL * ifwd * nnxyz + 1ULL * ix * nnyz + iy * nnz;
								if (index + ie >= stencil_ctx.fwd_size) {
									fprintf(stderr, "Thread %d: Out of bounds: index=%zu, fwd_size=%zu\n",
											omp_get_thread_num(),index,stencil_ctx.fwd_size);
									exit(1);
								}

								ux = &(v3[1ULL*ix*nnyz+iy*nnz]);
//								MSG("ifwd=%d",ifwd);
//								MSG("index=%d",1ULL*ifwd*nnxyz + 1ULL*ix*nnyz + iy*nnz);
								wx=&(data->fwd[1ULL*ifwd*nnxyz+1ULL*ix*nnyz+iy*nnz]);
#pragma ivdep
								for (int iz=ib; iz<ie; iz++) {
//									MSG("ix=%d,iy=%d,iz=%d,index=%d",ix,iy,iz,1ULL*ifwd*nnxyz + 1ULL*ix*nnyz + iy*nnz);
									wx[iz]=ux[iz];
								}
							}
							if( (gp->source_point_enabled==1)
								&& (gp->lsource_pt[2] >= ib ) //@KADIR
								&& (gp->lsource_pt[2] <  ie ) //@KADIR
								&& (gp->lsource_pt[1] >= yb ) //@KADIR
								&& (gp->lsource_pt[1] <  ye ) //@KADIR
								&& (gp->lsource_pt[0] == ix ) )	{
////									wx[data->src_x]-=data->source[t0+(t-tb)];	// delete source

//									gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] -= gp->src_exc_coef[isrc_exc2];
//									wx[  gp->lsource_pt[2] ] -= gp->src_exc_coef[isrc_exc2];
//								MSG("isrc_exc2=%d",isrc_exc2);
								data->fwd[1ULL*ifwd*nnxyz + 1ULL*ix*nnyz + gp->lsource_pt[1]*nnz+gp->lsource_pt[2]]-=gp->src_exc_coef[isrc_exc2];
								isrc_exc2++;
								//									gp->src_exc_coef[isrc_exc];
							}
						}
					}

////////////////////////
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
			} // diamond blocking in time (time loop)

        }
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
    }
}

void femwd_iso_ref_1st_grok(const int shape[3], const int zb, const int yb_r0, const int xb,
                       const int ze, const int ye_r0, const int xe, const real_t *coef,
                       hFloat *p11, hFloat *p12, hFloat *p13, hFloat *p21, hFloat *p22, hFloat *p23,
                       const hFloat *roc2, const hFloat *inv_rho, float *dampx, float *dampy, float *dampz,
                       int t_dim, int b_inc, int e_inc, int NHALO, int tb, int te, int t0, int ifwd,
                       stencil_ctx stencil_ctx, int mtid, tb_data_t *data)
{
    #pragma omp parallel num_threads(stencil_ctx.thread_group_size) \
    shared(shape, stencil_ctx, roc2, coef, mtid, tb, te, t_dim, NHALO, recv_rec, irecv_rec) \
    firstprivate(b_inc, e_inc)
    {
        const int nnx = shape[2], nny = shape[1], nnz = shape[0];
        const unsigned long nnyz = (unsigned long)nnz * nny;
        const int64_t nnxyz = (int64_t)nnx * nny * nnz;
        const int nnz_v = stencil_ctx.nz;
        const unsigned long nnyz_v = (unsigned long)stencil_ctx.nz * stencil_ctx.ny;

        int tgs = stencil_ctx.thread_group_size, nwf = stencil_ctx.num_wf;
        int tid = omp_get_thread_num(), gtid = tid + mtid * tgs;

        if (stencil_ctx.use_manual_cpu_bind) {
            if (sched_setaffinity(0, stencil_ctx.setsize, stencil_ctx.bind_masks[mtid * tgs + tid]) == -1) {
                printf("WARNING: Could not set CPU Affinity\n");
            }
        }

        hFloat *u1 = p11, *u2 = p12, *u3 = p13, *v1 = p21, *v2 = p22, *v3 = p23;
        int th_z = stencil_ctx.th_z, th_y = stencil_ctx.th_y, th_x = stencil_ctx.th_x;
        int tid_z = tid % th_z, tid_y = tid / th_z, tid_x = tid / (th_z * th_y);
        int yb_r = yb_r0, ye_r = ye_r0;

        if (th_y > 1) {
            if (b_inc && e_inc) {
                if (tid_y % 2 == 0) {
                    ye_r = (yb_r + ye_r) / 2;
                    e_inc = 0;
                } else {
                    yb_r = (yb_r + ye_r) / 2;
                    b_inc = 0;
                }
            } else {
                th_x *= th_y;
                tid_x = tid / th_z;
                if (nwf < th_x) nwf = th_x;
            }
        }

        int q = (ze - zb) / th_z, r = (ze - zb) % th_z;
        int ib = tid_z < r ? zb + tid_z * (q + 1) : zb + r * (q + 1) + (tid_z - r) * q;
        int ie = tid_z < r ? ib + (q + 1) : ib + q;
        int th_nwf = nwf / th_x, end = 0, iz_ = data->rcv_depth;

        const Myfloat dt_inv_dx = stencil_ctx.dt / stencil_ctx.dx;
        const Myfloat dt_inv_dy = stencil_ctx.dt / stencil_ctx.dy;
        const Myfloat dt_inv_dz = stencil_ctx.dt / stencil_ctx.dz;
        const Myfloat inv_dx = 1.0f / stencil_ctx.dx;
        const Myfloat inv_dy = 1.0f / stencil_ctx.dy;
        const Myfloat inv_dz = 1.0f / stencil_ctx.dz;

        // Backward phase
        if (data->flag_bwd && data->fwd && ifwd != -1) {
            int yb = yb_r, ye = ye_r, kt = xb, kte = xe;
            for (int t = tb; t < te; t++) {
                int mod = t % 2;
                if (mod == 0) {
                    #pragma omp simd
                    for (int ix = kt; ix < kte; ix++) {
                        if (((ix / th_nwf) % th_x) == tid_x) {
                            for (int iy = yb; iy < ye; iy++) {
                                float *restrict vx = &v3[ix * nnyz + iy * nnz];
                                float *restrict wx = &data->fwd[ifwd * nnxyz + ix * nnyz + iy * nnz];
                                float *restrict imgx = &data->img[(ix - NHALO) * nnyz_v + (iy - NHALO) * nnz_v];
                                float *restrict ilmx = &data->ilm[(ix - NHALO) * nnyz_v + (iy - NHALO) * nnz_v];
                                #pragma omp simd
                                for (int iz = ib; iz < ie; iz++) {
                                    imgx[iz] += vx[iz] * wx[iz];
                                    ilmx[iz] += wx[iz] * wx[iz];
                                }
                            }
                        }
                    }
                }
                yb += (t < t_dim) ? -b_inc : b_inc;
                ye += (t < t_dim) ? e_inc : -e_inc;
                kte = (end == 1) ? xe : max(kte - NHALO, xb);
                kt = max(kt - NHALO, xb);
            }
        }

        // Wavefield update
        for (int xi = xb; xi < xe; xi += nwf) {
            if (xe - xi <= nwf) {
                nwf = xe - xi;
                end = 1;
            }
            int yb = yb_r, ye = ye_r, kt = xi, kte = kt + nwf;
            int t_real = 0, tb_real = (tb / 2) + 1, t0_real = (t0 / 2) + 1;

            for (int t = tb; t < te; t++) {
                t_real = (t / 2) + 1;
                int mod = t % 2;
                double t_start = get_wall_time();

                if (mod) { // Velocity update
                    u1 = p11; u2 = p12; u3 = p13;
                    v1 = p21; v2 = p22; v3 = p23;
                    #pragma omp barrier
                    #pragma omp simd
                    for (int ix = kt; ix < kte; ix++) {
                        if (((ix / th_nwf) % th_x) == tid_x) {
                            for (int iy = yb; iy < ye; iy++) {
                                float *restrict v1_v = &v1[ix * nnyz + iy * nnz];
                                float *restrict v2_v = &v2[ix * nnyz + iy * nnz];
                                float *restrict v3_v = &v3[ix * nnyz + iy * nnz];
                                float *restrict u1_v = &u1[ix * nnyz + iy * nnz];
                                const float *restrict inv_rho_v = &inv_rho[(ix - NHALO) * nnyz_v + (iy - NHALO) * nnz_v];
                                #pragma omp simd
                                for (int iz = ib; iz < ie; iz++) {
                                    const Myfloat xum1 = u1_v[-nnyz + iz], xum2 = u1_v[-2 * nnyz + iz];
                                    const Myfloat xum3 = u1_v[-3 * nnyz + iz], xu0 = u1_v[iz];
                                    const Myfloat xup1 = u1_v[nnyz + iz], xup2 = u1_v[2 * nnyz + iz];
                                    const Myfloat xup3 = u1_v[3 * nnyz + iz];
                                    v1_v[iz] += inv_rho_v[iz] * dt_inv_dx * (
                                        FDM_O1_8_2_A1 * (xu0 - xum1) + FDM_O1_8_2_A2 * (xup1 - xum2) +
                                        FDM_O1_8_2_A3 * (xup2 - xum3) + FDM_O1_8_2_A4 * (u1_v[4 * nnyz + iz] - xum3));

                                    const Myfloat yum1 = u1_v[-nnz + iz], yum2 = u1_v[-2 * nnz + iz];
                                    const Myfloat yum3 = u1_v[-3 * nnz + iz], yu0 = u1_v[iz];
                                    const Myfloat yup1 = u1_v[nnz + iz], yup2 = u1_v[2 * nnz + iz];
                                    const Myfloat yup3 = u1_v[3 * nnz + iz];
                                    v2_v[iz] += inv_rho_v[iz] * dt_inv_dy * (
                                        FDM_O1_8_2_A1 * (yu0 - yum1) + FDM_O1_8_2_A2 * (yup1 - yum2) +
                                        FDM_O1_8_2_A3 * (yup2 - yum3) + FDM_O1_8_2_A4 * (u1_v[4 * nnz + iz] - yum3));

                                    const Myfloat zum1 = u1_v[-1 + iz], zum2 = u1_v[-2 + iz];
                                    const Myfloat zum3 = u1_v[-3 + iz], zu0 = u1_v[iz];
                                    const Myfloat zup1 = u1_v[1 + iz], zup2 = u1_v[2 + iz];
                                    const Myfloat zup3 = u1_v[3 + iz];
                                    v3_v[iz] += inv_rho_v[iz] * dt_inv_dz * (
                                        FDM_O1_8_2_A1 * (zu0 - zum1) + FDM_O1_8_2_A2 * (zup1 - zum2) +
                                        FDM_O1_8_2_A3 * (zup2 - zum3) + FDM_O1_8_2_A4 * (u1_v[4 + iz] - zum3));
                                }
                            }
                        }
                    }
                } else { // Pressure update
                    u1 = p21; u2 = p22; u3 = p23;
                    v1 = p11; v2 = p12; v3 = p13;
                    #pragma omp simd
                    for (int ix = kt; ix < kte; ix++) {
                        if (((ix / th_nwf) % th_x) == tid_x) {
                            for (int iy = yb; iy < ye; iy++) {
                                float *restrict v3_v = &v3[ix * nnyz + iy * nnz];
                                float *restrict u1_v = &u1[ix * nnyz + iy * nnz];
                                float *restrict u2_v = &u2[ix * nnyz + iy * nnz];
                                float *restrict u3_v = &u3[ix * nnyz + iy * nnz];
                                const float *restrict coef0_v = &roc2[(ix - NHALO) * nnyz_v + (iy - NHALO) * nnz_v];
                                #pragma omp simd
                                for (int iz = ib; iz < ie; iz++) {
                                    const Myfloat xum1 = u1_v[-nnyz + iz], xum2 = u1_v[-2 * nnyz + iz];
                                    const Myfloat xum3 = u1_v[-3 * nnyz + iz], xum4 = u1_v[-4 * nnyz + iz];
                                    const Myfloat xu0 = u1_v[iz], xup1 = u1_v[nnyz + iz];
                                    const Myfloat xup2 = u1_v[2 * nnyz + iz], xup3 = u1_v[3 * nnyz + iz];
                                    Myfloat d_vx_x = (FDM_O1_8_2_A1 * (xu0 - xum1) + FDM_O1_8_2_A2 * (xup1 - xum2) +
                                                     FDM_O1_8_2_A3 * (xup2 - xum3) + FDM_O1_8_2_A4 * (xup3 - xum4)) * inv_dx;

                                    const Myfloat yum1 = u2_v[-nnz + iz], yum2 = u2_v[-2 * nnz + iz];
                                    const Myfloat yum3 = u2_v[-3 * nnz + iz], yum4 = u2_v[-4 * nnz + iz];
                                    const Myfloat yu0 = u2_v[iz], yup1 = u2_v[nnz + iz];
                                    const Myfloat yup2 = u2_v[2 * nnz + iz], yup3 = u2_v[3 * nnz + iz];
                                    Myfloat d_vy_y = (FDM_O1_8_2_A1 * (yu0 - yum1) + FDM_O1_8_2_A2 * (yup1 - yum2) +
                                                     FDM_O1_8_2_A3 * (yup2 - yum3) + FDM_O1_8_2_A4 * (yup3 - yum4)) * inv_dy;

                                    const Myfloat zum1 = u3_v[-1 + iz], zum2 = u3_v[-2 + iz];
                                    const Myfloat zum3 = u3_v[-3 + iz], zum4 = u3_v[-4 + iz];
                                    const Myfloat zu0 = u3_v[iz], zup1 = u3_v[1 + iz];
                                    const Myfloat zup2 = u3_v[2 + iz], zup3 = u3_v[3 + iz];
                                    Myfloat d_vz_z = (FDM_O1_8_2_A1 * (zu0 - zum1) + FDM_O1_8_2_A2 * (zup1 - zum2) +
                                                     FDM_O1_8_2_A3 * (zup2 - zum3) + FDM_O1_8_2_A4 * (zup3 - zum4)) * inv_dz;

                                    v3_v[iz] += coef0_v[iz] * (d_vx_x + d_vy_y + d_vz_z) * dampx[ix] * dampy[iy] * dampz[iz];
                                }

                                if (data->flag_bwd) {
                                    double time_term = t0_real - (t_real - tb_real);
                                    int64_t sismos_ind = (int64_t)data->rcv_len * (int64_t)time_term +
                                                         (int64_t)(ix - 4) * (nny - 2 * NHALO) + (iy - 4);
                                    v3_v[iz_] += data->sismos[sismos_ind];
                                }
                                if (data->flag_fwd && data->src_depth != -1 && data->rec_sismos) {
                                    double time_term = t0_real + (t_real - tb_real);
                                    int64_t sismos_ind = (int64_t)data->rcv_len * (int64_t)time_term +
                                                         (int64_t)(ix - 4) * (nny - 2 * NHALO) + (iy - 4);
                                    data->sismos[sismos_ind] = v3_v[iz_];
                                }
                            }
                            if (data->flag_fwd && gp->source_point_enabled &&
                                gp->lsource_pt[2] >= ib && gp->lsource_pt[2] < ie &&
                                gp->lsource_pt[1] >= yb && gp->lsource_pt[1] < ye &&
                                gp->lsource_pt[0] == ix) {
                                gp->U1[((int64_t)gp->lsource_pt[0] * gp->ldomain_shape[1] +
                                        gp->lsource_pt[1]) * gp->ldomain_shape[0] + gp->lsource_pt[2]] += gp->src_exc_coef[isrc_exc++];
                            }
                        }
                    }
                }

                yb += (t < t_dim) ? -b_inc : b_inc;
                ye += (t < t_dim) ? e_inc : -e_inc;
                kte = (end == 1) ? xe : max(kte - NHALO, xb);
                kt = max(kt - NHALO, xb);
                stencil_ctx.t_wait[gtid] += get_wall_time() - t_start;
            }
        }

        // Forward phase
        if (data->flag_fwd && data->fwd && ifwd != -1) {
            int yb = yb_r, ye = ye_r, kt = xb, kte = xe;
            for (int t = tb; t < te; t++) {
                int mod = t % 2;
                if (mod == 0) {
                    #pragma omp simd
                    for (int ix = kt; ix < kte; ix++) {
                        if (((ix / th_nwf) % th_x) == tid_x) {
                            for (int iy = yb; iy < ye; iy++) {
                                size_t index = (size_t)ifwd * nnxyz + ix * nnyz + iy * nnz;
                                if (index + ie >= stencil_ctx.fwd_size) {
                                    fprintf(stderr, "Thread %d: Out of bounds: index=%zu, fwd_size=%zu\n",
                                            omp_get_thread_num(), index, stencil_ctx.fwd_size);
                                    exit(1);
                                }
                                float *restrict ux = &v3[ix * nnyz + iy * nnz];
                                float *restrict wx = &data->fwd[index];
                                #pragma omp simd
                                for (int iz = ib; iz < ie; iz++) {
                                    wx[iz] = ux[iz];
                                }
                            }
                            if (gp->source_point_enabled && gp->lsource_pt[2] >= ib &&
                                gp->lsource_pt[2] < ie && gp->lsource_pt[1] >= yb &&
                                gp->lsource_pt[1] < ye && gp->lsource_pt[0] == ix) {
                                data->fwd[(int64_t)ifwd * nnxyz + ix * nnyz + gp->lsource_pt[1] * nnz + gp->lsource_pt[2]] -=
                                    gp->src_exc_coef[isrc_exc2++];
                            }
                        }
                    }
                }
                yb += (t < t_dim) ? -b_inc : b_inc;
                ye += (t < t_dim) ? e_inc : -e_inc;
                kte = (end == 1) ? xe : max(kte - NHALO, xb);
                kt = max(kt - NHALO, xb);
            }
        }
    }
}

void intra_diamond_mwd_comp(Parameters *p, int yb_r, int ye_r, int b_inc,int e_inc, int tb, int te, int tid,int t0,int ifwd){
//	MSG("inside intra_diamond_mwd_comp_std t0=%d",t0);
    //@KADIR1 EXECUTED IN DIAMOND

    int t, x, xb[32], xe[32];
    int xb0,xe0;
    int yb, ye;
    int time_len = te-tb;
    double t1, t2, t3;

    // wavefront prologue
    t1 = get_wall_time();
    t2 = get_wall_time();
    // main wavefront loop
    yb = yb_r;
    ye = ye_r;
    //xb0 = (te-tb)*p->stencil.r;
    xb0 = p->stencil.r;
    xe0 = p->ldomain_shape[2]-p->stencil.r;

    int iifwd=ifwd;
    int fwd_steps=p->stencil_ctx.fwd_steps;

//    MSG("p->stencil_ctx.t_len=%d\n",p->stencil_ctx.t_len);
//    exit(1);

//    if (p->data->flag_fwd == 1) {
//		if (ifwd % fwd_steps == 0) iifwd = ifwd / fwd_steps;
//    }
//    if (p->data->flag_bwd == 1) {
//		if (ifwd % fwd_steps == 0) iifwd = ifwd / fwd_steps;
//    }
//    exit(1);

//	MSG("intra_diamond_mwd_comp, ifwd=%d\n",iifwd);
    p->stencil.mwd_func(p->ldomain_shape,p->stencil.r,yb,
                        xb0,p->lstencil_shape[0]+p->stencil.r, ye, xe0,
                        p->coef,p->U1,p->U1,p->U1,
                        p->U2, p->U3,p->U4,
						p->U5,p->U6,
                        p->dampx,p->dampy,p->dampz,
                        p->t_dim, b_inc, e_inc, p->stencil.r,
                        tb,te,t0,iifwd,p->stencil_ctx,tid,p->data);

    t3 = get_wall_time();
    p->stencil_ctx.t_wf_prologue[tid] += t2-t1;
    p->stencil_ctx.t_wf_main[tid]     += t3-t2;
    p->stencil_ctx.t_wf_epilogue[tid] += get_wall_time() - t3;
}

void intra_diamond_mwd_comp_std(Parameters *p, int yb_r, int ye_r, int b_inc,int e_inc, int tb, int te, int tid,int t_coord){
//	MSG("inside intra_diamond_mwd_comp_std t0=%d",t0);
    //@KADIR1 EXECUTED IN DIAMOND
    int t, x, xb[32], xe[32];
    int xb0,xe0;
    int yb, ye;
    int time_len = te-tb;
    double t1, t2, t3;

    // wavefront prologue
    t1 = get_wall_time();
    t2 = get_wall_time();
    // main wavefront loop
    yb = yb_r;
    ye = ye_r;
    //xb0 = (te-tb)*p->stencil.r;
    xb0 = p->stencil.r;
    xe0 = p->ldomain_shape[2]-p->stencil.r;

    int ifwd,iifwd,t0;
    iifwd = -1;
    int fwd_steps=p->stencil_ctx.fwd_steps;
//    int t_coord=st->t_pos[y_coord];
//    int t_coord=1;

//    MSG("p->stencil_ctx.t_len=%d\n",p->stencil_ctx.t_len);
//    exit(1);


    if (p->data->flag_fwd == 1) {
    	ifwd = t_coord;
		if (ifwd % fwd_steps == 0) iifwd = ifwd / fwd_steps;
		t0 = 1 + t_coord*(p->t_dim+1);
    }
    if (p->data->flag_bwd == 1) {
    	ifwd = p->t_len - 1 - t_coord;
		if (ifwd % fwd_steps == 0) iifwd = ifwd / fwd_steps;
//		t0 = ctx->time_steps - 2 - t_coord*(p->t_dim+1);original
		t0 = p->nt - 2 - t_coord*(p->t_dim+1);
    }
//    exit(1);
//    MSG("ifwd=%d,iifwd=%d\n",ifwd,iifwd);

    p->stencil.mwd_func(p->ldomain_shape,p->stencil.r,yb,
                        xb0,p->lstencil_shape[0]+p->stencil.r, ye, xe0,
                        p->coef,p->U1,p->U1,p->U1,
                        p->U2, p->U3,p->U4,
						p->U5,p->U6,
                        p->dampx,p->dampy,p->dampz,
                        p->t_dim, b_inc, e_inc, p->stencil.r,
                        tb,te,t0,iifwd,p->stencil_ctx,tid,p->data);

    t3 = get_wall_time();
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
    p->t_len=t_len;
    int num_thread_groups = get_ntg(*p);

    y_len_l = p->lstencil_shape[1] / (diam_width);
    y_len_r = y_len_l;
    if(p->is_last == 1) y_len_r++;

    int i, y, t;
    double t1,t2,t3,t4;
    int yb,ye;
    double db_t;
    int fwd_steps=p->stencil_ctx.fwd_steps;

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

    isrc_exc=0;
    isrc_exc2=0;
    MSG("isrc_exc=%d",isrc_exc);
    gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]=F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]) + gp->src_exc_coef[isrc_exc]);//@KADIR
    isrc_exc++;
    isrc_exc2++;

	MSG("p->data->flag_fwd=%d",p->data->flag_fwd);

	int ifwd,t0,t_coord;
    // Prologue
    MSG("ENTERING TB Prologue");
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
                        ifwd=-1;
                        if (p->data->flag_fwd==1){
                        	t0=1;
                        }
                        else{
                        	t0=p->nt-2;
                        }
                        if(p->stencil.type == REGULAR){
                            // we set t0=1, because it is the prologue. and we process diamonds with t0 equal to 1 there.
                        	intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, tb, te, tid,t0,ifwd);
                        }
                    }
                }
            }
        }
    }
    t2 = get_wall_time();

//	exit(1);

    // main loop
    //dynamic_intra_diamond_main_loop(p);
    MSG("ENTERING TB main loop");
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
                            t_coord = st.t_pos[y_coord];
//                            MSG("t_coord=%d\n",t_coord);
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
//                                        	MSG("t0=%d,t_coord=%d",t0,t_coord);
//                                        	printf("p->stencil.type= %s\n",p->stencil.type);
//                                        	MSG("Main loop. ifwd=%d",ifwd);
//                                        	t_coord=st->t_pos[y_coord];
//                                        	t_coord = get_t_coord(y_coord);
//                                        	MSG("t_coord=%d\n",t_coord);
//                                        	MSG("y_coord=%d,t_coord=%d\n",y_coord,t_coord);
                                            intra_diamond_mwd_comp_std(p,yb,ye,b_inc,e_inc,tb,te,tid,t_coord);
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
    MSG("ENTERING TB Epilogue");
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
//                    int t0 = t_len*(p->t_dim+1) + 1;

                    ifwd = -1;
                    if (p->data->flag_fwd == 1) {
                        t0 = p->t_len*(p->t_dim+1) + 1;}
                    if (p->data->flag_bwd == 1) {
						t0 = p->t_dim - 1;}
                    //@2
//					MSG("Epilogue. ifwd=%d",ifwd);
                    intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, tb, te, tid,t0,ifwd);
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

void dynamic_intra_diamond_ts_combined_backward(Parameters *p) {
    //@5
    MSG("inside dynamic_intra_diamond_ts_combined, backward phase.");
    int t_dim = p->t_dim;
    diam_width = (t_dim+1) * 2 *p->stencil.r;
    if(p->stencil.type == REGULAR){
        t_len = 2*( (p->nt-2)/((t_dim+1)*2) ) - 1;
    } else if(p->stencil.type == SOLAR){
        t_len = 2*( (p->nt)/((t_dim+1)*2) ) - 1;
    }
    p->t_len=t_len;
    int num_thread_groups = get_ntg(*p);

    y_len_l = p->lstencil_shape[1] / (diam_width);
    y_len_r = y_len_l;
    if(p->is_last == 1) y_len_r++;

    int i, y, t;
    double t1,t2,t3,t4;
    int yb,ye;
    double db_t;
    int fwd_steps=p->stencil_ctx.fwd_steps;

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

    isrc_exc=0;
    isrc_exc2=0;
    MSG("isrc_exc=%d",isrc_exc);
    gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]=F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]) + gp->src_exc_coef[isrc_exc]);//@KADIR
    isrc_exc++;
    isrc_exc2++;

	MSG("p->data->flag_fwd=%d",p->data->flag_fwd);

	int ifwd,t0,t_coord;
    // Prologue
    MSG("ENTERING TB Prologue, backward phase");
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
                        ifwd=-1;
                    	t0=1;
                        if(p->stencil.type == REGULAR){
                            // we set t0=1, because it is the prologue. and we process diamonds with t0 equal to 1 there.
                        	intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, tb, te, tid,t0,ifwd);
                        }
                    }
                }
            }
        }
    }
    t2 = get_wall_time();

//	exit(1);

    // main loop
    //dynamic_intra_diamond_main_loop(p);
    MSG("ENTERING TB main loop, backward phase");
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
                            t_coord = st.t_pos[y_coord];
//                            MSG("t_coord=%d\n",t_coord);
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
//                                        	MSG("t0=%d,t_coord=%d",t0,t_coord);
//                                        	printf("p->stencil.type= %s\n",p->stencil.type);
//                                        	MSG("Main loop. ifwd=%d",ifwd);
//                                        	t_coord=st->t_pos[y_coord];
//                                        	t_coord = get_t_coord(y_coord);
//                                        	MSG("t_coord=%d\n",t_coord);
//                                        	MSG("y_coord=%d,t_coord=%d\n",y_coord,t_coord);
                                            intra_diamond_mwd_comp_std(p,yb,ye,b_inc,e_inc,tb,te,tid,t_coord);
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
    MSG("ENTERING TB Epilogue, backward phase");
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
//                    int t0 = t_len*(p->t_dim+1) + 1;

                    ifwd = -1;
                    if (p->data->flag_fwd == 1) {
                        t0 = p->t_len*(p->t_dim+1) + 1;}
                    if (p->data->flag_bwd == 1) {
						t0 = p->t_dim - 1;}
                    //@2
//					MSG("Epilogue. ifwd=%d",ifwd);
                    intra_diamond_mwd_comp(p, yb, ye, b_inc, e_inc, tb, te, tid,t0,ifwd);
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


///////////////////////////////////////////////////////////////////////
void kernel_spatial_blocking_separate_mode_1st(const int nnx, const int nny, const int nnz,
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
                                               int b_inc,
                                               int e_inc,
                                               const int stencilr,
                                               const int tb,
                                               const int te,
                                               const int thread_group_size,
                                               const int groupid,
                                               const int setsize,
                                               cpu_set_t ** bind_masks,
                                               const tb_data_t* data,
                                               const int t0,
                                               const int ifwd,
                                               tb_timer_t* timer) {
    const int nnxy = nnx * nny;
    const int64_t nnxyz = 1ULL * nnx * nny * nnz;
    float *restrict u;
    float *restrict ux;

    float *restrict vx;
    float *restrict vx_x;

    float *restrict vy;
    float *restrict vy_x;

    float *restrict vz;
    float *restrict vz_x;

    const float *restrict rx;
    float *restrict wx;
    float *restrict imgx;
    float *restrict ilmx;
#pragma omp parallel default(shared) private(rx,vx_x,vy_x,vz_x,wx, imgx, ilmx) num_threads(thread_group_size)
    {
        int tid = 0;
        int gtid = 0;
#if defined(_OPENMP)
        tid = omp_get_thread_num();
        gtid = tid + groupid * thread_group_size;
#endif
        int err = sched_setaffinity(0, setsize, bind_masks[gtid]);
        if (err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n", groupid, tid, gtid);

        int yb, ye;
        // backward
        if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {
            yb = yb_r;
            ye = ye_r;
            for (int t = tb; t < te; t++) {
                u = u_r;
                if (t <= t_dim) {
#pragma omp for schedule(static)
                    for (int k = zb[t]; k < ze[t]; k++) {
                        for (int j = yb; j < ye; j++) {
                            // load fwd wavefield
                            // left most and right most
                            if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                                ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                                ux = &(u[1ULL * k * nnxy + j * nnx]);
                                wx = &(data->fwd[1ULL * ifwd * nnxyz + 1ULL * k * nnxy + j * nnx]);
                                imgx = &(data->img[1ULL * (k - 4) * (nnx - 8) * (nny - 8) + (j - 4) * (nnx - 8) - 4]);
                                ilmx = &(data->ilm[1ULL * (k - 4) * (nnx - 8) * (nny - 8) + (j - 4) * (nnx - 8) - 4]);
                                for (int i = xb; i < xe; i++) {
                                    imgx[i] += ux[i] * wx[i];
                                    ilmx[i] += wx[i] * wx[i];
                                }
                            }
                        }
                    }
#pragma omp barrier
                }
                if (t < t_dim) { // inverted trapezoid (or lower half of the diamond)
                    yb -= b_inc;
                    ye += e_inc;
                } else { // trapezoid  (or upper half of the diamond)
                    yb += b_inc;
                    ye -= e_inc;
                }
            }
        }

#pragma omp barrier
        yb = yb_r;
        ye = ye_r;
        for (int t = tb; t < te; t++) {
            u = u_r;
            vx = vx_r;
            vy = vy_r;
            vz = vz_r;
            int mod = (t) % 2;
#pragma omp for schedule(static)    // compute v from p
            for (int k = zb[t]; k < ze[t]; k++) {
                for (int j = yb; j < ye; j++) {
                    // compute
                    ux = &(u[1ULL * k * nnxy + j * nnx]);
                    vx_x = &(vx[1ULL * k * nnxy + j * nnx]);
                    vy_x = &(vy[1ULL * k * nnxy + j * nnx]);
                    vz_x = &(vz[1ULL * k * nnxy + j * nnx]);
//                    vx = &(   v[1ULL*k*nnxy + j*nnx]);
                    //          rx = &(roc2[1ULL*k*nnxy + j*nnx]);
                    rx = &(roc2[1ULL * (k - 4) * (nnx - 8) * (nny - 8) + (j - 4) * (nnx - 8) - 4]);
                    for (int i = xb; i < xe; i++) {
                        FUNC_BODY_1st_ord_Vsweep();     // compute v from p
                    }
                }
            }

#pragma omp for schedule(static)    // compute p from v
            for (int k = zb[t]; k < ze[t]; k++) {
                for (int j = yb; j < ye; j++) {
                    // compute
                    ux = &(u[1ULL * k * nnxy + j * nnx]);
                    vx_x = &(vx[1ULL * k * nnxy + j * nnx]);
                    vy_x = &(vy[1ULL * k * nnxy + j * nnx]);
                    vz_x = &(vz[1ULL * k * nnxy + j * nnx]);
//                        vx = &(   v[1ULL*k*nnxy + j*nnx]);
                    //          rx = &(roc2[1ULL*k*nnxy + j*nnx]);
                    rx = &(roc2[1ULL * (k - 4) * (nnx - 8) * (nny - 8) + (j - 4) * (nnx - 8) - 4]);
                    for (int i = xb; i < xe; i++) {
                        // compute p from v
                        FUNC_BODY_1st_ord_Psweep();
                    }
                    if (data->flag_bwd == 1) {
                        // add sismos
                        if (k == data->rcv_depth) {
                            for (unsigned int ir = 0; ir < data->rcv_len; ir++) {
                                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                                    ux[data->ix[ir]] += data->sismos[data->rcv_len * (t0 - (t - tb)) + ir];
                                }
                            }
                        }
                    }
                    if (data->flag_fwd == 1) {
                        // save sismos
                        if (k == data->rcv_depth) {
                            for (unsigned int ir = 0; ir < data->rcv_len; ir++) {
                                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                                    data->sismos[data->rcv_len * (t0 + (t - tb)) + ir] = ux[data->ix[ir]];
                                }
                            }
                        }
                        // add source
                        if (k == data->src_depth) {
                            if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) &&
                                (data->src_x < xe)) {
                                ux[data->src_x] += data->source[t0 + (t - tb)];
                            }
                        }
                    }
                }
            }

            if (t < t_dim) { // inverted trapezoid (or lower half of the diamond)
                yb -= b_inc;
                ye += e_inc;
            } else { // trapezoid  (or upper half of the diamond)
                yb += b_inc;
                ye -= e_inc;
            }
#pragma omp barrier
            /////
        }
        //////// forward
        if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {
            yb = yb_r;
            ye = ye_r;

            for (int t = tb; t < te; t++) {
                u = u_r;
                if (t >= t_dim) {
#pragma omp for schedule(static)
                    for (int k = zb[t]; k < ze[t]; k++) {
                        for (int j = yb; j < ye; j++) {
                            // save fwd
                            // left most and right most
                            if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                                ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                                ux = &(u[1ULL * k * nnxy + j * nnx]);
                                wx = &(data->fwd[1ULL * ifwd * nnxyz + 1ULL * k * nnxy + j * nnx]);
                                for (int i = xb; i < xe; i++) {
                                    wx[i]=ux[i];
                                }

                                if (k == data->src_depth) {
                                    if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) &&
                                        (data->src_x < xe)) {
                                        wx[data->src_x] -= data->source[t0 + (t - tb)];
                                    }
                                }
                            }
                        }
                    }
                }

                if (t < t_dim) { // inverted trapezoid (or lower half of the diamond)
                    yb -= b_inc;
                    ye += e_inc;
                } else { // trapezoid  (or upper half of the diamond)
                    yb += b_inc;
                    ye -= e_inc;
                }
#pragma omp barrier
            }
        }
    } // openmp
}

static void cpu_bind(tb_t *ctx) {
    int ncpus = get_nprocs();
    printf("ncpus : %d\n",ncpus);
    ctx->setsize = CPU_ALLOC_SIZE(ncpus);
    ctx->bind_masks = (cpu_set_t**) malloc(ctx->num_threads*sizeof(cpu_set_t*));

    for(int i = 0; i < ctx->num_threads; i++){
        ctx->bind_masks[i] = CPU_ALLOC( ncpus );
        if (ctx->bind_masks[i] == NULL) {
            ERR_MSG("bind mask error : %s",strerror(errno));
        }
        CPU_ZERO_S(ctx->setsize, ctx->bind_masks[i]);
//    CPU_SET_S(i, ctx->setsize, ctx->bind_masks[i]);
    }

    int *phys_cpu = (int*) calloc(ctx->num_threads,sizeof(int));

    printf("errno001 : %s\n",strerror(errno));
    //  omp_set_nested(1);
    int nb_level = 2;
    omp_set_max_active_levels (nb_level);
    printf("errno001 : %s\n",strerror(errno));

    printf("ntg %d, tgs %d \n",ctx->num_thread_groups,ctx->thread_group_size);
//  int* user_affinity = calloc(ctx->num_thread_groups*ctx->thread_group_size,sizeof(int));

    if (strcmp("NONE",ctx->affinity_file) == 0) {
        for(int i = 0; i < ctx->num_threads; i++){
            CPU_SET_S(i, ctx->setsize, ctx->bind_masks[i]);
        }
    } else {

        FILE * file = fopen(ctx->affinity_file,"r");
        CHK(file == NULL, "failed to open the affinity file");

        int value;

        for(int i = 0; i < ctx->num_threads; i++){
            fscanf(file,"%d",&value);
            CPU_SET_S(value, ctx->setsize, ctx->bind_masks[i]);
        }

        fclose(file);
    }

    // Set the affinity to reduce the cost of first run
#pragma omp parallel num_threads(ctx->num_thread_groups)
    {
        int mtid = omp_get_thread_num();
#pragma omp parallel shared(mtid) num_threads(ctx->thread_group_size)
        {
            int tid = omp_get_thread_num();
            int gtid = tid + mtid * ctx->thread_group_size;
            int err = sched_setaffinity(0, ctx->setsize, ctx->bind_masks[gtid]);
            if(err == -1) MSG("sched_setaffinity (th :%d) %s\n", gtid, strerror(errno));
            phys_cpu[gtid] = sched_getcpu();
        }
    }

    // checking
    printf("Threads binding (tid->OS tid):");
    for(int i = 0; i < ctx->num_threads; i++){
        printf(" %d->%d", i, phys_cpu[i]);
    }
    printf("\n");

    free(phys_cpu); phys_cpu = NULL;
}

void wave_tb_init(tb_t *ctx,
                  sismap_t *s,
                  parser *p) {
    // get thread info
    ctx->thread_group_size = parser_get_int(p, "tb_thread_group_size");
    ctx->num_thread_groups = parser_get_int(p, "tb_nb_thread_groups");
    ctx->num_threads = ctx->thread_group_size * ctx->num_thread_groups;

    ctx->th_x = parser_get_int(p, "tb_th_x");
    ctx->th_y = parser_get_int(p, "tb_th_y");
    ctx->th_z = parser_get_int(p, "tb_th_z");

    ctx->fwd_steps = parser_get_int(p,"fwd_steps");

    ctx->affinity_file = parser_get_string(p,"tb_affinity");

    ctx->mode = s->mode;
    switch (ctx->mode) {
        case 2:
            ctx->kernel_spatial_blocking_1st = kernel_spatial_blocking_separate_mode_1st;
//            ctx->kernel_tiling_blocking_1st = kernel_tiling_blocking_separate_mode_1st;
//            ctx->kernel_spatial_blocking = kernel_spatial_blocking_separate_mode;
//            ctx->kernel_tiling_blocking = kernel_tiling_blocking_separate_mode;
            break;
        default:
            ERR_MSG("tb_mode (%d) is unknown\n", ctx->mode);
    }


    // set thread affinity
    cpu_bind(ctx);

    // blocking setting
    ctx->t_dim = parser_get_int(p, "tb_t_dim");
    if (ctx->t_dim % 2 == 0) {
        ERR_MSG("t_dim (%d) should be odd\n", ctx->t_dim);
    }

    ctx->num_wf = parser_get_int(p, "tb_num_wf");


    // halo
    if ((s->sx != 4) || (s->sy != 4) || (s->sz != 4)) {
        ERR_MSG("sx (%d),sy (%d) and sz (%d) should be equal to 4\n", s->sx, s->sy, s->sz);
    }

    ctx->r = 4;

    // stride
    ctx->nnx = s->dimx + 2 * ctx->r;
    ctx->nny = s->dimy + 2 * ctx->r;
    ctx->nnz = s->dimz + 2 * ctx->r;

    // grid
    ctx->stencilx = s->dimx;
    ctx->stencily = s->dimy;
    ctx->stencilz = s->dimz;

    // coef
    ctx->coefx = s->coefx;
    ctx->coefy = s->coefy;
    ctx->coefz = s->coefz;

    ctx->diam_width = (ctx->t_dim + 1) * 2 * ctx->r;

    if (ctx->th_x * ctx->th_y * ctx->th_z != ctx->thread_group_size) {
        ERR_MSG("thread group (%d, %d, %d) is not match with its size (%d)\n",
                ctx->th_x, ctx->th_y, ctx->th_z, ctx->thread_group_size);
    }

    // check if thread group sizes are equal
    if (ctx->num_threads % ctx->thread_group_size != 0) {
        fprintf(stderr, "###ERROR: threads number must be multiples of thread group size\n");
        exit(1);
    }

    int diam_concurrency = ctx->stencily / ctx->diam_width; //t is mpi topology so its shape is 1x1x1
    if (ctx->num_thread_groups > diam_concurrency) {
        printf("###ERROR: the number of thread groups exceed the available concurrency. Consider using %d thread groups or less. \
        diam_concurrency:%d, num_thread_groups:%d, diam_width:%d, p->lstencil_shape[1]:%d, p->t.shape[1]:%d, p->lstencil_shape[1]/p->t.shape[1]:%d \
        \n", ((diam_concurrency > 1) ? diam_concurrency - 1 : 1), diam_concurrency, ctx->num_thread_groups,
               ctx->diam_width, ctx->stencily);
        exit(1);
    }

    if (ctx->stencily % ctx->diam_width != 0) {
        ERR_MSG("stencily (%d) should be multiple of diam_width (%d)\n", ctx->stencily, ctx->diam_width);
    }

    if( (ctx->num_wf%ctx->th_x != 0) && (ctx->thread_group_size != 1) ){
        fprintf(stderr,"ERROR: num_wavefronts must be multiples of thread groups size\n");
        exit(1);
    }
    if(ctx->t_dim < 1){
        fprintf(stderr,"ERROR: Diamond method does not support unrolling in time less than 1\n");
        exit(1);
    }

    if (ctx->stencily < ctx->r*(ctx->t_dim+1)*2){
        fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to fit at least one diamond: %d elements in Y [stencil_radius*2*(time_unrolls+1)]. Given %d elements\n"
                ,ctx->r*(ctx->t_dim+1)*2, ctx->stencily);
        exit(1);
    }
    if (floor(ctx->stencily/(ctx->r*(ctx->t_dim+1)*2.0)) != ctx->stencily / (ctx->r*(ctx->t_dim+1)*2.0)){
        fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to be multiples of the diamond width: %d elements [stencil_radius*2*(time_unrolls+1)]\n"
                ,ctx->r*(ctx->t_dim+1)*2);
        exit(1);
    }

    ////////// scheduling
    ////////// round the number of time steps to the nearest valid number
//    int nt_old=s->time_steps;
//    s->time_steps=s->time_steps*2;
//    int remain = (s->time_steps-2) % ((ctx->t_dim+1)*2);
//    if(remain != 0){
//        int nt2 = s->time_steps+(ctx->t_dim+1)*2 - remain;
//        if(nt2 != s->time_steps){
//            MSG("INFO: Modified nt from %03d to %03d for the intra-diamond method\n",nt_old,nt2/2);
//            s->time_steps=nt2;
//        }
//    }
    ////////////////////////////////////////
////  int remain=(s->time_steps) % ((ctx->t_dim+1)*2);
//    int remain = (s->time_steps-2) % ((ctx->t_dim+1)*2);  // pavel modification
//    if (remain != 0) {
//        int nt2 = s->time_steps+(ctx->t_dim+1)*2-remain;
//        if (nt2 != s->time_steps) {
//            MSG("INFO: Modified nt from %03d to %03d for the intra-diamond method\n",s->time_steps,nt2);
//            s->time_steps = nt2;
//        }
//    }
    ////////////////////////////////////////
    ctx->time_steps = s->time_steps;
    ctx->t_len = 2 * ((s->time_steps-2)/((ctx->t_dim + 1) * 2)) - 1;
    ctx->y_len_l = s->dimy/ctx->diam_width;
    ctx->y_len_r = ctx->y_len_l + 1;
    ////////////////////////////////////////

    ctx->t_pos = calloc(ctx->y_len_r,sizeof(int));
    ctx->avail_list = calloc(ctx->y_len_r,sizeof(int));

    // nb stencils main and total
    ctx->nb_stencils_main = ctx->t_len * (ctx->t_dim + 1)*ctx->stencilx * ctx->stencily * ctx->stencilz; //1LL *
//    ctx->nb_stencils_total_fwd =  ctx->time_steps * ctx->stencilx * ctx->stencily * ctx->stencilz; //1LL *
//    ctx->nb_stencils_total_bwd = (ctx->nb_stencils_total_fwd + ctx->nb_stencils_main) / 2;

//    MSG("ctx->nb_stencils_main=%d",ctx->nb_stencils_main);
//    MSG("ctx->nb_stencils_total_fwd=%d",ctx->nb_stencils_total_fwd);
//    MSG("ctx->nb_stencils_total_fwd=%" PRIu64 "\n", ctx.nb_stencils_total_fwd);
//    MSG("ctx->nb_stencils_total_bwd=%d",ctx->nb_stencils_total_bwd);

//    printf("ctx->nb_stencils_main=%d\n",ctx->nb_stencils_main);
//    printf("ctx->nb_stencils_total_fwd=%d\n",ctx->nb_stencils_total_fwd);
//    printf("ctx->nb_stencils_total_bwd=%d\n",ctx->nb_stencils_total_bwd);

//    MSG("ctx->nb_stencils_main=%llu\n", (unsigned long long)ctx->nb_stencils_main);
//    MSG("ctx->nb_stencils_total_fwd=%llu\n", (unsigned long long)ctx->nb_stencils_total_fwd);
//    MSG("ctx->nb_stencils_total_bwd=%llu\n", (unsigned long long)ctx->nb_stencils_total_bwd);

//    MSG("ctx->time_steps=%llu\n",ctx->time_steps);
//    MSG("ctx->nb_stencils_main=%llu\n", ctx->nb_stencils_main);
//    MSG("ctx->nb_stencils_total_fwd=%llu\n",ctx->nb_stencils_total_fwd);
//    MSG("ctx->nb_stencils_total_bwd=%llu\n",ctx->nb_stencils_total_bwd);

    // damping
    ctx->dampx = calloc(ctx->nnx, sizeof(float));
    ctx->dampy = calloc(ctx->nny, sizeof(float));
    ctx->dampz = calloc(ctx->nnz, sizeof(float));

    float alpha = 0.2;
    float tabdamp[NDAMP];


    for (int i = 1; i <= NDAMP; i++) {
        tabdamp[NDAMP - i] = exp(-alpha * (1.0 * i / NDAMP) * (1.0 * i / NDAMP));
    }

    for (int i = ctx->r; i < ctx->nnx - ctx->r; i++) {
        ctx->dampx[i] = 1.0;
    }
    for (int i = ctx->r; i < ctx->nny - ctx->r; i++) {
        ctx->dampy[i] = 1.0;
    }
    for (int i = ctx->r; i < ctx->nnz - ctx->r; i++) {
        ctx->dampz[i] = 1.0;
    }

    for (int i = 0; i < NDAMP; i++) {
        ctx->dampx[ctx->r + i] = tabdamp[i];
        ctx->dampy[ctx->r + i] = tabdamp[i];
        //ctx->dampz[ctx->r+i] = tabdamp[i];

        ctx->dampx[ctx->nnx - 1 - i - ctx->r] = tabdamp[i];
        ctx->dampy[ctx->nny - 1 - i - ctx->r] = tabdamp[i];
        ctx->dampz[ctx->nnz - 1 - i - ctx->r] = tabdamp[i];
    }
}

void wave_tb_init_p(tb_t *ctx,
                  sismap_t *s,
				  Parameters *p) {
	int enable_all_sizes = 0;
	if(enable_all_sizes == 0){
		p->nt=ctx->time_steps*2;
		// round the number of time steps to the nearest valid number
		int remain = (p->nt-2)%((p->t_dim+1)*2);
		if(remain != 0){
			int nt2 = p->nt + (p->t_dim+1)*2 - remain;
			if(nt2 != p->nt){
				if( (p->verbose ==1) ){
					printf("###INFO: Modified nt from %03d to %03d for the intra-diamond method to work properly\n",ctx->time_steps,nt2/2);
					fflush(stdout);
				}
				p->nt=nt2;
			}
		}
	}
}

void wave_tb_info(tb_t * ctx) {
    MSG(" ");
    MSG("-------------------------------------------");
    MSG("wave temporal blocking info");
    MSG("-------------------------------------------");
    MSG("thread group : (%d,%d,%d)",ctx->th_x,ctx->th_y,ctx->th_z);
    MSG("group size: %d num_group : %d",ctx->thread_group_size,ctx->num_thread_groups);
    MSG("-------------------------------------------");
    MSG("temporal blocking");
    MSG("t_dim : %d, num_wf : %d, diam_width : %d",ctx->t_dim, ctx->num_wf, ctx->diam_width);
    MSG("-------------------------------------------");
    MSG("scheduling info");
    MSG("t_len : %d, y_len_l : %d, y_lenr : %d",ctx->t_len, ctx->y_len_l, ctx->y_len_r);
    MSG("-------------------------------------------");
    MSG("grid info");
    MSG("stride       : (%d,%d,%d)",ctx->nnx, ctx->nny, ctx->nnz);
    MSG("stencil grid : (%d,%d,%d)",ctx->stencilx, ctx->stencily, ctx->stencilz);
    MSG("halo : %d",ctx->r);
    MSG("rtm info");
    MSG("fwd_steps : %d",ctx->fwd_steps);
    MSG("-------------------------------------------");
    MSG("mode info");
    switch (ctx->mode) {
        case 1:
            MSG("MODE (%d) : FUSED MODE",ctx->mode);
            break;

        case 2:
            MSG("MODE (%d) : SEPARATE MODE",ctx->mode);
            break;

        case 3:
            MSG("MODE (%d) : SEPARATE MODE with I/O",ctx->mode);
            break;

        default :
            MSG("MODE(%d) : UNKNOWN", ctx->mode);
    }
    MSG("-------------------------------------------");
}

void wave_tb_save_lastshot(sismap_t* s,
                           shot_t *shot,
                           float* u0,
                           float *u1) {
    MSG("inside wave_tb_save_lastshot");
    FILE * fd;

    char *snap_fd_name = (char*)malloc(20*sizeof(char));
    sprintf(snap_fd_name, "snapshot_u0_TB2nd_%u",s->time_steps);
    fd = fopen(snap_fd_name, "wb+");
    CHK(fd == NULL, "failed to open snapshot file");
    CHK(fwrite(u0, sizeof(float), s->size, fd) != s->size,
        "failed to write snapshot");
    CHK(fclose(fd)!=0,"failed to close snapshot file");
    if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);

    fd = 0;
    sprintf(snap_fd_name, "snapshot_TB2nd_%u",s->time_steps);
    fd = fopen(snap_fd_name, "wb+");
    CHK(fd == NULL, "failed to open snapshot file");
    CHK(fwrite(u1, sizeof(float), s->size, fd) != s->size,
        "failed to write snapshot");
    CHK(fclose(fd)!=0,"failed to close snapshot file");
    if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);
}

void wave_tb_save_lastshot_1st(sismap_t* s,
                               shot_t *shot,
                               float* u0) {
    MSG("inside wave_tb_save_lastshot_1st");
    FILE * fd;
    char *snap_fd_name = (char*)malloc(20*sizeof(char));
//    sprintf(snap_fd_name, "snapshot_TB1st_%u",s->time_steps/2);
    sprintf(snap_fd_name, "snapshot_TB1st_%u",s->time_steps);
    fd = fopen(snap_fd_name, "wb+");
    CHK(fd == NULL, "failed to open snapshot file");
    CHK(fwrite(u0, sizeof(float), s->size, fd) != s->size,"failed to write snapshot");
    CHK(fclose(fd)!=0,"failed to close snapshot file");
    if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);
    fd = 0;
}

void wave_sb_save_lastshot(sismap_t* s,
                           shot_t *shot,
                           float* u0,
                           float *u1,
                           unsigned int t) {
    FILE * fd;

    fd = fopen("SB_lastshot_u0", "wb+");
    CHK(fd == NULL, "failed to open snapshot file");
    CHK(fwrite(u0, sizeof(float), s->size, fd) != s->size,
        "failed to write snapshot");
    CHK(fclose(fd)!=0,"failed to close snapshot file");
    if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);

    fd = 0;

    fd = fopen("SB_lastshot_u1", "wb+");
    CHK(fd == NULL, "failed to open snapshot file");
    CHK(fwrite(u1, sizeof(float), s->size, fd) != s->size,
        "failed to write snapshot");
    CHK(fclose(fd)!=0,"failed to close snapshot file");
    if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);
    fd = 0;
}

void wave_tb_data_init(tb_data_t * data,
                       tb_t *tb,
                       sismap_t *s,
                       const int nb_thread_groups,
                       const int shotid,
                       size_t groupsize) {
    data->dx=s->dx;
    data->dy=s->dy;
    data->dz=s->dz;
    data->dt=s->dt;
    data->ix = calloc(s->rcv_len, sizeof(int));
    data->iy = calloc(s->rcv_len, sizeof(int));

    data->src_depth = s->src_depth;
    data->rcv_depth = s->rcv_depth;
    data->wave = NULL;

    data->gfwd_y0 = calloc(nb_thread_groups,sizeof(int));
    data->gfd     = calloc(nb_thread_groups,sizeof(FILE*));

    char tmp[512];
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    system(tmp);
    sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SNAP_BASE, shotid);


    for (int i = 0; i < nb_thread_groups; i++) {
        data->gfd[i] = fopen(tmp,"wb+");
        CHK(data->gfd[i] == NULL, "failed to open snapshot file");
    }

    strcpy(data->gfd_name,tmp);

    data->groupsize = groupsize;

//  mlbs_create(&(data->mlbs0));
//  mlbs_create(&(data->mlbs1));
//  mlbs_init(data->mlbs0,tb->nnx,tb->nny,tb->diam_width,tb->num_thread_groups,tmp, 63);
//  mlbs_init(data->mlbs1,tb->nnx,tb->nny,tb->diam_width,tb->num_thread_groups,tmp,127);

    bwriter_create(&(data->bwriter0));
    bwriter_create(&(data->bwriter1));
    bwriter_init(data->bwriter0,tb->nnx,tb->nny,tb->diam_width,tb->num_thread_groups,tmp);
    bwriter_init(data->bwriter1,tb->nnx,tb->nny,tb->diam_width,tb->num_thread_groups,tmp);
}

void wave_tb_data_set_src(tb_data_t * data,
                          sismap_t *s,
                          const unsigned int src_idx,
                          float * source) {
    data->source    = source;
    data->src_idx   = src_idx;
    data->src_depth = s->src_depth + 4;

//    data->src_x = data->src_idx % (s->dimx + 2* s->sx);
//    data->src_y = data->src_idx / (s->dimx + 2* s->sx);

    data->src_x = data->src_idx / (s->dimy + 2* s->sy);
    data->src_y = data->src_idx % (s->dimy + 2* s->sy);

    data->src_z = data->src_depth;
    MSG("SRC: ix=%d,iy=%d,iz=%d\n",data->src_x,data->src_y,data->src_z);
//    exit(0);
}

void wave_tb_data_unset_src(tb_data_t * data) {
    data->src_depth = -1;
}

void wave_tb_data_set_rcv(tb_data_t * data,
                          sismap_t *s,
                          float* sismos) {
    data->rcv       = s->rcv;
    data->rcv_len   = s->rcv_len;
    data->rcv_depth = s->rcv_depth + 4;
    MSG("data->rcv_depth=%d,s->rcv_depth=%d",data->rcv_depth,s->rcv_depth);
//    exit(1);
    data->sismos    = sismos;

//  printf("set rcv stride : %d\n",s->dimx + 2 * s->sx);

    for (int i = 0; i < s->rcv_len; i++) {
//    	data->ix[i] = s->rcv[i] % (s->dimx + 2 * s->sx);
//    	data->iy[i] = s->rcv[i] / (s->dimx + 2 * s->sx);

        data->ix[i] = s->rcv[i] / (s->dimy + 2 * s->sy);
        data->iy[i] = s->rcv[i] % (s->dimy + 2 * s->sy);
    }
}

void wave_tb_data_unset_rcv(tb_data_t * data) {
    data->rcv_depth = -1;
}

void wave_tb_data_set_wave(tb_data_t * data,
                           sismap_t *s) {
    data->wave = calloc(1ULL*(s->time_steps+1)
                        *(s->dimx + 2* s->sx)
                        *(s->dimy + 2* s->sy)
                        *(s->dimz + 2* s->sz),sizeof(float));
}

void wave_tb_data_unset_wave(tb_data_t* data) {
    free(data->wave);
    data->wave == NULL;
}

void wave_tb_data_dump_wave(tb_data_t *data,
                            sismap_t* s) {
    FILE * fd;
    size_t wavesize = 1ULL*(s->time_steps+1)
                      *(s->dimx + 2* s->sx)
                      *(s->dimy + 2* s->sy)
                      *(s->dimz + 2* s->sz);

    fd = fopen("forwardwave", "wb+");
    CHK(fd == NULL, "failed to open forwardwave file");
    CHK(fwrite(data->wave, sizeof(float), wavesize, fd) != wavesize,
        "failed to write forwardwave");
    CHK(fclose(fd)!=0,"failed to close forwardwave file");
    if (s->verbose) MSG("... saving last forwardwave (size %d)",wavesize);

    fd = 0;

}

void wave_tb_data_free(tb_data_t * data,
                       const int nb_thread_groups) {
//  mlbs_write(data->mlbs0,-1,NULL,0,0);
//  mlbs_write(data->mlbs1,-1,NULL,0,0);
//  mlbs_free(data->mlbs0);
//  mlbs_free(data->mlbs1);
//  mlbs_destroy(&data->mlbs0);
//  mlbs_destroy(&data->mlbs1);

    bwriter_free(data->bwriter0);
    bwriter_free(data->bwriter1);
    bwriter_destroy(&data->bwriter0);
    bwriter_destroy(&data->bwriter1);

    free(data->ix), data->ix = NULL;
    free(data->iy), data->iy = NULL;

    if (data->wave != NULL) {
        free(data->wave);
        data->wave = NULL;
    }

    for (int i = 0; i < nb_thread_groups; i ++) {
        fclose(data->gfd[i]); data->gfd[i] = NULL;
    }

    remove(data->gfd_name);
    memset(data->gfd_name,'\0',512);

    free(data->gfwd_y0); data->gfwd_y0 = NULL;
    free(data->gfd); data->gfd = NULL;

}

void wave_tb_data_info(tb_data_t* data) {
    MSG(" ");
    MSG("-------------------------------------------");
    MSG("wave data info");
    MSG("src_depth : %d",data->src_depth);
    if (data->src_depth != -1) {
        MSG("-------------------------------------------");
        MSG("src_depth : %d",data->src_depth);
        MSG("src_idx   : %d",data->src_idx);
        MSG("position  : (%d,%d,%d)",data->src_x,data->src_y,data->src_z);
    }

    if (data->rcv_depth != -1) {
        MSG("-------------------------------------------");
        MSG("rcv_depth : %d",data->rcv_depth);
        MSG("rcv_len   : %d",data->rcv_len);
/*
    for (int i = 0; i < data->rcv_len; i++) {
      MSG("%d : (%d,%d)",i,data->ix[i],data->iy[i]);
    }
*/
    }
    MSG("-------------------------------------------");
}

void wave_tb_free(tb_t* ctx) {
    int i;

    for (i = 0; i < ctx->num_threads; i++) {
        CPU_FREE(ctx->bind_masks[i]);
    }
    free(ctx->bind_masks); ctx->bind_masks = NULL;

    free((void*)(ctx->t_pos));      ctx->t_pos = NULL;
    free((void*)(ctx->avail_list)); ctx->avail_list = NULL;

    free(ctx->dampx); ctx->dampx = NULL;
    free(ctx->dampy); ctx->dampy = NULL;
    free(ctx->dampz); ctx->dampz = NULL;
}

static void intra_diamond_mwd_comp_1st(tb_t * ctx,
                                       tb_data_t *data,
                                       tb_timer_t* timer,
                                       float * restrict u0,
                                       float * restrict vx,
                                       float * restrict vy,
                                       float * restrict vz,
                                       const float * restrict roc2,
                                       int yb_r,
                                       int ye_r,
                                       int b_inc,
                                       int e_inc,
                                       int tb,
                                       int te,
                                       int t0,
                                       int ifwd,
                                       const int groupid) {
    int t, z, zb, ze;
    int yb, ye;
    int xb, xe;
    int zb_array[ctx->t_dim*2+1];
    int ze_array[ctx->t_dim*2+1];


    int time_len = te-tb;
    double t1,t2,t3,t4;

    t1 = wtime();

    // wavefront prologue
    yb = yb_r;
    ye = ye_r;
    xb = ctx->r;
    xe = ctx->stencilx + ctx->r;

    for(t = tb; t < te-1; t++) {
        zb_array[t] = ctx->r;
        ze_array[t] = ctx->r * (time_len - (t-tb));
    }

    ctx->kernel_spatial_blocking_1st(ctx->nnx, ctx->nny, ctx->nnz,
                                     ctx->r,                 yb, zb_array,
                                     ctx->r + ctx->stencilx, ye, ze_array,
                                     ctx->coefx, ctx->coefy, ctx->coefz,
                                     ctx->dampx, ctx->dampy, ctx->dampz,
                                     u0, vx,vy,vz, roc2,
                                     ctx->t_dim, b_inc, e_inc, ctx->r, tb, te - 1,
                                     ctx->thread_group_size,
                                     groupid,ctx->setsize,ctx->bind_masks,
                                     data,t0, ifwd, timer);
    t2 = wtime();

    // wavefront main loop
    yb = yb_r;
    ye = ye_r;
    zb = (te-tb) * ctx->r;
    ze = ctx->stencilz + ctx->r;

    ctx->kernel_tiling_blocking_1st(ctx->nnx, ctx->nny, ctx->nnz,
                                ctx->r,                 yb, zb,
                                ctx->r + ctx->stencilx, ye, ze,
                                ctx->coefx, ctx->coefy, ctx->coefz,
                                ctx->dampx, ctx->dampy, ctx->dampz,
                                u0, vx,vy,vz, roc2,
                                ctx->t_dim, b_inc, e_inc, ctx->r, tb, te,
                                ctx->num_wf, ctx->thread_group_size,
                                ctx->th_x,ctx->th_y,ctx->th_z,
                                groupid,ctx->setsize, ctx->bind_masks,
                                data, t0, ifwd, timer);

    t3 = wtime();

    // wavefront epilogue
    yb = yb_r;
    ye = ye_r;
    xb = ctx->r;
    xe = ctx->r + ctx->stencilx;

    if(tb < ctx->t_dim) { // lower half of the diamond
        yb -= b_inc;
        ye += e_inc;
    } else { // upper half of the diamond
        yb += b_inc;
        ye -= e_inc;
    }

    for(t = tb + 1; t < te; t++){
        ze_array[t] = ctx->stencilz + ctx->r;
        zb_array[t] = ctx->stencilz + ctx->r - (t-tb) * ctx->r;
    }

    ctx->kernel_spatial_blocking_1st(ctx->nnx, ctx->nny, ctx->nnz,
                                 ctx->r,                 yb, zb_array,
                                 ctx->r + ctx->stencilx, ye, ze_array,
                                 ctx->coefx, ctx->coefy, ctx->coefz,
                                 ctx->dampx, ctx->dampy, ctx->dampz,
                                 u0, vx,vy,vz, roc2,
                                 ctx->t_dim, b_inc, e_inc, ctx->r, tb + 1, te,
                                 ctx->thread_group_size,
                                 groupid,ctx->setsize,ctx->bind_masks,
                                 data,t0, ifwd, timer);

    t4 = wtime();

    timer->t_wf_prologue[groupid] += t2-t1;
    timer->t_wf_mainloop[groupid] += t3-t2;
    timer->t_wf_epilogue[groupid] += t4-t3;
}

static inline void intra_diamond_get_info(tb_t *ctx,
                                          tb_timer_t * timer,
                                          const int y_coord,
                                          const int t_coord,
                                          const int groupid,
                                          int *yb,
                                          int *ye,
                                          int *b_inc,
                                          int *e_inc){
    double diam_size;
    if( (y_coord == ctx->y_len_l - 1) && (t_coord % 2 == 0) ){ // right most process & left-half diamond
        // left half computations
        *yb = ctx->r + ctx->stencily - ctx->r;
        *ye = *yb + ctx->r;
        *b_inc = ctx->r;
        *e_inc = 0;
        diam_size = 0.5;
    }else if( (y_coord == ctx->y_len_r - 1) && (t_coord % 2 == 0) ){ // right most process & right-half diamond
        // right half computations
        *b_inc = 0;
        *e_inc = ctx->r;
        *yb = ctx->r;
        *ye = *yb + ctx->r;
        diam_size = 0.5;
    }else{ // full diamond computation
        if(t_coord % 2 == 0)// row shifted to the right
            *yb = ctx->r + ctx->diam_width   - ctx->r + y_coord * ctx->diam_width;
        else// row shifted to the left
            *yb = ctx->r + ctx->diam_width/2 - ctx->r + y_coord * ctx->diam_width;
        *ye = *yb + 2 * ctx->r;
        *b_inc = ctx->r;
        *e_inc = ctx->r;
        diam_size = 1;
    }

    timer->wf_num_resolved_diamonds[groupid] += diam_size;
}

static void fwd_setup(tb_t * ctx,
                      int * y0,
                      size_t * fwdsize,
                      size_t * offset,
                      const int y_coord,
                      const int t_coord,
                      const int iifwd) {

    if (t_coord % 2 == 0) {
        if (y_coord <= ctx->y_len_l-1) {
            *y0 = ctx->r + ctx->diam_width/2 + y_coord * ctx->diam_width;
        } else {
            *y0 = ctx->r;
        }
    } else {
        *y0 = ctx->r + ctx->diam_width * y_coord;
    }

    *fwdsize = 1LL* ctx->nnx * ctx->nny * ctx->diam_width;
    if ((t_coord % 2 == 0) && (y_coord >= ctx->y_len_l-1)) *fwdsize = (*fwdsize)/2;

    *offset = 1LL * (*y0) * ctx->nnx * ctx->nnz * sizeof(float) +
              1LL * iifwd   * ctx->nnx * ctx->nnz * ctx->nny * sizeof(float);

//  printf("%d %d (y_coord = %d, t_coord %d, y0 = %d,fwdsize = %ld)\n",ctx->y_len_r,ctx->y_len_l,y_coord,t_coord,*y0,*fwdsize);
}

static void fwd_load(tb_t * ctx,
                     FILE* fp,
                     float * fwd,
                     const int groupid,
                     const size_t offset,
                     const size_t fwdsize) {

    int rc;
    size_t size_done;

    rc = fseek(fp,offset, SEEK_SET);
    assert(rc == 0);

    size_done = fread(&(fwd[1LL* ctx->nnx * ctx->nnz * ctx->diam_width* groupid]),
                      sizeof(float),fwdsize, fp);
    assert(size_done != fwdsize);
}

static void fwd_save(tb_t * ctx,
                     FILE* fp,
                     float * fwd,
                     const int groupid,
                     const size_t offset,
                     const size_t fwdsize) {

    int rc;
    size_t size_done;

    rc = fseek(fp,offset, SEEK_SET);
    assert(rc == 0);

    size_done = fwrite(&(fwd[1LL* ctx->nnx * ctx->nnz * ctx->diam_width* groupid]),
                       sizeof(float),fwdsize, fp);
    assert(size_done != fwdsize);
}

void wave_tb_forward_1st(tb_t* ctx,
                         tb_data_t * data,
						 Parameters *p,
                         tb_timer_t* timer,
                         float * restrict u0,
                         float * restrict vx,
                         float * restrict vy,
                         float * restrict vz,
                         const float * restrict roc2,
						 const float *restrict inv_rho) {
    data->order=1;
    double t1,t2,t3,t4;
    if (data->src_depth !=-1) {
        u0[1ULL*data->src_depth* ctx->nnx * ctx->nny + data->src_idx] += data->source[0];
    }
    ////////////////////////////////////////////////////
//    Parameters *p = (Parameters*) calloc(1, sizeof(Parameters));
    ////////////////////////////////////////////////////
    if(p == NULL){
        printf("Error in allocating for Girih Parameters.\n");
        exit(0);
    }
    // girih print_param(p);
    p->lstencil_shape[0] = ctx->stencilz;
    p->lstencil_shape[1] = ctx->stencily;
    p->lstencil_shape[2] = ctx->stencilx;
    p->ldomain_shape[0] = ctx->nnz;
    p->ldomain_shape[1] = ctx->nny;
    p->ldomain_shape[2] = ctx->nnx;
    p->stencil_ctx.bs_y = ctx->stencily;
    p->stencil_ctx.dx = data->dx;
    p->stencil_ctx.dy = data->dy;
    p->stencil_ctx.dz = data->dz;
    p->stencil_ctx.nx = ctx->stencilx;
    p->stencil_ctx.ny = ctx->stencily;
    p->stencil_ctx.nz = ctx->stencilz;
    p->stencil_ctx.bs_y = ctx->stencily;

//    p->stencil_ctx.fwd_steps = ctx->fwd_steps;

    p->target_ts = 2;
    p->stencil.r = 4; // Stencil Kernel semi-bandwidth
    p->n_tests = 1;
    p->verbose = 1;
    p->t_dim=ctx->t_dim;
//    printf("data.dt:%f\n",data->dt);
//    exit(0);
    p->stencil_ctx.fwd_size=ctx->fwd_size;
    p->stencil_ctx.dt=data->dt;
    ///////////////////////////////////////////////////
    p->nt=ctx->time_steps*2;
    int enable_all_sizes = 0;
    if(enable_all_sizes == 0){
        // round the number of time steps to the nearest valid number
        int remain = (p->nt-2)%((p->t_dim+1)*2);
        if(remain != 0){
            int nt2 = p->nt + (p->t_dim+1)*2 - remain;
            if(nt2 != p->nt){
                if( (p->verbose ==1) ){
                    printf("###INFO: Modified nt from %03d to %03d for the intra-diamond method to work properly\n",ctx->time_steps,nt2/2);
                    fflush(stdout);
                }
                p->nt=nt2;
            }
        }
    }
    ///////////////////////////////////////////////////
    uint64_t nelm = (uint64_t) p->lstencil_shape[0] * p->lstencil_shape[1] * p->lstencil_shape[2];
    //    printf("nelm:%llu\n",nelm);
	uint64_t n_sample=0; n_sample=(uint64_t) p->nt * nelm;
//    printf("n_sample:%llu\n",n_sample);




	int t_len_tmp = 2*( (p->nt-2)/((ctx->t_dim+1)*2) ) - 1;

	ctx->nb_stencils_main =(uint64_t) t_len_tmp * (ctx->t_dim + 1)*ctx->stencilx * ctx->stencily * ctx->stencilz; //1LL *
//    ctx->nb_stencils_total_fwd =  p->nt * ctx->stencilx * ctx->stencily * ctx->stencilz; //1LL *
	ctx->nb_stencils_total_fwd =  n_sample; //1LL *

//	ctx->nb_stencils_total_bwd = (ctx->nb_stencils_total_fwd + ctx->nb_stencils_main) / 2;
	ctx->nb_stencils_total_bwd=n_sample;

	MSG("p->nt=%llu",p->nt);
	MSG("nelm=%llu", nelm);
//	MSG("stencil_nelm=%llu",ctx->stencilx * ctx->stencily * ctx->stencilz);
//	MSG("t_len_tmp=%llu", t_len_tmp);
//	MSG("(ctx->t_dim + 1)=%llu",  (ctx->t_dim + 1) );
	MSG("ctx->nb_stencils_main=%llu", ctx->nb_stencils_main);
	MSG("ctx->nb_stencils_total_fwd=%llu",ctx->nb_stencils_total_fwd);
	MSG("ctx->nb_stencils_total_bwd=%llu",ctx->nb_stencils_total_bwd);
//	exit(1);

    // same as SB
//  n_flop= p->nt * nelm * ((3 * 14) + 7+12) ;//pavel's proposed formula.
//	n_flop   = nt_corrected * nelm * ((6 * NB_OP_O2_8) + 10 + 4);  // TODO CHECK
    p->verify = 0;
    p->alignment = 8;
    p->source_point_enabled = 1;
    p->halo_concat = 0;
    int nthreads = 0;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    p->num_threads = nthreads;
    int tgs = ctx->th_x * ctx->th_y * ctx->th_z;
    if (nthreads % tgs){
        printf("nthreads %d cannot be dividable by thread group size %d. The remainder is %d\n. Exiting...\n", nthreads, tgs, nthreads%tgs);
        fflush(stdout);
        exit(0);
    }
    p->stencil_ctx.num_wf = ctx->num_wf;

    p->orig_thread_group_size = tgs;
    p->stencil_ctx.thread_group_size = tgs;

    // best on shaheen with tgs=4
    p->stencil_ctx.th_x = ctx->th_x;//4;//2;
    p->stencil_ctx.th_y = ctx->th_y;//2;
    p->stencil_ctx.th_z = ctx->th_z;//1;
    p->stencil_ctx.th_c = 1;
    p->stencil_ctx.fwd_steps = ctx->fwd_steps;
    p->t_len = ctx->t_len;

    p->th_block = 1;
    p->th_stride = 1;
    p->t.shape[0] = 1;
    p->t.shape[1] = 1;
    p->t.shape[2] = 1;

    // Kadir printed
    p->stencil.type = REGULAR;
    p->mwd_type = 1;
    p->wavefront = -1;
    p->is_last = 1;
    p->in_auto_tuning = 0;
    p->stencil_ctx.setsize = 8;
    p->stencil_ctx.use_manual_cpu_bind = 1;

    int diam_width = (p->t_dim+1) * 2 * p->stencil.r;
    int diam_concurrency = (p->lstencil_shape[1]/p->t.shape[1]) / diam_width; //t is mpi topology so its shape is 1x1x1
    int num_thread_groups = get_ntg(*p);

    if(num_thread_groups > diam_concurrency)
    {
        printf("###ERROR: the number of thread groups exceed the available concurrency. Consider using %d thread groups or less. \
                diam_concurrency:%d, num_thread_groups:%d, diam_width:%d, p->lstencil_shape[1]:%d, p->t.shape[1]:%d, p->lstencil_shape[1]/p->t.shape[1]:%d \
                \n", ((diam_concurrency>1)?diam_concurrency-1:1), diam_concurrency, num_thread_groups, diam_width, p->lstencil_shape[1], p->t.shape[1], (p->lstencil_shape[1]/p->t.shape[1]));
        exit(1);
    }
    // check for thread assignment validity
    if(p->stencil_ctx.thread_group_size > p->num_threads){
        printf("###WARNING: Requested thread group size is larger the total available threads \n");
    }

    // check thread group size validity
    if(p->stencil_ctx.thread_group_size != p->stencil_ctx.th_x * p->stencil_ctx.th_y*
                                           p->stencil_ctx.th_z * p->stencil_ctx.th_c){
        fprintf(stderr, "###ERROR: Thread group size must be consistent with parallelizm in all dimensions\n");
        exit(1);
    }

    // check number of threads along z-axis validity
    if(p->stencil_ctx.th_z > p->lstencil_shape[0]){
        fprintf(stderr, "###ERROR: no sufficient concurrency along the z-axis\n");
        exit(1);
    }
    // check if thread group sizes are equal
    if(p->num_threads%p->stencil_ctx.thread_group_size != 0){
        fprintf(stderr, "###ERROR: threads number must be multiples of thread group size\n");
        exit(1);
    }

    if( (p->stencil_ctx.num_wf%p->stencil_ctx.th_x != 0) && (p->stencil_ctx.thread_group_size != 1) ){
        fprintf(stderr,"ERROR: num_wavefronts must be multiples of thread groups size\n");
        exit(1);
    }
    if(p->t_dim < 1){
        fprintf(stderr,"ERROR: Diamond method does not support unrolling in time less than 1\n");
        exit(1);
    }
    if (p->lstencil_shape[1] < p->stencil.r*(p->t_dim+1)*2){
        fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to fit at least one diamond: %d elements in Y [stencil_radius*2*(time_unrolls+1)]. Given %d elements\n"
                ,p->stencil.r*(p->t_dim+1)*2, p->lstencil_shape[1]);
        exit(1);
    }
    if (floor(p->lstencil_shape[1] / (p->stencil.r*(p->t_dim+1)*2.0)) != p->lstencil_shape[1] / (p->stencil.r*(p->t_dim+1)*2.0)){
        fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to be multiples of the diamond width: %d elements [stencil_radius*2*(time_unrolls+1)]\n"
                ,p->stencil.r*(p->t_dim+1)*2);
        exit(1);
    }

    /// define source
    p->lsource_pt[0] = data->src_x;
    p->lsource_pt[1] = data->src_y;
    p->lsource_pt[2] = data->src_z;
//    MSG("!!!0=%d\n",data->src_x);
//    MSG("!!0=%d, 1=%d, 2=%d\n",p->lsource_pt[0],p->lsource_pt[1],p->lsource_pt[2]);
    //////////////////
    p->stencil_ctx.idz = ((real_t)1.)/((real_t)data->dz);
    p->stencil_ctx.idy = ((real_t)1.)/((real_t)data->dy);
    p->stencil_ctx.idx = ((real_t)1.)/((real_t)data->dx);
    p->stencil_ctx.idzyx_sum = p->stencil_ctx.idz + p->stencil_ctx.idy + p->stencil_ctx.idx;

//	p->stencil.stat_sched_func = stat_sched_iso_ref;
    p->stencil.mwd_func = femwd_iso_ref_1st;

    // Allocate the wavefront profiling timers
//	int num_thread_groups = get_ntg(*p);
    p->stencil_ctx.t_wait        = (double *) malloc(sizeof(double)*p->num_threads);
    p->stencil_ctx.t_wf_main     = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_comm = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_prologue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_epilogue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.wf_num_resolved_diamonds = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_group_wait = (double *) malloc(sizeof(double)*num_thread_groups);
    ////////////////////////////////////////////////////
    p->U1 = u0;
    p->U2 = vx;
    p->U3 = vy;
    p->U4 = vz;
    p->U5 = roc2;
    p->U6 = inv_rho;

    size_t size_src_exc_coef = p->nt * sizeof(real_t);
    p->src_exc_coef = (real_t*) malloc(size_src_exc_coef);
    int it=0;
    for (it=0;it<p->nt;it++)
    {
        p->src_exc_coef[it] = data->source[it];
    }
    p->dampx=ctx->dampx;
    p->dampy=ctx->dampy;
    p->dampz=ctx->dampz;
    gp = p;
    ////////////////////////////////////////////////////
    double elapse_time = 0.0;
    double *p_elapse_time;
    p_elapse_time = &elapse_time;

    //////////////////////////////////////////////////////////////////////
    p->data=data;
    //////////////////////

    reset_timers(&(p->prof));
    reset_wf_timers(p);
    p->stencil_ctx.use_manual_cpu_bind=1;
    cpu_bind_init(p);
    //@KADIR gcall
    double wall0 = get_wall_time();

    t1 = wtime();
    t2 = wtime();
    dynamic_intra_diamond_ts_combined(p);
    t3 = wtime();
    t4 = wtime();
    //////////////////////
    double wall1 = get_wall_time();
    *p_elapse_time = wall1 - wall0;
    //////////////////////
    timer->ts_main   += (t3-t2);
    timer->ts_others += (t2-t1) + (t4-t3);
    timer->total     += timer->ts_main + timer->ts_others;
    //////////////////////
//    free(p->coef);
//    free(p);
}

void wave_tb_backward_1st(tb_t* ctx,
                          tb_data_t * data,
						  Parameters *p,
                          tb_timer_t* timer,
                          float * restrict u0,
                          float * restrict vx,
                          float * restrict vy,
                          float * restrict vz,
                          const float * restrict roc2,
						  const float *restrict inv_rho) {
	MSG("inside wave_tb_backward_1st");
    data->order=1;
    double t1,t2,t3,t4;
    ////////////////////////////////////////////////////
    p->U1 = u0;
    p->U2 = vx;
    p->U3 = vy;
    p->U4 = vz;
    p->U5 = roc2;
    p->U6 = inv_rho;
    gp = p;
    ////////////////////////////////////////////////////
    double elapse_time = 0.0;
    double *p_elapse_time;
    p_elapse_time = &elapse_time;

    //////////////////////////////////////////////////////////////////////
//    p->data=data;
    //////////////////////

    reset_timers(&(p->prof));
    reset_wf_timers(p);
    p->stencil_ctx.use_manual_cpu_bind=1;
    cpu_bind_init(p);

    //@KADIR gcall
    double wall0 = get_wall_time();
    t1 = wtime();
    t2 = wtime();
    dynamic_intra_diamond_ts_combined(p);
//    dynamic_intra_diamond_ts_combined_backward(p);
    t3 = wtime();
    t4 = wtime();
    //////////////////////
    double wall1 = get_wall_time();
    *p_elapse_time = wall1 - wall0;
    //////////////////////
    timer->ts_main   += (t3-t2);
    timer->ts_others += (t2-t1) + (t4-t3);
    timer->total     += timer->ts_main + timer->ts_others;
    //////////////////////
}

/////////////////////////////////////////////////////////////
void wave_tb_timer_init(tb_timer_t * timer,
                        const int thread_group_size,
                        const int num_thread_groups) {
    timer->thread_group_size = thread_group_size;
    timer->num_thread_groups = num_thread_groups;
    timer->num_threads       = thread_group_size * num_thread_groups;

    timer->compute = 0.0;
    timer->total   = 0.0;
    timer->others  = 0.0;
    timer->ts_main = 0.0;
    timer->ts_others = 0.0;

    timer->t_wait = calloc(thread_group_size*num_thread_groups,sizeof(double));
    timer->t_wf_prologue = calloc(num_thread_groups,sizeof(double));
    timer->t_wf_mainloop = calloc(num_thread_groups,sizeof(double));
    timer->t_wf_epilogue = calloc(num_thread_groups,sizeof(double));
    timer->t_wf_snapio   = calloc(num_thread_groups,sizeof(double));
    timer->t_group_wait  = calloc(num_thread_groups,sizeof(double));
    timer->wf_num_resolved_diamonds  = calloc(num_thread_groups,sizeof(double));
}

void wave_tb_timer_clear(tb_timer_t * timer) {
    int i;

    timer->compute = 0.0;
    timer->total   = 0.0;
    timer->others  = 0.0;
    timer->ts_main = 0.0;
    timer->ts_others = 0.0;

    for (i = 0; i < timer->num_threads;       i++) {
        timer->t_wait[i]        = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->t_wf_prologue[i] = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->t_wf_mainloop[i] = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->t_wf_epilogue[i] = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->t_wf_snapio[i] = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->t_group_wait[i]  = 0.0;
    }
    for (i = 0; i < timer->num_thread_groups; i++) {
        timer->wf_num_resolved_diamonds[i]  = 0.0;
    }
}

void wave_tb_timer_free(tb_timer_t * timer) {
    free(timer->t_wait);        timer->t_wait = NULL;
    free(timer->t_wf_prologue); timer->t_wf_prologue = NULL;
    free(timer->t_wf_mainloop); timer->t_wf_mainloop = NULL;
    free(timer->t_wf_epilogue); timer->t_wf_epilogue = NULL;
    free(timer->t_wf_snapio);   timer->t_wf_snapio = NULL;
    free(timer->t_group_wait);  timer->t_group_wait = NULL;
    free(timer->wf_num_resolved_diamonds); timer->wf_num_resolved_diamonds = NULL;
}

void wave_tb_timer_info(tb_timer_t * timer,
                        const int64_t nb_stencils_total,
                        const int64_t nb_stencils_main) {
    int i;
    MSG(" ");
    MSG("-------------------------------------------");
    MSG("Global info:");
    MSG("Total:        %f (s) -%06.2f%%",timer->total,timer->total/timer->total*100.0);
    MSG("mainloop:     %f (s) -%06.2f%%",timer->ts_main,timer->ts_main/timer->total*100.0);
    MSG("pro/epilogue: %f (s) -%06.2f%%",timer->ts_others,timer->ts_others/timer->total*100.0);
    MSG("-------------------------------------------");

    MSG("Speed info:");
    MSG("Total: %f GStencils/s",nb_stencils_total/1e9/timer->total);
    MSG("Main:  %f GStencils/s",nb_stencils_main /1e9/timer->ts_main);
    MSG("nb_stencils_total: %llu",nb_stencils_total);
    MSG("nb_stencils_main:  %llu",nb_stencils_main);
    MSG("-------------------------------------------");

    MSG("Wavefront info:");
    printf("%-30s", "Metric \\ core:");
    for(i = 0; i < timer->num_threads; i++) {
        printf("  core %02d  ", i);
    }
    printf("\n");

    printf("%-27s", "Wavefront synchronization [s]:");
    for(i = 0; i < timer->num_threads; i++) {
        printf("  %4.3e", timer->t_wait[i]);
    }
    printf("\n");

    printf("%-27s", "Wavefront synchronization [%]:");
    for(i = 0; i < timer->num_threads; i++) {
        printf("  %05.2f    ", timer->t_wait[i]/timer->total*100.0);
    }
    printf("\n");

    printf("%-27s", "Metric \\ thread group:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  group %02d ", i);
    }
    printf("\n");

    printf("%-27s", "Wavefront steady state [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %4.3e", timer->t_wf_mainloop[i]);
    }
    printf("\n");

    printf("%-27s", "Wavefront steady state [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %05.2f    ", timer->t_wf_mainloop[i]/(timer->ts_main+timer->ts_others)*100);
    }
    printf("\n");

    printf("%-27s", "Wavefront startup/end [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %4.3e", timer->t_wf_prologue[i] + timer->t_wf_epilogue[i]);
    }
    printf("\n");

    printf("%-27s", "Wavefront startup/end [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %05.2f    ",(timer->t_wf_prologue[i] + timer->t_wf_epilogue[i])/(timer->ts_main+timer->ts_others)*100);
    }
    printf("\n");

    /*
    printf("%-27s", "Wavefront communication [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) printf("  %e", timer->t_wf_comm[i]);
    printf("\n");

    printf("%-27s", "Wavefront communication [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) printf("  %05.2f       ", timer->t_wf_comm[i]/(timer->ts_main+timer->ts_others)*100);
    printf("\n");
    */

    printf("%-27s", "Wavefront I/O [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %4.3e", timer->t_wf_snapio[i]);
    }
    printf("\n");

    printf("%-27s", "Wavefront I/O [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %05.2f    ", timer->t_wf_snapio[i]/(timer->ts_main+timer->ts_others)*100);
    }
    printf("\n");

    printf("%-27s", "Wavefront others [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %4.3e",timer->ts_main+timer->ts_others - (timer->t_wf_mainloop[i] + timer->t_wf_prologue[i] + timer->t_wf_epilogue[i] + timer->t_wf_snapio[i]));
    }
    printf("\n");

    printf("%-27s", "Wavefront others [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %05.2f    ",(timer->ts_main+ timer->ts_others - (timer->t_wf_mainloop[i] + timer->t_wf_prologue[i] + timer->t_wf_epilogue[i]))/(timer->ts_main+ timer->ts_others)*100);
    }
    printf("\n");

    printf("%-27s", "Group spin-wait [s]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %4.3e", timer->t_group_wait[i]);
    }
    printf("\n");

    printf("%-27s", "Group spin-wait [%]:");
    for(i = 0; i < timer->num_thread_groups; i++) {
        printf("  %05.2f    ", timer->t_group_wait[i]/(timer->total)*100);
    }
    printf("\n");

    printf("%-27s", "Resolved diamonds:");
    for(i = 0; i < timer->num_thread_groups; i++) printf("  %4.3e", timer->wf_num_resolved_diamonds[i]);
    printf("\n");

    MSG("-------------------------------------------");
}
