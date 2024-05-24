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
#include "simwave/macros.h"
#include "simwave/sismap.h"
#include "simwave/parser.h"
#include "simwave/wave_tb.h"
#include "simwave/wtime.h"

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
}
//ux[i] = dampx[i] * ux[i] + (1 - dampx[i]) * vx[i];                          \
//ux[i] = dampy[j] * ux[i] + (1 - dampy[j]) * vx[i];                          \
//ux[i] = dampz[k] * ux[i] + (1 - dampz[k]) * vx[i];                          \

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

void kernel_spatial_blocking_separate_mode(const int nnx, const int nny, const int nnz,
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
  const int nnxy  = nnx * nny;
  const int64_t nnxyz = 1ULL * nnx * nny * nnz;

  float *restrict u;
  float *restrict v;

  float *restrict ux;
  const float *restrict vx;
  const float *restrict rx;
  float *restrict wx;
  float *restrict imgx;
  float *restrict ilmx;

  float coef0 = coefx[0] + coefy[0] + coefz[0];

  #pragma omp parallel default(shared) private(rx,vx,ux,wx,imgx,ilmx) num_threads(thread_group_size)
  {
    int tid  = 0;
    int gtid = 0;
    #if defined(_OPENMP)
    tid = omp_get_thread_num();
    gtid = tid + groupid * thread_group_size;
    #endif
    int err = sched_setaffinity(0, setsize, bind_masks[gtid]);
    if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

    int yb,ye;

    if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {
      yb = yb_r;
      ye = ye_r;

      for(int t = tb; t < te; t++) {
        if(t % 2 == 1) {
          u = u_r; v = v_r;
        } else {
          u = v_r; v = u_r;
        }

        if (t <= t_dim) {
          #pragma omp for schedule(static)
          for(int k = zb[t]; k < ze[t]; k++) {
            for(int j = yb; j < ye; j++) {
              // load fwd wavefield
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                vx   = &(        v[1ULL*k*nnxy + j*nnx]);
                wx   = &(data->fwd[1ULL*ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                imgx = &(data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                ilmx = &(data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                for (int i = xb; i < xe; i++) {
                  imgx[i] += vx[i]*wx[i];
                  ilmx[i] += wx[i]*wx[i];
                }
              }
            }
          }
          #pragma omp barrier
        }

        if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
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

    for(int t = tb; t < te; t++) {

      if(t % 2 == 1) {
        u = u_r; v = v_r;
      } else {
        u = v_r; v = u_r;
      }

      #pragma omp for schedule(static)
      for(int k = zb[t]; k < ze[t]; k++) {
        for(int j = yb;j < ye; j++) {

          // compute
          ux = &(   u[1ULL*k*nnxy + j*nnx]);
          vx = &(   v[1ULL*k*nnxy + j*nnx]);
          //          rx = &(roc2[1ULL*k*nnxy + j*nnx]);
          rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

          for(int i = xb; i < xe; i++) {
            FUNC_BODY()
          }

          if (data->flag_bwd == 1) {
            // add sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                }
              }
            }
          }

          if (data->flag_fwd == 1) {
            // save sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                }
              }
            }
            // add source
            if (k == data->src_depth) {
              if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) && (data->src_x < xe)) {
                ux[data->src_x] += data->source[t0+(t-tb)];
              }
            }
          }
        }
      }

      if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
        yb -= b_inc;
        ye += e_inc;
      } else { // trapezoid  (or upper half of the diamond)
        yb += b_inc;
        ye -= e_inc;
      }

      #pragma omp barrier
    }

    if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {

      yb = yb_r;
      ye = ye_r;

      for(int t = tb; t < te; t++) {

        if(t % 2 == 1) {
          u = u_r; v = v_r;
        } else {
          u = v_r; v = u_r;
        }

        if (t >= t_dim) {
          #pragma omp for schedule(static)
          for(int k = zb[t]; k < ze[t]; k++) {
            for(int j = yb; j < ye; j++) {

              // save fwd
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                ux = &(   u[1ULL*k*nnxy + j*nnx]);
              wx = &(   data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
              for (int i = xb; i < xe; i++) {
                wx[i] = ux[i];
              }

              if (k == data->src_depth) {
                if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) && (data->src_x < xe)) {
                  wx[data->src_x] -= data->source[t0+(t-tb)];
                }
              }
            }
          }
        }
      }

        if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
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

void kernel_spatial_blocking_separate_mode_io(const int nnx, const int nny, const int nnz,
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
  const int nnxy  = nnx * nny;
  const int64_t nnxyz = 1ULL * nnx * nny * nnz;

  float *restrict u;
  float *restrict v;

  float *restrict ux;
  const float *restrict vx;
  const float *restrict rx;
  float *restrict wx;
  float *restrict imgx;
  float *restrict ilmx;

  float coef0 = coefx[0] + coefy[0] + coefz[0];

  #pragma omp parallel default(shared) private(rx,vx,ux,wx,imgx,ilmx) num_threads(thread_group_size)
  {
    int tid  = 0;
    int gtid = 0;
    #if defined(_OPENMP)
    tid = omp_get_thread_num();
    gtid = tid + groupid * thread_group_size;
    #endif
    int err = sched_setaffinity(0, setsize, bind_masks[gtid]);
    if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

    int yb,ye;

    if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {
      yb = yb_r;
      ye = ye_r;

      for(int t = tb; t < te; t++) {
        if(t % 2 == 1) {
          u = u_r; v = v_r;
        } else {
          u = v_r; v = u_r;
        }

        if (t <= t_dim) {
          #pragma omp for schedule(static)
          for(int k = zb[t]; k < ze[t]; k++) {
            for(int j = yb; j < ye; j++) {
              // load fwd wavefield
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                vx   = &(        v[1ULL*k*nnxy + j*nnx]);
                wx   = &(data->fwd[1ULL*data->groupsize*groupid + (j-data->gfwd_y0[groupid])*nnx*nnz+k*nnx]);
                imgx = &(data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                ilmx = &(data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                for (int i = xb; i < xe; i++) {
                  imgx[i] += vx[i]*wx[i];
                  ilmx[i] += wx[i]*wx[i];
                }
              }
            }
          }
          #pragma omp barrier
        }

        if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
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

    for(int t = tb; t < te; t++) {

      if(t % 2 == 1) {
        u = u_r; v = v_r;
      } else {
        u = v_r; v = u_r;
      }

      #pragma omp for schedule(static)
      for(int k = zb[t]; k < ze[t]; k++) {
        for(int j = yb;j < ye; j++) {

          // compute
          ux = &(   u[1ULL*k*nnxy + j*nnx]);
          vx = &(   v[1ULL*k*nnxy + j*nnx]);
          //          rx = &(roc2[1ULL*k*nnxy + j*nnx]);
          rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

          for(int i = xb; i < xe; i++) {
            FUNC_BODY()
          }

          if (data->flag_bwd == 1) {
            // add sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                }
              }
            }
          }

          if (data->flag_fwd == 1) {
            // save sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                }
              }
            }

            // add source
            if (k == data->src_depth) {
              if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) && (data->src_x < xe)) {
                ux[data->src_x] += data->source[t0+(t-tb)];
              }
            }
          }
        }
      }

      if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
        yb -= b_inc;
        ye += e_inc;
      } else { // trapezoid  (or upper half of the diamond)
        yb += b_inc;
        ye -= e_inc;
      }

      #pragma omp barrier
    }

    if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) {

      yb = yb_r;
      ye = ye_r;

      for(int t = tb; t < te; t++) {

        if(t % 2 == 1) {
          u = u_r; v = v_r;
        } else {
          u = v_r; v = u_r;
        }

        if (t >= t_dim) {
          #pragma omp for schedule(static)
          for(int k = zb[t]; k < ze[t]; k++) {
            for(int j = yb; j < ye; j++) {

              // save fwd
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                ux = &(   u[1ULL*k*nnxy + j*nnx]);
       //       printf("(j %d, y0 = %d)\n",j,data->gfwd_y0[groupid]);
                wx = &(data->fwd[1ULL*data->groupsize*groupid + (j-data->gfwd_y0[groupid])*nnx*nnz+k*nnx]);

                for (int i = xb; i < xe; i++) {
                  wx[i] = ux[i];
                }

                if (k == data->src_depth) {
                  if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) && (data->src_x < xe)) {
                    wx[data->src_x] -= data->source[t0+(t-tb)];
                  }
                }
              }
            }
          }
        }

        if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
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

void kernel_spatial_blocking_fuse_mode(const int nnx, const int nny, const int nnz,
                                       const int xb,  int yb, const int* zb,
                                       const int xe,  int ye, const int* ze,
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
  int i, j, k, jb, je, t;

  const int nnxy  = nnx * nny;
  const int64_t nnxyz = 1ULL * nnx * nny * nnz;

  float *restrict u;
  float *restrict v;

  float *restrict ux;
  const float *restrict vx;
  const float *restrict rx;
  float *restrict wx;
  float *restrict imgx;
  float *restrict ilmx;

  float coef0 = coefx[0] + coefy[0] + coefz[0];

  for(t = tb; t < te; t++) {

    if(t % 2 == 1) {
      u = u_r; v = v_r;
    } else {
      u = v_r; v = u_r;
    }

    #pragma omp parallel default(shared) num_threads(thread_group_size)
    {
      int tid  = 0;
      int gtid = 0;
      #if defined(_OPENMP)
      tid = omp_get_thread_num();
      gtid = tid + groupid * thread_group_size;
      #endif
      int err = sched_setaffinity(0, setsize, bind_masks[gtid]);
      if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

      #pragma omp for private(k,j,i,rx,vx,ux,wx,imgx,ilmx) schedule(static)
      for(k = zb[t]; k < ze[t]; k++) {
        for(j = yb;j < ye; j++) {
          if (data->flag_bwd == 1) {
            // load fwd wavefield
            if ((t <= t_dim) && (data->fwd != NULL) && (ifwd != -1)) {
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                vx = &(   v[1ULL*k*nnxy + j*nnx]);
                wx = &(   data->fwd[1ULL*ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                imgx = &(   data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                ilmx = &(   data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                for (i = xb; i < xe; i++) {
                  imgx[i] += vx[i]*wx[i];
                  ilmx[i] += wx[i]*wx[i];
                }
              }
            }
          }

          // compute
          ux = &(   u[1ULL*k*nnxy + j*nnx]);
          vx = &(   v[1ULL*k*nnxy + j*nnx]);
          //          rx = &(roc2[1ULL*k*nnxy + j*nnx]);
          rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

          for(i = xb; i < xe; i++) {
            FUNC_BODY()
          }

          if (data->flag_bwd == 1) {
            // add sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                }
              }
            }
          }

          if (data->flag_fwd == 1) {
            // save sismos
            if (k == data->rcv_depth) {
              for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                if ((data->iy[ir] == j) && (data->ix[ir] >= xb) && (data->ix[ir] < xe)) {
                  data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                }
              }
            }

            // save fwd
            if ((t >= t_dim) && (data->fwd != NULL) && (ifwd != -1)) {
              // left most and right most
              if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                  ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                ux = &(   u[1ULL*k*nnxy + j*nnx]);
                wx = &(   data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                for (i = xb; i < xe; i++) {
                  wx[i] = ux[i];
                }
              }
            }

            // add source
            if (k == data->src_depth) {
              if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= xb) && (data->src_x < xe)) {
                ux[data->src_x] += data->source[t0+(t-tb)];
              }
            }
          }
        }
      }
    }

    if(t < t_dim) { // inverted trapezoid (or lower half of the diamond)
      yb -= b_inc;
      ye += e_inc;
    } else { // trapezoid  (or upper half of the diamond)
      yb += b_inc;
      ye -= e_inc;
    }

  }
}

void kernel_tiling_blocking_fuse_mode(const int nnx, const int nny, const int nnz,
                                      const int xb, const int yb_r0, const int zb,
                                      const int xe, const int ye_r0, const int ze,
                                      const float * restrict coefx,
                                      const float * restrict coefy,
                                      const float * restrict coefz,
                                      const float * restrict dampx,
                                      const float * restrict dampy,
                                      const float * restrict dampz,
                                      float * restrict u,
                                      float * restrict v,
                                      const float * restrict roc2,
                                      const int t_dim,
                                      int b_inc,
                                      int e_inc,
                                      const int stencilr,
                                      const int tb,
                                      const int te,
                                      const int num_wf,
                                      const int thread_group_size,
                                      const int threadx,
                                      const int thready,
                                      const int threadz,
                                      const int groupid,
                                      const int setsize,
                                      cpu_set_t ** bind_masks,
                                      tb_data_t* data,
                                      const int t0,
                                      const int ifwd,
                                      tb_timer_t* timer) {

  #pragma omp parallel default(shared) firstprivate(u, v, b_inc, e_inc) \
  num_threads(thread_group_size)
  {
    int tgs, nwf, th_nwf, zi, yb, ye, ib, ie, kt, t, k, j, i, q, r, err;

    const int nnxy  = nnx * nny;
    const int64_t nnxyz = 1ULL * nnx * nny * nnz;

    int tid  = 0;
    int gtid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
    gtid = tid + groupid * thread_group_size;
#endif

    err = sched_setaffinity(0, setsize, bind_masks[gtid]);
    if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

    float * restrict u_r = u;
    float * restrict v_r = v;
    float * restrict ux;
    float * restrict vx;
    const float * restrict rx;
    float * restrict wx;
    float * restrict imgx;
    float * restrict ilmx;

    int th_x = threadx;
    int th_y = thready;
    int th_z = threadz;

    // tid = tid_z*(th_x*th_y) + tid_y*th_x + tid_x
    int tid_x = tid % th_x;
    int tid_y = tid / th_x;
    int tid_z = tid / (th_x*th_y);

    int yb_r = yb_r0;
    int ye_r = ye_r0;

    nwf = num_wf;

    if(thready > 1 ){
      if(b_inc != 0 && e_inc != 0){ // split only at full diamonds
        if (tid_y % 2 == 0){ // left thread
          ye_r = (yb_r + ye_r)/2;
          e_inc = 0;
        } else{
          yb_r = (yb_r + ye_r)/2;
          b_inc = 0;
        }
      }else{// use the y-threads along z-axis make sure to use sufficient number of frontlines
        th_z *= th_y;
        tid_z = tid/th_x;
        if (nwf < th_z) nwf = th_z;
      }
    }

    int nbx = (xe-xb)/th_x;
    q = (int)((xe-xb)/th_x);
    r = (xe-xb)%th_x;
    if(tid_x < r) {
      ib = xb + tid_x * (q+1);
      ie = ib + (q+1);
    } else {
      ib = xb + r * (q+1) + (tid_x - r) * q;
      ie = ib + q;
    }
    th_nwf = nwf/th_z;

    float coef0 = coefx[0] + coefy[0] + coefz[0];
    double t1;

    for(zi = zb; zi < ze; zi += nwf) { // wavefront loop (Z direction)

      if(ze-zi < nwf) {
        nwf = ze-zi;
      }

      yb = yb_r;
      ye = ye_r;

      kt = zi;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        for(j = yb; j < ye; j++) { // Y
          for(k = kt; k < kt + nwf; k++) { // Z
            if( ((k-stencilr)/th_nwf) % th_z == tid_z ) {

              if (data->flag_bwd == 1) { // load fwd wavefield and compute IC
                if ((t <= t_dim) && (data->fwd != NULL) && (ifwd != -1)) {
                  // left most and right most
                  if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                      ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                    vx   = &(        v[1ULL*k*nnxy + j*nnx]);
                    wx   = &(data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                    imgx = &(data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                    ilmx = &(data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                    for (i = ib; i < ie; i++) {
                      imgx[i] += vx[i]*wx[i];
                      ilmx[i] += wx[i]*wx[i];
                    }
                  }
                }
              }

              // compute propagation
              ux = &(   u[1ULL*k*nnxy + j*nnx]);
              vx = &(   v[1ULL*k*nnxy + j*nnx]);
              rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

              for (i = ib; i < ie; i++) {
                FUNC_BODY()
              }

              /*
              // save wavefield
              if (data->wave != NULL) {
                wx = &(data->wave[1ULL*(t0+t-tb)*nnz*nnxy + 1ULL*k*nnxy + j*nnx]);
                for (i = xb; i < xe; i++) {
                  wx[i] = ux[i];
                }
              }
              */

              if (data->flag_bwd == 1) {
                // add sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                    }
                  }
                }
              }

              if (data->flag_fwd == 1) {
                // save sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                    }
                  }
                }

                // save fwd
                if ((t >= t_dim) && (data->fwd != NULL) && (ifwd != -1)) {
                  // left most and right most
                  if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                      ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                    ux = &(   u[1ULL*k*nnxy + j*nnx]);
                    wx = &(   data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                    for (i = ib; i < ie; i++) {
                      wx[i] = ux[i];
                    }
                  }
                }

                // add source
                if (k == data->src_depth) {
                  if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= ib) && (data->src_x < ie)) {
                    ux[data->src_x] += data->source[t0+(t-tb)];
                  }
                }
              }

            }
          }
        }

        // Update block size in Y
        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;

        t1 = wtime();
        #pragma omp barrier
        timer->t_wait[gtid] += wtime() - t1;

      } // diamond blocking in time (time loop)
    } // wavefront loop

  } // parallel region
}

void kernel_tiling_blocking_separate_mode(const int nnx, const int nny, const int nnz,
                                          const int xb, const int yb_r0, const int zb,
                                          const int xe, const int ye_r0, const int ze,
                                          const float * restrict coefx,
                                          const float * restrict coefy,
                                          const float * restrict coefz,
                                          const float * restrict dampx,
                                          const float * restrict dampy,
                                          const float * restrict dampz,
                                          float * restrict u,
                                          float * restrict v,
                                          const float * restrict roc2,
                                          const int t_dim,
                                          int b_inc,
                                          int e_inc,
                                          const int stencilr,
                                          const int tb,
                                          const int te,
                                          const int num_wf,
                                          const int thread_group_size,
                                          const int threadx,
                                          const int thready,
                                          const int threadz,
                                          const int groupid,
                                          const int setsize,
                                          cpu_set_t ** bind_masks,
                                          tb_data_t* data,
                                          const int t0,
                                          const int ifwd,
                                          tb_timer_t* timer) {

  #pragma omp parallel default(shared) firstprivate(u, v, b_inc, e_inc) \
  num_threads(thread_group_size)
  {
    int tgs, nwf, th_nwf, zi, yb, ye, ib, ie, kt, t, k, j, i, q, r, err;

    const int nnxy  = nnx * nny;
    const int64_t nnxyz = 1ULL * nnx * nny * nnz;

    int tid  = 0;
    int gtid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
    gtid = tid + groupid * thread_group_size;
#endif

    err = sched_setaffinity(0, setsize, bind_masks[gtid]);
    if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

    float * restrict u_r = u;
    float * restrict v_r = v;
    float * restrict ux;
    float * restrict vx;
    const float * restrict rx;
    float * restrict wx;
    float * restrict imgx;
    float * restrict ilmx;

    int th_x = threadx;
    int th_y = thready;
    int th_z = threadz;

    // tid = tid_z*(th_x*th_y) + tid_y*th_x + tid_x
    int tid_x = tid % th_x;
    int tid_y = tid / th_x;
    int tid_z = tid / (th_x*th_y);

    int yb_r = yb_r0;
    int ye_r = ye_r0;

    nwf = num_wf;

    if(thready > 1){
      if(b_inc != 0 && e_inc != 0){ // split only at full diamonds
        if (tid_y % 2 == 0){ // left thread
          ye_r = (yb_r + ye_r)/2;
          e_inc = 0;
        } else{
          yb_r = (yb_r + ye_r)/2;
          b_inc = 0;
        }
      }else{// use the y-threads along z-axis make sure to use sufficient number of frontlines
        th_z *= th_y;
        tid_z = tid/th_x;
        if (nwf < th_z) nwf = th_z;
      }
    }

    int nbx = (xe-xb)/th_x;
    q = (int)((xe-xb)/th_x);
    r = (xe-xb)%th_x;
    if(tid_x < r) {
      ib = xb + tid_x * (q+1);
      ie = ib + (q+1);
    } else {
      ib = xb + r * (q+1) + (tid_x - r) * q;
      ie = ib + q;
    }
    th_nwf = nwf/th_z;

    float coef0 = coefx[0] + coefy[0] + coefz[0];
    double t1;

    int th_nz = (ze-zb)/th_z;

    // Load wavefield and imaging condition
    if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // load fwd wavefield and compute IC
      yb = yb_r;
      ye = ye_r;

      kt = zb;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        if (t <= t_dim) {
          for(k = kt; k < kt+ze-zb; k ++) { // wavefront loop (Z direction)
            if( ((k-stencilr)/th_nz) % th_z == tid_z ) {
              for(j = yb; j < ye; j++) { // Y
                // left most and right most
                if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                    ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                  vx   = &(        v[1ULL*k*nnxy + j*nnx]);
                wx   = &(data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                imgx = &(data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                ilmx = &(data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                for (i = ib; i < ie; i++) {
                  imgx[i] += vx[i]*wx[i];
                  ilmx[i] += wx[i]*wx[i];
                }
              }
            }
          }
        }
      }

        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;
      }
    }

    #pragma omp barrier

    for(zi = zb; zi < ze; zi += nwf) { // wavefront loop (Z direction)

      if(ze-zi < nwf) {
        nwf = ze-zi;
      }

      yb = yb_r;
      ye = ye_r;

      kt = zi;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        for(j = yb; j < ye; j++) { // Y
          for(k = kt; k < kt + nwf; k++) { // Z
            if( ((k-stencilr)/th_nwf) % th_z == tid_z ) {

              // compute propagation
              ux = &(   u[1ULL*k*nnxy + j*nnx]);
              vx = &(   v[1ULL*k*nnxy + j*nnx]);
              rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

              for (i = ib; i < ie; i++) {
                FUNC_BODY()
              }

              if (data->flag_bwd == 1) {
                // add sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                    }
                  }
                }
              }

              if (data->flag_fwd == 1) {
                // save sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                    }
                  }
                }

                // add source
                if (k == data->src_depth) {
                  if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= ib) && (data->src_x < ie)) {
                    ux[data->src_x] += data->source[t0+(t-tb)];
                  }
                }
              }

            }
          }
        }

        // Update block size in Y
        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;

        t1 = wtime();
        #pragma omp barrier
        timer->t_wait[gtid] += wtime() - t1;

      } // diamond blocking in time (time loop)
    } // wavefront loop

    #pragma omp barrier

    if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // save fwd wavefield
      yb = yb_r;
      ye = ye_r;

      kt = zb;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        if (t >= t_dim) {
          for(k = kt; k < kt+ze-zb; k ++) { // wavefront loop (Z direction)
            if( ((k-stencilr)/th_nz) % th_z == tid_z ) {
              for(j = yb; j < ye; j++) { // Y
                // left most and right most
                if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                    ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                  ux = &(   u[1ULL*k*nnxy + j*nnx]);
                wx = &(   data->fwd[1ULL * ifwd*nnxyz + 1ULL*k*nnxy + j*nnx]);
                for (i = ib; i < ie; i++) {
                  wx[i] = ux[i];
                }

                if (k == data->src_depth) {
                  if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= ib) && (data->src_x < ie)) {
                    wx[data->src_x] -= data->source[t0+(t-tb)];
                  }
                }
              }
            }
          }
        }
      }

        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;
      }
    }
  } // parallel region
}

void kernel_tiling_blocking_separate_mode_io(const int nnx, const int nny, const int nnz,
                                             const int xb, const int yb_r0, const int zb,
                                             const int xe, const int ye_r0, const int ze,
                                             const float * restrict coefx,
                                             const float * restrict coefy,
                                             const float * restrict coefz,
                                             const float * restrict dampx,
                                             const float * restrict dampy,
                                             const float * restrict dampz,
                                             float * restrict u,
                                             float * restrict v,
                                             const float * restrict roc2,
                                             const int t_dim,
                                             int b_inc,
                                             int e_inc,
                                             const int stencilr,
                                             const int tb,
                                             const int te,
                                             const int num_wf,
                                             const int thread_group_size,
                                             const int threadx,
                                             const int thready,
                                             const int threadz,
                                             const int groupid,
                                             const int setsize,
                                             cpu_set_t ** bind_masks,
                                             tb_data_t* data,
                                             const int t0,
                                             const int ifwd,
                                             tb_timer_t* timer) {

  #pragma omp parallel default(shared) firstprivate(u, v, b_inc, e_inc) \
  num_threads(thread_group_size)
  {
    int tgs, nwf, th_nwf, zi, yb, ye, ib, ie, kt, t, k, j, i, q, r, err;

    const int nnxy  = nnx * nny;
    const int64_t nnxyz = 1ULL * nnx * nny * nnz;

    int tid  = 0;
    int gtid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
    gtid = tid + groupid * thread_group_size;
#endif

    err = sched_setaffinity(0, setsize, bind_masks[gtid]);
    if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d %d)\n",groupid,tid,gtid);

    float * restrict u_r = u;
    float * restrict v_r = v;
    float * restrict ux;
    float * restrict vx;
    const float * restrict rx;
    float * restrict wx;
    float * restrict imgx;
    float * restrict ilmx;

    int th_x = threadx;
    int th_y = thready;
    int th_z = threadz;

    // tid = tid_z*(th_x*th_y) + tid_y*th_x + tid_x
    int tid_x = tid % th_x;
    int tid_y = tid / th_x;
    int tid_z = tid / (th_x*th_y);

    int yb_r = yb_r0;
    int ye_r = ye_r0;

    nwf = num_wf;

    if(thready > 1){
      if(b_inc != 0 && e_inc != 0){ // split only at full diamonds
        if (tid_y % 2 == 0){ // left thread
          ye_r = (yb_r + ye_r)/2;
          e_inc = 0;
        } else{
          yb_r = (yb_r + ye_r)/2;
          b_inc = 0;
        }
      }else{// use the y-threads along z-axis make sure to use sufficient number of frontlines
        th_z *= th_y;
        tid_z = tid/th_x;
        if (nwf < th_z) nwf = th_z;
      }
    }

    int nbx = (xe-xb)/th_x;
    q = (int)((xe-xb)/th_x);
    r = (xe-xb)%th_x;
    if(tid_x < r) {
      ib = xb + tid_x * (q+1);
      ie = ib + (q+1);
    } else {
      ib = xb + r * (q+1) + (tid_x - r) * q;
      ie = ib + q;
    }
    th_nwf = nwf/th_z;

    float coef0 = coefx[0] + coefy[0] + coefz[0];
    double t1;

    int th_nz = (ze-zb)/th_z;

    // Load wavefield and imaging condition
    if ((data->flag_bwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // load fwd wavefield and compute IC
      yb = yb_r;
      ye = ye_r;

      kt = zb;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        if (t <= t_dim) {
          for(k = kt; k < kt+ze-zb; k ++) { // wavefront loop (Z direction)
            if( ((k-stencilr)/th_nz) % th_z == tid_z ) {
              for(j = yb; j < ye; j++) { // Y
                // left most and right most
                if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                    ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                  vx   = &(        v[1ULL*k*nnxy + j*nnx]);
                  wx   = &(data->fwd[1ULL*data->groupsize*groupid + (j-data->gfwd_y0[groupid])*nnx*nnz+k*nnx]);
                  imgx = &(data->img[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                  ilmx = &(data->ilm[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);
                  for (i = ib; i < ie; i++) {
                    imgx[i] += vx[i]*wx[i];
                    ilmx[i] += wx[i]*wx[i];
                  }
                }
              }
            }
          }
        }

        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;
      }
    }

    #pragma omp barrier

    for(zi = zb; zi < ze; zi += nwf) { // wavefront loop (Z direction)

      if(ze-zi < nwf) {
        nwf = ze-zi;
      }

      yb = yb_r;
      ye = ye_r;

      kt = zi;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        for(j = yb; j < ye; j++) { // Y
          for(k = kt; k < kt + nwf; k++) { // Z
            if( ((k-stencilr)/th_nwf) % th_z == tid_z ) {

              // compute propagation
              ux = &(   u[1ULL*k*nnxy + j*nnx]);
              vx = &(   v[1ULL*k*nnxy + j*nnx]);
              rx = &(roc2[1ULL*(k-4)*(nnx-8)*(nny-8) + (j-4)*(nnx-8) - 4]);

              for (i = ib; i < ie; i++) {
                FUNC_BODY()
              }

              if (data->flag_bwd == 1) {
                // add sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      ux[data->ix[ir]] += data->sismos[data->rcv_len*(t0-(t-tb))+ir];
                    }
                  }
                }
              }

              if (data->flag_fwd == 1) {
                // save sismos
                if (k == data->rcv_depth) {
                  for (unsigned int ir = 0; ir < data->rcv_len; ir ++) {
                    if ((data->iy[ir] == j) && (data->ix[ir] >= ib) && (data->ix[ir] < ie)) {
                      data->sismos[data->rcv_len*(t0+(t-tb))+ir] = ux[data->ix[ir]];
                    }
                  }
                }

                // add source
                if (k == data->src_depth) {
                  if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= ib) && (data->src_x < ie)) {
                    ux[data->src_x] += data->source[t0+(t-tb)];
                  }
                }
              }

            }
          }
        }

        // Update block size in Y
        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;

        t1 = wtime();
        #pragma omp barrier
        timer->t_wait[gtid] += wtime() - t1;

      } // diamond blocking in time (time loop)
    } // wavefront loop

    #pragma omp barrier

//    MSG("... fwd=: %d\n", data->fwd );
//    MSG("... flag_fwd=: %d\n", data->flag_fwd );
//    MSG("... ifwd=: %d\n", ifwd );
    if ((data->flag_fwd == 1) && (data->fwd != NULL) && (ifwd != -1)) { // save fwd wavefield
      yb = yb_r;
      ye = ye_r;

      kt = zb;
      for(t = tb; t < te; t++){ // Diamond blocking in time
        if(t % 2 == 0){ //swap pointers
          u = v_r; v = u_r; // U->V
        } else{
          u = u_r; v = v_r; // V->U
        }

        if (t >= t_dim) {
          for(k = kt; k < kt+ze-zb; k ++) { // wavefront loop (Z direction)
            if( ((k-stencilr)/th_nz) % th_z == tid_z ) {
              for(j = yb; j < ye; j++) { // Y
                // left most and right most
                if (((b_inc != 0) && (yb <= j) && (j < yb + stencilr)) ||
                    ((e_inc != 0) && (ye - stencilr <= j) && (j < ye))) {
                  ux = &(   u[1ULL*k*nnxy + j*nnx]);
                  wx = &(data->fwd[1ULL*data->groupsize*groupid + (j-data->gfwd_y0[groupid])*nnx*nnz+k*nnx]);
                  for (i = ib; i < ie; i++) {
                    wx[i] = ux[i];
                  }

                  if (k == data->src_depth) {
                    if ((k == data->src_z) && (j == data->src_y) && (data->src_x >= ib) && (data->src_x < ie)) {
                      wx[data->src_x] -= data->source[t0+(t-tb)];
                    }
                  }
                }
              }
            }
          }
        }

        if(t < t_dim) { // lower half of the diamond
          yb -= b_inc;
          ye += e_inc;
        } else { // upper half of the diamond
          yb += b_inc;
          ye -= e_inc;
        }

        kt -= stencilr;
      }
    }
  } // parallel region
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

void wave_tb_init(tb_t* ctx,
                  sismap_t* s,
                  parser *p) {

  // get thread info
  ctx->thread_group_size = parser_get_int(p, "tb_thread_group_size");
  ctx->num_thread_groups = parser_get_int(p, "tb_nb_thread_groups");
  ctx->num_threads = ctx->thread_group_size * ctx->num_thread_groups;

  ctx->th_x = parser_get_int(p, "tb_th_x");
  ctx->th_y = parser_get_int(p, "tb_th_y");
  ctx->th_z = parser_get_int(p, "tb_th_z");

  ctx->fwd_steps = parser_get_int(p, "fwd_steps");

  if (ctx->th_x * ctx->th_y * ctx->th_z != ctx->thread_group_size) {
    ERR_MSG("thread group (%d, %d, %d) is not match with its size (%d)\n",
            ctx->th_x, ctx->th_y, ctx->th_z , ctx->thread_group_size);
  }

  ctx->affinity_file = parser_get_string(p, "tb_affinity");

  ctx->mode = s->mode;
  switch (ctx->mode) {
    case 1:
    ctx->kernel_spatial_blocking = kernel_spatial_blocking_fuse_mode;
    ctx->kernel_tiling_blocking  = kernel_tiling_blocking_fuse_mode;
    break;

    case 2:
    ctx->kernel_spatial_blocking = kernel_spatial_blocking_separate_mode;
    ctx->kernel_tiling_blocking  = kernel_tiling_blocking_separate_mode;
    break;

    case 3:
    ctx->kernel_spatial_blocking = kernel_spatial_blocking_separate_mode_io;
    ctx->kernel_tiling_blocking  = kernel_tiling_blocking_separate_mode_io;
    break;

    default:
    ERR_MSG("tb_mode (%d) is unknown\n",ctx->mode);
  }


  // set thread affinity
  cpu_bind(ctx);

  // blocking setting
  ctx->t_dim  = parser_get_int(p, "tb_t_dim");
  if (ctx->t_dim %2 == 0) {
    ERR_MSG("t_dim (%d) should be odd\n",ctx->t_dim);
  }

  ctx->num_wf = parser_get_int(p, "tb_num_wf");


  // halo
  if ((s->sx != 4) || (s->sy != 4) || (s->sz != 4)) {
    ERR_MSG("sx (%d),sy (%d) and sz (%d) should be equal to 4\n",s->sx,s->sy,s->sz);
  }

  ctx->r = 4;

  // stride
  ctx->nnx  = s->dimx + 2 * ctx->r;
  ctx->nny  = s->dimy + 2 * ctx->r;
  ctx->nnz  = s->dimz + 2 * ctx->r;

  // grid
  ctx->stencilx = s->dimx;
  ctx->stencily = s->dimy;
  ctx->stencilz = s->dimz;

  // coef
  ctx->coefx = s->coefx;
  ctx->coefy = s->coefy;
  ctx->coefz = s->coefz;

  ctx->diam_width = (ctx->t_dim+1) * 2 * ctx->r;
  if (ctx->stencily % ctx->diam_width != 0) {
    ERR_MSG("stencily (%d) should be multiple of diam_width (%d)\n",ctx->stencily,ctx->diam_width);
  }

  // scheduling
  // round the number of time steps to the nearest valid number
  int remain = (s->time_steps) % ((ctx->t_dim+1)*2);
  if(remain != 0){
    int nt2 = s->time_steps + (ctx->t_dim+1)*2 - remain;
    if(nt2 != s->time_steps){
      MSG("INFO: Modified nt from %03d to %03d for the intra-diamond method\n", s->time_steps,nt2);
      s->time_steps= nt2;
    }
  }

  ctx->time_steps = s->time_steps;
  ctx->t_len =  2*( (s->time_steps)/((ctx->t_dim+1)*2) ) - 1;
  ctx->y_len_l = s->dimy / ctx->diam_width;
  ctx->y_len_r = ctx->y_len_l + 1;

  ctx->t_pos = calloc(ctx->y_len_r,sizeof(int));
  ctx->avail_list = calloc(ctx->y_len_r,sizeof(int));

  // nb stencils main and total
  ctx->nb_stencils_main = 1LL * ctx->t_len * (ctx->t_dim+1)
  * ctx->stencilx * ctx->stencily * ctx->stencilz;

  ctx->nb_stencils_total_fwd = 1LL * ctx->time_steps
  * ctx->stencilx * ctx->stencily * ctx->stencilz;

  ctx->nb_stencils_total_bwd = (ctx->nb_stencils_total_fwd + ctx->nb_stencils_main) / 2;

  // damping
  ctx->dampx = calloc(ctx->nnx,sizeof(float));
  ctx->dampy = calloc(ctx->nny,sizeof(float));
  ctx->dampz = calloc(ctx->nnz,sizeof(float));

  float alpha = 0.2;
  float tabdamp[NDAMP];


  for (int i = 1; i <= NDAMP; i++) {
    tabdamp[NDAMP-i] = exp( -alpha * (1.0*i/NDAMP) * (1.0*i/NDAMP) );
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
    ctx->dampx[ctx->r+i] = tabdamp[i];
    ctx->dampy[ctx->r+i] = tabdamp[i];
    //ctx->dampz[ctx->r+i] = tabdamp[i];

    ctx->dampx[ctx->nnx-1-i-ctx->r] = tabdamp[i];
    ctx->dampy[ctx->nny-1-i-ctx->r] = tabdamp[i];
    ctx->dampz[ctx->nnz-1-i-ctx->r] = tabdamp[i];
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
                           float *u1,
                           unsigned int t) {
    MSG("inside wave_tb_save_lastshot");
    FILE * fd;

  fd = fopen("TB_lastshot_u0", "wb+");
  CHK(fd == NULL, "failed to open snapshot file");
  CHK(fwrite(u0, sizeof(float), s->size, fd) != s->size,
      "failed to write snapshot");
  CHK(fclose(fd)!=0,"failed to close snapshot file");
  if (s->verbose) MSG("... saving last snapshot (size %d)",s->size);

  fd = 0;

  fd = fopen("TB_lastshot_u1", "wb+");
  CHK(fd == NULL, "failed to open snapshot file");
  CHK(fwrite(u1, sizeof(float), s->size, fd) != s->size,
      "failed to write snapshot");
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
  data->ix = calloc(s->rcv_len, sizeof(int));
  data->iy = calloc(s->rcv_len, sizeof(int));

  data->src_depth = -1;
  data->rcv_depth = -1;
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

  data->src_x = data->src_idx % (s->dimx + 2* s->sx);
  data->src_y = data->src_idx / (s->dimx + 2* s->sx);
  data->src_z = data->src_depth;
}

void wave_tb_data_unset_src(tb_data_t * data) {
  data->src_depth = -1;
}

void wave_tb_data_set_rcv(tb_data_t * data,
                          sismap_t *s,
                          float* sismos) {
  data->rcv       = s->rcv;
  data->rcv_depth = s->rcv_depth + 4;
  data->rcv_len   = s->rcv_len;
  data->sismos    = sismos;

//  printf("set rcv stride : %d\n",s->dimx + 2 * s->sx);

  for (int i = 0; i < s->rcv_len; i++) {
    data->ix[i] = s->rcv[i] % (s->dimx + 2 * s->sx);
    data->iy[i] = s->rcv[i] / (s->dimx + 2 * s->sx);
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

static void intra_diamond_mwd_comp(tb_t * ctx,
                                   tb_data_t *data,
                                   tb_timer_t* timer,
                                   float * restrict u0,
                                   float * restrict v0,
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

  ctx->kernel_spatial_blocking(ctx->nnx, ctx->nny, ctx->nnz,
                               ctx->r,                 yb, zb_array,
                               ctx->r + ctx->stencilx, ye, ze_array,
                               ctx->coefx, ctx->coefy, ctx->coefz,
                               ctx->dampx, ctx->dampy, ctx->dampz,
                               u0, v0, roc2,
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

  ctx->kernel_tiling_blocking(ctx->nnx, ctx->nny, ctx->nnz,
                              ctx->r,                 yb, zb,
                              ctx->r + ctx->stencilx, ye, ze,
                              ctx->coefx, ctx->coefy, ctx->coefz,
                              ctx->dampx, ctx->dampy, ctx->dampz,
                              u0, v0, roc2,
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

  ctx->kernel_spatial_blocking(ctx->nnx, ctx->nny, ctx->nnz,
                               ctx->r,                 yb, zb_array,
                               ctx->r + ctx->stencilx, ye, ze_array,
                               ctx->coefx, ctx->coefy, ctx->coefz,
                               ctx->dampx, ctx->dampy, ctx->dampz,
                               u0, v0, roc2,
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

static inline void intra_diamond_resolve(tb_t *ctx,
                                         tb_data_t *data,
                                         tb_timer_t* timer,
                                         const int y_coord,
                                         const int groupid,
                                         float * restrict u0,
                                         float * restrict v0,
                                         const float * restrict roc2) {
  int yb, ye, b_inc, e_inc, ifwd, iifwd;
  int t0;
  size_t offset,fwdsize;

  double tm0;

  intra_diamond_get_info(ctx, timer, y_coord, ctx->t_pos[y_coord], groupid,
                         &yb, &ye, &b_inc, &e_inc);

  iifwd = -1;

  if (data->flag_fwd == 1) {
    ifwd = ctx->t_pos[y_coord];
    if (ifwd % ctx->fwd_steps == 0) iifwd = ifwd / ctx->fwd_steps;
    t0 = 1 + ctx->t_pos[y_coord]*(ctx->t_dim+1);
  }
  if (data->flag_bwd == 1) {
    ifwd = ctx->t_len - 1 - ctx->t_pos[y_coord];
    if (ifwd % ctx->fwd_steps == 0) iifwd = ifwd / ctx->fwd_steps;
    t0 = ctx->time_steps - 2 - ctx->t_pos[y_coord]*(ctx->t_dim+1);
  }

  if (ctx->mode == 3) {
    if (iifwd != -1) {
      fwd_setup(ctx,&(data->gfwd_y0[groupid]),&fwdsize,&offset,
                y_coord,ctx->t_pos[y_coord],iifwd);

      if ((data->flag_bwd == 1) &&(data->fwd != NULL)) {
        tm0 = wtime();
        fwd_load(ctx,data->gfd[groupid],data->fwd,groupid, offset,fwdsize);
        timer->t_wf_snapio[groupid] += wtime()-tm0;
      }
    }
  }

//  MSG("... t0=: %d\n", t0 );
  intra_diamond_mwd_comp(ctx, data, timer, u0, v0, roc2, yb, ye, b_inc, e_inc,
                         0, ctx->t_dim * 2 + 1, t0, iifwd, groupid);

  if (ctx->mode == 3) {
    if (iifwd != -1) {
      if ((data->flag_fwd == 1) &&(data->fwd != NULL)) {
        tm0 = wtime();
//        fwd_save(ctx,data->gfd[groupid],data->fwd,groupid,  offset,fwdsize);

        if (groupid < ctx->num_thread_groups/2) {
          bwriter_write(data->bwriter0,groupid,&data->fwd[1LL* ctx->nnx * ctx->nnz * ctx->diam_width* groupid],offset,fwdsize);
        } else {
          bwriter_write(data->bwriter1,groupid,&data->fwd[1LL* ctx->nnx * ctx->nnz * ctx->diam_width* groupid],offset,fwdsize);
        }

        timer->t_wf_snapio[groupid] += wtime()-tm0;
      }
    }
  }
}

static inline void dynamic_intra_diamond_prologue(tb_t * ctx,
                                                  tb_data_t *data,
                                                  tb_timer_t* timer,
                                                  float * restrict u0,
                                                  float * restrict v0,
                                                  const float * restrict roc2) {
  int i, yb, ye, ifwd;

  ifwd = -1;

  #pragma omp parallel num_threads(ctx->num_thread_groups)
  {
    int b_inc = ctx->r;
    int e_inc = ctx->r;

    int groupid = 0;
    #if defined(_OPENMP)
    groupid = omp_get_thread_num();
    #endif

    #pragma omp for schedule(dynamic) private(i,yb,ye)
    for(i = 0; i < ctx->y_len_l; i++){
      yb = ctx->r + i * ctx->diam_width + ctx->r;
      ye = yb + ctx->diam_width - 2*ctx->r;
      intra_diamond_mwd_comp(ctx, data, timer, u0, v0, roc2, yb, ye, b_inc, e_inc,
                             ctx->t_dim + 1, ctx->t_dim * 2 + 1, 1,
                             ifwd, groupid);
    }
  }
}

static inline void dynamic_intra_diamond_prologue_full(tb_t * ctx,
                                                       tb_data_t *data,
                                                       tb_timer_t* timer,
                                                       float * restrict u0,
                                                       float * restrict v0,
                                                       const float * restrict roc2) {
  int i, yb, ye;
  int ifwd;

  ifwd = -1;

  #pragma omp parallel num_threads(ctx->num_thread_groups)
  {
    int b_inc = ctx->r;
    int e_inc = ctx->r;

    int groupid = 0;
    #if defined(_OPENMP)
    groupid = omp_get_thread_num();
    #endif

    #pragma omp for schedule(dynamic) private(i,yb,ye)
    for(i = 0; i < ctx->y_len_l; i++){
      yb = ctx->r + i * ctx->diam_width;
      ye = yb + ctx->diam_width;
      intra_diamond_mwd_comp(ctx, data, timer, u0, v0, roc2, yb, ye, b_inc, e_inc,
                             ctx->t_dim, ctx->t_dim * 2 + 1, ctx->time_steps-1,
                             ifwd, groupid);
    }
  }
}


static inline void dynamic_intra_diamond_epilogue(tb_t* ctx,
                                                  tb_data_t *data,
                                                  tb_timer_t* timer,
                                                  float * restrict u0,
                                                  float * restrict v0,
                                                  const float * restrict roc2) {
  int i, yb, ye, ifwd;
  int t0;
  if (data->flag_fwd == 1) {
    ifwd = -1;
    t0 = ctx->t_len*(ctx->t_dim+1) + 1;
  }

  if (data->flag_bwd == 1) {
    ifwd = -1;
    t0 = ctx->t_dim - 1;
  }

  #pragma omp parallel num_threads(ctx->num_thread_groups)
  {
    int b_inc = ctx->r;
    int e_inc = ctx->r;
    int yb_r = ctx->r + ctx->diam_width/2 - ctx->r;
    int ye_r = yb_r + 2 * ctx->r;
    int groupid = 0;
    #if defined(_OPENMP)
    groupid = omp_get_thread_num();
    #endif

    #pragma omp for schedule(dynamic) private(i,yb,ye)
    for(i = 0; i < ctx->y_len_l; i++){
      yb = yb_r + i * ctx->diam_width;
      ye = ye_r + i * ctx->diam_width;
      intra_diamond_mwd_comp(ctx, data, timer, u0, v0, roc2, yb, ye, b_inc, e_inc,
                             0, ctx->t_dim + 1, t0, ifwd, groupid);
    }
  }
}

#define T_POS_L(y) (ctx->t_pos[(((y)+(ctx->y_len_l)) % (ctx->y_len_l))])
#define T_POS_R(y) (ctx->t_pos[(((y)+(ctx->y_len_r)) % (ctx->y_len_r))])

static inline void update_state(tb_t *ctx,
                                int y_coord){
  int sh;
  ctx->t_pos[y_coord]++; // advance the current tile in time

  if(ctx->t_pos[y_coord] % 2 == 0){ // right row case
    // add the current diamond to the ready queue if dependencies are satisfied
    if(ctx->t_pos[y_coord] < ctx->t_len){
      // if left-half diamond, no dependencies. Add to the list
      if(y_coord == ctx->y_len_l - 1){
        ctx->avail_list[wave_tb_head % ctx->y_len_r] = y_coord;
        wave_tb_head ++;
      } else if(T_POS_R(y_coord+1) >= ctx->t_pos[y_coord]) {
        //the reset have the same circular dependence (except the right-half diamond) if:
        // 1) the current tile did not reach the end of the temporal dimension
        // 2) the right neighbor is at least at the same time step
        ctx->avail_list[wave_tb_head % ctx->y_len_r] = y_coord;
        wave_tb_head ++;
      }
    } // check validity in range of temporal dimension

    // add the dependent diamond to the ready queue if other dependencies are satisfied:
    if (T_POS_R(y_coord-1) < ctx->t_len){
      // add the right-half diamond automatically when the left most diamond is updated
      if(y_coord == 0){ // no dependencies. Add to the list
        ctx->t_pos[ctx->y_len_r-1]++; // advance the right-half diamond in time
        ctx->avail_list[wave_tb_head % ctx->y_len_r] = ctx->y_len_r - 1;
        wave_tb_head ++;
      }
      else if(T_POS_R(y_coord-1) == ctx->t_pos[y_coord]) {
        // 1) the neighbor did not reach the end of the temporal dimension
        // 2) the left neighbor is at the same time step
        // 3) is not the right-half diamond
        ctx->avail_list[wave_tb_head % ctx->y_len_r] = (y_coord - 1 + ctx->y_len_r) % ctx->y_len_r; // add the dependent neighbor to the list if the dependency is satisfied
        wave_tb_head ++;
      }
    } // check validity in temporal dimension
  } //end right row case

  else if(ctx->t_pos[y_coord] % 2 == 1){ // left row
      // add the current diamond to the ready queue if dependencies are satisfied:
    if( (T_POS_R(y_coord-1) >= ctx->t_pos[y_coord]) && (ctx->t_pos[y_coord] < ctx->t_len)  && (y_coord != ctx->y_len_r-1) ) {
      // 1) the left neighbor is at least at the same time step
      // 2) the current diamond did not reach the end of the temporal dimension
      // 3) is not the right-half diamond
      ctx->avail_list[wave_tb_head % ctx->y_len_r] = y_coord;
      wave_tb_head ++;
    }

    // add the dependent diamond to the ready queue if other dependencies are satisfied:
    if( (T_POS_R(y_coord+1) == ctx->t_pos[y_coord]) && (T_POS_R(y_coord+1) < ctx->t_len) && (y_coord != ctx->y_len_l-1) ) {
      // 1) the right neighbor is at the same time step
      // 2) the neighbor did not reach the end of the temporal dimension
      // 3) is not the right most diamond in space
      ctx->avail_list[wave_tb_head % ctx->y_len_r] = (y_coord + 1 + ctx->y_len_r) % ctx->y_len_r; // add the dependent neighbor to the list if the dependency is satisfied
      wave_tb_head ++;
    }
  } // end left row case
}

void dynamic_intra_diamond_mainloop(tb_t* ctx,
                                    tb_data_t *data,
                                    tb_timer_t* timer,
                                    float * restrict u0,
                                    float * restrict v0,
                                    const float * restrict roc2) {
  int not_complete, th_y_coord, i, groupid;
  double t1;
  uint64_t il;
  uint64_t num_diamonds = ctx->y_len_l *  (ctx->t_len-1)/2
  + ctx->y_len_r * ((ctx->t_len-1)/2 +1);

  for(i = 0; i < ctx->y_len_r; i++){
    ctx->avail_list[i] = i;
  }

  wave_tb_head = ctx->y_len_r;
  wave_tb_tail = 0;

  for(i = 0; i < ctx->y_len_r; i++) {
    ctx->t_pos[i] = 0;
  }

#pragma omp parallel num_threads(ctx->num_thread_groups) \
  shared(wave_tb_head, wave_tb_tail) private(groupid)
  {

   groupid = 0;
#if defined(_OPENMP)
   groupid =  omp_get_thread_num();

   int gtid = groupid * ctx->thread_group_size;
   int err = sched_setaffinity(0, ctx->setsize, ctx->bind_masks[gtid]);
   if(err == -1) printf("WARNING: Could not set CPU Affinity (%d %d)\n",groupid,gtid);

#endif


#pragma omp for schedule(dynamic) private(il, th_y_coord, not_complete,t1)//shared(head,tail)
   for (il = 0; il < num_diamonds; il++){

    not_complete = 1;
    th_y_coord = -1;
    while(not_complete)
    {
      t1 = wtime();
        while(wave_tb_head - wave_tb_tail < 1); // spin-wait for available tasks
        timer->t_group_wait[groupid] += wtime() - t1;

#pragma omp critical// (consumer)
        {
#pragma omp flush (wave_tb_head, wave_tb_tail)
          if(wave_tb_head - wave_tb_tail > 0){ // make sure there is still available work
            th_y_coord = ctx->avail_list[wave_tb_tail % ctx->y_len_r]; //acquire task
            wave_tb_tail ++;
          }
        }
        if(th_y_coord >= 0){
          intra_diamond_resolve(ctx, data, timer, th_y_coord, groupid, u0, v0, roc2);
#pragma omp critical// (producer)
          {
#pragma omp flush (wave_tb_head)
            update_state(ctx, th_y_coord);
          }
          not_complete = 0;
        }
      }
    }

  }
}

void wave_tb_forward(tb_t* ctx,
                     tb_data_t * data,
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict v0,
                     const float * restrict roc2) {

  double t1,t2,t3,t4;

  if (data->src_depth !=-1) {
    u0[1ULL*data->src_depth* ctx->nnx * ctx->nny + data->src_idx] += data->source[0];
  }

  t1 = wtime();
  dynamic_intra_diamond_prologue(ctx, data, timer, u0, v0, roc2);
  t2 = wtime();
  dynamic_intra_diamond_mainloop(ctx, data, timer, u0, v0, roc2);
  t3 = wtime();
  dynamic_intra_diamond_epilogue(ctx, data, timer, u0, v0, roc2);
  t4 = wtime();

  timer->ts_main   += (t3-t2);
  timer->ts_others += (t2-t1) + (t4-t3);
  timer->total     += timer->ts_main + timer->ts_others;
}


void wave_tb_backward(tb_t* ctx,
                      tb_data_t * data,
                      tb_timer_t* timer,
                      float * restrict u0,
                      float * restrict v0,
                      const float * restrict roc2) {

  double t1,t2,t3,t4;

  t1 = wtime();
  dynamic_intra_diamond_prologue_full(ctx, data, timer, u0, v0, roc2);
  t2 = wtime();
  dynamic_intra_diamond_mainloop(ctx, data, timer, u0, v0, roc2);
  t3 = wtime();
//  dynamic_intra_diamond_epilogue(ctx, data, u0, v0, roc2);
  t4 = wtime();

  timer->ts_main   += (t3-t2);
  timer->ts_others += (t2-t1) + (t4-t3);
  timer->total     += timer->ts_main + timer->ts_others;
}


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
  MSG("Total:        %f (s) -%06.2f%%",timer->total,        timer->total/timer->total*100.0);
  MSG("mainloop:     %f (s) -%06.2f%%",timer->ts_main,      timer->ts_main/timer->total*100.0);
  MSG("pro/epilogue: %f (s) -%06.2f%%",timer->ts_others,    timer->ts_others/timer->total*100.0);
  MSG("-------------------------------------------");

  MSG("Speed info:");
  MSG("Total: %f Mstencils/s",nb_stencils_total/1e6/timer->total);
  MSG("Main:  %f Mstencils/s",nb_stencils_main /1e6/timer->ts_main);
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
