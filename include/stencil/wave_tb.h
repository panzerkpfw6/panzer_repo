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

typedef struct tb_s tb_t;
typedef struct tb_data_s tb_data_t;
typedef struct tb_timer_s tb_timer_t;
int get_ntg(Parameters p);

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

  // estimation
  int64_t nb_stencils_main;
  int64_t nb_stencils_total_fwd;
  int64_t nb_stencils_total_bwd;

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

void wave_tb_free(tb_t* ctx);

void wave_tb_save_lastshot(sismap_t* s,
                           shot_t *shot,
                           float* u0,
                           float * u1,
                           unsigned int t);

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
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict vx,
                     float * restrict vy,
                     float * restrict vz,
                     const float * restrict roc2);

void wave_tb_backward_1st(tb_t* ctx,
                     tb_data_t* data,
                     tb_timer_t* timer,
                     float * restrict u0,
                     float * restrict vx,
                     float * restrict vy,
                     float * restrict vz,
                     const float * restrict roc2);

#endif