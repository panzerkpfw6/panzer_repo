#ifndef __SIMWAVE_SISMAP_H_
#define __SIMWAVE_SISMAP_H_
///
/// @copyright Copyright 2017- Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// @b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// @b simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author Issam Said
/// @file simwave/sismap.h
///
/// Contains all the routines necessary for a wave simulation with boundary
/// conditions (using the <b>Perfectly Matched Layer</b> method) on CPU.
///
#include <stdio.h>
#include <stdbool.h>
//#include <cuda_runtime.h>
#include <simwave/config.h>
#include <simwave/shot.h>

/// @brief A data structure that contains fields related to the physics.
typedef struct __sismap_t {
  ///
  /// velocity meta-data:
  ///
  char *vel_file;
  unsigned int vel_dimx, vel_dimy, vel_dimz;
  unsigned int dcdp, dline, ddepth;
  ///
  /// acquisition meta-data:
  ///
  unsigned int drcv, dshot;
  unsigned int src_depth, rcv_depth;
  unsigned int rcv_len;
  unsigned int *rcv;
  int nb_shots;
  int    first;
  int     last;
  shot_t **shots;
  ///
  /// simulation meta-data:
  ///
  float courant_number, vmin, vmax, fmax;
  float hdx2, hdy2, hdz2;
  float dt;
  float lambda;
  float cfl;
  ///
  /// compute grids meta-data:
  ///
  unsigned int img_dimx, img_dimy, img_dimz;
  unsigned int dimx, dimy, dimz;
  unsigned int dx, dy, dz;
  unsigned int dtrpx, dtrpy, dtrpz;
  unsigned int sx, sy, sz;
  unsigned int pmlx, pmly, pmlz;
  unsigned int snap_idx;
  unsigned int nb_snap;
  unsigned int nyquist_sampling;
  unsigned int next_snap;
  int time_steps;
  size_t size;
  size_t size_eff;
  size_t size_img;
  ///
  /// read-only arrays that won't changed during the simulation:
  ///
  /// the stencil coefficients.
  float *coefx, *coefy, *coefz;

  // damping coefficients
  float *dampx, *dampy, *dampz;
  ///
  /// GPUs meta-data:
  ///
  int device;
  //dim3 g, g_img, o, l;
  float *d_coefx, *d_coefy, *d_coefz;
  /// receivers.
  unsigned int *d_rcv;
  ///
  /// control flags:
  ///
  bool modeling, verbose, cpu, check, dim2;

  int mode;
} sismap_t;

#endif // __SIMWAVE_WAVE_H_
