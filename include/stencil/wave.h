#ifndef __STENCIL_WAVE_H_
#define __STENCIL_WAVE_H_
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
/// @file stencil/wave.h
///
/// Contains all the routines necessary for a wave simulation with boundary
/// conditions (using the <b>Perfectly Matched Layer</b> method) on CPU.
///
#include <stdio.h>
#include <stdbool.h>

#include <stencil/config.h>
#include <stencil/sismap.h>
#include <stencil/shot.h>

/// SB cache blocking parameters.We prefer not to hardcode it. def:10,22,9999

void array_openmp_init(float* u,sismap_t*s);
void array_openmp_inner_init(float* u,sismap_t*s);
void wave_init_numerics(sismap_t *s);

/// @brief setup the velocity and compute geometries.
void wave_init_dimensions(sismap_t *s);
void wave_init_damp(sismap_t *s);

void wave_init_acquisition(sismap_t *s);

/// @brief Creates a wave descriptor for the CPU
/// @return a pointer to the created wave
void wave_init_buffers(sismap_t *s, float* phi, float* eta);

/// @brief Deletes a wave descriptor
/// @param w is a pointer to the wave descriptor to be deleted
void wave_release(sismap_t *s);

/// @brief Prints informations about a wave descriptor
/// @param w is a pointer to the wave descriptor to be printed
void wave_print(sismap_t *s);

/// @brief Computes the pressure on each element of the domain
/// @param w is a pointer to the wave descriptor to be updated
void wave_update_fields(sismap_t *s, float *,
	                      float *, float*, float*, float*);
void wave_update_fields_1st(sismap_t *s,
                            float *u0, float *vx,float *vy,float *vz,
                            float *roc2, float *phi, float *eta);
void wave_update_fields_block(sismap_t *s, float *,
                        float *, float*, float*, float*);
void wave_update_fields_block_bis_orig(sismap_t *s,
                                       float *restrict u0,
                                       float *restrict u1,
                                       float *restrict roc2,
                                       float *restrict phi,
                                       float *restrict eta);
void wave_update_fields_block_bis(sismap_t *s,
                                  float* restrict u0,
                                  float* restrict u1,
                                  float* restrict roc2,
                                  float* restrict phi,
                                  float* restrict eta);
void wave_update_fields_block_1st(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict vx,
                                  float *restrict vy,
                                  float *restrict vz,
                                  const float *restrict roc2,
								  const float *restrict inv_rho);

/// @brief Adds an impulse to the grid point situated at the
/// source location
/// @param w is a pointer to the wave descriptor to be updated
/// @param time_step is the current time step
void wave_update_source(sismap_t *s, shot_t *, float*, float);

/// @brief Extract the sismograms on the receivers positions.
/// @param w is a pointer to the wave descriptor to be updated
/// @param time_step is the current time step
void wave_extract_sismos(sismap_t *s, float*,
												 unsigned int time_step, float*);

/// @brief Inject the sismograms on the receivers positions.
/// @param w is a pointer to the wave descriptor to be updated
/// @param time_step is the current time step
void wave_inject_sismos(sismap_t *s, float*,
	                      unsigned int time_step, float*);

/// @brief Switches the @ref wave_t::u0 and @ref wave_t::u1 pointers
/// @param w is a pointer to the wave descriptor to be updated
#define WAVE_SWAP_POINTERS(u0, u1) \
  float* tmp = u0;                 \
  u0 = u1;                         \
  u1 = tmp

/// @brief Saves a snapshot of the wave fields
/// @param w is a pointer to the wave descriptor to be saved
void wave_save_snapshot(sismap_t *s, shot_t *shot, float*, unsigned int t);

/// @brief Reads a snapshot of the wave fields
/// @param w is a pointer to the wave descriptor to be saved
void wave_read_snapshot(sismap_t *s, shot_t *shot, float*, unsigned int t);

/// @brief Implement the imaging condition
/// @param w is a pointer to the wave descriptor to be saved
void wave_image_condition(sismap_t *s,
	                        float*, float *, float *, float*, unsigned int t);

void wave_image_condition_block(sismap_t *s,
                                float*, float *, float *, float*, unsigned int t);


/// @brief Saves sismos of the wave fields
/// @param w is a pointer to the wave descriptor to be saved
void wave_save_sismos(sismap_t *s, shot_t *shot, float*);

/// @brief Reads sismos of the wave fields
/// @param w is a pointer to the wave descriptor to be saved
void wave_read_sismos(sismap_t *s, shot_t *shot, float*);

/// @brief Saves the final image after migration.
/// @param w is a pointer to the wave descriptor to be saved
void wave_save_image(sismap_t *s, float*, char*);

void wave_save_fwd_dbg(sismap_t* s, shot_t *shot, float* u1, unsigned int t);
void wave_save_bwd_dbg(sismap_t* s, shot_t *shot, float* u1, unsigned int t);

void wave_save_img(sismap_t* s, shot_t *shot, float* , float*);

//int  wave_check_fields(float *tab, size_t len);

void wave_min_max(char* str, float *tab, size_t len);

#endif // __STENCIL_WAVE_H_
