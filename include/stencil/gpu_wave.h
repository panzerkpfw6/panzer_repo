#ifndef __STENCIL_GPU_WAVE_H_
#define __STENCIL_GPU_WAVE_H_
///
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
/// @file stencil/gpu_wave.h
/// @brief A \em reflection-less wave (\b wave) handler on the GPU.
///
/// Contains all the routines necessary for a wave simulation, with 
/// application of absorbing boundary conditions
/// (using the <b>Perfectly Matched Layer</b> method)
/// on the GPU side. We call the new wave \b gpu_wave for sake of clarity. 
///
#include <stencil/sismap.h>
#include <stencil/shot.h>

/// @brief Computes the pressure on each element of the domain
/// using an OpenCL kernel on the GPU
/// @param gpu_w is a pointer to the wave descriptor to be updated 
void gpu_wave_update_fields(sismap_t *s, float*, 
                            float*, float*, float*, float* );

/// @brief Adds an impulse to the grid point situated at the
/// source location using an OpenCL kernel on the GPUn 
/// @param gpu_w is a pointer to the wave descriptor to be updated 
/// @param sterm is the source term to add at the current time step.
void gpu_wave_update_source(sismap_t *s, shot_t *shot, float *, float);

void gpu_wave_image_condition(sismap_t* s, float *d_u1, 
                              float *d_fwd,
                              float *d_img, float *d_ilm, unsigned int t);
void gpu_wave_image_gather(sismap_t* s, float *d_img_shot, 
                           float *d_ilm_shot, float *d_img);
/// @brief Switches the @ref gpu_wave_t::d_u0 and @ref gpu_wave_t::d_u1 pointers 
/// @param gpu_w is a pointer to the wave descriptor to be updated 
#define GPU_WAVE_SWAP_POINTERS(u0, u1) \
  float* tmp = u0;                     \
  u0 = u1;                             \
  u1 = tmp       

/// @brief Retrieves data from the GPU and saves a snapshot of the 
/// wave fields. This can be improved by using asynchronous communications
/// @param gpu_w is a pointer to the wave descriptor to be save
void gpu_wave_save_snapshot(sismap_t *s, shot_t *shot,
                            float *, float *, unsigned int t);

void gpu_wave_read_snapshot(sismap_t* s, shot_t *shot, 
                            float *d_fwd, float *fwd, unsigned int t);

/// @brief Retrieves the sismos from the GPU and saves them on disk.  
void gpu_wave_save_sismos(sismap_t* s, shot_t *shot, 
	                        float *d_sismos, float *sismos);
/// @brief read the sismos and send to the GPU.  
void gpu_wave_read_sismos(sismap_t* s, shot_t *shot, 
	                        float *d_sismos, float *sismos);

void gpu_wave_extract_sismos(sismap_t* s, 
                             float *d_u1, unsigned int t, float *d_sismos);

void gpu_wave_inject_sismos(sismap_t* s, 
                            float *d_u1, unsigned int t, float *d_sismos);

/// @brief Creates a GPU resources descriptor
/// @param p is the parser that may contain a list of additional options
void gpu_wave_set(int);

void gpu_wave_unset();

void gpu_wave_init(sismap_t *s);

/// @brief Releases all the CUDA resources properly
void gpu_wave_release(sismap_t *s);

/// @brief Prints informations about the used OpenCL resources
void gpu_wave_info(sismap_t *s);

void gpu_wave_save_fwd_dbg(sismap_t* s, 
                           shot_t *shot, 
                           float *d_u1, float* u1, unsigned int t);
void gpu_wave_save_bwd_dbg(sismap_t* s, 
                           shot_t *shot, 
                           float *d_u1, float* u1, unsigned int t);
void gpu_wave_save_img(sismap_t* s, 
                       shot_t *shot, 
                       float *d_img, float *d_ilm, float *tmp);
#endif //  __STENCIL_GPU_WAVE_H_
