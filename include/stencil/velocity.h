#ifndef __STENCIL_VELOCITY_H_
#define __STENCIL_VELOCITY_H_
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
/// @file velocity.h
/// @brief A velocity model generator.
///
/// Generates a phase velocity model for the simulation.
///
#include <stencil/config.h>
#include <stencil/sismap.h>

/// @brief Generates a three dimensional velocity model in m/s
/// @param vtab is an array that will contain the model terms
/// @param dimx is the width of the domain
/// @param dimy is the height of the domain
/// @param dimz is the depth of the domain
/// @param vmin is the minimum velocity required
/// @param vmax is the maximum velocity tolerated
/// @param layers is the maximum velocity layers number
void velocity_generate_model(sismap_t *s, float* vtab, unsigned int layers);

void velocity_query_model(sismap_t *s);

void fill_coef_matrix(sismap_t *s,float *vtab,float *dens);

void dump_vel(sismap_t *s,float *vtab,float *dens);

void dump_coef(sismap_t *s, float *vtab);

void velocity_load_model(sismap_t *s, float *vtab);

void velocity_load_model_2d(sismap_t *s, float *vtab);

void velocity_load_model_3d(sismap_t *s, float *vtab);

void velocity_const_model2(sismap_t *s, float *vtab);

void velocity_2layer_model(sismap_t *s, float *vtab, unsigned int layers);

void velocity_load_model(sismap_t *s, float *vtab);

void velocity_load_salt3d(sismap_t *s, float *vtab);

void density_const_model(sismap_t *s,float *dens,float *inv_rho);

#endif // __STENCIL_VELOCITY_H_
