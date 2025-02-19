#ifndef __STENCIL_SOURCE_H_
#define __STENCIL_SOURCE_H_

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
/// @file stencil/source.h
/// @brief A wave source generator.
///
/// Create a periodic signal that simulates waves pressure impulses.
///
#include <stencil/config.h>
#include <stencil/sismap.h>

/// @brief Creates a Ricker Wavelet to be used as a wave source.
/// @see http://subsurfwiki.org/wiki/Ricker_wavelet
/// @param source an array that will contains the source discretized terms
/// @param dt is the duration of a time step
/// @param time_steps the number of the simulation time steps
/// @param fmax is the max frequency of the input signal
void source_ricker_wavelet(sismap_t *s, float* source);
void source_ricker_wavelet_1st(sismap_t *s, float *source);
void source_ricker_wavelet_2nd(sismap_t *s, float *source);

#endif // __STENCIL_SOURCE_H_
