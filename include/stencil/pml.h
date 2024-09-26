#ifndef __STENCIL_PML_H_
#define __STENCIL_PML_H_
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
/// @file stencil/pml.h
///
/// Contains all the routines necessary for the boundary
/// conditions (using the <b>Perfectly Matched Layer</b> method) computations.
///
#include <stencil/sismap.h>

#define PML_TAPER 0

/// @brief Calculate the damping factors.
void pml_compute_coefs(sismap_t *s, float* pml_array);

#endif // __STENCIL_PML_H_
