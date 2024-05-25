#ifndef __SIMWAVE_PML_H_
#define __SIMWAVE_PML_H_
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
/// @file simwave/pml.h
///
/// Contains all the routines necessary for the boundary
/// conditions (using the <b>Perfectly Matched Layer</b> method) computations.
///
#include <simwave/sismap.h>

#define PML_TAPER 0

/// @brief Calculate the damping factors.
void pml_compute_coefs(sismap_t *s, float* pml_array);

#endif // __SIMWAVE_PML_H_
