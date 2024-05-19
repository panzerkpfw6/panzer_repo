///
/// @copyright Copyright 2017- Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// \b simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author Issam Said
/// @file src/interp.c
/// @brief interpretation routines.
///
#include <stdio.h>
#include <string.h>
#include <simwave/interp.h>

float linearinterp(float a, float b, float t) {
	return a*(1-t)+b*t;
}

float bilinearinterp(float c00, float c01, 
	                          float c10, float c11, float tx, float ty) {
 	float a = linearinterp(c00, c01, tx); 
  float b = linearinterp(c10, c11, tx); 
  return linearinterp(a, b, ty); 
}

float trilinearinterp(float c000, float c001, 
	                           float c010, float c011, 
	                           float c100, float c101,
	                           float c110, float c111, 
	                           float tx, float ty, float tz) {
 	float e = bilinearinterp(c000, c001, c010, c011, tx, ty); 
  float f = bilinearinterp(c100, c101, c110, c111, tx, ty); 
  return linearinterp(e, f, tz); 
}

