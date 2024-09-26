/// @file src/interp.c
/// @brief interpretation routines.
///
#include <stdio.h>
#include <string.h>
#include <stencil/interp.h>

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

