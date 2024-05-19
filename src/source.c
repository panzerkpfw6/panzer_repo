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
/// @file src/source.c
/// @brief This file contains the implementation of the source simulator.
///
#include <math.h>
#include <stdio.h>
#include <simwave/config.h>
#include <simwave/macros.h>
#include <simwave/source.h>
#include <string.h>
#include <errno.h>

#define SOURCE_BASE "src"
//#define __DUMP_SOURCE

void source_ricker_wavelet(sismap_t *s, float* source) {
  char tmp[256];
  #ifdef __DUMP_SOURCE
  sprintf(tmp, "%s/%s.txt", OUTDIR, SOURCE_BASE);
  FILE *fd = fopen(tmp, "w");
  #endif // __DUMP_SOURCE

  #ifdef __USE_SIGMA
  int it;
  float sigma, tau, t, scale;
  sigma = 0.6*s->fmax;
  tau   = 1.0;
  scale = 8.0;
  for(it=0; it < s->time_steps; ++it) {
    t = s->dt*(float)(it);
    source[it] =
      -2.0*scale*sigma*
      (sigma-2.0*sigma*scale*(sigma*t-tau)*(sigma*t-tau))*
    exp(-scale*(sigma*t-tau)*(sigma*t-tau));
    #ifdef __DUMP_SOURCE
    fprintf(fd, "%d %f\n", it, source[it]);
    #endif // __DUMP_SOURCE
  }
  #else
  int it;
  float t1, t0;
  float PI = 4.0f * atan(1.0f);
  t0=1.0/s->fmax;

  for (it=0; it<s->time_steps; it++) {
    t1=it*s->dt;
    source[it] = exp(-PI*PI*s->fmax*s->fmax*(t1-t0)*(t1-t0))*
            (1.0-2.*PI*PI*s->fmax*s->fmax*(t1-t0)*(t1-t0));
    #ifdef __DUMP_SOURCE
    fprintf(fd, "%d %f\n", it, source[it]);
    #endif // __DUMP_SOURCE
  }
  #endif // __USE_SIGMA
  #ifdef __DUMP_SOURCE
  fclose(fd);
  sprintf(tmp, "%s/%s.raw", OUTDIR, SOURCE_BASE);
  fd = fopen(tmp, "wb");
  ERR_IF(fd == NULL, "failed to open the source file for dumping");
  CHK(fwrite(source, sizeof(float),
             s->time_steps, fd) != s->time_steps,
             "failed to write in source file");
  fclose(fd);
  #endif // __DUMP_SOURCE
}




  // if (Freq == 0.0) Freq = 30.0;
  // if (dt == 0.0) dt = 0.004;
  // Bpar = sqrt(6.0) / (PI * Freq);
  // N = ceil(1.35 * Bpar / dt);
  // Np1 = N;
  // *Npoint = 2 * N + 1;

  // Amp = alloc1float(*Npoint);

  // Amp[Np1] = 1.0;

  // for (i = 1; i <= N; i++) {
  //   t = dt * (float)i;
  //   u = 2.0 * sqrt(6.0) * t / Bpar;
  //   Amp[Np1 + i] = Amp[Np1 - i] = 0.5 * (2.0 - u * u) * exp(-u * u / 4.0);
  // }
  // return Amp;
