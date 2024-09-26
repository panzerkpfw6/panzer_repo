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
/// @file src/pml.c
/// @brief This file contains the PML coefficients calculations.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/pml.h>
#include <omp.h>

#define ETA(z,y,x) ( eta[(2+s->dimx)*((2+s->dimy)*(z+1) + (y+1)) + (x+1)])
//#define PML_SCALE_FACTOR(taper, f) 3.*f*log(1000.)/(2.*taper);
#define PML_SCALE_FACTOR(taper, f) 0
/// @brief compute PML parameters
void pml_compute_coefs(sismap_t *s, float* eta) {
  int    x,  y,  z;
  float vx, vy, vz;
  float scalex, scaley, scalez;
  scalex = PML_SCALE_FACTOR(PML_TAPER, s->fmax);
  scaley = PML_SCALE_FACTOR(PML_TAPER, s->fmax);
  scalez = PML_SCALE_FACTOR(PML_TAPER, s->fmax);
  for(z = 0; z < s->dimz+2; ++z) {
    for(y = 0; y < s->dimy+2; ++y) {
      for(x = 0; x < s->dimx+2; ++x) {
        /// free surface on the top.
        if(z >= ((s->dimz+2) - s->pmlz))
          vz = pow(1.*(z - ((s->dimz+2) - s->pmlz - 1))/s->pmlz, 2)*scalez;
        else
          vz = 0.;
        if(y < s->pmly)
          vy = pow(1.*(s->pmly - y)/s->pmly, 2)*scaley;
        else if(y >= ((s->dimy+2) - s->pmly))
          vy = pow(1.*(y - ((s->dimy+2) - s->pmly - 1))/s->pmly, 2)*scaley;
        else
          vy = 0.;
        if(x < s->pmlx)
          vx = pow(1.*(s->pmlx - x)/s->pmlx, 2)*scalex;
        else if(x >= ((s->dimx+2) - s->pmlx))
          vx = pow(1.*(x - ((s->dimx+2) - s->pmlx - 1))/s->pmlx, 2)*scalex;
        else
          vx = 0.;
        ETA(z-1,y-1,x-1) = (MAX(vx, MAX(vy, vz)))*s->dt;
      }
    }
  }
  #ifdef __DUMP_PML
  char tmp[512];
  sprintf(tmp, "mkdir -p %s", OUTDIR);
  CHK(system(tmp), "failed to create output directory for snapshots");
  sprintf(tmp, "%s/pml.raw", OUTDIR);
  FILE *fd  = fopen(tmp, "wb");
  CHK(fwrite(eta,
            (s->dimx+2)*(s->dimy+2)*(s->dimz+2)*sizeof(float), 1, fd)!=1,
      "failed to write the pml file");
  fclose(fd);
  #endif // __DUMP_PML
}
