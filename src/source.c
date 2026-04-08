///
/// @copyright Copyright . All rights reserved.
/// This file is part of \b stencil.
///
/// @b stencil is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// @b stencil is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b stencil.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author
/// @file src/source.c
/// @brief This file contains the implementation of the source simulator.
///
#include <math.h>
#include <stdio.h>
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/source.h>
#include <string.h>
#include <errno.h>

#define SOURCE_BASE "src"
#define __DUMP_SOURCE

void source_ricker_wavelet(sismap_t *s, float *source) {
    // 2nd order ricker wavelet
    char tmp[256];
#ifdef __DUMP_SOURCE
    sprintf(tmp, "%s/%s.txt", OUTDIR, SOURCE_BASE);
    FILE *fd = fopen(tmp,"w");
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
//    float PI=3.1415926535897 ;
//    t0 = 1.0 / s->fmax;
    t0 = 1.5 * sqrt(6.) / (PI *s->fmax);
    MSG("s->fmax=%s\n",s->fmax);

    for (it = 0; it < s->time_steps; it++) {
        t1 = it * s->dt;
        source[it] = exp(-PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0)) *
                     (1.0 - 2. * PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0));
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

void source_ricker_wavelet_1st(sismap_t *s, float *source) {
    // 1st order ricker wavelet
    char tmp[256];
#ifdef __DUMP_SOURCE
    sprintf(tmp, "%s/first%s.txt", OUTDIR, SOURCE_BASE);
    FILE *fd = fopen(tmp,"w");
#endif // __DUMP_SOURCE
    float PI = 4.0f * atan(1.0f);
    int  it;
    float t1,t0;
    float deltaT,deltaT2,deltaT3;
    float a  = PI*s->fmax;
    float a2 = a  * a;
    float a4 = a2 * a2;
//    t0 = 1.0 / (s->fmax);
    t0 = 1.5 * sqrt(6.) / (PI *s->fmax);

    for (it = 0; it < s->time_steps; it++) {
        t1 = it * s->dt;
        deltaT=(t1-t0);
        deltaT2=deltaT  * deltaT;
        deltaT3=deltaT2 * deltaT;
        source[it] = exp(-a2*deltaT2) *(1.0 - 2. * a2 * deltaT2);
//        source[it] = exp(-a2*deltaT2) *(1.0 - 2. * a2 * deltaT2)/s->dt;
        //        source[it] = exp(-a2 * deltaT2) *(- 6.*a2*deltaT +4.*a4*deltaT3 );
//        source[it] = exp(-PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0)) *
//                    (1.0 - 2. * PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0)); ///s->dt;
#ifdef __DUMP_SOURCE
        fprintf(fd, "%d %f\n", it, source[it]);
#endif // __DUMP_SOURCE
    }

#ifdef __DUMP_SOURCE
    fclose(fd);
    sprintf(tmp, "%s/first%s.raw", OUTDIR, SOURCE_BASE);
    fd = fopen(tmp, "wb");
    ERR_IF(fd == NULL, "failed to open the source file for dumping");
    CHK(fwrite(source, sizeof(float),
               s->time_steps, fd) != s->time_steps,
               "failed to write in source file");
    fclose(fd);
#endif // __DUMP_SOURCE
}

void source_ricker_wavelet_2nd(sismap_t *s, float *source) {
    // 1st order ricker wavelet
    char tmp[256];
#ifdef __DUMP_SOURCE
    sprintf(tmp, "%s/second%s.txt", OUTDIR, SOURCE_BASE);
    FILE *fd = fopen(tmp,"w");
#endif // __DUMP_SOURCE
    float PI = 4.0f * atan(1.0f);
    int  it;
    float t1,t0,deltaT,deltaT2,deltaT3;
    t0 = 1.5 * sqrt(6.) / (PI *s->fmax);
    float a  = PI*s->fmax;
    float a2 = a  * a;
    float a4 = a2 * a2;
    for (it = 0; it < s->time_steps; it++) {
        t1 = it * s->dt;
        deltaT=(t1-t0);
        deltaT2=deltaT  * deltaT;
        deltaT3=deltaT2 * deltaT;
        source[it] = (- 6.*a2*deltaT +4.*a4*deltaT3 )*exp(-a2 * deltaT2);
//        source[it] = exp(-PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0)) *
//                     (1.0 - 2. * PI * PI * s->fmax * s->fmax * (t1 - t0) * (t1 - t0)); ///s->dt;
#ifdef __DUMP_SOURCE
        fprintf(fd, "%d %f\n", it, source[it]);
#endif // __DUMP_SOURCE
    }

#ifdef __DUMP_SOURCE
    fclose(fd);
    sprintf(tmp, "%s/second%s.raw", OUTDIR, SOURCE_BASE);
    fd = fopen(tmp, "wb");
    ERR_IF(fd == NULL, "failed to open the source file for dumping");
    CHK(fwrite(source, sizeof(float),
               s->time_steps, fd) != s->time_steps,
        "failed to write in source file");
    fclose(fd);
#endif // __DUMP_SOURCE
}
