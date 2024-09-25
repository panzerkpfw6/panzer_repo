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
/// @file src/wave.c
/// @brief This file contains the implementation of the new CPU wave descriptor.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <simwave/config.h>
#include <simwave/macros.h>
#include <simwave/wave.h>
#include <simwave/velocity.h>
#include <simwave/source.h>
#include <simwave/shot.h>
#include <simwave/pml.h>
#include <omp.h>

#define U0(z,y,x)   (u0[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define DAMPX(z, y, x)   (s->dampx[(x)])
#define DAMPY(z, y, x)   (s->dampy[(y)])
#define DAMPZ(z, y, x)   (s->dampz[(z)])
#define VX(z, y, x)   (vx[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define VY(z, y, x)   (vy[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define VZ(z, y, x)   (vz[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define U1(z,y,x)   (u1[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define FWD(z,y,x)  (fwd[(x+s->sx) + (2*s->sx + s->dimx) * \
                    ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])
#define ROC2(z,y,x) (roc2[x + (s->dimx * (s->dimy*(z) + y))])
#define  PHI(z,y,x) ( phi[x + (s->dimx * (s->dimy*(z) + y))])
#define  ETA(z,y,x) ( eta[(2+s->dimx)*((2+s->dimy)*(z+1) + (y+1)) + (x+1)])

#define IMG_S(z,y,x)  (img_shot[x + (s->img_dimx * (s->img_dimy*(z) + y))])
#define ILM_S(z,y,x)  (ilm_shot[x + (s->img_dimx * (s->img_dimy*(z) + y))])
#define   IMG(z,y,x)  (img[x + (s->img_dimx * (s->img_dimy*(z) + y))])

#define NDAMP 20


#define SETUP_STENCIL_COEFFS(t) \
t[0] = -205. /  72.;            \
t[1] =    8. /   5.;            \
t[2] =   -1. /   5.;            \
t[3] =    8. / 315.;            \
t[4] =   -1. / 560.;

#define SETUP_STENCIL_COEFFS_1st(t) \
t[0] =    1225. /   1024.;            \
t[1] =   -245. /   3072.;            \
t[2] =    49. / 5120.;            \
t[3] =   -5. / 7168.;

#define SCALE_STENCIL_COEFFS(delta, t) \
t[0] /= pow(delta,2);                  \
t[1] /= pow(delta,2);                  \
t[2] /= pow(delta,2);                  \
t[3] /= pow(delta,2);                  \
t[4] /= pow(delta,2);

#define SCALE_STENCIL_COEFFS_1st(delta, t) \
t[0] /= pow(delta,1);                  \
t[1] /= pow(delta,1);                  \
t[2] /= pow(delta,1);                  \
t[3] /= pow(delta,1);

void array_openmp_init(float* u, sismap_t *s) {
    const int nnx = s->dimx + 2* s->sx;
    const int nny = s->dimy + 2* s->sy;
    const int nnxy = nnx*nny;
    float * restrict ux;
    #pragma omp parallel for collapse(2) private(ux)
    for (int zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            int zmax = zmin+BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            int ymax = ymin+BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for(int z = zmin; z < zmax; z++) {
                for(int y = ymin; y < ymax; y++) {
                    ux = &(u[1ULL* (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
//          #pragma simd
                    for (int x = 0; x < s->dimx; x++) {
                        ux[x] = 0.;
                    }
                }
            }
        }
    }
}

void array_openmp_inner_init(float *u, sismap_t *s) {
    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnxy = nnx * nny;
    float * restrict ux;
#pragma omp parallel for collapse(2) private(ux)
    for (int zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            int zmax = zmin + BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            int ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for (int z = zmin; z < zmax; z++) {
                for (int y = ymin; y < ymax; y++) {
                    ux = &(u[1ULL * z * s->dimx * s->dimy + y * s->dimx]);
//          #pragma simd
                    for (int x = 0; x < s->dimx; x++) {
                        ux[x] = 0.;
                    }
                }
            }
        }
    }
}

void wave_init_numerics(sismap_t *s) {
    unsigned int i;
    s->dim2 = (s->vel_dimy == 1);
    s->sx = 4;
    s->sy = 4;
    s->sz = 4;
    MSG("before velocity_query_model");
    velocity_query_model(s);
    s->lambda = ((unsigned int) (s->vmin / s->fmax));
//    s->dx = MIN(s->lambda / 5, s->dx);
//    s->dy = s->dim2 ? 1 : MIN(s->lambda / 5, s->dy);
//    s->dz = MIN(s->lambda / 5, s->dz);

    s->hdx2 = 1. / (4. * pow(s->dx, 2));
    s->hdy2 = s->dim2 ? 0.0 : 1. / (4. * pow(s->dy, 2));
    s->hdz2 = 1. / (4. * pow(s->dz, 2));
    MSG("wave_init_numerics. point1");
//    CREATE_BUFFER(s->coefx, 2 * s->sx + 1);
//    CREATE_BUFFER(s->coefy, 2 * s->sy + 1);
//    CREATE_BUFFER(s->coefz, 2 * s->sz + 1);
    CREATE_BUFFER(s->coefx, 4);
    CREATE_BUFFER(s->coefy, 4);
    CREATE_BUFFER(s->coefz, 4);
    if (s->order==1) {
        MSG("initialize 1st order scheme coefs");
        SETUP_STENCIL_COEFFS_1st(s->coefx);
        SETUP_STENCIL_COEFFS_1st(s->coefy);
        SETUP_STENCIL_COEFFS_1st(s->coefz);
//        MSG("!!!coefx:0=%f,1=%f,2=%f,3=%f,4=%f,5=%f,6=%f,7=%f\n",s->coefx[0],s->coefx[1],s->coefx[2],s->coefx[3],s->coefx[4],s->coefx[5],s->coefx[6]);
        MSG("!!!coefx:0=%f,1=%f,2=%f,3=%f\n",s->coefx[0],s->coefx[1],s->coefx[2],s->coefx[3]);
//        SCALE_STENCIL_COEFFS_1st(s->dx, s->coefx);
//        SCALE_STENCIL_COEFFS_1st(s->dy, s->coefy);
//        SCALE_STENCIL_COEFFS_1st(s->dz, s->coefz);
    } else{
        MSG("initialize 2nd order scheme coefs");
        SETUP_STENCIL_COEFFS(s->coefx);
        if (!s->dim2) {
            SETUP_STENCIL_COEFFS(s->coefy);
        } else {
            for (i = 0; i < (2 * s->sy + 1); i++) s->coefy[i] = 0.0;
        }
        SETUP_STENCIL_COEFFS(s->coefz);
        SCALE_STENCIL_COEFFS(s->dx, s->coefx);
        if (!s->dim2) SCALE_STENCIL_COEFFS(s->dy, s->coefy);
        SCALE_STENCIL_COEFFS(s->dz, s->coefz);
    }

    s->pmlx = PML_TAPER * s->lambda / s->dx;
    s->pmly = s->dim2 ? 0 : PML_TAPER * s->lambda / s->dy;
    s->pmlz = PML_TAPER * s->lambda / s->dz;
    s->courant_number =
            fabs(s->coefx[0]) + fabs(s->coefy[0]) + fabs(s->coefz[0]);
    for (i = 1; i < s->sx + 1; i++)
        s->courant_number += 2. * fabs(s->coefx[i]);
    for (i = 1; i < s->sy + 1; i++)
        s->courant_number += 2. * fabs(s->coefy[i]);
    for (i = 1; i < s->sz + 1; i++)
        s->courant_number += 2. * fabs(s->coefz[i]);
    s->courant_number = 2. / sqrt(s->courant_number);
    s->dt = s->cfl * s->courant_number / s->vmax;
//    s->dt =0.001;
    // check that nbsnap is < Nyquist limit:
    s->nyquist_sampling = floor(1.0 / (2.0 * s->fmax) / s->dt) + 1;
    if (s->nb_snap == -1) {
        MSG("... setting nb_snap = Nyquist sampling: %d\n", s->nyquist_sampling);
        s->nb_snap = s->nyquist_sampling;
    } else {
        if (s->nb_snap > s->nyquist_sampling) {
            MSG("... STENCIL information:");
            MSG("... Nyquist sampling=: %d\n", s->nyquist_sampling);
            ERR("... nb_snap cannot > Nyquist sampling");
        }
    }
}

void wave_init_dimensions(sismap_t *s) {
    s->dtrpx = s->dcdp / s->dx;
    s->dtrpy = s->dim2 ? 1 : s->dline / s->dy;
    s->dtrpz = s->ddepth / s->dz;
    s->img_dimx = (s->vel_dimx - 1) * s->dtrpx + 1;// ceil(((float)s->vel_dimx*s->dcdp)/(float)s->dx);
    s->img_dimy = (s->vel_dimy - 1) * s->dtrpy + 1;// ceil(((float)s->vel_dimy*s->dline)/(float)s->dy);
    s->img_dimz = (s->vel_dimz - 1) * s->dtrpz + 1;// ceil(((float)s->vel_dimz*s->ddepth)/(float)s->dz);
    s->size_img = 1ULL * (s->img_dimx) * (s->img_dimy) * (s->img_dimz);
    /// augment the dimensions to take into accounts the PML and free surface.
    s->dimx = s->img_dimx + 2 * s->pmlx;
    s->dimy = s->dim2 ? 1 : s->img_dimy + 2 * s->pmly;
    s->dimz = s->img_dimz + 1 * s->pmlz; /// free surface on the top.
    s->size = 1ULL * (2 * s->sx + s->dimx) * (2 * s->sy + s->dimy) * (2 * s->sz + s->dimz);
    s->size_eff = 1ULL * (s->dimx) * (s->dimy) * (s->dimz);

}

void wave_init_damp(sismap_t *s) {
    CREATE_BUFFER(s->dampx, 2 * s->sx + s->dimx);
    CREATE_BUFFER(s->dampy, 2 * s->sy + s->dimy);
    CREATE_BUFFER(s->dampz, 2 * s->sz + s->dimz);

    float alpha = 0.2;
    float tabdamp[NDAMP];

    for (int i = 1; i <= NDAMP; i++) {
        tabdamp[NDAMP - i] = exp(-alpha * (1.0 * i / NDAMP) * (1.0 * i / NDAMP));
    }

    for (int i = s->sx; i < s->sx + s->dimx; i++) {
        s->dampx[i] = 1.0;
    }
    for (int i = s->sy; i < s->sy + s->dimy; i++) {
        s->dampy[i] = 1.0;
    }
    for (int i = s->sz; i < s->sz + s->dimz; i++) {
        s->dampz[i] = 1.0;
    }

    for (int i = 0; i < NDAMP; i++) {
        s->dampx[s->sx + i] = tabdamp[i];
        s->dampy[s->sy + i] = tabdamp[i];
//    s->dampz[s->sz+i] = tabdamp[i];

        s->dampx[s->sx + s->dimx - 1 - i] = tabdamp[i];
        s->dampy[s->sy + s->dimy - 1 - i] = tabdamp[i];
        s->dampz[s->sz + s->dimz - 1 - i] = tabdamp[i];
    }
}

void wave_init_acquisition(sismap_t *s) {
    /// initialize shots geometry:
    int a = 0;
    s->rcv_len = 0;

    for (int y = (s->dim2 ? 0 : s->pmly + (s->drcv * s->dtrpy) - 1);
         y < s->dimy - s->pmly; y += (s->drcv * s->dtrpy)) {
        a++;
        for (int x = s->pmlx + (s->drcv * s->dtrpx) - 1;
             x < s->dimx - s->pmlx; x += (s->drcv * s->dtrpx)) {
            s->rcv_len = s->rcv_len + 1;
        }
    }
    MSG(" s->rcv_len, %d\n", s->rcv_len);
    MSG("number of receivers in y direction, %d\n",a);
    s->rcv = (unsigned int *) malloc(s->rcv_len * sizeof(unsigned int));
    int ir = -1;

    /// record acquisition file.
    char tmp[512];
    sprintf(tmp,"%s/acquisition.txt", OUTDIR);
    FILE *fd = fopen(tmp, "w");
    printf("%s\n", tmp);
    CHK(fd == NULL, "failed to open rcv file");

    for (int y = (s->dim2 ? 0 : s->pmly + (s->drcv * s->dtrpy) - 1);
         y < s->dimy - s->pmly; y += (s->drcv * s->dtrpy)) {
        for (int x = s->pmlx + (s->drcv * s->dtrpx) - 1;
             x < s->dimx - s->pmlx; x += (s->drcv * s->dtrpx)) {
            s->rcv[ir++] = (s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx);
//            fprintf(fd,"rec %3d in [%3d, %3d, %3d]\n", ir, x, y, s->rcv_depth);
            fprintf(fd,"rec %3d in [%3d, %3d, %3d] at %3d\n",ir,x,y,s->rcv_depth,(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx));
        }
        fprintf(fd, "\n");
    }
    fclose(fd);

    /// record shots file.
    char tmp2[512];
    sprintf(tmp2,"%s/shots.txt", OUTDIR);
    FILE *fd2 = fopen(tmp2, "w");
    printf("%s\n", tmp2);
    CHK(fd2 == NULL, "failed to open shots file");

    /// setup shots geometries.
    s->snap_idx = 0;
    s->nb_shots = 0;
    for (unsigned int ir = 0; ir < s->rcv_len; ir = ir + s->dshot) {
        s->nb_shots++;
    }
    if (s->first == -1) s->first = 0;
    if (s->last == -1) s->last = s->nb_shots - 1;
    s->shots = (shot_t **) malloc(sizeof(shot_t * ) * s->nb_shots);
    for (unsigned int idx = 0; idx < s->nb_shots; idx++)
        s->shots[idx] = (shot_t *) malloc(sizeof(shot_t));
    unsigned int idx = 0;
    for (unsigned int ir = 0; ir < s->rcv_len; ir = ir + s->dshot) {
        s->shots[idx]->srcidx = s->rcv[ir] + s->drcv / 2;
//      MSG("s->shots[idx]->srcidx=ir=, %d, %d, %d\n",s->shots[idx]->srcidx,ir);
        s->shots[idx]->id=idx;
        fprintf(fd2,"shot ir=%3d,srcidx=%3d,id=%3d,s->drcv/2=%d\n",ir,s->shots[idx]->srcidx,s->shots[idx]->id,s->drcv/2);
        idx++;
    }
    fclose(fd2);
}

void wave_release(sismap_t *s) {
    free(s->coefx);
    free(s->coefy);
    free(s->coefz);
    free(s->rcv);
    for (int idx = 0; idx < s->nb_shots; idx++) free(s->shots[idx]);
    free(s->shots);
    free(s->dampx);
    free(s->dampy);
    free(s->dampz);
}

void wave_print(sismap_t *s) {
    int i;
    MSG(" ");
    MSG("... simwave information:");
    MSG("... velocity size       = %u x %u x %u",
        s->vel_dimx, s->vel_dimy, s->vel_dimz);
    MSG("... velocity min/max    = %f - %f", s->vmin, s->vmax);
    MSG("... compute domain size = %u x %u x %u (%f MB)",
        s->dimx, s->dimy, s->dimz, s->size / 1024. / 1024.);
    MSG("... imaging domain size = %u x %u x %u (%f MB)",
        s->img_dimx, s->img_dimy, s->img_dimz,
        s->size_img / 1024. / 1024.);
    MSG("... cdp,line,depth      = %u x %u x %u", s->dcdp, s->dline, s->ddepth);
    MSG("... dS                  = %u x %u x %u", s->dx, s->dy, s->dz);
    MSG("... dtrp                = %u x %u x %u", s->dtrpx, s->dtrpy, s->dtrpz);
    MSG("... dt                  = %g", s->dt);
    MSG("... stencil size        = %u x %u x %u", 2 * s->sx + 1, 2 * s->sy + 1, 2 * s->sz + 1);
    MSG("... pml damping         = %u x %u x %u", s->pmlx, s->pmly, s->pmlz);
    MSG("... source depth        = %u ", s->src_depth);
    MSG_NLR("... coefx               = ");
    for (i = s->sx; i > 0; i--) MSG_INL("%g ", s->coefx[i]);
    MSG_INL("%g ", s->coefx[0]);
    for (i = 1; i < s->sx + 1; i++) MSG_INL("%g ", s->coefx[i]);
    MSG_JLR();
    MSG_NLR("... coefy               = ");
    for (i = s->sy; i > 0; i--) MSG_INL("%g ", s->coefy[i]);
    MSG_INL("%g ", s->coefy[0]);
    for (i = 1; i < s->sy + 1; i++) MSG_INL("%g ", s->coefy[i]);
    MSG_JLR();
    MSG_NLR("... coefz               = ");
    for (i = s->sz; i > 0; i--) MSG_INL("%g ", s->coefz[i]);
    MSG_INL("%g ", s->coefz[0]);
    for (i = 1; i < s->sz + 1; i++) MSG_INL("%g ", s->coefz[i]);
    MSG_JLR();
    MSG("... hdx2, hdy2, hdz2    = %f, %f, %f", s->hdx2, s->hdy2, s->hdz2);
    MSG("... time steps          = %u", s->time_steps);
    MSG("... acquisition geometry:");
    MSG("... number of receivers              = %u", s->rcv_len);
    MSG("... depth of the receivers           = %u", s->rcv_depth);
    MSG("... space between 2 receivers        = %u", s->drcv);
    MSG("... number of shots                  = %u", s->nb_shots);
    MSG(" ");
}

#define WAVE_COMPUTE_LAPLACIAN()                            \
  laplacian = coef0 * U0(z, y, x)                           \
    + s->coefx[1]*( U0(z,   y,   x+1) + U0(z,   y,   x-1))  \
    + s->coefy[1]*( U0(z,   y+1, x  ) + U0(z,   y-1, x  ))  \
    + s->coefz[1]*( U0(z+1, y,   x  ) + U0(z-1, y,   x  ))  \
    + s->coefx[2]*( U0(z,   y,   x+2) + U0(z,   y,   x-2))  \
    + s->coefy[2]*( U0(z,   y+2, x  ) + U0(z,   y-2, x  ))  \
    + s->coefz[2]*( U0(z+2, y,   x  ) + U0(z-2, y,   x  ))  \
    + s->coefx[3]*( U0(z,   y,   x+3) + U0(z,   y,   x-3))  \
    + s->coefy[3]*( U0(z,   y+3, x  ) + U0(z,   y-3, x  ))  \
    + s->coefz[3]*( U0(z+3, y,   x  ) + U0(z-3, y,   x  ))  \
    + s->coefx[4]*( U0(z,   y,   x+4) + U0(z,   y,   x-4))  \
    + s->coefy[4]*( U0(z,   y+4, x  ) + U0(z,   y-4, x  ))  \
    + s->coefz[4]*( U0(z+4, y,   x  ) + U0(z-4, y,   x  ))

#define WAVE_UPDATE_INNER_FIELDS() \
  U1(z,y,x) = 2.0f * U0(z,y,x) - U1(z,y,x) + ROC2(z,y,x) * laplacian

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD()                        \
ux[x] = 2.0f * vx[x] - ux[x]                                                   \
      + rx[x] * (coef0 * vx[x] + s->coefx[1] * (vx[x+1     ] + vx[x-1     ])   \
                               + s->coefy[1] * (vx[x+nnx   ] + vx[x-nnx   ])   \
                               + s->coefz[1] * (vx[x+nnxy  ] + vx[x-nnxy  ])   \
                               + s->coefx[2] * (vx[x+2     ] + vx[x-2     ])   \
                               + s->coefy[2] * (vx[x+2*nnx ] + vx[x-2*nnx ])   \
                               + s->coefz[2] * (vx[x+2*nnxy] + vx[x-2*nnxy])   \
                               + s->coefx[3] * (vx[x+3     ] + vx[x-3     ])   \
                               + s->coefy[3] * (vx[x+3*nnx ] + vx[x-3*nnx ])   \
                               + s->coefz[3] * (vx[x+3*nnxy] + vx[x-3*nnxy])   \
                               + s->coefx[4] * (vx[x+4     ] + vx[x-4     ])   \
                               + s->coefy[4] * (vx[x+4*nnx ] + vx[x-4*nnx ])   \
                               + s->coefz[4] * (vx[x+4*nnxy] + vx[x-4*nnxy]))  \

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v()                        \
vx0[x] = vx0[x]                                                                \
      + s->dt/s->dx*(    s->coefx[0] * (pr0[x+1] - pr0[x])   \
                   + s->coefx[1] * (pr0[x+2] - pr0[x-1])   \
                   + s->coefx[2] * (pr0[x+3] - pr0[x-2])   \
                   + s->coefx[3] * (pr0[x+4] - pr0[x-3])) ;  \
vy0[x] = vy0[x]                                                                \
        + s->dt/s->dy * (      s->coefy[0] * (pr0[x+1*nnx] - pr0[x])   \
                       + s->coefy[1] * (pr0[x+2*nnx] - pr0[x-1*nnx])   \
                       + s->coefy[2] * (pr0[x+3*nnx] - pr0[x-2*nnx])   \
                       + s->coefy[3] * (pr0[x+4*nnx] - pr0[x-3*nnx])) ;  \
vz0[x] = vz0[x]                                                         \
      + s->dt/s->dz * (                s->coefz[0] * (pr0[x+1*nnxy] - pr0[x])   \
                               + s->coefz[1] * (pr0[x+2*nnxy] - pr0[x-1*nnxy])   \
                               + s->coefz[2] * (pr0[x+3*nnxy] - pr0[x-2*nnxy])   \
                               + s->coefz[3] * (pr0[x+4*nnxy] - pr0[x-3*nnxy]))

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p()                        \
pr0[x] = pr0[x]                                                   \
      + rx[x] * (                                                             \
                    s->coefx[0]/s->dx * (vx0[x]           - vx0[x-1     ])   \
                   + s->coefy[0]/s->dy * (vy0[x]          - vy0[x-nnx   ])   \
                   + s->coefz[0]/s->dz * (vz0[x]          -vz0[x-nnxy  ])   \
                   + s->coefx[1]/s->dx * (vx0[x+1     ]   - vx0[x-2     ])   \
                   + s->coefy[1]/s->dy * (vy0[x+1*nnx ]   - vy0[x-2*nnx ])   \
                   + s->coefz[1]/s->dz * (vz0[x+1*nnxy]    -vz0[x-2*nnxy])   \
                   + s->coefx[2]/s->dx * (vx0[x+2     ]    -vx0[x-3     ])   \
                   + s->coefy[2]/s->dy * (vy0[x+2*nnx ]    -vy0[x-3*nnx ])   \
                   + s->coefz[2]/s->dz * (vz0[x +2*nnxy]   -vz0[x-3*nnxy])   \
                   + s->coefx[3]/s->dx * (vx0[x+3     ]   - vx0[x-4     ])   \
                   + s->coefy[3]/s->dy * (vy0[x+3*nnx ]   - vy0[x-4*nnx ])   \
                   + s->coefz[3]/s->dz * (vz0[x+ 3*nnxy]   -vz0[x-4*nnxy])                 \
                   )

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v_index()                            \
    VX(z, y, x) = VX(z, y, x)+ s->dt/s->dx * \
        (s->coefx[0]*( U0(z,   y,   x+1) - U0(z,y,x))  \
        + s->coefx[1]*( U0(z,   y,   x+2) - U0(z,   y,   x-1))  \
        + s->coefx[2]*( U0(z,   y,   x+3) - U0(z,   y,   x-2))  \
        + s->coefx[3]*( U0(z,   y,   x+4) - U0(z,   y,   x-3))  );\
    VY(z, y, x) = VY(z, y, x)+ s->dt/s->dy *\
        (s->coefy[0]*( U0(z,   y+1, x  ) - U0(z,   y, x  ))  \
        + s->coefy[1]*( U0(z,   y+2, x  ) - U0(z,   y-1, x  ))  \
        + s->coefy[2]*( U0(z,   y+3, x  ) - U0(z,   y-2, x  ))  \
        + s->coefy[3]*( U0(z,   y+4, x  ) - U0(z,   y-3, x  )));\
    VZ(z, y, x) = VZ(z, y, x)+  s->dt/s->dz *\
        (s->coefz[0]*( U0(z+1, y,   x  ) - U0(z, y,   x  ))  \
        + s->coefz[1]*( U0(z+2, y,   x  ) - U0(z-1, y,   x  ))  \
        + s->coefz[2]*( U0(z+3, y,   x  ) - U0(z-2, y,   x  ))  \
        + s->coefz[3]*( U0(z+4, y,   x  ) - U0(z-3, y,   x  )))

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p_index() \
    U0(z, y, x) = U0(z, y, x)+ROC2(z, y, x)*  \
        ( s->coefx[0]/s->dx*( VX(z,   y,   x) - VX(z,y,x-1))                \
        + s->coefy[0]/s->dy*( VY(z,   y,   x) - VY(z,y-1,x))                \
        + s->coefz[0]/s->dz*( VZ(z,   y,   x) - VZ(z-1,y,x))                \
        + s->coefx[1]/s->dx*( VX(z,   y,   x+1) - VX(z,   y,   x-2))  \
        + s->coefy[1]/s->dy*( VY(z,   y+1,   x) - VY(z,   y-2,   x))  \
        + s->coefz[1]/s->dz*( VZ(z+1,   y,   x) - VZ(z-2,   y,   x))      \
        + s->coefx[2]/s->dx*( VX(z,   y,   x+2) - VX(z,y,x-3))                \
        + s->coefy[2]/s->dy*( VY(z,   y+2,   x) - VY(z,y-3,x))                \
        + s->coefz[2]/s->dz*( VZ(z+2,   y,   x) - VZ(z-3,y,x))  \
        + s->coefx[3]/s->dx*( VX(z,   y,   x+3) - VX(z,   y,   x-4))  \
        + s->coefy[3]/s->dy*( VY(z,   y+3,   x) - VY(z,   y-4,   x))      \
        + s->coefz[3]/s->dz*( VZ(z+3,   y,   x) - VZ(z-4,y,x)));          \
    U0(z, y, x) = U0(z, y, x)*DAMPX(z, y, x)*DAMPY(z, y, x)*DAMPZ(z, y, x)

#define WAVE_UPDATE_PML_FIELDS()                                            \
  U1(z,y,x) =                                                               \
    ((2.-ETA(z,y,x)*ETA(z,y,x) + 2.*ETA(z,y,x))*U0(z,y,x)                   \
    - U1(z,y,x) + ROC2(z,y,x)*(laplacian + PHI(z,y,x)))/(1.+2.*ETA(z,y,x)); \
  PHI(z,y,x)= (PHI(z,y,x)-                                                  \
     (( ETA(z,   y,   x+1) - ETA(z,   y,   x-1))                            \
      *( U0(z,   y,   x+1) -  U0(z,   y,   x-1))*s->hdx2                    \
      +(ETA(z,   y+1, x  ) - ETA(z,   y-1, x  ))                            \
      *( U0(z,   y+1, x  ) -  U0(z,   y-1, x  ))*s->hdy2                    \
      +(ETA(z+1, y,   x  ) - ETA(z-1, y,   x  ))                            \
      *( U0(z+1, y,   x  ) -  U0(z-1, y,   x  ))*s->hdz2))/(1.+ETA(z,y,x))

void wave_update_fields(sismap_t *s,
                        float *u0, float *u1,
                        float *roc2, float *phi, float *eta) {
    unsigned int z, y, x;
    float laplacian;
    float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
#pragma omp parallel for private(laplacian, z, y, x)
    for (z = 0; z < s->dimz; z++) {
        for (y = 0; y < s->dimy; y++) {
            for (x = 0; x < s->dimx; x++) {
                WAVE_COMPUTE_LAPLACIAN();
                /// free surface on the top.
                if ((z >= 0) && (z < s->dimz - s->pmlz) &&
                    (y >= s->pmly) && (y < s->dimy - s->pmly) &&
                    (x >= s->pmlx) && (x < s->dimx - s->pmlx)) {
                    WAVE_UPDATE_INNER_FIELDS();
                } else {
                    WAVE_UPDATE_PML_FIELDS();
                }
            }
        }
    }
}

void wave_update_fields_1st(sismap_t *s,
                        float *u0, float *vx,float *vy,float *vz,
                        float *roc2, float *phi, float *eta) {
    unsigned int z, y, x;

    #pragma omp parallel for private(z, y, x)
        for (z = 0; z < s->dimz; z++) {     // v loop
            for (y = 0; y < s->dimy; y++) {
                for (x = 0; x < s->dimx; x++)   {
                    WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v_index();
                }
            }
        }
    #pragma omp parallel for private(z, y, x)
        for (z = 0; z < s->dimz; z++) {    // p loop
            for (y = 0; y < s->dimy; y++) {
                for (x = 0; x < s->dimx; x++) {
                    WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p_index();
                }
            }
        }
}
void wave_update_fields_block(sismap_t *s,
                              float *u0, float *u1,
                              float *roc2, float *phi, float *eta) {
    unsigned int z, y, x;
    float laplacian;
    float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
    unsigned int zmin, zmax, ymin, ymax;

#pragma omp parallel for collapse(2) private(laplacian, z, y, x, zmin, zmax, ymin, ymax)
    for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            zmax = zmin + BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for (z = zmin; z < zmax; z++) {
                for (y = ymin; y < ymax; y++) {
                    for (x = 0; x < s->dimx; x++) {
                        WAVE_COMPUTE_LAPLACIAN();
                        /// free surface on the top.
                        if ((z >= 0) && (z < s->dimz - s->pmlz) &&
                            (y >= s->pmly) && (y < s->dimy - s->pmly) &&
                            (x >= s->pmlx) && (x < s->dimx - s->pmlx)) {
                            WAVE_UPDATE_INNER_FIELDS();
                        } else {
                            WAVE_UPDATE_PML_FIELDS();
                        }
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_bis(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict u1,
                                  float *restrict roc2,
                                  float *restrict phi,
                                  float *restrict eta) {
    unsigned int z, y, x;
    float laplacian;
    float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
    unsigned int zmin, zmax, ymin, ymax;

    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnxy = nnx * nny;

    float *restrict ux;
    float *restrict vx;
    float *restrict rx;
#pragma omp parallel for collapse(2) private(laplacian, zmin, zmax, ymin, ymax, ux, vx, rx)
    for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            zmax = zmin + BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for (int z = zmin; z < zmax; z++) {
                for (int y = ymin; y < ymax; y++) {
                    ux = &(u1[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vx = &(u0[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    rx = &(roc2[1ULL * z * s->dimx * s->dimy + y * s->dimx]);
//           #pragma simd
                    for (int x = 0; x < s->dimx; x++) {
                        //            WAVE_COMPUTE_LAPLACIAN();
                        //            WAVE_UPDATE_INNER_FIELDS();
                        //            U1(z,y,x) = s->dampx[x+s->sx] * U1(z,y,x) + (1 - s->dampx[x+s->sx]) * U0(z,y,x);
                        //            U1(z,y,x) = s->dampy[y+s->sy] * U1(z,y,x) + (1 - s->dampy[y+s->sy]) * U0(z,y,x);
                        //            U1(z,y,x) = s->dampz[z+s->sz] * U1(z,y,x) + (1 - s->dampz[z+s->sz]) * U0(z,y,x);
                        WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD();
                        ux[x] = s->dampx[x + s->sx] * ux[x] + (1 - s->dampx[x + s->sx]) * vx[x];
                        ux[x] = s->dampy[y + s->sy] * ux[x] + (1 - s->dampy[y + s->sy]) * vx[x];
                        ux[x] = s->dampz[z + s->sz] * ux[x] + (1 - s->dampz[z + s->sz]) * vx[x];
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_1st(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict vx,
                                  float *restrict vy,
                                  float *restrict vz,
                                  float *restrict roc2,
                                  float *restrict phi,
                                  float *restrict eta) {
    unsigned int z, y, x;
    float laplacian;
    unsigned int zmin, zmax, ymin, ymax;

    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnxy = nnx * nny;
    float *restrict pr0;
    float *restrict vx0;
    float *restrict vy0;
    float *restrict vz0;
    float *restrict rx;

    // loop on the blocks .velocity
#pragma omp parallel for collapse(2) private(laplacian, zmin, zmax, ymin, ymax,pr0,vx0,vy0,vz0)
    for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            zmax = zmin + BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;

            for (int z = zmin; z < zmax; z++) {
                for (int y = ymin; y < ymax; y++) {
                    pr0 = &(u0[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vx0 = &(vx[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vy0 = &(vy[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vz0 = &(vz[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
//           #pragma simd
                    for (int x = 0; x < s->dimx; x++) {
                        WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v();
                    }
                }
            }
        }
    }

    // loop on the blocks.pressure
#pragma omp parallel for collapse(2) private(laplacian, zmin, zmax, ymin, ymax,pr0,vx0,vy0,vz0,rx)
    for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            zmax = zmin + BLOCKZ;
            if (zmax > s->dimz) zmax = s->dimz;
            ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;

            for (int z = zmin; z < zmax; z++) {
                for (int y = ymin; y < ymax; y++) {
                    pr0 = &(u0[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vx0 = &(vx[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vy0 = &(vy[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    vz0 = &(vz[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                    rx = &(roc2[1ULL * z * s->dimx * s->dimy + y * s->dimx]);
//           #pragma simd
                    for (int x = 0; x < s->dimx; x++) {
                        WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p();
                        pr0[x] = s->dampx[x + s->sx] * pr0[x];
                        pr0[x] = s->dampy[y + s->sy] * pr0[x];
                        pr0[x] = s->dampz[z + s->sz] * pr0[x];
                    }
                }
            }
        }
    }

}


void wave_update_source(sismap_t *s, shot_t *shot, float *u0, float sterm) {
//    MSG("shot->srcidx=%d\n",shot->srcidx);
    u0[(s->src_depth + s->sz) * (2 * s->sx + s->dimx)*(2 * s->sy + s->dimy)+shot->srcidx] += sterm;
}

void wave_extract_sismos(sismap_t *s, float *u1, unsigned int t, float *sismos) {
    for (unsigned int ir = 0; ir < s->rcv_len; ir++) {
        sismos[s->rcv_len * t + ir] = u1[(s->rcv_depth + s->sz) * (2 * s->sx + s->dimx) * (2 * s->sy + s->dimy)+s->rcv[ir]];
    }
}

void wave_inject_sismos(sismap_t *s, float *u0, unsigned int t, float *sismos) {
    for (unsigned int ir = 0; ir < s->rcv_len; ir++) {
        // u0[(s->rcv_depth + s->sz)*(2*s->sx + s->dimx)*(2*s->sy + s->dimy)
        //                   + s->rcv[ir]] += sismos[s->rcv_len*t + ir];
        u0[(s->rcv_depth + s->sz) * (2 * s->sx + s->dimx) * (2 * s->sy + s->dimy)
           + s->rcv[ir]] += sismos[s->rcv_len * t + ir];
    }
}

/// @brief save snapshot.
void wave_save_snapshot(sismap_t *s, shot_t *shot, float *u1, unsigned int t) {
    if (t % s->nb_snap == 0) {
        CHK(fwrite(u1, sizeof(float), s->size, shot->fd_snap) != s->size,
            "failed to write snapshot");
        if (s->verbose) MSG("... saving snapshot %d t=%u",s->snap_idx,t);
        s->snap_idx = s->snap_idx + 1;
    }
}

void wave_read_snapshot(sismap_t *s, shot_t *shot, float *fwd, unsigned int t) {
    if (t % s->nb_snap == 0) {
        if (s->verbose) MSG("... reading snapshot %d t=%u", s->snap_idx, t);
        CHK(fseek(shot->fd_snap, s->snap_idx * s->size * sizeof(float), SEEK_SET) != 0,
            "failed to fseek file");
        CHK(fread(fwd, sizeof(float), s->size, shot->fd_snap) != s->size,
            "failed to read snapshot");
        s->snap_idx = s->snap_idx - 1;
    }
}

void wave_image_condition(sismap_t *s, float *u1,
                          float *fwd,
                          float *img_shot, float *ilm_shot, unsigned int t) {
    if (t % s->nb_snap == 0) {
        unsigned int zfwd = 0;
        unsigned int zu1 = 0;
        unsigned int x, y, z;
        for (z = 0; z < s->img_dimz; z++) {
            for (y = 0; y < s->img_dimy; y++) {
                for (x = 0; x < s->img_dimx; x++) {
                    IMG_S(z, y, x) += U1(z, y + s->pmly, x + s->pmlx) * FWD(z, y + s->pmly, x + s->pmlx);
                    ILM_S(z, y, x) += FWD(z, y + s->pmly, x + s->pmlx) * FWD(z, y + s->pmly, x + s->pmlx);
                    if (U1(z, y + s->pmly, x + s->pmlx) == 0) zu1++;
                    if (FWD(z, y + s->pmly, x + s->pmlx) == 0) zfwd++;
                }
            }
        }
        wave_min_max("u1 (img)", u1, s->size);
        wave_min_max("fwd (img)", fwd, s->size);
        wave_min_max("img_shot", img_shot, s->size);
        wave_min_max("ilm_shot", ilm_shot, s->size);
        printf("nbz u1  %f\n", zu1 * 100.0 / s->size_eff);
        printf("nbz fwd %f\n", zfwd * 100.0 / s->size_eff);
    }
}

void wave_image_condition_block(sismap_t *s, float *u1,
                                float *fwd,
                                float *img_shot, float *ilm_shot,
                                unsigned int t) {
    if (t % s->nb_snap == 0) {
        const int nnx = s->dimx + 2 * s->sx;
        const int nny = s->dimy + 2 * s->sy;
        const int nnxy = nnx * nny;

        float *restrict
        ux;
        float *restrict
        wx;
        float *restrict
        imgx;
        float *restrict
        ilmx;

#pragma omp parallel for collapse(2) private(ux, imgx, ilmx, wx)
        for (int zmin = 0; zmin < s->img_dimz; zmin += BLOCKZ) {
            for (int ymin = 0; ymin < s->img_dimy; ymin += BLOCKY) {
                int zmax = zmin + BLOCKZ;
                if (zmax > s->img_dimz) zmax = s->img_dimz;
                int ymax = ymin + BLOCKY;
                if (ymax > s->img_dimy) ymax = s->img_dimy;

                for (int z = zmin; z < zmax; z++) {
                    for (int y = ymin; y < ymax; y++) {
                        ux = &(u1[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                        wx = &(fwd[1ULL * (z + s->sz) * nnxy + (y + s->sy) * nnx + s->sx]);
                        imgx = &(img_shot[1ULL * z * s->dimx * s->dimy + y * s->dimx]);
                        ilmx = &(ilm_shot[1ULL * z * s->dimx * s->dimy + y * s->dimx]);
//            #pragma simd
                        for (int x = 0; x < s->img_dimx; x++) {
//                            if (ux[x]>0.0001)
//                                {MSG("ux[x]=%f,wx[x]=%f,I=%f",ux[x],wx[x],ux[x]*wx[x]);}
                            imgx[x] += ux[x] * wx[x];
                            ilmx[x] += wx[x] * wx[x];
//              IMG_S(z, y, x)+=  U1(z,y+s->pmly,x+s->pmlx)*FWD(z,y+s->pmly,x+s->pmlx);
//              ILM_S(z, y, x)+= FWD(z,y+s->pmly,x+s->pmlx)*FWD(z,y+s->pmly,x+s->pmlx);
                        }
                    }
                }
            }
        }
    }
}

void wave_save_sismos(sismap_t *s, shot_t *shot, float *sismos) {
    size_t sz = s->rcv_len * s->time_steps * sizeof(float);
    printf("sismos size is %lu\n", sz);
    char tmp[512];
    sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SISMOS_BASE, shot->id);
    FILE *fd = fopen(tmp, "wb");
    CHK(fd == NULL, "failed to open sismos file");
    CHK(fwrite(sismos, sz, 1, fd) != 1,
        "failed to properly write sismos");
    fclose(fd);
#ifdef __DEBUG
    sprintf(tmp, "%s/sismos_%d.txt", OUTDIR, shot->id);
    fd = fopen(tmp, "w");
    for (unsigned int t = 0; t < s->time_steps; t++) {
      fprintf(fd, "%d", t);
      for (unsigned int r = 0; r < s->rcv_len; r++) {
        fprintf(fd, " %f", sismos[t*s->rcv_len+r]);
//        MSG("t=: %d,... s=: %d\n",t,sismos[t*s->rcv_len+r]);
      }
      fprintf(fd,"\n");
    }
    fclose(fd);
#endif // __DEBUG
}

void wave_read_sismos(sismap_t *s, shot_t *shot, float *sismos) {
    char tmp[512];
    FILE *fd;
    size_t sz = s->rcv_len * s->time_steps * sizeof(float);
    sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SISMOS_BASE, shot->id);

    MSG("reading sismos=%s",tmp);
    fd = fopen(tmp, "rb");
    CHK(fd == NULL,"sismos file not found, aborting (run modeling to generate it)");
    CHK(fread(sismos, sz, 1, fd) != 1,"failed to properly read sismos");
    fclose(fd);
}

/// @brief save the final image.
void wave_save_image(sismap_t *s, float *img, char *fname) {
    char tmp[512];
    FILE *fd;
    printf("saving %s\n", fname);
    fd = fopen(fname, "wb+");
    CHK(fd == NULL, "failed to open image file");
    CHK(fwrite(img, s->size_img * sizeof(float), 1, fd) != 1,
        "failed to properly write image");
    fclose(fd);
}

void wave_save_fwd_dbg(sismap_t *s, shot_t *shot, float *u1, unsigned int t) {
    MSG("... dimz=: %d\n", s->dimz);
    MSG("... dimx=: %d\n", s->dimx);
    if (t % s->nb_snap == 0) {
        for (unsigned int z = 0; z < s->dimz; z++)
            for (unsigned int x = 0; x < s->dimx; x++)
                fwrite(&U1(z, shot->srcidx / (s->dimx + 2 * s->sx), x),
                       sizeof(float), 1, shot->fd_fwd);
    }
}

void wave_save_bwd_dbg(sismap_t *s, shot_t *shot, float *u1, unsigned int t) {
    if (t % s->nb_snap == 0) {
        for (unsigned int z = 0; z < s->dimz; z++)
            for (unsigned int x = 0; x < s->dimx; x++)
                fwrite(&U1(z, shot->srcidx / (s->dimx + 2 * s->sx), x),
                       sizeof(float), 1, shot->fd_bwd);
    }
}

void wave_save_img(sismap_t *s, shot_t *shot,
                   float *img_shot, float *ilm_shot) {

    CHK(fwrite(ilm_shot, sizeof(float), s->size_img, shot->fd_ilm) != s->size_img,
        "failed to save ilm");

    CHK(fwrite(img_shot, sizeof(float), s->size_img, shot->fd_img) != s->size_img,
        "failed to save img");
}

int wave_check_fields(float *tab, size_t len) {
    int not_valid = 0;
    for (size_t i = 0; i < len; ++i) {
        if (isnan(tab[i])) not_valid++;
    }
    return not_valid;
}

void wave_min_max(char *str, float *tab, size_t len) {
    float max = 0.0;
    float min = 99999999999999999999.0;
    for (unsigned int i = 0; i < len; i++) {
        if (tab[i] < min) min = tab[i];
        if (tab[i] > max) max = tab[i];
    }
    printf("... %s min max %f %f\n", str, min, max);
}
