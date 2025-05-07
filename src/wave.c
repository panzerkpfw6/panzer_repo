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
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/wave.h>
#include <stencil/velocity.h>
#include <stencil/source.h>
#include <stencil/shot.h>
#include <stencil/pml.h>
#include <omp.h>
#include <immintrin.h> // For AVX2 intrinsics
#include <stdio.h>
#include <errno.h>


#ifdef _WIN32
#include <direct.h> // For Windows _mkdir
#else
#include <sys/stat.h> // For POSIX mkdir
#endif


#define ALIGNMENT 32 // For AVX2 (256-bit = 8 floats)
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
#define PHI(z,y,x) (phi[x + (s->dimx * (s->dimy*(z) + y))])
#define ETA(z,y,x) (eta[(2+s->dimx)*((2+s->dimy)*(z+1) + (y+1)) + (x+1)])

#define IMG_S(z,y,x)  (img_shot[x + (s->img_dimx * (s->img_dimy*(z) + y))])
#define ILM_S(z,y,x)  (ilm_shot[x + (s->img_dimx * (s->img_dimy*(z) + y))])
#define IMG(z,y,x)  (img[x + (s->img_dimx * (s->img_dimy*(z) + y))])

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

int create_folder(const char *path) {
#ifdef _WIN32
    int result = _mkdir(path); // Windows
#else
    int result = mkdir(path, 0755); // POSIX, permissions 0755
#endif
    if (result == 0) {
        printf("Folder '%s' created successfully.\n", path);
        return 0;
    } else {
        // Handle errors
        switch (errno) {
            case EEXIST:
                printf("Folder '%s' already exists.\n", path);
                return 0; // Not an error if folder already exists
            case EACCES:
                fprintf(stderr, "Error: Permission denied to create '%s'.\n", path);
                break;
            case ENOENT:
                fprintf(stderr, "Error: Path component in '%s' does not exist.\n", path);
                break;
            default:
                fprintf(stderr, "Error: Failed to create '%s' (errno: %d).\n", path, errno);
        }
        return -1;
    }
}

void array_openmp_init(float* u, sismap_t *s) {
    const int nnx = s->dimx + 2 * s->sx;  // Total x-size with halo
    const int nny = s->dimy + 2 * s->sy;  // Total y-size with halo
    const int nnz = s->dimz + 2 * s->sz;  // Total z-size with halo
    const int nnyz = nny * nnz;           // Size of yz-plane for XYZ order
    float *restrict ux;
    const int BLOCKX=s->blockx;
    const int BLOCKY=s->blocky;
    const int BLOCKZ=s->blockz;

#pragma omp parallel for collapse(2) private(ux)
    for (int xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
        for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            int xmax = xmin + BLOCKX;
            if (xmax > s->dimx) xmax = s->dimx;
            int ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for (int x = xmin; x < xmax; x++) {
                for (int y = ymin; y < ymax; y++) {
                    ux = &(u[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                    for (int z = 0; z < s->dimz; z++) {
                        ux[z] = 0.;
                    }
                }
            }
        }
    }
}

void array_openmp_inner_init(float *u, sismap_t *s) {
    const int nnyz = s->dimy * s->dimz;  // Physical yz-plane size (no halo)
    float *restrict ux;
    const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;

    #pragma omp parallel for collapse(2) private(ux)
    for (int xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
        for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            int xmax = xmin + BLOCKX;
            if (xmax > s->dimx) xmax = s->dimx;
            int ymax = ymin + BLOCKY;
            if (ymax > s->dimy) ymax = s->dimy;
            for (int x = xmin; x < xmax; x++) {
                for (int y = ymin; y < ymax; y++) {
                    ux = &(u[1ULL * x * nnyz + y * s->dimz]);
                    #pragma omp simd
                    for (int z = 0; z < s->dimz; z++) {
                        ux[z] = 0.;
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
        SCALE_STENCIL_COEFFS(s->dy, s->coefy);
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
//    s->dt = s->cfl * s->courant_number / s->vmax; // for 2nd order
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
//    s->dtrpx = s->dcdp / s->dx;
//    s->dtrpy = s->dim2 ? 1 : s->dline / s->dy;
//    s->dtrpz = s->ddepth / s->dz;
    /// set to 1 for simplicity
    s->dtrpx=1;s->dtrpy=1;s->dtrpz=1;

    s->img_dimx = (s->vel_dimx - 1) * s->dtrpx + 1;// ceil(((float)s->vel_dimx*s->dcdp)/(float)s->dx);
    s->img_dimy = (s->vel_dimy - 1) * s->dtrpy + 1;// ceil(((float)s->vel_dimy*s->dline)/(float)s->dy);
    s->img_dimz = (s->vel_dimz - 1) * s->dtrpz + 1;// ceil(((float)s->vel_dimz*s->ddepth)/(float)s->dz);
    s->size_img = 1ULL * (s->img_dimx) * (s->img_dimy) * (s->img_dimz);

    /// Commented out because I am nousing PML regions
    /// augment the dimensions to take into accounts the PML and free surface.
//    s->dimx = s->img_dimx + 2 * s->pmlx;
//    s->dimy = s->dim2 ? 1 : s->img_dimy + 2 * s->pmly;
//    s->dimz = s->img_dimz + 1 * s->pmlz; /// free surface on the top.

    s->dimx = s->img_dimx;
    s->dimy = s->img_dimy;
    s->dimz = s->img_dimz;

    s->size = 1ULL * (2 * s->sx + s->dimx) * (2 * s->sy + s->dimy) * (2 * s->sz + s->dimz);
    s->size_eff = 1ULL * (s->dimx) * (s->dimy) * (s->dimz);
}

void wave_init_damp(sismap_t *s) {
    CREATE_BUFFER(s->dampx, 2 * s->sx + s->dimx);
    CREATE_BUFFER(s->dampy, 2 * s->sy + s->dimy);
    CREATE_BUFFER(s->dampz, 2 * s->sz + s->dimz);

//    MSG("s->sx=%d\n",s->sx);
//    MSG("s->dimx=%d\n",s->dimx);

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
	const char *folder_path = "data"; // Folder to create
	create_folder(folder_path);
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
//    int ir = -1;
    int ir=0;

    /// record acquisition file.
    char tmp[512];
    sprintf(tmp,"%s/acquisition.txt", OUTDIR);
    FILE *fd = fopen(tmp, "w");
    printf("%s\n", tmp);
    CHK(fd == NULL, "failed to open rcv file");

	for (int x = s->pmlx + (s->drcv * s->dtrpx) - 1;
		 x < s->dimx - s->pmlx; x += (s->drcv * s->dtrpx)) {
		for (int y = (s->dim2 ? 0 : s->pmly + (s->drcv * s->dtrpy) - 1);
			 y < s->dimy - s->pmly; y += (s->drcv * s->dtrpy)) {
//            s->rcv[ir]=(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx);   //simwave:zyx
            s->rcv[ir]=(s->sx+x)*(s->dimy+2*s->sy)+(y+s->sy);   //stencil:xyz
//            fprintf(fd,"rec %3d in [%3d, %3d, %3d] at %3d\n",ir,x,y,s->rcv_depth,(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx));   //simwave:zyx
            fprintf(fd,"rec %3d in [%3d, %3d, %3d] at %3d\n",ir,x,y,s->rcv_depth,(s->sx+x)*(s->dimy+2*s->sy)+(y+s->sy));   //stencil:xyz
            ir++;
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
//        s->shots[idx]->srcidx = s->rcv[ir]+s->drcv/2;   // should I really use "s->drcv/2"?
        s->shots[idx]->srcidx=s->rcv[ir];
        s->shots[idx]->id=idx;
        int isx = (s->rcv[ir] / (s->dimy + 2 * s->sy)) - s->sx;
        int isy = (s->rcv[ir] % (s->dimy + 2 * s->sy)) - s->sy;
//        fprintf(fd2,"ir=%3d at %3d,shot id=%3d, isx=%d, isy=%d, isz=%d\n",ir,s->shots[idx]->srcidx,s->shots[idx]->id,isx,isy,isz);
        fprintf(fd2,"ir=%3d,shot id=%3d, isx=%d, isy=%d, isz=%d, isz_rcv=%d\n",ir,s->shots[idx]->id,isx,isy,s->src_depth,s->rcv_depth);
        ////////////////////////////////
////////////////////////////////
        idx++;
    }
    fclose(fd2);
}

void wave_init_acquisition_orig(sismap_t *s) {
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
//    int ir = -1;
    int ir=0;

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
            s->rcv[ir]=(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx);   //simwave:zyx
//            s->rcv[ir]=(s->sx+x)*(s->dimy+2*s->sy)+(y+s->sy);   //stencil:xyz
            fprintf(fd,"rec %3d in [%3d, %3d, %3d] at %3d\n",ir,x,y,s->rcv_depth,(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx));   //simwave:zyx
            ir++;
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
    if (s->verbose) wave_print(s);

    for (unsigned int idx = 0; idx < s->nb_shots; idx++)
        s->shots[idx] = (shot_t *) malloc(sizeof(shot_t));

    if (s->verbose) wave_print(s);
    unsigned int idx = 0;
    for (unsigned int ir = 0; ir < s->rcv_len; ir = ir + s->dshot) {
//        s->shots[idx]->srcidx = s->rcv[ir]+s->drcv/2;   // should I really use "s->drcv/2"?
        s->shots[idx]->srcidx=s->rcv[ir];
        s->shots[idx]->id=idx;
        fprintf(fd2,"ir=%3d at %3d,shot id=%3d\n",ir,s->shots[idx]->srcidx,s->shots[idx]->id);
        ////////////////////////////////
////////////////////////////////
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
    MSG("... stencil information:");
    MSG("... velocity size       = %u x %u x %u",
        s->vel_dimx, s->vel_dimy, s->vel_dimz);
    MSG("... velocity min/max    = %f - %f", s->vmin, s->vmax);
    MSG("... compute domain size = %u x %u x %u (%f MB)",
        s->dimx, s->dimy, s->dimz, s->size / 1024. / 1024.);
    MSG("... imaging domain size = %u x %u x %u (%f MB)",
        s->img_dimx, s->img_dimy, s->img_dimz,
        s->size_img / 1024. / 1024.);
    MSG("... cdp,line,depth      = %u x %u x %u", s->dcdp, s->dline, s->ddepth);
    MSG("... dx*dy*dz            = %f x %f x %f", s->dx, s->dy, s->dz);
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
ux[z] = 2.0f * vx[z] - ux[z]                                                   \
      + rx[z] * (coef0 * vx[z] + s->coefx[1] * (vx[z+1*nnyz] + vx[z-1*nnyz])   \
                               + s->coefy[1] * (vx[z+nnz] + vx[z-nnz])   \
                               + s->coefz[1] * (vx[z+1] + vx[z-1  ])   \
                               + s->coefx[2] * (vx[z+2*nnyz] + vx[z-2*nnyz])   \
                               + s->coefy[2] * (vx[z+2*nnz ] + vx[z-2*nnz ])   \
                               + s->coefz[2] * (vx[z+2] + vx[z-2])   \
                               + s->coefx[3] * (vx[z+3*nnyz] + vx[z-3*nnyz])   \
                               + s->coefy[3] * (vx[z+3*nnz ] + vx[z-3*nnz ])   \
                               + s->coefz[3] * (vx[z+3] + vx[z-3])   \
                               + s->coefx[4] * (vx[z+4*nnyz] + vx[z-4*nnyz])   \
                               + s->coefy[4] * (vx[z+4*nnz ] + vx[z-4*nnz ])   \
                               + s->coefz[4] * (vx[z+4] + vx[z-4]))  \

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v()                        \
vx0[z] = vx0[z]                                                                \
      + s->dt*inv_dx*(s->coefx[0] * (pr0[z+1*nnyz] - pr0[z])   \
                   + s->coefx[1] * (pr0[z+2*nnyz] - pr0[z-1*nnyz])   \
                   + s->coefx[2] * (pr0[z+3*nnyz] - pr0[z-2*nnyz])   \
                   + s->coefx[3] * (pr0[z+4*nnyz] - pr0[z-3*nnyz])) ;  \
vy0[z] = vy0[z]                                                                \
        + s->dt*inv_dy * (s->coefy[0] * (pr0[z+1*nnz] - pr0[z])   \
                       + s->coefy[1] * (pr0[z+2*nnz] - pr0[z-1*nnz])   \
                       + s->coefy[2] * (pr0[z+3*nnz] - pr0[z-2*nnz])   \
                       + s->coefy[3] * (pr0[z+4*nnz] - pr0[z-3*nnz])) ;  \
vz0[z] = vz0[z]                                                         \
      + s->dt*inv_dz * (          s->coefz[0] * (pr0[z+1] - pr0[z])   \
                               + s->coefz[1] * (pr0[z+2] - pr0[z-1])   \
                               + s->coefz[2] * (pr0[z+3] - pr0[z-2])   \
                               + s->coefz[3] * (pr0[z+4] - pr0[z-3]))

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p()                        \
pr0[z]=pr0[z]                                                   \
      + rx[z] * (                                                             \
                    s->coefx[0]*inv_dx * (vx0[z]           - vx0[z-1*nnyz     ])   \
                   + s->coefy[0]*inv_dy * (vy0[z]          - vy0[z-1*nnz   ])   \
                   + s->coefz[0]*inv_dz * (vz0[z]          - vz0[z-1  ])   \
                   + s->coefx[1]*inv_dx * (vx0[z+1*nnyz]   - vx0[z-2*nnyz])   \
                   + s->coefy[1]*inv_dy * (vy0[z+1*nnz ]   - vy0[z-2*nnz ])   \
                   + s->coefz[1]*inv_dz * (vz0[z+1]         -vz0[z-2])   \
                   + s->coefx[2]*inv_dx * (vx0[z+2*nnyz]    -vx0[z-3*nnyz])   \
                   + s->coefy[2]*inv_dy * (vy0[z+2*nnz ]    -vy0[z-3*nnz ])   \
                   + s->coefz[2]*inv_dz * (vz0[z+2]   -      vz0[z-3])   \
                   + s->coefx[3]*inv_dx * (vx0[z+3*nnyz]   - vx0[z-4*nnyz])   \
                   + s->coefy[3]*inv_dy * (vy0[z+3*nnz ]   - vy0[z-4*nnz ])   \
                   + s->coefz[3]*inv_dz * (vz0[z+ 3]   -     vz0[z-4])                 \
                   )

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v_index()                            \
    VX(z, y, x) = VX(z, y, x)+ s->dt*inv_dx * \
        (s->coefx[0]*( U0(z,   y,   x+1) - U0(z,y,x))  \
        + s->coefx[1]*( U0(z,   y,   x+2) - U0(z,   y,   x-1))  \
        + s->coefx[2]*( U0(z,   y,   x+3) - U0(z,   y,   x-2))  \
        + s->coefx[3]*( U0(z,   y,   x+4) - U0(z,   y,   x-3))  );\
    VY(z, y, x) = VY(z, y, x)+ s->dt*inv_dy *\
        (s->coefy[0]*( U0(z,   y+1, x  ) - U0(z,   y, x  ))  \
        + s->coefy[1]*( U0(z,   y+2, x  ) - U0(z,   y-1, x  ))  \
        + s->coefy[2]*( U0(z,   y+3, x  ) - U0(z,   y-2, x  ))  \
        + s->coefy[3]*( U0(z,   y+4, x  ) - U0(z,   y-3, x  )));\
    VZ(z, y, x) = VZ(z, y, x)+  s->dt*inv_dz *\
        (s->coefz[0]*( U0(z+1, y,   x  ) - U0(z, y,   x  ))  \
        + s->coefz[1]*( U0(z+2, y,   x  ) - U0(z-1, y,   x  ))  \
        + s->coefz[2]*( U0(z+3, y,   x  ) - U0(z-2, y,   x  ))  \
        + s->coefz[3]*( U0(z+4, y,   x  ) - U0(z-3, y,   x  )))

#define WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p_index() \
    U0(z, y, x) = U0(z, y, x)+ROC2(z, y, x)*  \
        ( s->coefx[0]*inv_dx*( VX(z,   y,   x) - VX(z,y,x-1))                \
        + s->coefy[0]*inv_dy*( VY(z,   y,   x) - VY(z,y-1,x))                \
        + s->coefz[0]*inv_dz*( VZ(z,   y,   x) - VZ(z-1,y,x))                \
        + s->coefx[1]*inv_dx*( VX(z,   y,   x+1) - VX(z,   y,   x-2))  \
        + s->coefy[1]*inv_dy*( VY(z,   y+1,   x) - VY(z,   y-2,   x))  \
        + s->coefz[1]*inv_dz*( VZ(z+1,   y,   x) - VZ(z-2,   y,   x))      \
        + s->coefx[2]*inv_dx*( VX(z,   y,   x+2) - VX(z,y,x-3))                \
        + s->coefy[2]*inv_dy*( VY(z,   y+2,   x) - VY(z,y-3,x))                \
        + s->coefz[2]*inv_dz*( VZ(z+2,   y,   x) - VZ(z-3,y,x))  \
        + s->coefx[3]*inv_dx*( VX(z,   y,   x+3) - VX(z,   y,   x-4))  \
        + s->coefy[3]*inv_dy*( VY(z,   y+3,   x) - VY(z,   y-4,   x))      \
        + s->coefz[3]*inv_dz*( VZ(z+3,   y,   x) - VZ(z-4,y,x)));          \
    U0(z, y, x) = U0(z, y, x)*DAMPX(z, y, x)*DAMPY(z, y, x)*DAMPZ(z, y, x)

#define WAVE_UPDATE_PML_FIELDS()                                            \
  U1(z,y,x) =                                                               \
    ((2.-ETA(z,y,x)*ETA(z,y,x) + 2.*ETA(z,y,x))*U0(z,y,x)                   \
    - U1(z,y,x) + ROC2(z,y,x)*(laplacian + PHI(z,y,x)))/(1.+2.*ETA(z,y,x)   ); \
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
    const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
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
    const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
    const float inv_dx=1. /(s->dx);
    const float inv_dy=1. /(s->dy);
    const float inv_dz=1. /(s->dz);

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
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
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
    const int BLOCKX=s->blockx;
    const int BLOCKY=s->blocky;
    const int BLOCKZ=s->blockz;
    unsigned int z, y, x;
    float laplacian;
    const float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
    unsigned int xmin,xmax,zmin,zmax,ymin,ymax;

    const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const int sx = s->sx;
	const int sy = s->sy;
	const int sz = s->sz;
    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnz = s->dimz + 2 * s->sz;
    const long int nnxy=(long int)nnx * nny; // XYZ order: x-slowest, z-fastest
    const long int nnyz=(long int)nny * nnz;

    const float inv_dx = 1. / (s->dx);
	const float inv_dy = 1. / (s->dy);
	const float inv_dz = 1. / (s->dz);

	// Precompute coefficients with dt for velocity updates
	const float dt_inv_dx = s->dt * inv_dx;
	const float dt_inv_dy = s->dt * inv_dy;
	const float dt_inv_dz = s->dt * inv_dz;

	// Host coefficient arrays outside the loop
	const float *restrict coefx = s->coefx;
	const float *restrict coefy = s->coefy;
	const float *restrict coefz = s->coefz;
	const float *restrict dampx = s->dampx;
	const float *restrict dampy = s->dampy;
	const float *restrict dampz = s->dampz;

    float *restrict ux;
    float *restrict vx;
    float *restrict rx;
#pragma omp parallel for collapse(3) schedule(dynamic) private(laplacian, xmin, xmax, zmin, zmax, ymin, ymax, ux, vx, rx)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
            	const Myint xmax = fmin(dimx, xmin + BLOCKX);
				const Myint ymax = fmin(dimy, ymin + BLOCKY);
				const Myint zmax = fmin(dimz, zmin + BLOCKZ);
                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
//                        s->rcv[ir]=(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx);   //simwave:zyx
                        ux = &(u1[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx = &(u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        rx = &(roc2[1ULL*x*dimy*dimz + y*dimz]);
						#pragma omp simd
                        for (int z = zmin; z < zmax; z++) {
///////                            WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD();

                            ux[z] = 2.0f * vx[z] - ux[z]                                                   \
                                  + rx[z] * (coef0 * vx[z] + \
                                		  + coefx[1] * (vx[z+1*nnyz] + vx[z-1*nnyz])   \
								   + coefy[1] * (vx[z+nnz] + vx[z-nnz])   \
								   + coefz[1] * (vx[z+1] + vx[z-1  ])   \
								   + coefx[2] * (vx[z+2*nnyz] + vx[z-2*nnyz])   \
								   + coefy[2] * (vx[z+2*nnz ] + vx[z-2*nnz ])   \
								   + coefz[2] * (vx[z+2] + vx[z-2])   \
								   + coefx[3] * (vx[z+3*nnyz] + vx[z-3*nnyz])   \
								   + coefy[3] * (vx[z+3*nnz ] + vx[z-3*nnz ])   \
								   + coefz[3] * (vx[z+3] + vx[z-3])   \
								   + coefx[4] * (vx[z+4*nnyz] + vx[z-4*nnyz])   \
								   + coefy[4] * (vx[z+4*nnz ] + vx[z-4*nnz ])   \
								   + coefz[4] * (vx[z+4] + vx[z-4]));
                            ux[z] = dampx[x + sx] * ux[z] + (1 - dampx[x + sx]) * vx[z];
                            ux[z] = dampy[y + sy] * ux[z] + (1 - dampy[y + sy]) * vx[z];
                            ux[z] = dampz[z + sz] * ux[z] + (1 - dampz[z + sz]) * vx[z];
                        }
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_bis_better(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict u1,
                                  float *restrict roc2) {
    const int BLOCKX = s->blockx;
    const int BLOCKY = s->blocky;
    const int BLOCKZ = s->blockz;
    const float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    const int sx = s->sx;
    const int sy = s->sy;
    const int sz = s->sz;
    const int nnx = dimx + 2 * sx;
    const int nny = dimy + 2 * sy;
    const int nnz = dimz + 2 * sz;
    const long int nnyz = (long int)nny * nnz;

    const float *restrict coefx = s->coefx;
    const float *restrict coefy = s->coefy;
    const float *restrict coefz = s->coefz;
    const float *restrict dampx = s->dampx;
    const float *restrict dampy = s->dampy;
    const float *restrict dampz = s->dampz;
    float *restrict ux;
    float *restrict vx;
    float *restrict rx;

#pragma omp parallel for collapse(3) schedule(static) private(ux, vx, rx)
    for (int xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (int ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (int zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);

                for (int x = xmin; x < xmax; x++) {
                    const long long base_x = (x + sx) * nnyz + sz;
                    float *restrict ux_base = u1 + base_x;
                    float *restrict vx_base = u0 + base_x;
                    float *restrict rx_base = roc2 + x * dimy * dimz;
                    const float damp_x = dampx[x + sx];

                    for (int y = ymin; y < ymax; y++) {
                        const long long offset_y = (y + sy) * nnz;
                        float *restrict ux = ux_base + offset_y;
                        float *restrict vx = vx_base + offset_y;
                        float *restrict rx = rx_base + y * dimz;
                        const float damp_xy = damp_x * dampy[y + sy];

#pragma omp simd
                        for (int z = zmin; z < zmax; z += 4) {
                            for (int i = 0; i < 4 && z + i < zmax; i++) {
                                const int zi = z + i;
                                ux[zi] = 2.0f * vx[zi] - ux[zi] + rx[zi] * (
                                    coef0 * vx[zi] +
                                    coefx[1] * (vx[zi + nnyz] + vx[zi - nnyz]) +
                                    coefy[1] * (vx[zi + nnz] + vx[zi - nnz]) +
                                    coefz[1] * (vx[zi + 1] + vx[zi - 1]) +
                                    coefx[2] * (vx[zi + 2 * nnyz] + vx[zi - 2 * nnyz]) +
                                    coefy[2] * (vx[zi + 2 * nnz] + vx[zi - 2 * nnz]) +
                                    coefz[2] * (vx[zi + 2] + vx[zi - 2]) +
                                    coefx[3] * (vx[zi + 3 * nnyz] + vx[zi - 3 * nnyz]) +
                                    coefy[3] * (vx[zi + 3 * nnz] + vx[zi - 3 * nnz]) +
                                    coefz[3] * (vx[zi + 3] + vx[zi - 3]) +
                                    coefx[4] * (vx[zi + 4 * nnyz] + vx[zi - 4 * nnyz]) +
                                    coefy[4] * (vx[zi + 4 * nnz] + vx[zi - 4 * nnz]) +
                                    coefz[4] * (vx[zi + 4] + vx[zi - 4])
                                );

                                // Optimized damping
                                const float damp_z = dampz[zi + sz];
                                const float damp_xyz = damp_xy * damp_z;
                                const float v_coeff = (1.0f - damp_x) + (1.0f - dampy[y + sy]) * damp_x + (1.0f - damp_z) * damp_xy;
                                ux[zi] = damp_xyz * ux[zi] + v_coeff * vx[zi];
                            }
                        }
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_bis_old(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict u1,
                                  float *restrict roc2) {
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	unsigned int z, y, x;
	float laplacian;
	float coef0 = s->coefx[0] + s->coefy[0] + s->coefz[0];
	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;

	const int nnx = s->dimx + 2 * s->sx;
	const int nny = s->dimy + 2 * s->sy;
	const int nnz = s->dimz + 2 * s->sz;
	const int nnxy = nnx * nny;
	const int nnyz = nny * nnz;

	float *restrict ux;
	float *restrict vx;
	float *restrict rx;
#pragma omp parallel for collapse(3) private(laplacian, xmin, xmax, zmin, zmax, ymin, ymax, ux, vx, rx)
	for (xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
		for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
			for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
				zmax = zmin + BLOCKZ;
				if (zmax > s->dimz) zmax = s->dimz;
				ymax = ymin + BLOCKY;
				if (ymax > s->dimy) ymax = s->dimy;
				xmax = xmin + BLOCKX;
				if (xmax > s->dimx) xmax = s->dimx;

				for (int x = xmin; x < xmax; x++) {
					for (int y = ymin; y < ymax; y++) {
//                        s->rcv[ir]=(s->sy+y)*(s->dimx+2*s->sx)+(x+s->sx);   //simwave:zyx
						ux = &(u1[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
						vx = &(u0[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
						rx = &(roc2[1ULL * x * s->dimy * s->dimz + y * s->dimz]);
						for (int z = zmin; z < zmax; z++) {
							WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD();
							ux[z] = s->dampx[x + s->sx] * ux[z] + (1 - s->dampx[x + s->sx]) * vx[z];
							ux[z] = s->dampy[y + s->sy] * ux[z] + (1 - s->dampy[y + s->sy]) * vx[z];
							ux[z] = s->dampz[z + s->sz] * ux[z] + (1 - s->dampz[z + s->sz]) * vx[z];
						}
					}
				}
			}
		}
	}
}

void wave_update_fields_block_1st_orig(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict vx,
                                  float *restrict vy,
                                  float *restrict vz,
                                  float *restrict roc2,
                                  float *restrict phi,
                                  float *restrict eta) {
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	unsigned int z, y, x;
    float laplacian;
    unsigned int xmin,xmax,zmin,zmax,ymin,ymax;

    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnz = s->dimz + 2 * s->sz;
    const float inv_dx=1. /(s->dx);
    const float inv_dy=1. /(s->dy);
    const float inv_dz=1. /(s->dz);
    long int nnxy = nnx * nny;
    long int nnyz = nny * nnz;
    float *restrict pr0;
    float *restrict vx0;
    float *restrict vy0;
    float *restrict vz0;
    float *restrict rx;
//    MSG("BLOCKX=%d, BLOCKY=%d, BLOCKZ=%d\n",BLOCKX,BLOCKY,BLOCKZ);
//    exit(0);
    // loop on the blocks .velocity
#pragma omp parallel for collapse(3) schedule(dynamic) private(xmin,xmax,zmin,zmax,ymin,ymax,pr0,vx0,vy0,vz0)
    for (xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(s->dimx, xmin + BLOCKX);
                const Myint ymax = fmin(s->dimy, ymin + BLOCKY);
                const Myint zmax = fmin(s->dimz, zmin + BLOCKZ);
                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
//                        MSG("x=%d, y=%d, z=%d\n",x,y,z);
                        pr0 = &(u0[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vx0 = &(vx[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vy0 = &(vy[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vz0 = &(vz[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
#pragma ivdep
                        for (int z = zmin; z < zmax; z++) {
                            // WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v();
                            vx0[z] = vx0[z]+ s->dt*inv_dx*(s->coefx[0] * (pr0[z+1*nnyz] - pr0[z])   \
                               + s->coefx[1] * (pr0[z+2*nnyz] - pr0[z-1*nnyz])   \
                               + s->coefx[2] * (pr0[z+3*nnyz] - pr0[z-2*nnyz])   \
                               + s->coefx[3] * (pr0[z+4*nnyz] - pr0[z-3*nnyz]));  \
                            vy0[z] = vy0[z]+ s->dt*inv_dy * (s->coefy[0] * (pr0[z+1*nnz] - pr0[z])   \
                                                   + s->coefy[1] * (pr0[z+2*nnz] - pr0[z-1*nnz])   \
                                                   + s->coefy[2] * (pr0[z+3*nnz] - pr0[z-2*nnz])   \
                                                   + s->coefy[3] * (pr0[z+4*nnz] - pr0[z-3*nnz]));  \
                            vz0[z] = vz0[z]+ s->dt*inv_dz * (s->coefz[0] * (pr0[z+1] - pr0[z])   \
                                                           + s->coefz[1] * (pr0[z+2] - pr0[z-1])   \
                                                           + s->coefz[2] * (pr0[z+3] - pr0[z-2])   \
                                                           + s->coefz[3] * (pr0[z+4] - pr0[z-3]));
                        }
                    }
                }
            }
        }
    }

    // loop on the blocks.pressure
#pragma omp parallel for collapse(3) schedule(dynamic) private(laplacian,xmin,xmax,zmin,zmax,ymin,ymax,pr0,vx0,vy0,vz0,rx)
    for (xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(s->dimx, xmin + BLOCKX);
                const Myint ymax = fmin(s->dimy, ymin + BLOCKY);
                const Myint zmax = fmin(s->dimz, zmin + BLOCKZ);

                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
                        pr0 = &(u0[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vx0 = &(vx[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vy0 = &(vy[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        vz0 = &(vz[1ULL * (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz]);
                        rx = &(roc2[1ULL * x * s->dimy * s->dimz + y * s->dimz]);
#pragma ivdep
                        for (int z = zmin; z < zmax; z++) {
                            // WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p();
                            pr0[z]=pr0[z]+ rx[z] * ( \
                                s->coefx[0]*inv_dx * (vx0[z]           - vx0[z-1*nnyz     ])   \
                               + s->coefy[0]*inv_dy * (vy0[z]          - vy0[z-1*nnz   ])   \
                               + s->coefz[0]*inv_dz * (vz0[z]          - vz0[z-1  ])   \
                               + s->coefx[1]*inv_dx * (vx0[z+1*nnyz]   - vx0[z-2*nnyz])   \
                               + s->coefy[1]*inv_dy * (vy0[z+1*nnz ]   - vy0[z-2*nnz ])   \
                               + s->coefz[1]*inv_dz * (vz0[z+1]         -vz0[z-2])   \
                               + s->coefx[2]*inv_dx * (vx0[z+2*nnyz]    -vx0[z-3*nnyz])   \
                               + s->coefy[2]*inv_dy * (vy0[z+2*nnz ]    -vy0[z-3*nnz ])   \
                               + s->coefz[2]*inv_dz * (vz0[z+2]   -      vz0[z-3])   \
                               + s->coefx[3]*inv_dx * (vx0[z+3*nnyz]   - vx0[z-4*nnyz])   \
                               + s->coefy[3]*inv_dy * (vy0[z+3*nnz ]   - vy0[z-4*nnz ])   \
                               + s->coefz[3]*inv_dz * (vz0[z+ 3]   -     vz0[z-4]) );
                            pr0[z] = s->dampx[x + s->sx]*s->dampy[y + s->sy]*s->dampz[z + s->sz] * pr0[z];
                        }
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_1st_(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict vx,
                                  float *restrict vy,
                                  float *restrict vz,
                                  float *restrict roc2,
                                  float *restrict phi __attribute__((unused)),
                                  float *restrict eta __attribute__((unused))) {
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnz = s->dimz + 2 * s->sz;
    const long int nnxy = (long int)nnx * nny;
    const long int nnyz = (long int)nny * nnz;

    // Precompute coefficients with dt and inverse distances
    const float dt_inv_dx = s->dt / s->dx;
    const float dt_inv_dy = s->dt / s->dy;
    const float dt_inv_dz = s->dt / s->dz;
    float coefx[4] __attribute__((aligned(ALIGNMENT))) = {
        s->coefx[0] * dt_inv_dx, s->coefx[1] * dt_inv_dx,
        s->coefx[2] * dt_inv_dx, s->coefx[3] * dt_inv_dx
    };
    float coefy[4] __attribute__((aligned(ALIGNMENT))) = {
        s->coefy[0] * dt_inv_dy, s->coefy[1] * dt_inv_dy,
        s->coefy[2] * dt_inv_dy, s->coefy[3] * dt_inv_dy
    };
    float coefz[4] __attribute__((aligned(ALIGNMENT))) = {
        s->coefz[0] * dt_inv_dz, s->coefz[1] * dt_inv_dz,
        s->coefz[2] * dt_inv_dz, s->coefz[3] * dt_inv_dz
    };

    // Ensure arrays are aligned (assumes caller allocates aligned memory)
    float *restrict pr0 = __builtin_assume_aligned(u0, ALIGNMENT);
    float *restrict vx0 = __builtin_assume_aligned(vx, ALIGNMENT);
    float *restrict vy0 = __builtin_assume_aligned(vy, ALIGNMENT);
    float *restrict vz0 = __builtin_assume_aligned(vz, ALIGNMENT);
    float *restrict rx = __builtin_assume_aligned(roc2, ALIGNMENT);

    #pragma omp parallel
    {
        // Single parallel region to reduce overhead
        #pragma omp for collapse(3) schedule(static) nowait
        for (int xmin = 0; xmin < s->dimx; xmin += BLOCKX) {
            for (int ymin = 0; ymin < s->dimy; ymin += BLOCKY) {
                for (int zmin = 0; zmin < s->dimz; zmin += BLOCKZ) {
                    const int xmax = xmin + BLOCKX < s->dimx ? xmin + BLOCKX : s->dimx;
                    const int ymax = ymin + BLOCKY < s->dimy ? ymin + BLOCKY : s->dimy;
                    const int zmax = zmin + BLOCKZ < s->dimz ? zmin + BLOCKZ : s->dimz;

                    for (int x = xmin; x < xmax; x++) {
                        for (int y = ymin; y < ymax; y++) {
                            const long int base_idx = (x + s->sx) * nnyz + (y + s->sy) * nnz + s->sz;
                            float *restrict p = pr0 + base_idx;
                            float *restrict vx_p = vx0 + base_idx;
                            float *restrict vy_p = vy0 + base_idx;
                            float *restrict vz_p = vz0 + base_idx;
                            float *restrict r = rx + (x * s->dimy * s->dimz + y * s->dimz);
                            const float damp_xy = s->dampx[x + s->sx] * s->dampy[y + s->sy];

                            // Vectorized inner loop (assuming zmax - zmin >= 8 for AVX2)
                            int z;
                            for (z = zmin; z <= zmax - 8; z += 8) {
                                // Load 8 elements at once
                                __m256 p_vec = _mm256_load_ps(p + z);
                                __m256 vx_vec = _mm256_load_ps(vx_p + z);
                                __m256 vy_vec = _mm256_load_ps(vy_p + z);
                                __m256 vz_vec = _mm256_load_ps(vz_p + z);
                                __m256 r_vec = _mm256_load_ps(r + z);
                                __m256 damp_z = _mm256_load_ps(&s->dampz[z + s->sz]);

                                // Velocity updates (x-direction)
                                __m256 px0 = _mm256_load_ps(p + z + 1 * nnyz);
                                __m256 px1 = _mm256_load_ps(p + z - 1 * nnyz);
                                __m256 px2 = _mm256_load_ps(p + z + 2 * nnyz);
                                __m256 px3 = _mm256_load_ps(p + z - 2 * nnyz);
                                __m256 px4 = _mm256_load_ps(p + z + 3 * nnyz);
                                __m256 px5 = _mm256_load_ps(p + z - 3 * nnyz);
                                __m256 px6 = _mm256_load_ps(p + z + 4 * nnyz);
                                __m256 px7 = _mm256_load_ps(p + z - 4 * nnyz);
                                vx_vec = _mm256_add_ps(vx_vec, _mm256_mul_ps(_mm256_set1_ps(coefx[0]), _mm256_sub_ps(px0, p_vec)));
                                vx_vec = _mm256_add_ps(vx_vec, _mm256_mul_ps(_mm256_set1_ps(coefx[1]), _mm256_sub_ps(px2, px1)));
                                vx_vec = _mm256_add_ps(vx_vec, _mm256_mul_ps(_mm256_set1_ps(coefx[2]), _mm256_sub_ps(px4, px3)));
                                vx_vec = _mm256_add_ps(vx_vec, _mm256_mul_ps(_mm256_set1_ps(coefx[3]), _mm256_sub_ps(px6, px7)));

                                // Velocity updates (y-direction)
                                __m256 py0 = _mm256_load_ps(p + z + 1 * nnz);
                                __m256 py1 = _mm256_load_ps(p + z - 1 * nnz);
                                __m256 py2 = _mm256_load_ps(p + z + 2 * nnz);
                                __m256 py3 = _mm256_load_ps(p + z - 2 * nnz);
                                __m256 py4 = _mm256_load_ps(p + z + 3 * nnz);
                                __m256 py5 = _mm256_load_ps(p + z - 3 * nnz);
                                __m256 py6 = _mm256_load_ps(p + z + 4 * nnz);
                                __m256 py7 = _mm256_load_ps(p + z - 4 * nnz);
                                vy_vec = _mm256_add_ps(vy_vec, _mm256_mul_ps(_mm256_set1_ps(coefy[0]), _mm256_sub_ps(py0, p_vec)));
                                vy_vec = _mm256_add_ps(vy_vec, _mm256_mul_ps(_mm256_set1_ps(coefy[1]), _mm256_sub_ps(py2, py1)));
                                vy_vec = _mm256_add_ps(vy_vec, _mm256_mul_ps(_mm256_set1_ps(coefy[2]), _mm256_sub_ps(py4, py3)));
                                vy_vec = _mm256_add_ps(vy_vec, _mm256_mul_ps(_mm256_set1_ps(coefy[3]), _mm256_sub_ps(py6, py7)));

                                // Velocity updates (z-direction)
                                __m256 pz0 = _mm256_load_ps(p + z + 1);
                                __m256 pz1 = _mm256_load_ps(p + z - 1);
                                __m256 pz2 = _mm256_load_ps(p + z + 2);
                                __m256 pz3 = _mm256_load_ps(p + z - 2);
                                __m256 pz4 = _mm256_load_ps(p + z + 3);
                                __m256 pz5 = _mm256_load_ps(p + z - 3);
                                __m256 pz6 = _mm256_load_ps(p + z + 4);
                                __m256 pz7 = _mm256_load_ps(p + z - 4);
                                vz_vec = _mm256_add_ps(vz_vec, _mm256_mul_ps(_mm256_set1_ps(coefz[0]), _mm256_sub_ps(pz0, p_vec)));
                                vz_vec = _mm256_add_ps(vz_vec, _mm256_mul_ps(_mm256_set1_ps(coefz[1]), _mm256_sub_ps(pz2, pz1)));
                                vz_vec = _mm256_add_ps(vz_vec, _mm256_mul_ps(_mm256_set1_ps(coefz[2]), _mm256_sub_ps(pz4, pz3)));
                                vz_vec = _mm256_add_ps(vz_vec, _mm256_mul_ps(_mm256_set1_ps(coefz[3]), _mm256_sub_ps(pz6, pz7)));

                                // Pressure updates
                                __m256 vx_diff0 = _mm256_sub_ps(vx_vec, _mm256_load_ps(vx_p + z - 1 * nnyz));
                                __m256 vy_diff0 = _mm256_sub_ps(vy_vec, _mm256_load_ps(vy_p + z - 1 * nnz));
                                __m256 vz_diff0 = _mm256_sub_ps(vz_vec, _mm256_load_ps(vz_p + z - 1));
                                __m256 p_update = _mm256_mul_ps(_mm256_set1_ps(coefx[0]), vx_diff0);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefy[0]), vy_diff0, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefz[0]), vz_diff0, p_update);

                                __m256 vx_diff1 = _mm256_sub_ps(_mm256_load_ps(vx_p + z + 1 * nnyz), _mm256_load_ps(vx_p + z - 2 * nnyz));
                                __m256 vy_diff1 = _mm256_sub_ps(_mm256_load_ps(vy_p + z + 1 * nnz), _mm256_load_ps(vy_p + z - 2 * nnz));
                                __m256 vz_diff1 = _mm256_sub_ps(_mm256_load_ps(vz_p + z + 1), _mm256_load_ps(vz_p + z - 2));
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefx[1]), vx_diff1, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefy[1]), vy_diff1, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefz[1]), vz_diff1, p_update);

                                __m256 vx_diff2 = _mm256_sub_ps(_mm256_load_ps(vx_p + z + 2 * nnyz), _mm256_load_ps(vx_p + z - 3 * nnyz));
                                __m256 vy_diff2 = _mm256_sub_ps(_mm256_load_ps(vy_p + z + 2 * nnz), _mm256_load_ps(vy_p + z - 3 * nnz));
                                __m256 vz_diff2 = _mm256_sub_ps(_mm256_load_ps(vz_p + z + 2), _mm256_load_ps(vz_p + z - 3));
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefx[2]), vx_diff2, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefy[2]), vy_diff2, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefz[2]), vz_diff2, p_update);

                                __m256 vx_diff3 = _mm256_sub_ps(_mm256_load_ps(vx_p + z + 3 * nnyz), _mm256_load_ps(vx_p + z - 4 * nnyz));
                                __m256 vy_diff3 = _mm256_sub_ps(_mm256_load_ps(vy_p + z + 3 * nnz), _mm256_load_ps(vy_p + z - 4 * nnz));
                                __m256 vz_diff3 = _mm256_sub_ps(_mm256_load_ps(vz_p + z + 3), _mm256_load_ps(vz_p + z - 4));
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefx[3]), vx_diff3, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefy[3]), vy_diff3, p_update);
                                p_update = _mm256_fmadd_ps(_mm256_set1_ps(coefz[3]), vz_diff3, p_update);

                                p_vec = _mm256_add_ps(p_vec, _mm256_mul_ps(r_vec, p_update));
                                p_vec = _mm256_mul_ps(p_vec, _mm256_mul_ps(_mm256_set1_ps(damp_xy), damp_z));

                                // Store results
                                _mm256_store_ps(vx_p + z, vx_vec);
                                _mm256_store_ps(vy_p + z, vy_vec);
                                _mm256_store_ps(vz_p + z, vz_vec);
                                _mm256_store_ps(p + z, p_vec);
                            }

                            // Handle remainder
                            for (; z < zmax; z++) {
                                vx_p[z] += coefx[0] * (p[z + 1 * nnyz] - p[z]) +
                                          coefx[1] * (p[z + 2 * nnyz] - p[z - 1 * nnyz]) +
                                          coefx[2] * (p[z + 3 * nnyz] - p[z - 2 * nnyz]) +
                                          coefx[3] * (p[z + 4 * nnyz] - p[z - 3 * nnyz]);
                                vy_p[z] += coefy[0] * (p[z + 1 * nnz] - p[z]) +
                                          coefy[1] * (p[z + 2 * nnz] - p[z - 1 * nnz]) +
                                          coefy[2] * (p[z + 3 * nnz] - p[z - 2 * nnz]) +
                                          coefy[3] * (p[z + 4 * nnz] - p[z - 3 * nnz]);
                                vz_p[z] += coefz[0] * (p[z + 1] - p[z]) +
                                          coefz[1] * (p[z + 2] - p[z - 1]) +
                                          coefz[2] * (p[z + 3] - p[z - 2]) +
                                          coefz[3] * (p[z + 4] - p[z - 3]);

                                float p_update = coefx[0] * (vx_p[z] - vx_p[z - 1 * nnyz]) +
                                                coefy[0] * (vy_p[z] - vy_p[z - 1 * nnz]) +
                                                coefz[0] * (vz_p[z] - vz_p[z - 1]) +
                                                coefx[1] * (vx_p[z + 1 * nnyz] - vx_p[z - 2 * nnyz]) +
                                                coefy[1] * (vy_p[z + 1 * nnz] - vy_p[z - 2 * nnz]) +
                                                coefz[1] * (vz_p[z + 1] - vz_p[z - 2]) +
                                                coefx[2] * (vx_p[z + 2 * nnyz] - vx_p[z - 3 * nnyz]) +
                                                coefy[2] * (vy_p[z + 2 * nnz] - vy_p[z - 3 * nnz]) +
                                                coefz[2] * (vz_p[z + 2] - vz_p[z - 3]) +
                                                coefx[3] * (vx_p[z + 3 * nnyz] - vx_p[z - 4 * nnyz]) +
                                                coefy[3] * (vy_p[z + 3 * nnz] - vy_p[z - 4 * nnz]) +
                                                coefz[3] * (vz_p[z + 3] - vz_p[z - 4]);
                                p[z] += r[z] * p_update;
                                p[z] *= damp_xy * s->dampz[z + s->sz];
                            }
                        }
                    }
                }
            }
        }
    }
}

void wave_update_fields_block_1st_v2(sismap_t *s,
                                  float *restrict u0,
                                  float *restrict vx,
                                  float *restrict vy,
                                  float *restrict vz,
                                  float *restrict roc2,
                                  float *restrict phi,
                                  float *restrict eta) {
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	unsigned int z, y, x;
    float laplacian;
    unsigned int xmin,xmax,zmin,zmax,ymin,ymax;

    const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const int sx = s->sx;
	const int sy = s->sy;
	const int sz = s->sz;
	const int nnx = dimx + 2 * sx;
	const int nny = dimy + 2 * sy;
	const int nnz = dimz + 2 * sz;

    const float inv_dx=1. /(s->dx);
    const float inv_dy=1. /(s->dy);
    const float inv_dz=1. /(s->dz);

    const long int nnyz = (long int)nnx * nny; // XYZ order: x-slowest, z-fastest
    const long int nnxy = (long int)nny * nnz;

    // Precompute coefficients with dt for velocity updates
	const float dt_inv_dx = s->dt * inv_dx;
	const float dt_inv_dy = s->dt * inv_dy;
	const float dt_inv_dz = s->dt * inv_dz;

	// Hoist coefficient arrays outside the loop
	const float *restrict coefx = s->coefx;
	const float *restrict coefy = s->coefy;
	const float *restrict coefz = s->coefz;
	const float *restrict dampx = s->dampx;
	const float *restrict dampy = s->dampy;
	const float *restrict dampz = s->dampz;

    float *restrict pr0;
    float *restrict vx0;
    float *restrict vy0;
    float *restrict vz0;
    float *restrict rx;
//    MSG("BLOCKX=%d, BLOCKY=%d, BLOCKZ=%d\n",BLOCKX,BLOCKY,BLOCKZ);
//    exit(0);
    // loop on the blocks .velocity
#pragma omp parallel for collapse(3) schedule(dynamic) private(laplacian,xmin,xmax,zmin,zmax,ymin,ymax,pr0,vx0,vy0,vz0)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);
                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
//                        MSG("x=%d, y=%d, z=%d\n",x,y,z);
                        pr0 = &(u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx0 = &(vx[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vy0 = &(vy[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vz0 = &(vz[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
#pragma ivdep
                        for (int z = zmin; z < zmax; z++) {
                            // WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_v();
                            vx0[z] = vx0[z]+ dt_inv_dx*(coefx[0] * (pr0[z+1*nnyz] - pr0[z])   \
                               + coefx[1] * (pr0[z+2*nnyz] - pr0[z-1*nnyz])   \
                               + coefx[2] * (pr0[z+3*nnyz] - pr0[z-2*nnyz])   \
                               + coefx[3] * (pr0[z+4*nnyz] - pr0[z-3*nnyz]));  \
                            vy0[z] = vy0[z]+ dt_inv_dy * (coefy[0] * (pr0[z+1*nnz] - pr0[z])   \
                                                   + coefy[1] * (pr0[z+2*nnz] - pr0[z-1*nnz])   \
                                                   + coefy[2] * (pr0[z+3*nnz] - pr0[z-2*nnz])   \
                                                   + coefy[3] * (pr0[z+4*nnz] - pr0[z-3*nnz]));  \
                            vz0[z] = vz0[z]+ dt_inv_dz * (coefz[0] * (pr0[z+1] - pr0[z])   \
                                                           + coefz[1] * (pr0[z+2] - pr0[z-1])   \
                                                           + coefz[2] * (pr0[z+3] - pr0[z-2])   \
                                                           + coefz[3] * (pr0[z+4] - pr0[z-3]));
                        }
                    }
                }
            }
        }
    }

    // loop on the blocks.pressure
#pragma omp parallel for collapse(3) schedule(dynamic) private(laplacian,xmin,xmax,zmin,zmax,ymin,ymax,pr0,vx0,vy0,vz0,rx)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);

                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
                        pr0 = &(u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx0 = &(vx[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vy0 = &(vy[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vz0 = &(vz[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        rx = &(roc2[1ULL * x * dimy * dimz + y * dimz]);
#pragma ivdep
                        for (int z = zmin; z < zmax; z++) {
                            // WAVE_COMPUTE_LAPLACIAN_AND_UPDATE_INNER_FIELD_1st_p();
                            pr0[z]=pr0[z]+ rx[z] * ( \
                                 coefx[0]*inv_dx * (vx0[z]           - vx0[z-1*nnyz     ])   \
                               + coefy[0]*inv_dy * (vy0[z]          - vy0[z-1*nnz   ])   \
                               + coefz[0]*inv_dz * (vz0[z]          - vz0[z-1  ])   \
                               + coefx[1]*inv_dx * (vx0[z+1*nnyz]   - vx0[z-2*nnyz])   \
                               + coefy[1]*inv_dy * (vy0[z+1*nnz ]   - vy0[z-2*nnz ])   \
                               + coefz[1]*inv_dz * (vz0[z+1]         -vz0[z-2])   \
                               + coefx[2]*inv_dx * (vx0[z+2*nnyz]    -vx0[z-3*nnyz])   \
                               + coefy[2]*inv_dy * (vy0[z+2*nnz ]    -vy0[z-3*nnz ])   \
                               + coefz[2]*inv_dz * (vz0[z+2]   -      vz0[z-3])   \
                               + coefx[3]*inv_dx * (vx0[z+3*nnyz]   - vx0[z-4*nnyz])   \
                               + coefy[3]*inv_dy * (vy0[z+3*nnz ]   - vy0[z-4*nnz ])   \
                               + coefz[3]*inv_dz * (vz0[z+ 3]   -     vz0[z-4]) );
                            pr0[z] = dampx[x + sx]*dampy[y + sy]*dampz[z+sz] * pr0[z];
                        }
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
								 float *restrict inv_rho) {
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	unsigned int z, y, x;

	unsigned int xmin, xmax, zmin, zmax, ymin, ymax;
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    const int sx = s->sx;
    const int sy = s->sy;
    const int sz = s->sz;
    const int nnx = dimx + 2 * sx;
    const int nny = dimy + 2 * sy;
    const int nnz = dimz + 2 * sz;
    const long int nnxy=(long int)nnx * nny; // XYZ order: x-slowest, z-fastest
    const long int nnyz=(long int)nny * nnz;

    const float inv_dx = 1. / (s->dx);
    const float inv_dy = 1. / (s->dy);
    const float inv_dz = 1. / (s->dz);
//    const unsigned long long  nnxy=nnx*nny; // XYZ order: x-slowest, z-fastest
//    const unsigned long long  nnyz=nny*nnz;

    // Precompute coefficients with dt for velocity updates
    const float dt_inv_dx = s->dt * inv_dx;
    const float dt_inv_dy = s->dt * inv_dy;
    const float dt_inv_dz = s->dt * inv_dz;

//    MSG("dt_inv_dx=%f\n",dt_inv_dx);
//    MSG("nnyz=%llu\n",nnyz);


    // Host coefficient arrays outside the loop
    const float *restrict coefx = s->coefx;
    const float *restrict coefy = s->coefy;
    const float *restrict coefz = s->coefz;
    const float *restrict dampx = s->dampx;
    const float *restrict dampy = s->dampy;
    const float *restrict dampz = s->dampz;

    float *restrict pr0;
    float *restrict vx0;
    float *restrict vy0;
    float *restrict vz0;
    float *restrict rx;
    float *restrict inv_r;

    // Velocity update loop
    #pragma omp parallel for collapse(3) schedule(dynamic) private(xmin, xmax, zmin, zmax, ymin, ymax, pr0, vx0, vy0, vz0)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);
                for (int x = xmin; x < xmax; x++) {
//                	MSG("vz0[z]=%f\n",vz0[z]);
                    for (int y = ymin; y < ymax; y++) {
                        pr0 = &(u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx0 = &(vx[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vy0 = &(vy[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vz0 = &(vz[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        inv_r = &(inv_rho[1ULL * x * dimy * dimz + y * dimz]);
                        #pragma omp simd
                        for (int z = zmin; z < zmax; z++) {
                            vx0[z] = vx0[z] + inv_r[z]*dt_inv_dx * (coefx[0] * (pr0[z + 1 * nnyz] - pr0[z]) +
                                                          coefx[1] * (pr0[z + 2 * nnyz] - pr0[z - 1 * nnyz]) +
                                                          coefx[2] * (pr0[z + 3 * nnyz] - pr0[z - 2 * nnyz]) +
                                                          coefx[3] * (pr0[z + 4 * nnyz] - pr0[z - 3 * nnyz]));
                            vy0[z] = vy0[z] + inv_r[z]*dt_inv_dy * (coefy[0] * (pr0[z + 1 * nnz] - pr0[z]) +
                                                          coefy[1] * (pr0[z + 2 * nnz] - pr0[z - 1 * nnz]) +
                                                          coefy[2] * (pr0[z + 3 * nnz] - pr0[z - 2 * nnz]) +
                                                          coefy[3] * (pr0[z + 4 * nnz] - pr0[z - 3 * nnz]));
                            vz0[z] = vz0[z] + inv_r[z]*dt_inv_dz * (coefz[0] * (pr0[z + 1] - pr0[z]) +
                                                          coefz[1] * (pr0[z + 2] - pr0[z - 1]) +
                                                          coefz[2] * (pr0[z + 3] - pr0[z - 2]) +
                                                          coefz[3] * (pr0[z + 4] - pr0[z - 3]));
//                            MSG("vz0[z]=%f\n,",vz0[z]);
                        }
                    }
                }
            }
        }
    }

    // Pressure update loop
    #pragma omp parallel for collapse(3) schedule(dynamic) private(xmin, xmax, zmin, zmax, ymin, ymax, pr0, vx0, vy0, vz0, rx)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
                const Myint xmax = fmin(dimx, xmin + BLOCKX);
                const Myint ymax = fmin(dimy, ymin + BLOCKY);
                const Myint zmax = fmin(dimz, zmin + BLOCKZ);

                for (int x = xmin; x < xmax; x++) {
                    for (int y = ymin; y < ymax; y++) {
                        pr0 = &(u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vx0 = &(vx[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vy0 = &(vy[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        vz0 = &(vz[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);
                        rx = &(roc2[1ULL * x * dimy * dimz + y * dimz]);
                        #pragma omp simd
                        for (int z = zmin; z < zmax; z++) {
                            pr0[z] = pr0[z] + rx[z] * (coefx[0] * inv_dx * (vx0[z] - vx0[z - 1 * nnyz]) +
                                                      coefy[0] * inv_dy * (vy0[z] - vy0[z - 1 * nnz]) +
                                                      coefz[0] * inv_dz * (vz0[z] - vz0[z - 1]) +
                                                      coefx[1] * inv_dx * (vx0[z + 1 * nnyz] - vx0[z - 2 * nnyz]) +
                                                      coefy[1] * inv_dy * (vy0[z + 1 * nnz] - vy0[z - 2 * nnz]) +
                                                      coefz[1] * inv_dz * (vz0[z + 1] - vz0[z - 2]) +
                                                      coefx[2] * inv_dx * (vx0[z + 2 * nnyz] - vx0[z - 3 * nnyz]) +
                                                      coefy[2] * inv_dy * (vy0[z + 2 * nnz] - vy0[z - 3 * nnz]) +
                                                      coefz[2] * inv_dz * (vz0[z + 2] - vz0[z - 3]) +
                                                      coefx[3] * inv_dx * (vx0[z + 3 * nnyz] - vx0[z - 4 * nnyz]) +
                                                      coefy[3] * inv_dy * (vy0[z + 3 * nnz] - vy0[z - 4 * nnz]) +
                                                      coefz[3] * inv_dz * (vz0[z + 3] - vz0[z - 4]));
                            pr0[z]*=dampx[x + sx] * dampy[y + sy] * dampz[z + sz];
//                            MSG("pr0[z]=%f\n,",pr0[z]);
                        }
                    }
                }
            }
        }
    }
}

void wave_update_source(sismap_t *s, shot_t *shot,float *u0,float sterm) {
//    MSG("shot->srcidx=%d\n",shot->srcidx);
    const int nnz = s->dimz + 2 * s->sz;
    const int nnx = s->dimx + 2 * s->sx;
    const int nny = s->dimy + 2 * s->sy;
    const int nnxy = nnx * nny;
    const int nnyz = nny * nnz;
//    u0[(s->src_depth + s->sz) * (2 * s->sx + s->dimx)*(2 * s->sy + s->dimy)+shot->srcidx] += sterm; // simwave:zyx
    u0[ shot->srcidx*nnz + (s->src_depth + s->sz) ] += sterm; // stencil:xyz
}

void wave_extract_sismos(sismap_t *s, float *u1, unsigned int t, float *sismos) {
	const int sz = s->sz;
	const int nnz=s->dimz + 2*sz;
    for (unsigned int ir = 0; ir < s->rcv_len; ir++) {
        // from stencil-dev
//        sismos[s->rcv_len * t + ir]=u1[ (rcvx+s->sx)*nnyz+(rcvy+s->sy)*nnz+(s->rcv_depth+s->sz) ];
        sismos[s->rcv_len * t + ir]=u1[ 1ULL*s->rcv[ir] * nnz + (s->rcv_depth+s->sz) ];
    }
}

void wave_inject_sismos(sismap_t *s, float *u0, unsigned int t, float *sismos) {
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const int sx = s->sx;
	const int sy = s->sy;
	const int sz = s->sz;
	const int nnx = dimx + 2 * sx;
	const int nny = dimy + 2 * sy;
	const int nnz = dimz + 2 * sz;
	const long int nnxy=(long int)nnx * nny; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)nny * nnz;

    for (unsigned int ir = 0; ir < s->rcv_len; ir++) {
    	//// u0[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]
    	//// s->rcv[ir]=(s->sx+x)*(s->dimy+2*s->sy)+(y+s->sy);
    	u0[ 1ULL*s->rcv[ir] * nnz + (s->rcv_depth+s->sz) ]+= sismos[s->rcv_len*t + ir];
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
    const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
    if (t % s->nb_snap == 0) {
        const int nnx = s->dimx + 2 * s->sx;
        const int nny = s->dimy + 2 * s->sy;
        const int nnxy = nnx * nny;

        float *restrict ux;
        float *restrict wx;
        float *restrict imgx;
        float *restrict ilmx;

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
                        for (int x = 0; x < s->img_dimx;x++) {
//                            if (ux[x]>0.0001)
//                                {MSG("ux[x]=%f,wx[x]=%f,I=%f",ux[x],wx[x],ux[x]*wx[x]);}
                            imgx[x] += ux[x]*wx[x];
                            ilmx[x] += wx[x]*wx[x];
//              IMG_S(z, y, x)+=  U1(z,y+s->pmly,x+s->pmlx)*FWD(z,y+s->pmly,x+s->pmlx);
//              ILM_S(z, y, x)+= FWD(z,y+s->pmly,x+s->pmlx)*FWD(z,y+s->pmly,x+s->pmlx);
                        }
                    }
                }
            }
        }
    }
}

void wave_image_condition_block_xyz(sismap_t *s, float *u1,
                                float *fwd,
                                float *img_shot, float *ilm_shot,
                                unsigned int t) {
    const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;

	const int img_dimx = s->img_dimx;
	const int img_dimy = s->img_dimy;
	const int img_dimz = s->img_dimz;

	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;

	const int nnx = s->dimx + 2 * s->sx;
	const int nny = s->dimy + 2 * s->sy;
	const int nnz = s->dimz + 2 * s->sz;
	const int sx = s->sx;
	const int sy = s->sy;
	const int sz = s->sz;

	const long int nnxy=(long int)nnx * nny; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)nny * nnz;

	unsigned int z, y, x;
	unsigned int xmin, xmax, zmin, zmax, ymin, ymax;

    if (t % s->nb_snap == 0) {
        float *restrict ux;
        float *restrict wx;
        float *restrict imgx;
        float *restrict ilmx;

#pragma omp parallel for collapse(3) schedule(dynamic) private(ux,imgx,ilmx,wx)
		for (xmin = 0; xmin < img_dimx; xmin += BLOCKX) {
			for (ymin = 0; ymin < img_dimy; ymin += BLOCKY) {
				for (zmin = 0; zmin < img_dimz; zmin += BLOCKZ) {
					const Myint xmax = fmin(img_dimx, xmin + BLOCKX);
					const Myint ymax = fmin(img_dimy, ymin + BLOCKY);
					const Myint zmax = fmin(img_dimz, zmin + BLOCKZ);
					for (int x = xmin; x < xmax; x++) {
		//                	MSG("vz0[z]=%f\n",vz0[z]);
						for (int y = ymin; y < ymax; y++) {
							ux = &(u1[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);	// backward. receiver wavefield.
							wx = &(fwd[1ULL * (x + sx) * nnyz + (y + sy) * nnz + sz]);	// forward. source wavefield.
							imgx = &(img_shot[1ULL * x * dimy * dimz + y * dimz]);
							ilmx = &(ilm_shot[1ULL * x * dimy * dimz + y * dimz]);
							#pragma omp simd
							for (int z = zmin; z < zmax; z++) {
								imgx[z] += ux[z]*wx[z];
								ilmx[z] += wx[z]*wx[z]; //source wavefield illumination
		//                            MSG("vz0[z]=%f\n,",vz0[z]);
							}
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
    // printf("Stubbed wave_save_img called\n");
    // return;
}

//int wave_check_fields(float *tab, size_t len) {
//    int not_valid = 0;
//    for (size_t i = 0; i < len; ++i) {
//        if (isnan(tab[i])) not_valid++;
//    }
//    return not_valid;
//}

void wave_min_max(char *str, float *tab, size_t len) {
    float max = 0.0;
    float min = 99999999999999999999.0;
    for (unsigned int i = 0; i < len; i++) {
        if (tab[i] < min) min = tab[i];
        if (tab[i] > max) max = tab[i];
    }
    printf("... %s min max %f %f\n", str, min, max);
}
