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
/// @file src/velocity.c
/// @brief This file contains the implementation of the velocity generator.
///
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/velocity.h>
#include <stencil/interp.h>

#define V(z,x) vtab[s->dimx*(z) + (x+s->pmlx)]
#define TMP(z, x) tmp[(s->vel_dimx+1)*(z) + x]

#define DOWNLOADSALT3D "\
#/bin/bash \n\
pwd \n\
mkdir ./velocity_models \n\
wget -q --content-disposition 'https://www.dropbox.com/scl/fi/05bc42ctyyhx1d7d10k51/salt3d_676x676x201_xyz.raw?rlkey=04646f0r2vil9ph5m9yjsiras&dl=0' \n\
rm ./velocity_models/* \n\
mv salt3d_676x676x201_xyz.raw ./velocity_models/salt3d_676x676x201_xyz.raw \n\
"

#define __DUMP_VEL

void velocity_generate_model(sismap_t *s, float *vtab, unsigned int layers) {
//    MSG("__________________________Inside function velocity_generate_model__________________________");
    unsigned int x, y, z, i;
    float delta = layers == 1 ? (s->vmax - s->vmin) :
                  (s->vmax - s->vmin) / (layers - 1);
    for (i = 0; i < layers; ++i) {
        for (z = (i * s->dimz) / layers; z < MIN(((i + 1) * s->dimz) / layers, s->dimz); z++) {
            for (y = 0; y < s->dimy; y++) {
                for (x = 0; x < s->dimx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] = s->vmin + (i * delta);

                }
            }
        }
    }

    for (z = 0; z < s->dimz; z++)
        for (y = 0; y < s->dimy; y++)
            for (x = 0; x < s->dimx; x++)
                if (s->order==1) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            pow(s->dt, 1) * pow(vtab[s->dimx * (s->dimy * z + y) + x], 2);}
                else{
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            pow(s->dt, 2) * pow(vtab[s->dimx * (s->dimy * z + y) + x], 2);}


//    MSG("vtab=%f\n", vtab[100]);
//    MSG("vtab=%f\n", vtab[1000]);
//    MSG("vtab=%f\n", vtab[5000]);
//    MSG("vtab=%f\n", vtab[10000]);
//    MSG("vtab=%f\n", vtab[20000]);
//    MSG("vtab=%f\n", vtab[30000]);

#ifdef __DUMP_VEL
    char tmp[512];
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    CHK(system(tmp), "failed to create output directory for snapshots");
    sprintf(tmp, "%s/augmented_vel.raw", OUTDIR);
    FILE *fd = fopen(tmp, "wb");
    CHK(fwrite(vtab, s->size_eff * sizeof(float), 1, fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
#endif //__DUMP_VEL
}

void velocity_query_model(sismap_t *s) {
//    MSG("__________________________Inside function velocity_query_model__________________________");
    if (strcmp("NONE", s->vel_file) == 0) {
        s->vmin = 1500.0;
        s->vmax = 4500.0;
    } else {
        size_t vel_size = 1LL * s->vel_dimx * s->vel_dimy * s->vel_dimz;
        float f, *tmp = (float *) malloc(vel_size * sizeof(float));

        FILE *fd = fopen(s->vel_file, "rb");
        MSG("velocity file name : %s\n", s->vel_file);
        MSG("velocity file dim : %d %d %d\n", s->vel_dimx, s->vel_dimy, s->vel_dimz);
        CHK(fd == NULL, "failed to open the velocity file");
        CHK(fread(tmp, sizeof(float), vel_size, fd) != vel_size,
            "failed to read from the velocity file");
        fclose(fd);
        /// search for vmin and vmax.
        s->vmin = FLT_MAX;
        s->vmax = FLT_MIN;
        for (unsigned int z = 0; z < s->vel_dimz; z++) {
            for (unsigned int y = 0; y < s->vel_dimy; y++) {
                for (unsigned int x = 0; x < s->vel_dimx; x++) {
                    f = tmp[1ULL * s->vel_dimx * (s->vel_dimy * z + y) + x];
                    if (f > s->vmax) s->vmax = f;
                    if (f < s->vmin) s->vmin = f;
                }
            }
        }
        free(tmp);
    }
}

void fill_coef_matrix(sismap_t *s,float *vtab,float *dens) {
    unsigned int x, y, z;

	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const long int nnxy=(long int)dimx*dimy; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)dimy*dimz;
	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;
	/// fill vtab matrix with rho*v^2*dt or rho*v^2*dt2 value

	MSG("entering fill_coef_matrix\n");
    if (s->order==1) {
        MSG("order=1,vtab\n");
		for (x = 0; x < dimx; x++) {
			for (y = 0; y < dimy; y++) {
				for (z = 0; z < dimz; z++) {
					unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
					if (idx >= s->size_eff) {
						fprintf(stderr, "Index out of bounds: idx=%llu, size_eff=%llu\n",idx,s->size_eff);
						exit(EXIT_FAILURE);
					}
					vtab[idx]=s->dt*dens[idx]*pow(vtab[idx],2);
				}
			}
		}

//#pragma omp parallel for collapse(3)
//		for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
//			for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
//				for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
//					const Myint xmax = fmin(dimx, xmin + BLOCKX);
//					const Myint ymax = fmin(dimy, ymin + BLOCKY);
//					const Myint zmax = fmin(dimz, zmin + BLOCKZ);
//					// loop on grid points inside block
//					for (int x = xmin; x < xmax; x++) {
//						for (int y = ymin; y < ymax; y++) {
//	#pragma omp simd
//							for (int z = zmin; z < zmax; z++) {
//								vtab[1ULL*(x)*nnyz+(y)*dimz+z]=
//										s->dt*dens[(x)*nnyz+(y)*dimz+z]*pow(vtab[(x)*nnyz+(y)*dimz+z],2);
//							}
//						}
//					}
//				}
//			}
//		}
    }else
    {
        MSG("order=2,vtab\n");
        for (x = 0; x < dimx; x++) {
			for (y = 0; y < dimy; y++) {
				for (z = 0; z < dimz; z++) {
					unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
					if (idx >= s->size_eff) {
						fprintf(stderr, "Index out of bounds: idx=%llu, size_eff=%llu\n",
								idx, s->size_eff);
						exit(EXIT_FAILURE);
					}
					vtab[idx]=pow(s->dt,2)*dens[idx]*pow(vtab[idx],2);
				}
			}
		}

//#pragma omp parallel for collapse(3)
//		for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
//			for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
//				for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
//					const Myint xmax = fmin(dimx, xmin + BLOCKX);
//					const Myint ymax = fmin(dimy, ymin + BLOCKY);
//					const Myint zmax = fmin(dimz, zmin + BLOCKZ);
//					// loop on grid points inside block
//					for (int x = xmin; x < xmax; x++) {
//						for (int y = ymin; y < ymax; y++) {
//	#pragma omp simd
//							for (int z = zmin; z < zmax; z++) {
//								vtab[1ULL*(x)*nnyz+(y)*dimz+z]=
//										pow(s->dt,2)*dens[(x)*nnyz+(y)*dimz+z]*pow(vtab[(x)*nnyz+(y)*dimz+z],2);
//							}
//						}
//					}
//				}
//			}
//		}
    }
}

void dump_vel(sismap_t *s,float *vtab,float *dens) {
#ifdef __DUMP_VEL
    MSG("... dump velocity ...");
    char tmp[512];
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    CHK(system(tmp), "failed to create output directory for snapshots");
    sprintf(tmp, "%s/velocity.raw", OUTDIR);
    FILE *fd = fopen(tmp, "wb");
    CHK(fwrite(vtab, s->size_eff * sizeof(float),1,fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
    ////////////////////////////////
    MSG("... dump density ...");
	char tmp2[512];
	sprintf(tmp2, "mkdir -p %s", OUTDIR);
	CHK(system(tmp2), "failed to create output directory for snapshots");
	sprintf(tmp2,"%s/density.raw", OUTDIR);
	FILE *fd2 = fopen(tmp2, "wb");
	CHK(fwrite(dens, s->size_eff * sizeof(float),1,fd2) != 1,
		"failed to write the density file");
	fclose(fd2);
#endif //__DUMP_VEL
}

void dump_coef(sismap_t *s, float *vtab) {
    MSG("... dump coefficient matrix ...");
    char tmp[512];
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    CHK(system(tmp), "failed to create output directory for snapshots");
    if (s->order==1) {
    	sprintf(tmp, "%s/coef_1st.raw", OUTDIR);
    }
    else{
    	sprintf(tmp, "%s/coef_2nd.raw", OUTDIR);
    }
    FILE *fd = fopen(tmp, "wb");
    CHK(fwrite(vtab, s->size_eff * sizeof(float),1,fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
}

void velocity_load_model_2d(sismap_t *s, float *vtab) {
//    MSG("__________________________Inside function velocity_load_model_2d__________________________");
    FILE *fd;
    float *vtmp, *tmp;
    if (strcmp("NONE", s->vel_file) == 0) {
        velocity_generate_model(s, vtab, 4);
    } else {
        size_t vel_size = s->vel_dimx * s->vel_dimz;
        vtmp = (float *) malloc(vel_size * sizeof(float));
        tmp = (float *) malloc((s->vel_dimx + 1) * (s->vel_dimz + 1) * sizeof(float));
        fd = fopen(s->vel_file, "rb");
        CHK(fd == NULL, "failed to open the velocity file");
        CHK(fread(vtmp, vel_size * sizeof(float), 1, fd) != 1,
            "failed to read from the velocity file");
        fclose(fd);

        for (unsigned int z = 0; z < s->vel_dimz; z++) {
            for (unsigned int x = 0; x < s->vel_dimx; x++) {
                tmp[(s->vel_dimx + 1) * z + x] =
                        pow(s->dt, 2) * pow(vtmp[s->vel_dimx * z + x], 2);
            }
        }
        for (unsigned int x = 0; x < s->vel_dimx; x++) {
            tmp[(s->vel_dimx + 1) * (s->vel_dimz) + x] =
                    pow(s->dt, 2) * pow(vtmp[s->vel_dimx * (s->vel_dimz - 1) + x], 2);
        }
        for (unsigned int z = 0; z < s->vel_dimz; z++) {
            tmp[(s->vel_dimx + 1) * z + (s->vel_dimx)] =
                    pow(s->dt, 2) * pow(vtmp[s->vel_dimx * z + (s->vel_dimx - 1)], 2);
        }
        free(vtmp);

        /// load the velocity in vtab after interpolation.
        for (unsigned int z = 0; z < s->vel_dimz; z++) {
            for (unsigned int x = 0; x < s->vel_dimx; x++) {
                vtab[s->dimx * (z * s->dtrpz) + ((x * s->dtrpx) + s->pmlx)] =
                        tmp[(s->vel_dimx + 1) * z + x];
                for (unsigned int iz = 0; iz < s->dtrpz; iz++) {
                    for (unsigned int ix = 0; ix < s->dtrpx; ix++) {
                        unsigned int zidx = (z * s->dtrpz) + iz;
                        unsigned int xidx = (x * s->dtrpx) + ix;
                        if ((zidx < s->dimz) && (xidx < s->dimx)) {
                            float c00 = TMP(z, x);
                            float c01 = TMP(z, x + 1);
                            float c10 = TMP(z + 1, x);
                            float c11 = TMP(z + 1, x + 1);
                            float tx = (float) ix / (float) s->dtrpx;
                            float tz = (float) iz / (float) s->dtrpz;
                            V(zidx, xidx) = bilinearinterp(c00, c01, c10, c11, tx, tz);
                        }
                    }
                }
            }
        }

        /// extend to the Z PML regions.
        for (unsigned int z = s->dimz - s->pmlz; z < s->dimz; z++) {
            for (unsigned int y = 0; y < s->dimy; y++) {
                for (unsigned int x = 0; x < s->dimx; x++) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            vtab[s->dimx * (s->dimy * (s->dimz - s->pmlz - 1) + y) + x];
                }
            }
        }
        /// extend to the Y PML regions.
        for (unsigned int z = 0; z < s->dimz - s->pmlz; z++) {
            for (unsigned int y = 0; y < s->pmly; y++) {
                for (unsigned int x = s->pmlx; x < s->dimx - s->pmlx; x++) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            vtab[s->dimx * (s->dimy * z + s->pmly) + x];
                }
            }
            for (unsigned int y = s->dimy - s->pmly; y < s->dimy; y++) {
                for (unsigned int x = s->pmlx; x < s->dimx - s->pmlx; x++) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            vtab[s->dimx * (s->dimy * z + (s->dimy - s->pmly - 1)) + x];
                }
            }
        }
        /// extend to the X PML regions.
        for (unsigned int z = 0; z < s->dimz - s->pmlz; z++) {
            for (unsigned int y = s->pmly; y < s->dimy - s->pmly; y++) {
                for (unsigned int x = 0; x < s->pmlx; x++) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            vtab[s->dimx * (s->dimy * (z) + y) + s->pmlx];
                }
                for (unsigned int x = s->dimx - s->pmlx; x < s->dimx; x++) {
                    vtab[s->dimx * (s->dimy * z + y) + x] =
                            vtab[s->dimx * (s->dimy * (z) + y) + s->dimx - s->pmlx - 1];
                }
            }
        }
        free(tmp);
    }
#ifdef __DUMP_VEL
    char stmp[512];
    sprintf(stmp, "mkdir -p %s", OUTDIR);
    CHK(system(stmp), "failed to create output directory for snapshots");
    sprintf(stmp, "%s/augmented_vel.raw", OUTDIR);
    fd = fopen(stmp, "wb");
    CHK(fwrite(vtab, s->size_eff * sizeof(float), 1, fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
#endif //__DUMP_VEL
}

#undef V
#undef TMP
#define V(z, y, x) vtab[s->dimx*(s->dimy*(z) + (y+s->pmly)) + (x+s->pmlx)]
#define TMP(z, y, x) tmp[s->vel_dimx*(s->vel_dimy*(z) + y) + x]

void velocity_load_model_3d(sismap_t *s, float *vtab) {
//    MSG("__________________________Inside function velocity_load_model_3d__________________________");
    FILE *fd;
    if (strcmp("NONE", s->vel_file) == 0) {
        velocity_generate_model(s, vtab, 4);
    } else {
        size_t vel_size = 1LL * s->vel_dimx * s->vel_dimy * s->vel_dimz;
        float *tmp = (float *) malloc(vel_size * sizeof(float));

        fd = fopen(s->vel_file, "rb");
        CHK(fd == NULL, "failed to open the velocity file");
        CHK(fread(tmp, vel_size * sizeof(float), 1, fd) != 1,
            "failed to read from the velocity file");
        fclose(fd);
        //#define __TRANSP_VEL
#ifdef __TRANSP_VEL
        char ktmp[512];
        sprintf(ktmp, "mkdir -p %s", OUTDIR);
        CHK(system(ktmp), "failed to create output directory for snapshots");
        sprintf(ktmp, "%s/transp_vel.raw", OUTDIR);
        fd  = fopen(ktmp, "wb");
        // CHK(fwrite(tab, s->size_eff*sizeof(float), 1, fd)!=1,
        //     "failed to write the velocity file");
        unsigned int vx = s->vel_dimx;
        unsigned int vy = s->vel_dimy;
        unsigned int vz = s->vel_dimz;
        printf("vx vy vz %u %u %u\n", vx, vy, vz);
        for (unsigned int x=0; x < vx ; x++)
          for (unsigned int y=0; y < vy ; y++)
            for (unsigned int z=0; z < vz ; z++)
              fwrite(&tmp[1ULL*s->vel_dimx*(s->vel_dimy*z + y) + x], sizeof(float), 1, fd);
        fclose(fd);
#endif // __TRANSP_VEL
        for (unsigned int z = 0; z < s->vel_dimz; z++) {
            for (unsigned int y = 0; y < s->vel_dimy; y++) {
                for (unsigned int x = 0; x < s->vel_dimx; x++) {
                    tmp[1ULL * s->vel_dimx * (s->vel_dimy * z + y) + x] =
                            pow(s->dt, 2) * pow(tmp[1ULL * s->vel_dimx * (s->vel_dimy * z + y) + x], 2);
                }
            }
        }

        /// load the velocity in vtab after interpolation (not done yet).
        for (unsigned int z = 0; z < s->vel_dimz - 1; z++) {
            for (unsigned int y = 0; y < s->vel_dimy - 1; y++) {
                for (unsigned int x = 0; x < s->vel_dimx - 1; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * (z * s->dtrpz) +
                                           ((y * s->dtrpy) + s->pmly)) +
                         ((x * s->dtrpx) + s->pmlx)] =
                            tmp[1ULL * s->vel_dimx * (s->vel_dimy * z + y) + x];
                    for (unsigned int iz = 0; iz < s->dtrpz; iz++) {
                        for (unsigned int iy = 0; iy < s->dtrpy; iy++) {
                            for (unsigned int ix = 0; ix < s->dtrpx; ix++) {
                                float c000 = TMP(z, y, x);
                                float c001 = TMP(z, y, x + 1);
                                float c010 = TMP(z, y + 1, x);
                                float c011 = TMP(z, y + 1, x + 1);
                                float c100 = TMP(z + 1, y, x);
                                float c101 = TMP(z + 1, y, x + 1);
                                float c110 = TMP(z + 1, y + 1, x);
                                float c111 = TMP(z + 1, y + 1, x + 1);
                                float tx = (float) ix / (float) s->dtrpx;
                                float ty = (float) iy / (float) s->dtrpy;
                                float tz = (float) iz / (float) s->dtrpz;
                                V((z * s->dtrpz) + iz,
                                  (y * s->dtrpy) + iy,
                                  (x * s->dtrpx) + ix) = trilinearinterp(c000, c001,
                                                                         c010, c011,
                                                                         c100, c101,
                                                                         c110, c111, tx, ty, tz);
                            }
                        }
                    }
                }
            }
        }

        /// extend to the Z PML regions.
        for (unsigned int z = s->dimz - s->pmlz; z < s->dimz; z++) {
            for (unsigned int y = 0; y < s->dimy; y++) {
                for (unsigned int x = 0; x < s->dimx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] =
                            vtab[1ULL * s->dimx * (s->dimy * (s->dimz - s->pmlz - 1) + y) + x];
                }
            }
        }
        /// extend to the Y PML regions.
        for (unsigned int z = 0; z < s->dimz - s->pmlz; z++) {
            for (unsigned int y = 0; y < s->pmly; y++) {
                for (unsigned int x = s->pmlx; x < s->dimx - s->pmlx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] =
                            vtab[1ULL * s->dimx * (s->dimy * z + s->pmly) + x];
                }
            }
            for (unsigned int y = s->dimy - s->pmly; y < s->dimy; y++) {
                for (unsigned int x = s->pmlx; x < s->dimx - s->pmlx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] =
                            vtab[1ULL * s->dimx * (s->dimy * z + (s->dimy - s->pmly - 1)) + x];
                }
            }
        }
        /// extend to the X PML regions.
        for (unsigned int z = 0; z < s->dimz - s->pmlz; z++) {
            for (unsigned int y = s->pmly; y < s->dimy - s->pmly; y++) {
                for (unsigned int x = 0; x < s->pmlx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] =
                            vtab[1ULL * s->dimx * (s->dimy * (z) + y) + s->pmlx];
                }
                for (unsigned int x = s->dimx - s->pmlx; x < s->dimx; x++) {
                    vtab[1ULL * s->dimx * (s->dimy * z + y) + x] =
                            vtab[1ULL * s->dimx * (s->dimy * (z) + y) + s->dimx - s->pmlx - 1];
                }
            }
        }
        free(tmp);
    }
#ifdef __DUMP_VEL
    char vtmp[512];
    sprintf(vtmp, "mkdir -p %s", OUTDIR);
    CHK(system(vtmp), "failed to create output directory for snapshots");
    sprintf(vtmp, "%s/augmented_vel.raw", OUTDIR);
    fd = fopen(vtmp, "wb");
    CHK(fwrite(vtab, 1ULL * s->size_eff * sizeof(float), 1, fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
#endif //__DUMP_VEL
}

void velocity_const_model2(sismap_t *s, float *vtab) {
//	MSG("... ! in velocity_const_model2! ...");
	unsigned int x, y, z, i;

	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const long int nnxy=(long int)dimx*dimy; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)dimy*dimz;
	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;

//	MSG("... ! before the loop! ...");
//	MSG("dimx=%d, dimy=%d, dimz=%d\n",dimx,dimy,dimz);
//	MSG("s->size_eff=%d\n",s->size_eff);
//#pragma omp parallel for collapse(3)
//    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
//    	MSG("xmin=%d",xmin);
//        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
//            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
//            	const Myint xmax = fmin(dimx, xmin + BLOCKX);
//				const Myint ymax = fmin(dimy, ymin + BLOCKY);
//				const Myint zmax = fmin(dimz, zmin + BLOCKZ);
////				MSG("xmax=%d,ymax=%d,zmax=%d",xmax,ymax,zmax);
//				// loop on grid points inside block
//				for (int x = xmin; x < xmax; x++) {
//					for (int y = ymin; y < ymax; y++) {
//#pragma omp simd
//						for (int z = zmin; z < zmax; z++) {
////							MSG("x=%d,y=%d,z=%d",x,y,z);
//							unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
//							vtab[idx]=1500;
//						}
//					}
//				}
//			}
//		}
//    }
    ////////////////////////////////
	for (x = 0; x < dimx; x++) {
		for (y = 0; y < dimy; y++) {
			for (z = 0; z < dimz; z++) {
				unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
				if (idx >= s->size_eff) {
					fprintf(stderr, "Index out of bounds: idx=%llu, size_eff=%llu\n",
							idx, s->size_eff);
					exit(EXIT_FAILURE);
				}
				vtab[idx] = 1500.0f;
			}
		}
	}
    ////////////////////////////////
    MSG("Exiting velocity_const_model2\n");
}

void velocity_2layer_model(sismap_t *s, float *vtab, unsigned int layers) {
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const long int nnxy=(long int)dimx*dimy; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)dimy*dimz;
    unsigned int x, y, z, i,val;
    for (z=0; z < s->dimz; z++) {
        if (z<120){
            val=1500;
        }else
        {
            val=3000;
        }
        for (y = 0; y < s->dimy; y++) {
            for (x = 0; x < s->dimx; x++) {
            	unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
                vtab[idx]=val;
            }
        }
    }

#ifdef __DUMP_VEL
    char tmp[512];
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    CHK(system(tmp), "failed to create output directory for snapshots");
    sprintf(tmp, "%s/augmented_vel.raw", OUTDIR);
    FILE *fd = fopen(tmp, "wb");
    CHK(fwrite(vtab, s->size_eff * sizeof(float), 1, fd) != 1,
        "failed to write the velocity file");
    fclose(fd);
#endif //__DUMP_VEL
    /////////////////////////////////////////////////////
}

void velocity_load_model(sismap_t *s, float *vtab) {
//    MSG("__________________________Inside function velocity_load_model__________________________");
    if (s->dim2) velocity_load_model_2d(s, vtab);
    else velocity_load_model_3d(s, vtab);
}

void velocity_load_salt3d(sismap_t *s, float *vtab) {
    /////////////////////////////////////////////////////
	// this salt3d model distribution is recorded in 676x676x201 grid.
	// this grid will be cropped according to specified (s->dimx ...)sizes
	/////////////////////////////////////////////////////
	// Load original model
	/////////////////////////////////////////////////////
    const char *file_namev = "./velocity_models/salt3d_676x676x201_xyz.raw";
    FILE *fd;
    MSG("velocity_load_salt3d");
    const int dimx_orig = 676;
	const int dimy_orig = 676;
	const int dimz_orig = 201;
    size_t vel_size = 1LL * dimx_orig * dimy_orig * dimz_orig;
    float *tmp=(float *) malloc(vel_size * sizeof(float));
    fd = fopen(file_namev, "rb");
    if (fd == NULL){
        MSG("reading from velocity file %s",file_namev);
        MSG("downloading velocity file from dropbox");
        system(DOWNLOADSALT3D);
        fd = fopen(file_namev,"rb");
    }
    CHK(fd == NULL,"failed to open the velocity file");
    CHK(fread(tmp,vel_size * sizeof(float), 1, fd) != 1,"failed to read from the velocity file");
    fclose(fd);
    /////////////////////////////////////////////////////
	// Crop it
	/////////////////////////////////////////////////////
	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const long int nnxy=(long int)dimx*dimy; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)dimy*dimz;
	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;
#pragma omp parallel for collapse(3)
    for (xmin = 0; xmin < dimx; xmin += BLOCKX) {
        for (ymin = 0; ymin < dimy; ymin += BLOCKY) {
            for (zmin = 0; zmin < dimz; zmin += BLOCKZ) {
            	const Myint xmax = fmin(dimx, xmin + BLOCKX);
				const Myint ymax = fmin(dimy, ymin + BLOCKY);
				const Myint zmax = fmin(dimz, zmin + BLOCKZ);
				// loop on grid points inside block
				for (int x = xmin; x < xmax; x++) {
					for (int y = ymin; y < ymax; y++) {
#pragma omp simd
						for (int z = zmin; z < zmax; z++) {
							// if velocity file is recorded in ZYX order
							vtab[ 1ULL*(x)*nnyz+(y)*dimz+z ]=tmp[ 1ULL*dimx_orig*(dimy_orig*(z)+(y))+(x) ];
							/////// if velocity file is recorded in XYZ order
//							vtab[1ULL*(x)*nnyz+(y)*dimz+z]=tmp[1ULL*(x)*nnyz+(y)*dimz+z];
						}
					}
				}
			}
		}
    }
    /////////////////////////////////
    MSG("custom velocity model read in memory");
}

////////////
void density_const_model(sismap_t *s,float *dens,float *inv_rho) {
	unsigned int x, y, z, i;

	const int BLOCKX=s->blockx;
	const int BLOCKY=s->blocky;
	const int BLOCKZ=s->blockz;
	const int dimx = s->dimx;
	const int dimy = s->dimy;
	const int dimz = s->dimz;
	const long int nnxy=(long int)dimx*dimy; // XYZ order: x-slowest, z-fastest
	const long int nnyz=(long int)dimy*dimz;
	unsigned int xmin,xmax,zmin,zmax,ymin,ymax;
	for (x = 0; x < dimx; x++) {
		for (y = 0; y < dimy; y++) {
			for (z = 0; z < dimz; z++) {
				unsigned long long idx = 1ULL * x * nnyz + 1ULL * y * dimz + z;
				if (idx >= s->size_eff) {
					fprintf(stderr, "Index out of bounds: idx=%llu, size_eff=%llu\n",
							idx, s->size_eff);
					exit(EXIT_FAILURE);
				}
				dens[idx]=2000.0f;		// 2000
				inv_rho[idx]=1/dens[idx];
			}
		}
	}
}

////////////


