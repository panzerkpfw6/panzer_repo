/// Function 'gather' stacks rtm images for different shots.
#include <stdio.h>
#include <unistd.h>  // For getcwd
#include <math.h>
#include <string.h>
#include <dirent.h>
//#include <cuda_runtime.h>
#include <stencil/parser.h>
#include <stencil/stencil.h>



void gather_img_div_ilm(unsigned int len,
                        float *img_shot, float *ilm_shot, float *img) {
	#pragma omp parallel for
  for (unsigned int i = 0; i < len; i++) {
    img[i] += img_shot[i]/ilm_shot[i];
  }
}

void gather_img_div_ilm_smart(unsigned int len,
                        float *img_shot, float *ilm_shot, float *img) {
    #pragma omp parallel for
    for (unsigned int i = 0; i < len; i++) {
        if (fabs(ilm_shot[i]) > 1e-6) {  // Avoid division by zero or near-zero
            img[i] += img_shot[i] / ilm_shot[i];
        } else {
            img[i] += 0.0;  // Skip contribution or use a small constant
        }
    }
}

void gather_img_ilm(unsigned int len,
                    float *img_shot, float *ilm_shot, float *img, float *ilm) {
  #pragma omp parallel for
  for (unsigned int i = 0; i < len; i++) {
    img[i] += img_shot[i];
    ilm[i] += ilm_shot[i];
  }
}

// New helper functions for Laplacian filtering
void apply_laplacian_filter(sismap_t *s, float *restrict img) {
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    float *temp = malloc(dimx * dimy * dimz * sizeof(float));
    for (int iter = 0; iter < 2; iter++) { // Apply filter twice for better suppression
        memcpy(temp, img, dimx * dimy * dimz * sizeof(float));
        #pragma omp parallel for collapse(3)
        for (int x = 1; x < dimx - 1; x++) {
            for (int y = 1; y < dimy - 1; y++) {
                for (int z = 4; z < dimz - 1; z++) {
                    long int idx = x * dimy * dimz + y * dimz + z;
                    img[idx] = 0.5f * (temp[idx + dimy * dimz] + temp[idx - dimy * dimz] +
                                       temp[idx + dimz] + temp[idx - dimz] +
                                       temp[idx + 1] + temp[idx - 1] -
                                       6 * temp[idx]);
                }
            }
        }
    }
    free(temp);
}

void apply_laplacian_filter2_small(sismap_t *s, float *restrict img) {
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    float *temp = malloc(dimx * dimy * dimz * sizeof(float));
    for (int iter = 0; iter < 3; iter++) { // Increased iterations for better suppression
        memcpy(temp, img, dimx * dimy * dimz * sizeof(float));
        #pragma omp parallel for collapse(3)
        for (int x = 1; x < dimx - 1; x++) {
            for (int y = 1; y < dimy - 1; y++) {
                for (int z = 1; z < dimz - 1; z++) {
                    long int idx = x * dimy * dimz + y * dimz + z;
                    float laplacian = 0.0f;
                    // Nearest neighbors (6-point)
                    laplacian += temp[idx + dimy * dimz] + temp[idx - dimy * dimz] +
                                 temp[idx + dimz] + temp[idx - dimz] +
                                 temp[idx + 1] + temp[idx - 1];
                    // Diagonal neighbors (additional 20 points, weighted less)
                    laplacian += 0.25f * (
                        temp[idx + dimy * dimz + dimz] + temp[idx + dimy * dimz - dimz] +
                        temp[idx - dimy * dimz + dimz] + temp[idx - dimy * dimz - dimz] +
                        temp[idx + dimy * dimz + 1] + temp[idx + dimy * dimz - 1] +
                        temp[idx - dimy * dimz + 1] + temp[idx - dimy * dimz - 1] +
                        temp[idx + dimz + 1] + temp[idx + dimz - 1] +
                        temp[idx - dimz + 1] + temp[idx - dimz - 1]
                    );
                    laplacian -= (6.0f + 0.25f * 12.0f) * temp[idx]; // Adjust center weight
                    img[idx] = 0.5f * laplacian; // Scale by 1/2 as per paper
                }
            }
        }
    }
    free(temp);
}

void apply_laplacian_filter2(sismap_t *s, float *restrict img) {
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    float *temp = malloc(dimx * dimy * dimz * sizeof(float));

    for (int iter = 0; iter < 3; iter++) { // Increased iterations for better suppression
        memcpy(temp, img, dimx * dimy * dimz * sizeof(float));
        #pragma omp parallel for collapse(3)
        for (int x = 2; x < dimx - 2; x++) { // Adjust bounds for 5x5x3 kernel
            for (int y = 2; y < dimy - 2; y++) {
                for (int z = 1; z < dimz - 1; z++) {
                    long int idx = x * dimy * dimz + y * dimz + z;
                    float laplacian = 0.0f;

                    // Nearest neighbors (X, Y at ±1, Z at ±1)
                    laplacian += 1.0f * (temp[idx + dimy * dimz] + temp[idx - dimy * dimz] +  // X±1
                                        temp[idx + dimz] + temp[idx - dimz] +                // Y±1
                                        temp[idx + 1] + temp[idx - 1]);                     // Z±1

                    // Farther neighbors in X, Y (at ±2)
                    laplacian += 0.5f * (temp[idx + 2 * dimy * dimz] + temp[idx - 2 * dimy * dimz] +  // X±2
                                        temp[idx + 2 * dimz] + temp[idx - 2 * dimz]);          // Y±2

                    // Diagonal neighbors in X, Y plane (at ±1, ±1; ±2, ±1; ±1, ±2; ±2, ±2)
                    laplacian += 0.25f * (
                        // X±1, Y±1, Z
                        temp[idx + dimy * dimz + dimz] + temp[idx + dimy * dimz - dimz] +
                        temp[idx - dimy * dimz + dimz] + temp[idx - dimy * dimz - dimz] +
                        // X±2, Y±1, Z
                        temp[idx + 2 * dimy * dimz + dimz] + temp[idx + 2 * dimy * dimz - dimz] +
                        temp[idx - 2 * dimy * dimz + dimz] + temp[idx - 2 * dimy * dimz - dimz] +
                        // X±1, Y±2, Z
                        temp[idx + dimy * dimz + 2 * dimz] + temp[idx + dimy * dimz - 2 * dimz] +
                        temp[idx - dimy * dimz + 2 * dimz] + temp[idx - dimy * dimz - 2 * dimz] +
                        // X±2, Y±2, Z
                        temp[idx + 2 * dimy * dimz + 2 * dimz] + temp[idx + 2 * dimy * dimz - 2 * dimz] +
                        temp[idx - 2 * dimy * dimz + 2 * dimz] + temp[idx - 2 * dimy * dimz - 2 * dimz]
                    );

                    // Diagonal neighbors involving Z (X±1, Y±1, Z±1; etc.)
                    laplacian += 0.1f * (
                        // X±1, Y, Z±1
                        temp[idx + dimy * dimz + 1] + temp[idx + dimy * dimz - 1] +
                        temp[idx - dimy * dimz + 1] + temp[idx - dimy * dimz - 1] +
                        // X, Y±1, Z±1
                        temp[idx + dimz + 1] + temp[idx + dimz - 1] +
                        temp[idx - dimz + 1] + temp[idx - dimz - 1] +
                        // X±1, Y±1, Z±1
                        temp[idx + dimy * dimz + dimz + 1] + temp[idx + dimy * dimz + dimz - 1] +
                        temp[idx + dimy * dimz - dimz + 1] + temp[idx + dimy * dimz - dimz - 1] +
                        temp[idx - dimy * dimz + dimz + 1] + temp[idx - dimy * dimz + dimz - 1] +
                        temp[idx - dimy * dimz - dimz + 1] + temp[idx - dimy * dimz - dimz - 1]
                    );

                    // Center weight: sum of all weights (excluding center) with negative sign
                    float center_weight = -(6.0f * 1.0f +        // Nearest neighbors
                                           4.0f * 0.5f +        // X, Y at ±2
                                           12.0f * 0.25f +      // Diagonals in X, Y
                                           12.0f * 0.1f);       // Diagonals involving Z
                    laplacian += center_weight * temp[idx];

                    img[idx] = 0.5f * laplacian; // Scale by 1/2 as per paper
                }
            }
        }
    }
    free(temp);
}

void apply_laplacian_xy_filter(sismap_t *s, float *restrict img) {
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    float *temp = malloc(dimx * dimy * dimz * sizeof(float));
    for (int iter = 0; iter < 3; iter++) { // Apply filter 3 times for better suppression
        memcpy(temp, img, dimx * dimy * dimz * sizeof(float));
        #pragma omp parallel for collapse(3)
        for (int x = 1; x < dimx - 1; x++) {
            for (int y = 1; y < dimy - 1; y++) {
                for (int z = 0; z < dimz; z++) { // No boundary checks in z-direction
                    long int idx = x * dimy * dimz + y * dimz + z;
                    float laplacian = 0.0f;
                    // 2D Laplacian in XY plane (5-point stencil)
                    laplacian += temp[idx + dimy * dimz] + temp[idx - dimy * dimz] + // x-direction
                                 temp[idx + dimz] + temp[idx - dimz];                 // y-direction
                    laplacian -= 4.0f * temp[idx]; // Center weight
                    img[idx] = 0.5f * laplacian; // Scale by 1/2 as per paper
                }
            }
        }
    }
    free(temp);
}

void smooth_illumination(sismap_t *s, float *restrict ilm) {
    const int dimx = s->dimx;
    const int dimy = s->dimy;
    const int dimz = s->dimz;
    float *temp = malloc(dimx * dimy * dimz * sizeof(float));
    memcpy(temp, ilm, dimx * dimy * dimz * sizeof(float));
    #pragma omp parallel for collapse(3)
    for (int x = 1; x < dimx - 1; x++) {
        for (int y = 1; y < dimy - 1; y++) {
            for (int z = 1; z < dimz - 1; z++) {
                long int idx = x * dimy * dimz + y * dimz + z;
                ilm[idx] = (temp[idx + dimy * dimz] + temp[idx - dimy * dimz] +
                            temp[idx + dimz] + temp[idx - dimz] +
                            temp[idx + 1] + temp[idx - 1] +
                            6 * temp[idx]) / 12.0f;
            }
        }
    }
    free(temp);
}

void normalize_image(sismap_t *s, float *restrict img_shot, float *restrict ilm_shot, float *restrict img_norm) {
    const unsigned int len = s->size_img;
    #pragma omp parallel for
    for (unsigned int i = 0; i < len; i++) {
        if (fabs(ilm_shot[i]) > 1e-6) {
            img_norm[i] = img_shot[i] / ilm_shot[i];
        } else {
            img_norm[i] = 0.0f;
        }
    }
}


int main_option1(int argc, char* argv[]) {
    sismap_t *s = (sismap_t*)malloc(sizeof(sismap_t));
    parser *p = parser_create("Reverse Time Migration using STENCIL");
    PARSER_BOOTSTRAP(p);
    parser_parse(p, argc, argv);
    s->verbose    = parser_get_bool(p, "verbose");
    s->cpu        = parser_get_bool(p, "cpu");
    s->time_steps = parser_get_int(p, "iter");
    s->cfl        = parser_get_float(p, "cfl");
    s->fmax       = parser_get_float(p, "fmax");
    s->vel_file   = parser_get_string(p, "in");
    s->vel_dimx   = parser_get_int(p, "n1");
    s->vel_dimy   = parser_get_int(p, "n2");
    s->vel_dimz   = parser_get_int(p, "n3");
    s->dx         = parser_get_int(p, "dx");
    s->dy         = parser_get_int(p, "dy");
    s->dz         = parser_get_int(p, "dz");
    s->dcdp       = parser_get_int(p, "dcdp");
    s->dline      = parser_get_int(p, "dline");
    s->drcv       = parser_get_int(p, "drcv");
    s->dshot      = parser_get_int(p, "dshot");
    s->ddepth     = parser_get_int(p, "ddepth");
    s->device     = parser_get_int(p, "device");
    s->first      = parser_get_int(p, "first");
    s->last       = parser_get_int(p, "last");
    s->src_depth  = parser_get_int(p, "src_depth");
    s->rcv_depth  = parser_get_int(p, "rcv_depth");
    s->modeling   = false;
    s->nb_snap    = parser_get_int(p, "nbsnap");
    char *dir     = parser_get_string(p, "dir");

    float* img;
    float* img_only;
    float* ilm_only;
    float *ilm_shot, *img_shot;

    wave_init_numerics(s);
    wave_init_dimensions(s);

    CREATE_BUFFER(img, s->size_img);
    NULIFY_BUFFER(img, s->size_img);
    CREATE_BUFFER(img_only, s->size_img);
    NULIFY_BUFFER(img_only, s->size_img);
    CREATE_BUFFER(ilm_only, s->size_img);
    NULIFY_BUFFER(ilm_only, s->size_img);
    CREATE_BUFFER(img_shot, s->size_img);
    CREATE_BUFFER(ilm_shot, s->size_img);

    if (s->verbose) {
        MSG(" ");
        MSG("... stencil information:");
        MSG("... compute domain size = %u x %u x %u (%f MB)",
            s->dimx, s->dimy, s->dimz, s->size/1024./1024.);
        MSG("... imaging domain size = %u x %u x %u (%f MB)",
            s->img_dimx, s->img_dimy, s->img_dimz,
            s->size_img/1024./1024.);
    }

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        printf("Current working directory: %s\n", cwd);
    } else {
        perror("getcwd");
        return EXIT_FAILURE;
    }

    struct dirent *de;
    DIR *dr = opendir(dir);
    CHK(dr == NULL, "failed to open img directory");

    char ilm_only_file[128];
    char img_only_file[128];
    char img_dilm_file[128];
    char img_file[128];
    char ilm_file[128];
    char img_pref[64];
    int idx;
    if (s->cpu) {
        sprintf(img_pref, "%s", "img_");
    } else {
        sprintf(img_pref, "%s", "gpu_img_");
    }

    while ((de = readdir(dr)) != NULL) {
        if (strstr(de->d_name, img_pref) != NULL) {
            if (s->cpu) {
                sscanf(de->d_name, "img_%d.raw", &idx);
                sprintf(img_file, "%s/img_%d.raw", dir, idx);
                sprintf(ilm_file, "%s/ilm_%d.raw", dir, idx);
            } else {
                sscanf(de->d_name, "gpu_img_%d.raw", &idx);
                sprintf(img_file, "%s/gpu_img_%d.raw", dir, idx);
                sprintf(ilm_file, "%s/gpu_ilm_%d.raw", dir, idx);
            }

            MSG("... stacking %s and %s", img_file, ilm_file);

            // Read img_shot
            FILE *fd = fopen(img_file, "rb");
            CHK(fd == NULL, "failed to open img file");
            CHK(fread(img_shot, sizeof(float), s->size_img, fd) != s->size_img,
                "failed to read img file");
            fclose(fd);

            // Read ilm_shot
            fd = fopen(ilm_file, "rb");
            CHK(fd == NULL, "failed to open ilm file");
            CHK(fread(ilm_shot, sizeof(float), s->size_img, fd) != s->size_img,
                "failed to read ilm file");
            fclose(fd);

            // Smooth ilm_shot to reduce low-frequency artifacts
            smooth_illumination(s, ilm_shot);

            // Normalize and stack into img
            gather_img_div_ilm_smart(s->size_img, img_shot, ilm_shot, img);

            // Also gather unnormalized images for img_only and ilm_only
            gather_img_ilm(s->size_img, img_shot, ilm_shot, img_only, ilm_only);
        }
    }
    closedir(dr);

    // Apply Laplacian filter to the final img
    apply_laplacian_filter2(s, img);

    // Save the final images
    if (s->cpu) {
        sprintf(img_dilm_file, "%s/img_filtered_opt1.raw", dir);
        sprintf(img_only_file, "%s/img_only.raw", dir);
        sprintf(ilm_only_file, "%s/ilm_only.raw", dir);
    } else {
        sprintf(img_dilm_file, "%s/gpu_img_filtered.raw", dir);
        sprintf(img_only_file, "%s/gpu_img_only.raw", dir);
        sprintf(ilm_only_file, "%s/gpu_ilm_only.raw", dir);
    }

    // Apply thresholding to ilm_only
    for (unsigned int i = 0; i < s->size_img; i++) {
        if (ilm_only[i] > 10) {
            ilm_only[i] = 0.15;
        }
    }

    wave_save_image(s, img, img_dilm_file);
    wave_save_image(s, img_only, img_only_file);
    wave_save_image(s, ilm_only, ilm_only_file);

    DELETE_BUFFER(img);
    DELETE_BUFFER(ilm_only);
    DELETE_BUFFER(img_only);
    DELETE_BUFFER(img_shot);
    DELETE_BUFFER(ilm_shot);

    printf("before free(s)\n");
    free(s);

    printf("before parser_delete(p)\n");
    parser_delete(p);
    return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  /// filtered_each_shot_individually
  /// structure to maintain the user choices.
  sismap_t *s = (sismap_t*)malloc(sizeof(sismap_t));
  /// create a parser.
  parser *p = parser_create("Reverse Time Migration using STENCIL");
  /// parse command line arguments.
  PARSER_BOOTSTRAP(p);
  parser_parse(p, argc, argv);
  s->verbose    = parser_get_bool(p, "verbose");
  s->cpu        = parser_get_bool(p, "cpu");
  s->time_steps = parser_get_int(p, "iter");
  s->cfl        = parser_get_float(p, "cfl");
  s->fmax       = parser_get_float(p, "fmax");
  s->vel_file   = parser_get_string(p, "in");
  s->vel_dimx   = parser_get_int(p, "n1");
  s->vel_dimy   = parser_get_int(p, "n2");
  s->vel_dimz   = parser_get_int(p, "n3");
  s->dx         = parser_get_int(p, "dx");
  s->dy         = parser_get_int(p, "dy");
  s->dz         = parser_get_int(p, "dz");
  s->dcdp       = parser_get_int(p, "dcdp");
  s->dline      = parser_get_int(p, "dline");
  s->drcv       = parser_get_int(p, "drcv");
  s->dshot      = parser_get_int(p, "dshot");
  s->ddepth     = parser_get_int(p, "ddepth");
  s->device     = parser_get_int(p, "device");
  s->first      = parser_get_int(p, "first");
  s->last       = parser_get_int(p, "last");
  s->src_depth  = parser_get_int(p, "src_depth");
  s->rcv_depth  = parser_get_int(p, "rcv_depth");
  s->modeling   = false;
  s->nb_snap    = parser_get_int(p, "nbsnap");
  char *dir     = parser_get_string(p, "dir");
//  char dir="./data";

  /// image for RTM.
  float* img;
  float* img_only;
  float* ilm_only;
  float *ilm_shot, *img_shot;
  float *img_norm; // Temporary buffer for normalized image per shot

  /// initialize the velocity and the compute sizes.
  wave_init_numerics(s);
  wave_init_dimensions(s);
  /// init the buffers:
  CREATE_BUFFER(img, s->size_img);
  NULIFY_BUFFER(img, s->size_img);
  CREATE_BUFFER(img_only, s->size_img);
  NULIFY_BUFFER(img_only, s->size_img);
  CREATE_BUFFER(ilm_only, s->size_img);
  NULIFY_BUFFER(ilm_only, s->size_img);
  CREATE_BUFFER(img_shot, s->size_img);
  CREATE_BUFFER(ilm_shot, s->size_img);
  CREATE_BUFFER(img_norm, s->size_img); // Buffer for normalized image

  /// print info:
  if (s->verbose) {
  	MSG(" ");
  	MSG("... stencil information:");
  	MSG("... compute domain size = %u x %u x %u (%f MB)",
      s->dimx, s->dimy, s->dimz, s->size/1024./1024.);
  	MSG("... imaging domain size = %u x %u x %u (%f MB)",
        s->img_dimx, s->img_dimy, s->img_dimz,
        s->size_img/1024./1024.);
  }
  /// browse the shot images/illuminations:
  char cwd[1024]; // Buffer to hold the current working directory
  // Get the current working directory
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
      printf("Current working directory: %s\n", cwd);
  } else {
      perror("getcwd"); // Print an error message if getcwd fails
      return EXIT_FAILURE;}

  struct dirent *de;
//  DIR *dr = opendir("./data");
  DIR *dr = opendir(dir);
  CHK(dr == NULL, "failed to open img directory");

  char ilm_only_file[128];
  char img_only_file[128];
  char img_dilm_file[128];
  char img_file[128];
  char ilm_file[128];
  char img_pref[64];
  int idx;
  if(s->cpu) {
  	sprintf(img_pref, "%s", "img_");
  } else {
  	sprintf(img_pref, "%s", "gpu_img_");
  }
  while((de = readdir(dr)) != NULL) {
		if (strstr(de->d_name, img_pref) != NULL) {
  		    if(s->cpu) {
				sscanf(de->d_name, "img_%d.raw", &idx);
				sprintf(img_file, "%s/img_%d.raw", dir, idx);
				sprintf(ilm_file, "%s/ilm_%d.raw", dir, idx);
			} else {
				sscanf(de->d_name, "gpu_img_%d.raw", &idx);
				sprintf(img_file, "%s/gpu_img_%d.raw", dir, idx);
				sprintf(ilm_file, "%s/gpu_ilm_%d.raw", dir, idx);
			}
//            MSG('img_file=%s',img_file);
//            MSG('img_file=%s',img_file);
//            MSG('dir=%s',dir);
//            MSG('idx=%s',idx);
			// read img:
  		    MSG("... stacking %s and %s", img_file,ilm_file);
			FILE *fd = fopen(img_file, "rb");
			CHK(fd == NULL, "failed to open img file");
			CHK(fread(img_shot, sizeof(float), s->size_img, fd) != s->size_img,
          "failed to read img file");
			fclose(fd);
			// read ilm:
			fd = fopen(ilm_file, "rb");
			CHK(fd == NULL, "failed to open ilm file");
			CHK(fread(ilm_shot, sizeof(float), s->size_img, fd) != s->size_img,
          "failed to read ilm file");
			fclose(fd);

			//////////////////////////////////////////////
			// Smooth ilm_shot to reduce low-frequency artifacts
//			apply_laplacian_filter2   smooth_illumination
//			apply_laplacian_filter2_small(s, ilm_shot);
//			apply_laplacian_filter2_small(s, img_shot);
//			apply_laplacian_xy_filter(s,img_shot);
			apply_laplacian_filter(s,img_shot);
//			apply_laplacian_filter(s,ilm_shot);

			// Normalize img_shot into img_norm
//			normalize_image(s, img_shot, ilm_shot, img_norm);

//            // Apply Laplacian filter to the normalized image
//            apply_laplacian_filter(s, img_norm);
//
//            // Stack the filtered image into the final img
//            #pragma omp parallel for
//            for (unsigned int i = 0; i < s->size_img; i++) {
//                img[i] += img_norm[i];
//            }
            //////////////////////////////////////////////

            // Also gather unnormalized images for img_only and ilm_only
			gather_img_div_ilm_smart(s->size_img,img_shot,ilm_shot,img);
			gather_img_ilm(s->size_img,img_shot,ilm_shot,img_only,ilm_only);
		}
  }
  closedir(dr);
  /// save the final image on disk.
	if (s->cpu) {
		sprintf(img_dilm_file, "%s/img_filtered_last.raw", dir);
		sprintf(img_only_file, "%s/img_only.raw",dir);
		sprintf(ilm_only_file, "%s/ilm_only.raw",dir);
  } else {
		sprintf(img_dilm_file, "%s/gpu_img.raw", "data");
		sprintf(img_only_file, "%s/gpu_img_only.raw", "data");
		sprintf(ilm_only_file, "%s/gpu_ilm_only.raw", "data");
  }
	for (unsigned int i=0; i<s->size_img; i++) {
  	if (ilm_only[i] > 10) {
//    	printf("%u %f\n", i, ilm_only[i]);
      ilm_only[i]=0.15;
    }
  }
  wave_save_image(s, img,img_dilm_file);
  wave_save_image(s, img_only, img_only_file);
  wave_save_image(s, ilm_only, ilm_only_file);
  DELETE_BUFFER(img);
  DELETE_BUFFER(ilm_only);
  DELETE_BUFFER(img_only);
  DELETE_BUFFER(img_shot);
  DELETE_BUFFER(ilm_shot);
  DELETE_BUFFER(img_norm);

//  printf("before wave_release\n");
  /// release stencil by each variable.not working.
//  wave_release(s);

  /// release the simulation structure.
  printf("before free(s)\n");
  free(s);

  /// delete the parser.
  printf("before parser_delete(p)\n");
  parser_delete(p);
  return EXIT_SUCCESS;
}

int main_original(int argc, char* argv[]) {
  /// original variant
  /// structure to maintain the user choices.
  sismap_t *s = (sismap_t*)malloc(sizeof(sismap_t));
  /// create a parser.
  parser *p = parser_create("Reverse Time Migration using STENCIL");
  /// parse command line arguments.
  PARSER_BOOTSTRAP(p);
  parser_parse(p, argc, argv);
  s->verbose    = parser_get_bool(p, "verbose");
  s->cpu        = parser_get_bool(p, "cpu");
  s->time_steps = parser_get_int(p, "iter");
  s->cfl        = parser_get_float(p, "cfl");
  s->fmax       = parser_get_float(p, "fmax");
  s->vel_file   = parser_get_string(p, "in");
  s->vel_dimx   = parser_get_int(p, "n1");
  s->vel_dimy   = parser_get_int(p, "n2");
  s->vel_dimz   = parser_get_int(p, "n3");
  s->dx         = parser_get_int(p, "dx");
  s->dy         = parser_get_int(p, "dy");
  s->dz         = parser_get_int(p, "dz");
  s->dcdp       = parser_get_int(p, "dcdp");
  s->dline      = parser_get_int(p, "dline");
  s->drcv       = parser_get_int(p, "drcv");
  s->dshot      = parser_get_int(p, "dshot");
  s->ddepth     = parser_get_int(p, "ddepth");
  s->device     = parser_get_int(p, "device");
  s->first      = parser_get_int(p, "first");
  s->last       = parser_get_int(p, "last");
  s->src_depth  = parser_get_int(p, "src_depth");
  s->rcv_depth  = parser_get_int(p, "rcv_depth");
  s->modeling   = false;
  s->nb_snap    = parser_get_int(p, "nbsnap");
  char *dir     = parser_get_string(p, "dir");
//  char dir="./data";

  /// image for RTM.
  float* img;
  float* img_only;
  float* ilm_only;
  float *ilm_shot, *img_shot;

  /// initialize the velocity and the compute sizes.
  wave_init_numerics(s);
  wave_init_dimensions(s);
  /// init the buffers:
  CREATE_BUFFER(img, s->size_img);
  NULIFY_BUFFER(img, s->size_img);
  CREATE_BUFFER(img_only, s->size_img);
  NULIFY_BUFFER(img_only, s->size_img);
  CREATE_BUFFER(ilm_only, s->size_img);
  NULIFY_BUFFER(ilm_only, s->size_img);
  CREATE_BUFFER(img_shot, s->size_img);
  CREATE_BUFFER(ilm_shot, s->size_img);
  /// print info:
  if (s->verbose) {
  	MSG(" ");
  	MSG("... stencil information:");
  	MSG("... compute domain size = %u x %u x %u (%f MB)",
      s->dimx, s->dimy, s->dimz, s->size/1024./1024.);
  	MSG("... imaging domain size = %u x %u x %u (%f MB)",
        s->img_dimx, s->img_dimy, s->img_dimz,
        s->size_img/1024./1024.);
  }
  /// browse the shot images/illuminations:
  char cwd[1024]; // Buffer to hold the current working directory
  // Get the current working directory
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
      printf("Current working directory: %s\n", cwd);
  } else {
      perror("getcwd"); // Print an error message if getcwd fails
      return EXIT_FAILURE;}

  struct dirent *de;
//  DIR *dr = opendir("./data");
  DIR *dr = opendir(dir);
  CHK(dr == NULL, "failed to open img directory");

  char ilm_only_file[128];
  char img_only_file[128];
  char img_dilm_file[128];
  char img_file[128];
  char ilm_file[128];
  char img_pref[64];
  int idx;
  if(s->cpu) {
  	sprintf(img_pref, "%s", "img_");
  } else {
  	sprintf(img_pref, "%s", "gpu_img_");
  }
  while((de = readdir(dr)) != NULL) {
		if (strstr(de->d_name, img_pref) != NULL) {
  		    if(s->cpu) {
				sscanf(de->d_name, "img_%d.raw", &idx);
				sprintf(img_file, "%s/img_%d.raw", dir, idx);
				sprintf(ilm_file, "%s/ilm_%d.raw", dir, idx);
			} else {
				sscanf(de->d_name, "gpu_img_%d.raw", &idx);
				sprintf(img_file, "%s/gpu_img_%d.raw", dir, idx);
				sprintf(ilm_file, "%s/gpu_ilm_%d.raw", dir, idx);
			}
//            MSG('img_file=%s',img_file);
//            MSG('img_file=%s',img_file);
//            MSG('dir=%s',dir);
//            MSG('idx=%s',idx);
			// read img:
  		    MSG("... stacking %s and %s", img_file,ilm_file);
			FILE *fd = fopen(img_file, "rb");
			CHK(fd == NULL, "failed to open img file");
			CHK(fread(img_shot, sizeof(float), s->size_img, fd) != s->size_img,
          "failed to read img file");
			fclose(fd);
			// read ilm:
			fd = fopen(ilm_file, "rb");
			CHK(fd == NULL, "failed to open ilm file");
			CHK(fread(ilm_shot, sizeof(float), s->size_img, fd) != s->size_img,
          "failed to read ilm file");
			fclose(fd);
			// do the shot gather:
			gather_img_div_ilm_smart(s->size_img,img_shot,ilm_shot,img);
			gather_img_ilm(s->size_img,img_shot,ilm_shot,img_only,ilm_only);
		}
  }
  closedir(dr);
  /// save the final image on disk.
	if (s->cpu) {
		sprintf(img_dilm_file, "%s/img.raw", dir);
		sprintf(img_only_file, "%s/img_only.raw",dir);
		sprintf(ilm_only_file, "%s/ilm_only.raw",dir);
  } else {
		sprintf(img_dilm_file, "%s/gpu_img.raw", "data");
		sprintf(img_only_file, "%s/gpu_img_only.raw", "data");
		sprintf(ilm_only_file, "%s/gpu_ilm_only.raw", "data");
  }
	for (unsigned int i=0; i<s->size_img; i++) {
  	if (ilm_only[i] > 10) {
//    	printf("%u %f\n", i, ilm_only[i]);
      ilm_only[i]=0.15;
    }
  }
  wave_save_image(s, img,img_dilm_file);
  wave_save_image(s, img_only, img_only_file);
  wave_save_image(s, ilm_only, ilm_only_file);
  DELETE_BUFFER(img);
  DELETE_BUFFER(ilm_only);
  DELETE_BUFFER(img_only);
  DELETE_BUFFER(img_shot);
  DELETE_BUFFER(ilm_shot);

//  printf("before wave_release\n");
  /// release stencil by each variable.not working.
//  wave_release(s);

  /// release the simulation structure.
  printf("before free(s)\n");
  free(s);

  /// delete the parser.
  printf("before parser_delete(p)\n");
  parser_delete(p);
  return EXIT_SUCCESS;
}
