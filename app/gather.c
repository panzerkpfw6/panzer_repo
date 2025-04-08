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

void gather_img_ilm(unsigned int len,
                    float *img_shot, float *ilm_shot, float *img, float *ilm) {
  #pragma omp parallel for
  for (unsigned int i = 0; i < len; i++) {
    img[i] += img_shot[i];
    ilm[i] += ilm_shot[i];
  }
}

int main(int argc, char* argv[]) {
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
  wave_init_acquisition(s);
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
			gather_img_div_ilm(s->size_img, img_shot, ilm_shot, img);
			gather_img_ilm(s->size_img, img_shot, ilm_shot, img_only, ilm_only);
		}
  }
  closedir(dr);
  /// save the final image on disk.
	if (s->cpu) {
		sprintf(img_dilm_file, "%s/img.raw", "data");
		sprintf(img_only_file, "%s/img_only.raw", "data");
		sprintf(ilm_only_file, "%s/ilm_only.raw", "data");
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
