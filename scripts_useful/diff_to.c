#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <libgen.h>

// structure to track MIN/MAX:
typedef struct __entry {
	size_t x;
	float  v;
} entry_t;

// main routine of the program:
int main(int argc, char **argv) {
  // file descriptor and names:
  FILE *fd;
  char dif_filename[256];
  char ref_filename[256];
  char gpu_filename[256];
  char tmp_filename[256];
  
  // tabs:
	float *gpu_tab;
	float *ref_tab;
	float *dif_tab;
	
	// variables for MIN/MAX:
	entry_t ref_min;
	entry_t ref_max;
	entry_t dif_min;
	entry_t dif_max;
	entry_t gpu_min;
	entry_t gpu_max;
	ref_min.v = FLT_MAX;
	ref_max.v = FLT_MIN;
	gpu_min.v = FLT_MAX;
	gpu_max.v = FLT_MIN;
	dif_min.v = FLT_MAX;
	dif_max.v = FLT_MIN;

	// RMS:
	float ref_rms=0.0;
	float gpu_rms=0.0;
	float dif_rms=0.0;

	// index and size:
	size_t tmp, sz, idx, n;
	
	// command line parsing:
	if (argc != 3 ) {
		fprintf(stderr, 
						"Usage: %s gpu_file ref_file\n", argv[0]);
    exit(EXIT_FAILURE);   
  }
  sprintf(gpu_filename, "%s", argv[1]);
  sprintf(tmp_filename, "%s", argv[1]);
  sprintf(ref_filename, "%s", argv[2]);
  sprintf(dif_filename, "%s/diff_%s",
          dirname(tmp_filename), basename(tmp_filename));
  // read from files and detect sizes:
  if ((fd=fopen(ref_filename, "rb")) == NULL) {
    fprintf(stderr, "... failed to open file %s\n", ref_filename);
    return EXIT_FAILURE;
  }
	// update size and count:
  fseek(fd, 0, SEEK_END);
  sz = ftell(fd);
  rewind(fd);
  n  = sz/sizeof(float);
  // allocate the tabs:
	dif_tab = (float*)malloc(n*sizeof(float));
	ref_tab = (float*)malloc(n*sizeof(float));
	gpu_tab = (float*)malloc(n*sizeof(float));
  
  if (fread(ref_tab, sizeof(float), n, fd) != n ) {
  	fprintf(stderr, "... failed to read from file %s\n", ref_filename);
  	return EXIT_FAILURE;
  }
  fclose(fd);
  // check files have the same size:
  if ((fd=fopen(gpu_filename, "rb")) == NULL) {
  	fprintf(stderr, "... failed to open file %s\n", gpu_filename);
  	return EXIT_FAILURE;
  }
  fseek(fd, 0, SEEK_END);
  tmp = ftell(fd);
  rewind(fd);
    
  if (tmp != sz) {
    fprintf(stderr, "... files should be of same size, aborting.\n");
    return EXIT_FAILURE;    
  }
  if (fread(gpu_tab, sizeof(float), n, fd) != n ) {
  	fprintf(stderr, "... failed to read from file %s\n", gpu_filename);
  	return EXIT_FAILURE;
  }
  fclose(fd);
  // calculate the difference, min, max and rms:
  for (idx=0; idx<n; ++idx) {
    dif_tab[idx]=ref_tab[idx]-gpu_tab[idx];
    // update ref MAX:
    if (ref_tab[idx] > ref_max.v) {
      ref_max.v = ref_tab[idx];
      ref_max.x = idx;
    }	
    // update gpu MAX:
    if (gpu_tab[idx] > gpu_max.v) {
      gpu_max.v = gpu_tab[idx];
      gpu_max.x = idx;
    }	
    // update diff MAX:
    if (dif_tab[idx] > dif_max.v) {
      dif_max.v = dif_tab[idx];
      dif_max.x = idx;
    }	
	  // update ref MIN:
    if (ref_tab[idx] < ref_min.v) {
      ref_min.v = ref_tab[idx];
      ref_min.x = idx;
    }	
	  // update gpu MIN:
    if (gpu_tab[idx] < gpu_min.v) {
      gpu_min.v = gpu_tab[idx];
      gpu_min.x = idx;
    }
    // update diff MIN:
    if (dif_tab[idx] < dif_min.v) {
      dif_min.v = dif_tab[idx];
      dif_min.x = idx;
    }		
    // RMS:
    ref_rms+=ref_tab[idx]*ref_tab[idx];
    gpu_rms+=gpu_tab[idx]*gpu_tab[idx];
    dif_rms+=dif_tab[idx]*dif_tab[idx];
  }
  // finalize RMS:
  ref_rms=sqrtf(ref_rms/n);
  gpu_rms=sqrtf(gpu_rms/n);
  dif_rms=sqrtf(dif_rms/n);
  // write the diff file:
  if ((fd=fopen(dif_filename, "wb")) == NULL) {
  	fprintf(stderr, "... failed to open file %s\n", dif_filename);
  	return EXIT_FAILURE;
  }
  if ( fwrite(dif_tab, sizeof(float), n, fd) != n ) {
  	fprintf(stderr, "... failed to write to file %s\n", dif_filename);
  	return EXIT_FAILURE;
  }
  fclose(fd);
  // clean temporary arrays:
  free(ref_tab);
  free(gpu_tab);
  free(dif_tab);
  // report results:
  fprintf(stdout, "...\n");
  fprintf(stdout, "... [DIFF TAB] generated is : %s (%lu elements)\n", 
          dif_filename, n);
  fprintf(stdout, "...\n");
  fprintf(stdout, 
  				"... [GPU_TAB] min: %20.7g @ [%20lu] | max: %20.7g @ [%20lu] | rms: %20.7g\n",
  				gpu_min.v, gpu_min.x,
  				gpu_max.v, gpu_max.x, gpu_rms);
  fprintf(stdout, 
  				"... [REF_TAB] min: %20.7g @ [%20lu] | max: %20.7g @ [%20lu] | rms: %20.7g\n",
  				ref_min.v, ref_min.x,
  				ref_max.v, ref_max.x, ref_rms);
  fprintf(stdout, 
  				"... [DIF_TAB] min: %20.7g @ [%20lu] | max: %20.7g @ [%20lu] | rms: %20.7g\n",
  				dif_min.v, dif_min.x,
  				dif_max.v, dif_max.x, dif_rms);
  fprintf(stdout, "...\n");
  // finish the program:
	return EXIT_SUCCESS;
}
