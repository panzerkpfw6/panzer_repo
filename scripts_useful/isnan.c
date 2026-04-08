#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

// structure to track MIN/MAX:
typedef struct __entry {
	size_t idx;
	float  v;
} entry_t;

// main routine of the program:
int main(int argc, char **argv) {
  // file descriptor and names:
  FILE *fd;
  char filename[256];
  
  // tabs:
	float *ref_tab;
	
	// variables for MIN/MAX:
	entry_t ref_min;
	entry_t ref_max;
	ref_min.v = FLT_MAX;
	ref_max.v = FLT_MIN;

	// RMS:
	float ref_rms=0.0;

	// index and size:
	size_t sz, N, not_valid=0;
	
	// command line parsing:
	if (argc != 2 ) {
		fprintf(stderr, 
						"Usage: %s file_name\n", argv[0]);
    exit(EXIT_FAILURE);   
  }
  sprintf(filename, "%s", argv[1]);
	// open the file:
  if ((fd=fopen(filename, "rb")) == NULL) {
  	fprintf(stderr, "... failed to open file %s\n", filename);
  	return EXIT_FAILURE;
  }
	// update size and count:
	fseek(fd, 0, SEEK_END);
	sz = ftell(fd);
	rewind(fd);
 	N  = sz/sizeof(float);
  // allocate the tabs:
	ref_tab = (float*)malloc(sz);
  // read from ref and gpu:
  if (fread(ref_tab, sizeof(float), N, fd) != N ) {
  	fprintf(stderr, "... failed to read from file %s\n", filename);
  	return EXIT_FAILURE;
  }
  fclose(fd);
  // calculate the difference, min, max and rms:
  for (size_t idx=0; idx<N; ++idx) {
  	// MAX:
		if (ref_tab[idx] > ref_max.v) {
  		ref_max.v   = ref_tab[idx];
  		ref_max.idx = idx;
  	}
		// MIN:
  	if (ref_tab[idx] < ref_min.v) {
  		ref_min.v   = ref_tab[idx];
  		ref_min.idx = idx;
  	}	
  	// RMS:
  	ref_rms+=ref_tab[idx]*ref_tab[idx];
  	// nan:
    if(isnan(ref_tab[idx])) not_valid++;
	}
  // finalize RMS:
  ref_rms=sqrtf(ref_rms/N);
  // clean temporary arrays:
  free(ref_tab);
  // report results:
  fprintf(stdout,"...\n"); 
  fprintf(stdout,"... file: %s\n", filename); 
  fprintf(stdout,"... size: %lu\n", sz); 
  fprintf(stdout,"... N   : %lu\n", N); 
  fprintf(stdout, 
  				"... min : %g @ %4lu\n", ref_min.v, ref_min.idx);
  fprintf(stdout, 
  				"... max : %g @ %4lu\n", ref_max.v, ref_max.idx);
  fprintf(stdout, 
  				"... rms : %g\n", ref_rms);
  fprintf(stdout, 
  				"... nan : %lu (%f %%)\n", not_valid, 
          (float)not_valid*100.0/(float)N);
  fprintf(stdout,"...\n"); 
  // finish the program:
	return EXIT_SUCCESS;
}
