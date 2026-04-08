
/*
#include <sys/time.h>
#include "wtime.h"

double wtime() {
  double t;
  static int sec = -1;
  struct timeval tv;

  gettimeofday(&tv, (void *)0);

  if (sec < 0) sec = tv.tv_sec;

  t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;

  return t;
}
*/
# include <stdlib.h>
# include <stdio.h>
#include <sys/time.h>
# include <time.h>
# include "stencil/wtime.h"
/******************************************************************************/

struct timespec start;

void wtime_init() {
   clock_gettime(CLOCK_REALTIME, &start);
}

double wtime ( ) {
  double t;
  struct timespec finish;
  clock_gettime(CLOCK_REALTIME,&finish);

  long seconds = finish.tv_sec - start.tv_sec;
  long ns = finish.tv_nsec - start.tv_nsec;

  if (start.tv_nsec > finish.tv_nsec) { // clock underflow
    --seconds;
    ns += 1000000000;
  }

  t = (double)seconds + (double)ns/(double)1000000000;
  return t;
}
