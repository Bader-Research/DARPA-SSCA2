#ifndef _ALG_2DFFT_H
#define _ALG_2DFFT_H

#include "simple.h"

#define EPSILON 0.00001       /* for comparing fp numbers */
#ifndef PI
#define PI 3.14159265358979   /* 4*atan(1.0) */
#endif

#define MAX_FFT_SIZE          (1<<10)

typedef struct {
  double r;             /* real      part */
  double i;             /* imaginary part */
} complex_t; 

/* swap a pair of complex numbers */
#define CMPLXSWAP(a,b) {double swap_temp=(a).r;(a).r=(b).r;(b).r=swap_temp;\
			       swap_temp=(a).i;(a).i=(b).i;(b).i=swap_temp;}


void
all_2dfft_mpi();

void
all_2dfft(int points, THREADED);

#endif
