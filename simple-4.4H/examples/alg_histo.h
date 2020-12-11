#ifndef _ALG_HISTO_H
#define _ALG_HISTO_H

#include <stdio.h>
#include <stdlib.h>
#include "simple.h"

int all_histo(int myL, int *A, int R, int *histo);
int all_histo_r(int myL, int *A, int R, int *histo, THREADED); 

int all_sel_sum(int myL, int *A);
int all_sel_sum_r(int myL, int *A, THREADED); 

#endif
