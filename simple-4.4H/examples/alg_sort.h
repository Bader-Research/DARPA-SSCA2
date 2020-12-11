#ifndef _ALG_SORT_H
#define _ALG_SORT_H

#include "umd.h"

int *seq_radixsort1(int *a, int *b, int n);
int *seq_radixsort2(int *a, int *b, int n);
int *seq_radixsort2_overlap(int *a, int *b, int n);
int *seq_radixsort_16(int *a, int *b, int n);
int *seq_radixsort_16_overlap(int *a, int *b, int n);

void check_seq_sort(int *a, int n, char *s);

void rsort_input_s(int n, int *x, int *y,
		   int *(*fun)(), char* s);

#define bits(x,k,j) ((x>>k) & ~(~0<<j))

#endif

