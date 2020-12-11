#ifndef _ALG_RADIX_H
#define _ALG_RADIX_H

#include "simple.h"

void all_countsort_node(int q,
			int *lKey,
			int *lSorted,
			int R,
			int bitOff, int m,
			THREADED);

#define all_radixsort_node(a,b,c,d)   all_radixsort_node_s3(a,b,c,d)
void all_radixsort_node_s3(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED);

void all_countsort_node_aux(int q,
			    int *lKey,
			    int *lSorted, int* auxKeys, int* auxSorted,
			    int R,
			    int bitOff, int m,
			    THREADED);

#define all_radixsort_node_aux(a,b,c,d)   all_radixsort_node_aux_s3(a,b,c,d)

void all_radixsort_node_aux_s3(int q,
			       int *lKeys,
			       int *lSorted, int* auxKeys, int* auxSorted,
			       THREADED);

#endif

