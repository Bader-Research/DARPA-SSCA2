#ifndef _ALG_RADIX_H
#define _ALG_RADIX_H

#include "simple.h"

void all_countsort_smp(int q,
			int *lKey,
			int *lSorted,
			int R,
			int bitOff, int m,
			THREADED);

#define all_radixsort_smp(a,b,c,d)   all_radixsort_smp_s3(a,b,c,d)
void all_radixsort_smp_s3(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED);
void all_radixsort_smp_s2(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED);
void all_radixsort20_smp_s1(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);
void all_radixsort20_smp_s2(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);

void all_radixsort_check(int q,
			 int *lSorted);

#endif

