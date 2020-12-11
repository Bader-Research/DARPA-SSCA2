#ifndef _ALG_RADIX_H
#define _ALG_RADIX_H

#include "simple.h"

void all_countsort(int q,
		   int *lKey,
		   int *lSorted,
		   int R,
		   int bitOff, int m);
void all_countsort_mpi(int q,
		       int *lKey,
		       int *lSorted,
		       int R,
		       int bitOff, int m);
void all_countsort_simple(int q,
			  int *lKey,
			  int *lSorted,
			  int R,
			  int bitOff, int m,
			  THREADED);
void all_countsort_node(int q,
			int *lKey,
			int *lSorted,
			int R,
			int bitOff, int m,
			THREADED);

#define all_radixsort(a,b,c)   all_radixsort_s3(a,b,c)
void all_radixsort_s3(int q,
		      int *lKeys,
		      int *lSorted);
void all_radixsort_s2(int q,
		      int *lKeys,
		      int *lSorted);
void all_radixsort20_s1(int q,
			int *lKeys,
			int *lSorted);
void all_radixsort20_s2(int q,
			int *lKeys,
			int *lSorted);

#define all_radixsort_mpi(a,b,c)   all_radixsort_mpi_s3(a,b,c)
void all_radixsort_mpi_s3(int q,
			  int *lKeys,
			  int *lSorted);
void all_radixsort_mpi_s2(int q,
			  int *lKeys,
			  int *lSorted);
void all_radixsort20_mpi_s1(int q,
			    int *lKeys,
			    int *lSorted);
void all_radixsort20_mpi_s2(int q,
			    int *lKeys,
			    int *lSorted);

#define all_radixsort_simple(a,b,c,d)   all_radixsort_simple_s3(a,b,c,d)
void all_radixsort_simple_s3(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);
void all_radixsort_simple_s2(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);
void all_radixsort20_simple_s1(int q,
			       int *lKeys,
			       int *lSorted,
			       THREADED);
void all_radixsort20_simple_s2(int q,
			       int *lKeys,
			       int *lSorted,
			       THREADED);

#define all_radixsort_node(a,b,c,d)   all_radixsort_node_s3(a,b,c,d)
void all_radixsort_node_s3(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED);
void all_radixsort_node_s2(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED);
void all_radixsort20_node_s1(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);
void all_radixsort20_node_s2(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED);

void all_radixsort_check(int q,
			 int *lSorted);

#endif

