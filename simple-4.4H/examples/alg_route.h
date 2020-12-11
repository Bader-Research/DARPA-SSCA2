#ifndef _ALG_ROUTE_H
#define _ALG_ROUTE_H

#include "simple.h"
#include "fastest.h"

int all_route_q(int q, 
		int *lKeys,
		int *lAddr,
		int *lRouted);

int all_route_q_mpi(int q, 
		    int *lKeys,
		    int *lAddr,
		    int *lRouted);

int all_route_q_simple(int q, 
		       int *lKeys,
		       int *lAddr,
		       int *lRouted,
		       THREADED);

extern int init_Alltoallv_param;
void init_Alltoallv();
void destroy_Alltoallv();

#endif
