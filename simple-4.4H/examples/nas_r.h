#ifndef _NAS_R_H
#define _NAS_R_H

#include "simple.h"

#define _NAS_BITS 19
#define _NAS_SEED 314159265.00
#define _NAS_MULT 1220703125.00

double   find_my_seed( long kn,       /* my processor rank, 0<=kn<=num procs */
                       long np,       /* np = num procs                      */
                       long nn,       /* total num of ran numbers, all procs */
                       double s,      /* Ran num seed, for ex.: 314159265.00 */
                       double a );    /* Ran num gen mult, try 1220703125.00 */

void	create_seq( double seed, double a , int q, int *arr);

void	create_seq_simple( double seed, double a , int q, int *arr, THREADED);

void	create_seq_random_simple( double seed, double a , int q, int *arr, THREADED);

#endif

#if 0

USAGE:
    create_seq( find_my_seed( my_rank, 
                              comm_size, 
                              4*TOTAL_KEYS,
                              314159265.00,      /* Random number gen seed */
                              1220703125.00 ),   /* Random number gen mult */
                              1220703125.00 );   /* Random number gen mult */

#endif
