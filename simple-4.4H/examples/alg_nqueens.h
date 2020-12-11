#ifndef _ALG_NQUEENS_H
#define _ALG_NQUEENS_H

#include <stdio.h>
#include <stdlib.h>
#include "simple.h"

/* R = recursive */
/* NR = non-recursive */

int nqueens_R_orig(int N);         /* sequential, simple DFS */
int nqueens_NR_orig(int N);         /* sequential, simple DFS */
int nqueens_R(int N, int k);       /* sequential, DFS from N^k positions */
int nqueens_NR(int N, int k);       /* sequential, DFS from N^k positions */


int all_nqueens_R_block(int N);
    /* message passing, first level partition */
int all_nqueens_NR_block(int N);
    /* message passing, first level partition */

int all_nqueens_R_cyclic(int N, int k);
    /* message passing, DFS from N^k positions, cyclicly distributed */
int all_nqueens_NR_cyclic(int N, int k);
    /* message passing, DFS from N^k positions, cyclicly distributed */

int all_nqueens_R_random(int N, int k);
    /* message passing, DFS from N^k positions, randomly distributed */
int all_nqueens_NR_random(int N, int k);
    /* message passing, DFS from N^k positions, randomly distributed */


int all_nqueens_R_block_r(int N, int k, THREADED);
    /* threaded MP, DFS from N^k positions, block distributed */
int all_nqueens_NR_block_r(int N, int k, THREADED);
    /* threaded MP, DFS from N^k positions, block distributed */

int all_nqueens_R_cyclic_r(int N, int k, THREADED);
    /* threaded MP, DFS from N^k positions, cyclicly distributed */
int all_nqueens_NR_cyclic_r(int N, int k, THREADED);
    /* threaded MP, DFS from N^k positions, cyclicly distributed */

int all_nqueens_R_random_r(int N, int k, THREADED); 
    /* threaded MP, DFS from N^k positions, randomly distributed */
int all_nqueens_NR_random_r(int N, int k, THREADED); 
    /* threaded MP, DFS from N^k positions, randomly distributed */

#endif
