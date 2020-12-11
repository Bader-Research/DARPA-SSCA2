/* 
   memread.h
   Test reading/writing of memory
*/

#ifndef MEMREAD_H
#define MEMREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "simple.h"
#include "alg_random.h"


/* Defaults */
#if RAND_MAX < 65536
#undef RAND_MAX
#define RAND_MAX 2147483647
#endif
#ifndef bool
typedef enum {false, true} bool;
#endif
/* #define INLINE inline */
#define INLINE

#endif /* MEMREAD_H */
