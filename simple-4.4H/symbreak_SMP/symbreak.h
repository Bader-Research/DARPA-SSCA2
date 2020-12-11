/* 
   graph_color.h
   pfh 9/29/99

   This is the header file for my implementation of JaJa's 
   3-coloring algorithm, under the tutelage of D Bader.

   Moving defines, etc, here so that we can reference this from the 
   sorting source code. 

*/

#ifndef GRAPH_COLOR_H
#define GRAPH_COLOR_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "simple.h"
#include "alg_random.h"


/* DEC unix hacks */
/* #define DEC_UNIX */

#ifdef DEC_UNIX

/* Incorrect RAND_MAX */
#if RAND_MAX < 65536
#undef RAND_MAX
#define RAND_MAX 2147483647
#endif

/* Doesn't know the new ANSI bool type, either */
#ifndef bool
typedef enum {false, true} bool;
#endif

/* Yet Another Hack for cruddy DEC C - compiler doesn't understand */
/* inlining! Ye gods. */
#define INLINE 
#else

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

#endif


/* Const - if a graph is bigger than this, then don't print it out. */
const int SIZE_TOOBIG = 16;

/* Struct used in sorting by color */
typedef struct 
{
  int c; /* color */
  int v; /* vertex index */
} s_struct;

/* Const - number of bits in histogram. */
/* For 64-bit ints or less, 8 is enough. */
const int NUM_HBITS = 8;
const int NUM_HBINS = 256;

/* ---------------------------------------------------------------------------- */
/* Routine to return the kth LSB of an integer. */
/* No error checking, very fast. */
/* New 10-19-99 -- per DAB, change this fn to a macro. */
#define kth_lsb(in_num, k) (((in_num) >> (k)) & 0x1)

/* ---------------------------------------------------------------------------- */
/* Bader code - macro for sorting routine. */
#define bits(x,k,j) (((x)>>(k)) & ~(~0<<(j)))

#endif /* Graph_color_H */
