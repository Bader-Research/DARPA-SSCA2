/******************************************************************************
* FILE: mpi_prime.c
* OTHER FILES: make.mpi_prime.c
* DESCRIPTION:
*   Generates prime numbers.  All tasks distribute the work evenly, taking
*   every nth number, where n is the stride computed as:  (rank *2) + 1
*   so that even numbers are automatically skipped.  The method of using
*   stride is preferred over contiguous blocks of numbers, since numbers
*   in the higher range require more work to compute and may result in 
*   load imbalance.  This program demonstrates embarrassing parallelism.
*   Collective communications calls are used to reduce the only two data
*   elements requiring communications: the number of primes found and
*   the largest prime.
* AUTHOR: Blaise Barney 11/25/95 - adapted from version contributed by 
*   Richard Ng &  Wong Sze Cheong during MHPCC Singapore Workshop (8/22/95).
* LAST REVISED: 
******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "alg_prime.h"

#define LIMIT     2500000     /* Increase this to find more primes */
#define FIRST     0           /* Rank of first task */

int isprime(int n) {
int i,squareroot;
if (n>10) {
   squareroot = (int) sqrt((double)n);
   for (i=3; i<=squareroot; i=i+2)
      if ((n%i)==0)
         return 0;
   return 1;
   }
/* Assume first four primes are counted elsewhere. Forget everything else */
else {
  if ((n==2)||(n==3)||(n==5)||(n==7))
    return 1;
  else
    return 0;
}
}


void
all_prime_mpi() {

int   ntasks,               /* total number of tasks in partitiion */
      rank,                 /* task identifier */
      n,                    /* loop variable */
      pc,                   /* prime counter */
      pcsum,                /* number of primes found by all tasks */
      foundone,             /* most recent prime found */
      maxprime,             /* largest prime found */
      mystart,              /* where to start calculating */
      stride;               /* calculate every nth number */

double start_time,end_time;

ntasks = NODES;
rank   = MYNODE;

if (((ntasks%2) !=0) || ((LIMIT%ntasks) !=0)) {
  fprintf(stderr,
	  "Sorry - this exercise requires an even number of processors\n");
  fprintf(stderr,
	  "evenly divisible into %d.  Try 4 or 8.\n",LIMIT);
  return;
}

MPI_Barrier(MPI_COMM_WORLD);
start_time = MPI_Wtime();	/* Initialize start time */
mystart = (rank*2)+1;       /* Find my starting point - must be odd number */
stride = ntasks*2;          /* Determine stride, skipping even numbers */
pc=0;                       /* Initialize prime counter */
foundone = 0;               /* Initialize */

if (rank == FIRST) {
   printf("Using %d tasks to scan %d numbers\n",ntasks,LIMIT);
}

for (n=mystart; n<=LIMIT; n=n+stride) {
  if (isprime(n)) {
    pc++;
    foundone = n;
#if 0
    /* Optional: print each prime as it is found */ 
    printf("PE%3d: %d  (pc: %d)\n",MYNODE,foundone,pc);
#endif
  }
}

MPI_Reduce(&pc,&pcsum,1,MPI_INT,MPI_SUM,FIRST,MPI_COMM_WORLD);
MPI_Reduce(&foundone,&maxprime,1,MPI_INT,MPI_MAX,FIRST,MPI_COMM_WORLD);

if (rank == FIRST) {
  end_time=MPI_Wtime();
  printf("Done. Largest prime is %d Total primes %d\n",maxprime,pcsum);
  printf("Wallclock time elapsed: %.3lf\n",end_time-start_time);
}

return;
}

void
all_prime(THREADED) {

int   ntasks,               /* total number of tasks in partitiion */
      rank,                 /* task identifier */
      n,                    /* loop variable */
      pc,                   /* prime counter */
      pcsum,                /* number of primes found by all tasks */
      foundone,             /* most recent prime found */
      maxprime,             /* largest prime found */
      mystart,              /* where to start calculating */
      stride;               /* calculate every nth number */

double start_time,end_time;

ntasks = TID;
rank   = ID;

if (((ntasks%2) !=0) || ((LIMIT%ntasks) !=0)) {
  fprintf(stderr,
	  "Sorry - this exercise requires an even number of processors\n");
  fprintf(stderr,
	  "evenly divisible into %d.  Try 4 or 8.\n",LIMIT);
  return;
}

all_Barrier(TH);

start_time = get_seconds();	/* Initialize start time */
mystart = (rank*2)+1;       /* Find my starting point - must be odd number */
stride = ntasks*2;          /* Determine stride, skipping even numbers */
pc=0;                       /* Initialize prime counter */
foundone = 0;               /* Initialize */

if (rank == FIRST) {
   printf("Using %d tasks to scan %d numbers\n",ntasks,LIMIT);
}

for (n=mystart; n<=LIMIT; n=n+stride) {
  if (isprime(n)) {
    pc++;
    foundone = n;
#if 0
    /* Optional: print each prime as it is found */ 
    printf("%d\n",foundone);
#endif
  }
}

pcsum    = all_Reduce_i(pc, SUM, TH);
maxprime = all_Reduce_i(foundone, MAX, TH);

if (rank == FIRST) {
  end_time=get_seconds();
  printf("Done. Largest prime is %d Total primes %d\n",maxprime,pcsum);
  printf("Wallclock time elapsed: %.3lf\n",end_time-start_time);
}

return;
}
