#include "simple.h"
#include "alg_radix.h"
#include "alg_create_input.h"

#define DEBUG              0
#define TIMING             1

#define TEST_WIDE          1

#if TEST_WIDE
#define ARR_SIZE_SORT   (1<<23)
#else
#define ARR_SIZE_SORT   (1<<14)
#endif
#define MIN_TIME       0.000001



void
test_radixsort_smp(int N1, THREADED) {
  int *inArr, *outArr;

#if TIMING
  double secs, tsec;
#endif

  inArr  = (int *)node_malloc(N1 * sizeof(int), TH);
  outArr = (int *)node_malloc(N1 * sizeof(int), TH);

  create_input_nas_smp(N1, inArr, TH);

#if TIMING
  node_Barrier();
  secs = get_seconds();
#endif

  all_radixsort_smp(N1, inArr, outArr, TH);

#if TIMING
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = node_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(stdout,"T: %3d n: %12d  SSort: %9.6lf\n",
	    THREADS, N1, tsec);
    fflush(stdout);
  }
#endif

  node_Barrier();

  on_one all_radixsort_check(N1,  outArr);

  node_Barrier();

  node_free(outArr, TH);
  node_free(inArr, TH);
}



void *SIMPLE_main(THREADED)
{
  int  i;

#if TEST_WIDE
#define TEST_INC (i<<1)
#else
#define TEST_INC (i+=8)
#endif
  
#if DEBUG
  fprintf(stdout,"PE%3d: SMP_main()\n",MYTHREAD);
  fflush(stdout);
#endif

  node_Barrier();
  
  for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    test_radixsort_smp(i, TH);

  node_Barrier();

  /*******************************/
  /* End of program              */
  /*******************************/
  
  SIMPLE_done(TH);
}

