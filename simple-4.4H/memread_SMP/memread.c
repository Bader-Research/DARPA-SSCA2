#include "memread.h"
#include <unistd.h>

#define DO_SERIAL 1

#define DO_R 1
#define DO_S 1
#define DO_C 1

#define MIN_SIZE (1 << 6) /* 16 */
#define MAX_SIZE (1 << 15) /* 25 */
#define NUM_ITER       1000

#define DEBUG 0

#define _INP_RANDOM   0
#define _INP_SEQ      1
#define _INP_CYCLIC02 2
#define _INP_CYCLIC04 3
#define _INP_CYCLIC08 4
#define _INP_CYCLIC16 5
#define _INP_CYCLIC32 6
#define _INP_CYCLIC64 7


char *inputLabel(const int input) {
  char *a;
  switch(input) {
  case _INP_RANDOM:   a = "  R"; break;
  case _INP_SEQ:      a = "  S"; break;
  case _INP_CYCLIC02: a = "C02"; break;
  case _INP_CYCLIC04: a = "C04"; break;
  case _INP_CYCLIC08: a = "C08"; break;
  case _INP_CYCLIC16: a = "C16"; break;
  case _INP_CYCLIC32: a = "C32"; break;
  case _INP_CYCLIC64: a = "C64"; break;
  default: a = "E";
  }
  return(a);
}

/* ---------------------------------------------------------------------------- */
/* Function from D Bader to permute a given array.  */
/* This does an in-place modification of the argument, and expects an array of n+1  */
/* points from 0 to n. */
void permute_array(int *S, const int num_pts)
{
    int		i, rnum, itmp;
    double  dtmp;
	
    /* 
       Seed the random-number generator with current time so that
       the numbers will be different every time we run.
       Code snipped from MSDN.

       Per DAB, use srandom and random instead of rand/srand. Not in 
       C standard, but better PRNG.
    */

    /* srandom should return zero if init'd ok */
#if 0
    assert((srandom(time(NULL)) == 0));
#else
    srandom(time(NULL));
#endif

    for(i = num_pts - 1; i >= 1; i--)
    {
      /* Get random number from 0 to 1; map to [0, i] */
	dtmp = (double) random() / (double)RAND_MAX;
	dtmp = dtmp * i;
	rnum = (int) floor(dtmp);

	assert(rnum >= 0);
	assert(rnum <  num_pts);

	/* Swap s[i], s[j] */
	itmp = S[i];
	S[i] = S[rnum];
	S[rnum] = itmp;
    }

    return;
}

/* ---------------------------------------------------------------------------- */
/* Generate a graph, in the format of Ja Ja's successor array. */
/* Caller responsible for freeing the returned vector. */
int *generate_array(const int num_pts, const int input, THREADED)
{
    int *result;
    int i;
    int	*list;
    
    result = (int *) node_malloc(sizeof(int) * num_pts , TH);

    switch(input) {
    case _INP_RANDOM : {
      on_one_thread {
	list = (int *) malloc(sizeof(int) * num_pts);
	assert(list != NULL);

	for (i=0 ; i<num_pts ; i++) {
	  list[i] = i;
	}
	
	permute_array(list, num_pts);
	
	for(i = 0; i < num_pts-1 ; i++) 
	  result[list[i]] = list[i + 1];
	result[list[num_pts-1]] = list[0];

	free(list);
      }
    }
    break;
    case _INP_SEQ : {
      on_one_thread {
	/* Algorithm, courtesy of Chris Hurlburt, to convert  */
	/* the permuted vertex list into a valid successor array. */
	for (i = 0; i < num_pts-1 ; i++)
	  result[i] = i + 1;
	result[num_pts-1] = 0;
      }
    }
    break;
    case _INP_CYCLIC02 : {
      on_one_thread {
	for (i=0; i < num_pts - 2 ; i++)
	  result[i] = i + 2;
	for (i=0 ; i<2 ; i++)
	  result[num_pts-2+i] = (i+1) % 2;
      }
    }
    break;
    case _INP_CYCLIC04 : {
      on_one_thread {
	for (i=0; i < num_pts - 4 ; i++)
	  result[i] = i + 4;
	for (i=0 ; i<4 ; i++)
	  result[num_pts-4+i] = (i+1) % 4;
      }
    }
    break;
    case _INP_CYCLIC08 : {
      on_one_thread {
	for (i=0; i < num_pts - 8 ; i++)
	  result[i] = i + 8;
	for (i=0 ; i<8 ; i++)
	  result[num_pts-8+i] = (i+1) % 8;
      }
    }
    break;
    case _INP_CYCLIC16 : {
      on_one_thread {
	for (i=0; i < num_pts - 16 ; i++)
	  result[i] = i + 16;
	for (i=0 ; i<16 ; i++)
	  result[num_pts-16+i] = (i+1) % 16;
      }
    }
    break;
    case _INP_CYCLIC32 : {
      on_one_thread {
	for (i=0; i < num_pts - 32 ; i++)
	  result[i] = i + 32;
	for (i=0 ; i<32 ; i++)
	  result[num_pts-32+i] = (i+1) % 32;
      }
    }
    break;
    case _INP_CYCLIC64 : {
      on_one_thread {
	for (i=0; i < num_pts - 64 ; i++)
	  result[i] = i + 64;
	for (i=0 ; i<64 ; i++)
	  result[num_pts-64+i] = (i+1) % 64;
      }
    }
    break;
    default: fprintf(stderr,"ERROR: No input selected\n"); exit(-1);
    }
    
    /* Return S vector */
    return(result);
}


/* ---------------------------------------------------------------------------- */
/* Script function - run the algorithm at a given size for a given
   number of iterations, and return results & stats for same.
*/
void run_the_code(const int num_pts, const int num_iterations, 
		  double *min_time,
		  double *max_time,
		  double *avg_time,
		  double *sigma,
		  const int input,
		  THREADED)
{
    int loop;
    int curidx;
    double start_time, run_time, ttl_time;
    int *data_set = NULL;
    double *run_times = NULL;
    double dtmp, sum;

    ttl_time = 0.0;

    /* Set min, max to values such that will be replaced on the first run */
    *min_time = (double) num_pts * 10.0;
    *max_time = -0.01;

    /* Allocate memory to hold runtimes - don't know how much */
    /* we will need until runtime. */
    on_one_thread
    {
	run_times = (double *) malloc(sizeof(double) * num_iterations);
	assert(run_times != NULL);
    }

    for(loop = 0; loop < num_iterations; loop++)
    {

      /* Generating the data set is not counted in the runtime. */
      data_set = generate_array(num_pts, input, TH);

      start_time = get_seconds();

      on_one_thread {
	curidx = data_set[0];

	while (curidx != 0) 
	  curidx = data_set[curidx];

      }
      
      run_time = get_seconds() - start_time;

      /* Save time for stats later */
      on_one_thread
	run_times[loop] = run_time;

      /* Update min, max, total times - do these on one thread. */
      on_one_thread {

	if(run_time < *min_time)
	  *min_time = run_time;
	    
	if(run_time > *max_time)
	  *max_time = run_time;
	    
	ttl_time += run_time;
      }

      node_free(data_set, TH);

    }

    /* Compute stats */
    on_one_thread
    {
      /* Mean */
	*avg_time = (double) (ttl_time / num_iterations);

	/* Compute std deviation */
	sum = 0.0;

	for(loop = 0; loop < num_iterations; loop++)
	{
	    dtmp = (run_times[loop] - *avg_time);
	    sum += (dtmp * dtmp);
	}

	/* Sigma undefined for N==1 */
	if(num_iterations > 1)
	    *sigma = sqrt(sum) * (1.0 / (num_iterations - 1.0));
	else
	    *sigma = *avg_time;

	free(run_times);
    }

    node_Barrier();
    
    return;
}


/* ---------------------------------------------------------------------------- */
/* Main - create data, run the algorithm, save the stats. */
void *symbreak(const int input, int cur_size, int num_iterations, int testrun, THREADED)
{
    double avgtime, mintime, maxtime, sigtime;

    run_the_code(cur_size,
		 num_iterations,
		 &mintime, &maxtime, &avgtime, &sigtime, input, TH);

    if (!testrun) {
      on_one_thread {
	fprintf(stdout,"%3s\t%12d\t%12.9f\t%12.9f\t%12.9f\t%12.9f\n",
		inputLabel(input),
		cur_size,
		mintime,
		maxtime,
		avgtime,
		sigtime);
	fflush(stdout);
      }
    }

    return;
}

void *SIMPLE_main(THREADED) {
    const int min_size = MIN_SIZE;
    const int max_size = MAX_SIZE;
    const int num_iterations = NUM_ITER;

    int cur_size;

    assert(THREADS==1);
    
    on_one_thread    {
      fprintf(stdout,"Memread test\n");
      fprintf(stdout,"Min size: %d Max size: %d Iterations per size: %d\n\n",
	      min_size, max_size, num_iterations);

	fprintf(stdout,"INP\tSize\tMinimum\tMaximum\tAverage\tSigma\n");
	fflush(stdout);
    }

    for (cur_size = min_size; cur_size <= max_size; cur_size *= 2) {
#if DO_R
      symbreak(_INP_RANDOM,   cur_size, num_iterations, 0, TH);
      symbreak(_INP_RANDOM,   cur_size, num_iterations, 1, TH);
#endif
#if DO_S
      symbreak(_INP_SEQ,      cur_size, num_iterations, 0, TH);
      symbreak(_INP_SEQ,      cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC02, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC02, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC04, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC04, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC08, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC08, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC16, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC16, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC32, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC32, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
      symbreak(_INP_CYCLIC64, cur_size, num_iterations, 0, TH);
      symbreak(_INP_CYCLIC64, cur_size, num_iterations, 1, TH);
#endif
    }
    
  return;
}
