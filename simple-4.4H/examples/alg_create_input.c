#include "alg_create_input.h"
#include "nas_r.h"

#define DEBUG  0
#define TIMING 1

void create_input_nas(int n, int *x) {

#if TIMING
  double secs, tsec;
  secs = get_seconds();
#endif

  create_seq( find_my_seed( MYNODE,
			    NODES,
			    (n >> 2),
			    _NAS_SEED,    /* Random number gen seed */
			    _NAS_MULT),   /* Random number gen mult */
	      _NAS_MULT,                  /* Random number gen mult */
	      n/NODES,
	      x);   

#if TIMING
  secs = get_seconds() - secs;

  MPI_Reduce(&secs, &tsec, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
  on_one_node {
    fprintf(outfile,"PE%3d: n: %12d  FTime: %9.6lf\n",MYNODE,n,tsec);
    fflush(outfile);
  }
#endif
}

void create_input_nas_simple(int n, int *x, THREADED) {
  register int tsize, mynum, thtot;
#if TIMING
  double secs, tsec;
  secs = get_seconds();
#endif

  tsize = (n/NODES) / THREADS;
  mynum = ID;
  thtot = TID;
  
  create_seq_simple( find_my_seed( mynum,
				   thtot,
				   (n >> 2),
				   _NAS_SEED,    /* Random number gen seed */
				   _NAS_MULT),   /* Random number gen mult */
		     _NAS_MULT,                  /* Random number gen mult */
		     tsize,
		     x+(tsize*MYTHREAD),
		     TH);   

#if TIMING
  secs = get_seconds() - secs;
  tsec = all_Reduce_d(secs, MAX, TH);
  on_one {
    fprintf(outfile,"PE%3d: n: %12d  STime: %9.6lf\n",MYNODE,n,tsec);
    fflush(outfile);
  }
#endif
}


void create_input_random(int myN, int *x) {

  register int i;

#if 1
  srandom(317*(MYNODE+17));
#endif
    
  for (i=0 ; i<myN ; i++)
    x[i] = random(); 
}

void create_input_random_simple_old(int myN, int *x, THREADED) {

  register int i;

#if 1
  srandom(317*(ID+17));
#endif

#if DEBUG
  fprintf(outfile,"PE%3d(%3d): create_input_random_simple(): %12d\n",
	  MYNODE,MYTHREAD,myN);
  fflush(outfile);
  all_Barrier(TH);
#endif
  
  for (i=0 ; i<myN ; i++)
    x[i] = random();

#if DEBUG
  all_Barrier(TH);
  on_one fprintf(outfile,"PE%3d(%3d): create_input_random_simple(): done\n",
	  MYNODE,MYTHREAD);
  fflush(outfile);
  all_Barrier(TH);
#endif
}


void create_input_random_simple(int myN, int *x, THREADED) {
  create_seq_random_simple( 317*(ID+17),
			    _NAS_MULT,   
			    myN,
			    x,
			    TH);   
}

