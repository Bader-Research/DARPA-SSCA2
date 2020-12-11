/* Block Based Approach*/

#include "simple.h"


void input_gen(int **result,int k,int *n,THREADED) {
  int i;
  
  *n=(1<<k);

  *result=(int *)node_malloc(*n*sizeof(int),TH);
  
  if (*result == NULL) {
    exit(-1);
  }
  
  if (MYTHREAD == 0) {
    for (i=0 ; i<*n ; i++) 
      (*result)[i] = i+1;
  }

#if DEBUG
  for(i=0 ; i<*n ; i++) {
    fprintf(outfile,"T%3d: result[%5d]: %12d\n",MYTHREAD,i,(*result)[i]);
  }
  fflush(outfile);
#endif
  
}


void prefix_sums(int *result,int k,int n,THREADED) {
  int j;
  int i;
  int *p;
  int r;
  int start;
  int end;
  int add_value;

  r = n / THREADS;

#if DEBUG
  fprintf(outfile,"T%3d: value of r: %5d",MYTHREAD, r);
#endif

  p = (int *)node_malloc(NOSHARE(THREADS)*sizeof(int),TH);

#if DEBUG   
  if (p == NULL) {
    exit(-1);
  }

  if (r * (THREADS) != n) {
    fprintf(stderr,"error\n");
  }
#endif

  start =  MYTHREAD*r + 1;
  end   = (MYTHREAD+1)*r;
  
  for (j=start ; j<end ; j++) 
    result[j] += result[j-1];
  
  p[NOSHARE(MYTHREAD)] = result[end-1];

  node_Barrier();
  
  on_one_thread {
    for (j=1 ; j<THREADS ; j++)
      p[NOSHARE(j)] += p[NOSHARE(j-1)];
  }
    
  node_Barrier();

  if (MYTHREAD>0) {
    add_value=p[NOSHARE(MYTHREAD-1)];
    
    for (j=start-1 ; j<end ; j++)
      result[j] += add_value;
  }
  
  node_Barrier();

#if DEBUG
  on_one {
    for (i=0 ; i<n ; i++)
      fprintf(outfile,"prefix_sums are: %12d\n",result[i]);
  }
#endif

}


void timing(int *result,int k,int n,THREADED) {
  double t;

  t = get_seconds();
  prefix_sums(result,k,n,TH);
  t = get_seconds() - t;

  t = node_Reduce_d(t, MAX, TH);
  on_one_thread {
    fprintf(outfile,"pBlock P: %3d  k: %12d n: %12d Time: %9.6f\n",
	    THREADS, k, n, t);
    fflush(outfile);
  }

  return;
}

void *SIMPLE_main(THREADED) {

  int n;
  int k;
  int *result;
  int *resultseq;
  int i;

#if DEBUG
  fprintf(outfile,"T%3d: THARGC= %d\n",MYTHREAD,THARGC);
  for (i=0 ; i<THARGC ; i++)
    fprintf(outfile,"T%3d: THARGV[%d]=%s\n",MYTHREAD,i,THARGV[i]);
  
  node_Barrier();
#endif

  if (THARGC != 1) {
    fprintf(stderr,"ERROR:call the problem with one argument\n");
    exit(-1);
  }

  k = atoi(THARGV[0]);

  input_gen(&result,k,&n,TH);

  timing(result,k,n,TH);

  node_Barrier();

  input_gen(&resultseq, k,&n,TH);
    
  node_Barrier();
  
  on_one_thread {
    
    fprintf(outfile,"n: %d\n",n);
    fflush(outfile);
    for (i=1 ; i<n ; i++)
      resultseq[i] += resultseq[i-1];
    fprintf(outfile,"Checking the result:\n");
    fflush(outfile);
    for (i=0 ; i<n ; i++)
      if (result[i] != resultseq[i]) {
	fprintf(outfile,"ERROR: i = %12d  result[i]: %6d resultseq[i]: %6d\n",
		i, result[i], resultseq[i]);
      }
    fprintf(outfile,"Done Checking.\n");
    fflush(outfile);
  }
  
  node_Barrier();

  node_free(result, TH);
  node_free(resultseq, TH);

  node_Barrier();
  
  SIMPLE_done(TH);
}















