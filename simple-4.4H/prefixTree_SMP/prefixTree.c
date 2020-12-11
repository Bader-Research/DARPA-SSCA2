/* Tree-based Approach*/

#include "simple.h"

#define DEBUG 0

void input_gen (int **result,int k,int *n,THREADED) {
  
  int i;
  
  *n = (1<<k);

#if DEBUG
  printf("T%3d: input is %d\n",MYTHREAD,*n);
  fflush(outfile);
  node_Barrier();
#endif
  
  *result = (int *) node_malloc(*n*sizeof(int),TH);
  
  if (*result == NULL) {
    exit(-1);
  }
  on_one_thread {
    for(i=0 ; i<*n ; i++) {
      (*result)[i] = i+1;
    }
  }

#if DEBUG
  for (i=0 ; i<*n ; i++)
    fprintf(outfile,"T%3d: result[%5d]: %12d\n",MYTHREAD,i,(*result)[i]);
  fflush(outfile);
#endif
  
  node_Barrier();
  return;
}


void prefix_sums(int *result, int k, int n, THREADED) {
  int h;
  int p;
  int j;
  int i;
  int x;
  int y;

#if DEBUG
  on_one_thread {
    fprintf(outfile,"k: %d\n", k); 
  }
#endif

  node_Barrier();

  for (h=1; h<=k ; h++) { 
      
    p = (1<<(k-h));
    j = (1<<h); 

    pardo (x,1,p+1,1) {
      result[(j*x) -1] += result[((1<<(h-1))*((x<<1)-1)) -1];
    } 
       
    node_Barrier();
  }       

  for (h=k-2 ; h>=0 ; h--) {

    p = (1<<(k-h));
    j=(1<<h); 

    pardo (y,3,p+1,2) {
      result[(j*y) -1] += result[(j*(y-1)) -1];
    }
      
    node_Barrier();
  }
  
#if DEBUG
 on_one_thread {
    for (i=0 ; i<n; i++) {
      fprintf(outfile,"prefix_sum: %12d\n", result[i]);
    }
 }
#endif

 return;
} 

void timing(int *result,int k,int n,THREADED) {
  int i;
  double t;

  t = get_seconds();
  prefix_sums(result, k, n, TH);
  t = get_seconds() - t;

  t = node_Reduce_d(t, MAX, TH);
  on_one_thread {
    fprintf(outfile,"pTree  P: %3d  k: %12d n: %12d Time: %9.6f\n",THREADS, k, n, t);
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
  fprintf(outfile,"THARGC = %d\n",THARGC);
  for (i=0;i<THARGC;i++)
    fprintf(outfile,"THARGV[%d]=%s\n",i,THARGV[i]);
  node_Barrier();
#endif
 
  if (THARGC != 1) {
    fprintf(stderr,"ERROR: you must call the program with one argument\n");
    exit(-1);
  }
  
  k = atoi(THARGV[0]);
  
  input_gen(&result, k,&n,TH);
  
  timing(result, k,n,TH);
  
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

