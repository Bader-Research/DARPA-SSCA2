#include "alg_nqueens.h"
#include "umd.h"
#include <math.h>
#ifdef RRANDOM
#include "alg_random.h"
#endif

#define K_DEFAULT 3

#define TRUE  1
#define FALSE 0

#define TIMING 1
#define TIME_RANDOM 1

#define PRINT_SOLN 0
#define LOOP_NQUEENS 1

#define MASK_INT  1

#define GCNT 0
#if GCNT
static int gcnt;
#endif

#define bitT(a,b) ((a) & (1<<(b)))
#define toggle(a,b) ((a) ^= (1<<(b)))

#if MASK_INT
#define NMAX 32
static int right0[2*(NMAX-1)+1];
static int *right;
static int left[2*(NMAX-1)+1];
static int height[NMAX];
#else
#define NMAX 16
static int right, left, height;
#endif

static int soln[NMAX];
static int count;

typedef struct {
#if MASK_INT  
  int right0[2*(NMAX-1)+1];
  int *right;
  int left[2*(NMAX-1)+1];
  int height[NMAX];
#else
  int right;
  int left;
  int height;
#endif
  int count;
  int soln[NMAX];
} nq_t;

/********************** COMMON FUNCTIONS **********************/

static void
printsoln(int N, int *sln) {
  int i;
  fprintf(outfile,"%3d: ",count);
  for (i=0 ; i<N ; i++)
    fprintf(outfile,"%2d ",sln[i]);
  fprintf(outfile,"\n");
}


static void
printsoln_r(int N, nq_t *nqp) {
  int i;
  fprintf(outfile,"%3d: ",nqp->count);
  for (i=0 ; i<N ; i++)
    fprintf(outfile,"%2d ",nqp->soln[i]);
  fprintf(outfile,"\n");
}

static int
randset(int a, int b) {
  int r;
  double r1;

  r1  = (double)random();
  r1 /= (double)2147483647;
  r1 *= (double)(b-a);
  r = (int)r1 + a;

  return r;
}

static int
randset_r(int a, int b, THREADED) {
  int r;
  double r1;

#ifdef RRANDOM  
  r1  = (double)rrandom_th(TH);
#else
  r1  = (double)random();
#endif  
  r1 /= (double)2147483647;
  r1 *= (double)(b-a);
  r = (int)r1 + a;

  return r;
}

#if MASK_INT
static void
invalid(int *sln, int d, int N) {
  int i,
    col, row;

  for (i=0 ; i<d ; i++) {
    col = i;
    row = sln[col];
    right[row - col] = FALSE;
    left[row + col]  = FALSE;
    height[row]      = FALSE;
  }
}
#else
static void
invalid(int *sln, int d, int N) {
  right  = 0;
  left   = 0;
  height = 0;
}
#endif

#if MASK_INT
static void
invalid_r(nq_t *nqp, int d, int N) {
  int i,
    col, row;

  for (i=0 ; i<d ; i++) {
    col = i;
    row = nqp->soln[col];
    nqp->right[row - col] = FALSE;
    nqp->left[row + col]  = FALSE;
    nqp->height[row]      = FALSE;
  }
}
#else
static void
invalid_r(nq_t *nqp, int d, int N) {
  nqp->right  = 0;
  nqp->left   = 0;
  nqp->height = 0;
}
#endif

static int
valid(int *sln, int d, int N) {
  int i,
    col, row,
    rx, lx;

  for (i=0 ; i<d ; i++) {
    col = i;
    row = sln[col];
#if MASK_INT
    rx = row-col;
    lx = row+col;
    if (!(right[rx] || left[lx] || height[row])) {
      right[rx]    = TRUE;
      left[lx]     = TRUE;
      height[row]  = TRUE;
    }
#else
    rx = row-col+N-1;
    lx = row+col;
    if (!(bitT(right,rx) || bitT(left,lx) || bitT(height,row))) {
      toggle(right,rx);
      toggle(left,lx);
      toggle(height,row);
    }
#endif
    else {
      invalid(sln,i,N);
      return(0);
    }
  }
  return(1);
}


static int
valid_r(nq_t *nqp, int d, int N) {
  int i,
    col, row,
    rx, lx;

  for (i=0 ; i<d ; i++) {
    col = i;
    row = nqp->soln[col];
#if MASK_INT
    rx = row-col;
    lx = row+col;
    if (!(nqp->right[rx] || nqp->left[lx] || nqp->height[row])) {
      nqp->right[rx]    = TRUE;
      nqp->left[lx]     = TRUE;
      nqp->height[row]  = TRUE;
    }
#else
    rx = row-col+N-1;
    lx = row+col;
    if (!(bitT(nqp->right,rx) || bitT(nqp->left,lx) || bitT(nqp->height,row))) {
      toggle(nqp->right,rx);
      toggle(nqp->left,lx);
      toggle(nqp->height,row);
    }
#endif
    else {
      invalid_r(nqp,i,N);
      return(0);
    }
  }
  return(1);
}



t(a,b,c){int d=0,e=a&~b&~c,f=1;if(a)for(f=0;d=(e-=d)&-e;f+=t(a-d,(b+d)*2,(c+d)/2));return f;}

/********************** RECURSIVE VERSIONS ********************/

static void
gensoln_R(int N, int col) {
  int row;
  register int rx, lx;

#if GCNT
  gcnt++;
#endif
  
  if (col==N) {
    count++;
#if PRINT_SOLN    
    printsoln(N, soln);
#endif
  }
  else {
    for (row=0 ; row<N ; row++) {
#if MASK_INT
      rx = row-col;
      lx = row+col;
      if (!(right[rx] || left[lx] || height[row])) {
	right[rx]        = TRUE;
	left[lx]         = TRUE;
	height[row]      = TRUE;
	soln[col]        = row;
	gensoln_R(N, col + 1);
	right[rx]        = FALSE;
	left[lx]         = FALSE;
	height[row]      = FALSE;
      }
#else
      rx = row-col+N-1;
      lx = row+col;
      if (!(bitT(right,rx) || bitT(left,lx) || bitT(height,row))) {
	toggle(right, rx);
	toggle(left,  lx);
	toggle(height,row);
	soln[col]        = row;
	gensoln_R(N, col + 1);
	toggle(right, rx);
	toggle(left,  lx);
	toggle(height,row);
      }
#endif
    }
  }
}

static void
gensoln_R_r(int N, int col, nq_t *nqp) {
  int row;
  register int rx, lx;

  if (col==N) {
    (nqp->count)++;
#if PRINT_SOLN    
    printsoln_r(N, nqp);
#endif
  }
  else {
    for (row=0 ; row<N ; row++) {
#if MASK_INT
      rx = row-col;
      lx = row+col;
      if (!(nqp->right[rx] || nqp->left[lx] || nqp->height[row])) {
	nqp->right[rx]        = TRUE;
	nqp->left[lx]         = TRUE;
	nqp->height[row]      = TRUE;
	nqp->soln[col]        = row;
	gensoln_R_r(N, col + 1, nqp);
	nqp->right[rx]        = FALSE;
	nqp->left[lx]         = FALSE;
	nqp->height[row]      = FALSE;
      }
#else
      rx = row-col+N-1;
      lx = row+col;
      if (!(bitT(nqp->right,rx) || bitT(nqp->left,lx) || bitT(nqp->height,row))) {
	toggle(nqp->right, rx);
	toggle(nqp->left,  lx);
	toggle(nqp->height,row);
	nqp->soln[col]        = row;
	gensoln_R_r(N, col + 1, nqp);
	toggle(nqp->right, rx);
	toggle(nqp->left,  lx);
	toggle(nqp->height,row);
      }
#endif
    }
  }
}

int nqueens_R_orig(int N) {
  int i, loop;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

#if GCNT
  gcnt=0;
#endif
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    gensoln_R(N, 0);
  }
  
  secs = get_seconds() - secs;
  
#if TIMING
  /*  fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0)); */
  fprintf(outfile,"\nN: %3d  queensR: %12d  Time: %9.6f\n",
	  N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif

#if GCNT
  fprintf(outfile,"nqueensR (N: %3d) gcnt: %d\n",N,gcnt);
  fflush(outfile);
#endif

  return (count);
}


int nqueens_R(int N, int k) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

  if (k < 0)
    k = K_DEFAULT;

#if GCNT
  gcnt=0;
#endif
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    for (i=0 ; i<Nk ; i++) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_R(N, k);
	invalid(soln,k,N);
      }
    }
  }
  
  secs = get_seconds() - secs;
  
#if TIMING
  /*  fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0)); */
  fprintf(outfile,"\nN: %3d  queensR: %12d  Time: %9.6f\n",
	  N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
#if GCNT
  fprintf(outfile,"nqueensR (N: %3d  K: %2d) gcnt: %d  N^k: %d\n",N,k,gcnt,Nk);
  fflush(outfile);
#endif
  
  return (count);
}

int all_nqueens_R_block(int N) {
  int i, loop;
  int s, mymin, mymax;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    s = N/NODES;
    if (((double)N / (double)NODES) > (double)s)
      s++;
    mymin = s*MYNODE;
    mymax = s*(MYNODE+1);
    if (mymax > N) mymax = N; 

    /* gensoln_R(N, 0, mymin, mymax); */

    for (i=mymin ; i<mymax ; i++) {
      soln[0] = i;
      if (valid(soln,1,N)) { 
	gensoln_R(N, 1);
	invalid(soln,1,N);
      }
    }
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}


int all_nqueens_R_random(int N, int k) {
  int i, loop;
  int s;
  int Nk, j;
  int Ni[NMAX+1];
  int *a, tmp;
  int *b;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

  if (k < 0)
    k = K_DEFAULT;
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    on_one_node {
      a = (int *)malloc(Nk*sizeof(int));
      assert_malloc(a);
  
      for (i=0; i < Nk; i++)
	a[i] = i;
      for (i=0; i < Nk; i++) {
	j = randset(i,Nk);
	tmp  = a[i];
	a[i] = a[j];
	a[j] = tmp;
      }
    }

    s = Nk/NODES;
    b = (int *)malloc(s*sizeof(int));
    assert_malloc(b);
    
    UMD_Scatter_i(a, s, b);
    
    for (i=0 ; i<s ; i++) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (b[i] % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_R(N, k);
	invalid(soln,k,N);
      }
    }

    /* Search remaining nodes */
    on_node(0) {
      for (i=s*NODES ; i<Nk ; i++) {
	for (j=0 ; j<k ;  j++) 
	  soln[j] = (a[i] % Ni[j+1]) / Ni[j];
	if (valid(soln,k,N)) {
	  gensoln_R(N, k);
	  invalid(soln,k,N);
	}
      }
    }
    
    free(b);
    on_one_node free(a);
    
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}

int all_nqueens_R_cyclic(int N, int k) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];
    
    for (i=MYNODE ; i<Nk ; i+=NODES) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_R(N, k);
	invalid(soln,k,N);
      }
    }

    
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}

int all_nqueens_R_random_r(int N, int k, THREADED) {
  int i, loop;
  int s;
  int Nk, j;
  int Ni[NMAX+1];

  int *a;
  int *b;
  double secs;
#if TIME_RANDOM
  double rsecs;
#endif
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    on_one {
      int tmp;
      a = (int *)malloc(Nk*sizeof(int));
      assert_malloc(a);

#if TIME_RANDOM
      rsecs = get_seconds();
#endif
      for (i=0 ; i<Nk ; i++)
	a[i] = i;
      for (i=0; i < Nk; i++) {
	j = randset_r(i,Nk, TH);
	tmp  = a[i];
	a[i] = a[j];
	a[j] = tmp;
      }
#if TIME_RANDOM
      rsecs = get_seconds() - rsecs;
      fprintf(outfile,"PE%3d(%3d): nqueens n: %2d k: %2d Nk: %6d randomization: %9.6f\n",
	      MYNODE,MYTHREAD,N,k,Nk,rsecs);
      fflush(outfile);
#endif
    }

    s = Nk/NODES + 1;

    b = (int *)node_malloc(s*sizeof(int), TH);
    a = node_Bcast_ip(a, TH);
    
    s = all_ScatterX_i(a, Nk, b, TH);

    all_Barrier(TH);
    
    pardo(i,0,s,1) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (b[i] % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_R_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }

    node_free(b, TH);

    on_one
      free(a);
    
  }

  /*  node_Barrier(); */
  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}

int all_nqueens_R_cyclic_r(int N, int k, THREADED) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    all_pardo_cyclic(i,0,Nk) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_R_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }
    
  }

  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}

int all_nqueens_R_block_r(int N, int k, THREADED) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    all_pardo_block(i,0,Nk) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_R_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }
    
  }

  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}

/*********** NON-RECURSIVE VERSIONS ****************************/

#if MASK_INT
static void
gensoln_NR(int N, int col0) {
  int row, col, last_row=-1, last_col=-1;
  register int rx, lx;

  col       =  col0;
  soln[col] = -1;

  while (col >= col0) {
    if (last_col == col) {
      right[last_row - last_col] = FALSE;
      left[last_row + last_col]  = FALSE;
      height[last_row]           = FALSE;
      last_col--;
      last_row = soln[last_col];
    }
    row = soln[col];
    if (row >= N-1) {
      col--;
    }
    else {
      row = (++soln[col]);
      rx = row-col;
      lx = row+col;
      if (!(right[rx] || left[lx] || height[row])) {
	last_row = row;
	last_col = col;
	right[rx]   = TRUE;
	left[lx]    = TRUE;
	height[row] = TRUE;
	if (col == N-1) {
	  count++;
#if PRINT_SOLN    
	  printsoln(soln, N);
#endif
	}
	else {
	  col++;
	  soln[col] = -1;
	}
      }
    }
  }
}
#else
static void
gensoln_NR(int N, int col0) {
  int row, col, last_row=-1, last_col=-1;
  register int rx, lx;

  col       =  col0;
  soln[col] = -1;

  while (col >= col0) {
    if (last_col == col) {
      toggle(right, last_row-last_col+N-1);
      toggle(left,  last_row+last_col);
      toggle(height,last_row);
      last_col--;
      last_row = soln[last_col];
    }
    row = soln[col];
    if (row >= N-1) {
      col--;
    }
    else {
      row = (++soln[col]);
      rx = row-col+N-1;
      lx = row+col;
      if (!(bitT(right,rx) || bitT(left,lx) || bitT(height,row))) {
	last_row = row;
	last_col = col;
	toggle(right, rx);
	toggle(left,  lx);
	toggle(height,row);
	if (col == N-1) {
	  count++;
#if PRINT_SOLN    
	  printsoln(soln, N);
#endif
	}
	else {
	  col++;
	  soln[col] = -1;
	}
      }
    }
  }
}
#endif

#if MASK_INT
static void
gensoln_NR_r(int N, int col0, nq_t *nqp) {
  int row, col, last_row=-1, last_col=-1;
  register int rx, lx;

  col            =  col0;
  nqp->soln[col] = -1;

  while (col >= col0) {
    if (last_col == col) {
      nqp->right[last_row - last_col] = FALSE;
      nqp->left[last_row + last_col]  = FALSE;
      nqp->height[last_row]           = FALSE;
      last_col--;
      last_row = nqp->soln[last_col];
    }
    row = nqp->soln[col];
    if (row >= N-1) {
      col--;
    }
    else {
      row = (++(nqp->soln[col]));
      rx = row-col;
      lx = row+col;
      if (!(nqp->right[rx] || nqp->left[lx] || nqp->height[row])) {
	last_row = row;
	last_col = col;
	nqp->right[rx]   = TRUE;
	nqp->left[lx]    = TRUE;
	nqp->height[row] = TRUE;
	if (col == N-1) {
	  nqp->count++;
#if PRINT_SOLN    
	  printsoln_r(N, nqp);
#endif
	}
	else {
	  col++;
	  nqp->soln[col] = -1;
	}
      }
    }
  }
}
#else
static void
gensoln_NR_r(int N, int col0, nq_t *nqp) {
  int row, col, last_row=-1, last_col=-1;
  register int rx, lx;

  col       =  col0;
  nqp->soln[col] = -1;

  while (col >= col0) {
    if (last_col == col) {
      toggle(nqp->right, last_row-last_col+N-1);
      toggle(nqp->left,  last_row+last_col);
      toggle(nqp->height,last_row);
      last_col--;
      last_row = nqp->soln[last_col];
    }
    row = nqp->soln[col];
    if (row >= N-1) {
      col--;
    }
    else {
      row = (++(nqp->soln[col]));
      rx = row-col+N-1;
      lx = row+col;
      if (!(bitT(nqp->right,rx) || bitT(nqp->left,lx) || bitT(nqp->height,row))) {
	last_row = row;
	last_col = col;
	toggle(nqp->right, rx);
	toggle(nqp->left,  lx);
	toggle(nqp->height,row);
	if (col == N-1) {
	  nqp->count++;
#if PRINT_SOLN    
	  printsoln_r(N, nqp);
#endif
	}
	else {
	  col++;
	  nqp->soln[col] = -1;
	}
      }
    }
  }
}
#endif

int nqueens_NR_orig(int N) {
  int i, loop;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

#if GCNT
  gcnt=0;
#endif
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    gensoln_NR(N, 0);
  }
  
  secs = get_seconds() - secs;
  
#if TIMING
  /*  fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0)); */
  fprintf(outfile,"\nN: %3d  queensNR: %12d  Time: %9.6f\n",
	  N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
#if GCNT
  fprintf(outfile,"nqueensNR (N: %3d) gcnt: %d\n",N,gcnt);
  fflush(outfile);
#endif

  return (count);
}


int nqueens_NR(int N, int k) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

  if (k < 0)
    k = K_DEFAULT;

#if GCNT
  gcnt=0;
#endif
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    for (i=0 ; i<Nk ; i++) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_NR(N, k);
	invalid(soln,k,N);
      }
    }
  }
  
  secs = get_seconds() - secs;
  
#if TIMING
  /*  fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0)); */
  fprintf(outfile,"\nN: %3d  queensNR: %12d  Time: %9.6f\n",
	  N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
#if GCNT
  fprintf(outfile,"nqueensNR (N: %3d  K: %2d) gcnt: %d  N^k: %d\n",
	  N,k,gcnt,Nk);
  fflush(outfile);
#endif
  
  return (count);
}

int all_nqueens_NR_block(int N) {
  int i, loop;
  int s, mymin, mymax;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    s = N/NODES;
    if (((double)N / (double)NODES) > (double)s)
      s++;
    mymin = s*MYNODE;
    mymax = s*(MYNODE+1);
    if (mymax > N) mymax = N; 

    /*    gensoln_NR(N, 0, mymin, mymax); */ 

    for (i=mymin ; i<mymax ; i++) {
      soln[0] = i;
      if (valid(soln,1,N)) { 
	gensoln_NR(N, 1);
	invalid(soln,1,N);
      }
    }
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}


int all_nqueens_NR_random(int N, int k) {
  int i, loop;
  int s;
  int Nk, j;
  int Ni[NMAX+1];
  int *a, tmp;
  int *b;
  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }

  if (k < 0)
    k = K_DEFAULT;
  
  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    on_one_node {
      a = (int *)malloc(Nk*sizeof(int));
      assert_malloc(a);
  
      for (i=0; i < Nk; i++)
	a[i] = i;
      for (i=0; i < Nk; i++) {
	j = randset(i,Nk);
	tmp  = a[i];
	a[i] = a[j];
	a[j] = tmp;
      }
    }

    s = Nk/NODES;
    b = (int *)malloc(s*sizeof(int));
    assert_malloc(b);
    
    UMD_Scatter_i(a, s, b);
    
    for (i=0 ; i<s ; i++) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (b[i] % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_NR(N, k);
	invalid(soln,k,N);
      }
    }

    /* Search remaining nodes */
    on_node(0) {
      for (i=s*NODES ; i<Nk ; i++) {
	for (j=0 ; j<k ;  j++) 
	  soln[j] = (a[i] % Ni[j+1]) / Ni[j];
	if (valid(soln,k,N)) {
	  gensoln_NR(N, k);
	  invalid(soln,k,N);
	}
      }
    }
    
    free(b);
    on_one_node free(a);
    
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}

int all_nqueens_NR_cyclic(int N, int k) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    right = right0 + N-1;
    for (i=1-N ; i<N ; i++)
      right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      height[i] = FALSE;
#else
    right  = 0;
    left   = 0;
    height = 0;
#endif

    count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];
    
    for (i=MYNODE ; i<Nk ; i+=NODES) {
      for (j=0 ; j<k ;  j++) 
	soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid(soln,k,N)) {
	gensoln_NR(N, k);
	invalid(soln,k,N);
      }
    }

    
  }
    
  secs = get_seconds() - secs;

#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d: N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,N,count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  count = UMD_Reduce_i(count,SUM);

  return (count);
}

int all_nqueens_NR_random_r(int N, int k, THREADED) {
  int i, loop;
  int s;
  int Nk, j;
  int Ni[NMAX+1];

  int *a;
  int *b;
  double secs;
#if TIME_RANDOM
  double rsecs;
#endif
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    on_one {
      int tmp;
      a = (int *)malloc(Nk*sizeof(int));
      assert_malloc(a);

#if TIME_RANDOM
      rsecs = get_seconds();
#endif
      for (i=0 ; i<Nk ; i++)
	a[i] = i;
      for (i=0; i < Nk; i++) {
	j = randset_r(i,Nk, TH);
	tmp  = a[i];
	a[i] = a[j];
	a[j] = tmp;
      }
#if TIME_RANDOM
      rsecs = get_seconds() - rsecs;
      fprintf(outfile,"PE%3d(%3d): nqueens n: %2d k: %2d Nk: %6d randomization: %9.6f\n",
	      MYNODE,MYTHREAD,N,k,Nk,rsecs);
      fflush(outfile);
#endif
    }

    s = Nk/NODES + 1;

    b = (int *)node_malloc(s*sizeof(int), TH);
    a = node_Bcast_ip(a, TH);
    
    s = all_ScatterX_i(a, Nk, b, TH);

    all_Barrier(TH);
    
    pardo(i,0,s,1) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (b[i] % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_NR_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }

    node_free(b, TH);
    on_one
      free(a);
    
  }

  /*  node_Barrier(); */
  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}

int all_nqueens_NR_cyclic_r(int N, int k, THREADED) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    all_pardo_cyclic(i,0,Nk) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_NR_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }
    
  }

  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}

int all_nqueens_NR_block_r(int N, int k, THREADED) {
  int i, loop;
  int Nk, j;
  int Ni[NMAX+1];

  double secs;
  nq_t nq;

  if (N > NMAX) {
    fprintf(stderr,"N (%d) greater than NMAX (%d)\n",N,NMAX);
    exit(-1);
  }
  
  if (k < 0)
    k = K_DEFAULT;

  secs = get_seconds();

  for (loop=0 ; loop<LOOP_NQUEENS ; loop++) {
    
#if MASK_INT
    nq.right = nq.right0 + N-1;
    for (i=1-N ; i<N ; i++)
      nq.right[i] = FALSE;
    for (i=0 ; i<(2*(N-1) + 1) ; i++)
      nq.left[i] = FALSE;
    for (i=0 ; i<N ; i++)
      nq.height[i] = FALSE;
#else
    nq.right  = 0;
    nq.left   = 0;
    nq.height = 0;
#endif

    nq.count = 0;

    Ni[0] = 1;
    for (i=1 ; i<=k ; i++)
      Ni[i] = Ni[i-1] * N;
    Nk = Ni[k];

    all_pardo_block(i,0,Nk) {
      for (j=0 ; j<k ;  j++) 
	nq.soln[j] = (i % Ni[j+1]) / Ni[j];
      if (valid_r(&nq,k,N)) {
	gensoln_NR_r(N, k, &nq);
	invalid_r(&nq,k,N);
      }
    }
    
  }

  secs = get_seconds() - secs;
  
#if TIMING
  /*fprintf(outfile,"\n queens:%d\n",t(~(~0<<N),0,0));*/
  fprintf(outfile,"PE%3d(%3d): N: %3d  queensNR: %12d  Time: %9.6f\n",
	  MYNODE,MYTHREAD,N,nq.count,secs/(double)LOOP_NQUEENS);
  fflush(outfile);
#endif
  
  i = all_Reduce_i(nq.count,SUM,TH);

  return (i);
}
