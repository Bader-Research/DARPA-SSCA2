#include "alg_histo.h"
#include "umd.h"
#include <math.h>

#define DSIZE 256

#define DEBUG 0

int all_histo(int myL, int *A, int R, int *histo) {
  int i, j;
  int *Ap, *As;
  int s;
  int *htran;

  htran = (int *)malloc(R*sizeof(int));
  assert_malloc(htran);

  for (i=0 ; i<R ; i++)
    histo[i] = 0;

  As = (Ap=A) + myL;
  while (Ap < As) {
    histo[*Ap]++;
    Ap++;
  }

  s = R/NODES;
  UMD_Alltoall_i(histo, s, htran);

  for (i=0 ; i<s ; i++) 
    for (j=1 ; j<NODES ; j++) 
      htran[i] += htran[j*s + i];

  /* Optional? */
  UMD_Gather_i(htran, s, histo);

#if DEBUG
  on_one_node {
    for (i=0 ; i<R ; i++)
      fprintf(outfile,"histo[%3d]: %12d\n",i,histo[i]);
    fflush(outfile);
  }
#endif
  free(htran);
  return 0;
}


int all_histo_r(int myL, int *A, int R, int *histo, THREADED) {
  int i, j;
  int s;
  int *allhisto;
  int *myhisto;
  int *htran;

  allhisto = (int *)node_malloc(R*THREADS*sizeof(int), TH);
  htran = (int *)node_malloc(R*sizeof(int), TH);

  myhisto   = allhisto + (R * MYTHREAD);
  
  for (i=0 ; i<R ; i++)
    myhisto[i] = 0;

  pardo(i, 0, myL, 1) 
    myhisto[A[i]]++;

  node_Barrier();

  pardo(i, 0, R, 1) 
    for (j=1 ; j<THREADS ; j++)
      allhisto[i] += allhisto[j*DSIZE +  i];

  node_Barrier();
  
  /* We now have one array / node */
  
  s = R/NODES;
  all_Alltoall_i(allhisto, s, htran, TH);

  node_Barrier();

  pardo(i, 0, s, 1)
    for (j=1 ; j<NODES ; j++) 
      htran[i] += htran[j*s + i];

  node_Barrier();
  
  /* Optional? */
  all_Gather_i(htran, s, histo, TH);

#if DEBUG
  on_one {
    for (i=0 ; i<R ; i++)
      fprintf(outfile,"r histo[%3d]: %12d\n",i,histo[i]);
    fflush(outfile);
  }
#endif

  node_free(htran, TH);
  node_free(allhisto, TH);

  return 0;
}

int all_sel_sum(int myL, int *A) {
  int i, j;
  int *Ap, *As;
  int s;
  long D[DSIZE];
  long Dt[DSIZE];

  for (i=0 ; i<DSIZE ; i++)
    D[i] = 0;

  As = (Ap=A) + myL;
  while (Ap < As) {
    D[(*Ap ^ (*Ap>>16)) & 255] += *Ap;
    Ap++;
  }

  s = DSIZE/NODES;
  UMD_Alltoall_l(D, s, Dt);

  for (i=0 ; i<s ; i++) 
    for (j=1 ; j<NODES ; j++) 
      Dt[i] += Dt[j*s + i];

  /* Optional? */
  UMD_Gather_l(Dt, s, D);

#if DEBUG
  on_one_node {
    for (i=0 ; i<DSIZE ; i++)
      fprintf(outfile,"D[%3d]: %12ld\n",i,D[i]);
    fflush(outfile);
  }
#endif
  return 0;
}


int all_sel_sum_r(int myL, int *A, THREADED) {
  int i, j;
  int s;
  long *D, *myD;
  long *Dt;

  D  = (long *)node_malloc(DSIZE*THREADS*sizeof(long), TH);
  Dt = (long *)node_malloc(DSIZE*sizeof(long), TH);

  myD = D + (DSIZE * MYTHREAD);
  
  for (i=0 ; i<DSIZE ; i++)
    myD[i] = 0;

  pardo(i, 0, myL, 1) 
    myD[(A[i] ^ (A[i]>>16)) & 255] += A[i];

  node_Barrier();

  pardo(i, 0, DSIZE, 1) 
    for (j=1 ; j<THREADS ; j++)
      D[i] += D[j*DSIZE +  i];

  node_Barrier();
  
  /* We now have one array / node */
  
  s = DSIZE/NODES;
  all_Alltoall_l(D, s, Dt, TH);

  node_Barrier();

  pardo(i, 0, s, 1)
    for (j=1 ; j<NODES ; j++) 
      Dt[i] += Dt[j*s + i];

  node_Barrier();
  
  /* Optional? */
  all_Gather_l(Dt, s, D, TH);

#if DEBUG
  on_one {
    for (i=0 ; i<DSIZE ; i++)
      fprintf(outfile,"r D[%3d]: %12ld\n",i,D[i]);
    fflush(outfile);
  }
#endif

  node_free(Dt, TH);
  node_free(D, TH);

  return 0;
}

