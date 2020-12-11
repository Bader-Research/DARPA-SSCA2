#include <stdio.h>

#include "alg_radix_smp.h"

#define DEBUG  0
#define TIMING 0

#define bits(x,k,j) ((x>>k) & ~(~0<<j))

/****************************************************/
void all_countsort_node(int q,
			int *lKey,
			int *lSorted,
			int R,
			int bitOff, int m,
			THREADED)
/****************************************************/
/* R (range)      must be a multiple of NODES */
/* q (elems/proc) must be a multiple of NODES */
{
    register int
	j,
	k,
        last, temp,
	offset;
    
    int *myHisto,
        *mhp,
        *mps,
        *psHisto,
        *allHisto;

#if TIMING
    INIT_STEP_TH();
#endif

    myHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);
    psHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);

#if TIMING
    START_STEP_TH();
#endif

    mhp = myHisto + MYTHREAD*R;

    for (k=0 ; k<R ; k++)
      mhp[k] = 0;
    
#if TIMING
    END_STEP1_TH("N1init ",q,bitOff);
    START_STEP_TH();
#endif

    pardo(k, 0, q, 1)
      mhp[bits(lKey[k],bitOff,m)]++;

    node_Barrier();

#if TIMING
    END_STEP1_TH("N2histo",q,bitOff);
    START_STEP_TH();
#endif

    pardo(k, 0, R, 1) {
      last = psHisto[k] = myHisto[k];
      for (j=1 ; j<THREADS ; j++) {
	temp = psHisto[j*R + k] = last + myHisto[j*R +  k];
	last = temp;
      }
    }

    allHisto = psHisto+(THREADS-1)*R;

    node_Barrier();
    
#if TIMING
    END_STEP1_TH("N3hcomb",q,bitOff);
    START_STEP_TH();
#endif

    node_Barrier();

    offset = 0;

    mps = psHisto + (MYTHREAD*R);
    for (k=0 ; k<R ; k++) {
      mhp[k]  = (mps[k] - mhp[k]) + offset;
      offset += allHisto[k];
    }

    node_Barrier();
    
#if TIMING
    END_STEP1_TH("N4rankp",q,bitOff);
    START_STEP_TH();
#endif

    pardo(k, 0, q, 1) {
      j = bits(lKey[k],bitOff,m);
      lSorted[mhp[j]] = lKey[k];
      mhp[j]++;
    }

    node_Barrier();

#if TIMING
    END_STEP1_TH("N5place",q,bitOff);
    REPORT_STEP1_TH("N9TOTAL",q,bitOff);
#endif

    node_free(psHisto, TH);
    node_free(myHisto, TH);
}

/****************************************************/
void all_countsort_node_aux(int q,
			    int *lKey,
			    int *lSorted, int* auxKey, int* auxSorted,
			    int R,
			    int bitOff, int m,
			    THREADED)
/****************************************************/
/* R (range)      must be a multiple of SMPS */
/* q (elems/proc) must be a multiple of SMPS */
{
    register int
	j,
	k,
        last, temp,
	offset;
    
    int *myHisto,
        *mhp,
        *mps,
        *psHisto,
        *allHisto;

    myHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);
    psHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);

    mhp = myHisto + MYTHREAD*R;

    for (k=0 ; k<R ; k++)
      mhp[k] = 0;
    
    pardo(k, 0, q, 1)
      mhp[bits(lKey[k],bitOff,m)]++;

    node_Barrier();

    pardo(k, 0, R, 1) {
      last = psHisto[k] = myHisto[k];
      for (j=1 ; j<THREADS ; j++) {
	temp = psHisto[j*R + k] = last + myHisto[j*R +  k];
	last = temp;
      }
    }

    allHisto = psHisto+(THREADS-1)*R;

    node_Barrier();
    
    offset = 0;

    mps = psHisto + (MYTHREAD*R);
    for (k=0 ; k<R ; k++) {
      mhp[k]  = (mps[k] - mhp[k]) + offset;
      offset += allHisto[k];
    }

    node_Barrier();
    
    pardo(k, 0, q, 1) {
      j = bits(lKey[k],bitOff,m);
      lSorted[mhp[j]] = lKey[k];
      auxSorted[mhp[j]] = auxKey[k];
      mhp[j]++;
    }

    node_Barrier();

    node_free(psHisto, TH);
    node_free(myHisto, TH);
}



/****************************************************/
void all_radixsort_node_s3(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED)
/****************************************************/
{
  int *lTemp;

    lTemp = (int *)node_malloc(q*sizeof(int), TH);
			
    all_countsort_node(q, lKeys,   lSorted, (1<<11),  0, 11, TH);
    all_countsort_node(q, lSorted, lTemp,   (1<<11), 11, 11, TH);
    all_countsort_node(q, lTemp,   lSorted, (1<<10), 22, 10, TH);

    node_free(lTemp, TH);
}

/****************************************************/
void all_radixsort_node_s2(int q,
			   int *lKeys,
			   int *lSorted,
			   THREADED)
/****************************************************/
{
  int *lTemp;

    lTemp = (int *)node_malloc(q*sizeof(int), TH);
			
    all_countsort_node(q, lKeys,   lTemp,   (1<<16),  0, 16, TH);
    all_countsort_node(q, lTemp,   lSorted, (1<<16), 16, 16, TH);

    node_free(lTemp, TH);
}

/****************************************************/
void all_radixsort_node_aux_s3(int q,
			       int *lKeys,
			       int *lSorted, int* auxKey, int* auxSorted,
			       THREADED)
/****************************************************/
{
  int *lTemp, *lTemp2;

    lTemp = (int *)node_malloc(q*sizeof(int), TH);
    lTemp2 = (int *)node_malloc(q*sizeof(int), TH);
		
    all_countsort_node_aux(q, lKeys, lSorted, auxKey, auxSorted, (1<<11),  0, 11, TH);
    all_countsort_node_aux(q, lSorted, lTemp, auxSorted, lTemp2, (1<<11), 11, 11, TH);
    all_countsort_node_aux(q, lTemp, lSorted, lTemp2, auxSorted, (1<<10), 22, 10, TH);

    node_free(lTemp, TH);
    node_free(lTemp2, TH);
}











