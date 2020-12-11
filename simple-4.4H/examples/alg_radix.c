#include <stdio.h>

#include "alg_radix.h"
#include "alg_route.h"
#include "alg_sort.h"

#define DEBUG  0
#define TIMING 0

/****************************************************/
void all_countsort(int q,
		   int *lKey,
		   int *lSorted,
		   int R,
		   int bitOff, int m)
/****************************************************/
/* R (range)      must be a multiple of NODES */
/* q (elems/proc) must be a multiple of NODES */
{

    register int
	j,
	k,
	s,
	offset;
    
    int *lAddr,
	*lIndex,
	*lScanTran;

    intpair_t
	*lIntScanTran,
	*lScans;

#if TIMING
    INIT_STEP();
#endif
    
    lAddr     = (int *)malloc(q*sizeof(int));
    assert_malloc(lAddr);
    lIndex    = (int *)malloc(R*sizeof(int));
    assert_malloc(lIndex);
    lScanTran = (int *)malloc(R*sizeof(int));
    assert_malloc(lScanTran);

    lIntScanTran = (intpair_t *)malloc(R*sizeof(intpair_t));
    assert_malloc(lIntScanTran);
    lScans       = (intpair_t *)malloc(R*sizeof(intpair_t));
    assert_malloc(lScans);

#if TIMING
    START_STEP();
#endif

    if (R < NODES)
	R = NODES;

    s = R/NODES;

    for (k=0 ; k<R ; k++)
	lIndex[k] = 0;

    for (k=0 ; k<q ; k++) 
      lIndex[bits(lKey[k],bitOff,m)]++;

#if TIMING
    END_STEP1("F1histo",q,bitOff);
    START_STEP();
#endif

/* Perform Matrix Transpose of Index -> ScanTran */    

    UMD_Alltoall_i(lIndex, s, lScanTran);

#if TIMING
    END_STEP1("F2trans",q,bitOff);
    START_STEP();
#endif

    /* Calculate local Prefix Sums */    
    for (j=0 ; j<s ; j++)
	for (k=1 ; k<NODES ; k++)
	    lScanTran[k*s + j] += lScanTran[(k-1)*s + j];

/* Create IntLeaveScan by interleaving ScanTran[.][*] with ScanTran[P-1][*]*/

    for (j=0 ; j<s ; j++)
	for (k=0 ; k<NODES ; k++) {
	    lIntScanTran[k*s + j].a = lScanTran[k*s + j];
	    lIntScanTran[k*s + j].b = lScanTran[(NODES-1)*s + j];
	}

    
#if TIMING
    END_STEP1("F3count",q,bitOff);
    START_STEP();
#endif

    /* Perform Inverse Matrix Transpose of IntScanTran -> Scans */

    UMD_Alltoall_i((int *)lIntScanTran, (s<<1), (int *)lScans);
    
#if TIMING
    END_STEP1("F4itran",q,bitOff);
    START_STEP();
#endif

    offset = 0;

    for (k=0 ; k<R ; k++) {
	lIndex[k] = (lScans[k].a - lIndex[k]) + offset;
	offset  += lScans[k].b;
    }

    for (k=0 ; k<q ; k++) {
      j = bits(lKey[k],bitOff,m);
      lAddr[k] = lIndex[j];
      lIndex[j]++;
    }

#if TIMING
    END_STEP1("F5addrp",q,bitOff);
    START_STEP();
#endif
    all_route_q(q, lKey, lAddr, lSorted); 
#if TIMING
    END_STEP1("F6route",q,bitOff);
    REPORT_STEP1("FTOTAL",q,bitOff);
#endif

    free(lScans);
    free(lIntScanTran);
    free(lScanTran);
    free(lIndex);
    free(lAddr);
}

/****************************************************/
void all_countsort_mpi(int q,
		       int *lKey,
		       int *lSorted,
		       int R,
		       int bitOff, int m)
/****************************************************/
/* R (range)      must be a multiple of NODES */
/* q (elems/proc) must be a multiple of NODES */
{

    register int
	j,
	k,
	s,
	offset;
    
    int *lAddr,
	*lIndex,
	*lScanTran;

    intpair_t
	*lIntScanTran,
	*lScans;

#if TIMING
    INIT_STEP();
#endif
    
    lAddr     = (int *)malloc(q*sizeof(int));
    assert_malloc(lAddr);
    lIndex    = (int *)malloc(R*sizeof(int));
    assert_malloc(lIndex);
    lScanTran = (int *)malloc(R*sizeof(int));
    assert_malloc(lScanTran);

    lIntScanTran = (intpair_t *)malloc(R*sizeof(intpair_t));
    assert_malloc(lIntScanTran);
    lScans       = (intpair_t *)malloc(R*sizeof(intpair_t));
    assert_malloc(lScans);

#if TIMING
    START_STEP();
#endif

    if (R < NODES)
	R = NODES;

    s = R/NODES;

    for (k=0 ; k<R ; k++)
	lIndex[k] = 0;

    for (k=0 ; k<q ; k++) 
      lIndex[bits(lKey[k],bitOff,m)]++;

#if TIMING
    END_STEP1("M1histo",q,bitOff);
    START_STEP();
#endif

/* Perform Matrix Transpose of Index -> ScanTran */    

    MPI_Alltoall(lIndex,   s,MPI_INT,
		 lScanTran,s,MPI_INT,
		 MPI_COMM_WORLD);

#if TIMING
    END_STEP1("M2trans",q,bitOff);
    START_STEP();
#endif

    /* Calculate local Prefix Sums */    
    for (j=0 ; j<s ; j++)
	for (k=1 ; k<NODES ; k++)
	    lScanTran[k*s + j] += lScanTran[(k-1)*s + j];

/* Create IntLeaveScan by interleaving ScanTran[.][*] with ScanTran[P-1][*]*/

    for (j=0 ; j<s ; j++)
	for (k=0 ; k<NODES ; k++) {
	    lIntScanTran[k*s + j].a = lScanTran[k*s + j];
	    lIntScanTran[k*s + j].b = lScanTran[(NODES-1)*s + j];
	}

    
#if TIMING
    END_STEP1("M3count",q,bitOff);
    START_STEP();
#endif

    /* Perform Inverse Matrix Transpose of IntScanTran -> Scans */

    MPI_Alltoall(lIntScanTran, s, MPI_INTPAIR,
		 lScans,       s, MPI_INTPAIR,
		 MPI_COMM_WORLD);
    
#if TIMING
    END_STEP1("M4itran",q,bitOff);
    START_STEP();
#endif

    offset = 0;

    for (k=0 ; k<R ; k++) {
	lIndex[k] = (lScans[k].a - lIndex[k]) + offset;
	offset  += lScans[k].b;
    }

    for (k=0 ; k<q ; k++) {
      j = bits(lKey[k],bitOff,m);
      lAddr[k] = lIndex[j];
      lIndex[j]++;
    }

#if TIMING
    END_STEP1("M5addrp",q,bitOff);
    START_STEP();
#endif
    all_route_q_mpi(q, lKey, lAddr, lSorted); 
#if TIMING
    END_STEP1("M6route",q,bitOff);
    REPORT_STEP1("MTOTAL",q,bitOff);
#endif

    free(lScans);
    free(lIntScanTran);
    free(lScanTran);
    free(lIndex);
    free(lAddr);
}

/****************************************************/
void all_countsort_simple(int q,
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
	s,
	offset;
    
    int *lAddr,
        *myHisto,
        *mhp,
        *allIndex,
        *myIndex,
	*lScanTran;

    intpair_t
	*lIntScanTran,
	*lScans;

#if TIMING
    INIT_STEP_TH();
#endif

    lAddr     = (int *)node_malloc(q*sizeof(int), TH);
    lScanTran = (int *)node_malloc(R*sizeof(int), TH);

    lIntScanTran = (intpair_t *)node_malloc(R*sizeof(intpair_t), TH);
    lScans       = (intpair_t *)node_malloc(R*sizeof(intpair_t), TH);

    myHisto = (int *)node_malloc(THREADS*R*sizeof(int), TH);

    allIndex = (int *)node_malloc(R*sizeof(int), TH);

    myIndex = (int *)malloc(R*sizeof(int));
    assert_malloc(myIndex);
    
#if TIMING
    START_STEP_TH();
#endif

    if (R < NODES)
	R = NODES;

    s = R/NODES;

    mhp = myHisto + MYTHREAD*R;

    for (k=0 ; k<R ; k++)
      mhp[k] = 0;
    
    pardo(k, 0, q, 1)
      mhp[bits(lKey[k],bitOff,m)]++;

    node_Barrier();

    pardo(k, 0, R, 1) {
      allIndex[k] = myHisto[k];
      for (j=1 ; j<THREADS ; j++)
	allIndex[k] += myHisto[j*R +  k];
    }

    
#if TIMING
    END_STEP1_TH("S1histo",q,bitOff);
    START_STEP_TH();
#endif

/* Perform Matrix Transpose of Index -> ScanTran */    

    node_Barrier();
    all_Alltoall_i(allIndex, s, lScanTran, TH);
    node_Barrier();

#if TIMING
    END_STEP1_TH("S2trans",q,bitOff);
    START_STEP_TH();
#endif

    /* Calculate local Prefix Sums */
    pardo(j, 0, s, 1)
      for (k=1 ; k<NODES ; k++)
	lScanTran[k*s + j] += lScanTran[(k-1)*s + j];

    node_Barrier();
    
/* Create IntLeaveScan by interleaving ScanTran[.][*] with ScanTran[P-1][*]*/

    pardo(k, 0, NODES, 1)
      for (j=0 ; j<s ; j++) {
	lIntScanTran[k*s + j].a = lScanTran[k*s + j];
	lIntScanTran[k*s + j].b = lScanTran[(NODES-1)*s + j];
      }

#if TIMING
    END_STEP1_TH("S3count",q,bitOff);
    START_STEP_TH();
#endif

    /* Perform Inverse Matrix Transpose of IntScanTran -> Scans */

    node_Barrier();
    all_Alltoall_i((int *)lIntScanTran, (s<<1), (int *)lScans, TH);
    node_Barrier();
    
#if TIMING
    END_STEP1_TH("S4itran",q,bitOff);
    START_STEP_TH();
#endif
    
    offset = 0;

    for (k=0 ; k<R ; k++) {
      myIndex[k] = (lScans[k].a - allIndex[k]) + offset;
      for (j=0 ; j<MYTHREAD ; j++)
	myIndex[k] += myHisto[j*R + k];
      offset  += lScans[k].b;
    }

    pardo(k, 0, q, 1) {
      j = bits(lKey[k],bitOff,m);
      lAddr[k] = myIndex[j];
      myIndex[j]++;
    }

#if TIMING
    END_STEP1_TH("S5addrp",q,bitOff);
    START_STEP_TH();
#endif

    node_Barrier();
#if 1
    ti->udata = bitOff;
#endif
    all_route_q_simple(q, lKey, lAddr, lSorted, TH);

#if TIMING
    END_STEP1_TH("S6route",q,bitOff);
    REPORT_STEP1_TH("STOTAL",q,bitOff);
#endif

    free(myIndex);
    
    node_free(allIndex, TH);
    node_free(myHisto, TH);
    node_free(lScans, TH);
    node_free(lIntScanTran, TH);
    node_free(lScanTran, TH);
    node_free(lAddr, TH);
}


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
void all_radixsort_check(int q,
			 int *lSorted)
/****************************************************/
{
  int i,
    last,
    *prev;

  prev = (int *)malloc(NODES*sizeof(int));
  assert_malloc(prev);

  for (i=1; i<q ; i++) 
    if (lSorted[i-1] > lSorted[i]) {
      fprintf(outfile,
	      "ERROR: PE%3d: q:%d lSorted[%6d] > lSorted[%6d] (%6d,%6d)\n",
	      MYNODE,q,i-1,i,lSorted[i-1],lSorted[i]);
      fflush(outfile);
    }

  last = lSorted[q-1];
#if 0
  MPI_Allgather(&last, 1, MPI_INT, prev, 1, MPI_INT, MPI_COMM_WORLD);
#else
  UMD_Gather_i(&last, 1, prev);
  UMD_Bcast_buf_i(prev,NODES,prev);  /* <------ was bad, now fixed */
#endif

  if (MYNODE>0)
    if (lSorted[0] < prev[MYNODE-1]) {
      fprintf(outfile,
	      "ERROR: PE%3d: lSorted[%6d] < prev[%3d] (%6d,%6d)\n",
	      MYNODE,0,MYNODE-1,lSorted[0],prev[MYNODE-1]);
      fflush(outfile);
    }

  free(prev);
}

/****************************************************/
void all_radixsort_s3(int q,
		      int *lKeys,
		      int *lSorted)
/****************************************************/
{
    all_countsort(q, lKeys,   lSorted, (1<<11),  0, 11);
    all_countsort(q, lSorted, lSorted, (1<<11), 11, 11);
    all_countsort(q, lSorted, lSorted, (1<<10), 22, 10);
}

/****************************************************/
void all_radixsort_s2(int q,
		      int *lKeys,
		      int *lSorted)
/****************************************************/
{
    all_countsort(q, lKeys,   lSorted, (1<<16),  0, 16);
    all_countsort(q, lSorted, lSorted, (1<<16), 16, 16);
}

/****************************************************/
void all_radixsort20_s1(int q,
			int *lKeys,
			int *lSorted)
/****************************************************/
{
    all_countsort(q, lKeys,   lSorted, (1<<20),  0, 20);
}

/****************************************************/
void all_radixsort20_s2(int q,
			int *lKeys,
			int *lSorted)
/****************************************************/
{
    all_countsort(q, lKeys,   lSorted, (1<<10),  0, 10);
    all_countsort(q, lSorted, lSorted, (1<<10), 10, 10);
}


/****************************************************/
void all_radixsort_mpi_s3(int q,
			  int *lKeys,
			  int *lSorted)
/****************************************************/
{
    all_countsort_mpi(q, lKeys,   lSorted, (1<<11),  0, 11);
    all_countsort_mpi(q, lSorted, lSorted, (1<<11), 11, 11);
    all_countsort_mpi(q, lSorted, lSorted, (1<<10), 22, 10);
}

/****************************************************/
void all_radixsort_mpi_s2(int q,
			  int *lKeys,
			  int *lSorted)
/****************************************************/
{
    all_countsort_mpi(q, lKeys,   lSorted, (1<<16),  0, 16);
    all_countsort_mpi(q, lSorted, lSorted, (1<<16), 16, 16);
}

/****************************************************/
void all_radixsort20_mpi_s1(int q,
			    int *lKeys,
			    int *lSorted)
/****************************************************/
{
    all_countsort_mpi(q, lKeys,   lSorted, (1<<20),  0, 20);
}

/****************************************************/
void all_radixsort20_mpi_s2(int q,
			    int *lKeys,
			    int *lSorted)
/****************************************************/
{
    all_countsort_mpi(q, lKeys,   lSorted, (1<<10),  0, 10);
    all_countsort_mpi(q, lSorted, lSorted, (1<<10), 10, 10);
}



/****************************************************/
void all_radixsort_simple_s3(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED)
/****************************************************/
{
    all_countsort_simple(q, lKeys,   lSorted, (1<<11),  0, 11, TH);
    all_countsort_simple(q, lSorted, lSorted, (1<<11), 11, 11, TH);
    all_countsort_simple(q, lSorted, lSorted, (1<<10), 22, 10, TH);
}

/****************************************************/
void all_radixsort_simple_s2(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED)
/****************************************************/
{
    all_countsort_simple(q, lKeys,   lSorted, (1<<16),  0, 16, TH);
    all_countsort_simple(q, lSorted, lSorted, (1<<16), 16, 16, TH);
}

/****************************************************/
void all_radixsort20_simple_s1(int q,
			       int *lKeys,
			       int *lSorted,
			       THREADED)
/****************************************************/
{
    all_countsort_simple(q, lKeys,   lSorted, (1<<20),  0, 20, TH);
}

/****************************************************/
void all_radixsort20_simple_s2(int q,
			       int *lKeys,
			       int *lSorted,
			       THREADED)
/****************************************************/
{
    all_countsort_simple(q, lKeys,   lSorted, (1<<10),  0, 10, TH);
    all_countsort_simple(q, lSorted, lSorted, (1<<10), 10, 10, TH);
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
void all_radixsort20_node_s1(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED)
/****************************************************/
{
    all_countsort_node(q, lKeys,   lSorted, (1<<20),  0, 20, TH);
}

/****************************************************/
void all_radixsort20_node_s2(int q,
			     int *lKeys,
			     int *lSorted,
			     THREADED)
/****************************************************/
{
  int *lTemp;

    lTemp = (int *)node_malloc(q*sizeof(int), TH);
			
    all_countsort_node(q, lKeys,   lTemp,   (1<<10),  0, 10, TH);
    all_countsort_node(q, lTemp,   lSorted, (1<<10), 10, 10, TH);

    node_free(lTemp, TH);
}










