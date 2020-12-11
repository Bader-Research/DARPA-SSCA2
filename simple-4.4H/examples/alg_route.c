#include "alg_route.h"

#define MPIR_ALLTOALLV_TAG 10

#define DEBUG       0
#define TIMING      0
#define TIMING_COMM 0

int init_Alltoallv_param = 0;
MPI_Request *reqarray;
MPI_Status  *starray;

void init_Alltoallv() {
  reqarray = (MPI_Request *)malloc(NODES*sizeof(MPI_Request));
  assert_malloc(reqarray);
  starray  = (MPI_Status  *)malloc(NODES*sizeof(MPI_Status));
  assert_malloc(starray);
  init_Alltoallv_param = 1;
}

void destroy_Alltoallv() {
  free(starray);
  free(reqarray);
  init_Alltoallv_param = 0;
}

int my_Alltoallv ( void *sendbuf, int *sendcnts,
		   int *sdispls, MPI_Datatype sendtype, 
		   void *recvbuf, int *recvcnts, 
		   int *rdispls, MPI_Datatype recvtype,
		   MPI_Comm comm )
{
  int        i, k;
  MPI_Aint   send_extent, recv_extent;
  
  MPI_Type_extent(sendtype, &send_extent);
  MPI_Type_extent(recvtype, &recv_extent);

  if (!init_Alltoallv_param)
    init_Alltoallv();
       
  for ( k=0; k<NODES; k++ ) {
    i = k ^ MYNODE;
    MPI_Irecv((void *)((char *)recvbuf+rdispls[i]*recv_extent), 
	      recvcnts[i], recvtype, i, MPIR_ALLTOALLV_TAG,
	      comm, reqarray+i);
    MPI_Send((void *)((char *)sendbuf+sdispls[i]*send_extent), 
	     sendcnts[i], sendtype, i, MPIR_ALLTOALLV_TAG,
	     comm);
  }
  
  MPI_Waitall(NODES,reqarray,starray);
  
  return (0);
}

int all_route_q(int q, 
		int *lKeys,
		int *lAddr,
		int *lRout)
{
    register int i, hq_mask;

    int
      j,
      LOGQ,
      *lsend_off,
      *lsend_cnt,
      *lrecv_off,
      *lrecv_cnt;
    
    intpair_t 
      *lsend,
      *lrecv,
      **bptr;

#if TIMING
    INIT_STEP();
#endif

    lsend_cnt = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lsend_cnt);

    lsend_off = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lsend_off);

    lrecv_cnt = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lrecv_cnt);

    lrecv_off = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lrecv_off);

    lsend = (intpair_t *)malloc(q*sizeof(intpair_t));
    assert_malloc(lsend);

    lrecv = (intpair_t *)malloc(q*sizeof(intpair_t));
    assert_malloc(lrecv);

    bptr = (intpair_t **)malloc(NODES*sizeof(intpair_t *));
    assert_malloc(bptr);

#if TIMING
    START_STEP();
#endif

    LOGQ = log_2(q);

    for (i=0 ; i<NODES ; i++) 
      lsend_cnt[i] = 0;

    for (i=0 ; i<q ; i++) 
      lsend_cnt[lAddr[i] >> LOGQ]++;

#if TIMING
    END_STEP1("Froute:count",q,0);
    START_STEP();
#endif

    UMD_Alltoall_i(lsend_cnt, 1, lrecv_cnt);

#if TIMING
    END_STEP1("Froute:trans",q,0);
    START_STEP();
#endif

    bptr[0] = lsend;
    for (i=1 ; i<NODES ; i++) 
	bptr[i] = bptr[i-1] + lsend_cnt[i-1];
    for (i=0 ; i<q ; i++) {
      j = (lAddr[i]>>LOGQ); 
      bptr[j]->a = lKeys[i]; 
      bptr[j]->b = lAddr[i]; 
      bptr[j]++;
    }
    
    lsend_off[0] = 0;
    lrecv_off[0] = 0;
    for (i=1 ; i<NODES ; i++) {
	lsend_off[i] = lsend_off[i-1] + lsend_cnt[i-1];
	lrecv_off[i] = lrecv_off[i-1] + lrecv_cnt[i-1];
    }

#if TIMING
    END_STEP1("Froute:offst",q,0);
    START_STEP();
#endif

    for (i=0 ; i<NODES ; i++) {
      lsend_cnt[i] *= 2;
      lrecv_cnt[i] *= 2;
      lsend_off[i] *= 2;
      lrecv_off[i] *= 2;
    }
    UMD_Alltoallv_i((int *)lsend,lsend_cnt,lsend_off,
		     (int *)lrecv,lrecv_cnt,lrecv_off);
		  
#if TIMING
    END_STEP1("Froute:alltv",q,0);
    START_STEP();
#endif

    hq_mask = q - 1;
    for (i=0 ; i<q ; i++) 
	*(lRout + (lrecv[i].b & hq_mask)) = lrecv[i].a;

#if TIMING
    END_STEP1("Froute:place",q,0);
    REPORT_STEP1("Froute:TOTAL",q,0);
    on_one_node {
      fprintf(outfile,"BW/node: %9.6f\n",
	      2.0 * (double)q*4.0 / tsec[0] / (double)(1<<20));
      fflush(outfile);
    }
#endif

    free(bptr);
    free(lrecv);
    free(lsend);
    free(lrecv_off);
    free(lrecv_cnt);
    free(lsend_off);
    free(lsend_cnt);

    return(q);
}

int all_route_q_mpi(int q, 
		    int *lKeys,
		    int *lAddr,
		    int *lRout)
{
    register int i, hq_mask;

    int
      j,
      LOGQ,
      *lsend_off,
      *lsend_cnt,
      *lrecv_off,
      *lrecv_cnt;
    
    intpair_t 
      *lsend,
      *lrecv,
      **bptr;

#if TIMING
    INIT_STEP();
#endif

    lsend_cnt = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lsend_cnt);

    lsend_off = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lsend_off);

    lrecv_cnt = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lrecv_cnt);

    lrecv_off = (int *)malloc(NODES*sizeof(int));
    assert_malloc(lrecv_off);

    lsend = (intpair_t *)malloc(q*sizeof(intpair_t));
    assert_malloc(lsend);

    lrecv = (intpair_t *)malloc(q*sizeof(intpair_t));
    assert_malloc(lrecv);

    bptr = (intpair_t **)malloc(NODES*sizeof(intpair_t *));
    assert_malloc(bptr);

#if TIMING
    START_STEP();
#endif

    LOGQ = log_2(q);

    for (i=0 ; i<NODES ; i++) 
      lsend_cnt[i] = 0;

    for (i=0 ; i<q ; i++) 
      lsend_cnt[lAddr[i] >> LOGQ]++;

#if TIMING
    END_STEP("Mroute:count");
    START_STEP();
#endif

    MPI_Alltoall(lsend_cnt, 1, MPI_INT,
		 lrecv_cnt, 1, MPI_INT,
		 MPI_COMM_WORLD);

#if TIMING
    END_STEP("Mroute:trans");
    START_STEP();
#endif

    bptr[0] = lsend;
    for (i=1 ; i<NODES ; i++) 
	bptr[i] = bptr[i-1] + lsend_cnt[i-1];
    for (i=0 ; i<q ; i++) {
      j = (lAddr[i]>>LOGQ); 
      bptr[j]->a = lKeys[i]; 
      bptr[j]->b = lAddr[i]; 
      bptr[j]++;
    }
    
    lsend_off[0] = 0;
    lrecv_off[0] = 0;
    for (i=1 ; i<NODES ; i++) {
	lsend_off[i] = lsend_off[i-1] + lsend_cnt[i-1];
	lrecv_off[i] = lrecv_off[i-1] + lrecv_cnt[i-1];
    }

#if TIMING
    END_STEP("Mroute:offst");
    START_STEP();
#endif

    my_Alltoallv(lsend,lsend_cnt,lsend_off,MPI_INTPAIR,
		  lrecv,lrecv_cnt,lrecv_off,MPI_INTPAIR,
		  MPI_COMM_WORLD);
		  
#if TIMING
    END_STEP("Mroute:alltv");
    START_STEP();
#endif

    hq_mask = q - 1;
    for (i=0 ; i<q ; i++) 
	*(lRout + (lrecv[i].b & hq_mask)) = lrecv[i].a;

#if TIMING
    END_STEP("Mroute:place");
    REPORT_STEP("Mroute:TOTAL");
    on_one_node {
      fprintf(outfile,"BW/node: %9.6f\n",
	      2.0 * (double)q*4.0 / tsec[0] / (double)(1<<20));
      fflush(outfile);
    }
#endif

    free(bptr);
    free(lrecv);
    free(lsend);
    free(lrecv_off);
    free(lrecv_cnt);
    free(lsend_off);
    free(lsend_cnt);

    return(q);
}

int all_route_q_simple(int q, 
		       int *lKeys,
		       int *lAddr,
		       int *lRout,
		       THREADED)
{
    register int i, hq_mask;

    int
      j,
      LOGQ,
      *lsend_off,
      *lsend_cnt,
      *lrecv_off,
      *lrecv_cnt,
      *msend_cnt,
      *ps_lsend_cnt,
      *msc;
    
    intpair_t 
      *lsend,
      *lrecv,
      **bptr;

#if TIMING_COMM
    double secs_comm, tsec_comm;
#endif
#if TIMING
    INIT_STEP_TH();
#endif

    lsend_cnt = (int *)node_malloc(NODES*sizeof(int), TH);
    msend_cnt = (int *)node_malloc(THREADS*NODES*sizeof(int), TH);
    lsend_off = (int *)node_malloc(NODES*sizeof(int), TH);
    lrecv_cnt = (int *)node_malloc(NODES*sizeof(int), TH);
    lrecv_off = (int *)node_malloc(NODES*sizeof(int), TH);
    lsend = (intpair_t *)node_malloc(q*sizeof(intpair_t), TH);
    lrecv = (intpair_t *)node_malloc(q*sizeof(intpair_t), TH);

    bptr = (intpair_t **)malloc(NODES*sizeof(intpair_t *));
    assert_malloc(bptr);
    
    ps_lsend_cnt = (int *)malloc(NODES*sizeof(int));
    assert_malloc(ps_lsend_cnt);

#if TIMING
    START_STEP_TH();
#endif
    
    LOGQ = log_2(q);

    msc = msend_cnt + MYTHREAD*NODES;
    for (i=0 ; i<NODES ; i++)
      msc[i] = 0;

    pardo(i, 0, q, 1)
#if DEBUG 
      {
	int dst;
	dst = lAddr[i]>>LOGQ;
	if (dst >= NODES) {
	  fprintf(outfile,
		  "PE%3d(%3d): ERROR: route dst>=NODES dst:%d i:%d LOGQ:%d\n",
		  MYNODE,MYTHREAD,dst,i,LOGQ);
	  fflush(outfile);
	}
	msc[dst]++;
      }
#else
      msc[lAddr[i] >> LOGQ]++;
#endif
    
    node_Barrier();

    pardo(i, 0, NODES, 1) {
      lsend_cnt[i] = msend_cnt[i];
      for (j=1 ; j<THREADS ; j++)
	lsend_cnt[i] += msend_cnt[j*NODES + i];
    }

#if TIMING
    END_STEP1_TH("S1route:count",q,ti->udata);
    START_STEP_TH();
#endif

    node_Barrier();
    all_Alltoall_i(lsend_cnt, 1, lrecv_cnt, TH);
    node_Barrier();

#if DEBUG
    on_one_thread {
      if (lsend_cnt[MYNODE] != lrecv_cnt[MYNODE]) {
	fprintf(outfile,
		"PE%3d: ERROR: sendcnt!=recvcnt (%d, %d)\n",
		MYNODE,lsend_cnt[MYNODE],lrecv_cnt[MYNODE]);
	fflush(outfile);
      }
    }
    node_Barrier();
#endif

#if TIMING
    END_STEP1_TH("S2route:trans",q,ti->udata);
    START_STEP_TH();
#endif

    /* Need to set up bptr correctly here! */
    /* lsend: send buffer */
    /* lsend_cnt[i]: # of el's this node is sending to node i */
    /* msend_cnt[j*NODES + i]: # of el's that thread j has for node i */
    /* To calculate where my thread's first element to node i goes:
       0: lsend ,  lsend+msend_cnt[0] , lsend+msend_cnt[NODES], ...
       1: lsend+lsend_cnt[0], lsend+lsend_cnt[0] + msend_cnt[1], */
    /* I need PS of lend_cnt */

    ps_lsend_cnt[0] = 0;
    for (i=1 ; i<NODES ; i++)
      ps_lsend_cnt[i] = ps_lsend_cnt[i-1] + lsend_cnt[i-1];

    for (j=0 ; j<NODES ; j++) {
      bptr[j] = lsend + ps_lsend_cnt[j];
      for (i=0 ; i<MYTHREAD ; i++)
	bptr[j] += msend_cnt[i*NODES + j];
    }
    
    pardo(i, 0, q, 1) {
      j = (lAddr[i]>>LOGQ); 
#if DEBUG
      if (j>=NODES) {
	fprintf(outfile,
		"PE%3d: ERROR: j>=NODES  j:%d i:%d addr[i]:%d\n",
		MYNODE,j,i,lAddr[i]);
	fflush(outfile);
      }
#endif
      bptr[j]->a = lKeys[i]; 
      bptr[j]->b = lAddr[i]; 
      bptr[j]++;
    }

    on_one_thread {
      lsend_off[0] = 0;
      lrecv_off[0] = 0;
      for (i=1 ; i<NODES ; i++) {
	lsend_off[i] = lsend_off[i-1] + lsend_cnt[i-1];
	lrecv_off[i] = lrecv_off[i-1] + lrecv_cnt[i-1];
      }
    }
    
#if TIMING
    END_STEP1_TH("S3route:offst",q,ti->udata);
    START_STEP_TH();
#endif

    node_Barrier();
    
    pardo(i, 0, NODES, 1) {
      lsend_cnt[i] *= 2;
      lrecv_cnt[i] *= 2;
      lsend_off[i] *= 2;
      lrecv_off[i] *= 2;
    }

    node_Barrier();
#if TIMING_COMM
    secs_comm = get_seconds();
#endif
    all_Alltoallv_i((int *)lsend,lsend_cnt,lsend_off,
		    (int *)lrecv,lrecv_cnt,lrecv_off, TH);
#if TIMING_COMM
    secs_comm = get_seconds() - secs_comm;
    tsec_comm = all_Reduce_d(secs_comm, MAX, TH); 
    on_one {
      fprintf(outfile,"(%3d %3d) n: %12d bo: %3d %s: %9.6f\n",
	      NODES,THREADS,q*NODES,ti->udata,"Sroute:alltv",tsec_comm); 
      fflush(outfile);
    }
#endif
    node_Barrier();

#if TIMING
    END_STEP1_TH("S4route:alltv",q,ti->udata);
    START_STEP_TH();
#endif

    hq_mask = q - 1;
    pardo(i, 0, q, 1)
      *(lRout + (lrecv[i].b & hq_mask)) = lrecv[i].a;

#if TIMING
    END_STEP1_TH("S5route:place",q,ti->udata);
#if 0
    REPORT_STEP1_TH("Sroute:TOTAL",q,ti->udata); 
    on_one {
      fprintf(outfile,"BW/node: %9.6f\n",
	      2.0 * (double)q*4.0 / tsec[0] / (double)(1<<20));
      fflush(outfile);
    }
#endif
#endif

    node_Barrier();

    free(ps_lsend_cnt);
    free(bptr);

    node_free(lrecv, TH);
    node_free(lsend, TH);
    node_free(lrecv_off, TH);
    node_free(lrecv_cnt, TH);
    node_free(lsend_off, TH);
    node_free(msend_cnt, TH);
    node_free(lsend_cnt, TH);

    return(q);
}

