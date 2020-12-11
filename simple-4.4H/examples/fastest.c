#include"simple.h"
#include "umd.h"
#ifdef RRANDOM
#include "alg_random.h"
#endif
#include "alg_2dfft.h"

#define DEBUG              0

#define TEST_WIDE          1

#if TEST_WIDE
#define ARR_SIZE        (1<<10)
#define ARR_SIZE_SORT   /*(1<<23)*/ (1<<10)
#define ARR_SIZE_SELSUM (1<<10)
#else
#define ARR_SIZE        (1<<11)
#define ARR_SIZE_SORT   (1<<14)
#define ARR_SIZE_SELSUM (1<<14)
#endif
#define LOOP_BARRIER        100
#define LOOP_SEND           100
#define LOOP_TRANSPOSE        5
#define MIN_TIME       0.000001

#define NQUEENS_N          15
#define NQUEENS_MAXK        3

#define TEST_TEST          0
#define TEST_TESTNODE      0
#define TEST_BARRIER       0
#define TEST_BARRIER_NODE  0
#define TEST_NODEBCAST     1
#define TEST_NODEREDUCE    0
#define TEST_ALLREDUCE     0
#define TEST_ALLBCAST      0
#define TEST_ALL_ALLREDUCE 0
#define TEST_ALL_ALLTOALL  1
#define TEST_NQUEENS       0
#define  TEST_NQUEENS_R    0
#define  TEST_NQUEENS_NR   0
#define  TEST_NQUEENS_X    0
#define TEST_HISTO         0

#define TEST_SELSUM        0
#define TEST_PING          0
#define TEST_PING_NB       0
#define TEST_TRANSPOSE     0
#define TEST_TRANSPOSE_S   0
#define TEST_GATHER        0
#define TEST_GATHER_S      0
#define TEST_SCATTER       0
#define TEST_SCATTER_S     0
#define TEST_MPI           0
#define TEST_SORT          0
#define TEST_SORT_M        0
#define TEST_SORT_S        0
#define TEST_SPLASH_OCEAN  0
#define TEST_SPLASH_RADIX  0
#define TEST_PRIME_M       0
#define TEST_PRIME         0
#define TEST_2DFFT_M       0
#define TEST_2DFFT_S       0

#if TEST_NQUEENS
#include "alg_nqueens.h"
#endif

#if TEST_HISTO
#include "alg_histo.h"
#endif

#if (TEST_SORT||TEST_SORT_M||TEST_SORT_S)
#define TIMING             1
#include "alg_route.h"
#include "alg_radix.h"
#endif
MPI_Datatype MPI_INTPAIR;

void
pattern2(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++)
    arr[i] = i;
}

void
check_pattern2(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++) 
    if (arr[i] != i)
      fprintf(stderr,"PE%3d: ERROR: arr[%12d] != %12d  (%12d)\n",
	      MYNODE,i,i,arr[i]);
}

void
pattern3(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++)
    arr[i] = 0;
}

void
pattern4(int *arr, int len) {
  int i, k;
  k = len/NODES;
  for (i=0 ; i<len ; i++)
    arr[i] = i/k;
}

void
check_pattern4(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++) 
    if (arr[i] != MYNODE)
      fprintf(stderr,"PE%3d: ERROR: arr[%12d] != %12d  (%12d)\n",
	      MYNODE,i,i,arr[i]);
}

void
pattern5(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++)
    arr[i] = MYNODE + i;
}

void
check_pattern5(int *arr, int len) {
  int i, j;
  for (j=0 ; j<NODES ; j++)
    for (i=0 ; i<len ; i++) 
      if (arr[j*len + i] != j+i)
	fprintf(stderr,"PE%3d: ERROR: arr[%12d] != %12d  (%12d)\n",
		MYNODE,j*len + i,j+i,arr[j*len+i]);
}

void
pattern6(int *arr, int len) {
  int i,j;
  for (j=0 ; j<NODES ; j++)
    for (i=0 ; i<len ; i++)
      arr[j*len+i] = j + i;
}

void
check_pattern6(int *arr, int len) {
  int i;
  for (i=0 ; i<len ; i++) 
    if (arr[i] != MYNODE + i)
	fprintf(stderr,"PE%3d: ERROR: arr[%12d] != %12d  (%12d)\n",
		MYNODE,i,MYNODE+i,arr[i]);
}

#if TEST_TEST

#define DEBUG_TEST 0

#if 0
void
test_test(THREADED) {
  return;
}
#else
void
test_test(THREADED) {
  int *arr1, *arr2, *arr3;
  int len;
  int i, j, src, dst;
  int tag;
  
  len = 10;

  arr1 = (int *)malloc(len*sizeof(int));
  assert_malloc(arr1);
  arr2 = (int *)malloc(len*sizeof(int));
  assert_malloc(arr2);
  arr3 = (int *)malloc(len*sizeof(int));
  assert_malloc(arr3);
  
  all_Barrier(TH);
  on_one_thread {
    for (j=0 ; j<10000 ; j++) {
      if (j%NODES != 0) {
	dst = (MYNODE+j) % NODES;
	src = (MYNODE-(j%NODES)+NODES) % NODES;
	fprintf(outfile,"PE%3d: j: %6d src: %3d dst: %3d\n",
		MYNODE, j, src, dst);
	fflush(outfile);
	for (i=0 ; i<len ; i++) 
	  arr1[i] = i + MYNODE;

#if DEBUG_TEST
	fprintf(outfile,"PE%3d: ISend: %3d (before)\n", MYNODE, dst);
	fflush(outfile);
#endif
	UMD_ISend(dst, arr1, len*sizeof(int), &tag);
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: ISend: %3d (after)\n", MYNODE, dst);
	fprintf(outfile,"PE%3d: tst0: %3d \n", MYNODE, UMD_Test(tag));
	fflush(outfile);
#endif
	
	UMD_Recv(src, arr2, len*sizeof(int));
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: Recv: %3d \n", MYNODE, src);
	fprintf(outfile,"PE%3d: tst1: %3d \n", MYNODE, UMD_Test(tag));
	fprintf(outfile,"PE%3d: waiting for %3d \n", MYNODE, tag);
	fflush(outfile);
#endif
	UMD_Wait(tag);
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: Wait DONE!\n", MYNODE);
	fflush(outfile);
#endif
	for (i=0 ; i<len ; i++) {
	  if (arr2[i] != i + src) {
	    fprintf(outfile,"PE%3d: ERROR: arr2[%6d]: %12d  src: %3d\n",
		    MYNODE,i,arr2[i],src);
	    fflush(outfile);
	  }
	}

#if DEBUG_TEST
	fprintf(outfile,"PE%3d: IRecv: %3d (before)\n", MYNODE, dst);
	fflush(outfile);
#endif
	UMD_IRecv(dst, arr3, len*sizeof(int), &tag);
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: IRecv: %3d (after)\n", MYNODE, dst);
	fflush(outfile);
#endif
	
	UMD_Send(src, arr2, len*sizeof(int));
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: Send: %3d \n", MYNODE, src);
	fprintf(outfile,"PE%3d: tst0: %3d \n", MYNODE, UMD_Test(tag));
	fflush(outfile);
#endif
	
	UMD_Wait(tag);
#if DEBUG_TEST
	fprintf(outfile,"PE%3d: Wait DONE!\n", MYNODE);
	fflush(outfile);
#endif
	
	for (i=0 ; i<len ; i++) {
	  if (arr1[i] != arr3[i]) {
	    fprintf(outfile,"PE%3d: ERROR: arr1[%6d]: %12d  arr3[%6d]: %12d\n",
		    MYNODE,i,arr1[i],i,arr3[i]);
	    fflush(outfile);
	  }
	}
	UMD_Barrier();
	fprintf(outfile,"PE%3d: After Barrier  j: %12d\n", MYNODE, j);
	fflush(outfile);
      }
    }
  }
  all_Barrier(TH);
  free(arr1);
  free(arr2);
  free(arr3);
  return;
}
#endif

#if 0

void
test_test(THREADED) {
  int ret, ret1, ret2, pol;
  struct sched_param  prm;
  struct sched_param  prm2;

  ret = pthread_getschedparam(pthread_self(), &pol, &prm);

  printf("A PE%3d:  ret: %d pol: %2d  BG: %2d FG: %2d FIFO: %2d  RR: %2d  OTH: %2d\n",
	 MYTHREAD, ret, pol, SCHED_BG_NP, SCHED_FG_NP,
	 SCHED_FIFO, SCHED_RR, SCHED_OTHER);
  printf("A PE%3d:  pri: %2d  FIFO(%d,%d) RR(%d,%d) OTH(%d,%d) FG(%d, %d) BG(%d,%d)\n",
	 MYTHREAD, prm.sched_priority,
	 PRI_FIFO_MAX,  PRI_FIFO_MIN,
	 PRI_RR_MAX,    PRI_RR_MIN,
	 PRI_OTHER_MAX, PRI_OTHER_MIN,
	 PRI_FG_MAX_NP,    PRI_FG_MIN_NP,
	 PRI_BG_MAX_NP,    PRI_BG_MIN_NP);

  prm2.sched_priority = PRI_FIFO_MAX;

  ret1 = pthread_setschedparam(pthread_self(), SCHED_FIFO, &prm2);
  ret2 = pthread_getschedparam(pthread_self(), &pol, &prm);

  printf("B PE%3d:  ret(%d,%d) pol: %2d  BG: %2d FG: %2d FIFO: %2d  RR: %2d  OTH: %2d\n",
	 MYTHREAD, ret1, ret2, pol, SCHED_BG_NP, SCHED_FG_NP,
	 SCHED_FIFO, SCHED_RR, SCHED_OTHER);
  printf("B PE%3d:  pri: %2d  FIFO(%d,%d) RR(%d,%d) OTH(%d,%d) FG(%d, %d) BG(%d,%d)\n",
	 MYTHREAD, prm.sched_priority,
	 PRI_FIFO_MAX,  PRI_FIFO_MIN,
	 PRI_RR_MAX,    PRI_RR_MIN,
	 PRI_OTHER_MAX, PRI_OTHER_MIN,
	 PRI_FG_MAX_NP,    PRI_FG_MIN_NP,
	 PRI_BG_MAX_NP,    PRI_BG_MIN_NP);
  return;
}
#endif

#if 0
void
test_test(THREADED) {

  int i, *arr, *arrPtr, *lastPtr;
  int j = (1<<20);
  int ret;
  unsigned seed;
  int size;
  double d;
  
  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE A\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  arr = (int *)malloc(j*sizeof(int));
  assert_malloc(arr);

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE B\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  seed = 317*(ID+17);

  fprintf(outfile,"PE%3d(%3d): seed: %u\n",MYNODE,MYTHREAD,seed);
  fflush(outfile);

  all_Barrier(TH);
  
  ret = srandom(seed);
  if (ret<0)
    perror("srandom()");

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE C\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  d = get_seconds();
  arrPtr = arr;
  lastPtr = arr+j;
  while (arrPtr < lastPtr) {
    *arrPtr++ = random();
  }
  d = get_seconds() - d;
  
  all_Barrier(TH);

  fprintf(outfile,"PE%3d(%3d): Aft test= %4d  time per K random()s: %9.6f\n",
	  MYNODE,MYTHREAD,j, d*(1024.0)/(double)j);
  fflush(outfile);

  free(arr);

{
  int i, *arr, *arrPtr, *lastPtr;
  int j = (1<<20);
  int ret;
  unsigned seed;
  int size;
  double d;
  
  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE A\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  arr = (int *)malloc(j*sizeof(int));
  assert_malloc(arr);

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE B\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  seed = 317*(ID+17);

  fprintf(outfile,"PE%3d(%3d): seed: %u\n",MYNODE,MYTHREAD,seed);
  fflush(outfile);

  all_Barrier(TH);
  
  srrandom(seed);

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE C\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  d = get_seconds();
  arrPtr = arr;
  lastPtr = arr+j;
  while (arrPtr < lastPtr) {
    *arrPtr++ = (int)rrandom();
  }
  d = get_seconds() - d;
  
  all_Barrier(TH);

  fprintf(outfile,"PE%3d(%3d): Aft test= %4d  time per K rrandom()s: %9.6f\n",
	  MYNODE,MYTHREAD,j, d*(1024.0)/(double)j);
  fflush(outfile);

  free(arr);
}
#ifdef RRANDOM
{
  int i, *arr, *arrPtr, *lastPtr;
  int j = (1<<20);
  int ret;
  unsigned seed;


  int size;
  double d;
  
  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE A\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  arr = (int *)malloc(j*sizeof(int));
  assert_malloc(arr);

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE B\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  seed = 317*(ID+17);

  fprintf(outfile,"PE%3d(%3d): seed: %u\n",MYNODE,MYTHREAD,seed);
  fflush(outfile);

  all_Barrier(TH);

  srrandom_th(seed,TH);

  all_Barrier(TH);
  on_one {
    fprintf(outfile,"HERE C\n");
    fflush(outfile);
  }
  all_Barrier(TH);

  d = get_seconds();
  arrPtr = arr;
  lastPtr = arr+j;
  while (arrPtr < lastPtr) {
    *arrPtr++ = (int)rrandom_th(TH);
  }
  d = get_seconds() - d;
  
  all_Barrier(TH);

  fprintf(outfile,"PE%3d(%3d): Aft test= %4d  time per K rrandom_th()s: %9.6f\n",
	  MYNODE,MYTHREAD,j, d*(1024.0)/(double)j);
  fflush(outfile);

  free(arr);
}
#endif /* RRANDOM */
}
#endif /* if */
#endif /* TEST_TEST */

#if TEST_TESTNODE
#define LOOPCNT 10
#define WARMCNT   1
void
test_testnode(THREADED) {
  double secs;
  int *nodearr, *arr,
    i, loop, ths, maxarr, sz, warm;
  
  maxarr = 1<<10;
  nodearr = (int *)node_malloc(THREADS*maxarr*sizeof(int),TH);
  arr   = nodearr+(MYTHREAD*maxarr);

  on_one_node {

    for (i=0 ; i<maxarr ; i++)
      arr[i] = i;

    for (sz=1 ; sz<=maxarr ; sz = (sz<<1)) {
      for (warm=0 ; warm<WARMCNT ; warm++) {
	for (ths = 1 ; ths <= THREADS ; ths++) {
	  node_Barrier();
	  if (MYTHREAD < ths) {

	    secs=get_seconds();
	    for (loop=0 ; loop<LOOPCNT ; loop++) {
	      for (i=0 ; i<sz ; i++) {
		arr[i] += i;
	      }
	    }
	    secs=get_seconds()-secs;
	    fprintf(outfile,
		    "%3d sz: %8d ths:%3d time: %9.6f  (us): %9.6f warm %2d (TH: %3d)\n",
		    THREADS,sz,ths,secs,(secs*1000000.0)/(double)(LOOPCNT*sz),warm,MYTHREAD);
	    fflush(outfile);
	  }
	}
      }
    }

  }  
  node_free(nodearr,TH);
}
#endif

#if TEST_NODEBCAST
#define LOOPCNT 100
void test_nodebcast(THREADED) {
double secs;
int *nodearr, *arr, i, maxarr, sz;
maxarr = 1<<20;
on_one_node {
   
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
     nodearr = (int *)node_malloc(THREADS*sz*sizeof(int),TH);
     arr = nodearr+(MYTHREAD*sz);
     for(i=0;i<sz;i++)
     arr[i]=i;
     node_Barrier();
     secs=get_seconds();
	for(i=0;i<LOOPCNT;i++)
	 node_Bcast_i(*arr,TH);
     secs=get_seconds()-secs;
     fprintf(outfile, "%3d th: %2d  sz: %8d time: %9.6f\n",THREADS,MYTHREAD,sz,secs/(double)LOOPCNT); 
     fflush(outfile);
     node_free(nodearr,TH);
   }
  }
}
#endif

#if TEST_NODEREDUCE
#define LOOPCNT 100
void
test_nodereduce(THREADED) {
double secs;
int *nodearr, *arr, i, maxarr, sz ;
maxarr = 1<<20;
on_one_node {
     
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
     nodearr = (int *)node_malloc(THREADS*sz*sizeof(int),TH);
     arr = nodearr+(MYTHREAD*sz);
     for(i=0;i<sz;i++)
     arr[i]=MYTHREAD+i;
     node_Barrier();
     secs=get_seconds();
        for(i=0;i<LOOPCNT;i++)
         node_Reduce_i(*arr,MAX,TH);
     secs=get_seconds()-secs;
     fprintf(outfile, "%3d th: %2d  sz: %8d time: %9.6f\n",THREADS,MYTHREAD,sz,secs/(double)LOOPCNT);
     fflush(outfile);
     node_free(nodearr,TH);
   }
  }
}
#endif

#if TEST_ALLREDUCE
#define LOOPCNT 100
void
test_allreduce(THREADED) {
double secs;   
int  *arr, i,bytes, maxarr, sz ;
maxarr = 1<<15;
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
      bytes=sz*sizeof(int);
     arr = (int *)node_malloc(bytes,TH);
    on_one_thread{ 
     for(i=0;i<sz;i++)
      arr[i]=i;
    }  
     all_Barrier(TH);
     secs=get_seconds();
        for(i=0;i<LOOPCNT;i++)
         all_Reduce_i(*arr,MAX,TH);
     secs=get_seconds()-secs;
  on_one{
     fprintf(outfile, " Bytes : %7d secs:%5.6f\n",bytes,secs/(double)LOOPCNT);
     fflush(outfile);
        }
     node_free(arr,TH);
   
  }
}
#endif


#if TEST_ALLBCAST
#define LOOPCNT 100
void
test_allbcast(THREADED) {
double secs;
int  *arr, i,bytes, maxarr, sz ;
maxarr = 1<<15;
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
      bytes=sz*sizeof(int);
     arr = (int *)node_malloc(bytes,TH);
    on_one_thread{
     for(i=0;i<sz;i++)
      arr[i]=i;
    }
     all_Barrier(TH);
     secs=get_seconds();
        for(i=0;i<LOOPCNT;i++) 
         all_Bcast_i(*arr,TH);
     secs=get_seconds()-secs;
  on_one{
     fprintf(outfile, " Bytes : %7d secs:%5.6f\n",bytes,secs/(double)LOOPCNT);
     fflush(outfile);
        }
     node_free(arr,TH);

  }
}
#endif

#if TEST_ALL_ALLREDUCE
#define LOOPCNT 100   
void
test_all_allreduce(THREADED) {
double secs;
int  *arr, i,bytes, maxarr, sz ;
maxarr = 1<<15;
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
      bytes=sz*sizeof(int);
     arr = (int *)node_malloc(bytes,TH);
    on_one_thread{
     for(i=0;i<sz;i++)
      arr[i]=i;
    }
     all_Barrier(TH);
     secs=get_seconds();
        for(i=0;i<LOOPCNT;i++)
         all_Allreduce_i(*arr,MAX,TH);
     secs=get_seconds()-secs;
  on_one{
     fprintf(outfile, " Bytes : %7d secs:%5.6f\n",bytes,secs/(double)LOOPCNT);
     fflush(outfile); 
        }
     node_free(arr,TH);

  }
}
#endif

#if TEST_ALL_ALLTOALL
#define LOOPCNT 100   
void
test_all_alltoall(THREADED) {
double secs;
int  *arr,*arr1, i,bytes, maxarr, sz ;
maxarr = 1<<15;
   for(sz=1;sz<=maxarr;sz=(sz<<1)) {
      bytes=sz*sizeof(int);
     arr = (int *)node_malloc(bytes,TH);
     arr1 = (int *)node_malloc(bytes,TH);
    on_one_thread{
     for(i=0;i<sz;i++)
      arr[i]=i;
    }
     all_Barrier(TH);
     secs=get_seconds();
        for(i=0;i<LOOPCNT;i++)
         all_Alltoall_i(arr,(sz/NODES),arr1,TH);
     secs=get_seconds()-secs;
  on_one{
     fprintf(outfile, " Bytes : %7d secs: %5.6f\n",bytes,secs/(double)LOOPCNT);
     fflush(outfile);
        }
     node_free(arr,TH);

  }
}
#endif

void
test_ping(int arr_sz) {
  int *a, *b;
  int j, bytes;
  double secs, tsec;
#if TEST_MPI
  MPI_Status stat;
#endif
  
  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(bytes);
  assert_malloc(a);
  b = (int *)malloc(bytes);
  assert_malloc(b);

  pattern2(a, arr_sz);
  pattern2(b, arr_sz);

  UMD_Barrier();

  /*******************************/
  /* TEST Send and Receives      */
  /*******************************/

  secs = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) {
    if ((MYNODE & 0x1)==0) {
      UMD_Send(MYNODE+1,a,bytes);
    }
    else {
      UMD_Recv(MYNODE-1,b,bytes);
      /* check_pattern2(b,arr_sz); */
    }
  }
  secs = get_seconds() - secs;

#if 0
  if ((MYNODE & 0x1)==0) {
    fprintf(outfile,"PE%3d: Send to PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE+1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
  else {
    fprintf(outfile,"PE%3d: Recv fr PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE-1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
#else
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,
	    "PE%3d: UMD  Send: %8d bytes, %9.6lf s: %7.3lf MB/s  %7.3lf Mbps\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    ((double)bytes*(double)LOOP_SEND)/(1000000.0*tsec),
	    (8.0*(double)bytes*(double)LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }
#endif
  
#if TEST_MPI
  /*******************************/
  /* TEST MPI Send and Receives  */
  /*******************************/

  UMD_Barrier();

  bytes = arr_sz*sizeof(int);
  
  secs = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) {
    if ((MYNODE & 0x1)==0) {
      MPI_Send(a, arr_sz, MPI_INT, MYNODE+1, 0, MPI_COMM_WORLD);
    }
    else {
      MPI_Recv(b, arr_sz, MPI_INT, MYNODE-1, 0, MPI_COMM_WORLD, &stat);
      /* check_pattern2(b,arr_sz); */
    }
  }
  secs = get_seconds() - secs;
#if 0  
  if ((MYNODE & 0x1)==0) {
    fprintf(outfile,"PE%3d: MPI  Send to PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE+1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
  else {
    fprintf(outfile,"PE%3d: MPI  Recv fr PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE-1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
#else
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,
	    "PE%3d: MPI  Send: %8d bytes, %9.6lf s: %7.3lf MB/s  %7.3lf Mbps\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    ((double)bytes*(double)LOOP_SEND)/(1000000.0*tsec),
	    (8.0*(double)bytes*(double)LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }
#endif
  
#endif
  
  UMD_Barrier();
  
#if 0
  /*******************************/
  /* TEST Send it back           */
  /*******************************/

  if ((MYNODE & 0x1)==1) 
    UMD_Send(MYNODE-1,a,bytes);
  else {
    UMD_Recv(MYNODE+1,b,bytes);
    check_pattern2(b,arr_sz); 
  }
#endif
  
  free(b);
  free(a);
}

#if _NB_COMM
void
test_ping_nb(int arr_sz) {
  int *a, *b;
  int j, bytes;
  double secs0, tsec0;
  double secs1, tsec1;
  double secs3;
  msg_tag mtag;
  int loop;
  double da, db, dc;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(bytes);
  assert_malloc(a);
  b = (int *)malloc(bytes);
  assert_malloc(b);

  pattern2(a, arr_sz);
  pattern2(b, arr_sz);

  UMD_Barrier();

  /*******************************/
  /* TEST Send and Receives      */
  /*******************************/

  secs0 = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) {
    if ((MYNODE & 0x1)==0) {
      UMD_ISend(MYNODE+1,a,bytes, &mtag);
      dc = 0.0;
      da = 0.9;
      for (loop=0 ; loop<4*arr_sz ; loop++)
	dc += da * (db = (double)loop);
      UMD_Wait(mtag);
    }
    else {
      UMD_Recv(MYNODE-1,b,bytes);
      /* check_pattern2(b,arr_sz); */
    }
  }
  secs0 = get_seconds() - secs0;

  /*******************************/
  /* TEST Send and Receives      */
  /*******************************/

  secs1 = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) {
    if ((MYNODE & 0x1)==0) {
      UMD_Send(MYNODE+1,a,bytes);
      dc = 0.0;
      da = 0.9;
      for (loop=0 ; loop<4*arr_sz ; loop++)
	dc += da * (db = (double)loop);
    }
    else {
      UMD_Recv(MYNODE-1,b,bytes);
      /* check_pattern2(b,arr_sz); */
    }
  }
  secs1 = get_seconds() - secs1;

  /*******************************/
  /* TEST Send and Receives      */
  /*******************************/

  secs3 = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) {
    if ((MYNODE & 0x1)==0) {
      dc = 0.0;
      da = 0.9;
      for (loop=0 ; loop<4*arr_sz ; loop++)
	dc += da * (db = (double)loop);
    }
  }
  secs3 = get_seconds() - secs3;

#if 1
  secs0 -= secs3;
  secs1 -= secs3;
#endif
  
#if 0
  if ((MYNODE & 0x1)==0) {
    fprintf(outfile,"PE%3d: Send to PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE+1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
  else {
    fprintf(outfile,"PE%3d: Recv fr PE%3d: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,MYNODE-1,bytes,
	    secs/(double)LOOP_SEND,
	    (double)(bytes*LOOP_SEND)/(1000000.0*secs));
    fflush(outfile);
  }
#else
  secs0 = max(secs0,MIN_TIME);
  tsec0 = UMD_Reduce_d(secs0,MAX);
  secs1 = max(secs1,MIN_TIME);
  tsec1 = UMD_Reduce_d(secs1,MAX);
  on_one_node {
    fprintf(outfile,
	    "PE%3d: UMD ISend: %8d bytes, %9.6lf s: %7.3lf MB/s  %7.3lf Mbps\n",
	    MYNODE,bytes,
	    tsec0/(double)LOOP_SEND,
	    ((double)bytes*(double)LOOP_SEND)/(1000000.0*tsec0),
	    (8.0*(double)bytes*(double)LOOP_SEND)/(1000000.0*tsec0));
    fprintf(outfile,
	    "PE%3d: UMD  Send: %8d bytes, %9.6lf s: %7.3lf MB/s  %7.3lf Mbps\n",
	    MYNODE,bytes,
	    tsec1/(double)LOOP_SEND,
	    ((double)bytes*(double)LOOP_SEND)/(1000000.0*tsec1),
	    (8.0*(double)bytes*(double)LOOP_SEND)/(1000000.0*tsec1));
    fflush(outfile);
  }
#endif

  UMD_Barrier();
  
#if 0
  /*******************************/
  /* TEST Send it back           */
  /*******************************/

  if ((MYNODE & 0x1)==1) 
    UMD_Send(MYNODE-1,a,bytes);
  else {
    UMD_Recv(MYNODE+1,b,bytes);
    check_pattern2(b,arr_sz); 
  }
#endif
  
  free(b);
  free(a);
}
#endif


void
test_transpose(int arr_sz) {
  int *a, *b;
  int j, k, bytes;
  double secs, tsec;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(bytes);
  assert_malloc(a);
  b = (int *)malloc(bytes);
  assert_malloc(b);

  k = arr_sz/NODES;
    
  pattern4(a, arr_sz);

  UMD_Barrier();

  /*******************************/
  /* TEST Alltoall               */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_TRANSPOSE ; j++)
    UMD_Alltoall_i(a,k,b);
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,
	    "PE%3d: p: %3d UMD  Transpose: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,NODES,bytes,
	    tsec/(double)LOOP_TRANSPOSE,
	    (double)(bytes*LOOP_TRANSPOSE)/(1000000.0*tsec));
    fflush(outfile);
  }

  check_pattern4(b,arr_sz);

  UMD_Barrier();
  
  free(b);
  free(a);
}

#if TEST_TRANSPOSE_S
void
test_transpose_simple(int arr_sz, THREADED) {
  int *a, *b;  
  int j, k, bytes;
  double secs, tsec;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)node_malloc(bytes, TH);
  b = (int *)node_malloc(bytes, TH);
  
  node_Barrier();

  k = arr_sz/NODES;
    
  on_one_thread pattern4(a, arr_sz);

  all_Barrier(TH);

  /*******************************/
  /* TEST Alltoall               */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_TRANSPOSE ; j++)
    all_Alltoall_i(a,k,b, TH);
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs, MAX, TH);
  on_one {
    fprintf(outfile,
	    "PE%3d: p: %3d SIMP Transpose: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,NODES,bytes,
	    tsec/(double)LOOP_TRANSPOSE,
	    (double)(bytes*LOOP_TRANSPOSE)/(1000000.0*tsec));
    fflush(outfile);
  }

  on_one_thread check_pattern4(b,arr_sz);

  all_Barrier(TH);

  node_free(a, TH);
  node_free(b, TH);
}
#endif

#if TEST_MPI
void
test_transpose_mpi(int arr_sz) {
  int *a, *b;
  int j, k, bytes;
  double secs, tsec;
  
#if DEBUG
  fprintf(outfile,"PE%3d: Here 1\n",MYNODE);
  fflush(outfile);
#endif

  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(bytes);
  assert_malloc(a);
  b = (int *)malloc(bytes);
  assert_malloc(b);

  k = arr_sz/NODES;
    
#if DEBUG
  fprintf(outfile,"PE%3d: Here 2\n",MYNODE);
  fflush(outfile);
#endif
  pattern4(a, arr_sz);
#if DEBUG
  fprintf(outfile,"PE%3d: Here 3\n",MYNODE);
  fflush(outfile);
#endif

  /*******************************/
  /* TEST MPI Alltoall           */
  /*******************************/

  MPI_Barrier(MPI_COMM_WORLD);
  
  secs = get_seconds();
#if DEBUG
  fprintf(outfile,"PE%3d: Before loop. a:%p  b:%p  k:%d  arr_sz:%d\n",
	  MYNODE,a,b,k,arr_sz);
  fflush(outfile);
#endif
  for (j=0 ; j<LOOP_TRANSPOSE ; j++) 
    MPI_Alltoall(a, k, MPI_INT, b, k, MPI_INT, MPI_COMM_WORLD);
#if DEBUG
  fprintf(outfile,"PE%3d: After loop. a:%p  b:%p  k:%d  arr_sz:%d\n",
	  MYNODE,a,b,k,arr_sz);
  fflush(outfile);
#endif
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
#if DEBUG
  fprintf(outfile,"PE%3d: Here 4\n",MYNODE);
  fflush(outfile);
#endif
  MPI_Reduce(&secs, &tsec, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  /*  tsec = UMD_Reduce_d(secs,MAX); */
#if DEBUG
  fprintf(outfile,"PE%3d: Here 5\n",MYNODE);
  fflush(outfile);
#endif
  on_one_node {
    fprintf(outfile,
	    "PE%3d: p: %3d MPI  Transpose: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,NODES,bytes,
	    tsec/(double)LOOP_TRANSPOSE,
	    (double)(bytes*LOOP_TRANSPOSE)/(1000000.0*tsec));
    fflush(outfile);
  }
#if DEBUG
  fprintf(outfile,"PE%3d: Here 6\n",MYNODE);
  fflush(outfile);
#endif
  check_pattern4(b,arr_sz);
#if DEBUG
  fprintf(outfile,"PE%3d: Here 7\n",MYNODE);
  fflush(outfile);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  
#if DEBUG
  fprintf(outfile,"PE%3d: Here 8\n",MYNODE);
  fflush(outfile);
#endif
  free(b);
  free(a);
}
#endif


void
test_gather(int arr_sz) {
  int *a, *b;
  int j, bytes;
  double secs, tsec;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(bytes);
  assert_malloc(a);
  b = (int *)malloc(NODES*bytes);
  assert_malloc(b);

  pattern5(a, arr_sz);

  UMD_Barrier();

  /*******************************/
  /* TEST Gather                 */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_SEND ; j++)
    UMD_Gather_i(a,arr_sz,b);
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: UMD  Gather: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }

  on_one_node check_pattern5(b,arr_sz);
  
#if TEST_MPI
  /*******************************/
  /* TEST MPI Gather             */
  /*******************************/

  UMD_Barrier();

  bytes = arr_sz*sizeof(int);
  
  secs = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) 
    MPI_Gather(a, arr_sz, MPI_INT, b, arr_sz, MPI_INT, 0, MPI_COMM_WORLD);
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: MPI  Gather: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }
  on_one_node check_pattern5(b,arr_sz);
#endif
  
  UMD_Barrier();
  
  free(b);
  free(a);
}

#if TEST_GATHER_S
void 
test_gather_simple(int arr_sz, THREADED) {
  int *a, *b;
  int j, bytes;
  double secs, tsec;


  bytes = arr_sz*sizeof(int);

  a = (int *)node_malloc(bytes, TH);
  b = (int *)node_malloc(NODES*bytes, TH);
  on_one_thread {
    pattern5(a, arr_sz);
  }
  
  all_Barrier(TH);

  /*******************************/
  /* TEST Gather                 */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_SEND ; j++)
    all_Gather_i(a,arr_sz,b, TH);
  secs = get_seconds() - secs;


  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs,MAX, TH);


  on_one {
    fprintf(outfile,"PE%3d: SIMP Gather: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }

  on_one check_pattern5(b,arr_sz);
  
  all_Barrier(TH);

  node_free(b, TH);
  node_free(a, TH);
  

}
#endif

void
test_scatter(int arr_sz) {
  int *a, *b;
  int j, bytes;
  double secs, tsec;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)malloc(NODES*bytes);
  assert_malloc(a);
  b = (int *)malloc(bytes);
  assert_malloc(b);

  on_one_node pattern6(a, arr_sz);

  UMD_Barrier();

  /*******************************/
  /* TEST Scatter                */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_SEND ; j++)
    UMD_Scatter_i(a,arr_sz,b);
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: UMD  Scatter: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }

  check_pattern6(b,arr_sz);
  
#if TEST_MPI
  /*******************************/
  /* TEST MPI Scatter            */
  /*******************************/

  UMD_Barrier();

  bytes = arr_sz*sizeof(int);
  
  secs = get_seconds();
  for (j=0 ; j<LOOP_SEND ; j++) 
    MPI_Scatter(a, arr_sz, MPI_INT, b, arr_sz, MPI_INT, 0, MPI_COMM_WORLD);
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: MPI  Scatter: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }
  check_pattern6(b,arr_sz);
#endif
  
  UMD_Barrier();
  
  free(b);
  free(a);
}

#if TEST_SCATTER_S
void
test_scatter_simple(int arr_sz, THREADED) {
  int *a, *b;
  int j, bytes;
  double secs, tsec;
  
  bytes = arr_sz*sizeof(int);

  a = (int *)node_malloc(NODES*bytes, TH);
  b = (int *)node_malloc(bytes, TH);

  on_one pattern6(a, arr_sz);

  all_Barrier(TH);

  /*******************************/
  /* TEST Scatter                */
  /*******************************/

  secs = get_seconds();
  
  for (j=0 ; j<LOOP_SEND ; j++) 
    all_Scatter_i(a,arr_sz,b, TH); 
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(outfile,"PE%3d: SIMP Scatter: (%8d bytes, %9.6lf s): %9.3f MB/s\n",
	    MYNODE,bytes,
	    tsec/(double)LOOP_SEND,
	    (double)((NODES-1)*bytes*LOOP_SEND)/(1000000.0*tsec));
    fflush(outfile);
  }

  on_one_thread check_pattern6(b,arr_sz);

  all_Barrier(TH);

  node_free(b, TH);
  node_free(a, TH);
}
#endif

void
test_barriers() {
  int i;
  double secs, tsec;

  /*******************************/
  /* TEST BARRIERS               */
  /*******************************/

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++)
    UMD_Barrier1();
  secs = get_seconds() - secs;

#if 0
  fprintf(outfile,"PE%3d: UMD  Barrier1: total: %9.6lf  each: %9.6lf\n",
	  MYNODE,secs, secs/(double)LOOP_BARRIER);
  fflush(outfile);
#endif
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs/(double)LOOP_BARRIER, MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: UMD  Barrier1 each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }
  
  /*******************************/
  /* TEST BARRIERS 2             */
  /*******************************/

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++)
    UMD_Barrier2();
  secs = get_seconds() - secs;

#if 0
  fprintf(outfile,"PE%3d: UMD  Barrier2: total: %9.6lf  each: %9.6lf\n",
	  MYNODE,secs, secs/(double)LOOP_BARRIER);
  fflush(outfile);
#endif
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs/(double)LOOP_BARRIER, MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: UMD  Barrier2 each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }

  /*******************************/
  /* TEST BARRIERS 3             */
  /*******************************/

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++)
    UMD_Barrier3();
  secs = get_seconds() - secs;

#if 0
  fprintf(outfile,"PE%3d: UMD  Barrier3: total: %9.6lf  each: %9.6lf\n",
	  MYNODE,secs, secs/(double)LOOP_BARRIER);
  fflush(outfile);
#endif
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs/(double)LOOP_BARRIER, MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: UMD  Barrier3 each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }

#if TEST_MPI
  /*******************************/
  /* TEST BARRIERS MPI           */
  /*******************************/

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++)
    MPI_Barrier(MPI_COMM_WORLD);
  secs = get_seconds() - secs;
#if 0
  fprintf(outfile,"PE%3d: MPI  Barrier:  total: %9.6lf  each: %9.6lf\n",
	  MYNODE,secs, secs/(double)LOOP_BARRIER);
  fflush(outfile);
#endif
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs/(double)LOOP_BARRIER, MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: MPI  Barrier  each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }
#endif

}

void
test_node_barriers(THREADED) {
  int i;
  double secs, tsec;

  /*******************************/
  /* TEST BARRIERS               */
  /*******************************/

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++) 
    node_Barrier_sync();
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = node_Reduce_d(secs/(double)LOOP_BARRIER, MAX, TH);
  on_one_thread {
    fprintf(outfile,"PE%3d: NODE Barrier_sync each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }

  node_Barrier();

  secs = get_seconds();
  for (i=0 ; i<LOOP_BARRIER ; i++) 
    node_Barrier_tree(TH);
  secs = get_seconds() - secs;

  secs = max(secs,MIN_TIME);
  tsec = node_Reduce_d(secs/(double)LOOP_BARRIER, MAX, TH);
  on_one_thread {
    fprintf(outfile,"PE%3d: NODE Barrier_tree each: %9.6lf\n",MYNODE,tsec);
    fflush(outfile);
  }
}  

#if TEST_SORT
void
test_radixsort(int N1) {
  int *inArr, *outArr;
  int q;
#if TIMING
  double secs, tsec;
#endif

  UMD_Barrier();

  q = N1/NODES;

  inArr  = (int *)malloc(q * sizeof(int));
  assert_malloc(inArr);
  outArr = (int *)malloc(q * sizeof(int));
  assert_malloc(outArr);

  create_input_nas(N1, inArr);
#if TIMING
  UMD_Barrier();
  secs = get_seconds();
#endif
  all_radixsort(N1/NODES, inArr, outArr);
#if TIMING
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: n: %12d  FSort: %9.6lf\n",MYNODE,N1,tsec);
    fflush(outfile);
  }
#endif
  all_radixsort_check(N1/NODES,  outArr); 

  UMD_Barrier();
    
  free(outArr);
  free(inArr);
}
#endif

#if TEST_SORT_M
void
test_radixsort_mpi(int N1) {
  int *inArr, *outArr;
  int q;
#if TIMING
  double secs, tsec;
#endif

  MPI_Type_contiguous(2,MPI_INT,&MPI_INTPAIR);
  MPI_Type_commit(&MPI_INTPAIR);
  if (!init_Alltoallv_param)
    init_Alltoallv();

  UMD_Barrier();

  q = N1/NODES;

  inArr  = (int *)malloc(q * sizeof(int));
  assert_malloc(inArr);
  outArr = (int *)malloc(q * sizeof(int));
  assert_malloc(outArr);

  create_input_nas(N1, inArr);
#if TIMING
  UMD_Barrier();
  secs = get_seconds();
#endif
  all_radixsort_mpi(N1/NODES, inArr, outArr);
#if TIMING
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: n: %12d  MSort: %9.6lf\n",MYNODE,N1,tsec);
    fflush(outfile);
  }
#endif
  all_radixsort_check(N1/NODES,  outArr);

  UMD_Barrier();
    
  free(outArr);
  free(inArr);
  destroy_Alltoallv();
  MPI_Type_free(&MPI_INTPAIR);
}
#endif


#if TEST_SORT_S
void
test_radixsort_simple(int N1, THREADED) {
  int *inArr, *outArr;
  int q;
#if TIMING
  double secs, tsec;
#endif

  q = N1/NODES;

  inArr  = (int *)node_malloc(q * sizeof(int), TH);
  outArr = (int *)node_malloc(q * sizeof(int), TH);

#if 1
  create_input_nas_simple(N1, inArr, TH);
#else
  on_one_thread create_input_nas(N1, inArr);
#endif

#if CHECK_NAS_INPUT_CREATION
  all_Barrier(TH);
  {
    int i;
    on_one_thread {
      create_input_nas(N1,outArr);
      for (i=0 ; i<q ; i++)
	if (inArr[i] != outArr[i]) {
	  fprintf(outfile,"PE%3d: ERROR: q: %d (%12d)  F: %12d  S: %12d\n",
		  MYNODE,q,i,outArr[i],inArr[i]);
	  fflush(outfile);
	}
    }
  }
  all_Barrier(TH);
  on_one {
    fprintf(outfile,"Input CHECKS!\n");
    fflush(outfile);
  }
  all_Barrier(TH);
#endif

#if TIMING
  all_Barrier(TH);
  secs = get_seconds();
#endif
  if (NODES==1)
    all_radixsort_node(N1, inArr, outArr, TH);
  else
    all_radixsort_simple(N1/NODES, inArr, outArr, TH);
#if TIMING
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(outfile,"(%3d)PE%3d: T: %3d n: %12d  SSort: %9.6lf\n",
	    NODES,MYNODE,THREADS, N1,tsec);
    fflush(outfile);
  }
#endif

  all_Barrier(TH);

  on_one_thread all_radixsort_check(N1/NODES,  outArr);

  all_Barrier(TH);

  node_free(outArr, TH);
  node_free(inArr, TH);
}
#endif

#if TEST_HISTO
void
test_histo(int arr_sz) {
  int *myA;
  int q;
  double secs, tsec;

  UMD_Barrier();

  q = arr_sz/NODES;
  myA = (int *)malloc(q*sizeof(int));
  assert_malloc(myA);

  create_input_random(q, myA);

  UMD_Barrier();
  secs = get_seconds();

  all_histo(q, myA);

  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: n: %12d  Histo: %9.6lf\n",MYNODE,arr_sz,tsec);
    fflush(outfile);
  }

  UMD_Barrier();
  
  free(myA);
}

void
test_histo_r(int arr_sz, THREADED) {
  int *myA;
  int q;
  double secs, tsec;

  all_Barrier(TH);

  q = arr_sz/TID;
  myA = (int *)malloc(q*sizeof(int));
  assert_malloc(myA);

  create_input_random(q, myA);

  all_Barrier(TH);
  secs = get_seconds();

  all_histo_r(q, myA, TH);

  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(outfile,"PE%3d: n: %12d  Histo_r: %9.6lf\n",MYNODE,arr_sz,tsec);
    fflush(outfile);
  }

  all_Barrier(TH);
  
  free(myA);
}
#endif

#if TEST_SELSUM
void
test_sel_sum(int arr_sz) {
  int *myA;
  int q;
  double secs, tsec;

  UMD_Barrier();

  q = arr_sz/NODES;
  myA = (int *)malloc(q*sizeof(int));
  assert_malloc(myA);

  create_input_random(q, myA);

  UMD_Barrier();
  secs = get_seconds();

  all_sel_sum(q, myA);

  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = UMD_Reduce_d(secs,MAX);
  on_one_node {
    fprintf(outfile,"PE%3d: n: %12d  SelSum: %9.6lf\n",MYNODE,arr_sz,tsec);
    fflush(outfile);
  }

  UMD_Barrier();
  
  free(myA);
}

void
test_sel_sum_r(int arr_sz, THREADED) {
  int *myA, *shA;
  int q;
  double secs, tsec;

  all_Barrier(TH);

  q = arr_sz/TID;

  shA = (int *)node_malloc(q*THREADS*sizeof(int), TH);
  on_one_thread {
    create_input_random(q*THREADS, shA);
  }

  myA = shA + q*MYTHREAD;

  all_Barrier(TH);
  secs = get_seconds();

  all_sel_sum_r(q, myA, TH);

  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = all_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(outfile,"PE%3d: n: %12d  SelSum_r: %9.6lf\n",MYNODE,arr_sz,tsec);
    fflush(outfile);
  }

  all_Barrier(TH);

  node_free(shA, TH);
}
#endif

void *SIMPLE_main(THREADED)
{
  int  i;

#if 0
  fprintf(outfile,"PE%3d(%3d): rint(4.7) -> %d\n",
	  MYNODE,MYTHREAD, (int)rint(4.7));
#endif
  
#if TEST_WIDE
#define TEST_INC (i<<1)
#else
#define TEST_INC (i+=8)
#endif
  
#if DEBUG
  fprintf(outfile,"PE%3d(%3d): SIMPLE_main()\n",MYNODE,MYTHREAD);
  fflush(outfile);
#endif
  
#if TEST_BARRIER
  on_one {
    fprintf(outfile,"PE%3d: LOOP_BARRIER: %d\n",MYNODE, LOOP_BARRIER);
    fflush(outfile);
  }
  on_one_thread test_barriers();
#endif

  all_Barrier(TH);
  
#if TEST_BARRIER_NODE
  on_one {
    fprintf(outfile,"PE%3d: LOOP_BARRIER_NODE: %d\n",MYNODE, LOOP_BARRIER);
    fflush(outfile);
  }
  on_one_node test_node_barriers(TH);
#endif

  all_Barrier(TH);

#if TEST_TEST
  test_test(TH);
#endif
  
  all_Barrier(TH);

#if TEST_TESTNODE
  test_testnode(TH);
#endif
  
  all_Barrier(TH);

#if TEST_NODEBCAST 
  test_nodebcast(TH);
#endif

  all_Barrier(TH);

#if TEST_NODEREDUCE
  test_nodereduce(TH);
#endif
      
  all_Barrier(TH);

#if TEST_ALLREDUCE
  test_allreduce(TH);
#endif

  all_Barrier(TH);

#if TEST_ALLBCAST
  test_allbcast(TH);
#endif
     
  all_Barrier(TH);

#if TEST_ALL_ALLREDUCE
  test_all_allreduce(TH);
#endif
     
  all_Barrier(TH);

#if TEST_ALL_ALLTOALL
  test_all_alltoall(TH);
#endif
     
  all_Barrier(TH);

#if TEST_PING
  if (NODES > 1) {
    on_one {
      fprintf(outfile,"PE%3d: LOOP_PING: %d\n", MYNODE, LOOP_SEND);
      fflush(outfile);
    }
    for (i = 4 ; i<=ARR_SIZE ; i = TEST_INC) 
      on_one_thread test_ping(i);
  }
#endif
  
  all_Barrier(TH);

#if TEST_PING_NB
#if _NB_COMM
  if (NODES > 1) {
    on_one {
      fprintf(outfile,"PE%3d: LOOP_PING: %d\n", MYNODE, LOOP_SEND);
      fflush(outfile);
    }
    for (i = 4 ; i<=ARR_SIZE ; i = TEST_INC) 
      on_one_thread test_ping_nb(i);
  }
#endif
#endif
  
  all_Barrier(TH);

#if TEST_NQUEENS
#if TEST_NQUEENS_R
#if 1
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_R_random_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_R_random_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 1
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_R_block_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_R_block_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 1
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_R_cyclic_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_R_cyclic_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 0
  on_one_thread {
    int x;
    on_one {
      fprintf(outfile,"all_nqueens_R_block()  N: %d\n",NQUEENS_N);
      fflush(outfile);
    }
    x = all_nqueens_R_block(NQUEENS_N);
    on_one {
      fprintf(outfile,"NQUEENS COUNT: %d\n",x);
      fflush(outfile);
    }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_R_cyclic()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_R_cyclic(NQUEENS_N,k); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
      }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_R_random()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_R_random(NQUEENS_N,k); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
      }
  }
  all_Barrier(TH);
#endif


#if 0
  on_one_thread on_node(1) {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) {
      fprintf(outfile,"nqueens_R()  N: %d  k: %d\n",NQUEENS_N,k);
      fflush(outfile);
      nqueens_R(NQUEENS_N,k);
    }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread on_node(1) {
    fprintf(outfile,"nqueens_R_orig()  N: %d\n",NQUEENS_N);
    fflush(outfile);

    nqueens_R_orig(NQUEENS_N);
  }
  all_Barrier(TH);
#endif
#endif

#if TEST_NQUEENS_NR
#if 1
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_NR_random_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_NR_random_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 1
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_NR_block_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_NR_block_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 0
  {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_NR_cyclic_r()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_NR_cyclic_r(NQUEENS_N,k,TH); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
	all_Barrier(TH);
      }
  }
#endif
#if 0
  on_one_thread {
    int x;
    on_one {
      fprintf(outfile,"all_nqueens_NR_block()  N: %d\n",NQUEENS_N);
      fflush(outfile);
    }
    x = all_nqueens_NR_block(NQUEENS_N);
    on_one {
      fprintf(outfile,"NQUEENS COUNT: %d\n",x);
      fflush(outfile);
    }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_NR_cyclic()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_NR_cyclic(NQUEENS_N,k); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
      }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) 
      {
	int x;
	on_one {
	  fprintf(outfile,"all_nqueens_NR_random()  N: %d  k: %d\n",
		  NQUEENS_N,k);
	  fflush(outfile);
	}
	x = all_nqueens_NR_random(NQUEENS_N,k); 
	on_one {
	  fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	  fflush(outfile);
	}
      }
  }
  all_Barrier(TH);
#endif


#if 0
  on_one_thread on_node(1) {
    int k;
    for (k=1 ; k<=NQUEENS_MAXK ; k++) {
      fprintf(outfile,"nqueens_NR()  N: %d  k: %d\n",NQUEENS_N,k);
      fflush(outfile);
      nqueens_NR(NQUEENS_N,k);
    }
  }
  all_Barrier(TH);
#endif
#if 0
  on_one_thread on_node(1) {
    fprintf(outfile,"nqueens_NR_orig()  N: %d\n",NQUEENS_N);
    fflush(outfile);

    nqueens_NR_orig(NQUEENS_N);
  }
  all_Barrier(TH);
#endif
#endif

#if TEST_NQUEENS_X
  {
    int n, x, k=3;
    for (n=8 ; n<20 ; n++) {
      on_one {
	fprintf(outfile,"all_nqueens_R_random_r()  N: %d  k: %d\n",
		n,k);
	fflush(outfile);
      }
      x = all_nqueens_R_random_r(n,k,TH); 
      on_one {
	fprintf(outfile,"NQUEENS COUNT: %d\n",x);
	fflush(outfile);
      }
      all_Barrier(TH);
    }
  }
#endif
  
#endif
  
  all_Barrier(TH);

#if TEST_HISTO
  on_one_thread {
    on_one {
      fprintf(outfile,"all_histo()\n");
      fflush(outfile);
    }
    for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
      test_histo(i);
  }

  all_Barrier(TH);

  on_one {
    fprintf(outfile,"all_histo_r()\n");
    fflush(outfile);
  }
  for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    test_histo_r(i,TH);
#endif

  all_Barrier(TH);

#if TEST_SELSUM
  on_one_thread {
    on_one {
      fprintf(outfile,"all_sel_sum()\n");
      fflush(outfile);
    }
    for (i = (1<<12) ; i<=ARR_SIZE_SELSUM ; i = TEST_INC)
      test_sel_sum(i);
    if (NODES >= 4)
      test_sel_sum((int)pow(10.0,8.0)); 
  }

  all_Barrier(TH);

  on_one {
    fprintf(outfile,"all_sel_sum_r()\n");
    fflush(outfile);
  }
  for (i = (1<<12) ; i<=ARR_SIZE_SELSUM ; i = TEST_INC)
      test_sel_sum_r(i,TH);
  if (NODES >= 4)
    test_sel_sum_r((int)pow(10.0,8.0),TH);
#endif
  
  all_Barrier(TH);

#if TEST_TRANSPOSE
  on_one {
    fprintf(outfile,"PE%3d: LOOP_TRANSPOSE: %d\n", MYNODE, LOOP_TRANSPOSE);
    fflush(outfile);
  }
  for (i = (1<<10) ; i<=ARR_SIZE ; i = TEST_INC) 
    on_one_thread test_transpose(i);
#if TEST_TRANSPOSE_S
  for (i = (1<<10) ; i<=ARR_SIZE ; i = TEST_INC) 
    test_transpose_simple(i, TH);
#endif
  for (i = (1<<10) ; i<=ARR_SIZE ; i = TEST_INC) 
    on_one_thread test_transpose_mpi(i);
#endif
  
  all_Barrier(TH);

#if TEST_GATHER
  on_one {
    fprintf(outfile,"PE%3d: LOOP_GATHER: %d\n", MYNODE, LOOP_SEND);
    fflush(outfile);
  }
  for (i = NODES ; i<=ARR_SIZE ; i = TEST_INC)
    on_one_thread test_gather(i);
#endif
  
  all_Barrier(TH);

#if TEST_GATHER_S
  on_one {
    fprintf(outfile,"PE%3d: LOOP_GATHER: %d\n", MYNODE, LOOP_SEND);
    fflush(outfile);
  }
  for (i = NODES ; i<=ARR_SIZE ; i = TEST_INC)
    test_gather_simple(i, TH);
#endif
  
  all_Barrier(TH);

#if TEST_SCATTER
  on_one {
    fprintf(outfile,"PE%3d: LOOP_SCATTER: %d\n", MYNODE, LOOP_SEND);
    fflush(outfile);
  }
  for (i = NODES ; i<=ARR_SIZE ; i = TEST_INC)
    on_one_thread test_scatter(i);
#endif

  all_Barrier(TH);

#if TEST_SCATTER_S
  on_one {
    fprintf(outfile,"PE%3d: LOOP_SCATTER: %d\n", MYNODE, LOOP_SEND);
    fflush(outfile);
  }
  for (i = NODES ; i<=ARR_SIZE ; i = TEST_INC)
    test_scatter_simple(i, TH);
#endif

  all_Barrier(TH);

#if TEST_SORT
  for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    on_one_thread test_radixsort(i);
#endif

  all_Barrier(TH);

#if TEST_SORT_M
  for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    on_one_thread test_radixsort_mpi(i);
#endif

  all_Barrier(TH);

#if TEST_SORT_S
  for (i = (1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    test_radixsort_simple(i, TH);
#endif

  all_Barrier(TH);

#if TEST_SPLASH_OCEAN
  on_one_node
    ocean_main(TH);
#endif

  all_Barrier(TH);

#if TEST_SPLASH_RADIX
  on_one_node
    all_splash_radix(TH);
#endif

  all_Barrier(TH);

#if TEST_PRIME_M
  on_one_thread
    all_prime_mpi();
#endif

  all_Barrier(TH);

#if TEST_PRIME
  all_prime(TH);
#endif

  all_Barrier(TH);

#if TEST_2DFFT_M
  on_one_thread
    all_2dfft_mpi();
#endif

  all_Barrier(TH);

#if TEST_2DFFT_S
  for (i = (1<<7) ; i<=(MAX_FFT_SIZE) ; i = TEST_INC)
    all_2dfft(i, TH);
#endif

  all_Barrier(TH);

  /*******************************/
  /* End of program              */
  /*******************************/
  
  SIMPLE_done(TH);
}

