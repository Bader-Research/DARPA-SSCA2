/* ---------------------------------------------------------------
 * MPI - Two-Dimensional Fast Fourier Transform - C Version
 * FILE: mpi_2dfft.c
 * OTHER FILES: mpi_2dfft.h
 * DESCRIPTION: The image originates on a single processor (SOURCE_PROCESSOR).
 * This image, a[], is distributed by rows to all other processors.  Each
 * processor then performs a one-dimensional FFT on the rows of the image
 * stored locally.  The image is then transposed using the MPI_Alltoall()
 * routine; this partitions the intermediate image by columns.  Each
 * processor then performs a one-dimensional FFT on the columns of the
 * image.  Finally, the columns of the image are collected back at the
 * destination processor and the output image is tested for correctness.
 *
 * Input is a 512x512 complex matrix. The input matrix is initialized with 
 * a point source. Output is a 512x512 complex matrix that overwrites 
 * the input matrix.  Timing and Mflop results are displayed following 
 * execution.
 *
 * A straightforward unsophisticated 1D FFT kernel is used.  It is
 * sufficient to convey the general idea, but be aware that there are
 * better 1D FFTs available on many systems.
 *
 *  AUTHOR: George Gusciora
 *  LAST REVISED:  06/08/96 Blaise Barney
 * --------------------------------------------------------------- */
#include <stdio.h>
#include <sys/utsname.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "alg_2dfft.h"

#define IMAGE_SIZE		MAX_FFT_SIZE
#define NUM_CELLS		4
#define IMAGE_SLICE		(IMAGE_SIZE / NUM_CELLS)
#define SOURCE_PROCESSOR	0
#define DEST_PROCESSOR		SOURCE_PROCESSOR
#define MAXTIME                 32

#define MPI_2DFFT 0

#define DEBUG   0
#define VERBOSE 0

#define FAST_FFT 1

void fft1(complex_t *, complex_t *, int, int);
void bit_reverse(complex_t *, int);

#if MPI_2DFFT
static complex_t a_mpi[IMAGE_SIZE][IMAGE_SIZE]; /* input matrix: complex numbers */
static complex_t b_mpi[IMAGE_SIZE][IMAGE_SIZE]; /* intermediate matrix */

complex_t a_slice_mpi[IMAGE_SLICE][IMAGE_SIZE];	
complex_t a_chunks_mpi[NUM_CELLS][IMAGE_SLICE][IMAGE_SLICE];
complex_t b_slice_mpi[IMAGE_SIZE][IMAGE_SLICE];
#if !FAST_FFT
complex_t w_common_mpi[IMAGE_SIZE/2];  /* twiddle factors */
#endif
#else
static complex_t a_mpi[1][1]; /* input matrix: complex numbers */
static complex_t b_mpi[1][1]; /* intermediate matrix */
complex_t a_slice_mpi[1][1];	
complex_t a_chunks_mpi[1][1][1];
complex_t b_slice_mpi[1][1];
#if !FAST_FFT
complex_t w_common_mpi[1];  /* twiddle factors */
#endif
#endif

static
void all_print_arr_on_one(char *name,
			  complex_t *arr, int r, int c, int src, THREADED) {

  int i, j;

  on_one fprintf(outfile,"\n");
  fflush(outfile);
  all_Barrier(TH);
  if (ID==src) {
    fprintf(outfile,"(%3d): printing %3d by %3d array (%s):\n",ID,r,c,name);
    for (i=0 ; i<r ; i++) {
      fprintf(outfile,"(%3d): row %3d: ",ID,i);
      for (j=0 ; j<c ; j++) {
	fprintf(outfile,"(%3.0f,%3.0f) ",(arr+i*c+j)->r, (arr+i*c+j)->i);
      }
      fprintf(outfile,"\n");
    }
    fflush(outfile);
  }
  all_Barrier(TH);
}

static
void all_print_arr_on_node(char *name,
			   complex_t *arr, int r, int c, THREADED) {

  int i, j, p;

  for (p=0 ; p<NODES ; p++) {
    on_one fprintf(outfile,"\n");
    fflush(outfile);
    all_Barrier(TH);
    on_one_thread on_node(p) {
      fprintf(outfile,"(%3d): printing %3d by %3d array (%s):\n",ID,r,c,name);
      for (i=0 ; i<r ; i++) {
	fprintf(outfile,"(%3d): row %3d: ",ID,i);
	for (j=0 ; j<c ; j++) {
	  fprintf(outfile,"(%3.0f,%3.0f) ",(arr+i*c+j)->r, (arr+i*c+j)->i);
	}
	fprintf(outfile,"\n");
      }
      fflush(outfile);
    }
    all_Barrier(TH);
  }
}

static
void all_print_arr(char *name,
		   complex_t *arr, int r, int c, THREADED) {

  int i, j, p;

  for (p=0 ; p<TID ; p++) {
    on_one fprintf(outfile,"\n");
    fflush(outfile);
    all_Barrier(TH);
    if (ID==p) {
      fprintf(outfile,"(%3d): printing %3d by %3d array (%s):\n",ID,r,c,name);
      for (i=0 ; i<r ; i++) {
	fprintf(outfile,"(%3d): row %3d: ",ID,i);
	for (j=0 ; j<c ; j++) {
	  fprintf(outfile,"(%3.0f,%3.0f) ",(arr+i*c+j)->r, (arr+i*c+j)->i);
	}
	fprintf(outfile,"\n");
      }
      fflush(outfile);
    }
    all_Barrier(TH);
  }
}

void
all_2dfft_mpi()
{
  int     numtasks;		/* Number of processors */
  int     taskid;	        /* ID number for each processor */
  double    etime[MAXTIME];
  int	    checkpoint;
  double    dt[MAXTIME], sum;

  int cell, i, j, n, nx, logn, errors, sign, flops;
  double mflops;

   checkpoint=0;

   /* Initialize MPI environment and get task's ID and number of tasks 
    * in the partition. */

   numtasks = NODES;
   taskid   = MYNODE;

   if (numtasks != NUM_CELLS)
   {
      fprintf(stderr,
	      "Error: NUM_CELLS is %d, MP_PROCS is %d\n", NUM_CELLS, numtasks);
      /*  exit(1); */
      return;
   }

   n = IMAGE_SIZE;
   /* compute logn and ensure that n is a power of two */
   nx = n;
   logn = 0;
   while(( nx >>= 1) > 0)
      logn++;
   nx = 1;
   for (i=0; i<logn; i++)
      nx = nx*2;
   if (nx != n)
   {
      fprintf(stderr, "%d: fft size must be a power of 2\n", n);
      exit(0);
   }

   if (taskid == SOURCE_PROCESSOR)
   {
      for (i=0; i<n; i++)
         for (j=0; j<n; j++)
            a_mpi[i][j].r = a_mpi[i][j].i = 0.0;
      a_mpi[n/2][n/2].r =  a_mpi[n/2][n/2].i = (double)n; 
   
#if VERBOSE
      /* print table headings in anticipation of timing results */ 
      fprintf(outfile,"512 x 512 2D FFT\n");
      fprintf(outfile,"                                 Timings(secs)\n");
      fprintf(outfile,
	      "          scatter   1D-FFT-row  transpose 1D-FFT-col  gather");
      fprintf(outfile,"        total\n");
#endif
   }

#if !FAST_FFT
   /* precompute the complex constants (twiddle factors) for the 1D FFTs */
   for (i=0;i<n/2;i++)
   {
      w_common_mpi[i].r = (double) cos((double)((2.0*PI*i)/(double)n));
      w_common_mpi[i].i = (double) -sin((double)((2.0*PI*i)/(double)n));
   }
#endif
/*****************************************************************************/
/*                      Distribute Input Matrix By Rows                      */
/*****************************************************************************/
   MPI_Barrier(MPI_COMM_WORLD);

   etime[checkpoint++] = get_seconds();
   MPI_Scatter((char *) a_mpi,       IMAGE_SLICE * n * 2, MPI_DOUBLE,
	       (char *) a_slice_mpi, IMAGE_SLICE * n * 2, MPI_DOUBLE,
	       SOURCE_PROCESSOR, MPI_COMM_WORLD);
   etime[checkpoint++] = get_seconds();

/*****************************************************************************/
/*                      Perform 1-D Row FFTs                                 */
/*****************************************************************************/
   for (i=0;i<IMAGE_SLICE;i++)
#if FAST_FFT
     fft_complex(n, &a_slice_mpi[i][0]);
#else
     fft1(&a_slice_mpi[i][0], w_common_mpi, n, logn);
#endif
   etime[checkpoint++] = get_seconds();
/*****************************************************************************/
/*			Transpose 2-D image				     */
/*****************************************************************************/
   for(cell=0;cell<NUM_CELLS;cell++)
   {
      for(i=0;i<IMAGE_SLICE;i++)
      {
         for(j=0;j<IMAGE_SLICE;j++)
         {
            a_chunks_mpi[cell][i][j].r =
               a_slice_mpi[i][j + (IMAGE_SLICE * cell)].r;
            a_chunks_mpi[cell][i][j].i =
               a_slice_mpi[i][j + (IMAGE_SLICE * cell)].i;
         }
      }
   }

   etime[checkpoint++] = get_seconds();

   MPI_Alltoall(a_chunks_mpi, IMAGE_SLICE * IMAGE_SLICE * 2, MPI_DOUBLE,
		b_slice_mpi, IMAGE_SLICE * IMAGE_SLICE * 2, MPI_DOUBLE,
		MPI_COMM_WORLD);
   etime[checkpoint++] = get_seconds();
/*****************************************************************************/
/*                      Perform 1-D Column FFTs                              */
/*****************************************************************************/
   for(i=0;i<IMAGE_SLICE;i++)
   {
      for(j=0;j<n;j++)
      {
         a_slice_mpi[i][j].r = b_slice_mpi[j][i].r;
         a_slice_mpi[i][j].i = b_slice_mpi[j][i].i;
      }
   }

   etime[checkpoint++] = get_seconds();

   for (i=0;i<IMAGE_SLICE;i++)
#if FAST_FFT
      fft_complex(n, &a_slice_mpi[i][0]);
#else
      fft1(&a_slice_mpi[i][0], w_common_mpi, n, logn);
#endif
   etime[checkpoint++] = get_seconds();
/*****************************************************************************/
/*                    Undistribute Output Matrix by Rows                     */
/*****************************************************************************/
   MPI_Gather(a_slice_mpi, IMAGE_SLICE * n * 2, MPI_DOUBLE,
	      a_mpi,       IMAGE_SLICE * n * 2, MPI_DOUBLE,
	      DEST_PROCESSOR, MPI_COMM_WORLD);

   etime[checkpoint++] = get_seconds();

   if (taskid == DEST_PROCESSOR)
   {
      for(i=0;i<n;i++)
      {
         for(j=0;j<n;j++)
         {
            b_mpi[i][j].r = a_mpi[j][i].r;
            b_mpi[i][j].i = a_mpi[j][i].i;
         }
      }
   }

   etime[checkpoint++] = get_seconds();
   fflush(stdout);

   for (j=0 ; j<numtasks ; j++) {
     MPI_Barrier(MPI_COMM_WORLD);
     if (taskid == j) {
       /* Calculate event timings and flops -  then print them */
       for(i=1;i<checkpoint;i++)
	 dt[i] = etime[i] - etime[i-1];

       sum=0;
       for(i=1;i<checkpoint;i++)
	 sum+=dt[i];
#if VERBOSE
       fprintf(outfile,"cell %d:   ", taskid);
       for(i=1;i<checkpoint;i++)
	 fprintf(outfile,"%2.6f   ", dt[i]);
       fprintf(outfile,"  %2.6f \n", sum);
       fflush(outfile);
#endif
     }
     MPI_Barrier(MPI_COMM_WORLD);
   }

   for(i=1 ; i<checkpoint ; i++) {
     MPI_Reduce(&dt[i], &dt[0], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     if (taskid == 0) {
       fprintf(outfile,"P: %2d n: %12d 2DMPIFFT Step %2d: %9.6f\n",
	       numtasks, n, i, dt[0]);
       fflush(outfile);
     }
   }
   MPI_Barrier(MPI_COMM_WORLD);
    

   if (taskid == DEST_PROCESSOR)
   {
     flops = (n*n*logn)*10; 
     mflops = ((double)flops/1000000.0);
     mflops = mflops/(double)sum;
#if VERBOSE
     fprintf(outfile,"Total Mflops= %3.4f\n", mflops);
     fflush(outfile);
#endif

      errors = 0;
      for (i=0;i<n;i++)
      {
         if (((i+1)/2)*2 == i)
            sign = 1;
         else
            sign = -1;
         for (j=0;j<n;j++)
         {
            if (b_mpi[i][j].r > n*sign+EPSILON ||
                b_mpi[i][j].r < n*sign-EPSILON ||
                b_mpi[i][j].i > n*sign+EPSILON ||
                b_mpi[i][j].i < n*sign-EPSILON)
            {
               fprintf(outfile,"[%d][%d] is %f,%f should be %f\n", i, j,
		       b_mpi[i][j].r, b_mpi[i][j].i, (double) n*sign);
               errors++;
            }
            sign *= -1;
         }
      }
      if (errors)
      {
	fprintf(outfile,"%d errors!!!!!\n", errors);
	exit(0);
      }
   }
#if 0
   if (taskid == DEST_PROCESSOR)
   {
      for(i=0;i<n;i++)
      {
         printf("Row %d\n",i);
         for(j=0;j<n;j++)
            printf("%f %f\n", a_mpi[i][j].r, a_mpi[i][j].i);
      }
   }
#endif

   return;
}


void
all_2dfft(int points, THREADED)
{
  int     numtasks;		/* Number of processors */
  int     taskid;	        /* ID number for each processor */

  complex_t *a_MEM;
  complex_t *b_MEM;
  complex_t *a_slice_MEM;
  complex_t *a_slice;
  complex_t *a_chunks_MEM;
  /*  complex_t *a_chunks; */
  complex_t *b_slice_MEM;
  complex_t *b_slice;
#if !FAST_FFT
  complex_t *w_common;
#endif
  int a_i0;
  int ac_i0;
  int ac_i1;
  int ac_off;
  int as_i0;
  int bs_i0;
  int num_cells, image_slice;
  int slicesize;

  double    etime[MAXTIME];
  int	    checkpoint;
  double    dt[MAXTIME], sum;

  int cell, i, j, n, nx, logn, errors, sign, flops;
  int x, ix, jx;
  int blksz;
  double mflops;

  checkpoint=0;

  if (points <= 0)
    n           = IMAGE_SIZE;
  else
    n           = points;
  num_cells   = TID;
  image_slice = n / num_cells;
  slicesize   = n * image_slice;

  a_MEM = (complex_t *)node_malloc(n * n * sizeof(complex_t), TH);
  b_MEM = (complex_t *)node_malloc(n * n * sizeof(complex_t), TH);

  a_slice_MEM  = (complex_t *)node_malloc(
    THREADS * slicesize * sizeof(complex_t), TH);
  a_chunks_MEM = (complex_t *)node_malloc(
    THREADS * slicesize * sizeof(complex_t), TH);
  b_slice_MEM  = (complex_t *)node_malloc(
    THREADS * slicesize * sizeof(complex_t), TH);

  a_slice  = a_slice_MEM  + (MYTHREAD * slicesize);
  /* a_chunks = a_chunks_MEM + (MYTHREAD * slicesize);  */
  b_slice  = b_slice_MEM  + (MYTHREAD * slicesize);

#if !FAST_FFT
  w_common = (complex_t *)malloc((n/2)*sizeof(complex_t));
  assert_malloc(w_common);
#endif

  a_i0   = n;
  ac_i0  = image_slice * image_slice * THREADS;
  ac_i1  = image_slice;
  ac_off = image_slice * image_slice * MYTHREAD;
  as_i0  = n;
  bs_i0  = image_slice;
				 
  numtasks = TID;
  taskid   = ID;

  /* compute logn and ensure that n is a power of two */
  nx = n;
  logn = 0;
  while(( nx >>= 1) > 0)
    logn++;
  nx = 1;
  for (i=0; i<logn; i++)
    nx = nx*2;
  if (nx != n) {
    fprintf(stderr, "%d: fft size must be a power of 2\n", n);
    exit(0);
  }

  blksz = ((n*n)/NODES) * sizeof(complex_t);

  if (taskid == SOURCE_PROCESSOR) {
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
	(a_MEM + i*a_i0 + j)->r = (a_MEM + i*a_i0 + j)->i = 0.0;
    (a_MEM + ((n/2)*a_i0) + n/2)->r =
      (a_MEM + ((n/2)*a_i0) + n/2)->i = (double)n; 

#if DEBUG
    fprintf(outfile,"points: %d\n", points);
    fprintf(outfile,"IMAGE_SIZE: %d\n", IMAGE_SIZE);
    fprintf(outfile,"image_slice: %d\n", image_slice);
    fprintf(outfile,"num_cells: %d\n", num_cells);
    fprintf(outfile,"n: %d\n", n);
    fprintf(outfile,"blksz: %d\n", blksz);
    fprintf(outfile,"logn: %d\n", logn);
    fprintf(outfile,"sizeof(complex_t): %d\n", sizeof(complex_t));
    fprintf(outfile,"TID: %d\n", TID);
    fprintf(outfile,"numtasks: %d\n", numtasks);
#endif
    
#if VERBOSE
    /* print table headings in anticipation of timing results */
    fprintf(outfile,"%3d x %3d 2D FFT\n",n,n);
    fprintf(outfile,"                                 Timings(secs)\n");
    fprintf(outfile,
	    "          scatter   1D-FFT-row  transpose 1D-FFT-col  gather");
    fprintf(outfile,"        total\n");
#endif

#if DEBUG
    for (i=0; i<n; i++)
      for (j=0; j<n; j++) {
	if ((((a_MEM + (i*a_i0) + j)->r != (a_MEM + (i*a_i0) + j)->i) ||
	     ((a_MEM + (i*a_i0) + j)->r != 0.0)) &&
	    (i != n/2) && (j != n/2))
	  fprintf(outfile,"(%3d): Error 0 i: %3d j: %3d  (%f, %f)\n",
		  ID, i, j,
		  (a_MEM + (i*a_i0) + j)->r,
		  (a_MEM + (i*a_i0) + j)->i);
      }
    if ((a_MEM + ((n/2)*a_i0) + n/2)->r != (double)n)
      fprintf(outfile,"(%3d): Error 1\n",ID);
    if ((a_MEM + ((n/2)*a_i0) + n/2)->i != (double)n)
      fprintf(outfile,"(%3d): Error 2\n",ID);
#endif

  }
  all_Barrier(TH);
    
#if !FAST_FFT  
  /* precompute the complex constants (twiddle factors) for the 1D FFTs */
  for (i=0;i<n/2;i++) {
    w_common[i].r = (double) cos((double)((2.0*PI*i)/(double)n));
    w_common[i].i = (double) -sin((double)((2.0*PI*i)/(double)n));
  }
#endif
/*****************************************************************************/
/*                      Distribute Input Matrix By Rows                      */
/*****************************************************************************/
#if DEBUG
  all_print_arr_on_one("a_MEM",a_MEM, n, n, 0, TH);
#endif
  all_Barrier(TH);
  etime[checkpoint++] = get_seconds();
  all_Scatter_c((char *)a_MEM, blksz, (char *)a_slice_MEM, TH); 
  all_Barrier(TH);
  etime[checkpoint++] = get_seconds();
#if DEBUG
  all_print_arr("a_slice after Scatter",a_slice, image_slice, n, TH);
#endif
/*****************************************************************************/
/*                      Perform 1-D Row FFTs                                 */
/*****************************************************************************/
  for (i=0 ; i<image_slice ; i++) {
#if FAST_FFT
    fft_complex(n, a_slice + (i*as_i0)); 
#else
    fft1(a_slice + (i*as_i0), w_common, n, logn);
#endif
  }
  etime[checkpoint++] = get_seconds();
  node_Barrier();
#if DEBUG
  all_print_arr("a_slice after FFT",a_slice, image_slice, n, TH);
#endif
/*****************************************************************************/
/*			Transpose 2-D image				     */
/*****************************************************************************/
  for (cell=0 ; cell<num_cells ; cell++) {
    for (i=0 ; i<image_slice ; i++) {
      memcpy(a_chunks_MEM + (cell*ac_i0) + (i*ac_i1) + ac_off,   /* dst */
	     a_slice  + (i*as_i0)    + (image_slice * cell), /* src */
	     image_slice*sizeof(complex_t));
    }
  }

  etime[checkpoint++] = get_seconds();
  node_Barrier();

#if DEBUG
  all_print_arr_on_node("a_chunks_MEM",a_chunks_MEM, num_cells, ac_i0, TH);
#endif

#if 0
  all_Alltoall_c((char *)a_chunks_MEM,
		 blksz/NODES,
		 (char *)b_slice_MEM, TH);
  all_Barrier(TH);
#else
  on_one_thread
    UMD_Alltoall_d((double *)a_chunks_MEM,
		   ((n/NODES)*(n/NODES)<<1),
		   (double *)b_slice_MEM);
  all_Barrier(TH);
#endif
  
  etime[checkpoint++] = get_seconds();
#if DEBUG
  all_print_arr("b_slice after Alltoall",b_slice, image_slice, n, TH);
#endif
/*****************************************************************************/
/*                      Perform 1-D Column FFTs                              */
/*****************************************************************************/
#if 0
  /* works when p=1 or r=1 */
  for (i=0 ; i<image_slice ; i++) {
    for (j=0 ; j<n ; j++) {
      (a_slice + (i*as_i0) + j)->r = (b_slice + (j*bs_i0) + i)->r;
      (a_slice + (i*as_i0) + j)->i = (b_slice + (j*bs_i0) + i)->i;
    }
  }
#endif
#if 0
  for (i=0 ; i<NODES ; i++) {
    memcpy(a_chunks + i*ac_i0,
	   b_slice_MEM + (i*THREADS + MYTHREAD)*ac_i0,
	   ac_i0*sizeof(complex_t));
  }
  for (i=0 ; i<image_slice ; i++) {
    for (j=0 ; j<n ; j++) {
      (a_slice + (i*as_i0) + j)->r = (a_chunks + (j*bs_i0) + i)->r;
      (a_slice + (i*as_i0) + j)->i = (a_chunks + (j*bs_i0) + i)->i;
    }
  }
#endif
#if 1
  for (i=0 ; i<image_slice ; i++) {
    for (j=0 ; j<n ; j++) {
      x  = (j*bs_i0)+i;
      ix = x / ac_i0;
      jx = x % ac_i0;
      (a_slice + (i*as_i0) + j)->r =
	(b_slice_MEM + (ix*THREADS + MYTHREAD)*ac_i0 + jx)->r;
      (a_slice + (i*as_i0) + j)->i = 
	(b_slice_MEM + (ix*THREADS + MYTHREAD)*ac_i0 + jx)->i;
    }
  }
#endif

  node_Barrier();
  etime[checkpoint++] = get_seconds();

#if DEBUG
  all_print_arr("a_slice after rearrange",a_slice, image_slice, n, TH);
#endif
  for (i=0 ; i<image_slice ; i++) {
#if FAST_FFT
    fft_complex(n, a_slice + (i*as_i0));
#else
    fft1(a_slice + (i*as_i0), w_common, n, logn); 
#endif
  }
  etime[checkpoint++] = get_seconds();
/*****************************************************************************/
/*                    Undistribute Output Matrix by Rows                     */
/*****************************************************************************/
  node_Barrier();
#if DEBUG
  all_print_arr("a_slice after FFT",a_slice, image_slice, n, TH);
#endif
  
  all_Gather_c((char *)a_slice_MEM, blksz, (char *)a_MEM, TH);
  all_Barrier(TH);

  etime[checkpoint++] = get_seconds();

  if (taskid == DEST_PROCESSOR) {
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	(b_MEM + (i*a_i0) + j)->r = (a_MEM + (j*a_i0) + i)->r;
	(b_MEM + (i*a_i0) + j)->i = (a_MEM + (j*a_i0) + i)->i;
      }
    }
  }

  etime[checkpoint++] = get_seconds();
  fflush(stdout);

  for (j=0 ; j < numtasks ; j++) {
    all_Barrier(TH);
    if (taskid == j) {
      /* Calculate event timings and flops -  then print them */
      for(i=1;i<checkpoint;i++)
	dt[i] = etime[i] - etime[i-1];

      sum=0;
      for(i=1;i<checkpoint;i++)
	sum+=dt[i];
#if VERBOSE
      fprintf(outfile,"cell %d:   ", taskid);
      for(i=1;i<checkpoint;i++)
	fprintf(outfile,"%2.6f   ", dt[i]);
      fprintf(outfile,"  %2.6f \n", sum);
      fflush(outfile);
#endif
      fflush(outfile);
    }
    all_Barrier(TH);
  }

  for(i=1 ; i<checkpoint ; i++) {
    dt[0] = all_Reduce_d(dt[i], MAX, TH);
    on_one {
      fprintf(outfile,"P: %2d R: %2d n: %12d 2DFFT Step %2d: %9.6f\n",
	      NODES, THREADS, n, i, dt[0]);
      fflush(outfile);
    }
  }
  all_Barrier(TH);
    

  if (taskid == DEST_PROCESSOR) {

    flops = (n*n*logn)*10; 
    mflops = ((double)flops/1000000.0);
    mflops = mflops/(double)sum;
#if VERBOSE
    fprintf(outfile,"Total Mflops= %3.4f\n", mflops);
#endif

    errors = 0;
    for (i=0;i<n;i++) {
      if (((i+1)/2)*2 == i)
	sign = 1;
      else
	sign = -1;
      for (j=0;j<n;j++) {
	if ((b_MEM + (i*a_i0) + j)->r > n*sign+EPSILON ||
	    (b_MEM + (i*a_i0) + j)->r < n*sign-EPSILON ||
	    (b_MEM + (i*a_i0) + j)->i > n*sign+EPSILON ||
	    (b_MEM + (i*a_i0) + j)->i < n*sign-EPSILON) {
	  fprintf(outfile,"[%d][%d] is %f,%f should be %f\n", i, j,
		  (b_MEM + (i*a_i0) + j)->r,
		  (b_MEM + (i*a_i0) + j)->i,
		  (double) n*sign);
#if DEBUG
	  exit(0);
#endif
	  errors++;
	}
	sign *= -1;
      }
    }
    if (errors) {
      fprintf(outfile,"%d errors!!!!!\n", errors);
      exit(0);
    }
  }

#if 0
  if (taskid == DEST_PROCESSOR) {
    for(i=0;i<n;i++) {
      fprintf(outfile,"Row %d\n",i);
      for(j=0;j<n;j++)
	fprintf(outfile,"%f %f\n",
		(a_MEM+(i*a_i0)+j)->r,
		(a_MEM+(i*a_i0)+j)->i);
    }
  }
#endif

#if !FAST_FFT
  free(w_common);
#endif
  node_free(b_slice_MEM,  TH);
  node_free(a_chunks_MEM, TH);
  node_free(a_slice_MEM,  TH);
  node_free(b_MEM,        TH);
  node_free(a_MEM,        TH);
    
  return;
}

void fft1(complex_t *data, complex_t *w_c, int n, int logn)
{
  int incrvec, i0, i1, i2;
  double f0, f1;

  /* bit-reverse the input vector */
  bit_reverse(data,n);

  /* do the first logn-1 stages of the fft */
  i2 = logn;
  for (incrvec=2;incrvec<n;incrvec<<=1) {
    i2--;
    for (i0 = 0; i0 < incrvec >> 1; i0++) {
      for (i1 = 0; i1 < n; i1 += incrvec) {
        f0 = data[i0+i1 + incrvec/2].r * w_c[i0<<i2].r - 
	  data[i0+i1 + incrvec/2].i * w_c[i0<<i2].i;
        f1 = data[i0+i1 + incrvec/2].r * w_c[i0<<i2].i + 
	  data[i0+i1 + incrvec/2].i * w_c[i0<<i2].r;
        data[i0+i1 + incrvec/2].r = data[i0+i1].r - f0;
        data[i0+i1 + incrvec/2].i = data[i0+i1].i - f1;
        data[i0+i1].r = data[i0+i1].r + f0;
        data[i0+i1].i = data[i0+i1].i + f1;
      }
    }
  }

  /* do the last stage of the fft */
  for (i0 = 0; i0 < n/2; i0++) {
    f0 = data[i0 + n/2].r * w_c[i0].r - 
      data[i0 + n/2].i * w_c[i0].i;
    f1 = data[i0 + n/2].r * w_c[i0].i + 
      data[i0 + n/2].i * w_c[i0].r;
    data[i0 + n/2].r = data[i0].r - f0;
    data[i0 + n/2].i = data[i0].i - f1;
    data[i0].r = data[i0].r + f0;
    data[i0].i = data[i0].i + f1;
  }
}

/* 
 * bit_reverse - simple (but somewhat inefficient) bit reverse 
 */
void bit_reverse(complex_t *a, int n)
{
  int i,j,k;

  j = 0;
  for (i=0; i<n-2; i++){
    if (i < j)
      CMPLXSWAP(a[j],a[i]);
    k = n>>1;
    while (k <= j) {
      j -= k; 
      k >>= 1;
    }
    j += k;
  }
}

