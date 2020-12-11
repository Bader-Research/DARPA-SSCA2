#include "im.h"
#include <math.h>
#include <string.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define DEBUG   0
#define RESULTS 0
#define NIL    -1
#define SNFLIMIT 200

#define border_pix(ngrid,i,j)                                  \
(                                                              \
   (((ngrid)->myRow==           0) && ((i)==            0)) || \
   (((ngrid)->myCol==           0) && ((j)==            0)) || \
   (((ngrid)->myCol==(ngrid)->c-1) && ((j)==(ngrid)->ny-1)) || \
   (((ngrid)->myRow==(ngrid)->r-1) && ((i)==(ngrid)->nx-1))    \
) 

void
alloc_image_queue(image_queue_t *im_queue) {
  *im_queue = (image_queue_t)malloc(sizeof(struct image_queue));
  assert_malloc(*im_queue);
  return;
}

void
free_image_queue(image_queue_t im_queue) {
  free(im_queue->desc);
  free(im_queue);
  return;
}

void
alloc_image_queue_desc(image_queue_t im_queue, int num) {
  im_queue->desc = (image_desc_t)malloc(num * sizeof(struct image_desc));
  assert_malloc(im_queue->desc);
  return;
}

void
all_GetParams(image_queue_t im_queue, THREADED) {
  int i, argc;
  argc = THARGC;
  if (argc < 1)
    fprintf(stderr,"ERROR: argc < 1 (%d)\n",argc);
#if DEBUG
  on_one {
    int i;
    fprintf(outfile,"argc: %d\n",argc);
    for (i=0 ; i<argc ; i++)
      fprintf(outfile,"argv[%2d]: %s\n",i, THARGV[i]);
  }
#endif
  im_queue->num  = argc;
  alloc_image_queue_desc(im_queue, im_queue->num);
  for (i=0 ; i<im_queue->num ; i++) {
    strcpy(im_queue->desc[i].fname,THARGV[i]);
  }
  return;
}

void
all_MapNodes(grid_t ngrid, THREADED) {
  ngrid->r     = (int)ceil(sqrt((double)NODES));
  ngrid->c     = NODES / ngrid->r;
  ngrid->myRow = MYNODE / ngrid->c;   
  ngrid->myCol = MYNODE % ngrid->c;

  ngrid->border_tp = (ngrid->myRow == 0);
  ngrid->border_bm = (ngrid->myRow == (ngrid->r-1));
  ngrid->border_lt = (ngrid->myCol == 0);
  ngrid->border_rt = (ngrid->myCol == (ngrid->c-1));

  ngrid->ul = ((ngrid->border_tp || ngrid->border_lt) &&
	       !(ngrid->border_bm || ngrid->border_rt));
  ngrid->ur = ((ngrid->border_tp || ngrid->border_rt) &&
	       !(ngrid->border_bm || ngrid->border_lt));
  ngrid->ll = ((ngrid->border_bm || ngrid->border_lt) &&
	       !(ngrid->border_tp || ngrid->border_rt));
  ngrid->lr = ((ngrid->border_bm || ngrid->border_rt) &&
	       !(ngrid->border_tp || ngrid->border_lt));
  ngrid->ce = !(ngrid->border_tp || ngrid->border_bm ||
		ngrid->border_lt || ngrid->border_rt);

#if DEBUG
  on_one {
    fprintf(outfile,"(MapNodes) r: %3d  c: %3d\n",
	    ngrid->r, ngrid->c);
    fflush(outfile);
  }
#endif
  return;
}

void
all_MapImage(image_desc_t im_desc, grid_t ngrid) {
  ngrid->x  = im_desc->x;
  ngrid->y  = im_desc->y;
  ngrid->nx = ngrid->x / ngrid->r;
  ngrid->ny = ngrid->y / ngrid->c;
#if DEBUG
  on_one {
    fprintf(outfile,"(MapImage) nx: %3d  ny: %3d\n",
	    ngrid->nx, ngrid->ny);
    fflush(outfile);
  }
#endif
  return;
}

void
all_AllocateImage(grid_t ngrid, image_t *X, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(AllocateImage)\n");
    fflush(outfile);
  }
#endif
  
  *X = (int *)node_malloc(ngrid->nx * ngrid->ny * sizeof(int), TH);
  return;
}

void
all_FreeImage(image_t X, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(FreeImage)\n");
    fflush(outfile);
  }
#endif
  node_free(X, TH);
  return;
}

void
all_ImageCopy(grid_t ngrid, image_t dst, image_t src, THREADED) {
  int i, blk;
#if DEBUG
  on_one {
    fprintf(outfile,"(ImageCopy)\n");
    fflush(outfile);
  }
#endif
  
  node_Barrier();

  blk= ngrid->ny * sizeof(int);
  
  node_pardo(i, 0, ngrid->nx, 1)
    memcpy(dst+(i*ngrid->ny), src+(i*ngrid->ny), blk);
  
  node_Barrier();
  return;
}

void
all_LoadImage_orig(image_desc_t im_info, grid_t ngrid, image_t *X, THREADED) {
  register int
    i, j,          /* tile position */
    ii, jj;        /* logical grid position */

  FILE* infile;
  unsigned char anum;
  int *rowbar, rowbar_sz;
  int src, dst;
  int *rowp;

#if DEBUG
  on_one {
    fprintf(outfile,"(LoadImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif

  on_one {

    /* OPEN INPUT IMAGE */
    infile = fopen(im_info->fname,"r");
    if (infile == NULL) errprnt("input file does not exist.");
    rewind(infile);

    /* GET IMAGE DESCRIPTION */
    fscanf(infile,"%s\n",im_info->format);
    fscanf(infile,"%s\n",im_info->comment);
    fscanf(infile,"%d",  &(im_info->x));
    fscanf(infile,"%d\n",&(im_info->y));
    fscanf(infile,"%d\n",&(im_info->k));
  }

  im_info->x = all_Bcast_i(im_info->x, TH);
  im_info->y = all_Bcast_i(im_info->y, TH);
  im_info->k = all_Bcast_i(im_info->k, TH);
  
  all_MapImage(im_info, ngrid);
  all_AllocateImage(ngrid, X, TH);

  rowbar_sz = ngrid->ny*sizeof(int);
  rowbar = (int *)malloc(rowbar_sz);
  assert_malloc(rowbar);

  all_Barrier(TH);

  rowp = *X;

  on_one {

    for (ii = 0 ; ii<ngrid->r ; ii++)
      for (i=0 ; i<ngrid->nx ; i++) 
	for (jj = 0 ; jj<ngrid->c ; jj++) {
	  for (j=0 ; j<ngrid->ny ; j++) {
	    fscanf(infile,"%c",&anum);
	    rowbar[j] = (int)anum;
	  }
	  dst = (ii*ngrid->c) + jj;
	  if (dst == MYNODE) {
	    memcpy(rowp,rowbar,rowbar_sz);
	    rowp += ngrid->ny;
	  }
	  else
	    UMD_Send(dst, rowbar, rowbar_sz);
	}

    fclose(infile);
  }

  if (MYNODE != 0)
    on_one_thread {
      src = 0;
      for (i=0 ; i<ngrid->nx ; i++) {
	UMD_Recv(src, rowp, rowbar_sz);
	rowp += ngrid->ny;
      }
  }

  all_Barrier(TH);
  free(rowbar);
  return;
}


void
all_LoadImage(image_desc_t im_info, grid_t ngrid, image_t *X, THREADED) {
  register int
    i, j;          /* tile position */

  FILE* infile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset;

#if DEBUG
  on_one {
    fprintf(outfile,"(LoadImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif

  /* OPEN INPUT IMAGE */
  infile = fopen(im_info->fname,"r");
  if (infile == NULL) errprnt("input file does not exist.");
  rewind(infile);

  /* GET IMAGE DESCRIPTION */
  fscanf(infile,"%s\n",im_info->format);
  fscanf(infile,"%s\n",im_info->comment);
  fscanf(infile,"%d",  &(im_info->x));
  fscanf(infile,"%d\n",&(im_info->y));
  fscanf(infile,"%d\n",&(im_info->k));
  pre_offset = (int)ftell(infile);

  all_MapImage(im_info, ngrid);
  all_AllocateImage(ngrid, X, TH);

  all_Barrier(TH);

  rowp = *X + MYTHREAD * (ngrid->nx * ngrid->ny / THREADS);

  /* base_offest: number of pixels before my first */
  base_offset = pre_offset + ngrid->myCol * ngrid->ny +
    (ngrid->myRow * ngrid->nx * ngrid->y);
  
  node_pardo(i, 0, ngrid->nx, 1) {
    fseek(infile, base_offset + (i*ngrid->y), SEEK_SET);
    for (j=0 ; j<ngrid->ny ; j++) {
      fscanf(infile,"%c",&anum);
      *rowp = (int)anum;
      rowp++;
    }
  }

  fclose(infile);

  all_Barrier(TH);
  return;
}


void
all_LoadImage_unroll(image_desc_t im_info, grid_t ngrid, image_t *X, THREADED) {
  register int
    i, j;          /* tile position */

  FILE* infile;
  unsigned char anum[8];
  int *rowp;
  int base_offset, pre_offset;

#if DEBUG
  on_one {
    fprintf(outfile,"(LoadImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif

  /* OPEN INPUT IMAGE */
  infile = fopen(im_info->fname,"r");
  if (infile == NULL) errprnt("input file does not exist.");
  rewind(infile);

  /* GET IMAGE DESCRIPTION */
  fscanf(infile,"%s\n",im_info->format);
  fscanf(infile,"%s\n",im_info->comment);
  fscanf(infile,"%d",  &(im_info->x));
  fscanf(infile,"%d\n",&(im_info->y));
  fscanf(infile,"%d\n",&(im_info->k));
  pre_offset = (int)ftell(infile);

  all_MapImage(im_info, ngrid);
  all_AllocateImage(ngrid, X, TH);

  all_Barrier(TH);

  rowp = *X + MYTHREAD * (ngrid->nx * ngrid->ny / THREADS);

  /* base_offest: number of pixels before my first */
  base_offset = pre_offset + ngrid->myCol * ngrid->ny +
    (ngrid->myRow * ngrid->nx * ngrid->y);
  
  node_pardo(i, 0, ngrid->nx, 1) {
    fseek(infile, base_offset + (i*ngrid->y), SEEK_SET);
    for (j=0 ; j<ngrid->ny ; j+=8) {
      fscanf(infile,"%c%c%c%c%c%c%c%c",
	     anum,   anum+1, anum+2, anum+3,
	     anum+4, anum+5, anum+6, anum+7);
      *rowp++ = (int)*anum;
      *rowp++ = (int)*(anum+1);
      *rowp++ = (int)*(anum+2);
      *rowp++ = (int)*(anum+3);
      *rowp++ = (int)*(anum+4);
      *rowp++ = (int)*(anum+5);
      *rowp++ = (int)*(anum+6);
      *rowp++ = (int)*(anum+7);
    }
  }

  fclose(infile);

  all_Barrier(TH);
  return;
}


void
all_LoadImage_flock(image_desc_t im_info, grid_t ngrid, image_t *X, THREADED) {
  register int
    i, j;          /* tile position */

  FILE* infile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset;
  struct flock alock;

#if DEBUG
  on_one {
    fprintf(outfile,"(LoadImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif

  /* OPEN INPUT IMAGE */
  infile = fopen(im_info->fname,"r");
  if (infile == NULL) errprnt("input file does not exist.");
  rewind(infile);

  /* GET IMAGE DESCRIPTION */
  fscanf(infile,"%s\n",im_info->format);
  fscanf(infile,"%s\n",im_info->comment);
  fscanf(infile,"%d",  &(im_info->x));
  fscanf(infile,"%d\n",&(im_info->y));
  fscanf(infile,"%d\n",&(im_info->k));
  pre_offset = (int)ftell(infile);

  all_MapImage(im_info, ngrid);
  all_AllocateImage(ngrid, X, TH);

  all_Barrier(TH);

  rowp = *X + MYTHREAD * (ngrid->nx * ngrid->ny / THREADS);

  /* base_offest: number of pixels before my first */
  base_offset = pre_offset + ngrid->myCol * ngrid->ny +
    (ngrid->myRow * ngrid->nx * ngrid->y);
  
  alock.l_type   = F_RDLCK;
  alock.l_whence = SEEK_SET;
  alock.l_start  = 0;
  alock.l_len    = 0;
  alock.l_pid    = getpid();
  fcntl(infile, F_SETLK, &alock);

  node_pardo(i, 0, ngrid->nx, 1) {
    fseek(infile, base_offset + (i*ngrid->y), SEEK_SET);
    for (j=0 ; j<ngrid->ny ; j++) {
      fscanf(infile,"%c",&anum);
      *rowp = (int)anum;
      rowp++;
    }
  }

  alock.l_type   = F_UNLCK;
  fcntl(infile, F_SETLK, &alock);
  fclose(infile);

  all_Barrier(TH);
  return;
}


void
all_SaveImage_orig(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j,          /* tile position */
    ii, jj;        /* logical grid position */

  FILE* savefile;
  unsigned char anum;
  int *rowbar, rowbar_sz;
  int src, dst;
  int *rowp;

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  rowbar_sz = ngrid->ny*sizeof(int);
  rowbar = (int *)malloc(rowbar_sz);
  assert_malloc(rowbar);

  all_Barrier(TH);

  rowp = X;

  on_one {

    savefile = fopen(im_info->fname,"w");
    if (savefile == NULL) errprnt("output file cannot be written.");
    rewind(savefile);

    /* HEADER */
    fprintf(savefile,"%s\n",im_info->format);
    fprintf(savefile,"%s\n",im_info->comment);
    fprintf(savefile,"%d %d\n",im_info->x, im_info->y);
    fprintf(savefile,"%d\n",im_info->k);

    for (ii = 0 ; ii<ngrid->r ; ii++)
      for (i=0 ; i<ngrid->nx ; i++) 
	for (jj = 0 ; jj<ngrid->c ; jj++) {
	  src = (ii*ngrid->c) + jj;
	  if (src == MYNODE) {
	    memcpy(rowbar,rowp,rowbar_sz);
	    rowp += ngrid->ny;
	  }
	  else
	    UMD_Recv(src, rowbar, rowbar_sz);
	  for (j=0 ; j<ngrid->ny ; j++) {
	    anum = (unsigned char)rowbar[j];
	    fprintf(savefile,"%c",anum);
	  }
	}

    fclose(savefile);
  }

  if (MYNODE != 0)
    on_one_thread {
      dst = 0;
      for (i=0 ; i<ngrid->nx ; i++) {
	UMD_Send(dst, rowp, rowbar_sz);
	rowp += ngrid->ny;
      }
  }

  all_Barrier(TH);
  free(rowbar);
  return;
}

void
all_SaveImage_orig_unroll(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j,          /* tile position */
    ii, jj;        /* logical grid position */

  FILE* savefile;
  unsigned char anum;
  int *rowbar, rowbar_sz;
  int src, dst;
  int *rowp;

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  rowbar_sz = ngrid->ny*sizeof(int);
  rowbar = (int *)malloc(rowbar_sz);
  assert_malloc(rowbar);

  all_Barrier(TH);

  rowp = X;

  on_one {

    savefile = fopen(im_info->fname,"w");
    if (savefile == NULL) errprnt("output file cannot be written.");
    rewind(savefile);

    /* HEADER */
    fprintf(savefile,"%s\n",im_info->format);
    fprintf(savefile,"%s\n",im_info->comment);
    fprintf(savefile,"%d %d\n",im_info->x, im_info->y);
    fprintf(savefile,"%d\n",im_info->k);

    for (ii = 0 ; ii<ngrid->r ; ii++)
      for (i=0 ; i<ngrid->nx ; i++) 
	for (jj = 0 ; jj<ngrid->c ; jj++) {
	  src = (ii*ngrid->c) + jj;
	  if (src == MYNODE) {
	    memcpy(rowbar,rowp,rowbar_sz);
	    rowp += ngrid->ny;
	  }
	  else
	    UMD_Recv(src, rowbar, rowbar_sz);
	  for (j=0 ; j<ngrid->ny ; j+=8) {
	    fprintf(savefile,"%c%c%c%c%c%c%c%c",
		    (unsigned char)rowbar[j  ],
		    (unsigned char)rowbar[j+1],
		    (unsigned char)rowbar[j+2],
		    (unsigned char)rowbar[j+3],
		    (unsigned char)rowbar[j+4],
		    (unsigned char)rowbar[j+5],
		    (unsigned char)rowbar[j+6],
		    (unsigned char)rowbar[j+7]);
	  }
	}

    fclose(savefile);
  }

  if (MYNODE != 0)
    on_one_thread {
      dst = 0;
      for (i=0 ; i<ngrid->nx ; i++) {
	UMD_Send(dst, rowp, rowbar_sz);
	rowp += ngrid->ny;
      }
  }

  all_Barrier(TH);
  free(rowbar);
  return;
}

void
all_SaveImage_one(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j;          /* tile position */

  FILE* savefile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset;

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  all_Barrier(TH);

  rowp = X;

  savefile = fopen(im_info->fname,"w");
  if (savefile == NULL) errprnt("output file cannot be written.");
  rewind(savefile);

  on_one {
    /* Write the HEADER */
    fprintf(savefile,"%s\n",im_info->format);
    fprintf(savefile,"%s\n",im_info->comment);
    fprintf(savefile,"%d %d\n",im_info->x, im_info->y);
    fprintf(savefile,"%d\n",im_info->k);
    pre_offset = (int)ftell(savefile);
  }

  pre_offset = all_Bcast_i(pre_offset, TH);

  on_one_thread {

    rowp = X;

    /* base_offest: number of pixels before my first */
    base_offset = pre_offset + ngrid->myCol * ngrid->ny +
      (ngrid->myRow * ngrid->nx * ngrid->y);
  
    for (i=0 ; i<ngrid->nx ; i++) { /* my row number */
      fseek(savefile, base_offset + (i*ngrid->y), SEEK_SET);
      for (j=0 ; j<ngrid->ny ; j++) {
	anum = (unsigned char)*rowp;
	rowp++;
	fprintf(savefile,"%c",anum);
      }
    }
  }

  fflush(savefile);
  all_Barrier(TH);
  
  fclose(savefile);

  all_Barrier(TH);
  return;
}

void
all_SaveImage_all_lockf(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j, rc, n, ls;          /* tile position */

  int savefile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset;
  unsigned char *rbuf;
  char wbuf[MAXLEN];

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  all_Barrier(TH);

  savefile = creat(im_info->fname,0666);
  if (savefile == NULL) errprnt("output file cannot be written.");
  lseek (savefile, 0L, SEEK_SET);

  on_one {
    /* Write the HEADER */
    pre_offset = 0;
    write(savefile, im_info->format, strlen(im_info->format));
    pre_offset += strlen(im_info->format);
    write(savefile,"\n",1);
    pre_offset ++;
    write(savefile, im_info->comment, strlen(im_info->comment));
    pre_offset += strlen(im_info->comment);
    write(savefile,"\n",1);
    pre_offset ++;
    sprintf(wbuf,"%d %d\n%d\n",im_info->x, im_info->y, im_info->k);
    write(savefile, wbuf, strlen(wbuf));
    pre_offset += strlen(wbuf);
  }

  pre_offset = all_Bcast_i(pre_offset, TH);

  /* base_offest: number of pixels before my first */
  base_offset = pre_offset + ngrid->myCol * ngrid->ny +
    (ngrid->myRow * ngrid->nx * ngrid->y);
  
  ls   = ngrid->ny * sizeof(unsigned char);
  rbuf = (unsigned char *)malloc(ls);
  assert_malloc(rbuf);
  node_pardo(i, 0, ngrid->nx, 1) { /* my row number */
    rowp = X + (i*ngrid->ny);
    lseek(savefile, base_offset + (i*ngrid->y), SEEK_SET);
    for (j=0 ; j<ngrid->ny ; j++) {
      rbuf[j] = (unsigned char)*rowp;
      rowp++;
    }
#if 1
    rc = lockf(savefile, F_LOCK, ls);
    if (rc != 0)
      perror("lockf");
#endif
    rc = write(savefile, rbuf, ls);
    if (rc != ngrid->ny)
      fprintf(stderr,"(%3d %3d)ERROR fwrite: rc: %d  ny: %d\n",
	      MYNODE, MYTHREAD, rc, ngrid->ny);
#if 0
    rc = lockf(savefile, F_ULOCK, -ls);
    if (rc != 0)
      perror("unlockf");
#endif
  }
  free(rbuf);

  all_Barrier(TH);

  close(savefile);

  all_Barrier(TH);
  return;
}

void
all_SaveImage_one_lockf(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j, rc, n, ls;          /* tile position */

  int savefile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset;
  unsigned char *rbuf;
  char wbuf[MAXLEN];

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  all_Barrier(TH);

  on_one_thread {
    savefile = creat(im_info->fname,0666);
    if (savefile == NULL) errprnt("output file cannot be written.");
    lseek (savefile, 0L, SEEK_SET);
  }

  on_one {
    /* Write the HEADER */
    pre_offset = 0;
    write(savefile, im_info->format, strlen(im_info->format));
    pre_offset += strlen(im_info->format);
    write(savefile,"\n",1);
    pre_offset ++;
    write(savefile, im_info->comment, strlen(im_info->comment));
    pre_offset += strlen(im_info->comment);
    write(savefile,"\n",1);
    pre_offset ++;
    sprintf(wbuf,"%d %d\n%d\n",im_info->x, im_info->y, im_info->k);
    write(savefile, wbuf, strlen(wbuf));
    pre_offset += strlen(wbuf);
  }

  pre_offset = all_Bcast_i(pre_offset, TH);
  
  on_one_thread {

    rowp = X;

    /* base_offest: number of pixels before my first */
    base_offset = pre_offset + ngrid->myCol * ngrid->ny +
      (ngrid->myRow * ngrid->nx * ngrid->y);

    ls   = ngrid->ny * sizeof(unsigned char);
    rbuf = (unsigned char *)malloc(ls);
    assert_malloc(rbuf);
    for (i=0 ; i<ngrid->nx ; i++) { /* my row number */
      lseek(savefile, base_offset + (i*ngrid->y), SEEK_SET);
      for (j=0 ; j<ngrid->ny ; j++) {
	rbuf[j] = (unsigned char)*rowp;
	rowp++;
      }
#if 1
      rc = lockf(savefile, F_LOCK, ls);
      if (rc != 0)
	perror("lockf");
#endif
      rc = write(savefile, rbuf, ls);
      if (rc != ngrid->ny)
	fprintf(stderr,"(%3d %3d)ERROR fwrite: rc: %d  ny: %d\n",
		MYNODE, MYTHREAD, rc, ngrid->ny);
#if 0
      rc = lockf(savefile, F_ULOCK, -ls);
      if (rc != 0)
	perror("unlockf");
#endif
    }
    free(rbuf);
    UMD_Barrier();
    close(savefile);
  }

  all_Barrier(TH);
  return;
}

void
all_SaveImage_flock(image_desc_t im_info, grid_t ngrid, image_t X, THREADED) {
  register int
    i, j;          /* tile position */

  FILE* savefile;
  unsigned char anum;
  int *rowp;
  int base_offset, pre_offset, sloc;
  struct flock alock;

#if DEBUG
  on_one {
    fprintf(outfile,"(SaveImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif
  
  all_Barrier(TH);

  on_one_thread {
    savefile = fopen(im_info->fname,"w");
    if (savefile == NULL) errprnt("output file cannot be written.");
#if 0
    setbuf(savefile, NULL);
#endif
    rewind(savefile);
  }

  alock.l_whence = SEEK_SET;
  alock.l_pid    = getpid();

  on_one {
    /* Write the HEADER */
#if 1
    fseek(savefile, 16 + (im_info->x * im_info->y), SEEK_SET);
    fprintf(savefile,"%c",X[0]);
    fseek(savefile, 0, SEEK_SET);
#endif
    alock.l_start  = 0;
    alock.l_len    = 0;
    alock.l_type   = F_WRLCK;
    fcntl(savefile, F_SETLK, &alock);
    fprintf(savefile,"%s\n",im_info->format);
    fprintf(savefile,"%s\n",im_info->comment);
    fprintf(savefile,"%d %d\n",im_info->x, im_info->y);
    fprintf(savefile,"%d\n",im_info->k);
    pre_offset = (int)ftell(savefile);
    alock.l_type   = F_UNLCK;
    fcntl(savefile, F_SETLK, &alock);
  }


  pre_offset = all_Bcast_i(pre_offset, TH);
  
  on_one_thread {

    rowp = X;

    /* base_offest: number of pixels before my first */
    base_offset = pre_offset + ngrid->myCol * ngrid->ny +
      (ngrid->myRow * ngrid->nx * ngrid->y);
  
    alock.l_len    = ngrid->y;
    for (i=0 ; i<ngrid->nx ; i++) { /* my row number */
      sloc = base_offset + (i*ngrid->y);
      fseek(savefile, sloc, SEEK_SET);
      alock.l_type   = F_WRLCK;
      alock.l_start  = sloc;
      fcntl(savefile, F_SETLK, &alock);
      for (j=0 ; j<ngrid->ny ; j++) {
	anum = (unsigned char)*rowp;
	rowp++;
	fprintf(savefile,"%c",anum);
      }
      alock.l_type  = F_UNLCK;
      fcntl(savefile, F_SETLK, &alock);
    }

    fflush(savefile);
    sleep(20);
    fclose(savefile);
  }

  all_Barrier(TH);
  return;
}



/*************************************************************/
void all_Get_pType(grid_t ngrid, image_t pType, THREADED)
/*************************************************************/
/* Fill pType with [1,9] for each pixel:

   1 | 2 | 3
   -----------
   4 | 5 | 6
   -----------
   7 | 8 | 9
*/
{
    register int
	i,
	j;

    int q, r;

#if DEBUG
    on_one {
      fprintf(outfile,"(Get_pType)\n");
      fflush(outfile);
    }
#endif
  
#if 0
    all_init_timer();
    all_start_timer();
#endif

    q = ngrid->nx;
    r = ngrid->ny;

    node_pardo(i,1,q-1,1)
      for (j=1 ; j<r-1 ; j++)
	pType[i*r + j] = 5;
    node_pardo(j,1,r-1,1) {
      pType[          j] = 2;
      pType[(q-1)*r + j] = 8;
    }
    node_pardo(i,1,q-1,1) {
      pType[i*r        ] = 4;
      pType[i*r + (r-1)] = 6;
    }
    on_one_thread {
      pType[              0] = 1;
      pType[            r-1] = 3;
      pType[(q-1)*r        ] = 7;
      pType[(q-1)*r + (r-1)] = 9;
    }
    
#if 0
    all_stop_timer("Get pType");
    all_print_timer(outfile);
#endif

    node_Barrier();
    return;
}


void
allocate_ghost4(grid_t ngrid, ghost4_t *gh4) {
#if DEBUG
  fprintf(outfile,"(allocate_ghost4)\n");
  fflush(outfile);
#endif
  
  *gh4 = (ghost4_t)malloc(sizeof(struct ghost4));
  assert_malloc(*gh4);

  (*gh4)->N = (image_t)malloc(ngrid->ny * sizeof(int));
  assert_malloc((*gh4)->N);
  (*gh4)->S = (image_t)malloc(ngrid->ny * sizeof(int));
  assert_malloc((*gh4)->S);
  (*gh4)->E = (image_t)malloc(ngrid->nx * sizeof(int));
  assert_malloc((*gh4)->E);
  (*gh4)->W = (image_t)malloc(ngrid->nx * sizeof(int));
  assert_malloc((*gh4)->W);
  return;
}

void
all_AllocateGhost4(grid_t ngrid, ghost4_t *gh4, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(AllocateGhost4)\n");
    fflush(outfile);
  }
#endif

  on_one_thread 
    allocate_ghost4(ngrid, gh4);
  *gh4 = (ghost4_t)node_Bcast_ip((int*)*gh4, TH);
  return;
}

void
free_ghost4(ghost4_t gh4) {
#if DEBUG
  fprintf(outfile,"(free_ghost4)\n");
  fflush(outfile);
#endif

  free(gh4->N);
  free(gh4->S);
  free(gh4->E);
  free(gh4->W);
  free(gh4);
  return;
}

void
all_FreeGhost4(ghost4_t gh4, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(FreeGhost4)\n");
    fflush(outfile);
  }
#endif
  
  on_one_thread
    free_ghost4(gh4);
  return;
}

void
allocate_ghost8(grid_t ngrid, ghost8_t *gh8) {
#if DEBUG
  fprintf(outfile,"(allocate_ghost8)\n");
  fflush(outfile);
#endif
  
  *gh8 = (ghost8_t)malloc(sizeof(struct ghost8));
  assert_malloc(*gh8);

  allocate_ghost4(ngrid, &((*gh8)->g4));
  
  (*gh8)->NW = (image_t)malloc(sizeof(int));
  assert_malloc((*gh8)->NW);
  (*gh8)->NE = (image_t)malloc(sizeof(int));
  assert_malloc((*gh8)->NE);
  (*gh8)->SW = (image_t)malloc(sizeof(int));
  assert_malloc((*gh8)->SW);
  (*gh8)->SE = (image_t)malloc(sizeof(int));
  assert_malloc((*gh8)->SE);
  return;
}

void
all_AllocateGhost8(grid_t ngrid, ghost8_t *gh8, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(AllocateGhost8)\n");
    fflush(outfile);
  }
#endif

  all_Barrier(TH);
  on_one_thread 
    allocate_ghost8(ngrid, gh8);
  node_Barrier();
  *gh8 = (ghost8_t)node_Bcast_ip((int*)*gh8, TH);
  all_Barrier(TH);
  return;
}

void
free_ghost8(ghost8_t gh8) {
#if DEBUG
  fprintf(outfile,"(free_ghost8)\n");
  fflush(outfile);
#endif
  
  free(gh8->SE);
  free(gh8->SW);
  free(gh8->NE);
  free(gh8->NE);
  free_ghost4(gh8->g4);
  free(gh8);
  return;
}

void
all_FreeGhost8(ghost8_t gh8, THREADED) {
#if DEBUG
  on_one {
    fprintf(outfile,"(FreeGhost8)\n");
    fflush(outfile);
  }
#endif
  
  on_one_thread
    free_ghost8(gh8);
  return;
}


void
all_get_ghost4(grid_t ngrid, image_t X, ghost4_t gh4,
	       int nil_elem, THREADED)
{
  int i, rowlen, collen;
  int src, dst;
  int *carr;

#if DEBUG
  on_one {
    fprintf(outfile,"(get_ghost4)\n");
    fflush(outfile);
  }
#endif

  all_Barrier(TH);
  
  rowlen = (ngrid->ny)*sizeof(int);
  collen = (ngrid->nx)*sizeof(int);

  carr = (int *)malloc(collen);
  assert_malloc(carr);
  
  if (task_do(0)) {
  /* ghostN */
    src = MYNODE-ngrid->c;
    dst = MYNODE+ngrid->c;
    if (ngrid->border_tp) {
      if (!ngrid->border_bm)
	UMD_Send(dst, X+(ngrid->ny)*(ngrid->nx-1), rowlen);
      for (i=0 ; i<ngrid->ny ; i++)
	gh4->N[i] = nil_elem;
    }
    else {
      if (ngrid->border_bm) {
	UMD_Recv(src, gh4->N, rowlen);
      }
      else {
	UMD_Sendrecv(dst, X+(ngrid->ny)*(ngrid->nx-1), rowlen,
		     src, gh4->N, rowlen);
      }
    }
  }

  if (task_do(1)) {
  /* ghostS */
    src = MYNODE+ngrid->c;
    dst = MYNODE-ngrid->c;
    if (ngrid->border_bm) {
      if (!ngrid->border_tp)
	UMD_Send(dst, X, rowlen);
      for (i=0 ; i<ngrid->ny ; i++)
	gh4->S[i] = nil_elem;
    }
    else {
      if (ngrid->border_tp) {
	UMD_Recv(src, gh4->S, rowlen);
      }
      else {
	UMD_Sendrecv(dst, X, rowlen,
		     src, gh4->S, rowlen);
      }
    }
  }

  if (task_do(2)) {
  /* ghostE */
    src = MYNODE+1;
    dst = MYNODE-1;
    for (i=0 ; i<ngrid->nx ; i++)
      carr[i] = X[i*ngrid->ny];
    if (ngrid->border_rt) {
      if (!ngrid->border_lt)
	UMD_Send(dst, carr, collen);
      for (i=0 ; i<ngrid->nx ; i++)
	gh4->E[i] = nil_elem;
    }
    else {
      if (ngrid->border_lt) {
	UMD_Recv(src, gh4->E, collen);
      }
      else {
	UMD_Sendrecv(dst, carr, collen,
		     src, gh4->E, collen);
      }
    }
  }

  if (task_do(3)) {
  /* ghostW */
    src = MYNODE-1;
    dst = MYNODE+1;
    for (i=0 ; i<ngrid->nx ; i++)
      carr[i] = X[(i+1)*ngrid->ny - 1];
    if (ngrid->border_lt) {
      if (!ngrid->border_rt)
	UMD_Send(dst, carr, collen);
      for (i=0 ; i<ngrid->nx ; i++)
	gh4->W[i] = nil_elem;
    }
    else {
      if (ngrid->border_rt) {
	UMD_Recv(src, gh4->W, collen);
      }
      else {
	UMD_Sendrecv(dst, carr, collen,
		     src, gh4->W, collen);
      }
    }
  }

  free(carr);

  node_Barrier();
  return;
}

void
all_get_ghost8(grid_t ngrid, image_t X, ghost8_t gh8,
	       int nil_elem, THREADED)
{
  int myNW, myNE, mySW, mySE;
  int src, dst;
  
#if DEBUG
  on_one {
    fprintf(outfile,"(get_ghost8)\n");
    fflush(outfile);
  }
#endif
  
  all_get_ghost4(ngrid, X, gh8->g4, nil_elem, TH);

  if (task_do(0)) {
  /* ghostNW */
    mySE = X[ngrid->nx*ngrid->ny - 1];
    src = MYNODE-ngrid->c - 1;
    dst = MYNODE+ngrid->c + 1;
    if (ngrid->ul) {
      if (!ngrid->lr)
	UMD_Send(dst, &mySE, sizeof(int));
      *(gh8->NW) = nil_elem;
    }
    else {
      if (ngrid->lr) {
	UMD_Recv(src, gh8->NW, sizeof(int));
      }
      else {
	if (ngrid->ce) {
	  UMD_Sendrecv(dst, &mySE, sizeof(int),
		       src, gh8->NW, sizeof(int));
	}
	else {
	  *(gh8->NW) = nil_elem;
	}
      }
    }
  }

  if (task_do(1)) {
  /* ghostSE */
    myNW = X[0];
    src = MYNODE+ngrid->c + 1;
    dst = MYNODE-ngrid->c - 1;
    if (ngrid->lr) {
      if (!ngrid->ul)
	UMD_Send(dst, &myNW, sizeof(int));
      *(gh8->SE) = nil_elem;
    }
    else {
      if (ngrid->ul) {
	UMD_Recv(src, gh8->SE, sizeof(int));
      }
      else {
	if (ngrid->ce) {
	  UMD_Sendrecv(dst, &myNW, sizeof(int),
		       src, gh8->SE, sizeof(int));
	}
	else {
	  *(gh8->SE) = nil_elem;
	}
      }
    }
  }

  if (task_do(2)) {
  /* ghostNE */
    mySW = X[(ngrid->nx-1)*ngrid->ny];
    src = MYNODE-ngrid->c + 1;
    dst = MYNODE+ngrid->c - 1;
    if (ngrid->ur) {
      if (!ngrid->ll)
	UMD_Send(dst, &mySW, sizeof(int));
      *(gh8->NE) = nil_elem;
    }
    else {
      if (ngrid->ll) {
	UMD_Recv(src, gh8->NE, sizeof(int));
      }
      else {
	if (ngrid->ce) {
	  UMD_Sendrecv(dst, &mySW, sizeof(int),
		       src, gh8->NE, sizeof(int));
	}
	else {
	  *(gh8->NE) = nil_elem;
	}
      }
    }
  }

  if (task_do(3)) {
  /* ghostSW */
    myNE = X[ngrid->ny - 1];
    src = MYNODE+ngrid->c - 1;
    dst = MYNODE-ngrid->c + 1;
    if (ngrid->ll) {
      if (!ngrid->ur)
	UMD_Send(dst, &myNE, sizeof(int));
      *(gh8->SW) = nil_elem;
    }
    else {
      if (ngrid->ur) {
	UMD_Recv(src, gh8->SW, sizeof(int));
      }
      else {
	if (ngrid->ce) {
	  UMD_Sendrecv(dst, &myNE, sizeof(int),
		       src, gh8->SW, sizeof(int));
	}
	else {
	  *(gh8->SW) = nil_elem;
	}
      }
    }
  }

  node_Barrier();

  return;
}

/*************************************************************/
void
my_proc_nbr(grid_t ngrid, image_t tile, image_t pType,
	    ghost8_t gh8, int nil_elem,
	    int i, int j, int *mask, int *nbrs) 
/*************************************************************/
{
  int q, r, idx;

  q = ngrid->nx;
  r = ngrid->ny;
  idx = i*r+j;
  
/* N  */
  switch (pType[idx]) {
  case 1:
  case 2:
  case 3:
    if (gh8->g4->N[j] == nil_elem)
      mask[0] = 0;
    else 
      nbrs[0] = gh8->g4->N[j];
    break;
  default: nbrs[0] = tile[(i-1)*r + j];
  }
/* S  */
  switch (pType[idx]) {
  case 7:
  case 8:
  case 9:
    if (gh8->g4->S[j] == nil_elem)
      mask[1] = 0;
    else
      nbrs[1] = gh8->g4->S[j];
    break;
  default: nbrs[1] = tile[(i+1)*r + j];
  }
/* W  */
  switch (pType[idx]) {
  case 1:
  case 4:
  case 7:
    if (gh8->g4->W[i] == nil_elem)
      mask[2] = 0;
    else
      nbrs[2] = gh8->g4->W[i];
    break;
  default: nbrs[2] = tile[i*r + j-1];
  }
/* E  */
  switch (pType[idx]) {
  case 3:
  case 6:
  case 9:
    if (gh8->g4->E[i] == nil_elem)
      mask[3] = 0;
    else
      nbrs[3] = gh8->g4->E[i];
    break;
  default: nbrs[3] = tile[i*r + j+1];
  }
/* NW */
  switch (pType[idx]) {
  case 1:
    if (*(gh8->NW) == nil_elem)
      mask[4] = 0;
    else
      nbrs[4] = *(gh8->NW);
    break;
  case 2:
  case 3:
    if (gh8->g4->N[j-1] == nil_elem)
      mask[4] = 0;
    else
      nbrs[4] = gh8->g4->N[j-1];
    break;
  case 4:
  case 7:
    if (gh8->g4->W[i-1] == nil_elem)
      mask[4] = 0;
    else
      nbrs[4] = gh8->g4->W[i-1];
    break;
  default: nbrs[4] = tile[(i-1)*r + j-1];
  }
/* NE */
  switch (pType[idx]) {
  case 3:
    if (*(gh8->NE) == nil_elem)
      mask[5] = 0;
    else
      nbrs[5] = *(gh8->NE);
    break;
  case 1:
  case 2:
    if (gh8->g4->N[j+1] == nil_elem)
      mask[5] = 0;
    else
      nbrs[5] = gh8->g4->N[j+1];
    break;
  case 6:
  case 9:
    if (gh8->g4->E[i-1] == nil_elem)
      mask[5] = 0;
    else
      nbrs[5] = gh8->g4->E[i-1];
    break;
  default: nbrs[5] = tile[(i-1)*r + j+1];
  }
/* SW */
  switch (pType[idx]) {
  case 7:
    if (*(gh8->SW) == nil_elem)
      mask[6] = 0;
    else
      nbrs[6] = *(gh8->SW);
    break;
  case 8:
  case 9:
    if (gh8->g4->S[j-1] == nil_elem)
      mask[6] = 0;
    else
      nbrs[6] = gh8->g4->S[j-1];
    break;
  case 1:
  case 4:
    if (gh8->g4->W[i+1] == nil_elem)
      mask[6] = 0;
    else
      nbrs[6] = gh8->g4->W[i+1];
    break;
  default: nbrs[6] = tile[(i+1)*r + j-1];
  }
/* SE */
  switch (pType[idx]) {
  case 9:
    if (*(gh8->SE) == nil_elem)
      mask[7] = 0;
    else
      nbrs[7] = *(gh8->SE);
    break;
  case 7:
  case 8:
    if (gh8->g4->S[j+1] == nil_elem)
      mask[7] = 0;
    else
      nbrs[7] = gh8->g4->S[j+1];
    break;
  case 3:
  case 6:
    if (gh8->g4->E[i+1] == nil_elem)
      mask[7] = 0;
    else
      nbrs[7] = gh8->g4->E[i+1];
    break;
  default: nbrs[7] = tile[(i+1)*r + j+1];
  }
  return;
}


/*************************************************************/
void all_filter_SNF(grid_t ngrid, image_t X, image_t Y, 
		    image_t pType, int steps, int epsilon, THREADED)
/*************************************************************/
/* Perform SNF on interior points (do not modify 1-pixel image
   border).
*/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	d0, d1, d2,
	idx,
	sum;

    int q, r, v, w;

    ghost8_t gh8;

    int
      temp[8],
      mask[8],
      tot,
      fp,
      myFP,
      tempFP,
      step,
      pfstate;
    
    image_t state;
    image_t SNFin;
    
    double pfixed,
	limit,
	lastpf;
    
    char buf[MAXLEN];

#if DEBUG
    on_one {
      fprintf(outfile,"(all_filter_SNF)\n");
      fflush(outfile);
    }
#endif

    all_AllocateImage(ngrid, &state, TH);
    all_AllocateImage(ngrid, &SNFin, TH);
    all_AllocateGhost8(ngrid, &gh8, TH);
    
#if 0
    all_init_timer();
    all_start_timer();
#endif

    all_ImageCopy(ngrid, SNFin, X, TH);

    q = ngrid->nx;
    r = ngrid->ny;
    v = ngrid->r;
    w = ngrid->c;
    
    tot = ((ngrid->x)-2) * ((ngrid->y)-2);

    node_pardo(i, 0, q, 1)
      for (j=0 ; j<r ; j++) 
	state[i*r + j] = !border_pix(ngrid,i,j);
    node_Barrier();
      
/* constants */

      limit = 100.0;
      pfixed = 0.0;

    /* pfstate is the number of iterations that pfixed has not changed.
       After 3 iterations of no change, stop iterating. */

      pfstate = 0;

#if 0
      all_stop_timer("SNF Initialization");
#endif
    

    for (step=1 ; ((step<=steps)&&(pfixed<limit)&&(pfstate<3)) ; step++) {

/* CHECK */
      myFP = 0;

/* Get copy of input image */

      if (step > 1)
	all_ImageCopy(ngrid, SNFin, Y, TH);

#if 0
      all_start_timer();
#endif

/******************************************************/

/* Create Ghost Cells for SNF input image */

      all_get_ghost8(ngrid, SNFin, gh8, NIL, TH);

/********************************************************/

      node_pardo(i,0,q,1)
	for (j=0 ; j<r ; j++) {		    
	  idx = i*r + j;
	  d0 = SNFin[idx];

	  if ((state[idx]>0) && (state[idx]<4)) {

	    for (t=0 ; t<8 ; t++)
	      mask[t] = 1;
	    
	    my_proc_nbr(ngrid, SNFin, pType, gh8, NIL, i, j, mask, temp);
	    
	    sum = (d0 << 2) + 2;

	    d1 = temp[2];
	    d2 = temp[3];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[1];
	    d2 = temp[0];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[7];
	    d2 = temp[4];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[6];
	    d2 = temp[5];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    Y[idx] = sum>>2;
	    tempFP = (Y[idx] == d0);
	    myFP += tempFP;
	    if (tempFP) {
	      switch (state[idx]) {
	      case 1:
	      case 2:
	      case 3:
		state[idx]++;
		break;
	      default: fprintf(stderr,"ERROR; default case reached\n");
	      }
	    }
	    else {
	      switch (state[idx]) {
	      case 1:
		break;
	      case 2:
	      case 3:
		state[idx] = 1;
		break;
	      default: fprintf(stderr,"ERROR; default case reached\n");
	      }
	    };
	  }
	  else {
	    Y[idx] = d0;
	    if (state[idx]==4) myFP++;
	  }
	}
      
#if 0
      sprintf(buf,"SNF iteration %3d",step);
      all_stop_timer(buf);
#endif

      fp = all_Allreduce_i(myFP, SUM, TH);

      lastpf = pfixed;

      pfixed = 100.0 * fp / tot;

      if (lastpf==pfixed)
	pfstate++;
      else
	pfstate=0;

      if (RESULTS)
	on_one 
	  fprintf(outfile,"SNF3x3  epsilon= %3d   step= %4d   fixed= %8.4f\n",
		  epsilon,step,pfixed);

    }

    all_FreeGhost8(gh8, TH);
    all_FreeImage(SNFin, TH);
    all_FreeImage(state, TH);
#if 0
    all_print_timer(outfile);
#endif
    return;
}

/*************************************************************/
void all_filter_SNF_nocopy(grid_t ngrid, image_t X, image_t Y, 
			   image_t pType, int steps, int epsilon, THREADED)
/*************************************************************/
/* Perform SNF on interior points (do not modify 1-pixel image
   border).
*/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	d0, d1, d2,
	idx,
	sum;

    int q, r, v, w;

    ghost8_t gh8;

    int
      temp[8],
      mask[8],
      tot,
      fp,
      myFP,
      tempFP,
      step,
      pfstate;
    
    image_t state;
    image_t imA, imB, imT;
    
    double pfixed,
	limit,
	lastpf;
    
    char buf[MAXLEN];

#if DEBUG
    on_one {
      fprintf(outfile,"(all_filter_SNF_nocopy)\n");
      fflush(outfile);
    }
#endif

    all_AllocateImage(ngrid, &state, TH);
    all_AllocateGhost8(ngrid, &gh8, TH);
    
#if 0
    all_init_timer();
    all_start_timer();
#endif

    imA = X;
    imB = Y;

    q = ngrid->nx;
    r = ngrid->ny;
    v = ngrid->r;
    w = ngrid->c;
    
    tot = ((ngrid->x)-2) * ((ngrid->y)-2);

    node_pardo(i, 0, q, 1)
      for (j=0 ; j<r ; j++) 
	state[i*r + j] = !border_pix(ngrid,i,j);
    node_Barrier();
      
/* constants */

    limit = 100.0;
    pfixed = 0.0;

    /* pfstate is the number of iterations that pfixed has not changed.
       After 3 iterations of no change, stop iterating. */

    pfstate = 0;

#if 0
    all_stop_timer("SNF Initialization");
#endif
    

    for (step=1 ; ((step<=steps)&&(pfixed<limit)&&(pfstate<3)) ; step++) {

/* CHECK */
      myFP = 0;

#if 0
      all_start_timer();
#endif

/******************************************************/

/* Create Ghost Cells for SNF input image */

      all_get_ghost8(ngrid, imA, gh8, NIL, TH);

/********************************************************/

      node_pardo(i,0,q,1)
	for (j=0 ; j<r ; j++) {		    
	  idx = i*r + j;
	  d0 = imA[idx];

	  if ((state[idx]>0) && (state[idx]<4)) {

	    for (t=0 ; t<8 ; t++)
	      mask[t] = 1;
	    
	    my_proc_nbr(ngrid, imA, pType, gh8, NIL, i, j, mask, temp);
	    
	    sum = (d0 << 2) + 2;

	    d1 = temp[2];
	    d2 = temp[3];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[1];
	    d2 = temp[0];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[7];
	    d2 = temp[4];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    d1 = temp[6];
	    d2 = temp[5];
	    if ((d1-=d0) && (d2-=d0) && (d1+d2))
	      if (d1+d2 > 0)
		if (d1 < d2) {
		  if (d1 <= epsilon)
		    sum += d1;
		}
		else {
		  if (d2 <= epsilon)
		    sum += d2;
		}
	      else
		if (d1 < d2) {
		  if (d2 >= -epsilon)
		    sum += d2;
		}
		else {
		  if (d1 >= -epsilon)
		    sum += d1;
		}

	    imB[idx] = sum>>2;
	    tempFP = (imB[idx] == d0);
	    myFP += tempFP;
	    if (tempFP) {
	      switch (state[idx]) {
	      case 1:
	      case 2:
	      case 3:
		state[idx]++;
		break;
	      default: fprintf(stderr,"ERROR; default case reached\n");
	      }
	    }
	    else {
	      switch (state[idx]) {
	      case 1:
		break;
	      case 2:
	      case 3:
		state[idx] = 1;
		break;
	      default: fprintf(stderr,"ERROR; default case reached\n");
	      }
	    };
	  }
	  else {
	    imB[idx] = d0;
	    if (state[idx]==4) myFP++;
	  }
	}
      
#if 0
      sprintf(buf,"SNF iteration %3d",step);
      all_stop_timer(buf);
#endif

      fp = all_Allreduce_i(myFP, SUM, TH);

      lastpf = pfixed;

      pfixed = 100.0 * fp / tot;

      if (lastpf==pfixed)
	pfstate++;
      else
	pfstate=0;

      if (RESULTS)
	on_one 
	  fprintf(outfile,"SNF3x3  epsilon= %3d   step= %4d   fixed= %8.4f\n",
		  epsilon,step,pfixed);

      imT = imA;
      imA = imB;
      imB = imT;

    }

    if (imA != Y)
      all_ImageCopy(ngrid, Y, X, TH);
    
    all_FreeGhost8(gh8, TH);
    all_FreeImage(state, TH);
#if 0
    all_print_timer(outfile);
#endif
    return;
}

	       
/*************************************************************/
void all_filter_1NN(grid_t ngrid, image_desc_t im_info,
		    image_t X, image_t Y, image_t pType, THREADED)
/*************************************************************/
/* Perform 1-NN filter. */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
        k,
        idx,
	t,
	pix;

    image_t T;

    int	temp[8],
	mask[8],
	diff,
	num,
	step, fp, lastfp, myFP;

    ghost8_t gh8;

    char buf[MAXLEN];
    
#if DEBUG
    on_one {
      fprintf(outfile,"(all_filter_1NN)\n");
      fflush(outfile);
    }
#endif

    all_Barrier(TH);

    all_AllocateImage(ngrid, &T, TH);
    all_AllocateGhost8(ngrid, &gh8, TH);

#if 0
    all_init_timer();
#endif

    all_ImageCopy(ngrid, T, X, TH);

    lastfp = 0;
    fp = 1;

    k = im_info->k;
    
    for (step=1 ; (lastfp != fp) ; step++) {

	myFP = 0;

/* Get copy of input image */

	if (step > 1)
	  all_ImageCopy(ngrid, T, Y, TH);
	
#if 0
	all_start_timer();
#endif
	
/******************************************************/

/* Create Ghost Cells for input image */

	all_get_ghost8(ngrid, T, gh8, NIL, TH);

/********************************************************/

	node_pardo(i, 0, ngrid->nx, 1)
	  for (j=0 ; j<ngrid->ny ; j++) {
	    idx = i*ngrid->ny + j;
		    
	    for (t=0 ; t<8 ; t++)
	      mask[t] = 1;
	    
	    my_proc_nbr(ngrid, T, pType, gh8, NIL, i, j, mask, temp);

	    pix = k;

	    for (t=0 ; t<8 ; t++) 
	      if (mask[t]) {
		diff = abs(T[idx] - temp[t]);
		if (diff < pix) {
		  pix = diff;
		  num = t;
		}
	      }

	    Y[idx] = (temp[num] + T[idx]) / 2 ;
	    if (Y[idx] == T[idx])
	      myFP++;
	  }

#if 0
	sprintf(buf,"1NN iteration %3d",step);
	all_stop_timer(buf);
#endif

	lastfp = fp;

	fp = all_Allreduce_i(myFP, SUM, TH);

	if (RESULTS)
	  on_one 
	    fprintf(outfile,"1NN   step= %4d   fixed= %8.4f\n", step,
		    (double)fp /
		    (double)(ngrid->r*ngrid->c*ngrid->nx*ngrid->ny));
	all_Barrier(TH);
    }

#if 0
    all_print_timer(outfile);
#endif
    all_FreeGhost8(gh8, TH);
    all_FreeImage(T, TH);
    return;
}


/*************************************************************/
void all_filter_1NN_nocopy(grid_t ngrid, image_desc_t im_info,
			   image_t X, image_t Y, image_t pType, THREADED)
/*************************************************************/
/* Perform 1-NN filter. */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
        k,
        idx,
	t,
	pix;

    image_t imA, imB, imT;
    
    int	temp[8],
	mask[8],
	diff,
	num,
	step, fp, lastfp, myFP;

    ghost8_t gh8;

    char buf[MAXLEN];
    
#if DEBUG
    on_one {
      fprintf(outfile,"(all_filter_1NN_nocopy)\n");
      fflush(outfile);
    }
#endif

    all_Barrier(TH);

    all_AllocateGhost8(ngrid, &gh8, TH);

#if 0
    all_init_timer();
#endif

    lastfp = 0;
    fp = 1;

    k = im_info->k;
    imA = X;
    imB = Y;
    
    for (step=1 ; (lastfp != fp) ; step++) {

	myFP = 0;

/* Get copy of input image */

#if 0
	all_start_timer();
#endif
	
/******************************************************/

/* Create Ghost Cells for input image */

	all_get_ghost8(ngrid, imA, gh8, NIL, TH);

/********************************************************/

	node_pardo(i, 0, ngrid->nx, 1)
	  for (j=0 ; j<ngrid->ny ; j++) {
	    idx = i*ngrid->ny + j;
		    
	    for (t=0 ; t<8 ; t++)
	      mask[t] = 1;
	    
	    my_proc_nbr(ngrid, imA, pType, gh8, NIL, i, j, mask, temp);

	    pix = k;

	    for (t=0 ; t<8 ; t++) 
	      if (mask[t]) {
		diff = abs(imA[idx] - temp[t]);
		if (diff < pix) {
		  pix = diff;
		  num = t;
		}
	      }

	    imB[idx] = (temp[num] + imA[idx]) / 2 ;
	    if (imB[idx] == imA[idx])
	      myFP++;
	  }

#if 0
	sprintf(buf,"1NN iteration %3d",step);
	all_stop_timer(buf);
#endif

	lastfp = fp;

	fp = all_Allreduce_i(myFP, SUM, TH);

	if (RESULTS)
	  on_one 
	    fprintf(outfile,"1NN   step= %4d   fixed= %8.4f\n", step,
		    (double)fp /
		    (double)(ngrid->r*ngrid->c*ngrid->nx*ngrid->ny));
	all_Barrier(TH);
	imT = imA;
	imA = imB;
	imB = imT;
    }

    if (imA != Y)
      all_ImageCopy(ngrid, Y, X, TH);

#if 0
    all_print_timer(outfile);
#endif
    all_FreeGhost8(gh8, TH);

    return;
}


/*************************************************************/
void all_crop_border(grid_t ngrid, image_t X, image_t Y, int crop,
		     THREADED)
/*************************************************************/
{
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	row,
	col;

    int q, r, v, w, ii, jj;
    int qv, rw;

    
#if DEBUG
    on_one {
      fprintf(outfile,"(all_crop_border)\n");
      fflush(outfile);
    }
#endif

    all_Barrier(TH);

#if 0
    all_init_timer();
    all_start_timer();
#endif

    q = ngrid->nx;
    r = ngrid->ny;
    v = ngrid->r;
    w = ngrid->c;
    ii = ngrid->myRow;
    jj = ngrid->myCol;
    qv = ngrid->x - crop;
    rw = ngrid->y - crop;

    node_pardo(i, 0, q, 1)
      for (j=0 ; j<r ; j++) {
	row = ii*q + i;
	col = jj*r + j;
	if ( (row < crop) || (row >= qv) ||
	     (col < crop) || (col >= rw) )
	  Y[i*r + j] = 0;
	else
	  Y[i*r + j] = X[i*r + j];
      }
    
    all_Barrier(TH);

#if 0
    all_stop_timer("Crop Border");
    all_print_timer(outfile);
#endif
    return;
}

void alloc_image_desc(image_desc_t *im_desc) {
  *im_desc = (image_desc_t)malloc(sizeof(struct image_desc));
  assert_malloc(*im_desc);
  return;
}

void alloc_grid(grid_t *gd) {
  *gd = (grid_t)malloc(sizeof(struct grid));
  assert_malloc(*gd);
  return;
}

/*************************************************************/
int all_decide_artificial(grid_t ngrid, image_t X, image_t pType,
			  THREADED)
/*************************************************************/
{
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	pix,
	range,
	ct,
	t,
	q, r;

    ghost8_t gh8;
    
    int	temp[8],
	mask[8];
    
    all_AllocateGhost8(ngrid, &gh8, TH);

    q = ngrid->nx;
    r = ngrid->ny;
    
#if 0
    all_init_timer();
    all_start_timer();
#endif
    
/******************************************************/

/* Create Ghost Cells for input image */

    all_get_ghost8(ngrid, X, gh8, NIL, TH);

/********************************************************/

    ct = 0;
    
    node_pardo(i, 0, q, 1)
      for (j=0 ; j<r ; j++) {

	/* Test for image border pixels */
	if ( !border_pix(ngrid, i,j) ) {
		
	  my_proc_nbr(ngrid, X, pType, gh8, NIL, i, j, mask, temp);

	  pix = X[i*r + j];
	  range = 0;
	  
	  for (t=0 ; t<8 ; t++)
	    if (temp[t] != pix)
	      range = 1;

	  if (!range)
	    ct++;
	}
      }


    range = all_Allreduce_i(ct, SUM, TH);
    ct = (ngrid->x - 2)*(ngrid->y - 2);

    all_FreeGhost8(gh8, TH);

#if 0
    all_stop_timer("Decide Artificial");
    all_print_timer(outfile);
#endif
    
    if (((double)range / (double)ct) >= 0.50)
      return(1);

    return(0);
}


void all_EnhanceImage(grid_t ngrid, image_desc_t im_info,
		      image_t X, image_t Y, image_t pType, THREADED)
{

  int art;
  all_timer_init();
  all_timer_reset();
  
  art = all_decide_artificial(ngrid, X, pType, TH);
  all_timer_mark("decide");
  on_one {
    fprintf(outfile,"(decide) Image is %s\n",art?"articial":"real");
    fflush(outfile);
  }

  all_timer_start();
  all_filter_SNF_nocopy(ngrid, X, Y, pType, 4,        0,        TH);
  all_timer_mark("SNF phase 1");
  all_filter_SNF_nocopy(ngrid, Y, Y, pType, SNFLIMIT, 6/*eps*/, TH);
  all_timer_mark("SNF phase 2");
  all_filter_SNF_nocopy(ngrid, Y, Y, pType, SNFLIMIT, 0,        TH);
  all_timer_mark("SNF phase 3");
  all_filter_1NN_nocopy(ngrid, im_info, Y, Y, pType, TH); 
  all_timer_mark("1NN");
  all_crop_border(ngrid, Y, Y, 3, TH); 
  all_timer_mark("crop border");

  all_timer_report(outfile,"Enhancing Image");
  return;
}

void *hello_world(void *parg) {
  int x;
  x = *(int *)parg;
  {
    int i;
    double d;
    d = 0.0;
    for (i=0 ; i<100 ; i++)
      d = (d + (double)i/2.0);
  }
  fprintf(outfile,"Hello world! I'm task %3d\n",x);
  fflush(outfile);
  return;
}

void *SIMPLE_main_orig(THREADED)
{
  image_desc_t im_desc;
#if 1
  image_desc_t im_desc_out;
#endif
  image_queue_t im_input;
  image_queue_t im_output;
  grid_t ngrid;
  image_t X, Y, pType;
  btask_t t1;
  int     t1_arg;
  
  int im_count;

  char buf[100];
  char *im_dir = "/fs/alf01-a/loc/dbader/imagery";

  all_timer_init();
#if DEBUG
  fprintf(outfile,"PE%3d(%3d): SIMPLE_main()\n",MYNODE,MYTHREAD);
  fflush(outfile);
#endif

  all_Barrier(TH);

#if 1
  alloc_image_desc(&im_desc_out);
#endif
  alloc_image_queue(&im_input);
  alloc_image_queue(&im_output);
  alloc_grid(&ngrid);
  all_GetParams(im_input, TH);
  all_MapNodes(ngrid, TH);

#if 0
  t1_arg = ID;
  t1 = create_btask(hello_world, &t1_arg);
#endif
#if 1
  on_one_thread {
    t1_arg = ID;
    t1 = create_btask(hello_world, &t1_arg);
  }
#endif

  for (im_count = 0 ; im_count < im_input->num ; im_count ++) {
#if 1
    on_one_thread waitfor_btask(t1);
    node_Barrier();
#endif
    im_desc = &(im_input->desc[im_count]);
    on_one {
      fprintf(outfile,"Image %3d: %s\n",im_count,im_desc->fname);
      fflush(outfile);
    }
    all_timer_reset();
#if 0
    all_LoadImage(im_desc, ngrid, &X, TH);
    sprintf(buf,"load image %4d x %4d",im_desc->x, im_desc->y);
    all_timer_mark(buf);
#endif
#if 0
    all_LoadImage_orig(im_desc, ngrid, &X, TH);
    sprintf(buf,"load image orig %4d x %4d",im_desc->x, im_desc->y);
    all_timer_mark(buf);
#endif
#if 1
    all_LoadImage_unroll(im_desc, ngrid, &X, TH);
    sprintf(buf,"load image u %4d x %4d",im_desc->x, im_desc->y);
    all_timer_mark(buf);
#endif
    sprintf(buf,"fname: %s",im_desc->fname);
    all_AllocateImage(ngrid, &pType, TH);
    all_Get_pType(ngrid, pType, TH);
  
    all_AllocateImage(ngrid, &Y, TH);

    all_timer_start();
    all_EnhanceImage(ngrid, im_desc, X, Y, pType, TH);
    all_timer_mark("enhance image");
  
    im_desc_out->x       = im_desc->x;
    im_desc_out->y       = im_desc->y;
    im_desc_out->k       = im_desc->k;
    strcpy(im_desc_out->comment, im_desc->comment);
    strcpy(im_desc_out->format, im_desc->format);

#if 0
    sprintf(im_desc_out->fname,"%s/test.raw.all.lockf.%d",im_dir,im_count);
    all_SaveImage_all_lockf(im_desc_out, ngrid, Y, TH);
    all_timer_mark("save image (all lockf)");
#endif
    
#if 0
    sprintf(im_desc_out->fname,"%s/test.raw.one.lockf.%d",im_dir,im_count);
    all_SaveImage_one_lockf(im_desc_out, ngrid, Y, TH);
    all_timer_mark("save image (one lockf)");
#endif

#if 0
    sprintf(im_desc_out->fname,"%s/test.raw.orig.%d",im_dir,im_count);
    all_SaveImage_orig(im_desc_out, ngrid, Y, TH);
    all_timer_mark("save image (orig)");
#endif
    
#if 1
    sprintf(im_desc_out->fname,"%s/test.raw.unroll.%d",im_dir,im_count);
    all_SaveImage_orig_unroll(im_desc_out, ngrid, Y, TH);
    all_timer_mark("save image (orig unr)");
#endif
    
    all_timer_report(outfile,"(main)");

    all_Barrier(TH);
    all_FreeImage(Y, TH);
    all_FreeImage(pType, TH);
    all_FreeImage(X, TH);
#if 1
    if (im_count < im_input->num - 1) 
      on_one_thread {
	t1_arg = im_count*100 + ID;
	t1 = create_btask(hello_world, &t1_arg);
      }
#endif
  }
#if 0
  waitfor_btask(t1);
#endif

  /*******************************/
  /* End of program              */
  /*******************************/
  
  free(ngrid);
  free(im_output);
  free(im_input);
  SIMPLE_done(TH);
}

typedef struct {
  image_t *imPtr;
  image_desc_t desc;
  grid_t ngrid;
} bg_load_image_t;

void *bg_load_image(void *parg) {
  bg_load_image_t *x;
  grid_t ngrid;
  image_desc_t im_info;
  image_t *imPtr;
  image_t X;
  
  register int
    i, j;          /* tile position */

  FILE* infile;
  unsigned char anum[8];
  int *rowp;
  int base_offset, pre_offset;

  x = (bg_load_image_t *)parg;
  ngrid   = x->ngrid;
  im_info = x->desc;
  imPtr   = x->imPtr;
  
#if DEBUG
  on_one {
    fprintf(outfile,"(LoadImage) x: %4d  y: %4d  fname: %s\n",
	    im_info->x,im_info->y,im_info->fname);
    fflush(outfile);
  }
#endif

  /* OPEN INPUT IMAGE */
  infile = fopen(im_info->fname,"r");
  if (infile == NULL) {
    fprintf(stderr,"input file [%s] does not exist.\n",im_info->fname);
    exit(-1);
  }
  rewind(infile);

  /* GET IMAGE DESCRIPTION */
  fscanf(infile,"%s\n",im_info->format);
  fscanf(infile,"%s\n",im_info->comment);
  fscanf(infile,"%d",  &(im_info->x));
  fscanf(infile,"%d\n",&(im_info->y));
  fscanf(infile,"%d\n",&(im_info->k));
  pre_offset = (int)ftell(infile);

  all_MapImage(im_info, ngrid);
  X = (int *)malloc(ngrid->nx * ngrid->ny * sizeof(int));
  assert_malloc(X);
  *imPtr = X;

  /*  all_Barrier(TH); */

  rowp = X;

  /* base_offest: number of pixels before my first */
  base_offset = pre_offset + ngrid->myCol * ngrid->ny +
    (ngrid->myRow * ngrid->nx * ngrid->y);
  
  for (i=0 ; i<ngrid->nx ; i++) {
    fseek(infile, base_offset + (i*ngrid->y), SEEK_SET);
    for (j=0 ; j<ngrid->ny ; j+=8) {
      fscanf(infile,"%c%c%c%c%c%c%c%c",
	     anum,   anum+1, anum+2, anum+3,
	     anum+4, anum+5, anum+6, anum+7);
      *rowp++ = (int)*anum;
      *rowp++ = (int)*(anum+1);
      *rowp++ = (int)*(anum+2);
      *rowp++ = (int)*(anum+3);
      *rowp++ = (int)*(anum+4);
      *rowp++ = (int)*(anum+5);
      *rowp++ = (int)*(anum+6);
      *rowp++ = (int)*(anum+7);
    }
  }

  fclose(infile);

  return;
}

void *SIMPLE_main(THREADED)
{
  image_desc_t im_desc;
#if 1
  image_desc_t im_desc_out;
#endif
  image_queue_t im_input;
  image_queue_t im_output;
  grid_t ngrid;
  image_t Y, pType;
  image_t imHead, imTail;
  btask_t t1;
  bg_load_image_t t1_arg;
  
  int im_count;

  char buf[100];
  char *im_dir = "/fs/alf01-a/loc/dbader/imagery";

  all_timer_init();
#if DEBUG
  fprintf(outfile,"PE%3d(%3d): SIMPLE_main()\n",MYNODE,MYTHREAD);
  fflush(outfile);
#endif

  all_Barrier(TH);

#if 1
  alloc_image_desc(&im_desc_out);
#endif
  alloc_image_queue(&im_input);
  alloc_image_queue(&im_output);
  alloc_grid(&ngrid);
  all_GetParams(im_input, TH);
  all_MapNodes(ngrid, TH);

#if 1
  all_timer_reset();
  im_desc = &(im_input->desc[0]);
  all_LoadImage_unroll(im_desc, ngrid, &imHead, TH);
  sprintf(buf,"load image u %4d x %4d",im_desc->x, im_desc->y);
  all_timer_mark(buf);
#endif

  for (im_count = 0 ; im_count < im_input->num ; im_count ++) {

    all_timer_reset();

    if (im_count < im_input->num - 1) {
      t1_arg.imPtr = &imTail;
      t1_arg.ngrid = ngrid;
      t1_arg.desc  = &(im_input->desc[im_count+1]);
      t1 = create_btask(bg_load_image, &t1_arg);
    }

    im_desc = &(im_input->desc[im_count]);
    on_one {
      fprintf(outfile,"Image %3d: %s\n",im_count,im_desc->fname);
      fflush(outfile);
    }
    all_timer_start(); 

    all_AllocateImage(ngrid, &pType, TH);
    all_Get_pType(ngrid, pType, TH);
  
    all_AllocateImage(ngrid, &Y, TH);
    all_timer_mark("allocate aux images");

    all_EnhanceImage(ngrid, im_desc, imHead, Y, pType, TH);
    all_timer_mark("enhance image");
  
    im_desc_out->x       = im_desc->x;
    im_desc_out->y       = im_desc->y;
    im_desc_out->k       = im_desc->k;
    strcpy(im_desc_out->comment, im_desc->comment);
    strcpy(im_desc_out->format, im_desc->format);

#if 1
    sprintf(im_desc_out->fname,"%s/test.raw.unroll.%d",im_dir,im_count);
    all_SaveImage_orig_unroll(im_desc_out, ngrid, Y, TH);
    all_timer_mark("save image (orig unr)");
#endif
    
    all_timer_report(outfile,"(loop)");

    all_Barrier(TH);
    all_FreeImage(Y, TH);
    all_FreeImage(pType, TH);
    all_FreeImage(imHead, TH);

    if (im_count < im_input->num - 1) {
      waitfor_btask(t1);
      all_Barrier(TH);
      imHead = imTail;
    }
  }

  all_timer_report(outfile,"(main)");
  
  /*******************************/
  /* End of program              */
  /*******************************/
  
  free(ngrid);
  free(im_output);
  free(im_input);
  SIMPLE_done(TH);
}

