#ifndef _IM_H
#define _IM_H

#include "simple.h"

#define MAXLEN 80

typedef int* image_t;

typedef struct image_desc {
  int x;                /* rows     */
  int y;                /* columns  */
  int k;                /* k */
  char fname[MAXLEN];   /* filename */
  char comment[MAXLEN]; /* comment */
  char format[MAXLEN];  /* image format */
} *image_desc_t;

typedef struct image_queue {
  int num;                 /* number of images */
  image_desc_t desc;       /* image descriptors */
} *image_queue_t;

typedef struct grid {
  int r;          /* node grid rows    */
  int c;          /* node grid columns */
  int nx;         /* rows of image pixels per node */
  int ny;         /* columns of image pixels per node */
  int myRow;      /* my row */
  int myCol;      /* my col */
  int x;          /* image rows */
  int y;          /* image cols */
  int border_tp; /* is my node a top    border? */
  int border_bm; /* is my node a bottom border? */
  int border_lt; /* is my node a left   border? */
  int border_rt; /* is my node a right  border? */
  int ul;        /* upper left */
  int ur;        /* upper right */
  int ll;        /* lower left */
  int lr;        /* lower right */
  int ce;        /* center */
} *grid_t;

typedef struct ghost4 {
  image_t N;
  image_t S;
  image_t E;
  image_t W;
} *ghost4_t;

typedef struct ghost8 {
  ghost4_t g4;
  image_t NW;
  image_t NE;
  image_t SW;
  image_t SE;
} *ghost8_t;


#endif

