#include <stdio.h>

#include "alg_sort.h"

#define LOOPS      10

void check_seq_sort(int *a, int n, char *s) {
  register int i;
  for (i=1 ; i<n ; i++)
    if (a[i] < a[i-1]) 
      fprintf(stderr,"PE%3d: ERROR: (%s) a[%3d] < a[%3d]  (%6d, %6d)\n",
	      MYNODE,s,i,i-1,a[i],a[i-1]);
}

int *seq_radixsort1(int *a, int *b, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	pass;

    int	nextpass,
	*count,
	*sbitArr;

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(a[i],pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) b[--count[sbitArr[i]]] = a[i]; 

	nextpass = pass+1;
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(b[i],nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) a[--count[sbitArr[i]]] = b[i];
    }

    free(sbitArr);
    free(count);

    return(a);
}

int *seq_radixsort2(int *a, int *b, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	M;

    int	*count,
	*sbitArr;

    m = 11;                 /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(a[i],0,m)]++;
    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[sbitArr[i]]] = a[i]; 

    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(b[i],m,m)]++;
    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) a[--count[sbitArr[i]]] = b[i];

    m = 10;
    M = (1<<m);
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(a[i],22,m)]++;
    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[sbitArr[i]]] = a[i];

    free(sbitArr);
    free(count);

    return(b);
}

int *seq_radixsort2_overlap(int *a, int *b, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	M;

    int	*count1,
	*sbitArr1,
	*count2,
	*sbitArr2;

    m = 11;                 /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    if ((count1 = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr1 = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

    if ((count2 = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr2 = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (j=0 ; j<M ; j++) {
	count1[j] = 0;
	count2[j] = 0;
    }
    for (i=0 ; i<n ; i++) count1[sbitArr1[i] = bits(a[i],0,m)]++;
    for (j=1 ; j<M ; j++) count1[j] += count1[j-1];
    for (i=n-1 ; i>=0 ; i--) {
	b[j = --count1[sbitArr1[i]]] = a[i];
	count2[sbitArr2[j] = bits(a[i],m,m)]++;
    }

    count1[0] = 0;
    for (j=1 ; j<M ; j++) {
	count1[j] = 0;
	count2[j] += count2[j-1];
    }
    for (i=n-1 ; i>=0 ; i--) {
	a[j = --count2[sbitArr2[i]]] = b[i];
	count1[sbitArr1[j] = bits(b[i],22,m)]++;
    }

    m = 10;
    M = (1<<m);
    for (j=1 ; j<M ; j++) count1[j] += count1[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count1[sbitArr1[i]]] = a[i];

    free(sbitArr2);
    free(count2);
    free(sbitArr1);
    free(count1);

    return(b);
}

int *seq_radixsort_16(int *a, int *b, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	M;

    int	*count,
	*sbitArr;

    m = 16;
    M = 1 << m;               /* The range of each pass */
    
    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(a[i],0,m)]++;
    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[sbitArr[i]]] = a[i]; 

    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[sbitArr[i] = bits(b[i],m,m)]++;
    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) a[--count[sbitArr[i]]] = b[i];

    free(sbitArr);
    free(count);

    return(a);
}


int *seq_radixsort_16_overlap(int *a, int *b, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	k,
	M;

    int	*count1,
	*count2,
	*sbitArr1,
	*sbitArr2;

    M = 1 << 16;             /* The range of each pass */
    
    if ((count1 = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr1 = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

    if ((count2 = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");
    if ((sbitArr2 = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (j=0 ; j<M ; j++) {
	count1[j] = 0;
	count2[j] = 0;
    }
    for (i=0 ; i<n ; i++) count1[sbitArr1[i] = bits(a[i],0,16)]++;
    for (j=1 ; j<M ; j++) count1[j] += count1[j-1];
    for (i=n-1 ; i>=0 ; i--) {
	b[j = --count1[sbitArr1[i]]] = a[i];
	count2[sbitArr2[j] = bits(a[i],16,16)]++;
    }
	    
    for (j=1 ; j<M ; j++) count2[j] += count2[j-1];
    for (i=n-1 ; i>=0 ; i--) a[--count2[sbitArr2[i]]] = b[i];

    free(sbitArr2);
    free(count2);
    free(sbitArr1);
    free(count1);

    return(a);
}

void rsort_input_s(int n, int *x, int *y,
		   int *(*fun)(), char* s) {
    double lsec, secs;
    int i, *xcpy, *ret;

    xcpy = (int *)malloc(n*sizeof(int));
    assert_malloc(xcpy);
    
    lsec = get_seconds();
    for (i=0 ; i<LOOPS ; i++) {
	memcpy(xcpy, x, n*sizeof(int));
    }
    lsec = get_seconds() - lsec;

#if DEBUG
    fprintf(outfile,"PE%3d: (%s) lsec overhead: %9.6f\n",MYNODE,s,lsec);
    fflush(outfile);
#endif
    
    secs = get_seconds();

    for (i=0 ; i<LOOPS ; i++) {
	memcpy(xcpy, x, n*sizeof(int));
	ret = fun(xcpy, y, n);
    }

    secs = get_seconds() - secs;
    secs = (secs/(double)LOOPS) - lsec;

    fprintf(outfile,"PE%3d: Time: %s [%8d]: %9.6f \n",MYNODE,s,n,secs);
    fflush(outfile);

    check_seq_sort(ret,n,s);

    free(xcpy);
}

