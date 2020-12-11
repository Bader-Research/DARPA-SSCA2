#ifndef _ALG_SPLASH_H
#define _ALG_SPLASH_H

#include <stdio.h>
#include <math.h>

/*************************************************************************/
/*                                                                       */
/*   MACROS for DEC Alpha                                                */
/*                                                                       */
/*************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>
#include "simple.h"

#define MAIN_ENV
#define MAIN_INITENV(a,b)
#define MAIN_END           

#define EXTERN_ENV

static pthread_t *_th;
static int _th_cnt;

#define G_MALLOC(n)        malloc(n)
     
#define CREATEINIT(n)    { _th = (pthread_t *)malloc((n)*sizeof(pthread_t)); \
                           _th_cnt = 0; }
                            
#define CREATE(f)          pthread_create(&_th[_th_cnt++], NULL, f, NULL);
#define WAIT_FOR_END(n)    { int i; for (i=0 ; i<(n) ; i++) \
                             pthread_join(_th[i], NULL); }

#define PAUSEDEC(a)        int (a);
#define PAUSEINIT(a)       (a)=0;
#define SETPAUSE(a)        (a)=1;
#define WAITPAUSE(a)       while (!(a)) sched_yield();
#define CLEARPAUSE(a)      (a)=0;

#define LOCKDEC(a)         pthread_mutex_t (a);
#define LOCKINIT(a)        pthread_mutex_init(&(a), NULL);
#define LOCK(a)            pthread_mutex_lock(&(a));
#define UNLOCK(a)          pthread_mutex_unlock(&(a));

#define ALOCKDEC(a,p)      
#define ALOCKINIT(a,p)

#define BARDEC(a)          smp_barrier_t (a);
#define BARINIT(a, n)      (a) = smp_barrier_init(n);
#define BARRIER(a, p)      smp_barrier_wait(a);

#define TIMEA              0

#if TIMEA
#define TIMEDEC(a)         unsigned int (a);
#define CLOCK(a)           (a) = time(0);
#else

#define TIMEDEC(a)         double (a);
#define CLOCK(a)           (a) = get_seconds();

#endif

#endif
