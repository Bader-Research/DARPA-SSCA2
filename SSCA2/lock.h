#ifndef _LOCK_H_
#define _LOCK_H_

#include <pthread.h>
#include "simple.h"

typedef pthread_mutex_t* LOCK_T;

LOCK_T lock_init(THREADED);
void lock_init_array(LOCK_T*, int, THREADED);
void lock_it(LOCK_T lock);
void unlock_it(LOCK_T lock);
void lock_destroy(LOCK_T lock,THREADED);
void lock_destroy_array(LOCK_T *, int, THREADED);

#endif


