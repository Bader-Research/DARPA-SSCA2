#include "lock.h"
#include "simple.h"

LOCK_T lock_init(THREADED)
{
  LOCK_T lock;

  lock=node_malloc(sizeof(pthread_mutex_t),TH);
  on_one_thread{
	  pthread_mutex_init(lock,NULL);
  }

  lock=(LOCK_T) node_Bcast_cp((char*) lock,TH);
  return (lock);
}

void lock_init_array(LOCK_T * lock_A, int n,THREADED)
{
  int i;

  pardo(i,0,n,1) {
    lock_A[i]=malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(lock_A[i],NULL);
  }
}

void lock_it(LOCK_T lock)
{
  pthread_mutex_lock(lock);
}

void unlock_it(LOCK_T lock)
{
  pthread_mutex_unlock(lock);
}


void lock_destroy(LOCK_T lock,THREADED)
{
  on_one_thread pthread_mutex_destroy(lock);
  node_Barrier();
  node_free(lock, TH);
}

void lock_destroy_array(LOCK_T * lock_A,int n,THREADED)
{
  int i;
  pardo(i, 0, n, 1) {
    pthread_mutex_destroy(lock_A[i]); 
    free(lock_A[i]);
  }
  node_Barrier();
  node_free(lock_A, TH);
}

