#include "simple-f.h"
      
      subroutine fsmain(fth)
      fthreaded

      parameter (n=16)
      common /fsmainc/ buf
      integer buf(n)
      integer i, i0, ierr
      integer step
      double precision d, d0
      integer mutex

c     TEST BUILT-IN VARIABLES
      print *, 'T',MYTHREAD,': I am THREAD ',MYTHREAD,' of ',THREADS,
&       ' on NODE ',MYNODE, 'of ',NODES 

c     TEST NODE BARRIER
      call fnbar(fth)

      buf(MYTHREAD) = MYTHREAD

      call fnbar(fth)
      
      on_one_thread
         do 20 i=1,THREADS
            print *, 'buf(',i,') = ', buf(i)
 20      continue 
      endif

c     TEST NODE PARDO
      step = 1
      call fnpardo(1, 10, step, fth)
      do 30 i=_M1,_M2, step
         print *, 'T',MYTHREAD,': fnpardo i= ',i
 30   continue

c     TEST NODE BROADCAST INTEGER and DOUBLE
      on_one_thread
         i = 47
         d = 8.3
      endif
      call fnbcasti(i, i0, fth)
      call fnbcastd(d, d0, fth)
      print *, 'T',MYTHREAD,': fnbcasti i0= ',i0,' fnbcastd d0= ',d0
      
c     TEST NODE REDUCE INTEGER
      i = MYTHREAD
      call fnreducei(i, OP_SUM, i0, fth)
      print *, 'T',MYTHREAD,': fnreducei i= ',i,' i0= ',i0

c     TEST NODE REDUCE DOUBLE
      d = MYTHREAD
      call fnreduced(d, OP_SUM, d0, fth)
      print *, 'T',MYTHREAD,': fnreduced d= ',d,' d0= ',d0
      
c     TEST NODE SCAN INTEGER
      i = MYTHREAD
      call fnscani(i, OP_SUM, i0, fth)
      print *, 'T',MYTHREAD,': fnscani i= ',i,' i0= ',i0

c     TEST NODE SCAN DOUBLE
      d = MYTHREAD
      call fnscand(d, OP_SUM, d0, fth)
      print *, 'T',MYTHREAD,': fnscand d= ',d,' d0= ',d0

      call fnbar(fth)
      
c     TEST MUTEX
      call fmutexinit(mutex, ierr, fth)
      on_one_thread
        buf(1) = 0
      endif
      call fnbar(fth)

      call fmutexlock(mutex, fth)
      buf(1) = buf(1) + 1
      call fmutexunlock(mutex, fth)

      call fnbar(fth)
      on_one_thread
        print *, 'MUTEX test, buf = ',buf(1)
      endif
      call fmutexdestroy(mutex, ierr, fth)
      
      
      return
      END
      
