#include "simple.h"
#include "alg_random.h"

#define SEED 985456376

#define DEBUG 1

void *SIMPLE_main(THREADED)
{
  int i;
#if DEBUG
  fprintf(outfile,"PE%3d(%3d): SIMPLE_main()\n",MYNODE,MYTHREAD);
  fflush(outfile);
#endif

  all_Barrier(TH);
  
  rrand_init_th(SEED+MYTHREAD, TH);

  for (i=0 ; i<10 ; i++)
    fprintf(outfile,"PE%3d(%3d): sample %2d: %12d\n",MYNODE,MYTHREAD,i,rrand_th(TH));

  rrand_destroy_th(TH);

  /*******************************/
  /* End of program              */
  /*******************************/
  
  all_Barrier(TH);
  SIMPLE_done(TH);
}


