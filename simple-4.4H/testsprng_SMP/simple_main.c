#include "simple.h"
#include "alg_random.h"

#define DEBUG 1

void *SIMPLE_main(THREADED)
{
  int i;

#if DEBUG
  fprintf(outfile,"T(%3d): SIMPLE_main()\n",MYTHREAD);
  fflush(outfile);
#endif

  node_Barrier();

  rrand_init_th(make_sprng_seed(), TH);

  for (i=0 ; i<10 ; i++)
    fprintf(outfile,"T%3d: sample %2d: %12d\n",MYTHREAD,i,rrand_th(TH));

  rrand_destroy_th(TH);

  /*******************************/
  /* End of program              */
  /*******************************/
  
  node_Barrier();
  SIMPLE_done(TH);
}

