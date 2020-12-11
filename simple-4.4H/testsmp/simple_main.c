#include "simple.h"

#define DEBUG 1

void *SIMPLE_main(THREADED)
{
#if DEBUG
  fprintf(outfile,"T(%3d): SIMPLE_main()\n",MYTHREAD);
  fflush(outfile);
#endif

  node_Barrier();
  

  /*******************************/
  /* End of program              */
  /*******************************/
  
  node_Barrier();
  SIMPLE_done(TH);
}

