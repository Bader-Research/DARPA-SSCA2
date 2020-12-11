#include "simple.h"

#define DEBUG 1

void *SIMPLE_main(THREADED)
{
#if DEBUG
  fprintf(outfile,"PE%3d(%3d): SIMPLE_main()\n",MYNODE,MYTHREAD);
  fflush(outfile);
#endif

  all_Barrier(TH);
  

  /*******************************/
  /* End of program              */
  /*******************************/
  
  all_Barrier(TH);
  SIMPLE_done(TH);
}

