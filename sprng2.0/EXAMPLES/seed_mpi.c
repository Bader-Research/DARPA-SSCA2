/***************************************************************************/
/*           ____Demonstrates the use of make_seed with MPI____            */
/* 'make_sprng_seed' is called to produce a new seed each time the program */
/* is run. The same seed is produced on each process.                      */
/***************************************************************************/

#include <stdio.h>
#include <mpi.h>                /* MPI header file                         */


/* Uncomment the following line to get the interface with pointer checking */
/*#define CHECK_POINTERS                                                   */
 
#define USE_MPI                 /* SPRNG makes MPI calls                   */
#include "sprng.h"              /* SPRNG header file                       */


main(int argc, char *argv[])
{
  int streamnum, nstreams, seed, *stream, i, myid, nprocs;
  double rn;
  int gtype;  /*---    */


  /*************************** MPI calls ***********************************/

  MPI_Init(&argc, &argv);	/* Initialize MPI                          */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);	/* find process id                 */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* find number of processes      */

  /************************** Initialization *******************************/

  streamnum = myid;
  nstreams = nprocs;		/* one stream per processor                */
  seed = make_sprng_seed();	/* make new seed each time program is run  */
  
  /*--- node 0 is reading in a generator type */
  if(myid == 0)
  {
#include "gen_types_menu.h"
    printf("Type in a generator type (integers: 0,1,2,3,4,5):  ");
    scanf("%d", &gtype);
  }
  MPI_Bcast(&gtype,1,MPI_INT,0,MPI_COMM_WORLD );

  /* Seed should be the same on all processes                              */
  printf("Process %d: seed = %16d\n", myid, seed);
  
  stream = init_sprng(gtype,streamnum,nstreams,seed,SPRNG_DEFAULT);	/*initialize stream*/
  printf("\n\nProcess %d: Print information about stream:\n",myid);
  print_sprng(stream);

  /************************ print random numbers ***************************/

  for (i=0;i<3;i++)
  {
    rn = sprng(stream);		/* generate double precision random number */
    printf("process %d, random number %d: %f\n", myid, i+1, rn);
  }

  free_sprng(stream);           /* free memory used to store stream state  */

  MPI_Finalize();		/* Terminate MPI                           */
}
