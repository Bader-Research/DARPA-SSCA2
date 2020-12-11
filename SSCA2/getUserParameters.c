#include "getUserParameters.h"
#include "globals.h"
#include "simple.h"

void getUserParameters()
{
  /* Scalable Data Generator parameters - defaults */
  SCALE                 =  16;           /* Binary Scaling Heuristic */
  TOT_VERTICES          =  (1<<SCALE);   /* Total number of vertices in directed multigraph. */
  MAX_CLIQUE_SIZE       =  (1<<(SCALE/3));  /* Maximum allowed clique size in directed multigraph. */
  MAX_PARAL_EDGES       =   3;           /* Max num of parallel edges allowed between two vertices. */
  PERC_INT_WEIGHTS      =  0.6;          /* Percentage of integer (vs. char string) edge weights. */
  MAX_INT_WEIGHT        =  (1<<SCALE);   /* Max allowed integer value in any integer edge weight. */
  PROB_UNIDIRECTIONAL   =  0.1;
  PROB_INTERCL_EDGES    =  0.5;          /* Initial probability of a link between two cliques. */
  
  MAX_STRLEN            =  SCALE;

  SOUGHT_STRING         =  "";           /* Kernel 2: Character string sought: specify it here, or else it is */
                                         /* picked from a randomly selected entry in genScalData.c */

  
  SUBGR_EDGE_LENGTH     =    3;          /* Kernel 3: max. path length, measured by the no. of edges  */
                                         /* in the subgraph generated from the end Vertex of SI and SC lists */
                                         
  MAX_CLUSTER_SIZE      =  MAX_CLIQUE_SIZE;   /* Kernel 4: Clustering search box size. */
  
  /* Some implementation-specific vars, nothing to do with the specs */
  K3_DS                 =  2;            /* 0 - Array, 1 - Linked List, 2 - Dynamic Array */
    
}

