#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "simple.h"
#include "alg_radix_smp.h"
#include "alg_random.h"
#include "sprng.h"
#include "defs.h"
#include "genScalData.h"
#include "globals.h"
#include "lock.h"


void genScalData(graphSDG* SDGdataPtr, THREADED)
{
  /* Temp. vars for clique generation */
  int* cliqueSizes;
  ULONGINT_T* firstVsInCliques;
  ULONGINT_T* lastVsInCliques;
  ULONGINT_T* startVertex;
  ULONGINT_T* endVertex;

  int estTotCliques, totCliques;
  
  /* Edgelist vars */
  ULONGINT_T estTotEdges, numEdgesPlacedInCliques, numEdgesPlacedOutside,
  numEdgesPlaced, numStrEdges, edgeNum;
 
  /* The vars associated with the graph tuple */
  ULONGINT_T* permV;

  /* For SPRNG random no. generation*/
  int* stream;
  
  /* Other temp vars */
  int i_clique, i_cliqueSize, i_paralEdge;
  ULONGINT_T i_firstVsInClique, i_edgePtr, *i_edgeStartCounter,*i_edgeEndCounter;

  int endVClique;
  ULONGINT_T randNum, randNumEdges;
  ULONGINT_T *startV, *endV;
  ULONGINT_T** tmpEdgeCounter;
  ULONGINT_T* tempIndex;
  
  ULONGINT_T numStrWtEdges;
  ULONGINT_T tempVertex1, tempVertex2;
  ULONGINT_T d, i, j, k, i0, i1, t, t1, t2, h, l, m;
  float p;
  float rn;
  double timed;
  LOCK_T permLock;
  
  FILE* outfp;

  /*--------------------------------------------------------------------------*/
  /* STEP 0: Create the permutations required to randomize the vertices       */
  /*--------------------------------------------------------------------------*/
       
  /* random number gen. */
  stream = init_sprng(SPRNG_LCG, MYTHREAD, THREADS, (MYTHREAD+1)*get_seconds(), SPRNG_DEFAULT);
  
  /*rrand_init_th(make_sprng_seed(), TH); */
  
  node_Barrier();
  
  permLock = (LOCK_T) lock_init(TH);
  

  permV = (ULONGINT_T *) node_malloc(TOT_VERTICES*sizeof(ULONGINT_T), TH);
  
  pardo(i, 0, TOT_VERTICES, 1) {
    /* Initialize the array */
    permV[i] = i;
    
  }

  node_Barrier();
  
  pardo(i, 0, TOT_VERTICES, 1) {
  
    t1 = isprng(stream);
    
    t = i + t1 % (TOT_VERTICES - i);
    
    if (t != i) {

      lock_it(permLock);  
      t2 = permV[t];
      permV[t]=permV[i];
      permV[i]=t2;  
      unlock_it(permLock);

    }

  }
 

  node_Barrier();
  
  lock_destroy(permLock, TH);
  
  node_Barrier();
  
 
  /*--------------------------------------------------------------------------*/
  /* STEP 1: Create Cliques.                                                  */
  /*--------------------------------------------------------------------------*/

  /* Estimate number of clique required and pad by 50% */
  estTotCliques = ceil(1.5 * TOT_VERTICES / ((1+MAX_CLIQUE_SIZE)/2) );

  /* Allocate mem for Clique array */
  cliqueSizes = (int *) node_malloc(estTotCliques*sizeof(int), TH);

  /* Generate random clique sizes. */
  pardo(i, 0, estTotCliques, 1) {
    cliqueSizes[i] = 1 + (int) (isprng(stream) % MAX_CLIQUE_SIZE);
    /* printf("%d ", cliqueSizes[i]); */
  }
  
  totCliques = 0;
  
  node_Barrier();
  /* Allocate memory for cliqueList */
  lastVsInCliques = (ULONGINT_T *) node_malloc(estTotCliques*sizeof(ULONGINT_T), TH);
  firstVsInCliques = (ULONGINT_T *) node_malloc(estTotCliques*sizeof(ULONGINT_T), TH);

  node_Barrier();

  /* Sum up vertices in each clique to determine the lastVsInCliques array */
  
  on_one_thread {
    lastVsInCliques[0] = cliqueSizes[0]-1;
    for (i=1; i<estTotCliques; i++) {
      lastVsInCliques[i] = cliqueSizes[i]+lastVsInCliques[i-1];
      if (lastVsInCliques[i] >= TOT_VERTICES-1)
        break;
    }
    totCliques = i+1;
    
    /* Fix the size of the last clique */
    cliqueSizes[totCliques-1] = TOT_VERTICES - lastVsInCliques[totCliques-2] - 1;
    lastVsInCliques[totCliques-1] = TOT_VERTICES-1;
  }
  
  node_Barrier();

  totCliques = node_Bcast_i(totCliques, TH);
      
  
    /* Compute start Vertices in cliques. */
  firstVsInCliques[0] = 0;
  pardo(i, 1, totCliques, 1) {
    firstVsInCliques[i] = lastVsInCliques[i-1] + 1;
  }  

  node_Barrier();

  /* Write the generated cliques to file for comparison with Kernel 4 */
  on_one_thread {
    outfp = fopen("cliques.txt", "w");

    fprintf(outfp, "No. of cliques - %d\n", totCliques);
    for (i=0; i<totCliques; i++) {
      fprintf(outfp, "Clq %d - ", i);
      for (j=firstVsInCliques[i]; j<=lastVsInCliques[i]; j++) {      
	      fprintf(outfp, "%d ", permV[j]);
      }
      	
      fprintf(outfp, "\n");
    }
    fclose(outfp);
    
  }
  
  node_Barrier();    

  /*--------------------------------------------------------------------------*/
  /* STEP 2: Create the edges within the cliques                              */
  /*--------------------------------------------------------------------------*/
  
  /* Estimate number of edges - using an empirical measure */
  if (SCALE >= 12)
    estTotEdges = ceil(((MAX_CLIQUE_SIZE-1)*TOT_VERTICES));
  else
    estTotEdges = ceil(1.2*(((MAX_CLIQUE_SIZE-1)*TOT_VERTICES)*((1 + MAX_PARAL_EDGES)/2)+TOT_VERTICES*2));

  node_Barrier();
  /*
  on_one_thread
    printf("Estimated total edges is %d \n", estTotEdges);
  */  
  /* Initialize edge counter */
  edgeNum = 0;
  i_edgePtr = 0;
  p = PROB_UNIDIRECTIONAL;

  /* Partial edgeLists */
  if (THREADS > 3) {
    startV = (ULONGINT_T *) malloc(1.5*(estTotEdges/THREADS)*sizeof(ULONGINT_T));
    endV = (ULONGINT_T *) malloc(1.5*(estTotEdges/THREADS)*sizeof(ULONGINT_T));
  } else  {
    startV = (ULONGINT_T *) malloc((estTotEdges/THREADS)*sizeof(ULONGINT_T));
    endV = (ULONGINT_T *) malloc((estTotEdges/THREADS)*sizeof(ULONGINT_T));
  }
    
  /* Tmp array to keep track of the no. of parallel edges in each direction  */
  tmpEdgeCounter = (ULONGINT_T **) malloc(MAX_CLIQUE_SIZE*sizeof(ULONGINT_T *));
  for (i=0; i<MAX_CLIQUE_SIZE; i++) {
    tmpEdgeCounter[i] = (ULONGINT_T *) malloc(MAX_CLIQUE_SIZE*sizeof(ULONGINT_T));
  }
  
  node_Barrier();
  
  /* Create edges in parallel */
  pardo (i_clique, 0, totCliques, 1) {
    /* Get current clique parameters */
    i_cliqueSize = cliqueSizes[i_clique];
    i_firstVsInClique = firstVsInCliques[i_clique];
    
    /* First create atleast one edge between two vetices in a clique */
    for (i = 0; i < i_cliqueSize; i++) {
      for (j = 0; j < i; j++) {
          if (sprng(stream) >= p) {

            startV[i_edgePtr] = i + i_firstVsInClique;
	          endV[i_edgePtr] = j + i_firstVsInClique;
	          i_edgePtr = i_edgePtr + 1;
	          tmpEdgeCounter[i][j] = 1;
	          
            startV[i_edgePtr] = j + i_firstVsInClique;
	          endV[i_edgePtr] = i + i_firstVsInClique;
	          i_edgePtr = i_edgePtr + 1;
	          tmpEdgeCounter[j][i] = 1;
	          
          } else {

            if (sprng(stream) >= 0.5) {

              startV[i_edgePtr] = i + i_firstVsInClique;
	            endV[i_edgePtr] = j + i_firstVsInClique;
	            i_edgePtr = i_edgePtr + 1;
	            tmpEdgeCounter[i][j] = 1;
              tmpEdgeCounter[j][i] = 0;
            } else {

              startV[i_edgePtr] = j + i_firstVsInClique;
	            endV[i_edgePtr] = i + i_firstVsInClique;
	            i_edgePtr = i_edgePtr + 1;
	            tmpEdgeCounter[j][i] = 1;
	            tmpEdgeCounter[i][j] = 0;
            }

         }

      }
    }
    
    if (i_cliqueSize != 1) {
      randNumEdges = (int) (isprng(stream) % (2*i_cliqueSize*MAX_PARAL_EDGES));
      for (i_paralEdge = 0; i_paralEdge < randNumEdges; i_paralEdge++) {
               
	      i = (ULONGINT_T) (isprng(stream) % i_cliqueSize);
	      j = (ULONGINT_T) (isprng(stream) % i_cliqueSize);
	
	      if ((i != j) && (tmpEdgeCounter[i][j] < MAX_PARAL_EDGES)) {
	        if (sprng(stream) >= p) {
            /* Copy to edge structure. */
            startV[i_edgePtr] = i + i_firstVsInClique;
            endV[i_edgePtr] = j + i_firstVsInClique;
            i_edgePtr = i_edgePtr + 1;
	          tmpEdgeCounter[i][j] = tmpEdgeCounter[i][j] + 1;
	        }
	      }
      } 
    }
  }
  
  node_Barrier(); 
   
  for (i=0; i<MAX_CLIQUE_SIZE; i++)
     free(tmpEdgeCounter[i]);
  
  free(tmpEdgeCounter);

  node_Barrier();    
  
  /* Merge partial edge lists */
  
  i_edgeStartCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  i_edgeEndCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
    
  on_thread(MYTHREAD) {
    i_edgeEndCounter[MYTHREAD] = i_edgePtr;     
    i_edgeStartCounter[MYTHREAD] = 0;
  }

  node_Barrier();
  
  on_one_thread {
    for (i=1; i<THREADS; i++) {
      i_edgeEndCounter[i] = i_edgeEndCounter[i-1] + i_edgeEndCounter[i];
      i_edgeStartCounter[i] = i_edgeEndCounter[i-1]; 
    }
  }
    
  edgeNum = node_Reduce_d(i_edgePtr, SUM, TH);
  
  node_Barrier();
  
    /* Initialize edge list arrays */
  if (SCALE < 10)
    startVertex = (ULONGINT_T *) node_malloc(2*edgeNum*sizeof(ULONGINT_T), TH);
  else
    startVertex = (ULONGINT_T *) node_malloc((edgeNum+MAX_PARAL_EDGES*TOT_VERTICES)*sizeof(ULONGINT_T), TH);

  if (SCALE < 10)
    endVertex = (ULONGINT_T *) node_malloc(2*edgeNum*sizeof(ULONGINT_T), TH);
  else
    endVertex = (ULONGINT_T *) node_malloc((edgeNum+MAX_PARAL_EDGES*TOT_VERTICES)*sizeof(ULONGINT_T), TH);

  node_Barrier();
    
  on_thread(MYTHREAD) {
    for (j = i_edgeStartCounter[MYTHREAD]; j<i_edgeEndCounter[MYTHREAD]; j++) {
      startVertex[j] = startV[j-i_edgeStartCounter[MYTHREAD]];
      endVertex[j] = endV[j-i_edgeStartCounter[MYTHREAD]];
    }   
  }

  numEdgesPlacedInCliques = edgeNum;
   
  node_Barrier();

  /*--------------------------------------------------------------------------*/
  /* STEP 3: Connect the cliques.                                             */
  /*--------------------------------------------------------------------------*/
  
  i_edgePtr = 0;  
  node_Barrier();
  
  p = PROB_INTERCL_EDGES;

#if 0
  pardo(i_clique, 0, totCliques, 1) {

    /* Get current clique parameters */
    i_cliqueSize = cliqueSizes[i_clique];
    i_firstVsInClique = firstVsInCliques[i_clique];
  
    /* printf("Cliquenum %d, CliqueSize %d\n", i_clique, i_cliqueSize, i_firstVsInClique); */
    for(i = 0; i < i_cliqueSize; i++) {
      randNum = 1 + (ULONGINT_T) (isprng(stream) % i_cliqueSize);
       
      for (t=0; t<randNum; t++) {

	      tempVertex1 = firstVsInCliques[i_clique] + (ULONGINT_T) (isprng(stream) % i_cliqueSize);
	    
	      endVClique =  (ULONGINT_T) (isprng(stream) % totCliques);
	      if (endVClique != i_clique) {
	        tempVertex2 = firstVsInCliques[endVClique] + (int) (isprng(stream) % cliqueSizes[endVClique]);
	    
	        if (sprng(stream) <= p) {
	          startV[i_edgePtr] = tempVertex1;
	          endV[i_edgePtr] = tempVertex2;
	          i_edgePtr = i_edgePtr + 1;
	        }
        }
      }
    }

  }

  node_Barrier();

  /* Merge partial edge lists */

  on_thread(MYTHREAD) {
    i_edgeEndCounter[MYTHREAD] = i_edgePtr;     
    i_edgeStartCounter[MYTHREAD] = 0;
  }
  
  node_Barrier();
  
  on_one_thread {
    for (i=1; i<THREADS; i++) {
      i_edgeEndCounter[i] = i_edgeEndCounter[i-1] + i_edgeEndCounter[i];
      i_edgeStartCounter[i] = i_edgeEndCounter[i-1]; 
    }
  }  
    
  edgeNum = node_Reduce_d(i_edgePtr, SUM, TH);
  numEdgesPlacedOutside = edgeNum;

  on_thread(MYTHREAD) {
    for (j = i_edgeStartCounter[MYTHREAD]; j<i_edgeEndCounter[MYTHREAD]; j++) {
      startVertex[j+numEdgesPlacedInCliques] = startV[j-i_edgeStartCounter[MYTHREAD]];
      endVertex[j+numEdgesPlacedInCliques] = endV[j-i_edgeStartCounter[MYTHREAD]];
    }   
  }
#endif

  /* Generating inter-clique edges as given in the specs */
  pardo (i, 0, TOT_VERTICES, 1) {

    tempVertex1 = i;
    h = totCliques;
    l = 0;
    t = -1;
    while (h - l > 1) {
      m = (h+l)/2;
      if (tempVertex1>=firstVsInCliques[m])
        l = m;
      else {
        if ((tempVertex1<firstVsInCliques[m]) && (m>0)) {
          if (tempVertex1>=firstVsInCliques[m-1]) {
            t = m-1;
            break;
          }
          else
            h = m;
        }
      }
    }

    if (t == -1) {
      for (m = l+1; m<h; m++) {
        if (tempVertex1<firstVsInCliques[m])
          break;
      }
      t = m-1;
    }

    t1 = firstVsInCliques[t];

    for (d=1, p=PROB_INTERCL_EDGES; d<TOT_VERTICES; d = d*2, p = p/2) {

      if (sprng(stream) <= p) {

        tempVertex2 = (i+d) % TOT_VERTICES;

        h = totCliques;
        l = 0;
        t = -1;
        while (h - l > 1) {
          m = (h+l)/2;
          if (tempVertex2>=firstVsInCliques[m])
            l = m;
          else {
            if ((tempVertex2<firstVsInCliques[m]) && (m>0)) {
              if (firstVsInCliques[m-1] <= tempVertex2) {
                t = m-1;
                break;
              }
              else
                h = m;
            }
          }
        }

        if (t == -1) {
          for (m = l+1; m<h; m++) {
            if (tempVertex2<firstVsInCliques[m])
              break;
          }
          t = m-1;
        }
    
        t2 = firstVsInCliques[t];

        if (t1 != t2) {
          randNumEdges = isprng(stream) % MAX_PARAL_EDGES + 1;
          for (j=0; j<randNumEdges; j++) {
            startV[i_edgePtr] = tempVertex1;
            endV[i_edgePtr] = tempVertex2;
            i_edgePtr = i_edgePtr + 1;
          }
        }

      }

      if ((sprng(stream) <= p) && (i-d>=0)) {

        tempVertex2 = (i-d) % TOT_VERTICES;

        h = totCliques;
        l = 0;
        t = -1;
        while (h - l > 1) {
          m = (h+l)/2;
          if (tempVertex2>=firstVsInCliques[m])
            l = m;
          else {
            if ((tempVertex2<firstVsInCliques[m]) && (m>0)) {
              if (firstVsInCliques[m-1] <= tempVertex2) {
                t = m-1;
                break;
              }
              else
                h = m;
            }
          }
        }

        if (t == -1) {
          for (m = l+1; m<h; m++) {
            if (tempVertex2<firstVsInCliques[m])
              break;
          }
          t = m-1;
        }

        t2 = firstVsInCliques[t];

        if (t1 != t2) {
          randNumEdges = isprng(stream) % MAX_PARAL_EDGES + 1;
          for (j=0; j<randNumEdges; j++) {
            startV[i_edgePtr] = tempVertex1;
            endV[i_edgePtr] = tempVertex2;
            // printf("%d ", *currentEdgePtr);
            i_edgePtr = i_edgePtr + 1;
          }

        }

      }

    }
  }
    
  node_Barrier();

  on_thread(MYTHREAD) {
    i_edgeEndCounter[MYTHREAD] = i_edgePtr;
    i_edgeStartCounter[MYTHREAD] = 0;
  }

  node_Barrier();

  on_one_thread {
    for (i=1; i<THREADS; i++) {
      i_edgeEndCounter[i] = i_edgeEndCounter[i-1] + i_edgeEndCounter[i];
      i_edgeStartCounter[i] = i_edgeEndCounter[i-1];
    }
  }

  edgeNum = node_Reduce_d(i_edgePtr, SUM, TH);
  numEdgesPlacedOutside = edgeNum;

  on_thread(MYTHREAD) {
    for (j = i_edgeStartCounter[MYTHREAD]; j<i_edgeEndCounter[MYTHREAD]; j++) {
      startVertex[j+numEdgesPlacedInCliques] = startV[j-i_edgeStartCounter[MYTHREAD]];
      endVertex[j+numEdgesPlacedInCliques] = endV[j-i_edgeStartCounter[MYTHREAD]];
    }
  }

  node_Barrier();


  numEdgesPlaced = numEdgesPlacedInCliques + numEdgesPlacedOutside;
  
  SDGdataPtr->numEdgesPlaced = numEdgesPlaced;
  
  node_Barrier();
  
  on_one_thread {
    printf("Finished generating edges\n");
    printf("No. of intra-clique edges - %d\n", numEdgesPlacedInCliques);
    printf("No. of inter-clique edges - %d\n", numEdgesPlacedOutside);
    printf("Total no. of edges        - %d\n", numEdgesPlaced); 
  }
  
  node_free(i_edgeStartCounter, TH);
  node_free(i_edgeEndCounter, TH);
  node_Barrier();
  free(startV);
  free(endV);

  node_Barrier();
  
  node_free(cliqueSizes, TH);  
  node_free(firstVsInCliques, TH);
  node_free(lastVsInCliques, TH);

  node_Barrier();

  /*--------------------------------------------------------------------------*/
  /* STEP 4: Generate edge weights                                                 */
  /*--------------------------------------------------------------------------*/

  numStrWtEdges  = 0;
  SDGdataPtr->intWeight = (LONGINT_T *) node_malloc(numEdgesPlaced*sizeof(LONGINT_T), TH);
  p = PERC_INT_WEIGHTS;
  
  node_Barrier();
  
  pardo(i, 0, numEdgesPlaced, 1) {
    if (sprng(stream) <= p) {
      SDGdataPtr->intWeight[i] = 1 + (isprng(stream) % (MAX_INT_WEIGHT-1));
    } else {
      SDGdataPtr->intWeight[i] = -1;
      numStrWtEdges = numStrWtEdges + 1;
    }
  }

  node_Barrier();
  t = 0;
  on_one_thread {
    for (i=0; i<numEdgesPlaced; i++) {
      if (SDGdataPtr->intWeight[i] < 0) {
        SDGdataPtr->intWeight[i] = -t;
	      t = t + 1;
      }
    }
  } 
  
  node_Barrier();

  numStrWtEdges = node_Reduce_d(numStrWtEdges, SUM, TH);
  
  node_Barrier();
  
  SDGdataPtr->strWeight =  (char *) node_malloc(numStrWtEdges*MAX_STRLEN*sizeof(char), TH);
  
  node_Barrier();

  pardo(i, 0, numEdgesPlaced, 1) {
    if (SDGdataPtr->intWeight[i] <= 0) {
      for (j=0; j<MAX_STRLEN; j++) {
        SDGdataPtr->strWeight[(-SDGdataPtr->intWeight[i])*MAX_STRLEN+j] = (char) (1 + isprng(stream) % 127);
      }
    }
  }
  
  /* Choose SOUGHT STRING randomly if not assigned */
  if (strlen(SOUGHT_STRING) != MAX_STRLEN) {
    SOUGHT_STRING = (char *) node_malloc(MAX_STRLEN * sizeof(char), TH);
  }
  
  on_one_thread {  
    t = isprng(stream) % numStrWtEdges;
    for (j=0; j<MAX_STRLEN; j++) {
      SOUGHT_STRING[j] = (char) ((int) SDGdataPtr->strWeight[t*MAX_STRLEN+j]);
    }
  }

  node_Barrier();

  /*--------------------------------------------------------------------------*/
  /* STEP 5: Permute Vertices                                                 */
  /*--------------------------------------------------------------------------*/

  pardo(i, 0, numEdgesPlaced, 1) {
    /* printf("O:[%d %d] ", startVertex[i], endVertex[i]); */
    startVertex[i] = permV[(startVertex[i])];
    endVertex[i] = permV[(endVertex[i])];
    /* printf("[%d %d] ", startVertex[i], endVertex[i]); */
  }

  node_Barrier();

  /*--------------------------------------------------------------------------*/
  /* STEP 6: Sort Vertices                                           */
  /*--------------------------------------------------------------------------*/

  /* Radix sort with StartVertex as primary key */
  SDGdataPtr->startVertex = (ULONGINT_T *) node_malloc(numEdgesPlaced*sizeof(ULONGINT_T), TH);
  SDGdataPtr->endVertex = (ULONGINT_T *) node_malloc(numEdgesPlaced*sizeof(ULONGINT_T), TH);

  node_Barrier();

  all_radixsort_node_aux_s3(numEdgesPlaced, startVertex, SDGdataPtr->startVertex, endVertex, SDGdataPtr->endVertex, TH);

  node_Barrier();
  
  node_free(startVertex, TH);
  node_free(endVertex, TH);
  
  node_Barrier();
  

  if (SCALE <12) {
    /* Sort with endVertex as secondary key */
    on_one_thread {
      i0 = 0;
      i1 = 0;
      i = 0;
      while (i<numEdgesPlaced) {
        for (i=i0; i<numEdgesPlaced; i++) {
          if (SDGdataPtr->startVertex[i] != SDGdataPtr->startVertex[i1]) {
            i1 = i;
  	        break;
          }
        }

        for (j=i0; j<i1; j++) {
          for (k=j+1; k<i1; k++) {
            if (SDGdataPtr->endVertex[k] < SDGdataPtr->endVertex[j]) {
              t = SDGdataPtr->endVertex[j];
              SDGdataPtr->endVertex[j] = SDGdataPtr->endVertex[k];
              SDGdataPtr->endVertex[k] = t;
            }
          }
        }

        if (SDGdataPtr->startVertex[i0] != TOT_VERTICES-1)
          i0 = i1;
        else {
          for (j=i0; j<numEdgesPlaced; j++) {
            for (k=j+1; k<numEdgesPlaced; k++) {
        	    if (SDGdataPtr->endVertex[k] < SDGdataPtr->endVertex[j]) {
        	      t = SDGdataPtr->endVertex[j];
        	      SDGdataPtr->endVertex[j] = SDGdataPtr->endVertex[k];
        	      SDGdataPtr->endVertex[k] = t;
        	    }
            }
          }
        }
      }

    }
    tempIndex = (ULONGINT_T *) node_malloc(sizeof(ULONGINT_T), TH);
    node_Barrier();
  }
  else {

    tempIndex = (ULONGINT_T *) node_malloc((TOT_VERTICES+1)*sizeof(ULONGINT_T), TH);
    node_Barrier();

    /* Update degree of each vertex */
    tempIndex[0] = 0;
    tempIndex[TOT_VERTICES] = numEdgesPlaced;
    i0 = 0;
    i = 0;
    node_Barrier();
    on_one_thread {
      for (i=0; i<TOT_VERTICES; i++) {
        tempIndex[i+1] = tempIndex[i];
        for (j=i0; j<numEdgesPlaced; j++) {
          if (SDGdataPtr->startVertex[j] != SDGdataPtr->startVertex[i0]) {
            if (SDGdataPtr->startVertex[i0] == i) {
              tempIndex[i+1] = j;
              i0 = j;
              break;
            }
          }
        }
      }
    }

    node_Barrier();
    
    /* Insertion sort for now, replace with something better later on */
    pardo (i, 0, TOT_VERTICES, 1) {
      for (j=tempIndex[i]; j<tempIndex[i+1]; j++) {
        for (k=j+1; k<tempIndex[i+1]; k++) {
          if (SDGdataPtr->endVertex[k] < SDGdataPtr->endVertex[j]) {
            t = SDGdataPtr->endVertex[j];
            SDGdataPtr->endVertex[j] = SDGdataPtr->endVertex[k];
            SDGdataPtr->endVertex[k] = t;
          }
        }
      }
    }
    
  }
  
  node_Barrier();
  free(stream);
  node_free(permV, TH);
  node_free(tempIndex, TH);
  node_Barrier();
}

