#include <stdlib.h>
#include "simple.h"
#include "globals.h"
#include "defs.h"
#include "getUserParameters.h"
#include "getStartLists.h"
#include "genScalData.h"
#include "computeGraph.h"
#include "findSubGraphs.h"
#include "cutClusters.h"
 
void *SIMPLE_main(THREADED)
{
  /* tuple for Scalable Data Generation  */
  /* stores startVertex, endVertex, int weight and other info  */
  graphSDG* SDGdata;  
  
  /* The graph data structure for this benchmark - see defs.h */
  graph* G;
  
  /* Kernel 2 */
  edge *maxIntWtList, *soughtStrWtList;
  int maxIntWtListSize, soughtStrWtListSize;
  
  /* Kernel 3 */
  V *intWtVList, *strWtVList;
  Vl **intWtVLList, **strWtVLList;
  Vd *intWtVDList, *strWtVDList;
  
  
  /* Temp vars */
  double time;

  LONGINT_T i, j, k, tmpVar;
  ULONGINT_T numEdges, edgeCounter;
  Vl *currV, *tempV;
  
  FILE* outfp;

  outfp = fopen("results.txt", "w");
  
  /*-------------------------------------------------------------------------- */
  /*       Preamble -- Untimed                                                 */
  /*-------------------------------------------------------------------------- */
   
  /* User Interface: Configurable parameters, and global program control. */
  on_one_thread {
    printf("\nHPCS SSCA #2 Graph Analysis Executable Specification:");
    printf("\nRunning...\n\n");
  }

  getUserParameters();

  on_one_thread {
    printf("\nNumber of processors: %d\n", THREADS);
    printf("Problem Scale: %d\n\n", SCALE);
    fprintf(outfp, "No. of processors - %d\n", THREADS);
    fprintf(outfp, "No. of vertices   - 2^%d\n", SCALE);
  }

  /*-------------------------------------------------------------------------- */
  /*       Scalable Data Generator                                             */
  /*-------------------------------------------------------------------------- */

  on_one_thread {
    printf("\nScalable Data Generator - genScalData() beginning execution...\n");  
  }
  
  node_Barrier();
  
  time = get_seconds();
  
  SDGdata  = (graphSDG*) node_malloc(sizeof(graphSDG), TH);

  /* Generate edge data as per the Written Specification. */
  genScalData(SDGdata, TH);
  
  time  = get_seconds()-time;

  on_one_thread {
    printf("\nTime taken for Scalable Data Generation is %9.6f sec.\n\n", time);
    fprintf(outfp, "\nScalable Data Generation - %9.6f sec.\n", time);
    printf("\n\tgenScalData() completed execution.\n");
  }
  

  /*-------------------------------------------------------------------------- */
  /* Kernel 1 - Graph Construction.                                            */
  /*-------------------------------------------------------------------------- */
  
  /* From the input edges, construct the graph 'G'.  */
  on_one_thread {
    printf("\nKernel 1 - computeGraph() beginning execution...\n");
  }
  
  node_Barrier();
  
  time = get_seconds(); 
  
  /* Memory allocated for the graph data structure */
  G = (graph *) node_malloc(sizeof(graph), TH);
     
  computeGraph(G, SDGdata, TH);
  
  time  = get_seconds()-time;
  
  on_one_thread {
    printf("\n\tcomputeGraph() completed execution.\n");
    printf("\nTime taken for kernel 1 is %9.6f sec.\n", time);
    fprintf(outfp, "\nKernel 1 - %9.6f sec.\n", time);
  }

  /*-------------------------------------------------------------------------- */
  /* Kernel 2 - Find Max weight and sought string                              */
  /*-------------------------------------------------------------------------- */
  
  on_one_thread {
    printf("\nKernel 2 - getStartLists() beginning execution...\n");
  }
  
  node_Barrier();

  time = get_seconds();

  /* Initialize vars */

  maxIntWtListSize = 0;
  soughtStrWtListSize = 0;
  
  /* Allocate mem. for the start lists */
  maxIntWtList = (edge *) node_malloc(sizeof(edge), TH);
  soughtStrWtList = (edge *) node_malloc(sizeof(edge), TH);

  node_Barrier();

  getStartLists(G, &maxIntWtList, &maxIntWtListSize, &soughtStrWtList, &soughtStrWtListSize, TH);
  
  time  = get_seconds()-time;
  
  on_one_thread {
    printf("\n\tgetStartLists() completed execution.\n");
    printf("\nTime taken for kernel 2 is %9.6f sec.\n\n", time);
    fprintf(outfp, "\nKernel 2 - %9.6f sec.\n", time);
  }
  
  node_Barrier();
  
  /*-------------------------------------------------------------------------- */
  /* Kernel 3 - Graph Extraction                                               */
  /*-------------------------------------------------------------------------- */
  
  on_one_thread {
    printf("\nKernel 3 - findSubGraphs() beginning execution...\n");
  }
  
  time = get_seconds();

  if (K3_DS == 0) {
    
    /* Allocate memory for subGraph data structures */

    intWtVList = (V *) node_malloc(G->numVertices*maxIntWtListSize*sizeof(V), TH);
  
    strWtVList = (V *) node_malloc(G->numVertices*soughtStrWtListSize*sizeof(V), TH);

    node_Barrier();
  
    findSubGraphs0(G, intWtVList, strWtVList, maxIntWtList, maxIntWtListSize, soughtStrWtList, soughtStrWtListSize, TH);
  
    node_Barrier();
    
    time  = get_seconds()-time;
  
  }

  if (K3_DS == 1) {

    /* Allocate memory for subGraph data structures */

    intWtVLList = (Vl **) node_malloc(maxIntWtListSize*sizeof(Vl*), TH);

    strWtVLList = (Vl **) node_malloc(soughtStrWtListSize*sizeof(Vl*), TH);

    node_Barrier();

    findSubGraphs1(G, intWtVLList, strWtVLList, maxIntWtList, maxIntWtListSize, soughtStrWtList, soughtStrWtListSize, TH);

    node_Barrier();

    time  = get_seconds()-time;
    
    /*  Verification
    on_one_thread {
      for (i=0; i<maxIntWtListSize; i++) {
        printf("%d -- ", i); 
        currV = intWtVLList[i];
        while (currV != NULL) {
          printf("[%d %d] ", currV->num, currV->depth);
          currV = currV->next;
        }
        printf("\n");  
      }

      for (i=0; i<soughtStrWtListSize; i++) {
        printf("%d -- ", i); 
        currV = strWtVLList[i];
        while (currV != NULL) {
          printf("[%d %d] ", currV->num, currV->depth);
          currV = currV->next;
        }
        printf("\n");  
      }

    }

   */

  }

  if (K3_DS == 2) {

    /* Allocate memory for subGraph data structures */

    intWtVDList = (Vd *) node_malloc(maxIntWtListSize*sizeof(Vd), TH);

    strWtVDList = (Vd *) node_malloc(soughtStrWtListSize*sizeof(Vd), TH);

    node_Barrier();

    findSubGraphs2(G, intWtVDList, strWtVDList, maxIntWtList, maxIntWtListSize, soughtStrWtList, soughtStrWtListSize, TH);

    node_Barrier();

    time  = get_seconds()-time; 
    
    /* Verification */
    /*
    on_one_thread {
      printf("\nInt weight sub-graphs \n");
      for (i=0; i<maxIntWtListSize; i++) {
        printf("%d -- ", i);
        for (j=0; j<intWtVDList[i].numArrays; j++) {
          printf("\n [Array %d] - \n", j);
          for (k=0; k<intWtVDList[i].arraySize[j]; k++) {
            printf("[%d %d] ", intWtVDList[i].vList[j][k].num, intWtVDList[i].vList[j][k].depth);          
          }           
        
        }
        printf("\n");  
      }
      
      printf("\nStr weight sub-graphs \n");
      for (i=0; i<soughtStrWtListSize; i++) {
        printf("%d -- ", i);
        for (j=0; j<strWtVDList[i].numArrays; j++) {
          printf("\n [Array %d] - \n", j);
          for (k=0; k<strWtVDList[i].arraySize[j]; k++) {
            printf("[%d %d] ", strWtVDList[i].vList[j][k].num, strWtVDList[i].vList[j][k].depth);          
          }           
        
        }
        printf("\n");  
      }

    }
   */
    
  }
  
  
  on_one_thread {
    printf("\n\tfindSubGraphs() completed execution.\n");
    printf("\nTime taken for kernel 3 is %9.6f sec.\n\n", time);
    fprintf(outfp, "\nKernel 3 - %9.6f sec.\n", time);
  }

  node_Barrier();
  
  /*-------------------------------------------------------------------------- */
  /* Kernel 4 - Graph Clustering                                               */
  /*-------------------------------------------------------------------------- */
  
  on_one_thread {
    printf("\nKernel 4 - cutClusters() beginning execution...\n");  
  }
    
  time = get_seconds();  
    
  cutClusters(G, TH);
  
  node_Barrier();
  
  time  = get_seconds()-time;
  
  on_one_thread {
    printf("\n\tcutClusters() completed execution.\n");
  }
    
  on_one_thread {
    printf("\nTime taken for Kernel 4 is %9.6f sec.\n\n", time);
    fprintf(outfp, "\nKernel 4 - %9.6f sec.\n", time);
  }
    
  node_Barrier();
  
  
  /**************************************************************************/

  fclose(outfp);
    
  /* Free shared memory vars */
  node_free(G->outDegree, TH);
  node_free(G->outVertexIndex, TH);
  node_free(G->outVertexList, TH);
  node_free(G->paralEdgeIndex, TH);
  node_free(G->inDegree, TH);
  node_free(G->inVertexIndex, TH);
  node_free(G->inVertexList, TH);
  node_free(G->intWeight, TH);
  node_free(G->strWeight, TH);

  if (K3_DS == 0) {
    node_free(strWtVList, TH);
    node_free(intWtVList, TH);
  }

  if (K3_DS == 1) {
    on_one_thread {
      for (i=0; i<maxIntWtListSize; i++) {
        currV = intWtVLList[i];
        while (currV != NULL) {
          tempV = currV->next;
          free(currV);
          currV = tempV;
        }  
      }
      for (i=0; i<soughtStrWtListSize; i++) {
        currV = strWtVLList[i];
        while (currV != NULL) {
          tempV = currV->next;
          free(currV);
          currV = tempV;
        }  
      }    
    }
        
    node_Barrier();
      
    node_free(strWtVLList, TH);
    node_free(intWtVLList, TH);
  }

  if (K3_DS == 2) {
    
    on_one_thread {
      for (i=0; i<maxIntWtListSize; i++) {
        for (j=0; j<intWtVDList[i].numArrays; j++) {
          free(intWtVDList[i].vList[j]);
        }
        free(intWtVDList[i].vList);  
        free(intWtVDList[i].arraySize);
      }  
      for (i=0; i<soughtStrWtListSize; i++) {
        for (j=0; j<strWtVDList[i].numArrays; j++) {
          free(strWtVDList[i].vList[j]);
        }
        free(strWtVDList[i].vList);
        free(strWtVDList[i].arraySize);  
      }   
    }
    
    node_Barrier();
    
    node_free(strWtVDList, TH);
    node_free(intWtVDList, TH);
  }

  node_free(soughtStrWtList, TH);
  node_free(maxIntWtList, TH);
  node_free(SOUGHT_STRING, TH);
  node_free(G, TH);
  node_free(SDGdata, TH);
      
  SIMPLE_done(TH);

}
