#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "simple.h"
#include "globals.h"
#include "getStartLists.h"


void getStartLists(graph* GPtr, edge** maxIntWtListPtr, int* maxIntWtListSize, edge** soughtStrWtListPtr, int* soughtStrWtListSize, THREADED)
{
  LONGINT_T maxWeight, i, j, t;
  edge* tmpEdgeList;
  edge *maxIntWtList, *soughtStrWtList;
  int i_edgeCounter, *i_edgeStartCounter, *i_edgeEndCounter;

  maxWeight = 0;

  /* Find Max Wt on each thread */
  pardo(i, 0, GPtr->numEdges, 1)
    if (GPtr->intWeight[i] > maxWeight)
      maxWeight = GPtr->intWeight[i];
  
  node_Barrier();

  maxWeight = node_Reduce_d(maxWeight, MAX, TH);
  node_Barrier();

  /* Create partial lists */
  
  /* Allocate mem. for temp edge list for each thread */
  tmpEdgeList = (edge *) malloc((5+ceil(1.5*(GPtr->numIntEdges)/MAX_INT_WEIGHT))*sizeof(edge));
      
  i_edgeCounter = 0;

  node_Barrier();
    
  pardo (i, 0, GPtr->numEdges, 1) {
    if (GPtr->intWeight[i] == maxWeight) {
      /* Find the corresponding endVertex */
      for (j=0; j<GPtr->numDirectedEdges; j++) {
        if (GPtr->paralEdgeIndex[j]>i)
          break;
      }
      tmpEdgeList[i_edgeCounter].endVertex = GPtr->outVertexList[j-1];
      tmpEdgeList[i_edgeCounter].edgeNum = j-1;

      for (t=0; t<GPtr->numVertices; t++) {
        if (GPtr->outVertexIndex[t]>j-1)
           break;
      }
      tmpEdgeList[i_edgeCounter].startVertex = t-1;
      i_edgeCounter = i_edgeCounter + 1;
    }
  }
  
  node_Barrier();

  /* Merge partial edge lists */
  
  i_edgeStartCounter = (int *) node_malloc(THREADS*sizeof(int), TH);
  i_edgeEndCounter = (int *) node_malloc(THREADS*sizeof(int), TH);
  
  on_thread(MYTHREAD) {
    i_edgeEndCounter[MYTHREAD] = i_edgeCounter;     
    i_edgeStartCounter[MYTHREAD] = 0;
  }
  
  node_Barrier();

  on_one_thread {
    for (i=1; i<THREADS; i++) {
      i_edgeEndCounter[i] = i_edgeEndCounter[i-1] + i_edgeEndCounter[i];
      i_edgeStartCounter[i] = i_edgeEndCounter[i-1]; 
    }
  }  
  
  node_Barrier();
  
  *maxIntWtListSize = node_Reduce_d(i_edgeCounter, SUM, TH);
  
  node_Barrier();

  node_free(*maxIntWtListPtr, TH);

  maxIntWtList = (edge *) node_malloc((*maxIntWtListSize)*sizeof(edge), TH);

  on_thread(MYTHREAD) {
    for (j = i_edgeStartCounter[MYTHREAD]; j<i_edgeEndCounter[MYTHREAD]; j++) {
      (maxIntWtList[j]).startVertex = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].startVertex;
      (maxIntWtList[j]).endVertex = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].endVertex;
      (maxIntWtList[j]).edgeNum = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].edgeNum;
    } 
  } 

  node_Barrier();

  *maxIntWtListPtr = maxIntWtList;

  node_Barrier();

  i_edgeCounter = 0;

  pardo (i, 0, GPtr->numStrEdges, 1) {
    if (strncmp(GPtr->strWeight+i*MAX_STRLEN, SOUGHT_STRING, MAX_STRLEN) == 0) {
      /* Find the corresponding endVertex */
      for (t=0; t<GPtr->numEdges; t++) {
        if (GPtr->intWeight[t]==-i)
	      break;
      }
      for (j=0; j<GPtr->numDirectedEdges; j++) {
        if (GPtr->paralEdgeIndex[j]>t)
          break;
      }
      tmpEdgeList[i_edgeCounter].endVertex = GPtr->outVertexList[j-1];
      tmpEdgeList[i_edgeCounter].edgeNum = j-1;

      for (t=0; t<GPtr->numVertices; t++) {
        if (GPtr->outVertexIndex[t]>j-1)
          break;
      }
      tmpEdgeList[i_edgeCounter].startVertex = t-1;
      i_edgeCounter = i_edgeCounter + 1;
    }
  }
  
 
  node_Barrier();

  on_thread(MYTHREAD) {
    i_edgeEndCounter[MYTHREAD] = i_edgeCounter;     
    i_edgeStartCounter[MYTHREAD] = 0;
  }
  
  node_Barrier();

  on_one_thread {
    for (i=1; i<THREADS; i++) {
      i_edgeEndCounter[i] = i_edgeEndCounter[i-1] + i_edgeEndCounter[i];
      i_edgeStartCounter[i] = i_edgeEndCounter[i-1]; 
    }
  }  
  
  node_Barrier();
  
  *soughtStrWtListSize = node_Reduce_d(i_edgeCounter, SUM, TH);
  
  node_Barrier();

  node_free(*soughtStrWtListPtr, TH);

  soughtStrWtList = (edge*) node_malloc((*soughtStrWtListSize)*sizeof(edge), TH);
    
  on_thread(MYTHREAD) {
    for (j = i_edgeStartCounter[MYTHREAD]; j<i_edgeEndCounter[MYTHREAD]; j++) {
      (soughtStrWtList[j]).startVertex = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].startVertex;
      (soughtStrWtList[j]).endVertex = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].endVertex;
      (soughtStrWtList[j]).edgeNum = tmpEdgeList[j-i_edgeStartCounter[MYTHREAD]].edgeNum;
    } 
  } 

  node_Barrier();

  *soughtStrWtListPtr = soughtStrWtList;
  
  /* Free temp vars */
  free(tmpEdgeList);
  node_free(i_edgeStartCounter, TH);
  node_free(i_edgeEndCounter, TH);
  
  node_Barrier();
}
