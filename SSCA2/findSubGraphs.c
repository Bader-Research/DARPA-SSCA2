#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "simple.h"
#include "globals.h"
#include "findSubGraphs.h"

void findSubGraphs0 (graph* GPtr, V* intWtVList, V* strWtVList, edge* maxIntWtList,
int maxIntWtListSize, edge* soughtStrWtList, int soughtStrWtListSize, THREADED) {
  ULONGINT_T i, j, k, t, count, depth, verticesVisited, currIndex;
  LONGINT_T *visited;
  
 
  pardo(i, 0, maxIntWtListSize+soughtStrWtListSize, 1) {
    if (i < maxIntWtListSize) {
      for (j=0; j<GPtr->numVertices; j++) {
        intWtVList[i*(GPtr->numVertices)+j].num = 0;
        intWtVList[i*(GPtr->numVertices)+j].depth = 0;
      }  
    }
    else {
      t = i - maxIntWtListSize;
      for (j=0; j<GPtr->numVertices; j++) {
        strWtVList[t*(GPtr->numVertices)+j].num = 0;
        strWtVList[t*(GPtr->numVertices)+j].depth = 0;
      }
    }
  }
   
  node_Barrier();
  
  visited = (LONGINT_T *) malloc(GPtr->numVertices*sizeof(LONGINT_T));

  node_Barrier();
  
  /* Each thread runs a BFS from endvertex of maxIntWtList edgeList */
  
  pardo(i, 0, maxIntWtListSize+soughtStrWtListSize, 1) {
    
    for (k=0; k<GPtr->numVertices; k++)
      visited[k] = 0;
    
    if (i<maxIntWtListSize) {
    
      intWtVList[i*(GPtr->numVertices)+0].num = maxIntWtList[i].startVertex;
      intWtVList[i*(GPtr->numVertices)+0].depth = -1;

      intWtVList[i*(GPtr->numVertices)+1].num = maxIntWtList[i].endVertex;
      intWtVList[i*(GPtr->numVertices)+1].depth = 1;

      visited[(intWtVList[i*(GPtr->numVertices)+0]).num] = 1;
      visited[(intWtVList[i*(GPtr->numVertices)+1]).num] = 1;

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {
        depth = intWtVList[i*(GPtr->numVertices)+currIndex].depth+1;
        for (j=GPtr->outVertexIndex[intWtVList[i*(GPtr->numVertices)+currIndex].num]; j<GPtr->outVertexIndex[intWtVList[i*(GPtr->numVertices)+currIndex].num]+GPtr->outDegree[intWtVList[i*(GPtr->numVertices)+currIndex].num]; j++) {
          if (visited[GPtr->outVertexList[j]] == 0) {
            visited[GPtr->outVertexList[j]] = 1;
            intWtVList[i*(GPtr->numVertices)+verticesVisited].num = GPtr->outVertexList[j];
            intWtVList[i*(GPtr->numVertices)+verticesVisited].depth = depth;
            verticesVisited = verticesVisited + 1;
          }
        }
        if ((currIndex < verticesVisited - 1) && (verticesVisited <
        GPtr->numVertices)){
          currIndex++;
	        depth = intWtVList[i*(GPtr->numVertices)+currIndex].depth;
        } else
          break;
      }
    }

    else { 
  
      t = i - maxIntWtListSize;  

      strWtVList[t*(GPtr->numVertices)+0].num = (soughtStrWtList[t]).startVertex;
      strWtVList[t*(GPtr->numVertices)+0].depth = -1;

      strWtVList[t*(GPtr->numVertices)+1].num = (soughtStrWtList[t]).endVertex;
      strWtVList[t*(GPtr->numVertices)+1].depth = 1;

      visited[(strWtVList[t*(GPtr->numVertices)+0]).num] = 1;
      visited[(strWtVList[t*(GPtr->numVertices)+1]).num] = 1;

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {
        depth = strWtVList[t*(GPtr->numVertices)+currIndex].depth+1;
        for (j=GPtr->outVertexIndex[strWtVList[t*(GPtr->numVertices)+currIndex].num];
        j<GPtr->outVertexIndex[strWtVList[t*(GPtr->numVertices)+currIndex].num]+GPtr->outDegree[strWtVList[t*(GPtr->numVertices)+currIndex].num]; j++) {
          if (visited[GPtr->outVertexList[j]] == 0) {
            visited[GPtr->outVertexList[j]] = 1;
            strWtVList[t*(GPtr->numVertices)+verticesVisited].num = GPtr->outVertexList[j];
            strWtVList[t*(GPtr->numVertices)+verticesVisited].depth = depth;
            verticesVisited = verticesVisited + 1;
          }
        }
        if (currIndex < verticesVisited - 1) {
          currIndex++;
          depth = strWtVList[t*(GPtr->numVertices)+currIndex].depth;
        } else

          break;
      }
    }
  }
  
  node_Barrier();
  
  /*
  on_one_thread {
  for(i=0; i<soughtStrWtListSize;i++) {
    for (j=0; j<GPtr->numVertices; j++) {
      printf("[%d %d] ", strWtVList[i][j].num, strWtVList[i][j].depth);
    }
    printf("\n");
  }
  }
  
  node_Barrier();
  */
  
  /*
  on_one_thread {
  for(i=0; i<maxIntWtListSize; i++) {
    for (j=0; j<GPtr->numVertices; j++) {
      printf("[%d %d] ", intWtVList[i][j].num, intWtVList[i][j].depth);
    }
    printf("\n");
  }
  }

  node_Barrier();
  */
    
  
  free(visited);
  
}

void findSubGraphs1 (graph* GPtr, Vl** intWtVLList, Vl** strWtVLList, edge* maxIntWtList,
int maxIntWtListSize, edge* soughtStrWtList, int soughtStrWtListSize, THREADED) {

  Vl *currV, *startV;

  ULONGINT_T i, j, k, t, depth, num, currIndex, verticesVisited;
  char* visited;

  visited = (char *) malloc(GPtr->numVertices*sizeof(char));


  pardo(i, 0, maxIntWtListSize+soughtStrWtListSize, 1) {

    for (k=0; k<GPtr->numVertices; k++)  {
      visited[k] = 'u';
    }
    
    if (i <maxIntWtListSize) {

      intWtVLList[i] = (Vl *) malloc(sizeof(Vl));
      (intWtVLList[i])->num = maxIntWtList[i].startVertex;
      (intWtVLList[i])->depth = 0;
      visited[(intWtVLList[i])->num] = 'v';

      (intWtVLList[i])->next = (Vl *) malloc(sizeof(Vl));
      ((intWtVLList[i])->next)->num = maxIntWtList[i].endVertex;
      ((intWtVLList[i])->next)->depth = 1;
      visited[((intWtVLList[i])->next)->num] = 'v';

      currV = (intWtVLList[i])->next;
      startV = (intWtVLList[i])->next;

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((startV->depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {

        depth = startV->depth+1;
        
        for (j=GPtr->outVertexIndex[startV->num];
        j<GPtr->outVertexIndex[startV->num]+GPtr->outDegree[startV->num]; j++) {
          if (visited[GPtr->outVertexList[j]] == 'u') {

            visited[GPtr->outVertexList[j]] = 'v';

            currV->next = (Vl *) malloc(sizeof(Vl));
            (currV->next)->num = GPtr->outVertexList[j];
            (currV->next)->depth = depth;

            verticesVisited = verticesVisited + 1;
            currV = currV->next;
            
          }
        }

        if ((currIndex < verticesVisited - 1) && (verticesVisited < GPtr->numVertices)) {
          currIndex++;
          startV = startV->next;

        } else {
          break;
        }
      }
      currV->next = NULL;

    }
    
    else {
      t = i - maxIntWtListSize;
      strWtVLList[t] = (Vl *) malloc(sizeof(Vl));
      (strWtVLList[t])->num = soughtStrWtList[t].startVertex;
      (strWtVLList[t])->depth = 0;
      visited[(strWtVLList[t])->num] = 'v';

      (strWtVLList[t])->next = (Vl *) malloc(sizeof(Vl));
      ((strWtVLList[t])->next)->num = soughtStrWtList[t].endVertex;
      ((strWtVLList[t])->next)->depth = 1;
      visited[((strWtVLList[t])->next)->num] = 'v';

      currV = (strWtVLList[t])->next;
      startV = (strWtVLList[t])->next;

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((startV->depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {

        depth = startV->depth+1;

        for (j=GPtr->outVertexIndex[startV->num];
        j<GPtr->outVertexIndex[startV->num]+GPtr->outDegree[startV->num]; j++) {
          if (visited[GPtr->outVertexList[j]] == 'u') {
            visited[GPtr->outVertexList[j]] = 'v';

            currV->next = (Vl *) malloc(sizeof(Vl));
            (currV->next)->num = GPtr->outVertexList[j];
            (currV->next)->depth = depth;

            verticesVisited = verticesVisited + 1;
            currV = currV->next;

          }
        }

        if ((currIndex < verticesVisited - 1) && (verticesVisited < GPtr->numVertices)) {
          currIndex++;
          startV = startV->next;
        } else {
         
          break;
        }
      }
      currV->next = NULL;

    }
  }

  free(visited);

}

void findSubGraphs2 (graph* GPtr, Vd* intWtVDList, Vd* strWtVDList, edge* maxIntWtList,
int maxIntWtListSize, edge* soughtStrWtList, int soughtStrWtListSize, THREADED) {

  ULONGINT_T i, j, k, t, count, depth, verticesVisited, currIndex, arraySize, vNum;
  char *visited;
  ULONGINT_T NS;
  
  /* No. of sub-arrays */
  NS = 30;
  arraySize = 5*MAX_CLUSTER_SIZE;
  
  pardo(i, 0, maxIntWtListSize+soughtStrWtListSize, 1) {
    if (i < maxIntWtListSize) {
      /* Initialize the DS and create one sub-array */
      intWtVDList[i].numArrays = 1;
      intWtVDList[i].arraySize = (ULONGINT_T *) malloc(NS*sizeof(ULONGINT_T));
      intWtVDList[i].arraySize[0] = 0;
      intWtVDList[i].vList = (V **) malloc(NS*sizeof(V *));
      intWtVDList[i].vList[0] = (V *) malloc(arraySize*sizeof(V));
      for (j=0; j< arraySize; j++) {
        intWtVDList[i].vList[0][j].num = 0;
        intWtVDList[i].vList[0][j].depth = 0;
      }

    }
    else {
      t = i - maxIntWtListSize;
      strWtVDList[t].numArrays = 1;
      strWtVDList[t].arraySize = (ULONGINT_T *) malloc(NS*sizeof(ULONGINT_T));
      strWtVDList[t].arraySize[0] = 0;
      strWtVDList[t].vList = (V **) malloc(NS*sizeof(V *));
      strWtVDList[t].vList[0] = (V *) malloc(arraySize*sizeof(V));
      for (j=0; j< arraySize; j++) {
        strWtVDList[t].vList[0][j].num = 0;
        strWtVDList[t].vList[0][j].depth = 0;
      }
    }
  }

  node_Barrier();

  visited = (char *) malloc(GPtr->numVertices*sizeof(char));

  node_Barrier();

  /* Each thread runs a BFS from endvertex of maxIntWtList edgeList */
  
  pardo(i, 0, maxIntWtListSize+soughtStrWtListSize, 1) {

    for (k=0; k<GPtr->numVertices; k++)
      visited[k] = 'u';

    if (i<maxIntWtListSize) {

      intWtVDList[i].vList[0][0].num = maxIntWtList[i].startVertex;
      intWtVDList[i].vList[0][0].depth = -1;

      intWtVDList[i].vList[0][1].num = maxIntWtList[i].endVertex;
      intWtVDList[i].vList[0][1].depth = 1;
      
      intWtVDList[i].arraySize[0] = 2;
      
      visited[intWtVDList[i].vList[0][0].num] = 'v';
      visited[intWtVDList[i].vList[0][1].num] = 'v';      

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {
        vNum = intWtVDList[i].vList[currIndex/arraySize][currIndex%arraySize].num;
        depth = intWtVDList[i].vList[currIndex/arraySize][currIndex%arraySize].depth + 1;
        for (j=GPtr->outVertexIndex[vNum]; j<GPtr->outVertexIndex[vNum]+GPtr->outDegree[vNum]; j++) {
          if (visited[GPtr->outVertexList[j]] == 'u') {
            visited[GPtr->outVertexList[j]] = 'v';
            intWtVDList[i].vList[verticesVisited/arraySize][verticesVisited%arraySize].num = GPtr->outVertexList[j];
            intWtVDList[i].vList[verticesVisited/arraySize][verticesVisited%arraySize].depth = depth;
            intWtVDList[i].arraySize[verticesVisited/arraySize] ++;  
            verticesVisited = verticesVisited + 1;
          }
        }
        
        /* Check if we need to create a new array */
        if (((float) verticesVisited / (float) arraySize) > 0.5) {
          /* create a new sub-array */
          if (intWtVDList[i].numArrays != (verticesVisited/arraySize + 2)) {
            intWtVDList[i].numArrays++;
            intWtVDList[i].vList[intWtVDList[i].numArrays-1] = (V *) malloc(arraySize*sizeof(V));
            intWtVDList[i].arraySize[intWtVDList[i].numArrays-1] = 0;
          }  
        }
        
        if ((currIndex < verticesVisited - 1) && (verticesVisited < GPtr->numVertices)){
          currIndex++;
	  depth = intWtVDList[i].vList[currIndex/arraySize][currIndex%arraySize].depth;
        } else
          break;
      }
    }

    else {
      
      t = i - maxIntWtListSize;
      
      strWtVDList[t].vList[0][0].num = soughtStrWtList[t].startVertex;
      strWtVDList[t].vList[0][0].depth = -1;

      strWtVDList[t].vList[0][1].num = soughtStrWtList[t].endVertex;
      strWtVDList[t].vList[0][1].depth = 1;
      
      strWtVDList[t].arraySize[0] = 2;
      
      visited[strWtVDList[t].vList[0][0].num] = 'v';
      visited[strWtVDList[t].vList[0][1].num] = 'v';      

      depth = 1;
      verticesVisited = 2;
      currIndex = 1;

      while ((depth < SUBGR_EDGE_LENGTH) || (verticesVisited == GPtr->numVertices)) {
        vNum = strWtVDList[t].vList[currIndex/arraySize][currIndex%arraySize].num;
        depth = strWtVDList[t].vList[currIndex/arraySize][currIndex%arraySize].depth + 1;
        for (j=GPtr->outVertexIndex[vNum]; j<GPtr->outVertexIndex[vNum]+GPtr->outDegree[vNum]; j++) {
          if (visited[GPtr->outVertexList[j]] == 'u') {
            visited[GPtr->outVertexList[j]] = 'v';
            strWtVDList[t].vList[verticesVisited/arraySize][verticesVisited%arraySize].num = GPtr->outVertexList[j];
            strWtVDList[t].vList[verticesVisited/arraySize][verticesVisited%arraySize].depth = depth;
            strWtVDList[t].arraySize[verticesVisited/arraySize] ++;  
            verticesVisited = verticesVisited + 1;
          }
        }
        
        /* Check if we need to create a new array */
        if (((float) verticesVisited / (float) arraySize) > 0.5) {
          /* create a new sub-array */
          if (strWtVDList[t].numArrays != (verticesVisited/arraySize + 2)) {
            strWtVDList[t].numArrays++;
            strWtVDList[t].vList[strWtVDList[t].numArrays-1] = (V *) malloc(arraySize*sizeof(V));
            strWtVDList[t].arraySize[strWtVDList[t].numArrays-1] = 0;
          }  
        }
        
        if ((currIndex < verticesVisited - 1) && (verticesVisited < GPtr->numVertices)){
          currIndex++;
	  depth = strWtVDList[t].vList[currIndex/arraySize][currIndex%arraySize].depth;
        } else
          break;
      }
      
    }
  }

  node_Barrier();

  free(visited);



}

