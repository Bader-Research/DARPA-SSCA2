#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "simple.h"
#include "globals.h"
#include "computeGraph.h"
#include "lock.h"

void computeGraph(graph* GPtr, graphSDG* SDGdataPtr, THREADED) {
  ULONGINT_T i, i0, i1, j, j0, k, r, t, v;
  ULONGINT_T startVertex, edgeCounter, strWeightEdgeNum, maxNumVertices, numEdgesPlaced;
  ULONGINT_T outVertexListSize, inVertexListSize, outVertexListNum, inVertexListNum;
  ULONGINT_T *impliedEdgeList;
  ULONGINT_T **auxArr;
  LOCK_T* vLock;
  
  maxNumVertices = 0;
  numEdgesPlaced = SDGdataPtr->numEdgesPlaced;

  /* First determine the number of vertices by scanning the tuple startVertex list */
  pardo (i, 0, numEdgesPlaced, 1) {
    if (SDGdataPtr->startVertex[i] > maxNumVertices) {
      maxNumVertices = SDGdataPtr->startVertex[i];
    }
  }
  
  node_Barrier();

  maxNumVertices = node_Reduce_d(maxNumVertices, MAX, TH) + 1;

  GPtr->numVertices = maxNumVertices;
  GPtr->numEdges    = numEdgesPlaced;

  GPtr->intWeight = SDGdataPtr->intWeight;
  
  GPtr->strWeight = SDGdataPtr->strWeight;
  
  strWeightEdgeNum = 0;
  
  on_one_thread {
    for (i=0; i<numEdgesPlaced; i++) {
      if (GPtr->intWeight[numEdgesPlaced-i-1] <0) {
        GPtr->numStrEdges = - (GPtr->intWeight[numEdgesPlaced-i-1]) + 1;
        GPtr->numIntEdges = numEdgesPlaced - GPtr->numStrEdges;
        break;
      }
    }
  }
  
  node_Barrier();
  

  GPtr->outDegree  = (SHORTINT_T *) node_malloc((GPtr->numVertices)*sizeof(SHORTINT_T), TH);

  GPtr->outVertexIndex = (ULONGINT_T *) node_malloc((GPtr->numVertices)*sizeof(ULONGINT_T), TH);

  pardo(i, 0, GPtr->numVertices, 1) {
    GPtr->outDegree[i] = 0;
    GPtr->outVertexIndex[i] = 0;
  }
  
  outVertexListSize = 0;
  
  node_Barrier();
  
  i0 = -1;
  pardo(i, 0, GPtr->numVertices, 1) {
    k = i;
    if ((outVertexListSize == 0) && (k != 0)) {
      while (i0 == -1) {
        for (j=0; j<numEdgesPlaced; j++) {
          if (k == SDGdataPtr->startVertex[j]) {
            i0 = j;
            break;
          }

        }
        k = k-1;
      }

    }

    if ((outVertexListSize == 0) && (k == 0)) {
      i0 = 0;
    }

    for (j=i0; j<numEdgesPlaced; j++) {

      if (i == GPtr->numVertices-1) {
        break;
      }

      
      if ((i != SDGdataPtr->startVertex[j])) {
        if ((j>0) && (i == SDGdataPtr->startVertex[j-1])) {
          if (j-i0 >= 1) {
            outVertexListSize = outVertexListSize + 1;
            GPtr->outDegree[i] = GPtr->outDegree[i] + 1;
            for (t=i0+1; t<j; t++) {
              if (SDGdataPtr->endVertex[t] != SDGdataPtr->endVertex[t-1]) {
                outVertexListSize = outVertexListSize + 1;
                GPtr->outDegree[i] = GPtr->outDegree[i]+1;
              }
            }
          }
        }
        i0 = j;
        break;
      }
    }

    if (i == GPtr->numVertices-1) {
      if (numEdgesPlaced-i0 >= 0) {
        outVertexListSize = outVertexListSize + 1;
        GPtr->outDegree[i] = GPtr->outDegree[i] + 1;
        for (t=i0+1; t<numEdgesPlaced; t++) {
          if (SDGdataPtr->endVertex[t] != SDGdataPtr->endVertex[t-1]) {
            outVertexListSize = outVertexListSize + 1;
            GPtr->outDegree[i] = GPtr->outDegree[i]+1;
          }
        }

      }

    }
    
  }

  node_Barrier();

  prefix_sums(GPtr->outVertexIndex, GPtr->outDegree, GPtr->numVertices, TH);

  node_Barrier();

  outVertexListSize = node_Reduce_d(outVertexListSize, SUM, TH);
  
  node_Barrier();
  
  GPtr->numDirectedEdges = outVertexListSize;
  
  GPtr->outVertexList = (ULONGINT_T *) node_malloc(outVertexListSize*sizeof(ULONGINT_T), TH);
  GPtr->paralEdgeIndex = (ULONGINT_T *) node_malloc(outVertexListSize*sizeof(ULONGINT_T), TH);
  outVertexListNum = 0;

  GPtr->outVertexList[0] = SDGdataPtr->endVertex[0];
  
  node_Barrier();
  
  /* Evaluate outVertexList */
  i0 = -1;
  pardo (i, 0, GPtr->numVertices, 1) {
    k = i;
    while ((i0 == -1) && (k != 0)) {
      for (j=0; j<numEdgesPlaced; j++) {
        if (k == SDGdataPtr->startVertex[j]) {
          i0 = j;
          break;
        }
      }
      k = k-1;
    }

    if ((i0 == -1) && (k == 0)) {
      i0 = 0;
    }
    
    for (j=i0; j<numEdgesPlaced; j++) {
      if (i == GPtr->numVertices-1) {
        break;
      }
      if (i != SDGdataPtr->startVertex[j]) {
        if ((j>0) && (i == SDGdataPtr->startVertex[j-1])) {
          if (j-i0 >= 1) {
            r = 0;
            GPtr->paralEdgeIndex[GPtr->outVertexIndex[i]+r] = i0;
            GPtr->outVertexList[GPtr->outVertexIndex[i]+r] = SDGdataPtr->endVertex[i0];
            r = r + 1;
            for (t=i0+1; t<j; t++) {
              if (SDGdataPtr->endVertex[t] != SDGdataPtr->endVertex[t-1]) {
                GPtr->paralEdgeIndex[GPtr->outVertexIndex[i]+r] = t;
                GPtr->outVertexList[GPtr->outVertexIndex[i]+r] = SDGdataPtr->endVertex[t];
                r = r + 1;
              }
            }

          }
        }
        i0 = j;
        break;
      }
    }
    
    
    if (i == GPtr->numVertices-1) {
      r = 0;
      if (numEdgesPlaced-i0 >= 0) {
        GPtr->paralEdgeIndex[GPtr->outVertexIndex[i]+r] = i0;
        GPtr->outVertexList[GPtr->outVertexIndex[i]+r] = SDGdataPtr->endVertex[i0];
        r = r + 1;
        for (t=i0+1; t<numEdgesPlaced; t++) {
          if (SDGdataPtr->endVertex[t] != SDGdataPtr->endVertex[t-1]) {
            GPtr->paralEdgeIndex[GPtr->outVertexIndex[i]+r] = t;
            GPtr->outVertexList[GPtr->outVertexIndex[i]+r] = SDGdataPtr->endVertex[t];
            r = r + 1;
          }
        }

      }

    }
    
  }

  
  node_Barrier();
  
  node_free(SDGdataPtr->startVertex, TH);
  node_free(SDGdataPtr->endVertex, TH);
    

  GPtr->inDegree  = (SHORTINT_T *) node_malloc(GPtr->numVertices*sizeof(SHORTINT_T), TH);
  GPtr->inVertexIndex = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);

  pardo(i, 0, GPtr->numVertices, 1) {
    GPtr->inDegree[i] = 0;
    GPtr->inVertexIndex[i] = 0;
  }
  
  /* A temp. array to store the inplied edges */
  impliedEdgeList = (ULONGINT_T *) node_malloc(GPtr->numVertices*MAX_CLUSTER_SIZE*sizeof(ULONGINT_T), TH);
  
  pardo(i, 0, GPtr->numVertices*MAX_CLUSTER_SIZE, 1)
    impliedEdgeList[i] = 0;
  
  /* An auxiliary array to store implied edges, in case we overshoot MAX_CLUSTER_SIZE */
  auxArr = (ULONGINT_T **) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T*), TH);
 
 
  vLock = (LOCK_T *) node_malloc(GPtr->numVertices*sizeof(LOCK_T),TH);


  lock_init_array(vLock, GPtr->numVertices, TH);
  
  node_Barrier();
  

  pardo (i, 0, GPtr->numVertices, 1) {
    /* Inspect adjacency list of vertex i */
    for (j=GPtr->outVertexIndex[i]; j<GPtr->outVertexIndex[i]+GPtr->outDegree[i]; j++) {
      v = GPtr->outVertexList[j];
      for (k = GPtr->outVertexIndex[v]; k<GPtr->outVertexIndex[v]+GPtr->outDegree[v]; k++) {
        if (GPtr->outVertexList[k] == i) {
          break;
        }
      }
      if (k == GPtr->outVertexIndex[v]+GPtr->outDegree[v]) {
        /* Add i to the impliedEdgeList of v */
        lock_it(vLock[v]);
        if (GPtr->inDegree[v] < MAX_CLUSTER_SIZE) {
          impliedEdgeList[v*MAX_CLUSTER_SIZE + GPtr->inDegree[v]] = i;
          GPtr->inDegree[v] = GPtr->inDegree[v] + 1;
        }
        else {
          /* Use auxiliary array to store the implied edge */
          /* Create an array if it's not present already */
          if ((GPtr->inDegree[v] % MAX_CLUSTER_SIZE) == 0) {
            auxArr[v] = (ULONGINT_T *) malloc(MAX_CLUSTER_SIZE*sizeof(ULONGINT_T));
          }
            auxArr[v][GPtr->inDegree[v] % MAX_CLUSTER_SIZE] = i;
            GPtr->inDegree[v] = GPtr->inDegree[v] + 1;
        } 
        unlock_it(vLock[v]);
        
      }
    }   
    
  }

  node_Barrier();
  
  lock_destroy_array(vLock, (int) GPtr->numVertices, TH);
  
  node_Barrier();
  
  prefix_sums(GPtr->inVertexIndex, GPtr->inDegree, GPtr->numVertices, TH);
  
  GPtr->numUndirectedEdges = GPtr->inVertexIndex[GPtr->numVertices-1] + GPtr->inDegree[GPtr->numVertices-1];

  node_Barrier();

  /* Create the inVertex List */
  GPtr->inVertexList = (ULONGINT_T *) node_malloc(GPtr->numUndirectedEdges*sizeof(ULONGINT_T), TH); 
  
  pardo(i, 0, GPtr->numVertices, 1) {
    for (j=GPtr->inVertexIndex[i]; j<GPtr->inVertexIndex[i]+GPtr->inDegree[i]; j++) {
      if (j-GPtr->inVertexIndex[i] < MAX_CLUSTER_SIZE) {
        GPtr->inVertexList[j] = impliedEdgeList[i*MAX_CLUSTER_SIZE+j-GPtr->inVertexIndex[i]];
      } else {
        GPtr->inVertexList[j] = auxArr[i][(j-GPtr->inVertexIndex[i]) % MAX_CLUSTER_SIZE];
      }

    }
  }
  
  node_Barrier();
  
  node_free(impliedEdgeList, TH);

  pardo(i, 0, GPtr->numVertices, 1) {
    if (GPtr->inDegree[i] > MAX_CLUSTER_SIZE) {
      free(auxArr[i]);
    }
  }

  node_Barrier();
  
  node_free(auxArr, TH);

}

void prefix_sums(ULONGINT_T *result, SHORTINT_T *input, ULONGINT_T arraySize,THREADED) {
  ULONGINT_T i, j, r, start, end, add_value;
  ULONGINT_T *p;

  r = arraySize / THREADS;

  p = (ULONGINT_T *) node_malloc(NOSHARE(THREADS)*sizeof(ULONGINT_T),TH);

  start =  MYTHREAD*r + 1;
  end   = (MYTHREAD+1)*r;

  if (MYTHREAD == THREADS-1)
    end = arraySize;

  for (j=start ; j<end ; j++)
    result[j] = input[j-1] + result[j-1];

  p[NOSHARE(MYTHREAD)] = result[end-1];

  node_Barrier();

  on_one_thread {
    for (j=1 ; j<THREADS ; j++)
      p[NOSHARE(j)] += p[NOSHARE(j-1)];
  }

  node_Barrier();

  if (MYTHREAD>0) {
    add_value=p[NOSHARE(MYTHREAD-1)];

    for (j=start-1 ; j<end ; j++)
      result[j] += add_value;
  }

  node_Barrier();
  node_free(p, TH);

}
