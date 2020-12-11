#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "simple.h"
#include "alg_radix_smp.h"
#include "globals.h"
#include "cutClusters.h"
 
void cutClusters(graph* GPtr, THREADED)
{
  ULONGINT_T i, j, k, r, t, v;

  /*
  Global array to keep track of vertex status:
  -1 if a vertex hasn't been assigned to a cluster yet
  t if it belongs to a cluster; t = iteration*THREADS + MYTHREAD
  */
  int* vStatus;
  
  int iter;
  
  ULONGINT_T currIndex, cliqueSize, verticesVisited, clusterCounter, cutSetCounter, cutSetIndex, cutSetIndexPrev;

  /* Data struct. for storing edgeCut */
  edge *pCutSet, *cutSet;
  
  ULONGINT_T* Index;
  ULONGINT_T* neighbourArray;
  ULONGINT_T* IndexSorted;
  ULONGINT_T* neighbourArraySorted;
  
  /* Temp vars for merging edge lists */
  ULONGINT_T *edgeStartCounter, *edgeEndCounter; 
  
  ULONGINT_T* startV;
  ULONGINT_T* clusterSize;
  
  FILE *outfp1, *outfp2;
  
  /* Sort the vertex list by their degree */
  Index = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  neighbourArray = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  IndexSorted = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  neighbourArraySorted = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);

  pardo(i, 0, GPtr->numVertices, 1) {
    neighbourArray[i] = GPtr->inDegree[i] + GPtr->outDegree[i];
    Index[i] = i;
  }
  
  node_Barrier();
  
  all_radixsort_node_aux_s3(GPtr->numVertices, neighbourArray, neighbourArraySorted, Index, IndexSorted, TH);

  node_Barrier();
 
  node_free(Index, TH);
  node_free(neighbourArray, TH);
  
  node_Barrier();
 
  /* Initialize vStatus Array */
  vStatus = (int *) node_malloc(GPtr->numVertices*sizeof(int), TH);
  
  node_Barrier();
  
  pardo(i, 0, GPtr->numVertices, 1)
    vStatus[i] = -1;
  
  node_Barrier();
  
  /* Allocate mem. for the cut set list */
  /* Maintain local arrays initially and merge them in the end */
  if (SCALE < 12)
    pCutSet = (edge *) malloc((1*(GPtr->numDirectedEdges)/THREADS)*sizeof(edge));
  else
    pCutSet = (edge *) malloc((0.2*(GPtr->numDirectedEdges)/THREADS)*sizeof(edge));
    
  /* Vertex to start from, on each thread */
  startV = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  clusterSize = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  
  verticesVisited = 0;
  
  node_Barrier();
  
  outfp1 = fopen("clusters.txt", "w");
  on_one_thread 
    fprintf(outfp1, "\nKernel 4 - Extracted Clusters\n");

  node_Barrier();
      
  iter = 0;
  currIndex = 0;      
  cutSetIndex = 0;

  while (verticesVisited < GPtr->numVertices) {
    
    /* Clear start vertex array */
    on_thread(MYTHREAD) {
       startV[MYTHREAD] = -1;
       clusterSize[MYTHREAD] = 0;
    }    
    

    if (currIndex == GPtr->numVertices) 
      currIndex = 0; 

    node_Barrier();

    /* Choose vertices to start from */
    /* Done sequentially right now, can be parallelized */
    on_one_thread {
      for (t=0; t<THREADS; t++) {
        for (r=currIndex; r<GPtr->numVertices; r++) {
          if (vStatus[IndexSorted[GPtr->numVertices - r - 1]] == -1) {
            startV[t] = IndexSorted[GPtr->numVertices - r - 1];
            vStatus[startV[t]] = iter*THREADS+t;
            for (j=0; j<GPtr->outDegree[startV[t]]; j++)
              if (vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[startV[t]]]] == -1) {
                vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[startV[t]]]] = iter*THREADS+t;
                clusterSize[t] = clusterSize[t] + 1;
              }
            for (j=0; j<GPtr->inDegree[startV[t]]; j++)
              if (vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[startV[t]]]] == -1) {
                vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[startV[t]]]] = iter*THREADS+t;
                clusterSize[t] = clusterSize[t] + 1;
              }
            currIndex = r+1;
            break;
          }
        }
      }
    }
    
    node_Barrier();
    
    /* Determine clusters and cut sets in parallel */
    on_thread(MYTHREAD) {
      
      i = startV[MYTHREAD];

      cliqueSize = 0;
      
      /* If the thread has some vertex to start from */
      if (i != -1)  {
         
	      cliqueSize = 1;
        /* clusterSize[MYTHREAD] gives the no. of 'unassigned' vertices adjacent to the current vertex */
        if ((clusterSize[MYTHREAD] >= 0.6*(GPtr->inDegree[i]+GPtr->outDegree[i])) || ((iter > (GPtr->numVertices)/(THREADS*MAX_CLUSTER_SIZE)) && (clusterSize[MYTHREAD]>0))) {
          
          /*Most of the adjacent vertices are unassigned, should be able to extract a cluster easily */
          

          /*Inspect adjacency list */
          for (j=0; j<GPtr->outDegree[i]; j++) {

            clusterCounter = 0;
            cutSetIndexPrev = cutSetIndex;
            cutSetCounter = 0;

            if (vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[i]]] == iter*THREADS+MYTHREAD) {

              v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
              /* Inspect vertices adjacent to v and determine if it belongs to a cluster or not */
              for (k=0; k<GPtr->outDegree[v]; k++) {
                if (vStatus[GPtr->outVertexList[k+GPtr->outVertexIndex[v]]] == iter*THREADS+MYTHREAD) {
                  clusterCounter = clusterCounter + 1;
                } else {
                  cutSetCounter = cutSetCounter + 1;
                  if (vStatus[GPtr->outVertexList[k+GPtr->outVertexIndex[v]]] == -1) {
                    /* To ensure that an edge is not added twice to the list */
                    pCutSet[cutSetIndex].startVertex = v;
                    pCutSet[cutSetIndex].endVertex = GPtr->outVertexList[k+GPtr->outVertexIndex[v]];
                    cutSetIndex = cutSetIndex + 1;
                  }
                }
              }

             if ( ((cutSetCounter >= clusterCounter) ) || ((SCALE<9) &&
	      (clusterCounter <= 2) && (GPtr->inDegree[v]+GPtr->outDegree[v] > clusterCounter +
	      cutSetCounter)  && (clusterSize[MYTHREAD] > clusterCounter + 2)) || ((SCALE>9) &&
	      (clusterCounter < 0.5*clusterSize[MYTHREAD]) )){
	      
                /* printf("[%d %d, %d %d %d, %d] ",i, v, clusterCounter, cutSetCounter, clusterSize[MYTHREAD],
		MYTHREAD); */
                /* v doesn't belong to this clique, free it */
                vStatus[v] = -1;
                
                /* Also add this edge to cutset list, removing previously added edges */
                cutSetIndex = cutSetIndexPrev;
                pCutSet[cutSetIndex].startVertex = i;
                pCutSet[cutSetIndex].endVertex = v;
                cutSetIndex = cutSetIndex + 1;

              } else {

                cliqueSize = cliqueSize + 1;
                 /* Add edges in inVertexList also to cut Set */
                for (k=0; k<GPtr->inDegree[v]; k++) {
                  if (vStatus[GPtr->inVertexList[k+GPtr->inVertexIndex[v]]] == -1) {
                    pCutSet[cutSetIndex].startVertex = v;
                    pCutSet[cutSetIndex].endVertex = GPtr->inVertexList[k+GPtr->inVertexIndex[v]];
                    cutSetIndex = cutSetIndex + 1;
                  }
                }

              }

            }

          }

          /* Do the same for the implied edges too */
          for (j=0; j<GPtr->inDegree[i]; j++) {

            clusterCounter = 0;
            cutSetIndexPrev = cutSetIndex;
            cutSetCounter = 0;

            if (vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[i]]] == iter*THREADS+MYTHREAD) {

              v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]];
              /* Inspect vertices adjacent to v and determine if it belongs to a cluster or not */
              for (k=0; k<GPtr->outDegree[v]; k++) {
                if (vStatus[GPtr->outVertexList[k+GPtr->outVertexIndex[v]]] == iter*THREADS+MYTHREAD) {
                  clusterCounter = clusterCounter + 1;
                } else {
                  cutSetCounter = cutSetCounter + 1;
                  if (vStatus[GPtr->outVertexList[k+GPtr->outVertexIndex[v]]] == -1) {
                    /* To ensure that an edge is not added twice to the list */
                    pCutSet[cutSetIndex].startVertex = v;
                    pCutSet[cutSetIndex].endVertex = GPtr->outVertexList[k+GPtr->outVertexIndex[v]];
                    cutSetIndex = cutSetIndex + 1;
                  }
                }
              }

              if ( ((cutSetCounter >= clusterCounter) ) || ((SCALE<9) &&
	      (clusterCounter <= 2) && (GPtr->inDegree[v]+GPtr->outDegree[v] > clusterCounter +
	      cutSetCounter)  && (clusterSize[MYTHREAD] > clusterCounter + 2)) || ((SCALE>9) &&
	      (clusterCounter < 0.5*clusterSize[MYTHREAD]) )){
                /* v doesn't belong to this clique, free it */
                /* printf("[%d %d, %d %d %d %d] ",i, v, clusterCounter, cutSetCounter, clusterSize[MYTHREAD],
		MYTHREAD); */
		vStatus[v] = -1;
                cutSetIndex = cutSetIndexPrev;
                pCutSet[cutSetIndex].startVertex = i;
                pCutSet[cutSetIndex].endVertex = v;
                cutSetIndex = cutSetIndex + 1;

              } else {

                cliqueSize = cliqueSize + 1;
                /* Add edges in inVertexList also to cut Set */
                for (k=0; k<GPtr->inDegree[v]; k++) {
                  if (vStatus[GPtr->inVertexList[k+GPtr->inVertexIndex[v]]] == -1) {
                    pCutSet[cutSetIndex].startVertex = v;
                    pCutSet[cutSetIndex].endVertex = GPtr->inVertexList[k+GPtr->inVertexIndex[v]];
                    cutSetIndex = cutSetIndex + 1;
                  }
                }

              }

            }

          }

        }

        if (clusterSize[MYTHREAD] == 0)  {
          /* Only one vertex in cluster */
          cliqueSize = 1;
	  
        } else {
          if ((clusterSize[MYTHREAD] < 0.6*(GPtr->inDegree[i]+GPtr->outDegree[i])) && (iter <= GPtr->numVertices/(THREADS*MAX_CLUSTER_SIZE))) {
            /* High perc. of intra-clique edges, do not commit clique */
            cliqueSize = 0;
            vStatus[i] = -1;
            
            for (j=0; j<GPtr->outDegree[i]; j++) {
              if (vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[i]]] == iter*THREADS+MYTHREAD)
                vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[i]]] = -1;
            }

            for (j=0; j<GPtr->inDegree[i]; j++) {
              if (vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[i]]] == iter*THREADS+MYTHREAD)
                vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[i]]] = -1;
            }
          }
	  
        }
      }
    }

    node_Barrier();

    /* Print to results.clq file */
    
    on_one_thread {
      for (t=0; t<THREADS; t++) {
        if (startV[t] != -1) {
          if (vStatus[startV[t]] == iter*THREADS+t) {
	  
            fprintf(outfp1, "%d ", startV[t]);
            for (j=0; j<GPtr->outDegree[startV[t]]; j++) {
              if (vStatus[GPtr->outVertexList[j+GPtr->outVertexIndex[startV[t]]]] == iter*THREADS+t)
        	    fprintf(outfp1, "%d ", GPtr->outVertexList[j+GPtr->outVertexIndex[startV[t]]]);
            }

            for (j=0; j<GPtr->inDegree[startV[t]]; j++) {
              if (vStatus[GPtr->inVertexList[j+GPtr->inVertexIndex[startV[t]]]] == iter*THREADS+t)
        	    fprintf(outfp1, "%d ", GPtr->inVertexList[j+GPtr->inVertexIndex[startV[t]]]);
            }

            fprintf(outfp1, "\n");
          }
	      }
      }
    }

    node_Barrier();
   
    
    iter = iter + 1;
    
    verticesVisited = verticesVisited + node_Reduce_d(cliqueSize, SUM, TH);
    
    node_Barrier();  
    
    if ((verticesVisited >= 0.95*GPtr->numVertices) || (iter>GPtr->numVertices/2))
       break;
    
  }
  
  node_Barrier();
  
  /* Take care of unmarked vertices */
  on_one_thread {
    if (verticesVisited < GPtr->numVertices) { 
      for(i=0; i<GPtr->numVertices; i++) {
        if (vStatus[i] == -1) {
          vStatus[i] = iter*THREADS+MYTHREAD;
	        fprintf(outfp1, "%d\n", i);
	        iter = iter + 1;
        }
      }
    }
  }
  
  node_Barrier();
       
  /* Merge partial Cutset Lists */
  edgeStartCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  edgeEndCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
   
  on_thread(MYTHREAD) {  
    edgeEndCounter[MYTHREAD] = cutSetIndex;     
    edgeStartCounter[MYTHREAD] = 0;
  }

  node_Barrier();

  on_one_thread {
    for (t=1; t<THREADS; t++) {
      edgeEndCounter[t] = edgeEndCounter[t-1] + edgeEndCounter[t];
      edgeStartCounter[t] = edgeEndCounter[t-1]; 
    }
  }  
  
  node_Barrier();
	    
  cutSetCounter = node_Reduce_d(cutSetIndex, SUM, TH);

  cutSet = (edge *) node_malloc(cutSetCounter*sizeof(edge), TH); 

  on_thread(MYTHREAD) {
    for (j = edgeStartCounter[MYTHREAD]; j<edgeEndCounter[MYTHREAD]; j++) {
      cutSet[j].startVertex = pCutSet[j-edgeStartCounter[MYTHREAD]].startVertex;
      cutSet[j].endVertex = pCutSet[j-edgeStartCounter[MYTHREAD]].endVertex;
    }   
  } 
  
  node_Barrier();
  
  outfp2 = fopen("edgeCut.txt", "w");
  
  on_one_thread {
    /* printf("Cut Set Counter  - %d\n", cutSetCounter);*/
    fprintf(outfp2, "\nEdges in Cut Set - \n");
    for (i=0; i<cutSetCounter; i++)
      fprintf(outfp2, "[%d %d] ", cutSet[i].startVertex, cutSet[i].endVertex);     
  }
  
  node_Barrier();
  
 
  fclose(outfp1);
  fclose(outfp2);
  node_free(edgeStartCounter, TH);
  node_free(edgeEndCounter, TH);
  free(pCutSet);
  node_free(IndexSorted, TH);      
  node_free(neighbourArraySorted, TH);
  node_free(startV, TH);
  node_free(clusterSize, TH);
  node_free(cutSet, TH);
  node_free(vStatus, TH);
}
