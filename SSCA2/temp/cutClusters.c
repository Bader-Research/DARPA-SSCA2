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

  char* vStatusArr;
  /*
    Char array to keep track of vertex status:
    u if unmarked
    a if vertex is adjacent to the current vertex being inspected
    p if vertex is inspected and probably in a cluster 
    t if vertex is taken 
  */
  double dtime;
  int iterCount;
  
  ULONGINT_T verticesVisited, clusterCounter, cutSetCounter, pCutSetCounter, adjCounter, pAdjCounter;

  edge *pCutSet, *cutSet;
  
  ULONGINT_T* Index;
  ULONGINT_T* neighbourArray;
  ULONGINT_T* IndexSorted;
  ULONGINT_T* neighbourArraySorted;
  
  ULONGINT_T *edgeStartCounter, *edgeEndCounter; 
  
  FILE* outfp;
  
  /* Sort the vertex list by their degree */
  Index = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  neighbourArray = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  IndexSorted = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);
  neighbourArraySorted = (ULONGINT_T *) node_malloc(GPtr->numVertices*sizeof(ULONGINT_T), TH);

  pardo(i, 0, GPtr->numVertices, 1) {
    neighbourArray[i] = GPtr->inDegree[i] + GPtr->outDegree[i];
    /* printf("%d ", neighbourArray[i]); */
    Index[i] = i;
  }
  
  node_Barrier();
  
  all_radixsort_node_aux_s3(GPtr->numVertices, neighbourArray, neighbourArraySorted, Index, IndexSorted, TH);

  node_Barrier();
 
  node_free(Index, TH);
  node_free(neighbourArray, TH);
  
  node_Barrier();
 
    /* Initialize vStatus Array */
  vStatusArr = (char *) node_malloc(GPtr->numVertices*sizeof(char), TH);
  
  node_Barrier();
  
  pardo(i, 0, GPtr->numVertices, 1)
      vStatusArr[i] = 'u';
  
  node_Barrier();
  
  cutSet = (edge *) node_malloc(2*GPtr->numUndirectedEdges*sizeof(edge), TH);
  
  edgeStartCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  edgeEndCounter = (ULONGINT_T *) node_malloc(THREADS*sizeof(ULONGINT_T), TH);
  
  node_Barrier();
  
  /* Allocate mem. for the temp lists */
  pCutSet = (edge *) malloc(3*MAX_CLUSTER_SIZE*sizeof(edge));
  
  verticesVisited = 0;
  
  node_Barrier();
  
  /* For the vertices with no neighbours, write to the results file and mark them off */
  outfp = fopen("results.clq", "a");
  on_one_thread {
    fprintf(outfp, "\n\n");
    fprintf(outfp, "Results of graph clustering kernel\n");
    
    for (i=0; i<=TOT_VERTICES; i++) {
      if (neighbourArraySorted[i] == 0) {
	      vStatusArr[IndexSorted[i]] = 't';
	      fprintf(outfp, "%d\n", IndexSorted[i]);
	      verticesVisited = verticesVisited + 1;
      }
      else {
	      break;
      }
    }
  }
  
  node_Barrier();
  
  
  verticesVisited = node_Bcast_i(verticesVisited, TH);
  
  /*
  on_one_thread 
    printf("Vertices visited initially - %d \n", verticesVisited);
  */

  
  iterCount = 0;
  
  
  while (verticesVisited < GPtr->numVertices) {
    
    iterCount = iterCount + 1;
    
    if (iterCount <= 2) {
    /* Most of the vertices should be marked in two iterations */

      cutSetCounter = 0;
      
      for (r=0; r<GPtr->numVertices; r++) {
      /* Start from vertex with max. degree */
        
	      i = IndexSorted[GPtr->numVertices - r - 1];
        adjCounter = 0;
		

        if (vStatusArr[i] == 'u') {
          
          node_Barrier();
          /*on_one_thread
	           vStatusArr[i] = 'a';*/
	  
      	  /* Determine eff. no. of vertices adjacent to this vertex */
	        pardo (j, 0, GPtr->inDegree[i]+GPtr->outDegree[i], 1) {
            
	          if (j<GPtr->outDegree[i]) {
	            /* Inspect outVertexList */
	            v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
	            if (vStatusArr[v] == 'u') {
		            adjCounter = adjCounter + 1;
		            vStatusArr[v] = 'a';
	            }
            }

	          else {
	            /* Inspect in vertex list */
	            v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]-GPtr->outDegree[i]];
	            if (vStatusArr[v] == 'u') {
		            adjCounter = adjCounter + 1;
		            vStatusArr[v] = 'a';
	            }
            }

	        }

	        node_Barrier();

	        adjCounter = node_Reduce_d(adjCounter, SUM, TH);
          
	        if (adjCounter > MAX_CLUSTER_SIZE/2) {
            
	          node_Barrier();
	
	          pardo(j, 0, GPtr->inDegree[i]+GPtr->outDegree[i], 1) {

	            pCutSetCounter = 0;
	            pAdjCounter = 0;

	            if (j<GPtr->outDegree[i]) {
		            /* Inspect outVertexList */
		            v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
		            if (vStatusArr[v] == 'a') {
		              for (t=0; t < GPtr->inDegree[v]+GPtr->outDegree[v]; t++) {
		                if (t<GPtr->outDegree[v]) {
		                  k = GPtr->outVertexList[t+GPtr->outVertexIndex[v]];
		                  if (vStatusArr[k] == 'a' || vStatusArr[k] == 'p')
			                  pAdjCounter = pAdjCounter + 1;
		                  else
			                  pCutSetCounter = pCutSetCounter + 1;

		                }
		                else {
		                  k = GPtr->inVertexList[t+GPtr->inVertexIndex[v]-GPtr->outDegree[v]];
		                  if (vStatusArr[k] == 'a' || vStatusArr[k] == 'p')
			                  pAdjCounter = pAdjCounter + 1;
		                  else
			                  pCutSetCounter = pCutSetCounter + 1;

		                }

		              }

		              if ((pAdjCounter >= 0.75 * adjCounter) && (pAdjCounter > pCutSetCounter))
		                vStatusArr[v] = 'p';

		            }

	            }

	            else {
		            /* Inspect in vertex list */
		            v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]-GPtr->outDegree[i]];
        	      if (vStatusArr[v] == 'a') {
		              for (t=0; t < GPtr->inDegree[v]+GPtr->outDegree[v]; t++) {
		                if (t<GPtr->outDegree[v]) {
		                  k = GPtr->outVertexList[t+GPtr->outVertexIndex[v]];
		                  if (vStatusArr[k] == 'a' || vStatusArr[k] == 'p')
			                  pAdjCounter = pAdjCounter + 1;
		                  else
			                  pCutSetCounter = pCutSetCounter + 1;

		                }
		                else {
		                  k = GPtr->inVertexList[t+GPtr->inVertexIndex[v]-GPtr->outDegree[v]];
		                  if (vStatusArr[k] == 'a' || vStatusArr[k] == 'p')
			                  pAdjCounter = pAdjCounter + 1;
		                  else
			                  pCutSetCounter = pCutSetCounter + 1;
                    }

		               }
		               if ((pAdjCounter >= 0.75 * adjCounter) && (pAdjCounter > pCutSetCounter))
		                  vStatusArr[v] = 'p';

		             }

	            }

	          }

	       }

    	  node_Barrier();

    	  clusterCounter = 0;

    	  pardo(j, 0, GPtr->inDegree[i]+GPtr->outDegree[i], 1) {

    	    if (j<GPtr->outDegree[i]) {
    	      /* Inspect outVertexList */
    	      v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
    	      if (vStatusArr[v] == 'p') {
    		      clusterCounter = clusterCounter + 1;
    	      }

    	    }

    	    else {
    	      /* Inspect in vertex list */
            v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]-GPtr->outDegree[i]];
    	      if (vStatusArr[v] == 'p') {
    		      clusterCounter = clusterCounter + 1;
    	      }

    	    }

    	  }

    	  node_Barrier();

    	  clusterCounter = node_Reduce_d(clusterCounter, SUM, TH);


    	  if (clusterCounter >= 0.75*adjCounter) {

           /* Accept the cluster */
          verticesVisited = verticesVisited + clusterCounter;
    	    pCutSetCounter = 0;
    	    vStatusArr[i] = 't';

    	    on_one_thread 
              fprintf(outfp, "\n%d ", i);
	      
    	    node_Barrier();

    	    pardo (j, 0, GPtr->inDegree[i]+GPtr->outDegree[i], 1) {
    	      if (j<GPtr->outDegree[i]) {
    		      /* Inspect outVertexList */
    		      v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
    		      if (vStatusArr[v] == 'p') {
    		        vStatusArr[v] = 't';
    		        fprintf(outfp, "%d ", v);
    		      } else {
    	      	        vStatusArr[v] = 'u';
    		        pCutSet[pCutSetCounter].startVertex = i;
    		        pCutSet[pCutSetCounter].endVertex = v;
    		        pCutSetCounter = pCutSetCounter + 1;
    		      }


    	      }

    	      else {
    		      /* Inspect inVertex list */
            	v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]-GPtr->outDegree[i]];
    		      if (vStatusArr[v] == 'p') {
    	              vStatusArr[v] = 't';
    		      fprintf(outfp, "%d ", v);

    		      } else {
    	      	  vStatusArr[v] = 'u';
    		        pCutSet[pCutSetCounter].startVertex = i;
    		        pCutSet[pCutSetCounter].endVertex = v;
    		        pCutSetCounter = pCutSetCounter + 1;

    		      }

	         }

	       }
	    
	    node_Barrier();
	    
	    /* Merge partial Cut Set Lists */

	    on_thread(MYTHREAD) {
              edgeEndCounter[MYTHREAD] = pCutSetCounter;     
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
	    
            cutSetCounter = cutSetCounter + node_Reduce_d(pCutSetCounter, SUM, TH);
            
	    node_Barrier(); 
	     
            on_thread(MYTHREAD) {
              for (j = edgeStartCounter[MYTHREAD]; j<edgeEndCounter[MYTHREAD]; j++) {
        	cutSet[j+cutSetCounter].startVertex = pCutSet[j-edgeStartCounter[MYTHREAD]].startVertex;
        	cutSet[j+cutSetCounter].endVertex = pCutSet[j-edgeStartCounter[MYTHREAD]].endVertex;
              }   
            } 

	    
	  }   

	  else {
            
	    /* Reject cluster */
	    vStatusArr[i] = 't';
	    
	    node_Barrier();
	    
	    pardo (j, 0, GPtr->inDegree[i]+GPtr->outDegree[i], 1) {
	      if (j<GPtr->outDegree[i]) {
		/* Inspect outVertexList */
		v = GPtr->outVertexList[j+GPtr->outVertexIndex[i]];
		vStatusArr[v] == 'u';

	      }  

	      else {
		/* Inspect in vertex list */
        	v = GPtr->inVertexList[j+GPtr->inVertexIndex[i]-GPtr->outDegree[i]];
		vStatusArr[v] = 'u';

	      }

	    }
	    
	    
	  }
	  	  
	}
      }  	
    }
    
    if (iterCount > 2) 
      break;
       
  }
  
  node_Barrier();
  fclose(outfp);
  node_free(edgeStartCounter, TH);
  node_free(edgeEndCounter, TH);
  free(pCutSet);
  node_free(IndexSorted, TH);      
  node_free(neighbourArraySorted, TH);
  node_free(cutSet, TH);
  node_free(vStatusArr, TH);
}
