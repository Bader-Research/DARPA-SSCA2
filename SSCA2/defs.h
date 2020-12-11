#ifndef _DEFS_H
#define _DEFS_H

/* #define WEIGHT_T unsigned long int */
#define ULONGINT_T unsigned long int
#define LONGINT_T long int
#define SHORTINT_T short int

typedef struct 
{
  ULONGINT_T* startVertex;
  ULONGINT_T* endVertex;

  LONGINT_T* intWeight;
  /* The idea is to store the index of the string weights (as a negative value)
  in the int Weight array. A negative value because we need to sort on
  the intWeights in Kernel 2. Hence the long int
  */

  char* strWeight;
  ULONGINT_T numEdgesPlaced;

} graphSDG;

typedef struct /*the graph data structure*/
{ 
  ULONGINT_T numVertices;
  ULONGINT_T numEdges;
  
  ULONGINT_T numDirectedEdges;
  ULONGINT_T numUndirectedEdges;

  ULONGINT_T numIntEdges;
  ULONGINT_T numStrEdges;

  SHORTINT_T* outDegree;
  ULONGINT_T* outVertexIndex;
  ULONGINT_T* outVertexList;
  ULONGINT_T* paralEdgeIndex;

  SHORTINT_T* inDegree;
  ULONGINT_T* inVertexIndex;
  ULONGINT_T* inVertexList;

  LONGINT_T* intWeight;
  char* strWeight;

} graph;


typedef struct /*edge structure for Kernel 2*/
{
  ULONGINT_T startVertex;
  ULONGINT_T endVertex;
  ULONGINT_T edgeNum;
} edge;

typedef struct /*Vertex list returned by Kernel 3*/
{
  ULONGINT_T num;
  int depth;  
} V;

typedef struct l /*A Linked list for Kernel 3*/
{
  ULONGINT_T num;
  SHORTINT_T depth;
  struct l* next;
} Vl;

typedef struct /*A dynamic array for Kernel 3*/
{
  ULONGINT_T numArrays;
  ULONGINT_T* arraySize;
  V** vList;
} Vd;

#endif

