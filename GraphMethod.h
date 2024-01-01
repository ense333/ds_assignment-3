#ifndef _GRAPHMETHOD_H_
#define _GRAPHMETHOD_H_

#include "ListGraph.h"
#include "MatrixGraph.h"

bool BFS(Graph* graph, char option, int vertex, ofstream* fouts);     
bool DFS(Graph* graph, char option,  int vertex, ofstream* fouts);     
bool KWANGWOON(Graph* graph, int vertex, ofstream* fouts);  
bool Kruskal(Graph* graph, ofstream* fouts);
bool Dijkstra(Graph* graph, char option, int vertex, ofstream* fouts);    //Dijkstra
bool Bellmanford(Graph* graph, char option, int s_vertex, int e_vertex, ofstream* fouts); //Bellman - Ford
bool FLOYD(Graph* graph, char option, ofstream* fouts);   //FLoyd

#endif
