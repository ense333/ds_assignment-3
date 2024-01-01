#include "MatrixGraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

MatrixGraph::MatrixGraph(bool type, int size) : Graph(type, size)
{
    // Allocate memory for the two-dimensional array.
    m_Mat = new int*[size + 1];
    for (int i = 0; i < size + 1; i++) {
        m_Mat[i] = new int[size + 1]();
    }
}

MatrixGraph::~MatrixGraph()
{
    // Deallocate the memory for the two-dimensional array.
    for (int i = 0; i < m_Size + 1; i++) {
        delete[] m_Mat[i];
    }
    delete[] m_Mat;
}

void MatrixGraph::getAdjacentEdges(int vertex, map<int, int>* m)
{
    // Retrieve adjacent edges for the given vertex.
    for (int i = 1; i <= m_Size; ++i) {
        if (m_Mat[vertex][i] != 0) {
            (*m)[i] = m_Mat[vertex][i];
        }
    }
}

void MatrixGraph::getAdjacentEdgesDirect(int vertex, map<int, int>* m)
{
    // Clear the map before populating it.
    m->clear();

    // Retrieve adjacent edges for the given vertex.
    for (int i = 1; i <= m_Size; ++i) {
        if (m_Mat[vertex][i] != 0) {
            (*m)[i] = m_Mat[vertex][i];
        }
    }

    // Retrieve incoming edges for the given vertex in a directed graph.
    for (int i = 1; i <= m_Size; ++i) {
        if (m_Mat[i][vertex] != 0) {
            (*m)[i] = m_Mat[i][vertex];
        }
    }
}

void MatrixGraph::insertEdge(int from, int to, int weight)
{
    // Insert an edge with a specified weight.
    m_Mat[from][to] = weight;
}

bool MatrixGraph::printGraph(ofstream* fouts)
{
    //check if graph is empty
    if (m_Size == 0) {
        return false;
    }
    // Print the matrix representation of the graph.
    *fouts << "=======PRINT========\n";
    *fouts << "   ";
    for (int i = 1; i <= m_Size; ++i) {
        *fouts << "[" << i << "]";      //For index values
    }
    *fouts << endl;
    for (int i = 1; i <= m_Size; ++i) {
        *fouts << "[" << i << "]";      //For index values
        for (int j = 1; j <= m_Size; ++j) {
            *fouts << " " << m_Mat[i][j] << " ";        //Print each graph values
        }
        *fouts << endl;
    }
    *fouts << "===================\n\n";

    return true;
}
