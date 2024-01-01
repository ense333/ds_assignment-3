#include "ListGraph.h"
#include <iostream>
#include <utility>
#include <fstream>

ListGraph::ListGraph(bool type, int size) : Graph(type, size)
{
    m_List = new map<int, int>[size + 1];
}

ListGraph::~ListGraph()    
{
    delete[] m_List;
}

void ListGraph::getAdjacentEdges(int vertex, map<int, int>* m)    // Definition of getAdjacentEdges for directed graph
{
    // Copy adjacent edges from the list to the provided map.
    for (auto i = m_List[vertex].begin(); i != m_List[vertex].end(); i++) {
        m->insert({i->first, i->second});
    }
    return;
}

void ListGraph::getAdjacentEdgesDirect(int vertex, map<int, int>* m) {		//For Undirected graph
    m->clear();

    if (!m_List[vertex].empty()) {
        *m = m_List[vertex];
    }

    // Add incoming edges.
    for (int i = 1; i < m_Size + 1; ++i) {
        if (i != vertex && !m_List[i].empty()) {
            auto it = m_List[i].find(vertex);
            if (it != m_List[i].end()) {
                // Add incoming edges if they exist.
                (*m)[i] = it->second;
            }
        }
    }
}

void ListGraph::insertEdge(int from, int to, int weight) // Definition of insertEdge
{
    m_List[from][to] = weight;
}

bool ListGraph::printGraph(ofstream* fouts) // Definition of print Graph
{
    //check if graph is empty
    if (m_Size == 0) {
        return false;
    }
    *fouts << "=======PRINT========\n";
    for (int i = 1; i <= m_Size; ++i) {
        int count = 0;
        int totalCount = m_List[i].size();
        *fouts << "[" << i << "] ->";		//Print vertex value
        if (!m_List[i].empty()) {		//Then print vertex that has edge with above vertex
            for (const auto& edge : m_List[i]) {
                count++;
                if (count == totalCount) {
                    *fouts << " (" << edge.first << "," << edge.second << ")";
                } else {
                    *fouts << " (" << edge.first << "," << edge.second << ") ->";
                }
            }
        }
        *fouts << endl;
    }
    *fouts << "===================\n\n";
    return true;
}
