#ifndef _LIST_H_
#define _LIST_H_

#include "Graph.h"

class ListGraph : public Graph{	
private:
	map < int, int >* m_List;
	vector<int>* kw_graph;
	
public:	
	ListGraph(bool type, int size);
	~ListGraph();

	// Getter for kw_graph
    vector<int>* getKWGraph() const {
        return kw_graph;
    }
		
	void getAdjacentEdges(int vertex, map<int, int>* m);	
	void getAdjacentEdgesDirect(int vertex, map<int, int>* m);
	void insertEdge(int from, int to, int weight);	
	bool printGraph(ofstream *fout);

};

#endif