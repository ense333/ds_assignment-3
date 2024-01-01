#include <iostream>
#include <vector>
#include "GraphMethod.h"
#include <stack>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <utility>
#include <limits>
#include <algorithm>
#include <fstream>

using namespace std;

const int INF = numeric_limits<int>::max();

// Definition of the Node structure with vertex and weight
struct Node {
    int vertex, weight;
    Node(int v, int w): vertex(v), weight(w) {}

    // Overloaded greater-than operator for comparison
    bool operator>(const Node& n) const {
        return weight > n.weight;
    }
};

// Function for insertion sort
void insertionSort(vector<pair<int, pair<int, int>>>& arr, int low, int high) {
    for (int i = low + 1; i <= high; i++) {
        auto key = arr[i];
        int j = i - 1;

        // Move elements of arr[0..i-1] that are greater than key,
        // to one position ahead of their current position
        while (j >= low && arr[j].first > key.first) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

// Function for the partition step in quicksort
int partition(vector<pair<int, pair<int, int>>>& arr, int low, int high) {
    auto pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j].first <= pivot.first) {
            i++;
            swap(arr[i], arr[j]);   //Swap values
        }
    }
    swap(arr[i + 1], arr[high]);       //Swap with pivot
    return (i + 1);
}


// Hybrid quicksort function with a specified segment size
void newQuickSort(vector<pair<int, pair<int, int>>>& arr, int low, int high, int segment_size) {
    if (low < high) {
        // If the size of the segment is less than or equal to the specified segment size,
        // perform insertion sort on the segment
        if (high - low + 1 <= segment_size) {
            insertionSort(arr, low, high);
        } else {
            // Otherwise, partition the array and recursively apply newQuickSort to the subarrays
            int pivot = partition(arr, low, high);
            newQuickSort(arr, low, pivot - 1, segment_size);
            newQuickSort(arr, pivot + 1, high, segment_size);
        }
    }
}

// Disjoint set data structure class
class DisjointSet {
private:
    vector<int> parent, rank;

public:
    // Constructor initializes the disjoint set with 'n' elements
    DisjointSet(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    // Find operation with path compression
    int find(int u) {
        if (u != parent[u])
            parent[u] = find(parent[u]); // Path compression
        return parent[u];
    }

    // Merge two sets based on their ranks
    void merge(int u, int v) {
        u = find(u);
        v = find(v);
        if (u == v) return;
        if (rank[u] > rank[v]) swap(u, v);
        parent[u] = v; // Always attach the smaller rank tree under the root of the higher rank tree
        if (rank[u] == rank[v]) rank[v]++;
    }
};

// Find the weight of the edge between vertices u and v
int findWeight(const vector<pair<int, pair<int, int>>>& edges, int u, int v) {
    for (const auto& edge : edges) {
        if ((edge.second.first == u && edge.second.second == v) || 
            (edge.second.first == v && edge.second.second == u)) {
            return edge.first; // Return the weight of the corresponding edge
        }
    }
    return -1; // Return -1 if the edge is not found
}


// Breadth-First Search (BFS) function on a graph
bool BFS(Graph* graph, char option, int vertex, ofstream* fouts) {
    // Check for valid option, graph, and vertex
    if (option != 'Y' && option != 'N' || graph == NULL || vertex < 1 || vertex > graph->getSize()) {       //For exception check
        return false;
    }

    // Array to keep track of visited vertices
    bool* visited = new bool[graph->getSize()];
    fill(visited, visited + graph->getSize(), false);

    int countA = 0;
    // Check if the vertex is within a valid range
    if (vertex > 0 && vertex <= graph->getSize()) {
        *fouts << "==========BFS========\n";

        if (option == 'Y') {  // Directed graph BFS
            *fouts << "Directed Graph BFS result" << endl;
            *fouts << "startvertex : " << vertex << endl;
            
            queue<int> q;   //Using queue
            q.push(vertex);     //push first vertex
            visited[vertex - 1] = true;     //Then set its value to true

            // Continue the loop until the queue is empty
            while (!q.empty()) {
                int v = q.front();  // Get the front element of the queue
                if (countA > 0) {
                    *fouts << " -> ";
                }
                *fouts << v;  // Print the current vertex
                countA++;
                q.pop();  // Remove the front element from the queue

                map<int, int> adjMap;
                graph->getAdjacentEdges(v, &adjMap);  // Get adjacent edges for the current vertex in a directed graph

                // Iterate through adjacent vertices
                for (const auto& adj : adjMap) {
                    int adjVertex = adj.first;
                    if (!visited[adjVertex - 1]) {
                        q.push(adjVertex);  // Add the adjacent vertex to the queue
                        visited[adjVertex - 1] = true;  // Mark the adjacent vertex as visited
                    }
                }
            }
            *fouts << endl;
            *fouts << "========================\n\n";
            delete[] visited;
            return true;
        } else if (option == 'N') {  // Undirected graph BFS
            *fouts << "Undirected Graph BFS result" << endl;
            *fouts << "startvertex : " << vertex << endl;

            queue<int> q;
            q.push(vertex);
            visited[vertex - 1] = true;

            // Perform breadth-first traversal using a while loop until the queue is empty
            while (!q.empty()) {
                // Retrieve the front vertex from the queue
                int v = q.front();

                // Output arrow separator (->) if it's not the first vertex in the path
                if (countA > 0) {
                    *fouts << " -> ";
                }
                // Output the current vertex to the file stream
                *fouts << v;
                // Dequeue the front vertex from the queue
                q.pop();
                // Increment the count of visited vertices
                countA++;
                // Create an empty adjacent map for the current vertex
                map<int, int> adjMap;
                graph->getAdjacentEdgesDirect(v, &adjMap);      //Since it is undirected graph

                for (const auto& adj : adjMap) {
                    int adjVertex = adj.first;
                    if (!visited[adjVertex - 1]) {
                        q.push(adjVertex);      //Add adjacent vertex
                        visited[adjVertex - 1] = true;      //Also change visted array
                    }
                }
            }
            *fouts << endl;
            *fouts << "==========================\n\n";
            delete[] visited;
            return true;
        } else {  // Invalid option input
            delete[] visited;
            return false;
        }
    } else {  // Given vertex is out of range
        delete[] visited;
        return false;
    }
}

// Depth-First Search (DFS) function on a graph with an option for directed or undirected output
bool DFS(Graph* graph, char option, int vertex, ofstream* fouts) {
    // Check for valid option, graph, and vertex
    if (option != 'Y' && option != 'N' || graph == NULL || vertex < 1 || vertex > graph->getSize()) {       //For exception check
        return false;
    }

    int countA = 0;
    // Array to keep track of visited vertices
    bool* visited = new bool[graph->getSize()];
    fill(visited, visited + graph->getSize(), false);

    *fouts << "=============DFS===========\n";
    if (option == 'Y') {        //For direction option
        *fouts << "Directed Graph DFS result" << endl;
    } else if (option == 'N') {
        *fouts << "Undirected Graph DFS result" << endl;
    }

    *fouts << "startvertex : " << vertex << endl;

    stack<int> s;       //Using stack 
    s.push(vertex);     //Put vertex value

    while (!s.empty()) {        //Loop until stack is empty
        int v = s.top();
        s.pop();

        if (!visited[v - 1]) {      //Check wether vertex is visited
            visited[v - 1] = true;      
            if (countA > 0) {
                *fouts << " -> ";
            }

            *fouts << v;
            countA++;
            map<int, int> adjMap;
            if (option == 'N') {
                graph->getAdjacentEdgesDirect(v, &adjMap);  // For undirected graphs
            } else {
                graph->getAdjacentEdges(v, &adjMap);  // For directed graphs
            }

            for (auto it = adjMap.rbegin(); it != adjMap.rend(); ++it) {
                // Retrieve the adjacent vertex and check if it has not been visited
                int adjVertex = it->first;
                if (!visited[adjVertex - 1]) {
                    // If the adjacent vertex is not visited, push it onto the stack
                    s.push(adjVertex);
                }
            }
        }
    }
    *fouts << endl;
    *fouts << "===========================\n\n";
    delete[] visited;
    return true;
}


// Kruskal's algorithm to find Minimum Spanning Tree (MST) in a graph
bool Kruskal(Graph* graph, ofstream* fouts) {
    // Check if the graph is valid
    if (graph == NULL) {
        return false;
    }

    const int V = graph->getSize();
    vector<pair<int, pair<int, int>>> edges;  // (weight, (u, v))
    vector<vector<pair<int, int>>> mst(V + 1);  // Vector to store edges connected to each vertex in MST

    // Store all edges in the edges vector
    for (int u = 1; u <= V; ++u) {
        map<int, int> adjMap;
        graph->getAdjacentEdgesDirect(u, &adjMap);
        for (const auto& edge : adjMap) {
            int v = edge.first, weight = edge.second;
            if (u < v) {
                edges.push_back({weight, {u, v}});
            }
        }
    }

    // Sort all edges based on weight
    newQuickSort(edges, 0, edges.size() - 1, 6);
    DisjointSet sets(V + 1);
    int total_weight = 0;
    int setCount = 0;
    // Iterate through sorted edges and add to MST
    for (const auto& e : edges) {
        int weight = e.first;
        int u = e.second.first, v = e.second.second;
        if (sets.find(u) != sets.find(v)) {
            sets.merge(u, v);
            total_weight += weight;
            mst[u].push_back({v, weight});  // Add edge connected to vertex u
            mst[v].push_back({u, weight});  // Add edge connected to vertex v (for undirected graph)
            if (++setCount == V - 1) break;
        }
    }

    // Check if MST is properly formed
    if (setCount != V - 1) {    //If not all vertices are connected
        return false;
    }

    // Print the result
    *fouts << "======Kruskal======" << endl;
    for (int u = 1; u <= V; ++u) {
        sort(mst[u].begin(), mst[u].end());  // Sort edges connected to each vertex
        *fouts << "[" << u << "] -> ";
        for (const auto& edge : mst[u]) {
            *fouts << " " << edge.first << "(" << edge.second << ")";
        }
        *fouts << endl;
    }
    *fouts << "cost: " << total_weight << endl;
    *fouts << "===================\n\n";

    return true;
}

bool Dijkstra(Graph* graph, char option, int vertex, ofstream* fouts) {
    // Check if the option, graph, and vertex are valid
    if (option != 'Y' && option != 'N' || graph == NULL || vertex < 1 || vertex > graph->getSize()) {   //For exception check
        return false;
    }

    const int V = graph->getSize();
    vector<int> dist(V + 1, INF);  // Initialize distance information to infinity
    vector<int> prev(V + 1, -1);   // Array to trace the path
    priority_queue<Node, vector<Node>, greater<Node>> pq;       //Using priority_queue

    dist[vertex] = 0;
    pq.push(Node(vertex, 0));

    while (!pq.empty()) {      
        Node node = pq.top();
        pq.pop();
        int u = node.vertex;
        map<int, int> adjMap;
        if (option == 'Y') {        //For direction option
            graph->getAdjacentEdges(u, &adjMap);
        } else if (option == 'N') {
            graph->getAdjacentEdgesDirect(u, &adjMap);  // Get adjacent edges for directed or undirected graph
        } else {        //Wrong input
            return false;
        }
        for (const auto& adj : adjMap) {
            int v = adj.first, weight = adj.second;

            // Check for negative weights
            if (weight < 0) {
                return false; // Return false if the graph has negative weights
            }

            // Actual distance value should be assigned to dist[u] instead of INF
            if (dist[u] != INF && dist[v] > dist[u] + weight) {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.push(Node(v, dist[v]));
            }
        }
    }
    // Print the result
    *fouts << "======Dijkstra======" << endl;
    if (option == 'Y') {        //For direction option
        *fouts << "Directed Graph Dijkstra result\n";
    } else if (option == 'N') {
        *fouts << "Undirected Graph Dijkstra result\n";
    }
    *fouts << "startvertex: " << vertex << endl;
    // Iterate through all vertices in the graph
    for (int v = 1; v <= V; ++v) {
        // Check if there is no path to vertex v from the starting vertex
        if (dist[v] == INF) {
            *fouts << "[" << v << "] " << vertex << " -> " << v << " X" << endl;
        } else {
            // Check if the vertex is not the starting vertex
            if (v != vertex) {
                vector<int> path;
                // Trace the path from vertex v back to the starting vertex using prev array
                for (int at = v; at != -1; at = prev[at]) {
                    path.push_back(at);
                }
                reverse(path.begin(), path.end());
                // Print the path if it has more than one vertex 
                if (path.size() > 1) {
                    *fouts << "[" << v << "] " << vertex;
                    for (size_t i = 1; i < path.size(); ++i) {
                        *fouts << " -> " << path[i];
                    }
                    *fouts << " (" << dist[v] << ")" << endl;
                }
            }
        }
    }
    *fouts << "====================\n\n";

    return true;
}

bool Bellmanford(Graph* graph, char option, int s_vertex, int e_vertex, ofstream* fouts) {
    // Check if the option and graph are valid
    if (option != 'Y' && option != 'N' || graph == NULL || s_vertex < 0 || e_vertex > graph->getSize() ) {      //For exception check
        return false;
    }

    const int V = graph->getSize();
    vector<int> dist(V + 1, INF);  // Initialize distance information to infinity
    vector<int> prev(V + 1, -1);   // Array to trace the path
    dist[s_vertex] = 0;            // Initialize the distance of the start vertex to 0

    // Relaxation process
    for (int i = 0; i < V - 1; ++i) {
        for (int u = 1; u <= V; ++u) {
            map<int, int> adjMap;

            // Get adjacent edges based on the option
            if (option == 'Y') {       //Directed graph
                graph->getAdjacentEdges(u, &adjMap);
            } else if (option == 'N') {     //Undirected graph
                graph->getAdjacentEdgesDirect(u, &adjMap);
            }

            for (const auto& edge : adjMap) {
                int v = edge.first, weight = edge.second;
                if (dist[u] != INF && dist[v] > dist[u] + weight) {     //Compare distance
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                }
            }
        }
    }

    // Check for negative cycles
    for (int u = 1; u <= V; ++u) {
        map<int, int> adjMap;
        graph->getAdjacentEdges(u, &adjMap);
        for (const auto& edge : adjMap) {
            int v = edge.first, weight = edge.second;
            if (dist[u] != INF && dist[v] > dist[u] + weight) {      // Negative cycle exists
                return false;
            }
        }
    }

    // Print the result
    *fouts << "======Bellman-Ford======" << endl;
    if (option == 'Y') {
        *fouts << "Directed Graph Bellman-Ford result\n";
    } else if (option == 'N') {
        *fouts << "Undirected Graph Bellman-Ford result\n";
    }

    // Check if the destination vertex is reachable
    if (dist[e_vertex] == INF) {
        *fouts << "x" << endl;  // Unreachable case
    } else {
        // Trace and print the path
        vector<int> path;
        for (int at = e_vertex; at != -1; at = prev[at]) {
            path.push_back(at);
        }
        reverse(path.begin(), path.end());

        // Check for invalid input
        if (path.front() != s_vertex || path.size() < 2) {
            *fouts << "Error" << endl;  // Incorrectly provided vertices
            return false;
        }

        for (size_t i = 0; i < path.size(); ++i) {
            // Print an arrow "->" between vertices (except for the first vertex)
            if (i > 0) *fouts << " -> ";
                *fouts << path[i];
        }

        // Move to the next line after printing the path
        *fouts << endl;
        // Output the total cost or distance from the source to the end vertex
        *fouts << "cost: " << dist[e_vertex] << endl;
    }
    *fouts << "========================\n\n";
    return true;
}


// Execute Floyd-Warshall algorithm to find all-pairs shortest paths
bool FLOYD(Graph* graph, char option, ofstream* fouts) {
    // Check if the option and graph are valid
    if (option != 'Y' && option != 'N' || graph == NULL) {      //For exception check
        return false;
    }

    const int V = graph->getSize();
    vector<vector<int>> dist(V, vector<int>(V, INF)); // 2D vector to store shortest distance information

    // Initialize the distance information from the graph
    for (int u = 0; u < V; ++u) {
        dist[u][u] = 0; // Distance to itself is 0
        map<int, int> adjMap;

        // Get adjacent edges based on the option
        if (option == 'Y') {
            graph->getAdjacentEdges(u + 1, &adjMap);
        } else if (option == 'N') {
            graph->getAdjacentEdgesDirect(u + 1, &adjMap);
        }

        for (const auto& edge : adjMap) {
            int v = edge.first - 1; // Adjust index
            int weight = edge.second;
            dist[u][v] = weight; // Set the distance from u to v as weight
            if (option == 'N') {
                dist[v][u] = weight; // For undirected graphs, the opposite direction is the same
            }
        }
    }

    // Execute the Floyd-Warshall algorithm
    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                // Check if there is a shorter path from vertex i to j through vertex k
                if (dist[i][k] != INF && dist[k][j] != INF) {
                    // Update the shortest distance between vertices i and j if a shorter path is found
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }

    // Check for negative cycles
    for (int i = 0; i < V; ++i) {
        if (dist[i][i] < 0) {
            // Negative cycle detected
            return false;  // Indicate failure due to negative cycle
        }
    }

    // Print the result
    *fouts << "==========Floyd=========" << endl;
    if (option == 'Y') { // Directed graph
        *fouts << "Directed Graph BFS result" << endl;
    } else if (option == 'N') { // Undirected graph
        *fouts << "Undirected Graph BFS result" << endl;
    }

    // Print column indices
    *fouts << "   ";
    for (int i = 0; i < V; ++i) {
        *fouts << "[" << i + 1 << "]";
    }
    *fouts << endl;

    // Print the distance matrix
    for (int i = 0; i < V; ++i) {
        *fouts << "[" << i + 1 << "]";
        for (int j = 0; j < V; ++j) {
            if (dist[i][j] == INF) {       //IF we can't go to that vertex
                *fouts << " " << "X ";
            } else {
                *fouts << " " << dist[i][j] << " ";
            }
        }
        *fouts << endl;
    }
    *fouts << "========================\n\n";

    return true;
}

class SegmentTree {
    vector<int> tree;
    int size;
    vector<int>* kw_graph;  // Added member to store kw_graph

    // Build the segment tree recursively
    void buildTree(int node, int start, int end, const vector<int>& edgesCount) {
        // If the start index equals the end index, it means we have reached a leaf node
        if (start == end) {
            // Set the value of the leaf node to the edge count of the corresponding vertex
            tree[node] = edgesCount[start];
            return;
        }

        // Calculate the middle index of the current range
        int mid = (start + end) / 2;

        // Recursively build the left and right subtrees
        buildTree(node * 2, start, mid, edgesCount);
        buildTree(node * 2 + 1, mid + 1, end, edgesCount);

        // Set the value of the current node to the sum of its left and right children
        tree[node] = tree[node * 2] + tree[node * 2 + 1];
    }

    // Update the segment tree with a new value at the specified index
    void update(int node, int start, int end, int idx, int value) {
        // If the specified index is outside the current range, return without updating
        if (idx < start || idx > end) return;

        // If the current node represents a leaf node (start equals end), update the value
        if (start == end) {
            tree[node] += value;
            return;
        }

        // Calculate the middle index of the current range
        int mid = (start + end) / 2;

        // Recursively update the left and right subtrees
        update(node * 2, start, mid, idx, value);
        update(node * 2 + 1, mid + 1, end, idx, value);

        // Update the current node with the sum of its left and right children
        tree[node] = tree[node * 2] + tree[node * 2 + 1];
    }

public:
    // Constructor to initialize the segment tree with the given size and kw_graph
    SegmentTree(int n, vector<int>* kw_graph) : size(n), kw_graph(kw_graph) {
        tree.resize(4 * n);
        buildTree(1, 0, n - 1, *kw_graph);
    }

    // Public method to update the value at the specified index
    void updateValue(int idx, int value) {
        update(1, 0, size - 1, idx, value);
    }
};

bool KWANGWOON(Graph* graph, int vertex, ofstream *fouts) {
    if(graph->getType() == 1){          //if graph type is matrix, we don't execute
        return false;
    }
    const int V = graph->getSize();
    vector<bool> visited(V + 1, false);
    vector<int> paths; // Used to store the traversal path
    vector<int> edgesCount(V + 1, 0);

    // Initialize edgesCount with the number of edges for each vertex
    for (int i = 1; i <= V; ++i) {
        map<int, int> adjMap;
        graph->getAdjacentEdgesDirect(i, &adjMap);
        edgesCount[i] = adjMap.size();
    }

    // Initialize the Segment Tree
    SegmentTree segTree(V, &edgesCount);

    int current = vertex;
    visited[current] = true;
    paths.push_back(current);
    segTree.updateValue(current, -edgesCount[current]); // Set the edge count of the starting vertex to 0

    while (true) {
        // Retrieve the adjacent edges of the current vertex using the 'getAdjacentEdgesDirect' method
        map<int, int> adjMap;
        graph->getAdjacentEdgesDirect(current, &adjMap);

        // Initialize a vector to store the next vertices to visit
        vector<int> nextVertices;

        // Iterate through the adjacent edges and check if the neighboring vertices are not visited
        for (const auto& adj : adjMap) {
            if (!visited[adj.first]) {
                // Add the neighboring vertex to the list of next vertices to visit
                nextVertices.push_back(adj.first);
            }
        }

        // Check if there are no unvisited neighboring vertices; if so, exit the loop
        if (nextVertices.empty()) break;

        // Determine the next vertex to visit based on the KWANGWOON algorithm rules
        sort(nextVertices.begin(), nextVertices.end());
        int next = (nextVertices.size() % 2 == 0) ? nextVertices.front() : nextVertices.back();

        visited[next] = true;
        paths.push_back(next);

        // Update the Segment Tree
        segTree.updateValue(current, -1); // Decrement the edge count for the current vertex
        segTree.updateValue(next, -1);    // Decrement the edge count for the next vertex

        current = next;
    }
    // Output the traversal path
    *fouts << "=========KWANGWOON========" << endl;
    *fouts << "startvertex: " << vertex << endl;
    for (size_t i = 0; i < paths.size(); ++i) {
        *fouts << paths[i];
        if (i < paths.size() - 1) {
            *fouts << " -> ";
        }
    }
    *fouts << endl;
    *fouts << "==========================\n\n";

    return true;
}



