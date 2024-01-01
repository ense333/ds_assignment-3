#include "Manager.h"
#include "GraphMethod.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

Manager::Manager()	
{
	graph = nullptr;	
	fout.open("log.txt", ios::app);
	load = 0;	//Anything is not loaded
}

Manager::~Manager()
{
	if(load)	//if graph is loaded, delete graph
		delete graph;	
	if(fout.is_open())	//if fout is opened, close file
		fout.close();	//close log.txt File
}

void Manager::run(const char* command_txt){
	ifstream fin;	//Command File Input File Stream
	fin.open(command_txt, ios_base::in);//File open with read mode

	if(!fin) { //If command File cannot be read, Print error
		fout<<"command file open error"<<endl;
		return;	//Return
	}
	string commandLine = "";    //Get command line
	while(getline(fin, commandLine)){
		istringstream iss(commandLine);
		string commandStart = "";
		iss >> commandStart;
		if(commandStart == "LOAD"){		//If command is LOAD command
			string fileType = "";
			iss >> fileType;
			const char* fileName = fileType.c_str();
			if(!LOAD(fileName)){
				printErrorCode(100);
			}else{
				fout << "========== LOAD ==========\n" << "Success\n" << "==========================\n\n";
			} 
		}else if(commandStart == "PRINT"){		//IF command is PRINT
			if(!PRINT()){
				printErrorCode(200);
			}
		}else if(commandStart == "BFS"){		//If command is BFS
			int startVertex = 0;
			char directed;
			iss >> directed >> startVertex;		//Need to get direction value, startvertex
			if(!mBFS(directed, startVertex)){	//If there is error in calling BFS
				printErrorCode(300);	//PRINT error code
			} 
		}else if(commandStart == "DFS"){		//IF command is DFS
			int startVertex = 0;
			char directed;
			iss >> directed >> startVertex;		//Need to get direction value, startvertex
			if(!mDFS(directed, startVertex)){	//If there is error in calling DFS
				printErrorCode(400);			//PRINT error code
			}
		}else if(commandStart == "KWANGWOON"){	//IF command is KWANGWOON
			if(!mKwoonWoon(1)){		//Since we start Kwangwoon alogorithm with vertex 1
				printErrorCode(500);		//If error
			}
		}else if(commandStart == "KRUSKAL"){	//IF command is KRUSKAL
			if(!mKRUSKAL()){		//Only need to call function
				printErrorCode(600);		//If there is error
			}
		}else if(commandStart == "DIJKSTRA"){	//IF command is DIJKSTRA
			int startVertex = 0;
			char directed;
			iss >> directed >> startVertex;		//Need direction option, start vertex
			if(!mDIJKSTRA(directed, startVertex)){	
				printErrorCode(700);		//If error
			}

		}else if(commandStart == "BELLMANFORD"){	//IF command is BELLMANFORD
			int startVertex = 0, endVertex = 0;	
			char directed;
			iss >> directed >> startVertex >> endVertex;		//Need direction option, start, end vertex
			if(!mBELLMANFORD(directed, startVertex, endVertex)){
				printErrorCode(800);
			}

		}else if(commandStart == "FLOYD"){	//IF command is FLOYD
			char directed;
			iss >> directed;	//Need direction option
			if(!mFLOYD(directed)){
				printErrorCode(900);
			}
		}else if(commandStart == "EXIT"){	//IF command is Exit
			fout << "========== EXIT ==========\n" << "Success\n" << "==========================\n\n";		//Print Success 
			fin.close();	//close file
			return;		//Then close function
		}else{		//IF command input is wrong
			printErrorCode(1000);
		}
	}
	
	fin.close();  //close file
	return;
}


bool Manager::LOAD(const char* filename)
{
	string type, eachString;
	int numVertices = 0;
	if(strcmp(filename, "graph_L.txt") == 0){	//If graph type is LIST
		ifstream fin(filename);
		if(!fin.is_open()){		//If error opening file
			cerr << "Couldn't open file" << endl;
		}
		fin >> type;		//For graph type
		fin >> numVertices;		//Number of vertex
		if(graph){		//If graph exists
			delete graph;	//Delete existing graph
		}
		ListGraph* listG = new ListGraph(0, numVertices);   //Create new List graph
		graph = listG;

		while(getline(fin, eachString)){	//Read graph data and store to graph 
			istringstream iss(eachString);
			int from, to, weight;
			if(!(iss >> to)){		//Since there will be vertex that goes out to nothing
				cout << "";
			}
			if(iss >> weight){		//If we get vertex, weight value
				listG->insertEdge(from, to, weight);		//using insertEdge function to create Graph
			}else{  //If there is vertex that have out-degrees, in-degrees
				//Only one value
				from = to;
			}
		}
		return true;		//Successful

	}else if(strcmp(filename, "graph_M.txt") == 0){		//If graph type is matrix

		ifstream fin(filename);
		if(!fin.is_open()){		//If error opening file
			cerr << "Couldn't open file" << endl;
		}
		fin >> type;	//For graph type
		fin >> numVertices;		//Number of vertex
		if(graph){		//If graph exists
			delete graph;	//Delete existing graph
		}
		MatrixGraph* matrixG = new MatrixGraph(1, numVertices);		//Create new matrix graph
		graph = matrixG;	//allocate graph
		int colIndex = 1;
		getline(fin, eachString);
		while(getline(fin, eachString)){		//Read graph data and create graph
			istringstream iss(eachString);
			int weight = 0;
			for(int i = 1; i < numVertices + 1; i++){
				if(iss >> weight){
					if(weight != 0){
						matrixG->insertEdge(colIndex, i, weight);
					}
				}else{    //Not enough input has given
					graph = NULL;
					return false;
				}
				
			}
			colIndex++;
		}
		return true;
	}else{
		return false;
	}
}

bool Manager::PRINT()	
{
	if(graph == NULL){
		//printErrorCode(200);
		return false;
	}else{
		graph->printGraph(&fout);
		return true;
	}
}

bool Manager::mBFS(char option, int vertex)	
{
	if(!BFS(graph, option, vertex, &fout)){		//If there is error when calling BFS
		return false;
	}
	return true;
}

bool Manager::mDFS(char option, int vertex)	
{
	if(!DFS(graph, option, vertex, &fout)){		//If there is error when calling DFS
		return false;
	}
	return true;
}

bool Manager::mDIJKSTRA(char option, int vertex)	
{
	if(!Dijkstra(graph, option, vertex, &fout)){		//If there is error when calling DIJKSTRA
		return false;
	}
	return true;
}

bool Manager::mKRUSKAL()
{
 	if(!Kruskal(graph, &fout)){			//If there is error when calling KRUSKAL
		return false;
	}
	return true;
}

bool Manager::mBELLMANFORD(char option, int s_vertex, int e_vertex) 
{
	if(!Bellmanford(graph, option, s_vertex, e_vertex, &fout)){		//If there is error when calling BELLMANFORD
		return false;
	}
	return true;
}

bool Manager::mFLOYD(char option)
{
	if(!FLOYD(graph, option, &fout)){	//If there is error when calling FLOYD
		return false;
	}
	return true;
}

bool Manager::mKwoonWoon(int vertex) {
	if(!KWANGWOON(graph, vertex, &fout)){		//If there is error when calling KWANGWOON
		return false;
	}
	return true;
}

void Manager::printErrorCode(int n)
{
	fout<<"========ERROR======="<<endl;
	fout<<n<<endl;
	fout<<"===================="<<endl << endl;
}


