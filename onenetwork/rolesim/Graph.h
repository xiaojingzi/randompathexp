#ifndef _GRAPH_H
#define _GRAPH_H
// Graph: a generic graph model, implemented as a vector of edge lists.

//#include "Util.h"
#include <iostream>
#include <vector>
using namespace std;
/****
Graph file format:
The first line must be the string "graph_for_greach" (for silly historical reasons)
The second line must contain a single integer, the number of vertices
The remaining lines are edgelists of the form: "n: d1 d2 d3 ... dm #"
	where n is a node id, and {d1, d2,..., dm} are the ids of its neighbors.
	Node ids must be integers, numbered from 0.

If directed:
(1) Edgelists are interpreted as lists of out-neighbors
(2) For each neighbor di of each node n,
	add di into n's outList (forwards) and add n to di's inList (reverse)
(3) Set vl[n].edgeList = union (n's outlist and inlist)

If undirected:
(1) Edgelists are interpreted as lists of undirected neighbors
(2) For each neighbor di of each node n,
	add di into n's outList (forwards) and add n to di's inList (reverse)
	(so far, same steps as for directed graph)
(3) Set vl[n].edgeList = union (n's outlist and inlist)
(4) Update n's outlist <-- edgeList (that is, union(outlist, inlist) )
*/

typedef vector<int> VectIntType;
typedef vector<float> VectFloatType;
typedef vector<VectIntType> VectVectIntType;
typedef vector<VectFloatType> VectVectFloatType;
//typedef map<VectIntType, int> MapVectIntToInt;

typedef vector<int> EdgeList;	// edge list represented by vertex id list

struct Vertex {
	bool visited;			// used in GraphRoleUtil.findComponents
	//bool fat;	// fat node
	//int topo_id;	// topological order
	//int path_id;	// path id
	//int dfs_order;
	//int pre_order;	
	//int post_order;
	//int first_visit; // for test
	float weight;
	EdgeList edgeList; // undirected
	int degree;
};

typedef vector<Vertex> VertexList;	// vertices list (store real vertex property) indexing by id

struct In_OutList {
	EdgeList inList;
	EdgeList outList;
};
typedef vector<In_OutList> GRA;	// index graph

class Graph {
	private:
		GRA graph;
		VertexList vl;
		int vsize;
		bool directed;
		string name;
		
	public:
		Graph();
		Graph(int);
		Graph(GRA&, VertexList&);
		Graph(istream&);
		~Graph();
		void readGraph(istream&);
		void writeGraph(ostream&);
		void printGraph();
		void addVertex(int);
		void addEdge(int, int, bool directed=false);
		int num_vertices();
		int num_edges();
		VertexList& vertices();
		EdgeList& out_edges(int);
		EdgeList& in_edges(int);
		int out_degree(int);
		int in_degree(int);
		vector<int> getRoots();
		bool hasEdge(int, int);	
		Graph& operator=(const Graph&);
		Vertex& operator[](const int);
		void setDirected(bool);
		bool isDirected();
		void sortAndSetUndirectedEdges();
		void sortDirectedEdges (); // victor
		Graph inducedSubgraph(VectIntType, bool reMap=false, bool directed=false); // victor
		void appendGraph(Graph& g2, bool directed=false); // victor
};
#endif
