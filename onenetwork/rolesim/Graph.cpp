#include "Graph.h"
// Graph: a generic graph model, implemented as a vector of edge lists.
#include <cstdio>

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <map>
#include <set>
#include <vector>
#include <iterator>
#include <list>
#include <deque>
#include <assert.h>
#include <utility>
#include <cmath>
#include <iomanip>
#ifndef _MSC_VER
	#include <sys/time.h>
#endif
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

Graph::Graph() {
	graph = GRA();
	vl = VertexList();
}

Graph::Graph(int size) {
	vsize = size;
	vl = VertexList(size);
	graph = GRA(size, In_OutList());
}

Graph::Graph(GRA& g, VertexList& vlist) {
	vsize = vlist.size();
	graph = g;
	vl = vlist;
}

Graph::Graph(istream& in) {
	readGraph(in);
}

Graph::~Graph() {}

void Graph::printGraph() {
	writeGraph(cout);
}

void Graph::setDirected(bool d) { directed = d; }

bool Graph::isDirected() { return directed; }

void Graph::readGraph(istream& in) {
	fprintf(stderr,"READ GRAPH STARTING ... \n"); fflush(NULL);
	string buf, token;
	getline(in, buf);
	istringstream ss(buf);
	ss >> token;


	if (token != "graph_for_greach") {
		cout << "BAD FILE FORMAT!" << endl;
		exit(0);
	}
	
	int n;
	getline(in, buf);
	istringstream(buf) >> n;
	// initialize
	vsize = n;
	vl = VertexList(n);
	graph = GRA(n, In_OutList());	

	for (int i = 0; i < n; i++)
		addVertex(i);

	string sub;
	int idx;
	int sid = 0;
	int tid = 0;

	while (getline(in, buf)) {
		idx = buf.find(":");
		sub = buf.substr(0, idx);
		istringstream(sub) >> sid;
		buf.erase(0, idx+2);
		while (buf.find(" ") != string::npos) {
			sub = buf.substr(0, buf.find(" "));
			istringstream(sub) >> tid;
			buf.erase(0, buf.find(" ")+1);
			if( sid == tid )
				cout << "sid: " << sid << ", tid: " << tid << endl;
			addEdge(sid, tid);
			//printf("source id = %d, target id = %d\n", sid, tid);
		}
		++sid;
	}
}	

void Graph::writeGraph(ostream& out) {
	cout << "Graph size = " << graph.size() << endl;
	out << "graph_for_greach" << endl;
	out << vl.size() << endl;

	GRA::iterator git;
	EdgeList el;
	EdgeList::iterator eit;
	VertexList::iterator vit;
	int i = 0;
	for (i = 0; i < vl.size(); i++) {
		out << i << ": ";
		el = graph[i].outList;
		for (eit = el.begin(); eit != el.end(); eit++)
			out << (*eit) << " ";
		out << "#" << endl;
	}
/*
	cout << "** In List for graph **" << endl;
	for (i = 0; i < vl.size(); i++) {
		out << i << ": ";
		el = graph[i].inList;
		for (eit = el.begin(); eit != el.end(); eit++)
			out << (*eit) << " ";
		out << "#" << endl;
	}
*/
}

void Graph::addVertex(int vid) {
	if (vid >= vl.size()) {
		int size = vl.size();
		for (int i = 0; i < (vid-size+1); i++) {
			graph.push_back(In_OutList());
			vl.push_back(Vertex());
		}
		vsize = vl.size();
	}

	Vertex v;
	//v.visited = false;
	v.weight = 1.0;
	vl[vid] = v;
	//cout << vid << ": " << v.weight << endl;

	EdgeList il = EdgeList();
	EdgeList ol = EdgeList();
	In_OutList oil = In_OutList();
	oil.inList = il;
	oil.outList = ol;
	graph[vid] = oil;	
}

void Graph::addEdge(int sid, int tid, bool directed) {
	if (sid >= vl.size())
		addVertex(sid);
	if (tid >= vl.size())
		addVertex(tid);
	// update edge list
	graph[tid].inList.push_back(sid);
	graph[sid].outList.push_back(tid);
}

int Graph::num_vertices() {
	return vl.size();
}

int Graph::num_edges() {
	EdgeList el;
	GRA::iterator git;
	int num = 0;
	for (git = graph.begin(); git != graph.end(); git++) {
		el = git->outList;
		num += el.size();
	}
	return num;
}

// return out edges of specified vertex
EdgeList& Graph::out_edges(int src) {
	return graph[src].outList;
}

// return in edges of specified vertex
EdgeList& Graph::in_edges(int trg) {
	return graph[trg].inList;
}

int Graph::out_degree(int src) {
	return graph[src].outList.size();
}

int Graph::in_degree(int trg) {
	return graph[trg].inList.size();
}

// get roots of graph (root is zero in_degree vertex)
vector<int> Graph::getRoots() {
	vector<int> roots;
	GRA::iterator git;
	int i = 0;
	for (git = graph.begin(), i = 0; git != graph.end(); git++, i++) {
		if (git->inList.size() == 0)
			roots.push_back(i);
	}
	
	return roots;
}

// check whether the edge (src, trg) is in the graph
bool Graph::hasEdge(int src, int trg) {
	EdgeList el = graph[src].outList;
	EdgeList::iterator ei;
	for (ei = el.begin(); ei != el.end(); ei++)
		if ((*ei) == trg)
			return true;
	return false;

}

// return vertex list of graph
VertexList& Graph::vertices() {
	return vl;
}

Graph& Graph::operator=(const Graph& g) {
	if (this != &g) {
		graph = g.graph;
		vl = g.vl;
		vsize = g.vsize;
	}
	return *this;
}

// get a specified vertex property
Vertex& Graph::operator[](const int vid) {
	return vl[vid];
}
///////////////////////////////////////////////////////////////////////////////////
void Graph::sortDirectedEdges() {
	VectIntType temp;
	for (int i = 0; i < vsize; i++ ) {
		sort(graph[i].inList.begin(), graph[i].inList.end());
		sort(graph[i].outList.begin(), graph[i].outList.end());
	}
}
///////////////////////////////////////////////////////////////////////////////////
void Graph::sortAndSetUndirectedEdges()
/* Perform two critical functions:
(1) Sort the in- and out- edge lists, which is a necessary pre-condition for
set operations (set_union, set_intersection) and speeds up other operations
(2) Sort the union of in- and out-neighbors in vl[i].edgeList
This provides a ready resource for undirected graph processing.

Even better than storing the undirected neighbors in vl[i].edgeList is to
store them in g.out_edges(i).  This allows us to write functions that can handle
both directed and undirected graphs with the same code (for undirected graphs,
set inWeight = 0 */
{
	VectIntType temp;
	for (int i = 0; i < vsize; i++ ) {
		sort(graph[i].inList.begin(), graph[i].inList.end());
		sort(graph[i].outList.begin(), graph[i].outList.end());

		set_union( graph[i].inList.begin(), graph[i].inList.end(), 
			graph[i].outList.begin(), graph[i].outList.end(), back_inserter( temp ) );
		vl[i].edgeList = temp;
		vl[i].degree = temp.size();
		temp.clear();

		out_edges(i) = vl[i].edgeList;
		in_edges(i).clear();
	}
}

///////////////////////////////////////////////////////////////////////////////////
Graph Graph::inducedSubgraph(VectIntType sublist, bool reMap, bool directed) {
	// Because the Graph class only supports graphs with consecutively numbered
	// vertices, we are constrained in how subgraphs are represented.
	// If reMap = true, we map the sublist's vertex ids to a list of consecutive
	// integers starting with 0.
	// If reMap = false, we actually return a graph which contains ALL the
	// original vertices, but includes only the induced edges.
	Graph subg;
	map<int,int> idMap;
	EdgeList elist;
	EdgeList::iterator eit;
	VectIntType::iterator vit;

	cout << "Beginning Graph::inducedSubgraph" << endl;
	int subSize = sublist.size();
	cout << "sublist.size() = " << subSize << endl;

	// Sort the subgraph vertices
	sort(sublist.begin(), sublist.end());

	if (reMap) {
		for (int i = 0; i < subSize; i++) {
			idMap[sublist[i]] = i;
		}
		subg = Graph(subSize);
	}
	else {
		subg = Graph(vsize);
	}

	// For each vertex in sublist
	int numEdges = 0;
	for (vit = sublist.begin(); vit != sublist.end(); vit++) {
		// Find its neighbors that are also in sublist
		sort(graph[*vit].outList.begin(), graph[*vit].outList.end());

		// Do out-edges only
		set_intersection(graph[*vit].outList.begin(), graph[*vit].outList.end(), 
				sublist.begin(), sublist.end(), back_inserter(elist));
		for (eit = elist.begin(); eit != elist.end(); eit++) {
			if (reMap) {
				subg.addEdge(idMap[*vit], idMap[*eit], directed);
			}
			else {
				subg.addEdge(*vit, *eit, directed);
			}
			numEdges++;
		}
		elist.clear();

	}
	cout << "Induced graph has " << sublist.size() << " vertices and " << numEdges << " edges" <<endl;
	cout << "subg.num_edges() = " << subg.num_edges() << endl;
	cout << "Ending Graph::inducedSubgraph" << endl;
	return subg;

}
///////////////////////////////////////////////////////////////////////
void Graph::appendGraph(Graph& g2, bool directed) {
	cout << "Original graph has " << num_vertices() << " vertices and " << num_edges() << " edges." << endl;
	vector<int> destList;
	vector<int>::iterator dest;
	int offset = vsize;
	int prevNumEdges;
	int count = 0;
	for (int src = 0; src < g2.vsize; src++) {
		destList = g2.out_edges(src);
		for (dest = destList.begin(); dest != destList.end(); dest++) {
			prevNumEdges = num_edges();
			addEdge(src + offset, *dest + offset, directed);
			if (num_edges() != prevNumEdges + 1) {
				cout << "Error: " << prevNumEdges << " + (" << src + offset << "," + *dest + offset << ") = " << num_edges() << endl;
			}
			count++;
		}
	}
	cout << "Second graph has " << g2.vsize << " vertices and " << g2.num_edges() << " edges." << endl;
	cout << "Appended " << count << " edges to first graph" << endl;
	if (directed) {
		sortDirectedEdges();
	}
	else {
		sortAndSetUndirectedEdges ();
		cout << "Undirected graph: copying edges in the reverse direction" << endl;
	} 
	cout << "The combined graph has " << vsize << " vertices and " << num_edges() << " edges." << endl;
}

