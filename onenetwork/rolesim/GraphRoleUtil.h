// GraphRoleUtil - functions for preprocessing graphs
#ifndef _GRAPH_ROLEUTIL_H
#define _GRAPH_ROLEUTIL_H

#include "Graph.h"
#include <vector>

using namespace std;

struct EdgeSupplement {
public:
	int		src;
	int		dst;
	float	weight;
	vector<int> cycleList;
};
typedef vector<EdgeSupplement> EdgeSuppList;
struct edgeComp {
	bool operator() (EdgeSupplement a, EdgeSupplement b) {
		return (a.weight < b.weight);
	}
};


struct VertexSupplement {
	bool			inCurrentPath;
	EdgeSuppList	cycleEdges;
	
};
typedef vector<VertexSupplement> VertexSuppList;

class GraphRoleUtil {
private:
	Graph g;
	int nCycles;
	VertexSuppList vxl;
	bool verbose;
	int nVerticesVisited;
	int nEdgesInCycles;

public:
	GraphRoleUtil();
	GraphRoleUtil(const Graph&);
	~GraphRoleUtil();

	vector< vector<int> > findComponents(Graph& g, bool verbose = false);
	vector<int> findKCoreAssignment(Graph& g, bool verbose = false);
	vector< vector<int> > findKShells(Graph& g, bool verbose = false);
	void writeComponents(vector< vector<int> >& shells, string fileprefix);
	vector< vector<int> > medusaDecomposition(Graph& g, int minCoreSize, bool verbose = false);
	void writeSubgraphIdMap(const VectIntType&, ostream& out); // victor

	Graph findGiantComponent(Graph& g, bool verbose = false);
	vector< vector<int> > makeBlockShells(char* blockFilename, bool verbose);	

};


#endif
