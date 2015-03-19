/* GraphRoleUtil: provides utility functions for various graph decomposition tasks, which are useful pre- or post-
   processing steps for RoleSim.  However, the tasks are actually generic and may be useful for many other
   applications dealing with graphs.  It is kept separate from the Graph class because it requires enlarging
   the amount and type of data stored within the Graph
*/
#include "GraphRoleUtil.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <queue>
#include <utility>
#include <map>

//////////////////////////////////////////////////////////////////////////////////
GraphRoleUtil::GraphRoleUtil() {}
//////////////////////////////////////////////////////////////////////////////////
GraphRoleUtil::GraphRoleUtil(const Graph& gr) {
	g = gr;
}
//////////////////////////////////////////////////////////////////////////////////
GraphRoleUtil::~GraphRoleUtil() {}
//////////////////////////////////////////////////////////////////////////////////
vector< vector<int> > GraphRoleUtil::findComponents(Graph& g, bool verbose) {
	queue<int> searchQueue;
	EdgeList eList;
	EdgeList::iterator eit;
	vector<int> comp;
	vector< vector<int> > components;
	int N = g.num_vertices();
	int numComp = 0;
	int numVisited = 0;
	int v;

	cout << "Beginning GraphRoleUtil::findComponents" << endl;
	if (verbose) cout << "verbose" << endl;
	// Clear vertex flags
	for (v = 0; v < N; v++) {
		g[v].visited = false;
	}
	
	while (numVisited < N) {
		// find an unvisited vertex
		v = 0;
		while (g[v].visited) {
			v++;
		}
		if (verbose) cout << "new component starting at vertex " << v << endl;

		// put the unvisited vertex in the (empty) queue
		g[v].visited = true;
		searchQueue.push(v);
		comp.push_back(v);
		numVisited++;

		if (searchQueue.size() != 1) cout << "ERROR: searchQueue.size() = " << searchQueue.size() << endl;
		// use BFS to find all connected unvisited vertices
		while (!searchQueue.empty()) {
			v = searchQueue.front();
			searchQueue.pop();
			//if (verbose) cout << "pop.";
			eList = g[v].edgeList;
			//if (verbose) cout << "el.";
			for (eit = eList.begin(); eit < eList.end(); eit++) {
				if (!g[*eit].visited) {
					//if (verbose) cout <<"+";
					g[*eit].visited = true;
					searchQueue.push(*eit);
					comp.push_back(*eit);
					numVisited++;
				}
			}
			//if (verbose) cout << endl;
		}
		components.push_back(comp);
		comp.clear();
		numComp++;
	}
	int nSingle = 0;
	cout << "There are " << components.size() << " components:" << endl;
	for (int i = 0; i < components.size(); i++) {
		if (components[i].size() == 1) {
			nSingle++;
		}
		else {
			cout << "\tComponent " << i << " is size " << components[i].size() << endl;
			if (verbose) {
				vector<int>::iterator vit;
				for (vit = components[i].begin(); vit != components[i].end(); vit++) {
					cout << *vit << " ";
				}
				cout << endl << endl;
			}
		}
	}
	cout << "\t" << nSingle << " remaining components are singletons." << endl;

	cout << "Finishing GraphRoleUtil::findComponents" << endl;
	return components;
}

///////////////////////////////////////////////////////////////////////////
/* findKCoreAssignment -         
	Assign each vertex to a K-core, where K is the largest integer such that
	vertex v has at least K links to other members of the K-core.
	Return a vector of K values
	This implements "An O(m) algorithm for Cores Decomposition of Networks"
	by Batagelj et al.
 */
////////////////////////////////////////////////////////////////////////////
vector<int> GraphRoleUtil::findKCoreAssignment(Graph& g, bool verbose) {
	vector<int> core;		// Initially, degree values. Finally, core values.
	vector<int> sortedVert;		// vertex list, sorted by increasing degree
	vector<int> pos;			// position of vertex v in sortedVert
	vector<int> binStart;		// starting position (in sortedVert) for each bin
	int d, v;
	int maxDeg;
	EdgeList::iterator uit;
	// Note: pos and sortedVert are inverse maps of one another, so
	// sortVert[pos[v]] = v,  and pos[sortVert[i]] = i

	cout << "Beginning GraphRoleUtil::findKCoreAssignment" << endl;
	if (verbose) cout << "verbose" << endl;
	int N = g.num_vertices();
	sortedVert = vector<int>(N);
	pos = vector<int>(N);

	if (verbose) cout << "Step 1: Initialize degree vector and find out the maximum degree" << endl;
	// 1. Initialize degree vector and find out the maximum degree
	maxDeg = 0;
	for (v = 0; v < N; v++) {
		core.push_back(g[v].edgeList.size());
		maxDeg = (core[v] > maxDeg) ? core[v] : maxDeg;
	}
	if (verbose) cout << "maxDeg = " << maxDeg << endl;

	if (verbose) cout << "Step 2. Count the size of each degree bin and initialize the bin pointers" << endl;
	// 2. Count the size of each degree bin and initialize the bin pointers
	//binStart = new int[maxDeg + 1];
	binStart = vector<int>(maxDeg + 1);

	for (d = 0; d <= maxDeg; d++) {
		binStart[d] = 0;
	}
	for (v = 0; v < N; v++) {	// Count bin sizes
		binStart[core[v]]++;
	}
	int offset = 0;
	int binSize;
	for (int d = 0; d <= maxDeg; d++) {	// binStart[d] = cumulative size of preceding bins
		binSize = binStart[d];
		cout << "Degree " <<d<< ": " <<binSize<< " vertices" << endl;
		binStart[d] = offset;
		offset += binSize;
	}


	if (verbose) cout << "Step 3. Initialize pos[v] and sortedVert[v]" << endl;
	// 3. Initialize pos[v] and sortedVert[v].  Note: we shift binStart values.
	for (v = 0; v < N; v++) {
		pos[v] = binStart[core[v]];	// pos[v] = next open slot for v's degree
		//if (verbose) cout << "pos["<<v<<"] = " << binStart[core[v]] << "\t";
		sortedVert[pos[v]] = v;		// put value v at pos[v] in sortedVert
		//if (verbose) cout << "sortedVert["<<pos[v]<<"] = " << v << endl;
		binStart[core[v]] += 1;		// update the position of the next open slot
	}
	for (d = maxDeg; d > 0; d--) {	// retore binStart to the correct values.
		binStart[d] = binStart[d-1];
		if (verbose) cout << "Bin " << d << " starts at index " << binStart[d] << endl;
	}
	binStart[0] = 0;
	if (verbose) cout << "Bin " << 0 << " starts at index " << binStart[0] << endl;

	if (verbose) {
		cout << "pos: ";
		for (int v = 0; v < N; v++) {
			cout << pos[v] << " ";
		}
		cout << endl;
		cout << "sortedVert: ";
		for (int v = 0; v < N; v++) {
			cout << sortedVert[v] << " ";
		}
		cout << endl;
	}

	if (verbose) cout << "Step 4. For each neighbor u of each vertex v, if core[u] > core[v], decrement core[u]." << endl;
	// 4. For each neighbor u of each vertex v, if core[u] > core[v],
	// decrement core[u]. Update pos[], sortedVert[], and bin[] as needed.
	int i, u, w, temp;
	for (i = 0; i < N; i++) {
		//if (verbose) cout << "i = " << i << endl;
		v = sortedVert[i];
		//if (verbose) cout << "v = sortedVert[" << i << "] = " << v << endl;
		for (uit = g[v].edgeList.begin(); uit != g[v].edgeList.end(); uit++) {
			u = *uit;
			//if (verbose) cout << "  u = *uit = " << u << endl;
			//if (verbose) cout << "  if (core[u]=" << core[u] << " > core[v]=" << core[v] << ")" << endl;
			if (core[u] > core[v]) {
				//if (verbose) cout << "    w = sortedVert[binStart[core[u]]] = " << sortedVert[binStart[core[u]]] << endl;
				w = sortedVert[binStart[core[u]]]; // w = first vertex in u's bin
				temp = pos[u];				// swap positions of u and w
				pos[u] = pos[w];
				pos[w] = temp;

				sortedVert[pos[u]] = u;		// update sortedVec values
				sortedVert[pos[w]] = w;

				binStart[core[u]] += 1;		// shift the start of u's (old) bin
				core[u] = core[u] - 1;		// move u to the next smaller bin
			}
		}
	}
	cout << "Finishing GraphRoleUtil::findKCoreAssignment" << endl;
	return core;
}

///////////////////////////////////////////////////////////////////////
/* makeBlockShells -
	(1) Read a file which contains the quantity and sizes of the blocks of the graph,
	whose format is a single line:
	'num_blocks size_0 size_1 ... size_(n-1)'
	(2) Output a vector of vectors, where each is the IDs assigned to the vertices of each block, i.e.
	vector0 = (0, 1, 2, ... , size0 - 1)
	vector1 = (size0, size0 + 1, ... size0 + size1 - 1)
	...
	vectorn-1 = (size0+size1+...+sizen-2, ... , size0+size1+...+sizen-1 - 1)
 */
vector< vector<int> > GraphRoleUtil::makeBlockShells(char* blockFilename, bool verbose) {

	ifstream in(blockFilename);
	vector< vector<int> > shells;

	cout << "Beginning GraphRoleUtil::makeBlockShells" << endl;
	int numShells, shellSize;
	in >> numShells;
	shells = vector< vector<int> >(numShells);

	int vertexId = 0;
	for (int i = 0; i < numShells; i++) {
		in >> shellSize;
		for (int j = 0; j < shellSize; j++) {
			shells[i].push_back(vertexId++);
		}
	}
	if (verbose) {
		for (int i = 0; i < shells.size(); i++){
			cout << i << ": " << *(shells[i].begin()) << " - " << *(shells[i].end()) << endl;
		}
	}
	cout << "Finishing GraphRoleUtil::makeBlockShells" << endl;
	return shells;
}
///////////////////////////////////////////////////////////////////////
/* findKshells -
	(1) calls findKCoreAssignment, (2) returns the vertex lists
	corresponding to the max K-core and its surrounding K-shells
 */
vector< vector<int> > GraphRoleUtil::findKShells(Graph& g, bool verbose) {

	int maxK = 0;
	int N = g.num_vertices();
	vector< vector<int> > Kshells;

	cout << "Beginning GraphRoleUtil::findKShells" << endl;
	if (verbose) cout << "verbose" << endl;
	vector<int> coreNum = findKCoreAssignment(g, verbose);
	
	// Find out maxK
	for (int v = 0; v < N; v++) {
		maxK = (coreNum[v] > maxK) ? coreNum[v]: maxK;
	}

	// Create empty Kshells list
	vector<int> emptyVecInt;
	for (int k = 0; k <= maxK; k++) {
		Kshells.push_back(emptyVecInt);
	}

	// Finally, assign vertices to shells
	for (int v = 0; v < N; v++) {
		Kshells[coreNum[v]].push_back(v);
	}

	for (int k = 0; k <= maxK; k++) {
		cout <<"Shell"<< k <<": "<< Kshells[k].size() <<" vertices"<< endl;

		if (verbose && Kshells[k].size() > 1) {
			vector<int>::iterator vit;
			for (vit = Kshells[k].begin(); vit != Kshells[k].end(); vit++) {
				cout << *vit << " ";
			}
			cout << endl;
		}
	}
	cout << "Finishing GraphRoleUtil::findKShells" << endl;
	return Kshells;
}

//////////////////////////////////////////////////////////////////////////
/* medusaDecomposition -
	(1) Finds K-shells
	(2) Merges the highest K-cores until their size is at least minCoreSize
	(3) Finds the largest component not in the core --> Body
	(4) Returns a vector of vertex lists for components, where
		component[0] = Core, component[1] = Body, followed by remaining
		peripheral components
 */
 /////////////////////////////////////////////////////////////////////////
vector< vector<int> > GraphRoleUtil::medusaDecomposition(Graph& g, int minCoreSize, bool verbose) {
	vector< vector<int> > comp;
	vector<int> core, nonCore, body;
	Graph periphery;
	int maxK, coreLevel, bodyIndex;
	
	cout << "Beginning GraphRoleUtil::medusaDecomposition" << endl;
	if (verbose) cout << "verbose" << endl;
	vector< vector<int> > Kshells = findKShells(g, verbose);

	// Find the layer that satisfies minCoreSize, then
	// construct vertex lists for core and periphery
	maxK = Kshells.size() - 1; // first Kshell is Kshell[0]
	
	coreLevel = maxK; 
	core = Kshells[maxK];
	while (core.size() < minCoreSize) {
		coreLevel--;
		core.insert(core.end(), Kshells[coreLevel].begin(), Kshells[coreLevel].end());
	}
	if (verbose) cout << "core has " << core.size() << " vertices" << endl;
	for (int i = 0; i < coreLevel; i++) {
		nonCore.insert(nonCore.end(), Kshells[i].begin(), Kshells[i].end());
	}
	if (verbose) cout << "nonCore has " << nonCore.size() << " vertices" << endl;

	// Find all the connected components in the periphery
	periphery = g.inducedSubgraph(nonCore); // victor
	periphery.sortAndSetUndirectedEdges();
	if (verbose) periphery.printGraph();

	comp = findComponents(periphery, verbose);

	// Finally, construct the medusa.  Find the largest periphery component,
	// and move it to the front of the component list.  Insert the core in front
	// of the largest component.
	bodyIndex = 0;
	for (int i = 0; i < comp.size(); i++) {
		bodyIndex = (comp[i].size() > comp[bodyIndex].size()) ? i : bodyIndex;
	}
	//if (verbose) cout << "The largest periph component is comp["<<bodyIndex<<"] with size " << comp[bodyIndex].size() << endl;
	body = comp[bodyIndex];
	comp[bodyIndex] = comp[0];
	comp[0] = body;
	//if (verbose) cout << "After swapping, comp[0] is now size " << comp[0].size() << endl;

	comp.insert(comp.begin(), core);
	
	cout << coreLevel << "-core has " << comp[0].size() << " vertices" << endl;
	cout << "Body component has " << comp[1].size() << " vertices" << endl;
	cout << "Finishing GraphRoleUtil::medusaDecomposition" << endl;
	return comp;	
}

/////////////////////////////////////////////////////////////////////////
/* given an array of integers, write a text file in tabular format,
 * where each line contains two integers:
 * 'array[i]    i'
 * If the array represents the sorted subset of the nodes in a subgraph,
 * then the table is a map of vertex IDs from the original graph to the subgraph.
 */
 /////////////////////////////////////////////////////////////////////////
void GraphRoleUtil::writeSubgraphIdMap(const VectIntType& sublist, ostream& out) {
	int subSize = sublist.size();
	for (int i = 0; i < subSize; i++) {
		out <<	sublist[i] << "\t" << i << endl;
	}
}

/////////////////////////////////////////////////////////////////////////
/* Write each component, formatted as a list of node ID numbers, to a text file.
 * Write one file for each component, plus one common file of Matlab commands
 * which load each file into a vector.
 */
/////////////////////////////////////////////////////////////////////////
void GraphRoleUtil::writeComponents(vector< vector<int> >& compVec, string filePrefix) {
	
	string compFilename, matlabFilename, matlabPrefix, matlabLine;
	ostringstream ss;
	vector<int>::iterator it;

	matlabFilename = filePrefix + "_compMatCmd.txt";
	ofstream mat(matlabFilename.c_str());

	// make a filename prefix excluding any directory path
	matlabPrefix = filePrefix.substr(filePrefix.find_last_of("/") + 1, filePrefix.size());

	for (int i = 0; i < compVec.size(); i++) {
		ss << i;
		compFilename = filePrefix + ss.str() + ".txt";
		ofstream out(compFilename.c_str());
		for (it = compVec[i].begin(); it != compVec[i].end(); it++) {
			out << *it << " ";
		}
		out.close();

		matlabLine = matlabPrefix + ss.str() + " = load(" + compFilename + ");" ;
		mat << matlabLine << endl;
		ss.str("");
		ss.clear();
	}
	mat.close();
}

//////////////////////////////////////////////////////////////////////////////////
Graph GraphRoleUtil::findGiantComponent(Graph& g, bool verbose) {
	vector< vector<int> > Components = findComponents(g, verbose);
	int giantId = 0;

	// Find the largest component
	for (int i = 1; i < Components.size(); i++) {
		if (Components[i].size() > Components[giantId].size()) {
			giantId  = i;
		}
	}

	vector<int> giant = Components[giantId];
	sort(giant.begin(), giant.end());

	// print a table or map of the giant components vertices
	char* subgraphMap = "subgraphMap.txt";
	ofstream map(subgraphMap);
	for (int i = 0; i < giant.size(); i++) {
		map << giant[i] << "\t" << i << endl;
	}
	
	Graph giantSubgraph = g.inducedSubgraph(giant, true);
	return giantSubgraph;

}
//////////////////////////////////////////////////////////////////////////////////

