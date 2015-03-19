/* IcebergSimMap: a class for scalable RoleSim computation.  While providing a more efficient data structure than SimMatrix
  (a nodePair->simValue map for only the most similar nodePairs, discarding nodePairs that are not very similar),
  its functions are analogous to those in SimMatrix.  In the future, this class may be re-implemented as a subclass
  of SimMatrix or some new generic similarity superclass
  */
#include "IcebergSimMap.h"

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
#include <queue>
#include <ctime>

/*
 * This class encapsulates a data structure map< <int,int>, float > which seeks to store
 * the similarity values for only the most similar vertex pairs, therefore using less memory
 * that SimMatrix, which stores the full (triangular) matrix.
 * The Initialize function uses the degrees of a given pair and the degrees of their neighbors
 * to compute an upper bound RoleSim value.  If the value is above a threshold 'waterline',
 * it is included in the Map.  The similarity for any vertex pair not in the map is defined
 * as ratioWt*degree(u)/degree(v) + beta, where ratioWt is between 0 and 1.
 */

IcebergSimMap::IcebergSimMap(void)
{
}

IcebergSimMap::~IcebergSimMap(void)
{
}
/////////////////////////////////////////////
IcebergSimMap::IcebergSimMap(float w, float r, float b, float i, Graph& g, bool verb)
{
	waterline = w;
	ratioWt = r;
	beta = b;
	inWeight =i;
	gPtr = &g;
	N = g.num_vertices();
	lowerTri = true;
	verbose = verb;
}
/////////////////////////////////////////////
float IcebergSimMap::getWaterline() { return waterline; }
float IcebergSimMap::getRatioWt() { return ratioWt; }
int   IcebergSimMap::getN() { return N; }
SimMapType& IcebergSimMap::getMap() { return simMap; }
/////////////////////////////////////////////
void IcebergSimMap::setParameters(float w, float r, float b, float i)
{
	waterline = w;
	ratioWt = r;
	beta = b;
	inWeight =i;
}
void IcebergSimMap::setGraph(Graph& g) {
	gPtr = &g;
	N = g.num_vertices();
	//cout << "ASSUME GRAPH IS UNDIRECTED" << endl;
	//gPtr->setOutEdgesToUndirectedEdges();
}
/////////////////////////////////////////////
/////////////////////////////////////////////
inline
void IcebergSimMap::insert(int x, int y, float val)
{
	if (x > y) {
		simMap[make_pair(x,y)] = val;
	}
	else {
		simMap[make_pair(y,x)] = val;
	}
}
/////////////////////////////////////////////
float IcebergSimMap::valUnsafe(int i, int j)
{
	if (i == j) {
		return 1;
	}
	SimMapType::iterator mPtr = simMap.find(make_pair(i,j));
	if ( mPtr != simMap.end() )
	{
		return mPtr->second;
	}
	float ratio = gPtr->out_degree(i) / (float) gPtr->out_degree(j);
	ratio = min (ratio, 1/ratio);
	return (1-beta)*ratio*ratioWt + beta;
}
/////////////////////////////////////////////
float IcebergSimMap::val(int x, int y)
/* This is the safe-but-slow version, where the caller does not
 * need to know anything about the underlying data struture.
 */
{
	int i, j;
	if (x == y) {
		return 1;
	}
	if (x > y) {
		i = x;
		j = y;
	}
	else {
		i = y;
		j = x;
	}
	SimMapType::iterator mPtr = simMap.find(make_pair(i,j));
	if ( mPtr != simMap.end() )
	{
		return mPtr->second;
	}
	float ratio = gPtr->out_degree(i) / (float) gPtr->out_degree(j);
	ratio = min (ratio, 1/ratio);
	return (1-beta)*ratio*ratioWt + beta;
}
/////////////////////////////////////////////
bool IcebergSimMap::Diff( IcebergSimMap& PrevIceberg, float threshold )
{
	float diff;
	float maxDiff = 0.0;
	double totDiff = 0.0;
	int jIncr = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ ) {
		for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
			diff = abs( valUnsafe(i,j) - PrevIceberg.valUnsafe(i,j) );
			totDiff += 2*diff;		// count each cell twice
			if (diff > maxDiff) {
				maxDiff = diff;
			}
		}
		// Difference of the diagonal cell will usually be 0, but not for the 1st iteration
		diff = abs( val(i,i) - PrevIceberg.val(i,i) );
		totDiff += diff;
		if (diff > maxDiff) {
			maxDiff = diff;
		}
	}
	cout << "maxDiff = " <<maxDiff<< ", totDiff = " <<totDiff<< ", threshold = " <<threshold<<endl;
	if( totDiff < threshold ) {
		return true;
	}
	else {
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////
// Change 12/29/10: Also check the diagonal for differences
bool IcebergSimMap::DiffRelative( IcebergSimMap& MatrixOdd, float threshold, int turn )
{
	// Compute the relative difference between the previous and current hash map values.
	// Also display the relative difference of the virtual full matrix, for comparison to SimMatrix
	float diff;
	float maxDiff = 0;
	double totDiff = 0;
	double Sum = 0;
	// If turn is odd, then MatrixOdd is the newer one, otherwise this matrix is the newer one.
	bool odd = (turn % 2 == 1);

	// Only check the "full" matrix if N < 10000
	if (N <= 10000) {
		int jIncr = (lowerTri) ? -1 : +1;
		for( int i = 0; i < N; i++ ) {
			for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
				diff = abs( val(i,j) - MatrixOdd.val(i,j) );
				totDiff += 2*diff;		// count each cell twice
				if (diff > maxDiff) {
					maxDiff = diff;
				}
				// count the cell of the older matrix twice
				Sum = (odd) ? (Sum + 2*val(i,j)) : (Sum + 2*MatrixOdd.val(i,j));
			}
			// Difference of the diagonal cell (should be non-zero only for the 1st iteration)
			diff = abs( val(i,i) - MatrixOdd.val(i,i) );
			totDiff += diff;
			if (diff > maxDiff) {
				maxDiff = diff;
			}		
			Sum = (odd) ? (Sum + val(i,i)) : (Sum + MatrixOdd.val(i,i));
		}
		cout << "Full matrix difference:" <<endl;
		cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;
	}

	////////////////////////////////////////
	SimMapType oddMap = MatrixOdd.getMap();
	SimMapType::iterator itCurr, itPrev;
	if (odd) {
		itCurr = oddMap.begin();
		itPrev = simMap.begin();
	}
	else {
		itCurr = simMap.begin();
		itPrev = oddMap.begin();
	}

	maxDiff = 0.0;
	totDiff = 0.0;
	Sum = 0.0;
	while (itCurr != simMap.end() && itCurr != oddMap.end()) {
		diff = abs( itCurr->second - itPrev->second );
		if (diff > maxDiff) {
			maxDiff = diff;
		}
		if (itCurr->first.first == itCurr->first.second) { // if a diagonal cell
			totDiff += diff;
			Sum += itPrev->second;
		}
		else {
			totDiff += 2 * diff;
			Sum += 2 * itPrev->second;
		}
		itCurr++;
		itPrev++;
	}
	cout << "Iceberg map difference:" <<endl;
	cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;

	if( totDiff < Sum*threshold ) {
		return true;
	}
	else {
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////
bool IcebergSimMap::DiffRelative( IcebergSimMap& prevIceberg, float threshold )
{
	// Version of DiffRelative that doesn't need to be told what 'turn' it is.
	// This iceberg is always the current one, and the prevMap parameter is the previous Map.
	// Compute the relative difference between the previous and current hash map values.
	// Also display the relative difference of the virtual full matrix, for comparison to SimMatrix
	float diff;
	float maxDiff = 0;
	double totDiff = 0;
	double Sum = 0;

	// Compute differences for the whole "virtual" matrix, if small enough
	if (N <= 5000) {
		int jIncr = (lowerTri) ? -1 : +1;
		for( int i = 0; i < N; i++ ) {
			for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
				diff = abs( val(i,j) - prevIceberg.val(i,j) );
				totDiff += 2*diff;		// count each cell twice
				if (diff > maxDiff) {
					maxDiff = diff;
				}
				// count the cell of the older matrix twice
				Sum = Sum + 2*prevIceberg.valUnsafe(i,j);
			}
			// Difference of the diagonal cell (should be non-zero only for the 1st iteration)
			diff = abs( val(i,i) - prevIceberg.val(i,i) );
			totDiff += diff;
			if (diff > maxDiff) {
				maxDiff = diff;
			}		
			Sum = Sum + prevIceberg.valUnsafe(i,i);
		}
		cout << "Full matrix difference:" <<endl;
		cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;
	}

	// Compute difference only for the iceberg
	// If we want to "melt" the newer icebergs, then we cannot assume that currMap
	// and prevMap contain the same keys!
	// Must iterate through currMap, then look-up in prevMap
	SimMapType prevMap = prevIceberg.getMap();
	SimMapType::iterator itCurr, itPrev;

	itCurr = simMap.begin();
	//itPrev = prevMap.begin();

	maxDiff = 0;
	totDiff = 0;
	Sum = 0;
	//while (itCurr != simMap.end() && itCurr != prevMap.end()) {
	while ( itCurr != simMap.end() ) {
		itPrev = prevMap.find(itCurr->first);		// new
		diff = abs( itCurr->second - itPrev->second);
		if (diff > maxDiff) {
			maxDiff = diff;
		}
		if (itCurr->first.first == itCurr->first.second) { // if a diagonal cell
			totDiff += diff;
			Sum += itPrev->second;
		}
		else {
			totDiff += 2 * diff;
			Sum += 2 * itPrev->second;
		}
		itCurr++;
		//itPrev++;
	}
	cout << "Iceberg map difference:" <<endl;
	cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;

	if( totDiff < Sum*threshold ) {
		return true;
	}
	else {
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::Print( ostream& out, int prec, bool header )
{
	float num;
	char star = '*';
	if (header) {
		for( int i = 0; i < N; i++ ) {
			out << "\t" << fixed << setprecision(prec) << i;
		}
		out << endl;
	}

	for( int i = 0; i < N; i++ ) {
		if (header) {
			out << i << ":\t";
		}
		for( int j = 0; j < N; j++ ) {
			num = val(i,j);
			if (header && simMap.find(make_pair(i,j)) != simMap.end()) {
				out << fixed << setprecision(prec) << num << star << "\t";
			}
			else {
			out << fixed << setprecision(prec) << num << "\t";
			}
		}
		out << endl;
	}
	out << "# Hash table size = " << simMap.size() << endl;
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::PrintMatlab( ostream& out, int prec, bool offset)
{
	float num;
	double maxSim = (offset) ? DENDROGRAM_BASE : 0;

	for( int i = 0; i < N; i++ ) {
		for( int j = 0; j < N; j++ ) {
			if( i==j ) {
				out << fixed << setprecision(prec) << 0.000 << "\t";
				continue;
			}
			num = val(i,j);

			if( (int)(num) == 1 ) {
				out << fixed << setprecision(prec) << maxSim << "\t";
			}
			else {
				out << fixed << setprecision(prec) << 1-num << "\t";
			}
		}
		out << endl;
	}
	out << "# Hash table size = " << simMap.size() << endl;
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::PrintCompact( ostream& out, int precision, bool compact, bool header )
{
	// To minimize bytes needed to store the matrix, the following format is used:
	// 1. Since the matrix is symmetric, only one triangle is stored.  It is stored
	// in Matlab squareform format, which is the following:
	// Matlab stores (1,2)...(1,n),(2,3)...(2,n),...,(n-1,n) (1st version)
	// = (2,1),(3,1)...(n,1),(3,2),(4,2)...(n,2),...,(n,n-1) (2nd version)
	// For the default RoleSim option, use the 1st mode.
	// For RoleSim -tri option, we should use the 2nd mode.
	// However, Matlab would expect the diagonal to be all zeroes.
	// Our matrix is a similarity matrix with all ones.  If the saved matrix is
	// restored to full form by matlab, remember to change the diagonal to ones.
	// 2. Floating point numbers are represented as having only 'p' bits to
	// the right of the decimal point, then shifted this many bits to the left
	// (multiplied by 10^p).  Why?  Because, for example, the number
	// "0.471..." becomes "471".  We eliminate the two characters "0.".

	if ( !compact ) {
		Print( out, precision, header );
	}
	else {
		int multInt = 1;
		for (int i = 0; i < precision; i++) { // Compute power of 10 multiplier
			multInt = multInt * 10;
		}
		float multF = (float)multInt;

		//  I thought this might be a faster way to cast as (int), but it's slower
		//out << fixed << setprecision(0) << noshowpoint;
		if (lowerTri) {
			for (int col = 0; col < N-1; col++) {
				for (int row = col + 1; row < N; row++) {
					out << (int)(val(row,col)*multF) << "\t";
				}
				out << "\t";
			}
		}
		else {
			for (int row = 0; row < N-1; row++) {
				for (int col = row + 1; col < N; col++) {
					out << (int)(val(row,col)*multF) << "\t";
				}
				out << "\t";
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::PrintMatlabCompact( ostream& out, int precision, bool compact ) 
{
	// See notes for PrintCompactVector
	// Differences: This is a distance matrix instead of a similarity matrix
	// If the similarity = 1 (max distance), use DENDROGRAM_BASE as the value.

	if ( !compact ) {
		PrintMatlab( out, precision, true );
	}
	else {
		int multInt = 1;
		for (int i = 0; i < precision; i++) {
			multInt = multInt * 10;
		}
		float multF = (float)multInt;
		int num;

		// I thought this might be a faster way to cast as (int), but it's slower
		//out << setprecision(0) << noshowpoint;
		if (lowerTri) {
			for (int col = 0; col < N-1; col++) {
				for (int row = col + 1; row < N; row++) {
					num = (int)( (1.0 - val(row,col))*multF );
					out << num << "\t";
				}
				out << "\t";
			}
		}
		else {
			for (int row = 0; row < N-1; row++) {
				for (int col = row + 1; col < N; col++) {
					num = (int)( (1.0 - val(row,col))*multF );
					out << num << "\t";
				}
				out << "\t";
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::PrintTable (ostream& out)
{
	SimMapType simMap = getMap();
	SimMapType::iterator mItr;
	cout << "Iceberg Sim Table" << endl;
	cout << "X \t Y \t Value" <<endl;
	for (mItr = simMap.begin(); mItr != simMap.end(); mItr++)
	{

		int x = mItr->first.first;
		int y = mItr->first.second;
		float val = mItr->second;
		cout << x << "\t" << y << "\t" << val << endl;
	}
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::save( ostream& out, char* graphFilename, int precision, bool compact )
{
	// write everything necessary to recreate the iceberg to a file.
	// 1. Write the header line
	out << "# iceberg";
	out << " graph "     << graphFilename  << " N "        << N;
	out << " waterline " << waterline      << " ratioWt "  << ratioWt;
	out << " beta "      << beta           << " inWeight " << inWeight;
	out << " #" << endl;

	// 2. Write the map in tabular format
	int multiplierI = 1;
	for (int i = 0; i < precision; i++) { // Compute power of 10 multiplier
		multiplierI = multiplierI * 10;
	}
	float multiplierF = (float)multiplierI;
	int intVal;
	SimMapType::iterator mItr;

	if (!compact)
	{
		for (mItr = simMap.begin(); mItr != simMap.end(); mItr++) {
			intVal = (int)(mItr->second * multiplierF);
			out << mItr->first.first << "\t" << mItr->first.second << "\t" << intVal << endl;
		}
	}
	else
	{
		/* Tabular format requires 3 numbers for each entry: x, y, val = 3m
		Use an even more compact form, which approaches 2m:
		Bin node-pairs by x value.  Then for each bin, display the x value once,
		the size of the y-list, followed by a list of (y,val) pairs.
		To make it EVEN MORE compact, instead of x and y values, display differential
		values (x - prevX) and (y - prevY), which are often single-digit numbers.
		*/
		// 1. Create a map of maps
		map< int, map<int,float> > xMap;
		map< int, map<int,float> >::iterator xItr;
		map<int,float> yMap;
		map<int,float>::iterator yItr;
		int x, y;
		for (mItr = simMap.begin(); mItr != simMap.end(); mItr++) {
			x = mItr->first.first;
			y = mItr->first.second;
			xItr = xMap.find(x);
			if (xItr == xMap.end()) {
				// x is not yet in the xMap
				yMap.clear();
				yMap[y] = mItr->second;
				xMap[x] = yMap;
			}
			else {
				// Add (y,val) to the existing x entry
				(xItr->second)[y] = mItr->second;
			}
		}
		int prevX, prevY;
		prevX = 0;
		// 2. Write the map of maps
		for (xItr = xMap.begin(); xItr != xMap.end(); xItr++) {
			out << (xItr->first - prevX);
			prevX = xItr->first;

			yMap = xItr->second;
			out << '\t' << xItr->second.size();
			prevY = 0;
			for (yItr = yMap.begin(); yItr != yMap.end(); yItr++) {
				intVal = (int)(yItr->second * multiplierF);
				out << '\t' << (yItr->first - prevY) << '\t' << intVal;
				prevY = yItr->first;
			}
			out << endl;
		}
	}
}
////////////////////////////////////////////////////////////////////////////
int IcebergSimMap::load( istream& in, int precision, Graph& g, bool compact )
{
	cout << "Starting IcebergSimMap.load" << endl;
	// 1. Read the header
	string line, token1, token2;
	string graphFilename;
	istringstream ss;

	getline(in, line);
	ss.clear();
	ss.str(line);
	ss >> token1 >> token2;
	if (token1.compare("#") != 0 || token2.compare("iceberg") != 0) {
		cout << "ERROR: file header does not begin with '# iceberg'" << endl;
		return -1;
	}
	ss >> token1;
	while (token1.compare("#") != 0) {
		if (token1.compare("graph") == 0) {
			ss >> graphFilename >> token1;
		}
		else if (token1.compare("N") == 0) {
			ss >> N >> token1;
		}
		else if (token1.compare("waterline") == 0) {
			ss >> waterline >> token1;
		}
		else if (token1.compare("ratioWt") == 0) {
			ss >> ratioWt >> token1;
		}
		else if (token1.compare("beta") == 0) {
			ss >> beta >> token1;
		}
		else if (token1.compare("inWeight") == 0) {
			ss >> inWeight >> token1;
		}
	}
	//setParameters(waterline, ratioWt, beta, inWeight);

	// 2. Read the map table.
	simMap.clear();
	int multiplierI = 1;
	for (int i = 0; i < precision; i++) { // Compute power of 10 multiplier
		multiplierI = multiplierI * 10;
	}
	float multiplierF = (float)multiplierI;

	if (compact) {
		/* Bin node-pairs by x value.  Then for each bin, display the x value once,
		the size of the y-list, followed by a list of (y,val) pairs.
		To make it EVEN MORE compact, instead of x and y values, display differential
		values (x - prevX) and (y - prevY), which are often single-digit numbers.
		*/
		// 1. Read the file
		int diffX, prevX, x;
		int diffY, prevY, y;
		int size, val;
		prevX = 0;
		while (getline(in,line)) {
			ss.clear();
			ss.str("");
			ss.str(line);
			ss >> diffX;
			x = diffX + prevX;
			prevX = x;

			ss >> size;
			prevY = 0;
			for (int n = 0; n < size; n++) {
				ss >> diffY >> val;
				y = diffY + prevY;
				prevY = y;
				simMap[make_pair(x,y)] = val/multiplierF;
			}
		}
	}
	else {
		// file is in regular tabular format: x y value
		int x, y, val;
		while (getline(in,line)) {
			ss.clear();
			ss.str(line);
			ss >> x >> y >> val;
			simMap[make_pair(x,y)] = (val/multiplierF);
		}
	}
	int mapSize = simMap.size();
	long long graphTriMatrixSize = g.num_vertices() * (g.num_vertices() - 1 )/2;
	cout << "Read in " << mapSize << " map values,";
	cout << " which is " << mapSize*100.0/graphTriMatrixSize << "% of the full similarity matrix" << endl;

	// 3. Set the graph
	setGraph(g);
	cout << endl << "Loaded iceberg:" << endl;
	cout << "Graphfile " << graphFilename << " with " ;
		cout << gPtr->num_vertices() << " vertices and " << gPtr->num_edges() << " edges" << endl;
	cout << " waterline = " << waterline;
	cout << " ratioWt = " << ratioWt;
	cout << " beta = " << beta;
	cout << " inWeight = " << inWeight;
	cout << endl << "Finished IcebergSimMap.load" << endl;
	return 0;		// successful completion
}
////////////////////////////////////////////////////////////////////////////
string IcebergSimMap::timeToString(double time) {
	ostringstream ss;
	int hr = (int)(time / 3600.0);
	double remainingSec = time - hr*3600.0;
	int min = (int)(remainingSec / 60.0);
	string minSpace = (min < 10) ? ":0" : ":";
	double sec = remainingSec - min*60.0;
	string secSpace = (sec < 10) ? ":0" : ":";
	ss << hr << minSpace << min << secSpace << sec << " (" << time << " sec)";
	return ss.str();
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::displayElapsedTime(time_t& previous, ostream& out, bool newline) {
	time_t current = time(NULL);
	double diff = difftime(current, previous);
	out << "elapsed time = " << timeToString(diff);
	if (newline) out << endl;
	previous = current;
}
////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::InitializeRatioDepth1( Graph& g, bool shortcutMatching )
{
#ifdef _TIMEVAL
		gettimeofday(&before_time, NULL);
#else
	clock_t start_time(clock());	
#endif
	time_t previous = time(NULL);

	// Group the nodes in buckets according to degree
	VectVectIntType degBucket;		// size maxDegree+1, for each degree, a list of nodes
	VectVectIntType neighDegrees;	// size N,           for each node, a list of its neighbor's degrees
	
	// Find out the maxDegree of the graph
	int maxDegree = 0;
	for (int i = 0; i < N; i++)
	{
		maxDegree = (g[i].degree > maxDegree) ? g[i].degree : maxDegree;
	}
	cout << "1. Found maxDegree" << endl;
	displayElapsedTime(previous, cout, true);
	// Create a vector of buckets, of length maxDegree + 1, counting bucket 0
	for ( int d = 0; d <= maxDegree; d++ )
	{
		degBucket.push_back(VectIntType());
	}

	// For each node, make a "signature" list of the degrees of its neighbors
	for ( int i = 0; i < N; i++ )
	{
		VectIntType temp;
		EdgeList::iterator jItr;
		degBucket[g[i].degree].push_back(i);	// Add i to the appropriate bucket
		temp.clear();							// Build sorted neighbor degree lists
		for (jItr = g[i].edgeList.begin(); jItr != g[i].edgeList.end(); jItr++)
		{
			temp.push_back( g[*jItr].degree ); // degree of the j_th neighbor of i
		}
		sort( temp.begin(), temp.end() );	// sort degrees to create the signature
		neighDegrees.push_back( temp );
	}
	cout << "2. Made degree signatures of each node" << endl;
	displayElapsedTime(previous, cout, true);

	// Verify degree buckets: print the size of each bucket:
	cout << "Degree D: Number of Nodes" << endl;
	for ( int d = 0; d <= maxDegree; d++ )
	{
		printf("%6d:%9d", d, degBucket[d].size());
		if (d % 5 == 4) cout << endl;
		//cout << d << ": " << degBucket[d].size() << endl;
	}
	cout << endl;

	// For each pair of vertices u and v, compare the i_th item from u's degreeList Du_i
	// with the j_th item from v's degreeList, Dv_j.
	VectIntType::iterator uItr, vItr;
	VectIntType D_u, D_v;
	VectIntType::const_iterator diItr, djItr;
	float theta = (waterline - beta) / (1 - beta);
	cout << "theta = " << theta << endl;
	float oneMinusBeta = 1 - beta;
	MatchItem matchItem;
	int fullAbortCount = 0;
	float weight;
	int insertCount = 0;
	int abortCount = 0;

	//// For the degree-1 nodes, use a simplified version
	//cout << "d1=" << 1 << ",size " << degBucket[1].size();
	//cout << "...d2=" << 1 << ",size" << degBucket[1].size() << endl;
	//for (uItr = degBucket[1].begin(); uItr != degBucket[1].end(); uItr++)
	//{
	//	D_u = neighDegrees[*uItr];
	//	for (vItr = degBucket[1].begin(); vItr != degBucket[1].end() && vItr != uItr; vItr++)
	//	{
	//		D_v = neighDegrees[*vItr];

	//		if (D_u[0] > D_v[0])
	//		{
	//			weight = oneMinusBeta * D_v[0] / (float) D_u[0] + beta;
	//		}
	//		else
	//		{
	//			weight = oneMinusBeta * D_u[0] / (float) D_v[0] + beta;
	//		}
	//		if (weight > theta)
	//		{
	//			insert(*uItr, *vItr, oneMinusBeta*weight + beta);
	//			insertCount++;
	//		}
	//	}
	//}
	//displayElapsedTime(previous, cout, false);
	//cout << " #inserts = " << insertCount << endl;


	// for each degree-bin d1 & degree-bin d2 < d1 (but not too much less; must be theta-similar)
	for (int d1 = 1; d1 <= maxDegree; d1++)
	{
		cout << "d1=" << d1 << ",size " << degBucket[d1].size();
		for (int d2 = (int)(theta*d1 + 0.5); d2 <= d1; d2++)
		{
			cout << "...d2=" << d2 << ",size" << degBucket[d2].size() << endl;
			abortCount = 0;
			insertCount = 0;


			 // For each node u of degree d1
			for (uItr = degBucket[d1].begin(); uItr != degBucket[d1].end(); uItr++)
			{
				D_u = neighDegrees[*uItr];

				// For each node v such that (degree d2 <= d1) and (if d2 = d1, then v < u)
				// If d2 < d1, then uItr and vItr point to different lists.  Go through the whole list.
				// If d2 = d1, vItr should stop when it reaches uItr
				// Since these lists are sorted, this will insure that ever node pairing is included once (u > v)
				for (vItr = degBucket[d2].begin(); vItr != degBucket[d2].end() && vItr != uItr; vItr++)
				{
					D_v = neighDegrees[*vItr];
					if (verbose) cout << "(" <<*uItr<< "," <<*vItr<< ")";
					// Iterate through the neighbor degree lists of u and v
					// WRONG: If degrees of neighbors are theta-similar, add their best-case similarity
					// RIGHT: Emulate regular RoleSim iteration, but instead of using the previous "matrix",
					//	compute the prev. similarity S(x,y) = (1-beta)d(x)/d(y) + beta, where d(x)<d(y)

					//////// Start of Update RoleSim
					priority_queue<MatchItem, vector<MatchItem>, MatchLessThan> matchQueue;
					int size_i = D_u.size();
					int size_j = D_v.size();
					matchItem.x = 0;
					bool firstPair = true;
					bool abortMatching = false;
					for (diItr = D_u.begin(); diItr != D_u.end(); diItr++)

					{
						matchItem.y = 0;
						for (djItr = D_v.begin(); djItr != D_v.end(); djItr++)
						{
							// Find the weight of the neighbor pair, depend on the case
							if (*diItr == *djItr)
							{
								matchItem.weight = 1.0;
							}
							else if (*diItr > *djItr)
							{
								matchItem.weight = oneMinusBeta * (*djItr / (float) *diItr) + beta;
							}
							else
							{
								matchItem.weight = oneMinusBeta * (*diItr / (float) *djItr) + beta;
							}

							//Test if it is possible for the remaining pairs to have enough weight
							if (firstPair && shortcutMatching)
							{
								firstPair = false;
								if ((*diItr > *djItr) && theta*d1 - matchItem.weight + 1 > d2)
								{
									abortMatching = true; // right
									abortCount++;
								}
							}
							if (abortMatching) break;

							matchQueue.push(matchItem);
							matchItem.y += 1;
						}
						if (abortMatching) break;

						matchItem.x += 1;
					}

					if (abortMatching)
					{
						if (verbose && N < 25) cout << "Abort matching" << endl;
						abortMatching = false;
					}
					else {
						// Find the top weighted match-pairs
						VectIntType Flag1 = VectIntType( size_i, FALSE );
						VectIntType Flag2 = VectIntType( size_j, FALSE );
						float numerator = 0;
						int matchSize = min(size_i, size_j);
						int matchCount = 0;

						while (!matchQueue.empty() && matchCount < matchSize)
						{
							matchItem = matchQueue.top();
							if ( Flag1[matchItem.x]==FALSE && Flag2[matchItem.y]==FALSE )
							{
								Flag1[matchItem.x] = TRUE;
								Flag2[matchItem.y]= TRUE;
								numerator += matchItem.weight;
								matchCount++;
							}
							matchQueue.pop();
						}
						weight = numerator / d1;  // d1 >= d2 guaranteed
						if (weight >= theta) 
						{
							if (verbose && N < 25)
							{
								cout << ": weight = " <<numerator<< " / " <<d1<< " = " <<weight<< " *" << endl;
							}
							insert(*uItr, *vItr, oneMinusBeta*weight + beta);
							insertCount++;
						}
						else if (verbose && N < 25)
						{
							cout << ": weight = " <<numerator<< " / " <<d1<< " = " <<weight<< endl;
						}
					}

				} // vItr
				//cout << endl;
			} // uItr
			if (abortCount > 0)
			{
				cout << abortCount << " ABORTS, ";
				fullAbortCount += abortCount;
			}
			displayElapsedTime(previous, cout, false);
			cout << " #inserts = " << insertCount << endl;
		} //d2
		//cout << endl;
		//cout << "simMap.size() = " << simMap.size() << endl;
	} // d1
	cout << "FULL Abort Count = " << fullAbortCount << endl;

	int hashMapSize = simMap.size();
	int matrixSize = N*(N-1)/2;
	double sizePct = 100.0 * hashMapSize / matrixSize;

	cout << "Sim hash table size " << hashMapSize << " = " << sizePct << "% of full tri. array size " << matrixSize << endl;
	if (verbose)
	{
		PrintTable(cout);
	}
		
}

////////////////////////////////////////////////////////////////////////////
void IcebergSimMap::updateRoleSim( IcebergSimMap& preIceberg, float beta, bool divMax, bool melt)
{
	float ratio;
	SimMapType::iterator mPtr;
	SimMapType prevMap = preIceberg.getMap();
	float oneMinusBeta = 1 - beta;
	float betaRatioWt = ratioWt * oneMinusBeta;

	MatchItem matchItem;

	// For each vertex pair in the map, compute its next RoleSim value and store it in the returned map.
	simMap.clear();		// CLEAR the old map!!!
	EdgeList N_i, N_j;
	SimMapType::iterator mItr;
	for (mItr = prevMap.begin(); mItr != prevMap.end(); mItr++)
	{
		int i = mItr->first.first;	// SimMatType is map< pair<int,int>, float >
		int j = mItr->first.second;	//
		N_i = gPtr->out_edges(i);
		N_j = gPtr->out_edges(j);

		priority_queue<MatchItem, vector<MatchItem>, MatchLessThan> matchQueue;

		int size_i = N_i.size();
		int size_j = N_j.size();
		if( size_i == 0 || size_j == 0 )
		{
			continue;
		}

		VectIntType::const_iterator iItr, jItr;
		matchItem.x = 0;
		for (iItr = N_i.begin(); iItr != N_i.end(); iItr++)
		{
			matchItem.y = 0;
			for (jItr = N_j.begin(); jItr != N_j.end(); jItr++)
			{
				// Find the weight of the neighbor pair, depend on the case
				// Case 1: the neighbors are the same node
				if (*iItr == *jItr)
				{
					matchItem.weight = 1.0;
				}
				else
				{
					// Else, Check the hash map
					if (*iItr > *jItr)
					{
						mPtr = prevMap.find(make_pair(*iItr, *jItr));
					}
					else
					{
						mPtr = prevMap.find(make_pair(*jItr, *iItr));
					}

					// Case 2. The pair is in the hash map
					if ( mPtr != prevMap.end() )
					{
						matchItem.weight = mPtr->second;
					}
					// Case 3. Not in the hash map; compute a weight on the fly
					else
					{
						int xDegree = gPtr->out_degree(*iItr);
						int yDegree = gPtr->out_degree(*jItr);
						if (xDegree < yDegree)
						{
							ratio = xDegree / (float) yDegree;
						}
						else
						{
							ratio = yDegree / (float) xDegree;
						}
						matchItem.weight = ratio*betaRatioWt + beta;
					}
				}
				matchQueue.push(matchItem);
				matchItem.y += 1;
			}
			matchItem.x += 1;
		}

		// Find the top weighted match-pairs
		VectIntType Flag1 = VectIntType( size_i, FALSE );
		VectIntType Flag2 = VectIntType( size_j, FALSE );
		float numerator = 0;
		int matchSize = min(size_i, size_j);
		int matchCount = 0;

		while (!matchQueue.empty() && matchCount < matchSize)
		{
			matchItem = matchQueue.top();
			if ( Flag1[matchItem.x]==FALSE && Flag2[matchItem.y]==FALSE )
			{
				Flag1[matchItem.x] = TRUE;
				Flag2[matchItem.y]= TRUE;
				numerator += matchItem.weight;
				matchCount++;
			}
			matchQueue.pop();
		}

		float denominatorF = (float)max(size_i, size_j);
		//float denominator = sqrt ( (float)(size_i * size_j) ); // OLD def'n of denominator

		if (melt)
		{
			float defaultWeight = ratioWt * matchSize;
			//Insert only if new weight > default non-iceberg weight
			if (numerator < defaultWeight)
			{
				cout <<"mapSize="<< simMap.size() <<", erasing ("<<i<<","<<j<<")";
				cout << ", value " << numerator/denominatorF << " too low" << endl;
				//simMap.erase(make_pair(i,j));
			}
			else
			{
				insert(i, j, oneMinusBeta*numerator/denominatorF + beta);
			}
		}
		else
		{
			//assert(weight<1.1);
			insert(i, j, oneMinusBeta*numerator/denominatorF + beta);
		}
	}
	cout << "Size of new iceberg = " << simMap.size() << endl;
}
/////////////////////////////////////////////
