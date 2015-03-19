// SimMatrix: a class for a node-node similarity matrix, with functions for RoleSim, SimRank, and related algorithms.
#include "SimMatrix.h"

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
#include "PM/PerfectMatching.h" // 12/21/11: moved to PM subdirectory
#include "PM/PMimplementation.h" // 12/21/11: moved to PM subdirectory

/* Conceptually, a SimMatrix is a square NxN matrix of floating point numbers.
* Since the values are node similarity values, where sim(x,y) = sim(y,x),
* it is actually a triangular matrix, implemented as a vector (size N)
* of variable-length float vectors. This should be invisible to the user, however.
* Originally, this could be configured as either an upper-triangular or
* lower-triangular matrix, but now it is only lowerTri.
* The outer vector represents the rows, numbered from 0 to n-1.
* Row r has length (r+1).

The SimMatrix class also includes several node-similarity measures an algorithms

There are two formats for printing or displaying a SimMatrix. In each case, the
maximum number bits of precision can be specified.
(1) Standard rows and columns, one row per line, and columns separated by tabs.
Obviously, this is only appropriate for small matrices.
(2) Linear sequence of values, separated by spaces, also called a simVec.
Also, each value is a string of digit characters, so that the memory requirements
scale with the precision. This is used when the -compact option is specified.
The sequence of values is documented with the PrintCompact() function 
*/
SimMatrix::SimMatrix(void)
{
	N = 0;
	lowerTri = true;
}
////////////////////////////////////////////////////////////////////////////
/** Create a size N matrix.  Because the matrix is symmetric, we create only the lower triangle */
SimMatrix::SimMatrix(int size)
{
	N = size;
	lowerTri = true;
	VectFloatType temp;
	for ( int i = 0; i < N; i++ ) {
		temp = VectFloatType(i+1, 0.0);
		Matrix.push_back(temp);
	}
}
////////////////////////////////////////////////////////////////////////////
SimMatrix::~SimMatrix(void) { }
////////////////////////////////////////////////////////////////////////////
/* Given any "virtual" location in the matrix, return the "physical" value,
 * based on the assumption that only the lower triangle is actually implemented.
 * This function separates the implementation details (only the lower triangle
 * physically exists) from the user. */
inline
float SimMatrix::val(int i, int j)
{
	return (i > j) ? Matrix[i][j] : Matrix[j][i];
}
////////////////////////////////////////////////////////////////////////////
SimMatrixType& SimMatrix::matrix() { return Matrix; }
////////////////////////////////////////////////////////////////////////////
/* Implementation of set() assumes the matrix is physically implemented as
* a symmetric lower triangular matrix, but this is invisible to the user */
void SimMatrix::set(int i, int j, float val)
{
	if (i > j) {
		Matrix[i][j] = val;
	}
	else {
		Matrix[j][i] = val;
	}
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::Clear()
{
	VectFloatType::iterator cell;
	for( int i = 0; i < N; i++ ) {
		for (cell = Matrix[i].begin(); cell != Matrix[i].end(); cell++) {
			*cell = 0.0;
		}
	}
}
////////////////////////////////////////////////////////////////////////////
bool SimMatrix::Diff( SimMatrix& Matrix2, float threshold )
{
	float diff;
	float maxDiff = 0;
	double totDiff = 0;
	int jIncr = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ ) {
		for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
			diff = abs( Matrix[i][j] - Matrix2.val(i,j) );
			totDiff += 2*diff;		// count each cell twice
			if (diff > maxDiff) {
				maxDiff = diff;
			}
		}
		// Difference of the diagonal cell (should be non-zero only for the 1st iteration)
		diff = abs( Matrix[i][i] - Matrix2.val(i,i) );
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
bool SimMatrix::DiffRelative( SimMatrix& Matrix2, float threshold, int turn )
{
	float diff;
	float maxDiff = 0;
	double totDiff = 0;
	double Sum = 0;
	// If turn is odd, then Matrix2 is the newer one, otherwise Matrix1 is the newer one.
	bool odd = (turn % 2 == 1);

	int jIncr = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ ) {
		for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
			diff = abs( Matrix[i][j] - Matrix2.val(i,j) );
			totDiff += 2*diff;		// count each cell twice
			if (diff > maxDiff) {
				maxDiff = diff;
			}
			// count the cell of the older matrix twice
			Sum = (odd) ? (Sum + 2*Matrix[i][j]) : (Sum + 2*Matrix2.val(i,j));
		}
		// Difference of the diagonal cell (should be non-zero only for the 1st iteration)
		diff = abs( Matrix[i][i] - Matrix2.val(i,i) );
		totDiff += diff;
		if (diff > maxDiff) {
			maxDiff = diff;
		}		
		Sum = (odd) ? (Sum + Matrix[i][i]) : (Sum + Matrix2.val(i,i));
	}

	cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;
	if( totDiff < Sum*threshold ) {
		return true;
	}
	else {
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////
bool SimMatrix::DiffRelative( SimMatrix& previous, float threshold )
{
	// New version of DiffRelative that doesn't need the 'turn' parameter to know if this is an even or odd iteration.
	// This matrix is always the current one, and the matrix passes as a parameter is the previous one.
	float diff;
	float maxDiff = 0;
	double totDiff = 0;
	double Sum = 0;

	SimMatrixType prevMat = previous.matrix();
	lowerTri = true;
	int jIncr = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ ) {
		for( int j = i + jIncr; j >= 0 && j < N; j += jIncr ) { // tri, excluding diagonal
			diff = abs( Matrix[i][j] - prevMat[i][j] );
			totDiff += 2*diff;		// count each cell twice
			if (diff > maxDiff) {
				maxDiff = diff;
			}
			// count the cell of the older matrix twice
			Sum = Sum + 2*prevMat[i][j];
		}
		// Difference of the diagonal cell (should be non-zero only for the 1st iteration)
		diff = abs( Matrix[i][i] - prevMat[i][i] );
		totDiff += diff;
		if (diff > maxDiff) {
			maxDiff = diff;
		}		
		Sum = Sum + prevMat[i][i];
	}

	cout << "maxDiff = " << maxDiff << ", totDiff = " << totDiff << ", Sum*threshold = " << Sum*threshold << endl;
	if( totDiff < Sum*threshold ) {
		return true;
	}
	else {
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::Print( ostream& out, int prec, bool header )
{
	float num;
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
			out << fixed << setprecision(prec) << num << "\t";
		}
		out << endl;
	}
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::PrintMatlab( ostream& out, int prec, bool offset)
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
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::PrintCompact( ostream& out, int precision, bool compact, bool header )
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
					out << (int)(Matrix[row][col]*multF) << "\t";
				}
				out << endl;
			}
		}
		else {
			for (int row = 0; row < N-1; row++) {
				for (int col = row + 1; col < N; col++) {
					out << (int)(Matrix[row][col]*multF) << "\t";
				}
				out << endl;
			}
		}
	}

}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::PrintMatlabCompact( ostream& out, int precision, bool compact )
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
					num = (int)( (1.0 - Matrix[row][col])*multF );
					out << num << "\t";
				}
				out << "\t";
			}
		}
		else {
			for (int row = 0; row < N-1; row++) {
				for (int col = row + 1; col < N; col++) {
					num = (int)( (1.0 - Matrix[row][col])*multF );
					out << num << "\t";
				}
				out << "\t";
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::loadSimVec ( char* infile, int N, int precision, bool distance )
{
	cout << "start LoadSimVec" << endl;
	ifstream in(infile);
	float rangeD = (float) pow(10.0, (double)precision);
	long long vecLen = N * (N - 1) / 2;

	int x, y;
	int value;
	if (distance) {
		for (y = 0; y < N; y++) {
			for (x = y+1; x < N; x++) {
				if (!in.eof() && in >> value) {
					Matrix[x][y] = (rangeD - value)/rangeD;
				}
			}
		}
	}
	else {
		for (y = 0; y < N; y++) {
			for (x = y+1; x < N; x++) {
				if (!in.eof() && in >> value) {
					Matrix[x][y] = value/rangeD;
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
void SimMatrix::InitializeUniform( ){
	VectIntType temp;

	VectVectFloatType::iterator row;
	VectFloatType::iterator cell;
	for (row = Matrix.begin(); row != Matrix.end(); row++) {
		for (cell = row->begin(); cell != row->end(); cell++) {
			*cell = 1;
		}
	}
}
/////////////////////////////////////////////////////////////////////////

void SimMatrix::InitializeDiagonal( ){

	//int N = vl.size();
	for( int u = 0; u < N; u++ )
	{
		Matrix[u][u] = 1.0;
	}
}

/////////////////////////////////////////////////////////////////////////

void SimMatrix::InitializeDepthZero( Graph& g, bool directed )
{
	//int N = num_vertices();
	int incr = (lowerTri) ? -1 : +1;
	if (directed) {
		int Iu, Iv, Ou, Ov;
		for( int u = 0; u < N; u++ )
		{
			Matrix[u][u] = 1.0;			// take care of the diagonal
			Iu = g.in_degree(u);
			Ou = g.out_degree(u);
			for( int v = u + incr; v >= 0 && v < N; v += incr )// triangle, excluding diagonal
			{
				Iv = g.in_degree(v);
				Ov = g.out_degree(v);
				if (Iu == Iv && Ou == Ov)
				{
					Matrix[u][v] = 1.0;
				}
			}
		}
	}
	else
	{
		int Nu, Nv;
		for( int u = 0; u < N; u++ )
		{
			Matrix[u][u] = 1.0;			// take care of the diagonal
			Nu = g.out_degree(u);
			for( int v = u + incr; v >= 0 && v < N; v += incr ) // triangle, excluding diagonal
			{
				Nv = g.out_degree(v);
				if (Nu == Nv)
				{
					Matrix[u][v] = 1.0;
				}

			}
		}
	}
}

// Note: This function treats edges as undirected.
// It actually sets the values of vl[i].edgeList = union (inList, outList).
// This is not fatal if the graph's edges are directed.  It just means that the initialization
// will be optimistic with respect to automorphism.
void SimMatrix::InitializeDepthOne( Graph& g, bool directed )
{
	//bool lowerTri = true;
	// For each vertex, make a "Neighbor Signature" = sorted list of the degrees of its neighbors
	// Map vertices with the same Neighbor Signature into the same bucket
	MapVectIntToVectInt DegreeMap;
	VectIntType temp;
	for( int i = 0; i < N; i++ )
	{
		temp.clear();
		for( int j = 0; j < g[i].degree; j++ )
		{
			temp.push_back( g[ g[i].edgeList[j] ].degree ); // degree of the j_th neighbor of i
		}
		sort( temp.begin(), temp.end() );			// sort degrees to create the signature
		DegreeMap[ temp ].push_back( i );
	}

	// For each DegreeMap bucket, for each pair of members, initialize similarity to 1.
	for( MapVectIntToVectInt::iterator it = DegreeMap.begin(); it != DegreeMap.end(); it++ )
	{
		for (VectIntType::iterator p1 = (it->second).begin(); p1 != (it->second).end(); p1++)
		{
			for (VectIntType::iterator p2 = (it->second).begin(); p2 != (it->second).end(); p2++)
			{
				if ((*p1 > *p2 && lowerTri) || (*p1 < *p2 && !lowerTri)) {
					Matrix[*p1][*p2] = 1.0;
				}
				else {
					Matrix[*p2][*p1] = 1.0;
				}
			}
		}
	}
	//cout << "ccc" << endl;
}
///////////////////////////////////////////////////////////////////////////////
void SimMatrix::InitializeRatioDepth0( Graph& g, float beta )
{
	// Initialize Matrix(i,j) = degree(i)/degree(j), where 
	int iDegree, jDegree;
	float ratio;
	float oneMinusBeta = 1 - beta;
	int jIncrement = (lowerTri) ? -1 : +1;
	for ( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0;			// self-similarity
		iDegree = g.out_degree(i);
		for ( int j = i + jIncrement; 0 <= j && j < N; j += jIncrement )
		{
			jDegree = g.out_degree(j);

			if (iDegree == 0 || jDegree == 0) {
				ratio = 0.0;
			}
			else if (iDegree < jDegree) {
				ratio = oneMinusBeta * iDegree / (float) jDegree + beta;
			}
			else {
				ratio = oneMinusBeta * jDegree / (float) iDegree + beta;
			}
			Matrix[i][j] = ratio;
		}
	}
}
////////////////////////////////////////////////////////////////////////
int SimMatrix::simCode(char* simMeasure)
{
	if (strcmp(simMeasure,"rs") == 0) {
		return ROLESIM;
	}
	if (strcmp(simMeasure,"sr") == 0) {
		return SIMRANK;
	}
	if (strcmp(simMeasure,"srpp") == 0) {
		return SIMRANK_PP;
	}
	if (strcmp(simMeasure,"psr") == 0) {
		return P_SIMRANK;
	}
	if (strcmp(simMeasure,"bc") == 0) {
		return BIBLIO_COUPLING;
	}
	if (strcmp(simMeasure,"cc") == 0) {
		return CO_CITATION;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////
////    SimRank														////
////////////////////////////////////////////////////////////////////////

void SimMatrix::UpdateSimRankUndirected( Graph& g, SimMatrix& PreMatrix, float beta )
{
	//int N = g.num_vertices();
	//float outWeight = 1 - inWeight;

	// Tricky way to iterate over either lower triangle columns or upper triangle columns:
	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0;
		int size1 = g[i].degree;
		for( int j = i + jIncrement; 0 <= j && j < N; j = j + jIncrement )
		{
			int size2 = g[j].degree;
			if( size1 == 0 || size2 == 0 )
			{
				Matrix[i][j] = 0.0;
				continue;
			}

			float numerator = 0.0;
			EdgeList::iterator N_i, N_j;
			EdgeList& iList = g[i].edgeList;
			EdgeList& jList = g[j].edgeList;
			for (N_i = iList.begin(); N_i != iList.end(); N_i++)
			{
				for (N_j = jList.begin(); N_j != jList.end(); N_j++)
				{
					if (*N_i > *N_j)
					{
						//numerator += (g[i].weight)*(g[j].weight)*PreMatrix.val(*N_i, *N_j);
						numerator += PreMatrix.val(*N_i, *N_j);
					}
					else {
						//numerator += (g[i].weight)*(g[j].weight)*PreMatrix.val(*N_j, *N_i);
						numerator += PreMatrix.val(*N_j, *N_i);
					}
					//int xCoord = (*N_i > *N_j) ? *N_i : *N_j;
					//int yCoord = (*N_i < *N_j) ? *N_i : *N_j;
					//numerator += (g[i].weight)*(g[j].weight)*PreMatrix.val(xCoord, yCoord);
				}
			}
			int denominator = size1 * size2;
			float weight_float = numerator / (float)denominator;
			assert(weight_float <= 1.05);
			Matrix[i][j] = (1-beta)*weight_float;
		}
	}
}

////////////////////////////////////////////////////////////////////////

void SimMatrix::UpdateSimRank( Graph& g, SimMatrix& PreMatrix, float beta, float inWeight )
// Lower tri is hardcoded in the heart (neighbor comparison)
{
	if (!g.isDirected()) {
		UpdateSimRankUndirected(g, PreMatrix, beta);
		return;
	}
	int sizeI_i, sizeI_j, sizeO_i, sizeO_j;
	float sim_I, sim_O;
	//int N = g.num_vertices();
	float outWeight = 1 - inWeight;

	// Tricky way to iterate over either lower triangle columns or upper triangle columns:
	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0;
		sizeI_i = g.in_degree(i);
		sizeO_i = g.out_degree(i);
		for( int j = i + jIncrement; 0 <= j && j < N; j = j + jIncrement )
		{
			sizeI_j = g.in_degree(j);
			sizeO_j = g.out_degree(j);

			// In neighbors
			if( sizeI_i == 0 || sizeI_j == 0 )
			{
				sim_I = 0.0;
			}
			else
			{
				float numerator = 0.0;
				EdgeList::iterator N_i, N_j;
				for (N_i = g.in_edges(i).begin(); N_i != g.in_edges(i).end(); N_i++)
				{
					for (N_j = g.in_edges(j).begin(); N_j != g.in_edges(j).end(); N_j++)
					{
						int xCoord = (*N_i > *N_j) ? *N_i : *N_j;
						int yCoord = (*N_i < *N_j) ? *N_i : *N_j;
						numerator += PreMatrix.val(xCoord, yCoord);
					}
				}
				sim_I = numerator / (sizeI_i * sizeI_j);
			}

			// Out neighbors
			if( sizeO_i == 0 || sizeO_j == 0 )
			{
				sim_O = 0.0;
			}
			else
			{
				float numerator = 0.0;
				EdgeList::iterator N_i, N_j;
				for (N_i = g.out_edges(i).begin(); N_i != g.out_edges(i).end(); N_i++)
				{
					for (N_j = g.out_edges(j).begin(); N_j != g.out_edges(j).end(); N_j++)
					{
						int xCoord = (*N_i > *N_j) ? *N_i : *N_j;
						int yCoord = (*N_i < *N_j) ? *N_i : *N_j;
						numerator += PreMatrix.val(xCoord, yCoord);
					}
				}
				sim_O = numerator / (sizeO_i * sizeO_j);
			}
			Matrix[i][j] = (1-beta)*(inWeight*sim_I + outWeight*sim_O);
		}
	}
}
////////////////////////////////////////////////////////////////////////
////    RoleSim         										    ////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Note: This function works only for undirected graphs
void SimMatrix::UpdateRoleSimExact( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax )
{
	//construct M_1
	//int N = g.num_vertices();
	VectIntType Flag1, Flag2;
	vector<pair<int, int> > Matching;
	vector<pair<int, int> > ExactMatching;
	pair<int, int> TempMatching;
	vector<float> Weight;
	//float TempWeight;
	int start, end;
	int weight_int;
	float weight_float;
	map<int, int> MapIntToInt;
	float numerator;
	VectIntType N_u, N_v;
	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		//cout << "i:" << i;
		Matrix[i][i] = 1.0;

		//for( int j = i; j < N; j++ )
		for( int j = i + jIncrement; 0 <= j && j < N; j += jIncrement )
		{
			int size1 = g[i].degree;
			int size2 = g[j].degree;
			//cout << ", j:" << j << ", size1: " << size1 << ", size2: " << size2 << endl;

			if( size1 == 0 || size2 == 0 )
			{
				Matrix[i][j] = 0.0;
				continue;
			}

			int node_num = 0;
			int edge_num = 0;
			vector<int> Edges_Vector;
			vector<int> Weights_Vector;

			N_u = g[i].edgeList;
			N_v = g[j].edgeList;
			int m, n;

			if( size1 >= size2 )//all nodes on the left are real nodes, and add necessary virtual nodes on the right
			{
				node_num = 2*size1;
				for( m = 0; m < size1; m++ )
				{//cout << "m: " << m << endl;
					for( n = 0; n < size2; n++ )
					{
						start = (g[i].edgeList[m]<g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];
						end   = (g[i].edgeList[m]>g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];
						weight_int = int(PreMatrix.val(start, end)*SCALE_FOR_EXACTMATCHING);
						//cout << "n: " << n << ", start: " << start << ", end:" << end << ", weight: " << fixed << setprecision(4) << weight_int << endl;
						//if( weight_int > 0 )
						//{
							//cout << "1. push " << m << ", " << n << endl;
							Edges_Vector.push_back( m );
							Edges_Vector.push_back( size1+n );
							Weights_Vector.push_back( 0-weight_int );
						//}
					}

					if( n==size2 )
						while( n < size1 )
						{
							//cout << "2. push " << m << ", " << n << endl;
							Edges_Vector.push_back( m );
							Edges_Vector.push_back( size1+n );
							Weights_Vector.push_back( SCALE_FOR_EXACTMATCHING );
							n++;
						}
				}
			}
			//all nodes on the right are real nodes, and add necessary virtual nodes on the left
			else
			{
				node_num = 2*size2;
				for( m = 0; m < size1; m++ )
				{
					for( n = 0; n < size2; n++ )
					{
						start = (g[i].edgeList[m]<g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];
						end   = (g[i].edgeList[m]>g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];
						weight_int = int(PreMatrix.val(start, end)*SCALE_FOR_EXACTMATCHING);
						//if( weight_int > 0 )
						//{
							//cout << "3. push " << m << ", " << n << endl;
							Edges_Vector.push_back( m );
							Edges_Vector.push_back( size2+n );
							Weights_Vector.push_back( 0-weight_int );
						//}
					}
				}

				if( m==size1 )
					while( m < size2 )
					{
						for( n = 0; n < size2; n++ )
						{
							//cout << "4. push " << m << ", " << n << endl;
							Edges_Vector.push_back( m );
							Edges_Vector.push_back( size2+n );
							Weights_Vector.push_back( SCALE_FOR_EXACTMATCHING );
						}
						m++;
					}
			}

			edge_num = Weights_Vector.size();
			assert( Weights_Vector.size()*2 == Edges_Vector.size() );
			int* EdgesList = new int[Edges_Vector.size()];
			int* WeightsList = new int[Weights_Vector.size()];
			//cout << "\tedge: " << edge_num << ", node: " << node_num << endl;
			//convert to use PerfectMatching Procedure
			for( m = 0; m < Weights_Vector.size(); m++ )
			{
				EdgesList[2*m] = Edges_Vector[2*m];
				EdgesList[2*m+1] = Edges_Vector[2*m+1];
				WeightsList[m] = Weights_Vector[m];

				//cout << endl << "\t" << EdgesList[2*m] << ", " << EdgesList[2*m+1] << ", " << WeightsList[m];
			}

			//cout << "aaa" << endl;

			//exact matching with blossom method
			struct PerfectMatching::Options options;
			PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
			pm->options = options;
			//cout << "jumpstart: " << pm->options.fractional_jumpstart << endl;
			//cout << "bbb" << endl;
			for (int e=0; e<edge_num; e++) pm->AddEdge(EdgesList[2*e], EdgesList[2*e+1], WeightsList[e]);
			pm->options.verbose = false;
			pm->Solve();
			//cout << "ccc" << endl;
			ComputePerfectMatchingCost(node_num, edge_num, EdgesList, WeightsList, pm);
			//cout << "ddd" << endl;
			numerator = 0.0;
			for (m=0; m<node_num; m++)
			{
				n = pm->GetMatch(m);
				if( m < n )
				{
					//cout << "Candidate: " << m << ", " << n << endl;
					if( size1>=size2 && n < size1+size2 )
					{
						n = pm->GetMatch(m)-size1;
					}
					else
						if( size1<size2 && m < size1 )
						{
							n = pm->GetMatch(m)-size2;
						}
						else
							continue;

					start = (g[i].edgeList[m]<g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];
					end   = (g[i].edgeList[m]>g[j].edgeList[n])?g[i].edgeList[m]:g[j].edgeList[n];

					//cout << "Matching: " << m << ", " << n << endl;

					float temp_weight = PreMatrix.val(start, end);
					numerator += (PreMatrix.val(start, end)>1.0f)?1.0f:PreMatrix.val(start, end);
				}
			}
			delete pm;

			if (divMax)
			{
				weight_float = numerator / max( g[i].degree, g[j].degree );
			}
			else
			{
				weight_float = numerator / sqrt( (float)( g[i].degree * g[j].degree ) );
			}

			//cout << "m[" << i << "," << j << "]: " << weight_float << endl;
			if (weight_float > 1.0) {
				cout << " ERROR: S["<<i<<"]["<<j<<"] = " << numerator << "/max(" <<g[i].degree<< "," <<g[j].degree<<")"<<endl;
			}
			//assert(weight_float<1.1);
			Matrix[i][j] = (1-beta)*weight_float + beta;

			numerator = 0.0;
			Matching.clear();
			Weight.clear();
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void SimMatrix::UpdateRoleSimGreedyUndirected( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax )
{
	cout << "UndateRoleSimGreedyUndirected" << endl;
	VectIntType N_i, N_j;
	int N = g.num_vertices();
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0;			// take care of the diagonal
		N_i = g[i].edgeList;

		for( int j = i-1; j >= 0; j-- )	// take care of the lower triangle, excluding the diagonal
		{
			N_j = g[j].edgeList;
			Matrix[i][j] = (1-beta)*calcRoleSimGreedy(PreMatrix, N_i, N_j, divMax ) + beta;
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void SimMatrix::UpdateRoleSimGreedy( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax, float inWeight )

/* Two options for weighting in-vs-out similarity:
* If inWeight is in the range [0,1], then use it as expected: inWeight*inSim + (1-inWeight)*outSim
* Otherwise, use the sizes of the the in-vs-out neighbors to compute an inherent ratioing.
* Use either the max or the geometric mean size of the neighborhoods, selected by divMax. */
{
	if (!g.isDirected()) {
		UpdateRoleSimGreedyUndirected(g, PreMatrix, beta, divMax);
		return;
	}

	cout << "UndateRoleSimGreedy" << endl;
	VectIntType I_i, I_j, O_i, O_j;
	int sizeIi, sizeIj, sizeOi, sizeOj;
	float outWeight = 1 - inWeight;
	float denomI, denomO;
	//int N = g.num_vertices();

	// Tricky way to iterate over either lower triangle columns or upper triangle columns:
	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0; 			// take care of the diagonal
		I_i = g.in_edges(i);
		O_i = g.out_edges(i);
		sizeIi = I_i.size();
		sizeOi = O_i.size();
		for( int j = i + jIncrement; 0 <= j && j < N; j += jIncrement ) 	// do triangle
		{
			I_j = g.in_edges(j);
			O_j = g.out_edges(j);
			float inSimilarity = calcRoleSimGreedy(PreMatrix, I_i, I_j, divMax);	
			float outSimilarity = calcRoleSimGreedy(PreMatrix, O_i, O_j, divMax);

			//cout << "R(" << i << "," << j << ") = (" << 1-beta << ")*";
			//cout << "(" << inWeight << " * " << inSimilarity << " + " << outWeight << " * " << outSimilarity << ") + " << beta;

			if (0 <= inWeight && inWeight <= 1)
			{
				//float outSimilarity = 1 - inSimilarity;
				Matrix[i][j] = (1-beta)*(inWeight*inSimilarity + outWeight*outSimilarity) + beta;
			}
			else
			{
				sizeIj = I_j.size();
				sizeOj = O_j.size();
				if (divMax)
				{
					denomI = (float)( (sizeIi > sizeIj) ? sizeIi : sizeIj );
					denomO = (float)( (sizeOi > sizeOj) ? sizeOi : sizeOj );
				}
				else
				{
					denomI =  sqrt ( (float)(sizeIi * sizeIj) );
					denomO =  sqrt ( (float)(sizeOi * sizeOj) );
				}
				Matrix[i][j] = (1-beta)*(denomI*inSimilarity + denomO*outSimilarity)/(denomI+denomO) + beta;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////
float SimMatrix::calcRoleSimGreedy(SimMatrix& PreMatrix, const VectIntType& N_i, const VectIntType& N_j, bool divMax )
{
	vector<float> Weight;
	float weight;

	MatchItem matchItem;
	priority_queue<MatchItem, vector<MatchItem>, MatchLessThan> matchQueue;

	int size_i = N_i.size();
	int size_j = N_j.size();
	if( size_i == 0 || size_j == 0 )
	{
		return 0.0;
	}

	SimMatrixType& mPrev = PreMatrix.matrix();

	VectIntType::const_iterator iItr, jItr;
	matchItem.x = 0;
	for (iItr = N_i.begin(); iItr != N_i.end(); iItr++)
	{
		matchItem.y = 0;
		for (jItr = N_j.begin(); jItr != N_j.end(); jItr++)
		{
			if (*iItr == *jItr)
			{
				matchItem.weight = 1.0;
			}
			else if (*iItr > *jItr)
			{
				matchItem.weight = mPrev[*iItr][*jItr];
			}
			else {
				matchItem.weight = mPrev[*jItr][*iItr];
			}

			if (matchItem.weight > 0.0)
			{
				matchQueue.push(matchItem);
			}

			matchItem.y += 1;
		}
		matchItem.x += 1;
	}

	// Find the top weighted match-pairs
	VectIntType Flag1 = VectIntType( size_i, FALSE );
	VectIntType Flag2 = VectIntType( size_j, FALSE );
	float numerator = 0.0;
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

	//if (divideByMax)
	//{
		weight = numerator / (size_i > size_j ? size_i : size_j);
	//}
	//else
	//{
	//	weight = numerator / sqrt ( (float)(size_i * size_j) );
	//}
	//assert(weight<1.1);
	return weight;
}
////////////////////////////////////////////////////////////////////////
////    Bibligraphic Coupling (code bc)                             ////
////////////////////////////////////////////////////////////////////////
void SimMatrix::BibliographicCoupling( Graph& g )
/* Matrix(i,j) = |O_i intersect O_j|/|O_i union O_j| */
{
	VectIntType O_i, O_j, top, bot;

	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0; 			// take care of the diagonal
		O_i = g.out_edges(i);

		for( int j = i + jIncrement; 0 <= j && j < N; j += jIncrement ) 	// do triangle
		{
			O_j = g.out_edges(j);
			if (g.out_degree(i) == 0 || g.out_degree(j) == 0 )
			{
				Matrix[i][j] = 0.0;
			}
			else
			{
				set_intersection( O_i.begin(), O_i.end(), O_j.begin(), O_j.end(), back_inserter( top ) );
				set_union       ( O_i.begin(), O_i.end(), O_j.begin(), O_j.end(), back_inserter( bot ) );
				Matrix[i][j] = ((float)top.size()) /  bot.size();
				top.clear();
				bot.clear();
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
////    Co-Citation (code cc)     									////
////////////////////////////////////////////////////////////////////////
void SimMatrix::CoCitation( Graph& g )
{
	/* Matrix(i,j) = |I_i intersect I_j|/|I_i union I_j| */
	VectIntType I_i, I_j, top, bot;

	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0; 			// take care of the diagonal
		I_i = g.in_edges(i);

		for( int j = i + jIncrement; 0 <= j && j < N; j += jIncrement ) 	// do triangle
		{
			if (g.in_degree(i) == 0 || g.in_degree(j) == 0 )
			{
				Matrix[i][j] = 0.0;
			}
			else
			{
				I_j = g.in_edges(j);	
				set_intersection( I_i.begin(), I_i.end(), I_j.begin(), I_j.end(), back_inserter( top ) );
				set_union       ( I_i.begin(), I_i.end(), I_j.begin(), I_j.end(), back_inserter( bot ) );
				Matrix[i][j] = ((float)top.size()) /  bot.size();
				top.clear();
				bot.clear();
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
////    SimRank++ (code srpp, Antontelli, VLDB'08)        			////
////////////////////////////////////////////////////////////////////////
void SimMatrix::addSimRankPPEvidence(Graph& g, float inWeight)
{
/* SimRank++ (Antonellis, VLDB'08) make two additions to SimRank:
(1) supports weighted edges, where there may also be weight variance
(2) increases the similarity for pairs with more directedly shared neighbors,
which they call the Evidence factor.
We implement only the second feature (Evidence) at this time.
Evidence(i,j) = sum_i=1 to k (1/2)^k = 1 - (1/2)^k, where k = |N_i ^ N_j|
Since Evidence is a constant for a given vertex pair, we can simply tack on
the Evidence factor after regular SimRank has been computed. */
	
	cout << "addSimRankPPEvidence" << endl;
	float inWt = (g.isDirected()) ? inWeight : 0;   // if undirected, neutralize in-edges

	// Find max degree in the graph
	int maxDeg = 0;
	for (int i = 0; i < N; i++) {
		maxDeg = (g.in_degree(i) > maxDeg) ? g.in_degree(i) : maxDeg;
		maxDeg = (g.out_degree(i) > maxDeg) ? g.out_degree(i) : maxDeg;
	}
	cout << "Max degree in graph = " << maxDeg << endl;

	// Pre-compute evidence factors
	float factor = 0.5f;
	VectFloatType Evid;
	Evid.push_back(1.0f - factor);			// Evid[0] = 0.5
	for (int i = 1; i <= maxDeg; i++) {
		Evid.push_back(1.0f - factor);
		factor = factor / 2.0f;
	}
	cout << "Evidence factors: ";
	for (int j = 0; j <= maxDeg; j++) {
		cout << Evid[j] << " ";
	}

	// Multiply each element by its evidence factor.
	// For directed graphs, this is the weighted sum of the in-evidence and out-evidence
	int jIncrement = (lowerTri) ? -1 : +1;

	VectIntType I_i, I_j, I_iANDj;
	VectIntType O_i, O_j, O_iANDj;
	float outWt = 1 - inWt;
	for( int i = 0; i < N; i++ )
	{
		I_i = g.in_edges(i);
		O_i = g.out_edges(i);
		for( int j = i + jIncrement; 0 <= j && j < N; j = j + jIncrement )
		{
			I_j = g.in_edges(j);
			O_j = g.out_edges(j);
			I_iANDj.clear();
			O_iANDj.clear();
			set_intersection( I_i.begin(), I_i.end(), I_j.begin(), I_j.end(), back_inserter( I_iANDj ) );
			set_intersection( O_i.begin(), O_i.end(), O_j.begin(), O_j.end(), back_inserter( O_iANDj ) );
		
			Matrix[i][j] = (inWt * Evid[I_iANDj.size()] + outWt * Evid[O_iANDj.size()]) * Matrix[i][j];
		}
	}
}
////////////////////////////////////////////////////////////////////////
////    PSimRank (code psr, Fogaras, WWW'05)   						////
////////////////////////////////////////////////////////////////////////
void SimMatrix::UpdatePSimRank( Graph& g, SimMatrix& PreMatrix, float beta, float inWeight )
{
	/* Compute three components for each node pairs.
	Separate computations for in-neighbors and out-neighbors
	*/
	float inWt = (g.isDirected()) ? inWeight : 0; // if undirected, neutralize in-edges

	VectIntType I_i, I_j, I_iANDj, I_iNOTj, I_jNOTi;
	VectIntType O_i, O_j, O_iANDj, O_iNOTj, O_jNOTi;
	EdgeList::iterator N_i, N_j;
	float sim_I, sim_O;
	float outWt = 1 - inWt;

	int jIncrement = (lowerTri) ? -1 : +1;
	for( int i = 0; i < N; i++ )
	{
		Matrix[i][i] = 1.0;
		I_i = g.in_edges(i);
		O_i = g.out_edges(i);
		int sizeI_i = g.in_degree(i);
		int sizeO_i = g.out_degree(i);
		for( int j = i + jIncrement; 0 <= j && j < N; j = j + jIncrement )
		{
			I_j = g.in_edges(j);
			O_j = g.out_edges(j);
			int sizeI_j = g.in_degree(j);
			int sizeO_j = g.out_degree(j);

			// In neighbors
			if( sizeI_i == 0 || sizeI_j == 0 )
			{
				sim_I = 0.0;
			}
			else
			{
				I_iANDj.clear();
				I_iNOTj.clear();
				I_jNOTi.clear();
				set_intersection( I_i.begin(), I_i.end(), I_j.begin(), I_j.end(), back_inserter( I_iANDj ) );
				set_difference  ( I_i.begin(), I_i.end(), I_j.begin(), I_j.end(), back_inserter( I_iNOTj ) );
				set_difference  ( I_j.begin(), I_j.end(), I_i.begin(), I_i.end(), back_inserter( I_jNOTi ) );
				float fsizeI_iORj =(float)( I_iANDj.size() + I_iNOTj.size() + I_jNOTi.size() );

				// Term 1: similarity between I_i and I_j
				float SI_iANDj = (float) I_iANDj.size();

				// Term 2: similarity between I_iNOTj and I_j
				float sum_iNOTj = 0.0;
				for (N_i = I_iNOTj.begin(); N_i != I_iNOTj.end(); N_i++)
				{
					for (N_j = I_j.begin(); N_j != I_j.end(); N_j++)
					{
						// SimMatrix.val takes care of whether xCoord > yCoord or not
						sum_iNOTj += PreMatrix.val(*N_i, *N_j);
					}
				}
				float SI_iNOTj = sum_iNOTj / (float)sizeI_j;

				// Term 3: similarity between I_i and I_jNOTi
				float sum_jNOTi = 0.0;
				for (N_i = I_i.begin(); N_i != I_i.end(); N_i++)
				{
					for (N_j = I_jNOTi.begin(); N_j != I_jNOTi.end(); N_j++)
					{
						// SimMatrix.val takes care of whether xCoord > yCoord or not
						sum_jNOTi += PreMatrix.val(*N_i, *N_j);
					}
				}
				float SI_jNOTi = sum_jNOTi / (float)sizeI_i;

				// Combine the 3 terms
				sim_I = (SI_iANDj + SI_iNOTj + SI_jNOTi) / fsizeI_iORj;
			}


			// Out neighbors
			if( sizeO_i == 0 || sizeO_j == 0 )
			{
				sim_O = 0.0;
			}
			else
			{
				O_iANDj.clear();
				O_iNOTj.clear();
				O_jNOTi.clear();
				set_intersection( O_i.begin(), O_i.end(), O_j.begin(), O_j.end(), back_inserter( O_iANDj ) );
				set_difference  ( O_i.begin(), O_i.end(), O_j.begin(), O_j.end(), back_inserter( O_iNOTj ) );
				set_difference  ( O_j.begin(), O_j.end(), O_i.begin(), O_i.end(), back_inserter( O_jNOTi ) );
				float fsizeO_iORj =(float)( O_iANDj.size() + O_iNOTj.size() + O_jNOTi.size() );

				// Term 1: similarity between O_i and O_j
				float SO_iANDj = (float) O_iANDj.size();

				// Term 2: similarity between O_iNOTj and O_j
				float sum_iNOTj = 0.0;
				for (N_i = O_iNOTj.begin(); N_i != O_iNOTj.end(); N_i++)
				{
					for (N_j = O_j.begin(); N_j != O_j.end(); N_j++)
					{
						// SimMatrix.val takes care of whether xCoord > yCoord or not
						sum_iNOTj += PreMatrix.val(*N_i, *N_j);
					}
				}
				float SO_iNOTj = sum_iNOTj / (float)sizeO_j;

				// Term 3: similarity between O_i and O_jNOTi
				float sum_jNOTi = 0.0;
				for (N_i = O_i.begin(); N_i != O_i.end(); N_i++)
				{
					for (N_j = O_jNOTi.begin(); N_j != O_jNOTi.end(); N_j++)
					{
						// SimMatrix.val takes care of whether xCoord > yCoord or not
						sum_jNOTi += PreMatrix.val(*N_i, *N_j);
					}
				}
				float SO_jNOTi = sum_jNOTi / (float)sizeO_i;

				// Combine the 3 terms
				sim_O = (SO_iANDj + SO_iNOTj + SO_jNOTi) / fsizeO_iORj;
			}
			Matrix[i][j] = (1-beta)*(inWt * sim_I + outWt * sim_O);
		}
	}
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

