// SimMatrix: a class for a node-node similarity matrix, with functions for RoleSim, SimRank, and related algorithms.
#ifndef _SIMMATRIX_H
#define _SIMMATRIX_H

//#include "Util.h"
#include "Graph.h"
#include <map>
#include <vector>
#include <iostream>
using namespace std;

#define DENDROGRAM_BASE 0.01
#define SCALE_FOR_EXACTMATCHING 10000 
#define TRUE 1
#define FALSE 0
typedef map<VectIntType, VectIntType> MapVectIntToVectInt;
typedef vector< vector<float> > SimMatrixType;

// Define some constants to make simMethod processing easier to read:
const int ROLESIM = 0;
const int SIMRANK = 1;
const int SIMRANK_PP = 2;
const int P_SIMRANK = 3;
const int BIBLIO_COUPLING = 4;
const int CO_CITATION = 5;

//MatchItem is used to record a matching pair (calcRoleSimGreedy, find*TopPairs)
struct MatchItem {
public:
	int		x;
	int		y;
	float	weight;
};
struct MatchLessThan {
	bool operator() (MatchItem a, MatchItem b) {
		return (a.weight < b.weight);
	}
};
struct MatchGreaterThan {
	bool operator() (MatchItem a, MatchItem b) {
		return (a.weight > b.weight);
	}
};
typedef vector<MatchItem> MatchItemVec;

class SimMatrix
{
private:
	SimMatrixType Matrix;
	int N;
	bool lowerTri;

public:
	SimMatrix(void);
	SimMatrix(int N);
	virtual ~SimMatrix(void);
	float val(int i, int j);
	SimMatrixType& matrix();
	void set(int i, int j, float val);
	void Clear();
	bool Diff( SimMatrix& Matrix2, float threshold);
	bool DiffRelative( SimMatrix& Matrix2, float threshold, int turn );
	bool DiffRelative( SimMatrix& previous, float threshold );

	void Print( ostream& out, int precision, bool header=true );
	void PrintCompact( ostream& out, int precision, bool compact, bool header=false );
	void PrintMatlab( ostream& out, int precision, bool offset=true );
	void PrintMatlabCompact( ostream& out, int precision, bool compact);
	void loadSimVec ( char* infile, int N, int precision, bool distance );

	void InitializeUniform( );
	void InitializeDiagonal( );
	void InitializeDepthZero( Graph& g, bool directed=false );
	void InitializeDepthOne(  Graph& g, bool directed=false );
	void InitializeRatioDepth0( Graph& g, float beta );

	int simCode(char* simMeasure);

	// SimRank (code sr)
	void UpdateSimRank          ( Graph& g, SimMatrix& PreMatrix, float beta, float inWeight );
	void UpdateSimRankUndirected( Graph& g, SimMatrix& PreMatrix, float beta );

	// RoleSim (code rs)
	void UpdateRoleSimExact           ( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax );
	void UpdateRoleSimExact           ( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax, float inWeight );
	void UpdateRoleSimGreedy          ( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax, float inWeight );
	void UpdateRoleSimGreedyUndirected( Graph& g, SimMatrix& PreMatrix, float beta, bool divMax );
	float calcRoleSimGreedy( SimMatrix& PreMatrix, const VectIntType& N_u, const VectIntType& N_v, bool divMax );
	

	void BibliographicCoupling( Graph& g );
	void CoCitation( Graph& g );
	void UpdatePSimRank          ( Graph& g, SimMatrix& PreMatrix, float beta, float inWeight );
	void addSimRankPPEvidence(Graph& g, float inWeight);
};

#endif
