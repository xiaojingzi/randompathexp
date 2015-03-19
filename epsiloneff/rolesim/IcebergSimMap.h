/* IcebergSimMap: a class for scalable RoleSim computation.  While providing a more efficient data structure than SimMatrix
  (a nodePair->simValue map for only the most similar nodePairs, discarding nodePairs that are not very similar),
  its functions are analogous to those in SimMatrix.  In the future, this class may be re-implemented as a subclass
  of SimMatrix or some new generic similarity superclass
  */
#ifndef _ICEBERGSIMMAP_H
#define _ICEBERGSIMMAP_H

#include "Graph.h"
#include "SimMatrix.h"
//#include "Util.h"
using namespace std;

/*
 * This class encapsulates a data structure map< <int,int>, float > which seeks to store
 * the similarity values for only the most similar vertex pairs, therefore using less memory
 * that SimMatrix, which stores the full (triangular) matrix.
 * The Initialize function uses the degrees of a given pair and the degrees of their neighbors
 * to compute an upper bound RoleSim value.  If the value is above a threshold 'waterline',
 * it is included in the Map.  The similarity for any vertex pair not in the map is defined
 * as ratioWt*degree(u)/degree(v) + beta, where ratioWt is between 0 and 1.
 */

struct FloatPair {
	float even;
	float odd;
};
typedef map< pair<int, int>, float> SimMapType;

class IcebergSimMap
{
private:
	SimMapType simMap;
	Graph *gPtr;
	int N;
	bool lowerTri;
	float waterline;
	float ratioWt;
	float beta;
	float inWeight;
	bool verbose;

public:
	IcebergSimMap(void);
	IcebergSimMap(float waterline, float ratioWt, float beta, float inWeight, Graph& g, bool verbose=false);
	~IcebergSimMap(void);
	float getWaterline();
	float getRatioWt();
	int getN();
	SimMapType& getMap();
	void setParameters(float w, float r, float b, float i);
	void setGraph(Graph& g);
	void insert (int x, int y, float val);
	float val(int x, int y);
	float valUnsafe(int x, int y);

	bool Diff( IcebergSimMap& MatrixOdd, float threshold);
	bool DiffRelative( IcebergSimMap& MatrixOdd, float threshold, int turn );
	bool DiffRelative( IcebergSimMap& PrevIceberg, float threshold );

	void Print( ostream& out, int precision, bool header=true );
	void PrintCompact( ostream& out, int precision, bool compact, bool header=true );
	void PrintMatlab( ostream& out, int precision, bool offset=true );
	void PrintMatlabCompact( ostream& out, int precision, bool compact);
	void PrintTable (ostream& out);

	void save( ostream& out, char* graphFilename, int precision );
	void save( ostream& out, char* graphFilename, int precision, bool compact );
	int load( istream& in, int precision, Graph& g, bool compact );
string timeToString(double time);
void displayElapsedTime(time_t& previous, ostream& out, bool newline);
	void InitializeRatioDepth1( Graph& g, bool shortcutMatching );
	void updateRoleSim(IcebergSimMap& preIceberg, float beta, bool divMax, bool melt);
};
#endif
