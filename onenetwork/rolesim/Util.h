/* Util: Collection of functions used for RoleSim/SimRank post-processing.  In general, the functions here do
  not need direct knowledge of the original graph, only the SimMatrix or IcebergSimMap containing the results
 */
#ifndef _UTIL_H
#define _UTIL_H

//#include "SimMatrix.h"
#include "IcebergSimMap.h"
#include <map>

#ifndef _MSC_VER
	#include <sys/time.h>
#endif
//#include <hash_map.h>

using namespace std;

// Block* and Author* structs are used by the coauthor analysis functions
struct BlockItem {
public:
	int				s1;
	int				s2;
	int				size;
	int				sampleSize;
	double			sumSim;
	double			sumRank;
	double			sumSimSq;
	double			sumRankSq;
};
struct BlockCompareLevel {
	bool operator() (BlockItem a, BlockItem b) {
		if ( (a.s1 == b.s1) ) {
			return (a.s2 < b.s2);
		}
		else {
			return (a.s1 < b.s1);
		}
	}
};
struct BlockCompareAvgRank {
	bool operator() (BlockItem a, BlockItem b) {
		return (a.sumRank/a.sampleSize < b.sumRank/b.sampleSize);
	}
};
struct BlockCompareAvgSim {
	bool operator() (BlockItem a, BlockItem b) {
		return (a.sumSim/a.sampleSize < b.sumSim/b.sampleSize);
	}
};
struct BlockCompareDiffThenRank {
	bool operator() (BlockItem a, BlockItem b) {
		if ( (a.s2 - a.s1) == (b.s2 - b.s1) ) {
			return (a.sumRank/a.sampleSize < b.sumRank/b.sampleSize);
		}
		else {
			return ( (a.s2 - a.s1) < (b.s2 - b.s1) );
		}
	}
};

struct AuthorItem {
	string name;
	int id;
	double Grank;
	double Hrank;
};


class Util {
	private:
		int vsize;
		map<int, int> pctRankHash; // findPercentileRank
		bool verbose;

	public:
		Util(int N, bool verbose=false);
		~Util();
// PostProcessing functions: probably belong in SimMatrix or a separate class

// General Purpose Analysis Functions:
// Comparing simVecs (a simVec is a compact reprsentation of a triangular similarity matrix)
	vector<int> ReloadSimilarityTriangleVector(const char* inputFilename, int N, int bitsPrecision, bool invert);
	void        histogramSimValues(const vector<int>& simVec, int numBins, int bitsPrecision);
	void		compareTwoSimVecs(const vector<int>& oldVec, const vector<int>& newVec, int bitsPrecision);

// PercentileRank: Read sim results, Compute/Write/Read percentileRankBins
	vector<int> computePercentileRankBins(vector<int>& simVec, int numBins);
	vector<int> computePercentileRankBins(IcebergSimMap& iceberg, int numBins, int bitsPrecision);
	void        writePercentileRankBins(ostream& out, vector<int>& pctileVec);
	vector<int> readPercentileRankBins(string rankFileName);
	int         findPercentileRank(const vector<int>& pctileVec, int value, bool verbose);

	void		rankMatrixValues(const vector<int>& simVec, int bitsPrecision, map<int,float>& matPctRank);
	void		compareTwoRankings(const vector<int>& simVec1, const vector<int>& simVec2,
				map<int,float>& rankMap1, map<int,float>& rankMap2, int bitsPrecision);

	void		findAndSortTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
				MatchItemVec& topPairs, double Kmax);
	int			getVal(const vector<int>& simVec, int x, int y, int N);


// Block Analysis //////////////////////////////////////////////////////

// Shell Analysis //////////////////////////////////////////////////////
	void doShellAnalysis(vector<int>& simVec, vector< vector<int> >& shells, vector<int>& pctileVec,
			int bitsPrecision, bool verbose);
	
		
		vector<int> convertShellListsToMap(const vector< vector<int> >& shells, int N);
		void findShellsOfTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
								const vector<int> &shellMap, int numShells, double K, int bitsPrecision);
		void computeWithinShellSimilarity(vector< vector<int> >& shells, const vector<int>& simVec,
				const vector<int>& pctileVec, int bitsPrecision, bool includeSelfSim, bool verbose);
		void computeCrossShellSimilarity(vector< vector<int> >& shells, const vector<int>& simVec,
				const vector<int>& pctileVec, int bitsPrecision, bool verbose);
// Co-author Analysis //////////////////////////////////////////////////////
	void doCoauthorAnalysis(vector<int>& simVec, vector<AuthorItem> authorVec, vector<int>& pctileVec,
		int bitsPrecision, int size1, bool verbose);
		
		vector<AuthorItem> readAuthorRanks(char* rankFileName);
		vector< vector<int> > binAuthorsSublistByRank(const vector<AuthorItem>& authorVec, int startIndex, int subLength);

		void findAuthorRankOfTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
								const vector<AuthorItem> authorVec, int bitsPrecision, int size1);
		void computeCrossComponentSimilarity(vector< vector<int> >& bins1, vector< vector<int> >& bins2,
			const vector<int>& simVec, const vector<int>& pctileVec, int bitsPrecision, bool verbose);

// Iceberg Accuracy Analysis //////////////////////////////////////////////////////
	MatchItemVec recordTopPairs(const vector<int>& simVec, const vector<int>& pctileVec, int bitsPrecision, int minVal, bool valIsRank);
	void rankIcebergValues( IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
		map<int,float>& icePctRank, map<int,float>& matPctRank);

	void compareIcebergToSimMatrix(IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
				 map<int,float>& iceRank,  map<int,float>& matRank);
	void checkIcebergCoverage(IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
				 map<int,float>& iceRank,  map<int,float>& matRank);
};

#endif
