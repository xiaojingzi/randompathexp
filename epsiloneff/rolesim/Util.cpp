/* Util: Collection of functions used for RoleSim/SimRank post-processing.  In general, the functions here do
  not need direct knowledge of the original graph, only the SimMatrix or IcebergSimMap containing the results
 */
#include "Util.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <map>
#include <set>
#include <vector>
#include <assert.h>
#include <iterator>
#include <list>
#include <deque>
#include <utility>
#include <cmath>
#include <iomanip>
#include <queue>

/////////////////////////////////////////////////////////////////////////
//// Constructor(s) and Deconstructor
/////////////////////////////////////////////////////////////////////////
Util::Util(int N, bool verb)
{
	vsize = N;
	verbose = verb;
}
Util::~Util() {}
/////////////////////////////////////////////////////////////////////////
int	Util::getVal(const vector<int>& simVec, int x, int y, int N)
{
	int nMax = N - 1;
	int colBack = (nMax - y)*(nMax - y - 1)/2;
	int rowBack = (nMax - x);
	return simVec.at(simVec.size() - 1 - colBack - rowBack);
}
/////////////////////////////////////////////////////////////////////////
//// Post-Processing
/////////////////////////////////////////////////////////////////////////
vector<int> Util::ReloadSimilarityTriangleVector(const char* inputFilename, int numVert, int bitsPrecision, bool invert) {
	// Read a 'compact' similarity matrix file, and store the values in a vector<int> structure.
	cout << "calling ReloadSimilarityTriangleVector(" <<inputFilename<< "," <<numVert<< "," <<bitsPrecision<< "," <<invert<< ")" << endl;

	ifstream in(inputFilename);
	int vecLen = numVert*(numVert-1)/2;		// vector length = (numVert choose 2)
	vector<int> simVec(vecLen);
	int range = (int) pow(10.0, (double)bitsPrecision);
	cout << "numVert = " << numVert << ", vecLen = " << vecLen << ", range = " << range << endl;

	int row = 1;
	int i = 0;
	int value;
	int count;
	while (!in.eof() && i < vecLen) {
		if (in >> value) {
			simVec[i] = (invert) ? (range - value) : value;
			count = i+1;
			if (count % 100000 == 0) { // to show progress, print a dot every 100,000 values
				cout << ".";
				if (count % 5000000 == 0) { // with up to 50 dots per line (5 million cells)
					cout << endl;
				}
			}
			if (count % (vecLen/10) == 0) { cout << (count / (vecLen/10)) << "0%" << endl; }
		}
		else {
			cout << "ERROR at i = " << i << ": input is not an integer" << endl;
			break;
		}
		i++;
		row++;
	}
	if (i != vecLen) {
		cout << "ERROR: Wrong number of values" << endl;
	}
	cout << "number of values read = " << i << endl;
	return simVec;
}
///////////////////////////////////////////////////////////////////////////
void Util::compareTwoSimVecs(const vector<int>& oldVec, const vector<int>& newVec, int bitsPrecision) {

	int diff;
	int maxDiff = 0;
	int totDiff = 0;

	double relDiff;
	double maxRelDiff = 0.0;
	double totRelDiff = 0.0;
	
	double avgDiff;
	double avgRelDiff;


	double threshold[] = {0.001, 0.005, 0.01, 0.05, 0.1};
	int count[] = {0, 0, 0, 0, 0};
	int numThres = 5;
	int i;

	for (int c = 0; c < oldVec.size(); c++) {
		diff = abs( newVec[c] - oldVec[c] );
		if (diff > maxDiff) {
			maxDiff = diff;
		}
		totDiff += diff;

		relDiff = diff / (double) oldVec[c];
		if (relDiff > maxRelDiff) {
			maxRelDiff = relDiff;
		}
		totRelDiff += relDiff;

		i = 0;
		while (relDiff >= threshold[i] && i < numThres) {
			count[i]++;
			i++;
		}
	}

	avgDiff = totDiff / (double)oldVec.size();
	avgRelDiff = totRelDiff / oldVec.size();

	double range = pow(10.0, (double)bitsPrecision);
	cout << "range = " << range << endl;

	cout << "maxDiff = " << maxDiff/range << ",\t totDiff = " << totDiff/range << ",\t avgDiff = " << avgDiff/range << endl;
	cout << "maxRelDiff = " << maxRelDiff << ",\t totRelDiff = " << totRelDiff << ",\t avgRelDiff = " << avgRelDiff << endl;
	cout << endl;

	double pctFactor = 100.0/oldVec.size();
	cout << "Number/fraction of cells with relative difference greater than T" << endl;
	    printf("%4s  %12s  %12s\n","Thres", "number", "% of nodes");
	//cout << "T\t" << "number\t"  << "fraction\n";
	for (i = 0; i < numThres; i++) {
		printf("%.4f  %12d  %12.5f%%\n",threshold[i], count[i], count[i]*pctFactor);
	}
}
///////////////////////////////////////////////////////////////////////////
void Util::histogramSimValues(const vector<int>& simVec, int numBins, int bitsPrecision) {
	// Report how many similarity values fall within each bin range of values.
	cout << "histogramSimValues" << endl;
	int vecLen = simVec.size();
	double range = pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	int binStep = (int) (range / numBins);	// range of values for one step or bin
	cout << " range=" <<range<< ", numBins=" <<numBins<< ", binStep=" <<binStep<< endl;

	int binIndex;
	vector<int> bin(numBins, 0);
	for (int i = 0; i < vecLen; i++) {
		binIndex = min (simVec[i] / binStep, numBins - 1);
		bin[binIndex]++;
	}
	for (int b = 0; b < numBins; b++) {
		cout << b << "\t" << bin[b] << endl;
	}
}

/////////////////////////////////////////////////////////////////////////
vector<int> Util::computePercentileRankBins(vector<int>& simVec, int numBins) {
/* WARNING: This function sorts the input vector simVec.
 * The vector is therefore modified by this function!!!
 */
	cout << "computePercentileRankBins: " << numBins << " bins" << endl;
	//vector<int> simVec = simVecOrig;
	int vecLen = simVec.size();			// length of the vectorized similarity matrix
	float binStepF = vecLen / (float) numBins;	// number of sim values in one percentile bin
	cout << "binStepF = " << ios::fixed << binStepF << endl;
	vector<int> pctileVec(numBins+1);

	sort(simVec.begin(), simVec.end());
	
	int binIndex;
	//pctileVec[0] = 0;
	for (int i = 0; i < numBins; i++) {
		binIndex = (int) (i * binStepF);
		pctileVec[i] = simVec[binIndex];			// similarity value at this point
		cout << i << "\t" << i/(float)numBins << "\t" << pctileVec[i] << endl;
		//cout << "pctileVec[" << i << "] = simVec[" << binIndex << "] = " << pctileVec[i] << endl;
	}
	// Set the last bin at the end of the vector
	pctileVec[numBins] = simVec[vecLen - 1];
	cout << numBins << "\t" << 1 << "\t" << pctileVec[numBins] << endl;

	return pctileVec;
}
/////////////////////////////////////////////////////////////////////////
vector<int> Util::computePercentileRankBins(IcebergSimMap& iceberg, int numBins, int bitsPrecision ) {
// The challenge is to generate a histogram of values w/o materializing the sim matrix.

	cout << "computePercentileRankBins(iceberg): " << numBins << " bins" << endl;

	// 1. Scan the whole virtual sim matrix. Count how many times each sim value occurs.
	// That is, build a map whose key->value are sim_value->count.
	map<float,int> valCount; // for ever
	vsize = iceberg.getN();
	float val;
	for (int y = 0; y < vsize; y++) {
		for (int x = y+1; x < vsize; x++) {
			val = iceberg.val(x,y);
			if (valCount.find(val) == valCount.end()) {
				valCount[val] = 1;
			}
			else {
				valCount[val]++;
			}
		}
	}
	int vecLen = vsize*(vsize-1)/2;		// length of the (virtual) vectorized similarity matrix
	cout << "Virtual simVec contains " << vecLen << " cells with " << valCount.size() << " unique values" << endl;
	int cumulCount = 0;
	map<float,int>::iterator vit;
	for (vit = valCount.begin(); vit != valCount.end(); vit++) {
		cumulCount += vit->second;
		cout << setprecision(5) << vit->first << "\t" << vit->second << "\t" << cumulCount<< endl;
	}

	// 2. Sort the map: The map's keys are already sorted, by definition of C++ STL map.
	
	// 3. Compute that number cells per percentile bin
	double binStepD = vecLen / (double) numBins;	// number of sim values in one percentile bin
	//cout << "binStepD = " << ios::fixed << binStepD << endl;
	cout << "binStepD = " << vecLen << " / " << numBins << " = " << binStepD << endl;
	vector<int> pctileVec(numBins+1);

	// 4. Compute the cumulative # cells that have a value less than or equal to each value.
	// Insert into pctileVec[i] the sim value when the cumul count reaches i*binStepF.
	
	double range = pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	int i = 0;
	val = 0;
	//pctileVec[i] = val;
	//cout << i << "\t" << i/(float)numBins << "\t" << pctileVec[i] << endl;
	//i++;
	map<float,int>::iterator mit = valCount.begin();
	cumulCount = 0;
	
	
	while (i <= numBins  && mit != valCount.end()) {
		// Keep adding to the count until reaching the next binStep threshold
		while (cumulCount < (int)(i * binStepD) && mit != valCount.end()) {
			val = mit->first;
			cumulCount += mit->second;
			mit++;
		}
		// If the last count increment was huge, it might encompass more than one percentile step
		while ( cumulCount >= (int)(i * binStepD) && i <= numBins) {
			//val = mit->first;
			pctileVec[i] = (int)(val * range);
			cout << i << "\t" << i/(float)numBins << "\t" << pctileVec[i] << endl;
			i++;
		}
	}
	// If somehow we reached the end of the valCount map without finishing, fill up the remaining
	// with the last (and highest) value from the map.
	while (i <= numBins) {
		pctileVec[i] = (int)(val * range);
		cout << i << "\t" << i/(float)numBins << "\t" << pctileVec[i] << endl;
		i++;
	}
	cout << numBins << "\t" << 1 << "\t" << pctileVec[numBins] << endl;
	return pctileVec;
}
///////////////////////////////////////////////////////////////////
void Util::writePercentileRankBins(ostream& out, vector<int>& pctileVec) {
	vector<int>::iterator it;
	for (it = pctileVec.begin(); it != pctileVec.end(); it++) {
		out << *it << endl;
	}
}
///////////////////////////////////////////////////////////////////
vector<int> Util::readPercentileRankBins(string rankFileName){

	cout << "starting readPercentileRankBins" << endl;
	ifstream in(rankFileName.c_str());
	int intValue;
	vector<int> pctileVec;
	int i = 0;
	while (!in.eof()) {
		if (in >> intValue) {
			pctileVec.push_back(intValue);
		}
		else {
			cout << "Successfully read " << i-1 << " integer values" << endl;
		}
		i++;
	}
	return pctileVec;
}
///////////////////////////////////////////////////////////////////
MatchItemVec Util::recordTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
			int bitsPrecision, int minVal, bool valIsRank) {
/*	Scan the simVec, tracking the virtual (x,y) coordinates.
	If valIsRank = false:
		if simVec[i] >= minVal, then add  (pair(x,y), val) to TopValueVec.
	If valIsRank = true:
		if simVec[i] >= pctileRank[minVal], then add (pair(x,y), val) to TopValues map.
*/
	cout << "Starting recordTopPairs" << endl;
	MatchItem matchItem;
	MatchItemVec topPairs;
	int threshold;

	// Compute some parameters
	float range = (float)pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	int N = (int)( sqrt(2.0*simVec.size()) + 1);
	cout << "simVec.size() = " << simVec.size() << ", num Vertices = " << N << endl;
	if (valIsRank) {
		threshold = pctileVec[minVal];
		cout << "valIsRank is true. minRank = " << minVal << ", val[minRank] = " << threshold << endl;
	}
	else {
		threshold = minVal;
	}

	// Refer to SimMatrix::PrintCompact to understand the (row, col) pattern
	vector<int>::const_iterator sPtr = simVec.begin();
	for (int col = 0; col < N-1; col++) {
		for (int row = col + 1; row < N; row++) {
			if (*sPtr >= threshold) {
				matchItem.x = row;
				matchItem.y = col;
				matchItem.weight = (*sPtr)/ range;
				topPairs.push_back(matchItem);
			}
			sPtr++;
		}
	}
	int topSize = topPairs.size();
	double topPct = (100.0 * topSize)/simVec.size();
	cout << "Size of topValueVec = " << topSize << " = " << topPct << "% of the full simVec" << endl;
	cout << "Finished with recordTopPairs" << endl;
	return topPairs;
}
///////////////////////////////////////////////////////////////////
void Util::findAndSortTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
				MatchItemVec& topPairs, double Kmax) {
	
	int vecLen = simVec.size();
	int thresholdValue;
	if (Kmax > 1) {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax/(double)vecLen ))];
	}
	else {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax ))];
	}
	cout << "findAndSortTopPairs: thresholdValue = " << thresholdValue << endl;

	// (3) For each cell whose value > threshold percentile (step 2), save the tuple (X, Y, value)
	
	MatchItem m;
	int i = 0;
	for (int y = 0; y < vsize; y++) {	// y,x iteration order MUST match the order of simVec
		for (int x = y + 1; x < vsize; x++) {
			if (simVec[i] > thresholdValue) {
				m.x = x;
				m.y = y;
				m.weight = (float) simVec[i];
				topPairs.push_back(m);
			}
			i++;
		}
	}
	cout << topPairs.size() << " out of " << vecLen << " items saved" << endl;

	// (4) Sort this list of top values.
	sort(topPairs.begin(), topPairs.end(), MatchLessThan());
	//cout << "topPairs sorted" << endl;
}
///////////////////////////////////////////////////////////////////
int Util::findPercentileRank(const vector<int>& pctileVec, int value, bool verbose) {
	// Do binary search and interpolation to find the closest integer percentile rank
	// for the input value.   This function needs to be computationally efficient,
	// because it will be called once for every matrix value.
	int vecSize = pctileVec.size();
	int rank;

	if (vecSize < 2) { // error check
		cout << "ERROR in findPercentilerank: pctileVec is too small" << endl;
		return -1;
	}

	if (verbose) cout << "findPercentileRank for value = " << value << endl;
	// 1. Check if the value is stored in the hash table
	map<int,int>::iterator mit = pctRankHash.find(value);
	if (mit != pctRankHash.end()) {
		rank = mit->second;
		if (verbose) cout << "RankHash["<<value<<"]=" << rank << endl;
		return rank;
	}

	// 2. Else, find the percentile bins immediately above and below the input value
	int lower = 0;
	int upper = vecSize - 1;
	while ( (upper - lower > 1) && (pctileVec[lower] != pctileVec[upper]) ) {
		if (verbose && upper < 27) {
			cout << "  pctileVec[upper] - pctileVec[lower] = " << pctileVec[upper] - pctileVec[lower] << endl;
		}
		int mid = (upper + lower)/ 2;
		if (value > pctileVec[mid]) {
			lower = mid;
		}
		else {
			upper = mid;
		}
	}
	if (verbose) cout << " step2: lower=" << lower << ", upper=" << upper << endl;

	// 3. If the value spans more than one rank, find "average" rank for this value
	if (pctileVec[lower] == pctileVec[upper]) {
		while (lower > 0 && pctileVec[lower - 1] == pctileVec[upper]) {
			lower--;
		}
		while (upper < vecSize-1 && pctileVec[lower] == pctileVec[upper + 1]) {
			upper++;
		}
		if (verbose) cout << " step3: lower=" << lower << ", upper=" << upper << endl;
		rank = (lower + upper)/2;
		pctRankHash[value] = rank;
		if (verbose) cout << "Average Rank["<<value<<"]=" << rank << endl;
	}
	else {
		// 4. Else, linearly interpolate between the upper and lower bins
		// to estimate the nearest integer percentile rank for the input value
		double ratio = max( (value - pctileVec[lower])/(double)(pctileVec[upper] - pctileVec[lower]), 0.0 );
		rank = lower + (int)( ratio*(upper - lower) + 0.49 );
		pctRankHash[value] = rank;
		if (verbose) {
			cout << " interpolate between pVec["<<lower<<"]="<<pctileVec[lower]<<" and pVec["<<upper<<"]="<<pctileVec[upper]
				<< "=> ratio="<<ratio<<", rank="<<rank<<endl;
		}
	}
	return rank;
}
///////////////////////////////////////////////////////////////////////////
//// Shell Analysis Functions
///////////////////////////////////////////////////////////////////////////
void Util::computeWithinShellSimilarity(vector< vector<int> >& shells, const vector<int>& simVec,
			const vector<int>& pctileVec, int bitsPrecision, bool includeSelfSim, bool verbose) {
	cout << endl << "STARTING computeWithinShellSimilarity" << endl;
	// Recall that simVec is the vectorization of the lower triangle similarity matrix.
	// matrix location: (1,0),(2,0)...(n-1,0),(2,1),(2,1)...(n-1,1),...,(n-1,n-2)
	// Rules: The yth column has (n-1 - y) rows.  Total #cells = (n choose 2) = n*(n-1)/2.
	// Remember: the diagonal cells are not in the vector.
	// To find any (x,y): (1) Flip x <-> y so that x > y.
	// (2) Count from the end of the vector backwards
	// (because an end-aligned subset of the triangle is just a smaller triangle).
	// For (x,y) = (n-1,n-2) ==> simVec(length - 1);
	// If y < n-2, count back a sub-triangle of size (n-2 - y)*(n-1 - y)/2
	// If x < n-1, count back (n-1 - x)
	// Tests: Let n = 10.  Length = (10 choose 2) = 45.
	// (1,0) = (45-1) - [(10-2 - 0)*(10-1 - 0)/2] - [10-1 - 1] = 44 - 36 - 8 = 0
	// (9,0) = (45-1) - [(10-2 - 0)*(10-1 - 0)/2] - [10-1 - 9] = 44 - 36 - 0 = 8
	// (2,1) = (45-1) - [(10-2 - 1)*(10-1 - 1)/2] - [10-1 - 2] = 44 - 28 - 7 = 9
	// (9,8) = (45-1) - [(10-2 - 8)*(10-1 - 8)/2] - [10-1 - 9] = 44 - 0 - 0 = 44

	// Define: some useful quantities
	int vecMax = simVec.size() - 1;
	int numShells = shells.size();	// skip shells[0]
	int range = (int)pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	double twoPow32 = (float)(1 << 16) * (float)(1 << 16);
	int nMax = (int)( sqrt(2.0*simVec.size()) + 1) - 1;
	cout << "nMax computed from simVec.size = " << nMax << endl;

	// print header
	cout << "\tShell\t#Nodes\t#Pairs\t\tSumSim\tAvgSim\tStdvSim\t\tSumRank\tAvgRank\tStdvRank" << endl;

	// iterate through each shell
	double grandSumSimD = 0.0;
	double grandSumRankD = 0.0;
	double grandSumSimSqD = 0.0;
	double grandSumRankSqD = 0.0;
	int grandNumPairs = 0;

	for (int s = numShells-1; s >= 0; s--) {

		int shellSize = shells[s].size();
		int numPairs = (includeSelfSim) ? shellSize*shellSize : shellSize*(shellSize - 1)/2;
		if (verbose) cout << "Shell " << s << " has " << shellSize << " nodes " << endl;

		unsigned int sumSim = 0;
		unsigned int sumRank = 0;
		unsigned int sumSimSq = 0;
		unsigned int sumRankSq = 0;
		unsigned int overflowSim = 0;
		unsigned int overflowRank = 0;
		unsigned int overflowSimSq = 0;
		unsigned int overflowRankSq = 0;

		for (int j = 0; j < shellSize; j++) {

			int y = shells[s].at(j);
			if (verbose) cout << "y=" << y << ": x=";
			int colBack = (nMax - y)*(nMax - y - 1)/2;
			int xStart = (includeSelfSim) ? 0 : j+1;

			for (int i = xStart; i < shellSize; i++) {

				int x = shells[s].at(i);
				if (verbose) cout << x << ",";
				int rowBack = (nMax - x);
				int sim;
				if (x == y) {
					sim = range;
				} else {
					sim = simVec.at(vecMax - colBack - rowBack);
					if (verbose) cout << "simVec["<<vecMax - colBack - rowBack<<"]; ";
				}
				if (sim < 0 || sim > range) cout << " error: sim = " << sim << " at ("<<x<<","<<y<<")" << endl;
				int rank = findPercentileRank(pctileVec, sim, false);
				if (rank < 0 || rank > 100) cout << " error: rank of sim["<<sim<<"]=" << rank << " at ("<<x<<","<<y<<")" << endl;
				
				// Compute sumSim, sumRank, sumSimSq, and sumRankSq.  Deal with overflow
				unsigned int prevSumSim = sumSim;
				sumSim += sim;
				if (sumSim < prevSumSim) {
					overflowSim++;
					if (verbose) cout << " sumSim overflow in shell " <<s<<" at ("<<x<<","<<y<<")" << endl;
				}
				unsigned int prevSumRank = sumRank;
				sumRank += rank;
				if (sumRank < prevSumRank) {
					overflowRank++;
					if (verbose) cout << " sumRank overflow in shell " <<s<<" at ("<<x<<","<<y<<")" << endl;
				}

				unsigned int prevSumSimSq = sumSimSq;
				sumSimSq += sim*sim;
				if (sumSimSq < prevSumSimSq) {
					overflowSimSq++;
					if (verbose) cout << " sumSimSq overflow in shell " <<s<<" at ("<<x<<","<<y<<")" << endl;
				}
				unsigned int prevSumRankSq = sumRankSq;
				sumRankSq += rank*rank;
				if (sumRankSq < prevSumRankSq) {
					overflowRankSq++;
					if (verbose) cout << " sumRankSq overflow in shell " <<s<<" at ("<<x<<","<<y<<")" << endl;
				}
			}
			if (verbose) cout << endl;
		}

		// Make corrections for overflow, then compute Average and Variance
		double sumSimD = overflowSim*twoPow32 + (double)sumSim;
		double sumRankD = overflowRank*twoPow32 + (double)sumRank;
		double sumSimSqD = overflowSimSq*twoPow32 + (double)sumSimSq;
		double sumRankSqD = overflowRankSq*twoPow32 + (double)sumRankSq;

		double avgSim = sumSimD / numPairs;
		double avgRank = sumRankD / numPairs;
		double stdvSim = sqrt( (sumSimSqD / numPairs) - (avgSim*avgSim) );
		double stdvRank = sqrt( (sumRankSqD / numPairs) - (avgRank*avgRank) );

		grandSumSimD += sumSimD;
		grandSumRankD += sumRankD;
		grandSumSimSqD += sumSimSqD;
		grandSumRankSqD += sumRankSqD;
		grandNumPairs += numPairs;

		// Display results for this shell
		cout << "WITHIN:\t" << s <<"\t"<< shellSize <<"\t"<< numPairs <<"\t";
		cout << setw(12) << sumSimD <<"\t" << avgSim <<"\t"<< stdvSim <<"\t";
		cout << setw(12) << sumRankD <<"\t" << avgRank <<"\t"<< stdvRank << endl;
	}

	// Compute average and variance for all shells together
	double grandAvgSim = grandSumSimD / grandNumPairs;
	double grandAvgRank = grandSumRankD / grandNumPairs;
	double grandStdvSim = sqrt( (grandSumSimSqD / grandNumPairs) - (grandAvgSim*grandAvgSim) );
	double grandStdvRank = sqrt( (grandSumRankSqD / grandNumPairs) - (grandAvgRank*grandAvgRank) );

	// Display results for all shells combined
	cout << "WITHIN:\t" << "ALL" <<"\t"<< "N/A" <<"\t"<< grandNumPairs <<"\t";
	cout << setw(12) << grandSumSimD <<"\t"<< grandAvgSim <<"\t"<< grandStdvSim <<"\t";
	cout << setw(12) << grandSumRankD <<"\t"<< grandAvgRank <<"\t"<< grandStdvRank << endl;
}
//////////////////////////////////////////////////////////////
void Util::computeCrossShellSimilarity(vector< vector<int> >& shells, const vector<int>& simVec,
			const vector<int>& pctileVec, int bitsPrecision, bool verbose) {
	// Recall that simVec is the vectorization of the lower triangle similarity matrix.
	// matrix location: (1,0),(2,0)...(n-1,0),(2,1),(2,1)...(n-1,1),...,(n-1,n-2)
	// Rules: The yth column has (n-1 - y) rows.  Total length = (n choose 2) = n*(n-1)/2.
	// Remember: the diagonal cells are not in the vector.
	// To find any (x,y): (1) Flip x <-> y so that x > y.
	// (2) Count from the end of the vector backwards
	// (because an end-aligned subset of the triangle is just a smaller triangle).
	// The last element (x,y) = (n-1,n-2) ==> simVec(length - 1);
	// If y < n-2, count back a sub-triangle of size (n-2 - y)*(n-1 - y)/2
	// If x < n-1, count back (n-1 - x)
	// Verify: Let n = 10.  Length = (10 choose 2) = 45.
	// (1,0) = (45-1) - [(10-2 - 0)*(10-1 - 0)/2] - [10-1 - 1] = 44 - 36 - 8 = 0
	// (9,0) = (45-1) - [(10-2 - 0)*(10-1 - 0)/2] - [10-1 - 9] = 44 - 36 - 0 = 8
	// (2,1) = (45-1) - [(10-2 - 1)*(10-1 - 1)/2] - [10-1 - 2] = 44 - 28 - 7 = 9
	// (9,8) = (45-1) - [(10-2 - 8)*(10-1 - 8)/2] - [10-1 - 9] = 44 - 0 - 0 = 44
	cout << endl << "STARTING computeCrossShellSimilarity"  << endl;

	// Define: some useful quantities
	int vecMax = simVec.size() - 1;
	int numShells = shells.size();	// skip shells[0]
	int range = (int)pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	double twoPow32 = (float)(1 << 16) * (float)(1 << 16);
	//int numComb = (numShells - 1) * (numShells - 2)/2;	// exclude shells[0]
	
	// Prep: compute nMax = n - 1.  Sort shells.
	int nMax = -1;
	vector< vector<int> >::iterator Sit;
	for (Sit = shells.begin(); Sit != shells.end(); Sit++) {
		nMax += Sit->size();
		sort (Sit->begin(), Sit->end());
	}
	cout << "nMax computed by adding shell sizes = " << nMax << endl;
	nMax = (int)( sqrt(2.0*simVec.size()) + 1) - 1;
	cout << "nMax computed from simVec.size = " << nMax << endl;

	// Define: Class c = set of all block-pairs whose difference of levels (level_y - level_x) = c
	// Special case: Class 0 = all block-pairs
	BlockItem blockItem;
	vector<BlockItem> BlockResults;
	int numClasses = numShells;
	vector<BlockItem> ClassResults(numClasses);
	for (vector<BlockItem>::iterator it = ClassResults.begin(); it != ClassResults.end(); it++) {
		it->s1 = it->s2 = it->sampleSize = it->size = 0;
		it->sumSim = it->sumRank = it->sumSimSq = it->sumRankSq = 0.0;
	}
	
	// Compute: for each pair of shells, compute the average similarity and rank.
	srand(time(NULL));
	int s = 0;

	// print header
	cout << "\tShell1 \t#Nodes1 \t#Shell2 \t#Nodes2 \t\tSumSim\tAvgSim\tStdvSim\t\tSumRank\tAvgRank\tStdvRank" << endl;

	// Iterate through each pair of Shells
	for (int s1 = 0; s1 < numShells; s1++) {

		int size1 = shells[s1].size();
		blockItem.s1 = s1;
		if (size1 == 0) {
			cout << s1 << "\t" << size1 << "\tSKIP" << endl;
			continue;
		}

		for (int s2 = s1 + 1; s2 < numShells; s2++) {

			int size2 = shells[s2].size();
			blockItem.s2 = s2;
			if (size2 == 0) {
				cout << s1 << "\t" << size1 << "\t" << s2 << "\t" << size2 << "\tSKIP" << endl;
				continue;
			}
			blockItem.size = size1*size2;

			// Decide if sampliing is needed
			if (verbose) cout << " = " << blockItem.size << " pairs" << endl;
			bool rand1 = false;
			bool rand2 = false;
			//if (size1 > 1000) {
			//	rand1 = true;
			//	cout << "\t Shell " << s1 << " is too large; sampling 1000" << endl;
			//	size1 = 1000;
			//}
			//if (size2 > 1000) {
			//	rand2 = true;
			//	cout << "\t Shell " << s2 << " is too large; sampling 1000" << endl;
			//	size2 = 2000;
			//}
			blockItem.sampleSize = size1*size2;


			// Add up all pairs of items with the pair of shells
			unsigned int sumSim = 0;
			unsigned int sumRank = 0;
			unsigned int sumSimSq = 0;
			unsigned int sumRankSq = 0;
			unsigned int overflowSim = 0;
			unsigned int overflowRank = 0;
			unsigned int overflowSimSq = 0;
			unsigned int overflowRankSq = 0;

			for (int j = 0; j < size1; j++) {


				int pre_y = (rand1) ? shells[s1].at( rand() % 1000 ) : shells[s1].at(j);

				for (int i = 0; i < size2; i++) {

					int pre_x = (rand2) ? shells[s2].at( rand() % 1000 ) : shells[s2].at(i);

					// Pick x to the the greater of (pre_x, pre_y)
					int x = (pre_x > pre_y) ? pre_x : pre_y;
					int y = (pre_x > pre_y) ? pre_y : pre_x;

					if (verbose) cout << "y=" << y << ": x=";
					int colBack = (nMax - y)*(nMax - 1 - y)/2;
					if (verbose) cout << x << ",";
					int rowBack = (nMax - x);
					int sim = simVec[vecMax - colBack - rowBack];
					if (verbose) cout << "simVec["<<vecMax - colBack - rowBack<<"]; " << endl;

					if (sim < 0 || sim > range) cout << " error: sim = " << sim << " at ("<<x<<","<<y<<")" << endl;
					int rank = findPercentileRank(pctileVec, sim, false);
					if (rank < 0 || rank > 100) cout << " error: rank of sim["<<sim<<"]=" << rank << " at ("<<x<<","<<y<<")" << endl;
					
					// Compute sumSim, sumRank, sumSimSq, and sumRankSq.  Deal with overflow
					unsigned int prevSumSim = sumSim;
					sumSim += sim;
					if (sumSim < prevSumSim) {
						if (verbose) cout << " sumSim overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
						overflowSim++;
					}
					unsigned int prevSumRank = sumRank;
					sumRank += rank;
					if (sumRank < prevSumRank) {
						overflowRank++;
						if (verbose) cout << " sumRank overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
					unsigned int prevSumSimSq = sumSimSq;
					sumSimSq += sim*sim;
					if (sumSimSq < prevSumSimSq) {
						overflowSimSq++;
						if (verbose) cout << " sumSimSq overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
					unsigned int prevSumRankSq = sumRankSq;
					sumRankSq += rank*rank;
					if (sumRankSq < prevSumRankSq) {
						overflowRankSq++;
						if (verbose) cout << " sumRankSq overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
				}
				if (verbose) cout << endl;
			}
			// Make corrections for overflow, then compute Average and Variance
			blockItem.sumSim = sumSim + overflowSim*twoPow32;
			blockItem.sumRank = sumRank + overflowRank*twoPow32;
			blockItem.sumSimSq = sumSimSq + overflowSimSq*twoPow32;
			blockItem.sumRankSq = sumRankSq + overflowRankSq*twoPow32;
			BlockResults.push_back(blockItem);

			double avgSim = blockItem.sumSim / blockItem.sampleSize;
			double avgRank = blockItem.sumRank / blockItem.sampleSize;
			double stdvSim = sqrt( (blockItem.sumSimSq / blockItem.sampleSize) - (avgSim*avgSim) );
			double stdvRank = sqrt( (blockItem.sumRankSq / blockItem.sampleSize) - (avgRank*avgRank) );
			
			// Display Results for this Block-Pair
			//cout << "Shell1\tShell2\t#Nodes1\t#Nodes2\t\tSumSim\tAvgSim\tStdvSim\t\tSumRank\tAvgRank\tStdvRank" << endl;
			cout << "CROSS:\t" << s1 <<"\t"<< size1 <<"\t"<< s2 <<"\t"<< size2 << "\t";
			cout << setw(12) << blockItem.sumSim <<"\t"<< avgSim <<"\t"<< stdvSim << "\t";
			cout << setw(12) << blockItem.sumRank <<"\t"<< avgRank <<"\t"<< stdvRank << endl;

			// Add this shell-pair to the appropriate class
			double sizeRatio = blockItem.size / (double) blockItem.sampleSize;
			ClassResults.at(s2-s1).s2++;					// count the number of shell-pairs
			ClassResults.at(s2-s1).size += blockItem.size;
			ClassResults.at(s2-s1).sumSim += blockItem.sumSim*sizeRatio;
			ClassResults.at(s2-s1).sumRank += blockItem.sumRank*sizeRatio;
			ClassResults.at(s2-s1).sumSimSq += blockItem.sumSimSq*sizeRatio;
			ClassResults.at(s2-s1).sumRankSq += blockItem.sumRankSq*sizeRatio;
			s++;
		}
	}

	// Compute the Global totals for all cross-block-pairs
	for (int c = 1; c < numClasses; c++) {
		ClassResults.at(0).s2 += ClassResults.at(c).s2;
		ClassResults.at(0).size += ClassResults.at(c).size;
		ClassResults.at(0).sumSim += ClassResults.at(c).sumSim;
		ClassResults.at(0).sumRank += ClassResults.at(c).sumRank;
		ClassResults.at(0).sumSimSq += ClassResults.at(c).sumSimSq;
		ClassResults.at(0).sumRankSq += ClassResults.at(c).sumRankSq;
	}

	// Compute: combine results from pairs that have the same level difference
	cout << "GROUPED INTO CLASSES ACCORDING TO SHELL LEVEL DIFFERENCE" << endl;
	cout << "\ts2 - s1\t#Blocks\t#NodePr\tsumSim\t\tAvgSim\tstdvSim\tSumRank\t\tAvgRank\tStdvRank" << endl;
	for (int c = 0; c < numClasses; c++) {

		double avgSim = ClassResults[c].sumSim/ClassResults[c].size;
		double avgRank = ClassResults[c].sumRank/ClassResults[c].size;
		double stdvSim = sqrt( (ClassResults[c].sumSimSq / ClassResults[c].size) - (avgSim*avgSim) );
		double stdvRank = sqrt( (ClassResults[c].sumRankSq / ClassResults[c].size) - (avgRank*avgRank) );

		if (c == 0) {
			cout << "CROSS:\t" << "ALL" << "\t";
		}
		else {
			cout << "CROSS:\t" << c << "\t";
		}
		cout << ClassResults[c].s2 << "\t" << ClassResults[c].size << "\t";
		cout << setw(12) << ClassResults[c].sumSim <<"\t"<< avgSim <<"\t"<< stdvSim << "\t";
		cout << setw(12) << ClassResults[c].sumRank <<"\t"<< avgRank <<"\t"<< stdvRank << endl;
	}

	//cout << "\r\n SORTED RESULTS" << endl;
	//vector<BlockItem>::reverse_iterator bit;

	//cout << "\r\n SORTED BY LEVEL" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareLevel());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}

	//cout << "\r\n SORTED BY AVERAGE RANK" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareAvgRank());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}

	//cout << "\r\n SORTED BY LEVEL DIFFERENCE THEN RANK" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareDiffThenRank());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}
}
//////////////////////////////////////////////////////////////
void Util::computeCrossComponentSimilarity(vector< vector<int> >& bins1, vector< vector<int> >& bins2,									
			const vector<int>& simVec, const vector<int>& pctileVec, int bitsPrecision, bool verbose) {
	// Assumes that the graph is composed of two disconnected components.  Component1 occupies
	// the lower size1 nodes.  This function computes the average similarity rank between each bin
	// of bins1 with each bin of bins 2.
	cout << endl << "STARTING computeCrossCompoenentSimilarity" << endl << endl;

	// Define: some useful quantities
	int vecMax = simVec.size() - 1;
	int numBins = bins1.size();
	int range = (int)pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	double twoPow32 = (float)(1 << 16) * (float)(1 << 16);
	
	// Prep: compute nMax = n - 1.
	int nMax = (int)( sqrt(2.0*simVec.size()) + 1) - 1;
	cout << "nMax computed from simVec.size = " << nMax << endl;

	// Define: Class c = set of all block-pairs whose difference of levels (level_y - level_x) = c
	// Special case: Class 'numBins' = all block-pairs
	BlockItem blockItem;
	vector<BlockItem> BlockResults;
	int numClasses = numBins + 1;
	vector<BlockItem> ClassResults(numClasses);
	for (vector<BlockItem>::iterator it = ClassResults.begin(); it != ClassResults.end(); it++) {
		it->s1 = it->s2 = it->sampleSize = it->size = 0;
		it->sumSim = it->sumRank = it->sumSimSq = it->sumRankSq = 0.0;
	}
	
	// Compute: for each pair of bins, compute the average similarity and rank.
	srand(time(NULL));
	int s = 0;

	// print header
	cout << "Bin1\t#Nodes1\t#Bin2\t#Nodes2\t\tSumSim\tAvgSim\tStdvSim\t\tSumRank\tAvgRank\tStdvRank" << endl;

	// Iterate through each pair of Bins
	for (int s1 = 0; s1 < numBins; s1++) {

		int size1 = bins1[s1].size();
		blockItem.s1 = s1;
		if (size1 == 0) {
			cout << s1 << "\t" << size1 << "\tSKIP" << endl;
			continue;
		}

		for (int s2 = s1; s2 < numBins; s2++) {

			int size2 = bins2[s2].size();
			blockItem.s2 = s2;
			if (size2 == 0) {
				cout << s1 << "\t" << size1 << "\t" << s2 << "\t" << size2 << "\tSKIP" << endl;
				continue;
			}
			blockItem.size = size1*size2;

			// Decide if sampliing is needed
			if (verbose) cout << " = " << blockItem.size << " pairs" << endl;
			bool rand1 = false;
			bool rand2 = false;
			//if (size1 > 1000) {
			//	rand1 = true;
			//	cout << "\t Shell " << s1 << " is too large; sampling 1000" << endl;
			//	size1 = 1000;
			//}
			//if (size2 > 1000) {
			//	rand2 = true;
			//	cout << "\t Shell " << s2 << " is too large; sampling 1000" << endl;
			//	size2 = 2000;
			//}
			blockItem.sampleSize = size1*size2;


			// Add up all pairs of items with the pair of bins
			unsigned int sumSim = 0;
			unsigned int sumRank = 0;
			unsigned int sumSimSq = 0;
			unsigned int sumRankSq = 0;
			unsigned int overflowSim = 0;
			unsigned int overflowRank = 0;
			unsigned int overflowSimSq = 0;
			unsigned int overflowRankSq = 0;

			for (int j = 0; j < size1; j++) {


				int pre_y = (rand1) ? bins1[s1].at( rand() % 1000 ) : bins1[s1].at(j);

				for (int i = 0; i < size2; i++) {

					int pre_x = (rand2) ? bins2[s2].at( rand() % 1000 ) : bins2[s2].at(i);

					// Pick x to be the greater of (pre_x, pre_y)
					int x = (pre_x > pre_y) ? pre_x : pre_y;
					int y = (pre_x > pre_y) ? pre_y : pre_x;

					if (verbose) cout << "y=" << y << ": x=";
					int colBack = (nMax - y)*(nMax - 1 - y)/2;
					if (verbose) cout << x << ",";
					int rowBack = (nMax - x);
					int sim = simVec[vecMax - colBack - rowBack];
					if (verbose) cout << "simVec["<<vecMax - colBack - rowBack<<"]; " << endl;

					if (sim < 0 || sim > range) cout << " error: sim = " << sim << " at ("<<x<<","<<y<<")" << endl;
					int rank = findPercentileRank(pctileVec, sim, false);
					if (rank < 0 || rank > 100) cout << " error: rank of sim["<<sim<<"]=" << rank << " at ("<<x<<","<<y<<")" << endl;
					
					// Compute sumSim, sumRank, sumSimSq, and sumRankSq.  Deal with overflow
					unsigned int prevSumSim = sumSim;
					sumSim += sim;
					if (sumSim < prevSumSim) {
						if (verbose) cout << " sumSim overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
						overflowSim++;
					}
					unsigned int prevSumRank = sumRank;
					sumRank += rank;
					if (sumRank < prevSumRank) {
						overflowRank++;
						if (verbose) cout << " sumRank overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
					unsigned int prevSumSimSq = sumSimSq;
					sumSimSq += sim*sim;
					if (sumSimSq < prevSumSimSq) {
						overflowSimSq++;
						if (verbose) cout << " sumSimSq overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
					unsigned int prevSumRankSq = sumRankSq;
					sumRankSq += rank*rank;
					if (sumRankSq < prevSumRankSq) {
						overflowRankSq++;
						if (verbose) cout << " sumRankSq overflow in shells ("<<s1<<","<<s2<<") at ("<<x<<","<<y<<")" << endl;
					}
				}
				if (verbose) cout << endl;
			}
			// Make corrections for overflow, then compute Average and Variance
			blockItem.sumSim = sumSim + overflowSim*twoPow32;
			blockItem.sumRank = sumRank + overflowRank*twoPow32;
			blockItem.sumSimSq = sumSimSq + overflowSimSq*twoPow32;
			blockItem.sumRankSq = sumRankSq + overflowRankSq*twoPow32;
			BlockResults.push_back(blockItem);

			double avgSim = blockItem.sumSim / blockItem.sampleSize;
			double avgRank = blockItem.sumRank / blockItem.sampleSize;
			double stdvSim = sqrt( (blockItem.sumSimSq / blockItem.sampleSize) - (avgSim*avgSim) );
			double stdvRank = sqrt( (blockItem.sumRankSq / blockItem.sampleSize) - (avgRank*avgRank) );
			
			// Display Results for this Block-Pair
			//cout << "Shell1\tShell2\t#Nodes1\t#Nodes2\t\tSumSim\tAvgSim\tStdvSim\t\tSumRank\tAvgRank\tStdvRank" << endl;
			cout << s1 <<"\t"<< size1 <<"\t"<< s2 <<"\t"<< size2 << "\t";
			cout << setw(12) << blockItem.sumSim <<"\t"<< avgSim <<"\t"<< stdvSim << "\t";
			cout << setw(12) << blockItem.sumRank <<"\t"<< avgRank <<"\t"<< stdvRank << endl;

			// Add this bin-pair to the appropriate class
			double sizeRatio = blockItem.size / (double) blockItem.sampleSize;
			ClassResults.at(s2-s1).s2++;					// count the number of bin-pairs
			ClassResults.at(s2-s1).size += blockItem.size;
			ClassResults.at(s2-s1).sumSim += blockItem.sumSim*sizeRatio;
			ClassResults.at(s2-s1).sumRank += blockItem.sumRank*sizeRatio;
			ClassResults.at(s2-s1).sumSimSq += blockItem.sumSimSq*sizeRatio;
			ClassResults.at(s2-s1).sumRankSq += blockItem.sumRankSq*sizeRatio;
			s++;
		}
	}

	// Compute the Global totals for all cross-block-pairs
	for (int c = 0; c < numBins; c++) {
		ClassResults.at(numBins).s2 += ClassResults.at(c).s2;
		ClassResults.at(numBins).size += ClassResults.at(c).size;
		ClassResults.at(numBins).sumSim += ClassResults.at(c).sumSim;
		ClassResults.at(numBins).sumRank += ClassResults.at(c).sumRank;
		ClassResults.at(numBins).sumSimSq += ClassResults.at(c).sumSimSq;
		ClassResults.at(numBins).sumRankSq += ClassResults.at(c).sumRankSq;
	}

	// Compute: combine results from pairs that have the same level difference
	cout << "GROUPED INTO CLASSES ACCORDING TO SHELL LEVEL DIFFERENCE" << endl;
	cout << "s2 - s1\t#Blocks\t#NodePr\tsumSim\t\tAvgSim\tstdvSim\tSumRank\t\tAvgRank\tStdvRank" << endl;
	for (int c = 0; c < numClasses; c++) {

		double avgSim = ClassResults[c].sumSim/ClassResults[c].size;
		double avgRank = ClassResults[c].sumRank/ClassResults[c].size;
		double stdvSim = sqrt( (ClassResults[c].sumSimSq / ClassResults[c].size) - (avgSim*avgSim) );
		double stdvRank = sqrt( (ClassResults[c].sumRankSq / ClassResults[c].size) - (avgRank*avgRank) );

		if (c == numBins) {
			cout << "ALL" << "\t";
		}
		else {
			cout << c << "\t";
		}
		cout << ClassResults[c].s2 << "\t" << ClassResults[c].size << "\t";
		cout << setw(12) << ClassResults[c].sumSim <<"\t"<< avgSim <<"\t"<< stdvSim << "\t";
		cout << setw(12) << ClassResults[c].sumRank <<"\t"<< avgRank <<"\t"<< stdvRank << endl;
	}

	//cout << "\r\n SORTED RESULTS" << endl;
	//vector<BlockItem>::reverse_iterator bit;

	//cout << "\r\n SORTED BY LEVEL" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareLevel());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}

	//cout << "\r\n SORTED BY AVERAGE RANK" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareAvgRank());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}

	//cout << "\r\n SORTED BY LEVEL DIFFERENCE THEN RANK" << endl;
	//sort(BlockResults.begin(), BlockResults.end(), BlockCompareDiffThenRank());
	//for (bit = BlockResults.rbegin(); bit != BlockResults.rend(); bit++) {
	//	cout <<bit->s1<< "\t" <<bit->s2<< "\t" <<(bit->s2 - bit->s1)<< "\t" <<bit->sumRank/bit->sampleSize << endl;
	//}
}
///////////////////////////////////////////////////////////////////////////
void Util::doShellAnalysis(vector<int>& simVec, vector< vector<int> >& shells, vector<int>& pctileVec, int bitsPrecision, bool verbose) {

	cout << "Starting doShellAnalysis(simVec, shells, " << bitsPrecision << ", " << verbose << ")" << endl;
	histogramSimValues(simVec, 100, bitsPrecision);

	const bool NOT_VERBOSE = false;
	const bool VERBOSE = true;

	// note: findPercentileRank displays its result each time that it is called
	int rank;
	cout << "test findPercentileRank" << endl;
	int range = (int)pow(10.0, (double)bitsPrecision); 
	for (int sim = 0; sim < range; sim = sim + range/500) {
		rank = findPercentileRank(pctileVec, sim, verbose);
	}

	bool includeSelfSim = false;

	vector<int> shellMap = convertShellListsToMap(shells, vsize);
	double K;
	//cout << "Lookup the shell membership of the top what(K) percentage of similary pairs? ";
	K = 0.10;
	cout << "Find the top " << K*100.0 << "% pairs" << endl;
	int numShells = shells.size();
	findShellsOfTopPairs(simVec, pctileVec, shellMap, numShells, K, bitsPrecision);
	computeWithinShellSimilarity(shells, simVec, pctileVec, bitsPrecision, includeSelfSim, verbose);
	computeCrossShellSimilarity(shells, simVec, pctileVec, bitsPrecision, verbose);
}
///////////////////////////////////////////////////////////////////////////
vector<int> Util::convertShellListsToMap(const vector< vector<int> >& shells, int N) {
	// Makes a vector like: { 0 0 0 0 1 1 1 1 1 1 2 2 2 2 3 3 3 ... s-1 s-1 s-1}
	cout << "Starting convertShellListsToMap" << endl;
	vector<int> nodeToShellMap(N);
	for (unsigned int s = 0; s < shells.size(); s++) {
		vector<int>::const_iterator it;
		for (it = shells[s].begin(); it != shells[s].end(); it++) {
			nodeToShellMap[*it] = s;
		}
	}
	return nodeToShellMap;
}
///////////////////////////////////////////////////////////////////////////
void Util::findShellsOfTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
								const vector<int> &shellMap, int numShells, double Kmax, int bitsPrecision) {
	cout << "START findShellsOfTopPairs" << endl;
	// Need to sort, while retaining the pre-sort index location of the top cells.
	// Answer:
	// (1) Read in previously computed percentile ranking (automatic, input parameter)
	// (2) Figure out K = what percentage of the total number of cells
	// (3) Do a pass through the whole matrix: For each cell whose value > threshold
	// percentile (step 2), save the tuple (X, Y, value)
	// (4) Sort this list of top values.
	// (5) For each item in the list, get the shell numbers corresponding to X and Y.
	// (6) Display the value, percentile, X, Y, shell(X), shell(Y) information.
	
	bool verbose = false;
	const int NUM_THRES = 5;
	double Kvalues[NUM_THRES] = {0.0001, 0.001, 0.01, 0.05, 0.10};
	int Kcount[NUM_THRES];					// number of values in the top K percentile
	int thresholdValue;						// threshold sim. value for the top Kmax items
	int vecLen = simVec.size();
	// (2) Figure out K = what percentage of the total number of cells
	// K describes the number of pairs to be reported, in excess of the perfect matches
	// If K >= 1, interpret it as a number of pairs
	// If 0 < K < 1, interpret it as the fraction of the total number of pairs.
	
	for (int i = 0; i < NUM_THRES; i++) {
		Kcount[i] = max( (int)(Kvalues[i]*vecLen), 1);
	}

	if (Kmax > 1) {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax/(double)vecLen ))];
	}
	else {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax ))];
	}
	cout << "thresholdValue for top " << Kmax*100.0 << "% = " << thresholdValue << endl;

	// (3) For each cell whose value > threshold percentile (step 2), save the tuple (X, Y, value)
	vector<MatchItem> topPairs, perfectPairs, nextBestPairs;
	MatchItem m;
	int i = 0;
	for (int y = 0; y < vsize; y++) {	// y,x iteration order MUST match the order of simVec
		for (int x = y + 1; x < vsize; x++) {
			if (simVec[i] >= thresholdValue) {
				m.x = x;
				m.y = y;
				m.weight = (float) simVec[i];
				topPairs.push_back(m);
			}
			i++;
		}
	}
	cout << topPairs.size() << " out of " << vecLen << " items saved" << endl;

	// (4) Sort this list of top values, ***HIGHEST VALUES FIRST***
	sort(topPairs.begin(), topPairs.end(), MatchGreaterThan());

	// (5) For each item in the list, get the shell numbers corresponding to X and Y.

	// Create matrices to store the count of top pairs, indexed by shell
	vector< vector<int> > numPerfectPairs(numShells);
	vector< vector<int> > numClosePairs(numShells);
	vector<int> resetVector(numShells, 0);
	for (int i = 0; i < numShells; i++) {
		numPerfectPairs.at(i) = resetVector;
		numClosePairs.at(i) = resetVector;
	}

	// A. Find the pairs that match perfectly
	int range = (int) pow(10.0, (double)bitsPrecision);
	float PERFECT_VALUE = (float)(0.999 * range);
	int numSameShell = 0;
	int numDiffShell = 0;
	vector<MatchItem>::iterator rit = topPairs.begin();
	while (rit->weight >= PERFECT_VALUE) {
		perfectPairs.push_back(*rit);
		numPerfectPairs[shellMap.at(rit->x)][shellMap.at(rit->y)]++;
		if (shellMap.at(rit->x) == shellMap.at(rit->y)) {
			numSameShell++;
		}
		else if (numDiffShell < 5) {
			numDiffShell++;
			cout << endl << "Unexpected match: x=" <<rit->x << " in shell " <<shellMap.at(rit->x);
			cout << ", y=" <<rit->y << " in shell " <<shellMap.at(rit->y)<< endl;
		}
		rit++;
	}
	cout << endl << perfectPairs.size() << " pairs match perfectly:" << endl;
	cout << "Count how many perfect pairs are within-shell:" << endl;
	cout << "\tValue\tShell\tCount\tFrac of Total" << endl;
	float numPerfPairsF = (float)perfectPairs.size();
	for (int i = 0; i < numShells; i++) {
		cout << "EQUIV:\t" <<PERFECT_VALUE << "\t" << i << "\t" << numPerfectPairs[i][i] << "\t" << numPerfectPairs[i][i]/numPerfPairsF << endl;
	}
	cout << "EQUIV:\t" << PERFECT_VALUE << "\t" << "ALL" << "\t" << numSameShell << "\t" << numSameShell/numPerfPairsF << endl;

	if (verbose) {
		cout << endl << "Value\tX\tY\tShellX\tShellY" << endl;
		for (unsigned int i = 0; i < perfectPairs.size(); i++) {
			m = perfectPairs[i];
			cout <<m.weight<< "\t" <<m.x<< "\t" <<m.y<< "\t" <<shellMap[m.x]<< "\t" <<shellMap[m.y] << endl;
		}
	}

	// B. Record the shell locations of the top K-percentile pairs
	numSameShell = 0;
	int pMax = topPairs.size();
	int Kidx = 0;
	int p = 0;
	cout << endl << "The K% best matching pairs that occur in the same shell:" << endl;
	cout << "\tK% \tKcount \tThrshd \t#SameSh \tFrac. of Total" << endl;
	//for (int p = 0; p < pMax; p++) {
	while (p < pMax && Kidx < NUM_THRES) {
		//nextBestPairs.push_back(topPairs[p]);
		numClosePairs[shellMap.at(topPairs[p].x)][shellMap.at(topPairs[p].y)]++;
		if ( shellMap.at(topPairs[p].x) == shellMap.at(topPairs[p].y) ) {
			numSameShell++;
		}
		if (p == Kcount[Kidx]-1) {
			cout << "TOPK:\t" << Kvalues[Kidx]*100.0 << "\t" << Kcount[Kidx] << "\t" << topPairs[Kcount[Kidx]].weight << "\t";
			cout << numSameShell << "\t" << numSameShell/(float)(p+1) << endl;
			Kidx++;
		}
		p++;
	}

	//while (rit != topPairs.rend()) {
	//	nextBestPairs.push_back(*rit);
	//	numClosePairs[shellMap.at(rit->x)][shellMap.at(rit->y)]++;
	//	if (shellMap.at(rit->x) == shellMap.at(rit->y)) {
	//		numSameShell++;
	//	}
	//	rit++;
	//}

	//cout << "The " << topPairs.size() - perfectPairs.size() << " next best matching pairs:" << endl;
	//cout << "Thrshd\tShell\tCount\t%Total" << endl;
	//float numClosePairsF = (float)nextBestPairs.size();
	//for (int i = 0; i < numShells; i++) {
	//	cout << thresholdValue << "\t" << i << "\t" << numClosePairs[i][i] << "\t" << numClosePairs[i][i]/numClosePairsF << endl;
	//}
	//cout << thresholdValue << "\t" << "ALL" << "\t" << numSameShell << "\t" << numSameShell/numClosePairsF << endl;
}

///////////////////////////////////////////////////////////////////////////
//// CoAuthor Analysis Functions
///////////////////////////////////////////////////////////////////////////
void Util::doCoauthorAnalysis(vector<int>& simVec, vector<AuthorItem> authorVec, vector<int>& pctileVec,
								int bitsPrecision, int size1, bool verbose) {

	cout << "Starting doCoauthorAnalysis(simVec, authorVec, pctile, " << bitsPrecision << ", " << verbose << ")" << endl;
	histogramSimValues(simVec, 100, bitsPrecision);

	const bool NOT_VERBOSE = false;
	const bool VERBOSE = true;
	bool includeSelfSim = false;

	//vector<int> shellMap = convertShellListsToMap(shells, vsize);
	vector< vector<int> > authBins1 = binAuthorsSublistByRank(authorVec, 0, size1);
	vector< vector<int> > authBins2 = binAuthorsSublistByRank(authorVec, size1, authorVec.size() - size1);

	findAuthorRankOfTopPairs(simVec, pctileVec, authorVec, bitsPrecision, size1);

	cout << endl << "*** Within-Bin: First Set of Authors ***" << endl;
	computeWithinShellSimilarity(authBins1, simVec, pctileVec, bitsPrecision, includeSelfSim, verbose);
	cout << endl << "*** Within-Bin: Second Set of Authors ***" << endl;
	computeWithinShellSimilarity(authBins2, simVec, pctileVec, bitsPrecision, includeSelfSim, verbose);

	cout << endl << "*** Cross-Bin: First Set of Authors ***" << endl;
	computeCrossShellSimilarity(authBins1, simVec, pctileVec, bitsPrecision, verbose);
	cout << endl << "*** Cross-Bin: Second Set of Authors ***" << endl;
	computeCrossShellSimilarity(authBins2, simVec, pctileVec, bitsPrecision, verbose);
	cout << endl << "*** Cross-Component: First Set of Authors' Bins vs Second Set of Authors' Bins ***" << endl;
	computeCrossComponentSimilarity(authBins1, authBins2, simVec, pctileVec, bitsPrecision, verbose);
}
///////////////////////////////////////////////////////////////////////////
vector< vector<int> > Util::binAuthorsSublistByRank(const vector<AuthorItem>& authorVec, int startIndex, int subLength) {
	vector< vector<int> > authorBins(10);
	cout << "Starting binAuthorsSublistByRank, startIndex = " << startIndex << ", subLength = " << subLength << endl;

	cout << "authorVec.size() = " << authorVec.size() << endl;
	
	for (int i = startIndex; i < startIndex+subLength; i++) {
		if (i >= authorVec.size()) {
			cout << "Out-of-range Error: i = " << i << endl;
		}
		int bin = min((int)( authorVec.at(i).Grank * 10.0 ), 9);
		if (bin > 9) {
			cout << "Out-of-range Error: authorVec.at("<<i<<").Grank = " << authorVec.at(i).Grank << endl;
		}
		authorBins.at(bin).push_back(authorVec.at(i).id);
	}
	// report sizes
	cout << "Bin\tSize" << endl;
	for (int b = 0; b < 10; b++) {
		cout << b << "\t" << authorBins[b].size() << endl;
	}
	return authorBins;
}
///////////////////////////////////////////////////////////////////////////
vector<AuthorItem> Util::readAuthorRanks(char* rankFileName){
	cout << "Starting readAuthorRanks(" << rankFileName << ")" << endl;
	AuthorItem item;
	vector<AuthorItem> itemList;
	ifstream in(rankFileName);
	while (!in.eof()) {
		if (in >> item.id) {
			in >> item.Grank >> item.Hrank;
			itemList.push_back(item);
		}
	}
	cout << itemList.size() << " items read" << endl;
	return itemList;
}
///////////////////////////////////////////////////////////////////////////
void Util::findAuthorRankOfTopPairs(const vector<int>& simVec, const vector<int>& pctileVec,
								const vector<AuthorItem> authorVec, int bitsPrecision, int size1) {
	cout << "START findAuthorRankOfTopPairs" << endl;
	// Give a value of K,
	// (1) Read in previously computed percentile ranking (automatic, input parameter)
	// (2) Figure out the threshold similarity value corresponds to the top K% of cells
	// (3) Do a pass through the whole matrix: For each cell whose value > threshold(s)
	// percentile (step 2), save the tuple (X, Y, value)
	// (4) Sort this list of top values.
	// (5) For each top cell, look up X.Grank, Y.Grank, and compute the difference deltaGrank.
	// Do the same for Hrank.
	// (6) Compute the average and std deviation for the group of top cells.
	
	bool verbose = false;
	const int NUM_THRES = 6;
	double Kvalues[NUM_THRES] = {0.00001, 0.0001, 0.001, 0.01, 0.05, 0.10};
	int Kcount[NUM_THRES];					// number of values in the top K percentile
	int thresholdValue;						// threshold sim. value for the top Kmax items
	int vecLen = simVec.size();
	double Kmax = Kvalues[NUM_THRES-1];		// max. K% value
	// (2) Figure out K = what percentage of the total number of cells
	// K describes the number of pairs to be reported, in excess of the perfect matches
	// If K >= 1, interpret it as a number of pairs
	// If 0 < K < 1, interpret it as the fraction of the total number of pairs.

	for (int i = 0; i < NUM_THRES; i++) {
		Kcount[i] = max( (int)(Kvalues[i]*vecLen), 1);
	}

	if (Kmax > 1) {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax/(double)vecLen ))];
	}
	else {
		thresholdValue = pctileVec[(int)(100.0*( 1.0 - Kmax ))];
	}
	cout << "thresholdValue for top " << Kmax*100.0 << "% = " << thresholdValue << endl;

	// (3) For each cell whose value > threshold percentile (step 2), save the tuple (X, Y, value)
	vector<MatchItem> topPairs, perfectPairs, nextBestPairs;
	MatchItem m;
	int i = 0;
	for (int y = 0; y < vsize; y++) {	// y,x iteration order MUST match the order of simVec
		for (int x = y + 1; x < vsize; x++) {
			if (simVec[i] >= thresholdValue) {
				m.x = x;
				m.y = y;
				m.weight = (float) simVec[i];
				topPairs.push_back(m);
			}
			i++;
		}
	}
	cout << topPairs.size() << " out of " << vecLen << " items saved" << endl;

	// (4) Sort this list of top values, ***HIGHEST VALUES FIRST***
	sort(topPairs.begin(), topPairs.end(), MatchGreaterThan());

	// (5) For each top cell, look up X.Grank, Y.Grank, and compute the difference deltaGrank.

	vector<double> XYmid;
	vector<double> XYdiff;
	vector<double> XYdiffSq;

	double sumXYmid = 0.0;
	double sumXYdiff = 0.0;
	double sumXYdiffSq = 0.0;
	int count = 0;
	int countDiffComp = 0;

	// A. Find the pairs that match perfectly
	int range = (int) pow(10.0, (double)bitsPrecision);
	float PERFECT_VALUE = (float)(0.999 * range);
	
	if (verbose) {
		cout << "Value\tMid\tDiff" << endl;
	}
	vector<MatchItem>::iterator rit = topPairs.begin();
	while (rit->weight >= PERFECT_VALUE) {
		perfectPairs.push_back(*rit);
		if (rit->x >= size1 && rit->y < size1 || rit->x < size1 && rit->y >= size1) {
			countDiffComp++;
		}
		double xGrank = authorVec.at(rit->x).Grank;
		double yGrank = authorVec.at(rit->y).Grank;
		sumXYmid += (xGrank + yGrank)/2;
		sumXYdiff += abs(xGrank - yGrank);
		sumXYdiffSq += (xGrank- yGrank)*(xGrank - yGrank);
		if (verbose) {
			cout << rit->weight << "\t" << (xGrank + yGrank)/2 << "\t" << abs(xGrank - yGrank);
		}
		rit++;
		count++;
	}
	//cout << count << "pairs match perfectly:" << endl;
	cout << "K% \t#Pairs \tThrshd \t#DifCmp\t%DifCmp\tAvgRank\tAvgDiff\tStdvDiff" << endl;

	// Print totals
	double avgXYdiff = sumXYdiff/count;
	cout << setprecision(5)
		<< "Perfect" << "\t"  << count << "\t" << rit->weight << "\t"
		<< countDiffComp << "\t" << countDiffComp*100.0/count << "\t"
		<< sumXYmid/count << "\t" << avgXYdiff << "\t"
		<< sqrt( (sumXYdiffSq/count) - (avgXYdiff*avgXYdiff) ) << endl;

	// B. Compute values for the top K pairs, for a set of values of K
	sumXYmid = 0.0;
	sumXYdiff = 0.0;
	sumXYdiffSq = 0.0;
	count = 0;
	countDiffComp = 0;

	int pMax = topPairs.size();
	int Kidx = 0;
	cout << endl << "G-index difference for the K% best matching pairs:" << endl;
	cout << "K% \tKcount \tThrshd \t#SameSh \tFrac. of Total" << endl;
	for (rit = topPairs.begin(); rit != topPairs.end(); rit++) {
		double xGrank = authorVec.at(rit->x).Grank;
		double yGrank = authorVec.at(rit->y).Grank;
		sumXYmid += (xGrank + yGrank)/2;
		sumXYdiff += abs(xGrank - yGrank);
		sumXYdiffSq += (xGrank- yGrank)*(xGrank - yGrank);
		count++;
		if (rit->x >= size1 && rit->y < size1 || rit->x < size1 && rit->y >= size1) {
			countDiffComp++;
		}
		//if (verbose) {
		//	cout << rit->weight << "\t" << (xGrank + yGrank)/2 << "\t" << abs(xGrank - yGrank);
		//}
		if (count == Kcount[Kidx]-1) {
			// Print totals
			double avgXYdiff = sumXYdiff/count;
			cout << setprecision(5)
				<< Kvalues[Kidx]*100 << "\t" << Kcount[Kidx] << "\t" << rit->weight << "\t"
				<< countDiffComp << "\t" << countDiffComp*100.0/count << "\t"
				<< sumXYmid/count << "\t" << avgXYdiff << "\t"
				<< sqrt( (sumXYdiffSq/count) - (avgXYdiff*avgXYdiff) ) << endl;
			Kidx++;
		}
	}
	cout << "count = " << count << ", Kcount["<<Kidx<<"]-1 = " << Kcount[NUM_THRES-1] -1 << endl;
}
///////////////////////////////////////////////////////////////////////////
void Util::rankMatrixValues(const vector<int>& simVec, int bitsPrecision, map<int,float>& matPctRank)
{
	cout << "rankMatrixValues" << endl;

	// 1. Scan the matrix's simVec, counting how many times each sim value occurs.
	// That is, build a map whose key->value are sim_value->count.
	if(verbose) cout << "  Scan simVec; count instances of each value" << endl;
	int matVal;
	int N = (int)( sqrt(2.0*simVec.size()) + 1);
	map<int,int> matCount;
	// Iterate through x and y in the particular ordering of the simVec
	int i = 0;
	for (int y = 0; y < N; y++) {
		for (int x = y+1; x < N; x++) {
			matVal = simVec[i];
			if (matCount.find(matVal) == matCount.end()) {
				if(verbose) cout << "matCount[" << matVal << "]=1" << endl;
				matCount[matVal] = 1;
			}
			else {
				if(verbose) cout << "matCount[" << matVal << "]=" << matCount[matVal] + 1 << endl;
				matCount[matVal]++;
			}
			i++;
		}
	}

	// 2. Sort the maps (unncessary; C++ STL map is always sorted by key)
	// 3. Assign a percentile rank to each individual map item. If multiple entries have the same value, then
	//		assign all of them the average rank = (rank.first_in_set + rank.last_in_set)/2
	int startIdx, endIdx;
	int vecSize = simVec.size();
	double matRankScale = 1.0/vecSize;
	double pctRank;
	int val;

	// Make the matPctRank table
	startIdx = 0;
	endIdx = 0;
	cout << "  Compute percentile rank of each matrix value" << endl;
	if (verbose) printf(" %7s %7s %7s %7s %10s\n", "Value", "Count", "StartIdx", "StartIdx", "AvgPctRank");
	map<int,int>::iterator cItr;
	for (cItr = matCount.begin(); cItr != matCount.end(); cItr++) {
		endIdx = startIdx + cItr->second;
		if (endIdx >= vecSize) {
			pctRank = 1.0; // (int) range;			// Cells with maximum similarit are assigned rank 1
		}
		else {
			pctRank = ( matRankScale * (endIdx + startIdx + 1) / 2.0);
		}
		val = cItr->first;
		matPctRank[val] =(float) pctRank;
		if (verbose) printf(" %7d %7d %7d %7d  %7.5f\n", cItr->first, cItr->second, startIdx, endIdx, matPctRank[val]);
		startIdx = endIdx;
	}
}
///////////////////////////////////////////////////////////////////////////
void Util::rankIcebergValues( IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
		map<int,float>& icePctRank, map<int,float>& matPctRank)
{
	// map<int,float> = <similarity_value, percentile_rank>
	cout << "rankIcebergValues" << endl;
	icePctRank.clear();
	matPctRank.clear();

	int N = iceberg.getN();
	double range = pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	SimMapType iceMap = iceberg.getMap();
	int mapSize = iceMap.size();

	// 1. Scan the iceberg Map, counting how many times each sim value occurs.
	// That is, build a map whose key->value are sim_value->count.
	if(verbose) cout << "  Scan iceberg table and matrix subset; count instances of each value" << endl;
	int iceVal, matVal;
	int x, y;
	int numDiagonal = 0;
	map<int,int> iceCount, matCount;
	SimMapType::iterator mItr;
	for (mItr = iceMap.begin(); mItr != iceMap.end(); mItr++) {
		x = mItr->first.first;
		y = mItr->first.second;

		// don't bother with the diagonal cells
		if (x == y) {
			numDiagonal++;
			continue;
		}
		if(verbose) cout << "x=" << x << ",y=" << y << ": ";

		// find the iceberg value at (x,y) and increment iceCount
		iceVal = (int)(mItr->second * range);
		if (iceCount.find(iceVal) == iceCount.end()) {
			if (verbose) cout << "iceCount[" << iceVal << "]=1, ";
			iceCount[iceVal] = 1;
		}
		else {
			if (verbose) cout << "iceCount[" << iceVal << "]=" << iceCount[iceVal] + 1 << ", ";
			iceCount[iceVal]++;
		}

		// find the simVec value for (x,y) and increment the matCount
		matVal = getVal(simVec, x, y, N);
		if (matCount.find(matVal) == matCount.end()) {
			if(verbose) cout << "matCount[" << matVal << "]=1" << endl;
			matCount[matVal] = 1;
		}
		else {
			if(verbose) cout << "matCount[" << matVal << "]=" << matCount[matVal] + 1 << endl;
			matCount[matVal]++;
		}
	}

	// 2. Sort the maps (unncessary; C++ STL map is always sorted by key)
	// 3. Assign a percentile rank to each individual map item. If multiple entries have the same value, then
	//		assign all of them the average rank = (rank.first_in_set + rank.last_in_set)/2
	int startIdx, endIdx;
	float pctRank;
	int val;
	map<int,int>::iterator vItr;
	// Make the icePctRank table
	startIdx = 0;
	endIdx = 0;
	float iceRankScale = 1 / (float)(mapSize - numDiagonal); // old: range / mapSize;

	cout << "  Compute percentile rank of each iceberg value" << endl;
	if (verbose) printf(" %7s %7s %7s %7s %10s\n", "Value", "Count", "StartIdx", "EndIdx", "AvgPctRank");
	for (vItr = iceCount.begin(); vItr != iceCount.end(); vItr++) {
		endIdx = startIdx + vItr->second;
		if (endIdx >= mapSize) {
			pctRank = 1.0f; //(int) range;			// Cells with maximum similarity are assigned maximum rank
		}
		else {
			pctRank = ( iceRankScale * (endIdx + startIdx + 1) / 2.0f );
		}
		val = vItr->first;
		icePctRank[val] = pctRank;
		if (verbose) printf(" %7d %7d %7d %7d  %7.5f\n", vItr->first, vItr->second, startIdx, endIdx, icePctRank[val]);
		startIdx = endIdx;
	}
	// Make the matPctRank table
	startIdx = 0;
	endIdx = 0;
	cout << "  Compute percentile rank of each corresponding matrix value" << endl;
	if (verbose) printf(" %7s %7s %7s %7s %10s\n", "Value", "Count", "StartIdx", "StartIdx", "AvgPctRank");
	for (vItr = matCount.begin(); vItr != matCount.end(); vItr++) {
		endIdx = startIdx + vItr->second;
		if (endIdx >= mapSize) {
			pctRank = 1.0; // (int) range;			// Cells with maximum similarit are assigned rank 1
		}
		else {
			pctRank = ( iceRankScale * (endIdx + startIdx + 1) / 2.0);
		}
		val = vItr->first;
		matPctRank[val] = pctRank;
		if (verbose) printf(" %7d %7d %7d %7d  %7.5f\n", vItr->first, vItr->second, startIdx, endIdx, matPctRank[val]);
		startIdx = endIdx;
	}
}
///////////////////////////////////////////////////////////////////
void Util:: compareTwoRankings(const vector<int>& simVec1, const vector<int>& simVec2,
			 map<int,float>& rank1, map<int,float>& rank2, int bitsPrecision)
{
	cout << "Starting compareTwoRankings" << endl;
	// Compute Pearson correlation to compare rankings
	// Also compute average, min, and max difference.
	double sumRank1 = 0;
	double sumRank2 = 0;
	double delta1, delta2;
	// To minimize precision errors from adding lots of small numbers, sort the deltas first,
	// then add the smallest ones together first.
	priority_queue<double> prodDeltaList, sqDelta1List, sqDelta2List;
	double sumProductDelta = 0;
	double sumSqDelta1 = 0;
	double sumSqDelta2 = 0;

	double diff;
	double maxDiff = 0;
	double sumDiff = 0;

	int N = simVec1.size();

	// Psss 0: compute average rankings (expected to be range/2
	//int meanRank =(int) pow(10.0, (double)bitsPrecision)/2 ;
	//for (int i = 0; i < N; i++) {
	//	sumRank1 += rank1[simVec1[i]];
	//	sumRank2 += rank2[simVec2[i]];
	//}
	//avgRank1 = sumRank1 / N;
	//avgRank2 = sumRank2 / N;
	//cout << "COMP: avgRank1 = " << avgRank1 << endl;
	//cout << "COMP: avgRank2 = " << avgRank2 << endl;

	// Pass 1: compute delta values and sq-delta values; store them in priority queue.
	
	for (int i = 0; i < N; i++) {
		delta1 = abs(rank1[simVec1[i]] - 0.5);
		delta2 = abs(rank2[simVec2[i]] - 0.5);
		prodDeltaList.push(delta1 * delta2);
		sqDelta1List.push(delta1 * delta1);
		sqDelta2List.push(delta2 * delta2);

		diff = abs(rank1[simVec1[i]] - rank2[simVec2[i]]);
		if (diff > maxDiff) {
			maxDiff = diff;
		}
		sumDiff += diff;
	}

	// Psss 2: sum up values in priority queues
	for (int i = 0; i < N; i++) {
		sumProductDelta += prodDeltaList.top();
		prodDeltaList.pop();
		sumSqDelta1 += sqDelta1List.top();
		sqDelta1List.pop();
		sumSqDelta2 += sqDelta2List.top();
		sqDelta2List.pop();
	}

	double pearson = sumProductDelta / sqrt(sumSqDelta1 * sumSqDelta2);
	double avgDiff = sumDiff / N;
	cout << "  (Percentile) Rank values range from 0.00 to 1.00" << endl;
	cout << "COMP: avg Rank difference = " << avgDiff << endl;
	cout << "COMP: max Rank difference = " << maxDiff << endl;
	cout << "COMP: Pearson correlation coefficient = " << pearson << endl;
}
///////////////////////////////////////////////////////////////////
void Util::compareIcebergToSimMatrix(IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
							  map<int,float>& iceRankTable,  map<int,float>& matRankTable) {
	cout << endl << "Start compareIcebergToSimMatrix()" << endl;
	cout << "1. Compare the full set of iceberg values to regular similarity matrix values" << endl;
	/*
	1. For the set of node-pairs P that are in the iceberg:
		Rank the iceberg's set of values, and separately rank the matrix's corresponding set of values.
		a. Avg. and Max. absolute difference (full matrix value - iceberg value)
		c. Avg. and Max. relative difference (full matrix value - iceberg value) / full matrix value
		b. Avg. and Max. absolute difference (full matrix rank - iceberg rank)
	*/

	int x, y;
	int xMaxValDiff, xMaxRelValDiff, xMaxRankDiff;
	int yMaxValDiff, yMaxRelValDiff, yMaxRankDiff;
	int iceVal, matVal;
	float iceRank, matRank;

	int valDiff;
	int maxValDiff = -1;

	double rankDiff;
	double maxRankDiff = -1;;

	double relValDiff, relRankDiff;
	double maxRelValDiff = -1;
	double maxRelRankDiff = -1;

	double totValDiff = 0;
	double totRelValDiff = 0;
	double totRankDiff = 0;
	double totRelRankDiff = 0;

	int eqBoth = 0;
	int eqMatNotIce = 0;
	int eqIceNotMat = 0;
	
	// Variables used by the Pearson correlation computation
	double sumRank1 = 0;
	double sumRank2 = 0;
	double delta1, delta2;
	//priority_queue<double> prodDeltaList, sqDelta1List, sqDelta2List;
	vector<double> prodDeltaList, sqDelta1List, sqDelta2List;
	double sumProductDelta = 0;
	double sumSqDelta1 = 0;
	double sumSqDelta2 = 0;

	double diff;
	double maxDiff = 0;
	double sumDiff = 0;
	// End - Variables used by the Pearson correlation computation

	double rangeD = pow(10.0, (double)bitsPrecision); // range of values stored in simVec
	int range = (int) rangeD;
	int meanRank = range / 2;
	SimMapType iceMap = iceberg.getMap();
	int mapSize = iceMap.size();
	int N = iceberg.getN();

	cout << "iceberg contains " << mapSize << " items = " << (100.0 * mapSize)/(simVec.size() + N) << "% of the full matrix" << endl;
	cout << "To interpret absolute value differences:" << endl;
	cout << "   The range for both similarity values and for percentile ranks is [0, 1]" << endl; 
	int numComparisons = 0;
	int numDiagonal = 0;
		
	vector<int> breakpoints;
	breakpoints.push_back(1000);
	breakpoints.push_back(10000);
	breakpoints.push_back(100000);
	int quintile = (int) (mapSize/5.0);
	breakpoints.push_back(quintile*1);
	breakpoints.push_back(quintile*2);
	breakpoints.push_back(quintile*3);
	breakpoints.push_back(quintile*4);
	breakpoints.push_back(mapSize);
	sort(breakpoints.begin(), breakpoints.end());
	int Q = 0;

	SimMapType::iterator mItr;
	for (mItr = iceMap.begin(); mItr != iceMap.end(); mItr++) {
		x = mItr->first.first;
		y = mItr->first.second;

		// don't bother with the diagonal cells
		if (x == y) {
			numDiagonal++;
			continue;
		}
		numComparisons++;

		// a. Get values and ranks
		iceVal = (int)(mItr->second * rangeD);
		iceRank = iceRankTable[iceVal];
		matVal = getVal(simVec, x, y, N);
		matRank = matRankTable[matVal];

		// b. Update various max values
		valDiff = abs(iceVal - matVal);
		if (valDiff > maxValDiff) {
			maxValDiff = valDiff;
			xMaxValDiff = x;
			yMaxValDiff = y;
		}

		relValDiff = valDiff / (double)matVal;
		if (relValDiff > maxRelValDiff) {
			maxRelValDiff = relValDiff;
			xMaxRelValDiff = x;
			yMaxRelValDiff = y;
		}

		rankDiff = abs(iceRank - matRank);
		if (rankDiff > maxRankDiff) {
			maxRankDiff = rankDiff;
			xMaxRankDiff = x;
			yMaxRankDiff = y;
		}

		relRankDiff = rankDiff / (double)matRank;
		if (relRankDiff > maxRelRankDiff) {
			maxRelRankDiff = relRankDiff;
		}

		totValDiff += valDiff;
		totRankDiff += rankDiff;
		totRelValDiff += relValDiff;
		totRelRankDiff += relRankDiff;

		// c. Count occurences of equiv. node-pairs which have the max. possible similarity value
		if (matVal == range && iceVal == range) {
			eqBoth++;
		}
		else if (matVal == range && iceVal != range) {
			eqMatNotIce++;
		}
		else if (matVal != range && iceVal == range) {
			eqIceNotMat++;
		}
		// d. Pearson correlation
		delta1 = abs(matRank - 0.5);
		delta2 = abs(iceRank - 0.5);
		prodDeltaList.push_back(delta1 * delta2);
		sqDelta1List.push_back(delta1 * delta1);
		sqDelta2List.push_back(delta2 * delta2);

		diff = abs(matRank - iceRank);
		if (diff > maxDiff) {
			maxDiff = diff;
		}
		sumDiff += diff;

		if (numComparisons == breakpoints[Q] ) {
			Q++;
			cout << endl;
			cout <<Q<< "\tICE: For_values_equal_or_above " << iceVal << endl;
			cout << "numComparions = " << numComparisons << endl;

			double numCompD = (double) numComparisons;
			double avgValDiff = totValDiff / numCompD;
			double avgRankDiff = totRankDiff / numCompD;
			double avgRelValDiff = totRelValDiff / numCompD;
			double avgRelRankDiff = totRelRankDiff / numCompD;
			     
			cout << endl;
			cout <<Q<< "\tICE: Avg_difference_of_values " << avgValDiff / rangeD << endl;
			cout <<Q<< "\tICE: Max_difference_of_values " << maxValDiff / rangeD << endl;
			//cout <<Q<< "\tICE: Avg_difference_of_ranks " << avgRankDiff << endl;
			//cout <<Q<< "\tICE: Max_difference_of_ranks " << maxRankDiff << endl;
			cout << endl;
			cout <<Q<< "\tICE: Avg_relative_difference_of_values = " << avgRelValDiff << endl;
			cout <<Q<< "\tICE: Max_relative_difference_of_values = " << maxRelValDiff << endl;
			//cout <<Q<< "\tICE: Avg_relative_difference_of_ranks = " << avgRelRankDiff << endl;
			//cout <<Q<< "\tICE: Max_relative_difference_of_ranks = " << maxRelRankDiff << endl;
			cout << endl;
			cout <<Q<< "\tICE: Max abs val diff is at (" << xMaxValDiff << ", " << yMaxValDiff << ") ";
			cout << "Matrix val = " << getVal(simVec, xMaxValDiff, yMaxValDiff, N)/rangeD;
			cout << ", iceberg val = " << iceberg.val(xMaxValDiff, yMaxValDiff) << endl;

			cout <<Q<< "\tICE: Max rel val diff is at (" << xMaxRelValDiff << ", " << yMaxRelValDiff << ") ";
			cout << "Matrix val = " << getVal(simVec, xMaxRelValDiff, yMaxRelValDiff, N)/rangeD;
			cout << ", iceberg val = " << iceberg.val(xMaxRelValDiff, yMaxRelValDiff) << endl;

			cout <<Q<< "\tICE: Max rank diff is at (" << xMaxRankDiff << ", " << yMaxRankDiff << ") ";
			cout << "Matrix rank = " << matRankTable[ getVal(simVec, xMaxRankDiff, yMaxRankDiff, N) ];
			cout << ", iceberg rank = " << iceRankTable[(int)( iceberg.val(xMaxRankDiff, yMaxRankDiff)*rangeD) ] << endl;
			cout << endl;
			cout <<Q<< "\tICE: " << eqBoth << " node-pairs are auto. equivalent in both iceberg and simMatrix (including self-similarity)" << endl;
			cout <<Q<< "\tICE: " << eqMatNotIce << " node-pairs are auto. equivalent in simMatrix but not in iceberg" << endl;
			cout <<Q<< "\tICE: " << eqIceNotMat << " node-pairs are auto. equivalent in iceberg but not in simMatrix" << endl;

			// Psss 2: Pearson correlation: sum up values in priority queues
			for (int i = 0; i < numComparisons; i++) {		// numComparisons, not N or mapSize
				sumProductDelta += prodDeltaList.at(i);
				//prodDeltaList.pop();
				sumSqDelta1 += sqDelta1List.at(i);
				//sqDelta1List.pop();
				sumSqDelta2 += sqDelta2List.at(i);
				//sqDelta2List.pop();
			}

			double pearson = sumProductDelta / sqrt(sumSqDelta1 * sumSqDelta2);
			double avgDiff = sumDiff / numComparisons;	// numComparisons, not N or mapSize
			cout << "  (Percentile) Rank values range from 0.00 to 1.00" << endl;
			cout <<Q<< "\tICE: avg Rank difference = " << avgDiff << endl;
			cout <<Q<< "\tICE: max Rank difference = " << maxDiff << endl;
			cout <<Q<< "\tICE: Pearson correlation coefficient = " << pearson << endl;
		}
	}
	cout << "  Number of comparisons = " << numComparisons << endl;
	if (numComparisons != (mapSize - numDiagonal)) {
		cout << "    but expected " << (mapSize - numDiagonal ) << " comparisons" << endl;
	}


}
///////////////////////////////////////////////////////////////////
void Util::checkIcebergCoverage(IcebergSimMap& iceberg, const vector<int>& simVec, int bitsPrecision,
				 map<int,float>& iceRank,  map<int,float>& matRank)
{
	cout << endl << "Start checkIcebergCoverage()" << endl;
	cout << "1. Pruning accuracy: what fraction of top matrix values are included in iceberg?" << endl;
	cout << "2. Computing accuracy: for top cells, how close are iceberg values to matrix values?" << endl;

	/*
	2. Sort the matrix node-pairs from high to low values (most to least similar)
		a. What fraction of the matrix's top node-pairs (value > theta) also appear in the iceberg?
		b. What fraction of the matrix's top node-pairs (top K%) also appear in the iceberg?
		c. For those top values that appear in both sets, avg. and max. absolute difference?
		d. For those top values that appear in both sets, avg. and max. relative difference?
	*/
	// Define the thresholds to be used. Note: simVec values are integers on a scale of 0 to 1000
	int vecSize = simVec.size();
	int N = (int)( sqrt(2.0*simVec.size()) + 1);
	double rangeD = pow(10.0, (double)bitsPrecision); // range of values stored in simVec

	const int MAX = 5;
	int theta[MAX]     = {1000, 900, 800, 700, 600};
	double topPct[MAX] = {0.0001, 0.0003, 0.001, 0.003, 0.01};
	int topVal[MAX];

	// Find the values that correspond to the topRank thresholds
	// Use the fact that C++ STL maps are automatically sorted by index value
	// topRank is the inverse of rank: topPct = 1 - rank
	map<int,float>::reverse_iterator mItr = matRank.rbegin();
	int k = 0;
	while (k < MAX) {
		double limit = 1 - topPct[k];
		while ( mItr->second > limit) {
			mItr++;
		}
		topVal[k]  = mItr->first;
		k++;
	}
	// Display the discovered topKvalues
	cout << endl << "Similarity values corresponding to the top K fraction" << endl;
	cout << "K \t Value" << endl;
	for (int i = 0; i < MAX; i++) {
		cout << topPct[i] << "\t" << topVal[i] << endl;
	}

	// Scan the simVec, and log any values that exceed the theta or topVal thresholds

	vector<int> thetaCountMat(MAX, 0);
	vector<int> thetaCountIce(MAX, 0);
	vector<int> kCountMat(MAX, 0);
	vector<int> kCountIce(MAX, 0);

	int matVal, iceVal;
	int valDiff;
	vector<int> maxDiffTheta(MAX, 0);
	vector<int> totDiffTheta(MAX, 0);
	vector<int> maxDiffK(MAX, 0);
	vector<int> totDiffK(MAX, 0);
	double relDiff;
	vector<double> maxRelDiffTheta(MAX, 0);
	vector<double> totRelDiffTheta(MAX, 0);
	vector<double> maxRelDiffK(MAX, 0);
	vector<double> totRelDiffK(MAX, 0);

	// Iterate through x and y in the particular ordering of the simVec
	bool iceFound = false;
	SimMapType iceMap = iceberg.getMap();
	SimMapType::iterator iceItr;

	for (int y = 0; y < N; y++) {
		for (int x = y+1; x < N; x++) {
			matVal = getVal(simVec, x, y, N);
			// If the value exceeds a lowest threshold
			if (matVal >= theta[MAX -1] || matVal >= topVal[MAX - 1]) {
				// See if this node-pair is in the iceberg. If so, calc the difference of values
				iceItr = iceMap.find(make_pair(x,y));
				if (iceItr != iceMap.end()) {
					iceFound = true;
					iceVal = (int)(iceItr->second * rangeD);
					valDiff = abs(matVal - iceVal);
					relDiff = valDiff / (double) matVal;
				}
				// Check all the thresholds and tally 
				// 1. Check the theta thresholds
				int t = MAX - 1;
				while (matVal >= theta[t] && t >= 0) {
					thetaCountMat[t]++;
					if (iceFound) {
						thetaCountIce[t]++;
						if (valDiff > maxDiffTheta[t]) {
							maxDiffTheta[t] = valDiff;
						}
						if (relDiff > maxRelDiffTheta[t]) {
							maxRelDiffTheta[t] = relDiff;
						}
						totDiffTheta[t] += valDiff;
						totRelDiffTheta[t] += relDiff;
					}
					t--;
				}
				// 2. Check the topVal thresholds
				int k = MAX - 1;
				while (matVal >= topVal[k] && k >= 0) {
					kCountMat[k]++;
					if (iceFound) {
						kCountIce[k]++;
						if (valDiff > maxDiffK[k]) {
							maxDiffK[k] = valDiff;
						}
						if (relDiff > maxRelDiffK[k]) {
							maxRelDiffK[k] = relDiff;
						}
						totDiffK[k] += valDiff;
						totRelDiffK[k] += relDiff;
					}
					k--;
				}
				iceFound = false;
			} // if matVal
		} // x
	} // y

	// Report statistics
	double thres;
	int numMat, numIce;
	float hitRate;
	double avgDiff, maxDiff, avgRelDiff, maxRelDiff;

	cout << endl << "### Iceberg Coverage & Accuracy for values above a Threshold ###" << endl;
	cout << "Theta   #Matrix #Ice    HitRate avgDiff maxDiff avgRelDiff maxRelDiff" << endl;
	for (int t = 0; t < MAX; t++) {
		thres = theta[t]/rangeD;
		numMat = thetaCountMat[t];
		numIce = thetaCountIce[t];
		hitRate = (numMat != 0) ? numIce/(float)numMat : 0;
		avgDiff = (numIce != 0) ? totDiffTheta[t]/(numIce * rangeD) : 0;
		maxDiff = maxDiffTheta[t] / rangeD;
		avgRelDiff = (numIce != 0) ? totRelDiffTheta[t]/(float)numIce : 0;
		maxRelDiff = maxRelDiffTheta[t];
		printf("%7.4f %7d %7d %7.4f %7.3f %7.3f %11.3f %11.3f\n", thres, numMat, numIce, hitRate, avgDiff, maxDiff, avgRelDiff, maxRelDiff);
	}

	cout << endl << "### Iceberg Coverage & Accuracy for Top K% values ###" << endl;
	cout << "K       topFrac #Matrix #Ice    HitRate avgDiff maxDiff avgRelDiff maxRelDiff" << endl;
	for (int k = 0; k < MAX; k++) {
		
		thres = topVal[k]/rangeD;
		numMat = kCountMat[k];
		numIce = kCountIce[k];
		hitRate = (numMat != 0) ? numIce/(float)numMat : 0;
		avgDiff = (numIce != 0) ? totDiffK[k]/(numIce * rangeD) : 0;
		maxDiff = maxDiffK[k] / rangeD;
		avgRelDiff = (numIce != 0) ? totRelDiffK[k]/(float)numIce : 0;
		maxRelDiff = maxRelDiffK[k];
		printf("%7.4f %7.4f %7d %7d %7.4f %7.3f %7.3f %11.3f %11.3f\n", topPct[k], thres, numMat, numIce, hitRate, avgDiff, maxDiff, avgRelDiff, maxRelDiff);
	}

}
///////////////////////////////////////////////////////////////////
