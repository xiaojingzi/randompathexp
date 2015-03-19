/* Main program for RoleSim.
Including Iceberg mode and some pre- and post-processing options.
*/
//#define _TIMEVAL
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Graph.h"
#include "GraphRoleUtil.h"
#include "SimMatrix.h"
#include "IcebergSimMap.h" /*ice*/
#include "Util.h" /*post*/
//#include "timer.h"
#include <ctime>

using namespace std;
int N;
const float DEF_BETA = 0.1;
const float DEF_THRESHOLD = 0.01;
const float DEF_IN_WEIGHT = 0.5;
const float DEF_ICE_ESTIMATE_WT = 0.5; //ice
const int GREEDY = 1;
const int EXACT = 2;
const int ICEBERG = 3;
/*

*/

template <class T>
inline std::string toString (const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

void usage() {
	cout << "\nUsage:\n"
		"main -h\n"
/*ice*/	"main -graph g [-directed][-init i][-sim s] [-matching m][-waterline w][-iceDefault e][-beta b][-inWeight w][-threshold t][-maxIter k]\n"
		" [-result rfile] [-out ofile] [-dist] [-verbose] [-compact] [-precision p]\n"
		" [-shell | -components | -giant | -append g2]\n"
/*ice*/		" [-loadvec vecFile [-loadvec2 vecFile2]][-loadice iceFile][-resume]\n" 
		" [-compare | -writerank rankfile | -readrank rankfile][-authors afile][-size1 size1][-block\n"
		">>>Description:\n"
		"	-h					Print the help message.\n"
		"	-g, -graph g		dataset name, REQUIRED\n"
		"	-d, -directed		treat the graph as directed (default: undirected)\n"
		">>>Computation options:\n"
		"	-s, -sim s			similarity measure: rs(default)=RoleSim, sr=SimRank, sr++=SimRank++, psr=PSimRank\n"
		"	-i, -init i			initialization: 1(default)=All 1's, r=degree ratio, b=binary(1 if degrees match, else 0), s=self-similarity\n" 
		"	-m, -matching m		matching method:1-Greedy(default), 2-Exact, 3-Iceberg.\n"
		"	-inw, -inWeight a	for directed graphs, fractional weight (default 0.5) for in-neighbors, i.e., sim = a*in_sim + (1-a)*out_sim.\n"
		"	-b, -beta b			beta(default:0.1).\n"
		"	-t, -threshold t	threshold value(default: 0.01).\n"
		"	-maxIter k			maximum number of iterations (default = 25).\n"
/*ice*/	"	-w, -waterline w	for Iceberg RoleSim, only store node-pairs whose estimated sim value > w, 0 < w <= 1 (default 0.75)\n" 
/*ice*/	"	-iceDefault e		for Iceberg RoleSim, estimated similarity of non-stored node-pair (u,v) = iceD*degree(u)/degree(v) + beta (default 0.5)\n" 
		"	-resume				resume simulation from a saved state.  ALL the regular sim options are needed, PLUS -loadvec and -resume\n"
		"Output options:\n"
		"	-result    rfile	time and rounds result file name(default:result.txt).\n"
		"	-o, -out   outfile	output matrix file name(default:auto-generated name describing input options).\n"
		"	-dist				invert the output to be a distance matrix: distance = 1 - similarity.\n"
		"	-v, -verbose		send lots of progress info (including every matrix interation) to stdout\n"
		"	-c, -compact		Output matrices as compact vectorized triangle matrices\n"
		"	-p, -precision p	bits of precision to write/read for compact result values (default 3)\n"
		"Pre- and Post-processing options:\n"
		"	-shell				decompose the graph into shell components (pre)\n"
		"	-components			decompose the graph into connected components (pre)\n"
		"	-giant				find the giant component subgraph\n"
		"	-append				append a second graph to the first graph (disjoint components)\n"
		"	-loadvec vecFile	reload a similarity matrix (in compact format)\n"
		"	-loadvec2 vecFile2	reload a 2nd similarity matrix (in compact format)\n"
/*ice*/	"	-loadice iceFile	reload an iceberg save file.\n" 
/*post*/"	-writerank rankFile	write a table of RoleSim/SimRank percentile rank values\n"
/*post*/"	-readrank  rankFile	read a table of RoleSim/SimRank percentile rank values\n" 
/*post*/"	-authors   autFile	read a table of author ids, G-index rank, H-index rank\n" 
/*post*/"	-blocks    blckFile	(for block graphs) read a file containing the block sizes\n" 
		"	-size1				size of the first graph (if -append was used)\n"
/*post*/"	-compare			Compares two result files, if two have been loaded (-loadvec, -loadvec2 or -loadice).\n" 
		"EXAMPLES:\n"
 		"To run RoleSim on a small graph, with all-1's initialization and full matrix output:\n"
 		"./main -g myGraphFile -s rs > logFile\n"
		"\n"
 		"To run SimRank on a small graph, with self-similarity (identity matrix) initialization, full matrix output sent to a named file:\n"
 		"./main -g myGraphFile -s sr -i s -out mySimMatrix> logFile\n"
		"\n"
 		"To run RoleSim on a large graph, with degree-ratio initialization and compact output:\n"
 		"./main -g myGraphFile -s rs -i r -c -out mySimVec >logFile\n"
		"\n"
	<< endl;
}

string timeToString(double time) {
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

void displayElapsedTime(time_t& previous, ostream& out) {
	time_t current = time(NULL);
	double diff = difftime(current, previous);
	out << "elapsed time = " << timeToString(diff) << endl;
	previous = current;
}

void keepResult( char* resultFileName, char* experimentName, double constructionTime, char* simMeasure, char init,
				float threshold, int matchingMethod, float beta, int round, float inWeight, bool isDirected,
/*ice*/			float waterline, float iceDefaultWt) {

	ofstream out(resultFileName, ios_base::out|ios_base::app);
	time_t now = time(NULL);
	out << ctime(&now);
	out << "*********************************************************" << endl;
	out << experimentName << ", simMeas:" << simMeasure << ", init:" << init << ", matchMeth:" << matchingMethod;
	out << ", Thresh:" << threshold << ", Beta:" << beta << endl;
	out << ", Directed:" << isDirected << ", inWt:" << inWeight;
	out << ", waterline:" << waterline << ", iceDef:" << iceDefaultWt;
	out << ", Time = " << timeToString(constructionTime) << "\tRounds = " << round << endl;
	out.close();
	cout << "Time = " << timeToString(constructionTime);
}

string autoName(char* GraphFileName, bool isDirected, float inWeight, char* simMeasure, char init,
		int matchingMethod, float threshold, float beta, float waterline, float iceDefaultWt, bool distance, bool compact, int bitsPrecision )
{
	ostringstream ss;
	string GraphBase(GraphFileName);
	GraphBase = GraphBase.substr(0, GraphBase.find_last_of("."));
	ss << GraphBase;
	if (isDirected) {
		ss << "D";
		if (inWeight != DEF_IN_WEIGHT) {
			ss << inWeight * 10;
		}
	}
	ss << "-" << simMeasure << "-I" << init << "-M" << matchingMethod;
	if (threshold != DEF_THRESHOLD) {
		ss << "-T" << (threshold * 1000);
	}
	if (beta != DEF_BETA) {
		ss << "-B" << (beta * 1000);
	}
	if (matchingMethod == ICEBERG) { /*ice*/	
		ss << "-W" << waterline;
		ss << "-E" << iceDefaultWt;
	}

	if (compact && bitsPrecision != 3) {
		ss << "-p" << bitsPrecision;
	}
	else if (compact) {
		ss << "-c";
	}
	if (distance) {
		ss << "-dst";
	}
	else {
		ss << "-sim";
	}
	ss << ".txt";
	return ss.str();
}


int main(int argc, char* argv[]) {
	char* GraphFileName = "";
	string MatrixOutputFileName = "SimMatrix.txt";
	char* ResultFileName = "result.txt";
	char* MatrixInputFileName = "0";
	char* MatrixInput2FileName = "0";
/*ice*/	string PercentrankFileName = "0";
	char* Graph2FileName = "0";
/*ice*/	char* AuthorFileName = "0";
/*ice*/	char* BlockSizeFileName = "0";
/*ice*/	char* IcebergFileName = "0";
	float threshold = DEF_THRESHOLD;
	float beta = DEF_BETA;
	int matchingMethod = 1;
	char init = '1';
	char* simMeasure = "rs";
	bool distance = false;
	bool verbose = false;
	bool compact = false;
	int	maxIter = 25;
	bool findShells = false;
	bool findComp = false;
	bool findGiant = false;
	bool isDirected = false;
	bool loadVec = false;
	bool loadVec2 = false;
/*ice*/	bool loadIce = false;
	int bitsPrecision = 3;
/*ice*/	bool testload = false;
/*ice*/	bool writerank = false;
/*ice*/	bool readrank = false;
	bool append = false;
/*ice*/	bool readAuthors = false;
/*ice*/	bool readBlockSizes = false;
	float inWeight = DEF_IN_WEIGHT;
	int size1 = 0;
/*ice*/	float waterline = 0.75;
/*ice*/	float iceDefaultWt = DEF_ICE_ESTIMATE_WT;
	bool userOutputFile = false;
/*ice*/	bool compare = false;
	bool shortcutMatching = true;
	bool resume = false;
/*ice*/	bool meltIceberg = false;
#ifdef _TIMEVAL
	struct timeval after_time, before_time;
#endif
	double construction_time = 0.0;

	int i = 1;
	//command input process
	cout << argc << " arguments" << endl;
	while (i < argc) {
		cout <<  argv[i] << " " << endl;
		if(strcmp("-h", argv[i]) == 0) {
			usage();
			exit(0);
		}

		if(strcmp("-graph", argv[i]) == 0 || strcmp("-g", argv[i]) == 0) {
			i++;
			GraphFileName = argv[i++];
			fprintf(stderr,"Grapgh :::: %s \n",GraphFileName);
			fflush(NULL);
		}
		else if(strcmp("-directed", argv[i]) == 0 || strcmp("-d", argv[i]) == 0) {
			i++;
			isDirected = true;
		}
		else if (strcmp("-result", argv[i]) == 0) {
			i++;
			ResultFileName = argv[i++];
		}
		else if (strcmp("-out", argv[i]) == 0 || strcmp("-o", argv[i]) == 0) {
			i++;
			userOutputFile = true;
			MatrixOutputFileName = argv[i++];
		}
		else if (strcmp("-dist", argv[i]) == 0) {
			i++;
			distance = true;
		}
		else if(strcmp("-sim", argv[i]) == 0 || strcmp("-s", argv[i]) == 0) {
			i++;
			simMeasure = argv[i++];
		}
		else if(strcmp("-init", argv[i]) == 0 || strcmp("-i", argv[i]) == 0) {
			i++;
			init = argv[i++][0];
		}
		else if(strcmp("-threshold", argv[i]) == 0 || strcmp("-t", argv[i]) == 0) {
			i++;
			threshold = atof(argv[i++]);
		}
		else if(strcmp("-matching", argv[i]) == 0 || strcmp("-m", argv[i]) == 0) {
			i++;
			matchingMethod = atoi(argv[i++]);
		}
/*ice*/
		else if(strcmp("-waterline", argv[i]) == 0 || strcmp("-w", argv[i]) == 0) {
			i++;
			waterline = atof(argv[i++]);
		}
		//
		else if(strcmp("-inWeight", argv[i]) == 0 || strcmp("-inw", argv[i]) == 0) {
			i++;
			inWeight = atof(argv[i++]);
		}
/*ice*/
		else if(strcmp("-iceDefault", argv[i]) == 0) {
			i++;
			iceDefaultWt = atof(argv[i++]);
		}
		else if(strcmp("-melt", argv[i]) == 0) {
			i++;
			meltIceberg = true;
		}
		//
		else if(strcmp("-beta", argv[i]) == 0 || strcmp("-b", argv[i]) == 0) {
			i++;
			beta = atof(argv[i++]);
		}
		else if(strcmp("-maxIter", argv[i]) == 0) {
			i++;
			maxIter = atoi(argv[i++]);
		}
		else if(strcmp("-verbose", argv[i]) == 0 || strcmp("-v", argv[i]) == 0) {
			i++;
			verbose = true;
		}
		else if(strcmp("-compact", argv[i]) == 0 || strcmp("-c", argv[i]) == 0) {
			i++;
			compact = true;
		}
		else if(strcmp("-precision", argv[i]) == 0 || strcmp("-p", argv[i]) == 0) {
			i++;
			bitsPrecision = atoi(argv[i++]);
		}
		else if(strcmp("-shell", argv[i]) == 0) {
			i++;
			findShells = true;
		}
		else if(strcmp("-components", argv[i]) == 0) {
			i++;
			findComp = true;
		}
		else if(strcmp("-giant", argv[i]) == 0) {
			i++;
			findGiant = true;
		}
		else if(strcmp("-append", argv[i]) == 0) {
			i++;
			append = true;
			Graph2FileName = argv[i++];
		}
		else if(strcmp("-loadvec", argv[i]) == 0) {
			i++;
			MatrixInputFileName = argv[i++];
			loadVec = true;
			compact = true;
		}
		else if(strcmp("-loadvec2", argv[i]) == 0) {
			i++;
			MatrixInput2FileName = argv[i++];
			loadVec2 = true;
			compact = true;
		}
/*ice*/
		else if(strcmp("-loadice", argv[i]) == 0) {
			i++;
			loadIce = true;
			IcebergFileName = argv[i++];
		}
		//
/*post*/
		else if(strcmp("-writerank", argv[i]) == 0) {
			i++;
			writerank = true;
			PercentrankFileName = argv[i++];
		}
		else if(strcmp("-readrank", argv[i]) == 0) {
			i++;
			readrank = true;
			PercentrankFileName = argv[i++];
		}
		else if(strcmp("-authors", argv[i]) == 0) {
			i++;
			readAuthors = true;
			AuthorFileName = argv[i++];
		}
		else if(strcmp("-size1", argv[i]) == 0) {
			i++;
			size1 = atoi(argv[i++]);
		}
		else if(strcmp("-blocks", argv[i]) == 0) {
			i++;
			readBlockSizes = true;
			BlockSizeFileName = argv[i++];
		}
		else if(strcmp("-fastice", argv[i]) == 0) {
			i++;
			shortcutMatching = true;
		}
		else if(strcmp("-compare", argv[i]) == 0) {
			i++;
			compare = true;
		}
		//
		else if(strcmp("-resume", argv[i]) == 0) {
			i++;
			resume = true;
		}
		else {
			cerr << "ERROR: Unrecognized argument: " << argv[i] << endl << "TERMINATING" << endl;
			exit(1);
		}
	}
	fprintf(stderr,"data inited\n"); fflush(NULL);
	if (!userOutputFile)  {
		/*ice*/ MatrixOutputFileName = autoName(GraphFileName, isDirected, inWeight, simMeasure, init,matchingMethod, threshold, beta, waterline, iceDefaultWt, distance, compact, bitsPrecision );
	}
	cout << endl;
	cout << "Program Configuration:" << endl 
		<< "\t dataset:" << GraphFileName << endl
		<< "\t directed:" << isDirected << endl
		<< "\t init:" << init << endl
		<< "\t sim:" << simMeasure << endl
		<< "\t threshold:" << threshold << endl
		<< "\t maxIter:" << maxIter << endl
		<< "\t matchingMethod:" << matchingMethod << endl
		<< "\t inWeight:" << inWeight << endl
		<< "\t beta:" << beta << endl
/*ice*/	<< "\t waterline:" << waterline << endl
/*ice*/	<< "\t iceDefault:" << iceDefaultWt << endl
/*ice*/	<< "\t meltIceberg:" << meltIceberg << endl
		<< "\t result:" << ResultFileName << endl 
		<< "\t output matrix:" << MatrixOutputFileName << endl 
		<< "\t distance:" << distance << endl
		<< "\t verbose:" << verbose << endl
		<< "\t compact:" << compact << endl
		<< "\t shell:" << findShells << endl
		<< "\t connect:" << findComp << endl
		<< "\t giant:" << findGiant << endl
		<< "\t append:" << Graph2FileName << endl
		<< "\t loadVec:" << MatrixInputFileName << endl
		<< "\t loadVec2:" << MatrixInput2FileName << endl
/*ice*/	<< "\t loadIce:" << IcebergFileName << endl
		<< "\t bitsPrecision:" << bitsPrecision << endl
/*ice*/	<< "\t writerank:" << PercentrankFileName << endl
/*ice*/	<< "\t readrank:" << PercentrankFileName << endl
/*ice*/	<< "\t -authors:" << AuthorFileName << endl
		<< "\t size1:" << size1 << endl
/*ice*/	<< "\t readBlockSizes:" << BlockSizeFileName << endl
/*ice*/	<< "\t compare:" << compare << endl
		<< "\t resume:" << resume << endl
		;

	//******************************************
	// * Step 1: Read in graph; perform any pre- or post-processing
	//******************************************

	ifstream infile1(GraphFileName);
	if (!infile1) {
		cerr << "Error: Cannot open " << GraphFileName << endl;
		return -1;
	} else {
		cerr << "File Opened : " << GraphFileName << endl;
	}

	Graph g(infile1);
	fprintf(stderr,"g1 loaded .. \n");
	fflush(NULL);
	N = g.num_vertices();
	g.setDirected(isDirected);
	if (isDirected) {
		g.sortDirectedEdges(); // victor
	} else {
		g.sortAndSetUndirectedEdges();
	}
	if (verbose && N < 100) {
		g.writeGraph(cout);
	}
	cout << "\t #vertex size:" << N << "\t#edges size:" << g.num_edges() << endl;
	GraphRoleUtil grUtil(g);
/*post*/ Util util(N, verbose);


	//////// Post Processing Options ////////
/*post*/
	if (loadVec && loadVec2 && compare) {
		vector<int> simVec1 =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);
		vector<int> simVec2 =  util.ReloadSimilarityTriangleVector(MatrixInput2FileName, N, bitsPrecision, distance);
		cout << "comparing two similarity matrices" << endl;
		util.compareTwoSimVecs(simVec1, simVec2, bitsPrecision);

		cout << "comparing two similarity rankings" << endl;
		map<int,float> rankMap1, rankMap2;
		util.rankMatrixValues(simVec1, bitsPrecision, rankMap1);
		util.rankMatrixValues(simVec2, bitsPrecision, rankMap2);
		util.compareTwoRankings(simVec1, simVec2, rankMap1, rankMap2, bitsPrecision);

		exit(0);
	}
	//
/*ice*/
	else if (loadVec && loadIce & compare) {
		cout << "trying to reload matrix file " << MatrixInputFileName << endl;
		vector<int> simVec =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);

		IcebergSimMap iceberg = IcebergSimMap();
		cout << "trying to reload iceberg file " << IcebergFileName << endl;
		ifstream iceIn(IcebergFileName);
		int loadStatus = iceberg.load(iceIn, bitsPrecision, g, compact);
		iceIn.close();

		map<int,float> icePctRank, matIcePctRank, matFullPctRank;
		util.rankIcebergValues( iceberg, simVec, bitsPrecision, icePctRank, matIcePctRank );
		util.rankMatrixValues( simVec, bitsPrecision, matFullPctRank );

		util.checkIcebergCoverage(iceberg, simVec, bitsPrecision, icePctRank, matFullPctRank);
		util.compareIcebergToSimMatrix( iceberg, simVec, bitsPrecision, icePctRank, matIcePctRank );
		exit(0);
	}
	//
/*post*/
	if (loadVec && writerank && matchingMethod != ICEBERG) {
		vector<int> simVec =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);
		cout << "read vectorized matrix, size = " << simVec.size() << endl;
		vector<int> pctileVec = util.computePercentileRankBins(simVec, 100);
		cout << "writing prctileVec to " << PercentrankFileName << endl;
		ofstream out_p(PercentrankFileName.c_str());
		util.writePercentileRankBins(out_p, pctileVec);
		exit(0);
	}
	else if (loadVec && writerank && matchingMethod == ICEBERG) {
		IcebergSimMap iceberg = IcebergSimMap();
		ifstream iceIn(MatrixInputFileName);
		iceberg.load(iceIn, bitsPrecision, g, compact);
		vector<int> pctileVec = util.computePercentileRankBins(iceberg, 100, bitsPrecision);
		ofstream out_p;
		cout << "writing prctileVec to " << PercentrankFileName << endl;
		out_p.open(PercentrankFileName.c_str());
		util.writePercentileRankBins(out_p, pctileVec);
		exit(0);
	}
	else if (loadVec && readrank && readBlockSizes) {
		vector<int> simVec =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);
		cout << "read vectorized matrix, size = " << simVec.size() << endl;
		vector< vector<int> > shells = grUtil.makeBlockShells(BlockSizeFileName, verbose);
		vector<int> pctileVec = util.readPercentileRankBins(PercentrankFileName.c_str());
		util.doShellAnalysis(simVec, shells, pctileVec, bitsPrecision, verbose);
		exit(0);
	}
	else if (loadVec && readrank && readAuthors) {
		vector<int> simVec =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);
		cout << "read vectorized matrix, size = " << simVec.size() << endl;
		vector<int> pctileVec = util.readPercentileRankBins(PercentrankFileName.c_str());
		vector<AuthorItem> authorVec = util.readAuthorRanks(AuthorFileName);
		// graph has two components.  size1 = size of the first components ( node ids 0:(size1 - 1) )
		util.doCoauthorAnalysis(simVec, authorVec, pctileVec, bitsPrecision, size1, verbose);
		exit(0);
	}
	else if (loadVec && readrank && findShells) {
		vector<int> simVec =  util.ReloadSimilarityTriangleVector(MatrixInputFileName, N, bitsPrecision, distance);
		cout << "read vectorized matrix, size = " << simVec.size() << endl;
		vector< vector<int> > shells = grUtil.findKShells(g, false);
		vector<int> pctileVec = util.readPercentileRankBins(PercentrankFileName.c_str());
		util.doShellAnalysis(simVec, shells, pctileVec, bitsPrecision, verbose);
		exit(0);
	}
//
	//////// Perform any sub-graph pre-processing ////////
	if (findGiant) {
		Graph giantSubgraph = grUtil.findGiantComponent(g, verbose);
		string gName(GraphFileName);
		string giantName = gName.substr(0,gName.find_last_of(".")) + "_giant.gra";
		ofstream gout(giantName.c_str());
		giantSubgraph.writeGraph(gout);
		exit(0);
	}
	if (append) {
		cout << "Reading graph " << Graph2FileName << endl;
		ifstream infile2(Graph2FileName);
		if (!infile2) {
			cerr << "Error: Cannot open " << Graph2FileName << endl;
			return -1;
		}
		Graph g2(infile2);
		g.appendGraph(g2, isDirected);
		if (isDirected) {
			g2.sortDirectedEdges();
		} else {
			g2.sortAndSetUndirectedEdges();
		}
		string g1Name(GraphFileName);
		g1Name = g1Name.substr(0,g1Name.find_last_of("."));
		string g2Name(Graph2FileName);
		g2Name = g1Name + "_" + g2Name.substr(g2Name.find_last_of("/") + 1);
		cout << "Writing to file " << g2Name << endl;
		ofstream gout(g2Name.c_str());
		g.writeGraph(gout);
		exit(0);
	}
	if (findComp) {
		vector< vector<int> > compVec = grUtil.findComponents(g, verbose);
		string gName(GraphFileName);
		string compPrefix = gName.substr(0,gName.find_last_of(".")) + "_c";
		grUtil.writeComponents(compVec, compPrefix);
		exit(0);
	}
	if (findShells) {
		vector< vector<int> > shellVec = grUtil.findKShells(g, verbose);
		string gName(GraphFileName);
		string shellPrefix = gName.substr(0,gName.find_last_of(".")) + "_s";
		grUtil.writeComponents(shellVec, shellPrefix);
		exit(0);
	}

	//******************************************
	// * Step 2: Create Similarity Matrix; initialize all cells = 0
	//******************************************
	//IcebergSimMap simHashEven, simHashOdd;
	//SimMatrix matrixEven, matrixOdd;
/*ice*/	IcebergSimMap iceberg[2];
	SimMatrix matrix[2];
/* begin ice*/
	if (matchingMethod == ICEBERG)
	{
		iceberg[0] = IcebergSimMap( waterline, iceDefaultWt, beta, inWeight, g, verbose );
		iceberg[1] = IcebergSimMap( waterline, iceDefaultWt, beta, inWeight, g, verbose );
	}
	else
	{ /* end ice*/
		//matrixEven = SimMatrix(N);
		//matrixOdd = SimMatrix(N);
		matrix[0] = SimMatrix(N);
		matrix[1] = SimMatrix(N);
	} /*ice*/
	cout << "Created blank similarity matrices" << endl;

	//******************************************
	// * Step 3: Initialize Similarity Matrix
	//******************************************

	int simCode = matrix[0].simCode(simMeasure);
	cout << "simCode = " << simCode << endl;

#ifdef _TIMEVAL
		gettimeofday(&before_time, NULL);
#else
	clock_t start_time(clock());	
#endif
	time_t previous = time(NULL);
/* begin ice*/
	if (matchingMethod == ICEBERG) {
		if (resume) {
			iceberg[0] = IcebergSimMap();
			ifstream iceIn(IcebergFileName);
			int loadStatus = iceberg[0].load(iceIn, bitsPrecision, g, compact);
			iceIn.close();
		}
		else {
			iceberg[0].InitializeRatioDepth1(g, shortcutMatching);
		}
		if (verbose && N < 100) {
			cout << " Initialization:" << endl;
			iceberg[0].PrintCompact( cout, bitsPrecision, compact );
		}
	}

	else { /*end ice*/
		if (resume) {
			matrix[0].loadSimVec(MatrixInputFileName, N, bitsPrecision, distance);
		}
		else if(simCode == ROLESIM)	// RoleSim
		{
			if(init == '1' || init == 'u' || init == 'U') {
				matrix[0].InitializeUniform( );
			}
			else if(init == 's' || init == 'S') {
				matrix[0].InitializeDiagonal( );
			}
			else if(init == 'b' || init == 'B') {
				matrix[0].InitializeDepthZero( g );
			}
			else if(init == 'r' || init == 'R') {
				matrix[0].InitializeRatioDepth0( g, beta );
			}
		}
		else //if (simCode == SIMRANK)
		{
			matrix[0].InitializeDiagonal( );
		}
		if (verbose && N < 100) {
			cout << " Initialization:" << endl;
			matrix[0].PrintCompact( cout, bitsPrecision, compact );
		}
	} /*ice*/
	cout << "Init option = " << init << ", ";
	displayElapsedTime(previous, cout);
		
	//******************************************
	// * Step 4: Perform Iterations
	//******************************************
	int turn = 0;
	int curr = 0;
	int prev = 1;
	
//#ifdef _TIMEVAL
//		gettimeofday(&before_time, NULL);
//#else
//	clock_t start_time(clock());	
//#endif
	previous = time(NULL);

	bool divMax = true;
	bool header = false;

/* begin ice*/
	if (simCode == ROLESIM && matchingMethod ==ICEBERG && maxIter > 0 )
	{
		///////////// RoleSim - iceberg //////////////////////////////
		// Before the first iteration, the prev iceberg is empty, so diffRelative() doesn't work
		do {
			prev = curr;
			curr = (++turn) % 2;

			cout << "Round " << turn << endl;
			iceberg[curr].updateRoleSim(iceberg[prev], beta, divMax, meltIceberg);
			displayElapsedTime(previous, cout);

			if (verbose && N < 100) {
				iceberg[curr].PrintCompact( cout, bitsPrecision, compact );
			}
		}
		while( !iceberg[curr].DiffRelative( iceberg[prev], threshold ) && turn < maxIter);

	}
	else if ( simCode == ROLESIM && matchingMethod != ICEBERG )
/*end ice*/
	if ( simCode == ROLESIM )
	{
		cout << "Start RoleSim simMatrix iterations" << endl;
		///////////// RoleSim //////////////////////////////
		do {
			prev = curr;
			curr = (++turn) % 2;

			cout << "Round " << turn << endl;
			if( matchingMethod == GREEDY ) {
				matrix[curr].UpdateRoleSimGreedy( g, matrix[prev], beta, divMax, inWeight);
			}
			else if ( matchingMethod == EXACT ) {
				matrix[curr].UpdateRoleSimExact( g, matrix[prev], beta, divMax );
			}
			displayElapsedTime(previous, cout);
			if (verbose && N < 100) {
				matrix[curr].PrintCompact( cout, bitsPrecision, compact );
			}
		}
		while( !matrix[curr].DiffRelative( matrix[prev], threshold )  && turn < maxIter);
	}
																				
	else if (simCode == SIMRANK || simCode == SIMRANK_PP || simCode == P_SIMRANK)
	{	
		cout << "Start SimRank simMatrix iterations" << endl;
		///////////// SimRank //////////////////////////////
		do {
			prev = curr;
			curr = (++turn) % 2;

			cout << "Round " << turn << endl;
			if (simCode == P_SIMRANK)  {
				matrix[curr].UpdatePSimRank( g, matrix[prev], beta, inWeight );
			}
			else { // simCode = SIMRANK or SIMRANK_PP
				matrix[curr].UpdateSimRank( g, matrix[prev], beta, inWeight );
			}
			displayElapsedTime(previous, cout);
		}
		while( !matrix[curr].DiffRelative( matrix[prev], threshold ) && turn < maxIter);

		if (simCode == SIMRANK_PP) {
			matrix[curr].addSimRankPPEvidence(g, inWeight);
		}
	}
	else if (simCode == BIBLIO_COUPLING)
	{	
		cout << "Start Bibliographic Coupling" << endl;
		///////////// Biblio Coupling  //////////////////////////////
		curr = turn = 1;
		matrix[curr].BibliographicCoupling(g);
		displayElapsedTime(previous, cout);
	}
	else if (simCode == CO_CITATION)
	{	
		cout << "Start Co-Citation" << endl;
		///////////// Co-Citation  //////////////////////////////
		curr = turn = 1;
		matrix[curr].CoCitation(g);
		displayElapsedTime(previous, cout);
	}

	//******************************************
	// * Step 5: Display Results
	//******************************************

#ifdef _TIMEVAL
		gettimeofday(&after_time, NULL);
		construction_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
#else
		clock_t finish_time(clock());
		construction_time = (finish_time - start_time)/CLOCKS_PER_SEC;
#endif
/*begin ice*/
	if (matchingMethod == ICEBERG)
	{
		string outName = MatrixOutputFileName;
		string mapName = outName.substr(0,outName.find_last_of(".")) + "-ice.txt";
		ofstream mapOut(mapName.c_str());
		if ( N > 50000) {
			iceberg[curr].save(mapOut, GraphFileName, bitsPrecision, true);
		}
		else{
			iceberg[curr].save(mapOut, GraphFileName, bitsPrecision, compact);
		}
		mapOut.close();

		// write the regular matrix only if small and verbose
		if (verbose && N < 25000) {
			ofstream out;
			out.open(MatrixOutputFileName.c_str());
			if (distance) {
				iceberg[curr].PrintMatlabCompact( out, bitsPrecision, compact );
			}
			else {
				iceberg[curr].PrintCompact( out, bitsPrecision, compact, false );
			}
			out.close();
		}
	}
	else
	{ /*end ice*/
		ofstream out;
		out.open(MatrixOutputFileName.c_str());
		if (distance) {
			matrix[curr].PrintMatlabCompact( out, bitsPrecision, compact );
		}
		else {
			matrix[curr].PrintCompact( out, bitsPrecision, compact, false );
		}
		out.close();
	} /*ice*/

	cout << " Finished, Rounds " << turn << endl;

	keepResult( ResultFileName, GraphFileName, construction_time, simMeasure, init, threshold, matchingMethod, beta,
		turn, inWeight, isDirected, waterline, iceDefaultWt);
}
