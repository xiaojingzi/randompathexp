//----------------------------------------------------------------------
//		File:			ann_sample.cpp
//		Programmer:		Sunil Arya and David Mount
//		Last modified:	03/04/98 (Release 0.1)
//		Description:	Sample program for ANN
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------

#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <ANN/ANN.h>					// ANN declarations
#include <map>
#include <string>
#include <cstring>
#include <cstdio>
#include <sys/time.h>
using namespace std;					// make std:: accessible

//----------------------------------------------------------------------
// ann_sample
//
// This is a simple sample program for the ANN library.	 After compiling,
// it can be run as follows.
// 
// ann_sample [-d dim] [-max mpts] [-nn k] [-e eps] [-df data] [-qf query]
//
// where
//		dim				is the dimension of the space (default = 2)
//		mpts			maximum number of data points (default = 1000)
//		k				number of nearest neighbors per query (default 1)
//		eps				is the error bound (default = 0.0)
//		data			file containing data points
//		query			file containing query points
//
// Results are sent to the standard output.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Parameters that are set in getArgs()
//----------------------------------------------------------------------
void getArgs(int argc, char **argv);			// get command-line arguments

int				k				= 1;			// number of nearest neighbors
int				dim				= 2;			// dimension
double			eps				= 0;			// error bound
int				maxPts			= 1000;			// maximum number of data points

istream*		dataIn			= NULL;			// input for data points
istream*		queryIn			= NULL;			// input for query points
ostream*		resultOut		= NULL;			// output for result



std::map<int, std::string> names;

std::string nameFile = ""; 

double timestamp (){
  struct timeval tv;
  gettimeofday (&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

bool readName(){
	FILE * f = fopen(nameFile.c_str() , "r");
	if(f == NULL ) {
		return false;
	}else {
		char line[1024];
		int counter = 0;
		while(fgets(line,1023,f) != NULL){
			counter ++ ;
			std::string name = "";
			std::string val = "";
			bool nameState = true;
			int len = strlen(line);
			for(int i = 0 ; i < len ; i ++ ){
				if(line[i] == '\t') {
					nameState = false;
				}else if(nameState ) {
					name += line[i];
				}else{
					val += line[i];
				}
			}
			int n_val = atoi(val.c_str());
			// printf("name: %s , val : %d \n",name.c_str(),n_val);
			names.insert(std::make_pair(n_val,name));
			if(n_val == 154959){
				printf("test: %d %s \n",n_val,names[n_val].c_str());				
			}
			// printf("test: %d %s \n",n_val,names[n_val].c_str());
			// printf("line :%s\n", line);
		}
		printf("read finished , line size : %d\n", counter);
	}
}


bool readPt(istream &in, ANNpoint p)			// read point (false on EOF)
{
	for (int i = 0; i < dim; i++) {
		if(!(in >> p[i])) return false;
	}
	return true;
}

void printPt(ostream &out, ANNpoint p)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}



int main(int argc, char **argv)
{
	double t00 = timestamp()	;

	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	getArgs(argc, argv);						// read command-line arguments

	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = 0;									// read data points

	// readName();
	double t01 = timestamp()	;
	printf("Time of allocating  memory = %lf s\n", t01-t00);

	double t1 = timestamp()	;
	// cout << "Data Points:\n";
	while (nPts < maxPts && readPt(*dataIn, dataPts[nPts])) {
		//printPt(cout, dataPts[nPts]);
		nPts++;
		//if (nPts % 10000 == 0 ){
		//	cout << nPts << endl;
		//}
	}
	printf("data points = %d\n", nPts);

	double t2 = timestamp();
	printf("Time of loading vectors = %lf s\n", t2-t1);

	kdTree = new ANNkd_tree(					// build search structure
					dataPts,					// the data points
					nPts,						// number of points
					dim);						// dimension of space

	double t3 = timestamp();
	printf("Time of building kd tree = %lf s\n", t3-t2);

	while (readPt(*queryIn, queryPt)) {			// read query points
		// cout << "Query point: ";				// echo query point
		// printPt(cout, queryPt);

		kdTree->annkSearch(						// search
				queryPt,						// query point
				k,								// number of near neighbors
				nnIdx,							// nearest neighbors (returned)
				dists,							// distance (returned)
				eps);							// error bound

		// cout << "\tNN:\tIndex\tDistance\n";
		for (int i = 0; i < k; i++) {			// print summary
			dists[i] = sqrt(dists[i]);			// unsquare distance
			// cout << "\t" << i << "\t" << nnIdx[i] << ":" << names[nnIdx[i]] << "\t" << dists[i] << "\n";
			*resultOut << nnIdx[i] << " ";
		}
		*resultOut << "\n";
	}
	double t4 = timestamp();
	printf("Time of searching top-K = %lf s\n", t4-t3);

    delete [] nnIdx;							// clean things up
    delete [] dists;
    delete kdTree;
	annClose();									// done with ANN

	return EXIT_SUCCESS;
}

//----------------------------------------------------------------------
//	getArgs - get command line arguments
//----------------------------------------------------------------------

void getArgs(int argc, char **argv)
{
	static ifstream dataStream;					// data file stream
	static ifstream queryStream;				// query file stream
	static ofstream resultStrem;      			 // result file stream

	if (argc <= 1) {							// no arguments
		cerr << "Usage:\n\n"
		<< "  ann_sample [-d dim] [-max m] [-nn k] [-e eps] [-df data] [-rf result]"
		   " [-qf query]\n\n"
		<< "  where:\n"
		<< "    dim      dimension of the space (default = 2)\n"
		<< "    m        maximum number of data points (default = 1000)\n"
		<< "    k        number of nearest neighbors per query (default 1)\n"
		<< "    eps      the error bound (default = 0.0)\n"
		<< "    data     name of file containing data points\n"
		<< "    query    name of file containing query points\n\n"
		<< "    result   name of file containing the top k data points of query points\n\n"
		<< " Results are sent to the standard output.\n"
		<< "\n"
		<< " To run this demo use:\n"
		<< "    ann_sample -df data.pts -qf query.pts\n";
		exit(0);
	}
	int i = 1;
	while (i < argc) {							// read arguments
		if (!strcmp(argv[i], "-d")) {			// -d option
			dim = atoi(argv[++i]);				// get dimension to dump
		}
		else if (!strcmp(argv[i], "-max")) {	// -max option
			maxPts = atoi(argv[++i]);			// get max number of points
		}
		else if (!strcmp(argv[i], "-nn")) {		// -nn option
			k = atoi(argv[++i]);				// get number of near neighbors
		}
		else if (!strcmp(argv[i], "-e")) {		// -e option
			sscanf(argv[++i], "%lf", &eps);		// get error bound
		}
		else if (!strcmp(argv[i], "-df")) {		// -df option
			dataStream.open(argv[++i], ios::in);// open data file
			if (!dataStream) {
				cerr << "Cannot open data file\n";
				exit(1);
			}
			dataIn = &dataStream;				// make this the data stream
		}
		else if (!strcmp(argv[i], "-qf")) {		// -qf option
			queryStream.open(argv[++i], ios::in);// open query file
			if (!queryStream) {
				cerr << "Cannot open query file\n";
				exit(1);
			}
			queryIn = &queryStream;			// make this query stream
		}
		else if (!strcmp(argv[i], "-rf")) {		// -rf option
			resultStrem.open(argv[++i], ios::out);// open data file
			if (!resultStrem) {
				cerr << "Cannot open result file\n";
				exit(1);
			}
			resultOut = &resultStrem;				// make this the data stream
		}
		else {									// illegal syntax
			cerr << "Unrecognized option.\n";
			exit(1);
		}
		i++;
	}
	if (dataIn == NULL || queryIn == NULL) {
		cerr << "-df and -qf options must be specified\n";
		exit(1);
	}
}
