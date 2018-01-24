#ifndef _SR_PROJ_PARTINFO_
#define _SR_PROJ_PARTINFO_

#include "reads.h"
#include "mylib.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <iomanip>

using namespace std;

typedef vector<int> intVector;

typedef pair<int,int> pairInt;

class PartInfo {
public:
	// information from the reads
	int* colStat; // col[i*5+j] : # of char j appears in column i
	int* invarCol; // invariable columns
	int* coverage;
	int numCol;
	char nuclMap[256];

	// constructor
	PartInfo();

	// destructor
	~PartInfo();

	// load from the reads
	// only for nucleotide reads
	void loadFromReads(Reads& reads);
	
	// print read information
	void printReadInfo();
	
	// print read information
	void printReadInfo(int i);

	// get the average coverage of the reads
	double getAvgCoverage();
	
	// normalize the coverage if it is too large
	// after normaliation, if the freq < min_freq, set the freq to zero
	void normalize(int max_thres, int min_freq);

	// reset the statistics
	// option = 1: using the seq
	// option = 2: using the updatedSeq
	void updateFromReads(Reads& reads, int option);

private:
	// initialize the array nuclMap
	void initNuclMap();

	// delete all the memory used in the array
	void deleteMem();
};

#endif
