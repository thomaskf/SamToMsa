#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <bitset>
#include <cstdlib>
#include <cmath>
// #include <ctgmath>
#include "reads.h"
#include "partInfo.h"
#include "mylib.h"

#define VERSION "1.11.1" // version of the program

#define MAX_HAP_NUM 3
#define SEQ_ERR_RATE 0.01
// #define PCR_ERR_RATE 0.25
// #define MIN_COVER 0.001

#define MAX_COVER 200 // if the coverage is too large, will normalize it
#define MIN_FREQ 2 // after normalization, if the frequency is smaller than MIN_FREQ, then set it to zero

#define WIN_SIZE_ERR_CORR 50 // window size for error correction = 50
#define WIN_SIZE 50 // window size = 50
#define MIN_COVER_FOR_WIN 20 // minimium coverage for each window
#define MIN_COVER_FOR_WIN_CHECKING 100 // minimium coverage for getting the best window distance
#define MIN_FREQ_FOR_SUBSEQ 2
#define MIN_FREQ_RATIO_FOR_SUBSEQ 0.05 // 0.005
#define ERR_SUBSEQ_RATE 0.01
#define DELTA_SUBSEQ_LOGL 3.0
#define MIN_SUBSEQ_LOGL_THRES -2000.0
#define COVER_THRES_NEAR_END 0.001 // 0.1 % of the average total coverage

#define MANY_ERR 3 // too many errors if the read aligned with any constructed haplotype with more than 2 mismatches

// #define OUT_FIRST_PRE_SEQ

// #define NO_REMOVE_GAP
// #define NO_REMOVE_LOW_COVERAGE

#define DO_BASIC_ERROR_CORR
#define DO_KMER_ERROR_CORR
// #define FILTER_MANY_ERR_RDS
// #define USE_TWO_ITERATIONS
// #define USE_OLD_METHOD1

#define SHOW_DEBUG_MESG
// #define SHOW_DEBUG2_MESG
// #define USE_WEIGHT_APPROACH

using namespace std;

// the following classes are designed for computation of the best combination
class Node {
public:
	Node* pre;
	int pattern;
	int combin;
	// constructor
	Node(Node* pre, int pattern, int combin);
};

class CombinSet {
public:
	Node* lastNode;
	int numItems;
	int numHaps;
	long double cumLogL;
	int numPattern;
	int numCombin;

	char* availAssign;
	// availAssign[i*numCombin + j] = 1 if pattern i can be assigned to combination j
	
	// constructor
	CombinSet(int pattern, int combin, long double logL, int numPattern, int numCombin, vector<Node*>& nodeCollection);
	CombinSet(CombinSet* aCombinSet, int pattern, int combin, long double logL, vector<Node*>& nodeCollection);
	
	// destructor
	~CombinSet();
	
	// show the content
	void show();
	
	// show the content of availAssign
	void showAvailAssign();

	// collect items and patterns
	void collect(vector<int>& items, vector<int>& patterns);

private:
	// identify all unavailable next assignments
	void identUnavailAssign(int pattern, int combin);
};

class CompareCombinSet {
public:
	bool operator() (const CombinSet* lhs, const CombinSet* rhs) const {
		return lhs->cumLogL < rhs->cumLogL;
	}
};

class BestCombin {
public:
	int numCombin;
	int numPattern;
	int numItem;
	int totSize;
	long double* logL;
	// logL[i*numCombin + j] = the log likelihood value of the pattern i for the combination j
	// i starts from 0, but j starts from 1
	
	// constructor
	BestCombin();

	// destructor
	~BestCombin();
	
	// update the dimension
	void resize(int numC, int numP, int numI);

	// compute the best combination with the maximum sum of log-likelihood
	bool compute(vector<int>& items, vector<int>& patterns, long double& bestLogL, long double& best2logL);
};

class ConstructHaplo {
public:

	char nucl_map[256];

	// partition information
	PartInfo partInfo;
	
	// the reads
	Reads rds;
	
	// number of columns
	int numCols;
	
	// number of reads
	int numRds;

	// number of haplotypes
	int numHap;
	
	// number of all possible combination of haplotypes
	int numCombin;
	
	// length of the estimated haplotypes
	int hapLen;
	
	// the corresponding ratios of the haplotype. The sum of the values should be 1.0
	long double* haploRatios;
	
	// haplo ratios for different combinitions
	long double* hapRateCombins;
	long double* logF;
	long double* logFF;
	
	// P(ci | Hj)
	// likelihood (probability) of the char ci, given the set of haplotypes Hi
	long double* baseCombinL;
	
	// number of hits for each read to each haplotype
	vector<vector<int> > hapReadHits;
	
	// mapping between mapPos and mapOrigPos
	// vector<vector<int> > mapping;
	
	// haplotype sequences (ignored the deletions)
	vector<string> haplos;
	
	// full haplotype sequences
	vector<string> fullHaplos;
	
	// set of reads belonging to each haplotype
	vector<vector<int> > hapReads;
	
	// constructor
	ConstructHaplo();

	// destructor
	~ConstructHaplo();

	// load the read file
	void loadReadFile(char* readFile);
	
	// initialise the values of haplotype frequency ratios
	// and allocate the memory according to the number of frequencies inside the file
	void loadFreqFile(char* freqFile);

	// calculate the haplo ratio for each combinition
	void calHapRateCombins();
	
	// calculate P(ci | Hj) for all bases of all columns and combinations
	void calBaseCombinL();
	
	// get the reads beloging to each haplotype
	void getHapReads();

	// get the hits of the read to each estimated haplotype
	void getHapReadHits();

	// construct the haplotypes based on ColCombinL
	void constructHaplos();

	// get the full haplotypes based on hapReadHits
	void getFullHaplos();

	// list out the haplotypes
	void listHaplos();

	// list out the full haplotypes
	void listFullHaplos();

	// show the sequences of the reads for a particular haplotype
	void showReads(int haploID);
	
	// update the reads according to the estimated haplotypes
	void updateReads();

	// remove the gaps inside the haplotype (first version)
	void rmGapInsideHaplo();

	// remove the gaps inside the full haplotype
	void rmGapInsideFullHaplo();

	// remove the gaps and the leading/ending N's inside the full haplotype
	// and remove the region with the too-low-coverage near the ends of the full haplotype
	void rmGapAndLowCoverInsideFullHaplo();

	// correct the reads
	// and normalize the frequencies
	void readCorrection();

	// identify and delete those reads have too many errors
	void deleteErrReads();

	// combine the overlap reads of the same pair into a long read
	// and update the frequencies
	void combineRdAndUpdateFreq();

	//---------------------------------------------------------
	// The methods considering windows
	//---------------------------------------------------------
	vector<vector<string> > hapSubSeqs;
	vector<vector<int> > hapSubSeqPos;
	vector<double> hapSubSeqWeight;

	// haplo ratios for different combinations for subsequences
	long double* sRateCombins;
	long double* slogF;
	long double* slogFF;

	// calculate the ratio for each combinations of subsequences
	void calRateForSubSeqCombins();
	
	// get the set of subsequences belonging to each haplotype
	void getSubSeqForHaplos();

	// update the haplotype sequences according to the set of subsequences
	void showHapSubSeqs();

	// update the haplotype sequences according to the set of subsequences
	void updateHaploBySubSeqs();	

	// get the underlining frequencies of the haplotypes (using windows)
	void getHaploFreq(int winSize, vector<double>& hapFreqs, int& winNum);

	//---------------------------------------------------------
	// The methods considering window-pairs
	//---------------------------------------------------------

	// get the set of pairs of subsequences belonging to each haplotype
	void getSubSeqPairForHaplos(int winDist, int winSize);

	// get the underlining frequencies of the haplotypes (using window-pairs)
	void getHaploFreq(int winDist, int winSize, vector<double>& hapFreqs, int& winNum);

private:

	// distance between haplo and the read
	int distHapRead(string& hap, string& read, int mpos);

	// return the max orig column index
	// int maxOrigColIndex(int hapID);
	
	// delete the memory of the array
	void deleteMem();
	
	// get the set of the read sequences for the specific region
	void getReadSeq(int frPos, int toPos, map<string,int>& seqs);

	// supporting function
	// check the best hits for the read
	// for single read, set rd2 = NULL
	void getBestHits(vector<int>& bestH, Read* rd, Read* rd2);
};
