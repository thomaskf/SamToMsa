#ifndef _SR_PROJ_READS_
#define _SR_PROJ_READS_

#include<string>
#include<vector>
#include<map>
#include<fstream>
#include<iostream>
#include<algorithm>
#include "fileHander.h"
#include "mylib.h"

#define MIN_OVERLAP_COMBINE_READ 5

// for window-pair approach
#define WIN_PAIR_LEN 50
#define MIN_WIN_PAIR_DIST 0
#define MAX_WIN_PAIR_DIST 500
#define DIST_INTERVAL 10
#define POS_INTERVAL 10

using namespace std;

class Insert {
public:
	int col;
	int readPos;
	string x;
	int readid;
	int len;
	Insert(int c, int p, string s, int id);
	void show();
};

class Read {
public:
	int id;
	string seq;
	int mapPos;
	string mapHap; // for debugging
	// int mapLen;
	int seqLen;
	string origSeq;
	string updateSeq;
	string desc;
	// string baseDeleted;
	// int mapOrigPos;
	bool isDiscard; // the read is discarded because it has too many errors
	
	// Read(string oSeq, string s, int mPos, string mHap, string baseDel); // constructor
	Read(string oSeq, string s, int mPos, string mHap, string des); // constructor

	// get the number of deleted bases in the beginning of the read
	// int numDelBegin();
	
	// combine this read with another read
	// return true if successful
	bool combineRead(Read* aRead, int& newMPos, string& newSeq);
};

class Reads {
public:
	int readNum;
	// a list of reads
	vector<Read*> readlist;
	vector<Insert*> inserts;
	
	vector<Read*> pairRead1; // paired-end reads
	vector<Read*> pairRead2;
	vector<Read*> unpairReads; // single-end reads
	
	// avg dist between the pair-ended reads
	double avgReadDist;
	
	// compute the avg distance between the starting positions of the pair-ended reads
	void computeAvgReadPosDist();
	
	// check the best pair-end window with the highest coverage
	// report the average coverage
	double computeBestWindowPair(int& bestDist);
	
	// check the best pair-end window with the highest diversity
	// report the average coverage
	double computeBestDivWindowPair(int numHaps, int minWinCover, int minFreq,
		double minFreqRatio, int& bestDist);

	// load reads from the SAM file
	void readSamFile(char* samFile);

	// return the max column index
	int maxColIndex();
	
	// sort the reads
	void sortReads();

	// update according to the insertions
	void updateForInserts();

	bool updateSeq(string& seq, string& cigarStr, int readId, int mapPos);

	// get the pairs of subsequences for a region (with freq >= minFreq)
	// return the total number of pairs of reads covering the region
	double getSubSeqPairs(int fr_pos1, int to_pos1, int fr_pos2, int to_pos2, int minFreq,
			double minFreqRatio, vector<string>& subSeq1, vector<string>& subSeq2, vector<double>& subSeqFreq);

	// get all the subsequences for a region
	// return the total number of fragments covering the region
	int getAllSubSeqs(int fr_pos, int to_pos, map<string,double>& subSeqMap, vector<Read*>* readArr);

	// get the subsequences for a region (with freq >= minFreq and freq >= tot * minFreqRatio)
	// return the total number of fragments above the freq threshold covering the region
	double getSubSeqs(int fr_pos, int to_pos, int minFreq, double minFreqRatio, vector<string>& subSeq, vector<double>& subSeqFreq);

	// print out the reads
	void printReads(int fr_pos, int to_pos);
	// print out the reads
	// void printOrigReads(int fr_pos, int to_pos);
	// print out the reads indicated inside the vector
	// void printOrigReads(int fr_pos, int to_pos, vector<int>& haploIDs);
	// print out all the original reads
	// void printAllOrigReads();

	// print out all the reads
	void printAllReads();

	// print out the reads indicated inside the vector
	void printReads(int fr_pos, int to_pos, vector<int>& haploIDs);
	
	// print out the read groups
	void printReadGrps(int fr_pos, int to_pos);
	
	// assign the IDs for the reads
	void assignReadIDs();
	
	// combine the overlapped read to a longer read
	void combineReads();

	// output the reads
	void outReads(char* pairRdFile, char* singleRdFile);
	
	// error correction by using relative k-mer approach
	void errorCorrection(int winSize, int minCover, int minFreq, double minFreqRatio);
	
	// get frequencies of the haplotypes
	void getFreqs(int winSize, int minCover, int minFreq, double minFreqRatio, int numHaps, vector<double>& freqs);

private:

	int readNumCoverWinPair(int winStart1, int winEnd1, int winStart2, int winEnd2);

	// check the best pair-end window with the highest diversity
	// report the average coverage
	double computeBestDivWindowPair(int numHaps, int minWinCover, int minFreq,
		double minFreqRatio, int& bestDist, int pos_interval, vector<int> d_candi,
		vector<int>& valid_candi);

};

#endif
