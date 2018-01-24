#include "partInfo.h"

// constructor
PartInfo::PartInfo() {
	colStat = NULL;
	invarCol = NULL;
	coverage = NULL;
	initNuclMap();
}

// destructor
PartInfo::~PartInfo() {
	deleteMem();
}

// delete all the memory used in the array
void PartInfo::deleteMem() {
	if (colStat != NULL)
		delete[] colStat;
	if (invarCol != NULL)
		delete[] invarCol;
	if (coverage != NULL)
		delete[] coverage;
}

void PartInfo::initNuclMap() {
	// A -> 1
	// C -> 2
	// G -> 3
	// T -> 4
	// others ( likely: '-' ) -> 0
	for (int i=0; i<256; i++)
		nuclMap[i] = 0;
	nuclMap['A'] = nuclMap['a'] = 1;
	nuclMap['C'] = nuclMap['c'] = 2;
	nuclMap['G'] = nuclMap['g'] = 3;
	nuclMap['T'] = nuclMap['t'] = 4;
}

// get the average coverage of the reads
double PartInfo::getAvgCoverage() {
	int i;
	int sumCover = 0;
	int numCover = 0;
	for (i=0; i<numCol; i++) {
		if (coverage[i] > 0) {
			numCover++;
			sumCover+= coverage[i];
		}
	}
	return (double) sumCover / numCover;
}

// load from the reads
// only for nucleotide reads
void PartInfo::loadFromReads(Reads& reads) {

	numCol = reads.maxColIndex()+1;
	int i,k;
	int c;
	int isInvar;

	// delete all the memory used in the array
	deleteMem();

	// create and initialize the array
	if (numCol > 0) {
		colStat = new int[numCol*5];
		invarCol = new int[numCol];
		coverage = new int[numCol];
		memset(colStat, 0, numCol*5*sizeof(int));
		memset(coverage, 0, numCol*sizeof(int));
		for (i=0; i<numCol; i++)
			invarCol[i] = -1;
	}
	
	// get the statistics
	// check the pair-end read
	Read *rd1, *rd2;
	int mpos1, mpos2;
	string s1, s2;
	int j1, j2;
	for (i=0; i<reads.pairRead1.size(); i++) {
		rd1 = reads.pairRead1[i];
		mpos1 = rd1->mapPos;
		s1 = rd1->seq;
		rd2 = reads.pairRead2[i];
		mpos2 = rd2->mapPos;
		s2 = rd2->seq;
		j1 = j2 = 0;
		while (j1 < s1.length() || j2 < s2.length()) {
			if (j1 >= s1.length() || j2+mpos2 < j1+mpos1) {
				// the region is only covered by read 2
				c = nuclMap[s2[j2]];
				k = (mpos2+j2)*5+c;
				if (k < numCol*5 && k >= 0)
					colStat[k]++;
				coverage[mpos2+j2]++;
				j2++;
			} else if (j2 >= s2.length() || j1+mpos1 < j2+mpos2) {
				// the region is only covered by read 1
				c = nuclMap[s1[j1]];
				k = (mpos1+j1)*5+c;
				if (k < numCol*5 && k >= 0)
					colStat[k]++;
				coverage[mpos1+j1]++;
				j1++;
			} else {
				// the region is covered by both read 1 and read 2
				if (s1[j1] == s2[j2]) {
					c = nuclMap[s1[j1]];
					k = (mpos1+j1)*5+c;
					if (k < numCol*5 && k >= 0)
						colStat[k]++;
					coverage[mpos1+j1]++;
				}
				j1++;j2++;
			}
		}
	}
	
	Read* rd;
	int mpos;
	string s;
	int j;
	// for the single-end read
	for (i=0; i<reads.unpairReads.size(); i++) {
		rd = reads.unpairReads[i];
		mpos = rd->mapPos;
		s = rd->seq;
		for (j=0; j<s.length(); j++) {
			c = nuclMap[s[j]];
			k = (mpos+j)*5+c;
			if (k < numCol*5 && k > 0) {
				colStat[k]++;
			}
			coverage[mpos+j]++;
		}
	}

	// compute the invarCol array
	for (i=0; i<numCol; i++) {
		c = -1;
		isInvar = 0;
		for (j=0; j<5; j++) {
			if (colStat[i*5+j]>0) {
				if (c==-1) {
					isInvar = 1;
					c = j;
				} else {
					isInvar = 0;
					break;
				}
			}
		}
		if (isInvar) {
			if (c == -1)
				c = 0;
			invarCol[i] = c;
		}
	}
}

// reset the statistics
// option = 1: using the seq
// option = 2: using the updatedSeq
void PartInfo::updateFromReads(Reads& reads, int option) {

	numCol = reads.maxColIndex()+1;
	int mpos;
	int i,j,k;
	int c;
	int isInvar;
	string s;

	// reset the array
	memset(colStat, 0, numCol*5*sizeof(int));
	memset(coverage, 0, numCol*sizeof(int));
	for (i=0; i<numCol; i++)
		invarCol[i] = -1;
	
	// get the statistics
	// for the first *updated* read
	for (i=0; i<reads.readNum; i++) {
		mpos = reads.readlist[i]->mapPos;
		if (option==1)
			s = reads.readlist[i]->seq;
		else
			s = reads.readlist[i]->updateSeq;
		for (j=0; j<s.length(); j++) {
			c = nuclMap[s[j]];
			k = (mpos+j)*5+c;
			if (k < numCol*5 && k > 0) {
				colStat[k]++;
			}
			coverage[mpos+j]++;
		}
	}
	
	// compute the invarCol array
	for (i=0; i<numCol; i++) {
		c = -1;
		isInvar = 0;
		for (j=0; j<5; j++) {
			if (colStat[i*5+j]>0) {
				if (c==-1) {
					isInvar = 1;
					c = j;
				} else {
					isInvar = 0;
					break;
				}
			}
		}
		if (isInvar) {
			if (c == -1)
				c = 0;
			invarCol[i] = c;
		}
	}
}

// print read information
void PartInfo::printReadInfo() {
	int i,j;
	for (i=0; i<numCol; i++) {
		for (j=0; j<5; j++) {
			cout << colStat[i*5+j] << " ";
		}
		cout << endl;
	}
}

// print read information
void PartInfo::printReadInfo(int i) {
	int j;
	for (j=0; j<5; j++) {
		cout << colStat[i*5+j] << " ";
	}
	cout << endl;
}

// normalize the coverage if it is too large
// after normaliation, if the freq < min_freq, set the freq to zero
void PartInfo::normalize(int max_thres, int min_freq) {
	int i,j;
	int max_j, max_cover;
	double ratio;
	int sum;
	for (i=0; i<numCol; i++) {
		if (coverage[i] > max_thres) {
			ratio = (double)max_thres / coverage[i];
			// check which is the largest
			max_cover = colStat[i*5];
			max_j = 0;
			for (j=1; j<5; j++) {
				if (colStat[i*5+j] > max_cover) {
					max_cover = colStat[i*5+j];
					max_j = j;
				}
			}
			// normalize the frequencies
			sum = 0;
			for (j=0; j<5; j++) {
				if (j==max_j)
					continue;
				colStat[i*5+j] = colStat[i*5+j] * ratio;
				if (colStat[i*5+j] < min_freq)
					colStat[i*5+j] = 0;
				sum += colStat[i*5+j];
			}
			colStat[i*5 + max_j] = max_thres - sum;
			coverage[i] = max_thres;
			if (sum==0) {
				invarCol[i] = max_j;
			} else {
				invarCol[i] = -1;
			}
		}
	}
}
