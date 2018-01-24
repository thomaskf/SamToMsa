#include "constructHaplo.h"

// constructor
Node::Node(Node* pre, int pattern, int combin) {
	this->pre = pre;
	this->pattern = pattern;
	this->combin = combin;
}

// identify all unavailable next assignments
// combin is zero-based
void CombinSet::identUnavailAssign(int pattern, int combin) {
	int i,j;
	// make the pattern to be unavailable for all combination
	for (j=0; j<numCombin; j++) {
		availAssign[pattern * numCombin + j] = 0;
	}
	// make the combinations which contain any item in "combin" to be unavailable
	// cout << "pattern = " << pattern << endl;
	// cout << "combin = " << combin << endl;
	for (j=0; j<numCombin; j++) {
		if ((j & combin) != 0) {
			for (i=0; i<numPattern; i++) {
				// cout << "set availAssign[" << i << "," << j << "] = 0" << endl;
				availAssign[i * numCombin + j] = 0;
			}
		}
	}
	// cout << "---------------" << endl;
}

// constructor
CombinSet::CombinSet(int pattern, int combin, long double logL, int numPattern, int numCombin, vector<Node*>& nodeCollection) {
	Node* pre = NULL;
	int item = 0;
	int i,j;
	
	// initialization
	this->numPattern = numPattern;
	this->numCombin = numCombin;
	numItems = 1;
	numHaps = bitset<32>(combin).count();
	cumLogL = logL;
	lastNode = new Node(NULL, pattern, combin);
	nodeCollection.push_back(lastNode);
	
	availAssign = new char[numPattern * numCombin];
	for (i=0; i<numPattern; i++) {
		for (j=0; j<numCombin; j++) {
			availAssign[i*numCombin + j] = 1;
		}
	}
	
	// set up the availAssign array
	identUnavailAssign(pattern, combin);
}

// constructor
// combin is zero-based
CombinSet::CombinSet(CombinSet* aCombinSet, int pattern, int combin, long double logL, vector<Node*>& nodeCollection) {
	numPattern = aCombinSet->numPattern;
	numCombin = aCombinSet->numCombin;
	availAssign = new char[numPattern * numCombin];
	memcpy(availAssign, aCombinSet->availAssign, numPattern * numCombin * sizeof(char));
	
	lastNode = new Node(aCombinSet->lastNode, pattern, combin);
	nodeCollection.push_back(lastNode);
	
	numItems = aCombinSet->numItems + 1;
	numHaps = aCombinSet->numHaps + bitset<32>(combin).count();
	cumLogL = aCombinSet->cumLogL + logL;
	identUnavailAssign(pattern, combin);
}

// destructor
CombinSet::~CombinSet() {
	delete[] availAssign;
}

// show the content of combination
// combin is zero-based
void showCombin(int combin) {
	int i = 0;
	int num = 0;
	int c = combin;
	cout << "{";
	while (c > 0) {
		if ((c & 1)==1) {
			if (num > 0)
				cout << ",";
			cout << i;
			num++;
		}
		c = (c >> 1);
		i++;
	}
	cout << "}";
}

void CombinSet::show() {
	Node* cNode = lastNode;
	
	while (cNode != NULL) {
		if (cNode->combin>0) {
			cout << "H";
			showCombin(cNode->combin);
		} else {
			cout << "NIL";
		}
		cout << "[S" << cNode->pattern << "] ";
		cNode = cNode->pre;
	}
	cout << " " << cumLogL;
	cout << endl;
}

// show the content of availAssign
void CombinSet::showAvailAssign() {
	int i,j,k;
	for (i=0; i<numPattern; i++) {
		k = i*numCombin;
		for (j=0; j<numCombin; j++) {
			cout << (int) availAssign[k+j] << " ";
		}
		cout << endl;
	}
}

// collect items and patterns
void CombinSet::collect(vector<int>& items, vector<int>& patterns) {
	int pattern;
	int item;
	int c;
	int i;
	Node* cNode = lastNode;
	items.clear();
	patterns.clear();
	while (cNode != NULL) {
		pattern = cNode->pattern;
		c = cNode->combin;
		i = 0;
		while (c > 0) {
			if ((c & 1) == 1) {
				items.push_back(i);
				patterns.push_back(pattern);
			}
			c = (c >> 1);
			i++;
		}
		cNode = cNode->pre;
	}
	/*
	// show the content of items and patterns
	for (i=0; i<items.size(); i++) {
		cout << items[i] << " -> " << patterns[i] << endl;
	}
	exit(1);
	 */
}

// constructor
BestCombin::BestCombin() {
	logL = NULL;
	totSize = 0;
}

// destructor
BestCombin::~BestCombin() {
	if (logL != NULL)
		delete[] logL;
}

// update dimension
void BestCombin::resize(int numC, int numP, int numI) {
	numCombin = numC;
	numPattern = numP;
	numItem = numI;
	int newSize = numC * numP;
	if (totSize < newSize) {
		if (logL != NULL)
			delete[] logL;
		logL = new long double[newSize];
	}
	totSize = newSize;
}

// compute the best combination with the maximum sum of log-likelihood
bool BestCombin::compute(vector<int>& items, vector<int>& patterns, long double& bestLogL, long double& best2LogL) {
	priority_queue<CombinSet*, vector<CombinSet*>, CompareCombinSet> allCombinSets;
	vector<Node*> nodeCollections;
	int i,j,k;
	CombinSet* currCombinSet;
	CombinSet* newCombinSet;
	long double maxLogL;
	int ansNum = 0;
	
	/*
	// show the array
	for (i=0; i<numPattern; i++) {
		for (j=0; j<numCombin; j++) {
			cout << logL[i*numCombin + j] << "\t";
		}
		cout << endl;
	}*/

	// initialization
	// only consider the first pattern
	i=0;
	for (j=0; j<numCombin; j++) {
		allCombinSets.push(new CombinSet(i, j, logL[i*numCombin + j], numPattern, numCombin, nodeCollections));
	}

	while (ansNum < 2) {
		currCombinSet = allCombinSets.top();
		allCombinSets.pop();
		// currCombinSet->show();
		// currCombinSet->showAvailAssign();
		i=currCombinSet->numItems;
		if (i == numPattern) {
			if (currCombinSet->numHaps == numItem) {
				// this is a valid answer
				if (ansNum==0) {
					// the first answer
					bestLogL = currCombinSet->cumLogL;
					// cout << "the above combination is the first answer (" << bestLogL << ")" << endl;
					currCombinSet->collect(items, patterns);
					if (numPattern == 1) {
						// second answer is not available
						best2LogL = bestLogL - (DELTA_SUBSEQ_LOGL*100);
						ansNum++;
					}
				} else {
					// the second answer
					best2LogL = currCombinSet->cumLogL;
					// cout << "the above combination is the second answer (" << best2LogL << ")" << endl;
				}
				ansNum++;
			}
		} else {
			k = i*numCombin;
			for (j=0; j<numCombin; j++) {
				if (currCombinSet->availAssign[k+j]==1) {
					// cout << "(" << i << "," << j << ") is selected" << endl;
					newCombinSet = new CombinSet(currCombinSet, i, j, logL[i*numCombin + j], nodeCollections);
					// newCombinSet->show();
					// exit(1);
					allCombinSets.push(newCombinSet);
				}
			}
		}
		delete(currCombinSet);
	}
	
	// delete all the items inside priority list
	while (!allCombinSets.empty()) {
		allCombinSets.pop();
	}
	
	// delete the node collections
	for (i=0; i<nodeCollections.size(); i++)
		delete(nodeCollections[i]);
	
	if (bestLogL - best2LogL >= DELTA_SUBSEQ_LOGL)
		return true;
	else
		return false;
}

ConstructHaplo::ConstructHaplo() {
	haploRatios=NULL;
	hapRateCombins=NULL;
	baseCombinL=NULL;
	logF=NULL;
	logFF=NULL;
	sRateCombins=NULL;
	slogF=NULL;
	slogFF=NULL;
	
	initializeCumLogI();

	// initialize nucl_map
	int i;
	for (i=0; i<256; i++)
		nucl_map[i] = 0;
	nucl_map['A'] = nucl_map['a'] = 1;
	nucl_map['C'] = nucl_map['c'] = 2;
	nucl_map['G'] = nucl_map['g'] = 3;
	nucl_map['T'] = nucl_map['t'] = 4;
}

// destructor
ConstructHaplo::~ConstructHaplo() {
	deleteMem();
}

// delete the memory of the array
void ConstructHaplo::deleteMem() {
	if (haploRatios!=NULL)
		delete[] haploRatios;
	if (hapRateCombins!=NULL)
		delete[] hapRateCombins;
	if (baseCombinL!=NULL)
		delete[] baseCombinL;
	if (logF!=NULL)
		delete[] logF;
	if (logFF!=NULL)
		delete[] logFF;
	if (sRateCombins!=NULL)
		delete[] sRateCombins;
	if (slogF!=NULL)
		delete[] slogF;
	if (slogFF!=NULL)
		delete[] slogFF;
}

/*
// return the max orig column index
int ConstructHaplo::maxOrigColIndex(int hapID) {
	int i;
	int m = 0;
	Read* rd;
	for (i=0; i<hapReads[hapID].size(); i++) {
		rd = rds.readlist[hapReads[hapID].at(i)];
		if (m < rd->mapOrigPos+rd->origSeq.length()) {
			m = rd->mapOrigPos+rd->origSeq.length()-1;
		}
	}
	return m;
}
*/

// load the read file
void ConstructHaplo::loadReadFile(char* readFile) {
	rds.readSamFile(readFile);
	rds.updateForInserts();
	
/*
	// for debugging
	// sort the reads
	rds.sortReads();
	// print out the reads
	int fr_pos = 0;
	int to_pos = 100;
	int k;
	for (k=0; k<20000; k+=100) {
		cout << "[" << fr_pos+k << " -- " << to_pos+k << "]" << endl;
		rds.printReads(fr_pos+k, to_pos+k);
	}
	exit(1);
*/	
	partInfo.loadFromReads(rds);
	numCols = partInfo.numCol;
	numRds = rds.readNum;
}

// initialise the values of haplotype frequency ratios
// and allocate the memory according to the number of frequencies inside the file
void ConstructHaplo::loadFreqFile(char* freqFile) {
	double sum;
	vector<double> freqs;
	string aline;
	ifstream fin;
	int i;
	// double f, err;
	// read the file
	fin.open(freqFile);
	while (getline(fin, aline)) {
		if (aline.length() > 0) {
			freqs.push_back(atof(aline.c_str()));
		}
	}
	fin.close();
	
	numHap = freqs.size();
	numCombin = pow(2,numHap);
	
	if (numHap > MAX_HAP_NUM) {
		cout << "Error! The number of haplotypes cannot be more than " << MAX_HAP_NUM << endl;
		exit(1);
	}

	// delete the memory
	deleteMem();

	// allocate the memory
	haploRatios=new long double[numHap];
	hapRateCombins = new long double[numCombin];
	baseCombinL = new long double[numCols*numCombin*5];
	logF = new long double[numCombin];
	logFF = new long double[numCombin];
	
	// compute the sum
	sum = 0.0;
	for (i=0; i<freqs.size(); i++)
		sum += freqs[i];
		
	// assign the haplotype freq ratios
	for (i=0; i<freqs.size(); i++) {
		haploRatios[i] = freqs[i] / sum;
	}
	/*
	// assign the haplotype freq ratios (with consideration of sequence error rate)
	cerr << "After considering the errors, the frequency of each haplotype:" << endl;
	for (i=0; i<freqs.size(); i++) {
		f = freqs[i] / sum;
		err = ((1.0-f)/3.0 - f) * SEQ_ERR_RATE;
		haploRatios[i] = f + err;
		cerr << "Haplotype " << i+1 << ": " << haploRatios[i] << endl;
	}
	exit(1);
	 */
}

// correct the reads
// and normalize the frequencies
void ConstructHaplo::readCorrection() {
	int i,k;
	int m;
	Read* rd;
	int mpos;
	char int2Nucl[] = {'-','A','C','G','T'};
	
	// first normalize the frequencies
	partInfo.normalize(MAX_COVER, MIN_FREQ);

	// update the reads
	for (i=0; i<numRds; i++) {
		rd = rds.readlist[i];
		mpos = rd->mapPos;
		m=0;
		for (k=0; k<rd->seq.length(); k++) {
			if (partInfo.invarCol[mpos+k] >= 0) {
				// it is an invariable site
				rd->seq[k] = int2Nucl[partInfo.invarCol[mpos+k]];
			}
			if (partInfo.invarCol[mpos+k] != 0) {
				// we cannot delete this char
				if (m < k) {
					rd->seq[m] = rd->seq[k];
				}
				m++;
			}
		}
		if (m==0)
			rd->seq = "";
		else if (m < rd->seq.length())
			rd->seq.resize(m);
	}
	
	// get the relationship between old mapping position and the new mapping position
	int* newPos = new int[numCols];
	m=0;
	for (i=0; i<numCols; i++) {
		if (partInfo.invarCol[i] != 0) {
			// we cannot delete this char
			newPos[i] = m;
			m++;
		} else {
			newPos[i] = -1;
		}
	}
	cerr << "Number of gaps deleted: " << numCols - m << endl;
	
	// update the mapping position
	for (i=0; i<numRds; i++) {
		rd = rds.readlist[i];
		rd->mapPos = newPos[rd->mapPos];
	}
}

// combine the overlap reads of the same pair into a long read
// and update the frequencies
void ConstructHaplo::combineRdAndUpdateFreq() {
	rds.combineReads();
	// update the frequencies
	partInfo.loadFromReads(rds);
	numCols = partInfo.numCol;
}

// calculate the haplo ratio for each combinition
void ConstructHaplo::calHapRateCombins() {
	int i,j;
	long double sumRatio;
	long double error;
	
	for (i=0; i<numCombin; i++) {
		sumRatio = 0.0;
		for (j=0; j<numHap; j++) {
			if ((i >> j) & 1) {
				// this combination has haplotype j
				// cerr << j << ",";
				sumRatio += haploRatios[j];
			}
		}
		// consider the sequencing errors
		error = ((1.0 - sumRatio) / 3.0 - sumRatio) * SEQ_ERR_RATE;
		hapRateCombins[i] = sumRatio + error;
		// cerr << "\t" << hapRateCombins[i] << endl;
		logF[i] = logl(hapRateCombins[i]);
		logFF[i] = logl(1.0-hapRateCombins[i]);
		// cerr << i << " " << logF[i] << " " << logFF[i] << endl;
	}
}

// calculate P(ci | Hj) for all bases of all columns and combinations
void ConstructHaplo::calBaseCombinL() {
	int c, i, x, m, n, totNum;
	long double logC;
	
	memset(baseCombinL, 0, numCols*numCombin*5*sizeof(long double));
	for (c=0; c<numCols; c++) {
		// c = 1883;
		totNum = partInfo.coverage[c];
		for (x=0; x<5; x++) {
			// if (partInfo.colStat[c*5+x] > 0) {
				m = partInfo.colStat[c*5+x];
				n = totNum - m;
				logC = logCombin(totNum, m);
				for (i=0; i<numCombin; i++) {
					// baseCombinL[c*5*numCombin + x*numCombin + i] = expl(logC + m * logF[i] + n * logFF[i]);
					baseCombinL[c*5*numCombin + x*numCombin + i] = logC + m * logF[i] + n * logFF[i];
				}
			// }
			// cout << partInfo.colStat[c*5+x] << "\t";
		}
		
		/*
		cout << endl;
		for (x=0; x<5; x++) {
			for (i=0; i<numCombin; i++) {
				cout << baseCombinL[c*5*numCombin + x*numCombin + i] << "\t";
			}
			cout << endl;
		}
		exit(1);
		*/
	}
}

#ifdef USE_OLD_METHOD1
// construct the haplotypes based on baseCombinL
// OLD VERSION
void ConstructHaplo::constructHaplos() {
	
	int i, k, x, h;
	int bestX;
	int bestI;
	long double currL;
	long double bestL;
	char int2Nucl[] = {'-','A','C','G','T'};
	
	haplos.resize(numHap);
	for (h=0; h<numHap; h++)
		haplos[h]="";
	for (k=0; k<numCols; k++) {
		if (partInfo.invarCol[k]>=0) {
			for (h=0; h<numHap; h++) {
				haplos[h].append(1, int2Nucl[partInfo.invarCol[k]]);
			}
		} else {
			for (h=0; h<numHap; h++) {
				bestI = -1;
				for (x=0; x<5; x++) {
					for (i=0; i<numCombin; i++) {
						if (((i+1) >> h) & 1) {
							currL = baseCombinL[k*5*numCombin + x*numCombin + i];
							if (bestI==-1 || currL > bestL) {
								bestL = currL;
								bestI = i;
								bestX = x;
							}
						}
					}
				}
				// cout << k << " " << h << " " << bestL << " " << logl(bestL) << endl;
				// cout << k << " " << h << " " << bestL << endl;
				if (bestI != -1) {
					haplos[h].append(1,int2Nucl[bestX]);
				} else {
					haplos[h].append(1,'N');
				}
			}
		}
	}
	hapLen = haplos[0].length();
}

#else

// construct the haplotypes based on baseCombinL
// NEW VERSION
void ConstructHaplo::constructHaplos() {
	
	int i, k, x, h;
	int bestX;
	int bestI;
	long double currL;
	long double bestL;
	long double best2L;
	char int2Nucl[] = {'-','A','C','G','T'};
	vector<int> items;
	vector<int> patterns;
	
	haplos.resize(numHap);
	BestCombin bestCombin;
	bestCombin.numCombin = numCombin;
	bestCombin.numPattern = 5;
	bestCombin.numItem = 3;
	
	// cout << "bestCombin.numCombin = " << bestCombin.numCombin << endl;
	
	for (h=0; h<numHap; h++)
		haplos[h]="";
	for (k=0; k<numCols; k++) {
		if (partInfo.invarCol[k]>=0) {
			for (h=0; h<numHap; h++) {
				haplos[h].append(1, int2Nucl[partInfo.invarCol[k]]);
			}
		} else {
			// set the log likelihood values
			bestCombin.logL = &(baseCombinL[k*5*numCombin]);
			bestCombin.compute(items, patterns, bestL, best2L);
			for (h=0; h<numHap; h++) {
				haplos[items[h]].append(1,int2Nucl[patterns[h]]);
			}
			// reset to NULL
			bestCombin.logL = NULL;
		}
	}
	hapLen = haplos[0].length();
}

#endif

// distance between haplo and the read
int ConstructHaplo::distHapRead(string& read, string& hap, int mpos) {
	int d = 0;
	int i;
	for (i=0; i<read.length(); i++) {
		d += (read[i] != hap[mpos+i]);
	}
	return d;
}

// identify and delete those reads have too many errors
void ConstructHaplo::deleteErrReads() {
	int i,j;
	int d;
	string readSeq;
	int mpos,bestDist;
	Read* rd;
	int numRdDiscarded = 0;
	// check the smallest distance for each read
	for (i=0; i<numRds; i++) {
		rd = rds.readlist[i];
		readSeq = rd->seq;;
		mpos = rd->mapPos;
		// check with the 1st haplotype
		bestDist = distHapRead(readSeq, haplos[0], mpos);
		for (j=1; j<numHap; j++) {
			// check with the other haplotypes
			d = distHapRead(readSeq, haplos[j], mpos);
			if (d < bestDist) {
				bestDist = d;
			}
		}
		// check whether the read has too many errors
		if (bestDist > MANY_ERR) {
			rd->isDiscard = true;
			numRdDiscarded++;
		}
	}
	// delete the discarded reads
	j=0;
	for (i=0; i<numRds; i++) {
		if (!rds.readlist[i]->isDiscard) {
			if (j < i)
				rds.readlist[j] = rds.readlist[i];
			j++;
		}
	}
	rds.readlist.resize(j);
	numRds = j;
	cerr << "Number of reads discarded due to too many errors: " << numRdDiscarded << endl;
}

// supporting function
// check the best hits for the read
// for single read, set rd2 = NULL
void ConstructHaplo::getBestHits(vector<int>& bestH, Read* rd, Read* rd2) {
	int mpos,bestDist,d;
	string readSeq;
	int mpos2,d2;
	string readSeq2;
	int j;

	readSeq = rd->seq;
	mpos = rd->mapPos;
	// check with the 1st haplotype
	bestDist = distHapRead(readSeq, haplos[0], mpos);
	if (rd2 != NULL) {
		readSeq2 = rd2->seq;
		mpos2 = rd2->mapPos;
		bestDist += distHapRead(readSeq2, haplos[0], mpos2);
	}
	// cout << i << "\t" << 0 << "\t" << bestDist << endl;
	bestH.clear();
	bestH.push_back(0);
	for (j=1; j<numHap; j++) {
		// check with the other haplotypes
		d = distHapRead(readSeq, haplos[j], mpos);
		if (rd2 != NULL) {
			d += distHapRead(readSeq2, haplos[j], mpos2);
		}
		// cout << i << "\t" << j << "\t" << d << endl;
		if (d < bestDist) {
			bestDist = d;
			bestH.clear();
			bestH.push_back(j);
		} else if (d == bestDist) {
			bestH.push_back(j);
		}
	}
}

// get the reads beloging to each haplotype
void ConstructHaplo::getHapReads() {
	int i,j;
	vector<int> bestH;

	// reset hapReads
	hapReads.resize(numHap);
	for (i=0; i<numHap; i++)
		hapReads[i].clear();
	
	// assign the IDs for the reads
	rds.assignReadIDs();
	
	// for pair-end reads
	for (i=0; i<rds.pairRead1.size(); i++) {
		getBestHits(bestH, rds.pairRead1[i], rds.pairRead2[i]);
		for (j=0; j<bestH.size(); j++) {
			hapReads[bestH[j]].push_back(rds.pairRead1[i]->id);
			hapReads[bestH[j]].push_back(rds.pairRead2[i]->id);
		}
	}
	
	// for single-end reads
	for (i=0; i<rds.unpairReads.size(); i++) {
		getBestHits(bestH, rds.unpairReads[i], NULL);
		for (j=0; j<bestH.size(); j++) {
			hapReads[bestH[j]].push_back(rds.unpairReads[i]->id);
		}
	}
	// exit(1);
}

// show the sequences of the reads for a particular haplotype
void ConstructHaplo::showReads(int haploID) {
	int i;
	Read* rd;
	for (i=0; i<hapReads[haploID].size(); i++) {
		rd = rds.readlist[hapReads[haploID].at(i)];
		if (rd->mapPos >0) {
			cout << string(rd->mapPos, ' ');
		}
		cout << rd->seq << " [" << i << "]" << endl;
	}
}

// get the hits of the read to each estimated haplotype
void ConstructHaplo::getHapReadHits() {
	int i, j, k,p;
	int mpos;
	string readSeq;
	Read* rd;

	// initMapping();
	hapReadHits.resize(numHap);
	for (i=0; i<numHap; i++) {
		hapReadHits[i].clear();
		hapReadHits[i].insert(hapReadHits[i].begin(),5*numCols,0);
	}

	for (i=0; i<numHap; i++) {
		// for haplotype i

		// update the hits
		for (j=0; j<hapReads[i].size(); j++) {
			rd = rds.readlist[hapReads[i].at(j)];
			mpos = rd->mapPos;
			readSeq = rd->seq;
			for (k=0; k<readSeq.length(); k++) {
				p = (mpos+k)*5 + nucl_map[readSeq[k]];
				hapReadHits[i].at(p)++;
			}
		}
	}
}

// get the full haplotypes based on hapReadHits
void ConstructHaplo::getFullHaplos() {
	int i,j,k;
	int bestChar, bestNumHit;
	char int2Nucl[] = {'-','A','C','G','T'};
	
	fullHaplos.resize(numHap);
	for (i=0; i<numHap; i++) {
		fullHaplos[i] = string(numCols, 'N');
	}

	for (i=0; i<numHap; i++) {
		for (j=0; j<numCols; j++) {
			if (fullHaplos[i].at(j)=='N') {
				bestChar = 'N';
				bestNumHit = 0;
				for (k=0; k<5; k++) {
					if (hapReadHits[i].at(j*5+k) > bestNumHit) {
						bestNumHit = hapReadHits[i].at(j*5+k);
						bestChar = int2Nucl[k];
					}
				}
				fullHaplos[i].at(j) = bestChar;
			}
		}
		
#ifndef STRICT
		// heuristic method
		for (j=0; j<numCols;j++) {
			if (fullHaplos[i].at(j)=='N') {
				fullHaplos[i].at(j) = haplos[i].at(j);
			}
		}
#endif		
	}
}

// remove the gaps and the leading/ending N's inside the first-version haplotype
void ConstructHaplo::rmGapInsideHaplo() {
	int i,j,k;
	bool toDelete;
	for (i=0; i<numHap; i++) {
		k=0;
		for (j=0; j<haplos[i].length(); j++) {
			toDelete = false;
			if (haplos[i].at(j)=='-')
				toDelete = true;
			else if (k==0 && haplos[i].at(j)=='N')
				toDelete = true;
			if (!toDelete) {
				if (k < j) {
					haplos[i].at(k) = haplos[i].at(j);
				}
				k++;
			}
		}
		while (k>0 && haplos[i].at(k-1)=='N')
			k--;
		haplos[i].resize(k);
	}
}

// remove the gaps and the leading/ending N's inside the full haplotype
void ConstructHaplo::rmGapInsideFullHaplo() {
	int i,j,k;
	bool toDelete;
	for (i=0; i<numHap; i++) {
		k=0;
		for (j=0; j<fullHaplos[i].length(); j++) {
			toDelete = false;
			if (fullHaplos[i].at(j)=='-')
				toDelete = true;
			else if (k==0 && fullHaplos[i].at(j)=='N')
				toDelete = true;
			if (!toDelete) {
				if (k < j) {
					fullHaplos[i].at(k) = fullHaplos[i].at(j);
				}
				k++;
			}
		}
		while (k>0 && fullHaplos[i].at(k-1)=='N')
			k--;
		fullHaplos[i].resize(k);
	}
}

// remove the gaps and the leading/ending N's inside the full haplotype
// and remove the region with the too-low-coverage near the ends of the full haplotype
void ConstructHaplo::rmGapAndLowCoverInsideFullHaplo() {
	int i,j,k;
	bool toDelete;
	double avgCoverage = partInfo.getAvgCoverage();
	string tooLowCover;
	for (i=0; i<numHap; i++) {
		k=0;
		tooLowCover = string(fullHaplos[i].length(),0);
		for (j=0; j<fullHaplos[i].length(); j++) {
			toDelete = false;
			if (fullHaplos[i].at(j)=='-')
				toDelete = true;
			else if (k==0 && fullHaplos[i].at(j)=='N')
				toDelete = true;
			else if (k==0 && (double) partInfo.coverage[j] < avgCoverage * COVER_THRES_NEAR_END)
				toDelete = true;
			if (!toDelete) {
				if (k < j) {
					fullHaplos[i].at(k) = fullHaplos[i].at(j);
				} 
				if ((double) partInfo.coverage[j] < avgCoverage * COVER_THRES_NEAR_END)
					tooLowCover[k]=1;
				k++;
			}
		}
		while (k>0 && (fullHaplos[i].at(k-1)=='N'
				|| tooLowCover[k]==1))
					k--;
		fullHaplos[i].resize(k);
	}
}

// update the reads according to the estimated haplotypes
void ConstructHaplo::updateReads() {
	int i,j;
	int d;
	int min_dist;
	vector<int> min_j;
	Read* currRd;
	string s;
	bool allSame;
	for (i=0; i<numRds; i++) {
		currRd = rds.readlist[i];
		min_j.clear();
		min_j.push_back(0);
		min_dist = distHapRead(currRd->seq, fullHaplos[0], currRd->mapPos);
		for (j=1; j<numHap; j++) {
			d = distHapRead(currRd->seq, fullHaplos[j], currRd->mapPos);
			if (d == min_dist) {
				min_j.push_back(j);
			} else if (d < min_dist) {
				min_j.clear();
				min_j.push_back(j);
				min_dist = d;
			}
		}
		// cout << i << "\t" << currRd->mapPos << "\t" << min_dist << "\t";
		/*
		for (j=0; j<min_j.size(); j++) {
			if (j>0)
				cout << ",";
			cout << min_j[j];
		}
		cout << endl;
		 */ 
		if (min_j.size()==1) {
			currRd->updateSeq = fullHaplos[min_j[0]].substr(currRd->mapPos, currRd->seq.length());
		} else {
			s = fullHaplos[min_j[0]].substr(currRd->mapPos, currRd->seq.length());
			allSame = true;
			for (j=1; j<min_j.size(); j++) {
				if (fullHaplos[min_j[j]].substr(currRd->mapPos, currRd->seq.length())!=s) {
					allSame = false;
					break;
				}
			}
			if (allSame)
				currRd->updateSeq = s;
			else
				currRd->updateSeq = currRd->seq;
		}
		
		currRd->seq = currRd->updateSeq;
		
		/*
		// for debugging
		if (currRd->seq != currRd->updateSeq) {
			cout << "0> " << currRd->seq << endl;
			cout << "1> " << currRd->updateSeq << endl;
		}
		 */
	}
}


// list out the haplotypes
void ConstructHaplo::listHaplos() {
	int i;
	for (i=0; i<numHap; i++) {
		cout << ">1_" << i+1 << endl;
		cout << haplos[i] << endl;
	}
}

// list out the full haplotypes
void ConstructHaplo::listFullHaplos() {
	int i;
	for (i=0; i<numHap; i++) {
		cout << ">" << VERSION << "_" << i+1 << endl;
		cout << fullHaplos[i] << endl;
	}
}

//---------------------------------------------------------
// The methods considering windows
//---------------------------------------------------------

// calculate the haplo ratio for each combinition
void ConstructHaplo::calRateForSubSeqCombins() {
	int i,j;
	long double sumRatio;
	long double error;
	
	// allocate the memories for the array
	if (sRateCombins==NULL)
		sRateCombins = new long double[numCombin];
	if (slogF==NULL)
		slogF = new long double[numCombin];
	if (slogFF==NULL)
		slogFF = new long double[numCombin];
	
	for (i=0; i<numCombin; i++) {
		sumRatio = 0.0;
		for (j=0; j<numHap; j++) {
			if (i >> j & 1) {
				// this combination has haplotype j
				// cerr << j << ",";
				sumRatio += haploRatios[j];
			}
		}
		if (sumRatio == 0.0)
			sRateCombins[i] = ERR_SUBSEQ_RATE; // / WIN_SIZE;
		else
			sRateCombins[i] = sumRatio * (1.0 - ERR_SUBSEQ_RATE);
		// cerr << "\t" << sRateCombins[i] << endl;
		slogF[i] = logl(sRateCombins[i]);
		slogFF[i] = logl(1.0-sRateCombins[i]);
	}
}

// get the underlining frequencies of the haplotypes (using window)
void ConstructHaplo::getHaploFreq(int winSize, vector<double>& hapFreqs, int& winNum) {
	int i,k;
	double tot;
	int fr_pos, to_pos;
	vector<string> subSeq;
	vector<double> subSeqFreq;
	hapFreqs.resize(numHap);
	for (k=0; k<numHap; k++)
		hapFreqs[k]=0.0;
	winNum=0;
	for (k=0; k<numCols-winSize+1; k++) {
		fr_pos = k; to_pos = fr_pos+winSize-1;
		tot = rds.getSubSeqs(fr_pos, to_pos, MIN_FREQ_FOR_SUBSEQ, MIN_FREQ_RATIO_FOR_SUBSEQ, subSeq, subSeqFreq);
		if (subSeq.size() == numHap) {
			sort(subSeqFreq.begin(), subSeqFreq.end());
			cout << fr_pos;
			for (i=0; i<numHap; i++) {
				hapFreqs[i] += subSeqFreq[i];
				cout << "\t" << subSeqFreq[i];
			}
			cout << endl;
			winNum++;
		}
	}
	for (i=0; i<numHap; i++)
		hapFreqs[i] = hapFreqs[i] / winNum;
}

// get the underlining frequencies of the haplotypes (using window-pairs)
void ConstructHaplo::getHaploFreq(int winDist, int winSize, vector<double>& hapFreqs, int& winNum) {
	int i,k;
	double tot;
	int fr_pos1, fr_pos2, to_pos1, to_pos2;
	vector<string> subSeq1;
	vector<string> subSeq2;
	vector<double> subSeqFreq;
	hapFreqs.resize(numHap);
	for (k=0; k<numHap; k++)
		hapFreqs[k]=0.0;
	winNum=0;
	for (k=0; k<numCols-winSize-winDist+1; k++) {
		fr_pos1 = k; to_pos1 = fr_pos1+winSize-1;
		fr_pos2 = to_pos1+winDist; to_pos2 = fr_pos2+winSize-1;
		tot = rds.getSubSeqPairs(fr_pos1, to_pos1, fr_pos2, to_pos2, MIN_FREQ_FOR_SUBSEQ,
				MIN_FREQ_RATIO_FOR_SUBSEQ, subSeq1, subSeq2, subSeqFreq);
		if (subSeq1.size() == numHap) {
			sort(subSeqFreq.begin(), subSeqFreq.end());
			cout << "[" << fr_pos1 << "-" << fr_pos1+winSize-1 << "],[" << fr_pos2 << "-" << fr_pos2+winSize-1 << "]";
			for (i=0; i<numHap; i++) {
				hapFreqs[i] += subSeqFreq[i];
				cout << "\t" << subSeqFreq[i];
			}
			cout << endl;
			winNum++;
		}
	}
	for (i=0; i<numHap; i++)
		hapFreqs[i] = hapFreqs[i] / winNum;
}

// get the set of pairs of subsequences belonging to each haplotype
void ConstructHaplo::getSubSeqPairForHaplos(int winDist, int winSize) {
	vector<string> subSeq1;
	vector<string> subSeq2;
	vector<double> subSeqFreq;
	
	double tot;
	int fr_pos1, fr_pos2, to_pos1, to_pos2;
	int i,j,k,l,m,n;
	vector<int> items;
	vector<int> patterns;
	BestCombin bestCombin;
	int numPatterns;

	long double logC;
	long double logL;
	long double maxLogL;
	long double max2LogL;

	// reset the arrays
	hapSubSeqs.resize(numHap);
	hapSubSeqPos.resize(numHap);
	hapSubSeqWeight.clear();
	for (i=0; i<numHap; i++) {
		hapSubSeqs[i].clear();
		hapSubSeqPos[i].clear();
	}
	
	for (k=0; k<numCols-winSize-winDist+1; k++) {
		#ifdef SHOW_DEBUG_MESG
			// for debugging
			fr_pos1 = 2755;
		#else
			fr_pos1 = k;
		#endif
		to_pos1 = fr_pos1+winSize-1;
		fr_pos2 = to_pos1+winDist; to_pos2 = fr_pos2+winSize-1;
		tot = rds.getSubSeqPairs(fr_pos1, to_pos1, fr_pos2, to_pos2, MIN_FREQ_FOR_SUBSEQ,
				MIN_FREQ_RATIO_FOR_SUBSEQ, subSeq1, subSeq2, subSeqFreq);
	
		#ifdef SHOW_DEBUG_MESG
		cout << "tot = " << tot << endl;
		for (i=0; i<subSeq1.size(); i++) {
			cout << subSeq1[i] << " " << subSeq2[i] << " {" << subSeqFreq[i] << "} " << (double) subSeqFreq[i] / tot << endl;
		}
		#endif

		numPatterns = subSeq1.size();
		bestCombin.resize(numCombin, numPatterns, numHap);
		
		// set up the logL array
		for (i=0; i<numPatterns; i++) {
			m = (int) (subSeqFreq[i]+0.5);
			n = tot - m;
			logC = logCombin(tot, m);
			for (j=0; j<numCombin; j++) {
				logL = logC + m * slogF[j] + n * slogFF[j];
				bestCombin.logL[i*numCombin + j] = logL;
			}
		}
		
		#ifdef SHOW_DEBUG_MESG
		// show the log-likelihood values
		for (i=0; i<numPatterns; i++) {
			for (j=0; j<numCombin; j++) {
				cout << bestCombin.logL[i*numCombin + j] << "\t";
			}
			cout << endl;
		}
		#endif

		// compute the best combination
		bestCombin.compute(items, patterns, maxLogL, max2LogL);
		#ifdef USE_WEIGHT_APPROACH
		if (maxLogL < MIN_SUBSEQ_LOGL_THRES) {
		#else
		if (maxLogL < MIN_SUBSEQ_LOGL_THRES || maxLogL - max2LogL < DELTA_SUBSEQ_LOGL) {
		#endif

			#ifdef SHOW_DEBUG_MESG
				cout << "The difference between maxLogL and max2LogL is too small: " << maxLogL - max2LogL << " tot=" << tot << endl;
				exit(1);
			#else
				continue;
			#endif
		}
		
		for (j=0; j<numHap; j++) {
			hapSubSeqs[items[j]].push_back(subSeq1[patterns[j]]);
			hapSubSeqPos[items[j]].push_back(fr_pos1);
			hapSubSeqs[items[j]].push_back(subSeq2[patterns[j]]);
			hapSubSeqPos[items[j]].push_back(fr_pos2);
		}
		
		#ifdef USE_WEIGHT_APPROACH
			hapSubSeqWeight.push_back(maxLogL - max2LogL);
			hapSubSeqWeight.push_back(maxLogL - max2LogL);
		#else
			hapSubSeqWeight.push_back(1.0);
			hapSubSeqWeight.push_back(1.0);
		#endif
		
		#ifdef SHOW_DEBUG_MESG
			cout << "Resulting logL = " << maxLogL << endl;
			cout << "Diff logL = " << maxLogL - max2LogL << endl;
			for (j=0; j<numHap; j++) {
				cout << "pos: " << hapSubSeqPos[j].at(hapSubSeqPos[j].size()-2);
				cout << "," << hapSubSeqPos[j].at(hapSubSeqPos[j].size()-1);
				cout << " H" << j+1 << ": " << hapSubSeqs[j].at(hapSubSeqs[j].size()-2);
				cout << "," << hapSubSeqs[j].at(hapSubSeqs[j].size()-1);
				cout << " tot=" << tot << endl;
			}
		#endif

		#ifdef SHOW_DEBUG_MESG
			exit(1);
		#endif
	}
}

// (new algorithm) get the set of subsequences belonging to each haplotype
void ConstructHaplo::getSubSeqForHaplos() {
	
	int i,j,k,l,m,n;
	int tot;
	int fr_pos, to_pos;
	long double logC;
	long double logL;
	long double maxLogL;
	long double max2LogL;

	vector<string> subSeq;
	vector<double> subSeqFreq;
	
	vector<int> items;
	vector<int> patterns;

	BestCombin bestCombin;
	int numPatterns;

	// reset the arrays
	hapSubSeqs.resize(numHap);
	hapSubSeqPos.resize(numHap);
	hapSubSeqWeight.clear();
	for (i=0; i<numHap; i++) {
		hapSubSeqs[i].clear();
		hapSubSeqPos[i].clear();
	}

	#ifdef SHOW_DEBUG2_MESG
		int debug_pos = 2491;
		int debug_pos2 = 2785;
	#endif

	for (k=0; k<numCols-WIN_SIZE+1; k++) {
	// for (k=2491; k<=2785; k++) {
		// get the subsequences for the window
		fr_pos = k;
		#ifdef SHOW_DEBUG_MESG
			fr_pos = 2755;
		#endif
		to_pos = fr_pos+WIN_SIZE-1;
		if (k > 0 && k % 1000 == 0) {
			cerr << k << "/" << numCols-WIN_SIZE << endl << flush;
		}
		tot = (int) (rds.getSubSeqs(fr_pos, to_pos, MIN_FREQ_FOR_SUBSEQ, MIN_FREQ_RATIO_FOR_SUBSEQ, subSeq, subSeqFreq) + 0.5);
		/*
		if (subSeq.size()==2 || subSeq.size()==3) {
			cout << k << "\t" << tot;
			for (i=0; i<subSeq.size(); i++) {
				cout << "\t" << subSeqFreq[i];
			}
			cout << endl;
		}*/

		#ifdef SHOW_DEBUG_MESG
		cout << "tot = " << tot << endl << flush;
		for (i=0; i<subSeq.size(); i++) {
			cout << subSeq[i] << " {" << subSeqFreq[i] << "} " << (double) subSeqFreq[i] / tot << endl;
		}
		if (tot < MIN_COVER_FOR_WIN)
			exit(1);
		#endif
	
		if (tot < MIN_COVER_FOR_WIN)
			continue;
		
		#ifdef SHOW_DEBUG_MESG
			cout << endl;
			cout << "[" << fr_pos << "," << to_pos << "]" << endl;
			cout << "----------------" << endl;
		#endif
		
		numPatterns = subSeq.size();
		bestCombin.resize(numCombin, numPatterns, numHap);
		
		// set up the logL array
		for (i=0; i<numPatterns; i++) {
			m = (int) (subSeqFreq[i]+0.5);
			n = tot - m;
			logC = logCombin(tot, m);
			for (j=0; j<numCombin; j++) {
				logL = logC + m * slogF[j] + n * slogFF[j];
				bestCombin.logL[i*numCombin + j] = logL;
			}
		}
		
	#ifdef SHOW_DEBUG_MESG
		// show the log-likelihood values
		for (i=0; i<numPatterns; i++) {
			for (j=0; j<numCombin; j++) {
				cout << bestCombin.logL[i*numCombin + j] << "\t";
			}
			cout << endl;
		}
	#endif

		// compute the best combination
		bestCombin.compute(items, patterns, maxLogL, max2LogL);
		#ifdef USE_WEIGHT_APPROACH
		if (maxLogL < MIN_SUBSEQ_LOGL_THRES) {
		#else
		if (maxLogL < MIN_SUBSEQ_LOGL_THRES || maxLogL - max2LogL < DELTA_SUBSEQ_LOGL) {
		#endif
			#ifdef SHOW_DEBUG2_MESG
				if (fr_pos >= debug_pos-WIN_SIZE && fr_pos <=debug_pos2) {
					cerr << fr_pos << " ";
					cerr << string(fr_pos - (debug_pos-WIN_SIZE), ' ');
					if (items[0]==0)
						cerr << subSeq[patterns[0]];
					else if (items[1]==0)
						cerr << subSeq[patterns[1]];
					else
						cerr << subSeq[patterns[2]];
					cerr << " #discard*";
					cerr << " diff=" << maxLogL - max2LogL;
					cerr << " tot=" << tot;
					cerr << " maxLogL=" << maxLogL;
					cerr << endl;
				}
			#endif

			#ifdef SHOW_DEBUG_MESG
				cout << "The difference between maxLogL and max2LogL is too small: " << maxLogL - max2LogL << " tot=" << tot << endl;
				exit(1);
			#else
				continue;
			#endif
		}
		
		for (j=0; j<numHap; j++) {
			hapSubSeqs[items[j]].push_back(subSeq[patterns[j]]);
			hapSubSeqPos[items[j]].push_back(fr_pos);
		}
		
		#ifdef USE_WEIGHT_APPROACH
			hapSubSeqWeight.push_back(maxLogL - max2LogL);
		#else
			hapSubSeqWeight.push_back(1.0);
		#endif
		
		#ifdef SHOW_DEBUG_MESG
			cout << "Resulting logL = " << maxLogL << endl;
			cout << "Diff logL = " << maxLogL - max2LogL << endl;
			for (j=0; j<numHap; j++) {
				cout << "pos: " << hapSubSeqPos[j].at(hapSubSeqPos[j].size()-1);
				cout << " H" << j+1 << ": " << hapSubSeqs[j].at(hapSubSeqs[j].size()-1);
				cout << " tot=" << tot << endl;
			}
		#endif
		#ifdef SHOW_DEBUG_MESG
			exit(1);
		#endif
		#ifdef SHOW_DEBUG2_MESG
		if (fr_pos >= debug_pos-WIN_SIZE && fr_pos <=debug_pos2) {
			cerr << fr_pos << " ";
			cerr << string(fr_pos - (debug_pos-WIN_SIZE), ' ');
			cerr << hapSubSeqs[0].at(hapSubSeqs[0].size()-1);
			cerr << " tot=" << tot;
			cerr << " diff=" << maxLogL - max2LogL;
			cerr << " maxLogL=" << maxLogL << endl;
		}
		#endif
	}
	#ifdef SHOW_DEBUG2_MESG
		exit(1);
	#endif
}

// update the haplotype sequences according to the set of subsequences
void ConstructHaplo::showHapSubSeqs() {
	int i,j;
	for (i=0; i<numHap; i++) {
		cout << "----------------------------------------" << endl;
		cout << "Subseqs for haplo " << i+1 << ":" << endl;
		for (j=0; j<hapSubSeqs[i].size(); j++) {
			cout << hapSubSeqPos[i].at(j) << "\t" << hapSubSeqs[i].at(j) << endl;
		}
		cout << "----------------------------------------" << endl;
	}
}

// update the haplotype sequences according to the set of subsequences
void ConstructHaplo::updateHaploBySubSeqs() {
	int i,j,k;
	int p;
	string s;
	double weight;
	double* stat = new double[numCols * 5];
	double* totCov = new double[numCols];
	int max_k;
	double max_freq;
	char int2Nucl[] = {'-','A','C','G','T'};
	for (i=0; i<numHap; i++) {
		// get the statistics for the haplotype
		memset(stat, 0, numCols*5*sizeof(double));
		memset(totCov, 0, numCols*sizeof(double));
		for (j=0; j<hapSubSeqPos[i].size(); j++) {
			p = hapSubSeqPos[i].at(j);
			s = hapSubSeqs[i].at(j);
			weight = hapSubSeqWeight[j];
			for (k=0; k<s.length(); k++) {
				// stat[(p+k)*5 + nucl_map[s[k]]]++;
				// totCov[p+k]++;
				stat[(p+k)*5 + nucl_map[s[k]]]+=weight;
				totCov[p+k]+=weight;
			}
		}
		// update the haplotype
		for (j=0; j<numCols; j++) {
			if (totCov[j] > 0) {
				max_k=0;
				max_freq=stat[j*5];
				for (k=1; k<5; k++) {
					if (stat[j*5+k] > max_freq) {
						max_freq = stat[j*5+k];
						max_k = k;
					}
				}
				haplos[i].at(j) = int2Nucl[max_k];
			}
		}
	}
	delete[] stat;
	delete[] totCov;
}

#ifndef NOMAIN

int main(int argc, char** argv) {

	cerr << "Welcome to ConstructHaplo ";
	cerr << "version " << VERSION;

#ifdef STRICT
	cerr << " - STRICT";
#endif

	cerr << endl;

	if (argc < 3) {
		cout << "Syntax: " << argv[0] << " [sam/bam file] [frequency file]  >  [estimated haplotype fasta file]" << endl;
		exit(1);
	}
	
	int c, i, j;
	double avgCoverage;

	ConstructHaplo read_likes;
	char* readFile = argv[1];
	char* freqFile = argv[2];
	
	// get the best window pair
	int bestDist;
	double bestCoverage;


	// show the window size
	cerr << "Window size = " << WIN_SIZE << endl;
	cerr << "Window-pair size = " << WIN_PAIR_LEN << endl;

	// load the read and frequency files
	read_likes.loadReadFile(readFile);
	
	read_likes.loadFreqFile(freqFile);

#ifdef DO_BASIC_ERROR_CORR
	// correct the reads
	read_likes.readCorrection();
#endif

#ifdef DO_KMER_ERROR_CORR
	// correct the reads by using relative kmer approach
	read_likes.rds.errorCorrection(WIN_SIZE_ERR_CORR, MIN_COVER_FOR_WIN, MIN_FREQ_FOR_SUBSEQ, MIN_FREQ_RATIO_FOR_SUBSEQ);
#endif

	// combine the overlap reads of the same pair into a long read
	// and update the frequencies
	read_likes.combineRdAndUpdateFreq();

	// get the average of the coverage
	avgCoverage = read_likes.partInfo.getAvgCoverage();
	cerr << "Average coverage = " << avgCoverage << endl << flush;
	
	bestCoverage = read_likes.rds.computeBestDivWindowPair(read_likes.numHap, MIN_COVER_FOR_WIN_CHECKING,
		MIN_FREQ_FOR_SUBSEQ, MIN_FREQ_RATIO_FOR_SUBSEQ, bestDist);
	
	cerr << "computing the best distance between the window pair..." << endl << flush;
	cerr << "Best distance between the window pair is: " << bestDist << endl << flush;

	// show the distance between pair-end reads
	// read_likes.rds.computeAvgReadPosDist();
	// exit(1);
	
	// bestCoverage = read_likes.rds.computeBestWindowPair(bestDist);

	// check the best pair-end window with the highest diversity
	// report the average coverage
	/*
	int minCheckingWinCover = (int) (0.1 * avgCoverage);
	if (minCheckingWinCover < MIN_COVER_FOR_WIN_CHECKING)
		minCheckingWinCover = MIN_COVER_FOR_WIN_CHECKING;*/

	// get the underlining frequencies of the haplotypes (using window)
	vector<double> hapFreqs;
	int winNum;
	// read_likes.getHaploFreq(WIN_SIZE, hapFreqs, winNum);
	
	// get the underlining frequencies of the haplotypes (using window-pairs)
	read_likes.getHaploFreq(bestDist, WIN_PAIR_LEN, hapFreqs, winNum);

	exit(1);
	
	// sort the reads
	read_likes.rds.sortReads();
/*
	// print out the reads
	int fr_pos = atoi(argv[3]);
	int to_pos = atoi(argv[4]);
	read_likes.rds.printReadGrps(fr_pos, to_pos);
	exit(1);
*/
	// calculate the haplo ratio for each combinition
	cerr << "calculate the haplo ratio for each combinition" << endl << flush;
	read_likes.calHapRateCombins();
	
	// calculate P(ci | Hj) for all bases of all columns and combinations
	cerr << "calculate P(ci | Hj) for all bases of all columns and combinations" << endl << flush;
	read_likes.calBaseCombinL();

	// construct the haplotypes by using site-by-site approach
	cerr << "construct the haplotypes" << endl << flush;
	read_likes.constructHaplos();
	
	/*
	// remove the gaps inside the haplotype (first version)
	cerr << "remove the gaps inside the haplotype (first version)" << endl << flush;
	read_likes.rmGapInsideHaplo();

	// list out the version 1 of haplotypes
	cerr << "list out the version 1 of haplotypes" << endl;
	read_likes.listHaplos();
	exit(1);
	*/

	//---------------------------------------------------------
	// using window approach
	//---------------------------------------------------------
	
	// calculate the ratio for each combinations of subsequences
	cerr << "calculate the ratio for each combinations of subsequences" << endl << flush;
	read_likes.calRateForSubSeqCombins();
	
	// get the set of subsequences belonging to each haplotype
	cerr << "get the set of subsequences belonging to each haplotype" << endl << flush;
	read_likes.getSubSeqForHaplos();

	// cerr << "show the subsequences for each haplotype" << endl;
	// read_likes.showHapSubSeqs();

	// update the haplotype sequences according to the set of subsequences
	cerr << "update the haplotype sequences according to the set of subsequences" << endl << flush;
	read_likes.updateHaploBySubSeqs();	

	//---------------------------------------------------------
	// using window-pair approach
	//---------------------------------------------------------
	// 
	
	// get the set of pairs of subsequences belonging to each haplotype
	cerr << "get the set of pairs of subsequences belonging to each haplotype" << endl << flush;
	read_likes.getSubSeqPairForHaplos(bestDist, WIN_PAIR_LEN);


	// update the haplotype sequences according to the set of subsequences
	cerr << "update the haplotype sequences according to the set of pairs of subsequences" << endl << flush;
	read_likes.updateHaploBySubSeqs();	




/*
#ifndef NO_REMOVE_GAP
	// remove the gaps inside the haplotype (first version)
	cerr << "remove the gaps inside the haplotype (first version)" << endl << flush;
	read_likes.rmGapInsideHaplo();
#endif
*/

#ifdef OUT_FIRST_PRE_SEQ
	// list out the first version of haplotypes
	cerr << "list out the first version of haplotypes" << endl;
	read_likes.listHaplos();
	exit(1);
#endif

	// get the reads beloging to each haplotype
	cerr << "get the reads beloging to each haplotype" << endl << flush;
	read_likes.getHapReads();

/*
	// print out the read
	int fr_pos = atoi(argv[3]);
	int to_pos = atoi(argv[4]);
	int which_haplo = atoi(argv[5]);
	// print out the reads indicated inside the vector
	read_likes.rds.printReads(fr_pos, to_pos, read_likes.hapReads[which_haplo]);
	exit(1);
*/

	// get the hits of the read to each estimated haplotype
	// cerr << "get the hits of the read to each estimated haplotype" << endl << flush;
	read_likes.getHapReadHits();

	// get the full haplotypes based on hapReadHits
	// cerr << "get the full haplotypes based on hapReadHits" << endl << flush;
	read_likes.getFullHaplos();
	
	/*
	// list out the region of the haplotypes
	for (i=0; i<read_likes.numHap; i++)
		cout << read_likes.fullHaplos[i].substr(1472, 50) << endl;
	*/

#ifndef NO_REMOVE_GAP
	#ifndef NO_REMOVE_LOW_COVERAGE
		// remove both gaps inside the full haplotype
		// and remove the too-low-coverage region around the ends of the haplotype
		read_likes.rmGapAndLowCoverInsideFullHaplo();
	#else
		// remove the gaps inside the full haplotype
		read_likes.rmGapInsideFullHaplo();
	#endif
#endif

	// list out the full haplotypes
	// cerr << "list out the full haplotypes" << endl << flush;
	read_likes.listFullHaplos();
}

#endif