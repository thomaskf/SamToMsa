#include "reads.h"

bool cmpInsert(Insert* i1, Insert* i2) {
	if (i1->col == i2->col) {
		if (i1->len == i2->len) {
			return i1->readid < i2->readid;
		} else {
			return i1->len > i2->len;
		}
	}
	return i1->col > i2->col;
}

Insert::Insert(int c, int p, string s, int id) {
	col = c;
	readPos = p;
	x = s;
	readid = id;
	len = s.length();
}

void Insert::show() {
	cout << col << "\t" << readPos << "\t" << x << "\t" << readid << endl;
}

// constructor
/*
Read::Read(string oSeq, string s, int mPos, string mHap, string baseDel) {
	origSeq = oSeq;
	seq = s;
	mapPos = mPos;
	mapHap = mHap;
	mapLen = seq.length();
	seqLen = oSeq.length();
	baseDeleted = baseDel;
}
*/

// constructor
Read::Read(string oSeq, string s, int mPos, string mHap, string des) {
	origSeq = oSeq;
	seq = s;
	mapPos = mPos;
	mapHap = mHap;
	desc = des;
	// mapLen = seq.length();
	seqLen = oSeq.length();
	
	isDiscard = false;
}

// combine this read with another read
// return true if successful
bool Read::combineRead(Read* aRead, int& newMPos, string& newSeq) {
	int lastMapPos;
	int overlapPosOnFirstRead;
	int i;
	int overlaplen;
	if (mapPos <= aRead->mapPos) {
		// this read is before aRead
		lastMapPos = mapPos + seq.length() - 1;
		overlaplen = lastMapPos - aRead->mapPos + 1;
		if (overlaplen > MIN_OVERLAP_COMBINE_READ) {
			// check whether the overlapped region is consistent
			overlapPosOnFirstRead = seq.length() - overlaplen;
			for (i=0; i<overlaplen; i++) {
				if (aRead->seq[i] != seq[overlapPosOnFirstRead+i])
					return false;
			}
			newMPos = mapPos;
			newSeq = seq;
			newSeq.append(aRead->seq.substr(overlaplen));
		} else
			return false;
	} else
		return aRead->combineRead(this, newMPos, newSeq);
	return true;
}

/*
// get the number of deleted bases in the beginning of the read
int Read::numDelBegin() {
	int i=0;
	while (i<baseDeleted.length() && baseDeleted[i]=='1')
		i++;
	return i;
}*/

// compute the avg distance between the starting positions of the pair-ended reads
void Reads::computeAvgReadPosDist() {
	vector<int> dist;
	Read* rd1;
	Read* rd2;
	int d,i;
	for (i=0; i<pairRead1.size(); i++) {
		rd1 = pairRead1[i];
		rd2 = pairRead2[i];
		d = abs(rd2->mapPos - rd1->mapPos);
		dist.push_back(d);
	}
	// show the array "dist"
	for (i=0; i<dist.size(); i++) {
		cout << dist[i] << endl;
	}
	exit(1);
}

int Reads::readNumCoverWinPair(int winStart1, int winEnd1, int winStart2, int winEnd2) {
	int num = 0;
	int posStart1, posStart2;
	int posEnd1, posEnd2;
	bool rd1CoverWin1, rd1CoverWin2;
	bool rd2CoverWin1, rd2CoverWin2;
	int i;
	// check pair-end reads
	for (i=0; i<pairRead1.size(); i++) {
		posStart1 = pairRead1[i]->mapPos;
		posStart2 = pairRead2[i]->mapPos;
		posEnd1 = posStart1 + pairRead1[i]->seq.length() - 1;
		posEnd2 = posStart2 + pairRead2[i]->seq.length() - 1;
		rd1CoverWin1 = (posStart1 <= winStart1 && posEnd1 >= winEnd1);
		rd2CoverWin1 = (posStart2 <= winStart1 && posEnd2 >= winEnd1);
		rd1CoverWin2 = (posStart1 <= winStart2 && posEnd1 >= winEnd2);
		rd2CoverWin2 = (posStart2 <= winStart2 && posEnd2 >= winEnd2);
		if ((rd1CoverWin1 || rd2CoverWin1) && (rd1CoverWin2 || rd2CoverWin2))
			num++;
	}
	// check single-end reads
	for (i=0; i<unpairReads.size(); i++) {
		posStart1 = unpairReads[i]->mapPos;
		posEnd1 = posStart1 + unpairReads[i]->seq.length() - 1;
		rd1CoverWin1 = (posStart1 <= winStart1 && posEnd1 >= winEnd1);
		rd1CoverWin2 = (posStart1 <= winStart2 && posEnd1 >= winEnd2);
		if (rd1CoverWin1 && rd1CoverWin2)
			num++;
	}
	return num;
}

// check the best pair-end window
// report the average coverage
double Reads::computeBestWindowPair(int& bestDist) {
	int d,p,toP,numP;
	int i;
	int winStart1, winStart2;
	int winEnd1, winEnd2;
	double currCoverage=0.0;
	double bestCoverage=-1.0;
	
	int maxCol = maxColIndex();
	
	// first, compute the best window distance
	// d = MIN_WIN_PAIR_DIST - DIST_INTERVAL;
	for (d=MIN_WIN_PAIR_DIST; d<=MAX_WIN_PAIR_DIST; d+= DIST_INTERVAL) {
	// while (currCoverage >= bestCoverage) {
	// while (d <= 400) {
		// d += DIST_INTERVAL;
		currCoverage = 0.0;
		toP = maxCol - d - 2 * WIN_PAIR_LEN;
		numP = 0;
		for (p=0; p<toP; p+=POS_INTERVAL) {
			winStart1 = p;
			winEnd1 = p + WIN_PAIR_LEN - 1;
			winStart2 = winEnd1 + d + 1;
			winEnd2 = winStart2 + WIN_PAIR_LEN - 1;
			currCoverage += readNumCoverWinPair(winStart1, winEnd1, winStart2, winEnd2);
			numP++;
		}
		currCoverage = currCoverage / numP;
		if (currCoverage > 0.0 && currCoverage > bestCoverage) {
			bestCoverage = currCoverage;
			bestDist = d;
		}
		// cout << d << "\t" << currCoverage << endl;
	}
	
	// second, compute the best window length
	cout << "bestDist = " << bestDist << endl;
	cout << "bestCoverage = " << bestCoverage << endl;
	
	return bestCoverage;
}

// check the best pair-end window with the highest diversity
// report the average coverage
double Reads::computeBestDivWindowPair(int numHaps, int minWinCover, int minFreq,
		double minFreqRatio, int& bestDist) {
			
	vector<int> d_candi;
	vector<int> d_valids;
	int i,d;
	int pos_interval;
	double coverage;
	for (d=MIN_WIN_PAIR_DIST; d<=MAX_WIN_PAIR_DIST; d+= DIST_INTERVAL) {
		d_candi.push_back(d);
	}
	// first test
	cerr << "step 1 / 3" << endl << flush;
	pos_interval = 50;
	coverage = computeBestDivWindowPair(numHaps, minWinCover, minFreq, minFreqRatio, bestDist, 
		pos_interval, d_candi, d_valids);
	if (d_valids.size() == 1)
		return coverage;
		/*
	cout << "after " << pos_interval << endl;
	for (i=0; i<d_valids.size(); i++)
		cout << d_valids[i] << " ";
	cout << endl;
	if (d_valids.size() == 1)
		return coverage;
		 */ 
	// second test
	cerr << "step 2 / 3" << endl << flush;
	pos_interval = 10;
	d_candi.clear();
	for (i=0; i<d_valids.size(); i++)
		d_candi.push_back(d_valids[i]);
	coverage = computeBestDivWindowPair(numHaps, minWinCover, minFreq, minFreqRatio, bestDist, 
		pos_interval, d_candi, d_valids);
	/*
	cout << "after " << pos_interval << endl;
	for (i=0; i<d_valids.size(); i++)
		cout << d_valids[i] << " ";
	cout << endl;
	 */
	if (d_valids.size() == 1)
		return coverage;
	// last test
	cerr << "step 3 / 3" << endl << flush;
	pos_interval = 5;
	d_candi.clear();
	for (i=0; i<d_valids.size(); i++)
		d_candi.push_back(d_valids[i]);
	return  computeBestDivWindowPair(numHaps, minWinCover, minFreq, minFreqRatio, bestDist, 
		pos_interval, d_candi, d_valids);
}

// check the best pair-end window with the highest diversity
// report the average coverage
double Reads::computeBestDivWindowPair(int numHaps, int minWinCover, int minFreq,
		double minFreqRatio, int& bestDist, int pos_interval, vector<int> d_candi,
		vector<int>& valid_candi) {
	int d,p,toP,numP;
	int i;
	int winStart1, winStart2;
	int winEnd1, winEnd2;
	vector<string> subSeq1;
	vector<string> subSeq2;
	vector<double> subSeqFreq;
	double currCoverage;
	double avgCoverage;

	int maxP = 0;
	double theCoverage = 0.0;
	
	int maxCol = maxColIndex();
	valid_candi.clear();
	
	// first, compute the best window distance
	// d = MIN_WIN_PAIR_DIST - DIST_INTERVAL;
	int x = d_candi.size() / 10;
	for (i=0; i<d_candi.size(); i++) {
		/*
		if (d_candi.size() > 10) {
			if (i>0 && i%x==0) {
				cout << i << "/" << d_candi.size() << endl;
			}
		} else {
			cout << i << "/" << d_candi.size() << endl;
		}*/
		d = d_candi[i];
		avgCoverage = 0.0;
		toP = maxCol - d - 2 * WIN_PAIR_LEN + 1;
		numP = 0;
		for (p=0; p<=toP; p+=pos_interval) {
			winStart1 = p;
			winEnd1 = p + WIN_PAIR_LEN - 1;
			winStart2 = winEnd1 + d + 1;
			winEnd2 = winStart2 + WIN_PAIR_LEN - 1;
			currCoverage = getSubSeqPairs(winStart1, winEnd1, winStart2, winEnd2, minFreq,
				minFreqRatio, subSeq1, subSeq2, subSeqFreq);
			// cout << d << "\t" << p << "\t" << subSeq1.size() << endl;
			if (subSeq1.size() == numHaps) {
				numP++;
				avgCoverage+=currCoverage;
			}
		}
		avgCoverage = avgCoverage / numP;
		cout << d << "\t" << numP << "\t" << avgCoverage << endl << flush;
		if (avgCoverage > (double) minWinCover && numP >= maxP) {
			if (numP == maxP) {
				valid_candi.push_back(d);
			} else {
				valid_candi.clear();
				valid_candi.push_back(d);
			}
			if (numP > maxP || avgCoverage > theCoverage) {
				theCoverage = avgCoverage;
				maxP = numP;
				bestDist = d;
			}
		}
		// cout << d << "\t" << currCoverage << endl;
	}
	
	cout << "bestDist = " << bestDist << endl;
	// cout << "theCoverage = " << theCoverage << endl;
	// cout << "maxP = " << maxP << endl;
	
	return theCoverage;
}

// update according to the insertions
void Reads::updateForInserts() {
	
	// sort the insertions
	sort(inserts.begin(), inserts.end(), cmpInsert);
	
	cerr << "Number of insertions: " << inserts.size() << endl;
	
	// update the positions and sequences
	int i,j;
	int cpos,rpos;
	string s;
	int rid;
	int l;
	int pre_cpos = -1;
	int rlen;
	Read* curRd;
	for (i=0; i<inserts.size(); i++) {
		cpos = inserts[i]->col;
		rpos = inserts[i]->readPos;
		s = inserts[i]->x;
		rid = inserts[i]->readid;
		l = inserts[i]->len;
		rlen = readlist[rid]->seq.length();
		if (cpos != pre_cpos) {
			// create l columns
			for (j=0; j<readlist.size(); j++) {
				curRd = readlist[j];
				if (curRd->mapPos >= cpos) {
					curRd->mapPos += l;
				} else if (curRd->mapPos < cpos && curRd->seq.length()+curRd->mapPos > cpos) {
					curRd->seq.insert(cpos - curRd->mapPos, l, '-');
				}
			}
			pre_cpos = cpos;
		}
		// put the characters into the corresponding positions of the read
		for (j=0; j<l; j++) {
			if (rpos+j < rlen)
				readlist[rid]->seq.at(rpos+j) = s[j];
			else
				readlist[rid]->seq.append(1, s[j]);
		}
	}
	
}

// update the seq according to the cigar string
// and store the insertions if exist
// return false if the cigarStr is not valid
bool Reads::updateSeq(string& seq, string& cigarStr, int readId, int mapPos) {
	int i,j,k;
	int seqI;
	string s;
	
	i=0;
	s="";
	seqI = 0;
	for (j=0; j<cigarStr.length(); j++) {
		// should begin with a number, and followed by a char
		if (isdigit(cigarStr[j]))
			continue;
		else {
			if (j==i) {
				// cout << "Warning! Invalid cigar string: " << cigarStr << endl;
				return false;
			}
			k = atoi(cigarStr.substr(i,j-i).c_str());
			switch (cigarStr[j]) {
				case 'I':
					inserts.push_back(new Insert(mapPos+seqI,seqI,seq.substr(seqI,k),readId));
				case 'S':
					seq.erase(seqI,k);
					break;
				case 'D':
				case 'N':
					seq.insert(seqI,k,'-');
					seqI += k;
					break;
				case 'M':
					seqI += k;
					break;
				case 'P':
				case 'H':
					// do nothing
					break;
			}
			i=j+1;
		}
	}
	return true;
}

/*
// update the seq according to the cigar string
bool updateSeq(string& seq, string& cigarStr, string& baseDel) {
	int i,j,k;
	int seqI;
	string s;
	
	i=0;
	s=seq;
	seqI = 0;
	baseDel = "";
	for (j=0; j<cigarStr.length(); j++) {
		// should begin with a number, and followed by a char
		if (isdigit(cigarStr[j]))
			continue;
		else {
			if (j==i) {
				// cout << "Warning! Invalid cigar string: " << cigarStr << endl;
				return false;
			}
			k = atoi(cigarStr.substr(i,j-i).c_str());
			switch (cigarStr[j]) {
				case 'I':
				case 'S':
					seq.erase(seqI,k);
					baseDel.append(k, (char)'1');
					break;
				case 'D':
				case 'N':
					seq.insert(seqI,k,'-');
					seqI += k;
					break;
				case 'M':
					seqI += k;
					baseDel.append(k, (char)'0');
					break;
				case 'P':
				case 'H':
					// do nothing
					break;
			}
			i=j+1;
		}
	}
	return true;
}
*/

// vector<Read*> pairRead1; // paired-end reads
// vector<Read*> pairRead2;
// vector<Read*> unpairReads; // single-end reads
// load reads from the SAM/BAM file
void Reads::readSamFile(char* samFile) {
	string aline;
	string mhap;
	int mpos;
	int flag;
	string origSeq,seq;
	string cigar;
	vector<string> token;
	SamBamFileHander fileHander;
	string preSeqName = "";
	string currSeqName;
	Read* preRead = NULL;
	Read* currRead;
	int indx = 0;
	// clear the arrays
	pairRead1.clear();
	pairRead2.clear();
	unpairReads.clear();
	readlist.clear();
	// open the file
	fileHander.openFile(samFile);
	while (fileHander.getNextSeq(aline)) {
		if (aline.length() > 0 && aline[0]!='@') {
			tokenizer(aline, "\t", &token);
			if (token.size() > 9) {
				currSeqName = token[0];
				flag = atoi(token[1].c_str());
				if ((flag&2048)!=0) {
					// this is a supplementry alignment
					// discard it
					continue;
				}
				mhap = token[2];
				mpos = atoi(token[3].c_str())-1;
				cigar = token[5];
				seq = token[9];
				origSeq = seq;
				if (updateSeq(seq, cigar, indx, mpos)) {
					// alignment is valid
					currRead = new Read(origSeq, seq, mpos, mhap, currSeqName);
					readlist.push_back(currRead);
					indx++;
					if (preRead!=NULL) {
						if (currSeqName == preSeqName) {
							// pair-ended reads
							if (preRead->mapPos <= currRead->mapPos) {
								pairRead1.push_back(preRead);
								pairRead2.push_back(currRead);
							} else {
								pairRead1.push_back(currRead);
								pairRead2.push_back(preRead);
							}
							preRead = NULL;
							preSeqName = "";
						} else {
							// single-end read
							unpairReads.push_back(preRead);
							preRead = currRead;
							preSeqName = currSeqName;
						}
					} else {
						// first read (of the pair)
						preRead = currRead;
						preSeqName = currSeqName;
					}
				}
			}
		}
	}
	if (preRead!=NULL) {
		// single-end read
		unpairReads.push_back(preRead);
	}
	readNum = indx;
	fileHander.closeFile();
}

// combine the overlapped read to a longer read
// should run this program after basic error correction
void Reads::combineReads() {
	int i;
	int k=0;
	int numCombineReads=0;
	int newMPos;
	string newSeq;
	Read* currRead;
	for (i=0; i<pairRead1.size(); i++) {
		if (pairRead1[i]->combineRead(pairRead2[i], newMPos, newSeq)) {
			// two reads overlap each other consistently
			currRead = new Read(newSeq, newSeq, newMPos, pairRead1[i]->mapHap, pairRead1[i]->desc);
			unpairReads.push_back(currRead);
			numCombineReads++;
		} else {
			if (i > k) {
				pairRead1[k] = pairRead1[i];
				pairRead2[k] = pairRead2[i];
			}
			k++;
		}
	}
	if (k < pairRead1.size()) {
		pairRead1.resize(k);
		pairRead2.resize(k);
	}
	cerr << "Number of pairs with both ends combining into a long read: " << numCombineReads << endl;
}

// sort the reads
bool compRds(Read* rd1, Read* rd2) {
	if (rd1->mapPos == rd2->mapPos)
		return rd1->seq.length() < rd2->seq.length();
	else
		return rd1->mapPos < rd2->mapPos;
}
void Reads::sortReads() {
	sort(readlist.begin(), readlist.end(), compRds);
}

// return the max column index
int Reads::maxColIndex() {
	int i,k;
	int m = 0;
	for (i=0; i<readNum; i++) {
		k = readlist[i]->mapPos + readlist[i]->seq.length() - 1;
		if (m < k)
			m = k;
	}
	return m;
}

/*
// print out the reads
void Reads::printReads(int fr_pos, int to_pos) {
	int i;
	int mpos;
	for (i=0; i<readNum; i++) {
		mpos =  readlist[i]->mapPos;
		if ( mpos >= fr_pos && mpos <= to_pos) {
			if ( mpos > fr_pos)
				cout << string(mpos-fr_pos,' ');
			cout << readlist[i]->seq << endl;
		}
		if (mpos > to_pos)
			break;
	}
}
*/

// print out the reads
void Reads::printReads(int fr_pos, int to_pos) {
	int i;
	int mpos;
	int len;
	for (i=0; i<readNum; i++) {
		mpos =  readlist[i]->mapPos;
		len = readlist[i]->seq.length();
		if ( mpos+len > to_pos && mpos <= fr_pos) {
			cout << readlist[i]->seq.substr(fr_pos - mpos, to_pos-fr_pos+1);
			// cout << " [" << i << "]";
			cout << endl;
		}
		if (mpos > to_pos)
			break;
	}
}

// print out the read groups
void Reads::printReadGrps(int fr_pos, int to_pos) {
	int i;
	int mpos;
	int len;
	string s;
	map<string, int> readGrp;
	map<string, int>::iterator itr;
	for (i=0; i<readNum; i++) {
		mpos =  readlist[i]->mapPos;
		len = readlist[i]->seq.length();
		if ( mpos+len > to_pos && mpos <= fr_pos) {
			s = readlist[i]->seq.substr(fr_pos - mpos, to_pos-fr_pos+1);
			itr = readGrp.find(s);
			if (itr == readGrp.end()) {
				readGrp.insert(pair<string,int>(s,1));
			} else {
				itr->second++;
			}
		}
		if (mpos > to_pos)
			break;
	}
	for (itr=readGrp.begin(); itr!=readGrp.end(); itr++) {
		cout << itr->first << "\t" << itr->second << endl;
	}
}

// get the pairs of subsequences for a region (with freq >= minFreq)
// return the total number of pairs of reads covering the region
double Reads::getSubSeqPairs(int fr_pos1, int to_pos1, int fr_pos2, int to_pos2, int minFreq,
		double minFreqRatio, vector<string>& subSeq1, vector<string>& subSeq2, vector<double>& subSeqFreq) {
	
	int i,j,k;
	
	// initialize
	subSeq1.clear();
	subSeq2.clear();
	subSeqFreq.clear();
	
	// for storing the pairs of subsequences and the corresponding frequencies
	map<string, double> subSeqPairs;
	map<string, double>::iterator itr;
	
	// for pair-ended reads
	Read *rd1, *rd2;
	int mpos1, mpos2;
	int len1, len2;
	int regionLen1, regionLen2;
	bool rd1cover1, rd2cover1, rd1cover2, rd2cover2;
	double eachFreq;
	double tot = 0.0;
	string s;
	vector<string> s1; // subsequences for the region 1
	vector<string> s2; // subsequences for the region 2
	double num;
	double ratio_thres;
	regionLen1 = to_pos1 - fr_pos1 + 1;
	regionLen2 = to_pos2 - fr_pos2 + 1;
	for (i=0; i<pairRead1.size(); i++) {
		rd1 = pairRead1[i];
		rd2 = pairRead2[i];
		mpos1 = rd1->mapPos;
		mpos2 = rd2->mapPos;
		len1 = rd1->seq.length();
		len2 = rd2->seq.length();
		rd1cover1 = (mpos1+len1 > to_pos1 && mpos1 <= fr_pos1);
		rd2cover1 = (mpos2+len2 > to_pos1 && mpos2 <= fr_pos1);
		rd1cover2 = (mpos1+len1 > to_pos2 && mpos1 <= fr_pos2);
		rd2cover2 = (mpos2+len2 > to_pos2 && mpos2 <= fr_pos2);
		if ((rd1cover1 || rd2cover1) && (rd1cover2 || rd2cover2)) {
			// this pair of read cover the region
			s1.clear();
			s2.clear();
			if (rd1cover1)
				s1.push_back(rd1->seq.substr(fr_pos1-mpos1, regionLen1));
			if (rd2cover1)
				s1.push_back(rd2->seq.substr(fr_pos1-mpos2, regionLen1));
			if (rd1cover2)
				s2.push_back(rd1->seq.substr(fr_pos2-mpos1, regionLen2));
			if (rd2cover2)
				s2.push_back(rd2->seq.substr(fr_pos2-mpos2, regionLen2));
			eachFreq = 1.0 / (double) (s1.size() * s2.size());
			for (j=0; j<s1.size(); j++) {
				for (k=0; k<s2.size(); k++) {
					s = s1[j] + s2[k];
					itr = subSeqPairs.find(s);
					if (itr != subSeqPairs.end()) {
						itr->second += eachFreq;
					} else {
						subSeqPairs.insert(pair<string,double>(s,eachFreq));
					}
					tot += eachFreq;
				}
			}
		}
	}
	
	// for single-ended reads
	for (i=0; i<unpairReads.size(); i++) {
		rd1 = unpairReads[i];
		mpos1 = rd1->mapPos;
		len1 = rd1->seq.length();
		rd1cover1 = (mpos1+len1 > to_pos1 && mpos1 <= fr_pos1);
		rd1cover2 = (mpos1+len1 > to_pos2 && mpos1 <= fr_pos2);
		if (rd1cover1 && rd1cover2) {
			// this read covers the region
			s = rd1->seq.substr(fr_pos1-mpos1, regionLen1) + rd1->seq.substr(fr_pos2-mpos1, regionLen2);
			itr = subSeqPairs.find(s);
			if (itr != subSeqPairs.end()) {
				itr->second += 1.0;
			} else {
				subSeqPairs.insert(pair<string,double>(s,1.0));
			}
			tot += 1.0;
		}
	}
	
	// get the results
	num = 0.0;
	ratio_thres = (double) tot * minFreqRatio;
	for (itr=subSeqPairs.begin(); itr!=subSeqPairs.end(); itr++) {
		if (itr->second >= minFreq && itr->second >= ratio_thres) {
			subSeq1.push_back(itr->first.substr(0,regionLen1));
			subSeq2.push_back(itr->first.substr(regionLen1));
			subSeqFreq.push_back(itr->second);
			num += itr->second;
		}
	}
	return num;
}

// get all the subsequences for a region
// return the total number of fragments covering the region
int Reads::getAllSubSeqs(int fr_pos, int to_pos, map<string,double>& subSeqMap, vector<Read*>* readArr) {
	int i;
	int tot = 0;
	map<string,double>::iterator itr;
	
	// clear the arrays
	subSeqMap.clear();
	if (readArr!=NULL)
		readArr->clear();

	// for pair-ended reads
	Read *rd1, *rd2;
	int mpos1, mpos2;
	int len1, len2;
	bool cover1, cover2;
	string s1, s2;
	for (i=0; i<pairRead1.size(); i++) {
		rd1 = pairRead1[i];
		rd2 = pairRead2[i];
		mpos1 = rd1->mapPos;
		mpos2 = rd2->mapPos;
		len1 = rd1->seq.length();
		len2 = rd2->seq.length();
		cover1 = (mpos1+len1 > to_pos && mpos1 <= fr_pos);
		cover2 = (mpos2+len2 > to_pos && mpos2 <= fr_pos);
		if (cover1 && cover2) {
			// since both ends cover this region, count 0.5 for each of them
			if (readArr != NULL) {
				readArr->push_back(rd1);
				readArr->push_back(rd2);
			}
			s1 = rd1->seq.substr(fr_pos-mpos1, to_pos-fr_pos+1);
			s2 = rd2->seq.substr(fr_pos-mpos2, to_pos-fr_pos+1);
			itr = subSeqMap.find(s1);
			if (itr == subSeqMap.end())
				subSeqMap.insert(pair<string,double>(s1,0.5));
			else
				itr->second+=0.5;
			itr = subSeqMap.find(s2);
			if (itr == subSeqMap.end())
				subSeqMap.insert(pair<string,double>(s2,0.5));
			else
				itr->second+=0.5;
			tot++;
		} else if (cover1) {
			// only the first end cover this region
			if (readArr != NULL)
				readArr->push_back(rd1);
			s1 = rd1->seq.substr(fr_pos-mpos1, to_pos-fr_pos+1);
			itr = subSeqMap.find(s1);
			if (itr == subSeqMap.end())
				subSeqMap.insert(pair<string,double>(s1,1.0));
			else
				itr->second+=1.0;
			tot++;
		} else if (cover2) {
			// only the second end cover this region
			if (readArr != NULL)
				readArr->push_back(rd2);
			s2 = rd2->seq.substr(fr_pos-mpos2, to_pos-fr_pos+1);
			itr = subSeqMap.find(s2);
			if (itr == subSeqMap.end())
				subSeqMap.insert(pair<string,double>(s2,1.0));
			else
				itr->second+=1.0;
			tot++;
		}
	}
	
	// for single-end reads
	Read *rd;
	int mpos;
	int len;
	string s;
	for (i=0; i<unpairReads.size(); i++) {
		rd = unpairReads[i];
		mpos =  rd->mapPos;
		len = rd->seq.length();
		if ( mpos+len > to_pos && mpos <= fr_pos) {
			// it covers this region
			if (readArr != NULL)
				readArr->push_back(rd);
			s = rd->seq.substr(fr_pos - mpos, to_pos-fr_pos+1);
			itr = subSeqMap.find(s);
			if (itr == subSeqMap.end())
				subSeqMap.insert(pair<string,double>(s,1.0));
			else
				itr->second+=1.0;
			tot++;
		}
	}
	return tot;
}


// get the subsequences for a region (with freq >= minFreq and freq >= tot * minFreqRatio)
// return the total number of fragments above the freq threshold covering the region
double Reads::getSubSeqs(int fr_pos, int to_pos, int minFreq, double minFreqRatio,
		vector<string>& subSeq, vector<double>& subSeqFreq) {
	
	map<string,double> subSeqMap;
	map<string,double>::iterator itr;
	
	// clear the arrays
	subSeq.clear();
	subSeqFreq.clear();
	
	int tot = getAllSubSeqs(fr_pos, to_pos, subSeqMap, NULL);
	double num = 0.0;
	double ratio_thres = (int) ((double) tot * minFreqRatio);
	
	for (itr=subSeqMap.begin(); itr!=subSeqMap.end(); itr++) {
		if (itr->second >= minFreq && itr->second >= ratio_thres) {
			subSeq.push_back(itr->first);
			subSeqFreq.push_back(itr->second);
			num += itr->second;
		}
	}
	return num;
}

/*
// print out the reads
void Reads::printAllOrigReads() {
	int i;
	int mpos;
	int len;
	int firstpos = -1;
	for (i=0; i<readNum; i++) {
		mpos =  readlist[i]->mapPos-readlist[i]->numDelBegin();
		if (mpos< 0)
			mpos = 0;
		if (firstpos == -1)
			firstpos = mpos;
		if (mpos > firstpos)
			cout << string(mpos-firstpos, ' ');
		cout << readlist[i]->origSeq;
		cout << " [" << i << ";mpos:" << mpos << "]" << endl;
	}
}
*/

// print out the reads
void Reads::printAllReads() {
	cout << "show the read sequences" << endl;
	int i;
	int mpos;
	int len;
	int firstpos = -1;
	for (i=0; i<readNum; i++) {
		mpos =  readlist[i]->mapPos;
		if (mpos< 0)
			mpos = 0;
		if (firstpos == -1)
			firstpos = mpos;
		if (mpos > firstpos)
			cout << string(mpos-firstpos, ' ');
		cout << readlist[i]->seq;
		cout << " [" << i << ";mpos:" << mpos << "]" << endl;
	}
}

/*
// print out the reads indicated inside the vector
void Reads::printOrigReads(int fr_pos, int to_pos, vector<int>& haploIDs) {
	int i,k;
	int mpos;
	int len;
	int firstpos = -1;
	for (k=0; k<haploIDs.size(); k++) {
		i = haploIDs[k];
		mpos =  readlist[i]->mapPos-readlist[i]->numDelBegin();
		if (mpos< 0)
			mpos = 0;
		len = readlist[i]->origSeq.length();
		if ( mpos+len >= to_pos && mpos <= fr_pos) {
			if (firstpos == -1)
				firstpos = mpos;
			if (mpos > firstpos)
				cout << string(mpos-firstpos,' ');
			cout << readlist[i]->origSeq;
			cout << " [" << i << "]" << endl;
		}
		if (mpos > to_pos)
			break;
	}
}
*/

// print out the reads indicated inside the vector
void Reads::printReads(int fr_pos, int to_pos, vector<int>& haploIDs) {
	int i,k;
	int mpos;
	int len;
	int firstpos = -1;
	for (k=0; k<haploIDs.size(); k++) {
		i = haploIDs[k];
		mpos =  readlist[i]->mapPos;
		len = readlist[i]->seq.length();
		if ( mpos+len > to_pos && mpos <= fr_pos) {
			/*
			if (firstpos == -1)
				firstpos = mpos;
			if (mpos > firstpos)
				cout << string(mpos-firstpos,' ');
			cout << readlist[i]->seq;
			cout << " [" << i << "]" << endl;
			 */
			cout << readlist[i]->seq.substr(fr_pos - mpos, to_pos-fr_pos+1);
			// cout << " [" << i << "]" << endl;
			cout << endl;
		}
		if (mpos > to_pos)
			break;
	}
}

// assign the IDs for the reads
void Reads::assignReadIDs() {
	int i;
	for (i=0; i<readlist.size(); i++) {
		readlist[i]->id = i;
	}
}

// output the reads
void Reads::outReads(char* pairRdFile, char* singleRdFile) {
	int i;
	Read* rd;
	ofstream fout1, fout2;
	
	// for pair-ended read file
	cerr << "Outputing pair-ended reads" << endl;
	fout1.open(pairRdFile);
	for (i=0; i<pairRead1.size(); i++) {
		// read 1
		rd = pairRead1[i];
		fout1 << "Rd1: " << rd->desc << endl;
		fout1 << "pos: " << rd->mapPos << endl;
		fout1 << "seq: " << rd->seq << endl;
		// read 2
		rd = pairRead2[i];
		fout1 << "Rd2: " << rd->desc << endl;
		fout1 << "pos: " << rd->mapPos << endl;
		fout1 << "seq: " << rd->seq << endl;
	}
	fout1.close();
	
	// for single-end read file
	cerr << "Outputing single-ended reads" << endl;
	fout2.open(singleRdFile);
	for (i=0; i<unpairReads.size(); i++) {
		rd = unpairReads[i];
		fout2 << "Rd:  " << rd->desc << endl;
		fout2 << "pos: " << rd->mapPos << endl;
		fout2 << "seq: " << rd->seq << endl;
	}
	fout2.close();
}

int hamming(string& s1, string& s2) {
	// assuming s1.length == s2.length
	int d=0;
	int i;
	for (i=0; i<s1.length(); i++) {
		if (s1[i]!=s2[i])
			d++;
	}
	return d;
}

// error correction by using relative k-mer approach
void Reads::errorCorrection(int winSize, int minCover, int minFreq, double minFreqRatio) {
	int p;
	int tot;
	int i,j,k;
	map<string, double> subSeqMap;
	map<string, double>::iterator itr;
	vector<Read*> readArr;
	int maxColIdx = maxColIndex();
	vector<string> validSubSeqs;
	vector<string> inValidSubSeqs;
	map<string,string> seqCorrMap;
	map<string,string>::iterator itr2;
	int thres;
	int curr_d,best_d;
	vector<string> best_s;
	string invalidS;
	string validS;
	string s;
	cerr << "performing error correction by using relative k-mer approach" << endl;
	for (p=0; p<=maxColIdx-winSize; p++) {
		if (p > 0 && p%500 == 0) {
			cerr << p << "/" << maxColIdx-winSize << endl <<flush;
		}
		tot = getAllSubSeqs(p, p+winSize-1, subSeqMap, &readArr);
		if (tot >= minCover) {

			// for debugging
			/*
			cout << "p=" << p << endl;
			cout << "before correction" << endl;
			for (itr=subSeqMap.begin(); itr!=subSeqMap.end(); itr++) {
				cout << itr->first << "\t" << itr->second << endl;
			}
			cout << "--------------------" << endl;
			*/
			thres = (int) ((double) tot * minFreqRatio);
			if (thres < minFreq)
				thres = minFreq;
			validSubSeqs.clear();
			inValidSubSeqs.clear();
			seqCorrMap.clear();
			for (itr=subSeqMap.begin(); itr!=subSeqMap.end(); itr++) {
				if (itr->second >= thres) {
					validSubSeqs.push_back(itr->first);
				} else {
					inValidSubSeqs.push_back(itr->first);
				}
			}
			// for each invalid subseq, find the valid one most similar to it
			for (i=0; i<inValidSubSeqs.size(); i++) {
				invalidS = inValidSubSeqs[i];
				if (validSubSeqs.size() > 0) {
					validS = validSubSeqs[0];
					best_d = hamming(invalidS, validS);
					best_s.clear();
					best_s.push_back(validS);
					for (j=1; j<validSubSeqs.size(); j++) {
						validS = validSubSeqs[j];
						curr_d = hamming(invalidS, validS);
						if (curr_d < best_d) {
							best_d = curr_d;
							best_s.clear();
							best_s.push_back(validS);
						} else if (curr_d == best_d) {
							best_s.push_back(validS);
						}
					}
					if (best_s.size() == 1) {
						seqCorrMap.insert(pair<string,string>(invalidS,best_s[0]));
					}
				}
			}
			// correct the reads accordingly
			for (i=0; i<readArr.size(); i++) {
				k = p-readArr[i]->mapPos;
				s = readArr[i]->seq.substr(k, winSize);
				itr2 = seqCorrMap.find(s);
				if (itr2!=seqCorrMap.end()) {
					// update the read sequence
					for (j=0; j<winSize; j++) {
						readArr[i]->seq[k+j] = itr2->second.at(j);
					}
				}
			}
			/*
			// for debugging
			cout << "after correction" << endl;
			getAllSubSeqs(p, p+winSize-1, subSeqMap, &readArr);
			for (itr=subSeqMap.begin(); itr!=subSeqMap.end(); itr++) {
				cout << itr->first << "\t" << itr->second << endl;
			}
			cout << "--------------------" << endl;*/
		}
	}
}

// get frequencies of the haplotypes
void Reads::getFreqs(int winSize, int minCover, int minFreq, double minFreqRatio, int numHaps, vector<double>& freqs) {
	int p;
	int maxColIdx = maxColIndex();
	for (p=0; p<=maxColIdx-winSize; p++) {
	}
}
