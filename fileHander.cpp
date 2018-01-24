#include "fileHander.h"


// for BAM file
string getNext(int& pos, char* rec_seq, int rec_len) {
	string nextStr = "";
	
	// skip the spaces
	while (pos < rec_len && rec_seq[pos] == '\t')
		pos++;
	
	// get the sequence
	while (pos < rec_len && rec_seq[pos] != '\t') {
		nextStr.append(1, rec_seq[pos]);
		pos++;
	}
	return nextStr;
}

void SamBamFileHander::openFile(char* fileName) {
	
	// first assume the file is in BAM format
	bamQueryFile = bam_open ( fileName, "r" );
	bamHeader = bam_header_init();
	int i;
	string aline;
	
	i = bgzf_check_EOF(bamQueryFile);
	if (i == 0) {
		// cout << "The format of inut file: SAM" << endl;
		type = 2; // SAM format
		bam_close(bamQueryFile);
		fin.open(fileName);
	} else {
		// cout << "The format of inut file: BAM" << endl;
		type = 1; // BAM format
		bamHeader = bam_header_read ( bamQueryFile );
		bam = bam_init1();
	}
}

bool SamBamFileHander::getNextSeq(string& line) {
	bool b;
	char* rec_seq;
	int rec_len;
	int i;
	if (type==1) {
		// BAM format
		if (bam_read1 (bamQueryFile, bam) > 0) {
			rec_seq = bam_format1 (bamHeader, bam);
			rec_len = strlen(rec_seq);
			if (rec_len == 0)
				return false;
			line.resize(rec_len);
			for (i=0; i<rec_len; i++)
				line[i] = rec_seq[i];
			free(rec_seq);
		} else
			return false;
	} else {
		// SAM format
		// skip all the sentences starting with char "@"
		b = getline(fin, line);
		while (b && (line.length()==0 || (line.length() > 0 && line[0]=='@'))) {
			b = getline(fin, line);
		}
		if (!b)
			return false;
	}
	return true;
}

void SamBamFileHander::closeFile() {
	if (type==1) {
		// BAM format
		bam_close(bamQueryFile);
	} else {
		// SAM format
		fin.close();
	}
}
