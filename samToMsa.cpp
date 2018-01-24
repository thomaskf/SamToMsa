#include "constructHaplo.h"

int main(int argc, char** argv) {
	if (argc < 4) {
		cout << "This program will update all the read sequences by inserting gaps to form a" << endl;
		cout << "multiple sequence alignment. Moreover, it will produce a longer read by" << endl;
		cout << "connecting two overlapped reads in the same pair." << endl;
		cout << "Syntax: " << argv[0] << " [sam/bam file] [pair-end out file] [single-end out file]" << endl;
		exit(1);
	}
	
	int i, j;

	ConstructHaplo read_likes;
	
	// load the read and frequency files
	read_likes.loadReadFile(argv[1]);

	// correct the reads
	read_likes.readCorrection();

	// combine the overlap reads of the same pair into a long read
	// and update the frequencies
	read_likes.combineRdAndUpdateFreq();
	
	// output all the reads
	read_likes.rds.outReads(argv[2], argv[3]);
}