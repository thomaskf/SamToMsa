# SamToMsa

SamToMsa updates all the read sequences by inserting gaps to form a multiple sequence alignment. SamToMsa also produces a longer read by connecting two overlapped reads of the same pair.

[Installation]

The software was written in C++, and it has been tested under linux and MacOS platform. You need to have C++ compiler installed in the machine in order to compile the source codes. The compilation steps are shown as follows:

$ cd SamToMsa
$ make

Then an executable file named SamToMsa will appear

[Usage]

Syntax: ./SamToMsa <sam/bam file> <pair-end out file> <single-end out file>
