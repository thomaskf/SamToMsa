CC = g++
SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o
CFLAGS = -O3
ZFLAGS = -lz

all : SamToMsa

$(SAMLIB) :	force_look
	cd $(SAMLIB); $(MAKE)

SamToMsa : mylib.cpp mylib.h reads.cpp reads.h partInfo.h partInfo.cpp constructHaplo.h constructHaplo.cpp fileHander.h fileHander.cpp samToMsa.cpp $(SAMOBJLIBS)
	$(CC) $(CFLAGS) -DNOMAIN -o SamToMsa samToMsa.cpp mylib.cpp reads.cpp partInfo.cpp constructHaplo.cpp fileHander.cpp $(SAMOBJLIBS) $(ZFLAGS)

clean :
	rm -f SamToMsa $(SAMOBJLIBS)
