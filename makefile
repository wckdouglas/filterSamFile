cpp=g++
mkdir_p=mkdir -p


all: directory filterSoftClipped parsePileup sam2fastq pileup2bed bedToJunction

directory:
	$(mkdir_p) bin

filterSoftClipped:
	$(cpp)  src/filterSam.cpp -o bin/filterSoftClipped

parsePileup:
	$(cpp) src/parsePileup.cpp -o bin/parsePileup
	$(cpp) src/tRNAparsePileup.cpp -o bin/tRNAparsePileup

sam2fastq:
	$(cpp) src/sam2fastq.cpp -o bin/sam2fastq

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed

bedToJunction:
	$(cpp) src/bedToJunction.cpp -o bin/bedToJunction -std=c++0x
