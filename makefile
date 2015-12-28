cpp=g++
sh=bash
mkdir_p=mkdir -p
include=./bamtools/include
lib=./bamtools/lib


all: directory bamtools filterSoftClipped parsePileup sam2fastq pileup2bed bedToJunction

directory:
	$(mkdir_p) bin

filterSoftClipped:
	$(cpp)  src/filterSam.cpp -o bin/filterSoftClipped  -I $(include) -L $(lib) -lbamtools

parsePileup:
	$(cpp) src/parsePileup.cpp -o bin/parsePileup
	$(cpp) src/tRNAparsePileup.cpp -o bin/tRNAparsePileup

sam2fastq:
	$(cpp) src/sam2fastq.cpp -o bin/sam2fastq

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed

bedToJunction:
	$(cpp) src/bedToJunction.cpp -o bin/bedToJunction -std=c++0x
