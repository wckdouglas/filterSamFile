cpp=g++
sh=bash
mkdir_p=mkdir -p
include=./bamtools/include
lib=./bamtools/lib

<<<<<<< HEAD

all: directory bamtools filterSoftClipped parsePileup sam2fastq pileup2bed bedToJunction
=======
all: directory filterSoftClipped parsePileup sam2fastq pileup2bed bedToJunction pileup2base
>>>>>>> ce5cd3ad52483f36d47d64ffdad87c6ea6f961e0

directory:
	$(mkdir_p) ./bin

filterSoftClipped:
	$(cpp)  src/filterSam.cpp -o bin/filterSoftClipped  -I $(include) -L $(lib) -lbamtools

parsePileup:
	$(cpp) src/parsePileup.cpp -o bin/parsePileup
	$(cpp) src/tRNAparsePileup.cpp -o bin/tRNAparsePileup

sam2fastq:
	$(cpp) src/sam2fastq.cpp -o bin/sam2fastq

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed

pileup2base:
	$(cpp) src/pileupBases.cpp -o bin/pileup2base

bedToJunction:
	$(cpp) src/bedToJunction.cpp -o bin/bedToJunction -std=c++0x
