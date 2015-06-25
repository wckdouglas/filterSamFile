cpp=g++
mkdir_p=mkdir -p


all: directory filterSoftClipped parsePileup sam2fastq pileup2bed

directory:
	$(mkdir_p) bin

filterSoftClipped:
	$(cpp)  src/filterSam.cpp -o bin/filterSoftClipped

parsePileup:
	$(cpp) src/parsePileup.cpp -o bin/parsePileup
	$(cpp) src/tRNAparsePileup.cpp -o bin/tRNAparsePileup
	$(cpp) src/tRNAparsePileupNew.cpp -o bin/tRNAparsePileupNew

sam2fastq:
	$(cpp) src/sam2fastq.cpp -o bin/sam2fastq

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed
