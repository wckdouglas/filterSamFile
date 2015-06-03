cppFlag=g++
all: filterSoftClipped

filterSoftClipped:
	$(cppFlag)  filterSam.cpp -o filterSoftClipped
	
