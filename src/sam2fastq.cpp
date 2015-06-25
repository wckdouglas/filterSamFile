#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <string.h>
#include "stringManipulation.h"

using namespace std;

//print usage 
int usage(char *argv[])
{
    cerr << "usage: "<< argv[0] << " <sam file>|<stdin> " << endl;
    cerr << endl;
    cerr << "Takes in samFile:"<<endl;
    cerr << "samtools view <bamFile> | ";
    cerr << argv[0] << " - "<<endl;
	return 1;
}


//process lines
int printSeq(vector<string> columns, int i)
{
    string id, sequence, qual, chrom;
    chrom = columns[2];
    if (chrom == "*")
    {
        id = columns[0];
        sequence = columns[9];
        qual = columns[10];
        if (sequence.length() != qual.length())
        {
            cerr << "Wrong columns!!" << endl;
            abort();
        }
        i += 1;
		printf("@%s\n%s\n+\n%s\n",id.c_str(),sequence.c_str(),qual.c_str());
    }
    return i;
}

int readFile(const char* filename)
{
    int i = 0;
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        if (line[0] != '@') 
        {
            stringList columns = split(line,'\t');
            i += printSeq(columns, i);
        }
    }
    return i;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream()
{
    int i = 0;
    for (string line; getline(cin, line);)
    {
        if (line[0] != '@') 
        {
            stringList columns = split(line,'\t');
            i += printSeq(columns, i);
        }
    }
    return i;
}


int main(int argc, char *argv[])
{
	ios::sync_with_stdio(false);
    // warnings
    if (argc != 2)
    {
        usage(argv);
		return 0;
    }

    // read lines
    int i;
    if (strcmp(argv[1],"-") == 0)
    {
        i = readStream();
    }
    else
    {
        const char* filename = argv[1];
        i = readFile(filename);
    }
    cerr << "Written "<< i <<" sequences." << endl;
    return 0;
}
