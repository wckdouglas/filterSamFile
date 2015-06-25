#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <cassert>
#include "stringManipulation.h"
#include "pileupFix.h"

typedef map<string, string> seq_map;

//print usage 
int usage(char *argv[])
{
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
    cerr << endl;
	return 0;
}



//printing all variables
void printTable(string transcriptID, string mispos, string ref, 
				int cov, int A, int C, int T, int G, 
				int a, int c, int t, int g,
				int insertion, int deletion, int refCount)
{
    if (A + C + T + G > 0)
	{
	    cout << transcriptID << "\t" ;
		cout << mispos << "\t"; 
		cout << atoi(mispos.c_str())+1 << '\t';
		cout << ref << "\t";
		cout << cov << "\t"; 
		cout << '+' << "\t";
		cout << A << "\t" << C << "\t";
		cout << T << "\t" << G  << "\t" ;
		cout << insertion << "\t" << deletion;
		cout << '\n';
	}
	else if (a + c + t + g > 0)
	{
	    cout << transcriptID << "\t" ;
		cout << mispos << "\t"; 
		cout << atoi(mispos.c_str())+1 << '\t';
		cout << ref << "\t";
		cout << cov << "\t"; 
		cout << '-' << "\t";
		cout << a << "\t" << c << "\t";
		cout << t << "\t" << g << "\t" ;
		cout << insertion << "\t" << deletion;
		cout << '\n';
	}
}

// processing lines with mismatches 
void extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, 
					string ref, int qualThreshold, int coverageThreshold)
{
    int start = 0, end = 0, i = 0;
	int A = 0, C = 0, T = 0, G = 0, N = 0; 
	int a = 0, c = 0, t = 0, g = 0, n = 0; 
	int qual;
	int insertion = 0, deletion = 0, current = 0;
	int refCount = 0;
	fixpileup(A, C, T, G, N,
			a, c, t, g, n,
			deletion, insertion, reads, baseQuals,
			qualThreshold, cov, refCount, start, end);
	cov = cov +  deletion - n - N;
	
	if (cov > coverageThreshold && refCount != cov)
	{
		printTable(transcriptID, mispos, ref, cov,  
					A, C, T, G, a, c, t, g, 
					insertion, deletion, refCount);
		assert (A + T + G + C + 
				a + c + t + g + 
				refCount + deletion == cov);
	}
		
}


// extract from each line different columns
// and give them to further processing
void processLine( stringList columns, int qualThreshold, int coverageThreshold) 
{
    if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
    {
        string transcriptID, pos, ref, reads, baseQuals;
        int cov;
        if (columns.size() == 6) 
        {
            cov = atoi(columns[3].c_str());
            if (cov > coverageThreshold)
            { 
                transcriptID = columns[0];
                pos = columns[1];
                ref = columns[2];
                reads = columns[4];
                baseQuals = columns[5];
                assert ( baseQuals.length() == cov ) ;
				extractMismatches(reads, baseQuals, cov, transcriptID, 
								pos, ref, qualThreshold, coverageThreshold);
            }
        }
    }
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
void readFile(const char* filename, int qualThreshold, int coverageThreshold)
{
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold);
    }
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
void readStream(int qualThreshold, int coverageThreshold)
{
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold);
    }
}

void printHeader()
{
    cout << "chrom" << "\t" ;
    cout << "start" << "\t";
    cout << "end" << "\t";
    cout << "ref" << "\t";
    cout << "cov" << "\t";
	cout << "strand\t" ;
	cout << "A\tC\tT\tG\tinsertion\tdeletion" << endl;
}

// main function
int main(int argc, char *argv[])
{
	ios::sync_with_stdio(false);
    // warnings
    if (argc != 4)
    {
        usage(argv);
		return 0;
    }

    // create modified RNA index
	int qualThreshold, coverageThreshold;
	qualThreshold = atoi(argv[2]);
	coverageThreshold = atoi(argv[3]);
    // read lines
    printHeader();
    if (strcmp(argv[1],"-") == 0)
    {
		cerr << "Reading from stdin" << endl;
        readStream(qualThreshold, coverageThreshold);
    }
    else
    {
        const char* filename = argv[1];
		cerr << "Reading from: " << filename << endl;
        readFile(filename, qualThreshold, coverageThreshold);
    }
    return 0;
}
