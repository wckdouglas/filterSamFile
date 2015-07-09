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
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> <modification reference fasta> ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - <modification reference fasta> ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
    cerr << endl;
	return 0;
}



//printing all variables
int printingTable(string transcriptID, string mispos, string ref, 
				int cov, string  modifiedBase, 
				int A, int C, int T, int G, 
				int insertion, int deletion, int refCount, int block)
{
    int i;
    // correct format for read data in R
    if (modifiedBase ==  "\\")
    {
        modifiedBase = "\\\\";
    }
    else if (modifiedBase =="\"")
    {
          modifiedBase = "\\\"";
    }
    else if (modifiedBase == "\'")
    {
          modifiedBase = "\\\'";
    }

    // print lines
    cout << transcriptID << "\t" << mispos << "\t" << ref << "\t";
    cout << cov + block<< "\t"; 
	cout << modifiedBase << "\t";
	cout << A << "\t" << C << "\t";
	cout << T << "\t" << G  << "\t" ;
	cout << insertion << "\t" << deletion;
	cout << "\t" << block;
	cout <<'\n';
    return 0;
}

// processing lines with mismatches 
int extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, 
					string ref, string modifiedBase, int qualThreshold, int coverageThreshold)
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
	cov += deletion;
	if (cov > coverageThreshold && refCount != cov)
	{
		assert (N + A + T + G + C + 
				a + c + t + g + n + 
				refCount + deletion == cov);
		printingTable(transcriptID, mispos, ref, cov, modifiedBase, 
					A, C, T, G, insertion, deletion, refCount,end);
	}	
    return 0;
}


// extract from each line different columns
// and give them to further processing
int processLine( stringList columns, seq_map &seqIndex, int qualThreshold, int coverageThreshold) 
{
    if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
    {
        string transcriptID, pos, ref, reads, baseQuals, modifiedBase;
        int cov;
        if (columns.size() == 6) 
        {
            cov = atoi(columns[3].c_str());
			transcriptID = columns[0];
			pos = columns[1];
            modifiedBase = seqIndex[transcriptID][atoi(pos.c_str()) - 1]  ;
			if (modifiedBase != "A" && modifiedBase != "C" && modifiedBase != "G" && modifiedBase != "T" && cov > coverageThreshold)
			{
				ref = columns[2];
				reads = columns[4];
				baseQuals = columns[5];
				assert ( baseQuals.length() == cov ) ;
				extractMismatches(reads, baseQuals, cov, transcriptID, 
							pos, ref, modifiedBase, qualThreshold, coverageThreshold);
			}
        }
    }
    return 0;
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename, seq_map &seqIndex, int qualThreshold, int coverageThreshold)
{
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, seqIndex, qualThreshold, coverageThreshold);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream(seq_map &seqIndex, int qualThreshold, int coverageThreshold)
{
	int prePos;
	int preEnd;
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, seqIndex, qualThreshold, coverageThreshold);
    }
    return 0;
}

void printHeader()
{
    cout << "transcriptID" << "\t" ;
    cout << "mispos" << "\t";
    cout << "ref" << "\t";
    cout << "cov" << "\t";
    cout << "modifiedBase" << '\t';
	cout << "A\tC\tT\tG\tinsertion" << '\t';
	cout << "deletion\tblock" << endl;
}

// main function
int main(int argc, char *argv[])
{
	ios::sync_with_stdio(false);
    // warnings
    if (argc != 5)
    {
        usage(argv);
		return 0;
    }

    // create modified RNA index
    const char* modifiedFa = argv[2];
    ifstream fastaFile (modifiedFa);
    string id, line;
    stringList seqList, idList;   
	int qualThreshold, coverageThreshold;

    while ( getline(fastaFile,line) )
    {
        if (line.at(0) == '>')
        {
            id = line.erase(0,1);
            stringList seqid = split(id,' ');
            idList.push_back(seqid[0]);
        }
        else
        {
            seqList.push_back(line);
        }
    }
    seq_map seqIndex;
    for (int i = 0 ; i < seqList.size(); i++)
    {
        seqIndex[idList[i]] = seqList[i];
    }

    printHeader();
	qualThreshold = atoi(argv[3]);
	coverageThreshold = atoi(argv[4]);
    // read lines
    if (strcmp(argv[1],"-") == 0)
    {
        readStream(seqIndex, qualThreshold, coverageThreshold);
    }
    else
    {
        const char* filename = argv[1];
        readFile(filename,seqIndex, qualThreshold, coverageThreshold);
    }
    return 0;
}
