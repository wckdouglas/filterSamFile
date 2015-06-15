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

using namespace std;
typedef map<string, string> seq_map;
typedef vector <string> lists;

//print usage 
int usage(char *argv[])
{
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> <modification reference fasta> ";
	cerr << "quality threshold> <coverage threshold>" << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - <modification reference fasta>"<<endl;
    cerr << endl;
    abort();
}


//split function to split line with desired deliminator
lists split(const string &s, char delim) 
{
        stringstream ss(s);
        string item;
		lists result;
        while (getline(ss, item, delim)) 
        {
                result.push_back(item);
        }
        return result;
}


//printing all variables
int printingTable(string transcriptID, string mispos, string ref, 
				int cov, string  modifiedBase, 
				int A, int C, int T, int G, 
				int insertion, int deletion)
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
	ios::sync_with_stdio(false);
    cout << transcriptID << "\t" << mispos << "\t" << ref << "\t";
    cout << cov << "\t"; 
	cout << modifiedBase << "\t";
	cout << A << "\t" << C << "\t";
	cout << T << "\t" << G  << "\t" ;
	cout << insertion << "\t" << deletion << "\t" << endl;
    return 0;
}

// processing lines with mismatches 
int extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, 
					string ref, string modifiedBase, int qualThreshold, int coverageThreshold)
{
    string correctedReads; 
	char readPos;
    int start = 0, end = 0, i = 0;
	int A = 0, C = 0, T = 0, G = 0, N = 0; 
	int qual, j = 0, k = 0, l = 0, trycount = 0;
	int insertion = 0, deletion = 0, current = 0;
	int refCount = 0;
    while (i < reads.length())
    {
		readPos = reads.at(i);
        if (readPos == '+')
		//insertion
        {
			i ++ ; 
			current = 0;
			insertion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			if (current > 10)
			{
				i += current - 2;
			}
			else
			{
				i += current - 1;
			}
        }
		else if (readPos == '-')
		// deletion
		{
			i ++ ; 
			current = 0;
			deletion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			if (current > 10)
			{
				i += current - 2;
			}
			else
			{
				i += current - 1;
			}
		}
        else if (readPos == '^')
        {
            i ++;
            start = 1;
        }
        else if (readPos == '$')
        {
	        //i ++;
            end = 1;
        }
		else 
		{
            qual = baseQuals[j] - 33 ;
			j++;
			if (qual < qualThreshold || readPos == '*')
			{
				cov = cov - 1;
			}
			else 
			{
				if (readPos == 'A')
				{	
					A ++;
				}
				else if (readPos == 'C')
				{
					C ++;
				}
				else if(readPos == 'G')
				{
					G ++;
				}
				else if (readPos == 'T')
				{
					T ++;
				}
				else if (readPos == 'N')
				{
					N ++;
				}
				else if (readPos == '.')
				{
					refCount ++;
				}
			}
		}
		i++;
    }
	cov += deletion;
	if (cov > coverageThreshold && refCount != cov)
	{
		printingTable(transcriptID, mispos, ref, cov, modifiedBase, 
					A, C, T, G, insertion, deletion);
		assert (N + A + T + G + C + refCount + deletion == cov);
	}	
    return 0;
}


// extract from each line different columns
// and give them to further processing
int processLine( lists columns, seq_map seqIndex, int qualThreshold, int coverageThreshold) 
{
    if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
    {
        string transcriptID, pos, ref, reads, baseQuals, modifiedBase;
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
                modifiedBase = seqIndex[transcriptID][atoi(pos.c_str()) - 1]  ;
                assert ( baseQuals.length() == cov ) ;
				if (modifiedBase != "A" && modifiedBase != "C" && modifiedBase != "G" && modifiedBase != "T")
				{
					extractMismatches(reads, baseQuals, cov, transcriptID, 
									pos, ref, modifiedBase, qualThreshold, coverageThreshold);
				}
            }
        }
    }
    return 0;
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename, seq_map seqIndex, int qualThreshold, int coverageThreshold)
{
	ios::sync_with_stdio(false);
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        lists columns = split(line,'\t');
        processLine(columns, seqIndex, qualThreshold, coverageThreshold);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream(seq_map seqIndex, int qualThreshold, int coverageThreshold)
{
	ios::sync_with_stdio(false);
    for (string line; getline(cin, line);)
    {
        lists columns = split(line,'\t');
        processLine(columns, seqIndex, qualThreshold, coverageThreshold);
    }
    return 0;
}

int printHeader()
{
	ios::sync_with_stdio(false);
    cout << "transcriptID" << "\t" ;
    cout << "mispos" << "\t";
    cout << "ref" << "\t";
    cout << "cov" << "\t";
    cout << "modifiedBase" << '\t';
	cout << "A\tC\tT\tG\tinsertion\tdeletion" << endl;
    return 0;
}

// main function
int main(int argc, char *argv[])
{
    // warnings
    if (argc != 5)
    {
        usage(argv);
    }

    // create modified RNA index
    const char* modifiedFa = argv[2];
    ifstream fastaFile (modifiedFa);
    string id, line;
    lists seqList, idList;   
	int qualThreshold, coverageThreshold;

    while ( getline(fastaFile,line) )
    {
        if (line.at(0) == '>')
        {
            id = line.erase(0,1);
            lists seqid = split(id,' ');
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
