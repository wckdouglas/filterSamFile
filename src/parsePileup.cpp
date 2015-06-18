#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>
#include <cstdlib>
#include <fstream>
#include <cassert>

using namespace std;
typedef vector<string> stringList;

//print usage 
int usage(char *argv[])
{
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> " << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - "<<endl;
    abort();
}


//split function to split line with desired deliminator
stringList split(const string &s, char delim) 
{
        stringstream ss(s);
        string item;
		stringList result;
        while (getline(ss, item, delim)) 
        {
                result.push_back(item);
        }
        return result;
}

int countDigits(int number) 
{
	if (number < 10) 
	{
		return 1;
	}
	int count = 0;
	while (number > 0) 
	{
		number /= 10;
		count++;
	}
	return count;
}

//printing all variables
int printingTable(string transcriptID, string mispos, string ref, string correctedReads,
                    int cov,string baseQuals,int start , int end)
{
    char read;
    char strand;
    int qual;
    int i;
    for (i = 0 ; i < correctedReads.length(); i++)
    {
        read = correctedReads.at(i);
        if (read == 'A' || read == 'C' || read == 'T' || read == 'G')
        {
            strand = '+';
        }
		else if (read == 'a' || read == 'c' || read == 't' || read == 'g')
		{
			strand = '-';
		}
        qual = baseQuals[i] - 33 ;
		ios::sync_with_stdio(false);
        cout << transcriptID << "\t" << mispos << "\t" << ref << "\t";
        cout << read << "\t" << cov <<  "\t" << qual << "\t" ;
        cout << strand << "\t" << start << "\t" << end << "\n";
    }
    return 0;
}

// processing lines with mismatches 
int extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, string ref)
{
    string correctedReads; 
    string skip;
    int start = 0, end = 0, i = 0;
	char readPos;
	int current = 0;
    while (i < reads.length())
    {
		readPos = reads.at(i);
        if (reads[i] == '+' || reads[i] == '-')
        {
			i ++ ; 
			current = 0;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			i += current - countDigits(current);
        }
        else if (reads[i] == '^')
        {
            i += 2;
            start = 1;
        }
        else if (reads[i] == '$')
        {
            i += 1;
            end = 1;
        }
        else 
        {
            correctedReads.push_back(reads[i]);
            i ++;
        }
    }
    if (correctedReads.size() != cov)
    {
        cout << " wrong program parse mismatch strings!!\n " << endl;
        abort();
    }
    printingTable(transcriptID, mispos, ref, correctedReads, cov, baseQuals,start,end);
    return 0;
}


// extract from each line different columns
// and give them to further processing
int processLine(stringList columns) 
{
    string transcriptID, pos, ref, reads, baseQuals;
    int cov;
    if (columns.size() == 6) 
    {
        cov = atoi(columns[3].c_str());
        if (cov > 0)
        { 
            transcriptID = columns[0];
            pos = columns[1];
            ref = columns[2];
            reads = columns[4];
            baseQuals = columns[5];
            assert (baseQuals.length() == cov) ;
			extractMismatches(reads,baseQuals,cov, transcriptID,pos,ref);
        }
    }
    return 0;
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename)
{
	ios::sync_with_stdio(false);
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream()
{
	ios::sync_with_stdio(false);
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns);
    }
    return 0;
}

int printHeader()
{
    cout << "transcriptID\t" ;
    cout << "mispos\t";
    cout << "ref\t";
    cout << "read\t";
    cout << "cov\t";
    cout << "baseQual\t";
    cout << "strand\t";
    cout <<"start\t";
    cout << "end" << endl;
    return 0;
}

// main function
int main(int argc, char *argv[])
{
    // warnings
    if (argc != 2)
    {
        usage(argv);
    }

    printHeader();
    // read lines
    if (strcmp(argv[1],"-") == 0)
    {
        readStream();
    }
    else
    {
        const char* filename = argv[1];
        readFile(filename);
    }
    return 0;
}
