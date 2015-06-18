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
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> " ;
	cerr << "<quality threshold> <covThreshold>" << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - "<<endl;
	return 0;
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

string basepair( string base)
{
	string pairedBase;
	assert (base.length() == 1);
	if (base == "A")
	{
		pairedBase = "T"; 
	}
	else if (base == "C")
	{
		pairedBase = "G";
	}
	else if (base == "T")
	{
		pairedBase = "A";
	}
	else if (base == "G")
	{
		pairedBase = "C";
	}
	return pairedBase;
}

//printing all variables
int printingTable(string transcriptID, string mispos, string refbase, string correctedReads,
                    int cov,string baseQuals,int start , int end, int qualThreshold)
{
    char readbase = ' ';
    char strand = ' ';
    int qual;
    int i;
	int flag;
	ios::sync_with_stdio(false);
    for (i = 0 ; i < correctedReads.length(); i++)
    {
		qual = baseQuals[i] - 33 ;
		if (qual > qualThreshold)
		{
			flag = 0;
			readbase = correctedReads.at(i);
			if (readbase == 'A' || readbase == 'C' || readbase == 'T' || readbase == 'G')
			{
				strand = '+';
				flag = 1;
			}
			else if (readbase == 'a' || readbase == 'c' || readbase == 't' || readbase == 'g')
			{
				strand = '-';
				flag = 1;
				refbase = basepair(refbase);
			}
			if (flag  == 1)
			{
				cout << transcriptID << "\t" << mispos << "\t" << refbase << "\t";
				cout << readbase << "\t" << cov <<  "\t" << qual << "\t" ;
				cout << strand << "\t" << start << "\t" << end << '\n';
			}
		}
    }
    return 0;
}

// processing lines with mismatches 
int extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, string refbase, int qualThreshold)
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
			i += current - countDigits(current) + 1;
        }
        else if (readPos == '^')
        {
            i += 2;
            start = 1;
        }
        else if (readPos == '$')
        {
            i ++ ;
            end = 1;
        }
        else if (readPos != '<' && readPos != '>')
        {
            correctedReads.push_back(reads[i]);
            i ++;
        }
    }
    assert (correctedReads.size() == cov);
    printingTable(transcriptID, mispos, refbase, correctedReads, cov, baseQuals,start,end, qualThreshold);
    return 0;
}


// extract from each line different columns
// and give them to further processing
int processLine(stringList columns, int qualThreshold, int covThreshold) 
{
    string transcriptID, pos, refbase, reads, baseQuals;
    int cov;
    if (columns.size() == 6) 
    {
        cov = atoi(columns[3].c_str());
        refbase = columns[2];
        if (cov > covThreshold && refbase != "N")
        { 
            transcriptID = columns[0];
            pos = columns[1];
            reads = columns[4];
            baseQuals = columns[5];
            assert (baseQuals.length() == cov) ;
			extractMismatches(reads,baseQuals,cov, transcriptID,pos,refbase, qualThreshold);
        }
    }
    return 0;
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename,int qualThreshold, int covThreshold)
{
	ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns,qualThreshold, covThreshold);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream(int qualThreshold, int covThreshold)
{
	//ios::sync_with_stdio(false);
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, covThreshold);
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
    if (argc != 4)
    {
        usage(argv);
		return 0;
    }
	int qualThreshold = atoi(argv[2]);
	int covThreshold = atoi(argv[3]);

    printHeader();

    // read lines
    if (strcmp(argv[1],"-") == 0)
    {
        readStream(qualThreshold, covThreshold);
    }
    else
    {
        const char* filename = argv[1];
        readFile(filename,qualThreshold, covThreshold);
    }
    return 0;
}
