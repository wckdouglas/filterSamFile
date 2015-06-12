#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>
#include <cstdlib>
#include <fstream>

using namespace std;

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
vector<string> &split(const string &s, char delim, vector<string> &result) 
{
        stringstream ss(s);
        string item;
        while (getline(ss, item, delim)) 
        {
                result.push_back(item);
        }
        return result;
}


// given 3 char ouput the concatenate of the mas integer
string insertionNum(char a, char b, char c)
{
    string number;
    if (isdigit(c))
    {
        number.push_back(a);
        number.push_back(b);
        number.push_back(c);
    }
    else if (isdigit(b))
    {
        number.push_back(a);
        number.push_back(b);
    }
    else
    {
        number.push_back(a);
    }
    return number;
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
        read = correctedReads[i];
        if (read == 'A' || read == 'C' || read == 'T' || read == 'G')
        {
            strand = '+';
            qual = baseQuals[i] - 33 ;
            cout << transcriptID << "\t" << mispos << "\t" << ref << "\t";
            cout << read << "\t" << cov <<  "\t" << qual << "\t" ;
            cout << strand << "\t" << start << "\t" << end << endl;
        }
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
    while (i < reads.length())
    {
        if (reads[i] == '+' || reads[i] == '-')
        {
            skip = insertionNum(reads[i+1],reads[i+2],reads[i+3]);
            i += (strtol(skip.c_str(),0,10) + skip.length()  + 1);
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
int processLine( vector<string> columns) 
{
    string transcriptID, pos,ref,reads,baseQuals;
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
            if ( baseQuals.length() != cov )
            {
                cout << "Wrongly parsed on quality!\n";
                abort();
            }
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
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        vector<string> columns;
        columns = split(line,'\t',columns);
        processLine(columns);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream()
{
    for (string line; getline(cin, line);)
    {
        vector<string> columns;
        columns = split(line,'\t',columns);
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
