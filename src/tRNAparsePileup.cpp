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
    cerr << "usage: "<< argv[0] << " <filename>|<stdin> <modification reference fasta>" << endl;
    cerr << endl;
    cerr << "Takes in mpileup result format:"<<endl;
    cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
    cerr << argv[0] << " - <modification reference fasta>"<<endl;
    cerr << endl;
	return 0;
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
                    int cov,string baseQuals,int start , int end, string modifiedBase)
{
    char read;
    char strand;
    int qual;
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
    for (i = 0 ; i < correctedReads.length(); i++)
    {
        read = correctedReads[i];
        if (read == 'A' || read == 'C' || read == 'T' || read == 'G')
        {
            strand = '+';
            qual = baseQuals[i] - 33 ;
            cout << transcriptID << "\t" << mispos << "\t" << ref << "\t";
            cout << read << "\t" << cov <<  "\t" << qual << "\t" ;
            cout << modifiedBase << endl;
        }
    }
    return 0;
}

// processing lines with mismatches 
int extractMismatches(string reads, string baseQuals, int cov, 
                    string transcriptID, string mispos, string ref,
                    string modifiedBase)
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
    assert (correctedReads.size() == cov);
    printingTable(transcriptID, mispos, ref, correctedReads, cov, baseQuals,start,end,modifiedBase);
    return 0;
}


// extract from each line different columns
// and give them to further processing
int processLine( lists columns, seq_map seqIndex) 
{
    if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
    {
        string transcriptID, pos,ref,reads,baseQuals, modifiedBase;
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
                modifiedBase = seqIndex[transcriptID][atoi(pos.c_str()) - 1]  ;
                assert ( baseQuals.length() == cov ) ;
                extractMismatches(reads,baseQuals,cov, transcriptID,pos,ref,modifiedBase);
            }
        }
    }
    return 0;
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename, seq_map seqIndex)
{
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        lists columns = split(line,'\t');
        processLine(columns, seqIndex);
    }
    return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream(seq_map seqIndex)
{
    for (string line; getline(cin, line);)
    {
        lists columns = split(line,'\t');
        processLine(columns, seqIndex);
    }
    return 0;
}

int printHeader()
{
    cout << "transcriptID" << "\t" ;
    cout << "mispos" << "\t";
    cout << "ref" << "\t";
    cout << "read" << "\t";
    cout << "cov" << "\t";
    cout << "baseQual" << "\t";
    cout << "modifiedBase" << '\n';
    return 0;
}

// main function
int main(int argc, char *argv[])
{
	ios::sync_with_stdio(false);
    // warnings
    if (argc != 3)
    {
        usage(argv);
		return 0;
    }

    // create modified RNA index
    const char* modifiedFa = argv[2];
    ifstream fastaFile (modifiedFa);
    string id, line;
    lists seqList, idList;   

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
    // read lines
    if (strcmp(argv[1],"-") == 0)
    {
        readStream(seqIndex);
    }
    else
    {
        const char* filename = argv[1];
        readFile(filename,seqIndex);
    }
    return 0;
}
