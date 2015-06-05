#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

using namespace std;

typedef vector<string> stringList;
typedef vector<int> numList;

//split function to split line with desired deliminator
stringList &split(const string &s, char delim, stringList &result) 
{
        stringstream ss(s);
        string item;
        while (getline(ss, item, delim)) 
        {
                result.push_back(item);
        }
        return result;
}


// separate number and string
int regexSeparate(string tag,numList &numberslists, stringList &letterlists)
{
    int i, current = 0, size = tag.length();
    string letters = "";
    char c;
    for ( i = 0; i < size; i ++)
    {
        c = tag.at(i);
        if (isdigit(c))
        {
            current = current * 10 + (c-'0');
            letters = "";
            if (i == size-1)
            {
                numberslists.push_back(current);
                letterlists.push_back("N");
            }
        }
        else
        {
            numberslists.push_back(current);
            letters = letters + string(1,c);
            letterlists.push_back(letters);
            current = 0;
        }
    }
    return 0;
}

int allClipped (numList cigarNumber, stringList cigarLetter)
{
    int totalClipped = 0;
    for (int i = 0; i < cigarNumber.size(); i++)
    {
        if (cigarLetter.at(i) == "S")
        {
            totalClipped += cigarNumber[i];
        }
    }
    return totalClipped;
}

int filter(stringList columns, double singleEndSoftclippedThreshold, double bothEndSoftclippedThreshold)
{

    // define variables
    numList cigarNumber;
    stringList cigarLetter;
    double thresholdBoth, thresholdSingle;
    int cigarSize, seqLength, totalClip, passFlag = 1;
    string cigar;

    //get elements from cigar string
    cigar = columns[5];
    regexSeparate(cigar, cigarNumber, cigarLetter);
    cigarSize = cigarNumber.size();
    seqLength = accumulate(cigarNumber.begin(),cigarNumber.end(),0);

    //set threshold
    thresholdSingle = seqLength * singleEndSoftclippedThreshold;
    thresholdBoth = seqLength * bothEndSoftclippedThreshold;
    
    totalClip = allClipped(cigarNumber, cigarLetter);

    //filtering, if passed, return line, else return 0
    if (cigarLetter.at(0) == "S" )
    {
        if (cigarNumber[0] > thresholdSingle)
        {
            passFlag = 0;
        }
        else if (cigarLetter.at(cigarSize-1) == "S" && totalClip > thresholdBoth)
        {
            passFlag = 0;
        }
    }
    else if (cigarLetter.at(cigarSize-1) == "S" )
    {
        if (cigarNumber[cigarSize-1] > thresholdSingle)
        {
            passFlag = 0;
        }
    }
    return passFlag;
}

int pairedEndProcessLine(string line, int lineno, int headerline, stringList &samlines,
                double singleEndSoftclippedThreshold, double bothEndSoftclippedThreshold,
                stringList &ids, int debugging, numList &passFlags)
{
    //headers 
    if (line.at(0) == '@')
    {
        cout << line << endl;
        headerline ++;
    }
    else
    {
        //filter each sam lines
        string id1, id2;
        int lineNum = lineno - headerline, passFlagTotal;
        stringList columns;
        split(line,'\t',columns);

        //first in pair
        if (remainder(lineNum, 2) != 0)
        {
            id1 = columns[0];
            ids[0] = id1;
            samlines[0] = line;
            passFlags[0] = filter(columns, singleEndSoftclippedThreshold, bothEndSoftclippedThreshold);
        }
        //second in pair
        else
        {
            id2 = columns[0];
            ids[1] = id2;
            samlines[1] = line;
            if ( ids[0] != ids[1] ) 
            {
                cerr << "sam input is not sorted!!" << endl;
                exit(EXIT_FAILURE);
            }
            passFlags[1] = filter(columns,singleEndSoftclippedThreshold, bothEndSoftclippedThreshold);

            //print if both of them passed the filter
            //do not contain "notpass"
            passFlagTotal = accumulate(passFlags.begin(),passFlags.end(),0);
            if (debugging == 0 && passFlagTotal == 2)
            {
                cout << samlines[0] << endl;
                cout << samlines[1] << endl;
            }
            else if (debugging == 1 && passFlagTotal < 2)
            {
                cout << samlines[0] << endl;
                cout << samlines[1] << endl;
            }
            stringList samlines(2);
            stringList ids(2);
            numList passFlags(2);
        }
    }
    return 0;
}

int pairedStreamFile(double singleEndSoftclippedThreshold, double bothEndSoftclippedThreshold, int debugging)
{
    int lineno = 0, headerline = 0;
    stringList samlines(2), ids(2);
    numList passFlags(2);
    for ( string line ; getline(cin, line);)
    {
        lineno ++;
        pairedEndProcessLine(line, lineno, headerline, samlines, singleEndSoftclippedThreshold, bothEndSoftclippedThreshold, ids, debugging, passFlags);
    }
	cerr << "Parsed " << lineno << " lines " << endl;
    return 0;
}

int singleStreamFile(double singleEndSoftclippedThreshold, double bothEndSoftclippedThreshold, int debugging)
{
    int pass;
	int lineno = 0;
    for ( string line ; getline(cin, line);)
    {
		lineno ++;
        stringList columns;
        split(line,'\t',columns);
        if (line.at(0) == '@')
        {
            cout << line << endl;
        }
        else
        {
            pass = filter(columns, singleEndSoftclippedThreshold, bothEndSoftclippedThreshold);
            if (pass != 0 && debugging == 0)
            {
                cout << line << endl;
            }
            else if (pass == 0 && debugging == 1)
            {
                cout << line << endl;
            }
        }
    }
	cerr << "Parse " << lineno  << " lines" << endl;
    return 0;
}

int usage(char *program)
{
	cerr << "****************************************************************" << '\n';
	cerr << "Filtering soft clipped reads from paired-end RNA-seq sam files" << '\n';
	cerr << "usage: cat <samFile> | " << program << " -s <oneSideSoftclipFractionThreshold> ";
	cerr << "-b <bothEndSoftclippedThreshold> [-vp]" << "\n" << endl;
	cerr << "<oneSideSoftclipFractionThreshold>" << "\t" << "Threshold for filtering one side softclip sequence. "<< endl;
    cerr << "                                  " << "\t" << "Must be between 0 and 1 [default: 0.3]" << endl;
	cerr << "<bothEndSoftclippedThreshold>     " << "\t" << "Threshold for filtering both side softclip sequence. "<< endl;
    cerr << "                                  " << "\t" << "Must be between 0 and 1 [default: 0.4]" << endl;
	cerr << "-v                                " << "\t" << "Debugging mode: print out all failed alignments" << endl;
	cerr << "-p                                " << "\t" << "paired-end mode [default = single end]" << endl;
	cerr << "If the soft clipped bases count > (threshold * [whole sequence length]), it will be filter out" << endl;
	cerr << "****************************************************************\n";
    cerr << endl;
    exit(EXIT_FAILURE);
	return 0;
}

int main(int argc, char **argv)
{
    int c;
    double singleEndSoftclippedThreshold = 0.3;
    double bothEndSoftclippedThreshold = 0.4;
    int debugging = 0, pairedEndFlag = 0;
    char *program = argv[0];
    if (argc < 2)
    {
        usage(program);
    }
    while ((c= getopt(argc, argv, "s:b:vp")) != -1)
    {
        switch(c)
        {
            case 's':
                singleEndSoftclippedThreshold = atof(optarg);
                break;
            case 'b':
                bothEndSoftclippedThreshold = atof(optarg);
                break;
            case 'v':
                debugging = 1;
                break;
            case 'p':
                pairedEndFlag = 1;
                break;
            case '?':
                usage(program);
                break;
            default:
                usage(program);
        }
    }
    if (pairedEndFlag == 1)
    {
        pairedStreamFile(singleEndSoftclippedThreshold,bothEndSoftclippedThreshold,debugging);
    }
    else
    {
        singleStreamFile(singleEndSoftclippedThreshold,bothEndSoftclippedThreshold,debugging);
    }
    return 0;
}
    
