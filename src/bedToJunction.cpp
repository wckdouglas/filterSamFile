#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <iomanip>
#include <cstring>
#include <map>

#include "stringManipulation.h"

using namespace std;

typedef map<string,int> stringHash;
typedef map<string,int>::iterator it_hash;


int processing(string line, stringHash &junctionList)
{
	stringList columns = split(line,'\t');
	if (columns.size() < 7)
	{
		cerr << "Input does not have 7 columns!! missing -cigar option?"<< endl;
		return 1;
	}
	string cigar = columns[6];
	if (cigar.find('N') != string::npos)
	{
		string chrom = columns[0];
		int readStart = atoi(columns[1].c_str());
		int readEnd = atoi(columns[2].c_str());
		string strand = columns[5];
		int initial = readStart;
		int junctionStart, junctionEnd, junctionLength;
		int junctionNumber;
		numList cigarNum(0);
		stringList cigarStr(0);
		regexSeparate(cigar,cigarNum,cigarStr);
		int fracs = cigarNum.size();
		numList cumCigarNum(fracs);
		for (int i = 0; i < fracs ; i ++)
		{
			initial = initial + cigarNum[i];
			cumCigarNum[i] = initial;
			if (cigarStr[i].compare("N")==0)
			{
				junctionNumber ++;
				junctionStart = cumCigarNum[i-1];
				junctionEnd = cumCigarNum[i];
				junctionLength = junctionEnd - junctionStart; 
				assert (junctionLength == cigarNum[i]);
				string junctionRecord = chrom + "_" + to_string(junctionStart) + "_" + to_string(junctionEnd)+'_'+strand;
				if (junctionList.find(junctionRecord) != junctionList.end())
				{
					junctionList[junctionRecord] += 1;
				}
				else
				{
					junctionList[junctionRecord] = 1;
				}
			}
		}
	}
	return 0;
}

// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
int readFile(const char* filename, stringHash &junctionList)
{
	int pass;
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        pass = processing(line,junctionList);
		if (pass == 1)
		{
			return 1;
		}
    }
	return 0;
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
int readStream(stringHash &junctionList)
{
	int pass;
    for (string line; getline(cin, line);)
    {
        pass = processing(line, junctionList);
		if (pass == 1)
		{
			return 1;
		}
    }
	return 0;
}

//usage
void usage(char *argv[])
{
	cerr << "usage " << argv[0] << " <filename>|<stdin>"<< "\n\n";
	cerr << "suggested usage: bamtobed -i <bamfile> -cigar | " << argv[0] << " - > junction.bed" << '\n';
	cerr << "**** use <-> when using stdin" << '\n';
	cerr << "Needed a bed file with cigar string" << '\n';
	cerr << "output file contains six columns: \n";
	cerr << "          column 1:        chromosome name\n";
	cerr << "          column 2:        junction start position (1-base start pos, same as bamfile)\n";
	cerr << "          column 3:        junction end position (1-base start pos, same as bamfile)\n";
	cerr << "          column 4:        name (ordered by chrom)\n";
	cerr << "          column 5:        number of reads supporting\n";
	cerr << "          column 6:        strand [+/-]\n";
	cerr << endl;
}


//main 
int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		usage(argv);
		return 0;
	}
	int j = 0,pass;
	string key,chrom,strand;
	int start,end;
	stringHash junctionList;

	//determine stdin or filename input
    if (strcmp(argv[1],"-") == 0)
    {
		cerr << "Reading from stdin" << endl;
        pass = readStream(junctionList);
    }
    else
    {
        const char* filename = argv[1];
		cerr << "Reading from: " << filename << endl;
        pass = readFile(filename,junctionList);
    }

	//printing out bed file
	if (pass == 1)
	{
		return 1;
	}
	for(it_hash iterator = junctionList.begin(); iterator != junctionList.end(); iterator++)
	{
		j ++;
		key = iterator->first; //chrom , start, end
		stringList info;
		info = split(key,'_');
		chrom = info[0];
		start = atoi(info[1].c_str());
		end = atoi(info[2].c_str());
		strand = info[3];
		cout << chrom << '\t' << start <<'\t' << end << '\t'; // chromosome name, start site, end site
		cout <<	"JUNC" << setfill('0') << setw(10) << j << '\t'; // junction name
		cout << iterator->second << '\t';	// supported by how many reads
		cout << strand << '\t'<< end - start <<'\n' ;// strand
	}
	cerr << "Finished!! " << endl;
	return 0;
}
