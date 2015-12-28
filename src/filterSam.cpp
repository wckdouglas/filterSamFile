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
#include "stringManipulation.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace BamTools;
using namespace std;

int processAlignment(BamAlignment aln, double singleThreshold, double bothThreshold)
{
	int headClipped = 0, tailClipped = 0, totalClipped = 0;
	int numberOfOperations = aln.CigarData.size();
	int passedFlag = 1;
	for (int i = 0; i < numberOfOperations ; i++)
	{
		const CigarOp& op = aln.CigarData.at(i);
		if (op.Type == 'S')
		{
			if (i == 0)
			{
				headClipped = op.Length;
			}
			if (i == (numberOfOperations-1))
			{
				tailClipped = op.Length;
			}
			totalClipped = totalClipped + op.Length;
		}
	}
	if (tailClipped < singleThreshold && headClipped < singleThreshold && totalClipped < bothThreshold)
	{
		passedFlag = 0;
	}
	return passedFlag;
}

int filterBam(double singleEndSoftclippedThreshold, double bothEndSoftclippedThreshold, int debugging, string inFile, string outFile)
{

	// opening inbam for reading bam
	BamReader reader;
	reader.Open(inFile);
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();	

	// open out bam for writing alignment
	BamWriter writer;
	writer.Open(outFile, header, references);

	//parse alignment
	
	// construct alignment object 
	BamAlignment aln;
	while(reader.GetNextAlignmentCore(aln)){
		int passedFlag = 1, seqLength = aln.Length;
		double singleThreshold = singleEndSoftclippedThreshold * seqLength;
		double bothThreshold = bothEndSoftclippedThreshold * seqLength;
		passedFlag = processAlignment(aln, singleThreshold, bothThreshold);
		
		if ((passedFlag == 0 && debugging ==0) || (passedFlag == 1 && debugging ==1))
		{
			writer.SaveAlignment(aln);
		}
	}
    return 0;
}

int usage(char *program)
{
	cerr << "****************************************************************" << '\n';
	cerr << "Filtering soft clipped reads from paired-end RNA-seq sam files" << '\n';
	cerr << program << "-i <inbam> -o <outbam> -s <oneSideSoftclipFractionThreshold> ";
	cerr << "-b <bothEndSoftclippedThreshold> [-v]" << "\n" << endl;
	cerr << "<inbam>" << "\t" << "Bam file to be filtered (can be: - if stdin) "<< endl;
	cerr << "<outbam>" << "\t" << "Bam file written to (can be: - if stdout) "<< endl;
	cerr << "<oneSideSoftclipFractionThreshold>" << "\t" << "Threshold for filtering one side softclip sequence. "<< endl;
    cerr << "                                  " << "\t" << "Must be between 0 and 1 [default: 0.3]" << endl;
	cerr << "<bothEndSoftclippedThreshold>     " << "\t" << "Threshold for filtering both side softclip sequence. "<< endl;
    cerr << "                                  " << "\t" << "Must be between 0 and 1 [default: 0.4]" << endl;
	cerr << "-v                                " << "\t" << "Debugging mode: print out all failed alignments" << endl;
	cerr << "If the soft clipped bases count > (threshold * [whole sequence length]), it will be filter out" << endl;
	cerr << "****************************************************************\n";
    cerr << endl;
    exit(EXIT_FAILURE);
	return 0;
}

int main(int argc, char **argv)
{
	ios::sync_with_stdio(false);
    int c;
    double singleEndSoftclippedThreshold = 0.3;
    double bothEndSoftclippedThreshold = 0.4;
    int debugging = 0, pairedEndFlag = 0;
    char *program = argv[0];
	string outFile = "stdin", inputFile = "stdin";
    if (argc < 2)
    {
        usage(program);
    }
    while ((c= getopt(argc, argv, "i:s:b:o:v")) != -1)
    {
        switch(c)
        {
            case 's':
                singleEndSoftclippedThreshold = atof(optarg);
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outFile = optarg;
                break;
            case 'b':
                bothEndSoftclippedThreshold = atof(optarg);
                break;
            case 'v':
                debugging = 1;
                break;
            case '?':
                usage(program);
                break;
            default:
                usage(program);
        }
    }
	if (strcmp("-",inputFile.c_str()) == 0)
	{
		inputFile = "stdin";
	}
	if (strcmp("-",outFile.c_str()) == 0)
	{
		outFile = "stdin";
	}

    filterBam(singleEndSoftclippedThreshold,bothEndSoftclippedThreshold,debugging, inputFile, outFile);
    return 0;
}
    
