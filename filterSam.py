#!/usr/bin/env python

import fileinput
import re
import sys



singleEndSoftclippedThreshold = 0.3
bothEndSoftclippedThreshold = 0.4


def processLine(line):
    id = line.split('\t')[0]
    cigar = line.split('\t')[5]

    #get elements from cigar string
    cigarNumber = re.findall('\d+',cigar)
    cigarLetter = re.findall('[A-Z]',cigar)
    seqLength = sum(int(num) for num in cigarNumber)

    #set threshold
    thresholdSingle = seqLength * singleEndSoftclippedThreshold
    thresholdBoth = seqLength * bothEndSoftclippedThreshold
    totalClip = sum(int(cigarNumber[i]) for i in range(len(cigarNumber)) if cigarLetter[i] == 'S')

    # filtering, if passed, return line, else return 0
    if (cigarLetter[0] == 'S' and \
            ((int(cigarNumber[0])) > thresholdSingle or \
            (cigarLetter[-1] == 'S' and \
            totalClip > thresholdBoth))):
            pass
    elif (cigarLetter[-1] == 'S' and \
            ((int(cigarNumber[-1])) > thresholdSingle or \
            (cigarLetter[0] == 'S' and \
            totalClip > thresholdBoth))):
            pass
    else:
        return line
    return 0

def streamFile():
    for line in fileinput.input():
        # print header
        line = line.strip()
        headerlines = 0
        if line[0] == '@':
            print line
            headerlines += 1
        else:
            #filter each sam lines
            lineNum = fileinput.lineno() - headerlines 
            if lineNum % 2 == 1:
                id1 = line.split('\t')[0]
                lines = [processLine(line)]
            else:
                id2 = line.split('\t')[0]
                assert id1==id2, 'sam input is not sorted!!'
                lines.append(processLine(line))
                if 0 not in lines:
                    print lines[0]
                    print lines[1]
    return 0

def usage(program):
	print '****************************************************************'
	print 'Filtering soft clipped reads from paired-end RNA-seq sam files\n'
	print 'usage: %s <samfile>|<->' %program
	print 'can be stdin when - supplied\n\n'
	print '****************************************************************\n'
	sys.exit()
	return 0

def main():
	if len(sys.argv) != 2:
		program = sys.argv[0]
		usage(program)
	else:
		streamFile()
	return 0

if __name__ =='__main__':
    main()
    
