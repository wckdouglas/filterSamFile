#!/bin/env python

import sys
import numpy as np
import os

if len(sys.argv) != 6:
    sys.exit('usage: python %s <bam file> <threads> <reads> <resultBam> <paired> \n' %(sys.argv[0]))
bamFile = sys.argv[1]
threads = int(sys.argv[2])
reads = int(sys.argv[3])
resultFile = sys.argv[4]
paired = sys.argv[5]
if paired == 'y' or paired == 'yes':
    reads = reads * 2 + 350
else:
    reads = reads + 350
reads = np.random.normal(reads,reads*0.1)

command = 'samtools sort -n -O sam -@ %i -T %s %s ' %(threads,resultFile,bamFile)   +\
        '| head -%i ' %(reads) +\
        '|samtools view -b@ %i - ' %(threads) +\
        '> %s' %(resultFile)
sys.stderr.write('Running: %s\n' %command)
os.system(command)
