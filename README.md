# filterSamFile

Softwares to filter sam files
filterSam.cpp and filterSam.py serves to remove alignments with length of soft clipped 

Installation:

``` bash
make
```

This will make:

---
1. filterSoftClipped

```
    Filtering soft clipped reads from paired-end RNA-seq sam files
    usage: cat \<samFile\> \| ./filterSoftClipped -s \<oneSideSoftclipFractionThreshold\> -b \<bothEndSoftclippedThreshold\> \[-vp\]
    
    \<oneSideSoftclipFractionThreshold\\>    oneSideSoftclipFractionThresholdThreshold for filtering one side softclip sequence. Must be between 0 and 1 \[default: 0.3\]
    \<bothEndSoftclippedThreshold\>         bothEndSoftclippedThresholdThreshold for filtering both side softclip sequence. Must be between 0 and 1 \[default: 0.4\]
    -v                                    Debugging mode: print out all failed alignments
    -p                                    paired-end mode \[default = single end\]
    If the soft clipped bases count \> (threshold * [whole sequence length]), the alignment will be filter out
```
---
