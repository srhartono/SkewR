#SkewR 1.00
---

## 1. Installation

`git clone https://github.com/srhartono/SkewR`

Install Bedtools version 2.16 or newer `https://code.google.com/p/bedtools/`  
Install StochHMM version 0.35 or newer `https://github.com/KorfLab/StochHMM`  
Install GitHub in your system if you haven't done so  

## 2. Synopsis

Running SkewR on human genome 19 hg19.fa:

`Perl bin/RunGC-Skew.pl -s hg19.fa -m model/GC_SKEW_7600.hmm -g hg19_gene.bed -b hg19_cpg.bed -o MyResult`

## 3. Usage

usage: Perl RunGC-SKEW.pl [options] -s Sequence -m Model -g Gene File -b CpG file -o Output Directory

Arguments:  
-s: Sequence file (Fasta format)  
-m: Model file (StochHMM HMM format)  
-g: Gene file (BED 6+ format)  
-b: CpG file (BED 3+ format)  
-o: Output Directory (string)  

Options:  

1. StochHMM Parameters:  
-t: Minimum posterior probability threshold (float [0-1]). Default: 0.9  
    Probability threshold is using \"equal or more than\" format (only take probability >= -t)  
-z: Thread number (integer)  

2. Skew Classes Parameters:  
-o: Project name, will be used as output directory (string)  
-x: Intersect with N,N bp of Transcription Start Site (default: -500 +1500 of TSS)  
    Example: -x -500,1500 is -500 +1500 of TSS  
-y: Intersect with Transcription Termintation Site (default: -1500 +500 of TTS)  
    Example: -y -1500,500  
-l: Minimum length of a SkewR peak to be recorded (integer). Default: 300  
Length is using \"equal or more than\" format (only take peak length >= -l)

## 4. Authors

Paul Ginno (http://www.fmi.ch/research/groupleader/website/schuebeler/lab_postdocs.php)  
Paul Lott (https://github.com/lottpaul)  
Stella Hartono (srhartono@ucdavis.edu)  

## 5. Copyright

Copyright 2012 Chedin Lab

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts. A copy of the license is included in the section entitled "GNU Free Documentation License".

## 6. Disclaimer

This script is provided "as is" without warranty of any kind.
