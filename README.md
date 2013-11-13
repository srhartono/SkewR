#SkewR 1.00b
---

## 0. Description

SkewR is a program that uses Hidden Markov Model (HMM) software by Paul Lott (StochHMM) to predict R-loop forming regions in any genome based on trained StochHMM models. The training phase of the HMM is very important for the accuracy of results; Therefore since the HMM models provided is human-trained, the result might be quite accurate for organisms with genome composition similar to human, but might not be very good otherwise. The program will produce R-loop peaks (G Skew or C Skew) and also cluster genes into 4 skew classes (Strong Skew, Weak Skew, No Skew, and Reverse Skew). Please see our publications for more informations.

Human SkewR peak and gene cluster is also available in our website (http://www.mcb.ucdavis.edu/faculty-labs/chedin/Resources.html)

Model files provided are human-trained:  
1. High threshold: GC_SKEW_1mil.hmm  
2. Med threshold: GC_SKEW_20k.hmm  
3. Low threshold: GC_SKEW_7600l.hmm  

## 1. Citation

If you use SkewR or our data for your publication, please cite:

Ginno, P. A., Lim, Y. W., Lott, P. L., Korf, I. & Chédin, F. GC skew at the 5‘ and 3’ ends of human genes links R-loop formation to epigenetic regulation and transcription termination. Genome research 23, 1590–1600 (2013).

Ginno, P. A. P., Lott, P. L. P., Christensen, H. C. H., Korf, I. I. & Chédin, F. F. R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters. Mol Cell 45, 814–825 (2012).

## 2. Installation

`git clone https://github.com/srhartono/SkewR`

Install Bedtools version 2.16 or newer `https://code.google.com/p/bedtools/`  
Install StochHMM version 0.35 or newer `https://github.com/KorfLab/StochHMM`  
Install GitHub in your system if you haven't done so  

## 3. Synopsis

Running SkewR on human genome 19 hg19.fa:

`Perl bin/RunGC-Skew.pl -s hg19.fa -m model/GC_SKEW_7600.hmm -g hg19_gene.bed -b hg19_cpg.bed -o MyResult`

## 4. Usage

usage: Perl RunGC-SKEW.pl [options] -s Sequence -m Model -g Gene File -b CpG file -o Output Directory

Arguments:  
-s: Sequence file (Fasta format)  
-m: Model file (StochHMM HMM format)  
-g: Gene file (BED 6+ format). Bed file containing gene coordinates and strand information (e.g. from UCSC/Refseq)
-b: CpG file (BED 3+ format). Bed file containing CpG island coordinate (e.g. from UCSC)
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

## 5. Contributors

Created SkewR:
Paul Ginno (http://www.fmi.ch/research/groupleader/website/schuebeler/lab_postdocs.php)  
Paul Lott (https://github.com/lottpaul)  

Manages and updates:
Stella Hartono (srhartono@ucdavis.edu)  

## 6. Copyright

Copyright 2012 Chedin Lab

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts. A copy of the license is included in the section entitled "GNU Free Documentation License".

## 7. Disclaimer

This script is provided "as is" without warranty of any kind.
