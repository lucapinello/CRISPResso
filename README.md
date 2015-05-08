![alt tag](https://github.com/lucapinello/CRISPResso/blob/master/CRISPRESSO_LOGO.png)



CRISPResso
========================

CRISPResso is a command line utility that implements a computational pipeline for the analysis of targeted CRISPR-Cas9 paired end sequence data. 
This algorithm allows the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences. 


CRISPResso automatizes and performs the following steps: 
1) filters low quality reads, 
2) trims adapters, 
3) aligns the reads to a reference amplicon, 
4) quantifies the proportion of HDR and NHEJ outcomes, 
5) produces a graphical report to visualize and quantify the indels distribution and location .

Requirements
------------
1) Python 2.7 Anaconda:  continuum.io/downloads
2) Java: http://java.com/download

Installation
------------

1) Download the setup file and decompress it
2) Run the command: python setup.py install

OR

Alternatively if you have pip installed just run:

pip install CRISPResso --verbose

The Setup will try to install these software for you:

1) Trimmomatic(tested with v0.32): http://www.usadellab.org/cms/?page=trimmomatic
2) Flash(tested with v1.2.11): http://ccb.jhu.edu/software/FLASH/
3) Needle from the EMBOSS suite(tested with 6.6.0): ftp://emboss.open-bio.org/pub/EMBOSS/

If the setup fails on your machine you have to install them manually and put these utilities/binary files in your path!

Usage
-----
CRISPResso requires as input two fastq files with paired end reads from a deep sequencing experiment, 
and a reference amplicon sequence to assess and quantify the efficiency of the targeted mutagenesis; 
optionally a donor template sequence for HDR can be provided and a sgRNA sequence can be provided to compare 
position of predicted cleavage to observed mutations. The reads are first filtered based on the quality score (phred33), 
allowing to remove potentially false positive indels. Subsequently the reads are automatically trimmed for adapters with Trimmomatic 
and  the paired ended sequences are merged with Flash.  The surviving reads are then aligned with needle from the EMBOSS suite, 
an optimal global sequence aligner, based on the Needleman-Wunsch algorithm, that can easily accounts for gaps. Finally, 
after analyzing the aligned reads, a set of informative graphs is generated, allowing the quantification and visualization of 
where and which types of outcomes are localized in the amplicon sequence.


NHEJ events:

In this case the required inputs are:
- two fastq files (pair-ended reads) in fastq format (fastaq.gz files are also accepted), 
- the reference amplicon seq

Example:

CRISPresso reads1.fq reads2.fq GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

HDR events:

In this case the required inputs are:
- two fastq files (pair-ended reads) in fastq format (fastaq.gz files are also accepted), 
- the reference amplicon seq
- the repaired seq

Example:
python CRISPresso.py reads1.fq reads2.fq GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA --repair_seq GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAAAAAAACGCCGTCACGACGTTATGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

Useful tips
-----------

- The log of the external utilities called are stored in the file CRISPResso_RUNNING_LOG.txt
- If you reads are not trimmed, you can use the option  --trim_sequences (trimmomatic is used in this case)
- Each of the command used: trimmomatic, flash and needle can be fully customize trough the options:
 	--trimmomatic_options_string 
        --flash_options_string FLASH_OPTIONS_STRING
        --needle_options_string NEEDLE_OPTIONS_STRING

- You can specificy the output folder with the option --output_folder 
- You can inspect intermediate files with the option --keep_intermediate


