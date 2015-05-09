.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso.png


CRISPResso is a command line utility that implements a computational pipeline for the analysis of targeted CRISPR-Cas9 paired end sequence data. 
This algorithm allows the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences. 


CRISPResso automatizes and performs the following steps summarized in the figure below: 
1) filters low quality reads, 
2) trims adapters, 
3) aligns the reads to a reference amplicon, 
4) quantifies the proportion of HDR and NHEJ outcomes, 
5) produces a graphical report to visualize and quantify the indels distribution and location.

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_pipeline.png


Requirements
------------
1) Python 2.7 Anaconda:  continuum.io/downloads
2) Java: http://java.com/download

Installation
------------

1) Download the setup file: https://github.com/lucapinello/CRISPResso/archive/master.zip and decompress it  
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

CRISPresso -r1 reads1.fq -r2 reads2.fq -a GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

HDR events:

In this case the required inputs are:
- two fastq files (paired end reads) in fastq format (fastq.gz files are also accepted), 
- the reference amplicon seq
- the amplicon seq with the donor sequence substituted

Example:
CRISPResso -r1 reads1.fq -r 2 reads2.fq -a GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA -d GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAAAAAAACGCCGTCACGACGTTATGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

OUTPUT
-----------
The output of CRISPResso consists in a set of informative graphs is generated, allowing the quantification and visualization of where and which types of outcomes are localized in the amplicon sequence
An example is shown for the determination of genome editing outcomes from human erythroid precursors transduced with Cas9 and sgRNA targeting BCL11A exon 2.

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_output.png

(A) Frequency distribution of sequence modifications (shown in blue) comprised of insertions, deletions, and substitutions. Reads with unmodified sequence are classified as unmodified (shown in red). (B) Quantification of editing frequency as determined by the percentage and number of sequence reads showing modified and unmodified alleles. (C, left panel) Frequency distribution of sequence modifications that increase read length with respect to the reference amplicon (positive indel size), which are classified as insertions. (C, middle panel) Frequency distribution of sequence modifications that reduce read length (negative indel size) with respect to the reference amplicon, which are classified as deletions. (C, right panel) Frequency distribution of sequence modifications that do not alter read length with respect to the reference amplicon, which are classified as substitutions. (D, left panel) Reads with insertions (red), deletions (purple), and substitution (green) mapped to position on the reference amplicon. The predicted cleavage site by CRISPR/Cas9 is indicated by a vertical dashed line. Only sequence positions directly adjacent to insertions or deletions, or those directly affected by substitution are plotted. (D, right panel)  Frequency distribution of sequence modification comprised of insertions, deletions, and substitutions mapped to position on the reference amplicon.


Useful tips
-----------

- The log of the external utilities called are stored in the file CRISPResso_RUNNING_LOG.txt
- If you reads are not trimmed, you can use the option  --trim_sequences (trimmomatic is used in this case)
- Each of the command used: trimmomatic, flash and needle can be fully customize trough the options:
 	  
 --trimmomatic_options_string 
 --flash_options_string 
 --needle_options_string 

- You can specificy the output folder with the option --output_folder 
- You can inspect intermediate files with the option --keep_intermediate


