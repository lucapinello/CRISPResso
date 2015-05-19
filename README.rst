.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso.png?raw=true


CRISPResso is a command line utility that implements a computational pipeline for the analysis of targeted CRISPR-Cas9 paired end sequence data. 
This algorithm allows the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences. 


CRISPResso automatizes and performs the following steps summarized in the figure below: 
1) filters low quality reads, 
2) trims adapters, 
3) aligns the reads to a reference amplicon, 
4) quantifies the proportion of HDR and NHEJ outcomes, 
5) produces a graphical report to visualize and quantify the indels distribution and location.

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_pipeline.png?raw=true


Requirements
------------
1) Python 2.7 Anaconda:  http://continuum.io/downloads
2) Java: http://java.com/download

Installation
------------

1) Open a terminal window
2) Type the command: pip install CRISPResso --verbose
3) Close the terminal window 

Alternatively if want to install the package without the PIP utility:

1) Download the setup file: https://github.com/lucapinello/CRISPResso/archive/master.zip and decompress it  
2) Open a terminal window
3) Type the command: python setup.py install
4) Close the terminal window 

The Setup will try to install these software for you:

1) Trimmomatic(tested with v0.33): http://www.usadellab.org/cms/?page=trimmomatic
2) Flash(tested with v1.2.11): http://ccb.jhu.edu/software/FLASH/
3) Needle from the EMBOSS suite(tested with 6.6.0): ftp://emboss.open-bio.org/pub/EMBOSS/

If the setup fails on your machine you have to install them manually and put these utilities/binary files in your path!

To check that the installation worked, open a terminal window and execute CRISPResso --help, you should see the help page.

The setup will create automatically a folder in your home folder called CRISPresso_dependencies, don't delete it or it will stop to work!

Usage
-----
CRISPResso requires as input two files for paired-end reads, or a single file for single-end reads, in fastq format (fastq.gz files are also accepted) from a deep sequencing experiment, 
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

- Two files for paired-end reads, or a single file for single-end reads, in fastq format (fastq.gz files are also accepted). The reads are assumed to be already trimmed for adapters, unless an option is specified to trim them
- The reference amplicon seq

Example:

CRISPresso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a GAATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT

HDR events:

In this case the required inputs are:

- Two files for paired-end reads, or a single file for single-end reads, in fastq format (fastq.gz files are also accepted). The reads are assumed to be already trimmed for adapters, unless an option is specified to trim them
- The reference amplicon seq
- The amplicon seq with the donor sequence substituted

Example:
CRISPResso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA -d GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAAAAAAACGCCGTCACGACGTTATGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

NOTES:
-----------

- It is important to check if your reads are trimmed or not. CRISPResso assumes that the reads ARE ALREADY TRIMMED. If not please use the option --trim_sequences. The default adapter file used is the Nextera. If you want to specify a custom adapter use the option  --trimmomatic_options_string. 
- It is possible to use CRISPResso with single end reads, in this case just omit the option -r2 to specify the second fastq file.
- It is possible to filter before the alignment the reads by the average quality using the option --min_bp_quality. A reasonable value for this parameter is 20.

OUTPUT
-----------
The output of CRISPResso consists in a set of informative graphs is generated, allowing the quantification and visualization of where and which types of outcomes are localized in the amplicon sequence
An example is shown for the determination of genome editing outcomes from human erythroid precursors transduced with Cas9 and sgRNA targeting BCL11A exon 2.

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_output.png?raw=true

(A) Frequency distribution of sequence modifications (shown in blue) comprised of insertions, deletions, and substitutions. Reads with unmodified sequence are classified as unmodified (shown in red). (B) Quantification of editing frequency as determined by the percentage and number of sequence reads showing modified and unmodified alleles. (C, left panel) Frequency distribution of sequence modifications that increase read length with respect to the reference amplicon (positive indel size), which are classified as insertions. (C, middle panel) Frequency distribution of sequence modifications that reduce read length (negative indel size) with respect to the reference amplicon, which are classified as deletions. (C, right panel) Frequency distribution of sequence modifications that do not alter read length with respect to the reference amplicon, which are classified as substitutions. (D, left panel) Reads with insertions (red), deletions (purple), and substitution (green) mapped to position on the reference amplicon. The predicted cleavage site by CRISPR/Cas9 is indicated by a vertical dashed line. Only sequence positions directly adjacent to insertions or deletions, or those directly affected by substitution are plotted. (D, right panel)  Frequency distribution of sequence modification comprised of insertions, deletions, and substitutions mapped to position on the reference amplicon.


TESTING CRISPResso
------------------

1) Download the two fastq files:

- http://bcb.dfci.harvard.edu/~lpinello/CRISPResso/reads1.fastq.gz 
- http://bcb.dfci.harvard.edu/~lpinello/CRISPResso/reads2.fastq.gz

2) Open a terminal and go to the folder where you have stored the files

3) Type CRISPResso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT

4) CRISPResso will create a folder with the processed data and the figures.

Useful tips
-----------

- The log of the external utilities called are stored in the file CRISPResso_RUNNING_LOG.txt
- You can specificy the output folder with the option --output_folder 
- You can inspect intermediate files with the option --keep_intermediate
- All the processed raw data used to generate the figures are available in the following plain text files:Quantification_of_editing_frequency.txt, effect_vector_combined.txt,effect_vector_deletion.txt,effect_vector_insertion.txt,effect_vector_substitution.txt


Acknowledgements
------------
- Daniel Bauer, Matthew Canver and Guo-Cheng Yuan contributed to the idea of CRISPResso
- Many people from Feng Zhang lab for the useful feedback and suggestions, in particular David Scott
