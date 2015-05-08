![alt tag](https://github.com/lucapinello/CRISPResso/blob/master/CRISPRESSO_LOGO.png)



REQUIREMENTS

1) Python 2.7 Anaconda:  continuum.io/downloads
2) Java: http://java.com/download

INSTALLATION

1) Download the setup file and decompress it
2) Run the command: python setup.py install

The Setup will try to install these software for you:

1) Trimmomatic(tested with v0.32): http://www.usadellab.org/cms/?page=trimmomatic
2) Flash(tested with v1.2.11): http://ccb.jhu.edu/software/FLASH/
3) Needle from the EMBOSS suite(tested with 6.6.0): ftp://emboss.open-bio.org/pub/EMBOSS/

If the setup fails on your machine you have to install them manually and put these utilities/binary files in your path!


USING CRISPResso

You can use CRISPResso to quantity the HR or the NHEJ events from PE sequencing experiment.

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

USEFUL TIPS
- The log of the external utilities called are stored in the file CRISPResso_RUNNING_LOG.txt
- If you reads are not trimmed, you can use the option  --trim_sequences (trimmomatic is used in this case)
- Each of the command used: trimmomatic, flash and needle can be fully customize trough the options:
 	--trimmomatic_options_string 
        --flash_options_string FLASH_OPTIONS_STRING
        --needle_options_string NEEDLE_OPTIONS_STRING

- You can specificy the output folder with the option --output_folder 
- You can inspect intermediate files with the option --keep_intermediate
