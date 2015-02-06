![alt tag](https://github.com/lucapinello/CRISPResso/blob/master/CRISPRESSO_LOGO.png)



REQUIREMENTS

1) Python 2.7: http://www.python.org
2) Python Libraries: 
	- pandas: http://pandas.pydata.org/, 
	- numpy: http://www.numpy.org/, 
	- matplotlib: http://matplotlib.org/


3) Java: http://java.com/download
4) Trimmomatic(tested with v0.32): http://www.usadellab.org/cms/?page=trimmomatic
5) Flash(tested with v1.2.11): http://ccb.jhu.edu/software/FLASH/
6) Needle from the EMBOSS suite(tested with 6.6.0): ftp://emboss.open-bio.org/pub/EMBOSS/

These utilities/binary files should be in your path!

INSTALLATION

1) Install all the requirements
2) Decompress the file Crispresso.tar.gz
3) Add the obtained folder to your PATH env variable

USING CRISPResso

You can use CRISPResso to quantity the HR or the NHEJ events from PE sequencing experiment.

NHEJ events:

In this case the required inputs are:
- two fastq files (pair-ended reads) in fastq format (fastaq.gz files are also accepted), 
- the reference amplicon seq

Example:

python CRISPresso.py reads1.fq reads2.fq GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA

HR events:

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
