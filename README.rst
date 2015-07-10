.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso.png?raw=true


CRISPResso is a software pipeline for the analysis of targeted CRISPR-Cas9 deep sequencing data. This algorithm allows for the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences.


CRISPResso automatizes and performs the following steps summarized in the figure below: 

1) filters low quality reads, 
2) trims adapters, 
3) aligns the reads to a reference amplicon, 
4) quantifies the proportion of HDR and NHEJ outcomes, 
5) quantifies frameshift/inframe mutations (if applicable) and identifies affected splice sites,
6) produces a graphical report to visualize and quantify the indels distribution and position.

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_pipeline.png?raw=true

TRY IT ONLINE! 
--------------
If you don't like command line tools you can also use CRISPResso online here:  http://crispresso.rocks


Installation and Requirements
------------
To install the command line version of CRISPresso some dependencies needs to be installed before running the setup:

1) Python 2.7 Anaconda:  http://continuum.io/downloads
2) Java: http://java.com/download
3) C compiler / make. If you have a Mac with a recente version of OSX just open the terminal app and type make, this will install the developer tools for you.

After checking that the required software is installed you can install CRISPResso from the official Python repository following these steps:

1) Open a terminal window
2) Type the command: 

.. code:: bash

        pip install CRISPResso --verbose
      
3) Close the terminal window 

Alternatively if want to install the package without the PIP utility:

1) Download the setup file: https://github.com/lucapinello/CRISPResso/archive/master.zip and decompress it  
2) Open a terminal window  and go to the folder where you have decompressed the zip file
3) Type the command: python setup.py install
4) Close the terminal window and open a new one  (this is important in order to setup correctly the PATH variable in your system).

The Setup will try to install these software for you:

1) Trimmomatic(tested with v0.33): http://www.usadellab.org/cms/?page=trimmomatic
2) Flash(tested with v1.2.11): http://ccb.jhu.edu/software/FLASH/
3) Needle from the EMBOSS suite(tested with 6.6.0): ftp://emboss.open-bio.org/pub/EMBOSS/

If the setup fails on your machine you have to install them manually and put these utilities/binary files in your path!

To check that the installation worked, open a terminal window and execute CRISPResso --help, you should see the help page.

The setup will create automatically a folder in your home folder called CRISPresso_dependencies, don't delete it or it will stop to work! If you want to put the folder in a different location, you need to set the environment variable: CRISPRESSO_DEPENDENCIES_FOLDER. For example to put the folder in /home/lpinello/other_stuff you can write in the terminal *BEFORE* the installation:

.. code:: bash
        
        export CRISPRESSO_DEPENDENCIES_FOLDER=/home/lpinello/other_stuff

Usage
-----
 CRISPResso requires two inputs: (1) single-end reads (single file) or paired-end reads (two files) in FASTQ format (fastq.gz files are also accepted)  from a deep sequencing experiment and (2) a reference amplicon sequence to assess and quantify the efficiency of the targeted mutagenesis. A donor template sequence to assess HDR frequency can be provided as an optional feature. An sgRNA sequence (without PAM sequence) can be provided to compare the predicted cleavage position to the position of the observed mutations. The reads are first filtered based on the quality score (phred33) in order to remove potentially false positive indels. The filtering based on the phred33 quality score can be modulated by adjusting the optimal parameters (see additional notes below). The adapters are trimmed from the reads using Trimmomatic and then sequences are merged with FLASha (if using paired-end data).The remaining reads are then aligned with needle from the EMBOSS suite, an optimal global sequence aligner based on the Needleman-Wunsch algorithm that can easily accounts for gaps. Finally, after analyzing the aligned reads, a set of informative graphs are generated, allowing for the quantification and visualization of the position and type of outcomes within the amplicon sequence.

NHEJ events:

The required inputs are: 

- Two files for paired-end reads or a single file for single-end reads in fastq format (fastq.gz files are also accepted). The reads are assumed to be already trimmed for adapters. If reads are not trimmed, please use the   --trim_sequences option and the   --trimmomatic_options_string  if you are using an adapter different than Nextera. 
- The reference amplicon sequence must also be provided.

Example:

.. code:: bash

                        CRISPResso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a GAATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT

HDR events:
The required inputs are: 

- Two files for paired-end reads or a single file for single-end reads in fastq format (fastq.gz files are also accepted). The reads are assumed to be already trimmed for adapters.
- The reference amplicon sequence.
- The expected amplicon sequence after HDR must also be provided.

Example:

.. code:: bash

                        CRISPResso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGAATGCCGTCACCACCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA -e GCTTACACTTGCTTCTGACACAACTGTGTTCACGAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAAAAAAACGCCGTCACGACGTTATGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGA
                        
IMPORTANT: You must input the entire reference amplicon sequence (’Expected HDR Amplicon sequence’ is the reference for the sequenced amplicon, not simply the donor sequence).  If only the donor sequence is provided, an error will result

Troubleshooting:
----------------

- It is important to check if your reads are trimmed or not. CRISPResso assumes that the reads are already trimmed! If reads are not trimmed, use the option --trim_sequences. The default adapter file used is the Nextera. If you want to specify a custom adapter use the option --trimmomatic_options_string.
- It is possible to use CRISPResso with single end reads. In this case, just omit the option -r2 to specify the second fastq file.
- It is possible to filter based on read quality before aligning reads using the option -q. A reasonable value for this parameter (phred33) is 30.
- The command line CRISPResso tool requires for use on Mac computers requires OS 10.7 or greater. It also requires that command line tools are installed on your machine. After the installation of Anaconda, open the Terminal app and type make, this should prompt you to install command line tools (requires internet connection).
- Once installed, simply typing CRISPResso into any new terminal should load CRISPResso (you will be greeted by the CRISPResso cup)
- Paired end sequencing files requires overlapping sequence from the paired sequencing data
- Use the following command to get to your folder (directory) with sequencing files, assuming that is /home/lpinello/Desktop/CRISPResso_Folder/Sequencing_Files_Folder: cd /home/lpinello/Desktop/CRISPResso_Folder/Sequencing_Files_Folder
- CRISPResso’s default setting is to output analysis files into your directory, otherwise use the --output parameter.

OUTPUT
-----------
The output of CRISPResso consists in of a set of informative graphs is generated, allowingthat allow for the quantification and visualization of where the position and  which types of outcomes are localized inwithin the an amplicon sequence. An example is shown below:

.. image:: https://github.com/lucapinello/CRISPResso/blob/master/CRISPResso_output.png?raw=true


TESTING CRISPResso
------------------

1) Download the two fastq files:

- http://bcb.dfci.harvard.edu/~lpinello/CRISPResso/reads1.fastq.gz 
- http://bcb.dfci.harvard.edu/~lpinello/CRISPResso/reads2.fastq.gz

2) Open a terminal and go to the folder where you have stored the files

3) Type: 

.. code:: bash

                        CRISPResso -r1 reads1.fastq.gz -r2 reads2.fastq.gz -a AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT -g TGAACCAGACCACGGCCCGT 

4) CRISPResso will create a folder with the processed data and the figures.

Useful tips
-----------

- The log of the external utilities called are stored in the file CRISPResso_RUNNING_LOG.txt
- You can specificy the output folder with the option --output_folder 
- You can inspect intermediate files with the option --keep_intermediate
- All the processed raw data used to generate the figures are available in the following plain text files:
        - Quantification_of_editing_frequency.txt, 
        - effect_vector_combined.txt,effect_vector_deletion.txt,
        - effect_vector_insertion.txt,
        - effect_vector_substitution.txt


Parameters of the command line
------------------------------

.. code-block:: bash

  -h, --help            show this help message and exit
  -r1 FASTQ_R1, --fastq_r1 FASTQ_R1
                        First fastq file (default: Fastq filename)
  -r2 FASTQ_R2, --fastq_r2 FASTQ_R2
                        Second fastq file for paired end reads (default: )
  -a AMPLICON_SEQ, --amplicon_seq AMPLICON_SEQ
                        Amplicon Sequence (default: None)
  -g GUIDE_SEQ, --guide_seq GUIDE_SEQ
                        sgRNA sequence, if more than one, please separate by
                        comma/s. Note that the sgRNA needs to be input as the
                        guide RNA sequence (usually 20 nt) immediately 5' of
                        the PAM sequence (usually NGG). If the PAM is found on
                        the opposite strand with respect to the Amplicon
                        Sequence, ensure the sgRNA sequence is also found on
                        the opposite strand. The CRISPResso convention is to
                        depict the expected cleavage position 3 nt 5' of the
                        PAM. (default: )
  -e EXPECTED_HDR_AMPLICON_SEQ, --expected_hdr_amplicon_seq EXPECTED_HDR_AMPLICON_SEQ
                        Amplicon sequence expected after HDR (default: )
  -d DONOR_SEQ, --donor_seq DONOR_SEQ
                        Donor Sequence. This optional input comprises a
                        subsequence of the expected HDR amplicon to be
                        highlighted in plots. (default: )
  -c CODING_SEQ, --coding_seq CODING_SEQ
                        Subsequence/s of the amplicon sequence covering one or
                        more coding sequences for the frameshift analysis.If
                        more than one (for example, split by intron/s), please
                        separate them by comma. (default: )
  -q MIN_AVERAGE_READ_QUALITY, --min_average_read_quality MIN_AVERAGE_READ_QUALITY
                        Minimum average quality score (phred33) to keep a read
                        (default: 0)
  -s MIN_SINGLE_BP_QUALITY, --min_single_bp_quality MIN_SINGLE_BP_QUALITY
                        Minimum single bp score (phred33) to keep a read
                        (default: 0)
  --min_identity_score MIN_IDENTITY_SCORE
                        Min identity score for the alignment (default: 50.0)
  -n NAME, --name NAME  Output name (default: )
  --max_insertion_size MAX_INSERTION_SIZE
                        Max insertion size tolerated for merging paired end
                        reads (default: 60)
  --hdr_perfect_alignment_threshold HDR_PERFECT_ALIGNMENT_THRESHOLD
                        Sequence homology % for an HDR occurrence (default:
                        98.0)
  --trim_sequences      Enable the trimming of Illumina adapters with
                        Trimmomatic (default: False)
  --trimmomatic_options_string TRIMMOMATIC_OPTIONS_STRING
                        Override options for Trimmomatic (default:
                        ILLUMINACLIP:/Users/luca/anaconda/lib/python2.7/site-
                        packages/CRISPResso-0.7.0-py2.7.egg/CRISPResso/data
                        /NexteraPE-PE.fa:0:90:10:0:true MINLEN:40)
  --needle_options_string NEEDLE_OPTIONS_STRING
                        Override options for the Needle aligner (default:
                        -gapopen=10 -gapextend=0.5 -awidth3=5000)
  --keep_intermediate   Keep all the intermediate files (default: False)
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
  --dump                Dump numpy arrays and pandas dataframes to file for
                        debugging purposes (default: False)
  --exclude_bp_from_sides EXCLUDE_BP_FROM_SIDES
                        Exclude bp from each side for the quantification of
                        the indels (default: 0)
  --save_also_png       Save also .png images additionally to .pdf files
                        (default: False)



Acknowledgements
----------------
- Daniel E. Bauer, Matthew C. Canver,Megan D Hoban and Guo-Cheng Yuan contributed to the idea of CRISPResso.
- Daniel E. Bauer, Matthew C., Megan D Hoban, Sorel Fitz-Gibbon and Donald B Kohn, for sharing the data used for the development.
- Many people from Guo-Cheng Yuan for testing CRISPResso.
- Many people from Feng Zhang's lab for the useful feedback and suggestions, in particular David Scott.
- The FAS Research Computing Team for hosting CRISPResso and for the great support, in particular Daniel Kelleher.
