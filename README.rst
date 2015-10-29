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
To install the command line version of CRISPResso, some dependencies must be installed before running the setup:

1) Python 2.7 Anaconda:  http://continuum.io/downloads
2) Java: http://java.com/download
3) C compiler / make. For Mac with OSX 10.7 or greater, open the terminal app and type and execute the command 'make', which will trigger the installation of OSX developer tools.Windows systems are not officially supported, although CRISPResso may work with Cygwin (https://www.cygwin.com/).

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

The setup will automatically create a folder in your home folder called CRISPResso_dependencies (if this folder is deleted, CRISPResso will not work!)! If you want to put the folder in a different location, you need to set the environment variable: CRISPRESSO_DEPENDENCIES_FOLDER. For example to put the folder in /home/lpinello/other_stuff you can write in the terminal *BEFORE* the installation:

.. code:: bash
        
        export CRISPRESSO_DEPENDENCIES_FOLDER=/home/lpinello/other_stuff

Usage
-----
CRISPResso requires two inputs: (1) paired-end reads (two files) or single-end reads (single file) in .fastq format (fastq.gz files are also accepted) from a deep sequencing experiment and (2) a reference amplicon sequence to assess and quantify the efficiency of the targeted mutagenesis. The amplicon sequence expected after HDR can be provided as an optional input to assess HDR frequency. One or more sgRNA sequences (without PAM sequences) can be provided to compare the predicted cleavage position/s to the position of the observed mutations. Coding sequence/s may be provided to quantify frameshift and potential splice site mutations. 

The reads are first filtered based on the quality score (phred33) in order to remove potentially false positive indels. The filtering based on the phred33 quality score can be modulated by adjusting the optimal parameters (see additional notes below). The adapters are trimmed from the reads using Trimmomatic and then sequences are merged with FLASha (if using paired-end data).The remaining reads are then aligned with needle from the EMBOSS suite, an optimal global sequence aligner based on the Needleman-Wunsch algorithm that can easily accounts for gaps. Finally, after analyzing the aligned reads, a set of informative graphs are generated, allowing for the quantification and visualization of the position and type of outcomes within the amplicon sequence.

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

Understanding the parameters of CRISPResso
------------------------------------------

Required parameters
 To run CRISPResso, only 2 parameters are required for single end reads, or 3 for paired end reads:

-r1 or --fastq_r1: This parameter allows for the specification of the first fastq file.
 
-r2 or  --fastq_r2 FASTQ_R2: This parameter allows for the specification of the second fastq file for paired end reads.

-a or --amplicon_seq: This parameter allows the user to enter the amplicon sequence used for the experiment.

Optional parameters
 In addition to the required parameters explained in the previous section, several optional parameters can be adjusted to tweak your analysis, and to ensure CRISPResso analyzes your data in the best possible way.
 
-g or --guide_seq or: This parameter allows for the specification of the sgRNA sequence. If more than one sequence are included, please separate by comma/s. If the guide RNA sequence is entered, then the position of the guide RNA and the cleavage site will be indicated on the output analysis plots. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20 nt) immediately 5' of the PAM sequence (usually NGG for SpCas9). If the PAM is found on the opposite strand with respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The CRISPResso convention is to depict the expected cleavage position using the value of the parameter cleavage_offset nt 3' from the end of the guide. In addition, the use of alternate nucleases to SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the cleavage_offset parameter to 1, since the default setting of -3 is suitable only for SpCas9. (default:None)

-e or --expected_hdr_amplicon_seq: This parameter allows for the specification of the amplicon sequence expected after HDR. If the data to be analyzed were derived from an experiment using a donor repair template for homology-directed repair (HDR for short), then you have the option to input the sequence of the expected HDR amplicon. This sequence is necessary for CRISPResso to be able to identify successful HDR events within the sequencing data.
 
--hdr_perfect_alignment_threshold: Sequence homology percentage for an HDR occurrence (default: 98.0). This parameter allows for the user to set a threshold for sequence homology for CRISPResso to count instances of successful HDR. This is useful to improve the analysis allowing some tolerance for technical artifacts present in the sequencing data such as sequencing errors or single nucleotide polymorphisms (SNPs) in the cells used in the experiment. Therefore, if you have a read that exhibits successful HDR but has a SNP or sequencing error within the amplicon, you can lower the sequence homology in order allow CRISPResso to count the read as a successful HDR event. If the data are completely free of sequencing errors or polymorphisms, then consider to set parameter to 100.

-d or -donor_seq:This parameter allows the user to highlight the critical subsequence of the expected HDR amplicon in plots. This parameter does not have any effect on the quantification of HDR events.
 
-c, --coding_seq:This parameter allows for the specification of the subsequence/s of the amplicon sequence covering one or more coding sequences for the frameshift analysis. If more than one (for example, split by intron/s), please separate by comma. (default: None)

-q, or --min_average_read_quality: This parameter allows for the specification of the minimum average quality score (phred33) to include a read for the analysis.(default: 0, minimum: 0, maximum: 40). This parameter is helpful to filter out low quality reads. If filtering based on average base quality is desired, a reasonable value for this parameter is greater than 30.

-s or --min_single_bp_quality: This parameter allows for the specification of the minimum single bp score (phred33) to include a read  for the analysis (default: 0, minimum: 0, maximum: 40). This parameter is helpful to filter out low quality reads. This filtering is more aggressive, since any read with a single bp below the threshold will be discarded. If you want to filter your reads based on single base quality to have very high quality reads, a reasonable value for this parameter is greater than 20.

--min_identity_score: This parameter allows for the specification of the min identity score for the alignment (default: 60.0). In order for a read to be considered properly aligned, it should pass this threshold. We suggest to lower this threshold only if really large insertions or deletions are expected in the experiment (>40% of the amplicon length).

-n or --name: This parameter allows for the specification of the output name of the report (default: the names is obtained from the filename of the fastq file/s used in input).

-o or --output_folder: This parameter allows for the specification of the output folder to use for the analysis (default: current folder).
 
--trim_sequences: This parameter enables the trimming of Illumina adapters with Trimmomatic (default: False)

--trimmomatic_options_string: This parameter allows the user the ability to override options for Trimmomatic (default: ILLUMINACLIP:/Users/luca/anaconda/lib/python2.7/site-packages/CRISPResso-0.8.0-py2.7.egg/CRISPResso/data/NexteraPE-PE.fa:0:90:10:0:true). This parameter is useful to specify different adaptor sequences used in the experiment if you need to trim them.

--min_paired_end_reads_overlap: This parameter allows for the specification of the minimum required overlap length between two reads to provide a confident overlap during the merging step. (default: 4, minimum: 1, max: read length)
  
-w ,--window_around_sgrna: This parameter allows for the specification of a window(s) in bp around each sgRNA to quantify the indels. Any indels outside this window are excluded. A value of -1 will disable this filter. (default: -1). This parameter is important since sequencing artifacts and/or SNPs can lead to false positives or false negatives in the quantification of indels and HDR occurrences. Therefore, the user can choose to create a window around the predicted double strand break site of the nuclease used in the experiment. This can help limit non-editing based alterations in an individual read from being inappropriately quantified in CRISPResso analysis.

--cleavage_offset: This parameter allows for the specification of the cleavage offset to use with respect to the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. The default is -3 and is suitable for the SpCas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 set this parameter to 1. (default: -3, minimum:1, max: reference amplicon length). Note: any large indel that partially overlap the window will be also fully quantified.

--exclude_bp_from_left: Exclude bp from the left side of the amplicon sequence for the quantification of the indels (default: 5). This parameter is helpful to avoid artifacts due to imperfect trimming of the reads.

--exclude_bp_from_right: Exclude bp from the right side of the amplicon sequence for the quantification of the indels (default: 5). This parameter is helpful to avoid artifacts due to imperfect trimming of the reads.

--needle_options_string: This parameter allows the user to override options for the Needle aligner (default: -gapopen=10 -gapextend=0.5 -awidth3=5000). More information on the meaning of these parameters can be found in the needle documentation (http://embossgui.sourceforge.net/demo/manual/needle.html). We suggest that only experienced users modify these values.

--keep_intermediate: This parameter allows the user to keep all the intermediate files (default: False). We suggest keeping this parameter disabled for most applications, since the intermediate files (processed reads and alignments) can be really large.

--dump: This parameter allows to dump numpy arrays and pandas dataframes to file for debugging purposes (default: False). 

--save_also_png: This  parameter allows the user to  also save.png images when creating the report., in addition to .pdf files.

Troubleshooting:
----------------

- It is important to check if your reads are trimmed or not. CRISPResso assumes that the reads are already trimmed! If reads are not trimmed, use the option --trim_sequences. The default adapter file used is the Nextera. If you want to specify a custom adapter use the option --trimmomatic_options_string.
- It is possible to use CRISPResso with single end reads. In this case, just omit the option -r2 to specify the second fastq file.
- It is possible to filter based on read quality before aligning reads using the option -q. A reasonable value for this parameter (phred33) is 30.
- The command line CRISPResso tool for use on Mac computers requires OS 10.7 or greater. It also requires that command line tools are installed on your machine. After the installation of Anaconda, open the Terminal app and type make, this should prompt you to install command line tools (requires internet connection).
- Once installed, simply typing CRISPResso into any new terminal should load CRISPResso (you will be greeted by the CRISPResso cup)
- Paired end sequencing files requires overlapping sequence from the paired sequencing data
- Use the following command to get to your folder (directory) with sequencing files, assuming that is /home/lpinello/Desktop/CRISPResso_Folder/Sequencing_Files_Folder: cd /home/lpinello/Desktop/CRISPResso_Folder/Sequencing_Files_Folder
- CRISPResso’s default setting is to output analysis files into your directory, otherwise use the --output parameter.

OUTPUT
-----------
The output of CRISPResso consists of a set of informative graphs that allow for the quantification and visualization of the position and type of outcomes within an amplicon sequence. An example is shown below:

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
- You can specify the output folder with the option --output_folder
- You can inspect intermediate files with the option --keep_intermediate
- All the processed raw data used to generate the figures are available in the following plain text files:
        - Mapping_statistics.txt: this file contains number of: reads in input, reads after preprocessing (merging or quality filtering) and reads properly aligned.
        - Quantification_of_editing_frequency.txt: quantification of editing frequency (number of reads aligned, reads with NHEJ, reads with HDR, and reads with mixed HDR-NHEJ);
        - Frameshift_analysis.txt: number of modified reads with frameshift, in-frame and noncoding mutations;
        - Splice_sites_analysis.txt: number of reads corresponding to potential affected splicing sites;
        - effect_vector_combined.txt: location of mutations (including deletions, insertions, and substitutions) with respect to the reference amplicon;
        - effect_vector_deletion.txt : location of deletions;
        - effect_vector_insertion.txt: location of insertions;
        - effect_vector_substitution.txt: location of substitutions. 
        - position_dependent_vector_avg_insertion_size.txt: average length of the insertions for each position.
        - position_dependent_vector_avg_deletion_size.txt: average length of the deletions for each position.

        



Acknowledgements
----------------
We are grateful to Feng Zhang and David Scott for useful feedback and suggestions; the FAS Research Computing Team, in particular Daniel Kelleher, for great support in hosting the web application of CRISPResso; and Sorel Fitz-Gibbon from UCLA for help in sharing data. Finally, we thank all members of the Guo-Cheng Yuan lab for testing the software.
