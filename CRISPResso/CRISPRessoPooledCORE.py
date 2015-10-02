# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:23:47 2015

@author: lpinello
"""


import os
import errno
import sys
import subprocess as sb
import glob
import argparse
import unicodedata
import string
import re

import pandas as pd
import numpy as np
import multiprocessing


import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info



_ROOT = os.path.abspath(os.path.dirname(__file__))

    
####Support functions###
def get_data(path):
        return os.path.join(_ROOT, 'data', path)

GENOME_LOCAL_FOLDER=get_data('genomes')

def force_symlink(src, dst):
    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_',})

def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))

def capitalize_sequence(x):
    return str(x).upper() if not pd.isnull(x) else x

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)

#the dependencies are bowtie2 and samtools
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def check_samtools():

    cmd_path=which('samtools')
    if cmd_path:
        sys.stdout.write('\n samtools is installed! (%s)' %cmd_path)
        return True
    else:
        sys.stdout.write('\nCRISPRessoPooled requires samtools')
        sys.stdout.write('\n\nPlease install it and add to your path following the instruction at: http://www.htslib.org/download/')
        return False

def check_bowtie2():

    cmd_path1=which('bowtie2')
    cmd_path2=which('bowtie2-inspect')

    if cmd_path1 and cmd_path2:
        sys.stdout.write('\n bowtie2 is installed! (%s)' %cmd_path1)
        return True
    else:
        sys.stdout.write('\nCRISPRessoPooled requires Bowtie2!')
        sys.stdout.write('\n\nPlease install it and add to your path following the instruction at: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2')
        return False

#this is overkilling to run for many sequences,
#but for few is fine and effective.
def get_align_sequence(seq,bowtie2_index):
    
    cmd='''bowtie2 -x  %s -c -U %s |\
    grep -v '@' | awk '{OFS="\t"; bpstart=$4; split ($6,a,"[MIDNSHP]"); n=0;  bpend=bpstart;\
    for (i=1; i<=length(a); i++){\
      n+=1+length(a[i]); \
      if (substr($6,n,1)=="S"){\
          bpstart-=a[i];\
          if (bpend==$4)\
            bpend=bpstart;\
      } else if( (substr($6,n,1)!="I")  && (substr($6,n,1)!="H") )\
          bpend+=a[i];\
    }if ( $2 AND 16)print $3,bpstart,bpend,"-",$1,$10,$11;else print $3,bpstart,bpend,"+",$1,$10,$11;}' ''' %(bowtie2_index,seq)
    p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
    return p.communicate()[0]

#if a reference index is provided aligne the reads to it
#extract region
def get_region_from_fa(chr_id,bpstart,bpend,uncompressed_reference):
    region='%s:%d-%d' % (chr_id,bpstart,bpend-1)
    p = sb.Popen("samtools faidx %s %s |   grep -v ^\> | tr -d '\n'" %(uncompressed_reference,region), shell=True,stdout=sb.PIPE)
    return p.communicate()[0]

def get_n_reads_compressed_fastq(compressed_fastq_filename):
     p = sb.Popen("zcat < %s | wc -l" % compressed_fastq_filename , shell=True,stdout=sb.PIPE)
     return float(p.communicate()[0])/4.0

#get a clean name that we can use for a filename
validFilenameChars = "+-_.() %s%s" % (string.ascii_letters, string.digits)

def clean_filename(filename):
    cleanedFilename = unicodedata.normalize('NFKD', unicode(filename)).encode('ASCII', 'ignore')
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)

def get_avg_read_lenght_fastq(fastq_filename):
     cmd=('z' if fastq_filename.endswith('.gz') else '' ) +('cat < %s' % fastq_filename)+\
                  r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
     p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
     return int(p.communicate()[0].strip())
    
    
def find_overlapping_genes(row):
    df_genes_overlapping=df_genes.ix[(df_genes.chrom==row.chr_id) &  
                                     (df_genes.txStart<=row.bpend) &  
                                     (row.bpstart<=df_genes.txEnd)]
    genes_overlapping=[]

    for idx_g,row_g in df_genes_overlapping.iterrows():
        genes_overlapping.append( '%s (%s)' % (row_g.name2,row_g['name']))

    row['gene_overlapping']=','.join(genes_overlapping)

    return row


###EXCEPTIONS############################
class FlashException(Exception):
    pass

class TrimmomaticException(Exception):
    pass

class Bowtie2Exception(Exception):
    pass

class AmpliconsNotUniqueException(Exception):
    pass

class AmpliconsNamesNotUniqueException(Exception):
    pass

class NoReadsAlignedException(Exception):
    pass

class DonorSequenceException(Exception):
    pass

class AmpliconEqualDonorException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass


def main():

    print '  \n~~~CRISPRessoPooled~~~'
    print '-Analysis of CRISPR/Cas9 outcomes from POOLED deep sequencing data-'
    print r'''
          )                                            )
         (           _______________________          (
        __)__       | __  __  __     __ __  |        __)__
     C\|     \      ||__)/  \/  \|  |_ |  \ |     C\|     \
       \     /      ||   \__/\__/|__|__|__/ |       \     /
        \___/       |_______________________|        \___/
    '''


    print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',

    __version__ = re.search(
        '^__version__\s*=\s*"(.*)"',
        open(os.path.join(_ROOT,'CRISPRessoCORE.py')).read(),
        re.M
        ).group(1)
    print 'Version %s\n' % __version__

    parser = argparse.ArgumentParser(description='CRISPRessoPooled Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r1','--fastq_r1', type=str,  help='First fastq file', required=True,default='Fastq filename' )
    parser.add_argument('-r2','--fastq_r2', type=str,  help='Second fastq file for paired end reads',default='')
    parser.add_argument('-f','--amplicons_file', type=str,  help='Amplicons description file', default='')
    parser.add_argument('-x','--bowtie2_index', type=str, help='Basename of Bowtie2 index for the reference genome', default='')

    #tool specific optional
    parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
    please select as table "knowGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')
    parser.add_argument('-p','--n_processes',help='Number of processes to use for the bowtie2 alignment',default=multiprocessing.cpu_count())
    parser.add_argument('--botwie2_options_string', type=str, help='Override options for the Bowtie2 alignment command',default=' -k 1 --end-to-end -N 0 --np 0 ')
    parser.add_argument('--min_perc_reads_to_use_region',  type=int, help='Minimum %% of reads that align to a region to perform the CRISPResso analysis', default=1.0)
    parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=1000)

    #general CRISPResso optional
    parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
    parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
    parser.add_argument('--min_identity_score', type=float, help='Min identity score for the alignment', default=50.0)
    parser.add_argument('-n','--name',  help='Output name', default='')
    parser.add_argument('-o','--output_folder',  help='', default='')
    parser.add_argument('--trim_sequences',help='Enable the trimming of Illumina adapters with Trimmomatic',action='store_true')
    parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic',default=' ILLUMINACLIP:%s:0:90:10:0:true MINLEN:40' % get_data('NexteraPE-PE.fa'))
    parser.add_argument('--min_paired_end_reads_overlap',  type=int, help='Minimum required overlap length between two reads to provide a confident overlap. ', default=4)
    parser.add_argument('-w','--window_around_sgrna', type=int, help='Window(s) in bp around each sgRNA to quantify the indels. Any indels outside this window is excluded. A value of -1 disable this filter.', default=50)
    parser.add_argument('--exclude_bp_from_left', type=int, help='Exclude bp from the left side of the amplicon sequence for the quantification of the indels', default=0)
    parser.add_argument('--exclude_bp_from_right', type=int, help='Exclude bp from the right side of the amplicon sequence for the quantification of the indels', default=0)
    parser.add_argument('--hdr_perfect_alignment_threshold',  type=float, help='Sequence homology %% for an HDR occurrence', default=98.0)
    parser.add_argument('--needle_options_string',type=str,help='Override options for the Needle aligner',default='-gapopen=10 -gapextend=0.5  -awidth3=5000')
    parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
    parser.add_argument('--dump',help='Dump numpy arrays and pandas dataframes to file for debugging purposes',action='store_true')
    parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')


    args = parser.parse_args()


    crispresso_options=['q','s','min_identity_score',
                               'w','exclude_bp_from_left',
                               'exclude_bp_from_right',
                               'hdr_perfect_alignment_threshold',
                              'needle_options_string',
                              'keep_intermediate',
                              'dump',
                              'save_also_png']
    
    def propagate_options(cmd,options,args):
    
        for option in options :
            if option:
                val=eval('args.%s' % option )
                if len(option)==1:
                    cmd+=' -%s %s' % (option,str(val))
                else:
                    if type(val)==str:
                        cmd+=' --%s "%s"' % (option,str(val))
                    else:
                        cmd+=' --%s %s' % (option,str(val))
            
        return cmd


    info('Checking dependencies...')

    if check_samtools() and check_bowtie2():
        print '\n\n All the required dependencies are present!'
    else:
        sys.exit(1)

    #check files
    check_file(args.fastq_r1)
    if args.fastq_r2:
        check_file(args.fastq_r2)

    if args.bowtie2_index:
        check_file(args.bowtie2_index+'.1.bt2')

    if args.amplicons_file:
        check_file(args.amplicons_file)

    if args.gene_annotations:
        check_file(args.gene_annotations)

    if args.amplicons_file and not args.bowtie2_index:
        RUNNING_MODE='ONLY_AMPLICONS'
        info('Only Amplicon description file was provided. The analysis will be perfomed using only the provided amplicons sequences.')

    elif args.bowtie2_index and not args.amplicons_file:
        RUNNING_MODE='ONLY_GENOME'
        info('Only bowtie2 reference genome index file provided. The analysis will be perfomed using only genomic regions where enough reads align.')
    elif args.bowtie2_index and args.amplicons_file:
        RUNNING_MODE='AMPLICONS_AND_GENOME'
        info('Amplicon description file and bowtie2 reference genome index files provided. The analysis will be perfomed using the reads that are aligned ony to the amplicons provided and not to other genomic regions.')
    else:
        error('Please provide the amplicons description file (-t or --amplicons_file option) or the bowtie2 reference genome index file (-x or --bowtie2_index option) or both.')
        sys.exit(1)


    ####TRIMMING AND MERGING
    get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

    if not args.name:
             if args.fastq_r2!='':
                     database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
             else:
                     database_id='%s' % get_name_from_fasta(args.fastq_r1)

    else:
             database_id=args.name
            


    OUTPUT_DIRECTORY='CRISPRessoPOOLED_on_%s' % database_id

    if args.output_folder:
             OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

    _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

    try:
             info('Creating Folder %s' % OUTPUT_DIRECTORY)
             os.makedirs(OUTPUT_DIRECTORY)
             info('Done!')
    except:
             warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

    log_filename=_jp('CRISPRessoPooled_RUNNING_LOG.txt')

    with open(log_filename,'w+') as outfile:
              outfile.write('[Command used]:\nCRISPRessoPooled %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

    if args.fastq_r2=='': #single end reads

         #check if we need to trim
         if not args.trim_sequences:
             #create a symbolic link
             symlink_filename=_jp(os.path.basename(args.fastq_r1))
             force_symlink(os.path.abspath(args.fastq_r1),symlink_filename)
             output_forward_filename=symlink_filename
         else:
             output_forward_filename=_jp('reads.trimmed.fq.gz')
             #Trimming with trimmomatic
             cmd='java -jar %s SE -phred33 %s  %s %s >>%s 2>&1'\
             % (get_data('trimmomatic-0.33.jar'),args.fastq_r1,
                output_forward_filename,
                args.trimmomatic_options_string.replace('NexteraPE-PE.fa','TruSeq3-SE.fa'),
                log_filename)
             #print cmd
             TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)

             if TRIMMOMATIC_STATUS:
                     raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')


         processed_output_filename=output_forward_filename

    else:#paired end reads case

         if not args.trim_sequences:
             output_forward_paired_filename=args.fastq_r1
             output_reverse_paired_filename=args.fastq_r2
         else:
             info('Trimming sequences with Trimmomatic...')
             output_forward_paired_filename=_jp('output_forward_paired.fq.gz')
             output_forward_unpaired_filename=_jp('output_forward_unpaired.fq.gz')
             output_reverse_paired_filename=_jp('output_reverse_paired.fq.gz')
             output_reverse_unpaired_filename=_jp('output_reverse_unpaired.fq.gz')

             #Trimming with trimmomatic
             cmd='java -jar %s PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
             % (get_data('trimmomatic-0.33.jar'),
                     args.fastq_r1,args.fastq_r2,output_forward_paired_filename,
                     output_forward_unpaired_filename,output_reverse_paired_filename,
                     output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
             #print cmd
             TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
             if TRIMMOMATIC_STATUS:
                     raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

             info('Done!')


         #Merging with Flash
         info('Merging paired sequences with Flash...')
         cmd='flash %s %s --min-overlap %d --max-overlap 80  -z -d %s >>%s 2>&1' %\
         (output_forward_paired_filename,
          output_reverse_paired_filename,
          args.min_paired_end_reads_overlap,
          OUTPUT_DIRECTORY,log_filename)

         FLASH_STATUS=sb.call(cmd,shell=True)
         if FLASH_STATUS:
             raise FlashException('Flash failed to run, please check the log file.')

         info('Done!')

         flash_hist_filename=_jp('out.hist')
         flash_histogram_filename=_jp('out.histogram')
         flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
         flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')

         processed_output_filename=_jp('out.extendedFrags.fastq.gz')





        
    #load gene annotation
    if args.gene_annotations:
        print 'Loading gene coordinates from annotation file: %s...' % args.gene_annotations
        try:
            df_genes=pd.read_table(args.gene_annotations,compression='gzip')
            df_genes.txEnd=df_genes.txEnd.astype(int)
            df_genes.txStart=df_genes.txStart.astype(int)
            df_genes.head()
        except:
            print 'Failed to load the gene annotations file.'
    

    if RUNNING_MODE=='ONLY_AMPLICONS' or  RUNNING_MODE=='AMPLICONS_AND_GENOME':

        #load and validate template file
        df_template=pd.read_csv(args.amplicons_file,names=[
                'Name','Amplicon_Sequence','sgRNA',
                'Expected_HDR','Coding_sequence'],comment='#',sep='\t')


        #remove empty amplicons/lines
        df_template.dropna(subset=['Amplicon_Sequence'],inplace=True)
        df_template.dropna(subset=['Name'],inplace=True)

        df_template.Amplicon_Sequence=df_template.Amplicon_Sequence.apply(capitalize_sequence)
        df_template.Expected_HDR=df_template.Expected_HDR.apply(capitalize_sequence)
        df_template.sgRNA=df_template.sgRNA.apply(capitalize_sequence)
        df_template.Coding_sequence=df_template.Coding_sequence.apply(capitalize_sequence)

        if not len(df_template.Amplicon_Sequence.unique())==df_template.shape[0]:
            raise Exception('The amplicons should be all distinct!')

        if not len(df_template.Name.unique())==df_template.shape[0]:
            raise Exception('The amplicon names should be all distinct!')

        df_template=df_template.set_index('Name')

        for idx,row in df_template.iterrows():

            wrong_nt=find_wrong_nt(row.Amplicon_Sequence)
            if wrong_nt:
                 raise NTException('The amplicon sequence %s contains wrong characters:%s' % (row.Name,' '.join(wrong_nt)))

            if not pd.isnull(row.sgRNA):
                wrong_nt=find_wrong_nt(row.sgRNA.strip().upper())
                if wrong_nt:
                    raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % ' '.join(wrong_nt))

                cut_points=[m.start() +len(row.sgRNA)-3 for m in re.finditer(row.sgRNA, row.Amplicon_Sequence)]+[m.start() +2 for m in re.finditer(reverse_complement(row.sgRNA), row.Amplicon_Sequence)]

                if not cut_points:
                    raise SgRNASequenceException('The guide sequence/s provided is(are) not present in the amplicon sequence! \n\nPlease check your input!')


    if RUNNING_MODE=='ONLY_AMPLICONS':
        #create a fasta file with all the amplicons
        amplicon_fa_filename=_jp('AMPLICONS.fa')
        fastq_gz_amplicon_filenames=[]
        with open(amplicon_fa_filename,'w+') as outfile:
            for idx,row in df_template.iterrows():
                if row['Amplicon_Sequence']:
                    outfile.write('>%s\n%s\n' %(clean_filename('AMPL_'+idx),row['Amplicon_Sequence']))

                    #create place-holder fastq files
                    fastq_gz_amplicon_filenames.append(_jp('%s.fastq.gz' % clean_filename('AMPL_'+idx)))
                    open(fastq_gz_amplicon_filenames[-1], 'w+').close()

        df_template['Demultiplexed_fastq.gz_filename']=fastq_gz_amplicon_filenames
        #create a custom index file with all the amplicons
        custom_index_filename=_jp('CUSTOM_BOWTIE2_INDEX')
        sb.call('bowtie2-build %s %s' %(amplicon_fa_filename,custom_index_filename), shell=True)


        #align the file to the amplicons (MODE 1)
        bam_filename_amplicons= _jp('CRISPResso_AMPLICONS_ALIGNED.bam')
        aligner_command= 'bowtie2 -x %s -p %s -k 1 --end-to-end -N 0 --np 0 -U %s | samtools view -bS - > %s' %(custom_index_filename,args.n_processes,args.fastq_r1,bam_filename_amplicons)

        sb.call(aligner_command,shell=True)

        s1=r"samtools view -F 4 %s | grep -v ^'@'" % bam_filename_amplicons
        s2=r'''|awk '{ gzip_filename=sprintf("gzip >> OUTPUTPATH%s.fastq.gz",$3);\
        print "@"$1"\n"$10"\n+\n"$11  | gzip_filename;}' '''

        cmd=s1+s2.replace('OUTPUTPATH',_jp(''))

        print cmd
        sb.call(cmd,shell=True)

        n_reads_aligned_amplicons=[]
        for idx,row in df_template.iterrows():
            n_reads_aligned_amplicons.append(get_n_reads_compressed_fastq(row['Demultiplexed_fastq.gz_filename']))
            crispresso_cmd='CRISPResso -r1 %s -a %s -o %s --name %s' % (row['Demultiplexed_fastq.gz_filename'],row['Amplicon_Sequence'],OUTPUT_DIRECTORY,idx)

            if n_reads_aligned_amplicons[-1]>args.min_reads_to_use_region:
                if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                    crispresso_cmd+=' -g %s' % row['sgRNA']

                if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                    crispresso_cmd+=' -e %s' % row['Expected_HDR']

                if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                    crispresso_cmd+=' -c %s' % row['Coding_sequence']
                
                crispresso_cmd=propagate_options(crispreso_cmd,crispresso_options,args)
                print crispresso_cmd
                sb.call(crispresso_cmd,shell=True)
            else:
                print '\nWARNING: Skipping amplicon [%s] since no reads are aligning to it\n'% idx

        df_template['n_reads']=n_reads_aligned_amplicons
        df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_AMPLICONS.txt'),sep='\t')

    if RUNNING_MODE=='AMPLICONS_AND_GENOME':
        print 'Mapping amplicons to the reference genome...'
        #find the locations of the amplicons on the genome and their strand and check if there are mutations in the reference genome
        additional_columns=[]
        for idx,row in df_template.iterrows():
            fields_to_append=list(np.take(get_align_sequence(row.Amplicon_Sequence, args.bowtie2_index).split('\t'),[0,1,2,3,5]))
            if fields_to_append[0]=='*':
                print 'The amplicon [%s] is not mappable to the reference genome provided!' % idx 
                additional_columns.append([idx,'NOT_ALIGNED',0,-1,'+',''])
            else:
                additional_columns.append([idx]+fields_to_append)
                print 'The amplicon [%s] was mapped to: %s ' % (idx,' '.join(fields_to_append[:3]) )
    
    
        df_template=df_template.join(pd.DataFrame(additional_columns,columns=['Name','chr_id','bpstart','bpend','strand','Reference_Sequence']).set_index('Name'))
        
        df_template.bpstart=df_template.bpstart.astype(int)
        df_template.bpend=df_template.bpend.astype(int)
        
        #Check reference is the same otherwise throw a warning
        for idx,row in df_template.iterrows():
            if row.Amplicon_Sequence != row.Reference_Sequence and row.Amplicon_Sequence != reverse_complement(row.Reference_Sequence):
                print 'Warning the amplicon sequence %s provided:\n%s\n\nis different from the reference sequence(both strand):\n\n%s\n\n%s\n' %(row.name,row.Amplicon_Sequence,row.Amplicon_Sequence,reverse_complement(row.Amplicon_Sequence))


    if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':

        ###HERE we recreate the uncompressed genome file if not available###

        #check you have all the files for the genome and create a fa idx for samtools
        
        uncompressed_reference=args.bowtie2_index+'.fa'
        
        if not os.path.exists(GENOME_LOCAL_FOLDER):
            os.mkdir(GENOME_LOCAL_FOLDER)

        if os.path.exists(uncompressed_reference):
            info('The uncompressed reference fasta file for %s is already present! Skipping generation.' % args.bowtie2_index)
        else:
            #uncompressed_reference=os.path.join(GENOME_LOCAL_FOLDER,'UNCOMPRESSED_REFERENCE_FROM_'+args.bowtie2_index.replace('/','_')+'.fa')
            info('Extracting uncompressed reference from the provided bowtie2 index since it is not available... Please be patient!')

            cmd_to_uncompress='bowtie2-inspect %s > %s' % (args.bowtie2_index,uncompressed_reference)
            print cmd_to_uncompress
            sb.call(cmd_to_uncompress,shell=True)

            info('Indexing fasta file with samtools...')
            #!samtools faidx {uncompressed_reference}
            sb.call('samtools faidx %s' % uncompressed_reference,shell=True)



    #####CORRECT ONE####
    #align in unbiased way the reads to the genome
    if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
        print 'Aligning reads to the provided genome index...'
        bam_filename_genome = _jp('%s_GENOME_ALIGNED.bam' % database_id)
        aligner_command= 'bowtie2 -x %s -p %s -k 1 --end-to-end -N 0 --np 0 -U %s | samtools view -bS - > %s' %(args.bowtie2_index,args.n_processes,processed_output_filename,bam_filename_genome)
        sb.call(aligner_command,shell=True)
        
        #REDISCOVER LOCATIONS and DEMULTIPLEX READS
        MAPPED_REGIONS=_jp('MAPPED_REGIONS/')
        if not os.path.exists(MAPPED_REGIONS):
            os.mkdir(MAPPED_REGIONS)

        s1=r'''samtools view -F 0x0004 %s |\
        awk '{OFS="\t"; bpstart=$4;  bpend=bpstart; split ($6,a,"[MIDNSHP]"); n=0;\
        for (i=1; i<=length(a); i++){\
            n+=1+length(a[i]);\
            if (substr($6,n,1)=="S"){\
                if (bpend==$4)\
                    bpstart-=a[i];\
                else
                    bpend+=a[i];
                }\
            else if( (substr($6,n,1)!="I")  && (substr($6,n,1)!="H") )\
                    bpend+=a[i];\
            }\
            if (and($2, 16))\
                print $3,bpstart,bpend,"-",$1,$10,$11;\
            else\
                print $3,bpstart,bpend,"+",$1,$10,$11;}' | ''' % (bam_filename_genome)
    
        s2=r'''  sort -k1,1 -k2,2n  | awk \
        'BEGIN{chr_id="NA";bpstart=-1;bpend=-1; fastq_filename="NA"}\
        { if ( (chr_id!=$1) || (bpstart!=$2) || (bpend!=$3) )\
            {\
            if (fastq_filename!="NA") {close(fastq_filename); system("gzip "fastq_filename)}\
            chr_id=$1; bpstart=$2; bpend=$3;\
            fastq_filename=sprintf("__OUTPUTPATH__REGION_%s_%s_%s.fastq",$1,$2,$3);\
            }\
        print "@"$5"\n"$6"\n+\n"$7 >> fastq_filename;\
        }' '''
        cmd=s1+s2.replace('__OUTPUTPATH__',MAPPED_REGIONS)
        print cmd
        print sb.call(cmd,shell=True)


    '''
    The most common use case, where many different target sites are pooled into a single 
    high-throughput sequencing library for quantification, is not directly addressed by this implementation. 
    Potential users of CRISPResso would need to write their own code to generate separate input files for processing. 
    Importantly, this preprocessing code would need to remove any PCR amplification artifacts 
    (such as amplification of sequences from a gene and a highly similar pseudogene ) 
    which may confound the interpretation of results. 
    This can be done by mapping of input sequences to a reference genome and removing 
    those that do not map to the expected genomic location, but is non-trivial for an end-user to implement.
    '''
    

    
    if RUNNING_MODE=='AMPLICONS_AND_GENOME':
        files_to_match=glob.glob(_jp('REGION*'))
        n_reads_aligned_genome=[]
        fastq_region_filenames=[]
    
        for idx,row in df_template.iterrows():
    
            #check if we have reads
            fastq_filename_region=_jp('REGION_%s_%s_%s.fastq.gz' % (row['chr_id'],row['bpstart'],row['bpend']))
    
            if os.path.exists(fastq_filename_region):
                
                N_READS=get_n_reads_compressed_fastq(fastq_filename_region)
                n_reads_aligned_genome.append(N_READS)
                fastq_region_filenames.append(fastq_filename_region)
                files_to_match.remove(fastq_filename_region)
                if N_READS>=args.min_reads_to_use_region:
                    print '\nThe amplicon [%s] has enough reads (%d) mapped to it! Running CRISPResso!\n' % (idx,N_READS)
    
                    crispresso_cmd='CRISPResso -r1 %s -a %s -o %s --name %s' % (fastq_filename_region,row['Amplicon_Sequence'],OUTPUT_DIRECTORY,idx)
    
                    if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                        crispresso_cmd+=' -g %s' % row['sgRNA']
    
                    if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                        crispresso_cmd+=' -e %s' % row['Expected_HDR']
    
                    if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                        crispresso_cmd+=' -c %s' % row['Coding_sequence']
                    
                    crispresso_cmd=propagate_options(crispreso_cmd,crispresso_options,args)
                    print crispresso_cmd
                    sb.call(crispresso_cmd,shell=True)
     
                else:
                    print '\nThe amplicon [%s] has not enough reads (%d) mapped to it! Skipping the running of CRISPResso!' % (idx,N_READS)
            else:
                fastq_region_filenames.append('')
                n_reads_aligned_genome.append(0)
                print "The amplicon %s don't have any read mapped to it!\n Please check your amplicon sequence." %  idx
    
        df_template['Amplicon_Specific_fastq.gz_filename']=fastq_region_filenames
        df_template['n_reads']=n_reads_aligned_genome
        
        if args.gene_annotations:
            df_template=df_template.apply(find_overlapping_genes,axis=1)
        
        df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_AND_AMPLICONS.txt'),sep='\t')
        
    
        files_to_match=glob.glob(_jp('REGION*'))  
        #Warn the user if we find reads mapped in other locations
        for fastq_filename_region in files_to_match:
            N_READS=get_n_reads_compressed_fastq(fastq_filename_region)
            if N_READS>=args.min_reads_to_use_region:
                print '\nWARNING: Region %s is not among your amplicons but contains %d reads!' %\
                ((os.path.basename(fastq_filename_region)\
                  .replace('REGION_','').replace('.fastq.gz','').\
                  replace('_',' ')), N_READS)


    if RUNNING_MODE=='ONLY_GENOME' :
        #Load regions and build REFERENCE TABLES 
        info('Parsing the demultiplexed files and extracting locations and reference sequences...')
        coordinates=[]
        for region in glob.glob(os.path.join(MAPPED_REGIONS,'REGION*.fastq.gz')):
            coordinates.append(os.path.basename(region).replace('.fastq.gz','').split('_')[1:4]+[region,get_n_reads_compressed_fastq(region)])
    
        df_regions=pd.DataFrame(coordinates,columns=['chr_id','bpstart','bpend','fastq_file','n_reads'])
        df_regions=df_regions.convert_objects(convert_numeric=True)
        
        df_regions.bpstart=df_regions.bpstart.astype(int)
        df_regions.bpend=df_regions.bpend.astype(int)

        info('Checking overlapping genes...')        
        if args.gene_annotations:
            df_regions=df_regions.apply(find_overlapping_genes,axis=1)
        
        df_regions['sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id,row.bpstart,row.bpend,uncompressed_reference),axis=1)
        df_regions.sort('n_reads',ascending=False,inplace=True)
        df_regions.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_ONLY.txt'),sep='\t',index=None)
        
        
        #run CRISPResso
        #demultiplex reads in the amplicons and call crispresso!
        info('Running CRISPResso on the regions discovered...')
        for idx,row in df_regions.iterrows():
    
            if row.n_reads > args.min_reads_to_use_region:
                info('\nRunning CRISPResso on: %s-%d-%d...'%(row.chr_id,row.bpstart,row.bpend ))
                crispresso_cmd='CRISPResso -r1 %s -a %s -o %s' %(row.fastq_file,row.sequence,OUTPUT_DIRECTORY)  
                crispresso_cmd=propagate_options(crispreso_cmd,crispresso_options,args)
                pcrispresso_cmd
                sb.call(crispresso_cmd,shell=True)
            else:
                info('Skipping region: %s-%d-%d , not enough reads (%d)' %(row.chr_id,row.bpstart,row.bpend, row.n_reads))



    #cleaaning up
    if not args.keep_intermediate:
         info('Removing Intermediate files...')
    
         if args.fastq_r2!='':
             files_to_remove=[processed_output_filename,flash_hist_filename,flash_histogram_filename,\
                          flash_not_combined_1_filename,flash_not_combined_2_filename] 
         else:
             files_to_remove=[processed_output_filename] 
    
         if args.trim_sequences and args.fastq_r2!='':
             files_to_remove+=[output_forward_paired_filename,output_reverse_paired_filename,\
                                               output_forward_unpaired_filename,output_reverse_unpaired_filename]
    
         if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
                 files_to_remove+=[bam_filename_genome]
             
         if RUNNING_MODE=='ONLY_AMPLICONS':  
            files_to_remove+=[bam_filename_amplicons,amplicon_fa_filename]
            for bowtie2_file in glob.glob(_jp('CUSTOM_BOWTIE2_INDEX.*')):
                files_to_remove.append(bowtie2_file)
    
         for file_to_remove in files_to_remove:
             try:
                     if os.path.islink(file_to_remove):
                         print 'LINK',file_to_remove
                         os.unlink(file_to_remove)
                     else:                             
                         os.remove(file_to_remove)
             except:
                     warn('Skipping:%s' %file_to_remove)


       
    info('All Done!')
    print'''     
              )             
             (              
            __)__           
         C\|     \          
           \     /          
            \___/
    '''
    sys.exit(0)


if __name__ == '__main__':
    main()