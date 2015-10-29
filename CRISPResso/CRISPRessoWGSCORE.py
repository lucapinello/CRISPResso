# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:23:47 2015

@author: lpinello
"""



import os
import sys
import subprocess as sb
import gzip
import argparse
import unicodedata
import string
import re

import pandas as pd


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

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})

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
        return True
    else:
        sys.stdout.write('\nCRISPRessoPooled requires samtools')
        sys.stdout.write('\n\nPlease install it and add to your path following the instruction at: http://www.htslib.org/download/')
        return False

#if a reference index is provided aligne the reads to it
#extract region
def get_region_from_fa(chr_id,bpstart,bpend,uncompressed_reference):
    region='%s:%d-%d' % (chr_id,bpstart,bpend-1)
    p = sb.Popen("samtools faidx %s %s |   grep -v ^\> | tr -d '\n'" %(uncompressed_reference,region), shell=True,stdout=sb.PIPE)
    return p.communicate()[0]


#get a clean name that we can use for a filename
validFilenameChars = "+-_.() %s%s" % (string.ascii_letters, string.digits)

def clean_filename(filename):
    cleanedFilename = unicodedata.normalize('NFKD', unicode(filename)).encode('ASCII', 'ignore')
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)

    
def find_overlapping_genes(row,df_genes):
    df_genes_overlapping=df_genes.ix[(df_genes.chrom==row.chr_id) &  
                                     (df_genes.txStart<=row.bpend) &  
                                     (row.bpstart<=df_genes.txEnd)]
    genes_overlapping=[]

    for idx_g,row_g in df_genes_overlapping.iterrows():
        genes_overlapping.append( '%s (%s)' % (row_g.name2,row_g['name']))

    row['gene_overlapping']=','.join(genes_overlapping)

    return row


def find_last(mylist,myvalue):
    return len(mylist) - mylist[::-1].index(myvalue) -1

def get_reference_positions( pos, cigar,full_length=True):
    positions = []
    
    ops = re.findall(r'(\d+)(\w)', cigar)
    
    for c in ops:
        l,op=c
        l=int(l)
        
        if op == 'S' or op == 'I':
            if full_length:
                for i in range(0,l):
                    positions.append(None)
        elif op == 'M':
            for i in range(pos,pos+l):
                positions.append(i)
            pos += l
        elif op == 'D' or op == 'N':
            pos += l

    return positions



def write_trimmed_fastq(in_bam_filename,bpstart,bpend,out_fastq_filename):
    p = sb.Popen(
                'samtools view %s | cut -f1,4,6,10,11' % in_bam_filename, 
                stdout = sb.PIPE,
                stderr = sb.STDOUT,
                shell=True
                )

    output=p.communicate()[0]
    n_reads=0
    
    with gzip.open(out_fastq_filename,'w+') as outfile:

        for line in output.split('\n'):
            if line:
                (name,pos,cigar,seq,qual)=line.split()
                #print name,pos,cigar,seq
                pos=int(pos)
                positions=get_reference_positions(pos,cigar)

                if positions[0]<=bpstart and  positions[-1]>=bpend:

                    st=positions.index(bpstart)
                    en=find_last(positions,bpend)
                    #print st,en,seq,seq[st:en]
                    n_reads+=1
                    #print '>%s\n%s\n+\n%s\n' %(name,seq[st:en],qual[st:en])
                    outfile.write('@%s\n%s\n+\n%s\n' %(name,seq[st:en],qual[st:en]))
    return n_reads

###EXCEPTIONS############################

class AmpliconsNamesNotUniqueException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass


def main():

    print '  \n~~~CRISPRessoWGS~~~'
    print '-Analysis of CRISPR/Cas9 outcomes from WGS data-'
    print r'''
       )                                 )
      (           ____________          (
     __)__       |     __  __ |        __)__
  C\|     \      ||  |/ _ (_  |     C\|     \
    \     /      ||/\|\__)__) |       \     /
     \___/       |____________|        \___/
    '''

    print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',

    __version__ = re.search(
        '^__version__\s*=\s*"(.*)"',
        open(os.path.join(_ROOT,'CRISPRessoCORE.py')).read(),
        re.M
        ).group(1)
    print 'Version %s\n' % __version__

    #tool specific optional
    parser = argparse.ArgumentParser(description='CRISPRessoWGS Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b','--bam_file', type=str,  help='WGS aligned bam file', required=True,default='bam filename' )
    parser.add_argument('-f','--region_file', type=str,  help='Regions description file. A BED format  file containing the regions to analyze, one per line. The REQUIRED\
    columns are: chr_id(chromosome name), bpstart(start position), bpend(end position), the optional columns are:name (an unique indentifier for the region), guide_seq, expected_hdr_amplicon_seq,coding_seq, see CRISPResso help for more details on these last 3 parameters)', default='')
    parser.add_argument('-r','--reference_file', type=str, help='A FASTA format reference file (for example hg19.fa for the human genome)', default='')
    parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=10)
    parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
    please select as table "knowGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')

    #general CRISPResso optional
    parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
    parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
    parser.add_argument('--min_identity_score', type=float, help='Min identity score for the alignment', default=60.0)
    parser.add_argument('-n','--name',  help='Output name', default='')
    parser.add_argument('-o','--output_folder',  help='', default='')
    parser.add_argument('--trim_sequences',help='Enable the trimming of Illumina adapters with Trimmomatic',action='store_true')
    parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic',default=' ILLUMINACLIP:%s:0:90:10:0:true MINLEN:40' % get_data('NexteraPE-PE.fa'))
    parser.add_argument('--min_paired_end_reads_overlap',  type=int, help='Minimum required overlap length between two reads to provide a confident overlap. ', default=4)
    parser.add_argument('-w','--window_around_sgrna', type=int, help='Window(s) in bp around the cleavage position (half on on each side) as determined by the provide guide RNA sequence to quantify the indels. Any indels outside this window are excluded. A value of -1 disables this filter.', default=-1)
    parser.add_argument('--cleavage_offset', type=int, help="Cleavage offset to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. The default is -3 and is suitable for the SpCas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1.", default=-3)        
    parser.add_argument('--exclude_bp_from_left', type=int, help='Exclude bp from the left side of the amplicon sequence for the quantification of the indels', default=5)
    parser.add_argument('--exclude_bp_from_right', type=int, help='Exclude bp from the right side of the amplicon sequence for the quantification of the indels', default=5)
    parser.add_argument('--hdr_perfect_alignment_threshold',  type=float, help='Sequence homology %% for an HDR occurrence', default=98.0)
    parser.add_argument('--needle_options_string',type=str,help='Override options for the Needle aligner',default=' -gapopen=10 -gapextend=0.5  -awidth3=5000')
    parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
    parser.add_argument('--dump',help='Dump numpy arrays and pandas dataframes to file for debugging purposes',action='store_true')
    parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')



    args = parser.parse_args()

    crispresso_options=['window_around_sgrna','cleavage_offset','min_average_read_quality','min_single_bp_quality','min_identity_score',
                               'min_single_bp_quality','exclude_bp_from_left',
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
  
                if type(val)==str:
                    cmd+=' --%s "%s"' % (option,str(val)) # this is for options with space like needle...
                elif type(val)==bool:
                    if val:
                        cmd+=' --%s' % option
                else:
                    cmd+=' --%s %s' % (option,str(val))
            
        return cmd


        
    info('Checking dependencies...')

    if check_samtools():
        info('\n All the required dependencies are present!')
    else:
        sys.exit(1)

    #check files
    check_file(args.bam_file)
    
    check_file(args.reference_file)

    if args.region_file:
        check_file(args.region_file)

    if args.gene_annotations:
        check_file(args.gene_annotations)
        

    #INIT
    get_name_from_bam=lambda  x: os.path.basename(x).replace('.bam','')

    if not args.name:
        database_id='%s' % get_name_from_bam(args.bam_file)
    else:
        database_id=args.name
            

    OUTPUT_DIRECTORY='CRISPRessoWGS_on_%s' % database_id

    if args.output_folder:
             OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

    _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

    try:
             info('Creating Folder %s' % OUTPUT_DIRECTORY)
             os.makedirs(OUTPUT_DIRECTORY)
             info('Done!')
    except:
             warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

    log_filename=_jp('CRISPRessoWGS_RUNNING_LOG.txt')
    logging.getLogger().addHandler(logging.FileHandler(log_filename))

    with open(log_filename,'w+') as outfile:
              outfile.write('[Command used]:\nCRISPRessoWGS %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
            
    #check if bam has the index already
    if os.path.exists(args.bam_file+'.bai'):
        info('Index file for input .bam file exists, skipping generation.')
    else:
        info('Creating index file for input .bam file...')
        sb.call('samtools index %s ' % (args.bam_file),shell=True)
        
        
        
    #load gene annotation
    if args.gene_annotations:
        info('Loading gene coordinates from annotation file: %s...' % args.gene_annotations)
        try:
            df_genes=pd.read_table(args.gene_annotations,compression='gzip')
            df_genes.txEnd=df_genes.txEnd.astype(int)
            df_genes.txStart=df_genes.txStart.astype(int)
            df_genes.head()
        except:
            info('Failed to load the gene annotations file.')
            
        
    #Load and validate the REGION FILE
    df_regions=pd.read_csv(args.region_file,names=[
            'chr_id','bpstart','bpend','Name','sgRNA',
            'Expected_HDR','Coding_sequence'],comment='#',sep='\t')


    #remove empty amplicons/lines
    df_regions.dropna(subset=['chr_id','bpstart','bpend'],inplace=True)
    
    df_regions.Expected_HDR=df_regions.Expected_HDR.apply(capitalize_sequence)
    df_regions.sgRNA=df_regions.sgRNA.apply(capitalize_sequence)
    df_regions.Coding_sequence=df_regions.Coding_sequence.apply(capitalize_sequence)

        
    #check or create names
    for idx,row in df_regions.iterrows():
        if pd.isnull(row.Name):
            df_regions.ix[idx,'Name']='_'.join(map(str,[row['chr_id'],row['bpstart'],row['bpend']]))


    if not len(df_regions.Name.unique())==df_regions.shape[0]:
        raise Exception('The amplicon names should be all distinct!')

    df_regions=df_regions.set_index('Name')
    df_regions.index=df_regions.index.str.replace(' ','_')


    #extract sequence for each region
    uncompressed_reference=args.reference_file
    
    if os.path.exists(uncompressed_reference+'.fai'):
        info('The index for the reference fasta file is already present! Skipping generation.')
    else:
        info('Indexing reference file... Please be patient!')
        sb.call('samtools faidx %s >>%s 2>&1' % (uncompressed_reference,log_filename),shell=True)
    
    df_regions['sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id,row.bpstart,row.bpend,uncompressed_reference),axis=1)

    for idx,row in df_regions.iterrows():

        if not pd.isnull(row.sgRNA):
            wrong_nt=find_wrong_nt(row.sgRNA.strip().upper())
            if wrong_nt:
                raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % ' '.join(wrong_nt))

            offset_fw=args.cleavage_offset+len(row.sgRNA)-1
            offset_rc=(-args.cleavage_offset)-1
            cut_points=[m.start() + offset_fw for m in re.finditer(row.sgRNA,  row.sequence)]+[m.start() + offset_rc for m in re.finditer(reverse_complement(row.sgRNA),  row.sequence)]

            #cut_points=[m.start() +len(row.sgRNA)-3 for m in re.finditer(row.sgRNA, row.sequence)]+[m.start() +2 for m in re.finditer(reverse_complement(row.sgRNA), row.sequence)]

            if not cut_points:
                #raise SgRNASequenceException('The guide sequence/s provided is(are) not present in the amplicon sequence! \n\nPlease check your input!')
                df_regions.ix[idx,'sgRNA']=''
    
    df_regions=df_regions.convert_objects(convert_numeric=True)
    
    df_regions.bpstart=df_regions.bpstart.astype(int)
    df_regions.bpend=df_regions.bpend.astype(int)       
    
    if args.gene_annotations:
        df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)
        

    #extract reads with samtools in that region and create a bam
    #create a fasta file with all the trimmed reads
    info('\nProcessing each regions...')
    
    ANALYZED_REGIONS=_jp('ANALYZED_REGIONS/')
    if not os.path.exists(ANALYZED_REGIONS):
        os.mkdir(ANALYZED_REGIONS)
    
    df_regions['n_reads']=0
    df_regions['bam_file_with_reads_in_region']=''
    df_regions['fastq.gz_file_trimmed_reads_in_region']=''
    
    for idx,row in df_regions.iterrows():
        
        if row['sequence']:
            
            fastq_gz_filename=os.path.join(ANALYZED_REGIONS,'%s.fastq.gz' % clean_filename('REGION_'+idx))
            bam_region_filename=os.path.join(ANALYZED_REGIONS,'%s.bam' % clean_filename('REGION_'+idx))
    
            #create place-holder fastq files
            open(fastq_gz_filename, 'w+').close()
    
            region='%s:%d-%d' % (row.chr_id,row.bpstart,row.bpend-1)
            info('\nExtracting reads in:%s and create the .bam file: %s' % (region,bam_region_filename))
            
            #extract reads in region
            cmd=r'''samtools view -b -F 4 %s %s > %s ''' % (args.bam_file, region, bam_region_filename)
            #print cmd
            sb.call(cmd,shell=True)
    
    
            #index bam file
            cmd=r'''samtools index %s ''' % (bam_region_filename)
            #print cmd
            sb.call(cmd,shell=True)
            
            info('Trim reads and create a fastq.gz file in: %s' % fastq_gz_filename)
            #trim reads in bam and convert in fastq
            n_reads=write_trimmed_fastq(bam_region_filename,row['bpstart'],row['bpend'],fastq_gz_filename)
            df_regions.ix[idx,'n_reads']=n_reads
            df_regions.ix[idx,'bam_file_with_reads_in_region']=bam_region_filename
            df_regions.ix[idx,'fastq.gz_file_trimmed_reads_in_region']=fastq_gz_filename
            
            
    df_regions.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_SELECTED_REGIONS_WGS.txt'),sep='\t')              

    #Run Crispresso    
    info('\nRunning CRISPResso on each regions...')
    for idx,row in df_regions.iterrows():
    
           if row['n_reads']>=args.min_reads_to_use_region:
                info('\nThe region [%s] has enough reads (%d) mapped to it!' % (idx,row['n_reads']))
    
                crispresso_cmd='CRISPResso -r1 %s -a %s -o %s --name %s' %\
                (row['fastq.gz_file_trimmed_reads_in_region'],row['sequence'],OUTPUT_DIRECTORY,idx)
    
                if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                    crispresso_cmd+=' -g %s' % row['sgRNA']
    
                if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                    crispresso_cmd+=' -e %s' % row['Expected_HDR']
    
                if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                    crispresso_cmd+=' -c %s' % row['Coding_sequence']
    
                crispresso_cmd=propagate_options(crispresso_cmd,crispresso_options,args)
                info('Running CRISPResso:%s' % crispresso_cmd)
                sb.call(crispresso_cmd,shell=True)
    
           else:
                info('\nThe region [%s] has not enough reads (%d) mapped to it! Skipping the running of CRISPResso!' % (idx,row['n_reads']))
        

    info('All Done!')
    print r'''
       )                                 )
      (           ____________          (
     __)__       |     __  __ |        __)__
  C\|     \      ||  |/ _ (_  |     C\|     \
    \     /      ||/\|\__)__) |       \     /
     \___/       |____________|        \___/
    '''
  
    sys.exit(0)


if __name__ == '__main__':
    main()