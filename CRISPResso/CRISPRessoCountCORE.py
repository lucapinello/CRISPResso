# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 13:17:59 2017

@author: Luca Pinello
"""

import os
import gzip
import argparse
import sys
import gzip
import subprocess as sb
from collections import defaultdict
import unicodedata
import re

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

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)
 


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCount!' % library_name)
                sys.exit(1)

def slugify(value): #adapted from the Django project
    
    value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '_', value).strip())
    value = unicode(re.sub('[-\s]+', '-', value))
    
    return str(value)

def filter_se_fastq_by_qual(fastq_filename,output_filename=None,min_bp_quality=20,min_single_bp_quality=0):

        if fastq_filename.endswith('.gz'):
                fastq_handle=gzip.open(fastq_filename)
        else:
                fastq_handle=open(fastq_filename)

        if not output_filename:
                output_filename=fastq_filename.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'

        try: 
            fastq_filtered_outfile=gzip.open(output_filename,'w+')

            for record in SeqIO.parse(fastq_handle, "fastq"):
                if np.array(record.letter_annotations["phred_quality"]).mean()>=min_bp_quality \
                and np.array(record.letter_annotations["phred_quality"]).min()>=min_single_bp_quality:
                    fastq_filtered_outfile.write(record.format('fastq'))
        except:
                raise Exception('Error handling the fastq_filtered_outfile')

 
        return output_filename

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))


def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return int(float(p.communicate()[0])/4.0)


pd=check_library('pandas')
np=check_library('numpy')
Bio=check_library('Bio')
from Bio import SeqIO


###EXCEPTIONS############################
class NTException(Exception):
    pass

class NoReadsAfterQualityFiltering(Exception):
    pass
#########################################

def main():
    try:
        print '  \n~~~CRISPRessoCount~~~'
        print '-Utility to perform sgRNA enumeration from deep sequencing data-'
        print r'''
              )                                             )
             (           ________________________          (
            __)__       | __   __            ___ |        __)__
         C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
           \     /      |\__, \__/ \__/ | \|  |  |       \     /
            \___/       |________________________|        \___/
        '''
        
        
        print'\n[Luca Pinello 2017, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',
        
        
        parser = argparse.ArgumentParser(description='CRISPRessoCount parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-r','--fastq', type=str,  help='fastq file', required=True,default='Fastq filename' )
        
        #optional
        parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
        parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
        parser.add_argument('-t','--tracrRNA',  help="tracr RNA sequence in each read, for single end reads it may necessary to change this parameter if the tracRNA is not fully sequenced, for example to GTTTTAGAG", default='GTTTTAGAGCTAGAAATAGC')
        parser.add_argument('-f','--sgRNA_file', type=str,  help='''sgRNA description file. The format required is one sgRNA per line, for example:AAAAAGATGATTTTTTTCTC\nAAAATATTTTTATCCCCTAA''')
        parser.add_argument('-n','--name',  help='Output name', default='')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('-l','--guide_length',  type=int,help='Lenght in bp to extract the sgRNA upstream of the tracrRNA sequence', default=20)
        parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
        
        args = parser.parse_args()
        
        
        #check files
        check_file(args.fastq)
        if args.sgRNA_file:
            check_file(args.sgRNA_file)
         
        #normalize name and remove not allowed characters
        if args.name:   
            clean_name=slugify(args.name)
            if args.name!= clean_name:
                   warn('The specified name %s contained characters not allowed and was changed to: %s' % (args.name,clean_name))
                   args.name=clean_name
                    
        
        if args.tracrRNA:
            #make evetything uppercase!
            args.tracrRNA=args.tracrRNA.strip().upper()
            wrong_nt=find_wrong_nt(args.tracrRNA)
            if wrong_nt:
                raise NTException('The amplicon sequence contains wrong characters:%s' % ' '.join(wrong_nt))
        
            info('Using specified tracrRNA sequence: %s' % args.tracrRNA)
        
        else:
            info('Using default tracrRNA sequence: %s' % args.tracrRNA)
            
                    
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')
        
        if not args.name:
                database_id='%s' % get_name_from_fasta(args.fastq)
        
        else:
                database_id=args.name
        
        
        OUTPUT_DIRECTORY='CRISPRessoCount_on_%s' % database_id
        
        if args.output_folder:
                OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
        
        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoCount_RUNNING_LOG.txt')
        
        
        try:
        
                 os.makedirs(OUTPUT_DIRECTORY)
                 logging.getLogger().addHandler(logging.FileHandler(log_filename))
        
                 with open(log_filename,'w+') as outfile:
                     outfile.write('[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)


        #filter reads by quality
        if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
            info('Filtering reads with average bp quality < %d and single bp quality < %d ...' % (args.min_average_read_quality,args.min_single_bp_quality))
            processed_output_filename=filter_se_fastq_by_qual(args.fastq,
                        output_filename=_jp(os.path.basename(args.fastq).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'),
                        min_bp_quality=args.min_average_read_quality,
                        min_single_bp_quality=args.min_single_bp_quality)
            info('Done!')
            
            #count reads 
            N_READS_INPUT=get_n_reads_fastq(args.fastq)
            N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
            if N_READS_AFTER_PREPROCESSING == 0:             
                raise NoReadsAfterQualityFiltering('No reads in input or no reads survived the average or single bp quality filtering.')
            else:
                info('Number of reads in input:%d\tNumber of reads after filtering:%d\n' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING))
        else:
            processed_output_filename=args.fastq
 
        if args.sgRNA_file:
            info('Using guides information in %s' % args.sgRNA_file)
            guides_count=dict()
            with open(args.sgRNA_file) as infile:
                for line in infile:
                    guides_count[line.strip()]=0
        else:
            info('No guide information file specified, counting all the guides')
            guides_count=defaultdict(lambda:0)
            
        
        if processed_output_filename.endswith('.gz'):
            infile=gzip.open(processed_output_filename)
        else:
            infile=open(processed_output_filename)
        
            
        info('Counting sgRNAs...')
        N_READS=0
        
        while infile.readline():
            read_seq=infile.readline().strip()
            infile.readline()
            infile.readline()
            
            N_READS+=1
            
            tracrRNA_idx=read_seq.find(args.tracrRNA)
            
            if tracrRNA_idx>=0:
                guide_seq=read_seq[tracrRNA_idx-args.guide_length: tracrRNA_idx]
        
        
                if args.sgRNA_file and not guide_seq in guides_count:
                    pass
                else:
                    guides_count[guide_seq]+=1
                        
        infile.close()
        info('Done!')    

        info('Writing output table...')
        df_guide_counts=pd.Series(guides_count,name='Read_Counts').to_frame()
        df_guide_counts.index.name='Guide_Sequence'
        df_guide_counts['Read_%']=df_guide_counts['Read_Counts']/N_READS*100
        df_guide_counts['RPM']=df_guide_counts['Read_Counts']/N_READS*1000000
        df_guide_counts.head()
        
        df_guide_counts.sort_values(by='Read_Counts',ascending=False).to_csv(_jp('CRISPRessoCount_%s_on_%s.txt' % \
            ('only_ref_guides' if args.sgRNA_file else 'no_ref_guides',args.fastq)),sep='\t')
        
        info('Done!')
        
        if not args.keep_intermediate:
            info('Removing Intermediate files...')
            
            files_to_remove=[]
                             
            if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
                files_to_remove.append(processed_output_filename)
                
            for file_to_remove in files_to_remove:
                         try:
                                 if os.path.islink(file_to_remove):
                                     os.unlink(file_to_remove)
                                 else:                             
                                     os.remove(file_to_remove)
                         except:
                                 warn('Skipping:%s' %file_to_remove)    
    
 
        info('All Done!')
        print r'''
              )                                             )
             (           ________________________          (
            __)__       | __   __            ___ |        __)__
         C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
           \     /      |\__, \__/ \__/ | \|  |  |       \     /
            \___/       |________________________|        \___/
        '''
        sys.exit(0)
    
    except Exception as e:
        error('\n\nERROR: %s' % e)
        sys.exit(-1)

if __name__ == '__main__':
    main()                       
