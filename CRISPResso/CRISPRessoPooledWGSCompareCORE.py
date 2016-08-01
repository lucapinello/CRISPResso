# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 13:53:26 2016

@author: lpinello
"""
import os
import errno
import sys
import subprocess as sb
import glob
import argparse
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


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoPooledWGSCompare!' % library_name)
                sys.exit(1)


def check_output_folder(output_folder):
    quantification_summary_file=os.path.join(output_folder,'SAMPLES_QUANTIFICATION_SUMMARY.txt')  

    if os.path.exists(quantification_summary_file):
        return quantification_summary_file
    else:
        raise OutputFolderIncompleteException('The folder %s  is not a valid CRISPRessoPooled or CRISPRessoWGS output folder.' % output_folder)
        



pd=check_library('pandas')


###EXCEPTIONS############################
class OutputFolderIncompleteException(Exception):
    pass




_ROOT = os.path.abspath(os.path.dirname(__file__))



def main():
    try:
    
        print '  \n~~~CRISPRessoPooledWGSCompare~~~'
        print '-Comparison of two CRISPRessoPooled or CRISPRessoWGS analysis-'
        print r'''
    
    
              )                                                                                     )
             (           ________________________________________________________________          (
            __)__       | __  __  __     __ __        __  __   __ __      __      __  __ |        __)__
         C\|     \      ||__)/  \/  \|  |_ |  \ /|  |/ _ (_   /  /  \|\/||__) /\ |__)|_  |     C\|     \
           \     /      ||   \__/\__/|__|__|__// |/\|\__)__)  \__\__/|  ||   /--\| \ |__ |       \     /
            \___/       |________________________________________________________________|        \___/
        '''
    
        print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',
    
        __version__ = re.search(
            '^__version__\s*=\s*"(.*)"',
            open(os.path.join(_ROOT,'CRISPRessoCORE.py')).read(),
            re.M
            ).group(1)
        print 'Version %s\n' % __version__
    
        parser = argparse.ArgumentParser(description='CRISPRessoPooledWGSCompare Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('crispresso_pooled_wgs_output_folder_1', type=str,  help='First output folder with CRISPRessoPooled or CRISPRessoWGS analysis')
        parser.add_argument('crispresso_pooled_wgs_output_folder_2', type=str,  help='Second output folder with CRISPRessoPooled or CRISPRessoWGS analysis')                           
    
        #OPTIONALS    
        parser.add_argument('-n','--name',  help='Output name', default='')    
        parser.add_argument('-n1','--sample_1_name',  help='Sample 1 name', default='Sample_1')
        parser.add_argument('-n2','--sample_2_name',  help='Sample 2 name', default='Sample_2')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')
    
        args = parser.parse_args()
        
        crispresso_compare_options=['save_also_png',]
    
        
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
    
        #check that the CRISPRessoPooled output is present
        quantification_summary_file_1=check_output_folder(args.crispresso_pooled_wgs_output_folder_1)
        quantification_summary_file_2=check_output_folder(args.crispresso_pooled_wgs_output_folder_2)  
      
        #create outputfolder and initialize the log
        get_name_from_folder=lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPRessoPooled_on_','').replace('CRISPRessoWGS_on_','')
    
        if not args.name:
                 database_id='%s_VS_%s' % (get_name_from_folder(args.crispresso_pooled_wgs_output_folder_1),get_name_from_folder(args.crispresso_pooled_wgs_output_folder_2))
        else:
                 database_id=args.name
        
        
        OUTPUT_DIRECTORY='CRISPRessoPooledWGSCompare_on_%s' % database_id
        
        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
        
        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')
        
        
        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)
        
        log_filename=_jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))
        
        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\nCRISPRessoPooledWGSCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
                  
                  
        #load data and calculate the difference
        df_quant_1=pd.read_table(quantification_summary_file_1)
        df_quant_2=pd.read_table(quantification_summary_file_2)
        df_comp=df_quant_1.set_index('Name').join(df_quant_2.set_index('Name'),lsuffix='_%s' % args.sample_1_name,rsuffix='_%s' % args.sample_2_name)
    
        df_comp['(%s-%s)_Unmodified%%' % (args.sample_1_name,args.sample_2_name)]=df_comp['Unmodified%%_%s' % args.sample_1_name]-df_comp['Unmodified%%_%s' % args.sample_2_name]
        df_comp['(%s-%s)_NHEJ%%' % (args.sample_1_name,args.sample_2_name)]=df_comp['NHEJ%%_%s' % args.sample_1_name]-df_comp['NHEJ%%_%s' % args.sample_2_name]
        df_comp['(%s-%s)_HDR%%' % (args.sample_1_name,args.sample_2_name)]=df_comp['HDR%%_%s' % args.sample_1_name]-df_comp['HDR%%_%s' % args.sample_2_name]
        df_comp['(%s-%s)_Mixed_HDR-NHEJ%%' % (args.sample_1_name,args.sample_2_name)]=df_comp['Mixed_HDR-NHEJ%%_%s' % args.sample_1_name]-df_comp['Mixed_HDR-NHEJ%%_%s' % args.sample_2_name]
    
        df_comp.fillna('NA').to_csv(_jp('COMPARISON_SAMPLES_QUANTIFICATION_SUMMARIES.txt'),sep='\t') 
                                                       
                                                       
        #now run CRISPRessoCompare for the pairs for wich we have data in both folders 
        for idx,row in df_comp.iterrows():
            if row.isnull().any():
                warn('Skipping sample %s since it was not processed in one or both conditions' % idx)
            else:
                crispresso_output_folder_1=os.path.join(args.crispresso_pooled_wgs_output_folder_1,'CRISPResso_on_%s' % idx)
                crispresso_output_folder_2=os.path.join(args.crispresso_pooled_wgs_output_folder_2,'CRISPResso_on_%s' % idx)
                crispresso_compare_cmd='CRISPRessoCompare "%s" "%s" -o "%s" -n1 "%s" -n2 "%s" ' % (crispresso_output_folder_1,
                                                                   crispresso_output_folder_2,
                                                                   OUTPUT_DIRECTORY,
                                                                   args.sample_1_name+'_%s' % idx,
                                                                   args.sample_2_name+'_%s' % idx,
                                                                  )
                
                cmd=propagate_options(crispresso_compare_cmd,crispresso_compare_options,args)
                info('Running CRISPRessoCompare:%s' % crispresso_compare_cmd)
                sb.call(crispresso_compare_cmd,shell=True)
            
            
        info('All Done!')
        print '''     
                      )             
                     (              
                    __)__           
                 C\|     \          
                   \     /          
                    \___/
             '''
        sys.exit(0)
    
    except Exception as e:
        error('\n\nERROR: %s' % e)
        sys.exit(-1)