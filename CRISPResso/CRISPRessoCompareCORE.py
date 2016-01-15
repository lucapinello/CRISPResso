# -*- coding: utf-8 -*-
import os
import errno
import sys
import argparse
import re
import cPickle as cp


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
                error('You need to install %s module to use CRISPRessoCompare!' % library_name)
                sys.exit(1)


def check_output_folder(output_folder):
    quantification_file=os.path.join(output_folder,'Quantification_of_editing_frequency.txt')  
    profile_file=os.path.join(output_folder,'effect_vector_combined.txt') 
    if os.path.exists(quantification_file) and profile_file:
        return quantification_file,profile_file
    else:
        raise OutputFolderIncompleteException('The folder %s  is not a valid CRISPResso output folder.' % output_folder)


def check_hdr_mode(output_folder_1,output_folder_2):
    hdr_1=os.path.exists(os.path.join(output_folder_1,'effect_vector_insertion_HDR.txt'))
    hdr_2=os.path.exists(os.path.join(output_folder_2,'effect_vector_insertion_HDR.txt'))
    
    if not hdr_1 ^ hdr_2:
        if hdr_1:
            return True
        else:
            return False 
    else:
        raise  MixedRunningModeException('You cannot mix outputs with different running modes (HDR/NHEJ/MIXED with only NHEJ')
    
    
    

def parse_quantification(quantification_file):
    with open(quantification_file) as infile:
        infile.readline()
        N_UNMODIFIED=float(re.findall("Unmodified:(\d+)",infile.readline())[0])
        N_MODIFIED=float(re.findall("NHEJ:(\d+)",infile.readline())[0])
        N_REPAIRED=float(re.findall("HDR:(\d+)", infile.readline())[0])
        N_MIXED_HDR_NHEJ=float(re.findall("Mixed HDR-NHEJ:(\d+)", infile.readline())[0])
        infile.readline()
        N_TOTAL=float(re.findall("Total Aligned:(\d+) reads",infile.readline())[0])
        return N_UNMODIFIED,N_MODIFIED,N_REPAIRED,N_MIXED_HDR_NHEJ,N_TOTAL


def parse_profile(profile_file):
    return np.loadtxt(profile_file,skiprows=1)


def load_cut_points_sgRNA_intervals(output_folder):
    cut_points_file=os.path.join(output_folder,'cut_points.pickle')
    sgRNA_intervals_file=os.path.join(output_folder,'sgRNA_intervals.pickle')
    if os.path.exists(cut_points_file):
        cut_points=cp.load(open(cut_points_file))
    else:
        cut_points=[]

        
    if os.path.exists(sgRNA_intervals_file):
        sgRNA_intervals=cp.load(open(sgRNA_intervals_file))
    else:
        sgRNA_intervals=[]
    
    return  cut_points,sgRNA_intervals
  
 
###EXCEPTIONS############################
class OutputFolderIncompleteException(Exception):
    pass

class MixedRunningModeException(Exception):
    pass

class DifferentAmpliconLengthException(Exception):
    pass
############################
    
    

matplotlib=check_library('matplotlib')
from matplotlib import font_manager as fm
font = {'size'   : 20}
matplotlib.rc('font', **font)
matplotlib.use('Agg')

plt=check_library('pylab')
np=check_library('numpy') 


_ROOT = os.path.abspath(os.path.dirname(__file__))



def main():

    try:
        print '  \n~~~CRISPRessoCompare~~~'
        print '-Comparison of two CRISPResso analysis-'
        print r'''
    
    
              )                                                )
             (           ___________________________          (
            __)__       | __ __      __      __  __ |        __)__
         C\|     \      |/  /  \|\/||__) /\ |__)|_  |     C\|     \
           \     /      |\__\__/|  ||   /--\| \ |__ |       \     /
            \___/       |___________________________|        \___/
        '''
    
        print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',
    
        __version__ = re.search(
            '^__version__\s*=\s*"(.*)"',
            open(os.path.join(_ROOT,'CRISPRessoCORE.py')).read(),
            re.M
            ).group(1)
        print 'Version %s\n' % __version__
    
        parser = argparse.ArgumentParser(description='CRISPRessoCompare Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('crispresso_output_folder_1', type=str,  help='First output folder with CRISPResso analysis')
        parser.add_argument('crispresso_output_folder_2', type=str,  help='Second output folder with CRISPResso analysis')                           
    
        #OPTIONALS    
        parser.add_argument('-n','--name',  help='Output name', default='')    
        parser.add_argument('-n1','--sample_1_name',  help='Sample 1 name', default='Sample_1')
        parser.add_argument('-n2','--sample_2_name',  help='Sample 2 name', default='Sample_2')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')
    
        args = parser.parse_args()
        
    
    
    
        
        #check that the CRISPResso output is present
        quantification_file_1,profile_file_1=check_output_folder(args.crispresso_output_folder_1)
        quantification_file_2,profile_file_2=check_output_folder(args.crispresso_output_folder_2)
        
        #check we are not mixing modes
        check_hdr_mode(args.crispresso_output_folder_1,args.crispresso_output_folder_2)
        
        
        get_name_from_folder=lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPResso_on_','')
    
        if not args.name:
                 database_id='%s_VS_%s' % (get_name_from_folder(args.crispresso_output_folder_1),get_name_from_folder(args.crispresso_output_folder_2))
        else:
                 database_id=args.name
        
        
        OUTPUT_DIRECTORY='CRISPRessoCompare_on_%s' % database_id
        
        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
        
        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoCompare_RUNNING_LOG.txt')
        
        
        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)
        
        log_filename=_jp('CRISPRessoCompare_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))
        
        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\nCRISPRessoCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
    
       
        #LOAD DATA
        N_UNMODIFIED_1,N_MODIFIED_1,N_REPAIRED_1,N_MIXED_HDR_NHEJ_1,N_TOTAL_1=parse_quantification(quantification_file_1)
        N_UNMODIFIED_2,N_MODIFIED_2,N_REPAIRED_2,N_MIXED_HDR_NHEJ_2,N_TOTAL_2=parse_quantification(quantification_file_2)
        
        profile_1=parse_profile(profile_file_1)
        profile_2=parse_profile(profile_file_2)
        
        try:
            assert np.all(profile_1[:,0]==profile_2[:,0])
        except:
            pass
            #raise DifferentAmpliconLengthException('Different amplicon lenghts for the two amplicons.')
        len_amplicon=profile_1.shape[0]
        effect_vector_any_1=profile_1[:,1]
        effect_vector_any_2=profile_2[:,1]
        cut_points,sgRNA_intervals=load_cut_points_sgRNA_intervals(args.crispresso_output_folder_1)
        
        
        #Quantification comparison barchart
        fig=plt.figure(figsize=(30,15))    
        n_groups = 4
        
        means_sample_1= np.array([N_UNMODIFIED_1,N_MODIFIED_1,N_REPAIRED_1,N_MIXED_HDR_NHEJ_1,])/N_TOTAL_1*100
        means_sample_2 = np.array([N_UNMODIFIED_2,N_MODIFIED_2,N_REPAIRED_2,N_MIXED_HDR_NHEJ_2])/N_TOTAL_2*100
        
        ax1=fig.add_subplot(1,2,1)
            
        index = np.arange(n_groups)
        bar_width = 0.35
        
        opacity = 0.4
        error_config = {'ecolor': '0.3'}
        
        ax1.bar(index, means_sample_1, bar_width,
                         alpha=opacity,
                         color=(0,0,1,0.4),
                         label=args.sample_1_name)
        
        ax1.bar(index + bar_width, means_sample_2, bar_width,
                         alpha=opacity,
                         color=(1,0,0,0.4),
                         label=args.sample_2_name)
        
        plt.ylabel('% Sequences')
        plt.title('%s VS %s' % (args.sample_1_name,args.sample_2_name))
        plt.xticks(index + bar_width, ('Unmodified', 'NHEJ', 'HDR', 'Mixed'))
        plt.legend()
        plt.xlim(index[0]-0.2,(index + bar_width)[-1]+bar_width+0.2)
        plt.tight_layout()
        
        ax2=fig.add_subplot(1,2,2)
        ax2.bar(index, means_sample_1- means_sample_2, bar_width+0.35,
                         alpha=opacity,
                         color=(0,1,1,0.4),
                         label='')
        
        
        plt.ylabel('% Sequences Difference')
        plt.title('%s - %s' % (args.sample_1_name,args.sample_2_name))
        plt.xticks(index + bar_width, ('Unmodified', 'NHEJ', 'HDR', 'Mixed'))

        
        plt.xlim(index[0]-bar_width, (index+bar_width)[-1]+2*bar_width)
        plt.tight_layout()
        plt.savefig(_jp('1.Comparison_Efficiency.pdf'), bbox_inches='tight')
        if args.save_also_png:
            plt.savefig(_jp('1.Comparison_Efficiency.png'), bbox_inches='tight')  
        
     
        #profile comparion  
        fig=plt.figure(figsize=(20,10))
        
        ax1=fig.add_subplot(1,2,1)
        plt.title('Mutation position distribution')
        y_max=max(effect_vector_any_1.max(),effect_vector_any_2.max())*1.2
         
        plt.plot(effect_vector_any_1,color=(0,0,1,0.3),lw=4,label='%s combined mutations' % args.sample_1_name)
        plt.hold(True)  
        plt.plot(effect_vector_any_2,color=(1,0,0,0.3),lw=4,label='%s combined mutations' % args.sample_2_name) 
        
        if cut_points:
            for idx,cut_point in enumerate(cut_points):
                if idx==0:    
                        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                else:
                        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2,label='_nolegend_')
             
                    
            for idx,sgRNA_int in enumerate(sgRNA_intervals):  
                if idx==0:    
                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA')
                else:
                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_')
                    
                   
        lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.3),ncol=1, fancybox=True, shadow=False)
         
     
        plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
        plt.xlabel('Reference amplicon position (bp)')
        plt.ylabel('Sequences %')
        plt.ylim(0,max(1,y_max))
        plt.xlim(xmax=len_amplicon-1)
        
        ax2=fig.add_subplot(1,2,2)
        
        effect_vector_any_diff=effect_vector_any_1-effect_vector_any_2
        
        y_max=effect_vector_any_diff.max()*1.2
        y_min=effect_vector_any_diff.min()*1.2
        
        plt.title('%s - %s' % (args.sample_1_name,args.sample_2_name))
        plt.plot(effect_vector_any_diff,color=(0,1,0,0.4),lw=3,label='Difference' )
    
            
        if cut_points:
            for idx,cut_point in enumerate(cut_points):
                if idx==0:    
                        plt.plot([cut_point,cut_point],[min(-1,y_min),max(1,y_max)],'--k',lw=2,label='Predicted cleavage position')
                else:
                        plt.plot([cut_point,cut_point],[min(-1,y_min),max(1,y_max)],'--k',lw=2,label='_nolegend_')
             
                    
            for idx,sgRNA_int in enumerate(sgRNA_intervals):  
                if idx==0:    
                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[min(-1,y_min),min(-1,y_min)],lw=10,c=(0,0,0,0.15),label='sgRNA')
                else:
                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[min(-1,y_min),min(-1,y_min)],lw=10,c=(0,0,0,0.15),label='_nolegend_')
        
        lgd2=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.2),ncol=1, fancybox=True, shadow=False)
        plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
        plt.xlabel('Reference amplicon position (bp)')
        plt.ylabel('Sequences Difference %')
        plt.xlim(xmax=len_amplicon-1)
        
        plt.ylim(min(-1,y_min),max(1,y_max))
        
        plt.savefig(_jp('2.Comparison_Combined_Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
        if args.save_also_png:
                plt.savefig(_jp('2.Comparison_Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')
      
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