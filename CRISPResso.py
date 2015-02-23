#!/usr/bin/env python 
print '  \n~~~CRISPResso~~~'
print '-Analysis of CRISPR/CAS9 events from pair ended (PE) sequencing data-'
print '''
        ..
      ..  ..
            ..
             ..
            ..
           ..
         ..
##       ..    ####
##.............##  ##
##.............##   ##
##.............## ##
##.............###
 ##...........##
  #############
  #############
#################'''
print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',

CRISPRESSO_VERSION=0.1

import sys
import os
import subprocess as sb
import argparse
import logging
import re

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
        error('You need to install %s to use Crispresso!' % library_name)
        sys.exit(1)

def which(program):
    import os
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

def check_program(binary_name,download_url=None):
    if not which(binary_name):
        error('You need to install and have the command #####%s##### in your PATH variable to use Crispresso!' % binary_name)
        if download_url:
            error('You can download it from here:%s' % download_url)
        sys.exit(1)


nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N'})

def reverse_complement(seq):
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

plt=check_library('pylab')
pd=check_library('pandas')
np=check_library('numpy')

check_program('java')
check_program('flash')
check_program('needle')

parser = argparse.ArgumentParser(description='CRISPRESSO Parameters')
parser.add_argument('fastq_r1', type=str,  help='' )
parser.add_argument('fastq_r2', type=str,  help='')
parser.add_argument('amplicon_seq', type=str,  help='')

#optional
parser.add_argument('--guide_seq',  help='', default='')
parser.add_argument('--repair_seq',  help='', default='')
parser.add_argument('--name',  help='', default='')
parser.add_argument('--MAX_INSERTION_SIZE',  type=int, help='', default=30)
parser.add_argument('--PERFECT_ALIGNMENT_THRESHOLD',  type=float, help='', default=98.0)
parser.add_argument('--trim_sequences',help='',action='store_true')
parser.add_argument('--trimmomatic_options_string', type=str, default=' ILLUMINACLIP:NexteraPE-PE.fa:1:30:10:5:true LEADING:1 TRAILING:1 SLIDINGWINDOW:2:28 MINLEN:60' )
parser.add_argument('--flash_options_string', type=str,default='')
parser.add_argument('--needle_options_string',type=str,default='-gapopen=10 -gapextend=0.5  -awidth3=300')
parser.add_argument('--keep_intermediate',help='',action='store_true')
parser.add_argument('--output_folder',  help='', default='')
parser.add_argument('--dump',help='',action='store_true')



args = parser.parse_args()

if args.guide_seq:
    cut_points=[m.start() +len(args.guide_seq)-3.5 for m in re.finditer(args.guide_seq, args.amplicon_seq)]+[m.start() +2.5 for m in re.finditer(reverse_complement(args.guide_seq), args.amplicon_seq)]

    if not cut_points:
        error('The guide sequences provided is not present in the amplicon sequence! \n\nPlease check your input!')
        sys.exit(1)
    info('Cut Points from guide seq:%s' % cut_points)
else:
    cut_points=[]



get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

if not args.name:
    database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
else:
    database_id=args.name

OUTPUT_DIRECTORY='Crispresso_on_%s' % database_id

if args.output_folder:
    OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
    
_jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

try:
    info('Creating Folder %s' % OUTPUT_DIRECTORY)
    os.makedirs(OUTPUT_DIRECTORY)
    info('Done!')
except:
    warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

log_filename=_jp('CRISPResso_RUNNING_LOG.txt')

with open(log_filename,'w+') as outfile:
    outfile.write('[Command used]:\nCRISPResso %s\n\n\n[Other tools log]:\n' % ' '.join(sys.argv))

if not args.trim_sequences:
    output_forward_paired_filename=args.fastq_r1
    output_reverse_paired_filename=args.fastq_r2
else:
    info('Trimming sequences with Trimmomatic...')
    output_forward_paired_filename=_jp('output_forward_paired.fq')
    output_forward_unpaired_filename=_jp('output_forward_unpaired.fq') 
    output_reverse_paired_filename=_jp('output_reverse_paired.fq') 
    output_reverse_unpaired_filename=_jp('output_reverse_unpaired.fq')

    #Trimming with trimmomatic
    cmd='java -jar /gcdata/gcproj/Luca/leo/programs/Trimmomatic-0.32/trimmomatic-0.32.jar PE   %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
    % (args.fastq_r1,args.fastq_r2,output_forward_paired_filename,output_forward_unpaired_filename,output_reverse_paired_filename,output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
    sb.call(cmd,shell=True)
    info('Done!')

#Merging with Flash
info('Merging paired sequences with Flash...')
if args.repair_seq:
    len_amplicon=len(args.repair_seq)+args.MAX_INSERTION_SIZE #considering some tolerance for new insertion
else:
    len_amplicon=len(args.amplicon_seq)+args.MAX_INSERTION_SIZE #considering some tolerance for new insertion

    
cmd='flash %s %s --max-overlap=%s -d %s %s >>%s 2>&1' %\
     (output_forward_paired_filename,output_reverse_paired_filename,len_amplicon,OUTPUT_DIRECTORY,args.flash_options_string,log_filename)
sb.call(cmd,shell=True)
info('Done!')

flash_output_filename=_jp('out.extendedFrags.fastq')
flash_hist_filename=_jp('out.hist')
flash_histogram_filename=_jp('out.histogram')
flash_not_combined_1_filename=_jp('out.notCombined_1.fastq')
flash_not_combined_2_filename=_jp('out.notCombined_2.fastq')


#plotting merged lenght reads
info('Calculating length distribution...')
df_hist=pd.read_table(flash_hist_filename,header=None,names=['overlap','density'])
plt.bar(df_hist['overlap'],df_hist['density'],align='center')
plt.xlim([len(args.amplicon_seq)-args.MAX_INSERTION_SIZE,len(args.amplicon_seq)+args.MAX_INSERTION_SIZE])
plt.ylabel('#Sequences')
plt.xlabel('Combined reads 1,2 lenght')
plt.hold(True)
plt.plot([len(args.amplicon_seq),len(args.amplicon_seq)],[0,df_hist['density'].max()*1.2],'--r',lw=3)
plt.ylim([0,df_hist['density'].max()*1.2])
plt.legend(['ref.amplicon (%dbp)' %len(args.amplicon_seq)])
plt.savefig(_jp('1.Merged_Reads_Lenghts_Distribution.pdf'))
info('Done!')


info('Preparing files for the alignment...')
#parsing flash output and prepare the files for alignment
data_to_parse=[]
with open(flash_output_filename) as r1_file:
    for idx,line in enumerate(r1_file):
        if (idx % 4) ==0:
            seq_id=line.split()[0]
        if (idx % 4) ==1:
            seq=line.strip()
            
        if (idx %4) == 3:
            qual=line.strip()
            data_to_parse.append((seq_id,seq,qual))
df_R1R2=pd.DataFrame(data_to_parse,columns=['ID','SEQ_R1R2','QUAL_R1R2']).set_index('ID')


database_fasta_filename=_jp('%s_database.fa' % database_id)
query_fasta_filename=_jp('%s_query.fa' % database_id)
needle_output_filename=_jp('needle_output_%s.txt' % database_id)
crispresso_output_filename=_jp('CRISPRresso_aligned_%s.txt' % database_id)

#write .fa files
with open(database_fasta_filename,'w+') as outfile:
    outfile.write('>%s\n%s\n' % (database_id,args.amplicon_seq))

with open(query_fasta_filename,'w+') as outfile:
    for seq_id,row in df_R1R2.iterrows():
        outfile.write('>%s\n%s\n' % (seq_id.replace(':','_')+'_R1R2',row['SEQ_R1R2']))

if args.repair_seq:
    database_repair_fasta_filename=_jp('%s_database_repair.fa' % database_id)
    needle_output_repair_filename=_jp('needle_output_repair_%s.txt' % database_id)
    
    with open(database_repair_fasta_filename,'w+') as outfile:
        outfile.write('>%s\n%s\n' % (database_id,args.repair_seq))
info('Done!')

info('Aligning sequences...')
#Alignment here
cmd='needle -asequence=%s -bsequence=%s -outfile=%s %s >>%s 2>&1' \
     %(database_fasta_filename,query_fasta_filename,needle_output_filename,args.needle_options_string,log_filename)
sb.call(cmd,shell=True)

#If we have a donor sequence we just compare the fq in the two cases and exit
if args.repair_seq:

    cmd='needle -asequence=%s -bsequence=%s -outfile=%s %s >>%s 2>&1'\
         %(database_repair_fasta_filename,query_fasta_filename,needle_output_repair_filename,args.needle_options_string,log_filename)
    sb.call(cmd,shell=True)
    info('Done!')

    info('Parsing aligned files and making plots...')
    def parse_needle_output_just_score(needle_filename,name=''):
        needle_data=[]
        with open(needle_filename) as needle_infile:

            line=needle_infile.readline()
            while line:

                while line and ('# Aligned_sequences' not  in line):
                    line=needle_infile.readline()

                if line:
                    #print line
                    needle_infile.readline() #skip another line

                    line=needle_infile.readline()
                    id_seq=line.split()[-1].replace('_',':')

                    for _ in range(5):
                        needle_infile.readline()

                    line=needle_infile.readline()
                    identity_seq=eval(line.strip().split(' ')[-1].replace('%',''))

                    needle_data.append([id_seq,identity_seq])

            return pd.DataFrame(needle_data,columns=['ID','score_'+name]).set_index('ID')


    df_database=parse_needle_output_just_score(needle_output_filename,'ref')
    df_database_repair=parse_needle_output_just_score(needle_output_repair_filename,'repaired')
    df_database_and_repair=df_database.join(df_database_repair) 


    df_database_and_repair['score_diff']=df_database_and_repair.score_ref-df_database_and_repair.score_repaired

    N_TOTAL=float(df_database_and_repair.shape[0])
    N_ALIGNED=sum((df_database_and_repair.score_ref>=args.PERFECT_ALIGNMENT_THRESHOLD))
    N_REPAIRED=sum((df_database_and_repair.score_diff<0) & (df_database_and_repair.score_repaired>=args.PERFECT_ALIGNMENT_THRESHOLD))
    N_FAILED=N_TOTAL-N_ALIGNED-N_REPAIRED

    fig=plt.figure(figsize=(8,3))
    ax1=fig.add_subplot(1,1,1)
    ax1.barh([0,1,2],[N_ALIGNED/N_TOTAL*100,N_REPAIRED/N_TOTAL*100,N_FAILED/N_TOTAL*100])
    ax1.set_yticks([0.5,1.5,2.5])
    ax1.set_yticklabels(['Not modified (%d)' % N_ALIGNED,'HR (%d)' % N_REPAIRED,'Others (%d)' %N_FAILED])
    plt.xlabel('% Sequences')
    fig.tight_layout()
    plt.savefig(_jp('2.Not_Modified_HR_others_bar_plot.pdf'))


    df_database_and_repair.to_csv(_jp('CRISPResso_SUMMARY_ALIGNMENT_IDENTITY_SCORE.txt'),header=['Identity_amplicon', 'Indentity_repaired_amplicon','Difference'],sep='\t')

    df_repaired=df_database_and_repair.ix[(df_database_and_repair.score_diff<0) & (df_database_and_repair.score_repaired>=args.PERFECT_ALIGNMENT_THRESHOLD)].sort('score_repaired',ascending=False)
    df_repaired.to_csv(_jp('CRISPResso_REPAIRED_ONLY_IDENTITY_SCORE.txt'),header=['Identity_amplicon', 'Indentity_repaired_amplicon','Difference'],sep='\t')

    info('Done!')
    
    if not args.keep_intermediate:
        info('Removing Intermediate files...')
        files_to_remove=[output_forward_paired_filename,output_reverse_paired_filename,\
                         flash_output_filename,flash_hist_filename,flash_histogram_filename,\
                         flash_not_combined_1_filename,flash_not_combined_2_filename,\
                         database_fasta_filename,query_fasta_filename,database_repair_fasta_filename] 

        if args.trim_sequences:
            files_to_remove+=[output_forward_unpaired_filename,output_reverse_unpaired_filename]

        
        for file_to_remove in files_to_remove:
            os.remove(file_to_remove)
   
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

info('Done!')

info('Parsing aligned files and making plots...')
#here we cover the case of the mutations plot instead..
df_needle_alignment=pd.DataFrame()
needle_data=[]
with open(needle_output_filename) as needle_infile,open(crispresso_output_filename,'w+') as outfile:
    
    for _ in range(18):
        needle_infile.readline()
        
    line=needle_infile.readline()
    id_seq=line.split()[-1].replace('_',':')
    
    for _ in range(13):
        needle_infile.readline()
    
    line=needle_infile.readline()
    aln_ref_seq=line.split()[2]
        
    aln_str=needle_infile.readline()[21:].rstrip()
    line=needle_infile.readline()
    aln_query_seq=line.split()[2]
    aln_query_len=line.split()[3]
    
    outfile.write('%s %s\n' %(id_seq,aln_query_len))
    outfile.write('\n'.join([aln_ref_seq,aln_str,aln_query_seq,'\n']))
    needle_data.append([id_seq,int(aln_query_len),aln_ref_seq,aln_str,aln_query_seq])
    
    end_of_file=False
    #idx=0
    N_TOTAL=0.0
    N_MODIFIED=0.0
    N_NON_MODIFIED=0.0
    while line and not end_of_file:
    
        N_TOTAL+=1
        try:
            for _ in range(6):
                needle_infile.readline()
       
            line=needle_infile.readline()
            id_seq=line.split()[-1].replace('_',':')
 
 
            for _ in range(13):
                needle_infile.readline()
 
 
            line=needle_infile.readline()
            aln_ref_seq=line.split()[2]
 
 
            aln_str=needle_infile.readline()[21:].rstrip()
            line=needle_infile.readline()
            aln_query_seq=line.split()[2]
            aln_query_len=line.split()[3]
            
            outfile.write('%s %s\n' %(id_seq,aln_query_len))
            outfile.write('\n'.join([aln_ref_seq,aln_str,aln_query_seq,'\n']))  
            needle_data.append([id_seq,int(aln_query_len),aln_ref_seq,aln_str,aln_query_seq])


            if len(set(aln_str))==1:
                N_NON_MODIFIED+=1
            else:
                N_MODIFIED+=1
 
        except:
            #print 'PROBLEMA:',line
            end_of_file=True

        


df_needle_alignment=pd.DataFrame(needle_data,columns=['ID','length','ref_seq','align_str','align_seq']).set_index('ID') 

fig=plt.figure(figsize=(8,3))
ax1=fig.add_subplot(1,1,1)
ax1.barh([0,1],[N_NON_MODIFIED/N_TOTAL*100,N_MODIFIED/N_TOTAL*100])
ax1.set_yticks([0.5,1.5])
ax1.set_yticklabels(['Not modified/perfect repair(%d)' % N_NON_MODIFIED,'NHEJ (%d)' % N_MODIFIED])
plt.xlabel('% Sequences')
fig.tight_layout()
plt.savefig(_jp('2.Not_Modified_NHEJ_bar_plot.pdf'))

      
#y_values_len,x_bins=plt.histogram(df_needle_alignment.ix[:,'length'],bins=range(80,180))
#plt.bar(x_bins[:-1],y_values_len,align='center')
#plt.title('Lenghts of the aligned fragments after filtering')
#plt.xlabel('size')
#plt.ylabel('# fragments')
#plt.xlim(80,180)
#plt.savefig(_jp('1b.Length_of_the_aligned_fragments_after_filtering_and_merging.pdf'))


df_needle_alignment['n_inserted']=df_needle_alignment['ref_seq'].apply(lambda x: x.count('-'))
df_needle_alignment['n_deleted']=df_needle_alignment['align_seq'].apply(lambda x: x.count('-'))
df_needle_alignment['n_mutated']=df_needle_alignment['align_str'].apply(lambda x: x.count('.'))


#(1) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);
y_values_mut,x_bins=plt.histogram(df_needle_alignment['n_mutated'],bins=range(0,60))
y_values_ins,x_bins=plt.histogram(df_needle_alignment['n_inserted'],bins=range(0,60))
y_values_del,x_bins=plt.histogram(df_needle_alignment['n_deleted'],bins=range(0,60))

fig=plt.figure(figsize=(20,10))

ax=fig.add_subplot(2,3,1)
ax.bar(x_bins[:-1],y_values_ins,align='center')
plt.title('Insertions')
plt.xlabel('size')
plt.ylabel('# sequences')

ax=fig.add_subplot(2,3,2)
ax.bar(x_bins[:-1],y_values_del,align='center')
plt.title('Deletions')
plt.xlabel('size')
plt.ylabel('# sequences')


ax=fig.add_subplot(2,3,3)
ax.bar(x_bins[:-1],y_values_mut,align='center')
plt.title('Mutations')
plt.xlabel('size')
plt.ylabel('# sequences')

ax=fig.add_subplot(2,3,4)
ax.bar(x_bins[:-1],y_values_ins/float(df_needle_alignment.shape[0])*100.0,align='center')
plt.xlabel('size')
plt.ylabel('% sequences')

ax=fig.add_subplot(2,3,5)
ax.bar(x_bins[:-1],y_values_del/float(df_needle_alignment.shape[0])*100.0,align='center')
plt.xlabel('size')
plt.ylabel('% sequences')

ax=fig.add_subplot(2,3,6)
ax.bar(x_bins[:-1],y_values_mut/float(df_needle_alignment.shape[0])*100.0,align='center')
plt.xlabel('size')
plt.ylabel('% sequences')
plt.savefig(_jp('3.Insertion_Deletion_Mutation_length_hist.pdf'))

#(2) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);
pos_idxs=[]
idx=0
for c in df_needle_alignment.iloc[3].ref_seq:
    if c in set(['A','T','C','G']):
        pos_idxs.append(idx)
        idx+=1
    else:
        pos_idxs.append(-idx)
pos_idxs=np.array(pos_idxs)

def compute_ref_positions(ref_seq):
    pos_idxs=[]
    idx=0
    for c in ref_seq:
        if c in set(['A','T','C','G']):
            pos_idxs.append(idx)
            idx+=1
        else:
            if idx==0:
                pos_idxs.append(-1)
            else:   
                pos_idxs.append(-idx)
    return np.array(pos_idxs)

#compute positions relative to alignmnet
df_needle_alignment['ref_positions']=df_needle_alignment['ref_seq'].apply(compute_ref_positions)

def find_largest_interval(row):
    st_alg=row.align_seq.find('-')
    en_alg=np.abs(row.align_seq.rfind('-'))
    
    st_ref=row.ref_seq.find('-')
    en_ref=np.abs(row.ref_seq.rfind('-'))
    
    int_len_alg=(en_alg-st_alg) if (st_alg!=-1) else -1
    int_len_ref=(en_ref-st_ref) if (st_ref!=-1) else -1
    
    #print  int_len_alg, int_len_ref
    

    if int_len_alg>=int_len_ref:
        return st_alg,en_alg,'DEL'
    else:
        return st_ref,en_ref,'INS'        

#make plot
effect_vector_insertion=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
effect_vector_deletion=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
effect_vector_mutation=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
problematic_seq=[]
for idx_row,row in df_needle_alignment.iterrows():
    mutated_positons=row.ref_positions[np.nonzero(np.array([1 if c=='.' else 0 for c in row.align_str]))[0]]
    effect_vector_mutation[mutated_positons]+=1

    if (row['n_inserted']==0) and (row['n_deleted']==0):
        pass
        #print 'perfect match'
    else:

        idx_start,idx_end,edit_type=find_largest_interval(row)
        #print row['align_seq'][idx_start],row['align_seq'][idx_end]

        if edit_type=='DEL':
            #print 'deletion'

            idx_start_ref=row.ref_positions[idx_start]
            idx_end_ref=row.ref_positions[idx_end]
            
            #print idx_start,idx_end, idx_start_ref,idx_end_ref, database_seq[idx_start_ref:idx_end_ref+1]
            
            effect_vector_deletion[idx_start_ref:idx_end_ref+1]+=1

        elif edit_type=='INS':
            #print 'insertion'

            try:
                idx_start_ref=row.ref_positions[idx_start-1]
                effect_vector_insertion[idx_start_ref]+=1
            except:
                #print row
                problematic_seq.append(row)
            
            try:
                idx_end_ref=row.ref_positions[idx_end+1]
                effect_vector_insertion[idx_end_ref]+=1
            except:
                #print row
                problematic_seq.append(row)
                
            #print idx_start,idx_end, idx_start_ref,idx_end_ref, database_seq[idx_start_ref],database_seq[idx_end_ref]

        else:
            raise Exception('Something wrong')

plt.figure()
plt.plot(effect_vector_insertion,'r',lw=2)
plt.hold(True)
plt.plot(effect_vector_deletion,'m',lw=2)
plt.plot(effect_vector_mutation,'g',lw=2)


y_max=max(max(effect_vector_insertion),max(effect_vector_deletion),max(effect_vector_mutation))*1.2


if cut_points:
    for cut_point in cut_points:
        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2)
    plt.legend(['Insertions','Deletions','Mutations','Cut point'])

else:
    plt.legend(['Insertions','Deletions','Mutations'])
    

plt.xlabel('bp')
plt.ylabel('# sequences')
plt.ylim(ymax=y_max)
plt.savefig(_jp('4.Insertion_Deletion_Mutation_Locations.pdf'))


plt.figure()

effect_vector_combined=100*(effect_vector_insertion+effect_vector_deletion+effect_vector_mutation)/float((df_needle_alignment.shape[0]-len(problematic_seq)))
y_max=max(effect_vector_combined)*1.2

if cut_points:
    for cut_point in cut_points:
        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2)
    plt.legend(['Cut point'])
    
plt.hold(True)    
plt.plot(effect_vector_combined,'r',lw=2)
plt.title('Combined Deletion, Insertion, Mutation')
plt.xlabel('bp')
plt.ylabel('% of sequences')
plt.ylim(ymax=y_max)
plt.savefig(_jp('5.Combined_Insertion_Deletion_Mutation_Locations.pdf'))

info('Done!')




if not args.keep_intermediate:
    info('Removing Intermediate files...')
    files_to_remove=[output_forward_paired_filename,output_reverse_paired_filename,\
                     flash_output_filename,flash_hist_filename,flash_histogram_filename,\
                     flash_not_combined_1_filename,flash_not_combined_2_filename,\
                     database_fasta_filename,query_fasta_filename,needle_output_filename] 

    if args.trim_sequences:
        files_to_remove+=[output_forward_unpaired_filename,output_reverse_unpaired_filename]

   
    for file_to_remove in files_to_remove:
        os.remove(file_to_remove)

if args.dump:
    info('Dumping all the processed data...')
    np.savez(_jp('effect_vector_insertion'),effect_vector_insertion)
    np.savez(_jp('effect_vector_deletion'),effect_vector_insertion)
    np.savez(_jp('effect_vector_mutation'),effect_vector_insertion)
    np.savez(_jp('effect_vector_combined'),(effect_vector_insertion+effect_vector_deletion+effect_vector_mutation)/float((df_needle_alignment.shape[0]-len(problematic_seq))))
    df_needle_alignment.to_pickle(_jp('df_needle_alignment'))        
    
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



