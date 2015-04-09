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

CRISPRESSO_VERSION=0.3
print 'Version %.2f\n' % CRISPRESSO_VERSION


import sys
import os
import subprocess as sb
import argparse
import logging
import re
import gzip


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


def filter_fastq_by_qual(fastq_filename,MIN_NT_QUALITY=10,output_filename=None):

    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

    if not output_filename:
        output_filename=fastq_filename.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'

    with gzip.open(output_filename,'w+') as fastq_filtered_outfile:

        for record in SeqIO.parse(fastq_handle, "fastq"):
            if (sum(np.array((record.letter_annotations["phred_quality"]))<MIN_NT_QUALITY)/float(len(record.letter_annotations["phred_quality"])))<1:
            #if min(record.letter_annotations["phred_quality"])>=MIN_NT_QUALITY:
                fastq_filtered_outfile.write(record.format('fastq'))

    return output_filename

trim_seq = lambda x: x[args.EXCLUDE_NT_FROM_SIDES:len(x)-args.EXCLUDE_NT_FROM_SIDES]


plt=check_library('pylab')
from matplotlib import font_manager as fm

pd=check_library('pandas')
np=check_library('numpy')
Bio=check_library('Bio')

check_program('java')
check_program('flash')
check_program('needle')

from Bio import SeqIO

parser = argparse.ArgumentParser(description='CRISPRESSO Parameters')
parser.add_argument('fastq_r1', type=str,  help='' )
parser.add_argument('fastq_r2', type=str,  help='')
parser.add_argument('amplicon_seq', type=str,  help='')

#optional
parser.add_argument('--guide_seq',  help='', default='')
parser.add_argument('--repair_seq',  help='', default='')
parser.add_argument('--MIN_NT_QUALITY', type=int, help='', default=30)
parser.add_argument('--min_identity_score', type=float, help='', default=50.0)
parser.add_argument('--name',  help='', default='')
parser.add_argument('--MAX_INSERTION_SIZE',  type=int, help='', default=30)
parser.add_argument('--PERFECT_ALIGNMENT_THRESHOLD',  type=float, help='', default=98.0)
parser.add_argument('--trim_sequences',help='',action='store_true')
parser.add_argument('--trimmomatic_options_string', type=str, default=' ILLUMINACLIP:NexteraPE-PE.fa:0:90:10:0:true MINLEN:40' )
parser.add_argument('--flash_options_string', type=str,default='')
parser.add_argument('--needle_options_string',type=str,default='-gapopen=10 -gapextend=0.5  -awidth3=5000')
parser.add_argument('--keep_intermediate',help='',action='store_true')
parser.add_argument('--output_folder',  help='', default='')
parser.add_argument('--dump',help='',action='store_true')
parser.add_argument('--EXLCUDE_MUTATIONS_COUNT',help='',action='store_true')
parser.add_argument('--EXCLUDE_NT_FROM_SIDES', type=int, help='', default=0)

args = parser.parse_args()

#make evetything uppercase!
if args.amplicon_seq:
    args.amplicon_seq=args.amplicon_seq.upper()
    
if args.repair_seq:
    args.repair_seq=args.repair_seq.upper()

if args.guide_seq:
    args.guide_seq=args.guide_seq.upper()
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

if args.MIN_NT_QUALITY>0:
    info('Filtering the reads with nt quality < %d ...' % args.MIN_NT_QUALITY)
    args.fastq_r1=filter_fastq_by_qual(args.fastq_r1,MIN_NT_QUALITY=args.MIN_NT_QUALITY,output_filename=_jp(os.path.basename(args.fastq_r1).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'))
    args.fastq_r2=filter_fastq_by_qual(args.fastq_r2,MIN_NT_QUALITY=args.MIN_NT_QUALITY,output_filename=_jp(os.path.basename(args.fastq_r2.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz'))

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
    cmd='java -jar /gcdata/gcproj/Luca/leo/programs/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
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


#plotting merged length reads
info('Calculating indel distribution based on the length of the merged sequences...')
df_hist=pd.read_table(flash_hist_filename,header=None,names=['overlap','density'])

plt.figure()
if args.guide_seq:
    min_cut=min(cut_points)
    max_cut=max(cut_points)
else:
    min_cut=len_amplicon/2
    max_cut=len_amplicon/2
    
len_amplicon=len(args.amplicon_seq)

plt.figure()
center_index=np.nonzero(df_hist.overlap==len_amplicon)[0]
plt.bar(0,df_hist['density'][center_index],color='red',linewidth=0)
plt.hold(True)
barlist=plt.bar(df_hist['overlap']-len_amplicon,df_hist['density'],align='center',linewidth=0)
barlist[center_index].set_color('r')
if args.guide_seq:
    plt.xlim([-min_cut,len_amplicon-max_cut])
else:
    plt.xlim([-min_cut,+max_cut])
    
plt.ylabel('# sequences')
plt.xlabel('Indel size (nt)')
plt.ylim([0,df_hist['density'].max()*1.2])
plt.title('Indel size distribution')
plt.legend(['unmodified','modified'])
plt.savefig(_jp('1a.Indel_size_distribution_n_sequences.pdf'))

plt.figure()
center_index=np.nonzero(df_hist.overlap==len_amplicon)[0]
plt.bar(0,df_hist['density'][center_index]/(df_hist['density'].sum())*100.0,color='red',linewidth=0)
plt.hold(True)
barlist=plt.bar(df_hist['overlap']-len_amplicon,df_hist['density']/(df_hist['density'].sum())*100.0,align='center',linewidth=0)
barlist[center_index].set_color('r')
plt.xlim([-min_cut,len_amplicon-max_cut])
plt.ylabel('% of sequences')
plt.xlabel('Indel size (nt)')
plt.title('Indel size distribution')
plt.legend(['unmodified','modified'])
plt.savefig(_jp('1b.Indel_size_distribution_percentage.pdf'))
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

def parse_needle_output(needle_filename,name='seq',just_score=False):
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
                
                identity_seq=eval(line.strip().split(' ')[-1].replace('%','').replace(')','').replace('(',''))
                    
                if just_score:
                    needle_data.append([id_seq,identity_seq])
                else:
                    for _ in range(7):
                        needle_infile.readline()

                    line=needle_infile.readline()
                    aln_ref_seq=line.split()[2]


                    aln_str=needle_infile.readline()[21:].rstrip()
                    line=needle_infile.readline()
                    aln_query_seq=line.split()[2]
                    aln_query_len=line.split()[3]
                    needle_data.append([id_seq,identity_seq,aln_query_len,aln_ref_seq,aln_str,aln_query_seq])

        if just_score:
            return pd.DataFrame(needle_data,columns=['ID','score_'+name]).set_index('ID')
        else:
            return pd.DataFrame(needle_data,columns=['ID','score_'+name,'length','ref_seq','align_str','align_seq']).set_index('ID')


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
    df_database=parse_needle_output(needle_output_filename,'ref')
    df_database_repair=parse_needle_output(needle_output_repair_filename,'repaired',just_score=True)
    df_database_and_repair=df_database.join(df_database_repair) 

    #filter bad alignments
    df_database_and_repair=df_database_and_repair.ix[(df_database_and_repair.score_ref>args.min_identity_score)|(df_database_and_repair.score_repaired>args.min_identity_score)]

    df_database_and_repair['score_diff']=df_database_and_repair.score_ref-df_database_and_repair.score_repaired

    N_REPAIRED=sum((df_database_and_repair.score_diff<0) & (df_database_and_repair.score_repaired>=args.PERFECT_ALIGNMENT_THRESHOLD))

    #df_database_and_repair.ix[:,['score_ref','score_repaired','score_diff']].to_csv(_jp('CRISPResso_SUMMARY_ALIGNMENT_IDENTITY_SCORE.txt'),header=['Identity_amplicon', 'Indentity_repaired_amplicon','Difference'],sep='\t')
    df_repaired=df_database_and_repair.ix[(df_database_and_repair.score_diff<0) & (df_database_and_repair.score_repaired>=args.PERFECT_ALIGNMENT_THRESHOLD)].sort('score_repaired',ascending=False)
    df_repaired.ix[:,['score_ref','score_repaired','score_diff']].to_csv(_jp('CRISPResso_REPAIRED_ONLY_IDENTITY_SCORE.txt'),header=['Identity_amplicon', 'Indentity_repaired_amplicon','Difference'],sep='\t')

#info('Parsing aligned files and making plots...')
#here we cover the case of the mutations plot instead..

#remove the HR events
if args.repair_seq:
    df_needle_alignment=df_database_and_repair.ix[df_database_and_repair.index.difference(df_repaired.index)]
    N_TOTAL=df_database_and_repair.shape[0]*1.0

else:
    df_needle_alignment=parse_needle_output(needle_output_filename,'ref')
    #filter out not aligned reads
    df_needle_alignment=df_needle_alignment.ix[df_needle_alignment.score_ref>args.min_identity_score]
    N_TOTAL=df_needle_alignment.shape[0]*1.0 #THIS SHOULD BE FIXED

if args.EXLCUDE_MUTATIONS_COUNT:
    check_seq_modified=lambda aln_str: 0 if (len(set(trim_seq(aln_str)).difference('.'))==1) else 1
else:
    check_seq_modified=lambda aln_str: 0 if (len(set(trim_seq(aln_str)))==1) else 1

df_needle_alignment['NHEJ']=df_needle_alignment.align_str.apply(check_seq_modified)

N_MODIFIED=df_needle_alignment['NHEJ'].sum()
N_UNMODIFIED=N_TOTAL-N_MODIFIED    

if args.repair_seq:
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(1,1,1)
    patches, texts, autotexts =ax.pie([N_UNMODIFIED,N_MODIFIED,N_REPAIRED],labels=['unmodified\n(%d)' %N_UNMODIFIED,'NHEJ\n(%d)' % N_MODIFIED, 'HR\n(%d)' %N_REPAIRED],explode=(0,0.05,0.1),colors=['w',(1,0,0,0.3),(0,0,1,0.3)],autopct='%1.1f%%')
    proptease = fm.FontProperties()
    proptease.set_size('xx-large')
    plt.setp(autotexts, fontproperties=proptease)
    plt.setp(texts, fontproperties=proptease)
    plt.savefig(_jp('2.Unmodified_NHEJ_HR_pie_chart.pdf'))

else:
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(1,1,1)
    patches, texts, autotexts =ax.pie([N_UNMODIFIED/N_TOTAL*100,N_MODIFIED/N_TOTAL*100],labels=['unmodified\n(%d)' %N_UNMODIFIED,'NHEJ\n(%d)' % N_MODIFIED],explode=(0,0.05),colors=['w',(1,0,0,0.2)],autopct='%1.1f%%')
    proptease = fm.FontProperties()
    proptease.set_size('xx-large')
    plt.setp(autotexts, fontproperties=proptease)
    plt.setp(texts, fontproperties=proptease)
    plt.savefig(_jp('2.Unmodified_NHEJ_pie_chart.pdf'))

df_needle_alignment['n_inserted']=df_needle_alignment['ref_seq'].apply(lambda x: trim_seq(x).count('-'))
df_needle_alignment['n_deleted']=df_needle_alignment['align_seq'].apply(lambda x: trim_seq(x).count('-'))
df_needle_alignment['n_mutated']=df_needle_alignment['align_str'].apply(lambda x: trim_seq(x).count('.'))


#(1) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);
y_values_mut,x_bins=plt.histogram(df_needle_alignment['n_mutated'],bins=range(0,60))
y_values_ins,x_bins=plt.histogram(df_needle_alignment['n_inserted'],bins=range(0,60))
y_values_del,x_bins=plt.histogram(df_needle_alignment['n_deleted'],bins=range(0,60))

fig=plt.figure(figsize=(20,10))

ax=fig.add_subplot(2,3,1)
ax.bar(x_bins[:-1],y_values_ins,align='center',linewidth=0)
barlist=ax.bar(x_bins[:-1],y_values_ins,align='center',linewidth=0)
barlist[0].set_color('r')
plt.title('Insertions')
plt.xlabel('size (nt)')
plt.ylabel('# sequences')
plt.legend(['Non-insertion','Insertion'][::-1])

ax=fig.add_subplot(2,3,2)
ax.bar(-x_bins[:-1],y_values_del,align='center',linewidth=0)
barlist=ax.bar(-x_bins[:-1],y_values_del,align='center',linewidth=0)
barlist[0].set_color('r')
plt.title('Deletions')
plt.xlabel('size (nt)')
plt.ylabel('# sequences')
plt.legend(['Non-deletion','Deletion'][::-1],loc=2)


ax=fig.add_subplot(2,3,3)
ax.bar(x_bins[:-1],y_values_mut,align='center',linewidth=0)
barlist=ax.bar(x_bins[:-1],y_values_mut,align='center',linewidth=0)
barlist[0].set_color('r')
plt.title('Mutations')
plt.xlabel('size (nt)')
plt.ylabel('# sequences')
plt.legend(['Non-mutation','Mutation'][::-1])

ax=fig.add_subplot(2,3,4)
ax.bar(x_bins[:-1],y_values_ins/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist=ax.bar(x_bins[:-1],y_values_ins/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist[0].set_color('r')
plt.xlabel('size (nt)')
plt.ylabel('% sequences')
plt.legend(['Non-insertion','Insertion'][::-1])

ax=fig.add_subplot(2,3,5)
ax.bar(-x_bins[:-1],y_values_del/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist=ax.bar(-x_bins[:-1],y_values_del/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist[0].set_color('r')
plt.xlabel('size (nt)')
plt.ylabel('% sequences')
plt.legend(['Non-deletion','Deletion'][::-1],loc=2)

ax=fig.add_subplot(2,3,6)
ax.bar(x_bins[:-1],y_values_mut/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist=ax.bar(x_bins[:-1],y_values_mut/float(df_needle_alignment.shape[0])*100.0,align='center',linewidth=0)
barlist[0].set_color('r')
plt.xlabel('size (nt)')
plt.ylabel('% sequences')
plt.legend(['Non-mutation','Mutation'][::-1])


plt.savefig(_jp('3.Insertion_Deletion_Mutation_size_hist.pdf'))



#(2) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);
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
    
    st_alg=trim_seq(row.align_seq).find('-')
    en_alg=trim_seq(row.align_seq).rfind('-')
    
    st_ref=trim_seq(row.ref_seq).find('-')
    en_ref=trim_seq(row.ref_seq).rfind('-')
    
    int_len_alg=(en_alg-st_alg) if (st_alg!=-1) else -1
    int_len_ref=(en_ref-st_ref) if (st_ref!=-1) else -1
    
    #print  int_len_alg, int_len_ref
    
    if ((int_len_alg==-1) and (int_len_ref==-1) )\
    or ((st_ref+args.EXCLUDE_NT_FROM_SIDES)==len(args.amplicon_seq)) or (st_ref==0) or (en_ref==(len(trim_seq(row.ref_seq))-1)): #here we remove the effect of the failed trimmed sequences 
        return None, None,'NOTHING'
    
    elif int_len_alg>=int_len_ref:
        return st_alg+args.EXCLUDE_NT_FROM_SIDES,en_alg+args.EXCLUDE_NT_FROM_SIDES,'DEL'
    
    else:
        return st_ref+args.EXCLUDE_NT_FROM_SIDES,en_ref+args.EXCLUDE_NT_FROM_SIDES,'INS'


#make plot
effect_vector_insertion=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
effect_vector_deletion=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
effect_vector_mutation=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))
effect_vector_any=np.zeros(len(df_needle_alignment.iloc[0].ref_seq))

problematic_seq=[]
exclude_idxs=range(args.EXCLUDE_NT_FROM_SIDES)+range(len(args.amplicon_seq)-args.EXCLUDE_NT_FROM_SIDES,len(args.amplicon_seq))
for idx_row,row in df_needle_alignment.iterrows():

    if not args.EXLCUDE_MUTATIONS_COUNT:
        mutated_positons=row.ref_positions[np.nonzero(np.array([1 if c=='.' else 0 for c in row.align_str]))[0]]
        effect_vector_mutation[mutated_positons]+=1
        effect_vector_mutation[-args.EXCLUDE_NT_FROM_SIDES:]=0
        effect_vector_mutation[:args.EXCLUDE_NT_FROM_SIDES]=0

        nt_modified=list(set(mutated_positons).difference(exclude_idxs))
    else:
        nt_modified=[]
        
    
    if (row['n_inserted']==0) and (row['n_deleted']==0):
        pass
        #print 'perfect match'
    else:

        idx_start,idx_end,edit_type=find_largest_interval(row)

        if edit_type=='DEL':
            #print 'deletion'
            
            idx_start_ref=row.ref_positions[idx_start]
            idx_end_ref=row.ref_positions[idx_end]
            
            if (not idx_start_ref > (len(args.amplicon_seq) - args.EXCLUDE_NT_FROM_SIDES)) or (not idx_end_ref > (len(args.amplicon_seq) - args.EXCLUDE_NT_FROM_SIDES)):
                effect_vector_deletion[idx_start_ref:idx_end_ref+1]+=1
                nt_modified+=range(idx_start_ref,idx_end_ref+1)
            

        elif edit_type=='INS':
            #print 'insertion'

                idx_start_ref=row.ref_positions[idx_start-1]
                if not idx_start_ref > (len(args.amplicon_seq) - args.EXCLUDE_NT_FROM_SIDES):
                    effect_vector_insertion[idx_start_ref]+=1
                    nt_modified.append(idx_start_ref)


                idx_end_ref=row.ref_positions[idx_end+1]
                if not idx_end_ref > (len(args.amplicon_seq) - args.EXCLUDE_NT_FROM_SIDES):
                    effect_vector_insertion[idx_end_ref]+=1
                    nt_modified.append(idx_end_ref)
                
        else:

            problematic_seq.append(row)


        if nt_modified!=[]:
            try:
                effect_vector_any[np.unique(nt_modified)]+=1
            except:
                print 'debug:',nt_modified
        

if len(problematic_seq)>0:
    warn('Skipped %d problematic sequences...' % len(problematic_seq) )
    pd.DataFrame(problematic_seq).to_csv(_jp('PROBLEMATIC_SEQUENCES.txt'),sep='\t')

plt.figure()
plt.plot(effect_vector_insertion,'r',lw=2)
plt.hold(True)
plt.plot(effect_vector_deletion,'m',lw=2)


if args.EXLCUDE_MUTATIONS_COUNT:
    labels_plot=['Insertions','Deletions']
else:
    plt.plot(effect_vector_mutation,'g',lw=2)
    labels_plot=['Insertions','Deletions','Mutations']

y_max=max(max(effect_vector_insertion),max(effect_vector_deletion),max(effect_vector_mutation))*1.2


if cut_points:
    for cut_point in cut_points:
        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2)
    lgd=plt.legend(labels_plot+['Predicted cleavage position'],loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)

else:
    lgd=plt.legend(labels_plot)
    

plt.xlabel('Amplicon position (nt)')
plt.ylabel('# sequences')
plt.ylim(ymax=y_max)
plt.xlim(xmax=len(args.amplicon_seq))
plt.title('Indel position distribution')
plt.savefig(_jp('4.Insertion_Deletion_Mutation_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')


plt.figure()

effect_vector_combined=100*effect_vector_any/float((N_TOTAL-len(problematic_seq)))
#effect_vector_combined=100*effect_vector_any/float((df_needle_alignment.shape[0]-len(problematic_seq)))

y_max=max(effect_vector_combined)*1.2

if cut_points:
    for cut_point in cut_points:
        plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2)
    lgd=plt.legend(['Predicted cleavage position'],loc='center', bbox_to_anchor=(0.5, -0.18),ncol=1, fancybox=True, shadow=True)
    
    
plt.hold(True)    
plt.plot(effect_vector_combined,'r',lw=2)
plt.title('Indel position distribution')
plt.xlabel('Amplicon position (nt)')
plt.ylabel('% of sequences')
plt.ylim(ymax=y_max)
plt.xlim(xmax=len(args.amplicon_seq))
plt.savefig(_jp('5.Combined_Insertion_Deletion_Mutation_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')


info('Done!')

if not args.keep_intermediate:
    info('Removing Intermediate files...')
    files_to_remove=[output_forward_paired_filename,output_reverse_paired_filename,\
                     flash_output_filename,flash_hist_filename,flash_histogram_filename,\
                     flash_not_combined_1_filename,flash_not_combined_2_filename,\
                     database_fasta_filename,query_fasta_filename] 

    if args.trim_sequences:
        files_to_remove+=[output_forward_unpaired_filename,output_reverse_unpaired_filename]

    if not args.dump:
        files_to_remove+=[needle_output_filename]
        if args.repair_seq:
            files_to_remove+=[needle_output_repair_filename]
        
    
    if args.repair_seq:
        files_to_remove+=[database_repair_fasta_filename,]

    if args.MIN_NT_QUALITY>0:
        files_to_remove+=[args.fastq_r1,args.fastq_r2]
        
   
    for file_to_remove in files_to_remove:
        os.remove(file_to_remove)

if args.dump:
    info('Dumping all the processed data...')
    np.savez(_jp('effect_vector_insertion'),effect_vector_insertion)
    np.savez(_jp('effect_vector_deletion'),effect_vector_deletion)
    np.savez(_jp('effect_vector_mutation'),effect_vector_mutation)
    np.savez(_jp('effect_vector_combined'),effect_vector_combined)
    #np.savez(_jp('effect_vector_combined'),(effect_vector_insertion+effect_vector_deletion+effect_vector_mutation)/float((df_needle_alignment.shape[0]-len(problematic_seq))))
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



