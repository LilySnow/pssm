#!/bin/bash
# 2018-02-26 14:08:54 CunliangGeng

if [[ $# -ne 4 && $# -ne 6 ]]; then
    echo -e  "Usage:\n${0}  <input_fasta_file>  <output_file_name> <output_format_index> [path_fasta] [path_output] [num_threads]"
    echo -e  "  output_format_index: recommend 7(Tabular) or 13(json)."
    echo -e  "  path_fasta, path_output: default is current working directory.\n"
    echo -e  "Example1:\n ${0} 1A22_A.fasta  1A22_chainA  7 /home/pssm/fasta  /home/pssm/out"
    echo -e  "Example2:\n ${0} 1A22_A.fasta  1A22_chainA  7\n"
    exit
fi


fseq=$1
fout=$2
index=$3
if [[ $# -eq 6 ]]; then
    dir_fin=$4
    dir_fout=$5
    num_threads=$6
else
    dir_fin=`pwd`
    dir_fout=`pwd`
    num_threads=1
fi


# Set parameters based on sequence length
## reference: https://www.ncbi.nlm.nih.gov/books/NBK279684/

seqlen=`grep -v "^>" ${dir_fin}/${fseq} | wc | awk '{print $3-$1}'`

if [[ $seqlen -lt 30 ]]; then
    wordSize=2
    gapOpen=9
    gapExtend=1
    scoringMatrix=PAM30
elif [[ $seqlen -ge 30 && $seqlen -lt 35 ]]; then
    wordSize=3
    gapOpen=9
    gapExtend=1
    scoringMatrix=PAM30
elif [[ $seqlen -ge 35 && $seqlen -lt 50 ]]; then
    wordSize=3
    gapOpen=10
    gapExtend=1
    scoringMatrix=PAM70
elif [[ $seqlen -ge 50 && $seqlen -le 85 ]]; then
    wordSize=3
    gapOpen=10
    gapExtend=1
    scoringMatrix=BLOSUM80
else
    wordSize=3
    gapOpen=11
    gapExtend=1
    scoringMatrix=BLOSUM62
fi


# Output psiblast command

# command="
# psiblast\
#  -query ${dir_fin}/${fseq}\
#  -db /projects/0/deeprank/software/ncbi-blast-2.7.1+/databases/nr/nr\
#  -out ${dir_fout}/${fout}.sa\
#  -evalue 0.0001\
#  -word_size ${wordSize}\
#  -gapopen ${gapOpen}\
#  -gapextend ${gapExtend}\
#  -matrix  ${scoringMatrix}\
#  -comp_based_stats 1\
#  -outfmt '${index} qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids stitle salltitles sstrand qcovs qcovhsp qcovus'\
#  -max_target_seqs 2000\
#  -export_search_strategy ${dir_fout}/${fout}.searchstrategy\
#  -num_iterations 3\
#  -out_pssm ${dir_fout}/${fout}.cptpssm\
#  -out_ascii_pssm ${dir_fout}/${fout}.pssm\
#  -save_pssm_after_last_round\
#  -save_each_pssm\
#  -num_threads 16
# "

command="
psiblast\
 -query ${dir_fin}/${fseq}\
 -db /projects/0/deeprank/software/ncbi-blast-2.7.1+/databases/nr/nr\
 -out ${dir_fout}/${fout}.sa\
 -evalue 0.0001\
 -word_size ${wordSize}\
 -gapopen ${gapOpen}\
 -gapextend ${gapExtend}\
 -matrix  ${scoringMatrix}\
 -comp_based_stats 1\
 -outfmt '${index} qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids stitle salltitles sstrand qcovs qcovhsp qcovus'\
 -max_target_seqs 2000\
 -num_iterations 3\
 -out_ascii_pssm ${dir_fout}/${fout}.pssm\
 -num_threads $num_threads
"


echo $command

##-- submit to cluster
#sbatch <<END
##!/bin/bash
##SBATCH -J pssm_${fseq}
##SBATCH -n 16
##SBATCH -t 08:00:00
##SBATCH -o slurm-%x-%j.out
#
#cd $dir_fin
#eval echo \$command
#$command
#
#END
#
