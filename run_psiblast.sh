#!/bin/bash
# 2018-02-26 14:08:54 CunliangGeng
# Li Xue

if [[ $# -ne 3 ]]; then
    echo -e  "Usage:\n${0}  <input_fasta_file>  <output_pssmfile> <output_format_index> [num_threads]"
    echo -e  "  output_format_index: recommend 7(Tabular) or 13(json)."
    echo -e  "  path_fasta, path_output: default is current working directory.\n"
    echo -e  "Example1:\n ${0} /home/fasta/1A22_A.fasta  1A22_chainA.pssm  7 16 "
    echo -e  "Example2:\n ${0} /home/fasta/1A22_A.fasta  1A22_chainA.pssm  7 "
    exit
fi

#blast_db=/projects/0/deeprank/software/ncbi-blast-2.7.1+/databases/nr/nr
blast_db=/data/lixue/DBs/blast_dbs/nr_v20180204/nr
psiblast='/home/clgeng/software/blast/ncbi-blast-2.7.1+/bin/psiblast'


fseq=$1
pssmFL=$2
index=$3
num_threads=$4

if [ -z $num_threads ];then
    num_threads=1
fi

outDIR=`dirname $pssmFL`

if [ ! -d $outDIR ];then
    mkdir -p $outDIR
fi

# Set parameters based on sequence length
## reference: https://www.ncbi.nlm.nih.gov/books/NBK279684/

seqlen=`grep -v "^>" $fseq | wc | awk '{print $3-$1}'`

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
$psiblast\
 -query ${fseq}\
 -db $blast_db\
 -out $pssmFL.sa\
 -evalue 0.0001\
 -word_size ${wordSize}\
 -gapopen ${gapOpen}\
 -gapextend ${gapExtend}\
 -matrix  ${scoringMatrix}\
 -comp_based_stats 1\
 -outfmt '${index} qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids stitle salltitles sstrand qcovs qcovhsp qcovus'\
 -max_target_seqs 2000\
 -num_iterations 3\
 -out_ascii_pssm $pssmFL\
 -num_threads $num_threads
"


echo $command


