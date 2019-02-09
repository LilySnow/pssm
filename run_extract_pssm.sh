#!/bin/bash
# 2017-01-20 14:38:40 CunliangGeng

dir_run="/home/clgeng/skempi/graph_kernel/features_test"
dir_out="${dir_run}"
list_mut="${dir_run}/list_mut_mall"
fout="pssm_ppro_itf.tsv"
dir_pssm="/home/clgeng/skempi/pssm/pssm_ori"
code_pssm="/home/clgeng/skempi/graph_kernel/features_test/extract_pssm.py"

# write graph list to array
arr=(`cat ${list_mut}`)
let num=${#arr[@]}-1
#echo ${arr[*]} # array index start from 0
#echo ${#arr[@]} # array length

for g1 in `seq 1 $num`; do
    # echo $g1,${arr[$g1]}
    mut=${arr[$g1]}
    read pdb chid reswt resi resmut <<< `echo $mut | awk 'BEGIN{}{split($1, arr, "_"); n=length(arr[3]); print arr[1], arr[2], substr(arr[3], 1, 3),  substr(arr[3], 4, n-6), substr(arr[3], n-2, 3)}'`
    python ${code_pssm} ${dir_pssm}/${pdb}_A_wild.pssm  ${dir_pssm}/${pdb}_B_wild.pssm ${mut} $chid $resi $reswt $resmut | sed -n '2p'
done > ${dir_out}/${fout}.tmp

python ${code_pssm} ${dir_pssm}/${pdb}_A_wild.pssm  ${dir_pssm}/${pdb}_B_wild.pssm ${mut} $chid $resi $reswt $resmut | sed -n '1p' > ${dir_out}/${fout}
cat ${dir_out}/${fout}.tmp >> ${dir_out}/${fout}
rm ${dir_out}/${fout}.tmp
