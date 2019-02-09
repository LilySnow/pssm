#!/bin/bash
# 2018-02-27 21:58:32 CunliangGeng

# Check the number of positions with IC=0 in the PSSM file.

[ $# -ne 1 ] && echo "Usage: $0 <PSSM path>" && echo "Example: $0 /home/test/1A22_A.pssm" && echo "Output: 1A22_A.pssm 15 (i.e. 15 positions have IC=0)" && exit
pssm=$1

nzero=`cat $pssm | tr -s " " | cut -d " "  -f 44 | grep -c "0.00"`

if [[ "nzero" -ge 1 ]]; then
    echo -e "$pssm\t$nzero"
fi
