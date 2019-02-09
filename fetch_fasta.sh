#!/bin/bash
# 2018-03-08 15:59:36 CunliangGeng

if [[ $# -ne 1 ]]; then
    echo
    echo "Usage: $0 <PDB ID>"
    echo "Output: PDBID.fasta"
    echo
    echo "Author: Cunliang Geng"
    echo "Email: gengcunliang AT gmail.com"
    exit
fi

wget -q https://www.rcsb.org/pdb/download/downloadFastaFiles.do\?structureIdList\=$1\&compressionType\=uncompressed -O $1.fasta

