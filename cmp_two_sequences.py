#!/usr/bin/env python

"""
Compare two sequences to see if they are same or not.

Usage: python cmp_two_sequences.py <FASTA file of 1st/query sequence>  <FASTA file of 2nd sequence>
Example: python cmp_two_sequences.py 7cei_A.fasta 7cei_B.fasta

Author: {0} ({1})
"""

import os
import sys
import glob
import numpy as np
from Bio import pairwise2

__author__ = "Cunliang Geng"
__email__ = "gengcunliang AT gmail.com"
USAGE = __doc__.format(__author__, __email__)

def check_input(args):
    """Validate user input
    
    Arguments:
        args {tuple} --  user input arguments
    """
    if len(args) != 2:
        sys.exit(USAGE)


def read_fasta_single_sequence(fastafile):
    """Read FASTA file with single sequence.

    Arguments:
        fastafile {str} -- FASTA file with single sequence.

    Returns:
        str -- one line sequence without description.
    """

    with open(fastafile) as f:
        lines =  f.readlines()
        if lines[0].startswith(">"):
            seq = "".join([i.strip() for i in lines[1:]])
        else:
            seq = "".join([i.strip() for i in lines[:]])

    return seq


def compare_2sequences(seqA, seqB):
    """Compare two protein sequences to see if they are same or not.

    Arguments:
        seqA {str} -- 1st protein sequence, or query sequence
        seqB {str} -- 2nd protein sequence, or template sequence

    Returns:
        str -- status
        str -- percentage of matched sequence in 1st sequence
        str -- percentage of matched sequence in 2nd sequence
        status:
            Same: two sequences are same
            Part: SeqA is a part of seqB
            Diff: two sequences might be different, need manual check
    """

    if seqA == seqB:
        status = "Same"
        identity1 = 1
        identity2 = 1
    else:
        len_seqA = len(seqA)
        len_seqB = len(seqB)

        ali = pairwise2.align.globalxs(seqA, seqB, -2, -1)
        ali_seqA = np.array([i for i in ali[0][0]])
        ali_seqB = np.array([i for i in ali[0][1]])
        # print(ali[0][0])
        # print(ali[0][1])
        n_match = np.count_nonzero(ali_seqA == ali_seqB)
        identity1 = n_match / len_seqA 
        identity2 = n_match / len_seqB

        # complexes are highly probably hetero when both identity values lower than 0.8
        if identity1 >= 0.8 or identity2 >= 0.8:
            status = "Part"
        else:
            status = "Diff"

    identity1 = '{:.0%}'.format(identity1)
    identity2 = '{:.0%}'.format(identity2)
    return status, identity1, identity2


def main(fseq1, fseq2):

    idseq1 = os.path.splitext(os.path.basename(fseq1))[0]
    idseq2 = os.path.splitext(os.path.basename(fseq2))[0]
    seq1 = read_fasta_single_sequence(fseq1)
    seq2 = read_fasta_single_sequence(fseq2)

    status, identity1, identity2 = compare_2sequences(seq1, seq2)

    header = ["#Sequence1", "Sequence2", "Status", "Identity1", "Identity2"]
    print("\t".join(header))
    print(idseq1, idseq2, status, identity1, identity2, sep="\t")

if __name__ == "__main__":

    check_input(sys.argv[1:])
    fseq1, fseq2 = sys.argv[1:]
    # fseq1 = "./test/4XC4.C.fasta ./test/4XC4.A.fasta"
    # fseq2 = "./test/4XC4.C.fasta ./test/5CNP.C.fasta"
    main(fseq1, fseq2)