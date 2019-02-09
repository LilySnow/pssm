#!/usr/bin/env python

"""
Map the sequence of a model/decoy to official PDB full sequences to get which chain it is from.
Note that only one of the most matched chains will be output.

Usage: python map_seq_model2pdb.py <FASTA file of model/decoy>  <PDBID>  <path of PDB chain sequences>
Example: python map_seq_model2pdb.py  model.A.fasta 7CEI ./test/fasta

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
    if len(args) != 3:
        sys.exit(USAGE)


def map_sequence_model2pdb(modelfastafile, pdbid, pdbfastadir):
    """Map model sequence to PDB full sequence to see which chain it is from.

    Arguments:
        modelfastafile {str} -- model FASTA file
        pdbid {str} -- PDB identifier
        pdbfastadir {str} -- path of PDB FASTA files

   Returns:
       list -- model sequence ID, PDB sequence ID and the status ("Same", "Part" and "Diff").
    """

    # read model seq
    modelseq = read_fasta_single_sequence(modelfastafile)
    modelseqID = os.path.splitext(os.path.basename(modelfastafile))[0]

    # read PDB sequences
    query = os.path.join(pdbfastadir, pdbid+".*.fasta")
    pdbfastafiles = glob.glob(query)
    if not pdbfastafiles:
        try:
            raise ValueError("FASTA files of PDB '{}' NOT exist".format(pdbid))
        except ValueError as e:
            sys.exit(e)

    pdbseq = {}
    for i in pdbfastafiles:
        pdbseqID = os.path.splitext(os.path.basename(i))[0]
        pdbseq[pdbseqID] = read_fasta_single_sequence(i)

    # compare sequence between model and PDB
    map = {}
    mapsame = []
    mappart = []
    mapdiff = []
    for i in pdbseq:
        status, selfidentity = compare_2sequences(modelseq, pdbseq[i])
        map[i] = [modelseqID, i, status, selfidentity]
        # map.append([modelseqID, i, status, selfidentity])

        if status == "Same":
            mapsame.append(i)
        elif status == "Part":
            mappart.append(i)
        elif status == "Diff":
            mapdiff.append(i)
        else:
            pass

    return map, mapsame, mappart, mapdiff



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
        str -- percentage of matched sequence in 1st or query sequence
        status:
            Same: two sequences are same
            Part: SeqA is a part of seqB
            Diff: two sequences might be different, need manual check
    """

    if seqA == seqB:
        status = "Same"
        selfidentity = 1  # percentage 100%
    else:
        len_seqA = len(seqA)

        ali = pairwise2.align.globalxs(seqA, seqB, -2, -1)
        ali_seqA = np.array([i for i in ali[0][0]])
        ali_seqB = np.array([i for i in ali[0][1]])

        n_match = np.count_nonzero(ali_seqA == ali_seqB)
        selfidentity = n_match / len_seqA

        if selfidentity >= 0.6:
            status = "Part"
        else:
            status = "Diff"

    selfidentity = '{:.0%}'.format(selfidentity)
    return status, selfidentity

if __name__ == "__main__":

    check_input(sys.argv[1:])
    modelfastafile, pdbid, pdbfastadir = sys.argv[1:]
    # modelfastafile, pdbid, pdbfastadir = ("5CNP.A.fasta", "5CNP", "./test")

    map, mapsame, mappart, mapdiff = map_sequence_model2pdb(modelfastafile, pdbid, pdbfastadir)
    modelseqID = os.path.splitext(os.path.basename(modelfastafile))[0]

    header = ["#Query", "PDBchain", "Status", "QueryIdentity"]
    print("\t".join(header))
    if mapsame:
        print(modelseqID, mapsame[0], "Same", map[mapsame[0]][3], sep="\t")
    elif mappart:
        print(modelseqID, mappart[0], "Part", map[mappart[0]][3],  sep="\t")
    else:
        print(modelseqID, mapdiff[0], "toCheck", map[mapdiff[0]][3], sep="\t")

    # for i in map:
    #     print("\t".join(i))
