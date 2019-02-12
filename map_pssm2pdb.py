#!/usr/bin/env python

"""
Get PSSM file and PDB file with consistent sequence between these two files (i.e. no gap, no residue X).

Usage: python map_pssm2pdb.py <input PSSM file> <input PDB file> <ChainID of PDB file> <Output path>
Example: python map_pssm2pdb.py  ./pssm/4CPA.A.pssm  ./pdb/4CPA.pdb  A  ./test
Output: 4CPA.A.pdb.pssm, 4CPA.A.pssm.pdb

Author: {0} ({1})
"""

import os
import sys
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
    if len(args) != 4:
        sys.exit(USAGE)


def get_pdb_seq(fpdb, chainID):
    """Get sequence residue names and IDs from PDB file.

    Arguments:
        fpdb {str} -- input pdb file
        chainID {str} -- target chain

    Raises:
        ValueError -- ChainID not exist in the pdb file

    Returns:
        numpy vector - seq_resn, sequence residue names
        numpy vector - seq_resi, sequence residue IDs
    """
    res_code_3to1 = dict([
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        ('ASX', 'B'), ('SEC', 'U'), ('GLX', 'Z'),
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ])

    records = set(['ATOM  ', 'HETATM'])
    chainID = chainID.upper()
    seq_resn = []
    seq_resi = []
    chains = set()
    read = set()
    with open(fpdb, "r") as f:
        for line in f:
            line = line.strip()
            if line[0:6] in records:
                resn = line[17:20]
                chain = line[21]
                resi = line[22:26]
                icode = line[26]
                r_uid = (resn, chain, resi, icode)
                chains.add(chain)
                if chain == chainID:
                    if r_uid not in read:
                        read.add(r_uid)
                    else:
                        continue
                    aa_resn = res_code_3to1.get(resn, 'X')
                    seq_resn.append(aa_resn)
                    seq_resi.append(resi.strip())
        if chainID not in chains:
            raise ValueError(
                "Chain `{}` NOT exist in PDB file '{}'".format(chainID, fpdb))

    return np.array(seq_resn), np.array(seq_resi)


def get_pssm(fpssm):
    """Get the content of PSSM file.

    Arguments:
        fpssm {str} -- input pssm file

    Raises:
        ValueError -- the line with number of columns not equal to 44

    Returns:
        [numpy 2D array] -- pssm content
    """
    rule = tuple(map(str, range(10)))
    pssm = []
    with open(fpssm, "r") as f:
        for line in f.readlines():
            line_raw = line
            line = line.strip()
            # only select lines that contain pssm values
            if line.startswith(rule):
                    # TODO parse pssm based on column index
                    # normal PSSM line have 44 columns.
                    # Abnormal <44 due to lakcing of gap between numbers.
                if len(line.split()) == 44:
                    pssm.append(line.split())
                else:
                    raise ValueError(
                        "Wrong format of the following line in PSSM file {}:\n{}".format(fpssm, line_raw))

    return np.array(pssm)


def get_aligned_sequences(seq1, seq2):
    """Align two sequnces using global alignment and return aligned sequences.
        Paramters of global alignment:
            match: 1
            mismtach: 0
            gap open: -2
            gap extend: -1

    Arguments:
        seq1 {str} -- 1st sequence.
        seq2 {str} -- 2nd sequence.

    Returns:
        [numpy vector] -- seq1_ali, aligned sequence for seq1
        [numpy vector] -- seq2_ali, aligned sequence for seq2
    """
    gap_open_cost = -2
    gap_extend_cost = -1
    ali = pairwise2.align.globalxs(seq1, seq2, gap_open_cost, gap_extend_cost)
    # ali[0][0] and alip[0][1] are strings of aligned sequences.
    seq1_ali = np.array([i for i in ali[0][0]])
    seq2_ali = np.array([i for i in ali[0][1]])

    return seq1_ali, seq2_ali


def get_consistent_sequence_indexes(pssm_seq_align, pdb_seq_align):
    """Get indexes to output PDB and PSSM that contain the consistent sequence
    and also the index of mutations.

    Arguments:
        pssm_seq_align {numpy vector} -- PSSM aligned sequence
        pdb_seq_align {numpy vector} -- PDB aligned sequence

    Returns:
        [numpy boolean vector] -- index to output PSSM with consistent sequence
        [numpy boolean vector] -- index to output PDB with consistent sequence
        [numpy boolean vector] -- index of mutations in consistent sequence
    """
    # get indexes for matched and mismatched residues.
    index_match = pdb_seq_align == pssm_seq_align
    index_mismatch = np.logical_not(index_match)

    # make a gap sequence (only "-") and X sequence (only "X")
    # that have same length as pdb/pssm_seq_align
    seqlen = len(pdb_seq_align)
    gap_seq = np.array(["-"] * seqlen)
    resX_seq = np.array(["X"] * seqlen)
    # get index of gap and residue X
    index_gappdb = gap_seq == pdb_seq_align
    index_resXpdb = resX_seq == pdb_seq_align
    index_gappssm = gap_seq == pssm_seq_align
    index_resXpssm = resX_seq == pssm_seq_align
    # get index of normal residues (not gap, not res X) for each sequence
    index_norm_pdb = np.logical_not(np.logical_or(index_gappdb, index_resXpdb))
    index_norm_pssm = np.logical_not(
        np.logical_or(index_gappssm, index_resXpssm))
    # get index of normal residues for both sequences
    index_norm_both = np.logical_and(index_norm_pdb, index_norm_pssm)

    # get index of match and mismatch with both norm
    index_match_norm = np.logical_and(index_match, index_norm_both)
    index_mismatch_norm = index_mut = np.logical_and(
        index_mismatch, index_norm_both)

    index_consistent = np.logical_or(index_match_norm, index_mismatch_norm)

    index_pssm_out = index_consistent[np.logical_not(index_gappssm)]
    index_pdb_out = index_consistent[np.logical_not(index_gappdb)]

    return index_pssm_out, index_pdb_out, index_mut


def display_mutations(pssmname, pdbname, chainID,
                      pssm_seq_align, pdb_seq_align, index_mut):
    """Display mutations warning in stdout if there are mutations.

    Arguments:
        pssmname {str} -- PSSM input file name
        pdbname {str} -- PDB input file name
        chainID {str} -- target chain ID
        pssm_seq_align {numpy vector} -- PSSM aligned sequence
        pdb_seq_align {numpy vector} -- PDB aligned sequence
        index_mut {numpy boolean vector} -- index of mutations in consistent sequence

    Raises:
        warning -- mutations exist in consistent sequence.
    """
    if len(set(index_mut)) > 1:
        mut_seq = []
        for i in index_mut:
            if i:
                mut_seq.append("^")
            else:
                mut_seq.append("_")
        try:
            raise Warning("Warning: Mutations exist in following sequences:\n"
                          ">{pdbname}_{chainid}:\n{pdbseq}\n"
                          ">{pssmname}:\n{pssmseq}\n{mutseq}\n".format(
                              pdbname=pdbname, chainid=chainID, pssmname=pssmname,
                              pdbseq="".join(pdb_seq_align), pssmseq="".join(pssm_seq_align),
                              mutseq="".join(mut_seq)
                          ))
        except Warning as e:
            print(e)


def write_consistent_pssm(pssm, pdbname, chainID, outdir,
                          index_pssm_out, index_pdb_out, pdb_seq_resn, pdb_seq_resi):
    """Write PSSM with consistent sequence to file.

    Arguments:
        pssm {numpy array} -- input PSSM content
        pdbname {str} -- PDB input file name
        chainID {str} -- target chain ID
        outdir {str} -- output path
        index_pssm_out {numpy boolean vector} -- index to output PSSM
        index_pdb_out {numpy boolean vector} -- index to output PDB
        pdb_seq_resn {numpy vector} -- input PDB sequence residue names
        pdb_seq_resi {numpy vector} -- input PDB sequence residue IDs
    """
    # Add the residue number and name of PDB file to the mapped pssm
    # For pssm,  only keep the score matrix and information content
    header = np.array(["pdbresi", "pdbresn", "seqresi", "seqresn", "A", "R",
                       "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
                       "F", "P", "S", "T", "W", "Y", "V", "IC"])

    resi_out = pdb_seq_resi[index_pdb_out]
    resn_out = pdb_seq_resn[index_pdb_out]

    pssm_out = pssm[index_pssm_out]
    pssm_out = np.concatenate((resi_out.reshape(resi_out.shape[0], 1), resn_out.reshape(
        resn_out.shape[0], 1), pssm_out[:, :22], pssm_out[:, -2:-1]), axis=1)
    pssm_out = np.concatenate((header.reshape(1, header.shape[0]), pssm_out))

    fopssm = os.path.join(outdir, os.path.splitext(
        pdbname)[0] + "." + chainID.upper() + ".pdb.pssm")
    with open(fopssm, "w") as f:
        for i in pssm_out:
            tmp1 = ["{:>7s}".format(j) for j in i[:4]]
            tmp2 = ["{:>4s}".format(j) for j in i[4:]]
            f.write(" ".join(tmp1+tmp2) + "\n")
    print (fopssm + " generated.")


def write_consistent_pdb(fipdb, pdbname, chainID, outdir, pdb_seq_resi, index_pdb_out):
    """Write PDB with consistent sequence of target chain to file.

    Arguments:
        fipdb {str} -- input PDB file
        pdbname {str} -- input PDB file name used for output
        chainID {str} -- target chain
        outdir {str} -- output path
        pdb_seq_resi {numpy vector} -- input PDB sequence residue IDs
        index_pdb_out {numpy boolean vector} -- index to output PDB with consistent sequence
    """
    fopdb = os.path.join(outdir, os.path.splitext(pdbname)[
                         0] + "." + chainID.upper() + ".pssm.pdb")
    fout = open(fopdb, "w")

    resi_out = pdb_seq_resi[index_pdb_out]
    _records = set(['ATOM  ', 'HETATM'])
    with open(fipdb, "r") as f:
        for line in f:
            if line[0:6] in _records:
                chain = line[21]
                resi = line[22:26].strip()
                if chain == chainID and resi in resi_out:
                    fout.write(line)
                else:
                    continue
            elif line[0:6] == "REMARK":
                fout.write(line)
            else:
                continue
    fout.close()
    print (fopdb + " generated.")


def main(fpssm, fpdb, chainID, outdir):
    """Map PDB sequence to PSSM sequence to get the consistent sequence,
    and output mapped PSSM and/or PDB file with consistent sequence.

    Arguments:
        fpssm {str} -- input PSSM file
        fpdb {str} -- input PDB file
        chainID {str} -- target chain of PDB
        outdir {str} -- output path

    Returns:
        mapped PSSM file:  pdbfilename.chainID.pdb.pssm
        mapped PDB file:   pdbfilename.chainID.pssm.pdb
    """
    # get pssm and pdb file name
    pdbname = os.path.basename(fpdb)
    pssmname = os.path.basename(fpssm)
    # get pdb sequence and residue numbers
    pdb_seq_resn, pdb_seq_resi = get_pdb_seq(fpdb, chainID)
    pdb_seq_str = "".join(pdb_seq_resn)
    # get pssm content and sequnce
    pssm = get_pssm(fpssm)
    pssm_seq_str = "".join(pssm[:, 1])

    # get aligned seqeuences
    pssm_seq_align, pdb_seq_align = get_aligned_sequences(
        pssm_seq_str, pdb_seq_str)

    # get consistent boolen index for output
    index_pssm_out, index_pdb_out, index_mut = get_consistent_sequence_indexes(
        pssm_seq_align, pdb_seq_align)

    # dispaly mutations
    display_mutations(pssmname, pdbname, chainID,
                      pssm_seq_align, pdb_seq_align, index_mut)

    # write consistent pssm and pdb
    write_consistent_pssm(pssm, pdbname, chainID, outdir,
                          index_pssm_out, index_pdb_out, pdb_seq_resn, pdb_seq_resi)
    write_consistent_pdb(fpdb, pdbname, chainID, outdir,
                         pdb_seq_resi, index_pdb_out)


if __name__ == "__main__":
    check_input(sys.argv[1:])
    fpssm, fpdb, chainID, outdir = sys.argv[1:]
    main(fpssm, fpdb, chainID, outdir)
