#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract PSSM values from PSSM file for specific position.

usage: python extract_pssm.py <pssmA file> <pssmB file> <output name> <chid> <resi> <wildtype resn> <mutant resn>
!NOTE: PSSM files must be in order of "<pssmA file> <pssmB file>".
example: python extract_pssm.py 4CPA_A_wild.pssm 4CPA_B_wild.pssm 4CPA_A_ARG2ALA A 2 ARG ALA
output:
    #ID    A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    PSSMwt    PSSMmut    PSSMmax    PSSMic
    4CPA_A_ARG2ALA    -3    8    -2    -3    -5    0    -2    -3    -2    -5    -3    1    -3    -5    -4    -2    -2    -2    -4    -4    8    -3    8    1.60

Author: {0} ({1})
Creation Time: 2016-09-23 12:06:07 CEST
Last Modified: 2017-03-24 15:14:27 CET
"""

import os
import sys
import re

__author__ = "Cunliang Geng"
__email__ = "gengcunliang@gmail.com"

USAGE = __doc__.format(__author__,  __email__)

def check_input(args):
    """
    check number of arguments
    """
    if len(args) != 7:
        sys.stderr.write(USAGE)
        sys.exit(1)

def read_pssm(fnpssm,  ncol=44):
    """read pssm values and InfoContent into a dict,
    the key is resi and the value is a dict contaning residue pssm and IC.
    Note that default number of columns ncol=44.
    e.g.
    dict_pssm['2'] = {'resi': '2', 'resn': 'R', 'max': '8', 'A': '-3', 'C': '-5', 'E': '-2', 'D': '-3', 'G': '-3', 'F': '-5',
    'I': '-5', 'H': '-2', 'K': '1', 'M': '-3', 'L': '-3', 'N': '-2', 'Q': '0', 'P': '-4', 'S': '-2', 'R': '8', 'T': '-2',
    'W': '-2', 'V': '-4', 'Y': '-4', 'IC': '1.60', 'pssm': ['-3', '8', '-2', '-3', '-5', '0', '-2', '-3', '-2', '-5', '-3', '1', '-3', '-5', '-4', '-2', '-2', '-2', '-4', '-4']}
    """
    # the resheader must match the residue order in pssm file header
    resheader = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    ncol = ncol
    dict_pssm = {}
    p = re.compile(r'^[0-9]')
    with open(fnpssm) as f:
        lines = [i.strip(" \n") for i in f]
    for i in lines:
        if p.match(i):
            lt_line = i.split()
            if len(lt_line) == ncol:
                pssml = {}
                for j in range(0, 20):
                    pssml[resheader[j]] = lt_line[j+2]
                pssml['pssm'] = lt_line[2:22]
                pssml['max'] = max(pssml.values())
                pssml['resi'] = lt_line[0]
                pssml['resn'] = lt_line[1]
                pssml['IC'] = lt_line[ncol-2]
                dict_pssm[lt_line[0]] = pssml
            else:
                print "\nThe number of columns is not {} in following lines:".format(ncol)
                print i
                sys.exit("Program quit!\n")
    return dict_pssm

def resn_1to3(resn):
    """convert 1 letters of residue name to 3 letters"""
    code1to3={
            "R":"ARG",
            "K":"LYS",
            "D":"ASP",
            "E":"GLU",
            "Q":"GLN",
            "N":"ASN",
            "H":"HIS",
            "S":"SER",
            "T":"THR",
            "Y":"TYR",
            "C":"CYS",
            "M":"MET",
            "W":"TRP",
            "A":"ALA",
            "I":"ILE",
            "L":"LEU",
            "F":"PHE",
            "V":"VAL",
            "P":"PRO",
            "G":"GLY"
          }
    resn = resn.upper()
    if len(resn) == 1:
        if resn in code1to3:
            resn = code1to3[resn]
        else:
            sys.exit("Error in residue name: {}".format(resn))
    elif len(resn) != 3 &  len(resn.upper()) != 1:
        sys.exit("Error in residue name: {}".format(resn))
    return resn

def resn_3to1(resn):
    """convert 3 letters of residue name to 1 letter"""
    code3to1={
            "ARG":"R",
            "LYS":"K",
            "ASP":"D",
            "GLU":"E",
            "GLN":"Q",
            "ASN":"N",
            "HIS":"H",
            "SER":"S",
            "THR":"T",
            "TYR":"Y",
            "CYS":"C",
            "MET":"M",
            "TRP":"W",
            "ALA":"A",
            "ILE":"I",
            "LEU":"L",
            "PHE":"F",
            "VAL":"V",
            "PRO":"P",
            "GLY":"G"
            }
    resn = resn.upper()
    if len(resn) == 3:
        if resn in code3to1:
            resn = code3to1[resn]
        else:
            sys.exit("Error in residue name: {}".format(resn))
    elif len(resn) != 1 & len(resn) != 3:
        sys.exit("Error in residue name: {}".format(resn))
    return resn

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def print_output(fpssmA,  fpssmB, nameout,  chid, resi, resnwt, resnmut=""):
    """
    formatted output
    #ID    A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    PSSMwt    PSSMmut    PSSMmax    PSSMic
    4CPA_A_ALA2ARG    -3    8    -2    -3    -5    0    -2    -3    -2    -5    -3    1    -3    -5    -4    -2    -2    -2    -4    -4    8    -3    8    1.60
    """
    pssm = {}
    pssm["A"] = read_pssm(fpssmA)
    pssm["B"] = read_pssm(fpssmB)
    resnwt = resn_3to1(resnwt)
    resnmut = resn_3to1(resnmut)

    #  header = ["#ID", 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V', "PSSMwt", "PSSMmut", "dPSSM", "PSSMmax", "PSSMic"]
    header = ["#ID", 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V', "PSSM_wt", "PSSM_mut", "PSSM_diff", "PSSMmax", "PSSMic"]
    print "\t".join(header)
    pssmwt = pssm[chid][resi][resnwt]
    pssmmut = pssm[chid][resi][resnmut]
    dpssm = int(pssmmut)-int(pssmwt)
    #  out = "{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}".format(nameout, "\t".join(pssm[chid][resi]['pssm']), pssm[chid][resi][resnwt], pssm[chid][resi][resnmut], pssm[chid][resi]['max'], float(pssm[chid][resi]['IC']))
    out = "{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}".format(nameout, "\t".join(pssm[chid][resi]['pssm']), pssmwt, pssmmut, dpssm, pssm[chid][resi]['max'], float(pssm[chid][resi]['IC']))
    print out

if __name__ == '__main__':
    check_input(sys.argv[1:])
    print_output(*sys.argv[1:])
