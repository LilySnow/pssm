#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Write sum and average information content for specific residues.

usage: python write_pssmic_residues.py <pssm file> <residueID list file>
example: python write_pssmic_residues.py 1FT4.pssm resID.list
output: the name of pssm file, the sum and average information content.
        e.g. 1FT4.pssm    5.02    1.67

Author: {0} ({1})
Creation Time: 2016-09-23 12:06:07 CEST
Last Modified: 2016-11-10 12:29:18 CET
"""

import os
import sys
import re
import numpy as np

try:
    import requests
except ImportError as e:
    print >> sys.stderr, "Can not import 'requests' Python module."
    sys.exit(1)

__author__ = "Cunliang Geng"
__email__ = "gengcunliang@gmail.com"

USAGE = __doc__.format(__author__, __email__)

def check_input(args):
    if len(args) != 2:
        sys.stderr.write(USAGE)
        sys.exit(1)


def read_pssm(fnpssm, ncol=44):
    """read pssm values and InfoContent into a dict,
    the key is resi and the value is a dict contaning residue pssm and IC.
    Note taht default number of columns ncol=44.
    e.g.
    dict_pssm['1'] = {'resi': '1', 'resn': 'F', 'A': '-3', 'C': '-3', 'E': '-4', 'D': '-5', 'G': '-4', 'F': '7',
    'I': '-1', 'H': '-2', 'K': '-4', 'M': '-1', 'L': '1', 'N': '-4', 'Q': '-4', 'P': '-5', 'S': '-3', 'R': '-4', 'T': '-3', 'W': '0', 'V': '-2', 'Y': '4', 'IC': '1.25'}
    """
    # the resheader must match the residue order in pssm file header
    resheader = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    ncol = ncol
    dict_pssm = {}
    p = re.compile(r'^[0-9]')
    with open(fnpssm) as f:
        lines = [ i.strip(" \n") for i in f ]
    for i in lines:
        if p.match(i):
            lt_line = i.split()
            if len(lt_line) == ncol:
                pssml = {}
                for j in range(0,20):
                    pssml[resheader[j]] = lt_line[j+2]
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

def read_iresi(firesi):
    """read interface residue ID into a list"""
    with open(firesi) as f:
        # lines = [ int(i.strip(" \n")) for i in f ]
        lines = [ i.strip(" \n") for i in f ]
    return lines

def is_non_zero_file(fpath):
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

if __name__ == '__main__':
    check_input(sys.argv[1:])
    fin = sys.argv[1]
    firesi = sys.argv[2]
    pssm = read_pssm(fin)
    # print pssm[resi], reswt, resmut
    iresi = read_iresi(firesi)
    # print iresi
    pssmic = []
    if is_non_zero_file(firesi):
        for i in iresi:
            if i!="":
                print pssm[i]
                pssmic.append(float(pssm[i]['IC']))
        out = "{}\t{:.2f}\t{:.2f}".format(fin, np.sum(pssmic), np.mean(pssmic))
    else:
        out = "{}\t0\t0".format(fin)

    # filename, ICsum, ICmean
    print out
