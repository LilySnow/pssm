#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Write node file with PSSM attributes.
For mutation positions, PSSM is the differnence of PSSM between mutant and wildtype residue.
For wildtype positions, PSSM is the PSSM of wildtype residue.

usage: python write_node_pssm.py <pssmA file> <pssmB file> <node file> <chid> <resi> <wildtype resn> <mutant resn>
!NOTE: PSSM files must be in order of "<pssmA file> <pssmB file>".
example: python write_node_pssm.py 1BRS_chain_A.pssm  1BRS_chain_B.pssm 1BRS_A_ARG57ALA_node.tsv A 57 ARG ALA
output:
    #index    node    PSSM    PSSMmax    PSSMic
    1    ALA:A:57:C    -6    5    0.43
    2    ARG:A:57:C    -6    5    0.43
    3    ARG:A:57:CA    -6    5    0.43

Author: {0} ({1})
Creation Time: 2016-09-23 12:06:07 CEST
Last Modified: 2016-11-15 16:00:20 CET
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
    if len(args) != 7:
        sys.stderr.write(USAGE)
        sys.exit(1)


def read_node(fnode, ncol=2):
    """read node flile into a dict.

    input:
    #index    node
    1    ALA:A:57:C
    2    ARG:A:57:C
    3    ARG:A:57:CA
    4    ALA:A:57:CA
    """
    dict_node = {}
    p = re.compile(r'^[0-9]')
    with open(fnode) as f:
        lines = [ i.strip(" \n") for i in f ]
    for i in lines:
        if p.match(i):
            line = i.split()
            if len(line) == ncol:
                node = {}
                nodeatt = line[1].split(":")
                # ['ALA', 'A', '57', 'C']
                node["node"] = line[1]
                node["resn"] = nodeatt[0]
                node["chid"] = nodeatt[1]
                node["resi"] = nodeatt[2]
                node["atom"] = nodeatt[3]
                index = int(line[0])
                dict_node[index] = node
            else:
                print "\nThe number of columns is not {} in following lines:".format(ncol)
                print i
                sys.exit("Program quit!\n")
    return dict_node


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

def print_output(fpssmA, fpssmB, fnode, chid, resi, resnwt, resnmut=""):
    """
    formatted output
    """
    pssm = {}
    pssm["A"] = read_pssm(fpssmA)
    pssm["B"] = read_pssm(fpssmB)
    nodes = read_node(fnode)
    resnwt = resn_3to1(resnwt)
    resnmut = resn_3to1(resnmut)

    header = ["#index", "node", "PSSM", "PSSMmax", "PSSMic" ]
    print "\t".join(header)
    for i,j in nodes.iteritems():
        # i is index
        # j is node, chid, resi, resi, atom
        if (j['chid'] == chid) and (j['resi'] == resi):
            # calcualte dPSSM = PSSMresmut - PSSMreswt
            dpssm = int(pssm[chid][resi][resnmut]) - int(pssm[chid][resi][resnwt])
            out = "{}\t{}\t{}\t{}\t{:.2f}".format(i, j["node"], dpssm, pssm[chid][resi]['max'], float(pssm[chid][resi]['IC']))
        elif (j['chid'] != chid) or (j['resi'] != resi):
            # calculate PSSM = PSSMreswt
            chid1 = j["chid"]
            resi1 = j["resi"]
            resn1 = resn_3to1(j["resn"])
            pssm1 = pssm[chid1][resi1][resn1]
            out = "{}\t{}\t{}\t{}\t{:.2f}".format(i, j["node"], pssm1, pssm[chid1][resi1]['max'], float(pssm[chid1][resi1]['IC']))
        print out


if __name__ == '__main__':
    check_input(sys.argv[1:])
    print_output(*sys.argv[1:])

    # fA = sys.argv[1]
    # fB = sys.argv[2]
    # fnode = sys.argv[3]
    # chid = sys.argv[4]
    # resi = sys.argv[5]
    # resnwt = sys.argv[6]
    # resnmut = sys.argv[7]
    # print_output(fA, fB, fnode, chid, resi, resnwt, resnmut)
