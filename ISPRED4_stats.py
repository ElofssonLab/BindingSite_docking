#!/usr/bin/env python3:q
import os
import sys
import math
import random
import argparse
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import multiprocessing as mp
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd

def mergesort_pred(linear):
    elem = 1
    while len(linear) > elem:
        joblist = []
        for idx in range(0, len(linear)+elem*2, elem*2):
            ida = idx+elem
            idb = idx+elem*2
            if len(linear) >= idb:
                a = linear[idx:ida]
                b = linear[ida:idb]
            elif len(linear) >= ida:
                a = linear[idx:ida]
                b = linear[ida:]
            elif len(linear) >= idx:
                a = linear[idx:]
                b = []
            else: continue
            joblist.append([a, b])

        pool = mp.Pool(processes=5)
        results = pool.map(merge, joblist)
        pool.close()
        pool.join()

        linear = [ el for res in results for el in res ]
        elem *= 2
    return linear

def merge(job):
    l = []
    l1 = job[0]
    l2 = job[1]
    p1 = p2 = 0
    while len(l) < len(l1)+len(l2):
        if p1 == len(l1): l.extend(l2[p2:])
        elif p2 == len(l2): l.extend(l1[p1:])
        elif l1[p1][1] >= l2[p2][1]:
            l.append(l1[p1])
            p1 += 1
        else:
            l.append(l2[p2])
            p2 += 1
    return l

def auc_stats(scores):

    TPRs = [0.0]
    FPRs = [0.0]
    PPVs = [0.0]
    TP = FP = FN = TN = 0
    for el in scores: 
        if el[2] == 1.0: FN += 1
        else: TN += 1

    prev = scores[0][1]
    for el in scores:
        if el[1] != prev:
            if  (TP != 0 or FN != 0)\
            and (FP != 0 or TN != 0)\
            and (TP != 0 or FP != 0): 
                TPRs.append(TP/(TP+FN))
                FPRs.append(FP/(FP+TN))
                PPVs.append(TP/(TP+FP))
        if el[2] == 1.0: 
            FN -= 1
            TP += 1
        else: 
            TN -= 1
            FP += 1
        prev = el[1]

    return TPRs, FPRs, PPVs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Analyze ISPRED4 predictions")
    parser.add_argument("-l", type=str, help="codes list")
    parser.add_argument("-p", type=str, help="ISPRED4 file folder")
    parser.add_argument("-t", type=str, help="ground truth ISPRED4 files folder")
    parser.add_argument("-s", type=str, help="structures files folder")
    ns = parser.parse_args()

    scores = []
    p = PDBParser(QUIET=True)
    codelist = [line.rstrip() for line in open(ns.l)]
    for code in codelist:
        isp = open(ns.p+code+'.ispred4.tsv')
        ist = open(ns.t+code+'.ispred4.tsv')
        pdb = p.get_structure('', ns.s+code+'.pdb')
        for c in pdb[0]: chain = c.get_id() 
        dssp = DSSP(pdb[0], ns.s+code+'.pdb')
        rasa = [line[3] for line in dssp]
        scores += [[lp.split('\t')[0].split('_')[0], float(lp.split('\t')[-1]),\
                    float(lt.split('\t')[-1])] for lp, lt, r in zip(isp,ist,rasa)\
                   if r >= 0.2]

    sorted_scores = mergesort_pred(scores)
    TPRs, FPRs, PPVs = auc_stats(sorted_scores)

    sb.lineplot(FPRs, TPRs)
    plt.show()

    sb.lineplot(TPRs, PPVs)
    plt.show()





