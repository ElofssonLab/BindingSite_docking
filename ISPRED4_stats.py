#!/usr/bin/env python3:q
import os
import sys
import math
import random
import argparse
import statistics
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import multiprocessing as mp
import numpy as np

from sklearn.metrics import auc
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
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
    PPVs = []
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

        if PPVs == []: PPVs.append(TP/(TP+FP))
        prev = el[1]

    if  (TP != 0 or FN != 0)\
    and (FP != 0 or TN != 0)\
    and (TP != 0 or FP != 0):
        TPRs.append(TP/(TP+FN))
        FPRs.append(FP/(FP+TN))
        PPVs.append(TP/(TP+FP))

    return TPRs, FPRs, PPVs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Analyze ISPRED4 predictions")
    parser.add_argument("-l", type=str, help="codes list")
    parser.add_argument("-p", type=str, help="ISPRED4 file folder")
    parser.add_argument("-t", type=str, help="ground truth ISPRED4 files folder")
    parser.add_argument("-s", type=str, help="structures files folder")
    ns = parser.parse_args()

    count = 0
    metrics = {}
    allscores = []
    codelist = [line.rstrip() for line in open(ns.l)]
    p = PDBParser(QUIET=True)
    for code in codelist:
        splitscores = []
        pdbcode = code.rstrip('_u12')
        isp = open(ns.p+code+'.ispred4.tsv')
        ist = open(ns.t+code+'.ispred4.tsv')
        pdb = p.get_structure('', ns.s+code+'.pdb')
        for c in pdb[0]: chain = c.get_id() 
        dssp = DSSP(pdb[0], ns.s+code+'.pdb')
        rasa = [line[3] for line in dssp]

        splitscores = [[lp.split('\t')[0].split('_')[0], float(lp.split('\t')[-1]),\
                        float(lt.split('\t')[-1])] for lp, lt, r in zip(isp,ist,rasa)\
                       if r >= 0.2]

        allscores += splitscores[:]
        sortedl = mergesort_pred(splitscores)
        TPRs, FPRs, PPVs = auc_stats(sortedl)
        if len(TPRs) > 2:
            PRauc = auc(TPRs, PPVs)
            ROCauc = auc(FPRs, TPRs)
            topPPV = [sum([1 for el in sortedl[:int(lbl)] if el[-1] == 1.0])/float(lbl)\
                      for lbl in ['5', '10', '20']]
            if not pdbcode in metrics: metrics[pdbcode] = [PRauc, ROCauc, TPRs, FPRs, PPVs] + topPPV
            elif metrics[pdbcode][0] > PRauc: 
                metrics[pdbcode] = [PRauc, ROCauc, TPRs, FPRs, PPVs] + topPPV

    sortedl = mergesort_pred(allscores)
    TPRs, FPRs, PPVs = auc_stats(sortedl)

    df = pd.DataFrame(metrics).T
    best20 = df.sort_values(1, ascending=False)[:10]
    worst20 = df.sort_values(1, ascending=True)[:10]
    b20 = list(best20.index)
    w20 = list(worst20.index)

    fs = 14
    bins = np.arange(0, 1.1, 0.1) 
    ##### all ROC curves
    fig, axes = plt.subplots(1, 2, figsize=(25,10))
    sb.lineplot(FPRs, TPRs, ci=None, linewidth=1.2, color='purple', ax=axes[0])
    for key in b20: sb.lineplot(metrics[key][3], metrics[key][2], ci=None, linewidth=0.6, color='blue', ax=axes[0])
    for key in w20: sb.lineplot(metrics[key][3], metrics[key][2], ci=None, linewidth=0.6, color='orange', ax=axes[0])
    axes[0].set_xlabel("1 - Specificity")
    axes[0].set_ylabel("Recall")
    axes[0].set_xlim(0,1)
    axes[0].set_ylim(0,1)

    rocs = [metrics[key][1] for key in metrics]
    sb.distplot(rocs, kde=False, bins=bins, ax=axes[1])
    mean = statistics.mean(rocs)
    median = statistics.median(rocs)
    #axes[1].vlines(mean, 0, 50)
    axes[1].vlines(median, 0, 70)
    axes[1].set_xlabel("ROC curve AUCs")
    axes[1].set_xlim(0,1)
    axes[1].set_ylim(0,70)
    fig.savefig('allROCs.png')

    ##### all PR curves

    best20 = df.sort_values(0, ascending=False)[:10]
    worst20 = df.sort_values(0, ascending=True)[:10]
    b20 = list(best20.index)
    w20 = list(worst20.index)

    fig, axes = plt.subplots(1, 2, figsize=(25,10)) 
    sb.lineplot(TPRs, PPVs, ci=None, linewidth=1.2, color='purple', ax=axes[0])
    for key in b20: sb.lineplot(metrics[key][2], metrics[key][4], ci=None, linewidth=0.6, color='blue', ax=axes[0])
    for key in w20: sb.lineplot(metrics[key][2], metrics[key][4], ci=None, linewidth=0.6, color='orange', ax=axes[0])
    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    axes[0].set_xlim(0,1)
    axes[0].set_ylim(0,1)

    prs = [metrics[key][0] for key in metrics]
    sb.distplot(prs, kde=False, bins=bins, ax=axes[1])
    mean = statistics.mean(prs)
    median = statistics.median(prs)
    #axes[1].vlines(mean, 0, 50)
    axes[1].vlines(median, 0, 100)
    axes[1].set_xlabel("PR curve AUCs")
    axes[1].set_xlim(0,1)
    axes[1].set_ylim(0,100)
    fig.savefig('allPRs.png')

#    ##### top PPVs
#    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(25,10))
#    top5 = [ metrics[key][5] for key in metrics ]
#    top10 = [ metrics[key][6] for key in metrics ]
#    top20 = [ metrics[key][7] for key in metrics ]
#    sb.distplot(top5, kde=False, ax=axes[0])
#    sb.distplot(top10, kde=False, ax=axes[1])
#    sb.distplot(top20, kde=False, ax=axes[2])
#    axes[0].set_xlabel('Top5')
#    axes[0].set_xlim(0,1)
#    axes[1].set_xlabel('Top10')
#    axes[1].set_xlim(0,1)
#    axes[2].set_xlabel('Top20')
#    axes[2].set_xlim(0,1)
#    fig.text(0.513, 0.02, 'Precision', ha='center')
#    fig.savefig('topPPVs.png')


