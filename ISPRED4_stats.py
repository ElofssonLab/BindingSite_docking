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
    TP = FP = FN = TN = 0
    TPRs = [0.0]
    FPRs = [0.0]
    PPVs = [0.0]
    thrs = []
    MCCs = []
    F1s = []

    for el in scores: 
        if el[2] == 1.0: FN += 1
        elif el[2] == 0.0: TN += 1

    prev = scores[0][1]
    for el in scores:
        if el[1] != prev:
            if  (TP != 0 or FN != 0)\
            and (FP != 0 or TN != 0)\
            and (TP != 0 or FP != 0):
                PPV = TP/(TP+FP)
                TPR = TP/(TP+FN)
                TPRs.append(TPR)
                PPVs.append(PPV)
                FPRs.append(FP/(FP+TN))
                if TPR != 0 or PPV != 0: F1s.append(2*(PPV*TPR)/(PPV+TPR))
                else: F1s.append(0.0)
                if TN != 0 or FN != 0:
                    num = (TP*TN)-(FP*FN)
                    den = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
                    MCCs.append(num/math.sqrt(den))
                else: MCCs.append(0.0)    
                thrs.append(el[1])

        if el[2] == 1.0: 
            FN -= 1
            TP += 1
        elif el[2] == 0.0: 
            TN -= 1
            FP += 1
        prev = el[1]

    if  (TP != 0 or FN != 0)\
    and (FP != 0 or TN != 0)\
    and (TP != 0 or FP != 0):
        PPV = TP/(TP+FP)
        TPR = TP/(TP+FN)
        TPRs.append(TPR)
        PPVs.append(PPV)
        FPRs.append(FP/(FP+TN))
        if PPV != 0 or TPR != 0: F1s.append(2*(PPV*TPR)/(PPV+TPR))
        else: F1s.append(0.0)
        if TN != 0 or FN != 0:
            num = (TP*TN)-(FP*FN)
            den = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
            MCCs.append(num/math.sqrt(den))
        else: MCCs.append(0.0)
        thrs.append(el[1])

    return TPRs, FPRs, PPVs, MCCs, F1s, thrs

def get_metrics(scorelist, thr):
    TP = FP = FN = TN = 0
    for el in scorelist: 
        if el[1] > thr:
            if el[2] == 1.0: TP += 1
            if el[2] == 0.0: FP += 1
        if el[1] <= thr:
            if el[2] == 1.0: FN += 1
            if el[2] == 0.0: TN += 1

    if TP != 0 or FP != 0: PPV = TP/(TP+FP) 
    else: PPV = 0
    if TP != 0 or FN != 0: TPR = TP/(TP+FN)
    else: TPR = 0
    if PPV != 0 or TPR != 0: F1 = 2*(PPV*TPR)/(PPV+TPR)
    else: F1 = 0

    if  (TP != 0 or FN != 0)\
    and (FP != 0 or TN != 0)\
    and (TP != 0 or FP != 0)\
    and (TN != 0 or FN != 0):
        num = (TP*TN)-(FP*FN)
        den = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
        MCC = num/math.sqrt(den)
    else: MCC = 0

    return F1, MCC, TPR, PPV

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
        pdbcode = code[:4]
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
        F1, MCC, TPR, PPV = get_metrics(sortedl, 0.75)
        TPRs, FPRs, PPVs, MCCs, F1s, thrs = auc_stats(sortedl)
        if len(TPRs) > 2:
            PRauc = auc(TPRs, PPVs)
            ROCauc = auc(FPRs, TPRs)

            if not pdbcode in metrics: 
                metrics[pdbcode] = [PRauc, ROCauc, TPRs, FPRs, PPVs, F1, MCC, TPR, PPV]
            if metrics[pdbcode][0] > PRauc: 
                metrics[pdbcode] = [PRauc, ROCauc, TPRs, FPRs, PPVs, F1, MCC, TPR, PPV]

    bins = np.arange(0, 1.1, 0.1)
    sortedl = mergesort_pred(allscores)
    TPRs, FPRs, PPVs, MCCs, F1s, thrs = auc_stats(sortedl)
    df = pd.DataFrame(metrics).T

    ##### ROC curves
    fig, axes = plt.subplots(1, 2, figsize=(25, 10))
    best20 = df.sort_values(1, ascending=False)[:10]
    worst20 = df.sort_values(1, ascending=True)[:10]
    b20 = list(best20.index)
    w20 = list(worst20.index)

    sb.lineplot(FPRs, TPRs, ci=None, linewidth=1.2, color='purple', ax=axes[0])
    for key in b20: sb.lineplot(metrics[key][3], metrics[key][2], ci=None, linewidth=0.6, color='blue', ax=axes[0])
    for key in w20: sb.lineplot(metrics[key][3], metrics[key][2], ci=None, linewidth=0.6, color='orange', ax=axes[0])
    axes[0].set_xlabel("1 - Specificity")
    axes[0].set_ylabel("Recall")
    axes[0].set_xlim(0,1)
    axes[0].set_ylim(0,1)

    ##### PR curves
    best20 = df.sort_values(0, ascending=False)[:10]
    worst20 = df.sort_values(0, ascending=True)[:10]
    b20 = list(best20.index)
    w20 = list(worst20.index)

    sb.lineplot(TPRs, PPVs, ci=None, linewidth=1.2, color='purple', ax=axes[1])
    for key in b20: sb.lineplot(metrics[key][2], metrics[key][4], ci=None, linewidth=0.6, color='blue', ax=axes[1])
    for key in w20: sb.lineplot(metrics[key][2], metrics[key][4], ci=None, linewidth=0.6, color='orange', ax=axes[1])
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_xlim(0,1)
    axes[1].set_ylim(0,1)
    fig.savefig('ROC-PR.png')

    ##### MCC - F1 scores
    fig, axes = plt.subplots(1, 2, figsize=(35, 10))
    sb.lineplot(thrs, TPRs[1:], ci=None, linewidth=1, color='purple', label='TPR', ax=axes[0])
    sb.lineplot(thrs, PPVs[1:], ci=None, linewidth=1, color='green', label='PPV', ax=axes[0])
    sb.lineplot(thrs, MCCs, ci=None, linewidth=1, color='blue', label='MCC', ax=axes[0])
    sb.lineplot(thrs, F1s, ci=None, linewidth=1, color='orange', label='F1', ax=axes[0])
    axes[0].set_xlabel("Cutoff")
    axes[0].set_xlim(0,1)
    axes[0].set_ylim(0,1)

    best = df.sort_values(8, ascending=False)
    dic = {'codes':[], 'score':[], 'metric':[]}
    for code in list(best.index):
        dic['codes'].append(code)
        dic['score'].append(metrics[code][7])
        dic['metric'].append('TPR')
        dic['codes'].append(code)
        dic['score'].append(metrics[code][8])
        dic['metric'].append('PPV')
    dfb = pd.DataFrame(dic)
    sb.barplot(x='codes', y='score', hue='metric', data=dfb, ax=axes[1])
    plt.xticks([])
    fig.savefig('MCC-F1-TPR-PPV.png')
    

