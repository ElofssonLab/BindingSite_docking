#!/usr/bin/env python3:q
import os
import sys
import math
import random
import argparse
import statistics
import numpy as np
import pandas as pd
import seaborn as sb
from sklearn.metrics import auc
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

def compute_metrics(TP, FP, FN, TN):
    if  (TP != 0 or FN != 0)\
    and (FP != 0 or TN != 0)\
    and (TP != 0 or FP != 0):
        PPV = TP/(TP+FP)
        TPR = TP/(TP+FN)
        FPR = FP/(FP+TN)
    
        if TPR != 0 or PPV != 0: F1 = 2*(PPV*TPR)/(PPV+TPR)
        else: F1 = 0.0
    
        if TN != 0 or FN != 0:
            num = (TP*TN)-(FP*FN)
            den = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
            MCC = num/math.sqrt(den)
        else: MCC = 0.0

        return TPR, FPR, PPV, MCC, F1
    else: 
        return 0.0, 0.0, 0.0, 0.0, 0.0

def auc_stats(df):
    TP = FP = FN = TN = 0
    TPRs = [0.0]
    FPRs = [0.0]
    PPVs = [0.0]
    MCCs = [0.0]
    F1s = [0.0]
    thrs = [1.0]

    for index, row in df.iterrows():
        if row[1] == 1.0: FN += 1
        elif row[1] == 0.0: TN += 1
    TPR, FPR, PPV, MCC, F1 = compute_metrics(TP, FP, FN, TN)
    TPRs.append(TPR)
    FPRs.append(FPR)
    PPVs.append(PPV)
    MCCs.append(MCC)
    F1s.append(F1)
    thrs.append(row[0])

    for index, row in df.iterrows():
        if row[1] == 1.0: 
            FN -= 1
            TP += 1
        elif row[1] == 0.0: 
            TN -= 1
            FP += 1
        TPR, FPR, PPV, MCC, F1 = compute_metrics(TP, FP, FN, TN)
        TPRs.append(TPR)
        FPRs.append(FPR)
        PPVs.append(PPV)
        MCCs.append(MCC)
        F1s.append(F1)
        thrs.append(row[0])

    return TPRs, FPRs, PPVs, MCCs, F1s, thrs

def cutoff_metrics(df, thr):
    TP = FP = FN = TN = 0
    
    for index, row in df.iterrows(): 
        if row[0] > thr:
            if row[1] == 1.0: TP += 1
            if row[1] == 0.0: FP += 1
        if row[0] <= thr:
            if row[1] == 1.0: FN += 1
            if row[1] == 0.0: TN += 1
    TPR, FPR, PPV, MCC, F1 = compute_metrics(TP, FP, FN, TN)

    return TP, TPR, PPV

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Analyze ISPRED4 predictions")
    parser.add_argument("-l", type=str, help="codes list")
    parser.add_argument("-p", type=str, help="ISPRED4 file folder")
    parser.add_argument("-t", type=str, help="ground truth ISPRED4 files folder")
    ns = parser.parse_args()

    allscores = []
    metrics_all = {}
    metrics_worst = {}
    codelist = [line.rstrip() for line in open(ns.l)]
    for code in codelist:
        pdbcode = code[:4]
        isp = open(ns.p+code+'.ispred4.tsv')
        ist = open(ns.t+code+'.ispred4.tsv')
        scores = [[float(lp.split('\t')[-1]),float(lt.split('\t')[-1])]\
                  for lp, lt in zip(isp,ist) if lt.split('\t')[-2] != 0.0\
                  and lt.split('\t')[-2] != 'None' and float(lt.split('\t')[-1]) != 0.5]

        scores = pd.DataFrame(scores)
        scores = scores.sort_values(0, ascending=False)
        TP, TPR, PPV = cutoff_metrics(scores, 0.9)
        metrics_all[code] = [PPV, TPR, TP]
        if not pdbcode in metrics_worst: 
            metrics_worst[pdbcode] = [PPV, TPR, TP]
        if metrics_worst[pdbcode][0] > PPV: 
            metrics_worst[pdbcode] = [PPV, TPR, TP]

    sb.set_style("whitegrid")
    fig, axes = plt.subplots(1, 1, figsize=(15, 15))
    dic = {'codes':[], 'PPV':[], 'TPR':[], 'TP':[], 'kind':[]}
    for code in metrics_all:
        dic['codes'].append(code)
        dic['PPV'].append(metrics_all[code][0])
        dic['TPR'].append(metrics_all[code][1])
        dic['kind'].append('single_chain')
    for code in metrics_worst:
        dic['codes'].append(code)
        dic['PPV'].append(metrics_worst[code][0])
        dic['TPR'].append(metrics_worst[code][1])
        dic['kind'].append('dimer_worst')
    sb.scatterplot(x='PPV', y='TPR', hue='kind', data=dic, palette=['orange', 'purple'], s=30, ax=axes)
    axes.set_xlim(0,1.01)
    axes.set_ylim(0,1.01)
    axes.set_xlabel("Precision")
    axes.set_ylabel("Recall")
    fig.savefig('MCC-F1-TPR-PPV.png')
    
    for code in metrics_worst: 
        if metrics_worst[code][0] > 0.25: 
            code1 = code+'_u1'
            code2 = code+'_u2'
            print (code1, metrics_all[code1], code2, metrics_all[code2])
