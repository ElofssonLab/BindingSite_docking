import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import sys

def main():
    codespr = {}
    codespp = {}
    codesgr = {}
    codesgp = {}
    for line in open('psite_top10_real'):
        codespr[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('psite_top10_pred'):
        codespp[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('gsite_top10_real'):
        codesgr[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('gsite_top10_pred'):
        codesgp[line.split()[0]] = float(line.split()[1].rstrip())
    
    
    set1 = set(codespr.keys())
    set2 = set(codespp.keys())
    set3 = set(codesgr.keys())
    set4 = set(codesgp.keys())
    common = list(set1.intersection(set2).intersection(set3))
    dic = {'pyros_real':[], 'pyros_pred':[], 'gramm_real':[], 'gramm_pred':[]}
    for code in common:
        dic['pyros_real'].append(codespr[code])
        dic['pyros_pred'].append(codespp[code])
        dic['gramm_real'].append(codesgr[code])
        dic['gramm_pred'].append(codesgp[code])
    
    dic = pd.DataFrame(dic)
    
    ax = sb.scatterplot(x='pyros_real', y='pyros_pred', data=dic)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.show()
    
    ax = sb.scatterplot(x='gramm_real', y='gramm_pred', data=dic)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.show()
    
    sb.scatterplot(x='pyros_real', y='gramm_real', data=dic)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.show()

def findcomm(list1, list2):
    dic1 = {}
    dic2 = {}
    for line in open(list1):
        dic1[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open(list2):
        dic2[line.split()[0]] = float(line.split()[1].rstrip())

    set1 = set(dic1.keys())
    set2 = set(dic2.keys())
    common = list(set1.intersection(set2))
    for code in common: print (code)

findcomm(sys.argv[1], sys.argv[2])
