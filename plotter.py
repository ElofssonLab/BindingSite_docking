import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
import seaborn as sb
import pandas as pd
import sys

def main():
    codespr = {}
    codespp = {}
    codesgr = {}
    codesgp = {}
    for line in open('pyros_real_all'):
        codespr[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('pyros_pred_0.1'):
        codespp[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('gramm_real_all'):
        codesgr[line.split()[0]] = float(line.split()[1].rstrip())
    for line in open('gramm_pred_0.1'):
        codesgp[line.split()[0]] = float(line.split()[1].rstrip())
    
    
    set1 = set(codespr.keys())
    set2 = set(codespp.keys())
    set3 = set(codesgr.keys())
    set4 = set(codesgp.keys())
    common = list(set1.intersection(set2).intersection(set3))
    dic = {'PyRosetta':[], 'pyros_pred':[], 'Gramm':[], 'gramm_pred':[]}
    for code in common:
        dic['PyRosetta'].append(codespr[code])
        dic['pyros_pred'].append(codespp[code])
        dic['Gramm'].append(codesgr[code])
        dic['gramm_pred'].append(codesgp[code])
    
    dic = pd.DataFrame(dic)
    
    #ax = sb.scatterplot(x='pyros_real', y='pyros_pred', data=dic, s=10)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    #plt.show()
    
    #ax = sb.scatterplot(x='gramm_real', y='gramm_pred', data=dic, s=10)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    #plt.show()
   
    fig, axes = plt.subplots(1, 1, figsize=(10, 10))
    sb.scatterplot(x='PyRosetta', y='Gramm', data=dic, s=20, ax=axes)
    plt.vlines(0.23, 0, 1, '#f0b326')
    plt.hlines(0.23, 0, 1, '#f0b326')
    axes.set_xlim(0,1)
    axes.set_ylim(0,1)
    fig.savefig('Gramm_PyRosetta.png')


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

main()
#findcomm(sys.argv[1], sys.argv[2])
