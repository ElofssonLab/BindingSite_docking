#!/usr/bin/env python3:q
import argparse
import random
import sys
import os

from Bio.PDB import *
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.DSSP import DSSP

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','GLX':'Q'}

def extract_structure(spath):
    seq = ''
    rasas = {}
    str_map = {}
    ##### load pdb file
    struc = p.get_structure('', spath)
    dssp = DSSP(struc[0], spath)
    ##### initialize a neighbour search toward partner chain
    if '_uc' in spath or '_bc' in spath:
        opp_chain = unfold_entities(struc[0][oppchain], 'A')
        ns = NeighborSearch(opp_chain)

    ##### join structure to relative accessibility values 
    ##### and neighbours presence in the thr. of 8A (1-yes, 0-no)
    keylist = []
    for res in struc[0][chain]:
        resnum = res.get_id()[1]
        resname = three2one[res.get_resname()]
        try: rasa = round(dssp[(chain, resnum)][3], 3)
        except: rasa = None
        seq += resname

        if '_uc' in spath or '_bc' in spath:
            neigh = 0
            for atom in res: 
                if ns.search(atom.get_coord(), 8) != []:
                    neigh = 1
                    break

            if type(rasa) == float: str_map[resnum] = [res, rasa, neigh]
            else: 
                print ('Pos {} missing from {}'.format(resnum, spath))
                str_map[resnum] = [res, 0.0, neigh]
        else: 
            if type(rasa) == float: str_map[resnum] = [res, rasa]
            else:
                print ('Pos {} missing from {}'.format(resnum, spath))
                str_map[resnum] = [res, 0.0]

        keylist.append(resnum)

    return str_map, seq, keylist

def interface_extract():
    ##### load structural data
    str_map = {}
    str_mapu, sequ, keysu = extract_structure(pathu)
    str_mapb, seqb, keysb = extract_structure(pathb)
    str_mapuc, sequc, keysuc = extract_structure(pathuc)
    str_mapbc, seqbc, keysbc = extract_structure(pathbc)

    ##### iterate over all aligned positions
    for pos in range(1, max(keysu)+1):
        ##### if position is mapped to bound structure
        if pos in str_mapb: 
            ##### mark residue as interface if rasa monomer > rasa dimer
            ##### or if its rasa is >= 0.2 and 
            ##### it is closer than thr to partner chain
            if str_mapb[pos][1] > str_mapbc[pos][1]:#\
            #or str_mapbc[pos][2] == 1 and str_mapbc[pos][1] >= 0.2:
                str_map[pos] = str_mapu[pos]+[1.0]
            else:
                str_map[pos] = str_mapu[pos]+[0.0]
        ##### same with unbound structure if position has no 
        ##### bound structure match
        else:
            if str_mapu[pos][1] > str_mapuc[pos][1]:#\
            #or str_mapuc[pos][2] == 1 and str_mapuc[pos][1] >= 0.2:
                str_map[pos] = str_mapu[pos]+[1.0]
            else:
                str_map[pos] = str_mapu[pos]+[0.0]

    ##### convert residue object to one letter aminoacid code
    for key in str_map:
        resname = three2one[str_map[key][0].get_resname()]
        str_map[key][0] = resname

    return str_map, sequ, keysu

def randomize_interfaces(tpr, ppv, intmap, keylist):
    randomized = [[intmap[pos][1], intmap[pos][2]] for pos in keylist]

    ##### compute number of real binding site (BS) residues
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false non-BS we need for desired TPR
    FNneeded = int(positives*(1-tpr))
    ##### mark random BS residues to be turned in false non-BS
    while FNneeded > 0:
        for p, (acc, score) in enumerate(randomized):
            if score == 1.0 and random.randint(0,1) == 1:
                randomized[p][1] = 0.75
                FNneeded -= 1
            if FNneeded <= 0: break

    ##### compute number of remaining residues correctly marked as BS
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false BS we need for desired PPV
    FPneeded = int(((positives/(ppv*100))*100)-positives)
    ##### mark random true non-BS residues to be turned in false BS
    while FPneeded > 0:
        nonegative = True
        for p, (acc, score) in enumerate(randomized):
            if score == 0.0 and acc > 0.2 : nonegative = False
            if score == 0.0 and acc > 0.2 and random.randint(0,1) == 1:
                randomized[p][1] = 0.25
                FPneeded -= 1
            if FPneeded <= 0: break
        if nonegative: break

    ##### turn all marked residues in the relative class
    for p, (acc, score) in enumerate(randomized):
        if score == 0.25: randomized[p][1] = 1.0
        if score == 0.75: randomized[p][1] = 0.0
    
    randomized = {pos:randomized[idx] for idx, pos in enumerate(keylist)}
    return randomized

def format_BIPSPI(bip):
    seq = ''
    predmap = {}
    maxval = 9.159
    struc = p.get_structure('', pathu)
    for line in open(bip):
        if line.startswith('chainId'): continue
        n = int(line.split()[1])
        res = struc[0][chain][n].get_resname()
        res = three2one[res]
        score = str(float(line.split()[-1])/maxval)
        predmap[n] = [res, score]
    resnumb = max(list(predmap.keys()))
    for res in struc[0][chain]:
        n = res.get_id()[1]
        resname = res.get_resname()
        resname = three2one[resname]
        if n in predmap: seq += predmap[n][0]
        else:
            predmap[n] = [resname, 0.0]
            seq += predmap[n][0]
    return predmap, seq

def format_ISPRED4(crf=False):
    seq = ''
    predmap = {}
    for n, line in enumerate(open(isp)):
        if crf:
            if line.split()[-2] == 'I': score = 1.0
            else: score = 0.0
            predmap[n+1] = [line.split()[3], score]
        else:
            predmap[n+1] = [line.split()[3], line.split()[-1]]
        seq += line.split()[3]
    return predmap, seq

def format_predus():
    seq = ''
    predmap = {}
    struc = p.get_structure('', pathu)
    for line in open(pre):
        for n in line.split(): 
            if n == 'PredUS': continue
            n = int(n)
            res = struc[0][chain][n].get_resname()
            res = three2one[res]
            predmap[n] = [res, '1.0']
    for res in struc[0][chain]:
        n = res.get_id()[1]
        if n not in predmap:
            res = struc[0][chain][n].get_resname()
            res = three2one[res]
            predmap[n] = [res, '0.0']
    for n in range(1, max(predmap.keys())+1): seq += predmap[n][0]
    return predmap, seq

def format_dynJET2(trace=False):
    seq = ''
    count = 1
    predmap = {}
    for line in open(dyn):
        if not line.startswith('AA'):
            resname = three2one[line.split()[0]]
            if not trace:
                score = max(float(line.split()[12]),
                            float(line.split()[13]),
                            float(line.split()[14]))
                if line.split()[15] == '1.0':
                    predmap[count] = [resname, score]
                else:
                    predmap[count] = [resname, 0.0]
            else:
                score = float(line.split()[11])
                predmap[count] = [resname, score]
            seq += resname
            count += 1
    return predmap, seq

def format_SPPIder():
    seq = ''
    count = 1
    predmap = {}
    previous = ''
    for line in open(spp):
        if not line.startswith('ATOM'): continue
        resnum = str(int(line[22:26]))
        resname = three2one[line[17:20]]
        score = round(float(line[61:66])/100, 2)
        if resnum != previous:
            predmap[count] = [resname, score]
            previous = resnum
            seq += resname
            count += 1
    return predmap, seq

def format_consurf():
    seq = ''
    predmap = {}
    main = False
    maxval = 3.666
    minval = -2.435
    outfile = open(con)
    for line in outfile:
        if line.startswith(' POS'): 
            outfile.readline()
            main=True
            continue
        if not main: continue
        if line.strip() == '': break
        n = int(line.split()[0])
        res = line.split()[1]
        score = (float(line.split()[3])-minval)/(maxval-minval)
        predmap[n] = [res, score]
        seq += res
    return predmap, seq

def align_predictions(pseq, rseq, pdic, rdic):
    alignments = pairwise2.align.globalxx(rseq, pseq)
    
    posr = posp = 1
    alignedr = alignments[0][0]
    alignedp = alignments[0][1]
    for r, p in zip(alignedr, alignedp):
        ##### if there is a match or mismatch
        if r != '-' and p != '-':
            if r != p: print ('Mismatch!', r, p)
            rdic[posr].append(pdic[posp][1])
            posr += 1
            posp += 1
        ##### if there is a gap in the prediction
        elif r != '-':
            rdic[posr].append(0.0)
            posr += 1
        ##### skip positions gapped in unbound sequence
        elif p != '-':
            posp += 1

def main(numb):
    print ('####################', pdb)
    ##### binding sites ground truth
    intmap, realseq, keylist = interface_extract()

    ##### ground truth controlled randomization
    for tpr, ppv in [[1,0.75], [1,0.5], [1,0.25], 
                     [0.75,1], [0.5,1], [0.25,1],
                     [0.5,0.5], [0.5,0.25], [0.25,0.5]]:
        rand = randomize_interfaces(tpr, ppv, intmap, keylist)
        for key in keylist:
            intmap[key].append(rand[key][1])

    ##### mapping predictions to unbound sequence

    pred, seq = format_BIPSPI(bip1)
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_BIPSPI(bip2)
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_ISPRED4()
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_ISPRED4(crf=True)
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_predus()
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_dynJET2()
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_dynJET2(trace=True)
    align_predictions(seq, realseq, pred, intmap)

    pred, seq = format_SPPIder()
    align_predictions(seq, realseq, pred, intmap)

    #pred, seq = format_consurf()
    #align_predictions(seq, realseq, pred, intmap)


    ##### write out
    outfile = open(ns.s+'/data/formatted_labels/'+pdb+'_u'+str(numb)+'.site', 'w')
    for pos in keylist:
        for n, val in enumerate(intmap[pos][1:]):
            if n != len(intmap[pos][1:])-1:
                outfile.write(str(round(float(val), 2))+'\t')
            else:
                outfile.write(str(round(float(val), 2))+'\n')
    outfile.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description = 'Converts PDB complex in ISPRED4 format')
    p.add_argument('-s', required= True, help='BindingSite_docking folder path')
    ns = p.parse_args()

    random.seed(42)
    p = PDBParser(QUIET=True)
    if not os.path.exists(ns.s+'/data/formatted_labels/'):
        os.mkdir(ns.s+'/data/formatted_labels/')

    for pdb in open(ns.s+'list_dimers'):
        pdb = pdb.rstrip()
        for n in [1,2]:
            pathu = '{}/data/processed_b4/{}_u{}.pdb'.format(ns.s, pdb, n)
            pathuc = '{}/data/processed_b4/{}_uc.pdb'.format(ns.s, pdb)
            pathb = '{}/data/processed_b4/{}_b{}.pdb'.format(ns.s, pdb, n)
            pathbc = '{}/data/processed_b4/{}_bc.pdb'.format(ns.s, pdb)
            bip1 = '{}/data/BIPSPI_predictions/mixed/{}_u{}.res'.format(ns.s, pdb, n)
            bip2 = '{}/data/BIPSPI_predictions/mixed_2/{}_u{}.res'.format(ns.s, pdb, n)
            isp = '{}/data/ISPRED_predictions/{}_u{}.ispred4.tsv'.format(ns.s, pdb, n)
            pre = '{}/data/predus_predictions/{}_u{}.int'.format(ns.s, pdb, n)
            dyn = '{}/data/dynJET2_predictions/{}_u{}_jet.res'.format(ns.s, pdb, n)
            spp = '{}/data/SPPIDER_predictions/{}_u{}.pdb'.format(ns.s, pdb, n)

            if n == 1:
                chain = 'A'
                oppchain = 'B'
            else:
                chain = 'B'
                oppchain = 'A'

            main(n)
            
    
    
