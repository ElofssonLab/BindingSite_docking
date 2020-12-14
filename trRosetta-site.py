import sys,os,json
import tempfile
import argparse
import numpy as np
import random as rnd
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
from sklearn.cluster import AgglomerativeClustering
from scipy.stats.stats import pearsonr

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.utility import vector1_int
from pyrosetta.rosetta.protocols.docking import DockingLowRes, DockMinMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose, renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import *

def ISPRED_thr_site(ispred, sep, thr=0.1):
    count = 0
    site = []
    for line in open(ispred):
        count += 1
        res = count + sep
        score = float(line.split()[-1])
        if seq[res-1] == 'G': continue
        if score > thr: site.append(res)
    return site

def ISPRED_top_site(ispred, sep, top=3):
    count = 0
    plist = []
    for line in open(ispred):
        count += 1
        res = count + sep
        score = float(line.split()[-1])
        if seq[res-1] == 'G': continue
        else: 
            for p, r in enumerate(plist):
                if score > r[1]: 
                    plist = plist[:p]+[[res,score]]+plist[p:]
                    break
            if [res,score] not in plist: plist.append([res,score])           
    site = [el[0] for el in plist]
    return site

def get_main_atom(res):
    atoms = [atom.get_id() for atom in res]
    if 'CB' in atoms: return res['CB']
    elif 'CA' in atoms: return res['CA']
    else: return None

def format_ISPRED_rst(siteA, siteB, strA, strB):
    siteatomsA = {}
    siteatomsB = {}
    for res in Selection.unfold_entities(strA, 'R'):
        if res.get_id[1] in siteA and get_main_coord(res) != None: 
            siteatomsA[res.get_id[1]] = get_main_atom(res)
    for res in Selection.unfold_entities(strB, 'R'):
        if res.get_id[1] in siteB and get_main_coord(res) != None:
            siteatomsB[res.get_id[1]] = get_main_atom(res)
    idxA = siteatomsA.keys()
    idxB = siteatomsB.keys()
    coordA = [siteatomsA[idx].get_coord() for idx in idxA]
    coordB = [siteatomsB[idx].get_coord() for idx in idxB]
    clust = AgglomerativeClustering(distance_threshold=12).fit(coordA)
    clusteredA = clust.labels_
    clust = AgglomerativeClustering(distance_threshold=12).fit(coordB)
    clusteredB = clust.labels_
    print (clusteredA, clusteredB)

    array = ['AtomPair CB {} CB {} FLAT_HARMONIC 8 1 4'.format(a, b) for a in siteA for b in siteB]
    return array

def apply_rst(pose, array, constraints_file):
    with open(constraints_file,'w') as f:
        for line in array: f.write(line+'\n')
    constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(constraints_file)
    constraints.add_constraints(True)
    constraints.apply(pose)
    os.remove(constraints_file)

def custom_docking(pose, file1, file2):
    
    ##### Scoring function #####
    score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    constraint = score_manager.score_type_from_name('atom_pair_constraint')

    sf1 = create_score_function('docking')
    sf1.set_weight(constraint, 1)
    SF = sf1

    ##### MoveMap #####
    mmap = MoveMap()
    mmap.set_bb(False)
    mmap.set_chi(True)
    mmap.set_jump(True)

    ##### Mover #####
    relax = rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(SF)
    relax.dualspace(True)
    relax.set_movemap(mmap)

    ##### Top Confidence constraint #####
    s1 = ISPRED_top_site(file1, 0)
    s2 = ISPRED_top_site(file2, len1)
    array = format_ISPRED_rst(s1, s2)
    print ('Extracted {} constraints'.format(len(array)))
    minimize(pose, relax, SF, array, 'top')

    ##### Confidence Thr. constraint #####
    #s1, ns1 = ISPRED_thr_site(file1, 0)
    #s2, ns2 = ISPRED_thr_site(file2, len1)
    #array = format_ISPRED_rst(s1, s2, ns1, ns2)
    #print ('Extracted {} constraints'.format(len(array)))
    #minimize(pose, relax, SF, array, 'thr')

def minimize(pose, mover, sf, array, tag):

    pose.remove_constraints()
    apply_rst(pose, array, tmpdir.name+'/minimize.cst')
    print ('True pose energy:'+str(sf(pose)))

    rotation = 30
    translation = 10
    dock_pert = RigidBodyPerturbMover(1, rotation, translation)

    for n in range(10):
        pose.remove_constraints()
        for _ in range(1000): dock_pert.apply(pose)
        apply_rst(pose, array, tmpdir.name+'/minimize.cst')
        print ('Pose energy before dock:'+str(sf(pose)))
        mover.apply(pose)
        print ('Pose energy after dock:'+str(sf(pose)))
        pose.dump_pdb(ns.out+'_'+tag+str(n+1)+'.pdb')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = "Class to implement trRosetta")
    parser.add_argument("str1", type=str, help="input structure 1")
    parser.add_argument("str2", type=str, help="input structure 2")
    parser.add_argument("out", type=str, help="output model path/prefix")
    parser.add_argument("-i1", type=str, help="ISPRED file for prot 1")
    parser.add_argument("-i2", type=str, help="ISPRED file for prot 2")
    parser.add_argument("-s", type=str, help="script directory")
    parser.add_argument("-d", type=str, help="data directory")
    ns = parser.parse_args()

    init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200\
          -rebuild_disulf false -detect_disulf false -out:level 100')
    scriptdir = ns.s.rstrip('/')+'/'
    datadir = ns.d.rstrip('/')+'/'
    tmpdir = tempfile.TemporaryDirectory(prefix=datadir)

    ##### biopython objects #####
    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', ns.str1)[0]
    str2 = p.get_structure('', ns.str2)[0]

    ##### Pose initialization #####
    p1 = pose_from_pdb(ns.str1)
    p2 = pose_from_pdb(ns.str2)
    len1 = len(p1.sequence())
    len2 = len(p2.sequence())
    append_pose_to_pose(p1,p2)

    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    to_centroid = SwitchResidueTypeSetMover('centroid')
    #to_centroid.apply(p1)

    renumber_pdbinfo_based_on_conf_chains(p1)
    seq = p1.sequence()
    lenseq = len(p1.sequence())

    custom_docking(p1, ns.i1, ns.i2)



