import sys,os,json
import math
import tempfile
import argparse
import numpy as np
import random as rnd
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats.stats import pearsonr

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.utility import vector1_int
from pyrosetta.rosetta.protocols.docking import DockingLowRes, DockMinMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose, renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def ISPRED_site(ispred, sequence, sep, thr=0.9):
    score_list = []
    for p, line in enumerate(open(ispred)):
        score = float(line.split()[-1])
        if score < thr or sequence[p] == 'G': continue
        score_list.append([p+sep+1, score]) 
    score_list = np.array(score_list, dtype=np.float32)
    score_list = score_list[score_list[:,1].argsort()[::-1]]
    score_list = np.array(score_list[:,0].flatten(), dtype=np.int)
    return score_list

def get_main_atom(res):
    atoms = [atom.get_id() for atom in res]
    if 'CB' in atoms: return 'CB'
    elif 'CA' in atoms: return 'CA'
    else: return None

def find_center(coord_list):
    acc = [0,0,0]
    for x, y, z in coord_list:
        acc[0]+=x
        acc[1]+=y
        acc[2]+=z
    center = [coord/len(coord_list) for coord in acc]
    return center

def get_closest_clusters(siteatoms, clusterdic, largest, thr):
    centers = {}
    for key in clusterdic:
        coords = [siteatoms[idx].get_coord() for idx in clusterdic[key]]
        centers[key] = find_center(coords)

    good_clust = []
    main_center = np.array(centers[largest])
    for key in centers:
        if key == largest: continue
        neigh = np.array(centers[key])
        dist = np.linalg.norm(main_center-neigh)
        if dist < thr: good_clust.append(key)
    return good_clust

def format_allvsall_rst(siteA, siteB):
    array = []
    for resA in siteA:
        resA = int(resA)
        atomA = get_main_atom(str1['A'][resA])
        for resB in siteB:
            resB = int(resB)
            atomB = get_main_atom(str2['B'][resB])
            if atomA is None or atomB is None: continue
            array += ['AtomPair {} {} {} {} FLAT_HARMONIC 8 10 4'\
                      .format(atomA, resA, atomB, resB)]
    return array

def format_cluster_rst(siteA, siteB):
    ##### dictionary to map indexes to CB/CA atoms
    siteatomsA = {}
    siteatomsB = {}
    for res in Selection.unfold_entities(str1, 'R'):
        if res.get_id()[1] in siteA and get_main_atom(res) != None: 
            siteatomsA[res.get_id()[1]] = get_main_atom(res)
    for res in Selection.unfold_entities(str2, 'R'):
        if res.get_id()[1] in siteB and get_main_atom(res) != None:
            siteatomsB[res.get_id()[1]] = get_main_atom(res)
    idxA = list(siteatomsA.keys())
    idxB = list(siteatomsB.keys())
    coordsA = np.array([siteatomsA[idx].get_coord() for idx in idxA])
    coordsB = np.array([siteatomsB[idx].get_coord() for idx in idxB])

    ##### clustering binding site predictions
    linkA = linkage(coordsA)
    linkB = linkage(coordsB)
    clustA = fcluster(linkA, 8, criterion='distance')
    clustB = fcluster(linkB, 8, criterion='distance')

    clustdicA = {}
    clustdicB = {}
    for idx, clust in enumerate(clustA):
        clustdicA[clust] = clustdicA.get(clust, []) + [idxA[idx]]
    for idx, clust in enumerate(clustB):
        clustdicB[clust] = clustdicB.get(clust, []) + [idxB[idx]]

    ##### find largest clusters
    maxlen = 0
    largestA = ''
    for idx in clustdicA: 
        if len(clustdicA[idx]) > maxlen: 
            largestA = idx
            maxlen = len(clustdicA[idx])
    maxlen = 0
    largestB = ''
    for idx in clustdicB: 
        if len(clustdicB[idx]) > maxlen:
            largestB = idx
            maxlen = len(clustdicB[idx])
    
    ##### max size of the smallest cluster
    dA = [siteatomsA[idx1]-siteatomsA[idx2] for idx1 in clustdicA[largestA] for idx2 in clustdicA[largestA]]
    dB = [siteatomsB[idx1]-siteatomsB[idx2] for idx1 in clustdicB[largestB] for idx2 in clustdicB[largestB]]
    dmax = min(max(dA), max(dB))

    neighbours = get_closest_clusters(siteatomsA, clustdicA, largestA, dmax)

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(coordsA[:,0], coordsA[:,1], coordsA[:,2], c=clustA, cmap='prism')  # plot points with cluster dependent colors
    #ax.scatter(coordsB[:,0], coordsB[:,1], coordsB[:,2], c=clustB, cmap='terrain')
    #plt.show()
    #sys.exit()

    siteA = clustdicA[largestA]
    siteB = clustdicB[largestB]
    arrays = [['AtomPair CB {} CB {} FLAT_HARMONIC 8 1 4'.format(a, b) for a in siteA for b in siteB]]
    for clust in neighbours:
        siteA = clustdicA[clust]
        subarray = ['AtomPair CB {} CB {} FLAT_HARMONIC 8 1 4'.format(a, b) for a in siteA for b in siteB]
        arrays.append(subarray)
    return arrays

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
    s1 = ISPRED_site(file1, seq1, 0)
    s2 = ISPRED_site(file2, seq2, len1)
    array = format_allvsall_rst(s1, s2)
    minimize(pose, relax, SF, array, 'top')

def minimize(pose, mover, sf, array, tag):
    print ('True pose energy:'+str(sf(pose)))
    rotation = 60
    translation = 10
    dock_pert = RigidBodyPerturbMover(1, rotation, translation)
    for n in range(10):
        for _ in range(1000): dock_pert.apply(pose)
        print ('Pose energy before dock:'+str(sf(pose)))
        apply_rst(pose, array, tmpdir.name+'/minimize.cst')
        mover.apply(pose)
        pose.remove_constraints()
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
    seq1 = p1.sequence()
    seq2 = p2.sequence()
    len1 = len(p1.sequence())
    len2 = len(p2.sequence())
    append_pose_to_pose(p1,p2)

    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    to_centroid = SwitchResidueTypeSetMover('centroid')
    #to_centroid.apply(p1)

    #renumber_pdbinfo_based_on_conf_chains(p1)
    seq = p1.sequence()
    lenseq = len(p1.sequence())

    custom_docking(p1, ns.i1, ns.i2)



