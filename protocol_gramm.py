#!/usr/bin/env python3
import math
import sys
import argparse
import numpy  as np
import tensorflow as tf
import multiprocessing as mp 
from functools import partial
from Bio.PDB import *


weights = {'C':12.011, 'N':14.007, 'O':15.999, 'S':32.06}

allowed_atoms = ['N','CA','C','O','OXT','CB',
                 'CZ','CZ1','CZ2','CE','CE1','CE2',
                 'CD','CD1','CD2','CG','CG1','CG2','CG3',
                 'CH2','OE','OE1','OE2','OD','OD1','OD2',
                 'OH','OG','OG1','OG2','NZ','NE','NE1','NE2',
                 'NH','NH1','NH2','ND','ND1','ND2',
                 'SG','SD','SE']

def get_sep(struc):

    three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
                 'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
                 'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                 'MET':'M','PHE':'F','PRO':'P','SER':'S',
                 'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
                 'MSE':'M'}
    seq = ''
    prv = ''
    chain = ''
    for line in open(struc):
        if line.startswith('ATOM'):
            if chain == '': chain = line[21]
            elif chain != line[21]: break

            if line[22:27].strip() != prv:
                seq += three2one[line[17:20]]
                prv = line[22:27].strip()
        if line.startswith('TER'): break

    return seq

def ligand_center(atom_list, geometric=False):
    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        if atom.get_id() in allowed_atoms: 
            masses.append(weights[atom.get_id()[0]])

        for i, coord in enumerate(atom.get_coord()): 
            positions[i].append(coord)

    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:       
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]

def get_rototranslation(data):
    count = 0
    result_batch = []
    with tf.Session(graph=get_rt, config=config) as sess:
        for pose, sc, m1, m2 in data:
            r_mat, t_vec = sess.run([r, t], feed_dict = {gr:m1, gt:m2, lcm:cm})
            result_batch.append([pose, sc, r_mat, t_vec])
            count += 1
            if count%1000 == 0: print ('Processed {} rt matrixes!'.format(count))
    return result_batch

def rototranslate_coord(data, c_only=True):
    count = 0
    result_batch = []
    with tf.Session(graph=rt_comp, config=config) as sess:
        for rt in data:
            rt_coord, p_score = sess.run([rtcoord, score], feed_dict = {shape:rt[1], r_mat:rt[2], t_vec:rt[3]})
            result_batch.append([rt[0], rt_coord, p_score])
            count += 1
            if count%1000 == 0: print ('Roto-translated {} structures!'.format(count))
    return result_batch

def get_main_coord(res):
    atoms = [atom.get_id() for atom in res]
    if 'CB' in atoms: return 'CB'
    elif 'CA' in atoms: return 'CA'
    else: return None

def split_jobs(job_list, cores):
    batch_list = []
    bs = math.floor(len(job_list)/cores)
    for n in range(0,cores-1): batch_list.append(job_list[n*bs:(n+1)*bs])
    batch_list.append(job_list[bs*(cores-1):])
    return batch_list

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

        pool = mp.Pool(processes=cores)
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
        elif l1[p1][2] >= l2[p2][2]:
            l.append(l1[p1])
            p1 += 1
        else:
            l.append(l2[p2])
            p2 += 1
    return l


if __name__ == "__main__":
    p = argparse.ArgumentParser(description = '- Plot PPV stats for a Cmap or a list of Cmaps')
    p.add_argument('-g', required= True, help='gramm output')
    p.add_argument('-i1', required= True, help='trRosetta predictions')
    p.add_argument('-i2', required= True, help='trRosetta predictions')
    p.add_argument('-s1', required= True, help='structure file 1')
    p.add_argument('-s2', required= True, help='structure file 2')
    p.add_argument('-o', required= False, default=None, help='output file path and prefix to write models (no format)')
    ns = p.parse_args()

    ##### parse contacts #####
    contactidxs = []
    scores1 = [float(line.rstrip().split('\t')[-1]) for line in open(ns.i1)]
    scores2 = [float(line.rstrip().split('\t')[-1]) for line in open(ns.i2)]
    scores1 = np.array(scores1, dtype=np.float32)
    scores2 = np.array(scores2, dtype=np.float32)
    for p in range(scores1.shape[0]): 
        if scores1[p] == 0.0: scores1[p] = 0.01
    for p in range(scores2.shape[0]): 
        if scores2[p] == 0.0: scores2[p] = 0.01
    scores1 = np.expand_dims(scores1, axis=1)
    scores2 = np.expand_dims(scores2, axis=0)
    cmap = scores1*scores2
    print (cmap, cmap.shape)

    sep = len(get_sep(ns.s1))
    contactids = [y+sep+1 for y in range(0, cmap.shape[1]) if np.any(cmap[:,y] >= 0.9)]

    ##### parse structures #####
    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', ns.s1)
    str2 = p.get_structure('', ns.s2)
    str3 = p.get_structure('', ns.s2)
    for c in str2[0]: lchainid = c.get_id()
    rec_res = Selection.unfold_entities(str1, 'R')
    lig_res = Selection.unfold_entities(str2, 'R')

    ##### ligand real interface CB/CA coordinates #####
    lcoordinates = []
    if (len(contactids)==0):
        print("No contacts found ")
        sys.exit(0)
    for idx in contactids: 
        atom = get_main_coord(str2[0][lchainid][idx])
        if atom is None: continue
        lcoordinates.append([c for c in str2[0][lchainid][idx][atom].get_coord()])
    lcoordinates = np.array(lcoordinates, dtype=np.float32)
    if (len(lcoordinates)==0):
        print("No ligand coordinates found ")
        sys.exit(0)
    
    ##### ligand CB/CA coordinates #####
    full_lig = []
    full_lig_id = []
    #print (lig_res)
    for res in lig_res:
        resid = res.get_id()
        for atom in res:
            atomid = atom.get_id()
            full_lig.append([c for c in atom.get_coord()])
            full_lig_id.append([resid, atomid])
    full_lig = np.array(full_lig, dtype=np.float32)

    ##### receptor CB/CA coordinates #####
    rcoordinates = []
    for res in rec_res: 
        atom = get_main_coord(res)
        if atom is None: continue
        rcoordinates.append([c for c in res[atom].get_coord()])
    rcoordinates = np.array(rcoordinates, dtype=np.float32)
    if (len(rcoordinates)==0):
        print("No receptor coordinates found")
        sys.exit(0)
    
    ##### get contact probabilities #####
    contactids = np.array(contactids, dtype=np.int)
    lrprobs = cmap[:,contactids-sep-1]
    print (lrprobs)

    ##### calculate lcm #####
    lig_atoms = Selection.unfold_entities(str2, 'A')
    cm = np.array(ligand_center(lig_atoms), dtype=np.float32)

    ##### graph to elaborate GRAMM output #######################################
    with tf.Graph().as_default() as get_rt:                                     #
        with tf.name_scope('input0'):                                           #
            pi = tf.constant(math.pi/180.0)                                     #
            gr = tf.placeholder(dtype=tf.float32, shape=(3))                    #
            gt = tf.placeholder(dtype=tf.float32, shape=(3))                    #
            lcm = tf.placeholder(dtype=tf.float32, shape=(3))                   #
                                                                                #
        gr1 = pi*gr                                                             #
        sina = tf.math.sin(gr1[0])                                              #
        sinb = tf.math.sin(gr1[1])                                              #
        sing = tf.math.sin(gr1[2])                                              #
        cosa = tf.math.cos(gr1[0])                                              #
        cosb = tf.math.cos(gr1[1])                                              #
        cosg = tf.math.cos(gr1[2])                                              #
                                                                                #
        r = tf.stack([tf.stack([cosg*cosa,                                      #
                                cosg*sina,                                      #
                                -sing]),                                        #
                      tf.stack([sinb*sing*cosa-cosb*sina,                       #
                                sinb*sing*sina+cosb*cosa,                       #
                                cosg*sinb]),                                    #
                      tf.stack([cosb*sing*cosa+sinb*sina,                       #
                                cosb*sing*sina-sinb*cosa,                       #
                                cosg*cosb])])                                   #
                                                                                #
        r_lcm = tf.linalg.matvec(r, lcm)                                        #
        t = gt+(lcm-r_lcm)                                                      #
        t = tf.expand_dims(t, axis=-1)                                          #
    #############################################################################

    ##### graph to roto-translate and score atom coordinates ####################
    with tf.Graph().as_default() as rt_comp:                                    #
        with tf.name_scope('input1'):                                           #
            pr = tf.constant(lrprobs)                                           #
            xyz = tf.constant(lcoordinates)                                     #
            rec = tf.constant(rcoordinates)                                     #
            shape = tf.placeholder(dtype=tf.float32, shape=())                  #
            t_vec = tf.placeholder(dtype=tf.float32, shape=(3, 1))              #
            r_mat = tf.placeholder(dtype=tf.float32, shape=(3, 3))              #
                                                                                #
        xyz = tf.transpose(xyz, perm=[1,0])                                     #
        rtcoord = tf.math.add(t_vec, tf.linalg.matmul(r_mat, xyz))              #
        rtcoord = tf.transpose(rtcoord, perm=[1,0])                             #
                                                                                #
        ##### scoring #####                                                     #
        d_thr = 6                                                               #
        pr = 1+tf.math.log(pr)                                                  #
        A = tf.expand_dims(rec, axis=1)                                         #
        B = tf.expand_dims(rtcoord, axis=0)                                     #
        distances = tf.math.sqrt(tf.reduce_sum((A-B)**2, axis=-1))              #
        zeros = tf.zeros(distances.shape, dtype=tf.float32)                     #
        scores = tf.where(tf.math.less(distances, d_thr), pr, zeros)            #
        score = tf.math.add(tf.math.reduce_sum(scores), shape)                  #
    #############################################################################

    ##### graph for full roto-translation #######################################
    with tf.Graph().as_default() as full_rt:                                    #
        with tf.name_scope('input2'):                                           #
            full_xyz = tf.constant(full_lig)                                    #
            full_t = tf.placeholder(dtype=tf.float32, shape=(3, 1))             #
            full_r = tf.placeholder(dtype=tf.float32, shape=(3, 3))             #
                                                                                #
        full_xyz = tf.transpose(full_xyz, perm=[1,0])                           #
        full_rtcoord = tf.math.add(full_t, tf.linalg.matmul(full_r, full_xyz))  #
        full_rtcoord = tf.transpose(full_rtcoord, perm=[1,0])                   #
    #############################################################################

    cores = 5 #mp.cpu_count()-1
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True

    ##### compute rototranslations #####
    mat_jobs = []
    poses = False
    for line in open(ns.g):
        if '[match]' in line: 
            poses = True
            continue
        if not poses: continue
        line = [float(el.strip()) for el in line.split() if el.strip() != '']
        grammr = np.array([line[2], line[3], line[4]], dtype=np.float32)
        grammt = np.array([line[5], line[6], line[7]], dtype=np.float32)
        mat_jobs.append([int(line[0]), int(line[1]), grammr, grammt])

    pool = mp.Pool(processes = cores)
    job_list = split_jobs(mat_jobs, cores)
    results = pool.map(get_rototranslation, job_list)
    pool.close()
    pool.join()

    ##### rototranslate and score poses #####
    rtlist = []
    rtdict = {}
    for batch in results: 
        rtlist.extend(batch) 
        for result in batch: rtdict[result[0]] = result[2:]

    pool = mp.Pool(processes = cores)
    job_list = split_jobs(rtlist, cores)
    results = pool.map(rototranslate_coord, job_list)
    pool.close()
    pool.join()

    ##### sort out results #####
    scorelist = []
    for batch in results: scorelist.extend(batch)
    sortedlist = mergesort_pred(scorelist)

    ##### print out #####
    io = PDBIO()
    with tf.Session(graph=full_rt, config=config) as sess:
        for n in range(1000):
            pose = sortedlist[n][0]
            score = sortedlist[n][2]
            print ('# {} - Pose {} - Score {}'.format(n+1, pose, score))
            if ns.o is None: continue
            if n+1 <= 10 or (n+1%10 == 0 and n+1 <= 100) or (n+1%100 == 0):
                rt_coord = sess.run(full_rtcoord, feed_dict = {full_r:rtdict[pose][0], full_t:rtdict[pose][1]})
                for coord, ids in zip(rt_coord, full_lig_id): str3[0][lchainid][ids[0]][ids[1]].set_coord(coord)
                io.set_structure(str3)
                io.save('{}_{}.pdb'.format(ns.o, n+1))

