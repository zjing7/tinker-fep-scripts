#!/usr/bin/env python
'''
Find the protein atoms of which the center of mass is close to the ligand
> get_grp_rest.py sample.xyz sample.key

Read any traj format and "ligand keyword" from Tinker key file
Print Tinker group definition


1. Instead of dicrete opt., this is a dirty method that optimizes the weights 
and then choose the atoms with larger weights.
2. It choose only heavy atoms and assumes equal mass when calculating CofM
'''

import numpy as np
import mdtraj as md
import sys
import os
import time
import scipy.optimize
from scipy.spatial import distance_matrix

def error_exit(msg):
    print("ERROR:", msg)
    sys.exit(1)
def write_tinker_idx(idxs):
    '''Convert list of numbers to ranges in Tinker format
    '''
    idx_out = []
    rs = []
    for i0 in sorted(idxs):
        if len(rs) == 0 or rs[-1][1]+1 < i0:
            rs.append([i0, i0])
        else: 
            rs[-1][1] = i0
    for r in rs:
        if r[0] == r[1]:
            idx_out.append(r[0])
        elif r[0] == r[1] - 1:
            idx_out.append(r[0])
            idx_out.append(r[1])
        else:
            idx_out.append(-r[0])
            idx_out.append(r[1])
    return idx_out
        
def read_tinker_idx(args):
    ''' Read tinker indices

    Example:
    ['5'] -> {5}
    ['-5', '7', '10', '-12', '15'] -> {5, 6, 7, 10, 12, 13, 14, 15}

    args: list of strings of integers
    return: a set of indices

    '''
    _range = []
    idxs = []
    for a in args:
        n = int(a)
        if n < 0 and len(_range) == 0:
            _range.append(-n)
        if n > 0:
            if len(_range) == 1:
                idxs.extend(list(range(_range.pop(), n+1)))
            else:
                idxs.append(n)
    return set(idxs)

def read_ligidx(fkey):
    '''Read tinker key file and return indices of the ligand
    '''
    ligand_idx = set() # ligand indices
    with open(fkey, 'r') as fh:
        for line in fh:
            w = line.split()
            if line.lower().startswith('ligand') and len(w) >= 2:
                idx_str = line.replace(',', ' ').split()[1:]
                ligand_idx |= read_tinker_idx(idx_str)
    return ligand_idx


def target_disp(weights, coord, alpha=0.1):
    '''target function for center of mass displacement and number of atoms in the group
    '''
    assert len(weights) == coord.shape[0]
    wts = np.array(weights).reshape((-1, 1))
    wts = np.maximum(wts, 0)
    assert sum(wts) > 0
    wts *= 1.0/np.mean(wts)

    loss = np.sum(np.abs(np.sum(wts*coord, axis=0)))
    loss += alpha*np.sum(np.abs(np.abs(wts*coord)))

    for m in np.arange(2, 20):
        mask0 = (wts > m)
        loss += alpha*np.sum((np.abs((wts*mask0))))

    return loss

def find_grp_idx(fxyz, ligidx0, atomnames='CA', rcutoff=1.2, alpha=1.0):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    ligidx = np.array(sorted(list(set(ligidx0) - set(t.topology.select('name H')))))
    if len(ligidx) == 0:
        error_exit("No ligand heavy atoms found")

    protidx0 = t.topology.select('name %s'%(atomnames))
    protidx0 = np.array(sorted(set(protidx0) - set(ligidx)))

    distmat = distance_matrix(t.xyz[0, ligidx, :], t.xyz[0, protidx0, :])
    distpro = np.min(distmat, axis=0)
    protidx = protidx0[distpro <= rcutoff]
    if len(protidx) == 0:
        return

    #print(distmat.shape)
    #print(np.min(distmat))
    imindist = np.argmin(distmat)
    iminlig = ligidx[imindist // len(protidx0)]

    com_lig = np.mean(t.xyz[0, list(ligidx), :], axis=0)
    com_lig = t.xyz[0, [iminlig], :]


    wts0 = np.ones(len(protidx))
    xyz_prot = t.xyz[0, protidx, :]
    xyz1_prot = xyz_prot - com_lig.reshape((1, -1))
    res = scipy.optimize.minimize(target_disp, wts0, args=(xyz1_prot, alpha))
    wts1 = np.array(res.x)
    wts1 = np.maximum(0, wts1)
    wts1 *= 1.0/np.sum(wts1)
    wtm = np.mean(wts1[wts1 > 0])
    mask1 = wts1 > 0.4*wtm
    #print(wts1)
    #print(np.mean(wts1))
    #print(mask1)
    com_p1 = (np.mean(t.xyz[0, protidx[mask1], :], axis=0)).reshape((1, -1))

    xyz_lig = t.xyz[0, ligidx, :]
    dists = np.linalg.norm(xyz_lig - com_p1, axis=1)
    #print('DIST', dists)
    imin = np.argmin(dists)

    dcom = np.linalg.norm(com_lig - com_p1)
    dmin = np.linalg.norm(xyz_lig[imin] - com_p1)


    idx_tinker = write_tinker_idx([_+1 for _ in protidx[mask1]])
    #print('#', (' '.join(['%5d'%(_) for _ in protidx[mask1]])))
    outp = ''
    sgrp = ''
    for n in idx_tinker:
        sgrp += ' %5d'%n
        if len(sgrp) > 50 and n > 0:
            #print('group 1 %s'%sgrp)
            outp += ('group 1 %s\n'%sgrp)
            sgrp = ''
            
    outp += ('group 2 %s\n'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    outp += ("#r_0=%.3f"%(dmin*10))
    #print('group 2 %s'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    #print("#r_0=%.3f"%(dmin*10))
    return dmin*10, outp

def calc_dg_rest(k0, r1, r2, V0=1, RT=0.59):
    '''
    V0: standard conc. in mol/L
    '''
    v0_ang = V0 * (1e-3/1e-30) / 6.02e23

    assert k0 > 0

    RMAX = np.sqrt(100/k0)
    DR = np.sqrt(0.0001/k0)
    xs = np.arange(0, r2+RMAX, DR)
    boltz = np.ones_like(xs)
    mask1 = xs < r1
    mask2 = (xs > r1) * (xs < r2)
    mask3 = xs > r2

    boltz[mask1] = np.exp(-k0*(xs[mask1]-r1)**2.0/RT)
    boltz[mask3] = np.exp(-k0*(xs[mask3]-r2)**2.0/RT)

    vols = (xs+DR*0.5)**3.0 - np.maximum(0, xs-DR*0.5)**3.0
    vols *= 4.0/3*np.pi
    #v1 = np.sum(boltz * 4*np.pi*xs**2.0*DR)
    v1 = np.sum(boltz * vols)
    return RT*np.log(v0_ang/v1)

def main():
    fxyz = sys.argv[1]
    fkey = sys.argv[2]
    ligidx = [_-1 for _ in sorted(read_ligidx(fkey))]
    if len(ligidx) == 0:
        error_exit("ligand keyword not found")
    #find_grp_idx(fxyz, ligidx, 'CA N')
    res = None
    r_thr = 0.5
    # try a few parameters and choose the best solution
    for r0 in (0.5, 0.6, 0.7):
        for a in (1.0, 0.5, 0.2):
            rmin, outp = find_grp_idx(fxyz, ligidx, 'CA C', r0, alpha=a)
            #rmin, outp = find_grp_idx(fxyz, ligidx, 'CA CB', r0, alpha=a)
            if res is None or rmin < res[0]:
                res = rmin, outp
            if rmin < r_thr:
                break
        if rmin < r_thr:
            break
    print(res[1])
    r0_rest = 2.0
    if res[0] >= r0_rest:
        print("#Warning: r0 = %.3f > %.3f Angstrom"%(res[0], r0_rest))
        print("# This means the centers of mass of the two groups are too far away,")
        print("# which will lead to poor convergence. Please adjust group definitions.")
    r0_recom = min(r0_rest, res[0]+0.5)
    dgrest = calc_dg_rest(15.0, 0, r0_recom)
    print("#Recommended restraint\n#restrain-groups 1 2 15.0 0.0 %.3f\n#dGrest(kcal/mol) %.5f"%(r0_recom, dgrest))

main()

