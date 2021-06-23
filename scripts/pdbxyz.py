#!/usr/bin/env python2
'''
Convert PDB to XYZ, or generate template file from PDB and XYZ
'''

from __future__ import print_function
import re
import numpy as np
import sys, os
import argparse

class tinkerxyz():
    def __init__(self):
        #tinkermd.__init__(self)
        self.top_natom = 0
        self.top_conn = {}
        self.top_type = {}
        self.top_name = {}

        self.nframes = 0
        self.frames = []
        self.usepbc = False

        self.fxyzout = None
        self.FORMAT = ' %5d %-4s %11.6f %11.6f %11.6f %5d'
        self.FORMAT1 = ' %5d'

    def read_xyz(self, xyzfile):
        '''
        Read topology information
        '''
        with open(xyzfile, 'r') as fin:
            iline = -1
            istart = 1
            for line in fin:
                iline = iline + 1
                w = line.split()
                if iline == 0:
                    self.top_natom = int(w[0])
                    continue
                elif iline == 1:
                    if len(w) == 6 and tuple(map(float, w[3:6])) == (90,90,90):
                        istart = 2
                        continue
                    else:
                        istart = 1
                if iline >= istart + self.top_natom:
                    break
                id = int(w[0])
                self.top_name[id] = w[1]
                self.top_type[id] = int(w[5])
                self.top_conn[id] = tuple(map(int, w[6:]))
            #self.get_residue()
            #self.top_protein = self.select(type=self.prm[0])
            return 0
        print("WARNING: Cannot open file %s"%xyzfile, file=sys.__stderr__)
        return 1
    def write_xyz(self, frame, fout):
        natom = len(self.top_name)
        outp = ' %5d\n'%(natom)
        if 'cell' in frame:
            cell = frame['cell']
            outp = outp + ' '+(' %11.6f'*6)%(cell[0], cell[1], cell[2], 90, 90, 90) + '\n'
        for i in range(1, natom+1):
            coord = frame[i]
            outp += self.FORMAT%(i, self.top_name[i], coord[0], coord[1], coord[2], self.top_type[i])
            conn = []
            if i in self.top_conn:
                conn = self.top_conn[i]
            outp += ''.join([self.FORMAT1%(K) for K in conn])
            outp += '\n'
        fout.write(outp)

    def append_frame(self, frame):
        self.frames.append(frame)

    def read_arc(self, inpfile, f):
        '''
        Read tinker traj file and apply function `f' to the current frame
        The frame is a dictionary that maps atom ids or `cell' to coordinates/cell parameters

        For example, `f' can be `append_frame'
        '''
        natom = self.top_natom
        if inpfile.endswith('.gz'):
            fin = gzip.open(inpfile, 'rb')
        else:
            fin = open(inpfile, 'r')
        if True:
        #with open(inpfile, 'r') as fin:
            iline = -1
            istart = 1
            iframe = -1
            currframe = {}
            for line in fin:
                iline = iline + 1
                w = line.split()
                if iline == 1:
                    if len(w) == 6 and tuple(map(float, w[3:6])) == (90,90,90):
                        istart = 2
                        currframe['cell'] = np.array(list(map(float, w[0:3])))
                        self.usepbc = True
                        continue
                    else:
                        istart = 1
                if iline >= 1:
                    #frame number
                    iatom = (iline - istart) % (natom + istart)
                    isec = (iline - istart)/(natom+istart)
                    if istart == 2 and iatom == natom+istart-1 and len(w) == 6:
                        currframe['cell'] = np.array(tuple(map(float, w[0:3])))
                    if iatom < natom and iline>=1:
                        id = int(w[0])
                        coord = np.array(tuple(map(float, w[2:5])))
                        currframe[id] = coord
                    if iatom == natom - 1:
                        f(currframe)
                        self.nframes = self.nframes + 1
                        currframe = {}
            return 0
        return 1


def splitline(line):
    '''
    Convert a line to fields separated by white spaces.
    Words inside "" are treated as one field.
    '''
    #result = list(re.finditer('"(.+)"|(\S+)\s', line))
    #result = list(re.finditer('(")?((?(1).+|\S+))(?(1)"|\s)', line))
    result = list(re.finditer('(")?((?(1)[^"]+|\S+))(?(1)"|\s)', line))
    words = [K.group(2) for K in result]
    return words

def convertpdb(pdbfile, res_def, verbose=True, bioclass=''):
    '''
    bioclass: pro or nuc
    '''
    with open(pdbfile) as fin:
        lines = fin.readlines()

        blocks = []
        last_res = [-999]
        curr_blk = []
        # convert PDB to residue blocks
        for line in lines:
            if len(line) >= 54 and line[:6] in ('ATOM  ', 'HETATM'):
                atid = int(line[6:11])
                atname = (line[12:16]).strip()
                resname = (line[17:20]).strip()
                if len(line) > 78:
                    elename = line[76:78]
                else:
                    elename = ''
                chain = (line[21])
                resid = int(line[22:26])
                r_x = float(line[30:38])
                r_y = float(line[38:46])
                r_z = float(line[46:54])

                # residue identifier
                res = [chain, resid, resname]

                if res != last_res:
                    # new residue
                    if len(curr_blk) > 0:
                        blocks.append(curr_blk)
                    curr_blk = [res]
                if res[0] != last_res[0]:
                    # first residue or after TER
                    blocks.append([' '])

                # name, xyz, atype
                curr_blk.append([atname, (r_x, r_y, r_z), 0])
                last_res = res

            elif line[:3] in ('TER', 'END'):
                last_res = [-999]
        if len(curr_blk) > 0:
            blocks.append(curr_blk)
    nblk = len(blocks)
    # assign type for each residue block
    for iblk, blk in enumerate(blocks):
        natom = len(blk) - 1
        if len(blk) < 2:
            continue
        # 0, N-term or 5'; 1, C_term or 3'; 2, non-term
        # Currently does handle cyclic
        resname = blk[0][2]
        resID = ' '.join([str(K) for K in blk[0]])
        if iblk == 0 or len(blocks[iblk-1]) == 1:
            flag_term = 0
            altnames = ('N'+resname, resname+'5', resname)
            if bioclass == 'nuc' or resname in 'DA DT DG DC A U G C'.split():
                blk[0][2] = resname+'5'
            else:
                blk[0][2] = 'N'+resname
        elif iblk == nblk-1 or len(blocks[iblk+1]) == 1:
            flag_term = 1
            altnames = ('C'+resname, resname+'3', resname)
            if bioclass == 'nuc' or resname in 'DA DT DG DC A U G C'.split():
                blk[0][2] = resname+'3'
            else:
                blk[0][2] = 'C'+resname
        else:
            flag_term = 2
            altnames = (resname,)
        for resff in altnames:
            resff_ID = (resff, natom)
            if resff_ID in res_def:
                curr_def = res_def[resff_ID]
                unique_name = curr_def[3]
                pdbnames = list(sorted([K[0] for K in blk[1:]]))
                ffnames = list(sorted(curr_def[0]))
                if unique_name and pdbnames == ffnames:
                    for curr_line in blk[1:]:
                        atype = curr_def[2][curr_line[0]]
                        curr_line[2] = atype
                elif len(pdbnames) == len(ffnames):
                    if verbose:
                        print(' Assign atom type based on the order of appearance...\n'\
                            ' Please make sure PDB has the same atom order as Template.')
                    for iat, curr_line in enumerate(blk[1:]):
                        atype = curr_def[1][iat]
                        curr_line[2] = atype
                    for iat in range(len(pdbnames)):
                        if pdbnames[iat] != ffnames[iat] and verbose:
                            print('WARNING: Conflicting atom name in %s: %5s %5s' \
                                  %(resID, pdbnames[iat], ffnames[iat]), file=sys.__stderr__)
                elif verbose:
                    print('WARNING: Numbers of atoms in PDB and FF do not match (%d, %d) for %s'\
                          %(len(pdbnames) , len(ffnames), resID), resff, file=sys.__stderr__)
                break
    return blocks


def read_template(inpfile):
    '''
    Read residue defination.
    Residue is identified by residue name and number of atoms in the residue.
    (str RESNAME, int NATOM):
            [list of NAME, 
             list of TYPE,
             dictionary that maps NAME to TYPE,
             boolean of whether NAME is unique]
    '''
    res_def = {}
    if inpfile == None:
        return res_def
    with open(inpfile) as fin:
        lines = fin.readlines()
        lastresidue = ''
        curr_blk = [[], [], {}, True]
        for line in lines:
            w = splitline(line.split('#')[0])
            nw = len(w)
            if nw >= 5 and w[0] == 'atom':
                atype = int(w[1])
                aname = w[3]
                residue = w[4]
                if (residue != lastresidue):
                    if len(curr_blk[0]) > 0:
                        res_ID = (lastresidue, len(curr_blk[0]))
                        res_def[res_ID] = curr_blk
                    curr_blk = [[], [], {}, True]
                curr_blk[0].append(aname)
                curr_blk[1].append(atype)

                lastresidue = residue
        if len(curr_blk[0]) > 0:
            res_ID = (residue, len(curr_blk[0]))
            res_def[res_ID] = curr_blk
    for residue in res_def:
        pair = zip(res_def[residue][0], res_def[residue][1])
        d_pair = dict(pair)
        if len(d_pair) < len(set(pair)):
            # Atom name not unique
            res_def[residue][3] = False
        else:
            res_def[residue][3] = True
            res_def[residue][2].update(d_pair)
    return res_def
                
def train_tplt(pdbfile, xyzfile, biocls='pro'):
    '''
    # Library for converting PDB to TINKER XYZ.

    # Residues are identified by resname and number of atoms. The TINKER atom 
    #   types are assigned by looking up the atom name in the library file. 
    #   If atom names are not unique or do not match, atom types will be
    #   assigned based on the order of appearance in the library.

    # Continuous atom entries with the same resname are considered one residue.
    # Two residues can be defined with same name but different numbers of atoms.
    #   In this case, these two residues must be separated by atom entries with
    #   a different resname. For example:

    #     atom    1    1 N     NH2
    #     atom    2    2 H1    NH2
    #     atom    3    3 H2    NH2

    #     atom    0    0 -     --------

    #     atom    4    4 N     NH2
    #     atom    5    5 H1    NH2
    #     atom    6    6 H2    NH2
    #     atom    7    7 H3    NH2


# atom; `atom type'; `atom class'; `atom name'; `residue name'; ...
    '''
    residues = convertpdb(pdbfile, {}, verbose=False, bioclass=biocls)
    t = tinkerxyz()
    t.read_xyz(xyzfile)
    t.read_arc(xyzfile, t.append_frame)

    frame = t.frames[0]

    pt_sym =  'C H O N P S'.split()
    pt_num =  [6,1,8,7,15,16,]
    pt_mass = [12,1,16,14,31,32,]
    crd2ff = {}
    for i in t.top_type:
        atype = t.top_type[i]
        crd = tuple(frame[i])
        aname = t.top_name[i]
        if aname[0] in pt_sym:
            ele_idx = pt_sym.index(aname[0])
            atnum = pt_num[ele_idx]
            atmass = pt_mass[ele_idx]
        else:
            atnum = 6
            atmass = 12
        crd2ff[crd] = [atype, atype, atnum, atmass, 0]

    reslib = {}
    r_thr = 0.5
    r_thr2 = 0.01
    used = set([])
    for res in residues:
        if len(res) <= 1:
            continue
        resname = res[0][2]
        missing = False
        natom = len(res) - 1
        curr_out = '\n# RESIDUE %s %d\n'%(resname, natom)
        for atom in res[1:]:
            crd = tuple(atom[1])
            aname = atom[0]
            if crd in crd2ff and crd not in used:
                ff = crd2ff[crd]
                used.add(crd)
                curr_out += 'atom %4d %4d %-5s %-20s %4d %7.4f %2d\n'%(ff[0], ff[1], aname, resname, ff[2], ff[3], ff[4])
            else:
                drmin = 1e5
                rmin = (0,0,0)
                for crd0 in crd2ff:
                    dr = np.linalg.norm(np.array(crd0)-np.array(crd))
                    if dr < drmin:
                        drmin = dr
                        rmin = crd0
                    if drmin < r_thr2:
                        break
                if drmin < r_thr and rmin not in used:
                    ff = crd2ff[rmin]
                    used.add(rmin)
                    curr_out += 'atom %4d %4d %-5s %-20s %4d %7.4f %2d\n'%(ff[0], ff[1], aname, resname, ff[2], ff[3], ff[4])
                else:
                    missing = True
            if missing:
                print("WARNING: Missing atom %s %7.4f %7.4f %7.4f in %s"%(aname, crd[0], crd[1], crd[2], pdbfile), file=sys.__stderr__)
                #break
        if not missing:
            reslib[(resname, natom)] = curr_out

    return reslib

def pdb_to_xyz(pdbfile, xyzfile, tpfile, use_element=False):
    d = read_template(tpfile)
    residues = convertpdb(pdbfile, d)

    t = tinkerxyz()

    iat = 0
    frame = {}
    for res in residues:
        for atom in res[1:]:
            iat += 1
            t.top_name[iat] = atom[0]
            t.top_type[iat] = atom[2]
            t.top_conn[iat] = []
            frame[iat] = list(atom[1])
    fout = open(xyzfile, 'w')
    t.write_xyz(frame, fout)
    fout.close()

def xyz_to_pdb(pdbfile, xyzfile, tpfile):
    d = read_template(tpfile)
    d2 = {}

    for res in d:
        tps = tuple(sorted(d[res][1]))
        tpdict = {}
        for tp, nm in zip(d[res][1], d[res][0]):
            if tp not in tpdict:
                tpdict[tp] = [nm]
            #elif nm not in tpdict[tp]:
            else:
                tpdict[tp].append(nm)
            pass
        #d2[tps] = [res, dict(zip(d[res][1], d[res][0]))]
        d2[tps] = [res, tpdict]

    t = tinkerxyz()
    t.read_xyz(xyzfile)
    t.read_arc(xyzfile, t.append_frame)
    residues = []
    reslengs = list(sorted(list(set([K[1] for K in d])), reverse=True))
    atypes = [t.top_type[K] for K in range(1, t.top_natom+1)]
    atom_found = np.zeros(t.top_natom)
    ilast = 0
    for i in range(t.top_natom):
        if i < ilast:
            continue
        for nres in reslengs:
            if i+nres > t.top_natom:
                continue
            cur_types = tuple(sorted(atypes[i:(i+nres)]))
            if cur_types in d2:
                residues.append((i, i+nres)) 
                ilast = i+nres
    pdb_i1 = 0
    pdb_i2 = 0
    pdb_out = 'AUTHOR    CONVERTED FROM TINKER XYZ\n'
    PDBFMT = 'ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f\n'
    for idxs in residues:
        atrange = list(range(idxs[0]+1, idxs[1]+1))
        cur_types = tuple(sorted(list(atypes[idxs[0]: idxs[1]])))
        curr_res = d2[cur_types]
        resname = curr_res[0][0]
        if len(resname) > 3:
            resname = resname[-3:]
        tpdict = curr_res[1]
        pdb_i2 += 1

        used_type = []
        for i in atrange:
            atype = t.top_type[i]
            atname = tpdict[t.top_type[i]][used_type.count(atype)]
            R = t.frames[0][i]
            pdb_i1 += 1
            pdb_out += PDBFMT%(pdb_i1%100000, atname, resname, ' ', pdb_i2%10000, R[0], R[1], R[2])
            used_type.append(atype)
    pdb_out += 'END\n'
    fout = open(pdbfile, 'w')
    fout.write(pdb_out)
    fout.close()

def read_prm(fprm):
    dprm = {}
    with open(fprm, 'r') as fh:
        for line in fh:
            w = line.split()
            nw = len(w)
            if line.startswith('atom') and nw >= 5:
                dprm[int(w[1])] = line
    return dprm

def xyz_to_lib(fprm, fxyz, flib, resname='UNK'):
    RES_SEP = '\natom    0    0 --    --------                0  0.0000  0\n'
    dprm = read_prm(fprm)
    t = tinkerxyz()
    t.read_xyz(fxyz)
    outl = RES_SEP
    for i in range(1, t.top_natom+1):
        atype = t.top_type[i]
        if atype in dprm:
            w = dprm[atype].split()
            i1 = dprm[atype].index(w[4])
            i2 = dprm[atype][i1:].index(w[5]) + i1
            idxs = [0]
            for iw in range(len(w)):
                i1 = idxs[-1]
                i2 = dprm[atype][i1:].index(w[iw]) + i1
                idxs.append(i2)
            i1 = idxs[5]
            i2 = idxs[-2]
            #print(i1, i2)
            outl += dprm[atype][:i1] + ' "%s" '%resname + dprm[atype][i2:]
        else:
            return
    with open(flib, 'a') as fh:
        fh.write(outl)

def pdbxyz_to_lib(pdb_list, xyz_list, biotype='pro', resname='UNK', outfile=None):
    reslib = {}
    for fxyz, fpdb in zip(xyz_list, pdb_list):
        reslib.update(train_tplt(fpdb, fxyz, biotype))
        pass
    outp = train_tplt.__doc__
    sorted_res = sorted(list(reslib), key=lambda t: (len(t[0]), t[0]))
    lastres = ('', 0)
    RES_SEP = '\natom    0    0 --    --------                0  0.0000  0\n'
    for res in (sorted_res):
        if res[0] == lastres[0]:
            outp += RES_SEP
        outp += (reslib[res])
        lastres = res
    if outfile is not None:
        fout = open(outfile, 'a')
    else:
        fout= sys.__stdout__
    print(outp, file=fout)
    if outfile is not None:
        fout.close()

def get_new_resname(inpfile, pref='R'):
    '''
    find new residue name not in inpfile
    '''
    res_def = read_template(inpfile)
    resnames = set([_[0] for _ in res_def])
    resname = 'UNK'
    for i in range(1, 100):
        r = '%s%02d'%(pref, i)
        if r not in resnames:
            resname = r
            break
    return resname

def main():
    parser = argparse.ArgumentParser(prog='pdbxyz.py', usage='%(prog)s [options]', description=__doc__)
    parser.add_argument('-m', nargs='?', metavar='mode', choices=['pdbxyz2lib', 'pdb2xyz', 'xyz2pdb', 'xyz2lib'], default='pdb2xyz', help='[pdbxyz2lib|pdb2xyz|xyz2pdb] pdbxyz2lib: Generate library file from PDB and XYZ; pdb2xyz: Convert PDB to XYZ using library file; xyz2pdb: Convert XYZ to PDB')
    parser.add_argument('-p', nargs='+', metavar='PDB', required=True, help='PDB file(s)')
    parser.add_argument('-x', nargs='+', metavar='XYZ', required=True, help='TINKER XYZ file(s)')
    parser.add_argument('-t', nargs='?', metavar='TEMPLATE', required=False, help='Library file')
    parser.add_argument('-k', nargs='?', metavar='KEY', required=False, help='Library file')
    parser.add_argument('-r', nargs='?', metavar='RESNAME', required=False, help='residue name for xyz2lib')
    parser.add_argument('--nuc', action='store_true', help='Add 5/3 to terminal residue names when generating template')
    parser.add_argument('--app', action='store_true', help='Append to template file; if not specified, will print to stdout')

    v = (parser.parse_args(sys.argv[1:]))
    biotype = 'pro'
    if v.nuc:
        biotype = 'nuc'
    if v.m == 'pdbxyz2lib':
        if v.app:
            outfile = v.t
        else:
            outfile = None
        pdbxyz_to_lib(v.p, v.x, biotype, outfile)
    elif v.m == 'pdb2xyz':
        for i in range(min(len(v.p), len(v.x))):
            pdb_to_xyz(v.p[i], v.x[i], v.t)
    elif v.m == 'xyz2pdb':
        for i in range(min(len(v.p), len(v.x))):
            xyz_to_pdb(v.p[i], v.x[i], v.t)
    elif v.m == 'xyz2lib':
        if v.r is None:
            resname = get_new_resname(v.t)
        else:
            resname = v.r
        xyz_to_lib(v.k, v.x[0], v.t, resname=resname)

main()

