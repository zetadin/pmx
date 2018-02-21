#!/usr/bin/env python

"""This module contains functions that can be used to setup the hybrid
structures and topologies needed for alchemical free energy calculations.
"""

from model import Model
from utils import get_ff_path, RangeCheckError, mtpError, UnknownResidueError
from geometry import Rotation, nuc_super, bb_super
from mutdb import read_mtp_entry
import library
import os

__all__ = ['apply_mutation']

# ==============
# Main Functions
# ==============
def apply_mutation(m, mut_resid, mut_resname, ff, infileB):

    # check selection is valid
    if not _check_residue_range(m, mut_resid):
        raise RangeCheckError(mut_resid)
    # get th residue
    residue = m.residues[mut_resid - 1]

    # Determine which mtp file to use
    ffpath = get_ff_path(ff=ff)
    # DNA mutation
    if residue.moltype == 'dna':
        mtp_file = os.path.join(ffpath, 'mutres_dna.mtp')
    # RNA mutation
    elif residue.moltype == 'rna':
        mtp_file = os.path.join(ffpath, 'mutres_rna.mtp')
    # Protein mutation
    elif residue.moltype == 'protein':
        mtp_file = os.path.join(ffpath, 'mutres.mtp')
    else:
        raise(ValueError, 'Cannot undertand mutation type needed '
                          'from the input strcture provided')

    # Mutation if Protein
    if residue.moltype == 'protein':
        new_aa_name = _convert_aa_name(mut_resname)
        apply_aa_mutation(m, residue, new_aa_name, mtp_file, infileB)
    # Mutation if DNA or RNA
    elif residue.moltype in ['dna', 'rna']:
        new_nuc_name = mut_resname.upper()
        apply_nuc_mutation(m, residue, new_nuc_name, mtp_file)


def apply_aa_mutation(m, residue, new_aa_name, mtp_file, infileB=None):

    if residue.resname == 'ILE':
        _rename_ile(residue)
    olkey = _convert_aa_name(residue.resname)

    # olkey should contain the correct one letter name of the WT residue
    # however, due to the different namings of the residues in the FFs
    # Lys needs to be checked once again: in OPLS Lys is non-protonated,
    # while in the other FFs it is protonated
    if ('opls' in mtp_file) and ('LYS' in residue.resname):
        olkey = _check_OPLS_LYS(residue)

    hybrid_residue_name = olkey+'2'+new_aa_name
    print('log_> Residue to mutate: %d | %s | %s ' % (residue.id, residue.resname, residue.chain_id))
    print('log_> Mutation to apply: %s->%s' % (olkey, new_aa_name))
    print('log_> Hybrid residue name: %s' % hybrid_residue_name)
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)
    bb_super(residue, hybrid_res)

    # VG rename residue atoms
    hash1 = _rename_to_match_library(residue)
    hash2 = _rename_to_match_library(hybrid_res)
    _set_conformation(residue, hybrid_res, rotdic)
    if infileB is not None:
        print("log_> Set Bstate geometry according to the provided structure")
        mB = Model(infileB, bPDBTER=True, for_gmx=True)
        residueB = mB.residues[residue.id-1]
        bb_super(residue, residueB)
        for atom in hybrid_res.atoms:
            if (atom.name[0] == 'D') or atom.name.startswith('HV'):
                for atomB in residueB.atoms:
                    if atomB.name == hybrid_res.morphes[atom.name]['n1']:
                        atom.x = atomB.x
    _rename_back(residue, hash1)
    _rename_back(hybrid_res, hash2)
    # VG rename residue atoms back

    m.replace_residue(residue, hybrid_res)
    print('log_> Inserted hybrid residue %s at position %d (chain %s)' %
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id))


def apply_nuc_mutation(m, residue, new_nuc_name, mtp_file):

    hybrid_residue_name, resname1, resname2 = get_nuc_hybrid_resname(residue, new_nuc_name)
    print 'log_> Residue to mutate: %d | %s | %s ' % (residue.id, residue.resname, residue.chain_id)
    print 'log_> Mutation to apply: %s->%s' % (residue.resname[1], new_nuc_name)
    print 'log_> Hybrid residue name: %s' % hybrid_residue_name
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)

    nuc_super(residue, hybrid_res, resname1, resname2)
    for atom in hybrid_res.atoms:
        if atom.name[0] != 'D':
            atom.x = residue[atom.name].x
    m.replace_residue(residue, hybrid_res)
    print('log_> Inserted hybrid residue %s at position %d (chain %s)' %
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id))


def get_hybrid_residue(residue_name, mtp_file='ffamber99sb.mtp'):
    print('log_> Scanning database for %s ' % residue_name)
    resi, bonds, imps, diheds, rotdic = read_mtp_entry(residue_name,
                                                       filename=mtp_file,
                                                       version='new')
    if len(resi.atoms) == 0:
        raise(mtpError("Hybrid residue %s not found in %s" %
                       (residue_name, mtp_file)))
    return resi, bonds, imps, diheds, rotdic


def get_nuc_hybrid_resname(residue, new_nuc_name):

    if residue.moltype == 'dna':
        firstLetter = 'D'
    if residue.moltype == 'rna':
        firstLetter = 'R'

    # identify if the nucleotide is terminal
    for a in residue.atoms:
        if a.name == 'H3T':
            r1 = firstLetter+residue.resname[1]+'3'
            r2 = firstLetter+new_nuc_name+'3'
            dict_key = r1+'_'+r2
            if residue.moltype == 'rna':
                hybrid_residue_name = _rna_names[dict_key]
            else:
                hybrid_residue_name = _dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
        elif a.name == 'H5T':
            r1 = firstLetter+residue.resname[1]+'5'
            r2 = firstLetter+new_nuc_name+'5'
            dict_key = r1+'_'+r2
            if residue.moltype == 'rna':
                hybrid_residue_name = _rna_names[dict_key]
            else:
                hybrid_residue_name = _dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
    hybrid_residue_name = residue.resname+new_nuc_name
    return(hybrid_residue_name, residue.resname[1], new_nuc_name)


# ===============
# HelperFunctions
# ===============
def _check_OPLS_LYS(res):
    if res.has_atom('HZ3'):
        return('K')
    else:
        return('O')


def _convert_aa_name(aa):
    # firstly, some special deprotonated cases
    if aa == 'CM':
        return(aa.upper())
    elif len(aa) == 1:
        return aa.upper()
    elif len(aa) in [3, 4]:
        return library._ext_one_letter[aa.upper()]
    else:
        raise UnknownResidueError(aa)


def _check_residue_range(m, idx):
    valid_ids = range(1, len(m.residues)+1)
    if idx not in valid_ids:
        return False
    return True


def _rename_ile(residue):
    dic = {'CD':  'CD1',
           'HD1': 'HD11',
           'HD2': 'HD12',
           'HD3': 'HD13'
           }
    for key, value in dic.items():
        try:
            atom = residue[key]
            atom.name = value
        except:
            pass


def _rename_to_match_library(res):
    name_hash = {}
    atoms = res.atoms
    for atom in atoms:
        foo = atom.name
        # for serine
        if (atom.resname == 'SER') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('S2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
        # for cysteine
        if (atom.resname == 'CYS') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('C2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'

        name_hash[atom.name] = foo
    return name_hash


def _rename_back(res, name_hash):
    for atom in res.atoms:
        atom.name = name_hash[atom.name]


def _set_conformation(old_res, new_res, rotdic):
    old_res.get_real_resname()
    dihedrals = library._aa_dihedrals[old_res.real_resname]
    for key, lst in rotdic.items():
        new = new_res.fetchm(lst)
        rotdic[key] = new
    chis = []
    for key in rotdic.keys():
        at1, at2 = key.split('-')
        for d in dihedrals:
            if d[1] == at1 and d[2] == at2 \
                   and d[-1] != -1:
                chis.append(d)
    for d in chis:
        atoms = old_res.fetchm(d[:4])
        phi = atoms[0].dihedral(atoms[1], atoms[2], atoms[3])
        atoms2 = new_res.fetchm(d[:4])
        phi2 = atoms2[0].dihedral(atoms2[1], atoms2[2], atoms2[3])
        diff = phi-phi2
        a1, a2 = new_res.fetchm(d[1:3])
        key = a1.name+'-'+a2.name
        atoms = rotdic[key]
        rot = Rotation(a1.x, a2.x)
        for atom in atoms:
            atom.x = rot.apply(atom.x, diff)
    for atom in new_res.atoms:
        if (atom.name[0] != 'D') and (not atom.name.startswith('HV')):
            atom.x = old_res[atom.name].x


# ==============
# Data/Libraries
# ==============
_dna_names = {
    'DA5_DT5': 'D5K',
    'DA5_DC5': 'D5L',
    'DA5_DG5': 'D5M',
    'DT5_DA5': 'D5N',
    'DT5_DC5': 'D5O',
    'DT5_DG5': 'D5P',
    'DC5_DA5': 'D5R',
    'DC5_DT5': 'D5S',
    'DC5_DG5': 'D5T',
    'DG5_DA5': 'D5X',
    'DG5_DT5': 'D5Y',
    'DG5_DC5': 'D5Z',
    'DA3_DT3': 'D3K',
    'DA3_DC3': 'D3L',
    'DA3_DG3': 'D3M',
    'DT3_DA3': 'D3N',
    'DT3_DC3': 'D3O',
    'DT3_DG3': 'D3P',
    'DC3_DA3': 'D3R',
    'DC3_DT3': 'D3S',
    'DC3_DG3': 'D3T',
    'DG3_DA3': 'D3X',
    'DG3_DT3': 'D3Y',
    'DG3_DC3': 'D3Z',
    # False names to avoid an error
    'DG3_DG3': 'FOO',
    'DC3_DC3': 'FOO',
    'DA3_DA3': 'FOO',
    'DT3_DT3': 'FOO',
    'DG5_DG5': 'FOO',
    'DC5_DC5': 'FOO',
    'DA5_DA5': 'FOO',
    'DT5_DT5': 'FOO',
    }

_rna_names = {
    'RA5_RU5': 'R5K',
    'RA5_RC5': 'R5L',
    'RA5_RG5': 'R5M',
    'RU5_RA5': 'R5N',
    'RU5_RC5': 'R5O',
    'RU5_RG5': 'R5P',
    'RC5_RA5': 'R5R',
    'RC5_RU5': 'R5S',
    'RC5_RG5': 'R5T',
    'RG5_RA5': 'R5X',
    'RG5_RU5': 'R5Y',
    'RG5_RC5': 'R5Z',
    'RA3_RU3': 'R3K',
    'RA3_RC3': 'R3L',
    'RA3_RG3': 'R3M',
    'RU3_RA3': 'R3N',
    'RU3_RC3': 'R3O',
    'RU3_RG3': 'R3P',
    'RC3_RA3': 'R3R',
    'RC3_RU3': 'R3S',
    'RC3_RG3': 'R3T',
    'RG3_RA3': 'R3X',
    'RG3_RU3': 'R3Y',
    'RG3_RC3': 'R3Z',
    # False names to avoid an error
    'RG3_RG3': 'FOO',
    'RC3_RC3': 'FOO',
    'RA3_RA3': 'FOO',
    'RU3_RU3': 'FOO',
    'RG5_RG5': 'FOO',
    'RC5_RC5': 'FOO',
    'RA5_RA5': 'FOO',
    'RU5_RU5': 'FOO',
    }
