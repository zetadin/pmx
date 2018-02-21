#!/usr/bin/env python

# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

"""Program to insert mutated residues in structure files for
free energy simulations.
"""

import sys
import os
import argparse
from pmx import library
from pmx.model import Model
from pmx.parser import read_and_format
from pmx.geometry import Rotation, nuc_super, bb_super
from pmx.mutdb import read_mtp_entry
from pmx.utils import get_ff_path, UnknownResidueError, RangeCheckError, mtpError

# set GMXLIB
path = os.path.abspath(library.__file__)
dir_path = os.path.dirname(path)
data_path = dir_path + '/data/mutff45dna'
os.environ['GMXLIB'] = data_path

# TODO: move dicts to library?
ext_one_letter = {
    'ALA':  'A',
    'ARG':  'R',
    'ASN':  'N',
    'ASP':  'D',
    'ASPH': 'B',
    'ASPP': 'B',
    'ASH':  'B',
    'CYS':  'C',
    'CYS2': 'C',
    'CYN':  'C',
    'CYX':  'C',
    'CYM':  'CM',
    'CYSH': 'C',
    'GLU':  'E',
    'GLUH': 'J',
    'GLUP': 'J',
    'GLH':  'J',
    'GLN':  'Q',
    'GLY':  'G',
    'HIS':  'H',
    'HIE':  'X',
    'HISE': 'X',
    'HSE':  'X',
    'HIP':  'Z',
    'HSP':  'Z',
    'HISH': 'Z',
    'HID':  'H',
    'HSD':  'H',
    'ILE':  'I',
    'LEU':  'L',
    'LYS':  'K',
    'LYSH': 'K',
    'LYP':  'K',
    'LYN':  'O',
    'LSN':  'O',
    'MET':  'M',
    'PHE':  'F',
    'PRO':  'P',
    'SER':  'S',
    'THR':  'T',
    'TRP':  'W',
    'TYR':  'Y',
    'VAL':  'V',
}

dna_names = {
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

rna_names = {
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

dna_one_letter = {'A': 'adenosine',
                  'C': 'cytosine',
                  'G': 'guanine',
                  'T': 'thymine'}

rna_one_letter = {'A': 'adenosine',
                  'C': 'cytosine',
                  'G': 'guanine',
                  'U': 'uracil'}

# ================
# Helper functions
# ================
def _print_sorted_dict(d):
    for key in sorted(d.iterkeys()):
        print("{0:>5}{1:>15}".format(key, d[key]))


def _int_input():
    inp = raw_input()
    try:
        inp = int(inp)
        return inp
    except:
        print 'You entered "%s" -> Try again' % inp
        return None


def _check_residue_name(res):
    if res.resname == 'LYS':
        if res.has_atom('HZ3'):
            res.set_resname('LYP')
    elif res.resname == 'HIS':
        if res.has_atom('HD1') and \
           res.has_atom('HE2'):
            res.set_resname('HIP')
        elif res.has_atom('HD1') and not res.has_atom('HE2'):
            res.set_resname('HID')
        elif not res.has_atom('HD1') and res.has_atom('HE2'):
            res.set_resname('HIE')
    elif res.resname == 'ASP':
        if res.has_atom('HD2'):
            res.set_resname('ASH')
    elif res.resname == 'GLU':
        if res.has_atom('HE2'):
            res.set_resname('GLH')
    elif res.resname == 'CYS':
        if not res.has_atom('HG'):
            print >>sys.stderr, ' Cannot mutate SS-bonded Cys %d' % res.id


def _check_OPLS_LYS(res):
    if res.has_atom('HZ3'):
        return('K')
    else:
        return('O')


def _check_residue_range(m, idx):
    valid_ids = range(1, len(m.residues)+1)
    if idx not in valid_ids:
        return False
    return True


def _ask_next():
    sys.stdout.write('\nApply another mutation [y/n]? ')
    res = raw_input().lower()
    if res == 'y':
        return True
    elif res == 'n':
        return False
    else:
        return _ask_next()


def _convert_aa_name(aa):
    # firstly, some special deprotonated cases
    if aa == 'CM':
        return(aa.upper())
    elif len(aa) == 1:
        return aa.upper()
    elif len(aa) in [3, 4]:
        return ext_one_letter[aa.upper()]
    else:
        raise UnknownResidueError(aa)


class InteractiveSelection:
    """Class containing fucntions related to the interactive selection of
    residues to be mutated.

    Parameters
    ----------
    m : Model object
        instance of pmx.model.Model
    ffpath : str
        path to forcefield files

    Attributes
    ----------
    mut_resid : int
        index of residue to be mutated
    mut_resname : str
        one-letter code of target residue

    """

    def __init__(self, m, ff):
        self.m = m
        self.ffpath = get_ff_path(ff)

        self.mut_resid = self.select_residue()
        self.mut_resname = self.select_mutation()

    def select_residue(self):
        """Ask for the residue id to mutate.
        """
        valid_ids = range(1, len(self.m.residues)+1)
        print '\nSelect residue to mutate:'
        for i, r in enumerate(self.m.residues):
            if r.moltype not in ['water', 'ion']:
                sys.stdout.write('%6d-%s-%s' % (r.id, r.resname, r.chain_id))
                if r.id % 6 == 0:
                    print("")
        print("")
        selected_residue_id = None
        while not selected_residue_id:
            sys.stdout.write('Enter residue number: ')
            selected_residue_id = _int_input()
            if selected_residue_id is not None and selected_residue_id not in valid_ids:
                print('Residue id %d not in range %d-%d -> Try again' %
                      (selected_residue_id, 1, len(self.m.residues)))
                selected_residue_id = None
        return selected_residue_id

    def select_mutation(self):
        """Ask which residue to mutate to.
        """

        residue = self.m.residues[self.mut_resid - 1]
        if residue.moltype == 'protein':
            aa = self.select_aa_mutation(residue)
        elif residue.moltype in ['dna', 'rna']:
            aa = self.select_nuc_mutation(residue)
        return aa

    def select_aa_mutation(self, residue):
        """Selection for protein residues.
        """

        _check_residue_name(residue)
        print '\nSelect new amino acid for %s-%s: ' % (residue.id, residue.resname)
        sys.stdout.write('Three- or one-letter code (or four-letter for ff specific residues): ')
        if residue.resname in ['HIE', 'HISE', 'HSE']:
            rol = 'X'
        elif residue.resname in ['HIP', 'HISH', 'HSP']:
            rol = 'Z'
        elif residue.resname in ['GLH', 'GLUH', 'GLUP']:
            rol = 'J'
        elif residue.resname in ['ASH', 'ASPH', 'ASPP']:
            rol = 'B'
        elif residue.resname in ['LYN', 'LYS', 'LSN']:
            rol = 'O'
        else:
            rol = library._one_letter[residue.resname]
        aa = None
        ol = library._aacids_dic.keys()
        tl = library._aacids_dic.values()
        ffpathlower = self.ffpath.lower()
        if 'amber' in ffpathlower:
                ol = library._aacids_ext_amber.keys()
                tl = library._aacids_ext_amber.values()
        if 'opls' in ffpathlower:
                ol = library._aacids_ext_oplsaa.keys()
                tl = library._aacids_ext_oplsaa.values()+['ASPP', 'GLUP', 'LSN']
        if 'charmm' in ffpathlower:
                ol = library._aacids_ext_charmm.keys()
                tl = library._aacids_ext_charmm.values()

        while aa is None:
            aa = raw_input().upper()
            # some special residues:
            #   CM - deprotonated cysteine
            #   YM - deprotonated tyrosine
            if aa == 'CM':
                sys.stdout.write('Special case for deprotonated residue')
            elif len(aa) != 1 and len(aa) != 3 and len(aa) != 4:
                sys.stdout.write('Nope!\nThree- or one-letter code (or four-letter for ff specific residues): ')
                aa = None
            elif (len(aa) == 1 and aa not in ol+['B', 'J', 'O', 'X', 'Z']) or \
                 (len(aa) == 3 and aa not in tl) or \
                 (len(aa) == 4 and aa not in tl):
                sys.stdout.write('Unknown aa "%s"!\nThree- or one-letter code (or four-letter for ff specific residues): ' % aa)
                aa = None
            if aa and (len(aa) == 3 or len(aa) == 4):
                aa = ext_one_letter[aa]
        print('Will apply mutation %s->%s on residue %s-%d'
              % (rol, aa, residue.resname, residue.id))
        return aa

    def select_nuc_mutation(self, residue):
        """Selection for nucleic acids.
        """
        aa = None
        print('\nSelect new base for %s-%s: ' % (residue.id, residue.resname))
        sys.stdout.write('One-letter code: ')
        while aa is None:
            aa = raw_input().upper()
            if residue.moltype == 'dna' and aa not in ['A', 'C', 'G', 'T']:
                sys.stdout.write('Unknown DNA residue "%s"!\nOne-letter code: ' % aa)
                aa = None
            elif residue.moltype == 'rna' and aa not in ['A', 'C', 'G', 'U']:
                sys.stdout.write('Unknown RNA residue "%s"!\nOne-letter code: ' % aa)
                aa = None
            if aa:
                print('Will apply mutation %s->%s on residue %s-%d'
                      % (residue.resname[1], aa, residue.resname, residue.id))
            return aa


def rename_to_match_library(res):
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


def rename_back(res, name_hash):
    for atom in res.atoms:
        atom.name = name_hash[atom.name]


def set_conformation(old_res, new_res, rotdic):
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
                hybrid_residue_name = rna_names[dict_key]
            else:
                hybrid_residue_name = dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
        elif a.name == 'H5T':
            r1 = firstLetter+residue.resname[1]+'5'
            r2 = firstLetter+new_nuc_name+'5'
            dict_key = r1+'_'+r2
            if residue.moltype == 'rna':
                hybrid_residue_name = rna_names[dict_key]
            else:
                hybrid_residue_name = dna_names[dict_key]
            return(hybrid_residue_name, residue.resname[1], new_nuc_name)
    hybrid_residue_name = residue.resname+new_nuc_name
    return(hybrid_residue_name, residue.resname[1], new_nuc_name)


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


def apply_aa_mutation(m, residue, new_aa_name, mtp_file, infileB=None):

    if residue.resname == 'ILE':
        rename_ile(residue)
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
    hash1 = rename_to_match_library(residue)
    hash2 = rename_to_match_library(hybrid_res)
    set_conformation(residue, hybrid_res, rotdic)
    if infileB is not None:
        print("log_> Set Bstate geometry according to the provided structure")
        mB = Model(infileB, bPDBTER=True)
        rename_atoms_to_gromacs(mB)
        mB.nm2a()
        residueB = mB.residues[residue.id-1]
        bb_super(residue, residueB)
        for atom in hybrid_res.atoms:
            if (atom.name[0] == 'D') or atom.name.startswith('HV'):
                for atomB in residueB.atoms:
                    if atomB.name == hybrid_res.morphes[atom.name]['n1']:
                        atom.x = atomB.x
    rename_back(residue, hash1)
    rename_back(hybrid_res, hash2)
    # VG rename residue atoms back

    m.replace_residue(residue, hybrid_res)
    print('log_> Inserted hybrid residue %s at position %d (chain %s)' %
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id))


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


def get_hybrid_residue(residue_name, mtp_file='ffamber99sb.mtp'):
    print('log_> Scanning database for %s ' % residue_name)
    resi, bonds, imps, diheds, rotdic = read_mtp_entry(residue_name,
                                                       filename=mtp_file,
                                                       version='new')
    if len(resi.atoms) == 0:
        raise(mtpError("Hybrid residue %s not found in %s" %
                       (residue_name, mtp_file)))
    return resi, bonds, imps, diheds, rotdic


def rename_ile(residue):
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


def rename_atoms_to_gromacs(m):
    for atom in m.atoms:
        if atom.name[0].isdigit():
            atom.name = atom.name[1:]+atom.name[0]





def parse_options():
    parser = argparse.ArgumentParser(description='''
This script applies mutations of residues in a structure file for subsequent
free energy calculations. It supports mutations to protein, DNA, and RNA
molecules.

The mutation information and dummy placements are taken from the hybrid residue
database "mutres.mtp". The best way to use this script is to take a pdb/gro file
that has been written with pdb2gmx with all hydrogen atoms present.

The program can either be executed interactively or via script. The script file
simply has to consist of "resi_number target_residue" pairs.

The script uses an extended one-letter code for amino acids to account for
different protonation states. Use the --resinfo flag to print the dictionary.

Currently available force fields:
    - amber99sb-star-ildn-mut (Best & Hummer, 2009; Lindorff-Larsen et al, 2010)
    - charmm36m-mut (Best et al, 2012)
''', formatter_class=argparse.RawTextHelpFormatter)

    ff_choices = ['amber99sb-star-ildn-mut', 'charmm36m-mut.ff',
                  'amber99sb-star-ildn-bsc1-mut', 'amber14sb-mut']

    parser.add_argument('-f',
                        metavar='infile',
                        dest='infile',
                        type=str,
                        help='Input structure file in PDB or GRO format. '
                        'Default is "protein.pdb"',
                        default='protein.pdb')
    parser.add_argument('-fB',
                        metavar='infileB',
                        dest='infileB',
                        type=str,
                        help='Input structure file of the B state in PDB '
                        'or GRO format (optional).',
                        default=None)
    parser.add_argument('-o',
                        metavar='outfile',
                        dest='outfile',
                        type=str,
                        help='Output structure file in PDB or GRO format. '
                        'Default is "mutant.pdb"',
                        default='mutant.pdb')
    parser.add_argument('-ff',
                        metavar='ff',
                        dest='ff',
                        type=str.lower,
                        help='Force field to use. '
                        'Default is "amber99sb-star-ildn-mut"',
                        choices=ff_choices,
                        default='amber99sb-star-ildn-mut')
    parser.add_argument('--script',
                        metavar='script',
                        dest='script',
                        type=str,
                        help='Text file with list of mutations (optional).',
                        default=None)
    parser.add_argument('--resinfo',
                        dest='resinfo',
                        help='Show the list of 3-letter -> 1-letter residues',
                        default=False,
                        action='store_true')

    args, unknown = parser.parse_known_args()

    # ------------------
    # residue dictionary
    # ------------------
    if args.resinfo is True:
        moltype = Model(args.infile).moltype
        if moltype == 'protein':
            print('\n ---------------------------')
            print(' Protein residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(ext_one_letter)
            print(' ---------------------------\n')
        elif moltype == 'dna':
            print('\n -----------------------')
            print(' DNA residues dictionary')
            print(' -----------------------')
            _print_sorted_dict(dna_one_letter)
            print(' -----------------------\n')
        elif moltype == 'rna':
            print('\n -----------------------')
            print(' RNA residues dictionary')
            print(' -----------------------')
            _print_sorted_dict(rna_one_letter)
            print(' -----------------------\n')
        else:
            print('\n ---------------------------')
            print(' Protein residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(ext_one_letter)
            print(' ---------------------------')
            print(' DNA residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(dna_one_letter)
            print(' ---------------------------')
            print(' RNA residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(rna_one_letter)
            print(' ---------------------------\n')
        exit()
    else:
        return args


def main(args):

    # input variables
    infile = args.infile
    infileB = args.infileB
    outfile = args.outfile
    ff = args.ff
    script = args.script

    # initialise Model
    m = Model(infile, bPDBTER=True)
    rename_atoms_to_gromacs(m)
    m.nm2a()

    # if script is provided, do the mutations in that file
    if script is not None:
        mutations_to_make = read_and_format(script, "is")
        for mut in mutations_to_make:
            _check_residue_name(m.residues[mut[0]-1])
            apply_mutation(m=m,
                           mut_resid=mut[0],
                           mut_resname=mut[1],
                           ff=ff,
                           infileB=infileB)  # FIXME: what about bDNA?
    # if not provided, interactive selection
    else:
        do_more = True
        while do_more:
            sele = InteractiveSelection(m, ff)
            print "m: ", m
            print "mut: ", sele
            print "infileB: ", infileB
            apply_mutation(m=m,
                           mut_resid=sele.mut_resid,
                           mut_resname=sele.mut_resname,
                           ff=ff,
                           infileB=infileB)
            if not _ask_next():
                do_more = False

    m.write(outfile)
    print
    print 'mutations done...........'
    print


if __name__ == '__main__':
    args = parse_options()
    main(args)
