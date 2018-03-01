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

from __future__ import print_function
import sys
import argparse
from pmx import library
from pmx.model import Model
from pmx.parser import read_and_format
from pmx.utils import get_ff_path, ff_selection
from pmx.alchemy import mutate

# resinfo
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
        print('You entered "%s" -> Try again' % inp)
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
            print(' Cannot mutate SS-bonded Cys %d' % res.id, file=sys.stderr)


def _ask_next():
    sys.stdout.write('\nApply another mutation [y/n]? ')
    res = raw_input().lower()
    if res == 'y':
        return True
    elif res == 'n':
        return False
    else:
        return _ask_next()


# ===============================
# Class for interactive selection
# ===============================
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
        print('\nSelect residue to mutate:')
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
        print('\nSelect new amino acid for %s-%s: ' % (residue.id, residue.resname))
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
                aa = library._ext_one_letter[aa]
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


# =============
# Input Options
# =============
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

''', formatter_class=argparse.RawTextHelpFormatter)

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
                        help='Force field to use. If none is provided, \n'
                        'a list of available ff will be shown.',
                        default=None)
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
            _print_sorted_dict(library._ext_one_letter)
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
            _print_sorted_dict(library._ext_one_letter)
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

    # ------------
    # ff selection
    # ------------
    if args.ff is None:
        args.ff = ff_selection()

    return args


def main(args):

    # input variables
    infile = args.infile
    infileB = args.infileB
    outfile = args.outfile
    ff = args.ff
    script = args.script

    # initialise Model
    m = Model(infile, bPDBTER=True, for_gmx=True)

    # if script is provided, do the mutations in that file
    if script is not None:
        mutations_to_make = read_and_format(script, "is")
        for mut in mutations_to_make:
            _check_residue_name(m.residues[mut[0]-1])
            mutate(m=m,
                   mut_resid=mut[0],
                   mut_resname=mut[1],
                   ff=ff,
                   refB=infileB,
                   inplace=True,
                   verbose=True)
    # if not provided, interactive selection
    else:
        do_more = True
        while do_more:
            sele = InteractiveSelection(m, ff)
            mutate(m=m,
                   mut_resid=sele.mut_resid,
                   mut_resname=sele.mut_resname,
                   ff=ff,
                   refB=infileB,
                   inplace=True,
                   verbose=True)
            if not _ask_next():
                do_more = False

    m.write(outfile)
    print('')
    print('mutations done...........')
    print('')


def entry_point():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    entry_point()
