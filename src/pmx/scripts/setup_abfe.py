#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import argparse
import os
import shutil
import numpy as np
from pmx.model import Model, merge_models, assign_masses_to_model, double_box
from pmx.alchemy import AbsRestraints, couple_mol, decouple_mol
from pmx.forcefield import Topology, merge_atomtypes
from pmx import gmx
from copy import deepcopy
from time import sleep
from cli import check_unknown_cmd


# TODO: allow providing indices for restraint
# TODO: build systems with pmx.gmx
# TODO: build folder structure and mdp files for equil or nonequil calcs
# TODO: make optional keeping intramolinteractions with singlebox setup
# TODO: singlebox allow using bLogestAxis

"""
Script to setup absolute binding free energy calculations for ligands binding to proteins.
"""


def parse_options():
    parser = argparse.ArgumentParser(description='''
describe...
''')
    parser.add_argument('-pt',
                        metavar='protop',
                        dest='pro_top',
                        type=str,
                        help='Input topology file for the protein. '
                        'Default is "protein.top".',
                        default='protein.top')
    parser.add_argument('-lt',
                        metavar='ligtop',
                        dest='lig_top',
                        type=str,
                        help='Input topology file for the ligand. It is '
                        'expected that all params needed for the ligand '
                        'are explicitly defined in this file. '
                        'Default is "ligand.itp".',
                        default='ligand.itp')
    parser.add_argument('-pc',
                        metavar='procrd',
                        dest='pro_crd',
                        type=str,
                        help='Input structure file in PDB or GRO format '
                        'for the protein. Default is "protein.gro".',
                        default='protein.gro')
    parser.add_argument('-lc',
                        metavar='ligcrd',
                        dest='lig_crd',
                        type=str,
                        help='Input structure file in PDB or GRO format '
                        'for the ligand. Default is "ligand.gro".',
                        default='ligand.gro')
    parser.add_argument('--build',
                        dest='build',
                        help='Whether to build the system (editconf, solvate, '
                        'genion) with a standard setup once the input files '
                        '(top, gro) are ready.',
                        default=False,
                        action='store_true')
    parser.add_argument('--singlebox',
                        dest='singlebox',
                        help='Whether to use the double-system single-box '
                        'setup.',
                        default=False,
                        action='store_true')
    parser.add_argument('--seed',
                        metavar='int',
                        dest='seed',
                        help='Random seed to use when picking atoms for the '
                        'restraints. The automated restraints selection is '
                        'stochastic, so if you want to have a reproducible '
                        'behaviour, provide a random seed.',
                        default=None,
                        type=int)

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args


def main(args):

    # Import GRO and TOP files
    lig = Model(args.lig_crd, renumber_residues=False)
    pro = Model(args.pro_crd, renumber_residues=False)
    ligtop = Topology(args.lig_top, is_itp=True, assign_types=True, self_contained=True)
    protop = Topology(args.pro_top, assign_types=False)

    # check lig has only 1 residue
    assert len(lig.residues) == 1

    # assign masses to ligand atoms
    assign_masses_to_model(lig, ligtop)

    # add position restraints to ligand topology if missing
    if not ligtop.footer and ligtop.has_posre is False:
        ligtop.make_posre()

    # ----------------------
    # merge gro into complex
    # ----------------------
    com = merge_models(lig, pro)
    com.box = pro.box
    com.renumber_atoms()

    # ------------------------
    # make topology of complex
    # ------------------------
    comtop = deepcopy(protop)
    comtop.include_itps.append(("ligand.itp", "top"))
    comtop.system = 'Complex'
    comtop.molecules.insert(0, [ligtop.name, 1])

    # merge atomtypes
    comtop.atomtypes = merge_atomtypes(ligtop.atomtypes, protop.atomtypes)

    # ------------------------------------------
    # figure out restraints (and save top entry)
    # ------------------------------------------
    restraints = AbsRestraints(pro, lig, seed=args.seed)

    # write restraints info
    restraints.write_summary()

    if args.build is False:
        # write restraints section only if we are not setting up the system
        # solvate/genion do not work well with topologies that contain
        # this restraints section
        comtop.ii = restraints.make_ii()

    # -------------------------------
    # save coordinates/topology files
    # -------------------------------
    com.write('complex.gro')
    ligtop.write('ligand.itp', stateBonded='A', write_atypes=False, posre_include=True)
    comtop.write('complex.top', stateBonded='A')

    # ---------------------------------------------------
    # gmx setup: editconf, solvate, genion, mdp files etc
    # ---------------------------------------------------
    if args.build is True:
        goodtogo = _check_topology_has_all_ff_info(protop)
        if goodtogo is False:
            print('Skipping system building ...')
            sleep(3)

        # ================
        # single-box setup
        # ================
        if args.singlebox is True and goodtogo is True:

            # create double-system single-box
            # the second ligand (lig) is inserted before the complex in the
            # output gro file (easier to modify topology this way atm)
            mout = double_box(m1=lig, m2=com, r=2.5, d=1.5,
                              bLongestAxis=False, verbose=False)
            mout.write('singlebox.gro')

            # identify an atom to restrain: this will be used in both ligands to
            # keep ligand and complex apart

            # reload because indices have changed after lig was merged into complex
            lig = Model(args.lig_crd, renumber_residues=False)
            assign_masses_to_model(lig, ligtop)
            posre_idx = _find_atom_near_com(lig)

            # create ligand topology with B-state
            # ------------------------------------
            ligtopAB = deepcopy(ligtop)
            decouple_mol(ligtopAB)
            _add_fixed_posre(top=ligtopAB, n=posre_idx)
            ligtopAB.write('ligandAB.itp', stateBonded='A', stateTypes='AB', stateQ='AB',
                           write_atypes=False, posre_include=True)

            # create second ligand topology with B-state
            # ------------------------------------------
            ligtopBA = deepcopy(ligtop)
            couple_mol(ligtopBA)
            _add_fixed_posre(top=ligtopBA, n=posre_idx)
            ligtopBA.name = '{}2'.format(ligtopBA.name)
            ligtopBA.write('ligandBA.itp', stateBonded='A', stateTypes='AB', stateQ='AB',
                           write_atypes=False, posre_include=True)

            # create system topology
            # ----------------------
            doubletop = deepcopy(comtop)
            doubletop.include_itps = [('ligandAB.itp', 'top'), ('ligandBA.itp', 'top')]
            doubletop.molecules.insert(1, [ligtopBA.name, 1])
            # merge atomtypes (basically add the dummies)
            doubletop.atomtypes = merge_atomtypes(doubletop.atomtypes, ligtopAB.atomtypes)
            # keep intramolecular interactions
            doubletop.make_nonbond_params(rule=2)
            # save topology
            doubletop.write('singlebox.top', stateBonded='A')

            # run gromacs setup
            # -----------------
            gmx.solvate(cp='singlebox.gro', cs='spc216.gro', p='singlebox.top', o='solvate.gro')
            gmx.write_mdp(mdp='enmin', fout='genion.mdp')
            gmx.grompp(f='genion.mdp', c='solvate.gro', p='singlebox.top', o='genion.tpr', maxwarn=1)
            gmx.genion(s='genion.tpr', p='singlebox.top', o='genion.gro', conc=0.15, neutral=True)

        # ==================================
        # standard setup with separate boxes
        # ==================================
        elif args.singlebox is False and goodtogo is True:

            # Setup complex
            # -------------
            os.mkdir('complex')
            os.chdir('complex')

            shutil.copy('../complex.gro', '.')
            shutil.copy('../complex.top', '.')
            shutil.copy('../ligand.itp', '.')

            gmx.editconf(f='complex.gro', o='editconf.gro', bt='cubic', d=1.2)
            gmx.solvate(cp='editconf.gro', cs='spc216.gro', p='complex.top', o='solvate.gro')
            gmx.write_mdp(mdp='enmin', fout='genion.mdp')
            gmx.grompp(f='genion.mdp', c='solvate.gro', p='complex.top', o='genion.tpr', maxwarn=1)
            gmx.genion(s='genion.tpr', p='complex.top', o='genion.gro', conc=0.15, neutral=True)

            # add restraints to topology
            comtop = Topology('complex.top', assign_types=False)
            comtop.ii = restraints.make_ii()
            comtop.write('complex.top', stateBonded='A')

            os.chdir('../')

            # Setup ligand
            # ------------
            os.mkdir('ligand')
            os.chdir('ligand')

            lig.write('ligand.gro')
            # setup topology
            # use the same protein ff for the water/ions in ligand sims
            ligtop.is_itp = False  # now we want to write it as top file
            ligtop.forcefield = protop.forcefield
            ligtop.footer = ['#include "{ff}.ff/tip3p.itp"'.format(ff=ligtop.forcefield),
                             '#ifdef POSRES_WATER',
                             '[ position_restraints ]',
                             '1    1       1000       1000       1000',
                             '#endif',
                             '#include "{ff}.ff/ions.itp"'.format(ff=ligtop.forcefield)]

            ligtop.write('ligand.top', stateBonded='A', write_atypes=True,
                         posre_include=True)

            # run gromacs setup
            gmx.editconf(f='ligand.gro', o='editconf.gro', bt='cubic', d=1.2)
            gmx.solvate(cp='editconf.gro', cs='spc216.gro', p='ligand.top', o='solvate.gro')
            gmx.write_mdp(mdp='enmin', fout='genion.mdp')
            gmx.grompp(f='genion.mdp', c='solvate.gro', p='ligand.top', o='genion.tpr', maxwarn=1)
            gmx.genion(s='genion.tpr', p='ligand.top', o='genion.gro', conc=0.15, neutral=True)

            os.chdir('../')


    print('\n\n          ********** Setup Completed **********\n\n')


def _check_topology_has_all_ff_info(top):
    '''This is needed in particular for host-guest systems where the host has
    not passed through pdb2gmx, and might not contain the ff info needed for
    solvate/genion/grompp (i.e. water/ion params to use).
    '''
    good = True

    if top.forcefield == '':
        good = False
        print('\nthe input topology file {} does not seem to '
              'contain information on the forcefield '
              'parameters to use.\nThese are usually defined via the following '
              'statement at the top of the topology file:\n'
              '#include "amber99sb-ildn.ff/forcefield.itp"\n'
              'Without this information it is not possible '
              'to setup the system in Gromacs.\n'.format(top.filename))

    if 'ions.itp' not in " ".join(top.footer):
        good = False
        print('\nthe input topology file {} does not seem to '
              'contain information on the ion '
              'parameters to use.\nThese are usualy defined '
              'just above [ system ], e.g.:\n'
              '#include "amber99sb-ildn.ff/ions.itp"\n'
              'Without this information it is not possible '
              'to setup the system in Gromacs.\n'.format(top.filename))

    if 'tip' not in " ".join(top.footer) and 'spc' not in " ".join(top.footer):
        good = False
        print('\nthe input topology file {} does not seem to '
              'contain information on the water '
              'parameters to use.\nThese are usualy defined '
              'just above [ system ], e.g.:\n'
              '#include "amber99sb-ildn.ff/tip3p.itp"\n'
              'Without this information it is not possible '
              'to setup the system in Gromacs.\n'.format(top.filename))

    return good


def _add_fixed_posre(top, n):
    ''' n : index of atom to restrain
    '''
    posre = ['[ position_restraints ]',
             ';  ai    funct            c0            c1          c2',
             '   {n}      1      1000   1000   1000'.format(n=n)]

    top.footer = posre + top.footer


def _find_atom_near_com(m):
    # get com
    com = np.array(m.com(vector_only=True))
    coords = [np.array(a.x) for a in m.atoms]
    # get all distances
    distances = [np.linalg.norm(c-com) for c in coords]
    # get index of atom closest to com. Add 1 since gromacs uses idx from 1
    com_idx = np.argmin(distances)+1
    return com_idx


def entry_point():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    entry_point()
