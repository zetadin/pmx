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

"""This module contains the Model class. It can use
GROMACS routines to read and write structure files. Moreover it
allows to modify structure files in various ways. E.g.:

1. Rename atoms, residues, chains
2. Delete or add atoms, residues and chains

By default, all residues are renumbered from 1. This option can be deactivated
by setting renumber_residues to False.

The Model instance contains:

* model.atoms       -> list of atoms
* model.residues    -> list of residues
* model.chains      -> list of chains
* model.chdic       -> chain dictionary (chdic['A'] returns chain A)

Examples
--------
Basic usage:

    >>> model = Model('input.pdb')

Some useful methods:

    >>> # returns all backbone atoms
    >>> model.fetch_atoms(['CA','N','C'])
    >>> # returns all carbon atoms
    >>> model.fetch_atoms('C',how='byelem')
    >>> # return all atoms except hydrogens
    >>> model.fetch_atoms('H',how='byelem',inv=True)

    >>> # return all ALA,TRP and CYS residues
    >>> model.fetch_residues(['ALA','TRP','CYS'])

    >>> # returns the first 10 residues
    >>> rl = model.residues[:10]
    >>> # return the last residue from chain A
    >>> rl = model.chdic['A'].residues[-1]
    >>> # returns a list with the first residues of each chain
    >>> rl = map(lamda m: m.residues[0], model.chains)
    >>> # remove chain A
    >>> del model['A']
    >>> # write new structure file
    >>> model.write('output.pdb')

"""

from __future__ import absolute_import, print_function, division
import sys
import copy
from . import _pmx as _p
from . import library
from . import chain
from .atomselection import Atomselection
from .molecule import Molecule
from .atom import Atom


__all__ = ['Model']


class Model(Atomselection):
    """Model Class.

    Parameters
    ----------
    filename : str
        filename of input structure
    pdbline : bool(?)
        what is pdbline?
    renumber_atoms : bool, optional
        renumber all atoms from 1. Default is True.
    renumber_residues : bool, optional
        renumber all residues from 1. In this way, each residue will have a
        unique ID, also across chains. Default is True.
    rename_atoms : bool, optional
        rename atoms so to conform to Gromacs format. Default is False.
    scale_coords : A|nm, optional
        whether to enforce the units of the coordinates to be in A or in nm.
        By default, PDB coordinates are assumed to be in Angstrom (A) and
        GRO coordinates in nanometers (nm). If you read a PDB file but would
        like operate on nm coordinates, then select "nm". Viceversa, if you
        read a GRO file but would like to work on A coordinates, select "A".
        Note that if you read a PDB file in A coordinates and select "A",
        nothing happens (same for GRO and "nm" selection).
    bPDBTER : bool(?)
        whether to recognize TER lines and other chain breaks, e.g.
        discontinuous residue indices(?). Default is True.
    bNoNewID : bool(?)
        whether to assign new chain IDs? If True, new chain IDs starting
        with 'pmx' will be assigned(?). Only relevant if bPDBTER is True.
        Default is True.

    Attributes
    ----------
    title : str
        title of model. Default is 'PMX MODEL'.
    filename : str
        filename from which the Model was imported, otherwise None.
    chains : list
        list of Chain instances
    chdic : dict
        dict with chain IDs as keys and Chain instances as values
    residues : list
        list of molecules/residues
    unity : str
        coordinates unit, either 'A' or 'nm'.
    box : 2d array
        3x3 array containing the box vectors. See Gromacs manual, Table 3.1
    moltype : str
        Type of system: protein, dna, rna, or unknown if organic molecule or
        a mix of molecules are in the system.
    """
    def __init__(self, filename=None, pdbline=None, renumber_atoms=True,
                 renumber_residues=True, rename_atoms=False, scale_coords=None,
                 bPDBTER=True, bNoNewID=True,
                 **kwargs):

        Atomselection.__init__(self)
        self.title = 'PMX MODEL'
        self.filename = filename
        self.chains = []
        self.chdic = {}
        self.residues = []
        self.name = None  # FIXME/QUESTION: self.name not used anywhere? remove?
        self.id = 0
        self.have_bonds = 0  # FIXME/QUESTION: same as above: never used -> remove?
        self.box = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.unity = 'A'
        for key, val in kwargs.items():
            setattr(self, key, val)

        if filename is not None:
            self.read(filename=filename, bPDBTER=bPDBTER, bNoNewID=bNoNewID)
        if pdbline is not None:
            self.__readPDB(pdbline=pdbline)
        if self.atoms:
            self.unity = self.atoms[0].unity
            self.make_chains()
            self.make_residues()
        if self.residues and not self.atoms:
            self.al_from_resl()
            self.make_chains()
            self.make_residues()
        if self.chains and not self.residues:
            self.resl_from_chains()
            self.al_from_resl()
            self.make_chains()
            self.make_residues()
        if self.chdic and not self.chains:
            for key, val in self.chdic.items():
                self.chains.append(val)
            if not self.atoms and not self.residues:
                self.resl_from_chains()
                self.al_from_resl()
        if renumber_atoms is True:
            self.renumber_atoms()
        if renumber_residues is True:
            self.renumber_residues()
        if rename_atoms is True:
            self.rename_atoms_to_gmx()
        if scale_coords is not None:
            if scale_coords == 'A':
                self.nm2a()
            elif scale_coords == 'nm':
                self.a2nm()
            else:
                raise ValueError('unknown unit %s for coordinates' % scale_coords)

        self.assign_moltype()

    def __str__(self):
        s = '< Model: moltype=%s, nchain=%d, nres=%d, natom=%d >' %\
            (self.moltype, len(self.chains), len(self.residues),
             len(self.atoms))
        return s

    def writePIR(self, filename, title=""):
        """Prints sequence to screen in PIR format"""
        fp = open(filename, "w")
        if not title:
            title = '_'.join(self.title.split())
        print('>P1;%s' % title, file=fp)
        print('sequence:::::::::', file=fp)
        for i in range(len(self.chains) - 1):
            print('%s/' % self.chains[i].get_sequence(), file=fp)
        print('%s*' % self.chains[-1].get_sequence(), file=fp)
        fp.close()

    def writeFASTA(self, filename, title=""):
        """Prints sequence to screen in FASTA format"""
        fp = open(filename, "w")
        if not title:
            title = '_'.join(self.title.split())
        if len(self.chains) == 1:
            print('> %s' % title, file=fp)
            print(self.chains[0].get_sequence(), file=fp)
        else:
            for chain in self.chains:
                print('> %s_chain_%s' % (title, chain.id), file=fp)
                print(chain.get_sequence(), file=fp)

    # FIXME/TODO: this function is overwriting the write function inherited from
    # Atomselection. Should we keep only one of the 2?
    # writePIR and writeFASTA could be moved to Atomselection so that all
    # write functions are in one place
    def write(self, fn, title='', nr=1, bPDBTER=False, bAssignChainIDs=False):
        """Writes Model to file. The format is deduced from the filename
        extension given, and available options are '.gro', '.pdb', '.fasta',
        '.pir'.

        Parameters
        ----------
        fn : str
            filename
        title : str, optional
            title of the Model to write in fn
        nr : int, optional
            what is nr?
        bPDBTER : bool, optional
            whether to separate chains with TER entries. Only relevant when
            writing PDB files. Default is False.
        bAssignChainIDs : bool, optional
            whether to write the chains IDs to file. Only relevant when
            writing PDB files. Default is False.
        """
        ext = fn.split('.')[-1]
        if ext == 'pdb':
            self.writePDB(fn, title, nr, bPDBTER, bAssignChainIDs)
        elif ext == 'gro':
            self.writeGRO(fn, title)
        elif ext == 'pir':
            self.writePIR(fn, title)
        elif ext == 'fasta':
            self.writeFASTA(fn, title)
        else:
            print('pmx_Error_> Can only write pdb or gro!', file=sys.stderr)
            sys.exit(1)

    def make_chains(self):
        """Initialises the Chain instances from the input"""
        self.chains = []
        self.chdic = {}
        cur_chain = None
        ch = None
        for atom in self.atoms:
            if atom.chain_id == cur_chain:
                if ch:
                    atom.chain = ch
                    ch.atoms.append(atom)
                else:
                    ch = chain.Chain()
                    cur_chain = atom.chain_id
                    ch.id = atom.chain_id
                    atom.chain = ch
                    ch.atoms.append(atom)
            else:
                if ch:
                    self.chains.append(ch)
                    ch = chain.Chain()
                    cur_chain = atom.chain_id
                    ch.id = cur_chain
                    atom.chain = ch
                    ch.atoms.append(atom)
                else:
                    ch = chain.Chain()
                    cur_chain = atom.chain_id
                    ch.id = atom.chain_id
                    atom.chain = ch
                    ch.atoms.append(atom)

        self.chains.append(ch)
        for ch in self.chains:
            ch.model = self
            idx = ch.id
            self.chdic[idx] = ch

    def make_residues(self):
        """Initialises the Molecule instances from the input"""
        self.residues = []
        for ch in self.chains:
            cur_mol = None
            mol = None
            for atom in ch.atoms:
                if atom.resnr == cur_mol:
                    if mol:
                        mol.atoms.append(atom)
                    else:
                        mol = Molecule()
                        cur_mol = atom.resnr
                        mol.resname = atom.resname
                        mol.id = cur_mol
                        mol.atoms.append(atom)
                else:
                    if mol:
                        mol.model = self
                        mol.chain = ch
                        ch.residues.append(mol)
                        self.residues.append(mol)
                        mol = Molecule()
                        cur_mol = atom.resnr
                        mol.resname = atom.resname
                        mol.id = cur_mol
                        mol.atoms.append(atom)
                    else:
                        mol = Molecule()
                        cur_mol = atom.resnr
                        mol.resname = atom.resname
                        mol.id = cur_mol
                        mol.atoms.append(atom)
            self.residues.append(mol)
            ch.residues.append(mol)
        for r in self.residues:
            r.assign_moltype()
            for atom in r.atoms:
                atom.molecule = r
                atom.model = self
                r.model = self
        for ch in self.chains:
            for r in ch.residues:
                r.chain = ch
                r.chain_id = ch.id

    def __readPDB(self, fname=None, pdbline=None):
        """Reads a PDB file"""
        if pdbline:
            lines = pdbline.split('\n')
        else:
            lines = open(fname, 'r').readlines()
        for line in lines:
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                a = Atom().readPDBString(line)
                self.atoms.append(a)
            if line[:6] == 'CRYST1':
                self.box = _p.box_from_cryst1(line)
        self.make_chains()
        self.make_residues()
        self.unity = 'A'
        return self

    # TODO: make readPDB and readPDBTER a single function. It seems like
    # readPDBTER is more general PDB reader?
    def __readPDBTER(self, fname=None, pdbline=None, bNoNewID=True):
        """Reads a PDB file with more options than __readPDB ?"""
        if pdbline:
            lines = pdbline.split('\n')
        else:
            lines = open(fname, 'r').readlines()

        chainIDstring = ('ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                         'abcdefghijklmnoprstuvwxyz'
                         '123456789')
        bNewChain = True
        chainID = ' '
        prevID = ' '
        prevResID = 0
        usedChainIDs = ''
        atomcount = 1

        for line in lines:
            if 'TER' in line:
                bNewChain = True
            if (line[:4] == 'ATOM') or (line[:6] == 'HETATM'):
                a = Atom().readPDBString(line, origID=atomcount)
                atomcount += 1
                # identify chain change by ID (when no TER is there)
                if (a.chain_id != prevID):
                    bNewChain = True
                if (a.resnr != prevResID):
                    try:
                        if a.resnr != prevResID+1:
                            bNewChain = True
                    except TypeError:
                        bNewChain = False
                prevID = a.chain_id
                prevResID = a.resnr
                if bNewChain is True:
                    if (a.chain_id == ' ') or (a.chain_id == chainID):
                        # find a new chain id
                        bFound = False
                        while bFound is False:
                            foo = chainIDstring[0]
                            chainIDstring = chainIDstring.lstrip(chainIDstring[0])
                            if foo not in usedChainIDs:
                                bFound = True
                                usedChainIDs = usedChainIDs+foo
                                chainID = foo
                                if bNoNewID is True:
                                    chainID = "pmx"+foo
                    else:
                        chainID = a.chain_id
                        usedChainIDs = usedChainIDs + chainID
                a.chain_id = chainID
                self.atoms.append(a)
                bNewChain = False
            if line[:6] == 'CRYST1':
                self.box = _p.box_from_cryst1(line)
        self.make_chains()
        self.make_residues()
        self.unity = 'A'
        return self

    def __readGRO(self, filename):
        """Reads a GRO file"""
        l = open(filename).readlines()
        # first line is name/comment
        name = l[0].rstrip()
        self.title = name
        # next line is number of atoms
        natoms = int(l[1])
        atoms_parsed = 0
        while atoms_parsed != natoms:
            line = l[atoms_parsed+2]
            resid = int(line[:5])
            resname = line[5:9].strip()
            name = line[10:15].strip()
            idx = int(line[15:20])
            rest = line[20:].split()
            assert len(rest) in [3, 6]
            x = float(rest[0])
            y = float(rest[1])
            z = float(rest[2])
            coords = [x, y, z]
            if len(rest) == 6:
                vx = float(rest[3])
                vy = float(rest[4])
                vz = float(rest[5])
                vel = [vx, vy, vz]
            else:
                vel = [0, 0, 0]
            a = Atom(id=idx, name=name, resname=resname,
                     resnr=resid, x=coords, v=vel, unity='nm')
            a.get_symbol()
            self.atoms.append(a)
            atoms_parsed += 1
        box_line = [float(i) for i in l[-1].split()]
        assert len(box_line) in [3, 9]
        box = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        box[0][0] = box_line[0]
        box[1][1] = box_line[1]
        box[2][2] = box_line[2]
        if len(box_line) == 3:
            box[0][1] = 0
            box[0][2] = 0
            box[1][0] = 0
            box[1][2] = 0
            box[2][0] = 0
            box[2][1] = 0
        else:
            box[0][1] = box_line[3]
            box[0][2] = box_line[4]
            box[1][0] = box_line[5]
            box[1][2] = box_line[6]
            box[2][0] = box_line[7]
            box[2][1] = box_line[8]
        self.box = box
        self.make_chains()
        self.make_residues()
        self.unity = 'nm'
        return self

    def assign_moltype(self):
        """Identifies what type of molecule the Model is: protein, dna, or rna.
        If it is a mix, or if it is an organic molecule, "unknown" is
        assigned to self.moltype.
        """
        residues = set([r.resname for r in self.residues])

        # do not consider water and ions
        residues -= library._water
        residues -= library._ions

        # determine type
        if residues.issubset(library._protein_residues):
            self.moltype = 'protein'
        elif residues.issubset(library._dna_residues):
            self.moltype = 'dna'
        elif residues.issubset(library._rna_residues):
            self.moltype = 'rna'
        else:
            self.moltype = 'unknown'

    def read(self, filename, bPDBTER=False, bNoNewID=True):
        """PDB/GRO file reader.

        Parameters
        ----------
        filename : str
            name of input file
        bPDBTER : bool, optional
            whether the file contains TER records?. Default is False.
        bNoNewID : bool, optional
            whether to assign new chain IDs. If True, new chain IDs starting
            with 'pmx' will be assigned(?). Only relevant if bPDBTER is also
            True. Default is True.
        """
        ext = filename.split('.')[-1]
        if ext == 'pdb':
            if bPDBTER is True:
                return self.__readPDBTER(fname=filename,
                                         pdbline=None,
                                         bNoNewID=bNoNewID)
            else:
                return self.__readPDB(fname=filename)
        elif ext == 'gro':
            return self.__readGRO(filename)
        else:
            raise IOError('ERROR: Can only read pdb or gro!')

    def renumber_residues(self):
        """Renumbers all residues from 1."""
        for i, res in enumerate(self.residues):
            res.set_orig_resid(res.id)
            res.set_resid(i+1)

    # TODO/FIXME: should add/remove/append atoms all be only once in
    # Atomselection? At the moment they are repeated in Molecule, Chain,
    # and Model
    def remove_atom(self, atom):
        """Removes an Atom instance.

        Parameters
        ----------
        atom : Atom
            Atom instance to remove
        """
        m = atom.molecule
        m.remove_atom(atom)

    def remove_residue(self, residue):
        """Removes a Molecule/residue instance.

        Parameters
        ----------
        residue : Molecule
            Molecule instance to remove
        """
        ch = residue.chain
        ch.remove_residue(residue)

    def remove_chain(self, key):
        """Removes a Chain instance given the chain ID.

        Parameters
        ----------
        key : str
            ID of the chain to remove
        """
        if key not in self.chdic:
            print('No chain %s to remove....' % key)
            print('No changes applied.')
            return
        for ch in self.chains:
            if ch.id == key:
                idx = self.chains.index(ch)
                del self.chains[idx]
                del self.chdic[key]
        self.resl_from_chains()
        self.al_from_resl()
        self.renumber_residues()
        self.renumber_atoms()

    def __delitem__(self, key):
        self.remove_chain(key)

    def insert_residue(self, pos, res, chain_id):
        """Inserts a residue in Model.

        Parameters
        ----------
        pos : int
            position/index where to insert the residue
        res : Molecule
            Molecule instance to insert containig the residue of interest
        chain_id : str
            ID of the chain where to insert the residue
        """
        ch = self.chdic[chain_id]
        ch.insert_residue(pos, res)

    def replace_residue(self, residue, new, bKeepResNum=False):
        """Replaces a residue.

        Parameters
        ----------
        residue : Molecule
            residue to replace
        new : Molecule
            residue to insert
        bKeepResNum : bool, optional
            whether to keep residue ID of the residue that is inserted.
            Default is False
        """
        ch = residue.chain
        ch.replace_residue(residue, new, bKeepResNum)

    def insert_chain(self, pos, new_chain):
        """Inserts a Chain in Model.

        Parameters
        ----------
        pos : int
            index where to insert the chain within the list of chains
        new_chain : Chain
            instance of the Chain to insert in Model
        """
        if new_chain.id in self.chdic:
            print('Chain identifier %s already in use!' % new_chain.id)
            print('Changing chain identifier to 0')
            new_chain.set_chain_id('0')
        self.chains.insert(pos, new_chain)
        self.resl_from_chains()
        self.al_from_resl()
        self.make_chains()
        self.make_residues()
        self.renumber_atoms()
        self.renumber_residues()

    def append(self, new_chain):
        """ we assume chain is a Chain"""
        idx = len(self.chains)
        self.insert_chain(idx, new_chain)

    def fetch_residue(self, idx, chain=None):
        """Get a residue based on its index, or index and chain.

        Parameters
        ----------
        idx : int
            ID of the residue to fetch.
        chain : str, optional
            chain ID of the residue. This is needed when you have multiple
            chains and you have not renumbered the residues.

        Returns
        -------
        residue : Molecule
            Molecule instance of the residue found.
        """

        # check idx is a valid selection
        if idx not in [r.id for r in self.residues]:
            raise ValueError('resid %s not found in Model residues' % idx)

        if chain is None:
            # check selection is unique
            if [r.id for r in self.residues].count(idx) != 1:
                raise ValueError('idx choice %s results in non-unique selection' % idx)
        else:
            # check chain is a valid selection
            if chain not in [c.id for c in self.chains]:
                raise ValueError('chain ID "%s" not found in Model chains' % chain)
            # check idx+chain is a valid selection
            if idx not in [r.id for r in self.chdic[chain].residues]:
                raise ValueError('resid %s not found in chain "%s"' % (idx, chain))
            # check selection is unique
            if [r.id for r in self.chdic[chain].residues].count(idx) != 1:
                raise ValueError('idx choice %s for chain "%s" results in non-unique selection' % (idx, chain))

        # then find and return the residue
        if chain is None:
            for r in self.residues:
                if r.id == idx:
                    return r
        elif chain is not None:
            for r in self.chdic[chain].residues:
                if r.id == idx:
                    return r

    def fetch_residues(self, key, inv=False):
        """Gets residues using a list of residue names.

        Parameters
        ----------
        key : str or list
            resname
        inv : bool
            invert selection. Default is False. If True, it finds all residues
            with resnames different from those in key.

        Returns
        -------
        residues : list
            list of Molecule instances with the residues
        """
        if not hasattr(key, "append"):
            key = [key]
        result = []
        if not inv:
            for r in self.residues:
                if r.resname in key:
                    result.append(r)
        else:
            for r in self.residues:
                if r.resname not in key:
                    result.append(r)
        return result

    def al_from_resl(self):
        self.atoms = []
        for r in self.residues:
            for atom in r.atoms:
                self.atoms.append(atom)

    def resl_from_chains(self):
        self.residues = []
        for ch in self.chains:
            for r in ch.residues:
                self.residues.append(r)

    def copy(self):
        return copy.deepcopy(self)

    def get_mol2_types(self):
        if self.atoms[0].symbol == '':
            self.get_symbol()
        for ch in self.chains:
            ch.get_mol2_types()

    def get_mol2_resname(self):
        for ch in self.chains:
            ch.get_mol2_resname()

    def get_nterms(self):
        nter = []
        for ch in self.chains:
            first = ch.residues[0]      # first residue
            if first.resname in library._one_letter.keys():
                nter.append(first)
        return nter

    def get_cterms(self):
        cter = []
        for ch in self.chains:
            last = ch.residues[-1]      # last residue
            if last.resname in library._one_letter.keys():
                cter.append(last)
        return last

    def rename_atoms(self):
        for c in self.chains:
            c.rename_atoms()

    def residue(self, idx):
        return self.residues[idx-1]

    def chain(self, iden):
        return self.chdic[iden]
