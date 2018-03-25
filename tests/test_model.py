#!/usr/bin/env python
import pytest
from pmx.model import Model


def test_title(empty_model):
    assert empty_model.title == 'PMX MODEL'


@pytest.mark.parametrize("inputpdb, expected_type", [
    ('peptide.pdb', 'protein'),
    ('dna.pdb', 'dna'),
    ('rna.pdb', 'rna')
])
def test_model_moltype(gf, inputpdb, expected_type):
    m = Model(gf(inputpdb))
    assert m.moltype == expected_type


def test_protein_residue_names(protein_model):
    residues = [r.resname for r in protein_model.residues]
    true_residues = ['NASN', 'LEU', 'TYR', 'ILE', 'GLN', 'TRP', 'LEU', 'LYS',
                     'ASP', 'GLY', 'GLY', 'PRO',  'SER', 'SER', 'GLY', 'ARG',
                     'PRO', 'PRO', 'PRO', 'CSER']
    assert residues == true_residues
