#!/usr/bin/env python

"""This module contains wrapper functions for some of the most used Gromacs
tools.
"""

from .utils import which
from subprocess import call


def get_gmx():
    """Gets path to gmx executable, and throws error if not found.
    """
    gmx = which('gmx')
    if gmx is not None:
        return gmx
    else:
        raise EnvironmentError('gmx executable not found')


def pdb2gmx(f='input.pdb', o='pdb2gmx.gro', p='topol.top', ff='amber99sb-mut',
            water='tip3p', other_flags=''):
    """Simple pdb2gmx wrapper.

    Parameters
    ----------
    f : str
        input structure
    o : str
        output structure file
    p : str
        name of topology file
    ff : str
        forcefield
    water : str
        water model
    other_flags : str
        additional flags as a string. These will be given to the
        ``gmx pdb2gmx`` command

    Returns
    -------
    None

    """
    gmx = get_gmx()

    call('{gmx} pdb2gmx -f {f} -o {o} -p {p} -ff {ff} -water {water} '
         '{other_flags}'.format(gmx=gmx, f=f, o=o, p=p, ff=ff, water=water, other_flags=other_flags),
         shell=True)


def solvate():
    pass


def genion():
    pass


def grompp():
    pass


def trjconv():
    pass
