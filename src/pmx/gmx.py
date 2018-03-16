#!/usr/bin/env python

"""This module contains wrapper functions for some of the most used Gromacs
tools.
"""

from __future__ import absolute_import, division, print_function
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


def editconf(f, o='editconf.gro', bt='cubic', d=1.2, other_flags=''):
    """Simple ``gmx editconf`` wrapper.

    Parameters
    ----------
    f : str
        input structure file
    o : str, optional
        name of output structure file. Default is "solvate.gro"
    bt : str
        box type:triclinic, cubic, dodecahedron, or octahedron
    d : float
        distance between the solute and the box (nm)
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    gmx = get_gmx()
    call('{gmx} editconf -f {f} -o {o} -bt {bt} -d {d} '
         '{other_flags}'.format(gmx=gmx, f=f, o=o, bt=bt, d=d, other_flags=other_flags),
         shell=True)


def pdb2gmx(f, o='pdb2gmx.gro', p='topol.top', ff='amber99sb-star-ildn-mut',
            water='tip3p', other_flags=''):
    """Simple ``gmx pdb2gmx`` wrapper.

    Parameters
    ----------
    f : str
        input structure file
    o : str, optional
        output structure file. Default is "pdb2gmx.gro".
    p : str, optional
        name of topology file. Default is "topol.top".
    ff : str, optional
        forcefield. Default is "amber99sb-star-ildn-mut".
    water : str, optional
        water model. Default is "tip3p".
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    gmx = get_gmx()
    call('{gmx} pdb2gmx -f {f} -o {o} -p {p} -ff {ff} -water {water} '
         '{other_flags}'.format(gmx=gmx, f=f, o=o, p=p, ff=ff, water=water, other_flags=other_flags),
         shell=True)


def solvate(cp, cs='spc216.gro', p='topol.top', o='solvate.gro',
            other_flags=''):
    """Simple ``gmx solvate`` wrapper.

    Parameters
    ----------
    cp : str
        input structure file
    cs : str, optional
        structure file of the solvent. Default is "spc216.gro".
    p : str, optional
        name of topology file. Default is "topol.top"
    o : str, optional
        name of output structure file. Default is "solvate.gro"
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    gmx = get_gmx()
    call('{gmx} solvate -cp {cp} -cs {cs} -p {p} -o {o} '
         '{other_flags}'.format(gmx=gmx, cp=cp, cs=cs, p=p, o=o, other_flags=other_flags),
         shell=True)


def grompp(f, c, p, o='grompp.tpr', maxwarn=0, other_flags=''):
    """Simple ``gmx grompp`` wrapper.

    Parameters
    ----------
    f : str
        input mdp file
    c : str
        input structure file
    p : str
        input topology file. Default is "topol.top"
    o : str, optional
        output tpr file. Default is "grompp.tpr"
    maxwarn : int, optional
        number of allowed warnings. Default is 0.
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    gmx = get_gmx()
    call('{gmx} grompp -f {f} -c {c} -p {p} -o {o} -maxwarn {maxwarn}'
         '{other_flags}'.format(gmx=gmx, f=f, c=c, p=p, o=o, maxwarn=maxwarn, other_flags=other_flags),
         shell=True)


def genion():
    pass


def trjconv():
    pass
