Mutate
------

example on the use of mutate script and api

example::

    >>> from pmx.model import Model
    >>> from pmx.alchemy import mutate

    >>> # load the PDB file
    >>> m = Model('peptide.pdb')
    >>> # perform mutation
    >>> mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut')
    >>> # save mutant PDB file
    >>> m.write('mutant.pdb')
