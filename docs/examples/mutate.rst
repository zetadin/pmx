Mutate
------

example on the use of mutate script and api

You can download an example PDB file from here:  :download:`peptide.pdb <peptide.pdb>`

example::

    >>> from pmx.model import Model
    >>> from pmx.alchemy import mutate

    >>> # load the PDB file
    >>> m = Model('peptide.pdb')
    >>> # perform mutation
    >>> m2 = mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut')
    >>> # save mutant PDB file
    >>> m2.write('mutant.pdb')


the model can also be modified in place: ::

    >>> mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut', inplace=True)
    >>> m.write('mutant.pdb')
