Analysis
--------

example on the use of analysis script and api


example::

    >>> from glob import glob
    >>> from pmx.analysis import read_dgdl_files
    >>> from pmx.estimators import BAR

    >>> # get dgdl.xvg files for the forward and reverse transitions
    >>> ff = glob('for*.xvg')  # forward
    >>> fr = glob('rev*.xvg')  # reverse

    >>> # now get integrated work values from them
    >>> wf = read_dgdl_files(ff, verbose=False)
    >>> wr = read_dgdl_files(fr, verbose=False)

    >>> # use BAR to estimate the free energy difference
    >>> bar = BAR(wf=wf, wr=wf, T=298)
    >>> print estimated dG and its uncertainty
    >>> print('%.2f +/- %.2f kJ/mol' % (bar.dg, bar.err))
    3.19 +/- 0.49 kJ/mol
