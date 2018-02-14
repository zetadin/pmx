from __future__ import print_function, division
from pmx import __version__
import sys
import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt
from copy import deepcopy
from scipy.special import erf
from utils import data2gauss, gauss_func, natural_sort
from estimators import Jarz, Crooks, BAR
from parser import read_and_format
import os
import time
import pickle
import warnings


__all__ = ['read_dgdl_files', 'integrate_dgdl',
           'ks_norm_test', 'make_cgi_plot', 'analyse_dgdl']


# ==============================================================================
# higher level function
# ==============================================================================
def analyse_dgdl(filesAB, filesBA, outfn=None, T=298.15, prec=4,
                 skip=1, sliceit=None, rand=None, index=None,
                 methods=['BAR'], nboots=0, nblocks=1, reverseB=False,
                 iA=None, iB=None, oA=None, oB=None,
                 do_ks_test=False, integ_only=False, quiet=True, units='kj',
                 do_pickle=False, cgi_plot=False, dpi=300, nbins=10):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments
    """

    # Constants
    kb = 0.00831447215   # kJ/(K*mol)

    # sort input
    filesAB = natural_sort(filesAB)
    filesBA = natural_sort(filesBA)

    # setup output
    if outfn is not None:
        out = open(outfn, 'w')
    else:
        out = None

    # start timing
    stime = time.time()

    # -------------------
    # Select output units
    # -------------------
    if units.lower() == 'kj':
        # kJ is the input from GMX
        unit_fact = 1.
        units = 'kJ/mol'
    elif units == 'kcal':
        unit_fact = 1./4.184
        units = 'kcal/mol'
    elif units.lower() == 'kt':
        unit_fact = 1./(kb*T)
        units = 'kT'
    else:
        exit('No unit type \'%s\' available' % units)

    print("# analyze_dgdl.py, pmx version = %s" % __version__, file=out)
    print("# pwd = %s" % os.getcwd(), file=out)
    print("# %s (%s)" % (time.asctime(), os.environ.get('USER')), file=out)
    print("# command = %s" % ' '.join(sys.argv), file=out)
    _tee(out, "\n", quiet=quiet)

    # ==========
    # Parse Data
    # ==========

    # If list of dgdl.xvg files are provided, parse dgdl
    if iA is None and iB is None:
        # If random selection is chosen, do this before reading files and
        # calculating the work values.
        if rand is not None:
            filesAB = np.random.choice(filesAB, size=rand, replace=False)
            filesBA = np.random.choice(filesBA, size=rand, replace=False)
            _tee(out, 'Selected random subset of %d trajectories.' % rand,
                 quiet=quiet)

        # If slice values provided, select the files needed. Again before
        # reading files so speed up the process
        if sliceit is not None:
            first = sliceit[0]
            last = sliceit[1]
            _tee(out, ' First trajectories read: %s and %s'
                 % (filesAB[first], filesBA[first]), quiet=quiet)
            _tee(out, ' Last trajectories  read: %s and %s'
                 % (filesAB[last-1], filesBA[last-1]), quiet=quiet)
            _tee(out, '', quiet=quiet)
            filesAB = filesAB[first:last]
            filesBA = filesBA[first:last]

        # If index values provided, select the files needed
        if index is not None:
            # Avoid index out of range error if "wrong" indices are provided
            filesAB = [filesAB[i] for i in index if i < len(filesAB)]
            filesBA = [filesBA[i] for i in index if i < len(filesBA)]
            # ...but warn if this happens
            if any(i > (len(filesAB) - 1) for i in index):
                warnings.warn('\nindex out of range for some of your chosen '
                              '\nindices for the forward work values. This means you are'
                              '\ntrying to select input files that are not present.')
            if any(i > (len(filesBA) - 1) for i in index):
                warnings.warn('\nindex out of range for some of your chosen'
                              '\nindices for the reverse work values. This means you are'
                              '\ntrying to select input files that are not present.')

        # when skipping start count from end: in this way the last frame is
        # always included, and what can change is the first one
        filesAB = list(reversed(filesAB[::-skip]))
        filesBA = list(reversed(filesBA[::-skip]))

        # --------------------
        # Now read in the data
        # --------------------
        if quiet is True:
            print('Processing the input data...')
        elif quiet is False:
            print(' ========================================================')
            print('                   PROCESSING THE DATA')
            print(' ========================================================')
            print('  Forward Data')
        res_ab = read_dgdl_files(filesAB, lambda0=0,
                                 invert_values=False, verbose=not quiet)
        if quiet is False:
            print('  Reverse Data')
        res_ba = read_dgdl_files(filesBA, lambda0=1,
                                 invert_values=reverseB, verbose=not quiet)

        if oA is not None:
            _dump_integ_file(oA, filesAB, res_ab)
        if oB is not None:
            _dump_integ_file(oB, filesBA, res_ba)

    # If work values are given as input instead, read those
    elif iA is not None and iB is not None:
        res_ab = []
        res_ba = []
        for fn in iA:
            print('\t\tReading integrated values (A->B) from ', fn)
            res_ab.extend(_data_from_file(fn))
        for fn in iB:
            print('\t\tReading integrated values (B->A) from ', fn)
            res_ba.extend(_data_from_file(fn))
    else:
        raise ValueError('you need to provide either none of both sets of '
                         'integrated work values.')

    # If asked to only do the integration of dhdl.xvg, exit
    if integ_only and quiet is False:
        print('\n    Integration done. Skipping analysis.')
        print('\n    ......done........\n')
        sys.exit(0)

    # ==============
    # Begin Analysis
    # ==============
    _tee(out, ' ========================================================', quiet=quiet)
    _tee(out, '                       ANALYSIS', quiet=quiet)
    _tee(out, ' ========================================================', quiet=quiet)
    _tee(out, '  Number of forward (0->1) trajectories: %d' % len(res_ab), quiet=quiet)
    _tee(out, '  Number of reverse (1->0) trajectories: %d' % len(res_ba), quiet=quiet)
    _tee(out, '  Temperature : %.2f K' % T, quiet=quiet)

    # ============================
    # Crooks Gaussian Intersection
    # ============================
    if 'cgi' in methods:
        _tee(out, '\n --------------------------------------------------------', quiet=quiet)
        _tee(out, '             Crooks Gaussian Intersection     ', quiet=quiet)
        _tee(out, ' --------------------------------------------------------', quiet=quiet)

        if quiet is True:
            print('Running CGI analysis...')
        elif quiet is False:
            print('  Calculating Intersection...')

        cgi = Crooks(wf=res_ab, wr=res_ba, nboots=nboots, nblocks=nblocks)
        if do_pickle is True:
            pickle.dump(cgi, open("cgi_results.pkl", "wb"))

        _tee(out, '  CGI: Forward Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mf*unit_fact,
                                                s=cgi.devf*unit_fact,
                                                p=prec, u=units),
             quiet=quiet)
        _tee(out, '  CGI: Reverse Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mr*unit_fact,
                                                s=cgi.devr*unit_fact,
                                                p=prec, u=units),
             quiet=quiet)

        if cgi.inters_bool is False:
            _tee(out, '\n  Gaussians too close for intersection calculation', quiet=quiet)
            _tee(out, '   --> Taking difference of mean values', quiet=quiet)

        _tee(out, '  CGI: dG = {dg:8.{p}f} {u}'.format(dg=cgi.dg*unit_fact,
                                                       p=prec, u=units), quiet=quiet)
        _tee(out, '  CGI: Std Err (bootstrap:parametric) = {e:8.{p}f} {u}'.format(e=cgi.err_boot1*unit_fact,
                                                                                  p=prec, u=units), quiet=quiet)

        if nboots > 0:
            _tee(out, '  CGI: Std Err (bootstrap) = {e:8.{p}f} {u}'.format(e=cgi.err_boot2*unit_fact,
                                                                           p=prec, u=units), quiet=quiet)

        if nblocks > 1:
            _tee(out, '  CGI: Std Err (blocks) = {e:8.{p}f} {u}'.format(e=cgi.err_blocks*unit_fact,
                                                                        p=prec, u=units), quiet=quiet)

    # --------------
    # Normality test
    # --------------
    if do_ks_test:
        if quiet is False:
            print('\n  Running KS-test...')
        q0, lam00, check0, bOk0 = ks_norm_test(res_ab)
        q1, lam01, check1, bOk1 = ks_norm_test(res_ba)

        _tee(out, '    Forward: gaussian quality = %3.2f' % q0, quiet=quiet)
        if bOk0:
            _tee(out, '             ---> KS-Test Ok', quiet=quiet)
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q0, check0), quiet=quiet)
        _tee(out, '    Reverse: gaussian quality = %3.2f' % q1, quiet=quiet)
        if bOk1:
            _tee(out, '             ---> KS-Test Ok', quiet=quiet)
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q1, check1), quiet=quiet)

    # ========================
    # Bennett Acceptance Ratio
    # ========================
    if 'bar' in methods:
        _tee(out, '\n --------------------------------------------------------', quiet=quiet)
        _tee(out, '             Bennett Acceptance Ratio     ', quiet=quiet)
        _tee(out, ' --------------------------------------------------------', quiet=quiet)

        if quiet is True:
            print('Running BAR analysis...')
        elif quiet is False:
            print('  Running Nelder-Mead Simplex algorithm... ')

        bar = BAR(res_ab, res_ba, T=T, nboots=nboots, nblocks=nblocks)
        if do_pickle:
            pickle.dump(bar, open("bar_results.pkl", "wb"))

        _tee(out, '  BAR: dG = {dg:8.{p}f} {u}'.format(dg=bar.dg*unit_fact, p=prec, u=units), quiet=quiet)
        _tee(out, '  BAR: Std Err (analytical) = {e:8.{p}f} {u}'.format(e=bar.err*unit_fact, p=prec, u=units), quiet=quiet)

        if nboots > 0:
            _tee(out, '  BAR: Std Err (bootstrap)  = {e:8.{p}f} {u}'.format(e=bar.err_boot*unit_fact, p=prec, u=units), quiet=quiet)
        if nblocks > 1:
            _tee(out, '  BAR: Std Err (blocks)  = {e:8.{p}f} {u}'.format(e=bar.err_blocks*unit_fact, p=prec, u=units), quiet=quiet)

        _tee(out, '  BAR: Conv = %8.2f' % bar.conv, quiet=quiet)

        if nboots > 0:
            _tee(out, '  BAR: Conv Std Err (bootstrap) = %8.2f' % bar.conv_err_boot, quiet=quiet)

    # =========
    # Jarzynski
    # =========
    if 'jarz' in methods:
        _tee(out, '\n --------------------------------------------------------', quiet=quiet)
        _tee(out, '             Jarzynski estimator     ', quiet=quiet)
        _tee(out, ' --------------------------------------------------------', quiet=quiet)

        if quiet is True:
            print('Running Jarz analysis...')
        jarz = Jarz(wf=res_ab, wr=res_ba, T=T, nboots=nboots, nblocks=nblocks, quiet=quiet)
        if do_pickle:
            pickle.dump(jarz, open("jarz_results.pkl", "wb"))

        _tee(out, '  JARZ: dG Forward = {dg:8.{p}f} {u}'.format(dg=jarz.dg_for*unit_fact,
                                                                p=prec, u=units), quiet=quiet)
        _tee(out, '  JARZ: dG Reverse = {dg:8.{p}f} {u}'.format(dg=jarz.dg_rev*unit_fact,
                                                                p=prec, u=units), quiet=quiet)
        _tee(out, '  JARZ: dG Mean    = {dg:8.{p}f} {u}'.format(dg=jarz.dg_mean*unit_fact,
                                                                p=prec, u=units), quiet=quiet)
        if nboots > 0:
            _tee(out, '  JARZ: Std Err Forward (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_for*unit_fact,
                                                                                    p=prec, u=units), quiet=quiet)
            _tee(out, '  JARZ: Std Err Reverse (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_rev*unit_fact,
                                                                                    p=prec, u=units), quiet=quiet)

        if nblocks > 1:
            _tee(out, '  JARZ: Std Err Forward (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_for*unit_fact,
                                                                                 p=prec, u=units), quiet=quiet)
            _tee(out, '  JARZ: Std Err Reverse (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_rev*unit_fact,
                                                                                 p=prec, u=units), quiet=quiet)


    _tee(out, ' ========================================================', quiet=quiet)

    if 'cgi' in methods and cgi_plot is not None:
        if quiet is False:
            print('\n   Plotting histograms......')
        make_cgi_plot(cgi_plot, res_ab, res_ba, cgi.dg, cgi.err_boot1,
                      nbins, dpi)

    if quiet is True:
        print('Done')
    if quiet is False:
        print('\n   ......done...........\n')

    if do_pickle and quiet is False:
        print('   NOTE: units of results in pickled files are as in the\n'
              '   provided dgdl.xvg or integ.dat files. These are typically\n'
              '   in kJ/mol when using dgdl.xvg files from Gromacs.\n')
    # execution time
    etime = time.time()
    h, m, s = _time_stats(etime-stime)
    if quiet is False:
        print("   Execution time = %02d:%02d:%02d\n" % (h, m, s))


# ==============================================================================
# Other functions
# ==============================================================================
def read_dgdl_files(lst, lambda0=0, invert_values=False, verbose=True):
    '''Takes a list of dgdl.xvg files and returns the integrated work values.

    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    w : array
        array of work values.
    lst : list
        sorted list of input dgdl.avg files corresponding to the work values
        in w.
    '''

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    _check_dgdl(lst[0], lambda0, verbose=verbose)
    first_w, ndata = integrate_dgdl(lst[0], lambda0=lambda0,
                                    invert_values=invert_values)
    w_list = [first_w]
    for idx, f in enumerate(lst[1:]):
        if verbose is True:
            sys.stdout.write('\r    Reading %s' % f)
            sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0,
                              invert_values=invert_values)
        if w is not None:
            w_list.append(w)

    if verbose is True:
        print('\n')
    return w_list


def integrate_dgdl(fn, ndata=-1, lambda0=0, invert_values=False):
    '''Integrates the data in a dgdl.xvg file.

    Parameters
    ----------
    fn : str
        the inpur dgdl.xvg file from Gromacs.
    ndata : int, optional
        number of datapoints in file. If -1, then ??? default is -1.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    integr : float
        result of the integration performed using Simpson's rule.
    ndata : int
        number of data points in the input file.
    '''

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    lines = open(fn).readlines()
    if not lines:
        return None, None

    # extract dgdl datapoints into r
    # TODO: we removed the check for file integrity. We could have an
    # optional files integrity check before calling this integration func

    lines = [l for l in lines if l[0] not in '#@&']
    r = map(lambda x: float(x.split()[1]), lines)

    if ndata != -1 and len(r) != ndata:
        try:
            print(' !! Skipping %s ( read %d data points, should be %d )'
                  % (fn, len(r), ndata))
        except:
            print(' !! Skipping %s ' % (fn))
        return None, None
    # convert time to lambda
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    # arrays for the integration
    # --------------------------
    # array of lambda values
    x = [lambda0+i*dlambda for i, dgdl in enumerate(r)]
    # array of dgdl
    y = r

    if lambda0 == 1:
        x.reverse()
        y.reverse()

    if invert_values is True:
        integr = simps(y, x) * (-1)
        return integr, ndata
    else:
        integr = simps(y, x)
        return integr, ndata


def ks_norm_test(data, alpha=0.05, refks=None):
    '''Performs a Kolmogorov-Smirnov test of normality.

    Parameters
    ----------
    data : array_like
        a one-dimensional array of values. This is the distribution tested
        for normality.
    alpha : float
        significance level of the statistics. Default if 0.05.
    refks : ???
        ???

    Returns
    -------
    Q : float
    lam0 : float
    check : float
    bOk : bool
    '''

    def ksref():
        f = 1
        potent = 10000
        lamb = np.arange(0.25, 2.5, 0.001)
        q = np.zeros(len(lamb), float)
        res = []
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f = -f
        for i in range(len(lamb)):
            res.append((lamb[i], q[i]))
        return res

    def ksfunc(lamb):
        f = 1
        potent = 10000
        q = 0
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f *= -1
        return q

    def edf(dg_data):
        edf_ = []
        ndata = []
        data = deepcopy(dg_data)
        data.sort()
        N = float(len(data))
        cnt = 0
        for item in data:
            cnt += 1
            edf_.append(cnt/N)
            ndata.append(item)
        ndata = np.array(ndata)
        edf_ = np.array(edf_)
        return ndata, edf_

    def cdf(dg_data):
        data = deepcopy(dg_data)
        data.sort()
        mean = np.average(data)
        sig = np.std(data)
        cdf = 0.5*(1+erf((data-mean)/float(sig*np.sqrt(2))))
        return cdf

    N = len(data)
    nd, ed = edf(data)
    cd = cdf(data)
    siglev = 1-alpha
    dval = []
    for i, val in enumerate(ed):
        d = abs(val-cd[i])
        dval.append(d)
        if i:
            d = abs(ed[i-1]-cd[i])
            dval.append(d)
    dmax = max(dval)
    check = np.sqrt(N)*dmax
    if not refks:
        refks = ksref()
    lst = filter(lambda x: x[1] > siglev, refks)
    lam0 = lst[0][0]
    if check >= lam0:
        bOk = False
    else:
        bOk = True

    q = ksfunc(check)
    return (1-q), lam0, check, bOk


def make_cgi_plot(fname, data1, data2, result, err, nbins, dpi=300):
    '''Plots work distributions and results for Crooks Gaussian Intersection'''

    def smooth(x, window_len=11, window='hanning'):

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than "
                             "window size.")
        if window_len < 3:
            return x
        if window not in ['flat', 'hanning', 'hamming',
                          'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                             "'bartlett', 'blackman'")
        s = np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
        # moving average
        if window == 'flat':
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')
        y = np.convolve(w/w.sum(), s, mode='same')
        return y[window_len-1:-window_len+1]

    plt.figure(figsize=(8, 6))
    x1 = range(len(data1))
    x2 = range(len(data2))
    if x1 > x2:
        x = x1
    else:
        x = x2
    mf, devf, Af = data2gauss(data1)
    mb, devb, Ab = data2gauss(data2)

    maxi = max(data1+data2)
    mini = min(data1+data2)

    sm1 = smooth(np.array(data1))
    sm2 = smooth(np.array(data2))
    plt.subplot(1, 2, 1)
    plt.plot(x1, data1, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
    plt.plot(x1, sm1, 'g-', linewidth=3)
    plt.plot(x2, data2, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
    plt.plot(x2, sm2, 'b-', linewidth=3)
    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    plt.ylabel(r'W [kJ/mol]', fontsize=20)
    plt.xlabel(r'# Snapshot', fontsize=20)
    plt.grid(lw=2)
    plt.xlim(0, x[-1]+1)
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplot(1, 2, 2)
    plt.hist(data1, bins=nbins, orientation='horizontal', facecolor='green',
             alpha=.75, normed=True)
    plt.hist(data2, bins=nbins, orientation='horizontal', facecolor='blue',
             alpha=.75, normed=True)

    x = np.arange(mini, maxi, .5)

    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plt.plot(y1, x, 'g--', linewidth=2)
    plt.plot(y2, x, 'b--', linewidth=2)
    size = max([max(y1), max(y2)])
    res_x = [result, result]
    res_y = [0, size*1.2]
    plt.plot(res_y, res_x, 'k--', linewidth=2,
             label=r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)


def _check_dgdl(fn, lambda0, verbose=True):
    '''Prints some info about a dgdl.xvg file.'''
    l = open(fn).readlines()
    if not l:
        return None
    r = []
    for line in l:
        if line[0] not in '#@&':
            r.append([float(x) for x in line.split()])
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    if verbose is True:
        print('    # data points: %d' % ndata)
        print('    Length of trajectory: %8.3f ps' % r[-1][0])
        print('    Delta lambda: %8.5f' % dlambda)


def _dump_integ_file(outfn, f_lst, w_lst):
    with open(outfn, 'w') as f:
        for fn, w in zip(f_lst, w_lst):
            f.write('{dhdl} {work}\n'.format(dhdl=fn, work=w))


def _data_from_file(fn):
    data = read_and_format(fn, 'sf')
    return map(lambda a: a[1], data)


def _tee(fp, s, quiet=False):
    print(s, file=fp)
    if quiet is False:
        print(s)


def _time_stats(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return h, m, s
