from __future__ import print_function, division
from scipy.integrate import simps
import sys
from matplotlib import pyplot as plt
import numpy as np
from copy import deepcopy
from scipy.special import erf

__all__ = ['read_dgdl_files', 'integrate_dgdl',
           'ks_norm_test', 'make_cgi_plot']


def read_dgdl_files(lst, lambda0=0, invert_values=False):
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

    _check_dgdl(lst[0], lambda0)
    first_w, ndata = integrate_dgdl(lst[0], lambda0=lambda0,
                                    invert_values=invert_values)
    w_list = [first_w]
    for idx, f in enumerate(lst[1:]):
        sys.stdout.write('\r    Reading %s' % f)
        sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0,
                              invert_values=invert_values)
        if w is not None:
            w_list.append(w)

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


def data2gauss(data):
    '''Takes a one dimensional array and fits a Gaussian.

    Parameters
    ----------
    data : list
        1D array of values

    Returns
    -------
    float
        mean of the distribution.
    float
        standard deviation of the distribution.
    float
        height of the curve's peak.
    '''
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A


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

    y1 = _gauss_func(Af, mf, devf, x)
    y2 = _gauss_func(Ab, mb, devb, x)

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


def _check_dgdl(fn, lambda0):
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

    print('    # data points: %d' % ndata)
    print('    Length of trajectory: %8.3f ps' % r[-1][0])
    print('    Delta lambda: %8.5f' % dlambda)


def _gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y
