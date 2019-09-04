import sys
import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt
from copy import deepcopy
from scipy.special import erf
from .utils import data2gauss, gauss_func

__all__ = ['read_dgdl_files', 'integrate_dgdl',
           'ks_norm_test', 'plot_work_dist']


def read_dgdl_files(lst, lambda0=0, lambda1=1, invert_values=False, verbose=True, sigmoid=0.0):
    '''Takes a list of dgdl.xvg files and returns the integrated work values.

    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    lambda0 : float
        lambda at which simulations started. Default is 0.
    lambda1 : float
        lambda at which simulations ended. Default is 1.
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

    # check lambda0 & 1 are between 0 and 1
    assert (lambda0>=0 and lambda0<=1), "Incorrect initial lambda "
    assert (lambda1>=0 and lambda1<=1), "Incorrect final lambda "

    #find the first good file
    good=False;
    idx=0;
    first_w=0;
    ndata=0;
    while(not good):
        try:
            _check_dgdl(lst[idx], lambda0=lambda0, lambda1=lambda1, verbose=verbose)
            first_w, ndata = integrate_dgdl(lst[idx], lambda0=lambda0, lambda1=lambda1,
                                        invert_values=invert_values, sigmoid=sigmoid)
            good=True
        except:
            print(' !! Error in checking %s' % (lst[idx]))
            good=False
            idx+=1


    w_list = [first_w]
    for idx, f in enumerate(lst[idx+1:]):
        if verbose is True:
            sys.stdout.write('\r    Reading %s' % f)
            sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0, lambda1=lambda1,
                              invert_values=invert_values, sigmoid=sigmoid)
        if w is not None:
            w_list.append(w)

    if verbose is True:
        print('\n')
    return w_list


def integrate_dgdl(fn, ndata=-1, lambda0=0, lambda1=1, invert_values=False, sigmoid=0.0):
    '''Integrates the data in a dgdl.xvg file.

    Parameters
    ----------
    fn : str
        the inpur dgdl.xvg file from Gromacs.
    ndata : int, optional
        number of datapoints in file. If -1, then ??? default is -1.
    lambda0 : float
        lambda at which simulations started. Default is 0.
    lambda1 : float
        lambda at which simulations ended. Default is 1.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    integr : float
        result of the integration performed using Simpson's rule.
    ndata : int
        number of data points in the input file.
    '''

    # check lambda0 is between 0 and 1
    assert (lambda0>=0 and lambda0<=1), "Incorrect initial lambda "
    assert (lambda1>=0 and lambda1<=1), "Incorrect final lambda "

    lines = open(fn, encoding="ISO-8859-1").readlines()
    if not lines:
        return None, None

    # extract dgdl datapoints into r
    # TODO: we removed the check for file integrity. We could have an
    # optional files integrity check before calling this integration func

    lines = [l for l in lines if l[0] not in '#@&']
    r=[]
    try:
        r = list(map(lambda x: float(x.split()[1]), lines))
    except:
        print(' !! Error in reading %s' % (fn))
        return None, None

    if ndata != -1 and len(r) != ndata:
        try:
            print(' !! Skipping %s ( read %d data points, should be %d )'
                  % (fn, len(r), ndata))
        except:
            print(' !! Skipping %s ' % (fn))
        return None, None
    # convert time to lambda
    ndata = len(r)
    dlambda = (lambda1-lambda0)/float(ndata)

    # arrays for the integration
    # --------------------------
    # array of lambda values
    x = [lambda0+i*dlambda for i, dgdl in enumerate(r)]
##### VG VG VG ####
    if sigmoid != 0.0:
        x = [1.0/(1.0+np.exp( -1.0*(lambda0+i*dlambda-0.5)*sigmoid ) ) for i, dgdl in enumerate(r)]
        x[0] = lambda0
        x[-1] = 1.0-lambda0
#    for i in x:
#        print(i)
##### VG VG VG ####
    # array of dgdl
    y = r

    if lambda0 > lambda1:
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
    lst = list(filter(lambda x: x[1] > siglev, refks))
    lam0 = lst[0][0]
    if check >= lam0:
        bOk = False
    else:
        bOk = True

    q = ksfunc(check)
    return (1-q), lam0, check, bOk


def plot_work_dist(wf, wr, fname='Wdist.png', nbins=20, dG=None, dGerr=None,
                   units='kJ/mol', dpi=300, window_len=11):
    '''Plots forward and reverse work distributions. Optionally, it adds the
    estimate of the free energy change and its uncertainty on the plot.

    Parameters
    ----------
    wf : list
        list of forward work values.
    wr : list
        list of reverse work values.
    fname : str, optional
        filename of the saved image. Default is 'Wdist.png'.
    nbins : int, optional
        number of bins to use for the histogram. Default is 20.
    dG : float, optional
        free energy estimate.
    dGerr : float, optional
        uncertainty of the free energy estimate.
    units : str, optional
        the units of dG and dGerr. Default is 'kJ/mol'.
    dpi : int
        resolution of the saved image file.
    window_len: int
        Width of smoothing window for plotting. Default is 11 samples.

    Returns
    -------
    None

    '''

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
    x1 = list(range(len(wf)))
    x2 = list(range(len(wr)))
    if x1 > x2:
        x = x1
    else:
        x = x2
    mf, devf, Af = data2gauss(wf)
    mb, devb, Ab = data2gauss(wr)

    maxi = max(wf+wr)
    mini = min(wf+wr)

    sm1 = smooth(np.array(wf), window_len=window_len)
    sm2 = smooth(np.array(wr), window_len=window_len)
    plt.subplot(1, 2, 1)
    plt.plot(x1, wf, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
    plt.plot(x1, sm1, 'g-', linewidth=3)
    plt.plot(x2, wr, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
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
    plt.hist(wf, bins=nbins, orientation='horizontal', facecolor='green',
             alpha=.75, density=True)
    plt.hist(wr, bins=nbins, orientation='horizontal', facecolor='blue',
             alpha=.75, density=True)

    x = np.arange(mini, maxi, .5)

    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plt.plot(y1, x, 'g--', linewidth=2)
    plt.plot(y2, x, 'b--', linewidth=2)
    size = max([max(y1), max(y2)])
    res_x = [dG, dG]
    res_y = [0, size*1.2]
    if dG is not None and dGerr is not None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f $\pm$ %.2f %s' % (dG, dGerr, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    elif dG is not None and dGerr is None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f %s' % (dG, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    else:
        plt.plot(res_y, res_x, 'k--', linewidth=2)

    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)


def _check_dgdl(fn, lambda0=0, lambda1=1, verbose=True):
    '''Prints some info about a dgdl.xvg file.'''
    lines = open(fn, encoding="ISO-8859-1").readlines()
    if not lines:
        return None
    r = []
    for line in lines:
        if line[0] not in '#@&':
            r.append([float(x) for x in line.split()])
    ndata = len(r)
    dlambda = (lambda1-lambda0)/float(ndata)

    if verbose is True:
        print('    # data points: %d' % ndata)
        print('    Length of trajectory: %8.3f ps' % r[-1][0])
        print('    Delta lambda: %8.5f' % dlambda)
