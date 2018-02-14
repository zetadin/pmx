"""Various utility functions
"""

import os
from glob import glob
import types
import numpy as np
from re import split as resplit


# =============
# File IO utils
# =============
def ffopen(filename, mode='r', backup=True):

    if mode == 'w':
        if os.path.isfile(filename):
            if backup:
                print 'Backing up %s to %s~' % (filename, filename)
                os.rename(filename, filename+'~')
        try:
            fp = open(filename, 'w')
            return fp
        except:
            print 'Error: Could not open file %s' % filename

    elif mode == 'r':
        try:
            fp = open(filename, 'r')
            return fp
        except:
            print 'No such file %s' % filename

    else:
        return open(filename, mode)


def listFiles(dir='./', ext=None, abs=True, backups=False):

    """ returns a list of files in
    directory dir, optionally only
    certain file types"""

    if dir[-1] != os.sep:
        dir += os.sep
    if dir[0] == '~':
        home = os.environ.get('HOME')
        dir = os.path.join(home, dir[2:])
    dir = os.path.abspath(dir)+os.sep

    l = os.listdir(dir)

    fl = []
    if not ext:
        for f in l:
            if os.path.isfile(dir+f):
                if backups:
                    if f[0] == '#' or \
                       f[-1] == '~':
                        fl.append(dir+f)
                else:
                    fl.append(dir+f)

    elif type(ext) in [types.ListType, types.TupleType]:
        for ex in ext:
            if backups:
                ff = glob(dir+'#*'+ex+'*')
                ff += glob(dir+'*.'+ex+'~')
            else:
                ff = glob(dir+'*.'+ex)
            fl.extend(ff)

    elif type(ext) == types.StringType:
        if backups:
            print dir+'*.'+ext+'~'
            fl = glob(dir+'#*.'+ext+'*')
            fl += glob(dir+'*.'+ext+'~')
        else:
            fl.extend(glob(dir+'*.'+ext))

    if not abs:
        new = []
        for f in fl:
            new.append(f.split('/')[-1])
        return new

    return fl


def listDirs(dir='./'):
    """ returns a list of directrories in
    directory dir"""

    if dir[-1] != os.sep:
        dir += os.sep
    if dir[0] == '~':
        home = os.environ.get('HOME')
        dir = os.path.join(home, dir[2:])
    dir = os.path.abspath(dir)+os.sep

    l = os.listdir(dir)

    dl = []

    for f in l:
        if os.path.isdir(dir+f):
            dl.append(dir+f)

    return dl


def killBackups(arg, dirname, fname):

    l = listFiles(dirname, arg[0], arg[1], arg[2])
    if arg[3]:
        for f in l:
            print '%s' % f


def removeBackups(dir, check=True):

    if dir[-1] != os.sep:
        dir += os.sep
    if dir[0] == '~':
        home = os.environ.get('HOME')
        dir = os.path.join(home, dir[2:])
    dir = os.path.abspath(dir)+os.sep

    os.path.walk(dir, killBackups, (False, True, True, check))


# ======
# Others
# ======
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


def gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in resplit('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
