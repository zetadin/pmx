#!/usr/bin/env python

from pmx.analysis import analyse_dgdl
import argparse

# =======
# Options
# =======
def parse_options():

    parser = argparse.ArgumentParser(description='Calculates free energies '
            'from fast growth thermodynamic integration simulations. '
            'Available methods for free energy estimation: '
            'Crooks Gaussian Intersection (CGI); '
            'Benett Acceptance Ratio (BAR); '
            'Jarzinski equality (JARZ).')

    exclus = parser.add_mutually_exclusive_group()

    parser.add_argument('-fA',
                        metavar='dgdl',
                        dest='filesAB',
                        type=str,
                        help='dgdl.xvg files for the A->B simulations. Use '
                        'wildcard to select multiple xvg files: e.g. "-fa '
                        './forward_results/dgdl*.xvg"',
                        required=True,
                        nargs='+')
    parser.add_argument('-fB',
                        metavar='dgdl',
                        dest='filesBA',
                        type=str,
                        help='dgdl.xvg files for the B->A simulations Use '
                        'wildcard to select multiple xvg files: e.g. "-fb '
                        './backward_results/dgdl*.xvg"',
                        required=True,
                        nargs='+')
    parser.add_argument('-m',
                        metavar='method',
                        type=str.lower,
                        dest='methods',
                        help='Choose one or more estimators to use from the '
                        'available ones: CGI, BAR, JARZ. Default is all.',
                        default=['cgi', 'bar', 'jarz'],
                        nargs='+')
    parser.add_argument('-t',
                        metavar='temperature',
                        dest='T',
                        type=float,
                        help='Temperature in Kelvin. Default is 298.15.',
                        default=298.15)
    parser.add_argument('-o',
                        metavar='result file',
                        dest='outfn',
                        type=str,
                        help='Filename of output result file. Default is '
                        '"results.txt."',
                        default='results.txt')
    parser.add_argument('-b',
                        metavar='nboots',
                        dest='nboots',
                        type=int,
                        help='Number of bootstrap samples to use for the '
                        'bootstrap estimate of the standard errors. Default '
                        'is 0 (no bootstrap).',
                        default=0)
    parser.add_argument('-n',
                        metavar='nblocks',
                        dest='nblocks',
                        type=int,
                        help='Number of blocks to divide the data into for '
                        'an estimate of the standard error. You can use this '
                        'when multiple independent equilibrium simulations'
                        'have been run so to estimate the error from the '
                        'repeats. Default is 1 (i.e. no repeats). It assumes '
                        'the dgdl files for each repeat are read in order and '
                        'are contiguous, e.g. dgdl_0 to dgdl_9 is the first '
                        'repeat, dgdl_10 to dgdl_19 is the second one, etc.',
                        default=1)
    parser.add_argument('--integ_only',
                        dest='integ_only',
                        help='Whether to do integration only; the integrated '
                        'values are computed and saved, and the program '
                        'terminated. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('-iA',
                        metavar='work input',
                        dest='iA',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation.')
    parser.add_argument('-iB',
                        metavar='work input',
                        dest='iB',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation.')
    parser.add_argument('-oA',
                        metavar='work output',
                        dest='oA',
                        type=str,
                        help='File where to save the list of input dgdl'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation. Default is '
                        '"integA.dat"',
                        default='integA.dat')
    parser.add_argument('-oB',
                        metavar='work output',
                        dest='oB',
                        type=str,
                        help='File where to save the list of input dgdl'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation. Default is '
                        '"integB.dat"',
                        default='integB.dat')
    parser.add_argument('--reverseB',
                        dest='reverseB',
                        help='Whether to reverse the work values for the '
                        'backward (B->A) transformation. This is useful '
                        'when in Gromacs both forward and reverse simulations '
                        'were run from lambda zero to one.'
                        'Default is False.',
                        default=False,
                        action='store_true')
    # The following are mutually exclusive options
    exclus.add_argument('--skip',
                        metavar='',
                        dest='skip',
                        type=int,
                        help='Skip files, i.e. pick every nth work value. '
                        'Default is 1 (all); with 2, every other work value '
                        'is discarded, etc.',
                        default=1)
    exclus.add_argument('--slice',
                        metavar='',
                        dest='sliceit',
                        type=int,
                        help='Subset of trajectories to analyze.'
                        'Provide list slice, e.g. "10 50" will'
                        ' result in selecting dgdl_files[10:50].'
                        ' Default is all.',
                        default=None,
                        nargs=2)
    exclus.add_argument('--rand',
                        metavar='',
                        dest='rand',
                        type=int,
                        help='Take a random subset of trajectories. '
                        'Default is None (do not take random subset)',
                        default=None)
    exclus.add_argument('--index',
                        metavar='',
                        dest='index',
                        type=int,
                        help='Zero-based index of files to analyze (e.g.'
                        ' 0 10 20 50 60). It keeps '
                        'the dgdl.xvg files according to their position in the'
                        ' list, sorted according to the filenames. Default '
                        'is None (i.e. all dgdl are used).',
                        default=None,
                        nargs='+')
    parser.add_argument('--prec',
                        metavar='',
                        dest='prec',
                        type=int,
                        help='The decimal precision of the screen/file output.'
                        ' Default is 2.',
                        default=2)
    parser.add_argument('--units',
                        metavar='',
                        dest='units',
                        type=str.lower,
                        help='The units of the output. Choose from "kJ", '
                        '"kcal", "kT". Default is "kJ."',
                        default='kJ',
                        choices=['kj', 'kcal', 'kt'])
    parser.add_argument('--pickle',
                        dest='do_pickle',
                        help='Whether to save the free energy results from '
                        'the estimators in pickled files. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('--no_ks',
                        dest='do_ks_test',
                        help='Whether to do a Kolmogorov-Smirnov test '
                        'to check whether the Gaussian assumption for CGI '
                        'holds. Default is True; this flag turns it to False.',
                        default=True,
                        action='store_false')
    parser.add_argument('--cgi_plot',
                        metavar='',
                        dest='cgi_plot',
                        type=str,
                        help='Whether to plot the work histograms along with '
                        'the CGI results. If the flag is used, you also need'
                        'to specify a filename.',
                        default=None)
    parser.add_argument('--nbins',
                        metavar='',
                        dest='nbins',
                        type=int,
                        help='Number of histograms bins for the plot. '
                        'Default is 10.',
                        default=10)
    parser.add_argument('--dpi',
                        metavar='',
                        dest='dpi',
                        type=int,
                        help='Resolution of the plot. Default is 300.',
                        default=300)
    parser.add_argument('--quiet',
                        dest='quiet',
                        help='Minimal screen output.',
                        default=False,
                        action='store_true')

    args, unknown = parser.parse_known_args()

    return args


# ====
# MAIN
# ====
def entry_point():
    args = parse_options()
    analyse_dgdl(**vars(args))

if __name__ == '__main__':
    entry_point()
