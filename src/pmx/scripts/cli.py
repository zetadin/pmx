#!/usr/bin/env python2

# simple interface to the available scripts
from __future__ import print_function, absolute_import
from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS
import sys


class PmxCli(object):

    def __init__(self):
        parser = ArgumentParser(
            description='''
    ------------------------
    pmx command line scripts
    ------------------------

    Available commands are:
        mutate     Mutate protein or DNA/RNA
        filltop    Fill hybrid topology with B states
        analyse    Estimate free energy from Gromacs xvg files

        genlib     Generate pmx ff library
        gmxlib     Show/set GMXLIB path''',
            formatter_class=RawTextHelpFormatter)

        parser.add_argument('command', help=SUPPRESS)
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def mutate(self):
        from . import mutate
        mutate.entry_point()

    def filltop(self):
        from . import fill_topology_bstates
        fill_topology_bstates.entry_point()

    def analyse(self):
        from . import analyze_dgdl
        analyze_dgdl.entry_point()

    def genlib(self):
        from . import generate_mutant_lib
        generate_mutant_lib.entry_point()

    def gmxlib(self):
        from . import set_gmxlib
        set_gmxlib.entry_point()


def entry_point():
    PmxCli()


if __name__ == '__main__':
    entry_point()
