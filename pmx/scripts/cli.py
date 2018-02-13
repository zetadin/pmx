#!/usr/bin/env python2

# simple interface to the available scripts

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
        analyse    Estimate free energy from Gromacs xvg files
        mutate     Mutate protein or DNA''',
            formatter_class=RawTextHelpFormatter)

        parser.add_argument('command', help=SUPPRESS)
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print 'Unrecognized command'
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def analyse(self):
        import analyze_dgdl
        analyze_dgdl.entry_point()

    def mutate(self):
        pass


def entry_point():
    PmxCli()

if __name__ == '__main__':
    entry_point()
