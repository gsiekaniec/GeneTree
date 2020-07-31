import argparse
import sys

import GeneTree

from . import clean_name_cli
from . import matrice_cli
from . import create_tree_cli

sub_commands = ['clean_name', 'matrice', 'create']

mains = {
    'clean_name':     GeneTree.clean_name.main,
    'matrice':   GeneTree.matrice.main,
    'create_tree':   GeneTree.create_tree.main
}

subparsers_funcs = {
	'clean_name': clean_name_cli.subparser,
	'matrice': matrice_cli.subparser,	
	'create_tree': create_tree_cli.subparser
}

def parser(description: str=None) -> argparse.ArgumentParser:
    sub_commands_str = '", "'.join(sorted(sub_commands))
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
        )
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'

    parser.add_argument('-l', '--logfile', metavar='FILE',
        help='log file for diagnostic messages warnings, and errors'
    )
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Verbose mode'
    )
    parser.add_argument('--version', action='version',
        version=f'GeneTree v{GeneTree.__version__}',
        help='Display SRcompressor version'
    )
    parser.add_argument('-h', '--help', action='help',
        help='Show this message and exit'
    )
    
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
        help='"' + sub_commands_str + '"'
    )

    for func in subparsers_funcs.values():
        func(subparsers)

    return parser

def parse_args(arglist=None):
    args = GeneTree.cli.parser().parse_args(arglist)

    GeneTree.logstream = sys.stderr

    if args.logfile:
        GeneTree.logfile = args.logfile

    GeneTree.setup_logger('logger', args.verbose, GeneTree.logfile)
    GeneTree.logger = GeneTree.logging.getLogger('logger')

    return args









