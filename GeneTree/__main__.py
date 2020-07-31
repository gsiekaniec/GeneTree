#!/usr/bin/env python3
# coding: utf-8

import GeneTree

def run(arglist=None):
    args = GeneTree.cli.parse_args(arglist)

    if args.cmd not in GeneTree.cli.mains:
        GeneTree.cli.parser("").parse_args(['-h'])

    GeneTree.logger.info(f'# GeneTree - {args.cmd} #')
    maincmd = GeneTree.cli.mains[args.cmd]
    GeneTree.logger.info(f'GeneTree v{GeneTree.__version__}')

    maincmd(args)

if __name__ == '__main__':
    run()
else:
    print (__name__)
