import argparse 
import GeneTree

def subparser(subparsers):

    desc = 'Create the files for the tree vizualisation in iTOL'

    subparser = subparsers.add_parser('create_tree', description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False
    )

    required = subparser.add_argument_group('required arguments')
    optional = subparser.add_argument_group('optional arguments')

    required.add_argument('-nc', '--ncorespe', metavar='IN_FILE1', type=str,
        required=True,
        help='Num/Name correspondance file from matrice subcommande.'
    )
    required.add_argument('-m', '--matrice', metavar='IN_FILE2', type=str,
        required=True,
        help='Strains X Genes matrices from matrice subcommande.'
    )

    required.add_argument('-n', '--newick_file', metavar='TREE',
        required=True,
        help='path to the clean Newick file.'
    )
    
    required.add_argument('-spn', '--speciesname', metavar='NAME_SPECIES',
        required=True,                  
        help='Species name in the newick file.'
    )
    
    required.add_argument('-clingo', '--clingo_path', metavar='CLINGO_PATH',
        required=True,                  
        help='Path to clingo.'
    )
    
    required.add_argument('-github', '--github_name', metavar='GITHUB_NAME',
        required=True,                  
        help='Name of the git owner of te repository.'
    )
    
    optional.add_argument('-outT', '--outtree', metavar='OUT_TREE', 
        default='<outfile>.nwk',
        help='Output newick file.'
    )
    
    optional.add_argument('-outF', '--outfile', metavar='OUT_FILE', 
        default='<outfile>.txt',
        help='Output popup file.'
    )

    optional.add_argument('-out_dir', '--output_directory', metavar='OUT_DIR',
        default= '<out_dir>',
        help='name of the output directory.'
    )

    optional.add_argument('-ngn', '--new_genes_name', metavar='NEW_NAME_GENES',
        default=None,
        help='Dictionary (pickle file) containing for every gene id a new name (pickle format). ex : {1: A, 2: B, 3: C}'
    )





