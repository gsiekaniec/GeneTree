import argparse 
import GeneTree

def subparser(subparsers):

    desc = 'Clean strains names from files'

    subparser = subparsers.add_parser('clean_name', description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False
    )

    required = subparser.add_argument_group('required arguments')
    optional = subparser.add_argument_group('optional arguments')

    required.add_argument('-c', '--core_file', metavar='IN_FILE', type=str,
        required=True,
        help='path to the core csv file from Microscope'
    )

    required.add_argument('-n', '--newick_file', metavar='TREE',
        required=True,
        help='path to the Newick file from Microscope.'
    )
    
    optional.add_argument('-o', '--out', metavar='OUT_TREE',
        default= '<out>.nwk',
        help='name of the new Newick file.'
    )
    
    optional.add_argument('-out_dir', '--output_directory', metavar='OUT_DIR',
        default= '<out_dir>',
        help='name of the output directory.'
    )
