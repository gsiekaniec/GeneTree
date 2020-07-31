import argparse 
import GeneTree

def subparser(subparsers):

    desc = 'Create the matrice genes against strains file'

    subparser = subparsers.add_parser('matrice', description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False
    )

    required = subparser.add_argument_group('required arguments')
    optional = subparser.add_argument_group('optional arguments')

    required.add_argument('-tsv', '--tsv_files', nargs='+', metavar='IN_FILE1', type=str,
	required=True,
	help='Tsv files from Microscope (without the pan files).'
    )
    
    optional.add_argument('-out_dir', '--output_directory', metavar='OUT_DIR',
        default= '<out_dir>',
        help='name of the output directory.'
    )


  



