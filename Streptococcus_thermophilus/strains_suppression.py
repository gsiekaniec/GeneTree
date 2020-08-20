from ete3 import Tree
import argparse

def cli_options():
    "Return parsed cli options"
    parser = argparse.ArgumentParser(description=__doc__)
    required = parser.add_argument_group('Required argument')
    required.add_argument('-f', '--files', nargs='+', metavar='IN_FILE', type=str,
        required=True,
        help='path to the tsv files from Microscope'
    )
    required.add_argument('-n', '--newick', metavar='TREE_FILE', type=str,
        required=True,
        help='path to the tree files from Microscope'
    )
    return parser.parse_args()


def main():
    args = cli_options()
    # Get the strains names
    t = Tree(args.newick, format=1)
    strains = set()
    for leaf in t:
        strains.add(leaf.name) 
    #Suppr from files
    for file in args.files:
        new_file = file.split('.tsv')[0]+'_new.tsv'
        with open(file,'r') as f:
            with open(new_file,'w') as o:
                line = f.readline()
                o.write(line)
                for line in f:
                    strain = line.split('\t')[2]
                    strain = ''.join(strain.split(':'))
                    if strain in strains:
                        o.write(line)
    
    
if __name__ == "__main__":
    main()