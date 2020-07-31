import GeneTree
from GeneTree import Timer
from ete3 import Tree
from difflib import SequenceMatcher
import os

D_OUT = '<out_dir>'
N_OUT = '<out>.nwk'

def similar(a, b):
    '''Take two str (a and b) and return the similarity pourcentage between them'''
    return SequenceMatcher(None, a, b).ratio()

def namesFromCore (core: str) -> set:
    '''Take the core file and returna set containing the strains names in this file'''
    setCleanName = set()
    with open(core,'r') as f:
        f.readline() 
        for line in f:
            name = line.split('\t')[2]
            if ':' in name:
                name = '_'.join(name.split(':'))
            setCleanName.add(name)
    return list(setCleanName)
    
def cleanNames(core,newick,out,out_dir):
    '''Write a new newick with corresponding strains names from the core file'''
    # Get the strains names from the core file 
    setCleanName = namesFromCore(core)
    # Get the Microscope tree
    t = Tree(newick, format=1)
    for leaf in t:
        # If the leaf name is not in names from core we search for the closest name
        if not leaf.name in setCleanName:
            sim = []
            for name in setCleanName:
                sim.append(similar(name,leaf.name))
            maxElement = max(sim)
            leaf.name = setCleanName[int(sim.index(maxElement))]
    # Write the new tree (newick format)
    t.write(format=1, outfile=out_dir+'/'+out+'temp')
    # Add un a root node for a better visualization on iTOL
    with open(out_dir+'/'+out,'w') as o:
        with open(out_dir+'/'+out+'temp','r') as f:
            line = f.readlines()
            o.write('(')  
            o.write(line[0][:-1])
            o.write('root:0.0001);')
    os.remove(out_dir+'/'+out+'temp')

def createDir(path : str):
    ''' Create the directory if it don't exist'''
    try:
        os.mkdir(path)
    except FileExistsError:
        GeneTree.crit(f'Output repository already exist: {path}', 1)

def main(args):

    # Check the inputs files
    GeneTree.info('Check inputs...')
    if not os.path.exists(args.core_file):
        GeneTree.crit(f'Failed to open {args.core_file}', 1)
    if not os.path.exists(args.newick_file):
        GeneTree.crit(f'Failed to open {args.newick_file}', 1)
    
    # if no out create, get the out name from the input name  
    if args.out == N_OUT:
        args.out = os.path.splitext(args.newick_file)[0]
        args.out = os.path.basename(args.out)+'_clean.nwk'
        
    # Create the output directory if it don't exist
    if args.output_directory == D_OUT:
        args.output_directory = 'output'
    createDir(args.output_directory)
        
    params = f'\tcore file -> {args.core_file}\n\
               \tnewick file -> {args.newick_file}\n\
               \tout newick file -> {args.out}\n\
               \tout directory -> {args.output_directory}\n\
               '
    GeneTree.info(f'Parameters:\n \t{params}')
    
    # Full time calculation
    with Timer() as full_time:
        cleanNames(args.core_file,args.newick_file,args.out,args.output_directory)
    GeneTree.info(f'Cleaning done in: {full_time.t}')

