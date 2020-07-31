import GeneTree
from GeneTree import Timer
from GeneTree.clean_name import createDir
import os

D_OUT = '<out_dir>'

def createVector(sts,strains):
    '''Create a vector of 1 and 0 corresponding to an absence or not of the gene in the strain'''
    st = []
    for s in strains:
        if s in sts:
            st.append(1)
        else:
            st.append(0)
    return st
        
def createMatrice(strains,dictRes,out_d):
        strains = sorted(list(strains))
        with open(out_d+'/matrice.tsv','w') as matrix:
            matrix.write('\t')
            for strain in strains:
                matrix.write(str(strain)+'\t')
            matrix.write('\n')
            for k in sorted(dictRes.keys()):
                vector = createVector(dictRes[k],strains)
                matrix.write(str(k)+'\t')
                for i in vector:
                    matrix.write(str(i)+'\t')
                matrix.write('\n')
           
def createName(numName,out_d):
    with open(out_d+'/correspgenesnames.txt','w') as cgn:
        for k,v in numName.items():
            cgn.write(f"{k}\t{' / '.join(list(v))}\n")
            
def createMatriceAndNumNamesFiles(files,out_d):
    ''' Create the matrice genesXstrains using Microscope files 
    (don't take into account singl from var file but take singl from spe file 
    and create new singl.tsv with corresponding singl new names)'''
    
    dictRes = {}
    numName = {}
    strains = set()
    nbsingl = 0
    with Timer() as read_file:
        with open(os.path.dirname(files[0])+'/singl.tsv','w') as singl:        
            for file in files:
                
                # get the file type
                if 'core' in file:
                    typ = 'core'
                elif 'var' in file:
                    typ = 'var'
                elif 'spe' in file:
                    typ = 'spe'
                    
                GeneTree.info(f'{typ} file treatment...')
                    
                with open(file,'r') as tsv:
                    
                    # Don't take first line into account
                    tsv.readline()
                    
                    for line in tsv:
                        treat = True
                        line = line.split('\t')
                        
                        if line[0] == 'singl':
                            if typ == 'var':    
                                treat = False
                        
                        if treat:
                            # gene name
                            gene = line[7]
                            if gene == '':
                                gene = 'Unknow'
                            
                            # strain
                            organism = line[2]
                            
                            # gene number (if single add a digit at the end)
                            
                            if line[0] != 'singl': 
                                num = line[0]
                            else:
                                num = line[0]+str(nbsingl)
                                singl.write(str(line[0]+str(nbsingl))+'\t'+str('\t'.join(line[1:])))
                                nbsingl+=1 
                            
                            # get all strains names
                            
                            if not organism in strains:
                                if ':' in organism:
                                    organism = '_'.join(organism.split(':'))
                                strains.add(organism) 
                            
                            # corresp num gene X organisms
                            
                            if num in dictRes.keys():
                                dictRes[num].add(organism)
                            else:
                                s = set()
                                s.add(organism)
                                dictRes[num]= s
                            
                            # corresp num gene X name gene
                            
                            if num in numName.keys():
                                numName[num].add(gene)
                            else: 
                                n = set()
                                n.add(gene)
                                numName[num]= n 
    GeneTree.info(f'Read files done in: {read_file.t}')
    with Timer() as create_num_name_file_time:
        GeneTree.info(f'nums X names file construction...')
        # Create num x name(s) file
        createName(numName,out_d)
    GeneTree.info(f'nums X names file creation done in: {create_num_name_file_time.t}')
    with Timer() as create_matrice_file_time:
        GeneTree.info(f'matrice construction...')
        # Create matrice      
        createMatrice(strains,dictRes,out_d)
    GeneTree.info(f'genes X strains matrice creation done in: {create_matrice_file_time.t}')
    
def main(args):
    
    # Check the inputs files
    GeneTree.info('Check inputs...')
    for file in args.tsv_files:
        if not os.path.exists(file):
            GeneTree.crit(f'Failed to open {file}', 1)
    
    # Test if the output directory exist and create one if not
    if args.output_directory == D_OUT:
        args.output_directory = 'output'
    if os.path.exists(args.output_directory):
        GeneTree.warn(f'Usage of an already existant output directory {args.output_directory}')
    else:
        createDir(args.output_directory)
        
    params = f'\ttsv files -> {args.tsv_files}\n\
               \tout directory -> {args.output_directory}\n\
               '
    GeneTree.info(f'Parameters:\n \t{params}')
    
    # Full time calculation
    with Timer() as full_time:
        createMatriceAndNumNamesFiles(args.tsv_files,args.output_directory)
    GeneTree.info(f'Matrice creation  done in: {full_time.t}')
    


