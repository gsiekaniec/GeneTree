import GeneTree,clyngor, os
from GeneTree import Timer
from GeneTree.clean_name import createDir
from clyngor import solve
from ete3 import Tree

D_OUT = '<out_dir>'
T_OUT = '<outfile>.nwk'
N_OUT = '<outfile>.txt'

ASP_CORE="""\
gene(N):- rel(_,N).
strain(S):- rel(S,_).

nb(gene,N):- N={gene(_)}.
nb(strain,N):- N={strain(_)}.

vargene(G):- gene(G); strain(S); not rel(S,G).
comment("Genes present in all strains (core genome)").
coregene(N,G):- gene(N,G); not vargene(N).
comment("Genes present in a single strain (specific)").
specgene(S,N,G):-  gene(N,G); rel(S,N), 1={rel(_,N)}.

comment("Statistics").
comment("vargene: Genes absent in at least one strain").
nb(vargene,N):- N={vargene(_)}.
comment("coregene: Genes present in all strains").
nb(coregene,N):- N={coregene(_,_)}.
comment("specgene: Genes present in a single strain (specific, signature gene)").
nb(specgene,N):- N={specgene(_,_,_)}.

#show coregene/2.
#show specgene/3.
#show nb/2.
#show comment/1.
"""

ASP_VAR="""\
strain(S):- straingene(S,_).
gene(N):- straingene(_,N).

vargene(N):- gene(N); strain(S); not straingene(S,N).
specgene(N):-  gene(N,G); straingene(S,N), 1={straingene(_,N)}.

%coregene(G):- gene(G); not vargene(G).
%nb(vargene,N):- N={vargene(_)}.

rel(S,N):- straingene(S,N); vargene(N); not specgene(N).

strainvar(S):-  rel(S,_).
genevar(N):-  rel(_,N).

nb(gene_var,N):- N={genevar(_)}.
nb(strain_var,N):- N={strainvar(_)}.

%Compute concepts

ext(X):- rel(X,_) ; rel(X,Y): int(Y).
int(Y):- rel(_,Y) ; rel(X,Y): ext(X).
        
% Avoid non-concept (no object or no attribute)
:- not ext(_).
:- not int(_).

%Compute AOC poset

% Is outsider any object or attribute that is linked to an attribute or object not in the concept.
ext_outsider(X):- ext(X) ; rel(X,Z) ; X!=Z ; not int(Z).
int_outsider(Y):- int(Y) ; rel(Z,Y) ; Y!=Z ; not ext(Z).

% Specific part of each concept.
%specext(X):- ext(X) ; not ext_outsider(X).
specint(Y):- int(Y) ; not int_outsider(Y).
%:- not specext(_).
:- not specint(_).

#show.
%#show coregene/1.
%#show nb/2.

%#show ext/1.
%#show int/1.

#show strain(S): ext(S).
#show spgene(N,G): specint(N); gene(N,G).
"""

##### EXTRACTION (ASP part)

def GenesCoreSpe (asp):
    ''' Get genes strain specific and core genome '''
    dictionary = {'core':set(),'spe':{}}
    toExec = asp+ASP_CORE
    answers = solve(inline=toExec)
    for answer in answers:
        for res in answer:
            if res[0] == 'coregene':
                num = res[1][0][1:-1]
                name = res[1][1][1:-1]
                dictionary['core'].add((num,name))
            elif res[0] == 'specgene':
                strain = res[1][0].split('"')[1]
                num = res[1][1]
                name = res[1][2]
                if strain in dictionary['spe']:
                    dictionary['spe'][strain].add((num,name))
                else:
                    dictionary['spe'][strain]=set()
                    dictionary['spe'][strain].add((num,name))
    return dictionary
        
def GenesVar (asp):
    ''' Get maximal biclusters {strains}X{genes} which contain at least 2 strains and not the core genome '''
    dictionary = {'var':{}}
    toExec = asp+ASP_VAR
    answers = solve(inline=toExec)
    for answer in answers:
        strains = set()
        genes = set()
        for res in answer:
            if res[0] == 'strain':
                strains.add(res[1][0][1:-1])
            elif res[0] == 'spgene':
                genes.add((res[1][0][1:-1],res[1][1][1:-1]))
        strains = frozenset(strains)
        if strains == frozenset(['Streptococcus thermophilus 1F8CT', 'Streptococcus thermophilus CIRM30']):
            print(strains,genes)
        if strains == frozenset(['Streptococcus thermophilus LMG 18311', 'Streptococcus thermophilus CIRM23']):
            print(strains,genes)
        if strains != set():
            dictionary['var'][strains]=genes
    return dictionary

def readgenenames(filename):
    '''Get genes names from file '''
    gene = ''
    l=[]
    with open(filename) as fd:
        for line in fd:
            sline=line.strip().split("\t")
            l.append('gene("'+str(sline[0].strip())+'\","'+str(sline[1].strip())+'\")')
    gene = str('. '.join(l))+'.'
    return(gene)

def readmatrix(filename):
    '''Get relation between genes and strains from matrice'''
    valueline=False
    l=[]
    strain = []
    with open(filename) as fd:
        for line in fd:
            sline=line.strip().split("\t")
            if valueline:     
                l.extend(['rel("'+str(strain[i])+'\","'+str(sline[0])+"\")" for i,x in enumerate(sline[1:]) if x=="1"])
            else : 
                valueline=True
                strain=[x.strip() for x in sline]
    rel = str('. '.join(l))+'.'
    return(rel)
      
########
##### CREATION    

def newickTreatment(newick,coreDict,varDict,out_tree,output,speciesname,out_dir,ngn=None):
    
    repositryname = out_dir+'/'+'_'.join(speciesname.split(' '))+'_data'
    if not os.path.exists(repositryname):
        createDir(repositryname)
    else:
        GeneTree.warn(f'Usage of an already existant directory {repositryname}')
    
    # read tree
    t = Tree(newick, format=1)

    GeneTree.info(f'the newick file contains: \n {t}')
    
    with open(output,'w') as popup:
        popup.write('POPUP_INFO\nSEPARATOR TAB\nDATA\n')
    
        number = 0
        # Tree walk
        #print(sorted(list(coreDict['spe'].keys())))
        for node in t.traverse("preorder"):
            itol_str = ''
            if node.is_root():
                pass
            elif node.name == 'root':
                with open(str(repositryname)+'/root.txt','w') as out:
                    itol_str += f"{node.name}\t {node.name} \t <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/GeneTree/master/Streptococcus_thermophilus/Streptococcus_thermophilus_data/root.txt'>&dagger;</a> </p>"
                    for gene in coreDict['core']:
                        if not ngn:
                            itol_str += f"id: {gene[0]} - name: {gene[1]} "
                            out.write(f"id: {gene[0]} - name: {gene[1]}\n")
                        else:
                            itol_str += f"id: {gene[0]} - name: {ngn[gene[0]]} ; "
                            out.write(f"id: {gene[0]} - name: {ngn[gene[0]]}\n")
                    popup.write(itol_str[:-1])
            elif node.is_leaf():
                filename = '_'.join(str(node.name).split(' '))+'.txt'
                with open(str(repositryname)+'/'+filename,'w') as out:
                    itol_str += f"{node.name}\t {node.name} \t  <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/GeneTree/master/Streptococcus_thermophilus/Streptococcus_thermophilus_data/{filename}'>&dagger;</a></p>"
                    name = node.name.strip()
                    if name in sorted(list(coreDict['spe'].keys())):
                        for i,gene in enumerate(coreDict['spe'][name]):
                            if not ngn:
                                itol_str +=  f"id: {gene[0]} - name: {gene[1]} "
                                out.write(f"id: {gene[0]} - name: {gene[1]}\n")
                            else:
                                itol_str += f"id: {gene[0]} - name: {ngn[gene[0]]} "
                                out.write(f"id: {gene[0]} - name: {ngn[gene[0]]}\n")
                            if i != int(len(coreDict['spe'][name])-1):
                                itol_str += " ; "
                    else:
                        GeneTree.info(f'Info : {name} do not contains specific genes')
                        itol_str +=  f"<p style='color:blue'>No specific genes </p>"
                    popup.write(itol_str)
            elif node.name == '':
                node.name = 'node '+str(number)
                filename = '_'.join(str(node.name).split(' '))+'.txt'
                with open(str(repositryname)+'/'+filename,'w') as out:
                    strains = set()
                    genes = ''
                    itol_str = f"node {str(number)}\tnode {str(number)}\t  <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/GeneTree/master/Streptococcus_thermophilus/Streptococcus_thermophilus_data/{filename}'>&dagger;</a></p>"
                    for leaf in node:
                        strains.add(str(leaf.name).strip())
                    try :
                        genes = varDict['var'][frozenset(strains)]
                        for i,g in enumerate(genes):
                            if not ngn:
                                itol_str +=  f"id: {g[0]} - name: {g[1]} "
                                out.write(f"id: {g[0]} - name: {g[1]}\n")
                            else:
                                itol_str += f"id: {g[0]} - name: {ngn[g[0]]} "
                                out.write(f"id: {g[0]} - name: {ngn[g[0]]}\n")
                            if i != int(len(genes)-1):
                                itol_str += " ; "
                    except KeyError:
                        itol_str += f"<p style='color:blue'>No specific genes </p> "
                        out.write("No specific genes\n")
                    number += 1
                    popup.write(itol_str)
            popup.write('\n')
            
    t.write(format=1, outfile=out_tree)

########

def main(args):
    
     # Check the inputs files
    GeneTree.info('Check inputs...')
    if not os.path.exists(args.ncorespe):
        GeneTree.crit(f'Failed to open {args.ncorespe}', 1)
    if not os.path.exists(args.matrice):
        GeneTree.crit(f'Failed to open {args.matrice}', 1)
    if not os.path.exists(args.newick_file):
        GeneTree.crit(f'Failed to open {args.newick_file}', 1)
    if not os.path.exists(args.clingo_path):
        GeneTree.crit(f'Failed to open {args.clingo_path}', 1)
    
    # Test if the output directory exist and create one if not
    if args.output_directory == D_OUT:
        args.output_directory = 'output'
    if os.path.exists(args.output_directory):
        GeneTree.warn(f'Usage of an already existant output directory {args.output_directory}')
    else:
        createDir(args.output_directory)
    
    # if no out tree create, get the out name from the input name  
    if args.outtree == T_OUT:
        args.outtree = os.path.splitext(args.newick_file)[0]
        args.outtree = args.output_directory+'/'+os.path.basename(args.outtree)+'_final.nwk'
        
    # if no out file create, get default the out name  
    if args.outfile == N_OUT:
        args.outfile = args.output_directory+'/popup.txt'
        
    params = f'\tnums X names file -> {args.ncorespe}\n\
               \tmatrice file -> {args.matrice}\n\
               \tnewick file -> {args.newick_file}\n\
               \tpath clingo -> {args.clingo_path}\n\
               \tout newick file -> {args.outtree}\n\
               \tout popup file -> {args.outfile}\n\
               \tout directory -> {args.output_directory}\n\
               '
    GeneTree.info(f'Parameters:\n \t{params}')
    
    # Full time calculation
    with Timer() as full_time:
    
        clyngor.CLINGO_BIN_PATH = args.clingo_path
        
        
        gene = readgenenames(args.ncorespe)
        rel = readmatrix(args.matrice)
        
        asp = gene+rel
        
        GeneTree.info(f'core/spe genes extraction...')
        with Timer() as coreSpe_time:
            # Dictionary Core/Spe
            coreSpe = GenesCoreSpe(asp)
        GeneTree.info(f'core/spe genes extraction done in: {coreSpe_time.t}')
        
        GeneTree.info(f'maximal biclusters (strains)X(genes) extraction...')
        with Timer() as var_time:
            # Dictionary Var
            var = GenesVar(asp)
        GeneTree.info(f'maximal biclusters (strains)X(genes) extraction done in: {var_time.t}')
        
        
       
        if args.new_genes_name:
            import pickle
            with open(args.new_genes_name,'rb') as infile:
                ngn = pickle.load(infile)
        else:
            ngn = None
            
        with Timer() as tree_creation_time:
            newickTreatment(args.newick_file,coreSpe,var,args.outtree,args.outfile,args.speciesname,args.output_directory,ngn)
        GeneTree.info(f'tree creation done in: {tree_creation_time.t}')
        
    GeneTree.info(f'done in: {full_time.t}')


