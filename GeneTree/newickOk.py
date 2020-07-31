from ete3 import Tree
import argparse, os


def cli_options():
    "Return parsed cli options"
    parser = argparse.ArgumentParser(description='')
    required = parser.add_argument_group('Required argument')
    optional = parser.add_argument_group('Optional argument')
    required.add_argument('--core', '-c', dest='core', help='ASP result for core/spefile.')
    required.add_argument('--var', '-v', dest='var', help='ASP result for varfile.')
    required.add_argument('--newick', '-n', dest='newick', help='Newick file.')
    optional.add_argument('--output', '-out', dest='output', help='Output newick file.')
    required.add_argument('--speciesname', '-spn', dest='speciesname', help='Species name in the newick file.')
    optional.add_argument('--new_genes_name', '-ngn', dest='ngn', help='Dictionary containing for every gene id a new name (pickle format).')
    return parser

def parseCoreFile(core):
    dictionary = {'core':set(),'spe':{}}
    with open(core,'r') as f:
        for line in f:
            line = line.split('\n')[0]
            if line.split('(')[0] == 'coregene':
                num = str(line.split('(')[1].split(')')[0].split(',')[0])[1:-1]
                name = str(line.split('(')[1].split(')')[0].split(',')[1])[1:-1]
                dictionary['core'].add(tuple([num,name]))
            elif line.split('(')[0] == 'specgene':
                souche = line.split('(')[1].split(')')[0].split(',')[0]
                souche = '_'.join(souche.split(':'))
                num = str(line.split('(')[1].split(')')[0].split(',')[1])[1:-1]
                name = str(line.split('(')[1].split(')')[0].split(',')[2])[1:-1]
                if souche[1:-1] in dictionary['spe'].keys():
                    dictionary['spe'][souche[1:-1]].add(tuple([num,name]))
                else:
                    dictionary['spe'][souche[1:-1]]=set()
                    dictionary['spe'][souche[1:-1]].add(tuple([num,name]))
    return dictionary

def parseVarFile(var):
    dictionary = {'var':{}}
    with open(var,'r') as f:
        for line in f:
            strains = set()
            genes = set()
            line = line.strip()
            line = line.split(')')
            for char in line:
                char = char.strip().split('(')
                if char[0] == 'strain':
                    souche = '_'.join(char[1][1:-1].split(':'))
                    strains.add(souche)
                elif char[0] == 'spgene':
                    gene = char[1].split(',')
                    genes.add((gene[0][1:-1],gene[1][1:-1]))
            if strains != set():
                dictionary['var'][frozenset(strains)]=genes
    return dictionary
         
def createDir(path : str):
    ''' Create the directory if it don't exist'''
    try:
        os.mkdir(path)
    except FileExistsError:
        print ('Warning repertory exist !!!\n')
        pass
       
def newickTreatment(newick,core,var,output,speciesname,ngn=None):
    
    repositryname = '_'.join(speciesname.split(' '))+'_data'
    createDir(repositryname)
    
    if not output:
        output = 'f.txt'
    
    t = Tree(newick, format=1)


    print('The newick file contains : \n', t)

    coreDict = parseCoreFile(core)
    varDict = parseVarFile(var)
    
    with open(output,'w') as popup:
        popup.write('POPUP_INFO\nSEPARATOR TAB\nDATA\n')
    
        number = 0
        for node in t.traverse("preorder"):
            itol_str = ''
            if node.is_root():
                pass
            elif node.name == 'root':
                with open(str(repositryname)+'/root.txt','w') as out:
                    itol_str += f"{node.name}\t {node.name} \t <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/Streptococcus_thermophilus_data/master/root.txt'>&dagger;</a> </p>"
                    for gene in coreDict['core']:
                        if not ngn:
                            itol_str += f"id: {gene[0]} - name: {gene[1]}"
                            out.write(f"id: {gene[0]} - name: {gene[1]}\n")
                        else:
                            itol_str += f"id: {gene[0]} - name: {ngn[gene[0]]} ; "
                            out.write(f"id: {gene[0]} - name: {ngn[gene[0]]}\n")
                    popup.write(itol_str[:-1])
            elif node.is_leaf():
                filename = '_'.join(str(node.name).split(' '))+'.txt'
                with open(str(repositryname)+'/'+filename,'w') as out:
                    itol_str += f"{node.name}\t {node.name} \t  <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/Streptococcus_thermophilus_data/master/{filename}'>&dagger;</a></p>"
                    name = node.name.split(speciesname)[-1].strip()
                    if name in sorted(list(coreDict['spe'].keys())):
                        for i,gene in enumerate(coreDict['spe'][name]):
                            if not ngn:
                                itol_str +=  f"id: {gene[0]} - name: {gene[1]}"
                                out.write(f"id: {gene[0]} - name: {gene[1]}\n")
                            else:
                                itol_str += f"id: {gene[0]} - name: {ngn[gene[0]]}"
                                out.write(f"id: {gene[0]} - name: {ngn[gene[0]]}\n")
                            if i != int(len(coreDict['spe'][name])-1):
                                itol_str += " ; "
                    else:
                        print('Info : ',name,'do not contains specific genes')
                        itol_str +=  f"<p style='color:blue'>No specific genes </p>"
                    popup.write(itol_str)
            elif node.name == '':
                node.name = 'node '+str(number)
                filename = '_'.join(str(node.name).split(' '))+'.txt'
                with open(str(repositryname)+'/'+filename,'w') as out:
                    strains = set()
                    genes = ''
                    itol_str = f"node {str(number)}\tnode {str(number)}\t  <p><FONT size='4pt'><b>Genes list</b></FONT>  <a href='https://raw.githubusercontent.com/gsiekaniec/Streptococcus_thermophilus_data/master/{filename}'>&dagger;</a></p>"
                    for leaf in node:
                        strains.add(str(leaf.name).split(speciesname)[-1].strip())
                    try :
                        genes = varDict['var'][frozenset(strains)]
                        for i,g in enumerate(genes):
                            if not ngn:
                                itol_str +=  f"id: {g[0]} - name: {g[1]}"
                                out.write(f"id: {g[0]} - name: {g[1]}\n")
                            else:
                                itol_str += f"id: {g[0]} - name: {ngn[g[0]]}"
                                out.write(f"id: {g[0]} - name: {ngn[g[0]]}\n")
                            if i != int(len(genes)-1):
                                itol_str += " ; "
                    except KeyError:
                        itol_str += f"<p style='color:blue'>No specific genes </p> "
                        out.write("No specific genes\n")
                    number += 1
                    popup.write(itol_str)
            popup.write('\n')
            
    t.write(format=1, outfile="new_tree.nw")


def main():
    parser = cli_options()
    options = parser.parse_args()
    if options.ngn:
        import pickle
        with open(options.ngn,'rb') as infile:
            ngn = pickle.load(infile)
    else:
        ngn = None
    try :
        newickTreatment(options.newick,options.core,options.var,options.output,options.speciesname,ngn)
    except AttributeError as error:
        print(parser.print_help())
        print('Caught this error: ' + repr(error))

if __name__ == "__main__":
    main()

