%*
Data:
 a tsv file representing a matrix genesXstrains of presence.   

*%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#script (python)

import clingo

def readgenenames(filename):
    l=[]
    with open(filename.string) as fd:
        for line in fd:
            sline=line.strip().split("\t")
            l.append((sline[0].strip(), sline[1].strip()))
    return(l)

def readmatrix(filename):
    valueline=False
    l=[]
    with open(filename.string) as fd:
        for line in fd:
            sline=line.strip().split("\t")
            if valueline:     
                l.extend([(strain[i],sline[0]) for i,x in enumerate(sline[1:]) if x=="1"])
            else : 
                valueline=True
                strain=[x.strip() for x in sline]
    return(l)

#end.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASP Program



gene(N,G):- (N,G) =@readgenenames("genenames.txt").

straingene(S,N):- (S,N)= @readmatrix("matrice.tsv").

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

%#show gene/2.
%#show strain/2.
