import networkx as nx
import pandas as pd
import numpy as np
import sys

def getnetwork(file):

    #Parse input
    fullnetwork = pd.read_csv(file, sep="\t")
    target = list(fullnetwork['Target'])
  
    regulator = list(fullnetwork['Regulator'])
    
    pairs = []
    for i in range(0,len(target)):
        pairs.append([regulator[i], target[i]])


    #Create a directional graph object
    g = nx.DiGraph()
    g.add_edges_from(pairs)

    rankfeatures = {}
    for i in list(g.nodes):
        targets = list(g.successors(i))
        rankfeatures[i] = targets

    return(rankfeatures)


def countfeatures(mappedfeatures, q):
    
    rankfeatures={}
    for i in mappedfeatures:
        rankfeatures[i] = len(mappedfeatures[i])
        
    node_order = list(rankfeatures.values())
    cut = np.quantile(list(node_order), q)
    
    selectedf = []
    for i in rankfeatures:
        if rankfeatures[i] > cut:
            selectedf.append(i)
    
    ret = (selectedf,cut)
    return(ret)


#main
f = sys.argv[1] #file name
q = float(sys.argv[2]) #Percentile
ranked = getnetwork(f)
nfeatures = countfeatures(ranked, q)
ffeatures = nfeatures[0]


with open("hubfeatures.tsv", 'w') as out:
	for i in ffeatures:
        	out.write(i+"\n")
            