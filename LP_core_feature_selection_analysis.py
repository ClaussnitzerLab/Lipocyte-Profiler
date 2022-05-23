import networkx as nx
import pandas as pd
import numpy as np
import sys

def getnetwork(fullnetwork):

    #Parse input
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
    
#Get mean MI per node
def MI_nodes(toprocess):
    df2 = toprocess.groupby(['Regulator', 'Category'], as_index = False).size()
    new = toprocess.merge(df2, on =['Regulator', 'Category'], how='left')
    out = new.groupby(['Regulator', 'Category','size'], as_index = False).agg(
        MI_score = ('MI', 'mean'))
    
    return(out)

#MI_average ranking cutoffs 
def MI_based_cutoffs(entire,q):
    entire = mibynode
    edge_filtered = entire[entire.Category == "primary"]
    tmp = entire[entire.Category == "all"]
    MI_filtered = tmp[tmp.MI_score < np.quantile(tmp['MI_score'],q)]
    percutoff = np.quantile(tmp['MI_score'],q)
    print("The cutoff for the 25 percentile is: %.2f" % percutoff)
    final_set = set(list(edge_filtered["Regulator"]) + list(MI_filtered["Regulator"]))
    size = len(final_set)
    print("The total size of core features is: ", size)
    toout = pd.DataFrame(final_set, columns = ['Features'])
    toout.to_csv('LipocyteProfiler_core_list.tsv')




#main
f = sys.argv[1] #file name
q = float(sys.argv[2]) #Percentile
p = float(sys.argv[3])
fullnetwork = pd.read_csv(f, sep="\t")
ranked = getnetwork(fullnetwork)
nfeatures = countfeatures(ranked, q)
ffeatures = nfeatures[0]

#with open("hubfeatures.tsv", 'w') as out:
#	for i in ffeatures:
#        	out.write(i+"\n")
            
c_all = fullnetwork
c_filter = pd.DataFrame(ffeatures, columns = ['Regulator'])
new = pd.merge(fullnetwork, c_filter, on = 'Regulator', how = 'inner')
core = ["primary"]*len(new.index)
new['Category'] = core
allf = ["all"]*len(fullnetwork.index)
c_all['Category'] = allf

toprocess = pd.concat([new, c_all])
mibynode = MI_nodes(toprocess)
MI_based_cutoffs(mibynode,p)  

