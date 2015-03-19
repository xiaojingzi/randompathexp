import sys
import time
import networkx as  nx
import numpy as np 
import re
import scipy.io 
if len(sys.argv) < 3:
    print "usage: python transfer2rwrformat.py dataset datasetpath"
    exit()

dataset = sys.argv[1]
datasetpath = sys.argv[2]
network_file = datasetpath+".graph"
map_file = datasetpath+".dict"


G = nx.Graph()

# read node map
for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    if len(elements) == 2:
        name = elements[0]
        id = int(elements[1])
    elif len(elements) == 1:
        name = "no name"
        id = int(elements[0])
    G.add_node(id, name=name, weight=0.0)
print "load nodes done."
node_num = G.number_of_nodes()
W = np.zeros((node_num, node_num))
# read network
i = 0
for line in open(network_file, "r"):
    if i == 0:
        i+=1
        continue
    i+=1
    elements = re.compile('\t').split(line.strip())
    node1 = int(elements[0])
    node2 = int(elements[1])
    w = float(elements[2])
    # G.add_edge(node1, node2, weight=w)
    # G.node[node1]['weight'] += w
    # G.add_edge(node2, node1, weight=w) 
    # G.node[node2]['weight'] += w
    W[node1][node2] = w
    W[node2][node1] = w

edge_num = G.number_of_edges()

scipy.io.savemat(dataset+'.mat',{'W': W})


