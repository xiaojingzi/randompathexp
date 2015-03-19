import sys
import networkx as nx
import time
import re
import numpy as np

if len(sys.argv) < 2:
    print "usage: dataset"
    print "output: dataset + '.RoleSimGraph'"
    exit()

dataset=sys.argv[1]
networkFile= "data/"+dataset+".graph"
mapFile = "data/"+dataset+".dict"

G = nx.DiGraph()

for line in open(mapFile, "r"):
    elements = re.compile('\t').split(line.strip())
    if len(elements) == 2:
        name = elements[0]
        id = int(elements[1])
    elif len(elements) == 1:
        name = "no name"
        id = int(elements[0])
    G.add_node(id, name=name, weight=0.0)
N = G.number_of_nodes()

i = 0
for line in open(networkFile, "r"):
    if i == 0:
        elements = re.compile('\t').split(line.strip())
    else:
        elements = re.compile('\t').split(line.strip())
        node1 = int(elements[0])
        node2 = int(elements[1])
        w = float(elements[2])
        G.add_edge(node1, node2, weight=w)
        # G.add_edge(node2, node1, weight=w)
    i +=1
    
M = G.number_of_edges()

print "N=",N,",M=",M
f = open("data/"+dataset +'.RoleSimGraph' , 'w')
f.write('graph_for_greach' + '\n')
f.write(str(N) + '\n')
edgeNum = 0
for u in G.nodes():
    f.write(str(u) + ': ')
    edgeNum += len(G.neighbors(u))
    for v in G.neighbors(u):
        f.write(str(v) + ' ')
    f.write('#\n')
f.close()
print "Edge number = ", edgeNum
    
