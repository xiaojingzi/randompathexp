import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict


if len(sys.argv) < 2:
    print "usage: dataset"
    print "output: dataset + '.cn'"
    exit()

dataset=sys.argv[1]
networkFile = "data/"+dataset+".graph"
mapFile = "data/"+dataset+".dict"

G = nx.Graph()

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


i =0
for line in open(networkFile, "r"):
    if i ==0:
        es = line.split('\t')
    else:
        elements = re.compile('\t').split(line.strip())
        node1 = int(elements[0])
        node2 = int(elements[1])
        w = float(elements[2])
        G.add_edge(node1, node2, weight=w)
        G.add_edge(node2, node1, weight=w)
    i +=1
M = G.number_of_edges()
print "node number = ", N , ", edge number = ", M

t1 = time.time()
comNb = defaultdict(int)
for u in G.nodes():
    u_neighbors = G.neighbors(u)
    for n1 in u_neighbors:
        for n2 in u_neighbors:
            comNb[str(n1)+"_"+str(n2)] +=1

print "pairs have cn = ",len(comNb)

# count = 0
# for u in G.nodes():
#     for v in G.nodes():
#         cn = sorted(nx.common_neighbors(G, u, v) )
#         if len(cn) >  0:
#             comNb[str(u)+"_"+str(v)] +=len(cn)
# print comNb["1121_1121"],
# print len(sorted(nx.common_neighbors(G, 1121,1121) ))

outputFile = "result/"+dataset+".cn"
f = open(outputFile, "w")
for k, v in comNb.items():
    keys = k.split("_")
    f.write(keys[0]+"\t"+keys[1]+"\t"+str(v)+"\n")
f.close()
# np.save(data + 'ComNb', comNb)
t2 = time.time()
print "Time of saving common_neighbors = ", (t2 - t1)
