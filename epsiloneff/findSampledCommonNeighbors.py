import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict
import random


if len(sys.argv) < 3:
    print "usage: dataset samplesize"
    print "output: dataset + '.cn'"
    exit()

dataset=sys.argv[1]
samplesize=int(sys.argv[2])
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

if N > samplesize:
    samplelist = random.sample(range(N), samplesize)
else:
    samplelist = range(N)

t1 = time.time()
comNb = defaultdict(int)
for u in samplelist:
    for v in G.neighbors(u):
        for w in G.neighbors(v):
            comNb[str(u)+"_"+str(w)] +=1

print "pairs have cn = ",len(comNb)

outputFile = "data/"+dataset+".cn"
f = open(outputFile, "w", 10000)
for k, v in comNb.items():
    keys = k.split("_")
    f.write(keys[0]+"\t"+keys[1]+"\t"+str(v)+"\n")
    f.flush
f.close()


# comNb = defaultdict(int)
# for u in G.nodes():
#     for v in G.nodes():
#         cn = sorted(nx.common_neighbors(G, u, v) )
#         if len(cn) >  0:
#             comNb[str(u)+"_"+str(v)] +=len(cn)

# outputFile = "data/"+dataset+".cntemp"
# f = open(outputFile, "w")
# for k, v in comNb.items():
#     keys = k.split("_")
#     f.write(keys[0]+"\t"+keys[1]+"\t"+str(v)+"\n")
# f.close()

outputFile = "data/"+dataset+".cnlist"
f = open(outputFile, "w")
for id in samplelist:
    f.write(str(id)+"\n")
f.close()


t2 = time.time()
print "Time of saving common_neighbors = ", (t2 - t1)
