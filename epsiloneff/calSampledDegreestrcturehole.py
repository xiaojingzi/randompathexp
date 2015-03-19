import sys
import time
import networkx as  nx
import numpy as np 
import re
import math
import random


def constraint(G,samplelist): 
    """
    Burt's constraint measure (equation 2.4, page 55 of Burt, 1992). Essentially a
    measure of the extent to which v is invested in people who are invested in
    other of v's alters (neighbors).  The "constraint" is characterized
    by a lack of primary holes around each neighbor.  Formally:
    constraint(v) = sum_{w in MP(v), w != v} localConstraint(v,w)
    where MP(v) is the subset of v's neighbors that are both predecessors and successors of v.
    """

    scores = [0 for i in range(len(samplelist))]
    
    i = 0
    for u in samplelist:
        if G.degree(u) <=1: 
            continue
        u_neighbors = G.neighbors(u)
        # print "u_neighbor:", u_neighbors
        for v in u_neighbors:
            # print "v_neighbor",G.neighbors(v)
            for w in G.neighbors(v):
                if w == u:
                    continue
                if w in u_neighbors:
                    # print "u=",u , "v=", v, "w=",w  
                    scores[i] +=  G[v][w]['weight'] * G[w][u]['weight'] 
        i+=1

    i = 0
    for u in samplelist:
        degree = G.degree(u)
        if degree > 1:
            scores[i] = math.pow(scores[i], 2) / degree
        else:
            scores[i] = 1.0
        i+=1
    return scores

    


if len(sys.argv) < 3:
    print "usage: python caldegreestructurehole.py dataset samplesize"
    exit()

dataset = sys.argv[1]
samplesize = int(sys.argv[2])
# t1 = int(sys.argv[2])
# t2 = float(sys.argv[3])
network_file = "data/"+dataset+".graph"
map_file = "data/"+dataset+".dict"



G = nx.Graph()

# read node map
for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    name = elements[0]
    id = int(elements[1])
    G.add_node(id, name=name,weight=0.0)

node_num = G.number_of_nodes()

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
    G.add_edge(node1, node2, weight=w)
    G.node[node1]['weight'] += w
    G.add_edge(node2, node1, weight=w) 
    G.node[node2]['weight'] += w

for edge in G.edges(data=True):
    edge[2]['weight']  = edge[2]['weight'] / G.node[edge[0]]['weight']


samplelist = random.sample(range(node_num), samplesize)

print "Calculate constaint values"
constraint = constraint(G, samplelist)
# clustercoefficient = nx.clustering(G)
# constraint = [ 0  for i in range(node_num)]
# print type(clustercoefficient)
# for k, v in clustercoefficient.items():
#     constraint[k] = v
print "Calculate degree values "
# degrees = nx.degree_centrality(G).items()
degrees = []
for n in samplelist:
    degrees.append( G.degree(n))

print "10% contraint = ", np.percentile(constraint, 10)
print "20% contraint = ", np.percentile(constraint, 20)
print "30% contraint = ", np.percentile(constraint, 30)
print "40% contraint = ", np.percentile(constraint, 40)
print "50% contraint = ", np.percentile(constraint, 50)
print "60% contraint = ", np.percentile(constraint, 60)
print "70% contraint = ", np.percentile(constraint, 70)
print "80% contraint = ", np.percentile(constraint, 80)
print "90% contraint = ", np.percentile(constraint, 90)
print "100% contraint = ", np.percentile(constraint, 100)

print "10% degrees = ", np.percentile(degrees, 10)
print "20% degrees = ", np.percentile(degrees, 20)
print "30% degrees = ", np.percentile(degrees, 30)
print "40% degrees = ", np.percentile(degrees, 40)
print "50% degrees = ", np.percentile(degrees, 50)
print "60% degrees = ", np.percentile(degrees, 60)
print "70% degrees = ", np.percentile(degrees, 70)
print "80% degrees = ", np.percentile(degrees, 80)
print "90% degrees = ", np.percentile(degrees, 90)
print "91% degrees = ", np.percentile(degrees, 91)
print "92% degrees = ", np.percentile(degrees, 92)
print "93% degrees = ", np.percentile(degrees, 93)
print "94% degrees = ", np.percentile(degrees, 94)
print "95% degrees = ", np.percentile(degrees, 95)
print "96% degrees = ", np.percentile(degrees, 96)
print "97% degrees = ", np.percentile(degrees, 97)
print "98% degrees = ", np.percentile(degrees, 98)
print "99% degrees = ", np.percentile(degrees, 99)
print "100% degrees = ", np.percentile(degrees, 100)

t1 = np.percentile(degrees, 90)
t2 = np.percentile(constraint, 20)

outputfile = "data/"+dataset+".spcid"
f = open(outputfile, "w", 10000)
for i in range(len(samplelist)):
    if degrees[i] >= t1 and constraint[i] <= t2:
        print degrees[i],constraint[i]
        f.write(str(samplelist[i])+ "\n")
f.close()


