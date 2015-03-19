import sys
import time
import networkx as  nx
import numpy as np 
import re
import math


def constraint(G, node_num): 
    """
    Burt's constraint measure (equation 2.4, page 55 of Burt, 1992). Essentially a
    measure of the extent to which v is invested in people who are invested in
    other of v's alters (neighbors).  The "constraint" is characterized
    by a lack of primary holes around each neighbor.  Formally:
    constraint(v) = sum_{w in MP(v), w != v} localConstraint(v,w)
    where MP(v) is the subset of v's neighbors that are both predecessors and successors of v.
    """

    scores = [0 for i in range(node_num)]
    
    for v in G:
        if G.degree(v) <=1: 
            continue
        v_neighbors = G.neighbors(v)  # v -> n1
        # print "v=", v, ", v'neighbors = ", v_neighbors
        for n1 in v_neighbors:
            n1_neighbors = G.neighbors(n1) # n1 -> n2
            # print "n1 = ", n1, ", n1's neighbors = ", n1_neighbors
            for n2 in n1_neighbors:
                # print "n2 = ", n2
                if n2 == v:
                    continue
                if n2 in v_neighbors:  # v -> n2
                    # print "++"
                    scores[n2] += G[v][n1]['weight'] * G[n1][n2]['weight'] 


    for v in G:
        degree = G.degree(v)
        if degree > 1:
            scores[v] = math.pow(scores[v], 2) / degree
        else:
            scores[v] = 1.0
    return scores

    


if len(sys.argv) < 4:
    print "usage: python caldegreestructurehole.py dataset t1(15) t2(0.0001) (t1 is the threshold of degree, t2 is the threshold of network constraint"
    exit()

dataset = sys.argv[1]
t1 = int(sys.argv[2])
t2 = float(sys.argv[3])
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

print "Calculate constaint values"
# constraint = constraint(G, node_num)
clustercoefficient = nx.clustering(G)
constraint = [ 0  for i in range(node_num)]
print type(clustercoefficient)
for k, v in clustercoefficient.items():
    constraint[k] = v
print "Calculate degree values "
# degrees = nx.degree_centrality(G).items()
degrees = []
for n in G:
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
print "100% degrees = ", np.percentile(degrees, 100)


outputfile = "data/"+dataset+".spcid"
f = open(outputfile, "w")
for n in G:
    if degrees[n] > t1 and constraint[n] < t2:
        print degrees[n],constraint[n]
        f.write(str(n)+ "\n")
f.close()


