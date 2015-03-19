import sys
import time
import networkx as  nx
import numpy as np 
import re

if len(sys.argv) < 2:
    print "usage: python refex.py dataset"
    exit()

dataset = sys.argv[1]
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
        i +=1
        continue
    i +=1
    elements = re.compile('\t').split(line.strip())
    node1 = int(elements[0])
    node2 = int(elements[1])
    w = float(elements[2])
    G.add_edge(node1, node2, weight=w)
    G.node[node1]['weight'] += w
    G.add_edge(node2, node1, weight=w) 
    G.node[node2]['weight'] += w

# print G.edges(data=True)


v = np.zeros((node_num, 6))
for i in range(node_num):
    if i % 5000 == 0:
        print "node :", i
    # print "======= Node ", i, "======"
    # local properties -------- degree / weighted_degree
    degree = G.degree(i)
    # print "degree", degree
    weightedDegree = G.degree(i, weight='weight')
    # print "weightedDegree", weightedDegree
   
    # egonetwork properties -------- edge number(weighted) in egonetwork, degree(weighted) outside egonetwork
    egonet = nx.ego_graph(G, i)
    egonet_in_edges = egonet.edges(data=True)
    # print "egonet_in_edges", egonet_in_edges
    egonet_in_edge_number = len(egonet_in_edges)
    # print "egonet_in_edge_number", egonet_in_edge_number
    egonet_in_weight = 0
    for edge in egonet_in_edges:
        egonet_in_weight += edge[2]['weight']

    # print "egonet_in_weight", egonet_in_weight
    egonet_degree_list = G.degree(egonet) 
    # print "egonet_degree_list", egonet_degree_list
    
    egonet_out_degree = 0
    # print egonet_degree_list
    for node, ndegree in egonet_degree_list.items():
        egonet_out_degree += ndegree
    # print "egonet_out_degree", egonet_out_degree
    egonet_out_degree -= egonet_in_edge_number*2

    # print "egonet_out_degree after removing egonet_in_edge_number", egonet_out_degree

    egonet_weight_list = G.degree(egonet, weight='weight')
    # print "egonet_weight_list", egonet_weight_list
    egonet_out_weight = 0
    for node, nweightdegree in egonet_weight_list.items():
        egonet_out_weight += nweightdegree
    # print "egonet_out_weight", egonet_out_weight
    egonet_out_weight -= egonet_in_weight*2
    # print "egonet_out_weight after removing egonet_out_weight", egonet_out_weight
    
    v[i][0] = degree
    v[i][1] = weightedDegree
    v[i][2] = egonet_in_edge_number
    v[i][3] = egonet_in_weight
    v[i][4] = egonet_out_degree
    v[i][5] = egonet_out_weight

v /= v.sum(axis = 0)

# print v

neighbors = G.neighbors(0)
meanvec = v[neighbors,:].mean(0)
sumvec = v[neighbors,:].sum(0)
newv = np.hstack((meanvec, sumvec))


for i in range(1,node_num):
    if i % 5000 == 0:
        print "neighbor feature: ", i

    if v[i][0] == 0:
        vec = np.zeros((1,12))
        newv = np.vstack((newv, vec))
    else:
        neighbors = G.neighbors(i)
        meanvec = v[neighbors,:].mean(0)
        sumvec = v[neighbors,:].sum(0)
        vec = np.hstack((meanvec, sumvec))
        newv = np.vstack((newv, vec))

v = np.hstack((v, newv))

# print v

output_file = "result/"+dataset+".ReFex"
f = open(output_file, "w")
d = v.shape[1]
for i in range(node_num):
    for j in range(d):
        f.write(str(v[i][j])+" ")
    f.write("\n")
f.close()

