import sys
import re
import numpy as np
import numpy.linalg as linalg
import networkx as nx

if len(sys.argv) < 2:
    print "usage: python basicfeatures2.py dataset"
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
        continue;
    i +=1
    elements = re.compile('\t').split(line.strip())
    node1 = int(elements[0])
    node2 = int(elements[1])
    w = float(elements[2])
    G.add_edge(node1, node2, weight=w)
    G.node[node1]['weight'] += w
    G.add_edge(node2, node1, weight=w)
    G.node[node2]['weight'] += w
 
print "degree"
degree = nx.degree_centrality(G)
print "clustering coefficient"
clustering = nx.clustering(G)
print "closeness"
closeness = nx.closeness_centrality(G)
print "pagerank"
pagerank = nx.pagerank(G)
print "triangle number"
triangle_number = nx.clustering(G)
print "betweenness"
betweenness = nx.betweenness_centrality(G,k=1000)


# v /= v.sum(axis = 0)

print "===output==="
output_file = "result/"+dataset+".Degree"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(degree[i])+"\n")
f.close()

output_file = "result/"+dataset+".ClusCoeff"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(clustering[i])+"\n")
f.close()

output_file = "result/"+dataset+".Closeness"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(closeness[i])+"\n")
f.close()

output_file = "result/"+dataset+".Pagerank"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(pagerank[i])+"\n")
f.close()

output_file = "result/"+dataset+".Trianglenumber"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(triangle_number[i])+"\n")
f.close()

output_file = "result/"+dataset+".Betweenness"
f = open(output_file, "w")
for i in range(node_num):
    f.write(str(betweenness[i])+"\n")
f.close()
