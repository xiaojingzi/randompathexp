import sys
import numpy as np
import re
import pickle
import time
from collections import defaultdict
import math
import networkx as nx


if len(sys.argv) < 4:
	print "usage: calculaterecall.py dataset algorithm P(only evaluate top-P nodes with highest degree)"
	exit()

dataset = sys.argv[1]
algorithm = sys.argv[2]
P = int(sys.argv[3])
print "dataset=", dataset,", algorithm=", algorithm
networkFile = "data/"+dataset+"1.graph"
topSimFile = "result/"+dataset+"1_"+dataset+"2."+algorithm

topNodes = []
G = nx.Graph()
id = 0
t1 = time.time()
for line in open(topSimFile, "r"):
	ids = line.strip().split(" ")
	topNodes.append(ids)	
	G.add_node(id)
	id +=1
t2 = time.time()
print "loading top-k node matrix = ", (t2-t1)
N = len(topNodes)
K = len(topNodes[0])

# read network
i = 0
for line in open(networkFile, "r"):
	if i == 0:
		i +=1
		continue
	i +=1
	elements = re.compile('\t').split(line.strip())
	node1 = int(elements[0])
	node2 = int(elements[1])
	w = float(elements[2])
	G.add_edge(node1, node2, weight=w)
	G.add_edge(node2, node1, weight=w) 

degrees = []
for n in G:
    degrees.append((n, G.degree(n)))
degreeList = sorted(degrees, key=lambda x:x[1], reverse=True)[:P]



recall = [0 for i in range(K+1)]
for id1, degree in degreeList:
	id2s = topNodes[id1]
	rank = 0 
	for id2 in id2s:
		if int(id2.strip()) == id1:
			# print commonName, rank
			break
		rank +=1
	recall[rank] +=1.0
	

recall[0] = recall[0] / P
for k in range(1,K):
	recall[k] = recall[k-1] + recall[k] / P
print recall

outputFile = dataset+"."+algorithm+".recall"
f = open(outputFile,"w")
for k in range(K):
	f.write(str(recall[k])+"\n")
f.close()



