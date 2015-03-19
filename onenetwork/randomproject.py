
import sys
import re
import numpy as np
import numpy.linalg as linalg
import networkx as nx
import time
import pickle 
import random
from collections import Counter
from collections import defaultdict
import math




def init( data_param=None, T_param=10, K_param=100, epsilon_param=0.01):
	global T, R, K, dataset, network_file, map_file
	T = T_param # path length
	K = K_param 
	dataset = data_param
	epsilon = epsilon_param

	network_file= "data/"+dataset+".graph"
	map_file= "data/"+dataset+".dict"

	c = 0.5
	delta=0.1
	C_T_2 = T*(T-1)/2
	R =  int(c/(epsilon*epsilon) * (math.log(C_T_2, 2) +1+ math.log(1/delta))) # sample size

	print "sample size = ", R


def loadNode():
	global G, map_file, N
	G = nx.Graph()
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
	N = G.number_of_nodes()


def loadLink():
	global G, newtork_file, N, M
	e = 0

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
		print node1, node2, w
		G.add_edge(node1, node2, weight=w)
		G.node[node1]['weight'] += w
		G.add_edge(node2, node1, weight=w)
		G.node[node2]['weight'] += w
		e +=1
		if e % 10000 == 0:
			print e, " edge"

	print "load edges dones."
	M = G.number_of_edges()
	print "#nodes = ", N, ", #edges = ", M

def generateRandomPath():
	start = time.time()
	global N, R, T, G, path2node, node2path
	path2node = []
	node2path = [[] for i in range(N)]

	# R = 1
	# randid = 0
	
	r = 0
	while r < R:
	# for i in G.nodes():
		# current_node = i
		current_node = random.randint(0, N-1)
		# current_node = randid
		neighbors = G.neighbors(current_node)
		# print "neighbors"
		# print neighbors
		if len(neighbors) == 0:
			continue
		
		path = []
		path.append(current_node)
		
		t = 0
		while t < T:
			# print "current node = ", current_node
			neighbors = G.neighbors(current_node)
			if len(neighbors) == 0:
				break

			# neighbor_index = random.randint(0,len(neighbors)-1)
			

			totalWeight = G.node[current_node]['weight']
			# print "totolweight = ", totalWeight
			# print "neighbor number = ", len(neighbors)

			# print "original weight"
			# for j in range(len(neighbors)):
			# 	print str(G.edge[current_node][neighbors[j]]['weight'])+" ",
			# print ""
			
			w = []
			w.append(G.edge[current_node][neighbors[0]]['weight']/totalWeight)
			for j in range(len(neighbors)-1):
				cWeight = G.edge[current_node][neighbors[j+1]]['weight']/totalWeight
				w.append(w[j] + cWeight)

			# print "accumulated weight"
			# print w
			rand = random.random()
			# print "random value = ", rand
			for j in range(len(neighbors)):
				if w[j] >= rand:
					break
			neighbor_index = j

			current_node = neighbors[neighbor_index]
			path.append(current_node)
			node2path[current_node].append(r) 
			t+=1
			
		r+=1
		path2node.append(path)
		if r % 10000 == 0 :
			print r
	end = time.time()
	print "Time of generating random paths = ", (end - start)
	
def calculatePathSim():
	start = time.time()
	global pathVecList, pathSimList, node2path, path2node, N, K
	pathVecList = [[] for i in range(N)]
	pathSimList = [[] for i in range(N)]

	for i in range(N):
		# print "======",i,"====="
		nodesInSamePath = defaultdict(int)
		for path in node2path[i]:
			nodes = path2node[path]
			# print "nodes in path = ", nodes
			nodeset = set(nodes)
			for node in nodeset:
				nodesInSamePath[node] += 1
		sortedNodes = sorted(nodesInSamePath.items(), key=lambda x:x[1], reverse=True)

		dim = min(K, len(sortedNodes))
		pathSimList[i] = sortedNodes[:dim]

		k = 0
		for key,value in sortedNodes:
			pathVecList[i].append(value)
			k +=1
			if k >= K:
				break
		
		if len(pathVecList[i]) < K:
			for j in range(k, K):
				pathVecList[i].append(0)
		if i % 10000 == 0 :
			print i
	end = time.time()
	print "Time of calcualating pathsim = ", (end - start)

def savePathSim():
	global dataset, pathSimList
	start = time.time()
	pathSimFile = "result/"+dataset+".pathsim"
	f = open(pathSimFile, "wb")
	for i in range(N):
		for j in range(len(pathSimList[i])):
			f.write(str(pathSimList[i][j][0])+" ")
		f.write("\n")
		if i % 10000 == 0 :
			print i
	f.close()
	end = time.time()
	print "Total time used for saving pathSim = ", (end - start)



def savePathVec():
	global dataset, pathVecList
	start = time.time()
	pathVecFile = "result/"+dataset+".pathvec"
	f = open(pathVecFile, "wb")
	for i in range(N):
		for j in range(len(pathVecList[i])):
			f.write(str(pathVecList[i][j])+" ")
		f.write("\n")
		if i % 10000 == 0 :
			print i
	f.close()
	end = time.time()
	print "Total time used for saving pathVec = ", (end - start)

def precompute(data=None,  T=10, K=100,epsilon=0.01):
	init(data, T, K, epsilon)
	t1 = time.time()
	loadNode()
	t2 = time.time()
	print "Time of loading nodes = ", (t2 - t1)
	loadLink()
	t3 = time.time()
	print "Time of loading edges = ", (t3 - t2)
	generateRandomPath()
	t4 = time.time()
	print "Time of generating random paths = ", (t4 - t3)
	calculatePathSim()
	t5 = time.time()
	print "Time of calcualating pathsim = ", (t5 - t4)
	savePathSim()
	savePathVec()
	t6 = time.time()
	print "Time of saving = ", (t6 - t5)


if len(sys.argv) < 5:
	print "usage: dataset T K epsilon"
	exit()
data=sys.argv[1]
T = int(sys.argv[2])
K = int(sys.argv[3])
epsilon = float(sys.argv[4])

print "dataset=", data, ",path length =", T, ",path vector length=", K, ", epsilon=", epsilon
precompute(data=data,T=T, K=K,epsilon=epsilon)

