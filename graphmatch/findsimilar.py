import sys
import numpy as np
import re
import pickle
import time
from collections import defaultdict
import math


def loadData(dataset=None, algorithm=None):
	start = time.time()
	# similarPath = pickle.load(open(pathSimFile, "rb"))
	# model = np.load(pathVecFile)
	# similarStructure = model['sim']
	global name2id, id2name, N, pathVecList, pathSimList
	
	nodemapfile = "data/"+dataset+".dict"
	
	name2id = {}
	id2name = {}
	i = 0
	for line in open(nodemapfile, "r"):
		elements = line.strip().split("\t")
		if len(elements) == 2:
			name = elements[0]
			id = int(elements[1])
		else:
			name = "no name"
			id = int(elements[0])
		name2id[name] = id
		id2name[id] = name
		
	N = len(id2name)
	print "node number = ", N

	pathVecFile = "result/"+dataset+"."+algorithm
	pathVecList = [[] for i in range(N)]
	
	i = 0
	for line in open(pathVecFile, "r"):
		elements = line.strip().split(" ")
		for element in elements:
			pathVecList[i].append(float(element))
		i +=1
	
	end = time.time()
	# print pathVecList
	print "Loading time = ", (end - start)

def queryTopNodesWithStructureSim( username=None, K=10):
	global pathVecList, name2id, id2name, N
	structureSim = []
	username = username.replace("_", " ")
	userid = name2id[username]
	vec1 = np.asarray(pathVecList[userid])
	
	for i in range(N):
		vec2 = np.asarray(pathVecList[i])
		# print vec2
		# sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1, ord=2) * np.linalg.norm(vec2,ord=2))
		sim = math.exp(-np.linalg.norm(vec1-vec2, ord=2))
		structureSim.append((id2name[i], sim))
		
	sortedstructureSimList = sorted(structureSim, key=lambda x:x[1], reverse=True)[:K]
	print "~~~Top ", K, " nodes with highest structureSim:", 
	print sortedstructureSimList


def queryAllPairs( K=10 ):
	global pathVecList, name2id, id2name, N
	start = time.time()
	for j in range(N):
		structureSim = []
		vec1 = np.asarray(pathVecList[j])
		for i in range(N):
			vec2 = np.asarray(pathVecList[i])
			# sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1, ord=2) * np.linalg.norm(vec2,ord=2))
			sim = math.exp(-np.linalg.norm(vec1-vec2, ord=2))
			structureSim.append((id2name[i], sim))
			
		if j % 100 == 0:
			print j
		sortedstructureSimList = sorted(structureSim, key=lambda x:x[1], reverse=True)[:K]
	end = time.time()
	print "Time of querying all pairs = ", (end - start)


if len(sys.argv) < 5:
	print "usage: findsimilar.py  dataset algorithm K username "
	exit()
dataset = "data-ma"
userid = 0
K = 10
if len(sys.argv) > 1:
	dataset = sys.argv[1]
if len(sys.argv) > 2:
	algorithm = sys.argv[2]
if len(sys.argv) > 3:
	K = int(sys.argv[3])
if len(sys.argv) > 4:
	username = sys.argv[4]

print dataset, algorithm, K, username

loadData(dataset=dataset, algorithm=algorithm)
queryTopNodesWithStructureSim(username=username, K=K )
# queryAllPairs(K)
