import sys
import numpy as np
import re
import pickle
import time
from collections import defaultdict
import math




def loadTwoNetwork(dataset1=None, dataset2=None, algorithm=None):
	global pathVecList1, pathVecList2, name2id1, name2id2, data1, data2,commonNames, suffix

	suffix = algorithm
	data1 = dataset1
	data2 = dataset2
	
	nodemapfile1 = "data/"+dataset1+".dict"
	name2id1 = {}
	
	for line in open(nodemapfile1, "r"):
		elements = line.strip().split("\t")
		name = elements[0]
		id = int(elements[1])
		name2id1[name] = id

	nodemapfile2 = "data/"+dataset2+".dict"
	name2id2 = {}
	
	for line in open(nodemapfile2, "r"):
		elements = line.strip().split("\t")
		name = elements[0]
		id = int(elements[1])
		name2id2[name] = id

	commonNames = set(name2id1.keys()).intersection(set(name2id2.keys()))


	pathVecFile1 = "result/"+dataset1+"."+suffix
	pathVecList1 = [[] for i in range(len(name2id1))]

	i = 0
	for line in open(pathVecFile1, "r"):
		elements = line.strip().split(" ")
		for element in elements:
			pathVecList1[i].append(float(element))
		i +=1

	pathVecFile2 = "result/"+dataset2+"."+suffix
	pathVecList2 = [[] for i in range(len(name2id2))]

	i = 0
	for line in open(pathVecFile2, "r"):
		elements = line.strip().split(" ")
		for element in elements:
			pathVecList2[i].append(float(element))
		i +=1




def queryAllPairs( K=100 ):
	global pathVecList1, pathVecList2, name2id1, name2id2,  data1, data2,commonNames, suffix
	recall = [0 for i in range(K+1)]

	start = time.time()
	j = 0
	for commonName in commonNames:
		id1 = name2id1[commonName]
		id2 = name2id2[commonName]
		structureSim = []
		vec1 = np.asarray(pathVecList1[id1])
		for commonName in commonNames:
			id_ = name2id2[commonName]
			vec2 = np.asarray(pathVecList2[id_])
			sim = math.exp(-np.linalg.norm(vec1-vec2, ord=2))
			structureSim.append((id_, sim))
		if j % 100 == 0:
			print j
		j +=1
		sortedstructureSimList = sorted(structureSim, key=lambda x:x[1], reverse=True)[:K]
		rank = 0
		for id_, sim in sortedstructureSimList:
			if rank == K:
				break
			if id_ == id2:
				# print commonName, rank
				break
			rank +=1
		recall[rank] +=1.0


	# print recall
	
	size = len(commonNames)
	recall[0] = recall[0] / size
	for k in range(1,K+1):
		recall[k] = recall[k-1] + recall[k] / size

	print recall
	outputFile = "result/"+data1+"_"+data2+"."+suffix+".recall"
	f = open(outputFile,"w")
	for k in range(K):
		f.write(str(recall[k])+"\n")
	f.close()

	end = time.time()
	print "Time of querying all pairs = ", (end - start)


if len(sys.argv) < 4:
	print "usage: calrecallCoauthor.py  dataset1  dataset2 algorthim K(top similar node number)"
	exit()
dataset1 = "kdd"
dataset2 = "icdm"
K = 10
if len(sys.argv) > 1:
	dataset1 = sys.argv[1]
if len(sys.argv) > 2: 
	dataset2 = sys.argv[2]
if len(sys.argv) > 3:
	algorithm = sys.argv[3]
if len(sys.argv) > 4:
	K = int(sys.argv[4])

print "data1=",dataset1, ", data2=",dataset2, ", algorithm=",algorithm, ", top K =", K

loadTwoNetwork(dataset1, dataset2, algorithm)
queryAllPairs(K)