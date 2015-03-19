import sys
import numpy as np
import re
import pickle
import time
from collections import defaultdict
import math
import random




def loadTwoNetwork(dataset1=None, dataset2=None):
	global name2id1, name2id2, data1, data2,commonNames

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


def queryAllPairs( K=100 ):
	global name2id1, name2id2,  data1, data2, commonNames
	recall = [0 for i in range(K+1)]

	start = time.time()
	j = 0
	for commonName in commonNames:
		if j % 100 == 0:
			print j
		j +=1
		randomList = random.sample(commonNames, K)
		rank = 0
		for name_ in randomList:
			if rank == K:
				break
			if name_ == commonName:
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
	outputFile = "result/"+data1+"_"+data2+".Random.recall"
	f = open(outputFile,"w")
	for k in range(K):
		f.write(str(recall[k])+"\n")
	f.close()

	end = time.time()
	print "Time of querying all pairs = ", (end - start)


if len(sys.argv) < 3:
	print "usage: calrecallRandom.py  dataset1  dataset2 K"
	exit()
dataset1 = "kdd"
dataset2 = "icdm"
K = 10
if len(sys.argv) > 1:
	dataset1 = sys.argv[1]
if len(sys.argv) > 2: 
	dataset2 = sys.argv[2]
if len(sys.argv) > 3:
	K = int(sys.argv[3])

print dataset1, dataset2, K

loadTwoNetwork(dataset1, dataset2)
queryAllPairs(K)
