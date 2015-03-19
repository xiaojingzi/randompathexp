import sys
import numpy as np
import re
import time
import math


def loadData(dataset=None):
	start = time.time()
	global pathVecList, name2id, data
	data = dataset
	
	nodemapfile = "data/"+dataset+".dict"
	name2id = {}
	
	for line in open(nodemapfile, "r"):
		elements = line.strip().split("\t")
		name = elements[0]
		id = int(elements[1])
		name2id[name] = id

	pathVecFile = "result/"+dataset+".pathvec"
	pathVecList = [[] for i in range(len(name2id))]

	i = 0
	for line in open(pathVecFile, "r"):
		elements = line.strip().split(" ")
		for element in elements:
			pathVecList[i].append(int(element))
		i +=1
	end = time.time()
	print "Loading time = ", (end - start)


def sigmod(x=None):
	print x
	return 2/(1+math.exp(-x))

def calculateStructureSimilarity( username=None):
	global pathVecList, name2id, data
	N = len(name2id)
	structureDistance = []
	simFile = username
	username = username.replace("_", " ")
	userid = name2id[username]
	vec1 = np.asarray(pathVecList[userid])
	print vec1
	f = open("result/"+simFile, "w")
	
	maxdistance = 0.0

	for i in range(N):
		vec2 = np.asarray(pathVecList[i])
		distance = np.linalg.norm(vec1-vec2, ord=2)
		if maxdistance < distance:
			maxdistance = distance
		structureDistance.append(distance)
		
	print "maxdistance=", maxdistance
	for distance in structureDistance:
		f.write(str(float(distance)/maxdistance)+"\n")
	f.close()


if len(sys.argv) < 3:
	print "usage: calculatesimilarity.py  dataset username "
	exit()
dataset = sys.argv[1]
username = sys.argv[2]

print dataset, username

loadData(dataset=dataset)
calculateStructureSimilarity(username=username)
