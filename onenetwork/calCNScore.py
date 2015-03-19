import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict

if len(sys.argv) < 5:
    print "usage: evaluateStructureSim.py data algorithm"
    exit()
data = sys.argv[1]
alg = sys.argv[4]

K = 50
print 'done with ' + data + ' with top '+ str(K) 
########### load ground truth ##########
groundTruthFile = dataset +'.structuregroundtruth'
groundTruth = []
f = open(groundTruthFile, 'r')
for line in f:
    groundTruth.append(int(line.strip()))
f.close()
groundTruthSet = set(groundTruth)

############### load calculated top similar nodes ############
topKsimFile = data + ".topKsim." + algorithm
f = open(topKsimFile, 'r')
topKsimList = defaultdict([])
currentNode = 0
for line in f:
    elements = line.strip().split(' ')
    for element in elements:
        topKsimList[currentNode].append(int(element.strip()))
    currentNode += 1
f.close()

#################  get the random top K similar nodes ##########
randomTopK = np.random.randint(numNodes, size = (numNodes, K))

#################  calculate average score of difference between algorithm A and random algorithm ##########
randomScore = defaultdict(float)
trueScore = defaultdict(float)
count = defaultdict(int)
score = [0 for i in range(K)]
for gt in groundTruth:
    cnt += 1

    l = len(topKsimList[gt])
    for k in range(l):
        count[k] +=1
        trueID = topKsimList[gt][k]
        if trueID in groundTruthSet:
            trueScore[k] += 1.0
        randomID = randomTopK[gt][k]
        if randomID in groundTruthSet:
           randomScore[k] += 1.0

################# output score ###############################
f = open(data + '.structuresim'+'.score', 'w')
for k in K:
    f.write(str((trueScore[k] - randomScore[k]) / (count[k] * k)) + '\n')
f.close()

