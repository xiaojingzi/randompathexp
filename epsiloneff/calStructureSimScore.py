import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict

if len(sys.argv) < 7:
    print "usage: evaluateStructureSim.py data T D Epsilon algorithm topK testParameterName datasize"
    exit()
dataset = sys.argv[1]
T = int(sys.argv[2])
D = int(sys.argv[3])
displayEpsilon = sys.argv[4]
print "displayEpsilon = ", displayEpsilon
algorithm = sys.argv[5]
K = int(sys.argv[6])
testParameterName=sys.argv[7]
datasize=int(sys.argv[8])
########### load ground truth ##########
groundTruthFile = "data/"+dataset +'.spcid'
groundTruth = []
f = open(groundTruthFile, 'r')
for line in f:
    groundTruth.append(int(line.strip()))
f.close()
groundTruthSet = set(groundTruth)

############### load calculated top similar nodes ############

topKFile = "result/"+dataset + "_"+str(T)+"_"+str(D)+"_"+str(displayEpsilon)+".topK." + algorithm
print topKFile
f = open(topKFile, 'r')
topKList = defaultdict(list)
currentNode = 0
for line in f:
    if line.strip() != "":
        elements = line.strip().split(' ')
        for element in elements:
            topKList[currentNode].append(int(element.strip()))
    else:
        topKList[currentNode] = []
    currentNode += 1
f.close()
print currentNode
N = len(topKList)
# print "data size = ", N
#################  get the random top K similar nodes ##########
print N, K
randomTopK = np.random.randint(N, size = (N, K))

#################  calculate average score of difference between algorithm A and random algorithm ##########
randomScore = 0
trueScore = 0
for gt in groundTruth:
    l = min(K,len(topKList[gt]))
    if l == 0:
        continue
    currentTrueScore = 0
    currentRandomScore = 0
    for k in range(l):
        trueID = topKList[gt][k]
        score = 0
        if trueID in groundTruthSet:
            score = 1
        currentTrueScore += score
        randomID = randomTopK[gt][k]
        score = 0
        if randomID in groundTruthSet:
            score = 1
        currentRandomScore += score
    trueScore += currentTrueScore
    randomScore += currentRandomScore

################# output score ###############################
f = open("result/"+testParameterName+'.structuresim.score', 'a')
f.write(str(datasize)+"\t"+str(float(trueScore - randomScore) / (len(groundTruth) * K)) + '\n')
print str(datasize), "\t", str(float(trueScore - randomScore) / (len(groundTruth) * K))
f.close()

