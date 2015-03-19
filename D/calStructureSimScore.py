import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict

if len(sys.argv) < 7:
    print "usage: evaluateStructureSim.py data T D Epsilon algorithm topK testParameterName"
    exit()
dataset = sys.argv[1]
T = int(sys.argv[2])
D = int(sys.argv[3])
Epsilon = float(sys.argv[4])
displayEpsilon = int(Epsilon*10000)
algorithm = sys.argv[5]
K = int(sys.argv[6])
testParameterName=sys.argv[7]

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

N = len(topKList)
# print "data size = ", N
#################  get the random top K similar nodes ##########
randomTopK = np.random.randint(N, size = (N, K))

#################  calculate average score of difference between algorithm A and random algorithm ##########
randomScore = defaultdict(float)
trueScore = defaultdict(float)
score = [0 for i in range(K)]
for gt in groundTruth:
    l = len(topKList[gt])
    currentTrueScore = defaultdict(float)
    currentRandomScore = defaultdict(float)
    for k in range(1,K+1):
        if l < k :
            currentTrueScore[k] = currentTrueScore[k-1] + 0
            currentRandomScore[k] = currentRandomScore[k-1] + 0
        else:
            trueID = topKList[gt][k-1]
            score = 0
            if trueID in groundTruthSet:
                score = 1
            currentTrueScore[k] = currentTrueScore[k-1] + score
            randomID = randomTopK[gt][k-1]
            score = 0
            if randomID in groundTruthSet:
                score = 1
            currentRandomScore[k] = currentRandomScore[k-1] + score
    for k in range(K):
        trueScore[k] += currentTrueScore[k+1]
        randomScore[k] += currentRandomScore[k+1]

################# output score ###############################
f = open("result/"+dataset + '_'+testParameterName+'.structuresim.score', 'a')
for k in range(K):
    # print count[k]
    f.write(str((trueScore[k] - randomScore[k]) / (len(groundTruth) * (k+1))) + '\n')
f.close()

