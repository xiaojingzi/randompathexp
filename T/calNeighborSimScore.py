import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict

if len(sys.argv) < 8:
    print "usage: evaluateNeighborSim.py data T D Epsilon algorithm topK testParameterName"
    exit()



dataset = sys.argv[1]
T = int(sys.argv[2])
D = int(sys.argv[3])
Epsilon = float(sys.argv[4])
displayEpsilon=int(Epsilon*10000)
algorithm = sys.argv[5]
K = int(sys.argv[6])
testParameterName=sys.argv[7]
########### load ground truth of common neighbors ##########
groundTruthFile = "result/"+dataset +'.cn'
gt = {}
f = open(groundTruthFile, 'r')
for line in f:
    elements = line.strip().split("\t")
    k = elements[0]+"_" +elements[1]
    v = int(elements[2])
    gt[k] =v
f.close()

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

#################  get the random top K similar nodes ##########
randomTopK = np.random.randint(N, size = (N, K))

#################  calculate average score of difference between algorithm A and random algorithm ##########
randomScore = defaultdict(float)
trueScore = defaultdict(float)
score = [0 for i in range(K)]
for i in range(N):
    l = len(topKList[i])
    currentTrueScore = defaultdict(float)
    currentRandomScore = defaultdict(float)
    for k in range(1,K+1):
        if l < k :
            currentTrueScore[k] = currentTrueScore[k-1] + 0
            currentRandomScore[k] = currentRandomScore[k-1] + 0
        else:
            trueID = topKList[i][k-1]
            trueKey = str(i)+"_"+str(trueID)
            cn = 0
            if trueKey in gt:
                cn = gt[trueKey]
            currentTrueScore[k] = currentTrueScore[k-1] + cn
            randomID = randomTopK[i][k-1]
            randomKey = str(i)+"_"+str(randomID)
            cn = 0
            if randomKey in gt:
                cn = gt[randomKey]
            currentRandomScore[k] = currentRandomScore[k-1] + cn
    for k in range(K):
        trueScore[k] += currentTrueScore[k+1]
        randomScore[k] += currentRandomScore[k+1]


################ output score ###############################
f = open("result/"+dataset + '_'+testParameterName+'.neighborsim.score', 'a')
for k in range(K):
    f.write(str((trueScore[k] - randomScore[k]) / (N * (k+1))) + '\n')
f.close()

