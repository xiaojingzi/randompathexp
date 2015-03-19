import sys
import networkx as nx
import time
import re
import numpy as np
from collections import defaultdict

if len(sys.argv) < 9:
    print "usage: evaluateNeighborSim.py data T D displayEpsilon algorithm topK testParameterName datasize"
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
########### load ground truth of common neighbors ##########
groundTruthFile = "data/"+dataset +'.cn'
print "ground truth file = ", groundTruthFile
gt = {}
f = open(groundTruthFile, 'r')
for line in f:
    elements = line.strip().split("\t")
    k = elements[0]+"_" +elements[1]
    v = int(elements[2])
    gt[k] =v
f.close()

cnlistFile = "data/"+dataset +'.cnlist'
print "cn file = ", cnlistFile
cnlist = []
f = open(cnlistFile, 'r')
for line in f:
    cnlist.append(int(line.strip()))
f.close()
print "load cn file over"
############### load calculated top similar nodes ############
topKFile = "result/"+dataset + "_"+str(T)+"_"+str(D)+"_"+str(displayEpsilon)+".topK." + algorithm
print "topfile=", topKFile
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

print "load over"
#################  get the random top K similar nodes ##########
randomTopK = np.random.randint(N, size = (N, K))

#################  calculate average score of difference between algorithm A and random algorithm ##########
randomScore = 0
trueScore = 0
for i in cnlist:
    l = min(K,len(topKList[i]))
    if l == 0:
        continue
    currentTrueScore = 0
    currentRandomScore = 0
    for k in range(l):
        trueID = topKList[i][k]
        trueKey = str(i)+"_"+str(trueID)
        cn = 0
        if trueKey in gt:
            cn = gt[trueKey]
        currentTrueScore += cn
        randomID = randomTopK[i][k]
        randomKey = str(i)+"_"+str(randomID)
        cn = 0
        if randomKey in gt:
            cn = gt[randomKey]
        currentRandomScore += cn
    trueScore += currentTrueScore
    randomScore += currentRandomScore

print "evaluate over"
# print trueScore, randomScore

################ output score ###############################
f = open("result/"+testParameterName+'.neighborsim.score', 'a')
f.write(str(datasize)+"\t"+str(float(trueScore - randomScore) / (len(cnlist) * K)) + '\n')
f.close()

