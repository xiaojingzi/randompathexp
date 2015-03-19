import sys
import networkx as nx
import time
import re
import numpy as np
import pickle

if len(sys.argv) < 4:
    print "usage: python getTopK.py dataset K algorithm"
    exit()
dataset = sys.argv[1]
K = int(sys.argv[2])
algorithm = sys.argv[3]

############## load similarity matrix ####################
f = open("result/" + dataset+'.topK.'+algorithm ,'w')

for line in open("result/"+dataset+"."+algorithm, "r"):
    elements = re.compile("\s+").split(line.strip())
    alist = []
    i = 0
    for element in elements:
        alist.append((i, float(element)))
        i +=1
    sortedList = sorted(alist, key=lambda x:x[1], reverse=True)[:K]
    for k,v in sortedList:
        f.write(str(k)+" ")
    f.write("\n")
f.close()
