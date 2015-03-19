import sys
if len(sys.argv) < 4:
    print "usage: data simtype testparameter"
    exit()
data = sys.argv[1]
simtype = sys.argv[2]
testparameter = sys.argv[3]

print 'start to deal with ' + data
f = open("result/"+data + '_'+testparameter+'.'+simtype+'.score', 'w')
f.close()
