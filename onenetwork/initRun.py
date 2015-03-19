import sys
if len(sys.argv) < 3:
    print "usage: data simtype"
    exit()
data = sys.argv[1]
simtype = sys.argv[2]


print 'start to deal with ' + data
f = open("result/"+data + '_T.'+simtype+'.score', 'w')
f.close()
