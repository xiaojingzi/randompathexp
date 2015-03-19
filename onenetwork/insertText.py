import sys
if len(sys.argv) < 4:
    print "usage: data T simtype"
    exit()
data = sys.argv[1]
simtype = sys.argv[3]
T = sys.argv[2]
f = open("result/"+data + '_T.'+simtype+".score", 'a')
f.write('#'+T + '\n')
f.close()

