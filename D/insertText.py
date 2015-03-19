import sys
if len(sys.argv) < 5:
    print "usage: data parametervalue simtype testparameter"
    exit()
data = sys.argv[1]
simtype = sys.argv[3]
parametervalue = sys.argv[2]
testparameter = sys.argv[4]
f = open("result/"+data + '_'+testparameter+'.'+simtype+".score", 'a')
f.write('#'+parametervalue + '\n')
f.close()

