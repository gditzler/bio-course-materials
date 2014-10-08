#!/usr/bin/env python
import sys

f = open(sys.argv[1], "rU")
total=0
for line in f :
    if(line[0] == ">") :
        if 'each' in locals():
            print("%s\t%i" % (h,each))
        h=line[1:-1]
#        h=h.split()[0]
        each=0
        continue
    each += len(line)-1
    total += len(line)-1
print("%s\t%i" % (h,each))
# print("Total %i nucleotides" % (total))
