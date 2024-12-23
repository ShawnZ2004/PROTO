#!/usr/bin/env python3

import os
import sys
sys.path.insert(0,"/usr/local/lib/xpython_modules_corey/")
from corey_python_functions import *
from xvariables import *
from xdefinitions import *

#with open("uniq.out","r") as fin:
#    _PROTOS_ALREADY_CONSIDERED=fin.read().splitlines()
#PROTOS_ALREADY_CONSIDERED=[]
#for p in _PROTOS_ALREADY_CONSIDERED:
#    PROTOS_ALREADY_CONSIDERED.append(p.split(".")[0])

with open("README_PROTO.TXT","r") as fin:
    lines=fin.read().splitlines()

lines_new=[]
for line in lines:
    if not line: continue
    if line[0]!="A": continue
    parts=line.split()
    if len(parts)<5: continue
    lparts=parts[0].split(".")
    if len(lparts)<2: continue
    label=lparts[0]
#    if label in PROTOS_ALREADY_CONSIDERED:
#        print(label+" already considered!")
#        continue
    dec=lparts[1]
    print("label="+str(label),"dec="+str(dec))
    dec_good=True
    for i in range(NARY):
        let=chr(ord("A")+i)
        print(let)
        if let not in dec: 
            dec_good=False
            break
    print("dec_good[1]="+str(dec_good))
    let=chr(ord("A")+NARY+1)
    if let in dec: 
        dec_good=False
    print("dec_good[2]="+str(dec_good))
    if not dec_good: 
        continue
    #if "B" not in dec: continue
    #if "C" not in dec: continue
    #if "D" in dec: continue
    compound=parts[4]
    if "-" in compound:
        compound=compound.split("-")[1]
    #skip pocc for now
    if "[" in compound or "(" in compound:
        continue
        #print(compound)
        #sys.exit()
    print("label="+str(label),"dec="+str(dec),"compound="+str(compound))
    elements=getElements(compound)
    print("elements="+str(elements))
    if len(elements)!=NARY: continue
    if not any([ e in elements for e in ANIONS2SEARCH  ]): continue
    lines_new.append(line)

with open("README_PROTO_TRUNCATED.TXT","w") as fout:
    fout.write("\n".join(lines_new))
