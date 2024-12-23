#!/usr/bin/env python3

import os
import sys
sys.path.insert(0,"/usr/local/lib/xpython_modules_corey/")
from corey_python_functions import *
from xvariables import *
from xdefinitions import *

COMMANDS=[]

################################################################
with open("README_PROTO_TRUNCATED.TXT","r") as fin:
    lines=fin.read().splitlines()

for line in lines:
    parts=line.split()
    lparts=parts[0].split(".")
    label=lparts[0]
    dec=lparts[1]
    compound=parts[4]
    if "-" in compound:
        print(compound)
        compound=compound.split("-")[1]
        print(compound)
        #sys.exit()
    print("label="+str(label),"dec="+str(dec),"compound="+str(compound))
    elements=getElements(compound)
    elements=sorted(elements)   #KC+TL20231009
    proto=label+":"+":".join(elements)  #you don't need dec here, it's always alphabetical
    #aflow-dev --proto=AB12C_cP14_221_a_h_b-001 | aflow-dev --get_cation_coordination_numbers --screen_only --print=json --quiet
    command="if [ ! -s "+DIR_OXIDATIONS+"/oxidations_"+proto+".json ]; then mkdir -p "+DIR_OXIDATIONS+"; "+AFLOW_BIN+" --proto="+proto+" | "+AFLOW_BIN+" --get_oxidation_number --screen_only --print=json --quiet > "+DIR_OXIDATIONS+"/oxidations_"+proto+".json; fi"
    COMMANDS.append(command)

with open("commands_oxidations.sh","w") as fout:
    fout.write("\n".join(COMMANDS))
