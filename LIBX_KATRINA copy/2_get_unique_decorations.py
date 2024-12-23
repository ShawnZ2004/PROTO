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
    print("label="+str(label))
    #aflow-dev --proto=AB12C_cP14_221_a_h_b-001 | aflow-dev --unique_atom_decorations --screen_only --print=json --quiet
    command="if [ ! -s "+DIR_DECORATIONS+"/decorations_"+label+".json ]; then mkdir -p "+DIR_DECORATIONS+"; "+AFLOW_BIN+" --proto="+label+" | "+AFLOW_BIN+" --unique_atom_decorations --screen_only --print=json --quiet > "+DIR_DECORATIONS+"/decorations_"+label+".json; fi"
    COMMANDS.append(command)

with open("commands_decorations.sh","w") as fout:
    fout.write("\n".join(COMMANDS))
