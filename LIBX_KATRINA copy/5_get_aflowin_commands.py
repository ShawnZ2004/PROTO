#!/usr/bin/env python3

import os,json,sys
sys.path.insert(0,"/usr/local/lib/xpython_modules_corey/")
from xvariables import *

if PERFORM_COMPARISONS==False:
    sys.exit("the script is not needed, run its predecessor to get commands_aflowin.sh")

JSONS=[ DIR_COMPARISONS+"/"+j for j in os.listdir(DIR_COMPARISONS) ]

COMMANDS_AFLOWIN=[]

for j in JSONS:
    with open(j,"r") as fin:
        data=json.load(fin)
    if len(data)!=1:
        sys.exit("unexpected json structure")
    data=data[0]
    already_run=False
    for d in data["structures_duplicate"]:
        if d["enthalpy_formation_atom"] is not None:
            already_run=True
            break
        #[GOOD CHECK]else:
        #[GOOD CHECK]    print(d)
    if already_run: continue
    proto_dec=j.replace(DIR_COMPARISONS+"/compare_","").replace(".json","")
    #print(proto_dec)
    #oxides MIGHT be magnetic, better to run static to get fully resolved properties
    aflowin_command=AFLOW_BIN+" --aflow_proto="+proto_dec+" --run_relax_static --bader --noldau2 --keep_nomix --aflow_db_run --force"
    print(aflowin_command)
    COMMANDS_AFLOWIN.append(aflowin_command)

with open("commands_aflowin.sh","w") as fout:
    fout.write("\n".join(COMMANDS_AFLOWIN))
