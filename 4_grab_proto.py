#!/usr/bin/env python3

import os
import sys
sys.path.insert(0,"/usr/local/lib/xpython_modules_corey/")
from corey_python_functions import *
from corey_VASP_crystal_functions import structure
from xvariables import *
from xdefinitions import *
import itertools
import copy
import json

debug_flag_set=False
if len(sys.argv)>1: #List of file and command line arguments 
    for arg in sys.argv:
        if arg=="--debug":
            debug_flag_set=True

DEBUG=(False or debug_flag_set)
TEST_PROTOS=(False or DEBUG)

COMMANDS=[]
################################################################
with open("README_PROTO_TRUNCATED.TXT","r") as fin: #Read prototype for NARY
    lines=fin.read().splitlines()

for line in lines:
    print("*"*30)
    parts=line.split()
    lparts=parts[0].split(".")
    label=lparts[0] #Prototype label
    dec=lparts[1] #Decoration elements
    compound=parts[4] #Natural compound stable 
    if "-" in compound:
        if DEBUG: print(compound)
        compound=compound.split("-")[1]
        if DEBUG: print(compound)
        #sys.exit()
    elements=getElements(compound) #List fo natural compound elements 
    elements_sorted = sorted(elements)   #KC+TL20231009
    print("label="+str(label),"dec="+str(dec),"compound="+str(compound),"elements="+str(elements))
    if len(elements)!=NARY: sys.exit("odd "+str(NARY)+"-nary prototype: "+label) #Print error


    #get oxidation states for this proto+chemistry, continue if not calculable or multi-valent
    proto=label+":"+":".join(elements)  #you don't need dec here, it's always alphabetical
    if DEBUG: print(proto)

    
    if DEBUG: print("trying to open file="+str(DIR_OXIDATIONS+"/oxidations_"+proto+".json"))
    #Get oxidation states for the prototype 
    try: 
        with open(DIR_OXIDATIONS+"/oxidations_"+proto+".json","r") as fin: 
            data=json.load(fin) #Creates a json object (dictionary)
    except:
        print("skipping proto="+str(proto)+", no json file found\n")
        continue
    #Keys are element and values are oxidation states as a dictionary 
    #{'Cl': {'0': 5, '1': 5}, 'K': {'2': 1, '3': 1}, 'O': {'4': -2, '5': -2, '6': -2, '7': -2, '8': -2, '9': -2}}
    if not data.keys(): #The element has no oxidation states -- Positions 
        print("skipping proto="+str(proto)+", no oxidation data available\n")
    has_consistent_oxstates=True
    old_oxstates=[]
    if DEBUG: print(data.keys())
    if DEBUG: print(data)
    for e in data.keys(): #Iterate through each element
        ox0=-999
        for ia in data[e].keys(): #For each oxidation string in the values("0")
            if ox0==-999:   #choose first oxdiation
                ox0=data[e][ia] #Integer value of the oxidation state (5)
                old_oxstates.append(ox0)
            ox=data[e][ia] #Changing oxidation
            #print("oxidation_state["+str(e)+"]="+str(ox))
            #Multi-valent element (multiple oxidation states) -- remove repetitive oxidation states 
            if ox0!=ox:
                #print("multi-valent proto="+str(proto)+", skipping")
                has_consistent_oxstates=False 
    if has_consistent_oxstates==False:
        print("multi-valent proto="+str(proto)+", skipping\n")
        continue
    #oxstates=[2,2,-2]   #REMOVE ME
    print("oxstates = ", old_oxstates)   
    
#----Oxdiation states of cations and anions 2 search 

    oxstates_cations2dec=[]
    for e in CATIONS2DECORATE:
        if e not in OXIDATION_STATES.keys():
            if DEBUG: print("no oxidation states available for e="+str(e)+", skipping")
            break
        #need to find ONLY positive oxstates for cation
        #here's a good example: aflow-dev --proto=ABC_hR3_160_a_a_a-001.BAC:Fe:Mn:O [decv=['BAC']; decv_dup=['BAC', 'BCA']] (new) vs. aflow-dev --proto=ABC_hR3_160_a_a_a-001.ABC:C:O:S (std)
        #this works if Fe is -2, we are not interested in this case
        #otherwise fix the test below
        #oxstates_possible_cations.append(OXIDATION_STATES[e])
        oxstates_cation=[ i for i in OXIDATION_STATES[e] if i>=0 ]
        if not oxstates_cation:
            if DEBUG: print("no positive oxidation states for e="+str(e)+", skipping")
            break
        oxstates_cations2dec.append(oxstates_cation)
    print("Cations oxstates = ", oxstates_cations2dec)
    
    oxstates_anions2dec=[]
    for e in ANIONS2DECORATE:
        if e not in OXIDATION_STATES.keys():
            if DEBUG: print("no oxidation states available for e="+str(e)+", skipping")
            break
        #need to find ONLY positive oxstates for cation
        #here's a good example: aflow-dev --proto=ABC_hR3_160_a_a_a-001.BAC:Fe:Mn:O [decv=['BAC']; decv_dup=['BAC', 'BCA']] (new) vs. aflow-dev --proto=ABC_hR3_160_a_a_a-001.ABC:C:O:S (std)
        #this works if Fe is -2, we are not interested in this case
        #otherwise fix the test below
        #oxstates_possible_cations.append(OXIDATION_STATES[e])
        oxstates_anion=[ i for i in OXIDATION_STATES[e] if i<=0 ]
        if not oxstates_anion:
            if DEBUG: print("no positive oxidation states for e="+str(e)+", skipping")
            break
        oxstates_anions2dec.append(oxstates_anion)
    print("Anions oxstates = ", oxstates_anions2dec)
    
#--Oldoxstates match new oxstates
    match_ox = []
    for i in range(NARY):
        empty = [] 
        match_ox.append(empty)
    for atomposition, oldox in enumerate(old_oxstates):    
        for cationindex, cationlist in enumerate(oxstates_cations2dec):
            for cat_ox in cationlist:
                if cat_ox == oldox:
                    match_ox[atomposition].append(CATIONS2DECORATE[cationindex])
        for anionindex, anionlist in enumerate(oxstates_anions2dec):
            for an_ox in anionlist:
                if an_ox == oldox:
                    match_ox[atomposition].append(ANIONS2DECORATE[anionindex])

#---Complete decoration 
    notempty = True
    for i in range(NARY):
        if len(match_ox[i]) == 0:
            notempty = False 
            print(f'No new ion for position {i}')
    print("Match oxdiation states = ", match_ox)

#---Combinations 
    proto_dec_set=set()
    if notempty:
        combinations = list(itertools.product(*match_ox))
        print("Combinations = ", combinations)
    
        if len(set(old_oxstates)) != len(old_oxstates):
            ox_atompositions = {}
            for atomposition, ox in enumerate(old_oxstates):
                if ox not in ox_atompositions:
                    ox_atompositions[ox] = []
                ox_atompositions[ox].append(atomposition)
            duplicate_positions = {ox: atompositions for ox, atompositions in ox_atompositions.items() if len(atompositions) > 1}
            print("duplicate_positions = ", duplicate_positions)
        
            deletionlist = [] 
            if duplicate_positions != {}:
                for candidate in combinations:
                    elements = []
                    for atompositions in duplicate_positions.values():
                        for i in atompositions:
                            elements.append(candidate[i]) 
                        if len(set(elements)) == 1:
                            deletionlist.append(candidate)
                            print("repeated element = ", candidate)
            for todel in deletionlist:
                combinations.remove(todel)
            print("remove repeated elements = ", combinations)

#---Decorations
        for combo in combinations:
  
            elements_new, sorted_dec =zip(*sorted(zip(list(combo), ["A", "B", "C"])))
            print("elements_new = ", elements_new)
            print("sorted_dec = ", sorted_dec)
            
   
            _proto_dec=label+"."+"".join(sorted_dec)+":"+":".join(elements_new)
            print(_proto_dec)
            proto_dec_set.add(_proto_dec)
            command1=AFLOW_BIN+" --proto="+_proto_dec
            print("RUN: "+command1)
            command2=AFLOW_BIN+" --proto="+str(label)+".ABC:"+":".join(elements_sorted)
            print("COMPARE TO: "+command2+"\n")
    
    
#--TEST Protos 
            if False:
                out2,err2=issue_command(command2)
                xstr2=structure(POSCAR_Lines=out2.splitlines())
                            #print(out2,err2,xstr2)
                for dec_dup in decv_dup:
                    _proto_dec_dup=label+"."+"".join(dec_dup)+":"+":".join(elements_new)
                    command1_dup=AFLOW_BIN+" --proto="+_proto_dec_dup
                                #if(command1!=command1_dup): 
                                #print("RUN(dup): "+command1_dup+"\n")
    
                    out1,err1=issue_command(command1_dup)
                    xstr1=structure(POSCAR_Lines=out1.splitlines())
                    if DEBUG: print(out1) #xstr1)
                    if DEBUG: print(out2) #xstr2)
                    _O_index=-1
                    for i,s in enumerate(xstr1.species):
                                        #if s==anion_new: #should be from ANIONS2DECORATE #"O":
                        if s in ANIONS2DECORATE:
                            _O_index=i
                            break
                    if _O_index==-1:
                        sys.exit("_O_index not found")
                    _anion_index=-1
                    for i,s in enumerate(xstr2.species):
                        if s in anion:  #should be from ANIONS2SEARCH
                            anion_index=i
                            break
                    if _anion_index==-1:
                        sys.exit("_anion_index not found")
                    if DEBUG: print("_O_index="+str(_O_index)+", num_each_type="+str(xstr1.num_each_type[_O_index]))
                    if DEBUG: print("_anion_index="+str(_anion_index)+", num_each_type="+str(xstr2.num_each_type[_anion_index]))
                    if xstr1.num_each_type[_O_index]!=xstr2.num_each_type[_anion_index]:
                        sys.exit("Found mismatch in num_each_type, check")
                    atoms_index_1=0
                    for i in range(_O_index): atoms_index_1+=xstr1.num_each_type[i]
                    if DEBUG: print("atoms_index_1="+str(atoms_index_1))
                    if DEBUG: print("atom["+str(atoms_index_1)+"].coordF["+str(0)+"]="+str(xstr1.atoms[atoms_index_1].coordF[0]),"atom["+str(atoms_index_1)+"].coordF["+str(1)+"]="+str(xstr1.atoms[atoms_index_1].coordF[1]),"atom["+str(atoms_index_1)+"].coordF["+str(2)+"]="+str(xstr1.atoms[atoms_index_1].coordF[2]))
                    atoms_index_2=0
                    for i in range(_anion_index): atoms_index_2+=xstr2.num_each_type[i]
                    if DEBUG: print("atoms_index_2="+str(atoms_index_2))
                    if DEBUG: print("atom["+str(atoms_index_2)+"].coordF["+str(0)+"]="+str(xstr2.atoms[atoms_index_2].coordF[0]),"atom["+str(atoms_index_2)+"].coordF["+str(1)+"]="+str(xstr2.atoms[atoms_index_2].coordF[1]),"atom["+str(atoms_index_2)+"].coordF["+str(2)+"]="+str(xstr2.atoms[atoms_index_2].coordF[2]))
                    if xstr1.atoms[atoms_index_1].coordF[0]!=xstr2.atoms[atoms_index_2].coordF[0]:
                        sys.exit("Found mismatch in atom sites[0], check")
                    if xstr1.atoms[atoms_index_1].coordF[1]!=xstr2.atoms[atoms_index_2].coordF[1]:
                        sys.exit("Found mismatch in atom sites[1], check")
                    if xstr1.atoms[atoms_index_1].coordF[2]!=xstr2.atoms[atoms_index_2].coordF[2]:
                        sys.exit("Found mismatch in atom sites[2], check")
                    if DEBUG: print("")
                    if DEBUG: print("")
                    if DEBUG: print("")
                    if DEBUG: print("")

    for proto_dec in sorted(proto_dec_set):
        if PERFORM_COMPARISONS:
            proto_command="if [ ! -s "+DIR_COMPARISONS+"/compare_"+proto_dec+".json ]; then mkdir -p "+DIR_COMPARISONS+"; "+AFLOW_BIN+" --proto="+proto_dec+" | "+AFLOW_BIN+" --compare2database --properties=enthalpy_formation_atom --screen_only --if DEBUG: print=json --quiet > "+DIR_COMPARISONS+"/compare_"+proto_dec+".json"+"; fi"
        else:
            proto_command=AFLOW_BIN+" --aflow_proto="+proto_dec+" --run_relax_static --bader --noldau2 --keep_nomix --aflow_db_run --force"
        if DEBUG: print(proto_command)
        COMMANDS.append(proto_command)
################################################################

filename="commands_comparison.sh"
if PERFORM_COMPARISONS==False:
    filename="commands_aflowin.sh"

with open(filename,"w") as fout:
    fout.write("\n".join(COMMANDS))


