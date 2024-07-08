# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 11:42:20 2024

@author: shawn
"""
#-----------------------------------------------
#collect data about old compound and check line matches target
import os
import sys
sys.path.insert(0,os.path.expanduser("~")+"/programs/python_modules/")
from corey_python_functions import *
from corey_VASP_crystal_functions import structure
from xvariables import *
from xdefinitions import *
import itertools
import copy
import json

#default no debug. if any line contains --debug, debug becomes true
debug_flag_set=False    
if len(sys.argv)>1:   
    for arg in sys.argv:
        if arg=="--debug":
            debug_flag_set=True

DEBUG=(False or debug_flag_set)  #set state to variable "DEBUG"
TEST_PROTOS=(False or DEBUG)     #set state of TEST_PROTOS to FALSE by default, TRUE if DEBUG is true (QUESTION: why not just DEBUG?)

COMMANDS=[]
################################################################
with open("README_PROTO_TRUNCATED.TXT","r") as fin:     #open file
    lines=fin.read().splitlines()

for line in lines:
    parts=line.split()
    lparts=parts[0].split(".")
    label=lparts[0]
    dec=lparts[1]
    compound=parts[4]
    if "-" in compound:
        if DEBUG: print(compound)
        compound=compound.split("-")[1]   #making sure all lines are split by alphebet-only words. separating lines by white space, and then if "-" then "-"
        if DEBUG: print(compound)
        #sys.exit()
    elements=getElements(compound)     
    elements_sorted=sorted(elements)   #KC+TL20231009                   #get the elements, if length is not NARY, say prototype wrong and then exit
    print("label="+str(label),"dec="+str(dec),"compound="+str(compound),"elements="+str(elements))
    if len(elements)!=NARY: sys.exit("odd "+str(NARY)+"-nary prototype: "+label)
    
    
    
#get old_oxstates    
    proto=label+":"+":".join(elements_sorted)  #you don't need dec here, it's always alphabetical
    if DEBUG: print("trying to open file="+str(DIR_OXIDATIONS+"/oxidations_"+proto+".json"))
    try: 
        with open(DIR_OXIDATIONS+"/oxidations_"+proto+".json","r") as fin:
            data=json.load(fin)
    except:
        if DEBUG: print("skipping proto="+str(proto)+", no json file found")
        continue
    #{'Cl': {'0': 5, '1': 5}, 'K': {'2': 1, '3': 1}, 'O': {'4': -2, '5': -2, '6': -2, '7': -2, '8': -2, '9': -2}}
    if not data.keys():
        if DEBUG: print("skipping proto="+str(proto)+", no oxidation data available")
    has_consistent_oxstates=True
    old_oxstates=[]
    for e in data.keys():
        ox0=-999
        for ia in data[e].keys():
            if ox0==-999:   #choose first one
                ox0=data[e][ia]
                old_oxstates.append(ox0)
            ox=data[e][ia]
            if DEBUG: print("oxidation_state["+str(e)+"]="+str(ox))
            if ox0!=ox:
                if DEBUG: print("multi-valent proto="+str(proto)+", skipping")
                has_consistent_oxstates=False
    if has_consistent_oxstates==False:   #QUESTION: why not in the for bloci on line 75? outside both for loops, only affects last round of loop... i guess it's not a "where?" but a "any at all?"
        continue
    #oxstates=[2,2,-2]   #REMOVE ME
    if DEBUG: print("oxstates="+str(oxstates))
    if len(old_oxstates)!=NARY:    #if not every atom's oxidation states is appended to oxstates, theres a bug so need rerun
        sys.exit("len(oxstates)!=NARY")


        
#get old anions       
    oldanion=[]
    for e in ANIONS2SEARCH:
        if e in elements:
            oldanion.append(e)
    if len(oldanion)==0:
        sys.exit("could not identify anion")
 
        
#get new cation oxstates (list of lists)
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
        #QUESTION: maybe add "if OXIDATION_STATES[e] >0" above?
        oxstates_cation=[ i for i in OXIDATION_STATES[e] if i>=0 ]
        if not oxstates_cation:
            if DEBUG: print("no positive oxidation states for e="+str(e)+", skipping")
            break
        oxstates_cations2dec.append(oxstates_cation)
    

#get new anion oxstates (list of lists)
    oxstates_anions2dec=[]
    for e in ANIONS2DECORATE:
        if e not in OXIDATION_STATES.keys():
            if DEBUG: print("no oxidation states available for e="+str(e)+", skipping")
            break
        oxstates_anion=[ i for i in OXIDATION_STATES[e] if i<=0 ]
        if not oxstates_anion:
            if DEBUG: print("no negative oxidation states for e="+str(e)+", skipping")
            break
        oxstates_anions2dec.append(oxstates_anion)
        
        

#see if any new oxstates match old oxstates, and in which position 
    matchox_list=[]
    for i in len(old_oxstates):
        empty=[]
        matchox_list.append(empty)
    for atomposition, oldox in enumerate(old_oxstates):
        for oxposition,cationox in enumerate(oxstates_cations2dec):
            for oxcat in cationox:
                if oxcat==oldox:
                    matchox_list[atomposition].append(CATIONS2DECORATE[oxposition])
        for oxposition,anionox in enumerate(oxstates_anions2dec):
            for oxan in anionox:
                if oxan==oldox:
                    matchox_list[atomposition].append(ANIONS2DECORATE[oxposition])
    
    ionpossible = True
    for i in range(len(matchox_list)):
        if len(matchox_list[i])==0:
            print(f'no new ion for position{i}') #for complete decoration
            ionpossible = False
            
            #is it complete redecoration for each prototype?
            #do we need a new section when oxygen is in old compound and need to remain (variable =True/False after getElements)
            

#combination of matching ox ions
    if ionpossible:
    


#decoration to see if new compound from dif dec exist (and it probably does)



#test proto (new)



#creating and writing commands into file
            
    
        
    
    