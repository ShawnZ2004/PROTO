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
    #elements=sorted(elements)   #KC+TL20231009
    print("label="+str(label),"dec="+str(dec),"compound="+str(compound),"elements="+str(elements))
    if len(elements)!=NARY: sys.exit("odd "+str(NARY)+"-nary prototype: "+label) #Print error

    #[GOOD IDEA, let's be more robust]count_metals=0
    #[GOOD IDEA, let's be more robust]for e in elements:
    #[GOOD IDEA, let's be more robust]    if e in GROUP_METALS:
    #[GOOD IDEA, let's be more robust]        count_metals+=1
    #[GOOD IDEA, let's be more robust]if DEBUG: print("count_metals="+str(count_metals))
    #[GOOD IDEA, let's be more robust]continue
    #[GOOD IDEA, let's be more robust]sys.exit()
    #[GOOD IDEA, let's be more robust]if count_metals!=2: continue
    """
    if len(ANIONS2SEARCH)>1 and NARY>2:
        sys.exit("this functionality is not available yet")
    if len(ANIONS2DECORATE)>1:
        sys.exit("this functionality is not available yet")
    """
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
    
#get old anions       
    oldanion=[]
    for e in ANIONS2SEARCH:
        if e in elements:
            oldanion.append(e)
    if len(oldanion)==0:
        sys.exit("could not identify anion")
 
    
#----Oxdiation states of cations and anions 2 search 

    oxstates_cations2dec=[]
    for e in CATIONS2DECORATE:
        if e not in OXIDATION_STATES.keys():
            if DEBUG: print("no oxidation states available for e="+str(e)+", skipping")
            break
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
    
#see if any new oxstates match old oxstates, and in which position 
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

#Complete decoration 
    notempty = True
    for i in range(NARY):
        if len(match_ox[i]) == 0:
            notempty = False 
            print(f'No new ion for position {i}')
    print("Match oxdiation states = ", match_ox)

#Combinations 
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
        
            if duplicate_positions != {}:
                deletionlist=[]
                for candidate in combinations:
                    for atompositions in duplicate_positions.values():
                        elements = []
                        for i in atompositions:
                            elements.append(candidate[i]) 
                        if len(set(elements)) == 1:
                            deletionlist.append(candidate)
                            print("repeated element = ", candidate)
            print("remove repeated elements = ", combinations)

        for todelete in deletionlist:
            combinations.remove(todelete)
            

    


#decoration to see if new compound from dif dec exist (and it probably does)



#test proto (new)



#creating and writing commands into file
            
    
        
    
    
