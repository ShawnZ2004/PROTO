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

#Complete decoration 
    empty = True
    for i in range(NARY):
        if len(match_ox[i]) == 0:
            empty = False 
            print(f'No new ion for position {i}')
    print("Match oxdiation states = ", match_ox)








    negs = 0
    pos = 0
    for i in oxstates: #Iterate through each oxidation state (only one - can't do multi valence)
        if i < 0: #Negative oxidation state = anion
            negs += 1
    pos = NARY - negs #Number or cations
    #if negs > 1:
        #print("two negs!!!") #Only for three elements (two anions) 
    #Make multi-valent ions 
    if len(oxstates)!=NARY: #Number of oxidation states -- one for each element
        sys.exit("len(oxstates)!=NARY")



    #get index of the anion, "O" or "I" in most cases
    anion=""
    cation=""
    cation_index=""
    for e in ANIONS2SEARCH: #Iterate through the anions needed to be searched 
        if e in elements: #Anion2search is part of the natural list of elements
            anion=e #The element in prototype compound that is being searched
            break 
    if negs > 1:
        for ox in oxstates:
            if ox > 0:
                index = oxstates.index(ox) #Index of the positive oxdiation state - cation 
        cation=elements[index]
        cation_index = index
    if not anion: #No anion in the search that matches the prototype
        print("uh oh")
        sys.exit("could not identify anion")
    if DEBUG: print("cation", cation)
    #print("anion " + anion)
    anion_index=elements.index(anion) #For only one anion 
    #print("index " + str(anion_index))
    anion_let=chr(ord("A")+anion_index) #Decoration letter 
    if DEBUG: print("anion="+str(anion),"anion_index="+str(anion_index),"anion_let="+str(anion_let))
    other_index=[ i for i in range(NARY) if i!=anion_index ] #Other indexes of elements that are not the anion of interest
    other_let=[ chr(ord("A")+i) for i in other_index ] #Other element letters (except for ABC on the anion)
    if DEBUG: print("other_index="+str(other_index),"other_let="+str(other_let))
    
    #get new anion for redecoration
    anion_new=ANIONS2DECORATE[0] #Replace with anion2decorate 
    if DEBUG: print("anion_new="+str(anion_new)) 

    proto_dec_set=set() #will save the unique part of the proto command
    #go through all possible N-ARY-1 combinations of cations
    if DEBUG: print("proto dec", proto_dec_set)
    
    #I need all possible combinations of anions of the right length (with oxygen always included)
    #anion_combos = list(itertools.combinations(ANIONS2DECORATE[:negs],negs))
    #print(anion_combos)

    cation_combos = list(itertools.combinations(CATIONS2DECORATE, NARY-negs)) #Number of cation combinations for elements in cations2decorate (3-number of anions)
    #print(cation_combos)
    
    anion_combos = list(itertools.combinations(ANIONS2DECORATE, negs)) #Number of anion combinations 
    combos = []
    for com in cation_combos:
        for cam in anion_combos:
            fin = com + cam #All cation + anion combinations from elements to decorate
            #fin = com + tuple(ANIONS2DECORATE[:negs])
            #fin = com + anion_combo
            combos.append(fin) 
    #print(combos)

    #combos = list(itertools.product(anion_combos,cation_combos))
    #print(combos)

    #for combo in list(itertools.combinations(CATIONS2DECORATE,NARY-negs)):
    for combo in combos:
        #print(combo)
        cations_new=list(combo[:(NARY-negs)]) #Cations in the combo
        #anions_new=list(ANIONS2DECORATE[:negs])
        anions_new = list(combo[NARY-negs:]) #Anions in the combo 
        #print(anions_new)
        oxstates_possible_anion=[] 
        
        
        #anions_new=list(itertools.combinations(ANIONS2DECORATE, negs))
        #print(anions_new)
        #for e in ANIONS2DECORATE[:negs]:
        for e in anions_new: #Iterate through new decorated anions 
            #print(e)
            oxstates_anion= [i for i in OXIDATION_STATES[e] if i < 0] #Dictionary of element's oxidation states -- gets oxidation states of the new anion 
            if DEBUG: print("oxstates_anion", oxstates_anion)
            oxstates_possible_anion.append(oxstates_anion)
        #print(oxstates_possible_anion)
        #cations_new=["N","S"]   #REMOVE ME
        
        if DEBUG: print("cations_new="+str(cations_new))
        
        #get all possible combinations of oxidation states between cations
        oxstates_possible_cations=[]
        for e in cations_new:
            if e not in OXIDATION_STATES.keys():
                print("no oxidation states available for e="+str(e)+", skipping\n")
                break

            #need to find ONLY positive oxstates for cation
            #here's a good example: aflow-dev --proto=ABC_hR3_160_a_a_a-001.BAC:Fe:Mn:O [decv=['BAC']; decv_dup=['BAC', 'BCA']] (new) vs. aflow-dev --proto=ABC_hR3_160_a_a_a-001.ABC:C:O:S (std)
            #this works if Fe is -2, we are not interested in this case
            #otherwise fix the test below
            #oxstates_possible_cations.append(OXIDATION_STATES[e])
            #oxstates_cation=[ i for i in OXIDATION_STATES[e] if i>=0 ]
            oxstates_cation=[i for i in OXIDATION_STATES[e] if i > 0] #Access the oxidation states of the new cations decorated
            if DEBUG: print("oxstates_cation", oxstates_cation)

            if not oxstates_cation: #No oxidation state for element that is positive
                if DEBUG: print("no positive oxidation states for e="+str(e)+", skipping")
                break
            oxstates_possible_cations.append(oxstates_cation)
        if len(oxstates_possible_cations)!=len(cations_new): continue #ISSUE: Only one possible oxidation state for the cation?
        if DEBUG: print("oxstates_possible_cations="+str(oxstates_possible_cations))
        #Cation oxidation state combinations 
        oxstates_cations_combos=list(itertools.product(*oxstates_possible_cations)) #Cartesian product, or what I call "enumerations"
        if DEBUG: print("oxstates_cations_combos="+str(oxstates_cations_combos))
        #Anion oxidation state combinations
        oxstates_anions_combos=list(itertools.product(*oxstates_possible_anion))
        #if negs > 1:
            #print("oxstates_anions_combos="+str(oxstates_anions_combos))
        #find an oxidation-state set that matches the original proto+chemistry, then get all relevant decorations (.ABC...)
        if DEBUG: print("oxstates="+str(oxstates),"anion_index="+str(anion_index))
        
        _oxstates_new=[oxstates[anion_index]] #List of anion oxidation states 
        #if negs >1:
            #_oxstates_new=[oxstates[cation_index]]
        if DEBUG: print("_oxstates_new="+str(_oxstates_new))

        oxstates_sorted=sorted(oxstates)
        oxstates_all = [] #Oxidation state combinations of cation and anion 
        for com in oxstates_cations_combos:
            for cam in oxstates_anions_combos:
                fin = com + cam
                oxstates_all.append(fin)

        #oxstates_all = list(itertools.product(oxstates_cations_combos,oxstates_anions_combos))
        #print("oxos", oxstates_all)
        
        #for combo in oxstates_cations_combos:
        for combo in oxstates_all:
            
            #oxstates_new=list(_oxstates_new)+list(combo)
            oxstates_new=combo
            #print(oxstates_new)
            """
            if negs >1:
                for combo2 in oxstates_anions_combos:
                    oxstates_new=list(combo)+list(combo2)
                    print("oxstates_new(possible)="+str(oxstates_new)+"\n")
            
                    if oxstates_sorted==sorted(oxstates_new):   #check sorted, this is a super easy way to see that we have the right oxidation states in the set
                    #ok, we have a match of oxidation states, now the tricky part, figure out the right decorations (ABC)
                        elements_new=anions_new+cations_new
                        print("elements_new(pre)="+str(elements_new))
                        elements_new,oxstates_new=zip(*sorted(zip(elements_new,oxstates_new)))    #sort two lists together by first list (elements_new)
                        print("elements_new(post)="+str(elements_new))
                        print("oxstates_new(works)="+str(oxstates_new))
                        _decv=[]
                        dec_all=list(itertools.permutations([ chr(ord("A")+i) for i in range(NARY) ]))
                        for dec in dec_all:
                            print("dec="+str(dec))

                            #THE MOST TRICKY PART
                            #if we have decoration BCA, we have two options for converting this into an integer
                            #we could translate the letters RAW: B->1; C->2; A->0: 120 #option 1
                            #or we could look at each position and ask where it's respective letter is:
                            #position 0: where is A? at position 2
                            #position 1: where is B? at position 0
                            #position 2: where is C? at position 1; 
                            #thus 201    #option 2
                            #option 2 is the one coded in aflow, very tricky...
                            dec_int1=[ ord(i)-65 for i in dec ] #option1
                            dec_int2=[]
                            for il,let in enumerate(dec):
                                dec_int2.append(dec.index(chr(ord("A")+il)))  #option2, amazingly tricky
                            if(dec_int1!=dec_int2):
                                if DEBUG: print("dec_int1="+str(dec_int1),"dec_int2="+str(dec_int2))
                                #sys.exit()
                            dec_int=list(dec_int2)  #the one coded in aflow
                            print("dec="+str(dec),"dec_int="+str(dec_int))
                            _dec_int,__oxstates_new=zip(*sorted(zip(dec_int,oxstates_new)))    #sort two lists together by first list (_dec_int), be careful: _oxstates_new is used above
                            print("__oxstates_new="+str(__oxstates_new))

                            if list(oxstates)==list(__oxstates_new):    #we get the required flipping, use this decoration (or a duplicate one that is more alphabetic)
                                print("WORKS")
                                _decv.append("".join(list(dec)))
                        print("_decv="+str(_decv))
                        with open(DIR_DECORATIONS+"/decorations_"+label+".json","r") as fin:
                            data=json.load(fin)
                        decorations=data["atom_decorations_equivalent"]
                        #_decv=['ACB', 'BCA']
                        #decorations=[['ABC', 'BAC'], ['CAB', 'ACB'], ['BCA', 'CBA']]
                        #decorations=[['ABC', 'BAC'], ['BCA', 'ACB'], ['CAB', 'CBA']] #REMOVE ME, this should only give one dec
                        print("decorations="+str(decorations))                
                        decv=set()
                        decv_dup=[]
                        for dec in _decv:
                            print("dec="+str(dec))
                            for _dset in decorations:
                                if DEBUG: print("_dset="+str(_dset))
                                dset=sorted(_dset)
                                if DEBUG: print("dset="+str(dset))
                                dec_canonical=dset[0]
                                if dec in dset:
                                    if DEBUG: print("FOUND")
                                    decv.add(dec_canonical)
                                    decv_dup.append(dec)
                        decv=list(decv)
                        print("decv="+str(decv))
                        print("decv_dup="+str(decv_dup))  #the duplicate is the non-alphabetic option that should best match the original proto+chemistry, keep so we can test the atom positions easily

                        for dec in decv:
                            _proto_dec=label+"."+"".join(dec)+":"+":".join(elements_new)
                            print(_proto_dec)
                            proto_dec_set.add(_proto_dec)
                            command1=AFLOW_BIN+" --proto="+_proto_dec
                            print("RUN: "+command1)
                            command2=AFLOW_BIN+" --proto="+str(label)+".ABC:"+":".join(elements)
                            print("COMPARE TO: "+command2+"\n")
 

            print("oxstates_new(possible)="+str(oxstates_new)+"\n")
            
            """
            if oxstates_sorted==sorted(oxstates_new):   #check sorted, this is a super easy way to see that we have the right oxidation states in the set
                #ok, we have a match of oxidation states, now the tricky part, figure out the right decorations (ABC)
                elements_new=cations_new + anions_new
                print("elements_new(pre)="+str(elements_new))
                elements_new,oxstates_new=zip(*sorted(zip(elements_new,oxstates_new)))    #sort two lists together by first list (elements_new)
                print("elements_new(post)="+str(elements_new))
                print("oxstates_new(works)="+str(oxstates_new))

                #[doesn't work if we have duplicates]#the tricky bit
                #[doesn't work if we have duplicates]dec_indices=[]
                #[doesn't work if we have duplicates]for i,oxnew in enumerate(oxstates_new):
                #[doesn't work if we have duplicates]    for j,ox in enumerate(oxstates):
                #[doesn't work if we have duplicates]        if oxstates_new[i]==oxstates[j]:
                #[doesn't work if we have duplicates]            dec_indices.append(j)
                #[doesn't work if we have duplicates]            break
                #[doesn't work if we have duplicates]dec_new=[ chr(ord("A")+i) for i in dec_indices ]
                #[doesn't work if we have duplicates]if DEBUG: print("dec_new="+str(dec_new))
                
                #go through all possible decorations, find the one that sorts the right way
                #must be done this way, otherwise you might miss duplicates that arise from cases like ["Ca","Mg"] which are both +2
                _decv=[]
                dec_all=list(itertools.permutations([ chr(ord("A")+i) for i in range(NARY) ])) #Get all permutations (order matters) for element decorations 
                for dec in dec_all:
                    #print("dec="+str(dec))

                    #THE MOST TRICKY PART
                    #if we have decoration BCA, we have two options for converting this into an integer
                    #we could translate the letters RAW: B->1; C->2; A->0: 120 #option 1
                    #or we could look at each position and ask where it's respective letter is:
                    #position 0: where is A? at position 2
                    #position 1: where is B? at position 0
                    #position 2: where is C? at position 1; 
                    #thus 201    #option 2
                    #option 2 is the one coded in aflow, very tricky...
                    dec_int1=[ ord(i)-65 for i in dec ] #option1 (order of ABC) #The letter as a number 
                    dec_int2=[]
                    for il,let in enumerate(dec):
                        dec_int2.append(dec.index(chr(ord("A")+il)))  #option2, amazingly tricky
                    if(dec_int1!=dec_int2):
                        print("dec_int1="+str(dec_int1),"dec_int2="+str(dec_int2))
                        #sys.exit()
                    dec_int=list(dec_int2)  #the one coded in aflow
                    #print("dec="+str(dec),"dec_int="+str(dec_int))
                    _dec_int,__oxstates_new=zip(*sorted(zip(dec_int,oxstates_new)))    #sort two lists together by first list (_dec_int), be careful: _oxstates_new is used above
                    print("__oxstates_new="+str(__oxstates_new))

                    if list(oxstates)==list(__oxstates_new):    #we get the required flipping, use this decoration (or a duplicate one that is more alphabetic)
                        print("WORKS")
                        _decv.append("".join(list(dec))) #Decoration order coirrect
                #print("_decv="+str(_decv))
        
                #we have our decorations, but we need to find the equivalent decoration that is the most "alphabetical", use aflow
                #decorations_AB2C_hP12_152_a_c_b-001.json
                #label="ABC2_hR4_166_a_b_c-003"  #REMOVE ME
                with open(DIR_DECORATIONS+"/decorations_"+label+".json","r") as fin: #Aflow decoration file
                    data=json.load(fin)
                decorations=data["atom_decorations_equivalent"]
                #_decv=['ACB', 'BCA']
                #decorations=[['ABC', 'BAC'], ['CAB', 'ACB'], ['BCA', 'CBA']]
                #decorations=[['ABC', 'BAC'], ['BCA', 'ACB'], ['CAB', 'CBA']] #REMOVE ME, this should only give one dec
                if DEBUG: print("decorations="+str(decorations))                
                decv=set()
                decv_dup=[]
                for dec in _decv: #Decoration order 
                    if DEBUG: print("dec="+str(dec))
                    for _dset in decorations: #All combinations
                        if DEBUG: print("_dset="+str(_dset))
                        dset=sorted(_dset)
                        if DEBUG: print("dset="+str(dset))
                        dec_canonical=dset[0]
                        if dec in dset:
                            if DEBUG: print("FOUND")
                            decv.add(dec_canonical)
                            decv_dup.append(dec)
                decv=list(decv)
                if DEBUG: print("decv="+str(decv))
                if DEBUG: print("decv_dup="+str(decv_dup))  #the duplicate is the non-alphabetic option that should best match the original proto+chemistry, keep so we can test the atom positions easily

                for dec in decv:
                    _proto_dec=label+"."+"".join(dec)+":"+":".join(elements_new)
                    print(_proto_dec)
                    proto_dec_set.add(_proto_dec)
                    command1=AFLOW_BIN+" --proto="+_proto_dec
                    print("RUN: "+command1)
                    command2=AFLOW_BIN+" --proto="+str(label)+".ABC:"+":".join(elements)
                    print("COMPARE TO: "+command2+"\n")

                    #create all of the non-alphabetical decorations so we can test to see if we get the right atom counts and atom positions for the anion (unit test)
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
                                    _anion_index=i
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

                    #[WRONG THINKING]v_new=copy.deepcopy(v)
                    #[WRONG THINKING]for i,c in enumerate(v_new):
                    #[WRONG THINKING]    if c[0]==other_let[0]: v_new[i][1]=ei
                    #[WRONG THINKING]    if c[0]==other_let[1]: v_new[i][1]=ej
                    #[WRONG THINKING]if DEBUG: print("v="+str(v),"v_new="+str(v_new))
                    #[WRONG THINKING]v_new_sorted=sorted(v_new,key=lambda x: x[1])
                    #[WRONG THINKING]if DEBUG: print("v_new_sorted="+str(v_new_sorted))
    for proto_dec in sorted(proto_dec_set): #Prototypes that work -- aflow commands 
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

