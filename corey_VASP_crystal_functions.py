#!/usr/bin/env python

import sys
from corey_python_functions import *
import numpy as np
import math #nan
import string
import copy

print_scaling_output=True

zeroTol=1e-8
index0=True    #index starts at 1 for getCoordsCHGCAR
testConvertIndex12Index0=False   #just index-1
verifyCPInternal=False
getCoordC=False
angle_DEGREES=True
rad2deg=180/np.pi
deg2rad=rad2deg**-1.0

#effective mass units
twopi=2.0*np.pi
mass_e_SI=9.10938356e-31
mass_e_u=5.485799090e-4
convert_eV2J=1.602176565e-19    #SI-units
convert_A2m=1e-10               #Angs to meters
h_bar_SI=1.054571800e-34   #SI-units

#getSIUnits=False
#if not getSIUnits:
#    A2m=eV2J=h_bar=mass_e=1
#else:
A2m=convert_A2m
eV2J=convert_eV2J
h_bar=h_bar_SI
mass_e=mass_e_SI

#definitions
NA="Nuclear Attractor"
BCP="Bond Critical Point"
RCP="Ring Critical Point"
CCP="Cage Critical Point"
CP={NA:-3,BCP:-1,RCP:1,CCP:3}
CP_inverse=dict((val,key) for key,val in CP.items())
#CP_inverse={val:key for key,val in CP.items()} #CP.items()}

def F2C(f2c,vec,scale=1.0):
    return np.dot(scale*f2c,vec)

def C2F(c2f,vec,scale=1.0):
    return np.dot(scale*c2f,vec)

def getF2C(lattice_vectors,scale=1.0):
    return np.transpose(scale*lattice_vectors)

def getC2F(lattice_vectors,scale=1.0):
    return np.linalg.inv(np.transpose(scale*lattice_vectors))

def diffFloat(f1,f2):
    return abs(f1-f2)

def sameFloat(f1,f2,eps=1e-8):
    return diffFloat(f1,f2)<eps

def sameVec(f1,f2,eps=1e-8):
    #print(len(f1),len(f2))
    #print([ [f1[i],"-",f2[i]] for i in range(len(f1)) ],[ abs(f1[i]-f2[i]) for i in range(len(f1)) ],[ abs(f1[i]-f2[i])<eps for i in range(len(f1)) ], all(abs(f1[i]-f2[i])<eps for i in range(len(f1))), all([abs(f1[i]-f2[i])<eps for i in range(len(f1))]))
    return len(f1)==len(f2) and all(sameFloat(f1[i],f2[i],eps) for i in range(len(f1)))

def countSignChanges(X):
    #0 is NOT a sign change
    sign_changes=0
    if len(X)<2:
        return sign_changes
    last_sign=np.sign(X[0])
    for x in X:
        if last_sign==0:
            last_sign=np.sign(x)
            continue
        new_sign=np.sign(x)
        if new_sign==0:
            continue
        if last_sign!=new_sign:
            sign_changes+=1
            last_sign=new_sign
    return sign_changes

def countNumberMonotonicRanges(X):
    dx=np.diff(X)
    return countSignChanges(dx)

def isMonotonicFirstDerivative(X):
    dx=np.diff(X)
    #ddx=np.diff(dx)
    #return isMonotonic(dx)
    return countSignChanges(dx)==0

def isMonotonicLTE(x,tending_toward_first_derivative):
    dx=np.diff(x)
    return isMonotonic(x) and \
            tendsToward(dx,direction=tending_toward_first_derivative)
            #isMonotonic(dx,direction=direction_first_derivative)

def isMonotonic(x,direction=None):
    dx = np.diff(x)
    if direction is None:
        return np.all(dx <= 0) or np.all(dx >= 0)
    elif direction=="+":
        return np.all(dx >= 0)
    elif direction=="-":
        return np.all(dx <= 0)
    else:
        sys.exit("ERROR: Incorrect direction input")

def tendsToward(x,direction="0"):
    if not isMonotonic(x):
        return False
    if direction=="0":
        return abs(x[0])>=abs(x[-1])
    elif direction=="inf":  #this is either plus or minus inf
        return abs(x[0])<=abs(x[-1])
    #elif direction=="inf+":
    #elif direction=="inf-":
    else:
        sys.exit("ERROR: Incorrect direction input")

AFLOW_bandgap_command="aflow --bandgap="
def getBandGapOutput(path,full_exit=True):
    out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
    output,errors=issue_command(AFLOW_bandgap_command+path)
    if not output or errors:
        if full_exit:
            sys.exit("try running: "+AFLOW_bandgap_command+path)
        else:
            out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
            return out
    output_lines=output.split("\n")
    if "Majority Spin" in output:
        if not "Minority Spin" in output:
            if full_exit:
                sys.exit("try running: "+AFLOW_bandgap_command+path)
            else:
                out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
                return out
        for line in output_lines:
            if "Majority Spin" in line:
                out[0]=line.split()[-1]
            if "Minority Spin" in line:
                out[1]=line.split()[-1]
        if not (out[0] and out[1]):
            if full_exit:
                sys.exit("issue with spin-polarized bandgap calculation, try running: "+AFLOW_bandgap_command+path)
            else:
                out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
                return out
    elif "Net Result" in output:
        for line in output_lines:
            if "Net Result" in line:
                out[0]=line.split()[-1]
        if not out[0]:
            if full_exit:
                sys.exit("issue with non-spin-polarized bandgap calculation, try running: "+AFLOW_bandgap_command+path)
            else:
                out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
                return out
    else:
        if full_exit:
            sys.exit("try running: "+AFLOW_bandgap_command+path)
        else:
            out={0:""}  #None gives problems with: "metal" not in bandgap_data[i_spin]
            return out
    #print(out)
    #sys.exit()
    return out

AFLOW_lattice_command="aflow --lattice"
def getLatticeOutput(path_poscar,full_exit=True):
    out={"type":None,"variation":None}
    if not os.path.exists(path_poscar):
        if full_exit:
            sys.exit(path_poscar+" does not exist!")
        else:
            return out
    
    command=(("bzcat" if ".bz2" in path_poscar else "cat")+" "+path_poscar+" | "+AFLOW_lattice_command)
    output,errors=issue_command(command)
    if not output:
        if full_exit:
            sys.exit("ERROR: no output - try running "+command)
        else:
            return out
    parts=output.split(",")
    if len(parts)!=2:
        if full_exit:
            sys.exit("ERROR: odd parts count - try running "+command)
        else:
            return out
    out["type"]=parts[0]
    out["variation"]=parts[1]
    return out

def findMostRelaxedPOSCAR(path,full_exit=True):
    poscar=os.path.join(path,"CONTCAR.relax")       #RAW
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.bands")        #RAW
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"CONTCAR.bands.bz2")   #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.bands.bz2")    #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"CONTCAR.static.bz2")   #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.static.bz2")    #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"CONTCAR.relax2.bz2")   #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.relax2.bz2")    #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"CONTCAR.relax1")    #RAW
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"CONTCAR.relax1.bz2")   #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.relax1.bz2")    #LIB
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.orig")    #RAW
    if os.path.exists(poscar):
        return poscar
    poscar=os.path.join(path,"POSCAR.orig.bz2")    #LIB
    if os.path.exists(poscar):
        return poscar
    if full_exit:
        sys.exit("ERROR: not POSCAR found in "+path)
    return None

def getFullXY_LTE(x_left,x_right,y_left,y_right):
    x_full=x_right
    x_full=np.insert(x_full,0,x_left[:-1])
    y_full=y_right
    y_full=np.insert(y_full,0,y_left[:-1])
    return x_full,y_full

def getMonotonicRegionLTE(x,y):
    if len(x)!=len(y):
        sys.exit("ERROR: x,y must be same size")
    if len(x)%2==0:
        sys.exit("ERROR: x,y must be ODD")
    #if size is 11
    halves=len(x)/2     #5
    quarters=len(x)/4   #2

    #print(halves,quarters)
    x_left=x[:halves+1]
    y_left=y[:halves+1]
    x_right=x[halves:]
    y_right=y[halves:]
    #redundant
    x_full,y_full=getFullXY_LTE(x_left,x_right,y_left,y_right)
    #print(x)
    #print(x_full)
    #print(y)
    #print(y_full)

    changed=False
    VERBOSE=False
    
    if VERBOSE:
        x_left_OLD=copy.deepcopy(x_left)
        y_left_OLD=copy.deepcopy(y_left)
        x_right_OLD=copy.deepcopy(x_right)
        y_right_OLD=copy.deepcopy(y_right)

    threshold_size=5

    while not (isMonotonicLTE(y_left,"0") and isMonotonicLTE(y_right,"inf") and isMonotonic(y_full)) \
            and (len(x_left)>=threshold_size and len(x_right)>=threshold_size):
        changed=True
        if VERBOSE:
            print("BEFORE")
            print(len(x_left),len(y_left))
            print(len(x_right),len(y_right))
        
        if not isMonotonicLTE(y_left,"0"):
            x_left=np.delete(x_left,0)
            y_left=np.delete(y_left,0)
        if not isMonotonicLTE(y_right,"inf"):
            x_right=np.delete(x_right,-1)
            y_right=np.delete(y_right,-1)
        
        if VERBOSE:
            print("AFTER")
            print(len(x_left),len(y_left))
            print(len(x_right),len(y_right))
        
        if len(x_left)<threshold_size or len(x_right)<threshold_size:
            break
    
        x_full,y_full=getFullXY_LTE(x_left,x_right,y_left,y_right)

        if not isMonotonic(y_full):
            x_left=np.delete(x_left,0)
            y_left=np.delete(y_left,0)
            x_right=np.delete(x_right,-1)
            y_right=np.delete(y_right,-1)

    if VERBOSE:
        if not changed:
            print("NO CHANGE")
        else:
            print(x_left_OLD,x_left)
            print(y_left_OLD,y_left)
            print(x_right_OLD,x_right)
            print(y_right_OLD,y_right)

    if not (len(x_left)>=threshold_size and len(x_right)>=threshold_size):
        return None,None,\
                None,None,\
                None,None

    return copy.deepcopy(x_full),copy.deepcopy(y_full),\
            copy.deepcopy(x_left),copy.deepcopy(y_left),\
            copy.deepcopy(x_right),copy.deepcopy(y_right)


def getMonotonicRegionLTE_OLD(x,y):
    if len(x)!=len(y):
        sys.exit("ERROR: x,y must be same size")
    if len(x)%2==0:
        sys.exit("ERROR: x,y must be ODD")
    #if size is 11
    halves=len(x)/2     #5
    quarters=len(x)/4   #2

    #print(halves,quarters)
    #print(x,y,len(x))
    x_left=x[:halves+1]
    y_left=y[:halves+1]
    x_right=x[halves:]
    y_right=y[halves:]

    changed=False
    VERBOSE=False
    
    if VERBOSE:
        x_left_OLD=copy.deepcopy(x_left)
        y_left_OLD=copy.deepcopy(y_left)
        x_right_OLD=copy.deepcopy(x_right)
        y_right_OLD=copy.deepcopy(y_right)

    threshold_size=3

    #while not (isMonotonic(y_left) and isMonotonic(y_right)):
    #while not (countNumberMonotonicRanges(y_left)<2 and countNumberMonotonicRanges(y_right)<2):
    #while not (isMonotonicFirstDerivative(y_left) and isMonotonicFirstDerivative(y_right)):
    while not (isMonotonicLTE(y_left) and isMonotonicLTE(y_right)):
        changed=True
        #print(x_left[len(x_left)-quarters:])
        #print(x[halves])
        #print(x_right[:quarters])
        if VERBOSE:
            print("BEFORE")
            print(len(x_left),len(y_left))
            print(len(x_right),len(y_right))
        
        #if not isMonotonic(y_left):
        #if not countNumberMonotonicRanges(y_left)<2:
        #if not isMonotonicFirstDerivative(y_left):
        if not isMonotonicLTE(y_left):
            x_left=np.delete(x_left,0)
            y_left=np.delete(y_left,0)
        #if not isMonotonic(y_right):
        #if not countNumberMonotonicRanges(y_right)<2:
        #if not isMonotonicFirstDerivative(y_right):
        if not isMonotonicLTE(y_right):
            x_right=np.delete(x_right,-1)
            y_right=np.delete(y_right,-1)
        
        if VERBOSE:
            print("AFTER")
            print(len(x_left),len(y_left))
            print(len(x_right),len(y_right))
        
        if len(x_left)<threshold_size or len(x_right)<threshold_size:
            break

    if VERBOSE:
        if not changed:
            print("NO CHANGE")
        else:
            print(x_left_OLD,x_left)
            print(y_left_OLD,y_left)
            print(x_right_OLD,x_right)
            print(y_right_OLD,y_right)

    if not (len(x_left)>=threshold_size and len(x_right)>=threshold_size):
        return None,None,\
                None,None,\
                None,None,\
                None,None

    threshold_size-=1   #using for middle region now, should be monotonic and small "central" region for first derivative

    #if changed:
    #    print("check")
    #    print(len(x_left),len(y_left))
    #    print(len(x_right),len(y_right))

    #do full x
    #print(x)
    x_full=x_right
    #x_full=np.insert(x_full,0,x[halves])   #already included middle
    x_full=np.insert(x_full,0,x_left[:-1])
    #print(x_full)
    #sys.exit()
    #do full y
    y_full=y_right
    #y_full=np.insert(y_full,0,y[halves])   #already included middle
    y_full=np.insert(y_full,0,y_left[:-1])

    #grab same count from both sizes, not different sizes
    grab_from_either_side=min(quarters+1,len(x_left),len(x_right))      #+1 because we already included middle
    #grab_from_left=min(quarters+1,len(x_left),len(x_right))      #+1 because we already included middle
    #grab_from_right=min(quarters+1,len(x_right))    #+1 because we already included middle

    #middle region should consist of monotonic left and monotonic right for best results
    while not (isMonotonic(y_left[len(y_left)-grab_from_either_side:]) and isMonotonic(y_right[:grab_from_either_side])):
        grab_from_either_side-=1
        if grab_from_either_side<threshold_size:
            break

    if grab_from_either_side<threshold_size:
        return None,None,\
                None,None,\
                None,None,\
                None,None

    #grab x_middle from x_left,x_right
    #x_middle=x_right[:grab_from_right]
    x_middle=x_right[:grab_from_either_side]
    #x_middle=np.insert(x_middle,0,x[halves])   #already included middle
    #x_middle=np.insert(x_middle,0,x_left[len(x_left)-grab_from_left:-1])
    x_middle=np.insert(x_middle,0,x_left[len(x_left)-grab_from_either_side:-1])
    #grab y_middle from y_left,y_right
    #y_middle=y_right[:grab_from_right]
    y_middle=y_right[:grab_from_either_side]
    #y_middle=np.insert(y_middle,0,y[halves])   #already included middle
    #y_middle=np.insert(y_middle,0,y_left[len(y_left)-grab_from_left:-1])
    y_middle=np.insert(y_middle,0,y_left[len(y_left)-grab_from_either_side:-1])
    
    #x_middle=x[halves-quarters:halves+quarters+1]
    #y_middle=y[halves-quarters:halves+quarters+1]
    #COMPARE
    #print(x,y)
    #print(x_middle,x[halves-quarters:halves+quarters+1])
    #print(y_middle,y[halves-quarters:halves+quarters+1])
    #print(len(x_middle),len(y_middle))
    #sys.exit()
    return copy.deepcopy(x_full),copy.deepcopy(y_full),\
            copy.deepcopy(x_middle),copy.deepcopy(y_middle),\
            copy.deepcopy(x_left),copy.deepcopy(y_left),\
            copy.deepcopy(x_right),copy.deepcopy(y_right)

class structure:
    #def __init__(self,name,POSCAR_lines,CHGCAR_lines,CPF_lines):
    def __init__(self,Name=None,POSCAR_Lines=None,CHGCAR_Lines=None,KPOINTS_Lines=None,EIGENVAL_Lines=None,CPF_Lines=None):
        #POSCAR
        #self.title=""
        #self.name=""
        #self.scale=0.0
        #self.volume=0.0
        #self.modeVolume=False
        #self.lattice_vectors=np.empty([3,3])
        #self.geometry=np.empty([6,1])  #a,b,c,alpha,beta,gamma
        #self.f2c=np.empty([3,3])
        #self.c2f=np.empty([3,3])
        #self.num_each_type=[]
        #self.species=[]
        #self.num_atoms=0
        #self.atoms=[]
        #CHGCAR
        #self.ngrid=np.empty(3)
        #self.num_charges=0
        #self.charges=np.empty()
        #CPF

        self.name=""
        self.created=True

        if Name is not None:
            self.name=Name
        if POSCAR_Lines is not None:
            self.POSCAR_lines=POSCAR_Lines
            if self.POSCAR_lines:
                structure.__parsePOSCAR(self)
        if CHGCAR_Lines is not None:
            self.CHGCAR_lines=CHGCAR_Lines
            if POSCAR_Lines is None and self.CHGCAR_lines:
                structure.__parsePOSCAR(self.CHGCAR_lines) #not as precise, but it works
            elif self.CHGCAR_lines:
                structure.__parseCHGCAR(self)
        if KPOINTS_Lines is not None:
            self.KPOINTS_Lines=KPOINTS_Lines
            if self.KPOINTS_Lines:
                structure.__parseKPOINTS(self)
        if EIGENVAL_Lines is not None:
            self.EIGENVAL_Lines=EIGENVAL_Lines
            if self.EIGENVAL_Lines:
                structure.__parseEIGENVAL(self)
        if CPF_Lines is not None:
            self.CPF_lines=CPF_Lines
            if self.CPF_lines:
                structure.__parseCPF(self)

    def __parsePOSCAR(self):
        #initialize local variables
        title=""
        #name=""
        scale=0.0
        volume=0.0
        modeVolume=False
        lattice_vectors=np.eye(3,3)
        reciprocal_lattice_vectors=np.eye(3,3)
        geometry=np.empty([6,1])
        f2c=np.eye(3,3)
        c2f=np.eye(3,3)
        num_each_type=[]
        species=[]
        num_atoms=0
        atoms=[]
        
        lines=self.POSCAR_lines
        
        #line 0 contains title
        line_count=0
        title=lines[line_count].strip()
        
        #line 1 contains scale
        line_count+=1
        try:
            scale_volume=float(lines[line_count].strip())
            if scale_volume<0:
                volume=abs(scale_volume)
                modeVolume=True
            else:
                scale=scale_volume
        except:
            print(": ".join(["ERROR","Issue extracting scale/volume"]))
            self.created=False
            return None

        #lines 2-4 contain lattice vectors, UNLESS it shows geometry
        #see aflow_xatom Getabc_angles() getclat()
        line_count+=1
        parts=lines[line_count].split()
        if len(parts)==3:
            for c in range(3):
                try:
                    #parts=[ float(i) for i in lines[line_count].split() ]
                    parts=lines[line_count].split()
                    lattice_vectors[c][0]=float(parts[0])
                    lattice_vectors[c][1]=float(parts[1])
                    lattice_vectors[c][2]=float(parts[2])
                    line_count+=1
                except:
                    print(lines[line_count])
                    print(lines[line_count].split())
                    print(": ".join(["ERROR","Issue extracting lattice vectors"]))
                    self.created=False
                    return None
            geometry[0]=scale*modV(lattice_vectors[0])
            geometry[1]=scale*modV(lattice_vectors[1])
            geometry[2]=scale*modV(lattice_vectors[2])
            geometry[3]=angleV(lattice_vectors[1],lattice_vectors[2])
            geometry[4]=angleV(lattice_vectors[2],lattice_vectors[0])
            geometry[5]=angleV(lattice_vectors[0],lattice_vectors[1])
            if angle_DEGREES:
                geometry[3]*=rad2deg
                geometry[4]*=rad2deg
                geometry[5]*=rad2deg
        elif len(parts)==6:
            geometry[0]=scale*float(parts[0])             #a
            geometry[1]=scale*float(parts[1])             #b
            geometry[2]=scale*float(parts[2])             #c
            geometry[3]=float(parts[3])*deg2rad     #bc
            geometry[4]=float(parts[4])*deg2rad     #ca
            geometry[5]=float(parts[5])*deg2rad     #ab
            if(geometry[5]<zeroTol):
                print(": ".join(["ERROR","The angle gamma from a to b is too small"]))
                created=False
                return None
            lattice_vectors[0][0]=geometry[0]   #a
            lattice_vectors[0][1]=0.0
            lattice_vectors[0][2]=0.0
            lattice_vectors[1][0]=geometry[1]*np.cos(geometry[5])
            lattice_vectors[1][1]=geometry[1]*np.sin(geometry[5])
            lattice_vectors[1][2]=0.0
            lattice_vectors[2][0]=geometry[2]*np.cos(geometry[4])
            lattice_vectors[2][1]=geometry[2]*(np.cos(geometry[3])-np.cos(geometry[5])*np.cos(geometry[4]))/np.sin(geometry[5])
            lattice_vectors[2][2]=np.sqrt(np.abs(geometry[2]**2.0-lattice_vectors[2][1]**2.0-lattice_vectors[2][0]**2.0))
            line_count+=1
            #print(geometry)
            #print(lattice_vectors)
            #if ab is too small, that's bad, so check if there's a problem
        else:
            print(": ".join(["ERROR","Issue extracting lattice vectors, unknown parts count (not 3 or 6)"]))
            created=False
            return None

        #get scale/volume opposites
        if modeVolume:
            scale=(volume/np.linalg.det(lattice_vectors))**(1.0/3.0)
        else:
            volume=scale**(3.0)*np.linalg.det(lattice_vectors)

        try:
            reciprocal_lattice_vectors=self.getReciprocalLattice(lattice_vectors,scale)
            #print(reciprocal_lattice_vectors)
        except:
            print(": ".join(["ERROR","Issue calculating reciprocal lattice vectors"]))
            created=False
            return None

        #get f2c and c2f
        try:
            f2c=self.getF2C(lattice_vectors) #np.transpose(lattice_vectors)
        except:
            print(": ".join(["ERROR","Issue calculating f2c"]))
            created=False
            return None

        try:
            c2f=self.getC2F(lattice_vectors) #np.linalg.inv(f2c)
        except:
            print(": ".join(["ERROR","Issue calculating c2f"]))
            self.created=False
            return None

        #print(lattice_vectors)
        #print(f2c)
        #print(c2f)

        #line 5 contains num_each_type
        poscar_4_format=True
        poscar_5_format=False
        try:
            num_each_type=[ int(i) for i in lines[line_count].split() ]
        except:
            #anticipate poscar_5_format
            poscar_4_format=False
            poscar_5_format=True
            species=lines[line_count].split()
            #now try num_each_type
            line_count+=1
            try:
                num_each_type=[ int(i) for i in lines[line_count].split() ]
            except:
                print(lines[line_count])
                print(lines[line_count].split())
                print(": ".join(["ERROR","Issue extracting num_each_type"]))
                self.created=False
                return None
            #print(species)
            #print(num_each_type)
        
        num_atoms=sum(num_each_type)
        
        #do some consistency checks with atoms
        if sum(num_each_type)!=num_atoms:
            print(": ".join(["ERROR","Mismatch between num_each_type ("+str(num_each_type)+") and num_atoms ("+str(num_atoms)+")"]))
            self.created=False
            return None

        #skip next line, doesn't really matter
        line_count+=1

        #next num_atoms lines contain atoms
        line_count+=1
        atomCount=0
        for specie_indx in range(len(num_each_type)):
            for at_indx in range(num_each_type[specie_indx]):
                specie=None
                coordF=np.empty([3])
                try:
                    #parts=[ i.strip() for i in lines[line_count].split() ]
                    parts=lines[line_count].split()
                except:
                    print(lines[line_count])
                    print(lines[line_count].split())
                    print(": ".join(["ERROR","Issue extracting coordF for atom %i" % atomCount]))
                    self.created=False
                    return None
                try:
                    coordF[0]=float(parts[0])
                    coordF[1]=float(parts[1])
                    coordF[2]=float(parts[2])
                    if len(parts)>3:
                        #only look at fourth entry, ignore the rest
                        specie=parts[3]
                    coordF=bringInCell(coordF)
                except:
                    print(": ".join(["ERROR","Issue assigning coordF for atom %i" % atomCount]))
                    self.created=False
                    return None
                try:
                    coordC=np.dot(f2c,coordF)
                except:
                    print(": ".join(["ERROR","Issue calculating coordC for atom %i" % atomCount]))
                    self.created=False
                    return None
                if not poscar_5_format and specie is not None:  #if poscar5, we already got species
                    #species.add(specie)  #could be vector of nones
                    if specie not in species: species.append(specie)    #could be a vector of nones
                a=atom(specie_indx,specie,coordF,coordC)
                #print(a.specie)
                atoms.append(a)
                line_count+=1
                atomCount+=1

        #do some consistency checks with atoms
        if num_atoms!=len(atoms):
            print(": ".join(["ERROR","Mismatch between num_atoms ("+str(num_atoms)+") and len(atoms) ("+str(len(atoms))+")"]))
            self.created=False
            return None
        #set of None doesn't work here
        #if len(num_each_type)!=len(species):
        #    print(": ".join(["ERROR","Mismatch between num_each_type ("+str(num_each_type)+") and species ("+str(species)+")"]))
        #    self.created=False
        #    return None

        #save local variables to object
        self.title=title
        #self.name=name
        self.scale=scale
        self.volume=volume
        self.modeVolume=modeVolume
        self.lattice_vectors=lattice_vectors
        self.reciprocal_lattice_vectors=reciprocal_lattice_vectors
        self.geometry=geometry
        self.f2c=f2c
        self.c2f=c2f
        self.poscar_4_format=poscar_4_format
        self.poscar_5_format=poscar_5_format
        self.num_each_type=num_each_type
        self.species=list(species)
        self.num_atoms=num_atoms
        self.atoms=atoms

        return None

    def getF2C(self,lattice_vectors=None,scale=1.0):
        if lattice_vectors is None:
            lattice_vectors=self.lattice_vectors
        return getF2C(lattice_vectors,scale) #np.transpose(scale*lattice_vectors)

    def getC2F(self,lattice_vectors=None,scale=1.0):
        if lattice_vectors is None:
            lattice_vectors=self.lattice_vectors
        return getC2F(lattice_vectors,scale) #np.linalg.inv(np.transpose(scale*lattice_vectors))

    def fixLattice(self):
        self.f2c=self.getF2C(self.lattice_vectors)
        self.c2f=self.getC2F(self.lattice_vectors)
        self.reciprocal_lattice_vectors=self.getReciprocalLattice(self.lattice_vectors,self.scale)

    def updateNewScale(self,items,scale):
        for i,item in enumerate(items):
            items[i].coordC*=self.scale/scale

    #simply changes how poscar (lattice) looks
    def reScale(self,scale=1.0,lattice_vectors=None):
        if lattice_vectors is None:
            lattice_vectors=self.lattice_vectors
        self.updateNewScale(self.atoms,scale)
        if print_scaling_output:
            print("")
            print("start reScale:",scale)
            print("scale",scale)
            print("self.scale ",self.scale)
            print("np.linalg.det(lattice_vectors) ",np.linalg.det(lattice_vectors))
            print("self.volume ",self.volume)
        self.lattice_vectors*=self.scale/scale
        self.scale=scale
        #self.volume=self.scale**(3.0)*np.linalg.det(self.lattice_vectors)   #DOES NOT CHANGE, feel free to check
        self.fixLattice()
        if print_scaling_output:
            print("change")
            print("self.scale ",self.scale)
            print("np.linalg.det(lattice_vectors) ",np.linalg.det(lattice_vectors))
            print("self.volume ",self.volume)
            print("end reScale:",scale)
            print("")

    #actually stretches/shrinks cell
    def setVolume(self,in_vol,lattice_vectors=None):
        if lattice_vectors is None:
            lattice_vectors=self.lattice_vectors
        if print_scaling_output:
            print("")
            print("start setVol:",in_vol)
            print("in_vol",in_vol)
            print("self.scale ",self.scale)
            print("np.linalg.det(lattice_vectors) ",np.linalg.det(lattice_vectors))
            print("self.volume ",self.volume)
        scale=(in_vol/np.linalg.det(lattice_vectors))**(1.0/3.0)
        if print_scaling_output:
            print("new scale ",scale)
        self.updateNewScale(self.atoms,scale)   #volume stretches/shrinks cell, changing cpos but not fpos
        self.scale=scale
        self.fixLattice()
        #self.reScale(scale=scale)              #WRONG, reScaling simply changes the format, doesn't actually change volume
        self.volume=in_vol
        if print_scaling_output:
            print("change")
            print("self.scale ",self.scale)
            print("np.linalg.det(lattice_vectors) ",np.linalg.det(lattice_vectors))
            print("self.volume ",self.volume)
            print("end setVol:",in_vol)
            print("")

    def normalizeLattice(self):
        self.reScale(1.0)
        self.setVolume(np.linalg.det(self.lattice_vectors))

    def getReciprocalLattice(self,lattice_vectors=None,scale=1.0):
        if lattice_vectors is None:
            lattice_vectors=self.lattice_vectors
        #lattice_vectors=scale*lattice_vectors
        #scale=1.0
        reciprocal_lattice_vectors=np.eye(3,3)
        norm=twopi/(np.linalg.det(lattice_vectors)*scale)   #units twopi/SCALE (lattice constant, first number in POSCAR)
        b1=norm*np.cross(lattice_vectors[1],lattice_vectors[2])
        b2=norm*np.cross(lattice_vectors[2],lattice_vectors[0])
        b3=norm*np.cross(lattice_vectors[0],lattice_vectors[1])
        #print("b1:",b1)
        #print("b2:",b2)
        #print("b3:",b3)
        #sys.exit(0)
        reciprocal_lattice_vectors[0]=b1
        reciprocal_lattice_vectors[1]=b2
        reciprocal_lattice_vectors[2]=b3
        return reciprocal_lattice_vectors

    def modifyComposition(self,new_num_each_type,new_species):
        if not new_num_each_type:
            new_num_each_type=self.num_each_type
        if not new_species:
            new_species=self.species
        
        #do some consistency checks with atoms
        if len(new_num_each_type)!=len(self.num_each_type):
            sys.exit(": ".join(["ERROR","Mismatch between new_num_each_type ("+str(new_num_each_type)+") and structure.num_each_type ("+str(self.num_each_type)+")"]))
        if sum(new_num_each_type)!=self.num_atoms:
            sys.exit(": ".join(["ERROR","Mismatch between new_new_each_type ("+str(new_species)+") and structure.num_atoms ("+str(self.num_atoms)+")"]))
        if len(new_species)!=len(self.species):
            sys.exit(": ".join(["ERROR","Mismatch between new_species ("+str(new_species)+") and structure.species ("+str(self.species)+")"]))

        atomCount=0
        for i in range(len(new_num_each_type)):
            for j in range(new_num_each_type[i]):
                self.atoms[atomCount].atom_type=i
                self.atoms[atomCount].species=new_species[i]
                atomCount+=1
        
        self.num_each_type=new_num_each_type
        self.species=new_species

        #sort
        self.species,self.num_each_type=zip(*sorted(zip(self.species,self.num_each_type)))
        self.atoms.sort(key=lambda x: x.species)

        return None

    def printPOSCAR(self,geometry=False,toFile=False,File="POSCAR"):
        output=[]
        output.append(self.title)
        if self.modeVolume:
            output.append("-%.6f" % self.volume)
        else:
            output.append("%.6f" % self.scale)
        if geometry:
            #output.append(" ".join(map(str,self.geometry)))
            output.append(" ".join([ "%.2f" % g for g in self.geomtry ]))
        else:
            for i in range(3):
                #output.append("  "+" ".join([ "{:>15}".format("%.14g" % lat) for lat in self.lattice_vectors[i] ]))
                output.append("  "+" ".join([ "{:>17}".format(to_precision(lat,15)) for lat in self.lattice_vectors[i] ]))
        output.append(" ".join(map(str,self.num_each_type)))
        pseudo_species=string.ascii_uppercase[:len(self.num_each_type)]
        output.append("Direct("+str(self.num_atoms)+") ["+"".join([ str(sp)+str(num) for sp,num in zip(pseudo_species,self.num_each_type) ])+"]")
        for i in range(self.num_atoms):
            a=self.atoms[i]
            out="   "+" ".join([ "%.14f" % cF for cF in a.coordF ])
            if a.species:
                out+=" "+a.species
            output.append(out)
        if toFile:
            with open(File,"w") as fout:
                fout.write("\n".join(output))
        else:
            print("\n".join(output))

    def __parseCHGCAR(self):
        #initialize local variables
        ngrid=np.empty(3)
        num_charges=0
        #charges=np.empty()

        lines=self.CHGCAR_lines

        line_count=len(self.POSCAR_lines) #skip POSCAR in CHGCAR

        #read in ngrid
        try:
            parts=[ int(i) for i in lines[line_count].split() ]
            ngrid[0]=parts[0]
            ngrid[1]=parts[1]
            ngrid[2]=parts[2]
            num_charges=ngrid[0]*ngrid[1]*ngrid[2]
        except:
            print(": ".join(["ERROR","Issue extracting ngrid"]))
            self.created=False
            return None

        charges=np.empty([ngrid[0],ngrid[1],ngrid[2]]) #declare here
        #now start reading in chgs
        line_count+=1
        i_charge=0      #for charge indexing
        indexX=0;indexY=0;indexZ=0  #for grid indexing
        while i_charge<num_charges:  #because i_charge starts at 1
            try:
                parts=[ float(i) for i in lines[line_count].split() ]
            except:
                print(": ".join(["ERROR","Issue extracting chgs on line %i" % line_count]))
                self.created=False
                return None
            for c in parts:
                if indexX and indexX%ngrid[0]==0:
                    indexX=0
                    indexY+=1
                if indexY and indexY%ngrid[1]==0:
                    indexY=0
                    indexZ+=1
                #print(indexX,indexY,indexZ)
                charges[indexX][indexY][indexZ]=c
                indexX+=1
                i_charge+=1
            line_count+=1

        #save local variables to object
        self.ngrid=ngrid
        self.num_charges=num_charges
        self.charges=charges

    def __parseKPOINTS(self):
        lines=self.KPOINTS_Lines
        self.kpoints={}
        line_count=1
        self.kpoints["k_grid"]=int(lines[line_count].split("!")[0])
        line_count+=2
        if lines[line_count].strip()!="reciprocal":
            sys.exit("Corey, write reader for non-reciprocal coords in KPOINTS.bands")
        line_count+=1
        kpoints=[]
        while line_count<len(lines):
            line=lines[line_count].strip()
            if not line:
                line_count+=1
                continue
            kpoint={}
            parts=line.split("!")
            label=parts[1].strip()
            #make everything italics except Gamma
            if not "Gamma" in label:
                label="\it{"+label+"}"
            label="$"+label+"$"
            #label=("$"+label+"$" if "\\" in label else label)
            kpoint["label"]=label #("$"+label+"$" if "\\" in label else label)
            #kpoint["label"]=parts[1].replace("\\","").strip()
            kpoint["coordF"]=np.asarray([ float(x) for x in parts[0].split() ])
            kpoints.append(kpoint)
            line_count+=1
        
        self.kpoints["kpoints"]=kpoints
        #print(self.kpoints["kpoints"])
        #sys.exit()

    def __parseEIGENVAL(self):
        lines=self.EIGENVAL_Lines
        line_count=5
        tokens=lines[line_count].split()
        n_kpoints=int(tokens[1])    #keep
        n_energy=int(tokens[2])
        spin_polarized=(len(lines[8].split())==3)
        spin_count=(2 if spin_polarized else 1)

        kpoints=[]  #contain {} in order of path, "kpoint" and "weight"
        BS_points=np.zeros(shape=(n_kpoints,spin_count,n_energy),dtype=float)
        line_count+=2
        
        klattice=self.reciprocal_lattice_vectors
        kf2c=getF2C(klattice)

        debugkunits=False

        for k in range(n_kpoints):
            parts=lines[line_count].split()
            kpoint=np.asarray([ float(x) for x in parts[:3] ],dtype=float)
            #print(parts)
            weight=parts[3]
            #print("kpoint",kpoint)
            #if len(kpoints):
            #    print("kpoints[-1]",kpoints[-1])
            #print(kpoints,kpoint)
            if debugkunits:
                #we did a small test here to see whether F2C gives us some radian units or actual 1/A (1/distance)
                #see here:  https://cms.mpi.univie.ac.at/vasp/vasp/Entering_all_k_points_explicitly.html
                #keep in mind here that "a" is defined by conventional cell, but we can just think of it as some 
                #convenient way of defining the matrices
                #"a" in this case is simply 2*scale*element (element = lattice_vectors[0,1]
                #remember that lattice vector length (geometry) and "a" differ by sqrt(2)
                #with this convenient scaling ("a"), we can compare to what we see on the vasp page
                #you see that in order to get the (2pi/a) units, we need to divide out 2pi/a from answer
                #therefore, answer is already given in units of 1/A, true cartesian
                scaling_factor=(twopi/(2*self.lattice_vectors[0,1]*self.scale))
                print(klattice/scaling_factor)
                test=np.asarray([0.5,0.75,0.25],dtype=float)
                print(np.dot(kf2c,test)) #*(self.scale/twopi)**2
                print(np.dot(kf2c,test)/scaling_factor) #*(self.scale/twopi)**2
                sys.exit(0)
            #do NOT remove until you collect everything at the end, because otherwise there will be 0's at the end of each array
            #if len(kpoints) and sameVec(kpoint,kpoints[-1]["coordF"]):
            #    #print(sameVec(kpoint,kpoints[-1]))
            #    #print(lines[line_count])
            #    line_count+=n_energy+2
            #    #print(lines[line_count])
            #    continue
            k_dict={"coordF":kpoint,"coordC":F2C(kf2c,kpoint),"weight":weight}
            kpoints.append(k_dict)
            line_count+=1
            
            for indx,e in enumerate(range(n_energy)):
                tokens=lines[line_count].split()
                for i_spin in range(spin_count):
                    BS_points[len(kpoints)-1,i_spin,indx]=float(tokens[i_spin+1])
                    #print("BS_points["+str(len(kpoints)-1)+","+str(i_spin)+","+str(indx)+"]",BS_points[len(kpoints)-1,i_spin,indx])
                line_count+=1
            line_count+=1
        
        self.eigenval={}
        self.eigenval["n_kpoints"]=n_kpoints
        self.eigenval["n_energy"]=n_energy
        self.eigenval["spin_polarized"]=spin_polarized
        self.eigenval["kpoints"]=kpoints
        self.eigenval["BS_points"]=BS_points
        return None

    def getHighSymmetryIndices(self,assume_unique_only=True):
        k_grid=self.kpoints["k_grid"]
        high_symmetry_indices=[]
        high_symmetry_indices.append(0)
        high_symmetry_indices.append(high_symmetry_indices[-1]+k_grid-1)
        for ki in range(1,len(self.kpoints["kpoints"])-1,2):
            if sameVec(self.kpoints["kpoints"][ki]["coordF"],self.kpoints["kpoints"][ki+1]["coordF"]):
                high_symmetry_indices.append(high_symmetry_indices[-1]+k_grid-(1 if assume_unique_only else 0))
            else:
                high_symmetry_indices.append(high_symmetry_indices[-1]+1)
        high_symmetry_indices.append(high_symmetry_indices[-1]+k_grid-1)
        #for ki in range(1,len(self.kpoints["kpoints"])-1,2):
        #    if sameVec(self.kpoints["kpoints"][ki]["coordF"],self.kpoints["kpoints"][ki+1]["coordF"]):
        #        ranges_of_continuity[-1][1]+=k_grid-(1 if assume_unique_only else 0)
        #    else:
        #        ranges_of_continuity.append([ranges_of_continuity[-1][1]+1,ranges_of_continuity[-1][1]+1+(k_grid-1)])
        #print(high_symmetry_indices)
        #sys.exit()
        self.kpoints["high_symmetry_indices"]=high_symmetry_indices

    def getKPathContinuity(self,assume_unique_only=True):
        #if assume_unique_only, then duplicate points at boundaries are not counted
        k_grid=self.kpoints["k_grid"]
        ranges_of_continuity=[]
        ranges_of_continuity.append([0,k_grid-1])
        for ki in range(1,len(self.kpoints["kpoints"])-1,2):
            if sameVec(self.kpoints["kpoints"][ki]["coordF"],self.kpoints["kpoints"][ki+1]["coordF"]):
                ranges_of_continuity[-1][1]+=k_grid-(1 if assume_unique_only else 0)
            else:
                ranges_of_continuity.append([ranges_of_continuity[-1][1]+1,ranges_of_continuity[-1][1]+1+(k_grid-1)])
        #print(ranges_of_continuity)
        #sys.exit()
        self.kpoints["ranges_of_continuity"]=ranges_of_continuity

    def getUniqueKPathEigenval(self):
        #print(self.eigenval["BS_points"].shape)
        indices_to_delete=[]
        for ki in range(len(self.eigenval["kpoints"])-1):
            if sameVec(self.eigenval["kpoints"][ki+1]["coordF"],self.eigenval["kpoints"][ki]["coordF"]):
                indices_to_delete.append(ki+1)
        for i in reversed(indices_to_delete):
            del self.eigenval["kpoints"][i]
            self.eigenval["BS_points"]=np.delete(self.eigenval["BS_points"],i,0)
        #print(len(self.eigenval["kpoints"]))
        #print(self.eigenval["BS_points"].shape)
        #print(self.eigenval["BS_points"][:,:,9])
        #sys.exit()

    def adjustEIGENVAL_fermi(self,fermi=None,fermi_file=None):
        if fermi is None:
            if fermi_file is None:
                sys.exit("Please build automatic reader for DOSCAR to extract fermi")
            else:
                lines=read_file_lines(fermi_file)
                fermi=float(lines[5].split()[3])
        fermi_energy=float(fermi)
        self.eigenval["BS_points"]-=fermi_energy
        return None

    def getKpointGrid(self,k_grid=None,kpoints_file=None):
        if not hasattr(self,"kpoints"):
            self.kpoints={}
        if k_grid is None:
            if kpoints_file is None:
                sys.exit("Please build automatic reader for KPOINTS to extract k_grid")
            else:
                lines=read_file_lines(kpoints_file)
                k_grid=int(lines[1].split("!")[0])
        self.kpoints["k_grid"]=k_grid
        return None

    def getKpointPathVector(self,getLabels=True):
        if getLabels and not hasattr(self,"kpoints"):
            sys.exit("Need to parse KPOINTS file first before we can get labels")

        k_grid=self.kpoints["k_grid"]
        kpath=[0.0] #better to start with list and append, then convert to array
        
        kpoints=copy.deepcopy(self.eigenval["kpoints"])
        BS_points=copy.deepcopy(self.eigenval["BS_points"])
        special_hs_kpoints_coords=[]
        labels=[]

        ki=0
        count_hs_kpoints=0  #count of high symmetry kpoints defining grid
        #for i,k in enumerate(self.eigenval["kpoints"]):
        #    print(i,k["coordF"])
        #for i,k in enumerate(self.kpoints["kpoints"]):
        #    print(i,k["coordF"])

        indices_to_delete=[]

        while ki<len(self.eigenval["kpoints"])-1:
            if not ki:
                #first hs kpoint
                if getLabels:
                    #check of stupidity
                    if not sameVec(self.eigenval["kpoints"][ki]["coordF"],self.kpoints["kpoints"][count_hs_kpoints]["coordF"]):
                        print("SHIT")
                        for i,k in enumerate(self.eigenval["kpoints"]):
                            print(i,k["coordF"])
                        for i,k in enumerate(self.kpoints["kpoints"]):
                            print(i,k["coordF"])
                        sys.exit("Check KPOINTS file vs. EIGENVAL for mismatch in kpoints ("+str(ki)+","+str(count_hs_kpoints)+"): "+str(self.eigenval["kpoints"][ki]["coordF"])+" vs. "+str(self.kpoints["kpoints"][count_hs_kpoints]["coordF"]))
                    label={}
                    label["label"]=self.kpoints["kpoints"][count_hs_kpoints]["label"]
                    label["coord"]=kpath[-1]
                    labels.append(label)
                    count_hs_kpoints+=1
                special_hs_kpoints_coords.append(kpath[-1])
            else:
                #all other hs points except last
                if getLabels:
                    #check of stupidity
                    if not sameVec(self.eigenval["kpoints"][ki]["coordF"],self.kpoints["kpoints"][count_hs_kpoints]["coordF"]):
                        print("SHIT")
                        for i,k in enumerate(self.eigenval["kpoints"]):
                            print(i,k["coordF"])
                        for i,k in enumerate(self.kpoints["kpoints"]):
                            print(i,k["coordF"])
                        sys.exit("Check KPOINTS file vs. EIGENVAL for mismatch in kpoints ("+str(ki)+","+str(count_hs_kpoints)+"): "+str(self.eigenval["kpoints"][ki]["coordF"])+" vs. "+str(self.kpoints["kpoints"][count_hs_kpoints]["coordF"]))
                    if not sameVec(self.eigenval["kpoints"][ki+1]["coordF"],self.kpoints["kpoints"][count_hs_kpoints+1]["coordF"]):
                        print("SHIT")
                        for i,k in enumerate(self.eigenval["kpoints"]):
                            print(i,k["coordF"])
                        for i,k in enumerate(self.kpoints["kpoints"]):
                            print(i,k["coordF"])
                        sys.exit("Check KPOINTS file vs. EIGENVAL for mismatch in kpoints ("+str(ki+1)+","+str(count_hs_kpoints+1)+"): "+str(self.eigenval["kpoints"][ki+1]["coordF"])+" vs. "+str(self.kpoints["kpoints"][count_hs_kpoints+1]["coordF"]))
                if sameVec(self.eigenval["kpoints"][ki+1]["coordF"],self.eigenval["kpoints"][ki]["coordF"]):
                    indices_to_delete.append(ki+1)
                    if getLabels:
                        label={}
                        label["label"]=self.kpoints["kpoints"][count_hs_kpoints]["label"]
                        label["coord"]=kpath[-1]
                        labels.append(label)
                        count_hs_kpoints+=2
                else:
                    kpath.append(kpath[-1])
                    if getLabels:
                        label={}
                        label["label"]=self.kpoints["kpoints"][count_hs_kpoints]["label"]+"|"+self.kpoints["kpoints"][count_hs_kpoints+1]["label"]
                        label["coord"]=kpath[-1]
                        labels.append(label)
                        count_hs_kpoints+=2
                ki+=1
                special_hs_kpoints_coords.append(kpath[-1])
            #print("START")
            for kii in range(k_grid-1):
                #print(ki,len(self.eigenval["kpoints"]),range((1 if ki else 0),k_grid))
                kcdiff=self.eigenval["kpoints"][ki+1]["coordC"]-self.eigenval["kpoints"][ki]["coordC"]
                kpath.append(np.linalg.norm(kcdiff)+kpath[-1])
                ki+=1
                #print(kii)
        #get the last label
        if getLabels:
            label={}
            label["label"]=self.kpoints["kpoints"][count_hs_kpoints]["label"]
            label["coord"]=kpath[-1]
            labels.append(label)
        special_hs_kpoints_coords.append(kpath[-1])
        for i in reversed(indices_to_delete):
            del kpoints[i]
            BS_points=np.delete(BS_points,i,0)
        self.eigenval["kpoints_plot"]=kpoints
        self.eigenval["BS_points_plot"]=BS_points
        self.eigenval["kpath"]=np.asarray(kpath)
        self.eigenval["labels"]=labels
        self.eigenval["special_hs_kpoints_coords"]=np.asarray(special_hs_kpoints_coords)

    def getBandExtremaIndex(self,kpoints,BS_points,regional_index_count=5,band_to_look="valence_max",away_from_edge=False):
        if not (band_to_look=="valence_max" or band_to_look=="conduction_min"):
            sys.exit("getBandExtremaIndex(): not sure which band I'm looking at.")
        k_index=0
        e_index=0
        dist_from_zero=1e9
        compare_VM_zero=False   #if true, we find band closest to fermi energy, otherwise, find one closest to cushion above 0
        cushion_above_zero=0.75 #0.1 #eV
        cushion_away_edge=3 #index
        cushion_same_band=0.01 #eV
        #regional average will either be close to 0 (more negative) if valence_max, or strongly positive if conduction min
        #heavy and light same for both, since charges flip, outer band is lighter (smaller concavity) than inner band
        #we want regional average to approach zero ALWAYS (assume lighter variant always), therefore look for abs
        #for valence max, we want regional average to go from very negative to 0, so we look for largest (closest to 0) value here
        #for conduction min, we want regional average to go from very large to 0, so we look for smallest (closest to 0) value here
        #we could take abs and just look for smallest to 0, but leave like this in case we want heavy variant, etc.
        regional_average_keep=(-1e9 if band_to_look=="valence_max" else 1e9)
        for ei in range(BS_points.shape[1]):
            y=( BS_points[cushion_away_edge:-cushion_away_edge+1,ei] if away_from_edge else BS_points[:,ei] )
            index_interest=( np.argmax(y) if band_to_look=="valence_max" else np.argmin(y) )
            val_interest=y[index_interest]
            #new, try to get widest one
            _k_index=( index_interest+cushion_away_edge if away_from_edge else index_interest )
            krange=self.getRegionalIndices(kpoints,_k_index,regional_index_count)
            BS_regional=BS_points[krange[0]:krange[1],ei]
            regional_average=np.average(BS_regional)
            relevant_distance=( (abs(val_interest) if compare_VM_zero  else cushion_above_zero-val_interest) if band_to_look=="valence_max" else val_interest )
            if band_to_look=="valence_max":
                if val_interest<cushion_above_zero and relevant_distance<dist_from_zero:
                    dist_from_zero=relevant_distance
            elif band_to_look=="conduction_min":
                if val_interest>cushion_above_zero and relevant_distance<dist_from_zero:
                    dist_from_zero=relevant_distance
            else:
                sys.exit("What band are we looking at?")        
            #print(regional_average,regional_average_keep,dist_from_zero)
            #print(val_interest,cushion_above_zero,"<")
            #print(abs(val_interest),dist_from_zero,"<0.01 (abs diff)")
            #print(regional_average,regional_average_keep,">")
            #print("")
            #print("index_interest:",index_interest)
            #print("val_interest:",val_interest)
            #print("krange:",krange)
            #print("BS_regional:",BS_regional)
            #print("regional_average:",regional_average)
            #sys.exit()
            #if ( (val_interest<cushion_above_zero and abs(val_interest)<dist_from_zero) if band_to_look=="valence_max" else \
            #        (val_interest>cushion_above_zero and val_interest<dist_from_zero) ):
            if ( (val_interest<cushion_above_zero and abs(diffFloat(relevant_distance,dist_from_zero))<cushion_same_band and regional_average>regional_average_keep) if band_to_look=="valence_max" else \
                    (val_interest>cushion_above_zero and abs(diffFloat(relevant_distance,dist_from_zero))<cushion_same_band and regional_average<regional_average_keep) ):
                #dist_from_zero=(abs(val_interest) if band_to_look=="valence_max" else val_interest)
                #print("NEW DIST")
                #print(dist_from_zero)
                k_index=( index_interest+cushion_away_edge if away_from_edge else index_interest )
                e_index=ei #(index_interest+cushion_away_edge if away_from_edge else index_interest)
                #print("NEW regional average")
                #print(regional_average)
                regional_average_keep=regional_average
            #print(val_interest)
        #    print(y)
        #    print(y)
        #print(k_index,e_index)
        #sys.exit()

        #for ki,kpoint in enumerate(kpoints):
        #    if away_from_edge and ki<cushion_away_edge or ki>(len(kpoints)-1)-cushion_away_edge:
        #        #print(ki,"HERE")
        #        continue
        #    y=BS_points[ki,:]
        #    index_interest=( np.argmax(y) if band_to_look=="valence_max" else np.argmin(y) )
        #    val_interest=y[index_interest]
        #    print(y)
        #    print(val_interest)
        #    #sys.exit()

        #    if ( (val_interest<cushion_above_zero and abs(val_interest)<dist_from_zero) if band_to_look=="valence_max" else \
        #            (val_interest>cushion_above_zero and val_interest<dist_from_zero) ):
        #        dist_from_zero=(abs(val_interest) if band_to_look=="valence_max" else val_interest)
        #        print("NEW DIST")
        #        print(dist_from_zero)
        #        k_index=ki
        #        e_index=index_interest #(index_interest+cushion_away_edge if away_from_edge else index_interest)
            #for ei,energy in enumerate(BS_points[ki,:]):
            #    if away_from_edge and ki<3 or ki>len(kpoints)-3:
            #        continue
            #    if (abs(energy)<dist_from_zero if band_to_look=="valence_max" else (energy>cushion_above_zero and energy<dist_from_zero)):
            #        dist_from_zero=(abs(energy) if band_to_look=="valence_max" else energy)
            #        k_index=ki
            #        e_index=ei
        #print(k_index,e_index)
        #sys.exit()
        return (k_index,e_index)

    def getRegionalIndices(self,vec,index,count):
        #make sure odd so middle point is always minimum or maximum
        if count%2==0:
            count-=1
        start=0
        stop=len(vec)-1
        if count>=len(vec):
            return (start,stop)
        start=index-count/2
        #print(index)
        #print(start)
        while start>=0:
            #assume continuous
            if start<len(vec) and (start+count)<len(vec):
                return (start,start+count)
            start-=1
        #reset
        stop=index-count/2
        while stop<len(vec):
            #assume continuous
            if stop-count>=0 and stop-count<len(vec) and stop<len(vec):
                print(stop-count,stop)
                return (stop-count,stop)
            stop+=1

    def getEffectiveMass(self,index_region=5,plot=True,i_spin=0):
        output={}
        include_labels=True
        #index_region=5    #better if odd, min/max is in the middle
        if plot:
            if not hasattr(self,"kpoints"):
                print("If you want labels, parse KPOINTS file")
                include_labels=False
        #print(len(self.eigenval["BS_points"]))
        #self.getUniqueKPathEigenval()  #bad idea, sometimes the path has a hard break, just deal with it in kpointpathvector
        #print(len(self.eigenval["BS_points"]))
        #sys.exit()
        self.getKpointPathVector()
        kpoints=self.eigenval["kpoints_plot"]
        #kpoints_path=self.getKpointPathVector()
        #kpoints_path/=A2m  #1/A to 1/m
        #print(kpoints_path)
        #kpoints_path=[ x*3.0*(twopi) for x in kpoints_path ]
        spin_count=(2 if self.eigenval["spin_polarized"] else 1)
        if i_spin>spin_count-1:
            print(i_spin,spin_count-1)
            sys.exit("ERROR: requesting spin-polarized of a non-spin-polarized calculation")
        #for i_spin in range(spin_count):
        print("for spin",i_spin)
        BS_points=self.eigenval["BS_points_plot"][:,i_spin,:]   #kpoints,spin,energy
        for band in ["conduction_min","valence_max"]: #["valence_max"]: #["conduction_min","valence_max"]:
            output[ "hole" if band=="valence_max" else "electron" ]={}
            print(("\t"+"hole effective mass (band=="+band+"):" if band=="valence_max" else "\t"+"electron effective mass (band=="+band+"):"))
            #print(band)
            (k_index,e_index)=self.getBandExtremaIndex(kpoints,BS_points,regional_index_count=index_region,band_to_look=band,away_from_edge=True) #np.zeros(shape=(n_kpoints,spin_count,n_energy),dtype=float)
            ##################
            #david's stuff
            effective_mass_david=h_bar**2.0*(modV(kpoints[k_index+1]["coordC"]-kpoints[k_index]["coordC"])/A2m)**2.0/(2.0*(BS_points[k_index+1,e_index]-BS_points[k_index,e_index])*eV2J)/mass_e
            output[ "hole" if band=="valence_max" else "electron" ]["linear"]=abs(effective_mass_david)
            if band=="valence_max" and effective_mass_david>0:
                sys.exit("hole effective mass >0: "+str(effective_mass_david))
            elif band=="conduction_min" and effective_mass_david<0:
                sys.exit("electron effective mass <0: "+str(effective_mass_david))
            print("\t","\t","linear approximation:",abs(effective_mass_david))
            ##################
            #print((k_index,e_index))
            krange=self.getRegionalIndices(kpoints,k_index,index_region)
            #print(BS_points[:,e_index])
            #print(BS_points[k_index,e_index])
            #sys.exit(0)
            
            #print(krange)
            x=copy.deepcopy(self.eigenval["kpath"][krange[0]:krange[1]])
            y=copy.deepcopy(BS_points[krange[0]:krange[1],e_index])
            
            x/=A2m #*(twopi/self.scale) #/(twopi/self.scale)
            y*=eV2J
            #do fits on x,y
            #print("x",x)
            #print("y",y)
            z=np.polyfit(x,y,2)
            p=np.poly1d(z)
            xp=np.linspace(x[0],x[-1],100)
            sop=z[0]    #second order poly
            #print("second_order_poly:",sop)
            second_derivative=2.0*sop
            effective_mass=(h_bar**2)*1.0/(second_derivative)/mass_e
            output[ "hole" if band=="valence_max" else "electron" ]["curvature"]=abs(effective_mass)
            if band=="valence_max" and effective_mass>0:
                sys.exit("hole effective mass >0: "+str(effective_mass))
            elif band=="conduction_min" and effective_mass<0:
                sys.exit("electron effective mass <0: "+str(effective_mass))
            #print("second_derivative:",second_derivative)
            #print("eff_mass:",(h_bar**2/modV(kex))*1.0/(second_derivative) #/(2*np.pi))
            #print("eff_mass:",(h_bar**2)*1.0/(second_derivative)/mass_e #/(2*np.pi))
            #print("eff_mass_new:",1.0/(second_derivative)/mass_e_u #/(2*np.pi))
            print("\t","\t","full curvature:",abs(effective_mass))
            #print(len(x_full))
            #print(len(y_full))
            #sys.exit()
            if plot:
                import matplotlib.pyplot as plt
                from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
                x_full=copy.deepcopy(self.eigenval["kpath"])
                x_full/=A2m
                #print(x_full)
                for ie in range(self.eigenval["n_energy"]):
                    y_full=copy.deepcopy(BS_points[:,ie])
                    y_full*=eV2J
                    #print(y_full)
                    plt.plot(x_full,y_full,"-",linewidth=0.5)
                if band=="conduction_min":
                    plt.plot(xp,p(xp),"-",linewidth=4,color=[0,0,1])
                elif band=="valence_max":
                    plt.plot(xp,p(xp),"-",linewidth=4,color=[1,0,0])
                else:
                    sys.exit("Not sure what band this is.")
                if include_labels:
                    labels=[ label["label"] for label in self.eigenval["labels"] ]
                    coords=np.asarray([ label["coord"] for label in self.eigenval["labels"] ])
                    coords/=A2m
                    #print(coords)
                    plt.xticks(coords,labels)
                #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                #plt.rc('text', usetex=True)
                axes = plt.gca()
                axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                #axes.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                #y_formatter = ScalarFormatter(useOffset=False)
                #axes.yaxis.set_major_formatter(y_formatter)
                y_min=-10
                y_max=10
                y_min*=eV2J
                y_max*=eV2J
                axes.set_xlim([x_full[0],x_full[-1]])
                axes.set_ylim([y_min,y_max])
                #plt.getp(axes.spines.values()) #default linewidth is 0.8, can see by uncommenting here
                default_lw=0.8
                #ymin, ymax = axes.get_ylim()
                for kval in self.eigenval["special_hs_kpoints_coords"][1:-1]:
                    kval/=A2m
                    #print(kval)
                    plt.axvline(x=kval,color='k',linewidth=default_lw,snap=False)
                    #plt.plot((kval, kval), (ymin, ymax), 'k-')
                plt.axhline(y=0,color='k',linewidth=default_lw,snap=False)
                #if getSIUnits:
                plt.ylabel("Energy (J)")
                #else:
                #    plt.ylabel("Energy (eV)")
                #axes.spines.values()
                #print(plt.getp(axes.spines.values()))
                #print(plt.getp(axes.spines.left,'linewidth'))
                #print(axes.spines.values() #rcParams['axes.linewidth'])
                if spin_count==1:
                    plt.savefig(self.name+".pdf")
                else:
                    plt.savefig(self.name+"_spin_"+str(i_spin)+".pdf")
                #plt.show()
                plt.clf()
                plt.cla()
                plt.close()

            #print("eff_mass:",h_bar**2*1.0/(second_der) #/(2*np.pi))
            #print("eff_mass_electron:",h_bar**2*1.0/(second_der)/mass_e #/(2*np.pi))
        
        return output

	#start here!
    def showsLTEDescriptor(self,index_region=9,plot=True,plot_dir=".",i_spin=0):
        if index_region%2==0:
            sys.exit("ERROR: index_region must be odd!")
        if index_region<=3:
            sys.exit("ERROR: index_region must be >3!")
        isLTE=False
        include_labels=True
        #index_region=5    #better if odd, min/max is in the middle
        if plot:
            if not hasattr(self,"kpoints"):
                print("If you want labels, parse KPOINTS file")
                include_labels=False
        self.getKpointPathVector()
        kpoints=self.eigenval["kpoints_plot"]
        spin_count=(2 if self.eigenval["spin_polarized"] else 1)
        if i_spin>spin_count-1:
            print(i_spin,spin_count-1)
            sys.exit("ERROR: requesting spin-polarized of a non-spin-polarized calculation")
        print("for spin",i_spin)
        BS_points=self.eigenval["BS_points_plot"][:,i_spin,:]   #kpoints,spin,energy
        halves=index_region/2
        quarters=halves/2

        VERBOSE=True
        restrict_hs_points=True
        getSIUnits=False

        threshold_first_derivative=0.1  #close to 0, must be SMALLER than this
        threshold_second_derivative=10 #1.0 #5 #15.0 #10.0 #1.5 #0.0 #1.5 #3.0 #1.0 #significant curves, must be BIGGER than this
        #orginally 15, 5 works better to catch Natalio's request:  RHL/As2Na1Sn2_ICSD_82366
        #originally 5, 1 works better to catch Natalio's requeset: MCLC/Re4Si7_ICSD_151529

        check_first_derivative=False

        ##########################################
        #COREY NEW MODS 180409
        check_if_within_1eV=False           #old procedure
        check_if_within_300meV_vbe=True     #new procedure, check if within 300meV of the valence band edge
        energies_found=[]
        second_derivatives_found=[]
        ##########################################

        for ei in range(BS_points.shape[1]):    #through each band
            for continuous_range in self.kpoints["ranges_of_continuity"]:
                range_start=continuous_range[0]+halves
                range_end=continuous_range[1]-halves
                #print(continuous_range,halves,"range",range_start,range_end)
                for ki in range(range_start,range_end):
                    if restrict_hs_points and ki not in self.kpoints["high_symmetry_indices"]:
                        continue
                    
                    ##############################################################
                    #COREY NEW MODS 180409
                    #OLD PROCEDURE
                    if check_if_within_1eV and abs(BS_points[ki,ei])>1.0: #np.any(np.abs(y_full_range)>1.0):	#only check point for within -1:1, not full range (only used for curvature)
                        continue
                    
                    #NEW PROCEDURE
                    if check_if_within_300meV_vbe and (BS_points[ki,ei]>0 or BS_points[ki,ei]<-0.3):
                        continue
                    ##############################################################

                    if False:
                        x_full_range=copy.deepcopy(self.eigenval["kpath"][ki-halves:ki+halves+1])
                        y_full_range=copy.deepcopy(BS_points[ki-halves:ki+halves+1,ei])
                    else:
                        x_full_range=self.eigenval["kpath"][ki-halves:ki+halves+1]
                        y_full_range=BS_points[ki-halves:ki+halves+1,ei]
                    #print(ki-halves,ki+halves+1)
                    #print(y_full_range)
                    #print(self.kpoints["ranges_of_continuity"])
                    
                    if False:
                        x_full_range_OLD=copy.deepcopy(x_full_range)
                        y_full_range_OLD=copy.deepcopy(y_full_range)

                    x_full_range,y_full_range,\
                            x_left,y_left,\
                            x_right,y_right=getMonotonicRegionLTE(x_full_range,y_full_range)
                    
                    #x_full_range,y_full_range,\
                    #        x_middle,y_middle,\
                    #        x_left,y_left,\
                    #        x_right,y_right=getMonotonicRegionLTE(x_full_range,y_full_range)

                    if x_full_range is None:
                        continue

                    if False:
                        x_full_range_NEW=copy.deepcopy(x_full_range)
                        y_full_range_NEW=copy.deepcopy(y_full_range)

                    #print(BS_points[ki-halves:ki+halves+1,ei],ki-halves,ki+halves+1)
                    #x_point_of_interest=copy.deepcopy(self.eigenval["kpath"][ki])
                    #y_point_of_interest=copy.deepcopy(BS_points[ki,ei])
                    #x_middle=self.eigenval["kpath"][ki-quarters:ki+quarters+1] #copy.deepcopy(self.eigenval["kpath"][ki-quarters:ki+quarters+1])
                    #y_middle=BS_points[ki-quarters:ki+quarters+1,ei] #copy.deepcopy(BS_points[ki-quarters:ki+quarters+1,ei])
                    
                    #if check_first_derivative:
                    #    z=np.polyfit(x_middle,y_middle,1)
                    #    fop=z[0]    #first order poly
                    #    #print(fop)
                    #    if abs(fop)>threshold_first_derivative:
                    #        continue
                    #    #print(y_middle,fop)
                    
                    #x_left=copy.deepcopy(self.eigenval["kpath"][ki-halves:ki])
                    #y_left=copy.deepcopy(BS_points[ki-halves:ki,ei])
                    #x_right=copy.deepcopy(self.eigenval["kpath"][ki+1:ki+halves+1])
                    #y_right=copy.deepcopy(BS_points[ki+1:ki+halves+1,ei])
                    
                    #convert units
                    if getSIUnits:
                        x_left/=A2m
                        x_right/=A2m
                        y_left*=eV2J
                        y_right*=eV2J
                    
                    z_left=np.polyfit(x_left,y_left,2)
                    sop_left=z_left[0]    #second order poly
                    z_right=np.polyfit(x_right,y_right,2)
                    sop_right=z_right[0]    #second order poly
                    second_derivative_left=2.0*sop_left
                    second_derivative_right=2.0*sop_right

                    if abs(second_derivative_left)<threshold_second_derivative or abs(second_derivative_right)<threshold_second_derivative:
                        continue

                    if np.sign(second_derivative_left)!=np.sign(second_derivative_right):
                        #print(y_left,y_right)
                        #print(second_derivative_left,second_derivative_right)
                        #print("we found one!")
                        if VERBOSE:
                            #print("OLD",x_full_range_OLD,y_full_range_OLD)
                            #print("NEW",x_full_range_NEW,y_full_range_NEW)
                            #print("possible LTE:",y_full_range,"fop:",fop,"sop-:",second_derivative_left,"sop+:",second_derivative_right)
                            print("left",x_left,y_left)
                            print("right",x_right,y_right)
                            print("2ndDerivative-:",second_derivative_left,"2ndDerivative+:",second_derivative_right)
                            #if len(x_left)<7:
                            #    print(x_left,y_left)
                            #    print(x_right,y_right)
                        isLTE=True
                        energies_found.append(BS_points[ki,ei])
                        second_derivatives_found.append(second_derivative_left)
                        second_derivatives_found.append(second_derivative_right)
                        if plot:
                            if getSIUnits:
                                x_full_range/=A2m
                                y_full_range*=eV2J
                            
                            p_left=np.poly1d(z_left)
                            p_right=np.poly1d(z_right)
                            xp_left=np.linspace(x_left[0],x_left[-1],100)
                            xp_right=np.linspace(x_right[0],x_right[-1],100)
                            plt.plot(xp_left,p_left(xp_left),"-",linewidth=4,color=[0,0,1])
                            plt.plot(xp_right,p_right(xp_right),"-",linewidth=4,color=[0,0,1])
                            

        if isLTE:
            if plot:
                if getSIUnits:
                    x_full=copy.deepcopy(self.eigenval["kpath"])
                    x_full/=A2m
                else:
                    x_full=self.eigenval["kpath"]
                #print(x_full)
                for ie in range(self.eigenval["n_energy"]):
                    if getSIUnits:
                        y_full=copy.deepcopy(BS_points[:,ie])
                        y_full*=eV2J
                    else:
                        y_full=BS_points[:,ie]
                    #print(y_full)
                    plt.plot(x_full,y_full,"-",linewidth=0.5)
                if include_labels:
                    labels=[ label["label"] for label in self.eigenval["labels"] ]
                    coords=np.asarray([ label["coord"] for label in self.eigenval["labels"] ])
                    if getSIUnits:
                        coords/=A2m
                    #print(coords)
                    plt.xticks(coords,labels)
                #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                #plt.rc('text', usetex=True)
                axes = plt.gca()
                axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                #axes.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                #y_formatter = ScalarFormatter(useOffset=False)
                #axes.yaxis.set_major_formatter(y_formatter)
                y_min=-4
                y_max=4
                if getSIUnits:
                    y_min*=eV2J
                    y_max*=eV2J
                axes.set_xlim([x_full[0],x_full[-1]])
                axes.set_ylim([y_min,y_max])
                #plt.getp(axes.spines.values()) #default linewidth is 0.8, can see by uncommenting here
                default_lw=0.8
                #ymin, ymax = axes.get_ylim()
                for kval in self.eigenval["special_hs_kpoints_coords"][1:-1]:
                    if getSIUnits:
                        kval/=A2m
                    #print(kval)
                    plt.axvline(x=kval,color='k',linewidth=default_lw,snap=False)
                    #plt.plot((kval, kval), (ymin, ymax), 'k-')
                plt.axhline(y=0,color='k',linewidth=default_lw,snap=False)
                if getSIUnits:
                    plt.ylabel("Energy (J)")
                else:
                    plt.ylabel("Energy (eV)")
                plt.title(self.name)
                #axes.spines.values()
                #print(plt.getp(axes.spines.values()))
                #print(plt.getp(axes.spines.left,'linewidth'))
                #print(axes.spines.values() #rcParams['axes.linewidth'])
                mkdir_p(plot_dir)
                file_name=""
                if spin_count==1:
                    file_name=os.path.join(plot_dir,self.name.replace("/","_")+"_band_structure.pdf")
                    #plt.savefig("test.pdf")#self.name+".pdf")
                else:
                    file_name=os.path.join(plot_dir,self.name.replace("/","_")+"_spin_"+str(i_spin)+"_band_structure.pdf")
                    #plt.savefig("test.pdf")#self.name+"_spin_"+str(i_spin)+".pdf")
                plt.savefig(file_name)
                #plt.show()
                plt.clf()
                plt.cla()
                plt.close()
                return (True,max(energies_found))
        return (False,None)

    def plotBands(self,i_spin=0):
        include_labels=True
        self.getKpointPathVector()
        if not hasattr(self,"kpoints"):
            print("If you want labels, parse KPOINTS file")
            include_labels=False
        kpoints=self.eigenval["kpoints_plot"]
        spin_count=(2 if self.eigenval["spin_polarized"] else 1)
        if i_spin>spin_count-1:
            sys.exit("ERROR: requesting spin-polarized of a non-spin-polarized calculation")
        #for i_spin in range(spin_count):
        BS_points=self.eigenval["BS_points_plot"][:,i_spin,:]
        x_full=copy.deepcopy(self.eigenval["kpath"])
        x_full/=A2m
        #print(x_full)
        for ie in range(self.eigenval["n_energy"]):
            y_full=copy.deepcopy(BS_points[:,ie])
            y_full*=eV2J
            #print(y_full)
            plt.plot(x_full,y_full,"-",linewidth=0.5)
        if include_labels:
            labels=[ label["label"] for label in self.eigenval["labels"] ]
            coords=np.asarray([ label["coord"] for label in self.eigenval["labels"] ])
            coords/=A2m
            #print(coords)
            plt.xticks(coords,labels)
        #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.rc('text', usetex=True)
        axes = plt.gca()
        axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        #axes.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        #y_formatter = ScalarFormatter(useOffset=False)
        #axes.yaxis.set_major_formatter(y_formatter)
        y_min=-10
        y_max=10
        y_min*=eV2J
        y_max*=eV2J
        axes.set_xlim([x_full[0],x_full[-1]])
        axes.set_ylim([y_min,y_max])
        #plt.getp(axes.spines.values()) #default linewidth is 0.8, can see by uncommenting here
        default_lw=0.8
        #ymin, ymax = axes.get_ylim()
        for kval in self.eigenval["special_hs_kpoints_coords"][1:-1]:
            kval/=A2m
            #print(kval)
            plt.axvline(x=kval,color='k',linewidth=default_lw,snap=False)
            #plt.plot((kval, kval), (ymin, ymax), 'k-')
        plt.axhline(y=0,color='k',linewidth=default_lw,snap=False)
        plt.ylabel("Energy (J)")
        #axes.spines.values()
        #print(plt.getp(axes.spines.values()))
        #print(plt.getp(axes.spines.left,'linewidth'))
        #print(axes.spines.values() #rcParams['axes.linewidth'])
        if spin_count==1:
            plt.savefig(self.name+".pdf")
        else:
            plt.savefig(self.name+"_spin_"+str(i_spin)+".pdf")
        #plt.show()
        plt.clf()
        plt.cla()
        plt.close()

    def __parseCPF(self):
        #initialize local variables
        critical_points=[]

        lines=self.CPF_lines
        
        #set search booleans and default vals
        search4CoordI=False
        foundCoordI=False
        search4EigenV=False
        foundEigenV=False
        setCP=False
        coordI=np.empty(3)
        coordF=np.empty(3)
        coordC=np.empty(3)
        eigenV=np.empty(3)

        for index,line in enumerate(lines):
            if "A NEW ENTRY" in line:
                if index and not setCP:
                    print(": ".join(["ERROR","Issue extracting full critical point information on line " % index]))
                    self.created=False
                    return None
                coordI=np.empty(3)
                coordF=np.empty(3)
                coordC=np.empty(3)
                eigenV=np.empty(3)
                search4CoordI=False
                foundCoordI=False
                search4EigenV=False
                foundEigenV=False
                setCP=False
            if foundCoordI and foundEigenV:
                coords=getCoordsCHGCAR(coordI,self)
                try:
                    coordF=coords["coordF"]
                    coordC=coords["coordC"]
                except:
                    print(": ".join(["ERROR","Issue extracting coordF/coordC near line %i" % index]))
                    self.created=False
                    return None
                cp=criticalPoint(coordI,coordF,coordC,eigenV)
                #print("coordI=",cp.coordI," ngrid=",self.ngrid)
                #print("coordF=",cp.coordF)
                #print("coordC=",cp.coordC)
                #print("eigenV=",cp.eigenV)
                #print("type:",cp.point_type)
                #print("")
                critical_points.append(cp)
                foundCoordI=False
                foundEigenV=False
                setCP=True
            if not foundCoordI:
                if search4CoordI:
                    try:
                        parts=[ int(i) for i in line.split() ]
                        coordI[0]=parts[0]
                        coordI[1]=parts[1]
                        coordI[2]=parts[2]
                        #print("pre:",coordI)
                        coordI=convertIndex12Index0(coordI,self.ngrid)  #get correct index
                        #print("post:",coordI)
                        #print(self.charges[coordI[0]][coordI[1]][coordI[2]])
                    except:
                        print(": ".join(["ERROR","Issue extracting coordI on line %i" % index]))
                        self.created=False
                        return None
                    search4CoordI=False
                    foundCoordI=True
                else:
                    if "Critical point number:" in line:
                        search4CoordI=True
            if not foundEigenV:
                if search4EigenV:
                    try:
                        parts=[ float(i) for i in line.split() ]
                        if any([ math.isnan(i) for i in parts ]):
                            print(": ".join(["WARNING","Found NaN on line %i" % index]))
                            print(line)
                            search4EigenV=False     #sabotage the adding of this CP
                            setCP=True              #allow it to pass like nothing
                            continue
                        eigenV[0]=parts[0]
                        eigenV[1]=parts[1]
                        eigenV[2]=parts[2]
                    except:
                        print(": ".join(["ERROR","Issue extracting eigenV on line %i" % index]))
                        self.created=False
                        return None
                    search4EigenV=False
                    foundEigenV=True
                else:
                    if "Eigenvalues:" in line:
                        search4EigenV=True

        if verifyCPInternal and not verifyCP(self,critical_points):
            print(": ".join(["ERROR","Issue with critical points"]))
            self.created=False
            return None

        #save local variables to object
        self.critical_points=critical_points

    def plotStruct(self):
        x=np.empty(self.num_atoms)
        y=np.empty(self.num_atoms)
        z=np.empty(self.num_atoms)
        for i,atom in enumerate(self.atoms):
            print(atom.coordF,atom.coordC)
            x[i]=atom.coordC[0]
            y[i]=atom.coordC[2]
            z[i]=atom.coordC[2]
        #print(x)
        #print(y)
        #print(z)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig=plt.figure()
        ax=fig.gca(projection='3d')
        ax.scatter(x,y,z,s=40)
        plt.show()

    def __repr__(self):
        string="<structure: "
        string+=str(vars(self))
        string+=">"
        return string
    def __str__(self):
        if self.name:
            return "%s" % (self.name)
        #elif False and self.title:
        elif self.title:
            return "%s" % (self.title)
        else:
            string="<structure: "
            string+=str(vars(self))
            string+=">"
            return string

class atom:
    #def __init__(self,species,atom_type,coordF,coordC):
    def __init__(self,atom_type,species,coordF,coordC):
        self.species=""
        self.atom_type=-1
        self.coordF=np.empty(3)
        self.coordC=np.empty(3)

        #self.species=species
        self.atom_type=atom_type
        self.species=species
        self.coordF=coordF
        self.coordC=coordC
    def __repr__(self):
        string="<atom: "
        string+=str(vars(self))
        string+=">"
        return string
    def __str__(self):
        return "%s" % (self.atom_type)

class criticalPoint:
    #analyze type in here
    def __init__(self,coordI,coordF,coordC,eigenV):
        self.coordI=np.empty(3)
        self.coordF=np.empty(3)
        self.coordC=np.empty(3)
        self.eigenV=np.empty(3)
        self.point_type=""

        self.coordI=coordI
        self.coordF=coordF
        self.coordC=coordC
        self.eigenV=eigenV
        self.point_type=CP_inverse[np.sum(np.sign(eigenV))]
    def __repr__(self):
        string="<critical point: "
        string+=str(vars(self))
        string+=">"
        return string
    def __str__(self):
        return "%s" % (self.point_type)
        

#make methods for these instead
#def chg:
#    #def __init__(self,coordI,value,coordF,coordC):
#    def __init__(self,value,coordF,coordC):
#        #self.coordI=np.empty(3)
#        self.value=0.0
#        self.coordF=np.empty(3)
#        self.coordC=np.empty(3)
#
#        #self.coordI=coordI
#        self.value=value
#        self.coordF=coordF
#        self.coordC=coordC
#    def __repr__(self):
#        string="<chg: "
#        string+=str(vars(self))
#        string+=">"
#        return string
#    def __str__(self):
#        return "%s" % (self.value)

def bringInCell(coordF):
    newCoordF=np.empty(coordF.size)
    inCellEPS=1.0-zeroTol
    for i in range(coordF.size):
        c=coordF[i]
        while c>inCellEPS:
            c-=1.0
        while c<0.0:
            c+=1.0
        if abs(c)<zeroTol:
            c=0.0
        if c>inCellEPS:
            c=0.0
        newCoordF[i]=c
    return newCoordF

if index0:
    def getCoordCCHGCAR(coordI,ngrid,lattice_vectors):
        #starts at 0!
        #print("1:",lattice_vectors[0,:],0.5*lattice_vectors[0,:])
        #print("2:",lattice_vectors[1,:])
        #print("3:",lattice_vectors[2,:])
        coordC=float(coordI[0])/float(ngrid[0]-1)*lattice_vectors[0,:]+\
                float(coordI[1])/float(ngrid[1]-1)*lattice_vectors[1,:]+\
                float(coordI[2])/float(ngrid[1]-1)*lattice_vectors[2,:]
        return coordC
else:
    def getCoordCCHGCAR(coordI,ngrid,lattice_vectors):
        #starts at 1!
        coordC=float(coordI[0])/float(ngrid[0])*lattice_vectors[0,:]+\
                float(coordI[1])/float(ngrid[1])*lattice_vectors[1,:]+\
                float(coordI[2])/float(ngrid[1])*lattice_vectors[2,:]
        return coordC

if index0:
    def getCoordFCHGCAR(coordI,ngrid):
        #starts at 0!
        coordF=np.array([float(coordI[0])/float(ngrid[0]-1),float(coordI[1])/float(ngrid[1]-1),float(coordI[2])/float(ngrid[2]-1)])
        return bringInCell(coordF)
else:
    def getCoordFCHGCAR(coordI,ngrid):
        #starts at 1!
        coordF=np.array([float(coordI[0])/float(ngrid[0]),float(coordI[1])/float(ngrid[1]),float(coordI[2])/float(ngrid[2])])
        return bringInCell(coordF)

if getCoordC:
    def getCoordsCHGCAR(coordI,struct):
        coordC=getCoordCCHGCAR(coordI,struct.ngrid,struct.lattice_vectors)
        coordF=bringInCell(np.dot(struct.c2f,coordC))
        coordC=np.dot(struct.f2c,coordF)    #adjust for bringInCell
        return {"coordF":coordF,"coordC":coordC}
else:
    def getCoordsCHGCAR(coordI,struct):
        coordF=getCoordFCHGCAR(coordI,struct.ngrid)
        coordC=np.dot(struct.f2c,coordF)
        return {"coordF":coordF,"coordC":coordC}


if testConvertIndex12Index0:
    def convertIndex12Index0(coordI,ngrid):
        nChg=(coordI[0]-1)+(coordI[1]-1)*ngrid[0]+(coordI[2]-1)*ngrid[1]*ngrid[2]
        print("nChg:",nChg)
        i_charge=0      #for charge indexing
        indexX=0;indexY=0;indexZ=0  #for grid indexing
        while i_charge<nChg+1:
            if indexX and indexX%ngrid[0]==0:
                indexX=0
                indexY+=1
            if indexY and indexY%ngrid[1]==0:
                indexY=0
                indexZ+=1
            #print(indexX,indexY,indexZ)
            indexX+=1
            i_charge+=1
        indexX-=1
        return np.array([indexX,indexY,indexZ])
else:
    def convertIndex12Index0(coordI,ngrid):
        return np.array([coordI[0]-1,coordI[1]-1,coordI[2]-1])

def getThreshold(ngrid,lattice_vectors):
    return 1.5*max( [ modV(1.0/ngrid[i]*lattice_vectors[i]) for i in range(3) ] )

def getThresholdS(struct):
    return getThreshold(struct.ngrid,struct.lattice_vectors)

def samePositions(coord1,coord2,tol=0.2):
    #tol=0.2
    try:
        return distCoords(coord1,coord2)<tol 
        #return all([ abs(coord1[i]-coord2[i])<tol for i in range(coord1.size) ])
    except:
        return False

def modV(vec):
    return np.linalg.norm(vec,ord=2)

def angleV(v1,v2):
    return np.arccos(np.dot(v1,v2)/(modV(v1)*modV(v2)))

def distCoords(coord1,coord2):
    return abs(modV(coord1-coord2))

def distCoordsPBC(lattice_vectors,coord1,coord2):
    minDist=distCoords(coord1,coord2)
    #print(minDist)
    for i in range(-2,2):
        for j in range(-2,2):
            for k in range(-2,2):
                #print(abs(modV(coord1-coord2+float(i)*lattice_vectors[0,:]+float(j)*lattice_vectors[1,:]+float(k)*lattice_vectors[2,:])))
                minDist=min(minDist,abs(modV(coord1-coord2+float(i)*lattice_vectors[0,:]+float(j)*lattice_vectors[1,:]+float(k)*lattice_vectors[2,:])) )
    return minDist

def getTolRatio(struct):
    tol=getThresholdS(struct)
    tolRatios=[]
    indicesAssociated=[]
    if struct.critical_points:
        for indexCP,cp in enumerate(struct.critical_points):
            if cp.point_type!=NA:
                continue
            minCPIndex=0
            minDist=0
            minSet=False
            print("coordI=",cp.coordI," ngrid=",struct.ngrid)
            print("coordF=",cp.coordF)
            print("coordC=",cp.coordC)
            print("eigenV=",cp.eigenV)
            print("type:",cp.point_type)
            print("")
            for indexAtom,atom in enumerate(struct.atoms):
                print("atom=",indexAtom)
                print("coordF=",atom.coordF)
                print("coordC=",atom.coordC)
                newDist=distCoordsPBC(struct.lattice_vectors,atom.coordC,cp.coordC)
                if not minSet or newDist<minDist:
                    minCPIndex=indexCP
                    minDist=newDist
                    minSet=True
            if minSet:
                print("minCPIndex=",minCPIndex,"minDist=",minDist,"tol=",tol,"minDist/tol=",minDist/tol)
                tolRatios.append(minDist/tol)

        #for indexAtom,atom in enumerate(struct.atoms):
        #    minCPIndex=0
        #    minDist=0
        #    minSet=False
        #    for indexCP,cp in enumerate(struct.critical_points):
        #        if indexCP in indicesAssociated:
        #            continue
        #        ####
        #        #MUST FIGURE OUT WRAP AROUND DISTANCES
        #        ####
        #        newDist=distCoordsPBC(struct.lattice_vectors,atom.coordC,cp.coordC)
        #        #print("atom:",indexAtom,"cp:",indexCP,"dist:",newDist)
        #        if not minSet or newDist<minDist:
        #            minCPIndex=indexCP
        #            minDist=newDist
        #            minSet=True
        #    if minSet and struct.critical_points[minCPIndex].point_type==NA:
        #        #print("fract:",atom.coordF,cp.coordF)
        #        #print(minCPIndex,minDist,tol,minDist/tol)
        #        tolRatios.append(minDist/tol)
    else:
        print(": ".join(["WARNING","Found no critical points"]) )
    if tolRatios:
        #print(tolRatios)
        return max(tolRatios)
    else:
        print(": ".join(["WARNING"," ".join(["Found no",NA])]) )
        return 0.0


def verifyCP(struct,critical_points):
    tol=getThresholdS(struct)
    #indicesAssociated=[]
    print(len([ cp for cp in critical_points if cp.point_type==NA ]))
    for indexAtom,atom in enumerate(struct.atoms):
        outString="".join(["atom=",str(indexAtom)])
        minCPIndex=0
        minDist=0
        minSet=False
        for indexCP,cp in enumerate(critical_points):
            #if cp.point_type!=NA or indexCP in indicesAssociated:
            #if cp.point_type!=NA:
            #    continue
            newDist=distCoordsPBC(struct.lattice_vectors,atom.coordC,cp.coordC)
            if not minSet or newDist<minDist:
                minCPIndex=indexCP
                minDist=newDist
            if not indexAtom:
                print(cp)
        #indicesAssociated.append(minCPIndex)
        outString="".join([outString,", cp=",str(minCPIndex),", minDist=",str(minDist),", point_type=",critical_points[minCPIndex].point_type,", eigenV=",str(critical_points[minCPIndex].eigenV)])
        print(outString)

            #if samePositions(cp.coordF,atom.coordF):
            #    print(cp.coordF,atom.coordF)
            #    foundCP=True
            #    break
        #if not foundCP:
        #    print(": ".join(["ERROR","".join(["Did not find ",NA," for atom ",str(index)])]))
        #    for cp in critical_points:
        #        print(cp.coordF," // ",atom.coordF,", mod=",str(abs(np.linalg.norm(cp.coordF-atom.coordF,ord=2))))
        #    return False

    return True

def checkConvergenceOSZICAR(oszicar,steps=60):
    lines=read_file_lines(oszicar)
    if lines is None:
        sys.exit(oszicar+" file is empty")
    #check file finished writing
    if not lines:
        sys.exit(oszicar+" file is empty")
    last_line=lines[-1]
    if ("F=" not in last_line) or ("E0=" not in last_line) or ("d E =" not in last_line):
        return (False,"not finished writing")
    last_iteration=lines[-2]
    last_step=last_iteration.split(":")[1].split()[0].strip()
    try:
        last_step=int(last_step)
    except:
        if last_step=="***":
            last_step=0
            for line in lines:
                line=line.strip()
                if not line:
                    continue
                if "N" in line and "E" in line and "dE" in line and "d eps" in line and "ncg" in line and "rms" in line and "rms(c)" in line:
                    last_step=0
                    continue
                if "F=" in line and "E0=" in line and "d E =" in line:
                    pass
                else:
                    last_step+=1
            try:
                last_step=int(last_step)
            except:
                sys.exit(oszicar+" odd format, cannot get last_iteration_step")
    if last_step>=steps:
        return (False,"unconverged")
    return (True,"converged")

def getEnergyOUTCAR(outcar_file=None,outcar_lines=None):
    lines=None
    if lines is None and outcar_lines is not None:
        lines=outcar_lines
    if lines is None and outcar_file is not None:
        lines=read_file_lines(outcar_file)
    if not lines:
        sys.exit(outcar+" file is empty")
    energy_line=None
    for line in lines:
        if "energy without entropy" in line:
            energy_line=line
    if energy_line is None:
        sys.exit(outcar+" could not find energy line")
    energy_wo_entropy=energy_line.split("=")[1].split()[0]
    try:
        energy_wo_entropy=float(energy_wo_entropy)
    except:
        sys.exit(outcar+" odd format, cannot get energy_wo_entropy")
    return energy_wo_entropy

def checkConvergenceOUTCAR(outcar,check_neg_energy=True,check_full_occupancy=False,check_pressure=True):
    lines=read_file_lines(outcar)
    if lines is None:
        sys.exit(outcar+" file is empty")
    #check file finished writing
    if not lines:
        sys.exit(outcar+" file is empty")
    last_line=lines[-1]
    #print(last_line)
    if "Voluntary context switches" not in last_line:
        return (False,"not finished writing")
    #EENTRO must be less than 1 meV/atom
    #https://sites.google.com/site/notesinmaterialsscience/vasp
    eentro_line=None
    #energy_line=None
    pressure_line=None
    for line in lines:
        if "EENTRO" in line:
            eentro_line=line
        #if "energy without entropy" in line:
        #    energy_line=line
        if "in kB" in line:
            pressure_line=line
    if eentro_line is None:
        sys.exit(outcar+" could not find eentro line")
    eentro_energy=eentro_line.split()[-1]
    try:
        eentro_energy=float(eentro_energy)
    except:
        sys.exit(outcar+" odd format, cannot get eentro_energy")
    #print(eentro_energy)
    if eentro_energy>0.001:
        return (False,"unphysical EENTRO")
    energy_wo_entropy=getEnergyOUTCAR(outcar_file=outcar,outcar_lines=lines)
    #if energy_line is None:
    #    sys.exit(outcar+" could not find energy line")
    #energy_wo_entropy=energy_line.split("=")[1].split()[0]
    #try:
    #    energy_wo_entropy=float(energy_wo_entropy)
    #except:
    #    sys.exit(outcar+" odd format, cannot get energy_wo_entropy")
    if check_neg_energy and energy_wo_entropy>0:
        return (False,"positive energy")
    if check_full_occupancy:  #atoms only
        found_band_occs=False
        for line in lines:
            line=line.strip()
            #print(line,(not line))
            if found_band_occs:
                if not line:
                    found_band_occs=False
                    continue
                parts=line.split()
                try:
                    index=int(parts[0])
                    energy=float(parts[1])
                    occ=float(parts[2])
                except:
                    print(outcar)
                    sys.exit("odd format for bands")
                #print(energy,occ)
                if (not floatsEqual(occ,0.0)) and (not floatsEqual(occ,1.0)):
                    return (False,"partial occupancy found")
            else:
                if "band No.  band energies     occupation" in line:
                    found_band_occs=True
                    continue
    if check_pressure:
        pressures=pressure_line.replace("in kB","").strip().split()
        try:
            pressures=[ abs(float(p)) for p in pressures ]
        except:
            sys.exit(outcar+" odd format, cannot get pressures")
        if len(pressures)!=6:
            sys.exit(outcar+" odd format, not all stress tensor components provided")
        #print(pressures)
        if max(pressures)>10:
            return (False,"pressure too high max(P)="+str(max(pressures)))
    return (True,"converged")

def getNELM(incar):
    lines=read_file_lines(incar)
    if lines is None:
        sys.exit(incar+" file is empty")
    #check file finished writing
    if not lines:
        sys.exit(incar+" file is empty")
    for line in lines:
        if "NELM=" in line:
            nelm=line.strip().split("=")[1]
            try:
                nelm=int(nelm)
            except:
                sys.exit(incar+" badly formatted NELM")
            return nelm
    return 60
