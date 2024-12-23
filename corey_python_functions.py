import os
import sys
import shutil
import errno
import subprocess
#import bz2
#import gzip
#import zipfile
import math
import re
import json

ELEMENTS=["H" ,"He" ,"Li" ,"Be" ,"B" ,"C" ,"N" ,"O" ,"F" ,"Ne" ,"Na" ,"Mg" ,"Al"\
,"Si" ,"P" ,"S" ,"Cl" ,"Ar" ,"K" ,"Ca" ,"Sc" ,"Ti" ,"V" ,"Cr" ,"Mn" ,"Fe"\
,"Co" ,"Ni" ,"Cu" ,"Zn" ,"Ga" ,"Ge" ,"As" ,"Se" ,"Br" ,"Kr" ,"Rb" ,"Sr"\
,"Y" ,"Zr" ,"Nb" ,"Mo" ,"Tc" ,"Ru" ,"Rh" ,"Pd" ,"Ag" ,"Cd" ,"In" ,"Sn"\
,"Sb" ,"Te" ,"I" ,"Xe" ,"Cs" ,"Ba" ,"La" ,"Ce" ,"Pr" ,"Nd" ,"Pm" ,"Sm"\
,"Eu" ,"Gd" ,"Tb" ,"Dy" ,"Ho" ,"Er" ,"Tm" ,"Yb" ,"Lu" ,"Hf" ,"Ta" ,"W"\
,"Re" ,"Os" ,"Ir" ,"Pt" ,"Au" ,"Hg" ,"Tl" ,"Pb" ,"Bi" ,"Po" ,"At" ,"Rn"\
,"Fr" ,"Ra" ,"Ac" ,"Th" ,"Pa" ,"U" ,"Np" ,"Pu" ,"Am" ,"Cm" ,"Bk" ,"Cf"\
,"Es" ,"Fm" ,"Md" ,"No" ,"Lr" ,"Rf" ,"Db" ,"Sg" ,"Bh" ,"Hs" ,"Mt" ,"Ds"\
,"Rg" ,"Uub" ,"Uut" ,"Uuq" ,"Uup" ,"Uuh" ,"Uus" ,"Uuo"]

def replace(_file,pattern,substitution):
	#replaces pattern with substitution in file
	with open(_file,"r") as fin:
		text = fin.read()
	new_text = text.replace(pattern,substitution)
	with open(_file,"w") as fout:
		fout.write(new_text)

#def getElements(compound):
def getElements(compound,no_composition=True,remove_pp=True):
    compound=compound.split(":")[0]
    if no_composition and remove_pp:
        parts=re.findall('[A-Z][a-z]*', compound)
        return parts
    parts=re.findall('[A-Z][^A-Z]*', compound)
    #print(parts)
    #sys.exit()
    if no_composition:
        new_parts=[]
        for part in parts:
            if remove_pp:
                part=part.split(":")[0]
                part=part.split("_")[0]
            new_parts.append("".join([ i for i in part if not (i.isdigit() or i==".") ]))
        parts=new_parts
    #print(compound,parts)
    return parts

def getComposition(compound):
    return re.findall(r"[-+]?\d*\.\d+|\d+", compound)

def parseCompound(compound,sort=True):
    elements=[]
    composition=[]
    stoichiometry=[]
    if not compound:
        print("empty compound")
        return {"elements":elements,"composition":composition,"stoichiometry":stoichiometry}
    if not compound[0].isupper():
        print("elements not capitalized: "+compound)
        return {"elements":elements,"composition":composition,"stoichiometry":stoichiometry}
    current_element=""
    current_composition=""
    adding2element=True
    adding2composition=False
    for c in compound:
        if c.isalpha():
            if adding2composition:
                composition.append(float(current_composition))
                current_composition=""
                adding2composition=False
            if c.isupper() and current_element:
                elements.append(current_element)
                composition.append(1.0)
                current_element=""
            current_element+=c
            add2element=True
        else:
            if c.isdigit() or c==".":
                if add2element:
                    elements.append(current_element)
                    current_element=""
                    add2element=False
                current_composition+=c
                adding2composition=True
    if add2element and current_element:
        elements.append(current_element)
        composition.append(1.0)
        current_element=""
        add2element=False
    if adding2composition and current_composition:
        composition.append(float(current_composition))
        current_composition=""
        adding2composition=False
    if len(elements)!=len(composition):
        print("len(elements)!=len(composition): "+compound)
    if elements!=getElements(compound):
        print("elements!=getElements(compound): "+compound)
    total_comp=sum(composition)
    stoichiometry=[ float(i)/float(total_comp) for i in composition ]
    if sort:
        zipped=zip(elements,composition,stoichiometry)
        zipped.sort()
        elements,composition,stoichiometry=zip(*zipped)
    return {"elements":elements,"composition":composition,"stoichiometry":stoichiometry}

def content_empty(content):
    content=content.replace(" ","")
    content=content.replace("\t","")
    content=content.replace("\n","")
    return len(content)==0

def content_not_empty(content):
    return not content_empty

def issue_command(command):
    p=subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out,err=p.communicate()
    return out.decode('UTF-8'),err.decode('UTF-8')

def list_files(path):
    files=[]
    try:
        for _file in os.listdir(path):
            try:
                if os.path.isfile(os.path.join(path,_file)):
                    files.append(_file)
            except OSError as e:
                pass
    except OSError as e:
        pass
    return files
    #return [ _file for _file in os.listdir(path) if os.path.isfile(os.path.join(path,_file)) ]

def list_files_path(path):
    return [ os.path.join(path,_file) for _file in list_files(path) ]

def list_directories(path):
    directories=[]
    try:
        for directory in os.listdir(path):
            try:
                if os.path.isdir(os.path.join(path,directory)):
                    directories.append(directory)
            except OSError as e:
                pass
    except OSError as e:
        pass
    return directories
    #return [ directory for directory in os.listdir(path) if os.path.isdir(os.path.join(path,directory)) ]

def list_directories_path(path):
    return [ os.path.join(path,directory) for directory in list_directories(path) ]

#MATH
def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a if abs(a)>=1.0 else 1

def gcdm(*args):
    """Return gcd of args."""
    return reduce(gcd,args)

def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // gcd(a, b)

def lcmm(*args):
    """Return lcm of args."""
    return reduce(lcm, args)

def reduce_stoich(vec):
    return [ i/gcdm(*vec) for i in vec ]
#MATH

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

def split_path_prelim(p):
    # WARNING:  will return "" as first entry if path starts with /
    #           doesn't handle "/" at the end of the path well (os.path.split)
    a,b = os.path.split(p)
    return (split_path(a) if len(a) and len(b) else []) + [b]

def split_path(p):
    path=split_path_prelim(p)
    if not path[0]:
        path[0]=os.path.sep
    return path

def read_file(_filePath):
    if not os.path.exists(_filePath):
        return None
    if os.path.sep in _filePath:
        _file=split_path(_filePath)[-1]
    else:
        _file=_filePath
    vext=[".bz2",".xz",".gz",".zip",""] #stroke of genius
    vcommand=["bzcat","xzcat","zcat","unzip -p","cat"]
    for ie in range(len(vext)):
        if vext[ie] in _file:
            command=vcommand[ie]+" "+_filePath
            output,errors=issue_command(command)
            if content_not_empty(errors):
                return None
            fileContent=output
            break
    #if ".bz2" in _file:
    #    fin=bz2.BZ2File(_filePath,"r")
    #elif ".gz" in _file:
    #    fin=gzip.GzipFile(_filePath,"r")
    #elif ".zip" in _file:
    #    fin=zipfile.ZipFile(_filePath,"r")
    #else:
    #    try:
    #        fin=open(_filePath,"r")
    #    except:
    #        return None
    #fileContent=fin.read()
    if ".json" in _file:
        return json.loads(fileContent)
    return fileContent
    #if ".bz2" in _file:
    #    with bz2.BZ2File(_filePath,"r") as fin:
    #        return fin.read()
    #elif ".gz" in _file:
    #    with gzip.GzipFile(_filePath,"r") as fin:
    #        return fin.read()
    #elif ".zip" in _file:
    #    with zipfile.ZipFile(_filePath,"r") as fin:
    #        return fin.read()
    #else:
    #    try:
    #        with open(_filePath,"r") as fin:
    #            return fin.read()
    #    except:
    #        return None

def read_file_lines(_filePath):
    _fileContents=read_file(_filePath)
    if _fileContents is None:
        return None
    if ".json" in _filePath:
        return _fileContents
    return _fileContents.splitlines()

def read_file_lines_CLEAN(_filePath):
    _fileContentLines=read_file_lines(_filePath)
    if _fileContentLines is None:
        return None
    if ".json" in _filePath:
        return _fileContentLines
    _fileContentLines=_fileContentLines.strip()
    return [ content.strip() for content in _fileContentLines ]

def add_line_to_file(file,pattern,input,count,DEBUG=True):
    #adds input to file after count number of lines below pattern
    #count 1 is right after line
    #we automatically add endl
    #with open(file,"r") as fin:
    #	text = fin.read()
    #lines = text.strip().splitlines()
    lines=read_file_lines(file)
    for index,line in enumerate(lines):
        if index == len(lines) + 1 - count:
            break
        if pattern in line:
            lines.insert(index+count,input)
            break
    output="\n".join(map(str,lines))
    if DEBUG:
        print(output)
    else:
        with open(file,"w") as fout:
            fout.write(output)

#def check_files( file_path, outcars_to_check, oszicars_to_check ):
#    finishStringOUTCAR = "Voluntary context switches:"
#    for outcar in outcars_to_check:
#        outcar_path = file_path + "/" + outcar
#        #with bz2.BZ2File(outcar_path,"r") as fin:
#        #    lines = fin.readlines()
#        lines=read_file_lines(outcar_path)
#        if lines is None:
#            return("error","cannot open "+outcar)
#        count = 0
#        for index, line in enumerate(lines):
#            if line == "\n":
#                del lines[index]
#            if finishStringOUTCAR in line:
#                count += 1
#        if count == 0:
#            return("error",outcar + " did not finish writing")
#        elif count > 1:
#            return("error",outcar + " is corrupted--there is more than one end-line printed")
#        else:
#            if finishStringOUTCAR not in lines[-1]:
#                return("error",outcar + " is corrupted--last line looks suspicious")
#            continue
#
#    for oszicar in oszicars_to_check:
#        oszicar_path = file_path + "/" + oszicar
#        #with bz2.BZ2File(oszicar_path,"r") as fin:
#        #    lines = fin.readlines()
#        lines=read_file_lines(oszicar_path)
#        if lines is None:
#            return("error","cannot open "+oszicar)
#        for index, line in enumerate(lines):
#            if line == "\n":
#                del lines[index]
#        if ("F=" not in lines[-1]) or ("E0=" not in lines[-1]) or ("d E =" not in lines[-1]):
#            return("error",oszicar + " did not finish writing")
#        elif ": 120" in lines[-2]:
#            return("error",oszicar + " did not converge")
#        else:
#            continue
#    return("good","good")

#http://randlet.com/blog/python-significant-figures-format/
#sig figs with trailing 0s
def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

def generateError(soliloquy,message,Quit=True):
    full_message=" ".join(["ERROR in",soliloquy,"-",message])
    if Quit:
        sys.exit(full_message)
    else:
        print(full_message)

def generateMessage(soliloquy,message):
    full_message=" ".join(["MESSAGE from",soliloquy,"-",message])
    print(full_message)

class Options:
    """
    """
    def __init__(self,Argv=[],Prefix="--",Delimiter="="):
        self.argv=[]
        #one dict, {flag_name:{flag:True,scheme:""}}
        #self.flag_names[]   #actual name
        #self.flags=[]       #on/off
        #self.schemes=[]     #content associated with it
        self.flag_dict={}
        self.flag_template={"flag":True,"scheme":None}
        self.prefix=Prefix
        self.delimiter=Delimiter

        self.argv=Argv
        Options.__argv2Flag(self)

    def __argv2Flag(self):
        for arg in self.argv:
            Options.__arg2Flag(self,arg,flag=True)

    def __arg2Flag(self,arg,flag=True):
        soliloquy="Options.__args2Flag()"
        #arg="--no_abstract=True"
        if not self.prefix in arg:
            return None
        arg=arg.replace(self.prefix,"")
        if self.delimiter in arg:
            parts=arg.split(self.delimiter)
            if len(parts)!=2:
                generateError(soliloquy,"too many delimiters in flag input: "+arg)
            Options.addOption(self,Name=parts[0],Flag=flag,Scheme=parts[1])
        else:
            Options.addOption(self,Name=arg,Flag=flag)
        return None

    def addOption(self,Name="",Flag=True,Scheme=None):
        self.flag_dict[Name]=dict(self.flag_template)
        self.flag_dict[Name]["flag"]=Flag
        #be smart with scheme, maybe later
        self.flag_dict[Name]["scheme"]=Scheme

        return None

    def printOptions(self):
        for option in self.flag_dict.keys():
            print("{:<30}{:<30}{:<30}".format(option,self.flag_dict[option]["flag"],self.flag_dict[option]["scheme"]))

        return None

    def toggleOption(self,Name):
        soliloquy="Options.toggleOption()"
        try:
            self.flag_dict[Name]["flag"]=(not self.flag_dict[Name]["flag"])
        except:
            generateError(soliloquy,"flag not set yet: "+Name)

    def setFlag(self,Name,Flag):
        soliloquy="Options.setFlag()"
        try:
            self.flag_dict[Name]["flag"]=Flag
        except:
            generateError(soliloquy,"flag not set yet: "+Name)

    def getOptions(self):
        return self.flag_dict.keys()

    def optionOn(self,Name):
        return Options.setFlag(self,Name,True)

    def optionOff(self,Name):
        return Options.setFlag(self,Name,False)

    def setScheme(self,Name,Scheme):
        self.flag_dict[Name]["scheme"]=Scheme

    def detachScheme(self,Name):
        soliloquy="Options.detachScheme()"
        try:
            self.flag_dict[Name]["scheme"]=None
        except:
            generateError(soliloquy,"flag not set yet: "+Name)

    def removeScheme(self,Name):
        return detachScheme(self,Name)

    def flag(self,Name):
        try:
            return self.flag_dict[Name]["flag"]
        except:
            return False

    def scheme(self,Name):
        try:
            return self.flag_dict[Name]["scheme"]
        except:
            return None

def floatsEqual(a,b,tol=1e-8):
    return abs(a-b)<tol

def floatVectorsEqual(va,vb,tol=1e-8):
    if len(va)!=len(vb):
        return False
    for i in range(len(va)):
        if not floatsEqual(va[i],vb[i],tol):
            return False
    return True

#CRC STUFF
# https://gist.github.com/larsimmisch/9915766
# Portions Copyright (c) 1996-2001, PostgreSQL Global Development Group
# (Any use permitted, subject to terms of PostgreSQL license; see.)

# If we have a 64-bit integer type, then a 64-bit CRC looks just like the
# usual sort of implementation. (See Ross Williams' excellent introduction
# A PAINLESS GUIDE TO CRC ERROR DETECTION ALGORITHMS, available from
# ftp://ftp.rocksoft.com/papers/crc_v3.txt or several other net sites.)
# If we have no working 64-bit type, then fake it with two 32-bit registers.
#
# The present implementation is a normal (not "reflected", in Williams'
# terms) 64-bit CRC, using initial all-ones register contents and a final
# bit inversion. The chosen polynomial is borrowed from the DLT1 spec
# (ECMA-182, available from http://www.ecma.ch/ecma1/STAND/ECMA-182.HTM):
#
# x^64 + x^62 + x^57 + x^55 + x^54 + x^53 + x^52 + x^47 + x^46 + x^45 +
# x^40 + x^39 + x^38 + x^37 + x^35 + x^33 + x^32 + x^31 + x^29 + x^27 +
# x^24 + x^23 + x^22 + x^21 + x^19 + x^17 + x^13 + x^12 + x^10 + x^9 +

crc_table = [
    0x0000000000000000, 0x42F0E1EBA9EA3693,
    0x85E1C3D753D46D26, 0xC711223CFA3E5BB5,
    0x493366450E42ECDF, 0x0BC387AEA7A8DA4C,
    0xCCD2A5925D9681F9, 0x8E224479F47CB76A,
    0x9266CC8A1C85D9BE, 0xD0962D61B56FEF2D,
    0x17870F5D4F51B498, 0x5577EEB6E6BB820B,
    0xDB55AACF12C73561, 0x99A54B24BB2D03F2,
    0x5EB4691841135847, 0x1C4488F3E8F96ED4,
    0x663D78FF90E185EF, 0x24CD9914390BB37C,
    0xE3DCBB28C335E8C9, 0xA12C5AC36ADFDE5A,
    0x2F0E1EBA9EA36930, 0x6DFEFF5137495FA3,
    0xAAEFDD6DCD770416, 0xE81F3C86649D3285,
    0xF45BB4758C645C51, 0xB6AB559E258E6AC2,
    0x71BA77A2DFB03177, 0x334A9649765A07E4,
    0xBD68D2308226B08E, 0xFF9833DB2BCC861D,
    0x388911E7D1F2DDA8, 0x7A79F00C7818EB3B,
    0xCC7AF1FF21C30BDE, 0x8E8A101488293D4D,
    0x499B3228721766F8, 0x0B6BD3C3DBFD506B,
    0x854997BA2F81E701, 0xC7B97651866BD192,
    0x00A8546D7C558A27, 0x4258B586D5BFBCB4,
    0x5E1C3D753D46D260, 0x1CECDC9E94ACE4F3,
    0xDBFDFEA26E92BF46, 0x990D1F49C77889D5,
    0x172F5B3033043EBF, 0x55DFBADB9AEE082C,
    0x92CE98E760D05399, 0xD03E790CC93A650A,
    0xAA478900B1228E31, 0xE8B768EB18C8B8A2,
    0x2FA64AD7E2F6E317, 0x6D56AB3C4B1CD584,
    0xE374EF45BF6062EE, 0xA1840EAE168A547D,
    0x66952C92ECB40FC8, 0x2465CD79455E395B,
    0x3821458AADA7578F, 0x7AD1A461044D611C,
    0xBDC0865DFE733AA9, 0xFF3067B657990C3A,
    0x711223CFA3E5BB50, 0x33E2C2240A0F8DC3,
    0xF4F3E018F031D676, 0xB60301F359DBE0E5,
    0xDA050215EA6C212F, 0x98F5E3FE438617BC,
    0x5FE4C1C2B9B84C09, 0x1D14202910527A9A,
    0x93366450E42ECDF0, 0xD1C685BB4DC4FB63,
    0x16D7A787B7FAA0D6, 0x5427466C1E109645,
    0x4863CE9FF6E9F891, 0x0A932F745F03CE02,
    0xCD820D48A53D95B7, 0x8F72ECA30CD7A324,
    0x0150A8DAF8AB144E, 0x43A04931514122DD,
    0x84B16B0DAB7F7968, 0xC6418AE602954FFB,
    0xBC387AEA7A8DA4C0, 0xFEC89B01D3679253,
    0x39D9B93D2959C9E6, 0x7B2958D680B3FF75,
    0xF50B1CAF74CF481F, 0xB7FBFD44DD257E8C,
    0x70EADF78271B2539, 0x321A3E938EF113AA,
    0x2E5EB66066087D7E, 0x6CAE578BCFE24BED,
    0xABBF75B735DC1058, 0xE94F945C9C3626CB,
    0x676DD025684A91A1, 0x259D31CEC1A0A732,
    0xE28C13F23B9EFC87, 0xA07CF2199274CA14,
    0x167FF3EACBAF2AF1, 0x548F120162451C62,
    0x939E303D987B47D7, 0xD16ED1D631917144,
    0x5F4C95AFC5EDC62E, 0x1DBC74446C07F0BD,
    0xDAAD56789639AB08, 0x985DB7933FD39D9B,
    0x84193F60D72AF34F, 0xC6E9DE8B7EC0C5DC,
    0x01F8FCB784FE9E69, 0x43081D5C2D14A8FA,
    0xCD2A5925D9681F90, 0x8FDAB8CE70822903,
    0x48CB9AF28ABC72B6, 0x0A3B7B1923564425,
    0x70428B155B4EAF1E, 0x32B26AFEF2A4998D,
    0xF5A348C2089AC238, 0xB753A929A170F4AB,
    0x3971ED50550C43C1, 0x7B810CBBFCE67552,
    0xBC902E8706D82EE7, 0xFE60CF6CAF321874,
    0xE224479F47CB76A0, 0xA0D4A674EE214033,
    0x67C58448141F1B86, 0x253565A3BDF52D15,
    0xAB1721DA49899A7F, 0xE9E7C031E063ACEC,
    0x2EF6E20D1A5DF759, 0x6C0603E6B3B7C1CA,
    0xF6FAE5C07D3274CD, 0xB40A042BD4D8425E,
    0x731B26172EE619EB, 0x31EBC7FC870C2F78,
    0xBFC9838573709812, 0xFD39626EDA9AAE81,
    0x3A28405220A4F534, 0x78D8A1B9894EC3A7,
    0x649C294A61B7AD73, 0x266CC8A1C85D9BE0,
    0xE17DEA9D3263C055, 0xA38D0B769B89F6C6,
    0x2DAF4F0F6FF541AC, 0x6F5FAEE4C61F773F,
    0xA84E8CD83C212C8A, 0xEABE6D3395CB1A19,
    0x90C79D3FEDD3F122, 0xD2377CD44439C7B1,
    0x15265EE8BE079C04, 0x57D6BF0317EDAA97,
    0xD9F4FB7AE3911DFD, 0x9B041A914A7B2B6E,
    0x5C1538ADB04570DB, 0x1EE5D94619AF4648,
    0x02A151B5F156289C, 0x4051B05E58BC1E0F,
    0x87409262A28245BA, 0xC5B073890B687329,
    0x4B9237F0FF14C443, 0x0962D61B56FEF2D0,
    0xCE73F427ACC0A965, 0x8C8315CC052A9FF6,
    0x3A80143F5CF17F13, 0x7870F5D4F51B4980,
    0xBF61D7E80F251235, 0xFD913603A6CF24A6,
    0x73B3727A52B393CC, 0x31439391FB59A55F,
    0xF652B1AD0167FEEA, 0xB4A25046A88DC879,
    0xA8E6D8B54074A6AD, 0xEA16395EE99E903E,
    0x2D071B6213A0CB8B, 0x6FF7FA89BA4AFD18,
    0xE1D5BEF04E364A72, 0xA3255F1BE7DC7CE1,
    0x64347D271DE22754, 0x26C49CCCB40811C7,
    0x5CBD6CC0CC10FAFC, 0x1E4D8D2B65FACC6F,
    0xD95CAF179FC497DA, 0x9BAC4EFC362EA149,
    0x158E0A85C2521623, 0x577EEB6E6BB820B0,
    0x906FC95291867B05, 0xD29F28B9386C4D96,
    0xCEDBA04AD0952342, 0x8C2B41A1797F15D1,
    0x4B3A639D83414E64, 0x09CA82762AAB78F7,
    0x87E8C60FDED7CF9D, 0xC51827E4773DF90E,
    0x020905D88D03A2BB, 0x40F9E43324E99428,
    0x2CFFE7D5975E55E2, 0x6E0F063E3EB46371,
    0xA91E2402C48A38C4, 0xEBEEC5E96D600E57,
    0x65CC8190991CB93D, 0x273C607B30F68FAE,
    0xE02D4247CAC8D41B, 0xA2DDA3AC6322E288, 
    0xBE992B5F8BDB8C5C, 0xFC69CAB42231BACF, 
    0x3B78E888D80FE17A, 0x7988096371E5D7E9, 
    0xF7AA4D1A85996083, 0xB55AACF12C735610, 
    0x724B8ECDD64D0DA5, 0x30BB6F267FA73B36, 
    0x4AC29F2A07BFD00D, 0x08327EC1AE55E69E, 
    0xCF235CFD546BBD2B, 0x8DD3BD16FD818BB8, 
    0x03F1F96F09FD3CD2, 0x41011884A0170A41, 
    0x86103AB85A2951F4, 0xC4E0DB53F3C36767, 
    0xD8A453A01B3A09B3, 0x9A54B24BB2D03F20, 
    0x5D45907748EE6495, 0x1FB5719CE1045206, 
    0x919735E51578E56C, 0xD367D40EBC92D3FF, 
    0x1476F63246AC884A, 0x568617D9EF46BED9, 
    0xE085162AB69D5E3C, 0xA275F7C11F7768AF, 
    0x6564D5FDE549331A, 0x279434164CA30589, 
    0xA9B6706FB8DFB2E3, 0xEB46918411358470, 
    0x2C57B3B8EB0BDFC5, 0x6EA7525342E1E956, 
    0x72E3DAA0AA188782, 0x30133B4B03F2B111, 
    0xF7021977F9CCEAA4, 0xB5F2F89C5026DC37, 
    0x3BD0BCE5A45A6B5D, 0x79205D0E0DB05DCE, 
    0xBE317F32F78E067B, 0xFCC19ED95E6430E8, 
    0x86B86ED5267CDBD3, 0xC4488F3E8F96ED40, 
    0x0359AD0275A8B6F5, 0x41A94CE9DC428066, 
    0xCF8B0890283E370C, 0x8D7BE97B81D4019F, 
    0x4A6ACB477BEA5A2A, 0x089A2AACD2006CB9, 
    0x14DEA25F3AF9026D, 0x562E43B4931334FE, 
    0x913F6188692D6F4B, 0xD3CF8063C0C759D8, 
    0x5DEDC41A34BBEEB2, 0x1F1D25F19D51D821, 
    0xD80C07CD676F8394, 0x9AFCE626CE85B507 
] 

class CRC64(object):

    def __init__(self):
        self.crc = 0xffffffffffffffff

    def append(self, buffer):
        for c in buffer:
            tab_index = ((self.crc >> 56) ^ ord(c)) & 0xFF
            self.crc = crc_table[tab_index] ^ ((self.crc << 8) &
                                               0xffffffffffffffff)

    def fini(self):
        tmp=self.crc ^ 0 #L #doesn't work in python3
        return tmp


def crc64(buffer):
    crc = CRC64()
    crc.append(buffer)
    
    return crc.fini()
