#Partition an alignment in Phylip format to disjoint subsets (a partition scheme)
#Author: Thulek@gmail.com
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time
import array as arr
#import runCommand

import argparse

text = "=======mPartitions======="

parser = argparse.ArgumentParser(description = text)  
#parser.parse_args()  

jobIds = arr.array('i')

parser.add_argument("-f", "--filex", help="The alignment file")
parser.add_argument("-o", "--output", help="Output directory")
parser.add_argument("-t", "--tper", help="Minimum length")
parser.add_argument("-mset", "--mset", help="mset model")

# read arguments from the command line
args = parser.parse_args()

filex = args.filex
output = "Output"
if args.output:
	output = args.output

####Check DNA or AA###
prot = 0	
cdna = open(filex,"r")
zz = 0
for ll in cdna:
	if zz == 1:
		ll = ll.upper()
		if " " in ll:
			if len(ll.split(" ",1)[1].strip()) == ll.split(" ",1)[1].strip().count("-") or len(ll.split(" ",1)[1].strip()) == ll.split(" ",1)[1].strip().count("N"):
				zz = zz-1
			else:
				seq = ll.split(" ",1)[1].strip().replace("A","").replace("T","").replace("G","").replace("C","").replace("X","").replace("-","").replace("N","")
				if len(seq) > 0:
					prot = 0
				else:
					prot = 1
		elif "\t" in ll:
			if len(ll.split("\t",1)[1].strip()) == ll.split("\t",1)[1].strip().count("-") or len(ll.split("\t",1)[1].strip()) == ll.split("\t",1)[1].strip().count("N"):
				zz = zz-1
			else:
				seq = ll.split("\t",1)[1].strip().replace("A","").replace("T","").replace("G","").replace("C","").replace("X","").replace("-","").replace("N","")
				#print(seq)
				if len(seq) > 0:
					prot = 0
				else:
					prot = 1
	zz += 1
cdna.close()

mset = "JC69,F81,HKY,GTR"
#if prot == 0:
#	mset = "LG,WAG,JTT"
#else:
#	mset = "JC69,F81,HKY,GTR"
if args.mset:
	mset = args.mset

tper = 50.0
if args.tper:
	tper = float(args.tper)


f = open(filex,"r")
line = f.readline()
f.close()
if (" " in line):
	maxlength = int(line.split(" ",1)[1].strip())
else:
	maxlength = int(line.split("\t",1)[1].strip())
if tper >= 50:
	tper = float(tper*100/maxlength)


if not os.path.isdir(output):
	cmd = 'mkdir '+output
	os.system(cmd)
else:
	print(output+" is exists")
	
if not os.path.isdir("logs"):
	cmd = 'mkdir logs'
	os.system(cmd)
else:
	print("logs is exists")
	

###CREATE Results Folder####
if not os.path.isdir(output+"/Results"):
	cmd = 'mkdir '+output+"/Results"
	os.system(cmd)
else:
	print("Don't create "+output+"/Results.")

if not os.path.isdir(output+"/siteLH"):
	cmd = 'mkdir '+output+"/siteLH"
	os.system(cmd)
else:
	print("Don't create "+output+"/siteLH.")

def getFilename(filex):
	treefn = filex
	if "/" in filex:
		treef = filex.split("/")
		treefn = treef[len(treef)-1].strip()
	else:
		treefn = filex
	return treefn
	

alignList = []
alignList.append(filex)
num_run = 1
treefn = ""
if "/" in filex:
	treef = filex.split("/")
	treefn = treef[len(treef)-1].strip()
else:
	treefn = filex
treefile = output+"/"+treefn+".tree"
print("mset: "+mset)
while len(alignList) > 0:
	for align in alignList:
		if num_run == 1 and len(alignList) == 1:
			command = "python mPartition_3part.py -f "+align+" -m "+str(maxlength)+" -tper "+str(tper)+" -mset "+mset+" -o "+output
			os.system(command)

			extS = 0
			if(os.path.isfile(output+"/"+getFilename(align)+"P1")):
				alignList.append(output+"/"+getFilename(align)+"P1")
				extS = 1
			if(os.path.isfile(output+"/"+getFilename(align)+"P2")):
				alignList.append(output+"/"+getFilename(align)+"P2")
				extS = 1
			if(os.path.isfile(output+"/"+getFilename(align)+"P3")):
				alignList.append(output+"/"+getFilename(align)+"P3")
				extS = 1
			alignList.remove(align)
			if extS == 0:
				os.system("cp "+align+" "+output+"/Results/")
		else:
			parfile = output + "/par."+treefn
			tiger_file = "rate_"+treefn
			inv_file = "inv_"+treefn
			command = "python mPartition_3part.py -f "+align+" -m "+str(maxlength)+" -tiger "+tiger_file+" -inv"+inv_file+" -tper "+str(tper)+" -mset "+mset+" -t "+treefile+" -p "+parfile+" -o "+output

			os.system(command)
			extS = 0
			if(os.path.isfile(output+"/"+getFilename(align)+"P1")):
				alignList.append(output+"/"+getFilename(align)+"P1")
				extS = 1
			if(os.path.isfile(output+"/"+getFilename(align)+"P2")):
				alignList.append(output+"/"+getFilename(align)+"P2")
				extS = 1
			if(os.path.isfile(output+"/"+getFilename(align)+"P3")):
				alignList.append(output+"/"+getFilename(align)+"P3")
				extS = 1
			alignList.remove(align)
			if extS == 0:
				os.system("cp "+align+" "+output+"/Results/")
	print("Number of Executes: "+str(num_run)) 
	num_run += 1
print("Mission Completed.")

if(os.path.isfile(output + "/par."+treefn)):
	par = open(output + "/par."+treefn,"r")
	for line in par:
		v = line.split(";")
		vfile = open(output+"/par."+treefn+"_parf_"+v[1].replace("-","_").strip(),"a+")
		vfile.write(v[0].strip()+" ")
		vfile.close()
	par.close()
	if(os.path.isfile(output + "/Results/par."+treefn)):
		os.system("rm "+output + "/Results/par."+treefn)
	finishParfile = open(output + "/Results/par."+treefn,"a+")
	finishParfile.write("#nexus\nbegin sets;\n")
	for file in os.listdir(output):
		if "par."+treefn+"_parf_" in file:
			finishParfile.write("\tcharset "+file.replace("par."+treefn+"_parf_","")+" = ")
			wp = open(output+"/"+file,"r")
			for line in wp:
				finishParfile.write(line.strip()+";\n")
	finishParfile.write("end;\n")
	finishParfile.close()	
	os.system("rm "+output + "/par."+treefn+"_parf_*")
else:
	if(os.path.isfile(output + "/Results/"+treefn)):
		if(os.path.isfile(output + "/Results/par."+treefn)):
			os.system("rm "+output + "/Results/par."+treefn)
		finishParfile = open(output + "/Results/par."+treefn,"a+")
		finishParfile.write("#nexus\nbegin sets;\n")
		finishParfile.write("\tcharset Par1 = 1-"+str(maxlength)+";\n")
		finishParfile.write("end;\n")
		finishParfile.close()	

os.system("rm *"+treefn+"*.bsub")
os.system("rm "+output+"/"+treefn+"*")

if(not os.path.isdir("Results")):
	os.system("mkdir Results")
os.system("cp "+output+"/Results/*.* Results/")
os.system("rm -rf "+output+"/Results")
