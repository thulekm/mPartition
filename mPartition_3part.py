#Partition an alignment in Phylip format to 3 subsets 
#Author: Thulek@gmail.com
#using these tools: IQ-TREE, TIGER
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time
import array as arr
import config
#import runCommand

import argparse

text = "===================\nSplit partitions from an alignments\n===================="

parser = argparse.ArgumentParser(description = text)  
#parser.parse_args()  

parser.add_argument("-f", "--filex", help="The alignment file")
parser.add_argument("-t", "--treefile", help="The tree file")
parser.add_argument("-tiger", "--tiger", help="The tiger file")
parser.add_argument("-m", "--maxlength", help="Max length of sequences")
parser.add_argument("-p", "--parfile", help="Partition file")
parser.add_argument("-tper", "--tper", help="Percent - minimum length")
parser.add_argument("-mset", "--mset", help="mset model, separated by commas")
parser.add_argument("-o", "--output", help="The output directory")
# read arguments from the command line
args = parser.parse_args()

filex = args.filex
output = args.output
if args.treefile:
	treefile = args.treefile
else:
	treefile = "NOTTREE"
	

tig_rate = "0"
if args.tiger:
	tig_rate = args.tiger

#iqtree_path = "$IQTREE/" #IQ-TREE path including splash
#tiger_path = "/home/lkthu/Partitions/tiger_original/"

iqtree_path = config.iqtree_path
tiger_path = config.tiger_path

mset = "LG,WAG"
if args.mset:
	mset = args.mset

tper = 10.0
if args.tper:
	tper = float(args.tper)

if tper >= 50: 
	sys.exit("The minimum percent parameter must be less than 50.")

if args.maxlength:
	maxlength = int(args.maxlength)
else:
	f = open(filex)
	line = f.readline()
	f.close()
	maxlength = int(line.split(" ",1)[1].strip())

if not os.path.isdir(output):
	cmd = 'mkdir '+output
	os.system(cmd)
else:
	print(output+" is exists")
	
	
def countLines(files):
	if(os.path.isfile(files)):
		with open(files) as f:
			return len(f.readlines())
	else:
		return 0

def convertPhylip2Fasta(filex1,tofile):
	infile = open(filex1,"r")
	outfile = open(tofile,"w")
	i = 1
	for line in infile:
		if i > 1:
			if len(line.strip()) >= 1:
				if "\t" in line:
					tx = line.split("\t",1)
					outfile.write(">"+tx[0].strip()+"\n")
					outfile.write(tx[1].strip()+"\n")
				elif " " in line:
					tx = line.split(" ",1)
					outfile.write(">"+tx[0].strip()+"\n")
					outfile.write(tx[1].strip()+"\n")
		i = i + 1
	outfile.close()
	infile.close()
	if(countLines(tofile)>0):
		print "Finish to create the fasta file."
	else:
		print "Create fasta file error."

def line_prepender(filename, line):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)

def removeGapsPhylip(filex):
	cfile = open(filex,"r")
	i = 0
	t = 0
	treefn = ""
	if "/" in filex:
		treef = filex.split("/")
		treefn = "ck"+treef[len(treef)-1].strip()
	else:
		treefn = "ck"+filex
	if os.path.isfile(output+"/"+treefn):
		os.system("rm "+output+"/"+treefn)
		
	rfile = open(output+"/"+treefn,"w")
	for line in cfile:
		if i == 0:
			if " " in line:
				tax = int(line.strip().split(" ",1)[0].strip())
				sitex = int(line.strip().split(" ",1)[1].strip())
			elif "\t" in line:
				tax = int(line.strip().split("\t",1)[0].strip())
				sitex = int(line.strip().split("\t",1)[1].strip())
		elif i >= 1:
			if " " in line:
				checkLine = line.split(" ",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
					
			elif "\t" in line:
				checkLine = line.split("\t",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
		i += 1
	cfile.close()
	rfile.close()
	if t == 0:
		print("This file is ok.")
		if os.path.isfile(output+"/"+treefn):
			os.system("rm "+output+"/"+treefn)
		#newfile = filex
	else:
		print("This file has "+str(t)+" seqs only gaps or missing data.")
		tax = tax - t
		line_prepender(output+"/"+treefn,str(tax)+" "+str(sitex))
		

min_rate = 1.0
max_rate = 0.0
#tg_rate = []

treefn = ""

if "/" in filex:
	treef = filex.split("/")
	treefn = treef[len(treef)-1].strip()
else:
	treefn = filex
path = filex.strip(treefn)

#print treefn

parfile = ""
if args.parfile:
	parfile = args.parfile
else:
	parfile = output+"/par."+treefn+""
	if(os.path.isfile(parfile)):
		os.system("rm "+parfile)

removeGapsPhylip(filex) #remove gaps from phylip
newfile = filex
if os.path.isfile(output+"/ck"+treefn):
	newfile = output+"/ck"+treefn
print(newfile)
checkstate = 0
rvalue = []
min_rate = 0.0
max_rate = 0.0
if tig_rate == "0":
	convertPhylip2Fasta(newfile,output+"/"+treefn+".FASTA")
	if not os.path.isfile(output+"/rate_"+treefn):
		os.system(tiger_path+"tiger -in "+output+"/"+treefn+".FASTA -f s,r -rl "+output+"/rate_"+treefn) 

	os.system("rm "+output+"/"+treefn+".FASTA")
	
	if os.path.isfile(output+"/rate_"+treefn):
		ratefile = open(output+"/rate_"+treefn,"r")
		for line in ratefile:
			rvalue.append(float(line.strip()))
		ratefile.close()
		checkstate = 1
	else:
		print "Rate file's not exists, please check the input file. \n"
else:
	#if tiger_rate file's not exists.
	if not os.path.isfile(output+"/"+tig_rate) or not os.path.isfile(parfile):
		convertPhylip2Fasta(newfile,output+"/"+treefn+".FASTA")
		os.system(tiger_path+"tiger -in "+output+"/"+treefn+".FASTA -f s,r -rl "+output+"/rate_"+treefn) 
		os.system("rm "+output+"/"+treefn+".FASTA")
		if os.path.isfile(output+"/rate_"+treefn):
			ratefile = open(output+"/rate_"+treefn,"r")
			for line in ratefile:
				rvalue.append(float(line.strip()))
			ratefile.close()
		else:
			print "Rate file's not exists, please check the input file. \n"
	else:
		ftiger = open(output+"/"+tig_rate,"r")
		lines=ftiger.readlines()
		ftiger.close()	
		opar = open(parfile,"r")
		i = 0
		for ln in opar:
			if ln.split(";")[1].strip() in treefn:
				rvalue.append(float(lines[i].strip()))
			i+=1	
		opar.close()
	checkstate = 1
	
par1 = []
par2 = []
par3 = []
if checkstate == 1:
	
	newList = sorted(rvalue)
	if int(round(len(rvalue)*1/100))-1 >= 50:
		min_rate = float(newList[int(round(len(rvalue)*1/100))-1])
		max_rate = float(newList[int(round(len(rvalue)*99/100))-1])
	else:
		min_rate = float(newList[0])
		max_rate = float(newList[int(round(len(rvalue)))-1])
	avg_rate = min_rate+(max_rate-min_rate)/3
	a2vg_rate = min_rate+(max_rate-min_rate)*2/3 
	
	print "avg rate: "+str(avg_rate)
	print("2 check max rate: "+str(max_rate) + "; min rate: "+ str(min_rate)+"; 1/3 rate: "+str(avg_rate)+" 2/3 rate: "+str(a2vg_rate)+"\n")
	par1 = []
	par2 = []
	par3 = []
	i=1
	#ratefile = open(output+"/rate_"+treefn,"r")
	for line in rvalue:
		if float(line) > a2vg_rate:
			par3.append(i)
			#print("par2 added seq: "+str(i))
		elif float(line) < avg_rate:
			par1.append(i)
			#print("par1 added seq: "+str(i))
		else:
			par2.append(i)
		i = i + 1
	#ratefile.close()
	if (len(par1)==0 or len(par2)==0 or len(par3)==0):
		par1 = []
		par2 = []
		par3 = []
		i=1
		#ratefile = open(output+"/rate_"+treefn,"r")
		for line in rvalue:
			if float(line) >= a2vg_rate and len(par3) <= len(rvalue)/3:
				par3.append(i)
				#print("par2 added seq: "+str(i))
			elif float(line) <= avg_rate and len(par1) <= len(rvalue)/3:
				par1.append(i)
				#print("par1 added seq: "+str(i))
			else:
				par2.append(i)
			i = i + 1
		ratefile.close()
	
	parFile = open(output+"/S1_Par_"+treefn,"w")
	parFile.write("#nexus\nbegin sets;\n")
	line = "charset Par1 ="
	for p in par1:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\n")
	line = "charset Par2 ="
	for p in par2:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\n")
	line = "charset Par3 ="
	for p in par3:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\nend;")
	parFile.close()
	command = iqtree_path+"iqtree -s "+newfile+" -m MFP -fast -mset "+mset+" -spp "+output+"/S1_Par_"+treefn+" -pre "+output+"/"+treefn+"\n"
	os.system(command)
	if not os.path.isfile(output+"/"+treefn+"_AICc.iqtree"):
		command = iqtree_path+"iqtree -s "+newfile+" -m MFP -fast -mset "+mset+" -pre "+output+"/"+treefn+"_AICc\n"
		os.system(command)
	
model =""
first_model = ""
second_model = ""
third_model = ""
if(os.path.isfile(output+"/"+treefn+".iqtree")):
	gtModel = open(output+"/"+treefn+".iqtree","r")
	for l in gtModel:
		if "Best-fit model according to BIC: " in l:
			model = l.split("Best-fit model according to BIC: ")[1].strip()
			fmodel = model.split(":")
			first_model = fmodel[0].strip()
			second_model = fmodel[1].split(",")[1].strip()
			third_model = fmodel[2].split(",")[1].strip()
	gtModel.close()
first_model = first_model.replace("+ASC","")
second_model = second_model.replace("+ASC","")
third_model = third_model.replace("+ASC","")
print("First Best Model: "+first_model+" | Treefile: "+output+"/"+treefn+".treefile\n")
print("Second Best Model: "+second_model+" | Treefile: "+output+"/"+treefn+".treefile\n")
print("Third Best Model: "+third_model+" | Treefile: "+output+"/"+treefn+".treefile\n")

if treefile == "NOTTREE":
	os.system("cp "+output+"/"+treefn+"_AICc.treefile "+output+"/"+treefn+".tree")
	treefile = output+"/"+treefn+".tree"


def getAIC(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Akaike information criterion (AIC) score" in line.strip():
				text += line.split(" ")[5]
				return text

#get n
def getN(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Input data" in line.strip():
				text += line.split(" ")[2]
				print(text)
				return text

#get free para
def getK(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "free parameters" in line.strip():
				text += line.split(" ")[8]
				return text

command = iqtree_path+"iqtree -s "+newfile+" -m "+first_model+" -te "+treefile+" -wsl -pre "+output+"/"+treefn+"_G1\n"
os.system(command)
command = iqtree_path+"iqtree -s "+newfile+" -m "+second_model+" -te "+treefile+" -wsl -pre "+output+"/"+treefn+"_G2\n"
os.system(command)
command = iqtree_path+"iqtree -s "+newfile+" -m "+third_model+" -te "+treefile+" -wsl -pre "+output+"/"+treefn+"_G3\n"
os.system(command)

g1 = []
g2 = []
g3 = []
if(os.path.isfile(output+"/"+treefn+"_G1.sitelh") and os.path.isfile(output+"/"+treefn+"_G2.sitelh") and os.path.isfile(output+"/"+treefn+"_G3.sitelh")):
	getG1 = open(output+"/"+treefn+"_G1.sitelh","r")
	i = 1
	for l in getG1:
		if i > 1:
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g1.append(float(x))
					#print(x+"\n")
				z+=1
		i = i + 1
	getG1.close()
	getG2 = open(output+"/"+treefn+"_G2.sitelh","r")
	i = 1
	for l in getG2:
		if i > 1:
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g2.append(float(x))
				z+=1
		i = i + 1
	getG2.close()
	getG3 = open(output+"/"+treefn+"_G3.sitelh","r")
	i = 1
	for l in getG3:
		if i > 1:
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g3.append(float(x))
				z+=1
		i = i + 1
	getG3.close()

	par1 = []
	par2 = []
	par3 = []
	i = 0
	while i < len(g1):
		if g1[i] > g2[i] and g1[i]>g3[i]:
			par1.append(i+1)
		elif g2[i]>g3[i]:
			par2.append(i+1)
		else:
			par3.append(i+1)
		i += 1

	print("par1: "+str(len(par1))+" | par2: "+str(len(par2))+" | par3: "+str(len(par3)))
	
	if (len(par1) < int((maxlength)*tper/100) and len(par1)>0) or (len(par2)>0 and len(par2) < int((maxlength)*tper/100)) or (len(par3)>0 and len(par3) < int((maxlength)*tper/100)):
		print "Not good. Partition's length is very small."
	elif (len(par1) + len(par2) ==0) or (len(par1) + len(par3) ==0) or (len(par3) + len(par2) ==0):
		print "Don't need split."
	else:
		parFile = open(output+"/F1_Par_"+treefn,"w")
		if len(par1)>0:
			line = "P1 ="
			for p in par1:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		if len(par2)>0:
			line = "P2 ="
			for p in par2:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		if len(par3)>0:
			line = "P3 ="
			for p in par3:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		parFile.close()
		#os.system("python extractPartitions.py -f "+filex+" -p "+output+"/F1_Par_"+treefn+" -o "+output)
		os.system("python splitPartition.py -f "+newfile+" -p "+output+"/F1_Par_"+treefn+"")
		
		if(os.path.isfile(path+""+treefn+"P1") and not os.path.isfile(output+"/"+treefn+"P1")):
			os.system("cp "+path+""+treefn+"P1 "+output+"/"+treefn+"P1")
			os.system("rm "+path+""+treefn+"P1")
		if(os.path.isfile(path+""+treefn+"P2") and not os.path.isfile(output+"/"+treefn+"P2")):
			os.system("cp "+path+""+treefn+"P2 "+output+"/"+treefn+"P2")
			os.system("rm "+path+""+treefn+"P2")
		if(os.path.isfile(path+""+treefn+"P3") and not os.path.isfile(output+"/"+treefn+"P3")):
			os.system("cp "+path+""+treefn+"P3 "+output+"/"+treefn+"P3")
			os.system("rm "+path+""+treefn+"P3")
		
		if(os.path.isfile(output+"/ck"+treefn+"P1") and not os.path.isfile(output+"/"+treefn+"P1")):
			os.system("cp "+output+"/ck"+treefn+"P1 "+output+"/"+treefn+"P1")
			os.system("rm "+output+"/ck"+treefn+"P1")
		if(os.path.isfile(output+"/ck"+treefn+"P2") and not os.path.isfile(output+"/"+treefn+"P2")):
			os.system("cp "+output+"/ck"+treefn+"P2 "+output+"/"+treefn+"P2")
			os.system("rm "+output+"/ck"+treefn+"P2")
		if(os.path.isfile(output+"/ck"+treefn+"P3") and not os.path.isfile(output+"/"+treefn+"P3")):
			os.system("cp "+output+"/ck"+treefn+"P3 "+output+"/"+treefn+"P3")
			os.system("rm "+output+"/ck"+treefn+"P3")
		
		
		#os.system("mv "+output+"/"+treefn+"_Par1 "+output+"/"+treefn+"P1")
		#os.system("mv "+output+"/"+treefn+"_Par2 "+output+"/"+treefn+"P2")
		command = ""
		if(os.path.isfile(output+"/"+treefn+"P1")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P1 -m MFP -fast -mset "+mset+" -t "+treefile+"  -pre "+output+"/"+treefn+"P1_AICc\n"
			os.system(command)
		if(os.path.isfile(output+"/"+treefn+"P2")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P2 -m MFP -fast -mset "+mset+" -t "+treefile+"  -pre "+output+"/"+treefn+"P2_AICc\n"
			os.system(command)
		if(os.path.isfile(output+"/"+treefn+"P3")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P3 -m MFP -fast -mset "+mset+" -t "+treefile+"  -pre "+output+"/"+treefn+"P3_AICc\n"
			os.system(command)
		
		aic = 0.0
		k = 0.0
		n = 0.0
		files = [f for f in listdir(output)]
		for f in files:
			if (treefn+"P1_AICc.iqtree" in f) or (treefn+"P2_AICc.iqtree" in f)  or (treefn+"P3_AICc.iqtree" in f):
				n += float(getN(output+"/"+f))
				k += float(getK(output+"/"+f))
				aic += float(getAIC(output+"/"+f))
		AICc = aic + 2*k*(k+1)/(n-k-1)
		FirstAICc = float(getAIC(output+"/"+treefn+"_AICc.iqtree"))
		print "FirstAICc: "+str(FirstAICc)+" | New AICc: "+str(AICc)
		if AICc > FirstAICc:
			print "Worser."
			os.system("rm "+output+"/"+treefn+"P*")
			#os.system("rm "+output+"/"+treefn+"P2*")
		else:
			if(os.path.isfile(output+"/"+treefn+"P3") and os.path.isfile(output+"/"+treefn+"P2") and os.path.isfile(output+"/"+treefn+"P1")):
				print "Better. Extract "+ filex+ " to 3 partitions: "+output+"/"+treefn+"P[1,2,3]."
			else:
				print "Better. Extract "+ filex+ " to 2 partitions: "+output+"/"+treefn+"P[1,2]."
			if(os.path.isfile(parfile)):
				opar = open(parfile,"r")
				ipar = open(parfile+"_1","w")
				i = 1
				z = 1
				for line in opar:
					if line.split(";")[1].strip() in treefn:
						if g1[z-1] > g2[z-1] and g1[z-1]>g3[z-1]:
						#if g1[z-1] > g2[z-1]:
							ipar.write(str(i)+";"+treefn+"P1\n")
						#else:
						elif g2[z-1]>g3[z-1]:
							ipar.write(str(i)+";"+treefn+"P2\n")
						else:
							ipar.write(str(i)+";"+treefn+"P3\n")
						z = z+1
					else:
						ipar.write(line.strip()+"\n")
					i += 1
				ipar.close()
				opar.close()
				os.system("mv "+parfile+"_1 "+parfile)
			else:
				opar = open(parfile,"w")
				i = 0
				while i < len(g1):
					if g1[i] > g2[i] and g1[i]>g3[i]:
						opar.write(str(i+1)+";"+treefn+"P1\n")
					elif g2[i]>g3[i]:
						opar.write(str(i+1)+";"+treefn+"P2\n")
					else:
						opar.write(str(i+1)+";"+treefn+"P3\n")
					i += 1
				opar.close()
			#os.system("rm "+output+"/"+treefn+"P*AICc*")
			#os.system("rm "+output+"/"+treefn+"P2*")
		
