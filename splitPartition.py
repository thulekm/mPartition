#Author: Thulek@gmail.com
#Split an alignment into sub-alignments using a predefined partitioning scheme
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time

import array as arr
import random

import argparse
parser = argparse.ArgumentParser() 

parser.add_argument("-f", "--filex", help="File")
parser.add_argument("-p", "--part", help="par file")
args = parser.parse_args()

parFile = args.part
filex = args.filex

def line_prepender(filename, line):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)


parF = open(parFile,"r")
parNum = 1
for line in parF:
	if "=" in line:
		siteArr = line.split("=")[1].strip().replace(";","").replace(",","")
		siteList = siteArr.split(" ")
		if os.path.isfile(filex+"P"+str(z)):
			os.system("rm "+filex+"P"+str(z))
		par = open(filex+"P"+str(parNum),"w")
		fil = open(filex,"r")
		noOfTax = 0
		for lx in fil:
			if noOfTax > 0:
				if " " in lx:
					seq1 = lx.split(" ",1)
					par.write(str(seq1[0].strip())+" ")
					seqContent = seq1[1].strip()
					for x in siteList:
						par.write(seqContent[int(x.strip())-1])
					par.write("\n")
				elif "\t" in lx:
					seq1 = lx.split("\t",1)
					par.write(str(seq1[0].strip())+" ")
					seqContent = seq1[1].strip()
					for x in siteList:
						par.write(seqContent[int(x.strip())-1])
					par.write("\n")
			noOfTax += 1
		fil.close()
		par.close()
		line_prepender(filex+"P"+str(z), str(noOfTax-1)+" "+str(len(siteList)))
		parNum = parNum+1
parF.close()
