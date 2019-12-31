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

#def count_str(filename, char):
#	with open(filename, 'r+') as f:
#		content = f.read()
#		return content.count("=")

wp = open(parFile,"r")
z = 1
for line in wp:
	if "=" in line:
		#par1 = []
		dt = line.split("=")[1].strip().replace(";","").replace(",","")
		dt1 = dt.split(" ")
		if os.path.isfile(filex+"P"+str(z)):
			os.system("rm "+filex+"P"+str(z))
		par = open(filex+"P"+str(z),"w")
		fil = open(filex,"r")
		seq = 0
		i = 0
		for lx in fil:
			if i > 0:
				if " " in lx:
					seq1 = lx.split(" ",1)
					par.write(str(seq1[0].strip())+" ")
					for x in dt1:
						par.write(seq1[1].strip()[int(x.strip())-1])
					par.write("\n")
				elif "\t" in lx:
					seq1 = lx.split("\t",1)
					par.write(str(seq1[0].strip())+" ")
					for x in dt1:
						par.write(seq1[1].strip()[int(x.strip())-1])
					par.write("\n")
			i += 1
		fil.close()
		par.close()
		line_prepender(filex+"P"+str(z), str(i-1)+" "+str(len(dt1)))
		z = z+1
wp.close()