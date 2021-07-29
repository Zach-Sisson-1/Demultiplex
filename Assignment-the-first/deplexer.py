#!/usr/bin/env python

import gzip
import argparse
import Bioinfo as bf
import numpy as np
import matplotlib.pyplot as plt

def get_args():
	parser = argparse.ArgumentParser(description="A script to generate a distribution plot of mean Qscore per nucleotide")
	parser.add_argument("-f", "--file_name", help = "Path to file to run script on", type=str, required=True )
	parser.add_argument("-r", "--read_length", help = "Length of read (bp)", type=int, required=True )
	parser.add_argument("-o", "--output", help = "Title of graph to output", type=str, required=False, default="Distribution of mean qual scores per nucleotide" )
	return parser.parse_args()

args = get_args()

#Initialize list of zeros:
def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 'x' values of 0.0. where x = read_length passed in'''
    for i in range(args.read_length):
        lst.append(value)
    return lst
  
Qscore_list: list = []
Qscore_list = init_list(Qscore_list)


#Loops through file, keeping a running sum of the quality score at each base pair position, then calculating the mean after the final line. 
with gzip.open(args.file_name, "rt") as fh:
	index = 1
	line_count = 0
	for line in fh:
		line = line.strip('\n')
		if index % 4 == 0:
			for i in range(args.read_length):
				Qscore_list[i] += (bf.convert_phred(line[i]))
			line_count += 1
		index +=1
	#Updates for long runs
	if index % 4000000 ==0:
		print("Reading line "+str(index))


#Calculates mean Qscore at each base-pair position and generates distribution plot. 
mean = []
for i in Qscore_list:
	mean.append(i/line_count)

base_pos = []
for i in range(args.read_length):
	base_pos.append(i)

plt.bar(base_pos,mean)
plt.title(args.output)
plt.xlabel('Nucleotide Position')
plt.ylabel('Mean Quality Score')
plt.savefig(args.output)

