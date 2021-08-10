#!/usr/bin/env python

import argparse
import Bioinfo
import gzip 
import re
import os

'''The utility of this script is to demultiplex Illumina paired-end reads given a set of dual-matched indexes and the four output files (Read 1 Insert 1, Insert 2, Read 2) produced by the sequencing run. The script will output Read1 and Read2 files for each sample group matched via the dual index as well as a R1/R2 file for both unknown/low quality and index hopped reads. The script will also output a stats.txt file containing information on the number of reads in each file.'''

#Define all functions that will be used in this Script
def get_args():
	parser = argparse.ArgumentParser(description="A script to demultiplex Illumina reads given a set of dual-matched indexes and the four output files (Read 1 Insert 1, Insert 2, Read 2) produced by the sequencing run. The script will output Read1 and Read2 files for each sample group matched via the dual index as well as a R1/R2 file for both unknown/low quality and index hopped reads.")
	parser.add_argument("-indexes", "--indexes_file", help = "path to tsv file containing column information: sample, group, treatment, index, index-sequence. NOTE: first line of file must be header information, all subsequent lines sample information", type=str, required=True)
	parser.add_argument("-r1", "--read1_file_name", help = "path to FASTQ for Read1 input file", type=str, required=True)
	parser.add_argument("-r2", "--read2_file_name", help = "path to FASTQ for Read1 input file", type=str, required=True)
	parser.add_argument("-i1", "--insert1_file_name", help = "path to FASTQ for Read1 input file", type=str, required=True)
	parser.add_argument("-i2", "--insert2_file_name", help = "path to FASTQ for Read1 input file", type=str, required=True)
	parser.add_argument("-cutoff", "--cutoff", help = "index Q-score cutoff value for 'unkown/lowqual binning. Default = 30", type=int, required=False, default=30)
	return parser.parse_args()

args = get_args()

#Define all functions that will be used in script
def index_dictionary(path_to_tsv_file: str) -> dict:
	'''function inputs a tsv file containing column information: sample, group, treatment, index, index-sequence and outputs a dictionary where the keys = unique index sequences and values = string of all other information (i.e sample, group, etc.). Note: the first line of the indexes file must be a header line, all subsequent lines sample information.'''
	index_dictionary: dict[str,str] = {}
	with open(path_to_tsv_file, "r") as indexes_file:
		next(indexes_file)
		for line in indexes_file.readlines():
			line = line.strip('\n')
			working_list = re.split("\t", line)
			index_dictionary[working_list[4]] = str(str(working_list[0])+"_"+str(working_list[1])+"_"+str(working_list[2])+"_"+str(working_list[3]))
	return index_dictionary

def reverse_complement(DNA_string: str) -> str:
	'''Function inputs a DNA sequence as a string and outputs the reverse complement of that sequence and returns a string of the reverse_complement sequence'''	
	string = DNA_string.upper().strip('\n')
	stringlength = len(string)
	reverse_string = string[stringlength::-1] 
	reverse_complement: str = ""
	for i in reverse_string:
		if i == "A":
			reverse_complement = str(reverse_complement) + "T"
		if i == "T":
			reverse_complement = str(reverse_complement) + "A"
		if i == "C":
			reverse_complement = str(reverse_complement) + "G"
		if i == "G":
			reverse_complement = str(reverse_complement) + "C"
	return reverse_complement

def Qscore_eval(Qscore_line: str,cutoff: int):
	'''Function inputs an Illumina Qscore sequence line (stripped of newline) and returns TRUE if the average Q-score of the line is at or above the cuttoff value AND there are no 'N's in the sequence and FALSE otherwise.'''
	Qscore_counter = 0
	if 'N' in Qscore_line:
		return False
	for i in Qscore_line:
		Qscore_counter += Bioinfo.convert_phred(i)
	if Qscore_counter//8 >= cutoff:
		return True
	else:
		return False
		
def Read1opener(filename,L1,L2,L3,L4):
	'''Function takes a file name and 4 lines, and writes those four lines to the Read 1 file'''	
	file_name = str('/projects/bgmp/zsisson2/bioinformatics/Bi622/Demultiplex/Assignment-the-third/output_files/')+str(filename)+str('_L008__R1_001.fq.gz')
	f = gzip.open(file_name,"at")
	f.write(L1)
	f.write(L2)
	f.write(L3)
	f.write(L4)
	return

def Read2opener(filename,L1,L2,L3,L4):
	'''Function takes a file name and 4 lines, and writes those four lines to the Read 2 file'''	
	file_name = str('/projects/bgmp/zsisson2/bioinformatics/Bi622/Demultiplex/Assignment-the-third/output_files/')+str(filename)+str('_L008__R2_001.fq.gz')
	f = gzip.open(file_name,"at")
	f.write(L1)
	f.write(L2)
	f.write(L3)
	f.write(L4)
	return

def Read1closer(filename):
	'''Function closes name of file'''	
	file_name = str('/projects/bgmp/zsisson2/bioinformatics/Bi622/Demultiplex/Assignment-the-third/output_fastq_files/' +str(filename)+ '_L008__R1_001.fq.gz')
	file_name.close()
	return

def Read2closer(filename):
	'''Function closes name of file'''	
	file_name = str('/projects/bgmp/zsisson2/bioinformatics/Bi622/Demultiplex/Assignment-the-third/output_fastq_files/'+str(filename)+ '_L008__R2_001.fq.gz')
	file_name.close()
	return


#generate index dictionary for reference. Note keys = unique index barcode and values = string of all other sample information(sample,group,etc)
index_dictionary = index_dictionary(args.indexes_file)

#Initialize counters for the number of read-pairs observed in each file:
matched_counter = 0
hopped_counter = 0
unknown_counter = 0
total_counter = 0

#Create a list of all possible index-matches that will be used as a counter
index_pair_list = []
Index_list = list(index_dictionary.keys())
for num in range(len(Index_list)):
	for i in range(len(Index_list)):
		index_pair_list.append(str(Index_list[i])+'-'+str(Index_list[num]))
		
#Initalize dictionary where the keys = unique index-pairs and the values = occurences, initialized to zero 
index_pair_dict = {}
for i in index_pair_list:
	index_pair_dict[i] = 0


#Open four input FASTQ files, and grab record from each file
Read1 = gzip.open(args.read1_file_name, "rt")
Read2 = gzip.open(args.read2_file_name, "rt")
Index1 = gzip.open(args.insert1_file_name, "rt")
Index2 = gzip.open(args.insert2_file_name, "rt")

R1line1 = Read1.readline()
R1line2 = Read1.readline()
R1line3 = Read1.readline()
R1line4 = Read1.readline()

I1line1 = Index1.readline()
I1line2 = Index1.readline()
I1line3 = Index1.readline()
I1line4 = Index1.readline()

I2line1 = Index2.readline()
I2line2 = Index2.readline()
I2line3 = Index2.readline()
I2line4 = Index2.readline()

R2line1 = Read2.readline()
R2line2 = Read2.readline()
R2line3 = Read2.readline()
R2line4 = Read2.readline()

#Perform this loop so long as there is a record present in memory:
while R1line4:
	#append the index1 and index2 sequence to the header lines for the reads sequences
	R1_header = str(R1line1.strip('\n') + ' ' + I1line2.strip('\n')+ '-' + I2line2)
	R2_header = str(R2line1.strip('\n') + ' ' + I1line2.strip('\n')+ '-' + I2line2)
	#Determine avg quality score of index lines. If either of the indexes in the pair contain an unknown base call (N) or an avg Q-score below the cutoff then bin to unknown. 
	if Qscore_eval(I1line4,30) == False or Qscore_eval(I2line4,30) == False:
		Read1opener("Unknown",R1_header,R1line2,R1line3,R1line4)
		Read2opener("Unknown",R2_header,R2line2,R2line3,R2line4)
		unknown_counter += 1
		total_counter += 1

	else:
		#Check to see if Index 2 is a perfect reverse complement of index 1 and if the index is in the reference list, If so, bin the reads to the proper file.
		rev_comp = reverse_complement(I2line2)
		if I1line2 == rev_comp and I1line2 in index_dictionary.keys():
			Read1opener(index_dictionary[I1line2],R1_header,R1line2,R1line3,R1line4)
			Read2opener(index_dictionary[I1line2],R2_header,R2line2,R2line3,R2line4)
			matched_counter += 1
			total_counter += 1
			#Increments specific index-pair counter
			key_to_search = str(I1line2)+'-'+str(rev_comp)
			index_pair_dict[key_to_search] += 1

		
		#If the indexes are perfect complements but are not found in the reference list, bin to unknown.
		if I1line2 == rev_comp and I1line2 not in index_dictionary.keys():
			Read1opener("Unknown",R1_header,R1line2,R1line3,R1line4)
			Read2opener("Unknown",R2_header,R2line2,R2line3,R2line4)
			unknown_counter += 1
			total_counter += 1
		
		#If the indexes are NOT perfect reverse complements, BUT they are both found in the reference list, assign to index-hopped file.
		if  I1line2 != rev_comp and I1line2 in index_dictionary.keys() and rev_comp in index_dictionary.keys():
			Read1opener("Swapped",R1_header,R1line2,R1line3,R1line4)
			Read2opener("Swapped",R2_header,R2line2,R2line3,R2line4)
			hopped_counter += 1
			total_counter += 1
			#Increments specific index-pair counter
			key_to_search = str(I1line2)+'-'+str(rev_comp)
			index_pair_dict[key_to_search] += 1
		
		#If none of the above then bin to unknown.
		else:
			Read1opener("Unknown",R1_header,R1line2,R1line3,R1line4)
			Read2opener("Unknown",R2_header,R2line2,R2line3,R2line4)
			unknown_counter += 1
			total_counter += 1

#Grab next record for each file and continue
	R1line1 = Read1.readline()
	R1line2 = Read1.readline()
	R1line3 = Read1.readline()
	R1line4 = Read1.readline()

	I1line1 = Index1.readline()
	I1line2 = Index1.readline()
	I1line3 = Index1.readline()
	I1line4 = Index1.readline()

	I2line1 = Index2.readline()
	I2line2 = Index2.readline()
	I2line3 = Index2.readline()
	I2line4 = Index2.readline()

	R2line1 = Read2.readline()
	R2line2 = Read2.readline()
	R2line3 = Read2.readline()
	R2line4 = Read2.readline()

#Close all files
# for key, value in index_dictionary.items():
# 	Read1closer(index_dictionary[key])
# 	Read2closer(index_dictionary[key])

#sort index-paired dictionary numerically by value and calculate total read-pairs for matches (hopped and dual-matched)
total_matched_hopped = hopped_counter+matched_counter
sorted(index_pair_dict, key=lambda i: int(index_pair_dict[i]))

#Generate report of file contents
with open('/home/zsisson2/bgmp/bioinformatics/Bi622/Demultiplex/Assignment-the-third/output_files/stats.md',"w") as fto:
	fto.write(str('Statistic' + '\t' + 'Value' + '\t' + 'Percentage of Whole'+'\n'))
	fto.write(str('Total Number of Read-Pairs' + '\t' + str(total_counter) + '\t' + str((total_counter/total_counter)*100) +'%'+'\n'))
	fto.write(str('Number of Unknown Read-Pairs' + '\t' + str(unknown_counter) + '\t' + str((unknown_counter/total_counter)*100) +'%'+'\n'))
	fto.write(str('Number of Index-Hopped Read-Pairs' + '\t' + str(hopped_counter) + '\t' + str((hopped_counter/total_counter)*100) +'%'+'\n'))
	fto.write(str('Number of Dual-Matched Read-Pairs' + '\t' + str(matched_counter) + '\t' + str((matched_counter/total_counter)*100) +'%'+'\n'))
	fto.write(str('\n'))
	fto.write(str('Breakdown of Matched Pairs:'))
	fto.write(str('\n'))
	fto.write(str('Index Pair'+'\t'+'Number of Read-Pairs'+'\t'+'Percentage of Total Matched/Hopped Pair'+ 'Percentage of Total'))
	for key, value in index_pair_dict.items():
		fto.write(key)
		fto.write('\t')
		fto.write(str(value))
		fto.write('\t')
		fto.write(str((value/total_matched_hopped)*100))
		fto.write('%')
		fto.write('\t')
		fto.write(str((value/total_counter)*100))
		fto.write('%')
		fto.write('\n')




