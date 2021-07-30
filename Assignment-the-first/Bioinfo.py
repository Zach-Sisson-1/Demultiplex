#!/usr/bin/env python

def validate_base_seq(seq: str,RNAflag: bool = False)-> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

DNA = set('ATCGatcg')

#will only run for testing purposes when THIS script is run. Will NOT run when imported as a module. 
if __name__ == "__main__":
	assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
	assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
	assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
	assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
	print("Passed DNA and RNA tests")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    gc_counter = 0
    for nucleotide in str(DNA):
        if nucleotide == "G":
            gc_counter += 1
        if nucleotide == "C":
            gc_counter += 1
        gc_content = gc_counter/len(DNA)
    return gc_content

#will only run for testing purposes when THIS script is run. Will NOT run when imported as a module.
if __name__ == "__main__":
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")


def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred = (ord(letter))-33
    return phred

#will only run for testing purposes when THIS script is run. Will NOT run when imported as a module.
if __name__ == "__main__":
    assert convert_phred('A') == 32
    assert convert_phred('G') == 38
    assert convert_phred('V') == 53
    print("correctly converted phred score")

def qual_score(phred_score):
    """Takes the original, unmodified phred_score string as a parameter. This function should calculate the average quality score of the whole phred string."""
    total = 0
    for letter in phred_score:
        total += convert_phred(letter)
    return total/(len(phred_score))

#will only run for testing purposes when THIS script is run. Will NOT run when imported as a module.
if __name__ == "__main__":
    assert qual_score('KKKK') == 42
    assert qual_score('kkkk') == 74
    assert qual_score('KKkk') == 58
    assert qual_score('?vs&6') == 44.6
    print("correctly converted string into a mean phred score")



#Combines the upper two
def phred_score(string, output ='none'):
    """Takes a string of Qscores and returns a list of phred scores for each index. if 'avg' is passed in as the second argument, function returns an average quality score for the entire string."""
    phred_list=[]
    for char in string:
        phred_list.append(convert_phred(char))
    
    #Calculates average phred score of entire string if specified
    output = output.lower()
    if output == 'none':
         return phred_list
    if output == 'avg':
        return qual_score(string)


def barcode_counter(file_name, output = 'pts'):
    """Inputs a fastq file as a string (located in same directory) in which the sequence line is the 4th line of every entry and returns a dictionary count of each of the 5 letter barcodes to the screen. Note -- if "ptf" is passed in as second argument, the function will output the barcode as a .tsv file in the current directory."""  
    barcodes = {}
    with open('%s' % file_name, "r") as file:
        index = 5
        for line in file:
            i = line[0:5]
            if index % 4 == 2:
                if i in barcodes:
                    barcodes[i] += 1
                if i not in barcodes:
                    barcodes[i] = 1   
            index += 1
        dictlist = sorted(barcodes.items(), key=lambda x:x[1], reverse=True)
        barcodes = dict(dictlist)
   
    #Prints to screen or saves to file
    output = output.lower()
    if output == 'pts':
        return barcodes
    if output == 'ptf':
        with open('./%s.counted.csv' % file_name, "w") as file2:
            for key,value in barcodes.items():
                file2.write(key)
                file2.write("\t")
                file2.write(str(value))
                file2.write("\n")

def unwrap_fasta(file):
	"""Function inputs a fasta file as a string and outputs the file unwrapped (i.e one line for header and one line for entire sequence) in same directory with name file.unwrapped"""
	with open(file, "r") as fh:
		with open(str(file+'.unwrapped'),"w") as fto:
			for line in fh:
				line = line.strip('\n')
				if '>' in line:
					fto.writelines('\n'+line+'\n')
				else:
					fto.writelines(line)


def kspec(kmer_size: int, read_length_size:int, file_name:str, file_path:str):
    """Script takes input of a FASTQ file and path and creates a k-mer spectra of the frequency and distribution of k-mers of length 'kmer_size' and read length size 'read_length_size' found in the file."""
    #Initialize dictionary that will contain each k-mer {keys} of size k, and the count of how many times {value} that kmer is found in the sequence reads.
    kmer_dict:dict[str,int] = {}

    #loop through file grabbing each sequence line
    with open(str(file_path+file_name), "r") as fh:
        index:int = 3
        for line in fh:
            if index % 4 == 0:
                i:int = 0
                j:int = kmer_size
                for i in range(read_length_size-kmer_size+1):		
                    #Start of line by line command - checks if kmer already exists, if so, increment tally, if not, set new key equal to 1
                    if line[i:j] in kmer_dict.keys(): 
                        kmer_dict[line[i:j]] += 1
                    else:  
                        kmer_dict[line[i:j]] = 1
                    i += 1
                    j += 1
            index += 1
            #Check to see if script is still running
            if index % 4000000 ==0:
                print("...still running...")


    #Initialize empty dictionary to hold summary of k-mer counts. {Keys} = number of times a unique kmer was found throughout the sequence reads. {Values} represent the count of how many unique kmers were found {key} times throughout the sequence reads. 
    kmer_freq_dict:dict[int,int] = {}

    for key,value in kmer_dict.items():
        if kmer_dict[key] in kmer_freq_dict.keys():
            kmer_freq_dict[value] += 1
        else:
            kmer_freq_dict[value] = 1
        

    #assigns keys in freq_dict to lst1 and values to lst2.
    lst1:list[int]=[]
    lst2:list[int]=[]
    for key,value in kmer_freq_dict.items():
        lst1.append(key)
        lst2.append(value)


    #graph the frequency distributions
    import matplotlib.pyplot as plt 


    # plt.gca().set_xlim(right=10000)
    plt.bar(lst1,lst2)    
    plt.xlabel('k-mer frequency')
    plt.ylabel('Number of k-mers in this category')
    title=str('K-mer spectrum for length of k ='+str(kmer_size))
    plt.title(title)
    plt.yscale("log")

    plt.show()
    return


def knorm(kmer_size:int, coverage_limit:int, read_length_size:int, file_name:str, file_path:str, outputfile_name:str):
    """This script inputs a FASTQ file from a specified path, along with a k-mer size variable \
    and coverage limit variable. With these variables, the script will normalize the k-mer coverage of the FASTQ file and \
    output a new FASTQ file with only the records that fall at or below the coverage limit."""
    #Initialize empty dictionary. Keys = k-mers, values = number of occurences of each k-mer.
    kmer_norm_dict: dict[str,int] = {}

    #Opens output file to write and input file to read
    with open(outputfile_name, "w") as fto:
        intcounter: int = 0
        running_counter: int = 0
        with open(str(file_path+file_name), "r") as fh:
            #Creates strings out of the first four lines of the input file for use
            for line in fh:
                line =line.strip('\n')
                if intcounter==0:
                    str1: str = str(line[0:read_length_size])
                if intcounter==1:
                    str2: str = str(line[0:read_length_size])
                if intcounter==2:
                    str3: str = str(line[0:read_length_size])
                if intcounter==3:
                    str4: str = str(line[0:read_length_size])

                #Do stuff with strings created from each loop:
                if intcounter ==3:
                    #slides down sequence line by length of k, k-merizing the sequence and incrementing the kmer_norm_dict for each new kmer matched.
                    slider: int =0
                    for i in range(len(str2)-kmer_size):
                        if str2[slider:(slider+kmer_size)] in kmer_norm_dict.keys():
                            kmer_norm_dict[str2[slider:(slider+kmer_size)]] += 1
                        else:
                            kmer_norm_dict[str2[slider:(slider+kmer_size)]] = 1
                        #increments slider
                        slider += 1
                    
                    #initiliaze empty list to store k-mer coverages for this read
                    kmer_lst: list[int] = []

                    #k-merizing the sequence read again, and for each kmer, appending the coverage value (from kmer_norm_dict) to kmer list, generating list of coverages
                    slider: int = 0
                    for i in range(len(str2)-kmer_size):
                        kmer_lst.append(kmer_norm_dict[str2[slider:(slider+kmer_size)]])
                        slider += 1

                    #Calculate median k-mer coverage of this read - stored as median_kmer
                    kmer_lst.sort()
                    kmer_lst_length:int = len(kmer_lst)
                    if kmer_lst_length % 2 == 0:
                        med1: int = kmer_lst[kmer_lst_length//2]
                        med2: int = kmer_lst[(kmer_lst_length//2)-1]
                        median_kmer: int = (med1+med2)//2
                    else:
                        median_kmer: int = kmer_lst[kmer_lst_length//2]

                    #Checks median k-mer coverage of read. If equal to or below the coverage limit, writes the sequence line to output file.
                    if median_kmer <= coverage_limit:
                        fto.writelines(str1)
                        fto.writelines('\n')
                        fto.writelines(str2)
                        fto.writelines('\n')
                        fto.writelines(str3)
                        fto.writelines('\n')
                        fto.writelines(str4)
                        fto.writelines('\n')

                        
                #increments counter to fill all four strings
                intcounter +=1
                running_counter +=1 
                

                #Restarts the intcounter to get fresh set of four lines in the fastq file
                if intcounter == 4:
                    intcounter = 0
        
                #Checking to make sure program is still running
                if running_counter % 4000000 ==0:
                    print('still running')
    return


def Velvetparser(input_file:str, kmer_size:int):
    """Function takes an input fatsa file generated by Velvetg (as a string) and a kmer-size (as an int) and outputs two files. One file contains stats run on the contigs.fa file produced by Velvetg and the other is a plot of the distribution of contigs."""
    import re

    #Initialize an empty dictionary where the {keys} will be the node number of each contig an the {values} will be lists containing both the contig lenth, and k-mer coverage of each contig extracted from each header.
    contig_dict = {}

    with open(input_file, "r") as fh:
        for line in fh:
            line = line.strip('\n')

            #Selects only the header lines
            if re.match(">",line):
                #Pulls out the Node number, contig length and contig coverage from each header
                node_number = re.search('NODE_[0-9]+',line).group(0)
                #cleans node number variable to filter out text
                node_number = re.search('[0-9]+',node_number).group(0)
                length_contig = re.search('length_[0-9]+',line).group(0)
                #cleans length_contig variable to filter out text
                length_contig = re.search('[0-9]+',length_contig).group(0)
                contig_coverage = re.search('[0-9]+\.[0-9]+',line).group(0)

                contig_dict[node_number]=[length_contig,contig_coverage]

    #print(contig_dict)

    #Convert each length and coverage in the dictionary to integers and floats, respectively (from strings)	
    for key in contig_dict:
        contig_dict[key][0] = int(contig_dict[key][0])
        contig_dict[key][1] = float(contig_dict[key][1])

    #Converts each length into the physical contig length by adding (kmer size - 1) to each length value (taken from manual pg 15, first paragraph)
    for key in contig_dict:
        contig_dict[key][0] = (contig_dict[key][0] + kmer_size - 1)

    #print(contig_dict) 

    #Calculate number of contigs
    number_of_contigs: int = 0
    for key in contig_dict:
        number_of_contigs += 1


    #Calculate maximum contig length
    maximum_contig_length: int = 0
    for key in contig_dict:
        if contig_dict[key][0] > maximum_contig_length:
            maximum_contig_length = contig_dict[key][0]

    #Calculate mean contig length
    contig_length_sum: int = 0
    for key in contig_dict:
        contig_length_sum += contig_dict[key][0]
    mean_contig_length:float = contig_length_sum/len(contig_dict) 

    #Calculate total length of genome assembly across the contigs
    total_length_of_genome_assembly = contig_length_sum

    #Calculate weighted mean depth of coverage for the contigs
    contig_coverage_sum : float = 0
    contig_total_length:int = 0
    for key in contig_dict:
        contig_total_length += contig_dict[key][0]
        contig_coverage_sum += (contig_dict[key][1] * contig_dict[key][0])
    mean_depth_of_contig_coverage:float = contig_coverage_sum/contig_total_length



    #Calculate the N50 of assembly
    #part 1 - extract all contig lengths into a list and sort it in descending numeric order
    contig_length_list:list = []
    for key in contig_dict:
        contig_length_list.append(contig_dict[key][0])
    contig_length_list.sort(reverse=True)


    #part 2 - continue summing the lenths of contigs in the list from largest to smallest until that sum reaches 50% or more of the total assembly length. 
    #the contig length associated with bumping the sum to, or past, 50% of the total assembly is the N50 value.
    running_sum: int = 0
    for i in contig_length_list:
        running_sum += i
        if running_sum >= total_length_of_genome_assembly/2:
            N50_value = i
            break
        else:
            continue

    #Calculate the distribution of contig lengths
    #part 1 -  - create and initialize a dictionary where the {keys} will be the values (100)(200)...(max contig length rounded down to hundreds) the {values} will be a count of how many contigs fall between the range(key + 99). The number in each key is inclusive so (0) contains contigs of length (0-99).
    contig_distribution_dict = {}

    bin_integer_counter:int = 0
    for i in range((contig_length_list[0]//100)+1):
        contig_distribution_dict[bin_integer_counter] = 0
        bin_integer_counter += 100

    #part 2 - loop through list of contig lengths, assigning them to their appropriate bin.
    for contig in contig_length_list:
        for keybin in contig_distribution_dict:
            if contig in range(keybin,keybin+100):
                contig_distribution_dict[keybin] += 1
                break
            else:
                continue


    #Plot the distribution 
    import matplotlib.pyplot as plt 

    #create list for bins used in plot. Bins will take the form [0-100),[100-200),[200-300)...up to bin containing the maximum contig length(rounded up).
    bin_integer_plot:int = 100
    bin_list =[0]
    for i in range((contig_length_list[0]//100) +1):
        bin_list.append(bin_integer_plot)
        bin_integer_plot += 100

    title_str = str('Distribution of contig frequency at k = '+ str(kmer_size))
    plt.hist(contig_length_list,bins=bin_list, rwidth=0.4)  
    # plt.gca().set_xlim(left=0,right=10000)
    plt.xlabel('Contig length (in bp)')
    plt.ylabel('Frequency')
    plt.yscale("log")
    plt.title(title_str)
    plt.savefig(str(title_str+'.png'))


    #prints results in organized table
    # def print_table():
    # 	print('Metric:','\t','Value')
    # 	print('Number of contigs:','\t', number_of_contigs)
    # 	print('Max contig length:','\t',maximum_contig_length)
    # 	print('Mean contig length:','\t',mean_contig_length)
    # 	print('Total length of genome assembly:','\t',total_length_of_genome_assembly)
    # 	print('Mean depth of contig coverage:','\t',mean_depth_of_contig_coverage)
    # 	print('N50 value:','\t',N50_value)
    # 	return(print)

    #prints results to file in organized table
    with open('Velvetparser_output',"w") as fto:
        fto.write('Metric:'+'\t'+'Value'+'\n')
        fto.write('Number of contigs:'+'\t'+ str(number_of_contigs)+'\n')
        fto.write('Max contig length:'+'\t'+str(maximum_contig_length)+'\n')
        fto.write('Mean contig length:'+'\t'+str(mean_contig_length)+'\n')
        fto.write('Total length of genome assembly:'+'\t'+str(total_length_of_genome_assembly)+'\n')
        fto.write('Mean depth of contig coverage:'+'\t'+str(mean_depth_of_contig_coverage)+'\n')
        fto.write('N50 value:'+'\t'+str(N50_value)+'\n')
        fto.write('\n')
        fto.write('# Contig length'+'\t'+'Number of contigs in this category'+'\n')
        for key,value in contig_distribution_dict.items():
            fto.write(str(key)+'\t'+str(value)+'\n')

    return
