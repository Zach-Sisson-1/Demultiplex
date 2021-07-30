# Assignment the First

## Part 1
1. Be sure to upload your Python script.
'deplexer.py'

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read1 |
| 1294_S1_L008_R2_001.fastq.gz | Index1 |
| 1294_S1_L008_R3_001.fastq.gz | Index2 |
| 1294_S1_L008_R4_001.fastq.gz |  Read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
        
        ![](Read_1.png)

        ![](Read_2.png)

        ![](Index_1.png)

        ![](Index_2.png)

    2. ```Looking at the mean distributions, the majority of base calls within the indexes are at or above a Qscore of 30, indicating an confidene call of 99.9%. Given this distribution, A good quality score cutoff for index reads should be Q30, that is to say any base call within an index sequence that is below Q30 should be counted as unknown. The Qscore cutoff for biological read pairs can be slightly lower, at Q20 (99% confidence). Since there are very few base-pairs (8) in the indexes, it takes fewer errors to have a large impact (i.e match the entire sequence to the incorrect group) and so the Qscore cutoff should be more stringent to ensure the data that is grouped together, belongs together, even if that means the data is slightly more sparse. On the otherhand, A Qscore of 20 indicates a 1 in 100 (99%) chance that the base call is accurate. Since downstream sequence analysis involves larger reads than the indexes, the cutoff can be relaxed because the impact of any one base being incorrect drops with a longer sequence since there are more bases on either side to validate a sequence match than in the index reads.```
    3. ```Total indexes with (N) base calls = 7,304,664 ```
        
        ```Bash command:```
            
            (base) [zsisson2@n226 2017_sequencing]$ zcat 1294_S1_L008_R2_001.fastq.gz  1294_S1_L008_R3_001.fastq.gz | grep -v -e '^@' -e '^+' | grep 'N' | wc -l
            >>7304664
## Part 2
1. Define the problem:

    Illumina output FASTQ files contain multiplexed sequence data from different biolgical samples. Here, we need to sort each record within the FASTQ files to be reassigned to files contaning only the FASTQ records for that biological sample group. To do this, paired-end sequence reads must be sorted using their unique dual-matched index barcode, however, the code needs to account for barcodes having undetermined base calls (N) as well as index hopping, which results in some reads containing two different indexes on either side of the sequence. 
2. Describe output:

    Given an input of four files, two containing sequence reads(R1/R2) and two containing index reads(I1/I2), the most informative output file structure would contain 52 total files with the following content:
    
    - 48 FASTQ record-containing files, each containing records with the following header structure: (original header) + (I1-I2); where I1/I2 is the index sequence. Ex = @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1 GTAGCGTA-TACGCTAC *note I2 should be reverse compliment of I1  
	
	    - 24 for Read 1; file name structure = sample_group_treatment_index_L008_R1_001.fastq
	
	    - 24 for Read 2; file name structure = sample_group_treatment_index_L008_R2_001.fastq
 
     - 2 FASTQ record-containing files for sequences with unknown or low quality, each containing records with the following header structure: (original header) + (I1-I2); where I1/I2 is the index sequence. Ex = @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1 GTANCGTA-NAGCCATG    

	    - 1 for Read 1; file name structure = Unknown_L008_R1_001.fastq

	    - 1 for Read 2; file name structure = Unknown_L008_R2_001.fastq

    - 2 FASTQ record-containing files for sequences with non-matching indexes (index hopped), each containing records with the following header structure: (original header) + (I1-I2); where I1/I2 is the index sequence. Ex = @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1 GAACGTA-GAGCGTA
	
	    - 1 for Read 1; file name structure = Swapped_L008_R1_001.fastq

	    - 1 for Read 2; file name structure = Swapped_L008_R2_001.fastq

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

    done

4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement



Pseudocode:

0) using given index file, loop through and create a dictionary where keys = indexes and the values = field info (i.e sample_group_treatment_index)

1) Open all 52 files for writing at once:
    Open 48 by looping through dictionary, assining R1/R2 file for each key/value pair
    Open 2 by hardcoding Unknown Read1/Read2
    Open 2 by hardcoding Swapped Read1/Read2

2) Open all four  FASTQ files for reading at once as Read1,Index1,Index2,Read2

Loop through each input file, storing first four lines (stripped of '\n') in memory for each file:
    
    R1line1
    R1line2
    R1line3
    R1line4

    I1line1
    I1line2
    I1line3
    I1line4

    I2line1
    I2line2
    I2line3
    I2line4

    R2line1
    R2line2
    R2line3
    R2line4
    

3) Carry out the following instructions and evaluate conditions for each record set of 16 lines:
     
    Count the number of N's in I1line2 and I2line2 and keep a record of the number (this will be used to filter out any indexes with 'N' base calls)

    For the pair, if this number >= 1: (if there are any N's)
    assign record to unknown file via writing:

        string(R1line1 + I1line2-I2line2) >> 'Unknown_L008_R1_001.fastq'  *note that this writes the original header + the index pair to new file
        R1line2 >> 'Unknown_L008_R1_001.fastq'
        R1line3 >> 'Unknown_L008_R1_001.fastq'
        R1line4 >> 'Unknown_L008_R1_001.fastq'
        and
        string(R2line1 + I1line2-I2line2) >> 'Unknown_L008_R2_001.fastq'  *note that this writes the original header + the index pair to new file
        R2line2 >> 'Unknown_L008_R2_001.fastq'
        R2line3 >> 'Unknown_L008_R2_001.fastq'
        R2line4 >> 'Unknown_L008_R2_001.fastq'


If this number == 0:

- calculate/store reverse compliment of Index1 (I1line2) (see FXNS)
- using a regex expression, search for the match of the reverse compliment of Index1 in Index2 (I2line2), matching only if each character is identical to the pattern. 

    If match = True AND Index1 is in dictionary.keys:
        Indexes match and should be assigned to proper group
        assign that R1/R2 record by writing to appropriate file via:

        string(R1line1 + I1line2-I2line2) >> 'dictionary[index1 match]_L008_R1_001.fastq'  *note that this writes the original header + the index pair to new file
        R1line2 >> 'dictionary[index1 match]_L008_R1_001.fastq'
        R1line3 >> 'dictionary[index1 match]_L008_R1_001.fastq'
        R1line4 >> 'dictionary[index1 match]_L008_R1_001.fastq'
        and
        string(R2line1 + I1line2-I2line2) >> 'dictionary[index1 match]_L008_R2_001.fastq'  *note that this writes the original header + the index pair to new file
        R2line2 >> dictionary[index1 match]_L008_R2_001.fastq'
        R2line3 >> dictionary[index1 match]_L008_R2_001.fastq'
        R2line4 >> dictionary[index1 match]_L008_R2_001.fastq'
        
    If match = False AND Index1 is in dictionary.keys() AND reverse compliment of Index2 is in dictionary.keys()
        Indexes must be swapped, assign R1/R2 record to swapped file via:

        string(R1line1 + I1line2-I2line2) >> 'Swapped_L008_R1_001.fastq'  *note that this writes the original header + the index pair to new file
        R1line2 >> 'Swapped_L008_R1_001.fastq'
        R1line3 >> 'Swapped_L008_R1_001.fastq'
        R1line4 >> 'Swapped_L008_R1_001.fastq'
        and
        string(R1line1 + I1line2-I2line2) >> 'Swapped_L008_R2_001.fastq'  *note that this writes the original header + the index pair to new file
        R2line2 >> 'Swapped_L008_R2_001.fastq'
        R2line3 >> 'Swapped_L008_R2_001.fastq'
        R2line4 >> 'Swapped_L008_R2_001.fastq'
    
    Else: Assign records to unknown file using code from above
        

4) Refresh 16 line record memory, reset counters for Qscore, increment line-grabbing counter and grab new record for each file
5) Iterate until end of file


Functions: 

def Index_dictionary(file):

```Function inputs an index file name/path, with the following file structure:```

    sample  group   treatment       index   index sequence
    1       2A      control B1      GTAGCGTA
    2       2B      control A5      CGATCGAT
    3       2B      control C1      GATCAAGG
   ```and loops through each line, assigning the index sequences to a index_dictionary.keys() and the field information as the values.  returns the dictionary```

return Index_dict
    
    Input: 'indexes.txt'
    Expected output: {'GTAGCGTA': '1_2A_control_B1', 'CGATCGAT': '2_2B_control_A5', 'GATCAAGG': '3_2B_control_C1'}


def reverse_compliment(DNAstring):

```Function inputs a DNA sequence as a string and outputs the reverse compliment of that sequence and returns a string of the reverse_compliment sequence```

return reverse_compliment
    
    Input: 'ATCG'
    Expected output: 'CGAT'


def N_counter(string1,string2): 

```Counts the total number of N's in string1, and string2, which will represent Index lines (I1line2 and I2line2) and returns an integer value```

return N_count

    Input:(ATCGN,NNTCG)
    Expected output: 3

    Input:(ATCGC,CCTCG)
    Expected output: 0


def rc_patternmatch(string1,string2):

```Function will input two strings, which will represent Index1 and Index 2 sequence lines and compare the reverse compliment of Index1 with Index 2, returning True/False if the pattern match is identical```

return boolean

    Input:(ATCG,CGAT)
    Expected output: True

    Input:(ATCG,CCCT)
    Expected output: False
