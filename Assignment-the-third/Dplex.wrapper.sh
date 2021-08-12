#!/usr/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-20:00:00
#SBATCH --output=DeplexerOUTPUT.%j
#SBATCH --error=Deplexer.OUTPUT.err

conda activate bgmp_py39
/usr/bin/time -v ./Dplexer.py -indexes '/home/zsisson2/bgmp/bioinformatics/Bi622/Demultiplex/Assignment-the-third/Testinput/indexes.txt' \
-r1 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' \
-r2 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz' \
-i1 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz' \
-i2 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' \
-cutoff 30
