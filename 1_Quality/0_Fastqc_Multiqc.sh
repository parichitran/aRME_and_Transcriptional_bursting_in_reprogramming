#! /bin/bash

#Running fastqc for all files
fastqc -t 8 *.fastq.gz -o 01_Fastqc/

#Change to to fastqc outputs
cd 01_Fastqc/

#Merging with multiaqc
multiqc .
 
