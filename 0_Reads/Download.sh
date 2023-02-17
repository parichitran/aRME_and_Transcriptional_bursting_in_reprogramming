#! /bin/bash
#IFS-Input field separator to split the data by comma and read mode to read data from columns
while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column5 rec_column7  
do
#prints the read details  
  echo "$rec_column1"
#Downloads the ftp file by curl and outputs as desired file name
  curl -C - "$rec_column2" -o "$rec_column1".fastq.gz
  #echo "$rec_column1"_2
  #curl -C - "$rec_column3" -o "$rec_column1"_2.fastq.gz
  echo "$rec_column3"
  echo "$rec_column5"
  echo "$rec_column7"
#Tail to skip the header in csv file
done < <(tail -n +2 Data.csv)
