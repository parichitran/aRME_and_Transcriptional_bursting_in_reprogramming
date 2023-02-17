library(dplyr)
#Getting allelic count files

AlleleA=read.table('../day0/129S1_filtered_grt_10_f25cells_exp.tsv',sep = "\t",header=TRUE)
AlleleB=read.table('../day0/CAST_filtered_grt_10_f25cells_exp.tsv', sep="\t",header=TRUE)

AlleleA=data.frame(AlleleA, row.names="ID")
AlleleB=data.frame(AlleleB, row.names="ID")

genename <- rownames(AlleleA)
sampname <- colnames(AlleleA)

N = apply(AlleleA + AlleleB, 2, sum)
lib.size = N/mean(N)


#spikein_read=spike_in
#remove cells with extreme library size factors
lib.size.filter=(lib.size>0.85 & lib.size < 1.3)
AlleleA=AlleleA[,lib.size.filter]; AlleleB=AlleleB[,lib.size.filter]

total <- AlleleA+AlleleB

#REmoving non expressing genes
total <- total[rowSums(total[])>0,]
AlleleA <- AlleleA [row.names(total),]
AlleleB <- AlleleB [row.names(total),]

#calculating ratios

AlleleA <- (AlleleA/total)*100
AlleleB <- (AlleleB/total)*100

write.table(AlleleA,file="AlleleAratio.txt",sep="\t")
write.table(AlleleB,file="AlleleBratio.txt",sep="\t")
rm(list=ls())

#Catogorise with the next python script
