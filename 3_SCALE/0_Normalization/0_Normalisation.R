library(GenomicFeatures)
library(data.table)
rm(list=ls())
txdb <- makeTxDbFromGFF('/full/path/to/Mus_musculus.GRCm38.102.gtf')

#whole gene length
#all.genes <- genes(txdb)
#my.genes.lengths <- width(all.genes[all.genes$gene_id])

#names(my.genes.lengths) <-  all.genes$gene_id
#h <- data.frame(my.genes.lengths)

#Method1
#https://www.biostars.org/p/83901/
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

leng <- data.frame(exonic.gene.sizes)
leng <- data.frame(row.names(leng),leng)
#Method2
#https://www.biostars.org/p/185665/
#exonic <- exonsBy(txdb, by="gene")
#red.exonic <- reduce(exonic)
#exon.lengths <- vapply(width(red.exonic), sum, numeric(1))

#Getting this file contains previously processed cells.So the consistency is maintained(2 or 3 cells came extra when whole data processed in this system).

#allelic <- read.table('/full/path/to//01_Rawdata/Test/mergeallelicnewipsc/burstipsc/129S1_allele.txt', header=T,sep = '\t',row.names=1)

#Getiing input of the Allelic counts that should be normalised
AlleleA <- read.table(file="/2_Allelic/1_AllelicMappingandCounting/129S1_SvImJ.txt",sep="\t",header=T,row.names=1)
AlleleB <- read.table(file="/2_Allelic/1_AllelicMappingandCounting/CAST.txt",sep="\t",header=T,row.names=1)

#AlleleB<- AlleleB[row.names(AlleleA),colnames(AlleleA)]
total <- AlleleA+AlleleB

total <- total[,colSums(total) >= 200]


#AlleleA <- AlleleA[rowSums(is.na(AlleleA))!= ncol(AlleleA), ]
#AlleleA[is.na(AlleleA)] = 0

#only subsample required cells
AlleleA <- AlleleA[,colnames(total)]

#order the genelength accordingly
l <- leng[row.names(AlleleA),]

#Making the vector useful for calculation
l <- l[,2] 

#Total mapped reads per sample
cSA <- colSums(AlleleA) 
#https://www.biostars.org/p/273537/ #Formula for RPKM
#https://www.biostars.org/p/96084/ #script for normalisation
#RPKM for AlleleA
rpkmA <- (10^9)*t(t(AlleleA/l)/cSA)

#Similarly done for AlleleB
AlleleB <- AlleleB[,colnames(total)]

#Again gene length is assigned for alleleB
lb <- leng[row.names(AlleleB),]

lb <- lb[,2] 

cSB <- colSums(AlleleB) #Total mapped reads per sample

rpkmB <- (10^9)*t(t(AlleleB/lb)/cSB)

#trp <- rpkmA+rpkmB

#v <- data.frame(rowSums(trp))


write.table(rpkmA, file="/full/path/to/rpkmA.tsv", sep="\t")
write.table(rpkmB, file="/full/path/to/rpkmB.tsv", sep="\t")

#Get cell data file to sort cells
Cell_info=read.table("/0_Reads/meta_NA626.tsv", sep="\t", header=T)

#check
head(Cell_info)

#subsample whole required
Cell_info <- Cell_info[colnames(total),]

#subsample to day point
sub <- Cell_info[which(Cell_info$Timepoint=='iPSCs'),]


norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

#Read geneid to chrom map file to remove allosomes
gtf2 <- read.table('/full/path/to//mouse.txt', header = T, sep = '\t',row.names = 1)

#Getting X & Y chromosome genes
x <- gtf2[gtf2$Chromosome %like% "X:", ]
y <- gtf2[gtf2$Chromosome %like% "Y:", ]

#Getting Autsomes
b <-gtf2[!(row.names(gtf2) %in% row.names(x)),]
a <-b[!(row.names(b) %in% row.names(y)),]

#Getting only Autosomal genes
nix <- intersect(row.names(a),row.names(norm_AlleleA))

#Only sampling autosomes
xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nix),]

write.table(xAi, file="/full/path/to/ipsc/ipsc129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/ipsc/ipsccastrpkm.tsv", sep="\t")

############################################################ DAY0 ################################################
sub <- Cell_info[which(Cell_info$Timepoint=='Day_0'),]
norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

nix <- intersect(row.names(a),row.names(norm_AlleleA))

xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nix),]

write.table(xAi, file="/full/path/to/day0/Day0129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/day0/Day0castrpkm.tsv", sep="\t")

######################################################  Day8 ########################################################
sub <- Cell_info[which(Cell_info$Timepoint=='Day_8'),]
norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

nix <- intersect(row.names(a),row.names(norm_AlleleA))
nixb <- intersect(row.names(a),row.names(norm_AlleleB))

xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nixb),]
#xAi[is.na(xAi)] = 0
#xBi[is.na(xBi)] = 0



write.table(xAi, file="/full/path/to/day8/Day8129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/day8/Day8castrpkm.tsv", sep="\t")
######################################################## Day9
sub <- Cell_info[which(Cell_info$Timepoint=='Day_9'),]
norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

nix <- intersect(row.names(a),row.names(norm_AlleleA))

xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nix),]


write.table(xAi, file="/full/path/to/day9/Day9129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/day9/Day9castrpkm.tsv", sep="\t")
######################################################## Day 10
sub <- Cell_info[which(Cell_info$Timepoint=='Day_10'),]
norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

nix <- intersect(row.names(a),row.names(norm_AlleleA))

xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nix),]

write.table(xAi, file="/full/path/to/day10/Day10129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/day10/Day10castrpkm.tsv", sep="\t")
######################################################## Day 12
sub <- Cell_info[which(Cell_info$Timepoint=='Day_12'),]
norm_AlleleA <- rpkmA[,row.names(sub)]
norm_AlleleB <- rpkmB[,row.names(sub)]

nix <- intersect(row.names(a),row.names(norm_AlleleA))

xAi <-  norm_AlleleA[c(nix),]
xBi <-  norm_AlleleB[c(nix),]

write.table(xAi, file="/full/path/to/day12/Day12129rpkm.tsv", sep="\t")

write.table(xBi, file="/full/path/to/day12/Day12castrpkm.tsv", sep="\t")
########################################################
