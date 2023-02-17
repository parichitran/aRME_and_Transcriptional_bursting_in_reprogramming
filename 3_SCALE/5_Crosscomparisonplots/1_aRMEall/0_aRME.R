library(UpSetR)
#library(clusterProfiler)
#library(org.Mm.eg.db)
#Getting ratio file for genes
day0<- read.table(file="../Day0/AlleleAratio.txt",sep="\t",header=T)
day8<- read.table(file="../Day8/AlleleAratio.txt",sep="\t",header=T)
day9<- read.table(file="../Day9/AlleleAratio.txt",sep="\t",header=T)
day10<- read.table(file="../Day10/AlleleAratio.txt",sep="\t",header=T)
day12<- read.table(file="../Day12/AlleleAratio.txt",sep="\t",header=T)
ipsc<- read.table(file="../ipsc/AlleleAratio.txt",sep="\t",header=T)


#Common gene across day points
commongene <- Reduce(intersect, list(row.names(day0),row.names(day8),row.names(day9),row.names(day10),row.names(day12),row.names(ipsc)))


#Getting the genes wise
day0cat1 <- read.table(file="../Day0/Cat1",sep="\t",header=T)
day0cat1<- c(intersect(commongene,day0cat1$ID))
day0cat2 <- read.table(file="../Day0/Cat2",sep="\t",header=T)
day0cat2<- c(intersect(commongene,day0cat2$ID))
day0cat3 <- read.table(file="../Day0/Cat3",sep="\t",header=T)
day0cat3<- c(intersect(commongene,day0cat3$ID))
day0cat4 <- read.table(file="../Day0/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day0cat4<- c(intersect(commongene,day0cat4$ID))

day8cat1 <- read.table(file="../Day8/Cat1",sep="\t",header=T)
day8cat1<- c(intersect(commongene,day8cat1$ID))
day8cat2 <- read.table(file="../Day8/Cat2",sep="\t",header=T)
day8cat2<- c(intersect(commongene,day8cat2$ID))
day8cat3 <- read.table(file="../Day8/Cat3",sep="\t",header=T)
day8cat3<- c(intersect(commongene,day8cat3$ID))
day8cat4 <- read.table(file="../Day8/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day8cat4<- c(intersect(commongene,day8cat4$ID))

day9cat1 <- read.table(file="../Day9/Cat1",sep="\t",header=T)
day9cat1<- c(intersect(commongene,day9cat1$ID))
day9cat2 <- read.table(file="../Day9/Cat2",sep="\t",header=T)
day9cat2<- c(intersect(commongene,day9cat2$ID ))
day9cat3 <- read.table(file="../Day9/Cat3",sep="\t",header=T)
day9cat3<- c(intersect(commongene,day9cat3$ID))
day9cat4 <- read.table(file="../Day9/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day9cat4<- c(intersect(commongene,day9cat4$ID))

day10cat1 <- read.table(file="../Day10/Cat1",sep="\t",header=T)
day10cat1<- c(intersect(commongene,day10cat1$ID))
day10cat2 <- read.table(file="../Day10/Cat2",sep="\t",header=T)
day10cat2<- c(intersect(commongene,day10cat2$ID))
day10cat3 <- read.table(file="../Day10/Cat3",sep="\t",header=T)
day10cat3<- c(intersect(commongene,day10cat3$ID))
day10cat4 <- read.table(file="../Day10/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day10cat4<- c(intersect(commongene,day10cat4$ID))


day12cat1 <- read.table(file="../Day12/Cat1",sep="\t",header=T)
day12cat1<- c(intersect(commongene,day12cat1$ID))
day12cat2 <- read.table(file="../Day12/Cat2",sep="\t",header=T)
day12cat2<- c(intersect(commongene,day12cat2$ID))
day12cat3 <- read.table(file="../Day12/Cat3",sep="\t",header=T)
day12cat3<- c(intersect(commongene,day12cat3$ID))
day12cat4 <- read.table(file="../Day12/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day12cat4<- c(intersect(commongene,day12cat4$ID))

ipsccat1 <- read.table(file="../ipsc/Cat1",sep="\t",header=T)
ipsccat1<- c(intersect(commongene,ipsccat1$ID))
ipsccat2 <- read.table(file="../ipsc/Cat2",sep="\t",header=T)
ipsccat2<- c(intersect(commongene,ipsccat2$ID))
ipsccat3 <- read.table(file="../ipsc/Cat3",sep="\t",header=T)
ipsccat3<- c(intersect(commongene,ipsccat3$ID))
ipsccat4 <- read.table(file="../ipsc/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
ipsccat4<- c(intersect(commongene,ipsccat4$ID))


#Creating list for upset plot
listInput <- list(day0cat1 = day0cat1, 
                day8cat1=day8cat1,
                day9cat1=day9cat1,
                day10cat1=day10cat1,
                day12cat1=day12cat1,
                iPSCcat1=ipsccat1,
                
                day0cat2=day0cat2,
                day8cat2=day8cat2,
                day9cat2=day9cat2,
                day10cat2=day10cat2,
                day12cat2=day12cat2,
                iPSCcat2=ipsccat2,
                
                day0cat3=day0cat3,
                day8cat3=day8cat3,
                day9cat3=day9cat3,
                day10cat3=day10cat3,
                day12cat3=day12cat3,
                iPSCcat3=ipsccat3,
                
                day0cat4=day0cat4,
                day8cat4=day8cat4,
                day9cat4=day9cat4,
                day10cat4=day10cat4,
                day12cat4=day12cat4,
                iPSCcat4=ipsccat4)




tiff("01_armecommon", units="in", width=28, height=25, res=300)




names(listInput) <-  c("Day0-cat1","Day8-cat1","Day9-cat1","Day10-cat1","Day12-cat1","iPSCs-cat1","Day0-cat2","Day8-cat2","Day9-cat2","Day10-cat2","Day12-cat2","iPSCs-cat2","Day0-cat3","Day8-cat3","Day9-cat3","Day10-cat3","Day12-cat3","iPSCs-cat3","Day0-cat4","Day8-cat4","Day9-cat4","Day10-cat4","Day12-cat4","iPSCs-cat4")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

#k <- k[7:24]
upset(fromList(listInput),keep.order = T,sets = c(k),order.by = "freq", sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =5,point.size=9 ,mb.ratio = c(0.55, 0.45),sets.bar.color=c("chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","red","red","red","red","red","red"))


dev.off()

#..

#Writing intersect genes for downstream analysis
#allcat3
allcat3 <- Reduce(intersect, list(day0cat3,day8cat3,day9cat3,day10cat3,day12cat3,ipsccat3))

write.table(allcat3,file="..01_Allcat3.txt",sep="\t")

#allcat3ipsccat4
allcat3ipsccat4 <- Reduce(intersect, list(day0cat3,day8cat3,day9cat3,day10cat3,day12cat3,ipsccat4))

write.table(allcat3ipsccat4,file="..02_allcat3ipsccat4.txt",sep="\t")


#allcat2
allcat2 <- Reduce(intersect, list(day0cat2,day8cat2,day9cat2,day10cat2,day12cat2,ipsccat2))

write.table(allcat2,file="..03_allcat2.txt",sep="\t")

#allcat3ipsccat2
allcat3ipsccat2 <- Reduce(intersect, list(day0cat3,day8cat3,day9cat3,day10cat3,day12cat3,ipsccat2))

write.table(allcat3ipsccat2,file="..04_allcat3ipsccat2.txt",sep="\t")

#allcat4
allcat4 <- Reduce(intersect, list(day0cat4,day8cat4,day9cat4,day10cat4,day12cat4,ipsccat4))
write.table(allcat4,file="..05_allcat4.txt",sep="\t")

#allcat3day0cat4
allcat3day0cat4 <- Reduce(intersect, list(day0cat4,day8cat3,day9cat3,day10cat3,day12cat3,ipsccat3))

write.table(allcat3day0cat4,file="..06_allcat3day0cat4.txt",sep="\t")

#allcat1
allcat1 <- Reduce(intersect, list(day0cat1,day8cat1,day9cat1,day10cat1,day12cat1,ipsccat1))

write.table(allcat1,file="..07_allcat1.txt",sep="\t")


rm(list=ls())

#################################################### No common genes
#Getting coordinated
day0<- read.table(file="../Day0/AlleleAratio.txt",sep="\t",header=T)
day8<- read.table(file="../Day8/AlleleAratio.txt",sep="\t",header=T)
day9<- read.table(file="../Day9/AlleleAratio.txt",sep="\t",header=T)
day10<- read.table(file="../Day10/AlleleAratio.txt",sep="\t",header=T)
day12<- read.table(file="../Day12/AlleleAratio.txt",sep="\t",header=T)
ipsc<- read.table(file="../ipsc/AlleleAratio.txt",sep="\t",header=T)


day0cat1 <- read.table(file="../Day0/Cat1",sep="\t",header=T)
day0cat2 <- read.table(file="../Day0/Cat2",sep="\t",header=T)
day0cat3 <- read.table(file="../Day0/Cat3",sep="\t",header=T)
day0cat4 <- read.table(file="../Day0/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

day8cat1 <- read.table(file="../Day8/Cat1",sep="\t",header=T)
day8cat2 <- read.table(file="../Day8/Cat2",sep="\t",header=T)
day8cat3 <- read.table(file="../Day8/Cat3",sep="\t",header=T)
day8cat4 <- read.table(file="../Day8/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)


day9cat1 <- read.table(file="../Day9/Cat1",sep="\t",header=T)
day9cat2 <- read.table(file="../Day9/Cat2",sep="\t",header=T)
day9cat3 <- read.table(file="../Day9/Cat3",sep="\t",header=T)
day9cat4 <- read.table(file="../Day9/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)


day10cat1 <- read.table(file="../Day10/Cat1",sep="\t",header=T)
day10cat2 <- read.table(file="../Day10/Cat2",sep="\t",header=T)
day10cat3 <- read.table(file="../Day10/Cat3",sep="\t",header=T)
day10cat4 <- read.table(file="../Day10/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

day12cat1 <- read.table(file="../Day12/Cat1",sep="\t",header=T)
day12cat2 <- read.table(file="../Day12/Cat2",sep="\t",header=T)
day12cat3 <- read.table(file="../Day12/Cat3",sep="\t",header=T)
day12cat4 <- read.table(file="../Day12/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

ipsccat1 <- read.table(file="../ipsc/Cat1",sep="\t",header=T)
ipsccat2 <- read.table(file="../ipsc/Cat2",sep="\t",header=T)
ipsccat3 <- read.table(file="../ipsc/Cat3",sep="\t",header=T)
ipsccat4 <- read.table(file="../ipsc/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

listInput<- list(day0cat1$ID, 
                   day8cat1$ID,
                   day9cat1$ID,
                   day10cat1$ID,
                   day12cat1$ID,
                   ipsccat1$ID,
                   
                   day0cat2$ID,
                   day8cat2$ID,
                   day9cat2$ID,
                   day10cat2$ID,
                   day12cat2$ID,
                   ipsccat2$ID,
                   
                   day0cat3$ID,
                   day8cat3$ID,
                   day9cat3$ID,
                   day10cat3$ID,
                   day12cat3$ID,
                   ipsccat3$ID,
                   
                   day0cat4$ID,
                   day8cat4$ID,
                   day9cat4$ID,
                   day10cat4$ID,
                   day12cat4$ID,
                   ipsccat4$ID)
names(listInput) <-  c("Day0-cat1","Day8-cat1","Day9-cat1","Day10-cat1","Day12-cat1","iPSCs-cat1","Day0-cat2","Day8-cat2","Day9-cat2","Day10-cat2","Day12-cat2","iPSCs-cat2","Day0-cat3","Day8-cat3","Day9-cat3","Day10-cat3","Day12-cat3","iPSCs-cat3","Day0-cat4","Day8-cat4","Day9-cat4","Day10-cat4","Day12-cat4","iPSCs-cat4")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

tiff("02_upsetarmenocommon.tiff", units="in", width=30, height=25, res=300)

#k <- k[7:24]
upset(fromList(listInput),keep.order = T,sets = c(k),order.by = "freq", sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =5,point.size=7 ,mb.ratio = c(0.55, 0.45),nintersects=50,sets.bar.color=c("chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","red","red","red","red","red","red"))

dev.off()

rm(list=ls())

###################################################################### day0vs ipsc
#Getting coordinated
day0<- read.table(file="../Day0/AlleleAratio.txt",sep="\t",header=T)
day8<- read.table(file="../Day8/AlleleAratio.txt",sep="\t",header=T)
day9<- read.table(file="../Day9/AlleleAratio.txt",sep="\t",header=T)
day10<- read.table(file="../Day10/AlleleAratio.txt",sep="\t",header=T)
day12<- read.table(file="../Day12/AlleleAratio.txt",sep="\t",header=T)
ipsc<- read.table(file="../ipsc/AlleleAratio.txt",sep="\t",header=T)


day0cat1 <- read.table(file="../Day0/Cat1",sep="\t",header=T)
day0cat2 <- read.table(file="../Day0/Cat2",sep="\t",header=T)
day0cat3 <- read.table(file="../Day0/Cat3",sep="\t",header=T)
day0cat4 <- read.table(file="../Day0/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)


ipsccat1 <- read.table(file="../ipsc/Cat1",sep="\t",header=T)
ipsccat2 <- read.table(file="../ipsc/Cat2",sep="\t",header=T)
ipsccat3 <- read.table(file="../ipsc/Cat3",sep="\t",header=T)
ipsccat4 <- read.table(file="../ipsc/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

listInput <- list(day0cat1 = day0cat1$ID, 
                  
                  iPSCcat1=ipsccat1$ID,
                  
                  day0cat2=day0cat2$ID,
                  
                  iPSCcat2=ipsccat2$ID,
                  
                  day0cat3=day0cat3$ID,
                  
                  iPSCcat3=ipsccat3$ID,
                  
                  day0cat4=day0cat4$ID,
                  
                  iPSCcat4=ipsccat4$ID)

names(listInput) <-  c("Day0-cat1","iPSCs-cat1","Day0-cat2","iPSCs-cat2","Day0-cat3","iPSCs-cat3","Day0-cat4","iPSCs-cat4")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

tiff("03_upsetday0vsipsc.tiff", units="in", width=25, height=17, res=300)

#k <- k[7:24]
upset(fromList(listInput),order.by = "freq",keep.order = T,sets = c(k), sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =5,point.size=9 ,mb.ratio = c(0.55, 0.45),sets.bar.color=c("chocolate4","chocolate4","lightseagreen","lightseagreen","darkslategray4","darkslategray4","red","red"))

dev.off()

#ipsccat1day0cat2
ipsccat1day0cat2 <- Reduce(intersect, list(day0cat2$ID,ipsccat1$ID))

write.table(ipsccat1day0cat2,file="../aRME/Day0vsipsc/01_ipsccat1day0cat2",sep="\t")


#ipsccat2day0cat1
ipsccat2day0cat1 <- Reduce(intersect, list(day0cat1$ID,ipsccat2$ID))

write.table(ipsccat2day0cat1,file="../aRME/Day0vsipsc/02_ipsccat2day0cat1",sep="\t")



#ipsccat1day0cat3
ipsccat1day0cat3 <- Reduce(intersect, list(day0cat3$ID,ipsccat1$ID))

write.table(ipsccat1day0cat3,file="../aRME/Day0vsipsc/03_ipsccat1day0cat3",sep="\t")

#ipsccat3day0cat1
ipsccat3day0cat1 <- Reduce(intersect, list(day0cat1$ID,ipsccat3$ID))

write.table(ipsccat3day0cat1,file="../aRME/Day0vsipsc/04_ipsccat3day0cat1",sep="\t")
rm(list=ls())


#################################################### maternal vs paternal
#Getting coordinated
day0<- read.table(file="../Day0/AlleleAratio.txt",sep="\t",header=T)
day8<- read.table(file="../Day8/AlleleAratio.txt",sep="\t",header=T)
day9<- read.table(file="../Day9/AlleleAratio.txt",sep="\t",header=T)
day10<- read.table(file="../Day10/AlleleAratio.txt",sep="\t",header=T)
day12<- read.table(file="../Day12/AlleleAratio.txt",sep="\t",header=T)
ipsc<- read.table(file="../ipsc/AlleleAratio.txt",sep="\t",header=T)


day0cat1 <- read.table(file="../Day0/Cat1",sep="\t",header=T)
#paternal sampling
day0cat1a <- subset(day0cat1,day0cat1$CAST_perc_y==100)
#Maternal sampling
day0cat1b <- subset(day0cat1,day0cat1$CAST_perc_y==0)
day0cat2 <- read.table(file="../Day0/Cat2",sep="\t",header=T)
day0cat3 <- read.table(file="../Day0/Cat3",sep="\t",header=T)
day0cat4 <- read.table(file="../Day0/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

#Same as day0
day8cat1 <- read.table(file="../Day8/Cat1",sep="\t",header=T)
day8cat1a <- subset(day8cat1,day8cat1$CAST_perc_y==100)
day8cat1b <- subset(day8cat1,day8cat1$CAST_perc_y==0)
day8cat2 <- read.table(file="../Day8/Cat2",sep="\t",header=T)
day8cat3 <- read.table(file="../Day8/Cat3",sep="\t",header=T)
day8cat4 <- read.table(file="../Day8/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)


day9cat1 <- read.table(file="../Day9/Cat1",sep="\t",header=T)
day9cat1a <- subset(day9cat1,day9cat1$CAST_perc_y==100)
day9cat1b <- subset(day9cat1,day9cat1$CAST_perc_y==0)
day9cat2 <- read.table(file="../Day9/Cat2",sep="\t",header=T)
day9cat3 <- read.table(file="../Day9/Cat3",sep="\t",header=T)
day9cat4 <- read.table(file="../Day9/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)


day10cat1 <- read.table(file="../Day10/Cat1",sep="\t",header=T)
day10cat1a <- subset(day10cat1,day10cat1$CAST_perc_y==100)
day10cat1b <- subset(day10cat1,day10cat1$CAST_perc_y==0)
day10cat2 <- read.table(file="../Day10/Cat2",sep="\t",header=T)
day10cat3 <- read.table(file="../Day10/Cat3",sep="\t",header=T)
day10cat4 <- read.table(file="../Day10/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

day12cat1 <- read.table(file="../Day12/Cat1",sep="\t",header=T)
day12cat1a <- subset(day12cat1,day12cat1$CAST_perc_y==100)
day12cat1b <- subset(day12cat1,day12cat1$CAST_perc_y==0)
day12cat2 <- read.table(file="../Day12/Cat2",sep="\t",header=T)
day12cat3 <- read.table(file="../Day12/Cat3",sep="\t",header=T)
day12cat4 <- read.table(file="../Day12/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

ipsccat1 <- read.table(file="../ipsc/Cat1",sep="\t",header=T)
ipsccat1a <- subset(ipsccat1,ipsccat1$CAST_perc_y==100)
ipsccat1b <- subset(ipsccat1,ipsccat1$CAST_perc_y==0)
ipsccat2 <- read.table(file="../ipsc/Cat2",sep="\t",header=T)
ipsccat3 <- read.table(file="../ipsc/Cat3",sep="\t",header=T)
ipsccat4 <- read.table(file="../ipsc/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)

listInput <- list(day0cat1a = day0cat1a$ID, 
                  day8cat1a = day8cat1a$ID,
                  day9cat1a = day9cat1a$ID, 
                  day10cat1a = day10cat1a$ID,
                  day12cat1a = day12cat1a$ID,
                  ipsccat1a = ipsccat1a$ID, 
                  day0cat1b = day0cat1b$ID,
                   
                  day8cat1b = day8cat1b$ID,
                  
                  day9cat1b = day9cat1b$ID,
                   
                  day10cat1b = day10cat1b$ID,
                  
                  day12cat1b = day12cat1b$ID,
                  
                  ipsccat1b = ipsccat1b$ID)
                  
                  

names(listInput) <-  c("Day0-Pat","Day8-Pat","Day9-Pat","Day10-Pat","Day12-Pat","iPSCs-Pat","Day0-Mat","Day8-Mat","Day9-Mat","Day10-Mat","Day12-Mat","iPSCs-Mat")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

tiff("04_upsetnarmematvspat.tiff", units="in", width=30, height=25, res=300)

#k <- k[7:24]
upset(fromList(listInput),keep.order = T,sets = c(k),order.by = "freq", sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =5,point.size=7 ,mb.ratio = c(0.55, 0.45),nintersects=120,sets.bar.color=c("chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen"))

dev.off()
rm(list=ls())

###################################################################### imprinted upset after filteration table
library(readxl)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)

imprintedgene <- read.table(file="..7_genelist/imprintedmouse.txt",sep="\t",header=T)
imprintedgene <- subset(imprintedgene,imprintedgene$Status=="Imprinted")

#Reading gene id to gene name file which generated form gtf
gtf2 <- read.table('/media/gayenlab/GayenLab/Mousegenome_102/mouse.txt', header = T, sep = '\t')

GeneSymbol <- gtf2[duplicated(gtf2$GeneSymbol),]

GeneSymbol <- setdiff(gtf2$GeneSymbol,GeneSymbol$GeneSymbol)

GeneSymbol <- data.frame(GeneSymbol)

gtf2 <- merge(gtf2,GeneSymbol,by="GeneSymbol")
dev <- merge(imprintedgene ,gtf2,by.x="Gene",by.y="GeneSymbol")
#Getting coordinated
day0<- read.table(file="..day0/SCALE.output.txt",sep="\t",header=T)
day8<- read.table(file="..day8/SCALE.output.txt",sep="\t",header=T)
day9<- read.table(file="..day9/SCALE.output.txt",sep="\t",header=T)
day10<- read.table(file="..day10/SCALE.output.txt",sep="\t",header=T)
day12<- read.table(file="..day12/SCALE.output.txt",sep="\t",header=T)
ipsc<- read.table(file="..ipsc/SCALE.output.txt",sep="\t",header=T)

commongene <- dev$gene_id

day0imp <- intersect(day0$genename,commongene)
day8imp <- intersect(day8$genename,commongene)
day9imp <- intersect(day9$genename,commongene)
day10imp <- intersect(day10$genename,commongene)
day12imp <- intersect(day12$genename,commongene)
ipscimp <- intersect(ipsc$genename,commongene)
time <- c("day0","day8","day9","day10","day12","ipsc")
imprinted <- c(length(day0imp),length(day8imp),length(day9imp),length(day10imp),length(day12imp),length(ipscimp))
out <- data.frame(time,imprinted)
tiff("05_imprintedafterfilteration.tiff", units="in", width=2, height=2, res=300)
grid.table(out)
dev.off()
rm(list=ls())

############################################## imprint genes in upset plot
library(readxl)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)
#imprint genes are downloaded from Geneimprint database
imprintedgene <- read.table(file="..7_genelist/imprintedmouse.txt",sep="\t",header=T)
imprintedgene <- subset(imprintedgene,imprintedgene$Status=="Imprinted")

#Reading gene id to gene name file which generated form gtf
gtf2 <- read.table('/media/gayenlab/GayenLab/Mousegenome_102/mouse.txt', header = T, sep = '\t')

GeneSymbol <- gtf2[duplicated(gtf2$GeneSymbol),]

GeneSymbol <- setdiff(gtf2$GeneSymbol,GeneSymbol$GeneSymbol)

GeneSymbol <- data.frame(GeneSymbol)

gtf2 <- merge(gtf2,GeneSymbol,by="GeneSymbol")
dev <- merge(imprintedgene ,gtf2,by.x="Gene",by.y="GeneSymbol")
#Getting coordinated
day0<- read.table(file="../Day0/AlleleAratio.txt",sep="\t",header=T)
day8<- read.table(file="../Day8/AlleleAratio.txt",sep="\t",header=T)
day9<- read.table(file="../Day9/AlleleAratio.txt",sep="\t",header=T)
day10<- read.table(file="../Day10/AlleleAratio.txt",sep="\t",header=T)
day12<- read.table(file="../Day12/AlleleAratio.txt",sep="\t",header=T)
ipsc<- read.table(file="../ipsc/AlleleAratio.txt",sep="\t",header=T)

commongene <- dev$gene_id

day0cat1 <- read.table(file="../Day0/Cat1",sep="\t",header=T)
day0cat1<- c(intersect(commongene,day0cat1$ID))
day0cat2 <- read.table(file="../Day0/Cat2",sep="\t",header=T)
day0cat2<- c(intersect(commongene,day0cat2$ID))
day0cat3 <- read.table(file="../Day0/Cat3",sep="\t",header=T)
day0cat3<- c(intersect(commongene,day0cat3$ID))
day0cat4 <- read.table(file="../Day0/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day0cat4<- c(intersect(commongene,day0cat4$ID))

day8cat1 <- read.table(file="../Day8/Cat1",sep="\t",header=T)
day8cat1<- c(intersect(commongene,day8cat1$ID))
day8cat2 <- read.table(file="../Day8/Cat2",sep="\t",header=T)
day8cat2<- c(intersect(commongene,day8cat2$ID))
day8cat3 <- read.table(file="../Day8/Cat3",sep="\t",header=T)
day8cat3<- c(intersect(commongene,day8cat3$ID))
day8cat4 <- read.table(file="../Day8/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day8cat4<- c(intersect(commongene,day8cat4$ID))

day9cat1 <- read.table(file="../Day9/Cat1",sep="\t",header=T)
day9cat1<- c(intersect(commongene,day9cat1$ID))
day9cat2 <- read.table(file="../Day9/Cat2",sep="\t",header=T)
day9cat2<- c(intersect(commongene,day9cat2$ID))
day9cat3 <- read.table(file="../Day9/Cat3",sep="\t",header=T)
day9cat3<- c(intersect(commongene,day9cat3$ID))
day9cat4 <- read.table(file="../Day9/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day9cat4<- c(intersect(commongene,day9cat4$ID))

day10cat1 <- read.table(file="../Day10/Cat1",sep="\t",header=T)
day10cat1<- c(intersect(commongene,day10cat1$ID))
day10cat2 <- read.table(file="../Day10/Cat2",sep="\t",header=T)
day10cat2<- c(intersect(commongene,day10cat2$ID))
day10cat3 <- read.table(file="../Day10/Cat3",sep="\t",header=T)
day10cat3<- c(intersect(commongene,day10cat3$ID))
day10cat4 <- read.table(file="../Day10/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day10cat4<- c(intersect(commongene,day10cat4$ID))


day12cat1 <- read.table(file="../Day12/Cat1",sep="\t",header=T)
day12cat1<- c(intersect(commongene,day12cat1$ID))
day12cat2 <- read.table(file="../Day12/Cat2",sep="\t",header=T)
day12cat2<- c(intersect(commongene,day12cat2$ID))
day12cat3 <- read.table(file="../Day12/Cat3",sep="\t",header=T)
day12cat3<- c(intersect(commongene,day12cat3$ID))
day12cat4 <- read.table(file="../Day12/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
day12cat4<- c(intersect(commongene,day12cat4$ID))

ipsccat1 <- read.table(file="../ipsc/Cat1",sep="\t",header=T)
ipsccat1<- c(intersect(commongene,ipsccat1$ID))
ipsccat2 <- read.table(file="../ipsc/Cat2",sep="\t",header=T)
ipsccat2<- c(intersect(commongene,ipsccat2$ID))
ipsccat3 <- read.table(file="../ipsc/Cat3",sep="\t",header=T)
ipsccat3<- c(intersect(commongene,ipsccat3$ID))
ipsccat4 <- read.table(file="../ipsc/Biallelic_Cat4_above_95_perc.txt",sep="\t",header=T)
ipsccat4<- c(intersect(commongene,ipsccat4$ID))


listInput <- list(day0cat1 = day0cat1, 
                  day8cat1=day8cat1,
                  day9cat1=day9cat1,
                  day10cat1=day10cat1,
                  day12cat1=day12cat1,
                  iPSCcat1=ipsccat1,
                  
                  day0cat2=day0cat2,
                  day8cat2=day8cat2,
                  day9cat2=day9cat2,
                  day10cat2=day10cat2,
                  day12cat2=day12cat2,
                  iPSCcat2=ipsccat2,
                  
                  day0cat3=day0cat3,
                  day8cat3=day8cat3,
                  day9cat3=day9cat3,
                  day10cat3=day10cat3,
                  day12cat3=day12cat3,
                  iPSCcat3=ipsccat3,
                  
                  day0cat4=day0cat4,
                  day8cat4=day8cat4,
                  day9cat4=day9cat4,
                  day10cat4=day10cat4,
                  day12cat4=day12cat4,
                  iPSCcat4=ipsccat4)




tiff("06_imprint", units="in", width=28, height=25, res=300)




names(listInput) <-  c("Day0-cat1","Day8-cat1","Day9-cat1","Day10-cat1","Day12-cat1","iPSCs-cat1","Day0-cat2","Day8-cat2","Day9-cat2","Day10-cat2","Day12-cat2","iPSCs-cat2","Day0-cat3","Day8-cat3","Day9-cat3","Day10-cat3","Day12-cat3","iPSCs-cat3","Day0-cat4","Day8-cat4","Day9-cat4","Day10-cat4","Day12-cat4","iPSCs-cat4")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

#k <- k[7:24]
upset(fromList(listInput),keep.order = T,sets = c(k),order.by = "freq", sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =5,point.size=9 ,mb.ratio = c(0.55, 0.45),sets.bar.color=c("chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","chocolate4","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","red","red","red","red","red","red"),nintersects = 100)


dev.off()
