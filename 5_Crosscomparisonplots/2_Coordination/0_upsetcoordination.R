library(ggplot2)
day0 <- read.table(file="../day0/SCALE.output.txt",header=T)


day0 <- data.frame(table(day0$gene.category),row.names = 1)

day0 <- (day0/(day0[1,]+day0[2,])*100)


day0 <-as.matrix (day0 [1:2,])
colnames(day0) <- c("D")
day0 <- data.frame(day0)
day0$type<- c("Bursty","Non-bursty")
Day="D0"
day0 <- cbind(day0,Day)
day8 <- read.table(file="../day8/SCALE.output.txt",header=T)


day8 <- data.frame(table(day8$gene.category),row.names = 1)

day8 <- (day8/(day8[1,]+day8[2,])*100)

day8 <-as.matrix (day8 [1:2,])
colnames(day8) <- c("D")
day8 <- data.frame(day8)
day8$type <- c("Bursty","Non-bursty")
Day="D8"
day8 <- cbind(day8,Day)
day9 <- read.table(file="../day9/SCALE.output.txt",header=T)


day9 <- data.frame(table(day9$gene.category),row.names = 1)

day9 <- (day9/(day9[1,]+day9[2,])*100)


day9 <-as.matrix (day9 [1:2,])
colnames(day9) <- c("D")
day9 <- data.frame(day9)
day9$type<- c("Bursty","Non-bursty")
Day="D9"
day9 <- cbind(day9,Day)
day10 <- read.table(file="../day10/SCALE.output.txt",header=T)


day10 <- data.frame(table(day10$gene.category),row.names = 1)

day10 <- (day10/(day10[1,]+day10[2,])*100)


day10 <-as.matrix (day10 [1:2,])
colnames(day10) <- c("D")
day10 <- data.frame(day10)
day10$type<- c("Bursty","Non-bursty")
Day="D10"
day10 <- cbind(day10,Day)

day12 <- read.table(file="../day12/SCALE.output.txt",header=T)


day12 <- data.frame(table(day12$gene.category),row.names = 1)

day12 <- (day12/(day12[1,]+day12[2,])*100)


day12 <-as.matrix (day12 [1:2,])
colnames(day12) <- c("D")
day12 <- data.frame(day12)
day12$type<- c("Bursty","Non-bursty")
Day="D12"
day12 <- cbind(day12,Day)

ipsc <- read.table(file="../ipsc/SCALE.output.txt",header=T)
ipsc <- data.frame(table(ipsc$gene.category),row.names = 1)

ipsc <- (ipsc/(ipsc[1,]+ipsc[2,])*100)


ipsc <-as.matrix (ipsc [1:2,])
colnames(ipsc) <- c("D")
ipsc <- data.frame(ipsc)
ipsc$type<- c("Bursty","Non-bursty")
Day="iPSCs"
ipsc <- cbind(ipsc,Day)







tiff("01_Burstyvsnonbursty.tiff", units="in", width=16, height=10, res=300)


datm <- rbind(day0,day8,day9,day10,day12,ipsc)

level_order <- c('D0','D8','D9','D10','D12','iPSCs') 
u <- ggplot(datm, aes(x = factor(Day, level = level_order),y= D, fill = type)) +
geom_bar(position = "stack", stat = "identity")+
  xlab("") + ylab("% of Genes")+
  theme_classic()+
  scale_fill_discrete(name=NULL)+
  geom_col(width = 0.05)+
  
  theme(legend.text=element_text(size=rel(4)))
u  + theme(axis.text = element_text(size = 50)) + theme(axis.title = element_text(size = 50)) + theme(legend.key.size = unit(3,"line"))+ theme(axis.text=element_text(colour="black"))+theme(axis.text.x = element_text(angle = 45, vjust=0.5))
dev.off()
rm(list=ls())
######################################### common genes #########################################################
day0 <- read.table(file="../day0/SCALE.output.txt",sep="\t",header=T)

day0 <- day0[which(day0$gene.category=='Biallelic.bursty'),]

day8 <- read.table(file="../day8/SCALE.output.txt",sep="\t",header=T)

day8 <- day8[which(day8$gene.category=='Biallelic.bursty'),]

day9 <- read.table(file="../day9/SCALE.output.txt",sep="\t",header=T)

day9 <- day9[which(day9$gene.category=='Biallelic.bursty'),]

day10 <- read.table(file="../day10/SCALE.output.txt",sep="\t",header=T)

day10 <- day10[which(day10$gene.category=='Biallelic.bursty'),]

day12 <- read.table(file="../day12/SCALE.output.txt",sep="\t",header=T)

day12 <- day12[which(day12$gene.category=='Biallelic.bursty'),]

ipsc <- read.table(file="../ipsc/SCALE.output.txt",sep="\t",header=T)

ipsc <- ipsc[which(ipsc$gene.category=='Biallelic.bursty'),]

commongene <- Reduce(intersect, list(day0$genename,day8$genename,day9$genename,day10$genename,day12$genename,ipsc$genename))

write.table(commongene,file="../upsetplots/1_Coordinationchanges/02_commongenes.txt",sep="\t")

H <- list(day0=day0$genename,day8=day8$genename,day9=day9$genename,day10=day10$genename,day12=day12$genename,iPSC=ipsc$genename)

library(ggVennDiagram)
library(ggplot2)
tiff("02_commongenesbursty.tiff", units="in", width=16, height=10, res=300)

ggVennDiagram(H, label = "count", edge_size = 2) + scale_fill_distiller(palette = "RdBu")
dev.off()

rm(list=ls())
#########################################   bursty vs non bursty in biallelic genes
library(UpSetR)
#library(clusterProfiler)
#library(org.Mm.eg.db)

#Getting coordinated

day0 <- read.table(file="../day0/SCALE.output.txt",sep="\t",header=T)
day0bursty <- subset(day0,day0$gene.category=="Biallelic.bursty")
day0nonbursty <- subset(day0,day0$gene.category=="Biallelic.nonbursty")

day8 <- read.table(file="../day8/SCALE.output.txt",sep="\t",header=T)
day8bursty <- subset(day8,day8$gene.category=="Biallelic.bursty")
day8nonbursty <- subset(day8,day8$gene.category=="Biallelic.nonbursty")
day9 <- read.table(file="../day9/SCALE.output.txt",sep="\t",header=T)
day9bursty <- subset(day9,day9$gene.category=="Biallelic.bursty")
day9nonbursty <- subset(day9,day9$gene.category=="Biallelic.nonbursty")
day10 <- read.table(file="../day10/SCALE.output.txt",sep="\t",header=T)
day10bursty <- subset(day10,day10$gene.category=="Biallelic.bursty")
day10nonbursty <- subset(day10,day10$gene.category=="Biallelic.nonbursty")
day12 <- read.table(file="../day12/SCALE.output.txt",sep="\t",header=T)
day12bursty <- subset(day12,day12$gene.category=="Biallelic.bursty")
day12nonbursty <- subset(day12,day12$gene.category=="Biallelic.nonbursty")
dayi <- read.table(file="../ipsc/SCALE.output.txt",sep="\t",header=T)
dayibursty <- subset(dayi,dayi$gene.category=="Biallelic.bursty")
dayinonbursty <- subset(dayi,dayi$gene.category=="Biallelic.nonbursty")

test <- do.call("rbind", list(day0,day8,day9,day10,day12,dayi))

dday0 <- append(day0bursty$genename,day0nonbursty$genename)
dday8 <- append(day8bursty$genename,day8nonbursty$genename)
dday9 <- append(day9bursty$genename,day9nonbursty$genename)
dday10 <- append(day10bursty$genename,day10nonbursty$genename)
dday12 <- append(day12bursty$genename,day12nonbursty$genename)
ddayi <- append(dayibursty$genename,dayinonbursty$genename)
commongene <- Reduce(intersect, list(dday0,dday8,dday9,dday10,dday12,ddayi))
day0bursty <- day0bursty[day0bursty$genename %in% commongene ,]
day0nonbursty <- day0nonbursty[day0nonbursty$genename %in% commongene ,]
day8bursty <- day8bursty[day8bursty$genename %in% commongene ,]
day8nonbursty <- day8nonbursty[day8nonbursty$genename %in% commongene ,]

day9bursty <- day9bursty[day9bursty$genename %in% commongene ,]
day9nonbursty <- day9nonbursty[day9nonbursty$genename %in% commongene ,]

day10bursty <- day10bursty[day10bursty$genename %in% commongene ,]
day10nonbursty <- day10nonbursty[day10nonbursty$genename %in% commongene ,]

day12bursty <- day12bursty[day12bursty$genename %in% commongene ,]
day12nonbursty <- day12nonbursty[day12nonbursty$genename %in% commongene ,]

dayibursty <- dayibursty[dayibursty$genename %in% commongene ,]
dayinonbursty <- dayinonbursty[dayinonbursty$genename %in% commongene ,]


listInput <- list(iPSCsnonbursty=dayinonbursty$genename,
                  day12nonbursty=day12nonbursty$genename,
                  day10nonbursty=day10nonbursty$genename,
                  day9nonbursty=day9nonbursty$genename,
                  day8nonbursty=day8nonbursty$genename,
                  day0nonbursty=day0nonbursty$genename,
                  iPSCsbursty = dayibursty$genename,
                  day12bursty = day12bursty$genename, 
                  day10bursty = day10bursty$genename, 
                  day9bursty = day9bursty$genename, 
                  day8bursty = day8bursty$genename, 
                  day0bursty = day0bursty$genename)

tiff("03_upsetnonburstyvsbursty.tiff", units="in", width=25, height=17, res=300)

#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
names(listInput) <-c("iPSCs-Non-bursty" ,"Day12-Non-bursty" ,"Day10-Non-bursty" ,"Day9-Non-bursty","Day8-Non-bursty"  ,"Day0-Non-bursty","iPSCs-bursty" ,"Day12-bursty" ,"Day10-bursty" ,"Day9-bursty" , "Day8-bursty",  "Day0-bursty") 
k <- c("iPSCs-bursty" ,"Day12-bursty" ,"Day10-bursty" ,"Day9-bursty" , "Day8-bursty",  "Day0-bursty","iPSCs-Non-bursty" ,"Day12-Non-bursty" ,"Day10-Non-bursty" ,"Day9-Non-bursty","Day8-Non-bursty"  ,"Day0-Non-bursty")





upset(fromList(listInput),keep.order = T,sets =k,order.by = "freq" , sets.x.label = "Gene set size",matrix.color="dodgerblue4",, 
      text.scale = c(5, 5, 4,1.9, 5, 5),point.size=9,line.size = 2 ,mb.ratio =c(0.55,0.45),nsets=30,sets.bar.color=c("lightcoral","lightcoral","lightcoral","lightcoral","lightcoral","lightcoral","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen"),queries = list(
        list(
          query = intersects,
          params = list("Day12-bursty" ,"Day10-bursty" ,"Day9-bursty" , "Day8-bursty",  "Day0-bursty","iPSCs-Non-bursty" ), 
          color = "red", 
          active = T,
          query.name = "Shift"
        )
      ))


dev.off()
#1allbursty
Allco <- Reduce(intersect, list(dayibursty$genename,day0bursty$genename,
                                day8bursty$genename,
                                day9bursty$genename,
                                day10bursty$genename,day12bursty$genename))
write.table(Allco,file="..4_Ontology/intersectgenes/nonbursty/01_allbursty",sep="\t")

#2ipsc nonbursty
a <-Reduce(intersect, list(dayinonbursty$genename,day0bursty$genename,
                           day8bursty$genename,
                           day9bursty$genename,
                           day10bursty$genename,day12bursty$genename))

write.table(a,file="..4_Ontology/intersectgenes/nonbursty/02_ipscnonburstyabursty",sep="\t")



#3day0nonburstyallbursty.txt"
b<-Reduce(intersect, list(dayibursty$genename,day0nonbursty$genename,
                          day8bursty$genename,
                          day9bursty$genename,
                          day10bursty$genename,day12bursty$genename))
write.table(b,file="..4_Ontology/intersectgenes/nonbursty/03_day0nonburstyallbursty.txt",sep="\t")
#4Allco_40genes
#nonburstyipscday0
m <- Reduce(intersect, list(dayinonbursty$genename,day0nonbursty$genename,
                            day8bursty$genename,
                            day9bursty$genename,
                            day10bursty$genename,day12bursty$genename))
write.table(m,file="..4_Ontology/intersectgenes/nonbursty/04_13_nonburstyinipscaday0.txt",sep="\t")
#5Allnonb
day9ico <-Reduce(intersect, list(dayinonbursty$genename,day0nonbursty$genename,
                                 day8nonbursty$genename,
                                 day9nonbursty$genename,
                                 day10nonbursty$genename,day12nonbursty$genename))
write.table(day9ico,file="..4_Ontology/intersectgenes/nonbursty/05_allnonbursty",sep="\t")
rm(list=ls())

########################################################### upset semivs high co #####################################################
library(UpSetR)
#library(clusterProfiler)
#library(org.Mm.eg.db)

#Getting coordinated

day0 <- read.table(file="../day0/Metadata.txt",sep="\t",header=T)

day8 <- read.table(file="../day8/Metadata.txt",sep="\t",header=T)

day9 <- read.table(file="../day9/Metadata.txt",sep="\t",header=T)

day10 <- read.table(file="../day10/Metadata.txt",sep="\t",header=T)

day12 <- read.table(file="../day12/Metadata.txt",sep="\t",header=T)

dayi <- read.table(file="../ipsc/Metadata.txt",sep="\t",header=T)


test <- do.call("rbind", list(day0,day8,day9,day10,day12,dayi))

commongene <- read.table(file="../upsetplots/1_Coordinationchanges/02_commongenes.txt",header=T)

test = test[test$genename %in% c(commongene$x) ,]

day0co <- subset(test,test$time=="day0"&test$category=="Coordinated")

day8co <- subset(test,test$time=="day8"&test$category=="Coordinated")
day9co <- subset(test,test$time=="day9"&test$category=="Coordinated")
day10co <- subset(test,test$time=="day10"&test$category=="Coordinated")
day12co <- subset(test,test$time=="day12"&test$category=="Coordinated")
dayico <- subset(test,test$time=="ipsc"&test$category=="Coordinated")

day0Semico <- subset(test,test$time=="day0"&test$category=="Semicoordinated")
day8Semico <- subset(test,test$time=="day8"&test$category=="Semicoordinated")
day9Semico <- subset(test,test$time=="day9"&test$category=="Semicoordinated")
day10Semico <- subset(test,test$time=="day10"&test$category=="Semicoordinated")
day12Semico <- subset(test,test$time=="day12"&test$category=="Semicoordinated")
dayiSemico <- subset(test,test$time=="ipsc"&test$category=="Semicoordinated")

day0Ind <- subset(test,test$time=="day0"&test$category=="Independent")
day8Ind <- subset(test,test$time=="day8"&test$category=="Independent")
day9Ind <- subset(test,test$time=="day9"&test$category=="Independent")
day10Ind <- subset(test,test$time=="day10"&test$category=="Independent")
day12Ind <- subset(test,test$time=="day12"&test$category=="Independent")
dayiInd <- subset(test,test$time=="ispc"&test$category=="Independent")

day0Highp2 <- subset(test,test$time=="day0"&test$category=="Highp2")
day8Highp2 <- subset(test,test$time=="day8"&test$category=="Highp2")
day9Highp2 <- subset(test,test$time=="day9"&test$category=="Highp2")
day10Highp2 <- subset(test,test$time=="day10"&test$category=="Highp2")
day12Highp2 <- subset(test,test$time=="day12"&test$category=="Highp2")
dayiHighp2 <- subset(test,test$time=="ipsc"&test$category=="Highp2")
listInput <- list(Day0Coordinated = day0co$genename, 
                  Day8Coordinated=day8co$genename,
                  Day9Coordinated=day9co$genename,
                  Day10Coordinated=day10co$genename,
                  Day12Coordinated=day12co$genename,
                  iPSCsCoordinated=dayico$genename,
                  
                  Day0Semicoordinated=day0Semico$genename,
                  Day8Semicoordinated=day8Semico$genename,
                  Day9Semicoordinated=day9Semico$genename,
                  Day10Semicoordinated=day10Semico$genename,
                  Day12Semicoordinated=day12Semico$genename,
                  iPSCsSemicoordinated=dayiSemico$genename,
                  
                  Day0Independant=day0Ind$genename,
                  Day8Independant=day8Ind$genename,
                  Day9Independant=day9Ind$genename,
                  Day10Independant=day10Ind$genename,
                  Day12Independant=day12Ind$genename,
                  iPSCsIndependant=dayiInd$genename,
                  
                  Day0Lowpohighp2=day0Highp2$genename,
                  Day8Lowpohighp2=day8Highp2$genename,
                  Day9Lowpohighp2=day9Highp2$genename,
                  Day10Lowpohighp2=day10Highp2$genename,
                  Day12Lowpohighp2=day12Highp2$genename,
                  iPSCsLowpohighp2=dayiHighp2$genename)

tiff("04_upsetall.tiff", units="in", width=40, height=25, res=300)

names(listInput) <-  c("Day0-Highly-Coordinted","Day8-Highly-Coordinted","Day9-Highly-Coordinted","Day10-Highly-Coordinted","Day12-Highly-Coordinted","iPSCs-Highly-Coordinted","Day0-Semi-Coordinted","Day8-Semi-Coordinted","Day9-Semi-Coordinted","Day10-Semi-Coordinted","Day12-Semi-Coordinted","iPSCs-Semi-Coordinted","Day0-Independant","Day8-Independant","Day9-Independant","Day10-Independant","Day12-Independant","iPSCs-Independant","Day0-Low p0 high p2","Day8-Low p0 high p2","Day9-Low p0 high p2","Day10-Low p0 high p2","Day12-Low p0 high p2","iPSCs-Low p0 high p2")
#For ordering
#https://github.com/hms-dbmi/UpSetR/issues/169
k <- rev(names(listInput))

k <- k[7:24]
#,show.numbers=F
upset(fromList(listInput),keep.order = T,sets = c(k),order.by = "freq", sets.x.label = "Gene set Size",matrix.color="dodgerblue4",text.scale =6.5,point.size=12 ,mb.ratio = c(0.55, 0.45),sets.bar.color=c("brown","brown","brown","brown","brown","brown","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","lightseagreen","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4","darkslategray4"),queries = list(
  list(
    query = intersects,
    params = list("Day0-Semi-Coordinted","Day8-Semi-Coordinted","Day9-Semi-Coordinted","Day10-Semi-Coordinted","Day12-Semi-Coordinted","iPSCs-Highly-Coordinted" ), 
    color = "red", 
    active = T,
    query.name = "Shift"
  ),list(
    query = intersects,
    params = list("Day0-Highly-Coordinted","Day8-Semi-Coordinted","Day9-Semi-Coordinted","Day10-Semi-Coordinted","Day12-Semi-Coordinted","iPSCs-Semi-Coordinted" ), 
    color = "blue", 
    active = T,
    query.name = "Shift"
  )
))


dev.off()

#1Allsemico
allsemico <- Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                                    day9Semico$genename,
                                    day10Semico$genename,
                                    day12Semico$genename,dayiSemico$genename))

write.table(allsemico,file="..4_Ontology/intersectgenes/01_Allsemico.txt",sep="\t")

#2ipscco -> all semico
b<- Reduce(intersect, list(day0Semico$genename,
                           day8Semico$genename,
                           day9Semico$genename,
                           day10Semico$genename,
                           day12Semico$genename,
                           dayico$genename))
write.table(b,file="..4_Ontology/intersectgenes/02_60ipsco.txt",sep="\t")

#3allsemico-> day0co 40 genes
a <-Reduce(intersect, list(day0co$genename,day8Semico$genename,
                           day9Semico$genename,
                           day10Semico$genename,
                           day12Semico$genename,dayiSemico$genename))

write.table(a,file="..4_Ontology/intersectgenes/03_40day0co.txt",sep="\t")


#4Allco_38genes
Allco <- Reduce(intersect, list(day0co$genename,day8co$genename,
                                                         day9co$genename,
                                                        day10co$genename,
                                                       day12co$genename,dayico$genename))
write.table(Allco,file="..4_Ontology/intersectgenes/04_Allcoordinated.txt",sep="\t")


   #5 day0ind->all semi  13 genes
v <-Reduce(intersect, list(day0Ind$genename,day8Semico$genename,
                           day9Semico$genename,
                           day10Semico$genename,
                           day12Semico$genename,dayiSemico$genename))
write.table(v,file="..4_Ontology/intersectgenes/05_13day0ind_semi.txt",sep="\t")

#6ipscind to all semi 11 genes
m <- Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                            day9Semico$genename,
                            day10Semico$genename,
                            day12Semico$genename,dayiInd$genename))
write.table(m,file="..4_Ontology/intersectgenes/06_11ipscind.txt",sep="\t")

#7 day9ipscco ->all semico 13 genes

day9ico <- Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                            day9co$genename,
                            day10Semico$genename,
                            day12Semico$genename,dayico$genename))
write.table(day9ico,file="..4_Ontology/intersectgenes/07_day9ipsc0allsemi.txt",sep="\t")

#8 day0day9co ->all semico  9 genes

day0day9co <-  Reduce(intersect, list(day0co$genename,day8Semico$genename,
                                      day9co$genename,
                                      day10Semico$genename,
                                      day12Semico$genename,dayiSemico$genename))

write.table(day0day9co,file="..4_Ontology/intersectgenes/08_day0day9c0allsemi.txt",sep="\t")

#9 day9co ->all semi co 12 genes
day9coasemi <-  Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                                      day9co$genename,
                                      day10Semico$genename,
                                      day12Semico$genename,dayiSemico$genename))
write.table(day9coasemi,file="..4_Ontology/intersectgenes/09_day9c0allsemi.txt",sep="\t")


#10 day8semico <- all co 9 genes
day8semiaco <- Reduce(intersect, list(day0co$genename,day8Semico$genename,
                       day9co$genename,
                       day10co$genename,
                       day12co$genename,dayico$genename))

write.table(day8semiaco,file="..4_Ontology/intersectgenes/10_day8semic0allco.txt",sep="\t")

#11 day10co <- all semico 21 genes
day10coasemi <- Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                                          day9Semico$genename,
                                          day10co$genename,
                                          day12Semico$genename,dayiSemico$genename))
write.table(day10coasemi,file="..4_Ontology/intersectgenes/11_day10coasemi.txt",sep="\t")

#12 day0&ipsc to all semico   6 genes
day0ipscoasemi <-  Reduce(intersect, list(day0co$genename,day8Semico$genename,
                                           day9Semico$genename,
                                           day10Semico$genename,
                                           day12Semico$genename,dayico$genename))
write.table(day0ipscoasemi,file="..4_Ontology/intersectgenes/12_day0ipscoasemi.txt",sep="\t")

#13 day8ind all semico 14 genes

day8indallsemi  <- Reduce(intersect, list(day0Semico$genename,day8Ind$genename,
                                          day9Semico$genename,
                                          day10Semico$genename,
                                          day12Semico$genename,dayiSemico$genename))
write.table(day8indallsemi,file="..4_Ontology/intersectgenes/13_day8indallsemi.txt",sep="\t")

#14 day0semi all coordinated 6 genes

day0semiaco <- Reduce(intersect, list(day0Semico$genename,day8co$genename,
                                      day9co$genename,
                                      day10co$genename,
                                      day12co$genename,dayico$genename))
write.table(day0semiaco,file="..4_Ontology/intersectgenes/14_day0semiaco.txt",sep="\t")

#15 day1oico all semi 8genes
day10icoasemi <-Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                                       day9Semico$genename,
                                       day10co$genename,
                                       day12Semico$genename,dayico$genename))

write.table(day10icoasemi,file="..4_Ontology/intersectgenes/15_day10icoasemi.txt",sep="\t")

#16 allsemi upto day12 8 genes
uptoday12asemi <- Reduce(intersect, list(day0Semico$genename,day8Semico$genename,
                                     day9Semico$genename,
                                     day10Semico$genename,
                                     day12Semico$genename))
write.table(uptoday12asemi,file="..4_Ontology/intersectgenes/16_allsemiuptod12.txt",sep="\t")

#17 day010co all semi 7 genes
day010coas<- Reduce(intersect, list(day0co$genename,day8Semico$genename,
                                    day9Semico$genename,
                                    day10co$genename,
                                    day12Semico$genename,dayiSemico$genename))

write.table(day010coas,file="..4_Ontology/intersectgenes/17_day010coas.txt",sep="\t")


rm(list=ls())
############################################ inutero #######################################
library(readxl)
library(tidyr)
library(ggplot2)
library(data.table)
#Reading gene id to gene name file which generated form gtf
gtf2 <- read.table('/media/gayenlab/GayenLab/Mousegenome_102/mouse.txt', header = T, sep = '\t')

GeneSymbol <- gtf2[duplicated(gtf2$GeneSymbol),]

GeneSymbol <- setdiff(gtf2$GeneSymbol,GeneSymbol$GeneSymbol)

GeneSymbol <- data.frame(GeneSymbol)

gtf2 <- merge(gtf2,GeneSymbol,by="GeneSymbol")

rm(GeneSymbol)

#Get input of which genes should be looked

#Got genes from gene ontology of ipsc highco
dev <- read.table(file="..7_genelist/uteroembrypdev",sep="\t")

dev <- merge(dev,gtf2,by.x="V1",by.y="gene_id")


#All metadata category
alldata <- read.table(file="../03_totalmetadatacat",sep="\t",header=T)

colnames(alldata)[1] <- "gene_id"

#Remvoing those duplicated gene names 
alldata <- merge(alldata,gtf2,by="gene_id")

#
t <- alldata[alldata$GeneSymbol %in% dev$GeneSymbol, ]

#get genes in id that to be catergory table
uteroembryo <-unique(t$gene_id)
out <- data.frame(uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,row.names = 1)
colnames(out) <- c("day0","day8","day9","day10","day12","ipsc")
out[] <- "absent"
for(gene in uteroembryo)
{
  print(gene)
  day <-  c("day0","day8","day9","day10","day12","ipsc")
  for (time in day)
  {
    date <- paste("../",time,"/Highlycoordinated.txt",sep="")
    
    dates <- paste("../",time,"/Semicoordinated.txt",sep="")
    
    datei <- paste("../",time,"/Independent.txt",sep="")
    
    dateh <- paste("../",time,"/Highp2.txt",sep="")
    
    dayhighco <- read.table(file=date)
    com <- intersect(gene,row.names(dayhighco))
    n <- as.numeric(length(com))
    daysemco <- read.table(file=dates )
    com <- intersect(gene,row.names(daysemco))
    m <- as.numeric(length(com))
    dayind <- read.table(file= datei )
    com <- intersect(gene,row.names(dayind))
    k <- as.numeric(length(com))
    dayhip <- read.table(file= dateh )
    com <- intersect(gene,row.names(dayhip))
    h <- as.numeric(length(com))
    
    if (n==1) {
      out[gene,time] <- "High-co"
    } else if (m==1) {
      out[gene,time] <- "Semi-co"
    } else if (k==1) {
      out[gene,time] <- "Indep"
    } else if (h==1) {
      out[gene,time] <- "Highp2"
    }else {C <- NULL}
  }
  
}
out$gene_id <- rownames(out)
out <- merge(gtf2,out,by="gene_id")
out <- out[,-(3:5)]

alldata$category1 <- as.numeric(c("Semicoordinated" = "0.75", "Coordinated" = "1","uncategorised" ="NA","Independent" ="0.1" ,"Highp2"="0.5")[alldata$category])

t3 <- alldata[alldata$GeneSymbol %in% t$GeneSymbol, ]
t3$time <- factor(t3$time, levels = c("day0","day8","day9","day10","day12","iPSC"))
wb <- t3[order(t3$time),]

ggplot(data=wb, aes(x=time, y=category1,colour=GeneSymbol,group=GeneSymbol)) + geom_point() +geom_line()+   ylim(0, 1) + ylab("Degree of coordination")+facet_wrap (~GeneSymbol)+
  theme_bw()

rm(list=ls())

#write.table(out,file="inutero")
######################################################## Glucose metabolism ###############33
gtf2 <- read.table('/media/gayenlab/GayenLab/Mousegenome_102/mouse.txt', header = T, sep = '\t')

GeneSymbol <- gtf2[duplicated(gtf2$GeneSymbol),]

GeneSymbol <- setdiff(gtf2$GeneSymbol,GeneSymbol$GeneSymbol)

GeneSymbol <- data.frame(GeneSymbol)

gtf2 <- merge(gtf2,GeneSymbol,by="GeneSymbol")

rm(GeneSymbol)

#Get input of which genes should be looked

#Got genes from gene ontology of ipsc highco
dev <- read.table(file="..7_genelist/05_glucose",sep="\t")

dev <- merge(dev,gtf2,by.x="V1",by.y="gene_id")


#All metadata category
alldata <- read.table(file="../03_totalmetadatacat",sep="\t",header=T)

colnames(alldata)[1] <- "gene_id"

#Remvoing those duplicated gene names 
alldata <- merge(alldata,gtf2,by="gene_id")

#
t <- alldata[alldata$GeneSymbol %in% dev$GeneSymbol, ]

#get genes in id that to be catergory table
uteroembryo <-unique(t$gene_id)
out <- data.frame(uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,uteroembryo,row.names = 1)
colnames(out) <- c("day0","day8","day9","day10","day12","ipsc")
out[] <- "absent"
for(gene in uteroembryo)
{
  print(gene)
  day <-  c("day0","day8","day9","day10","day12","ipsc")
  for (time in day)
  {
    date <- paste("../",time,"/Highlycoordinated.txt",sep="")
    
    dates <- paste("../",time,"/Semicoordinated.txt",sep="")
    
    datei <- paste("../",time,"/Independent.txt",sep="")
    
    dateh <- paste("../",time,"/Highp2.txt",sep="")
    
    dayhighco <- read.table(file=date)
    com <- intersect(gene,row.names(dayhighco))
    n <- as.numeric(length(com))
    daysemco <- read.table(file=dates )
    com <- intersect(gene,row.names(daysemco))
    m <- as.numeric(length(com))
    dayind <- read.table(file= datei )
    com <- intersect(gene,row.names(dayind))
    k <- as.numeric(length(com))
    dayhip <- read.table(file= dateh )
    com <- intersect(gene,row.names(dayhip))
    h <- as.numeric(length(com))
    
    if (n==1) {
      out[gene,time] <- "High-co"
    } else if (m==1) {
      out[gene,time] <- "Semi-co"
    } else if (k==1) {
      out[gene,time] <- "Indep"
    } else if (h==1) {
      out[gene,time] <- "Highp2"
    }else {C <- NULL}
  }
  
}
out$gene_id <- rownames(out)
out <- merge(gtf2,out,by="gene_id")
out <- out[,-(3:5)]

alldata$category1 <- as.numeric(c("Semicoordinated" = "0.75", "Coordinated" = "1","uncategorised" ="NA","Independent" ="0.1" ,"Highp2"="0.5")[alldata$category])

t3 <- alldata[alldata$GeneSymbol %in% t$GeneSymbol, ]
t3$time <- factor(t3$time, levels = c("day0","day8","day9","day10","day12","iPSC"))
wb <- t3[order(t3$time),]

ggplot(data=wb, aes(x=time, y=category1,colour=GeneSymbol,group=GeneSymbol)) + geom_point() +geom_line()+   ylim(0, 1) + ylab("Degree of coordination")+facet_wrap (~GeneSymbol)+
  theme_bw()

#write.table(out,file="../glucosemetabolismreactome")
