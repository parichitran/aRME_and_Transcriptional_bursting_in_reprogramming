##################################################################### Cat 2##################################################################################################
library(ggplot2)
day0 <- read.table(file="..day0/SCALE.output.txt",header=T)
cat2<- read.csv(file="..Category/Day0/Cat2",sep = "\t")

day0 <- day0[day0$genename %in% cat2$ID , ]

day0 <- data.frame(table(day0$gene.category),row.names = 1)

day0 <- (day0/(day0[1,]+day0[2,])*100)


day0 <-as.matrix (day0 [1:2,])
colnames(day0) <- c("D")
day0 <- data.frame(day0)
day0$type<- c("Bursty","Non-bursty")
Day="D0"
day0 <- cbind(day0,Day)
day8 <- read.table(file="..day8/SCALE.output.txt",header=T)
cat28<- read.csv(file="..Category/Day8/Cat2",sep = "\t")

day8 <- day8[day8$genename %in% cat28$ID , ]

day8 <- data.frame(table(day8$gene.category),row.names = 1)

day8 <- (day8/(day8[1,]+day8[2,])*100)

day8 <-as.matrix (day8 [1:2,])
colnames(day8) <- c("D")
day8 <- data.frame(day8)
day8$type <- c("Bursty","Non-bursty")
Day="D8"
day8 <- cbind(day8,Day)
day9 <- read.table(file="..day9/SCALE.output.txt",header=T)
cat29<- read.csv(file="..Category/Day9/Cat2",sep = "\t")

day9 <- day9[day9$genename %in% cat29$ID , ]

day9 <- data.frame(table(day9$gene.category),row.names = 1)

day9 <- (day9/(day9[1,]+day9[2,])*100)


day9 <-as.matrix (day9 [1:2,])
colnames(day9) <- c("D")
day9 <- data.frame(day9)
day9$type<- c("Bursty","Non-bursty")
Day="D9"
day9 <- cbind(day9,Day)
day10 <- read.table(file="..day10/SCALE.output.txt",header=T)
cat210<- read.csv(file="..Category/Day10/Cat2",sep = "\t")

day10 <- day10[day10$genename %in% cat210$ID , ]

day10 <- data.frame(table(day10$gene.category),row.names = 1)

day10 <- (day10/(day10[1,]+day10[2,])*100)


day10 <-as.matrix (day10 [1:2,])
colnames(day10) <- c("D")
day10 <- data.frame(day10)
day10$type<- c("Bursty","Non-bursty")
Day="D10"
day10 <- cbind(day10,Day)

day12 <- read.table(file="..day12/SCALE.output.txt",header=T)

cat212<- read.csv(file="..Category/Day12/Cat2",sep = "\t")

day12 <- day12[day12$genename %in% cat212$ID , ]
day12 <- data.frame(table(day12$gene.category),row.names = 1)

day12 <- (day12/(day12[1,]+day12[2,])*100)


day12 <-as.matrix (day12 [1:2,])
colnames(day12) <- c("D")
day12 <- data.frame(day12)
day12$type<- c("Bursty","Non-bursty")
Day="D12"
day12 <- cbind(day12,Day)

ipsc <- read.table(file="..ipsc/SCALE.output.txt",header=T)

cat2i<- read.csv(file="..Category/ipsc/Cat2",sep = "\t")

ipsc <- ipsc[ipsc$genename %in% cat2i$ID , ]


ipsc <- data.frame(table(ipsc$gene.category),row.names = 1)

ipsc <- (ipsc/(ipsc[1,]+ipsc[2,])*100)


ipsc <-as.matrix (ipsc [1:2,])
colnames(ipsc) <- c("D")
ipsc <- data.frame(ipsc)
ipsc$type<- c("Bursty","Non-bursty")
Day="iPSCs"
ipsc <- cbind(ipsc,Day)







tiff("0_Burstycat2.tiff", units="in", width=16, height=10, res=300)


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


rm(list=ls())

################################################# cat3 ###################################################################################################################################
day0 <- read.table(file="..day0/SCALE.output.txt",header=T)
cat2<- read.csv(file="..Category/Day0/Cat3",sep = "\t")

day0 <- day0[day0$genename %in% cat2$ID , ]

day0 <- data.frame(table(day0$gene.category),row.names = 1)

day0 <- (day0/(day0[1,]+day0[2,])*100)


day0 <-as.matrix (day0 [1:2,])
colnames(day0) <- c("D")
day0 <- data.frame(day0)
day0$type<- c("Bursty","Non-bursty")
Day="D0"
day0 <- cbind(day0,Day)
day8 <- read.table(file="..day8/SCALE.output.txt",header=T)
cat28<- read.csv(file="..Category/Day8/Cat3",sep = "\t")

day8 <- day8[day8$genename %in% cat28$ID , ]

day8 <- data.frame(table(day8$gene.category),row.names = 1)

day8 <- (day8/(day8[1,]+day8[2,])*100)

day8 <-as.matrix (day8 [1:2,])
colnames(day8) <- c("D")
day8 <- data.frame(day8)
day8$type <- c("Bursty","Non-bursty")
Day="D8"
day8 <- cbind(day8,Day)
day9 <- read.table(file="..day9/SCALE.output.txt",header=T)
cat29<- read.csv(file="..Category/Day9/Cat3",sep = "\t")

day9 <- day9[day9$genename %in% cat29$ID , ]

day9 <- data.frame(table(day9$gene.category),row.names = 1)

day9 <- (day9/(day9[1,]+day9[2,])*100)


day9 <-as.matrix (day9 [1:2,])
colnames(day9) <- c("D")
day9 <- data.frame(day9)
day9$type<- c("Bursty","Non-bursty")
Day="D9"
day9 <- cbind(day9,Day)
day10 <- read.table(file="..day10/SCALE.output.txt",header=T)
cat210<- read.csv(file="..Category/Day10/Cat3",sep = "\t")

day10 <- day10[day10$genename %in% cat210$ID , ]

day10 <- data.frame(table(day10$gene.category),row.names = 1)

day10 <- (day10/(day10[1,]+day10[2,])*100)


day10 <-as.matrix (day10 [1:2,])
colnames(day10) <- c("D")
day10 <- data.frame(day10)
day10$type<- c("Bursty","Non-bursty")
Day="D10"
day10 <- cbind(day10,Day)

day12 <- read.table(file="..day12/SCALE.output.txt",header=T)

cat212<- read.csv(file="..Category/Day12/Cat3",sep = "\t")

day12 <- day12[day12$genename %in% cat212$ID , ]
day12 <- data.frame(table(day12$gene.category),row.names = 1)

day12 <- (day12/(day12[1,]+day12[2,])*100)


day12 <-as.matrix (day12 [1:2,])
colnames(day12) <- c("D")
day12 <- data.frame(day12)
day12$type<- c("Bursty","Non-bursty")
Day="D12"
day12 <- cbind(day12,Day)

ipsc <- read.table(file="..ipsc/SCALE.output.txt",header=T)

cat2i<- read.csv(file="..Category/ipsc/Cat3",sep = "\t")

ipsc <- ipsc[ipsc$genename %in% cat2i$ID , ]


ipsc <- data.frame(table(ipsc$gene.category),row.names = 1)

ipsc <- (ipsc/(ipsc[1,]+ipsc[2,])*100)


ipsc <-as.matrix (ipsc [1:2,])
colnames(ipsc) <- c("D")
ipsc <- data.frame(ipsc)
ipsc$type<- c("Bursty","Non-bursty")
Day="iPSCs"
ipsc <- cbind(ipsc,Day)







tiff("1_Burstycat3.tiff", units="in", width=16, height=10, res=300)


datm <- rbind(day0,day8,day9,day10,day12,ipsc)

level_order <- c('D0','D8','D9','D10','D12','iPSCs') 
u <- ggplot(datm, aes(x = factor(Day, level = level_order),y= D, fill = type)) +
  geom_bar(position = "stack", stat = "identity")+
  xlab("") + ylab("% of Genes")+
  theme_classic()+
  scale_fill_discrete(name=NULL)+
  geom_col(width = 0.05)+
  
  theme(legend.text=element_text(size=rel(4)))
u  + theme(axis.text = element_text(size = 50)) + theme(axis.title = element_text(size = 50)) + theme(legend.key.size = unit(3,"line"))+ theme(axis.text=element_text(colour="black"))+theme(axis.text.x = element_text(angle = 45, vjust=0.5))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
rm(list=ls())

########################################### cat 4 #######################################3333
day0 <- read.table(file="..day0/SCALE.output.txt",header=T)
cat2<- read.csv(file="..Category/Day0/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

day0 <- day0[day0$genename %in% cat2$ID , ]

day0 <- data.frame(table(day0$gene.category),row.names = 1)

day0 <- (day0/(day0[1,]+day0[2,])*100)


day0 <-as.matrix (day0 [1:2,])
colnames(day0) <- c("D")
day0 <- data.frame(day0)
day0$type<- c("Bursty","Non-bursty")
Day="D0"
day0 <- cbind(day0,Day)
day8 <- read.table(file="..day8/SCALE.output.txt",header=T)
cat28<- read.csv(file="..Category/Day8/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

day8 <- day8[day8$genename %in% cat28$ID , ]

day8 <- data.frame(table(day8$gene.category),row.names = 1)

day8 <- (day8/(day8[1,]+day8[2,])*100)

day8 <-as.matrix (day8 [1:2,])
colnames(day8) <- c("D")
day8 <- data.frame(day8)
day8$type <- c("Bursty","Non-bursty")
Day="D8"
day8 <- cbind(day8,Day)
day9 <- read.table(file="..day9/SCALE.output.txt",header=T)
cat29<- read.csv(file="..Category/Day9/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

day9 <- day9[day9$genename %in% cat29$ID , ]

day9 <- data.frame(table(day9$gene.category),row.names = 1)

day9 <- (day9/(day9[1,]+day9[2,])*100)


day9 <-as.matrix (day9 [1:2,])
colnames(day9) <- c("D")
day9 <- data.frame(day9)
day9$type<- c("Bursty","Non-bursty")
Day="D9"
day9 <- cbind(day9,Day)
day10 <- read.table(file="..day10/SCALE.output.txt",header=T)
cat210<- read.csv(file="..Category/Day10/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

day10 <- day10[day10$genename %in% cat210$ID , ]

day10 <- data.frame(table(day10$gene.category),row.names = 1)

day10 <- (day10/(day10[1,]+day10[2,])*100)


day10 <-as.matrix (day10 [1:2,])
colnames(day10) <- c("D")
day10 <- data.frame(day10)
day10$type<- c("Bursty","Non-bursty")
Day="D10"
day10 <- cbind(day10,Day)

day12 <- read.table(file="..day12/SCALE.output.txt",header=T)

cat212<- read.csv(file="..Category/Day12/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

day12 <- day12[day12$genename %in% cat212$ID , ]
day12 <- data.frame(table(day12$gene.category),row.names = 1)

day12 <- (day12/(day12[1,]+day12[2,])*100)


day12 <-as.matrix (day12 [1:2,])
colnames(day12) <- c("D")
day12 <- data.frame(day12)
day12$type<- c("Bursty","Non-bursty")
Day="D12"
day12 <- cbind(day12,Day)

ipsc <- read.table(file="..ipsc/SCALE.output.txt",header=T)

cat2i<- read.csv(file="..Category/ipsc/Biallelic_Cat4_above_95_perc.txt",sep = "\t")

ipsc <- ipsc[ipsc$genename %in% cat2i$ID , ]


ipsc <- data.frame(table(ipsc$gene.category),row.names = 1)

ipsc <- (ipsc/(ipsc[1,]+ipsc[2,])*100)


ipsc <-as.matrix (ipsc [1:2,])
colnames(ipsc) <- c("D")
ipsc <- data.frame(ipsc)
ipsc$type<- c("Bursty","Non-bursty")
Day="iPSCs"
ipsc <- cbind(ipsc,Day)







tiff("02_Burstycat4.tiff", units="in", width=16, height=10, res=300)


datm <- rbind(day0,day8,day9,day10,day12,ipsc)

level_order <- c('D0','D8','D9','D10','D12','iPSCs') 
u <- ggplot(datm, aes(x = factor(Day, level = level_order),y= D, fill = type)) +
  geom_bar(position = "stack", stat = "identity")+
  xlab("") + ylab("% of Genes")+
  theme_classic()+
  scale_fill_discrete(name=NULL)+
  geom_col(width = 0.05)+
  
  theme(legend.text=element_text(size=rel(4)))
u  + theme(axis.text = element_text(size = 50)) + theme(axis.title = element_text(size = 50)) + theme(legend.key.size = unit(3,"line"))+ theme(axis.text=element_text(colour="black"))+theme(axis.text.x = element_text(angle = 45, vjust=0.5))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
