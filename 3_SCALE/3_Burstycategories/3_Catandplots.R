########################################### preliminary Category and plots from scale output#######################################
library("PerformanceAnalytics")
library(dplyr)
df=read.table("SCALE.output.txt", sep='\t', header=TRUE)
df = df[df$gene.category=="Biallelic.bursty",]
pdf('plot_for_filter.pdf',width=6,height=6)
smoothScatter(df$prop_Off,df$prop_ab,main = "Day0",xlab='P0 (% of cells expressing neither alleles)',ylab='P2 (% of cells expressing both alleles )',nrpoints=Inf,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.5)

#Coordinated line
points(c(0,1),c(1,0),type='l',col='deepskyblue4',lty=2,lwd=2)
p=0:1000/1000

#Independant lines
points((1-p)^2,p^2,type='l',col='red',lty=2,lwd=2)
points(((1-p)^2)+0.05,(p^2)+0.05,type='l',col='red',lty=2,lwd=2)
points(((1-p)^2)-0.05,(p^2)-0.05,type='l',col='red',lty=2,lwd=2)
legend(0.30, 1, c("Coordination", "Independence with shared kinetics"), 
       col = c("deepskyblue4", "red"), bg = "#FFFFFFAA", 
       lty = 2, lwd = 4, box.col = "#FFFFFF00")

#Getting genes near the curved lines(Independent)
library(sp)
p=0:1000/1000
middle <- cbind((1 - p)^2, p^2)
lower <- cbind(rev((1 - p)^2 - 0.05), rev((p^2) - 0.05))
uper <- cbind(rev((1 - p)^2 + 0.05), rev((p^2) + 0.05))
poly <- rbind(middle, uper)
polygon(poly)
inside <- point.in.polygon(df$prop_Off,df$prop_ab, poly[, 1], poly[, 2])
data1=df[inside==1, ]


library(sp)
p=0:1000/1000
middle <- cbind((1 - p)^2, p^2)
lower <- cbind(rev((1 - p)^2 - 0.05), rev((p^2) - 0.05))
uper <- cbind(rev((1 - p)^2 + 0.05), rev((p^2) + 0.05))
poly <- rbind(middle, lower)
polygon(poly)
inside <- point.in.polygon(df$prop_Off,df$prop_ab, poly[, 1], poly[, 2])
data2=df[inside==1, ]
dev.off()

data <- rbind(data1, data2)

############################Independent genes
write.csv(data, "Curve_line_genes")
##############################################
#Indpendent line plot
pdf('Gene_near_the_independent_curve.pdf',width=6,height=6)
smoothScatter(data$prop_Off,data$prop_ab,main = "Day0",xlab='P0 (% of cells expressing neither alleles)',ylab='P2 (% of cells expressing both alleles )',nrpoints=Inf,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.5)
points(c(0,1),c(1,0),type='l',col='deepskyblue4',lty=2,lwd=2)
p=0:1000/1000
points((1-p)^2,p^2,type='l',col='red',lty=2,lwd=2)
#points(((1-p)^2)+0.03,(p^2)+0.03,type='l',col='red',lty=2,lwd=2)
#points(((1-p)^2)-0.03,(p^2)-0.03,type='l',col='red',lty=2,lwd=2)
legend(0.30, 1, c("Coordination", "Independence with shared kinetics"), 
       col = c("deepskyblue4", "red"), bg = "#FFFFFFAA", 
       lty = 2, lwd = 4, box.col = "#FFFFFF00")

dev.off()

#Diagonal line cordinated genes plot and genes
df$sum=df$prop_Off+df$prop_ab
df2=df[df$sum >= 0.90,]

######################################
write.csv(df2, "Diagonal_line_genes")
######################################

pdf('Gene_near_the_cordinated_line.pdf',width=6,height=6)
smoothScatter(df2$prop_Off,df2$prop_ab,main = "Day0",xlab='P0 (% of cells expressing neither alleles)',ylab='P2 (% of cells expressing both alleles )',nrpoints=Inf,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.5)
points(c(0,1),c(1,0),type='l',col='deepskyblue4',lty=2,lwd=2)
p=0:1000/1000
points((1-p)^2,p^2,type='l',col='red',lty=2,lwd=2)
#points(((1-p)^2)+0.03,(p^2)+0.03,type='l',col='red',lty=2,lwd=2)
#points(((1-p)^2)-0.03,(p^2)-0.03,type='l',col='red',lty=2,lwd=2)
legend(0.30, 1, c("Coordination", "Independence with shared kinetics"), 
       col = c("deepskyblue4", "red"), bg = "#FFFFFFAA", 
       lty = 2, lwd = 4, box.col = "#FFFFFF00")

dev.off()
rm(list = ls())

############################################### uncategorised genes 

library("PerformanceAnalytics")
library(dplyr)
df=read.table("SCALE.output.txt", sep='\t', header=TRUE)
EPI = df[df$gene.category=="Biallelic.bursty",]
EPI = data.frame(EPI)
#pdf('plot_for_filter.pdf',width=6,height=6)
smoothScatter(EPI$prop_Off,EPI$prop_ab,main = "EPI E6.5",xlab='P0 (% of cells expressing neither alleles)',ylab='P2 (% of cells expressing both alleles )',nrpoints=Inf,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.5)
points(c(0,1),c(1,0),type='l',col='deepskyblue4',lty=2,lwd=2)
p=0:1000/1000
points((1-p)^2,p^2,type='l',col='red',lty=2,lwd=2)
points(((1-p)^2)+0.05,(p^2)+0.05,type='l',col='red',lty=2,lwd=2)
points(((1-p)^2)-0.05,(p^2)-0.05,type='l',col='red',lty=2,lwd=2)
legend(0.30, 1, c("Coordination", "Independence with shared kinetics"), 
       col = c("deepskyblue4", "red"), bg = "#FFFFFFAA", 
       lty = 2, lwd = 4, box.col = "#FFFFFF00")




#Below function identify will create an interactive session to select the uncats
pts=identify(EPI$prop_Off, EPI$prop_ab,tolerance = 0.5)
#Select manually
v <- pts
z <- c("1427")
pts <- append(pt,z)
pts <- as.integer(pts)
no_cat=EPI[pts,]
#write.csv(no_cat, "/media/gayenlab/GayenLab10TB1/01_Project2/03_SCALE/0_Filtered/day0/Uncat.csv")
rm(list = ls())
################################################ Final category ################################################
#Getting scale output after filteration
day <- read.table(file="SCALE.output.txt",sep="\t",header=T)

#
day <- day[which(day$gene.category=='Biallelic.bursty'),]

#Genes in dependant
dayd <- read.table(file="Diagonal_line_genes",sep=",",header=T)

#Genes indpependant
dayin <- read.table(file="Curve_line_genes",sep=",",header=T)

#Genes at high p2
commongenes <- intersect(dayd$genename,dayin$genename)


#Genes in uncategory
unc <- read.csv("Uncat.csv")
#Genes in semicoordinated
genetoremove <- unique(append(dayd$genename,dayin$genename))
b <- append(genetoremove,unc$genename)

semico <- setdiff(day$genename,b)


day <- data.frame(day,row.names = 1)
semico <- day[semico,]
highp2 <- day[commongenes,]

#High c0ordinated genes alone
daydf <- setdiff(dayd$genename,dayin$genename)

daydf <- day[daydf,]

#Independent
dayinf <- setdiff(dayin$genename,dayd$genename)

dayinf <- day[dayinf,]

#Write output
write.table(daydf,file="Highlycoordinated.txt",sep="\t")

write.table(dayinf,file="Independent.txt",sep="\t")

write.table(semico,file="Semicoordinated.txt",sep="\t")

write.table(highp2,file="Highp2.txt",sep = "\t")


time <- "day0"


category <- "Coordinated"
dayA <- cbind(row.names(daydf),time,category)

category <- "Independent"

dayB <- cbind(row.names(dayinf),time,category)

tempA <- rbind(dayA,dayB)

category <- "Semicoordinated"

dayC <- cbind(row.names(semico),time,category)

category <- "Highp2"

dayD <- cbind(row.names(highp2),time,category)

tempB <- rbind(dayC,dayD)

final <- rbind(tempA,tempB)

category <- "uncategorised"
dayunc <- cbind(unc$genename,time,category)
dayunc <- data.frame(dayunc)


final <- rbind(final,dayunc)

colnames(final)[1] <- "genename"



write.table(final,file="Metadata.txt",sep = "\t")
rm(list = ls())

################################################# Set for plot ########################################3333
day <- read.table(file="SCALE.output.txt",sep="\t",header=T)

#
day <- day[which(day$gene.category=='Biallelic.bursty'),]

day <- data.frame(day,row.names = 1)
meta <- read.table(file="Metadata.txt",sep="\t",header=T)

meta <- data.frame(meta,row.names = 1)

m <- merge(day,meta,by='row.names',all=T)

t <- m[,c("Row.names","prop_ab","prop_Off","category")]

write.table(t,file="catergoryplotfile.txt",sep="\t")


rm(list = ls())

################################################# category plots final ##########################################

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra) 

#EPI
df=read.table("catergoryplotfile.txt", header=T, sep="\t")
tiff("categoryplot.tiff", units="in", width=6, height=5, res=300)

curvy <- function(x) { ((1-x)^2)^x^2 }
straight <- function(x) 0.90-x
p=0:100/100
color <- c(Highp2="olivedrab2", Coordinated="lightblue4", Semicoordinated="lightseagreen", Independent="salmon3", uncategorised="black")
straight <- function(x) 0.90-x
q=ggplot(df, aes(prop_Off, prop_ab , col=category)) +geom_point(aes(shape=category),size = 1.5)+
  scale_shape_manual(values=c(8, 15, 17, 20, 22))+
  geom_line(data = data.frame(A = (1-p)^2, B = (p^2)),mapping = aes(A, B),color = "red",linetype = "dashed") +
  geom_line(data = data.frame(A = 0:1, B = 1:0),mapping = aes(A, B),color = "deepskyblue4",linetype = "dashed")+
  geom_line(data = data.frame(A = (1-p)^2-0.05, B = (p^2)-0.05), mapping = aes(A, B), color = "red", linetype = "dashed", size=0.8)+
  geom_line(data = data.frame(A = (1-p)^2+0.05, B = (p^2)+0.05), mapping = aes(A, B), color = "red", linetype = "dashed", size=0.8)+
  geom_line(stat='function', fun=straight, color='deepskyblue4', linetype='dashed', size=0.8)+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  scale_color_manual(values = color)+theme_cowplot()+theme(panel.border = element_rect(fill=NA,color="black", size=0.60, linetype="solid"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title=element_text(size=17))+
  scale_y_continuous(breaks=seq(0.0 , 1.0, 0.20),limits=c(0,1))+
  scale_x_continuous(breaks=seq(0.0, 1.0, 0.20),limits=c(0,1))+xlab('') + ylab('')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
q+theme(legend.position = "none")
#ggsave(plot = q, width = 8, height = 6, filename = "scatterplot.pdf")

dev.off()

tiff("legend.tiff", units="in", width=7, height=5, res=300)

legend <- cowplot::get_legend(q)

grid.newpage()
grid.draw(legend)
dev.off()

rm(list=ls())

################################################ freq and size #############################################
library(qvalue)
library("ggplot2")
day <- read.table(file="SCALE.output.txt",header=T)
day <- day[which(day$gene.category=="Biallelic.bursty"),]

#get only with true values
days <- subset(day,day$pval.size!="-")
day <- subset(day,day$pval.kon!="-")

days <- transform(days, pval.size = as.numeric(pval.size))
day <- transform(day, pval.kon = as.numeric(pval.kon))
#padjusted values for kon
kon <- qvalue(as.numeric(day$pval.kon)  )
day$padjust.on<- kon$qvalues

#Assigning significance
groupon <- as.factor(ifelse(day$padjust.on<=0.05, "Sig", "Nonsig"))
day$sig <- groupon

#Simlarly done for ksize
ksiz <- qvalue(as.numeric(days$pval.size)  ) 
days$padjust.size<-ksiz$qvalues

groups <- as.factor(ifelse(days$padjust.size<=0.05, "Sig", "Nonsig"))
days$sig <- groups


#Burst kinetics correlation
cor1.1=cor(log(as.numeric(day$konA)),log(as.numeric(day$konB)))

#Burst size correlation
g <- as.numeric(days$sA)
j <- as.numeric(days$koffA)
gb <- as.numeric(days$sB)
jb<- as.numeric(days$koffB)
cor1.2=cor(log(g/j),log(gb/jb))

tiff("freq.tiff", units="in", width=5.5, height=5, res=300)
my_ggp <- ggplot(day,aes(x=log(as.numeric(day$konB)),y=log(as.numeric(day$konA)),group=sig))+
  
  geom_point(aes(color=sig,shape=sig),size = 2.5)+
  xlab("") + ylab("")+
  coord_cartesian(xlim =c(min(log(as.numeric(day$konA)),log(as.numeric(day$konB))),max(log(as.numeric(day$konA)),log(as.numeric(day$konB)))), ylim = c(min(log(as.numeric(day$konA)),log(as.numeric(day$konB))),max(log(as.numeric(day$konA)),log(as.numeric(day$konB)))))+
   geom_abline( linetype = "dotted")+
  theme_bw()+
 
  scale_color_manual(values=c('grey','red'))
my_ggp + theme(axis.text = element_text(size = 30)) + theme(legend.position = "none")
dev.off()

tiff("size.tiff", units="in", width=5, height=5, res=300)

my_ggps <- ggplot(days,aes(x=log(gb/jb),y=log(g/j),group=sig))+
  
  geom_point(aes(color=sig,shape=sig),size = 2.5)+
  xlab("") + ylab("")+
  coord_cartesian(xlim =c(min(log(g/j),log(gb/jb)),max(log(g/j),log(gb/jb))), ylim = c(min(log(g/j),log(gb/jb)),max(log(g/j),log(gb/jb))))+
  geom_abline( linetype = "dotted")+
  theme_bw()+
  
  scale_color_manual(values=c('grey','red'))
my_ggps + theme(axis.text = element_text(size = 30)) + theme(legend.position = "none")
dev.off()

cor <- append(cor1.1,cor1.2)
write(cor,file="correlation")
rm(list=ls())
