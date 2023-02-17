#ploting the piechart
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra)    


# get cat data.
cat1<- read.csv(file="Cat1",sep = "\t")
cat2<- read.csv(file="Cat2",sep = "\t")
cat3<- read.csv(file="Cat3",sep = "\t")
cat4<- read.csv(file="Biallelic_Cat4_above_95_perc.txt",sep = "\t")
cat1genes<- as.numeric(length(row.names(cat1)))
cat2genes<- as.numeric(length(row.names(cat2)))
cat3genes<- as.numeric(length(row.names(cat3)))
cat4genes<- as.numeric(length(row.names(cat4)))



total <- cat1genes+cat2genes+cat3genes+cat4genes

print(total)
cat1percent <- cat1genes/total*100
cat2percent <- cat2genes/total*100
cat3percent <- cat3genes/total*100
cat4percent <- cat4genes/total*100

#Making dataframe for  plotiing
count.data<- data.frame(
  prop = c(2, round(cat2percent) , round(cat3percent), 8),
  class = c("Cat1","Cat2","Cat3","Cat4"),
  n = c(as.numeric(length(row.names(cat1))),as.numeric(length(row.names(cat2))), as.numeric(length(row.names(cat3))), as.numeric(length(row.names(cat4)))))


#https://www.datanovia.com/en/blog/how-to-create-a-pie-chart-in-r-using-ggplot2/

#Making label position
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data

#Colours input
mycols <- c("#0073C2FF", "#EFC000FF", "snow4", "#CD534CFF")

#donut piechart x=2 and xlim creates the hole in the piechart
t <- ggplot(count.data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  xlim(0.5, 2.5)

count.data$prop <- as.character(paste0(count.data$prop,"%"))
tiff("genwiseplot.tiff", units="in", width=10, height=10, res=300)


t+ geom_label_repel(aes(y = lab.ypos, label = count.data$prop), color = "white",size=17,fontface='bold')+
  guides(fill = guide_legend(title = "Group"))+
  theme_void()+ theme(legend.position="none")

dev.off()

legend<- cowplot::get_legend(t)
tiff("legendwhite.tiff", units="in", width=2, height=2, res=300)

grid.draw(legend)
dev.off()

rm(list=ls())
############################################### Mat vs pat similar to above plots
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra) 
# Create test data.
cat1<- read.csv(file="..Day0/Cat1",sep = "\t")
cat2<- read.csv(file="..Day0/Cat2",sep = "\t")
cat3<- read.csv(file="..Day0/Cat3",sep = "\t")
cat4<- read.csv(file="..Day0/Biallelic_Cat4_above_95_perc.txt",sep = "\t")
cat1a <- subset(cat1,cat1$CAST_perc_y==100)
cat1b <- subset(cat1,cat1$CAST_perc_y==0)

cat1agenes<- as.numeric(length(row.names(cat1a)))
cat1bgenes<- as.numeric(length(row.names(cat1b)))

cat2genes<- as.numeric(length(row.names(cat2)))
cat3genes<- as.numeric(length(row.names(cat3)))
cat4genes<- as.numeric(length(row.names(cat4)))



total <- cat1agenes+cat1bgenes+cat2genes+cat3genes+cat4genes

print(total)
cat1apercent <- cat1agenes/total*100
cat1bpercent <- cat1bgenes/total*100

cat2percent <- cat2genes/total*100
cat3percent <- cat3genes/total*100
cat4percent <- cat4genes/total*100

count.data<- data.frame(
  prop = c(2, 0.4,12 , 78, 7.6),
  class = c("Cat1P","Cat1M","Cat2","Cat3","Cat4"),
  n = c(as.numeric(length(row.names(cat1a))),as.numeric(length(row.names(cat1b))),as.numeric(length(row.names(cat2))), as.numeric(length(row.names(cat3))), as.numeric(length(row.names(cat4)))))



count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data


mycols <- c("#0073C2FF", "#EFC000FF", "snow4", "#CD534CFF","darkgreen")


t <- ggplot(count.data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  xlim(0.5, 2.5)

count.data$prop <- as.character(paste0(count.data$prop,"%"))
tiff("matvspat.tiff", units="in", width=16, height=10, res=300)
t+ geom_label_repel(aes(y = lab.ypos, label = count.data$prop), color = "white",size=17)+
  guides(fill = guide_legend(title = "Group"))+
  theme_void()+ theme(legend.position="none")
dev.off()

legend <- cowplot::get_legend(t)
tiff("legendmatvspat.tiff", units="in", width=16, height=10, res=300)

grid.draw(legend)
dev.off()
