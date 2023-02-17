library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
data <- read.csv(file="GF.csv",header=T)

data <- data[data$source %like% "BP", ]

data <- data[,c("term_id","adjusted_p_value","term_name","negative_log10_of_adjusted_p_value")]

write.table(data,file="BP.txt",sep="\t")

rm(data)
data <- read.table(file="BP.txt",header=T)


data <- data.frame(data$term_name,data$negative_log10_of_adjusted_p_value)

data <- data[order(data$data.negative_log10_of_adjusted_p_value,decreasing = T),]

#data <- data[1:30,]

colnames(data) <- c("x_axis","log10Padj")

tiff("BP.tiff", units="in", width=9, height=1.5, res=300)
# Reorder following the value of another column:
m <- data %>%
  mutate(x_axis = fct_reorder(x_axis, (log10Padj)))%>%
  ggplot( aes(x=x_axis, y=(log10Padj))) +
  geom_bar(stat="identity", fill="dodgerblue4", alpha=.6, width=.4,colour="dodgerblue4") +
  coord_flip() +
  ylab("-log10(Padj)") +
  xlab("") +
  theme_bw()
m+theme_classic() + theme(text = element_text(size = 15)) +
  coord_flip(ylim=c(min(data$log10Padj),max(data$log10Padj))) +
  theme(text = element_text(face="bold"))+
  theme(text=element_text(family="sans"))         
