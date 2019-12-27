#!/usr/bin/env Rscript
library(ggplot2)
library(grid)
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

setwd(args[1])
stat <- read.table(args[2],header =TRUE)
grob <- grobTree(textGrob(as.character(args[7]), x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))
png(args[3], width = 6, height = 2, units = "in", res = 300)
ggplot() + geom_line(data=stat, aes(x=start, y=log10(thetaW)),colour="orange") + 
  geom_line(data=stat, aes(x=start, y=log10(pi)),colour="red") +
  scale_y_continuous(limits = c(-5, 0),breaks = c(-5,-4,-3,-2,-1,0),labels=c("10-5","10-4",0.001,0.01,0.1,1)) +
  scale_x_continuous(breaks=seq(0,ceiling(max(stat$start)/1000000)*1000000,by = 1000000),labels=seq(0,ceiling(max(stat$start)/1000000),by = 1)) +
  labs(x = "Position on chromosome (Mb)", y = NULL, title = NULL)+
  annotation_custom(grob)
dev.off()

png(args[4], width = 6, height = 2, units = "in", res = 300)
ggplot() + geom_line(data=stat, aes(x=start, y=tajD)) + 
  geom_hline(aes(yintercept = 0),colour="black", linetype="dashed") + 
  scale_y_continuous(limits = c(-3, 4),breaks = c(-3,-2,-1,0,1,2,3,4),labels=c(-3,-2,-1,0,1,2,3,4)) +
  scale_x_continuous(breaks=seq(0,ceiling(max(stat$start)/1000000)*1000000,by = 1000000),labels=seq(0,ceiling(max(stat$start)/1000000),by = 1)) +
  labs(x = "Position on chromosome (Mb)", y = NULL, title = NULL)+
  annotation_custom(grob)
dev.off()