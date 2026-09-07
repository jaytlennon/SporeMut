library(ggplot2)
library(ggpmisc)
library(hm)
library(tibble)
library(grid)
library(gridExtra)
library(gtable)
require(cowplot)
library(reshape2)
library(plyr)
library(growthcurve)
library(ggpubr)

BacillusCurveData<-read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Fitness/WeightedAverage/BacillusFitnessWA_Trun10h.txt",sep="\t", header=TRUE)
BacillusCurveMelt<-melt(BacillusCurveData,id="Time")
names(BacillusCurveMelt)<-c("Time","Sample","OD600")
NonSpore <- BacillusCurveMelt$Sample %in% c("M4", "M13", "M17" , "M19", "M21", "M41", "M54", "M79")
Spore<-BacillusCurveMelt$Sample %in% c( "M23", "M26", "S1", "S6", "S11", "S22", "S51", "S95")
BacillusCurveMelt$CellType[NonSpore] <- "Vegetative"
BacillusCurveMelt$CellType[Spore] <- "Spore"
BacillusCurveMelt
Growth_sinR<-BacillusCurveMelt$Sample %in% c("M17" , "M19", "M21", "M41", "M54")
Growth_ywcC<-BacillusCurveMelt$Sample %in% c("M4", "M13", "M79")

BacillusCurveMelt$BiofilmGene[Spore] <- "Spore"
BacillusCurveMelt$BiofilmGene[Growth_sinR] <- "Vegetative:sinR"
BacillusCurveMelt$BiofilmGene[Growth_ywcC] <- "Vegetative:ywcC/slrR"

BacillusRepData<-read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Fitness/WeightedAverage/BacillusReplicateStats_Trunc10_5.txt",sep="\t", header=TRUE)
BacillusRepPlotData<-subset(BacillusRepData, , -c(b0, A, z, lag,umax.lw,umax.up, umax.lw.FI, umax.up.FI))
BacillusRepPlotData
BacillusRepMelt<-melt(BacillusRepPlotData, id=c("Curve","Sample","CellType","Run"))
names(BacillusRepMelt)<-c("Curve","Sample","CellType","Run", "Parameter","Value")
BacillusRepMelt$Par2 <- factor(BacillusRepMelt$Parameter, labels = c( "Lag","u[max]", "K"))
BacillusRepMelt
GrowthRep_sinR<-BacillusRepMelt$Sample %in% c("M17" , "M19", "M21", "M41", "M54")
GrowthRep_ywcC<-BacillusRepMelt$Sample %in% c("M4", "M13", "M79")
RepSpore<-BacillusRepMelt$Sample %in% c( "M23", "M26", "S1", "S6", "S11", "S22", "S51", "S95")
BacillusRepMelt$BiofilmGene[RepSpore] <- "Spore"
BacillusRepMelt$BiofilmGene[GrowthRep_sinR] <- "Vegetative:sinR"
BacillusRepMelt$BiofilmGene[GrowthRep_ywcC] <- "Vegetative:ywcC/slrR"
BacillusRepMelt

BacillusRepMelt$Parameter <- factor(BacillusRepMelt$Parameter, levels = c("Lag","u[max]","K"))


BacGrowthPlot<-ggplot()+
  geom_point(data=BacillusCurveMelt, aes(x=Time,y=OD600, color=CellType))+
  stat_growthcurve(data=BacillusCurveMelt, aes(x=Time,y=OD600, group=CellType, color=CellType), model="logistic4p",position = "identity", geom="line")+
  ylab("Absorbance (OD600)")+xlab("Time (hours)")+scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values=c("#d55e00BB","#009e73BB"))+
  theme_classic()+theme(axis.text.x = element_text(size=12),text = element_text(size=12),axis.text.y = element_text(size=12), axis.line=element_line(color="#888888"), legend.position="top", legend.title=element_blank())
BacGrowthPlot
BacGrowthStatPlot<-ggplot(data=BacillusRepMelt, aes(x=CellType, y=Value, color=CellType))+
  geom_jitter()+
  stat_summary(fun.y = mean, geom="point", pch=3, colour="#888888") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", colour="#888888")+
  facet_wrap(~Par2,scales = "free_y",strip.position="left", labeller=label_parsed)+
  #stat_compare_means(method="t.test",label="p.format")+
  ylab(NULL)+xlab(NULL)+
  scale_color_manual(values=c("#d55e00BB","#009e73BB"))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle=45,hjust=1, size=12),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.background = element_blank(),strip.placement ="outside",strip.text=element_text(size=12),legend.position="none")
  theme(axis.text.x = element_blank(),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.background = element_blank(),strip.placement ="outside",strip.text=element_text(size=12),legend.position="none")
BacGrowthStatPlot
compare_means(Value~CellType, data=BacillusRepMelt, group.by = "Parameter", method = "t.test")
compare_means(Value~BiofilmGene, data=BacillusRepMelt, group.by = "Parameter", method = "t.test")

TwoGrowthGroups<-plot_grid(BacGrowthPlot,BacGrowthStatPlot, nrow=2, labels=c("A","B"))

BacGrowthPlotGenes<-ggplot()+
  geom_point(data=BacillusCurveMelt, aes(x=Time,y=OD600, color=BiofilmGene))+
  stat_growthcurve(data=BacillusCurveMelt, aes(x=Time,y=OD600, group=BiofilmGene, color=BiofilmGene), model="logistic4p",position = "identity", geom="line")+
  ylab("Absorbance (600nm)")+xlab("Time (hours)")+scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values=c("#d55e00BB","#009e73BB","#000080BB"),labels=c("Spore",expression(paste(Vegetative~(italic('sinR')))),expression(paste(Vegetative~(italic('ywcC')~or~italic('espA-slrR'))))))+
  theme_bw()+
  theme(axis.text.x = element_blank(),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position=c(0.25,0.80),legend.title=element_text(size =10,face="bold"),legend.title.align=0.5,legend.text = element_text(size=9),legend.text.align = 0.5,legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  labs(col="Ecotype")

BacGrowthPlotGenes
BacGrowthStatPlotGenes<-ggplot(data=BacillusRepMelt, aes(x=BiofilmGene, y=Value, color=BiofilmGene))+
  geom_jitter()+
  stat_summary(fun.y = mean, geom="point", pch=3) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  facet_wrap(~Par2,scales = "free_y",strip.position="left", labeller=label_parsed)+
  #stat_compare_means(method="t.test",label="p.format")+
  ylab(NULL)+xlab(NULL)+
  scale_color_manual(values=c("#d55e00BB","#009e73BB","#000080BB"),labels=c("Spore",expression(paste(Vegetative~(italic('sinR')))),expression(paste(Vegetative~(italic('ywcC')~or~italic('slrR'))))))+
  theme_bw()+
  #theme(axis.text.x = element_text(angle=45,hjust=1, size=12),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.background = element_blank(),strip.placement ="outside",strip.text=element_text(size=12),legend.position="none")
  theme(axis.text.x = element_blank(),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.background = element_blank(),strip.placement ="outside",strip.text=element_text(size=12),legend.position="none",,panel.grid.major = element_blank(), panel.grid.minor = element_blank())
BacGrowthStatPlotGenes

BSBiofilmGenePlot
TopRepGroups<-plot_grid(BSBiofilmGenePlot,BacGrowthPlotGenes, nrow=1, labels=c("A","B"),rel_widths = c(0.5,1))
BottomRepGroups<-plot_grid(BacGrowthStatPlotGenes,labels=c("C"))
plot_grid(TopRepGroups,BottomRepGroups,nrow=2,hjust=1)

TwoRepGroups
expression(paste((Abs[550])[Sample]))

plot_grid(TwoGrowthGroups,TwoRepGroups)


VegLag <- BacillusRepMelt[which(BacillusRepMelt$CellType=='Vegetative'& BacillusRepMelt$Par2=='Lag'),]  
VegLag
SporeLag <- BacillusRepMelt[which(BacillusRepMelt$CellType=='Spore'& BacillusRepMelt$Par2=='Lag'),]  
(mean(VegLag$Value)-mean(SporeLag$Value))*60
