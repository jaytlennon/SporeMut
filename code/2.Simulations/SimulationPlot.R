library(ggplot2)
library(ggpmisc)
library(tibble)
library(grid)
library(gridExtra)
library(gtable)
require(cowplot)
library(TeachingDemos)

FilmTest<-read.table("/Users/megbehri/Desktop/LynchLab/JohnMuri/Biofilm/TestBiofilmDataset.txt", sep="\t", header=TRUE)

FilmTest.aov<-aov(Biofilm ~Treatment, data=FilmTest)
summary(FilmTest.aov)
TukeyHSD(FilmTest.aov)

ggplot(data=FilmTest, aes(x=Treatment, y=Biofilm, color=Treatment))+
  geom_jitter()+
  stat_summary(fun.y = mean, geom="point", pch=3, colour="black") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", colour="black")


BsubMut_Sim<-read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Simulations/Rplots/FinalData/Bacillus_MutCount_Simulation.txt", sep="\t", header=TRUE, check.names=FALSE)
BsubMut_Sim_Melt<-melt(BsubMut_Sim, id=c("Sample","Treatment","DataType"))
BsubMut_Sim_Melt
names(BsubMut_Sim_Melt)<-c("Sample","Treatment","DataType","Bin","Polymorphisms")
BsubMut_Sim_Melt$Bin<-as.numeric(as.character(BsubMut_Sim_Melt$Bin))
BsubMut_Sim_Melt$Polymorphisms<-as.numeric(as.character(BsubMut_Sim_Melt$Polymorphisms))
#BsubMut_Sim_Melt$Treatment <- factor(BsubMut_Sim_Melt$Treatment, c("STS-S","LTS-S","LTS-M"))
BsubMut_Obs<-read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Simulations/Rplots/FinalData/Bacillus_MutCount_Observed.txt", sep="\t", header=TRUE)
BsubMut_Obs
names(BsubMut_Obs)<-c("DataType","Treatment","Bin","Polymorphisms")
BsubMut_Obs$Bin<-as.numeric(BsubMut_Obs$Bin)
BsubMut_Obs$Polymorphisms<-as.numeric(BsubMut_Obs$Polymorphisms)

#BsubMut_Obs$Treatment <- factor(BsubMut_Obs$Treatment, c("STS-S","LTS-S","LTS-M"))
Expected<-BsubMut_Sim_Melt[which(BsubMut_Sim_Melt$Bin==0.06& BsubMut_Sim_Melt$Treatment=='S+NS[1000]'),]
sd(Expected$Polymorphisms)
TestValue<-c(5)

t.test(Expected$Polymorphisms,mu=4)


#S1000
z.test(5,mu=0.02,0.14,alternative="two.sided",n=1)

BsubSimulation<-ggplot()+
  stat_summary(data=BsubMut_Sim_Melt, aes(x=Bin, y=Polymorphisms), fun.y = mean, geom="point", pch=3, colour="#888888")+
  stat_summary(data=BsubMut_Sim_Melt, aes(x=Bin, y=Polymorphisms),fun.data = "median_hilow", geom = "errorbar", colour="#888888")+
  geom_point(data=BsubMut_Obs, aes(x=Bin, y=Polymorphisms, color=Treatment),size=2)+
  xlim(c(0,0.5))+
  scale_color_manual(values=c("#CC79A7BB","#56B4E9BB","#0072B2BB"))+
  facet_wrap(.~Treatment, nrow=3,labeller=label_parsed, strip.position="top")+  
  theme_bw()+ylab("Number of Polymorphisms")+xlab("Unfolded Allele Frequency")+
  #theme(axis.text.x = element_text(angle=45,hjust=1, size=8),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),strip.background =  element_rect(fill = NA, colour = NA),legend.position="none")
  theme(axis.text.x = element_text(angle=45,hjust=1, size=8),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")


BsubCircos<-ggdraw()+
  draw_image("/Users/megbehri/Desktop/LennonLab/Bacillus/Circos_Plot_Figure_2/LT_NotTreated/SNP_Figure2A_Edited.gene-1.pdf")

plot_grid(BsubCircos,BsubSimulation,labels=c("A","B"),rel_widths = c(1,0.5))       
