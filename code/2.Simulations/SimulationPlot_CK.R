library(ggplot2)
library(ggpmisc)
library(tibble)
library(grid)
library(gridExtra)
library(gtable)
require(cowplot)
library(TeachingDemos)
# modified 2023-11-29 Canan 
# libraries 
library(tidyverse)

mytheme <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=12))+
  theme(axis.text = element_text(size = 14), axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1))+
  theme(strip.text.x = element_text(size = 14))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1)) +
  theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
  theme(axis.title.x = element_text(margin=margin(10,0,0)),
        axis.title.y = element_text(margin=margin(0,10,0,0)),
        axis.text.x = element_text(margin=margin(10,0,0,0)),
        axis.text.y = element_text(margin=margin(0,10,0,0)))

# load data
BsubMut_Sim <- read.table("~/GitHub/SporMut/data/2.Simulations/RPlots/FinalData/Bacillus_MutCount_Simulation.txt", sep="\t", header=TRUE, check.names=FALSE)
BsubMut_Obs <- read.table("~/GitHub/SporMut/data/2.Simulations/RPlots/FinalData/Bacillus_MutCount_Observed.txt", sep="\t", header=TRUE, check.names=FALSE)

# tidy 
BsubMut_SimTidy <- BsubMut_Sim %>%
  pivot_longer(cols = -c("Sample","Treatment","DataType"), names_to = "Bin", values_to = "Polymorphisms") %>%
  mutate_at(c('Bin', 'Polymorphisms'), as.numeric) %>%
  dplyr::group_by(Treatment, Bin) %>%
  dplyr::summarise(meanPoly = mean(Polymorphisms), 
                   sdPoly = sd(Polymorphisms), 
                   sePoly = sd(Polymorphisms)/sqrt(length(Polymorphisms)))


# New facet label names 
pop.labs <- c("S[10]", "S[1000]", "S+NS[1000]")
names(pop.labs) <- c("Early spore fraction", "Endpoint spore fraction", "Endpoint total fraction")

BsubMut_SimTidy$Treatment <- factor(BsubMut_SimTidy$Treatment, levels = c("S[10]", "S[1000]", "S+NS[1000]"), 
                  labels =  c("Early spore fraction", "Endpoint spore fraction", "Endpoint total fraction"))


BsubSimulation<-ggplot()+
  geom_point(data =BsubMut_SimTidy, aes(x = Bin, y = meanPoly), color = "grey25", alpha = .5, size =3)+
  geom_errorbar(data =BsubMut_SimTidy,aes(x = Bin, ymax = meanPoly+sdPoly, ymin = meanPoly-sdPoly))+
  #geom_point(data=BsubMut_Obs, aes(x=Bin, y=Count, color=Treatment),size=2)+
  facet_wrap(.~Treatment, nrow=3, strip.position="top", scales = "free_y")+
  mytheme+
  ylab("Number of polymorphysims")+
  xlab("Unfolded allele frequency")+
  theme(strip.background = element_blank())+
    scale_y_continuous(sec.axis=dup_axis())+
  scale_x_continuous(breaks = c(0, 0.5, 1))

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
