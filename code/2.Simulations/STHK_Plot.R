library(ggplot2)
require(cowplot)
STHKTable <- read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Simulations/Rplots/Data_LTMv2/Sim_Rdata.txt", header=T)
STHKObs <- read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Simulations/Rplots/Data_LTMv2/Obs_Rdata.txt", header=T)
STHKTable$Treatment1 <- factor(STHKTable$Treatment, c("STS-S","LTS-S","LTS-M2"))
STHKObs$Treatment1 <- factor(STHKObs$Treatment, c("STS-S","LTS-S","LTS-M2"))
ggplot(data=STHKTable, aes(x=Bin, y=Count,group=Bin))+
  geom_boxplot(aes(size=3))+
  geom_point(data=STHKObs,aes(x=Bin, y=Count, col=Treatment1, size=3))+
  theme_classic()+theme(text = element_text(size=24),axis.text.x = element_text(angle=45, hjust=1),legend.position = "none")+facet_grid(~ Treatment1,scale="free")+xlab("Derived Allele Frequency")+ylab("Number of Polymorphisms")

