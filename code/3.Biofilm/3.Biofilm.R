library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)
library(dplyr)
library(scales)

BacillusBiofilm<-read.table("/Users/megbehri/Desktop/LennonLab/Bacillus/BacillusDataFinal/Biofilm/Biofilm_06_16_20.txt", sep= "\t", header=FALSE)
BacillusBiofilm


###Format Datasets By Plate#######
###100 Lines
Plate1BS<- BacillusBiofilm[which(BacillusBiofilm$V1=='plate1'),] 
Plate2BS<- BacillusBiofilm[which(BacillusBiofilm$V1=='plate2'),] 

names(Plate1BS)<- as.matrix(Plate1BS[1,])
Plate1BS<- Plate1BS[-1,]
Plate1BS[]<-lapply(Plate1BS,function(x) type.convert(as.character(x)))
Plate1BS
names(Plate2BS)<- as.matrix(Plate2BS[1,])
Plate2BS<- Plate2BS[-1,]
Plate2BS[]<-lapply(Plate2BS,function(x) type.convert(as.character(x)))
Plate2BS

Plate1BSMelt<-melt(Plate1BS,id=c("plate1", "Replicate"))
names(Plate1BSMelt)<-c("Plate","Replicate","Sample","OD550")
Plate1BSMelt

Plate2BSMelt<-melt(Plate2BS,id=c("plate2", "Replicate"))
names(Plate2BSMelt)<-c("Plate","Replicate","Sample","OD550")
Plate2BSMelt

####Identify Samples####

#NonSpore <- BacillusCurveMelt$Sample %in% c("M4", "M13", "M17" , "M19", "M21", "M41", "M54", "M79")
#Spore<-BacillusCurveMelt$Sample %in% c( "M23", "M26", "S1", "S6", "S11", "S22", "S51", "S95")

Blank <- Plate1BSMelt$Sample %in% c("Blank1")
Spore<-Plate1BSMelt$Sample %in% c( "M23", "M26")
Vegetative<-Plate1BSMelt$Sample %in% c("M4", "M13", "M17" , "M19", "M21", "M41", "M54", "M79")
Ancestor<-Plate1BSMelt$Sample %in% c( "Ancestor")
sinR<-Plate1BSMelt$Sample %in% c("M17" , "M19", "M21", "M41", "M54")
ywcC<-Plate1BSMelt$Sample %in% c("M4", "M13", "M79")


Plate1BSMelt$Treatment[Blank] <- "Blank"
Plate1BSMelt$Treatment[Spore] <- "Spore"
Plate1BSMelt$Treatment[Vegetative] <- "Vegetative"
Plate1BSMelt$Treatment[Ancestor] <- "Ancestor"
Plate1BSMelt$BiofilmMut[Ancestor] <- "Ancestor"
Plate1BSMelt$BiofilmMut[sinR] <- "Vegetative:sinR"
Plate1BSMelt$BiofilmMut[ywcC] <- "Vegetative:ywcC/slrR"
Plate1BSMelt$BiofilmMut[Spore] <- "Spore"
Plate1BSMelt$BiofilmMut[Blank] <- "Blank"

Plate1BSMelt

Blank <- Plate2BSMelt$Sample %in% c("Blank1", "Blank2", "Blank3", "Blank4", "Blank5", "Blank6")
Spore<-Plate2BSMelt$Sample %in% c( "S1", "S6","S22","S52","S95")

Plate2BSMelt$Treatment[Blank] <- "Blank"
Plate2BSMelt$Treatment[Spore] <- "Spore"
Plate2BSMelt$BiofilmMut[Blank] <- "Blank"
Plate2BSMelt$BiofilmMut[Spore] <- "Spore"

Plate2BSMelt

############Correct OD for Blank######
Plate1BSBlankOnly<- Plate10Melt[which(Plate1BSMelt$Treatment=='Blank'),] 
Plate2BSBlankOnly<- Plate11Melt[which(Plate2BSMelt$Treatment=='Blank'),] 

Plate1BSMelt$BlankAvg<-median(Plate1BSBlankOnly$OD550)
Plate2BSMelt$BlankAvg<-median(Plate2BSBlankOnly$OD550)

Plate1BSMelt$OD550_Corrected<-Plate1BSMelt$OD550-Plate1BSMelt$BlankAvg
Plate2BSMelt$OD550_Corrected<-Plate2BSMelt$OD550-Plate2BSMelt$BlankAvg

Plate1BSMelt$OD550_Corrected[Plate1BSMelt$OD550_Corrected<=0] <- 0.00001
Plate2BSMelt$OD550_Corrected[Plate2BSMelt$OD550_Corrected<=0] <- 0.00001


##Combine Plates and Normalize by Ancestor##
BacillusBiofilmPlates_withBlanks<-rbind(Plate1BSMelt,Plate2BSMelt)
BacillusBiofilmPlates_withBlanks
#BacillusBiofilmPlates<-BacillusBiofilmPlates_withBlanks[which(c(BacillusBiofilmPlates_withBlanks$Treatment=='Spore'|BacillusBiofilmPlates_withBlanks$Treatment=='Vegetative'| BacillusBiofilmPlates_withBlanks$Treatment=='Ancestor')),]
BacillusBiofilmPlates<-BacillusBiofilmPlates_withBlanks[which(c(BacillusBiofilmPlates_withBlanks$Treatment=='Spore'|BacillusBiofilmPlates_withBlanks$Treatment=='Vegetative')),]
BacillusBiofilmPlates
BacillusBiofilmAnsOnly<-BacillusBiofilmPlates_withBlanks[which(BacillusBiofilmPlates_withBlanks$Treatment=='Ancestor'),]
BacillusBiofilmAnsOnly
BSAncestor_BFAverage<-mean(BacillusBiofilmAnsOnly$OD550_Corrected)
sem(BacillusBiofilmAnsOnly$OD550_Corrected)
BacillusBiofilmPlates$AnsAvg<-mean(BacillusBiofilmAnsOnly$OD550_Corrected)
BacillusBiofilmPlates$OD_Ratio<-BacillusBiofilmPlates$OD550_Corrected/BacillusBiofilmPlates$AnsAvg

BacillusBiofilmPlates


###Variables to Plot Y axis in Scientific Notation###
BFbreaks=c(1e-1,1e0,1e1,1e2,1e3)
scientificNoOne <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("1e", "10^", scientific_format()(x)))))
}

##Make plot##

BSBiofilmPlot<-ggplot(data=BacillusBiofilmPlates, aes(x=Treatment, y=OD550_Corrected,color=Treatment))+
  geom_jitter(aes(color=Treatment))+
  stat_summary(fun.y=mean, geom="point", pch=3) +
  #stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text")+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  #stat_compare_means(comparisons = PhageBFCompare,method="t.test",aes(label="..p.adj.."), size=4)+
  #scale_color_manual(values=c('#A4A4A4FF','#A4A4A488','#e69f00FF','#e69f0088',"#D55E00FF","#D55E0088","#CC79A7FF","#CC79A788",'#56B4E9FF','#56B4E988',"#009E73FF","#009E7388"), labels = c("Ancestor","Ancestor:\u03d5 80+", "Large:1-day","Large:1-day:\u03d5 80+","Medium:1-day","Medium:1-day:\u03d5 80+","Small:1-day","Small:1-day:\u03d5 80+","Large:10-day","Large:10-day:\u03d5 80+","Large:100-day","Large:100-day:\u03d5 80+"))+
  scale_color_manual(values=c("#A4A4A4BB","#d55e00BB","#009e73BB"))+
  geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)),linetype="dashed")+
  geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)+(2*sem(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)-(2*sem(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  #scale_y_log10()+
  #coord_trans(y = "log10",limy = c(1e-1, 1e1))+ 
  coord_trans(y = "log10",limy = c(1e-3, 5))+ 
  scale_y_continuous(breaks = BFbreaks,label=scientificNoOne,)+
  xlab("Treatment")+ylab(expression(paste((Abs[550])[Sample])))+
  #ylab("Abs(550) Sample /Mean Abs(OD550) WT Ancestor ")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=12),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


BSBiofilmGenePlot<-ggplot(data=BacillusBiofilmPlates, aes(x=BiofilmMut, y=OD550_Corrected,color=BiofilmMut))+
  geom_jitter(aes(color=BiofilmMut))+
  stat_summary(fun.y=mean, geom="point", pch=3) +
  #stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text")+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  #stat_compare_means(comparisons = PhageBFCompare,method="t.test",aes(label="..p.adj.."), size=4)+
  #scale_color_manual(values=c('#A4A4A4FF','#A4A4A488','#e69f00FF','#e69f0088',"#D55E00FF","#D55E0088","#CC79A7FF","#CC79A788",'#56B4E9FF','#56B4E988',"#009E73FF","#009E7388"), labels = c("Ancestor","Ancestor:\u03d5 80+", "Large:1-day","Large:1-day:\u03d5 80+","Medium:1-day","Medium:1-day:\u03d5 80+","Small:1-day","Small:1-day:\u03d5 80+","Large:10-day","Large:10-day:\u03d5 80+","Large:100-day","Large:100-day:\u03d5 80+"))+
  scale_color_manual(values=c("#d55e00BB","#009e73BB","#000080BB"),labels=c("Spore",expression(paste(Vegetative~(italic('sinR')))),expression(paste(Vegetative~(italic('ywcC')~or~italic('espA-slrR'))))))+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)),linetype="dashed")+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)+(2*se(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)-(2*se(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  #scale_y_log10()+
  #coord_trans(y = "log10",limy = c(1e-1, 1e1))+ 
  coord_trans(y = "log10",limy = c(1e-3, 5))+ 
  scale_y_continuous(breaks = BFbreaks,label=scientificNoOne,)+
  xlab("")+ylab("Absorbance (550nm)")+
  #ylab("Abs(550) Sample /Mean Abs(OD550) WT Ancestor ")+
  theme_bw()+theme(legend.position="none")+
  #theme(axis.text.x = element_text(angle=19,hjust=1, size=12),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position=c(0.1,0.22))
  theme(axis.text.x = element_blank(),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #theme(legend.position=c(0.18,0.15),legend.title=element_text(size =10,face="bold"),legend.title.align=0.5,legend.text = element_text(size=9),legend.text.align = 0.5,legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  #labs(col="Ecotype")
  
BSBiofilmGenePlot
  
plot_grid(BSBiofilmPlot,BSBiofilmGenePlot,labels=c("A","B"),rel_widths = c(.75,1), rel_heights =c(.75,1))



BSBiofilmPlot
BacillusBiofilmPlates

ggplot(data=BacillusBiofilmPlates, aes(x=Treatment, y=OD550_Corrected,fill=Treatment))+
  geom_column(aes(color=Treatment))+
  #stat_summary(fun.y=mean, geom="point", pch=3) +
  #stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text")+
  #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  #stat_compare_means(comparisons = PhageBFCompare,method="t.test",aes(label="..p.adj.."), size=4)+
  #scale_color_manual(values=c('#A4A4A4FF','#A4A4A488','#e69f00FF','#e69f0088',"#D55E00FF","#D55E0088","#CC79A7FF","#CC79A788",'#56B4E9FF','#56B4E988',"#009E73FF","#009E7388"), labels = c("Ancestor","Ancestor:\u03d5 80+", "Large:1-day","Large:1-day:\u03d5 80+","Medium:1-day","Medium:1-day:\u03d5 80+","Small:1-day","Small:1-day:\u03d5 80+","Large:10-day","Large:10-day:\u03d5 80+","Large:100-day","Large:100-day:\u03d5 80+"))+
  scale_fill_manual(values=c("#A4A4A4BB","#d55e00BB","#009e73BB"))+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)),linetype="dashed")+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)+(2*sem(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  #geom_hline(data=BacillusBiofilmAnsOnly,aes(yintercept=mean(BacillusBiofilmAnsOnly$OD550_Corrected)-(2*sem(BacillusBiofilmAnsOnly$OD550_Corrected))),linetype="solid")+
  #scale_y_log10()+
  #coord_trans(y = "log10",limy = c(1e-1, 1e1))+ 
  #coord_trans(y = "log10",limy = c(1e-3, 5))+ 
  scale_y_continuous(breaks = BFbreaks,label=scientificNoOne,)+
  xlab("Treatment")+ylab(expression(paste((Abs[550])[Sample])))+
  #ylab("Abs(550) Sample /Mean Abs(OD550) WT Ancestor ")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=12),text = element_text(size=12),axis.text.y = element_text(size=12),axis.line=element_line(color="#888888"),strip.placement.x = "outside",strip.text.x=element_text(size=12, color="#000000",angle=360),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")


BacillusBiofilmPlates %>% group_by(Treatment) %>% summarize(Biofilm.var=var(OD550_Corrected))


BacillusBiofilmPlates.Spore<-BacillusBiofilmPlates[which(BacillusBiofilmPlates$Treatment=='Spore'),]
BacillusBiofilmPlates.Vegetative<-BacillusBiofilmPlates[which(BacillusBiofilmPlates$Treatment=='Vegetative'),]
BacillusBiofilmPlates.Vegetative$OD550_Corrected

var.test(BacillusBiofilmPlates.Vegetative$OD550_Corrected,BacillusBiofilmPlates.Spore$OD550_Corrected,ratio=1,alternative="two.sided")
compare_means(OD550_Corrected~Treatment,BacillusBiofilmPlates)
t.test(BacillusBiofilmPlates.Vegetative$OD550_Corrected,BacillusBiofilmPlates.Spore$OD550_Corrected,alternative="two.sided")

compare_means(OD550_Corrected~BiofilmMut,BacillusBiofilmPlates, method="t.test")
