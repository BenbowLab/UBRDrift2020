#Figures May 2020
############
##Import Step
###########
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ggpubr)
library(multcompView)
library(lmerTest)
library(lme4)
library(tiff)
library(glmmTMB)
library(bbmle)
library(GGally)
library(nlme)
library(DHARMa)
library(patchwork)

set.seed(45682)
#DRIFT FILES all taxa
otufull=read.table("DataClean\\SturgeonDriftInvertAbundanceMatrix8.15.19.txt",header=TRUE)
#head(otufull)
metadata=read.csv("DataClean\\SturgeonDriftMetadata5.5.2020.csv",header=TRUE)
metadata<-subset(metadata,Ninverts!=0) #Remove sample dates where  inverts were not sampled
metadata$PercentRiverDischargeSampled<-metadata$DischargeSampledByNight/metadata$Q*100
metadata$InvertsByDischargeSampled<-metadata$Ninverts100/metadata$DischargeSampledByNight
metadata$InvertsByRiverDischarge<-metadata$Ninverts100/metadata$Q
metadata$BiomassByRiverDischarge<- metadata$InvertBiomass100/metadata$Q
metadata$BiomassByDischargeSampled<-metadata$InvertBiomass100/metadata$DischargeSampledByNight
metadata$DriftInvertConc<-((metadata$Ninverts100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
metadata$DriftBiomassConc<-((metadata$InvertBiomass100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water
metadata$SturgeonConc<-((metadata$Nsturgeon*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate sturgeon larvae/ 100 m3 water
metadata$SuckerConc<-((metadata$Nsuckers100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate sturgeon larvae/ 100 m3 water


# 
# HistDischarge<-ggplot(metadata,aes(x=DischargeSampledByNight))+geom_histogram(fill="grey",color="black")+xlab(expression(Discharge~(m^3/sec)~Sampled))+ylab("Frequency")#+facet_wrap(~Year)
# HistDischarge
# theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
# mean(metadata$DischargeSampledByNight,na.rm=T)
# sd(metadata$DischargeSampledByNight,na.rm=T)
# 
# dev.off()
# tiff("Figures/Discharge_Sampled.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
# HistDischarge
# dev.off()
# 
# 
# HistQ<-ggplot(metadata,aes(x=Q))+geom_histogram(fill="grey",color="black")+xlab(expression(River~Discharge~(m^3/sec)))+ylab("Frequency")#+facet_wrap(~Year)
# HistQ
# theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
# 
# dev.off()
# tiff("Figures/Discharge_Q.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
# HistQ
# dev.off()


#discharge sampled sd = 0.407179 mean = 0.9361
#metadataOutliers

taxmatrixfull=as.matrix(read.table("DataClean\\SturgeonDriftInvertTaxNames8.15.19.txt"))
#head(taxmatrixfull)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
#OTU
TAX=tax_table(taxmatrixfull)
colnames(TAX)=("Family")

sampdat=sample_data(metadata)
sample_names(sampdat)=sampdat$SampleID
#head(sampdat)
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)#joins together OTU,TAX, and metadata 
physeq
levels(sample_data(physeq)$MoonPhase)
sample_data(physeq)$MoonPhase = factor(sample_data(physeq)$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))
#levels(sample_data(physeq)$CODE)=c("NM","WXC","FQ","WXG","FM","WAG","LQ","WNC")

Richness=plot_richness(physeq, x="DPFS", measures=c("Shannon"))#+geom_boxplot(aes(x=DPFS, y=value, color=DPFS), alpha=0.05)
#Richness$data

write.csv(Richness$data, "SturgeonMetadataWDiversity.csv")
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)
AllData<-metadata
AllData$Date2<-as.Date(AllData$Date,format= "%d-%B")

####################
#Figure 1
########################

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

############
#Fig 1A Invert By day and year
############
Trtdata <- ddply(AllData, c("Date2","Year"), summarise,
                 N    = length(Ninverts100),
                 meanInvert = mean(Ninverts100)
)
#Trtdata
InvertByDayOfYear<-ggplot(Trtdata, aes(x=Date2,y=meanInvert))+geom_bar(colour="black", stat="identity")+xlab("Date")+ylab("Macroinvertebrate Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
InvertByDayOfYear
##########
#Figure 1B Sucker abu by day and year
##########

Trtdata <- ddply(AllData, c("Date2","Year"), summarise,
                 N    = length(Nsuckers100),
                 meanSuckers = mean(Nsuckers100)
)

SuckersByDayOfYear<-ggplot(Trtdata, aes(x=Date2,y=meanSuckers))+geom_bar(colour="black", stat="identity")+xlab("Date")+ylab("Larval Catostomidae Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
SuckersByDayOfYear

##########
# Figure 1C sturgeon larvae abu by day and year 
##########
head(AllData)
Trtdata <- ddply(AllData, c("Date2","Year"), summarise,
                 N    = length(Nsturgeon),
                 meanSturgeon = mean(Nsturgeon)
)
#Trtdata
SturgeonByDayOfYear<-ggplot(Trtdata, aes(x=Date2,y=meanSturgeon))+geom_bar(colour="black", stat="identity")+xlab("Date")+ylab("Larval Sturgeon Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
SturgeonByDayOfYear


dev.off()
tiff("Figures/SturgeonByDayAndYear.tiff", width = 90, height = 110, units = 'mm', res = 1200)
SturgeonByDayOfYear
dev.off()


#########
#Old Figure 1C Sturgeon larvae CTU from first spawn
#########

# SpawningDischarge<-read.csv("DataClean\\SpawningSiteDischarge2012-2018JR.csv",header=T)
# head(SpawningDischarge)
# 
# SturgeonCTU<-ggplot(AllData,aes(x=CTUSturgeon))+geom_point(size=0.35,aes(y=Nsturgeon),shape=19)+geom_line(aes(y=DischargeSpawning/0.002,x=SturgeonCTU),data=SpawningDischarge,color="blue")+scale_y_continuous(sec.axis = sec_axis(~.*0.002,name=expression(River~Discharge~In~Spawning~Area~(m^3/sec))))+xlab("Cumulative Temperature Units")+facet_wrap(~Year)+ylab("Larval Sturgeon Abundance")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# SturgeonCTU
Trtdata <- ddply(AllData, c("CTUSturgeon","Year"), summarise,
                 N    = length(Nsturgeon),
                 meanSturgeon = mean(Nsturgeon)
)
ggplot(Trtdata, aes(x=CTUSturgeon,y=meanSturgeon))+geom_bar(colour="black", stat="identity",width = 10)+xlab("Date")+ylab("Larval Sturgeon Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_wrap(Year~.)#+scale_fill_manual(values=cbPalette)


SturgeonCTU<-ggplot(AllData,aes(x=CTUSturgeon,y=Nsturgeon))+geom_point(size=0.5)+xlab("Cumulative Temperature Units")+facet_wrap(~Year)+ylab("Larval Sturgeon Abundance")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
SturgeonCTU
head(AllData)

dev.off()
tiff("Figures/SturgeonCTU.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonCTU
dev.off()




###########
#Fig 1D Densities by CTU
##########

Sturgeon<-ggplot(ShannonSubset,aes(CTUSturgeon,SturgeonConc))+geom_point(color="#000000",shape=15,size=1)+theme(axis.title.y=element_blank())+xlab("Cumulative Temperature Units (CTU)")+ scale_y_continuous(breaks=c(0,25,50,75,100),limits=c(0,100))
Invert<-ggplot(ShannonSubset,aes(CTUSturgeon,DriftInvertConc))+geom_point(color="#E69F00",shape=17,size=1)+theme(axis.title.x=element_blank(),axis.title.y=element_blank())
Cato<-ggplot(ShannonSubset,aes(CTUSturgeon,SuckerConc))+geom_point(color="#0072B2",shape=16,size=1)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_text(size = 7))


dev.off()
tiff("Figures/Figure1D DensityByCTU.tiff", width = 87, height = 95, units = 'mm', res = 1200)
ggarrange(Invert,Cato,Sturgeon,
          labels = c("", "","",""),
          ncol = 1, nrow = 3)
dev.off()

############
#Join together Figure 1
###########

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/Fig1July2020.tiff", width = 174, height = 190, units = 'mm', res = 1200)
ggarrange(InvertByDayOfYear,SuckersByDayOfYear,SturgeonByDayOfYear,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()


#Fig1D Added in paint




##############
#Figure  Sucker GAM
##############
citation("gratia")
citation(gratia)
library(itsadug)
library(gratia)
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, DischargeSampledByNight!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

head(ShannonSubset)

SuckerAbuLog0 <- gam(log(SuckerConc+10e-5)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SuckerAbuLog1 <- gam(log(SuckerConc+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog2 <- gam(log(SuckerConc+10e-5)~ s(DayOfYear)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog3 <- gam(log(SuckerConc+10e-5)~ s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog4 <- gam(log(SuckerConc+10e-5)~ s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SuckerAbuLog5 <- gam(log(SuckerConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog6 <- gam(log(SuckerConc+10e-5)~ s(DayOfYear)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


SuckerAbuLog7 <- gam(log(SuckerConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog8 <- gam(log(SuckerConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SuckerAbuLogDischarge <- gam(log(Nsuckers100+10e-5)~ s(DischargeCentered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

ShannonSubset2<-subset(ShannonRichness, Q!= "NA")

summary(SuckerAbuLogDischarge)
SuckerAbuLogDischarge2 <- gam(log(Nsuckers100+10e-5)~ s(Q)+s(Year,bs="re"), data = ShannonSubset2, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
summary(SuckerAbuLogDischarge2)
#Nsuckers100
#SuckerConc
SuckerAbuGamma <- gam((SuckerConc+10e-5)~ s(CTUSturgeon)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SuckerAbuGamma$aic
summary(SuckerAbuGamma)
AICctab(SuckerAbuLog0,SuckerAbuLog1,SuckerAbuLog2,SuckerAbuLog3,SuckerAbuLog4,SuckerAbuLog5,SuckerAbuLog6,SuckerAbuLog7,SuckerAbuLog8,
        weights=T)
AICctab(SuckerAbuGamma,SuckerAbuLog5,weights=T)
SuckerAbuLog5$aic
SuckerAbuLog5$aic
summary(SuckerAbuLog5)
summary(SuckerAbuLog1)

# SuckerAbuLog6 <- gam(log(Nsuckers100+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# SuckerAbuLogConc <- gam(log(SuckerConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# AICctab(SuckerAbuLog6,SuckerAbuLogConc,weights=T)

#Models using drift conc performing better


##########33
#Top model
#########3
SuckerAbuLog <- gam(log(SuckerConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

dev.off()
tiff("Figures/SuckerGAMAppraisal.tiff", width = 174, height = 174, units = 'mm', res = 1200)
appraise(SuckerAbuLog5)
dev.off()



summary(SuckerAbuLog5)
plot<-plot_smooth(SuckerAbuLog, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
plot_smooth(SuckerAbuLog, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=SuckerConc))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")


SuckerGAMPlot<-ggplot(ShannonSubset,aes(CTUSturgeon,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Larval~Catostomidae~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))#+ylim(NA, 1000)
SuckerGAMPlot

PredictedDataFramelog<-data.frame((plot$fv$fit),plot$fv$CTUSturgeon,(plot$fv$ul),(plot$fv$ll))
colnames(PredictedDataFramelog)<-c("SuckerConc","CTUSturgeon","UpperCI","LowerCI")

ggplot(ShannonSubset,aes(CTUSturgeon,log(SuckerConc)))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5,aes(y=log(SuckerConc)))+  
  geom_line(data = PredictedDataFrame, aes(y = log(LowerCI)), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=log(UpperCI)),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Log(Larval~Catostomidae~Per~100~m^3~Drift)))



dev.off()
tiff("Figures/SuckerGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SuckerGAMPlot
dev.off()

#Percent Illumination
plot<-plot_smooth(SuckerAbuLog, view="percillum",rm.ranef=F,sim.ci = T)
plot_smooth(SuckerAbuLog, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$percillum,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","percillum","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=percillum,y=SuckerConc))+geom_point()+geom_line(aes(percillum, UpperCI),color="red")+geom_line(aes(percillum,LowerCI),color="red")


SuckerGAMPlot2<-ggplot(ShannonSubset,aes(percillum,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Percent Illumination (%)")+ylab(expression(Larval~Catostomidae~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))#+ylim(NA, 1000)
SuckerGAMPlot2

PredictedDataFramelog<-data.frame((plot$fv$fit),plot$fv$CTUSturgeon,(plot$fv$ul),(plot$fv$ll))
colnames(PredictedDataFramelog)<-c("SuckerConc","CTUSturgeon","UpperCI","LowerCI")

ggplot(ShannonSubset,aes(CTUSturgeon,log(SuckerConc)))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5,aes(y=log(SuckerConc)))+  
  geom_line(data = PredictedDataFrame, aes(y = log(LowerCI)), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=log(UpperCI)),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Log(Larval~Catostomidae~Per~100~m^3~Drift)))



dev.off()
tiff("Figures/SuckerGAM2.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SuckerGAMPlot2
dev.off()





##############
#Figure Sturgeon GAM
##############

library(itsadug)
library(gratia)

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, DischargeSampledByNight!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

head(ShannonSubset)

SturgeonAbuLog0 <- gam(log(SturgeonConc+10e-5)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SturgeonAbuLog1 <- gam(log(SturgeonConc+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog2 <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog3 <- gam(log(SturgeonConc+10e-5)~ s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog4 <- gam(log(SturgeonConc+10e-5)~ s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SturgeonAbuLog5 <- gam(log(SturgeonConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog6 <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


SturgeonAbuLog7 <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog7_ <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(DriftInvertConc)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog8 <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog9 <- gam(log(SturgeonConc+10e-5)~ s(DayOfYear)+s(CTUSturgeon)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SturgeonAbuLogDischarge <- gam(log(Nsturgeon+10e-5)~ s(DischargeCentered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
summary(SturgeonAbuLogDischarge)
AICctab(SturgeonAbuLog0,SturgeonAbuLog1,SturgeonAbuLog2,SturgeonAbuLog3,SturgeonAbuLog4,SturgeonAbuLog5,SturgeonAbuLog6,SturgeonAbuLog7,SturgeonAbuLog8,SturgeonAbuLog9,
        weights=T)
summary(SturgeonAbuLog7)
SturgeonAbuLog7$aic
#Gamma Models
SturgeonAbuGamma0 <- gam((SturgeonConc+10e-5)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))

SturgeonAbuGamma1 <- gam((SturgeonConc+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma2 <- gam((SturgeonConc+10e-5)~ s(DayOfYear)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma3 <- gam((SturgeonConc+10e-5)~ s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma4 <- gam((SturgeonConc+10e-5)~ s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))

SturgeonAbuGamma5 <- gam((SturgeonConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma6 <- gam((SturgeonConc+10e-5)~ s(DayOfYear)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))


SturgeonAbuGamma7 <- gam((SturgeonConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma8 <- gam((SturgeonConc+10e-5)~ s(DayOfYear)+s(temp_centered)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
#SturgeonAbuGamma9 <- gam((SturgeonConc+10e-5)~ s(CTUSturgeon)+s(percillum)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
SturgeonAbuGamma10 <- gam((SturgeonConc+10e-5)~ s(CTUSturgeon)+s(temp_centered)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))

SturgeonAbuGammaDischarge <- gam((Nsturgeon+10e-5)~ s(Q)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
summary(SturgeonAbuGammaDischarge)

summary(SturgeonAbuGamma7)
summary(SturgeonAbuGamma8)

cor(ShannonSubset$DayOfYear,ShannonSubset$CTUSturgeon)


AICctab(SturgeonAbuGamma0,SturgeonAbuGamma1,SturgeonAbuGamma2,SturgeonAbuGamma3,SturgeonAbuGamma4,SturgeonAbuGamma5,SturgeonAbuGamma6,SturgeonAbuGamma7,SturgeonAbuGamma8,
                SturgeonAbuGamma10,weights=T)
AICctab(SturgeonAbuGamma7,SturgeonAbuLog7,weights=T) #Gamma model performing much better

SturgeonAbuGamma7$aic
##########33
#Top model
#########3
summary(SturgeonAbuGamma7)
gam.check(SturgeonAbuGamma7)
plot(SturgeonAbuGamma7)
appraise(SturgeonAbuGamma7) #Model ok at best, not great


plot<-plot_smooth(SturgeonAbuGamma7, view="DayOfYear",rm.ranef=F,sim.ci = T)
plot_smooth(SturgeonAbuLog, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$DayOfYear,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SturgeonConc","DayOfYear","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=DayOfYear,y=SturgeonConc))+geom_point()+geom_line(aes(DayOfYear, UpperCI),color="red")+geom_line(aes(DayOfYear,LowerCI),color="red")

ggplot(ShannonSubset,aes(DayOfYear,SturgeonConc))+geom_point()

SturgeonGAMPlot<-ggplot(ShannonSubset,aes(DayOfYear,SturgeonConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Day of Year")+ylab(expression(Larval~Sturgeon~Per~100~m^3~Drift))#+ylim(NA, 1000)
SturgeonGAMPlot

PredictedDataFramelog<-data.frame((plot$fv$fit),plot$fv$CTUSturgeon,(plot$fv$ul),(plot$fv$ll))
colnames(PredictedDataFramelog)<-c("SturgeonConc","CTUSturgeon","UpperCI","LowerCI")

#ggplot(ShannonSubset,aes(CTUSturgeon,log(SturgeonConc)))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5,aes(y=log(SturgeonConc)))+  
  geom_line(data = PredictedDataFrame, aes(y = log(LowerCI)), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=log(UpperCI)),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Log(Larval~Catostomidae~Per~100~m^3~Drift)))



dev.off()
tiff("Figures/SturgeonGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonGAMPlot
dev.off()


dev.off()
tiff("Figures/SturgeonGAMAppraisal.tiff", width = 174, height = 174, units = 'mm', res = 1200)
appraise(SturgeonAbuGamma7)
dev.off()


#Temperature

plot<-plot_smooth(SturgeonAbuGamma7, view="temp_centered",rm.ranef=F,sim.ci = T)
#plot_smooth(SturgeonAbuLog, view="temp_centered",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$temp_centered,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SturgeonConc","temp_centered","UpperCI","LowerCI")
PredictedDataFrame$Temp<-PredictedDataFrame$temp_centered+mean(ShannonSubset$Temp) #Remove centering
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=Temp,y=SturgeonConc))+geom_point()+geom_line(aes(Temp, UpperCI),color="red")+geom_line(aes(Temp,LowerCI),color="red")

ggplot(ShannonSubset,aes(temp_centered,SturgeonConc))+geom_point()

SturgeonGAMPlot2<-ggplot(ShannonSubset,aes(Temp,SturgeonConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Temperature (C)")+ylab(expression(Larval~Sturgeon~Per~100~m^3~Drift))#+ylim(NA, 1000)
SturgeonGAMPlot2


dev.off()
tiff("Figures/SturgeonGAMTemp.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonGAMPlot2
dev.off()


#Sturgeon River Discharge

plot<-plot_smooth(SturgeonAbuGammaDischarge, view="Q",rm.ranef=F,sim.ci = T)
#plot_smooth(SturgeonAbuLog, view="temp_centered",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$Q,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("Nsturgeon","Q","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=Q,y=Nsturgeon))+geom_point()+geom_line(aes(Q, UpperCI),color="red")+geom_line(aes(Q,LowerCI),color="red")

ggplot(ShannonSubset,aes(Q,Nsturgeon))+geom_point()

SturgeonGAMPlot3<-ggplot(ShannonSubset,aes(Q,Nsturgeon))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("River Discharge (m3/sec)")+ylab(expression(Nightly~Larval~Sturgeon))#+ylim(NA, 1000)
SturgeonGAMPlot3

appraise(SturgeonAbuGammaDischarge)
dev.off()
tiff("Figures/SturgeonGAMDischarge.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonGAMPlot3
dev.off()


##########
#LMM Invertebrate Abundance Model
###########

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

min(ShannonSubset$DischargeSampledByNight,na.rm=T)
max(ShannonSubset$DischargeSampledByNight,na.rm=T)
mean(ShannonSubset$DischargeSampledByNight,na.rm=T)

min(ShannonSubset$Q,na.rm=T)
max(ShannonSubset$Q,na.rm=T)
mean(ShannonSubset$Q,na.rm=T)




InvertsByDischargeSampledDL0 = lme(log(DriftInvertConc)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL1 = lme(log(DriftInvertConc)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL2 = lme(log(DriftInvertConc)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL3 = lme(log(DriftInvertConc)~DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

InvertsByDischargeSampledDL4 = lme(log(DriftInvertConc)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL5 = lme(log(DriftInvertConc)~percillum+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL6 = lme(log(DriftInvertConc)~temp_centered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

InvertsByDischargeSampledDL7 = lme(log(DriftInvertConc)~percillum+DoYCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

InvertsByDischargeSampledDL8 = lme(log(DriftInvertConc)~CTUSturgeon,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL9 = lme(log(DriftInvertConc)~CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))





AICctab(InvertsByDischargeSampledDL0,InvertsByDischargeSampledDL1,InvertsByDischargeSampledDL2,InvertsByDischargeSampledDL3,
        InvertsByDischargeSampledDL4,InvertsByDischargeSampledDL5,InvertsByDischargeSampledDL6,InvertsByDischargeSampledDL7,
        InvertsByDischargeSampledDL8,InvertsByDischargeSampledDL9, weights=TRUE)

summary(InvertsByDischargeSampledDL3) #AICc weight = 0.786, AiCc = 437.8
summary(InvertsByDischargeSampledDL2) #AICc weight = 0.114 AICc = 441.7
summary(InvertsByDischargeSampledDL6) #
summary(InvertsByDischargeSampledDL0) #
summary(InvertsByDischargeSampledDL5) #
summary(InvertsByDischargeSampledDL4) #
summary(InvertsByDischargeSampledDL1) #
summary(InvertsByDischargeSampledDL7) #



#Model no AR1 structure
InvertsByDischargeSampledDL3NoCor = lme(log(DriftInvertConc)~DoYCentered,random =~1|Year,data=ShannonSubset)
InvertsByDischargeSampledDL3NoCor
AICctab(InvertsByDischargeSampledDL3NoCor,InvertsByDischargeSampledDL3)
#######
#Best model residual plots
#######

plot(resid(InvertsByDischargeSampledDL3) ~ temp_centered, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$temp_centered), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ percillum, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$percillum), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ DayOfYear, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$DayOfYear), col=2)


hist(resid(InvertsByDischargeSampledDL3))
plot(resid(InvertsByDischargeSampledDL3) ~ DischargeSampledByNight, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$DischargeSampledByNight), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ predict(InvertsByDischargeSampledDL3))
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ predict(InvertsByDischargeSampledDL3)), col=2)

summary(InvertsByDischargeSampledDL3)
confint(InvertsByDischargeSampledDL3)
fixef(InvertsByDischargeSampledDL3)
ranef(InvertsByDischargeSampledDL3)
coef(InvertsByDischargeSampledDL3)



#######
#Invert abu model parameter estimates
#######
summary(InvertsByDischargeSampledDL3)
intervals(InvertsByDischargeSampledDL3,which=c("fixed"))

#Intercept
exp(2.29942397) #9.97 estimate Intercept (centered)
exp(2.00731098) #7.44 CI
exp(2.591537) #13.35 CI


#DayOfYear
(exp(-0.03322902)-1)*100 #-3.27 estimate Intercept (centered)
(exp(-0.04684404)-1)*100 #-4.58 CI
(exp(-0.019614)-1)*100 #-1.94 CI


mean(ShannonSubset$DayOfYear) #158.4 june 6th or 7th


#Second best model
summary(InvertsByDischargeSampledDL2)
intervals(InvertsByDischargeSampledDL2,which=c("fixed"))

#Intercept
exp(2.2518052) #9.50 estimate Intercept (centered)
exp(1.9494452) #7.02 CI
exp(2.55416524) #12.86 CI

mean(ShannonSubset$Temp) #19.2

(exp(-0.1212527)-1)*100 #-11.41899% estimate temp centered
(exp(-0.1924751)-1)*100 #-17.51% CI
(exp(-0.05003044)-1)*100 #-4.88% CI


################
#Figure 2C Invert Abu Lmm DayOfYear plot
################

#Bootstrapped CIs rom https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)


ShannonSubset$Year<-as.factor(ShannonSubset$Year)
ShannonSubset$DriftLog<-log(ShannonSubset$DriftInvertConc)
InvertsByDischargeSampledDL3 = lme((DriftLog)~DoYCentered,random=~1|Year,data=ShannonSubset)#correlation=corAR1(form = ~DayOfYear|Year)

xvals <-  with(ShannonSubset,seq(min(DoYCentered),max(DoYCentered),length.out=100))
nresamp<-1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked <- mvrnorm(nresamp, mu = fixef(InvertsByDischargeSampledDL3), Sigma = vcov(InvertsByDischargeSampledDL3))


## predicted values: useful below
pframe <- with(ShannonSubset,data.frame(DoYCentered))
pframe$DriftLog <- predict(InvertsByDischargeSampledDL3,newdata=pframe,level=0)
head(pframe)
## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

## bootstrapping
sampfun <- function(fitted,data,idvar="Year") {
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted)
  dd <- data.frame(data,pred=pp,res=rr)
  ## sample groups with replacement
  iv <- levels(data[[idvar]])
  bsamp1 <- sample(iv,size=length(iv),replace=TRUE)
  bsamp2 <- lapply(bsamp1,
                   function(x) {
                     ## within groups, sample *residuals* with replacement
                     ddb <- dd[dd[[idvar]]==x,]
                     ## bootstrapped response = pred + bootstrapped residual
                     ddb$DriftLog <- ddb$pred +
                       sample(ddb$res,size=nrow(ddb),replace=TRUE)
                     return(ddb)
                   })
  res <- do.call(rbind,bsamp2)  ## collect results
  if (is(data,"groupedData"))
    res <- groupedData(res,formula=formula(data))
  return(res)
}

pfun <- function(fm) {
  predict(fm,newdata=pframe,level=0)
}

set.seed(101)
yvals2 <- replicate(nresamp,
                    pfun(update(InvertsByDischargeSampledDL3,data=sampfun(InvertsByDischargeSampledDL3,ShannonSubset,"Year"))))
c2 <- get_CI(yvals2,"boot_")
pframe2 <- data.frame(pframe,c2)
head(pframe2)
pframe2$DayOfYear<-pframe2$DoYCentered+mean(ShannonSubset$DayOfYear) #Remove centering on DOY for plotting
pframe2$DriftInvertConc<-exp(pframe2$DriftLog)

pframe2$boot_lwr<-exp(pframe2$boot_lwr)
pframe2$boot_upr<-exp(pframe2$boot_upr)




DOYInvertModelBoot<-ggplot(ShannonSubset,aes(DayOfYear,DriftInvertConc))+geom_point(color="darkgrey")+geom_line(data=pframe2,size=1.5)+  
  geom_line(data = pframe2, aes(y = boot_lwr), size = .75,linetype="dashed")+geom_line(data = pframe2,aes(y=boot_upr),size=0.75,linetype="dashed")+
  xlab("Day of Year")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift))
DOYInvertModelBoot

dev.off()
tiff("Figures/DOYInvertModelBoot.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DOYInvertModelBoot
dev.off()

#####################
#Combine Abundance Model Figure
#####################
DOYInvertModelBoot
SturgeonGAMPlot
SuckerGAMPlot

dev.off()
tiff("Figures/Fig2AbundanceModels.tiff", width = 174, height = 84, units = 'mm', res = 1200)

ggarrange(SturgeonGAMPlot,SuckerGAMPlot,DOYInvertModelBoot,
          labels = c("a", "b","c"),
          ncol = 3, nrow = 1)

dev.off()




###########
#Total InvertN By Moon Phase Normalized 
############

# Model<-aov(log(DriftInvertConc)~MoonPhase,data=ShannonRichness)
# summary(Model)
# hist(resid(Model))
# 
# kruskal.test(DriftInvertConc~MoonPhase, data=ShannonRichness)
# Tukey<-TukeyHSD(Model,"MoonPhase")
# 
# Tukey
# difference<-Tukey$MoonPhase[,"p adj"]
# Letters<-multcompLetters(difference)
# 
# Letters


kruskal.test(DriftInvertConc~MoonPhase, data=ShannonRichness)

compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Means=compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
Letters<-multcompLetters(difference)
Letters
Letters$Letters
#manually renamed due to some weirdness with plotting where new moon donsn't start with a

LettersRearranged<-c("bc","abc","abc","a","ab","ab","ab","c")
DriftPlotNAsRemoved<-subset(ShannonRichness, DriftInvertConc!="NA")
DriftPlotNAsRemoved
Trtdata <- ddply(DriftPlotNAsRemoved, c("MoonPhase"), summarise,
                 N    = length(DriftInvertConc),
                 meanSturgeon = mean(DriftInvertConc),
                 sd   = sd(DriftInvertConc),
                 se   = sd / sqrt(N),na.rm =T
                 
)
Trtdata

Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

TotalInvertAbuMoonPhase<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(fill=MoonPhase),stat="identity")+xlab("Moon Phase")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift~(SE)))+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 10))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+1,label=LettersRearranged))+theme(legend.position = "none")+geom_text(aes(x=5,y=25,label= "KW, Chi-sq = 39.98, P < 0.001"),size=4)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
TotalInvertAbuMoonPhase
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TotalInvertebrateAbuByPhaseNormalized.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TotalInvertAbuMoonPhase
dev.off()




compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

###############
#Hourly Results
###############


Hourly<-read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\SturgeonDrift\\SturgeonDrift\\DataClean\\SturgeonDriftDataHourly12.21.19.csv",header=T)
head(Hourly)
Hourly$Time = factor(Hourly$Time, levels = c("Ten","Eleven","Twelve","One","Two"))
Hourly$DriftInvertConc<-((Hourly$InvertSample100*100)/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
Hourly$DriftSuckerConc<-(Hourly$Suckers100*100/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
Hourly$DriftSturgeonConc<-(Hourly$LiveSturgeon*100/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)

head(Hourly$Suckers100)




StatSubset<-subset(Hourly, InvertSample100 > 0&DriftInvertConc!= "NA") #Remove sample dates where Inverts were not collected (no dates where drift collections occured which had no drifting invertebrates)

length(StatSubset$ï..SampleID)
length(Hourly$ï..SampleID)

hist(log(Hourly$InvertSample100))
plot(DriftInvertConc~Time,data= StatSubset)

# 
# 
# Model<-aov(log(DriftInvertConc)~Time,data=StatSubset)
# 
# summary(Model)
# plot(Model)
# Tukey<-TukeyHSD(Model,"Time")
# plot(Tukey)
# Tukey
# (exp(1.510104612)-1)*100 #352.7%Difference between 10 pm collection and midnight
# (exp(1.1965172)-1)*100 #230.8574% lower CI
# (exp(1.8236920)-1)*100 #519.4687% upper CI
# 
# difference<-Tukey$Time[,"p adj"]
# Letters<-multcompLetters(difference)
# Letters                         
# #rearrange letters so a is at ten pm, c=a,a=b,c=c
# LetterRearranged<-c("a","b","c","c","c")
# hist(resid(Model))

kruskal.test((DriftInvertConc)~Time,data=StatSubset)

compare_means(DriftInvertConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")


Means=compare_means(DriftInvertConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersInvert<-multcompLetters(difference)
LettersInvert<-LettersInvert$Letters


Trtdata <- ddply(StatSubset, c("Time"), summarise,
                 N    = length(DriftInvertConc),
                 meanInverts = mean(DriftInvertConc),
                 sd   = sd(DriftInvertConc),
                 se   = sd / sqrt(N)
)
Trtdata

HourlyInvertAbu<-ggplot(Trtdata,aes(x=Time,y=meanInverts))+geom_point()+xlab("Collection Time")+
  geom_errorbar(aes(ymin=meanInverts-se,ymax=meanInverts+se))+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanInverts+se+0.5,label=LettersInvert))+
  annotate("text", label = "Kruskal-Wallis,\n Chi-sq = 188.4,\n P < 0.001", size = 3.5, x = 1.5, y = 15)+ theme(axis.title.y = element_text(size = 10))
HourlyInvertAbu


dev.off()
tiff("Figures/Hourly_InvertAbundance.tiff", width = 84, height = 84, units = 'mm', res = 1200)
HourlyInvertAbu
dev.off()


#Sucker by collection
TrtdataSucker <- ddply(StatSubset, c("Time"), summarise,
                       N    = length(DriftSuckerConc),
                       meanSucker = mean(DriftSuckerConc),
                       sd   = sd(DriftSuckerConc),
                       se   = sd / sqrt(N)
)
TrtdataSucker



plot(DriftSuckerConc~Time,data= StatSubset)
kruskal.test((DriftSuckerConc)~Time,data=StatSubset)

compare_means(DriftSuckerConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")


Means=compare_means(DriftSuckerConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersSucker<-multcompLetters(difference)
LettersSucker

HourlySuckerAbu<-ggplot(TrtdataSucker,aes(x=Time,y=meanSucker))+geom_point()+xlab("Collection Time")+
  geom_errorbar(aes(ymin=meanSucker-se,ymax=meanSucker+se))+ylab(expression(Catastomidae~Larvae~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanSucker+se+5,label=LettersSucker$Letters))+
  annotate("text", label = "Kruskal-Wallis,\n Chi-sq = 126.1,\n P < 0.001", size = 3.5, x = 4.5, y = 90)+ theme(axis.title.y = element_text(size = 9))
HourlySuckerAbu


theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/Hourly_SuckerAbundance.tiff", width = 84, height = 84, units = 'mm', res = 1200)
HourlySuckerAbu
dev.off()


#Sturgeon
TrtdataSturgeon <- ddply(StatSubset, c("Time"), summarise,
                         N    = length(DriftSturgeonConc),
                         meanSturgeon = mean(DriftSturgeonConc),
                         sd   = sd(DriftSturgeonConc),
                         se   = sd / sqrt(N)
)
TrtdataSturgeon



plot(DriftSturgeonConc~Time,data= StatSubset)
kruskal.test((DriftSturgeonConc)~Time,data=StatSubset)

compare_means(DriftSturgeonConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")
plot(Model)


Means=compare_means(DriftSturgeonConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")
Means
Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersSturgeon<-multcompLetters(difference)


LettersSturgeon<-c("a","b","c","bc","bc")
TrtdataSturgeon
HourlySturgeonAbu<-ggplot(TrtdataSturgeon,aes(x=Time,y=meanSturgeon))+geom_point()+xlab("Collection Time")+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ylab(expression(Sturgeon~Larvae~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanSturgeon+se+0.5,label=LettersSturgeon))+
  annotate("text", label = "Kruskal-Wallis,\n Chi-sq = 107.7,\n P < 0.001", size = 3.5, x = 2, y = 7.5)+ theme(axis.title.y = element_text(size = 10))
HourlySturgeonAbu


theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/Hourly_SturgeonAbundance.tiff", width = 84, height = 84, units = 'mm', res = 1200)
HourlySturgeonAbu
dev.off()

###############
#Join together fig 3
##############
DOYInvertModelBoot
TotalInvertAbuMoonPhase
HourlyInvertAbu
HourlySuckerAbu
DriftTotalInvertBiomass #From below
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/Fig3July2020.tiff", width = 174, height = 174, units = 'mm', res = 1200)

ggarrange(DriftTotalInvertBiomass,HourlyInvertAbu,HourlySturgeonAbu,HourlySuckerAbu,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)

dev.off()



#########
#Models Invert Biomass standardized by Discharge Sampled
#########

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")
ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)



length(ShannonSubset$SampleID)

head(ShannonSubset)
hist(ShannonSubset$DriftBiomassConc)
hist(log(ShannonSubset$DriftBiomassConc))

ShannonSubset$percillumScale<-scale(ShannonSubset$percillum)
ShannonSubset$TempScale<-scale(ShannonSubset$Temp)
ShannonSubset$DischargeByNightScale<-scale(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSScale<-scale(ShannonSubset$DPFS)


# 
# Biomass0 = lmer(DriftBiomassConc~1+(1|Year),data=ShannonSubset,REML=F)
# Biomass1 = lmer(DriftBiomassConc~percillum+(1|Year),data=ShannonSubset,REML=F)
# Biomass2 = lmer(DriftBiomassConc~temp_centered+(1|Year),data=ShannonSubset,REML=F)
# Biomass4 = lmer(DriftBiomassConc~DPFS+(1|Year),data=ShannonSubset,REML=F)
# 
# Biomass5 = lmer(DriftBiomassConc~percillum+temp_centered+(1|Year),data=ShannonSubset,REML=F)
# Biomass6 = lmer(DriftBiomassConc~percillum+DPFS+(1|Year),data=ShannonSubset,REML=F)
# Biomass10 = lmer(DriftBiomassConc~temp_centered+DPFS+(1|Year),data=ShannonSubset,REML=F)
# 
# Biomass11 = lmer(DriftBiomassConc~temp_centered+DPFS+percillum+(1|Year),data=ShannonSubset,REML=F)
# AICctab(Biomass0,Biomass1,Biomass2,Biomass4,Biomass5,Biomass6,Biomass10,Biomass11, weights=TRUE)
# Log transforming doing a better job for biomass compared to linear model


BiomassLog0 = lme(log(DriftBiomassConc)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

BiomassLog1 = lme(log(DriftBiomassConc)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog2 = lme(log(DriftBiomassConc)~DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog3 = lme(log(DriftBiomassConc)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

BiomassLog4 = lme(log(DriftBiomassConc)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog5 = lme(log(DriftBiomassConc)~percillum+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog6 = lme(log(DriftBiomassConc)~temp_centered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

BiomassLog7 = lme(log(DriftBiomassConc)~percillum+DoYCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog8 = lme(log(DriftBiomassConc)~CTUSturgeon,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
BiomassLog9 = lme(log(DriftBiomassConc)~CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))


AICctab(BiomassLog0,BiomassLog1,BiomassLog2,BiomassLog3,BiomassLog4,BiomassLog5,BiomassLog6,BiomassLog7,BiomassLog8,BiomassLog9, weights=TRUE) 
summary(BiomassLog2)
summary(BiomassLog6)


mean(ShannonSubset$DayOfYear)
intervals(BiomassLog2,which=c("fixed"))
(exp(-0.0462115)-1)*100 #-4.516% estimate DoY
(exp(-0.06238553)-1)*100 #-6.05 CI
(exp(-0.03003747)-1)*100 #-2.96% CI

#Intercept
exp(-0.7472522) #0.47 g/100 m3
exp(-1.08021716) #0.33 g/100 m3
exp(-0.41428726) #0.66 g/ 100 m3




summary(BiomassLog6)

intervals(BiomassLog6,which=c("fixed"))

#Intercept
exp(-0.76390840) #0.4658 estimate intercept
exp(-1.10841480) # 0.330 CI
exp(-0.41940201) # 0.657 CI

(exp(-0.03869282)-1)*100 #-3.795 estimate DoY
(exp(-0.05861655)-1)*100 #-5.69% CI
(exp(-0.01876908)-1)*100 #-1.86% CI

intervals(BiomassLog6,which=c("fixed"))
(exp(-0.06113747)-1)*100 #-5.93 estimate temp
(exp(-0.15521462)-1)*100 #-14.38% CI
(exp(0.03293969)-1)*100 #3.35% CI

summary(BiomassLog3) #t=-3.40
intervals(BiomassLog3,which=c("fixed"))
(exp(-0.1502938)-1)*100 #-13.954% estimate temp only
(exp(-0.06363162)-1)*100 #-6.165 CI
(exp(-0.236956)-1)*100 #-21.097% CI


cor(ShannonSubset$DriftInvertConc,ShannonSubset$DriftBiomassConc) #Drift biomass and total inverts are highly correlated 0.9647
cor(ShannonSubset$temp_centered,ShannonSubset$DoYCentered) #Temp and discharge are relatively highlty correlated (0.4933)
######################
#Biomass LMM plot
#####################

#Bootstrapped CIs rom https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)


ShannonSubset$Year<-as.factor(ShannonSubset$Year)
ShannonSubset$DriftLog<-log(ShannonSubset$DriftBiomassConc)
BiomassLog2NoCor = lme((DriftLog)~DoYCentered,random=~1|Year,data=ShannonSubset)#correlation=corAR1(form = ~DayOfYear|Year)

xvals <-  with(ShannonSubset,seq(min(DoYCentered),max(DoYCentered),length.out=100))
nresamp<-1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked <- mvrnorm(nresamp, mu = fixef(BiomassLog2NoCor), Sigma = vcov(BiomassLog2NoCor))


## predicted values: useful below
pframe <- with(ShannonSubset,data.frame(DoYCentered))
pframe$DriftLog <- predict(BiomassLog2NoCor,newdata=pframe,level=0)
head(pframe)
## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

## bootstrapping
sampfun <- function(fitted,data,idvar="Year") {
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted)
  dd <- data.frame(data,pred=pp,res=rr)
  ## sample groups with replacement
  iv <- levels(data[[idvar]])
  bsamp1 <- sample(iv,size=length(iv),replace=TRUE)
  bsamp2 <- lapply(bsamp1,
                   function(x) {
                     ## within groups, sample *residuals* with replacement
                     ddb <- dd[dd[[idvar]]==x,]
                     ## bootstrapped response = pred + bootstrapped residual
                     ddb$DriftLog <- ddb$pred +
                       sample(ddb$res,size=nrow(ddb),replace=TRUE)
                     return(ddb)
                   })
  res <- do.call(rbind,bsamp2)  ## collect results
  if (is(data,"groupedData"))
    res <- groupedData(res,formula=formula(data))
  return(res)
}

pfun <- function(fm) {
  predict(fm,newdata=pframe,level=0)
}

set.seed(101)
yvals2 <- replicate(nresamp,
                    pfun(update(BiomassLog2NoCor,data=sampfun(BiomassLog2NoCor,ShannonSubset,"Year"))))
c2 <- get_CI(yvals2,"boot_")
pframe2 <- data.frame(pframe,c2)
head(pframe2)
pframe2$DayOfYear<-pframe2$DoYCentered+mean(ShannonSubset$DayOfYear) #Remove centering on DOY for plotting
pframe2$DriftBiomassConc<-exp(pframe2$DriftLog)

pframe2$boot_lwr<-exp(pframe2$boot_lwr)
pframe2$boot_upr<-exp(pframe2$boot_upr)




DOYBiomassModelBoot<-ggplot(ShannonSubset,aes(DayOfYear,DriftBiomassConc))+geom_point(color="darkgrey")+geom_line(data=pframe2,size=1.5)+  
  geom_line(data = pframe2, aes(y = boot_lwr), size = .75,linetype="dashed")+geom_line(data = pframe2,aes(y=boot_upr),size=0.75,linetype="dashed")+
  xlab("Day of Year")+ylab(expression(Macroinvertebrate~Biomass~(g)~Per~100~m^3~Drift))+ theme(axis.title.y = element_text(size = 8))
DOYBiomassModelBoot

dev.off()
tiff("Figures/DOYBiomassModelBoot.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DOYBiomassModelBoot
dev.off()












####################
#Total invert biomass By Moon Phase
#####################
head(ShannonRichness)
TotalInvertBiomass<-subset(ShannonRichness, DriftBiomassConc!="NA")
head(TotalInvertBiomass)
Trtdata <- ddply(TotalInvertBiomass, c("MoonPhase"), summarise,
                 N    = length(DriftBiomassConc),
                 meanSturgeon = mean(DriftBiomassConc),
                 sd   = sd(DriftBiomassConc),
                 se   = sd / sqrt(N)
)
head(Trtdata)
Trtdata
Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))


kruskal.test(DriftBiomassConc~MoonPhase, data=ShannonRichness)
ShannonRichness$MoonPhase = factor(ShannonRichness$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

compare_means(DriftBiomassConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Means=compare_means(DriftBiomassConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
Letters<-multcompLetters(difference)
Letters<-Letters$Letters
Trtdata
vector<-c("bc","bc","abc","a","ab","b","ab","c") #Relabeled due to new moon not being a
DriftTotalInvertBiomass<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(fill=MoonPhase),stat="identity")+xlab("Moon Phase")+ylab(expression(Macroinvertebrate~Biomass~(g)~Per~100~m^3~Drift~(SE)))+#Invertebrate Biomass (g) per 100 m3 drift (SEM)
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 8))+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+.08,label=vector))+geom_text(aes(x=5,y=1.5,label= "KW, chi-squared = 49.3, P < 0.001"),size=3)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
DriftTotalInvertBiomass
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/BiomassByPhase100m3.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
DriftTotalInvertBiomass
dev.off()








####################
#Invertebrate family richness
####################

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")
ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

RichnessLog0A = lme(log(Nfamilies)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

RichnessLog1A = lme(log(Nfamilies)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog2A = lme(log(Nfamilies)~DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog3A = lme(log(Nfamilies)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

RichnessLog4A = lme(log(Nfamilies)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog5A = lme(log(Nfamilies)~percillum+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog6A = lme(log(Nfamilies)~temp_centered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

RichnessLog7A = lme(log(Nfamilies)~percillum+DoYCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog8A = lme(log(Nfamilies)~DischargeCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog9A = lme(log(Nfamilies)~DischargeCentered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog10A = lme(log(Nfamilies)~DischargeCentered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog11A = lme(log(Nfamilies)~DischargeCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog12A = lme(log(Nfamilies)~DischargeCentered+CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog13A = lme(log(Nfamilies)~CTUSturgeon,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
RichnessLog14A = lme(log(Nfamilies)~CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))


AICctab(RichnessLog0A,RichnessLog1A,RichnessLog2A,RichnessLog3A,RichnessLog4A,RichnessLog5A,RichnessLog6A,RichnessLog7A,RichnessLog8A,
        RichnessLog9A,RichnessLog10A,RichnessLog11A,RichnessLog12A,RichnessLog13A,RichnessLog14A,weights=TRUE)


summary(RichnessLog2A)
summary(RichnessLog9A)
intervals(RichnessLog2A,which=c("fixed"))

#DOY estimates
(exp(-0.01001227)-1)*100 #-0.9962% estimate 
(exp(-0.01488792)-1)*100 #-1.477% estimate 
(exp(-0.00513661)-1)*100 #-0.512% estimate 



#Day 140
min(ShannonSubset$DayOfYear)
max(ShannonSubset$DayOfYear)
mean(ShannonSubset$DayOfYear) #mean =158.4 to account for centering 

exp(2.1725554+(-18.4*-0.01001227)) #N families Estimate at 140 days 10.556
exp(2.1725554+(-18.4*-0.00513661)) #Upper CI at 15 days 9.65
exp(2.1725554+(-18.4*-0.01488792)) #Lower CI at 15 days 11.55


#Day 170
170-158.4
exp(2.1725554+(11.6*-0.01001227)) #N families Estimate at 170 days 7.82
exp(2.1725554+(11.6*-0.00513661)) #Upper CI at 170 days 8.27
exp(2.1725554+(11.6*-0.01488792)) #Lower CI at 170 days 7.39

###############
#LMM Richness plot Day of Year
#################
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)
ShannonSubset$Qcentered = ShannonSubset$Q - mean(ShannonSubset$Q)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)


ShannonSubset$Year<-as.factor(ShannonSubset$Year)
ShannonSubset$LogFamilies<-log(ShannonSubset$Nfamilies)
RichnessLog2A = lme(LogFamilies~DoYCentered,random=~1|Year,data=ShannonSubset)


xvals <-  with(ShannonSubset,seq(min(DoYCentered),max(DoYCentered),length.out=100))
nresamp<-1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked <- mvrnorm(nresamp, mu = fixef(RichnessLog2A), Sigma = vcov(RichnessLog2A))


## predicted values: useful below
pframe <- with(ShannonSubset,data.frame(DoYCentered))
pframe$LogFamilies <- predict(RichnessLog2A,newdata=pframe,level=0)
head(pframe)
## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

## bootstrapping
sampfun <- function(fitted,data,idvar="Year") {
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted)
  dd <- data.frame(data,pred=pp,res=rr)
  ## sample groups with replacement
  iv <- levels(data[[idvar]])
  bsamp1 <- sample(iv,size=length(iv),replace=TRUE)
  bsamp2 <- lapply(bsamp1,
                   function(x) {
                     ## within groups, sample *residuals* with replacement
                     ddb <- dd[dd[[idvar]]==x,]
                     ## bootstrapped response = pred + bootstrapped residual
                     ddb$LogFamilies <- ddb$pred +
                       sample(ddb$res,size=nrow(ddb),replace=TRUE)
                     return(ddb)
                   })
  res <- do.call(rbind,bsamp2)  ## collect results
  if (is(data,"groupedData"))
    res <- groupedData(res,formula=formula(data))
  return(res)
}

pfun <- function(fm) {
  predict(fm,newdata=pframe,level=0)
}

set.seed(101)
yvals2 <- replicate(nresamp,
                    pfun(update(RichnessLog2A,data=sampfun(RichnessLog2A,ShannonSubset,"Year"))))
c2 <- get_CI(yvals2,"boot_")
pframe2 <- data.frame(pframe,c2)
head(pframe2)
pframe2$DayOfYear<-pframe2$DoYCentered+mean(ShannonSubset$DayOfYear) #Remove centering on DOY for plotting
pframe2$Nfamilies<-exp(pframe2$LogFamilies)

pframe2$boot_lwr<-exp(pframe2$boot_lwr)
pframe2$boot_upr<-exp(pframe2$boot_upr)




DOYFamilyModelBoot<-ggplot(ShannonSubset,aes(DayOfYear,Nfamilies))+geom_point(color="darkgrey")+geom_line(data=pframe2,size=1.5)+  
  geom_line(data = pframe2, aes(y = boot_lwr), size = .75,linetype="dashed")+geom_line(data = pframe2,aes(y=boot_upr),size=0.75,linetype="dashed")+
  xlab("Day of Year")+ylab("Macroinvertebrate Family Richness")+ theme(axis.title.y = element_text(size = 8))
DOYFamilyModelBoot

dev.off()
tiff("Figures/DOYFamilyModelBoot.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DOYFamilyModelBoot
dev.off()



############
#LMM models Shannon
#############
head(ShannonSubset)



ShannonLog0A = lme(log(value)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

ShannonLog1A = lme(log(value)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog2A = lme(log(value)~DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog3A = lme(log(value)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

ShannonLog4A = lme(log(value)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog5A = lme(log(value)~percillum+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog6A = lme(log(value)~temp_centered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

ShannonLog7A = lme(log(value)~percillum+DoYCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog8A = lme(log(value)~DischargeCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog9A = lme(log(value)~DischargeCentered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog10A = lme(log(value)~DischargeCentered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog11A = lme(log(value)~DischargeCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog12A = lme(log(value)~DischargeCentered+temp_centered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog13A = lme(log(value)~DischargeCentered+CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog14A = lme(log(value)~CTUSturgeon,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
ShannonLog15A = lme(log(value)~CTUSturgeon+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))


AICctab(ShannonLog0A,ShannonLog1A,ShannonLog2A,ShannonLog3A,
        ShannonLog4A,ShannonLog5A,ShannonLog6A,ShannonLog7A,ShannonLog8A,
        ShannonLog9A,ShannonLog10A,ShannonLog11A,ShannonLog12A,ShannonLog13A,ShannonLog14A,ShannonLog15A, weights=TRUE)


summary(ShannonLog8A)
summary(ShannonLog0A)


intervals(ShannonLog8A,which=c("fixed"))
(exp(0.4891075)-1)*100 #63.086% estimate 
(exp(0.15330047)-1)*100 #16.57% CI 
(exp(0.8249146)-1)*100 #128.17% CI

#############
#Diagnostic plots Shannon LMM
############




plot(resid(ShannonLog8A) ~ temp_centered, data=ShannonSubset)
lines(lowess(resid(ShannonLog8A) ~ ShannonSubset$temp_centered), col=2)


hist(resid(ShannonLog8A))
plot(resid(ShannonLog8A) ~ DischargeCentered, data=ShannonSubset)
lines(lowess(resid(ShannonLog8A) ~ ShannonSubset$DischargeCentered), col=2)

plot(resid(ShannonLog8A) ~ predict(ShannonLog8A))
lines(lowess(resid(ShannonLog8A) ~ predict(ShannonLog8A)), col=2)

summary(ShannonLog8A)
confint(ShannonLog8A)
fixef(ShannonLog8A)
ranef(ShannonLog8A)
confint(ShannonLog8A)
coef(ShannonLog8A)




#############
#Envfit plot
#############

physeqSubset<- subset_samples(physeq, percillum!="NA")

physeqSubset
physeqSubset<- subset_samples(physeqSubset, DischargeSampledByNight !="NA") #Remove samples without discharge info
physeqSubset<- subset_samples(physeqSubset, Q!="NA")


physeqSubset<- subset_samples(physeqSubset, DischargeSampledByNight < 2) #Remove discharge outliers

physeqSubset
GPdist=phyloseq::distance(physeqSubset, "bray")

vare.mds= ordinate(physeqSubset, "NMDS",GPdist)
#vare.mds <- metaMDS(VeganDist, trace = FALSE)
#vare.mds
cor(ShannonSubset$DayOfYear,ShannonSubset$CTUSturgeon)
metadataEnvfitSubset<-subset(sample_data(physeqSubset),DischargeSampledByNight !="NA"&DischargeSampledByNight < 2&Q!="NA",percillum!="NA")
EnvFitMeta=data.frame(metadataEnvfitSubset$percillum,metadataEnvfitSubset$DayOfYear,metadataEnvfitSubset$DischargeSampledByNight,metadataEnvfitSubset$Q) #metadataEnvfitSubset$CTUSturgeon to show ctu as well

head(EnvFitMeta)
colnames(EnvFitMeta)<-c("Lunar\n Illumination (%)","DOY","\nDischarge\n Sampled"," River\n  Discharge")#"CTU"
#EnvFitMeta <- EnvFitMeta[!is.na(EnvFitMeta)]
#EnvFitMeta
ef =envfit(vare.mds, EnvFitMeta, na.rm=TRUE, permu=999)

#cor(metadataEnvfitSubset$DischargeSampledByNight,metadataEnvfitSubset$Q)
ef
plot(vare.mds,display="sites")
envplot=plot(ef, p.max = 0.05)
#ordisurf(ord ~ A1, data = dune.env, add = TRUE, knots = 1)
envplot

dev.off()
tiff("Figures/EnvfitPlot.tiff", width = 84, height = 84, units = 'mm', res = 1200)
par(mar = c(4, 4, 0.5, 0.5))
plot(vare.mds,display="sites",cex=0.5)
plot(ef, p.max = 0.05,cex=0.6)

#mtext("NMDS2", side = 2, line = 1, cex = 1)

envplot=plot(ef, p.max = 0.05,cex=0.6)
dev.off()

#Save the plot as a single object 
colnames(EnvFitMeta)<-c("","","","")
ef =envfit(vare.mds, EnvFitMeta, na.rm=TRUE, permu=999)

dev.off()
pdf(NULL)
dev.control(displaylist = "enable")
par(mar = c(4, 4, 0.5, 0.5))
plot(vare.mds,display="sites",cex=0.5)
plot(ef, p.max = 0.05,cex=0.7)
p1<-recordPlot()
dev.off()
p1

############
#Join Together Figure 4
############

DOYFamilyModelBoot
DriftTotalInvertBiomass
DOYBiomassModelBoot
p1


library(gridGraphics)

dev.off()
tiff("Figures/Figure4July.tiff", width = 84, height = 174, units = 'mm', res = 1200)
ggarrange(DOYFamilyModelBoot,p1,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
dev.off()


############
#Total Abundance Family by Moon phase
############
DriftDataCombined<-read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\SturgeonDrift\\SturgeonDrift\\DataClean\\AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)
head(DriftDataCombined)



#Isonychiidae
kruskal.test(Isonychiidae~MoonPhase,data=DriftDataCombined)
Means<-compare_means(Isonychiidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means


Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersIso<-multcompLetters(difference)
LettersIso
#Rearrange letters so New moon letter = a instead of d
#d=a, c=b,a=c, b=d


#Heptageniidae
kruskal.test(Heptageniidae~MoonPhase,data=DriftDataCombined)
Means<-compare_means(Heptageniidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means


Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHep<-multcompLetters(difference)
LettersHep


#Crayfish
kruskal.test(Crayfish~MoonPhase,data=DriftDataCombined)
Means<-compare_means(Crayfish ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means


Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersCray<-multcompLetters(difference)
LettersCray

#Hydropsychiidae
kruskal.test(Hydropsychiidae~MoonPhase,data=DriftDataCombined)
Means<-compare_means(Hydropsychiidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means


Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHyd<-multcompLetters(difference)
LettersHyd

#Graph
physeqSubset<-subset_taxa(physeq,Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")
physeqSubset
df <- psmelt(physeqSubset)

ddply(df, c("Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)




LettersCray
LettersHep
LettersIso
#Rearrange iso letters so New moon letter = a instead of d
#d=a, c=b,a=c, b=d
TrtdataSorted<-Trtdata[order(Trtdata$Family),]
head(TrtdataSorted)




df$Abundance<-df$Abundance*20 #To get total number based on 5% subsample ID
head(df)
df$AbuPer100<-((df$Abundance*100)/(60*4*60*df$AreaSampled.m2.*df$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
df<-subset(df,AbuPer100!="NA")

head(df)
levels(df$Family)
FamiliesPer100m3<-compare_means(AbuPer100 ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")
FamiliesPer100m3

Trtdata <- ddply(df, c("MoonPhase","Family"), summarise,
                 N    = length(AbuPer100),
                 mean = mean(AbuPer100),
                 sd   = sd(AbuPer100),
                 se   = sd / sqrt(N)
)

Trtdata

Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

LettersCray
LettersHep
LettersIso
#Rearrange iso letters so New moon letter = a instead of d
#d=a, c=b,a=c, b=d
TrtdataSorted<-Trtdata[order(Trtdata$Family),]
head(TrtdataSorted)
vector<-c("","","","","","","","",
          "","","","","","","","",
          "a","a","ab","b","ab","a","a","a",
          "","","","","","","","",
          "a","ab","cd","d","bcd ","bc","abc   ","abc ",
          "","","","","","","","")


length(vector)
length(TrtdataSorted$N)
KruskalLabel<- c("Kruskal-Wallis,\n P-adj = 0.053","KW, P-adj = 0.89", "KW, P-adj < 0.001","KW, P-adj = 0.01","     P-adj < 0.001","KW, P-adj = 0.17")

dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.053","KW, P-adj = 0.89", "  KW, P-adj < 0.001"," KW, P-adj = 0.01","     P-adj < 0.001","KW, P-adj = 0.17"),
  Family   = c("Crayfish","Ephemerillidae","Heptageniidae","Hydropsychidae","Isonychiidae","Leptoceridae")
)
#TrtdataSorted

TopFamilyAbuPer100=ggplot(TrtdataSorted, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = MoonPhase),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Macroinvertebrates~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+geom_text(aes(x=MoonPhase,y= mean+se+1.5),label=vector,size=2)
TopFamilyAbuPer100
TopFamilyAbuPer100<-TopFamilyAbuPer100+ geom_text(data=dat_text,size=2.5,mapping = aes(x = 4, y = 15, label = label))
TopFamilyAbuPer100

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FamilyLevelAbuPer100m3.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TopFamilyAbuPer100
dev.off()

###############
#All families biomass
##############

BiomassAverages<-read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\SturgeonDrift\\SturgeonDrift\\DataClean\\BiomassFamilyAveragesR.txt",header=T)
head(BiomassAverages)
DriftDataCombined<-read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\SturgeonDrift\\SturgeonDrift\\DataClean\\AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)

otubiomass<-otufull

row.names(otubiomass)<-taxmatrixfull
#Top famlies Rel Abu above 3% <-c("Chironomidae","Crayfish","Ephemerillidae","Heptageniidae","Hydropsychidae","Isonychiidae","Leptoceridae")

otubiomasst<-t(otubiomass)
otubiomasst<-as.data.frame(otubiomasst)
otubiomasstAbove20 <- otubiomasst[,colSums(otubiomasst) > 20] #Remove families with Total N collected < 20
head(otubiomasstAbove20)
colnames(otubiomasstAbove20)[1]


for (i in colnames(otubiomasstAbove20)){
  Name<-factor(i)
  otubiomasstAbove20[i]<-otubiomasstAbove20[i]*sum(BiomassAverages[i])
  
}
otubiomasstAbove20[i]<-otubiomasstAbove20[i]*sum(BiomassAverages[i])

head(otubiomasstAbove20)


OTUBiomassObject<-as.matrix(t(otubiomasstAbove20))
OTUBiomassObject

TaxMatrixAbove20N<-as.matrix(colnames(otubiomasstAbove20))
TaxMatrixAbove20N


OTUBiomass=otu_table(OTUBiomassObject, taxa_are_rows=TRUE)
#OTU
TAXBiomass=tax_table(TaxMatrixAbove20N)
colnames(TAXBiomass)=("Family")

sampdat=sample_data(metadata)
sample_names(sampdat)=sampdat$SampleID
#head(sampdat)
taxa_names(TAXBiomass)=TaxMatrixAbove20N
physeqBiomass=phyloseq(OTUBiomass,TAXBiomass,sampdat)#joins together OTU,TAX, and metadata 
physeqBiomass
sample_data(physeqBiomass)$MoonPhase = factor(sample_data(physeqBiomass)$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))


df <- psmelt(physeqBiomass)
df$Abundance<-df$Abundance*20# To account for 5% subsample used for ID
df$BiomassPer100<-((df$Abundance*100)/(60*4*60*df$AreaSampled.m2.*df$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
head(df)
df$PercentageTotalBiomass<-df$BiomassPer100/df$DriftBiomassConc*100 #Percentage of total biomass for each family

dfSubset<-subset(df,BiomassPer100!="NA")

BiomassPer100m3<-compare_means(BiomassPer100 ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")
BiomassPer100m3

subset(BiomassPer100m3, p.adj < 0.05) #Isonychiidae, Heptageniidae, Hydropsychidae



Trtdata <- ddply(dfSubset, c("MoonPhase","Family"), summarise,
                 N    = length(BiomassPer100),
                 mean = mean(BiomassPer100),
                 sd   = sd(BiomassPer100),
                 se   = sd / sqrt(N),
                 meanPercent = mean(PercentageTotalBiomass),
                 sdPercent   = sd(PercentageTotalBiomass),
                 sePercent   = sdPercent / sqrt(N)
)
Trtdata
options(scipen = 999)
Trtdata



AllFamilyBiomass=ggplot(Trtdata, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Biomass~(g)~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.,scales="free_y")+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))#+ annotate("text", label = KruskalLabel, size = 2, x = 4, y = .7)+scale_fill_manual(values=cbPalette)
AllFamilyBiomass
theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/AllFamilyBiomass.tiff", width = 174, height = 174, units = 'mm', res = 1200)
AllFamilyBiomass
dev.off()


########
#Biomass top families only
###########
#See section above for data/table formatting

#
Trtdata <- ddply(dfSubset, c("Family"), summarise,
                 N    = length(BiomassPer100),
                 mean = mean(BiomassPer100),
                 sd   = sd(BiomassPer100),
                 se   = sd / sqrt(N),
                 meanPercent = mean(PercentageTotalBiomass),
                 sdPercent   = sd(PercentageTotalBiomass),
                 sePercent   = sdPercent / sqrt(N)
)

TrtdataSorted<-Trtdata[order(-Trtdata$meanPercent),]
TrtdataSorted

Top6BiomassFam<-(TrtdataSorted)[1:6,]
sum(Top6BiomassFam$meanPercent) #81% of total biomass comes from top six families
Top6BiomassFam
BiomassPer100m3<-compare_means(BiomassPer100 ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")
#BiomassPer100m3
subset(BiomassPer100m3,Family %in% Top6BiomassFam$Family)
subset(BiomassPer100m3, p.adj < 0.05) #Isonychiidae, Heptageniidae, Hydropsychidae


Trtdata <- ddply(dfSubset, c("Family","MoonPhase"), summarise,
                 N    = length(BiomassPer100),
                 mean = mean(BiomassPer100),
                 sd   = sd(BiomassPer100),
                 se   = sd / sqrt(N))
TrtdataSubset<-subset(Trtdata,Family %in% Top6BiomassFam$Family)
TrtdataSubset


#Heptageniidae
kruskal.test(Heptageniidae~MoonPhase,data=DriftDataCombined)# Use adj-p values
Means<-compare_means(Heptageniidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHep<-multcompLetters(difference)
LettersHep

#Isonychiidae
kruskal.test(Isonychiidae~MoonPhase,data=DriftDataCombined) #Use adj-p values
Means<-compare_means(Isonychiidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersIso<-multcompLetters(difference)
LettersIso

TrtdataSubset<-TrtdataSubset[order(Trtdata$Family),]
TrtdataSubset<-subset(TrtdataSubset,Family!="NA")


LettersHep
LettersIso
# vector<-c("","","","","","","","",
#           "","","","","","","","",
#           "","","","","","","","",
#           "a","a","ab","b","ab","a","a","a",
#           "a","ab","cd","d","bcd","bc","abc","abc",
#           "","","","","","","","")
# 
vector<-c("a","a","ab","b","ab","a","a","a","a","ab","cd","d","bcd","bc","abc","abc")
dat_text <- data.frame(
  label = c("KW, Chi-sq = 50.19,\n adj-P < 0.001",  "KW, Chi-sq = 49.56,\n      adj-P < 0.001"),
  Family   = c("Heptageniidae","Isonychiidae")
)

TrtdataSubsetHI<-subset(TrtdataSubset, Family=="Heptageniidae"|Family=="Isonychiidae")
TrtdataSubsetHI
HepIsoBiomass=ggplot(TrtdataSubsetHI, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = MoonPhase),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Biomass~(g)~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+ geom_text(data=dat_text,size=3,mapping = aes(x = 5, y = 0.9, label = label))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase,y= mean+se+.05),label=vector,size=2.5)
HepIsoBiomass


vector<-c("","","","","","","","",
          "","","","","","","","",
          "","","","","","","","",
          "a","a","ab","b","ab","a","a","a",
          "a","ab","cd","d","bcd","bc","abc    ","abc ",
          "","","","","","","","")


dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.2","KW, P-adj = 0.89","KW, P-adj = 0.55", "    KW, P-adj < 0.001",  "      P-adj < 0.001","KW, P-adj = 0.89"),
  Family   = c("Crayfish","Ephemerillidae","Gomphidae","Heptageniidae","Isonychiidae","Lepidostomatidae")
)
TopFamilyBiomass=ggplot(TrtdataSubset, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Biomass~(g)~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+ geom_text(data=dat_text,size=2.4,mapping = aes(x = 4, y = 0.9, label = label))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase,y= mean+se+.07),label=vector,size=2)
TopFamilyBiomass


theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TopFamiliesFamiliesBiomass.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TopFamilyBiomass
dev.off()

cor(ShannonSubset$DPFS,ShannonSubset$percillum)

#################
#LMM models top families
############
FamilyLMMData<-psmelt(physeq)
head(FamilyLMMData)

FamilyLMMData<-subset(FamilyLMMData, Temp!= "NA")
FamilyLMMData<-subset(FamilyLMMData, AverageNetFlowByNight!= "NA")
FamilyLMMData<-subset(FamilyLMMData,FamilyLMMData$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
FamilyLMMData$temp_centered = FamilyLMMData$Temp - mean(FamilyLMMData$Temp)
FamilyLMMData$Qcentered = FamilyLMMData$Q - mean(FamilyLMMData$Q)
FamilyLMMData$AverageFlowCentered = FamilyLMMData$AverageNetFlowByNight - mean(FamilyLMMData$AverageNetFlowByNight)
FamilyLMMData$DischargeCentered = FamilyLMMData$DischargeSampledByNight - mean(FamilyLMMData$DischargeSampledByNight)
FamilyLMMData$DPFSCentered = FamilyLMMData$DPFS - mean(FamilyLMMData$DPFS)
FamilyLMMData$DOYCentered = FamilyLMMData$DayOfYear - mean(FamilyLMMData$DayOfYear)

head(FamilyLMMData)


#########
#Crayfish model
#########
CrayfishAbu<-subset(FamilyLMMData,Family=="Crayfish")
head(CrayfishAbu)
CrayfishAbu$Abundance<-CrayfishAbu$Abundance*20 #Account for 5% sampling
CrayfishAbu$CrayfishPer100<-((CrayfishAbu$Abundance*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
CrayfishAbu$SturgeonPer100<-((CrayfishAbu$Nsturgeon*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
CrayfishAbu$SuckerPer100<-((CrayfishAbu$Nsuckers100*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)


hist(CrayfishAbu$CrayfishPer100)
hist(log(CrayfishAbu$CrayfishPer100))
CrayfishAbu$percillumScale<-scale(CrayfishAbu$percillum)
CrayfishAbu$TempScale<-scale(CrayfishAbu$Temp)
CrayfishAbu$DischargeByNightScale<-scale(CrayfishAbu$DischargeSampledByNight)
CrayfishAbu$DPFSScale<-scale(CrayfishAbu$DPFS)


head(CrayfishAbu)
CrayfishLog0 = lme(log(CrayfishPer100+10e-5)~1,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))

CrayfishLog1 = lme(log(CrayfishPer100+10e-5)~percillum,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))

CrayfishLog2 = lme(log(CrayfishPer100+10e-5)~DOYCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))
# CrayfishLog2Poly = lme((CrayfishPer100)~poly(DPFS,2),random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
# summary(CrayfishLog2Poly)



CrayfishLog3 = lme(log(CrayfishPer100+10e-5)~temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))

CrayfishLog4 = lme(log(CrayfishPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))
CrayfishLog5 = lme(log(CrayfishPer100+10e-5)~percillum+DOYCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))
CrayfishLog6 = lme(log(CrayfishPer100+10e-5)~temp_centered+DOYCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))

CrayfishLog7 = lme(log(CrayfishPer100+10e-5)~percillum+DOYCentered+temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DayOfYear|Year))
AICctab(CrayfishLog0,CrayfishLog1,CrayfishLog2,CrayfishLog3,CrayfishLog4,CrayfishLog5,CrayfishLog6,CrayfishLog7, weights=TRUE)
hist(resid(CrayfishLog2))
plot(resid(CrayfishLog2) ~ DOYCentered, data=CrayfishAbu)
lines(lowess(resid(CrayfishLog2) ~ CrayfishAbu$DOYCentered), col=2)
summary(CrayfishLog2)


#


CrayfishAbuLog0 <- gam(log(CrayfishPer100+10e-5)~ 1+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

CrayfishAbuLog1 <- gam(log(CrayfishPer100+10e-5)~ s(DayOfYear)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog2 <- gam((CrayfishPer100+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog3 <- gam(log(CrayfishPer100+10e-5)~ s(DayOfYear)+s(Year,bs="re")+s(temp_centered), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog4 <- gam(log(CrayfishPer100+10e-5)~ s(DayOfYear)+s(Year,bs="re")+s(percillum), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog5 <- gam(log(CrayfishPer100+10e-5)~ s(CTUSturgeon)+s(Year,bs="re")+s(percillum), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
AICctab(CrayfishAbuLog0,CrayfishAbuLog1,CrayfishAbuLog2,CrayfishAbuLog3,CrayfishAbuLog4,CrayfishAbuLog5,weights=T)
# 
appraise(CrayfishAbuLog2)

#Gamma models performing better by AIC but residual fit still off
# #NCrayfishs100
# #CrayfishConc
CrayfishAbuGamma0 <- gam((CrayfishPer100+10e-5)~ 1+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
CrayfishAbuGamma1 <- gam((CrayfishPer100+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
CrayfishAbuGamma2 <- gam((CrayfishPer100+10e-5)~ s(CTUSturgeon)+s(percillum)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
CrayfishAbuGamma3 <- gam((CrayfishPer100+10e-5)~ s(DayOfYear)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
CrayfishAbuGamma4 <- gam((CrayfishPer100+10e-5)~ s(DayOfYear)+s(temp_centered)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
CrayfishAbuGamma5 <- gam((CrayfishPer100+10e-5)~ s(DayOfYear)+s(Year,bs="re")+s(percillum), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
# 
AICctab(CrayfishAbuGamma0,CrayfishAbuGamma1,CrayfishAbuGamma2,CrayfishAbuGamma3,CrayfishAbuGamma4,CrayfishAbuGamma5,weights=T)
AICctab(CrayfishAbuGamma2,CrayfishAbuLog2,CrayfishLog2, weights=T)
# summary(CrayfishAbuGamma2)
CrayfishAbuLog2$aic
CrayfishAbuGamma2$aic
appraise(CrayfishAbuGamma2)
summary(CrayfishAbuGamma2)
summary(CrayfishAbuGamma1)
##########33
#Top model
#########3
#Note: Both Gamma and log based GAM models performed poorly, but models show a very strong effect of CTU on crayfish density
plot<-plot_smooth(CrayfishAbuGamma2, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
hist(resid(CrayfishAbuLog2))
hist(resid(CrayfishLog2))
summary(CrayfishLog2)
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("CrayfishPer100","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=CrayfishPer100))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")
CrayfishNoModelPlot<-ggplot(CrayfishAbu, aes(x=CTUSturgeon,y=CrayfishPer100))+geom_point(color="black")+ylab(expression(Crayfish~Per~100~m^3~Drift))+xlab("Cumulative Thermal Units")#+facet_wrap(~Year)
ggplot(CrayfishAbu, aes(x=DayOfYear,y=CrayfishPer100))+geom_point(color="black")+ylab(expression(Crayfish~Per~100~m^3~Drift))+xlab("Day of Year")#+facet_wrap(~Year)


CrayfishGAMPlot<-ggplot(CrayfishAbu,aes(CTUSturgeon,CrayfishPer100))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Crayfish~Per~100~m^3~Drift))
CrayfishGAMPlot
appraise(CrayfishAbuGamma2)
summary(CrayfishAbuLog2)
appraise(CrayfishAbuLog2)
summary(CrayfishAbuLog5)

dev.off()
tiff("Figures/CrayfishGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
CrayfishGAMPlot
dev.off()

dev.off()
tiff("Figures/CrayfishNoModel.tiff", width = 84, height = 84, units = 'mm', res = 1200)
CrayfishNoModelPlot
dev.off()

######
#Isonychiidae abundance model
########
FamilyLMMData<-psmelt(physeq)
head(FamilyLMMData)

FamilyLMMData<-subset(FamilyLMMData, Temp!= "NA")
FamilyLMMData<-subset(FamilyLMMData, AverageNetFlowByNight!= "NA")
FamilyLMMData<-subset(FamilyLMMData,FamilyLMMData$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
FamilyLMMData$temp_centered = FamilyLMMData$Temp - mean(FamilyLMMData$Temp)
FamilyLMMData$Qcentered = FamilyLMMData$Q - mean(FamilyLMMData$Q)
FamilyLMMData$AverageFlowCentered = FamilyLMMData$AverageNetFlowByNight - mean(FamilyLMMData$AverageNetFlowByNight)
FamilyLMMData$DischargeCentered = FamilyLMMData$DischargeSampledByNight - mean(FamilyLMMData$DischargeSampledByNight)
FamilyLMMData$DPFSCentered = FamilyLMMData$DPFS - mean(FamilyLMMData$DPFS)
FamilyLMMData$DOYCentered = FamilyLMMData$DayOfYear - mean(FamilyLMMData$DayOfYear)

head(FamilyLMMData)
IsoAbu<-subset(FamilyLMMData,Family=="Isonychiidae")
head(IsoAbu)
IsoAbu$Abundance<-IsoAbu$Abundance*20 #Account for 5% subsample

hist(IsoAbu$Abundance)
hist(log(IsoAbu$Abundance))
IsoAbu$IsoPer100<-((IsoAbu$Abundance*100)/(60*4*60*IsoAbu$AreaSampled.m2.*IsoAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
hist(IsoAbu$IsoPer100)
hist(log(IsoAbu$IsoPer100))


IsoLog0 = lme(log(IsoPer100+10e-5)~1,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))

IsoLog1 = lme(log(IsoPer100+10e-5)~percillum,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog2 = lme(log(IsoPer100+10e-5)~DOYCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog3 = lme(log(IsoPer100+10e-5)~temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))

IsoLog4 = lme(log(IsoPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog5 = lme(log(IsoPer100+10e-5)~percillum+DOYCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog6 = lme(log(IsoPer100+10e-5)~temp_centered+DOYCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))

IsoLog7 = lme(log(IsoPer100+10e-5)~percillum+DOYCentered+temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog8 = lme(log(IsoPer100+10e-5)~CTUSturgeon,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog9 = lme(log(IsoPer100+10e-5)~CTUSturgeon+percillum,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))
IsoLog10 = lme(log(IsoPer100+10e-5)~CTUSturgeon,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DayOfYear|Year))


AICctab(IsoLog0,IsoLog1,IsoLog2,IsoLog3,IsoLog4,IsoLog5,IsoLog6,IsoLog7,IsoLog8,IsoLog9,IsoLog10, weights=TRUE)


summary(IsoLog2)
summary(IsoLog1)



# summary(L6)
# summary(n6)

plot(resid(IsoLog2) ~ DOYCentered, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$DOYCentered), col=2)



plot(resid(IsoLog2) ~ temp_centered, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$temp_centered), col=2)

plot(resid(IsoLog2) ~ percillum, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$percillum), col=2)


intervals(IsoLog2)

summary(IsoLog2)
(exp(-0.09597084)-1)*100 #-9.15% estimate  

(exp(-0.1388636)-1)*100 #-12.97% CI
(exp(-0.05307808)-1)*100 #-5.17% CI 



IsoAbu$LogIso<-log(IsoAbu$IsoPer100+10e-5)
IsoLog2 = lme(LogIso~DOYCentered,random=~1|Year,data=IsoAbu)


xvals <-  with(IsoAbu,seq(min(DOYCentered),max(DOYCentered),length.out=100))
nresamp<-1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked <- mvrnorm(nresamp, mu = fixef(IsoLog2), Sigma = vcov(IsoLog2))


## predicted values: useful below
pframe <- with(IsoAbu,data.frame(DOYCentered))
pframe$LogIso <- predict(IsoLog2,newdata=pframe,level=0)
head(pframe)
## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

## bootstrapping
sampfun <- function(fitted,data,idvar="Year") {
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted)
  dd <- data.frame(data,pred=pp,res=rr)
  ## sample groups with replacement
  iv <- levels(data[[idvar]])
  bsamp1 <- sample(iv,size=length(iv),replace=TRUE)
  bsamp2 <- lapply(bsamp1,
                   function(x) {
                     ## within groups, sample *residuals* with replacement
                     ddb <- dd[dd[[idvar]]==x,]
                     ## bootstrapped response = pred + bootstrapped residual
                     ddb$LogIso <- ddb$pred +
                       sample(ddb$res,size=nrow(ddb),replace=TRUE)
                     return(ddb)
                   })
  res <- do.call(rbind,bsamp2)  ## collect results
  if (is(data,"groupedData"))
    res <- groupedData(res,formula=formula(data))
  return(res)
}

pfun <- function(fm) {
  predict(fm,newdata=pframe,level=0)
}

set.seed(101)
yvals2 <- replicate(nresamp,
                    pfun(update(IsoLog2,data=sampfun(IsoLog2,IsoAbu,"Year"))))
c2 <- get_CI(yvals2,"boot_")
pframe2 <- data.frame(pframe,c2)
head(pframe2)
pframe2$DayOfYear<-pframe2$DOYCentered+mean(ShannonSubset$DayOfYear) #Remove centering on DOY for plotting
pframe2$IsoPer100<-exp(pframe2$LogIso)

pframe2$boot_lwr<-exp(pframe2$boot_lwr)
pframe2$boot_upr<-exp(pframe2$boot_upr)




DOYFamilyModelBoot<-ggplot(IsoAbu,aes(DayOfYear,IsoPer100))+geom_point(color="darkgrey")+geom_line(data=pframe2,size=1.5)+
  geom_line(data = pframe2, aes(y = boot_lwr), size = .75,linetype="dashed")+geom_line(data = pframe2,aes(y=boot_upr),size=0.75,linetype="dashed")+
  xlab("Day of Year")+ylab("Macroinvertebrate Family Richness")+ theme(axis.title.y = element_text(size = 8))
DOYFamilyModelBoot

dev.off()
tiff("Figures/DOYFamilyModelBoot.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DOYFamilyModelBoot
dev.off()




##############
#Heptageniidae model
############
HepAbu<-subset(FamilyLMMData,Family=="Heptageniidae")
head(HepAbu)
hist(HepAbu$Abundance)
HepAbu$HepPer100<-((HepAbu$Abundance*100)/(60*4*60*HepAbu$AreaSampled.m2.*HepAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
hist(log(HepAbu$Abundance))


HepLog0 = lme(log(HepPer100+10e-5)~1,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))

HepLog1 = lme(log(HepPer100+10e-5)~percillum,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog2 = lme(log(HepPer100+10e-5)~DOYCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog3 = lme(log(HepPer100+10e-5)~temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))

HepLog4 = lme(log(HepPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog5 = lme(log(HepPer100+10e-5)~percillum+DOYCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog6 = lme(log(HepPer100+10e-5)~temp_centered+DOYCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))

HepLog7 = lme(log(HepPer100+10e-5)~percillum+DOYCentered+temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog8 = lme(log(HepPer100+10e-5)~CTUSturgeon,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog9 = lme(log(HepPer100+10e-5)~CTUSturgeon+percillum,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))
HepLog10 = lme(log(HepPer100+10e-5)~CTUSturgeon,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DayOfYear|Year))


AICctab(HepLog0,HepLog1,HepLog2,HepLog3,HepLog4,HepLog5,HepLog6,HepLog7,HepLog8,HepLog9,HepLog10, weights=TRUE)

summary(HepLog2)


plot(resid(HepLog2) ~ DOYCentered, data=HepAbu)
lines(lowess(resid(HepLog2) ~ HepAbu$DOYCentered), col=2)


plot(resid(HepLog2) ~ temp_centered, data=HepAbu)
lines(lowess(resid(HepLog2) ~ HepAbu$temp_centered), col=2)


plot(resid(HepLog2) ~ percillum, data=HepAbu)
lines(lowess(resid(HepLog2) ~ HepAbu$percillum), col=2)


hist(resid(HepLog2))
plot(resid(HepLog2) ~ DischargeSampledByNight, data=ShannonSubset)
lines(lowess(resid(HepLog2) ~ ShannonSubset$DischargeSampledByNight), col=2)

plot(resid(HepLog2) ~ predict(HepLog2))
lines(lowess(resid(HepLog2) ~ predict(HepLog2)), col=2)

intervals(HepLog2,which=c("fixed"))

summary(HepLog2)
(exp(-0.05580188)-1)*100 #-5.43% estimate  Doy

(exp(-0.07709384)-1)*100 #-7.42% CI
(exp(-0.03450992)-1)*100 #-3.39% CI 




################
#Percentage Nightly Biomass
################
#Percentage Nightly Biomass
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonRichness$SturgeonBiomassNightly<-ShannonRichness$Nsturgeon*0.005
ShannonRichness$SuckerBiomassNightly<-ShannonRichness$Nsuckers100*0.00204
ShannonRichness$CombinedNightlyBiomass<-ShannonRichness$SturgeonBiomassNightly+ShannonRichness$SuckerBiomassNightly+ShannonRichness$InvertBiomass100


ShannonRichness$PercentNightlyBiomassSturgeon<- ShannonRichness$SturgeonBiomassNightly/ShannonRichness$CombinedNightlyBiomass*100
ShannonRichness$PercentNightlyBiomassSucker<- ShannonRichness$SuckerBiomassNightly/ShannonRichness$CombinedNightlyBiomass*100
ShannonRichness$PercentNightlyBiomassInvert<- ShannonRichness$InvertBiomass100/ShannonRichness$CombinedNightlyBiomass*100


SturgeonDataFrame<-data.frame("Sturgeon",ShannonRichness$DPFS, ShannonRichness$DayOfYear,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassSturgeon)
colnames(SturgeonDataFrame)<-c("Taxa","DPFS","DayOfYear","Year","PercentNightly")

SuckerDataFrame<-data.frame("Catostomidae",ShannonRichness$DPFS, ShannonRichness$DayOfYear,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassSucker)
colnames(SuckerDataFrame)<-c("Taxa","DPFS","DayOfYear","Year","PercentNightly")

InvertDataFrame<-data.frame("Invertebrate",ShannonRichness$DPFS,ShannonRichness$DayOfYear,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassInvert)
colnames(InvertDataFrame)<-c("Taxa","DPFS","DayOfYear","Year","PercentNightly")

PlottingDataframe<-rbind(SturgeonDataFrame,SuckerDataFrame,InvertDataFrame)
PlottingDataframe[is.na(PlottingDataframe)] <- 0 #Change NA values to 0 as no taxa were observed on those dates but sampling occured

head(PlottingDataframe)
PlottingDataframe$Taxa = factor(PlottingDataframe$Taxa, levels = c("Invertebrate","Sturgeon","Catostomidae"))

cbPalette2 <- c("#E69F00", "#000000", "#0072B2")
cbrewer<-c("#66c2a5","#fc8d62","#8da0cb")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")

PercentBiomassPlot<-ggplot(PlottingDataframe, aes(x=DayOfYear, y= PercentNightly,fill=Taxa))+facet_wrap(~Year)+geom_bar(stat="identity",lwd=0.1)+
  theme(legend.justification=c(1,0), legend.position=c(1,-0.04))+ylab("Nightly Biomass (%)")+xlab("Day of Year")+
  scale_fill_manual(values=cbPalette2)+scale_color_manual(values=cbPalette2)+theme(legend.text = element_text(size = 5),legend.title = element_blank())+ 
  theme(legend.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))#+ guides(shape = guide_legend(override.aes = list(size=2)))#+geom_line()
PercentBiomassPlot

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/PercentBiomassPlot.tiff", width = 84, height = 84, units = 'mm', res = 1200)
PercentBiomassPlot
dev.off()


##########
#Join Together Fig 4
#########
TopFamilyAbuPer100
TotalInvertAbuMoonPhase
CrayfishNoModelPlot
TopFamilyBiomass
PercentBiomassPlot

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/Fig5July2020.tiff", width = 174, height = 190, units = 'mm', res = 1200)
ggarrange(TopFamilyAbuPer100,CrayfishGAMPlot,PercentBiomassPlot,HepIsoBiomass,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()

###########
#Proportion collected by night fig
###########

head(AllData)
"#999999", "#E69F00", "#56B4E9", "#009E73"
Invert<-ggplot(AllData, aes(x=Date2,y=PerInvert))+geom_point(shape=17,color="#E69F00")+xlab("Date")+ylab("Proportion of Yearly Macroinvertebrates (%)")+theme(axis.title.x=element_blank())#theme(axis.title.y = element_text(size = 10))
Invert
Cat<-ggplot(AllData, aes(x=Date2,y=PerCat))+geom_point(shape=16,color="#0072B2")+xlab("Date")+ylab("Proportion of Yearly Catostomidae (%)")
Cat
#+geom_point(aes(x=CTUSturgeon,y= PerCat,group=Year),shape=17)+facet_wrap(~Year)+xlab("Date")+ylab("Proportion of Yearly Total (%)")+geom_point(aes(x=CTUSturgeon,y=PerStur),color="#E69F00",shape=15)
#+geom_line()
Stur<-ggplot(AllData, aes(x=Date2,y=PerStur))+geom_point(shape=15,color="#000000")+xlab("Date")+ylab("Proportion of Yearly Sturgeon (%)")#+geom_line()
Stur

dev.off()
tiff("Figures/ProportionYearlyTotal.tiff", width = 174, height = 190, units = 'mm', res = 1200)
ggarrange(Invert,Cat,Stur,
          labels = c("a", "b","c",""),
          ncol = 2, nrow = 2)
dev.off()

