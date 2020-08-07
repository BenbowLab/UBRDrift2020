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
metadata<-subset(metadata,Ninverts!=0) #Remove sample dates where  inverts were not collected
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
#####################
#General Result info
####################
AllData<-metadata
AllData$Date2<-as.Date(AllData$Date,format= "%d-%B")
head(AllData)
AllData$Year
SamplesByYear <- ddply(AllData, c("Year"), summarise,
                 N    = length(SampleID),
                 Sturgeon = sum(Nsturgeon),
                 Inverts5 = sum(Ninverts),
                 Inverts100 = sum(Ninverts100),
                 Catostomidae = sum(Nsuckers100)
)



SamplesByYear
mean(SamplesByYear$N) #Average number of days per year= 30
sum(SamplesByYear$N) #Total number of sampling days= 240
sum(AllData$Ninverts) #Number of invertebrates IDed = 21356
sum(AllData$Ninverts100) #100% numbers for inverts = 427,120
sum(AllData$Nsturgeon) #Total number of sturgeon larvae collected = 108,674
sum(AllData$Nsuckers100)#number of catostomidae larvae = 1,619,320


mean(AllData$percillum) #48.96%
min(AllData$percillum) #0
max(AllData$percillum) #100


mean(AllData$DPFS) #36.425
median(AllData$DPFS) #36
min(AllData$DPFS) #min 15
max(AllData$DPFS) #max 68


mean(AllData$Ninverts100) #1779.66
min(AllData$Ninverts100) #40
max(AllData$Ninverts100) #12340
sd   = sd(AllData$Ninverts100) #
sd
sd(AllData$Ninverts100) #1574.3
sd / sqrt(length(AllData$Ninverts100))
se #101.62


mean(AllData$Nsturgeon) #452.8
min(AllData$Nsturgeon) #0
max(AllData$Nsturgeon) #16,897
sd   = sd(AllData$Nsturgeon) #
sd(AllData$Nsturgeon)
 sd / sqrt(length(AllData$Nsturgeon)) #101,6

mean(AllData$Nsuckers100) #6747.16
min(AllData$Nsuckers100) #0
max(AllData$Nsuckers100)# 185,380

sd   = sd(AllData$Nsuckers100) #
sd / sqrt(length(AllData$Nsuckers100)) #1261.7

head(AllData)

mean(AllData$InvertBiomass100) #96.4
min(AllData$InvertBiomass100) #0.59
max(AllData$InvertBiomass100)#85.52


#RiverQ
AllDataQ<-subset(AllData,Q!= "NA") #Remove sample dates where no flow data was collected (to calc se correctly using length)
sd(metadata$PercentRiverDischargeSampled,na.rm=T)
mean(AllDataQ$Q,na.rm=T) #8.165 m3/sec
min(AllDataQ$Q,na.rm=T) #4.035 m3/sec
max(AllDataQ$Q,na.rm=T) #16.83 m3/sec
sd   = sd(AllDataQ$Q,na.rm=T)
se   = sd / sqrt(length(AllDataQ$Q))
sd #2.22
se # 0.152
hist(AllDataQ$Q,xlab = "River Discharge (m3/sec)")
RiverDischargePlot<-ggplot(AllDataQ,aes(x=Q))+geom_histogram(binwidth = 1,fill="grey",color="black")+xlab(expression(River~Discharge~(m^3/~sec)))+ylab("Frequency")
dev.off()
tiff("Figures/RiverDischarge.tiff", width = 74, height = 74, units = 'mm', res = 1000)
RiverDischargePlot
dev.off()

AllDataTemp<-subset(AllData,Temp!= "NA") #Remove sample dates where no temp data was collected (to calc se correctly using length)

mean(AllDataTemp$Temp,na.rm=T) #19.03
min(AllDataTemp$Temp,na.rm=T) #12.22
max(AllDataTemp$Temp,na.rm=T) #23.34
sd   = sd(AllDataTemp$Temp,na.rm=T)
se   = sd / sqrt(length(AllDataTemp$Temp))
sd #2.53
se #0.165

#Day of year sampled
min(AllData$DayOfYear) #130 April 17
max(AllData$DayOfYear) #186 May 9th 

#Percent River discharge sampled
AllDataPercent<-subset(AllData,PercentRiverDischargeSampled!= "NA") #Remove sample dates where no flow data was collected (to calc se correctly using length)

mean(AllDataPercent$PercentRiverDischargeSampled,na.rm=T) #12.44
min(AllDataPercent$PercentRiverDischargeSampled,na.rm=T) #4.89
max(AllDataPercent$PercentRiverDischargeSampled,na.rm=T) #45.31
sd= sd(AllDataPercent$PercentRiverDischargeSampled, na.rm=T)
se   = sd / sqrt(length(AllDataPercent$PercentRiverDischargeSampled))
sd #6.16
se #0.433

head(AllData)
ggplot(AllData,aes(x=Date,y=InvertBiomass100))+geom_point()+geom_smooth()+ylab("Inverts Biomass (g)")+xlab("Percent Illumination")
ggplot(AllData,aes(x=Q,y=Ninverts100))+geom_point()+facet_wrap(~Year)+xlab("River Discharge (m3/sec)")
ggplot(AllData,aes(x=Q,y=Ninverts100))+geom_point()+facet_wrap(~Year)+geom_smooth()+xlab("River Discharge (m3/sec)")
ggplot(AllData,aes(x=Q,y=Ninverts100))+geom_point()+xlab("River Discharge (m3/sec)")

ggplot(AllData,aes(x=DischargeSampledByNight,y=Ninverts100))+geom_point()+xlab("Discharge Sampled (m3/sec)")+facet_wrap(~Year)
ggplot(AllData,aes(x=DischargeSampledByNight,y=InvertBiomass100))+geom_point()+xlab("Discharge Sampled (m3/sec)")+geom_smooth()

ggplot(AllData,aes(x=DPFS,y=DischargeSampledByNight))+geom_point()+facet_wrap(~Year)

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

AllDataGaps=read.csv("DataClean\\SturgeonDriftMetadata5.6.2020WGAPS.csv",header=TRUE)#USE THIS FILE INSTEAD OF METADATA TO ADD BREAKS TO DISCHARGE LINE FOR GAPS IN DRIFT STURGEON COLLECTIONS

SturgeonCTU<-ggplot(AllDataGaps,aes(x=CTUSturgeon))+geom_point(size=0.75,aes(y=Nsturgeon),shape=19)+geom_line(aes(y=Q/0.002),color="blue")+scale_y_continuous(sec.axis = sec_axis(~.*0.002,name=expression(River~Discharge~(m^3/sec))))+xlab("Cumulative Thermal Units from First Spawn")+facet_wrap(~Year,scale="free_y")+ylab("Lake Sturgeon Larvae")
SturgeonCTU
dev.off()
tiff("Figures/SturgeonCTU.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonCTU
dev.off()


ggplot(AllData,aes(x=CTUSturgeon,y=Nsuckers100))+geom_point()+xlab("Cumulative Thermal Units")+facet_wrap(~Year)+ylab("Sucker Larvae")

min(AllData$Nsturgeon,na.rm=T)
max(AllData$Nsturgeon,na.rm=T)

min(AllData$Q,na.rm=T)
max(AllData$Q,na.rm=T)

cor(AllData$DayOfYear,AllData$percillum,method="pearson")
cor(AllData$DayOfYear,AllData$CTUSturgeon,method="pearson")

ggplot(AllData,aes(x=DPFS,y=Nsturgeon))+geom_point()+xlab("Days Post First Spawning (DPFS)")+facet_wrap(~Year)

ggplot(AllData,aes(x=DayOfYear,y=log(Ninverts100)))+geom_point()+xlab("Day of Year")#+facet_wrap(~Year)

ggplot(AllData,aes(x=Date2,y=DayOfYear))+geom_point()+ facet_wrap(~Year)#geom_jitter()

ggplot(AllData,aes(x=DayOfYear,y=LunarDay))+geom_point()+ facet_wrap(~Year)#geom_jitter()
AllData$MoonPhase = factor(AllData$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

ggplot(AllData,aes(x=MoonPhase,y=percillum))+geom_point()+ facet_wrap(~Year)#geom_jitter()
ggplot(AllData,aes(x=Date2,y=percillum))+geom_point()+ facet_wrap(~Year)#geom_jitter()

ggplot(AllData,aes(x=Date2,y=CTUSturgeon))+geom_point()+facet_wrap(~Year)

################
#Total Abu graphs
#################

SturgeonDataFrame<-data.frame("Sturgeon",ShannonRichness$DayOfYear, ShannonRichness$CTUSturgeon,ShannonRichness$Year,ShannonRichness$Nsturgeon,ShannonRichness$InvertsByDischargeSampled)
colnames(SturgeonDataFrame)<-c("Taxa","DayOfYear","CTUSturgeon","Year","Abundance","Conc")

SuckerDataFrame<-data.frame("Catostomidae",ShannonRichness$DayOfYear, ShannonRichness$CTUSturgeon,ShannonRichness$Year,ShannonRichness$Nsuckers100,ShannonRichness$SuckerConc)
colnames(SuckerDataFrame)<-c("Taxa","DayOfYear","CTUSturgeon","Year","Abundance","Conc")

InvertDataFrame<-data.frame("Invertebrate",ShannonRichness$DayOfYear,ShannonRichness$CTUSturgeon,ShannonRichness$Year,ShannonRichness$Ninverts100, ShannonRichness$SturgeonConc)
colnames(InvertDataFrame)<-c("Taxa","DayOfYear","CTUSturgeon","Year","Abundance","Conc")

PlottingDataframe<-rbind(SturgeonDataFrame,SuckerDataFrame,InvertDataFrame)
PlottingDataframe[is.na(PlottingDataframe)] <- 0 #Change NA values to 0 as no taxa were observed on those dates but sampling occured

head(PlottingDataframe)
PlottingDataframe$Taxa = factor(PlottingDataframe$Taxa, levels = c("Invertebrate","Sturgeon","Catostomidae"))

theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

TotalAbuAllPlot<-ggplot(PlottingDataframe, aes(x=DayOfYear, y= Abundance))+facet_grid(Taxa~.,scales = "free_y")+geom_point(size=1)+
  theme(legend.justification=c(1,0), legend.position=c(1,-0.01))+ylab(expression(Abundance~Per~100~m^3~Drift))+xlab("Calendar Date")+
  scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+theme(legend.text = element_text(size = 7),legend.title = element_blank())+ 
  theme(legend.background=element_blank())+ guides(shape = guide_legend(override.aes = list(size=2)))#+geom_line()
TotalAbuAllPlot

#######
########Data by quantile ranks
#######

AllDataDPFS= mutate(AllData, quantile_rank = ntile(AllData$percillum,4))
head(AllDataDPFS)
Trtdata <- ddply(AllDataDPFS, c("quantile_rank"), summarise,
                 N    = length(Ninverts100),
                 meanInverts = mean(Ninverts100),
                 sd   = sd(Ninverts100),
                 se   = sd / sqrt(N)
)
ggplot(Trtdata,aes(x=quantile_rank,y=meanInverts))+geom_point()+xlab("Percent Illumination Quartiles")+
  geom_errorbar(aes(ymin=meanInverts-se,ymax=meanInverts+se))+ylab("Invertebrates Per Night (+/- SE)")




ShannonRichnessDischargeSampled = mutate(AllData, quantile_rank = ntile(AllData$DischargeSampledByNight,4))
head(ShannonRichnessDischargeSampled)
Trtdata <- ddply(ShannonRichnessDischargeSampled, c("quantile_rank"), summarise,
                 N    = length(Ninverts100),
                 meanInverts = mean(Ninverts100),
                 sd   = sd(Ninverts100),
                 se   = sd / sqrt(N)
)
ggplot(Trtdata,aes(x=quantile_rank,y=meanInverts))+geom_point()+xlab("Discharge Sampled Quartiles")+
  geom_errorbar(aes(ymin=meanInverts-sd,ymax=meanInverts+sd))+ylab("Invertebrates Per Night (+/- 1 SD)")

Trtdata <- ddply(ShannonRichnessDischargeSampled, c("quantile_rank"), summarise,
                 N    = length(InvertBiomass100),
                 meanInverts = mean(InvertBiomass100),
                 sd   = sd(InvertBiomass100),
                 se   = sd / sqrt(N)
)
ggplot(Trtdata,aes(x=quantile_rank,y=meanInverts))+geom_point()+xlab("Discharge Sampled Quartiles")+
  geom_errorbar(aes(ymin=meanInverts-sd,ymax=meanInverts+sd))+ylab("Invertebrate Biomass Per Night (+/- 1 SD)")
############
#Combined Summary plots
###########
AllData<-metadata
SamplesByYear <- ddply(AllData, c("Year"), summarise,
                       N    = length(SampleID),
                       Sturgeon = sum(Nsturgeon),
                       SturgeonByNight = mean(Nsturgeon),
                       sdSturgeonByNight = sd(Nsturgeon),
                       seSturgeonByNight = sdSturgeonByNight/sqrt(N),
                       Inverts5 = sum(Ninverts),
                       Inverts100 = sum(Ninverts100),
                       InvertsByNight = mean(Ninverts100),
                       sdInvertsByNight = sd(Ninverts100),
                       seInvertsByNight = sdInvertsByNight/sqrt(N),
                       CatostomidaeByNight = mean(Nsuckers100),
                       sdSuckers = sd(Nsuckers100),
                       seSuckers = sdSuckers/sqrt(N)
)

SamplesByYear
SturgeonByYearByNight<- ggplot(SamplesByYear, aes(x=Year,y=SturgeonByNight))+geom_bar(stat="identity",color="black",fill = "grey")+ylab("Sturgeon Per Night (SE)")+
  geom_errorbar(aes(ymin=SturgeonByNight-seSturgeonByNight,ymax=SturgeonByNight+seSturgeonByNight))
InvertsByYearByNight<-ggplot(SamplesByYear, aes(x=Year,y=InvertsByNight))+geom_bar(stat="identity",color="black",fill = "grey")+ylab("Invertebrates Per Night (SE)")+
  geom_errorbar(aes(ymin=InvertsByNight-seInvertsByNight,ymax=InvertsByNight+seInvertsByNight))
SuckersByYearByNight<-ggplot(SamplesByYear, aes(x=Year,y=CatostomidaeByNight))+geom_bar(stat="identity",color="black",fill = "grey")+ylab("Catostomidae larvae Per Night (SE)")+
  geom_errorbar(aes(ymin=CatostomidaeByNight-seSuckers,ymax=CatostomidaeByNight+seSuckers))
SturgeonByYearByNight
InvertsByYearByNight
SuckersByYearByNight

SturgeonByYear<- ggplot(SamplesByYear, aes(x=Year,y=Sturgeon))+geom_bar(stat="identity",color="black",fill = "grey")+ylab("Sturgeon larvae")
SturgeonByYear
InvertsByYear<-ggplot(SamplesByYear, aes(x=Year,y=Inverts100))+geom_bar(stat="identity",color="black",fill = "grey")+ylab("Invertebrates collected")




theme_set(theme_bw(base_size = 14)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TotalsByYearCombined.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
ggarrange(SturgeonByYear,InvertsByYear,SturgeonByYearByNight,InvertsByYearByNight,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Sucker abu by day and year
Trtdata <- ddply(AllData, c("DPFS","Year"), summarise,
                 N    = length(Nsuckers100),
                 meanSuckers = mean(Nsuckers100)
)
#Trtdata





SuckersByDPFS<-ggplot(Trtdata, aes(x=DPFS,y=meanSuckers))+geom_bar(colour="black", stat="identity")+xlab("Days Post First Spawning")+ylab("Catostomidae Abundance")+
 theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
SuckersByDPFS

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/SuckerAbundanceSupplementalFig.tiff", width = 6.85, height = 4, units = 'in', res = 800)
ggarrange(SuckersByYearByNight,SuckersByDPFS,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)
dev.off()
#Sturgeon abu by day and year
Trtdata <- ddply(AllData, c("DPFS","Year"), summarise,
                 N    = length(Nsuckers100),
                 meanSturgeon = mean(Nsturgeon)
)
#Trtdata
SturgeonByDPFS<-ggplot(Trtdata, aes(x=DPFS,y=meanSturgeon))+geom_bar(colour="black", stat="identity")+xlab("Days Post First Spawning")+ylab("Larval Sturgeon Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
SturgeonByDPFS

#Invertebrate abu by day and year
Trtdata <- ddply(AllData, c("DPFS","Year"), summarise,
                 N    = length(Ninverts100),
                 meanInvert = mean(Ninverts100)
)
#Trtdata
InvertByDPFS<-ggplot(Trtdata, aes(x=DPFS,y=meanInvert))+geom_bar(colour="black", stat="identity")+xlab("Days Post First Spawning")+ylab("Invertebrate Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
InvertByDPFS

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/SturgeonInvertByDPFSByYear.tiff", width = 174, height = 174, units = 'mm', res = 1000)
ggarrange(SturgeonByYearByNight,InvertsByYearByNight,SturgeonByDPFS, InvertByDPFS,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()
#Sucker supplemental figure
dev.off()
tiff("Figures/SuckersCombinedByYear.tiff", width = 174, height = 84, units = 'mm', res = 1000)
ggarrange(SuckersByYearByNight,SuckersByDPFS,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)
dev.off()

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))




#####################
#Alpha diversity
####################
#Shannon Richness
head(sample_data(physeq))
Richness=plot_richness(physeq, x="DPFS", measures=c("Shannon"))#+geom_boxplot(aes(x=DPFS, y=value, color=DPFS), alpha=0.05)
Richness+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("DPFS")
#Richness$data

write.csv(Richness$data, "SturgeonMetadataWDiversity.csv")
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)
head(ShannonRichness)
levels(ShannonRichness$MoonPhase)
levels(metadata$MoonPhase)
#Days post first spawn
Trtdata <- ddply(ShannonRichness, c("DPFS"), summarise,
                 N    = length(value),
                 meanShannon = mean(value),
                 sd   = sd(value),
                 se   = sd / sqrt(N)
)
#Trtdata
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
ggplot(Trtdata, aes(x=DPFS,y=meanShannon))+geom_bar(aes(),colour="black", stat="identity")+xlab("Days post first spawn")+ylab("Shannon (SEM)")+
  geom_errorbar(aes(ymin=meanShannon-se,ymax=meanShannon+se))

hist(ShannonRichness$percillum)
hist(log(ShannonRichness$Ninverts100))
hist((ShannonRichness$Ninverts100))


#ShannonSubset$percillum
hist(ShannonRichness$Ninverts100)
m1 = glmer(Ninverts100~1+(1|Year),data=ShannonRichness,family = poisson)
m2 = glmer(Ninverts100~percillum+(1|Year),data=ShannonRichness,family = poisson)
#m3 = glmer.nb(Ninverts100~DPFS+percillum+(1|Year),data=ShannonSubset)

m1=glm(Ninverts100~percillum,data=ShannonRichness,family=poisson)
summary(m1)


head(ShannonRichness)
AIC(m1,m2)
summary(m2)
plot(Ninverts100~percillum,data=ShannonRichness)
lines(predict(m2) ~ ShannonRichness$percillum)

ShannonRichness$MoonPhase = factor(ShannonRichness$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))
levels(ShannonRichness$MoonPhase)
AbuBoxplot<-ggplot(ShannonRichness, aes(x=MoonPhase,y=Ninverts100))+geom_boxplot()+scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Total Invertebrates Per Night")+xlab("Moon Phase")

dev.off()
tiff("Figures/AbuBoxplotByMoonPhase.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
AbuBoxplot
dev.off()

AbuDotplotDPFS<-ggplot(ShannonRichness, aes(x=DPFS,y=Ninverts100))+geom_point()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+ylab("Total Invertebrates Per Night")+xlab("Days Post First Spawning")
AbuDotplotDPFS
dev.off()
tiff("Figures/AbuDotplotByMoonPhase.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
AbuDotplotDPFS
dev.off()







#By calendar date
Trtdata <- ddply(ShannonRichness, c("DayOfYear"), summarise,
                 N    = length(value),
                 meanShannon = mean(value),
                 sd   = sd(value),
                 se   = sd / sqrt(N)
)

ggplot(Trtdata, aes(x=DayOfYear,y=meanShannon))+geom_bar(colour="black", stat="identity")+xlab("Calendar date")+ylab("Shannon (SEM)")+
  geom_errorbar(aes(ymin=meanShannon-se,ymax=meanShannon+se))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Moon phase
Trtdata <- ddply(ShannonRichness, c("MoonPhase"), summarise,
                 N    = length(value),
                 meanShannon = mean(value),
                 sd   = sd(value),
                 se   = sd / sqrt(N)
)
Trtdata
ggplot(Trtdata, aes(x=MoonPhase,y=meanShannon))+geom_bar(aes(fill = MoonPhase),colour="black", stat="identity")
ggplot(Trtdata, aes(x=MoonPhase,y=meanShannon))+geom_bar(aes(fill = MoonPhase),colour="black", stat="identity")+xlab("Moon Phase")+ylab("Shannon Diversity (SEM)")+
  geom_errorbar(aes(ymin=meanShannon-se,ymax=meanShannon+se))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=cbPalette)+ theme(legend.position = "none")+
  geom_text(x=6.5,y=2,label="KW, Chi2=10.5, p= 0.158")

kruskal.test(value~MoonPhase,data=ShannonRichness)

hist(ShannonRichness$value)


ggplot(ShannonRichness,aes(x=LunarDay,y=value))+geom_point()

#Temperature 
ggplot(ShannonRichness, aes(x=Temp,y=value))+geom_point()+xlab("Calendar date")+ylab("Shannon (SEM)")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_smooth()
Temps<-subset(ShannonRichness, Temp!= "NA")
Temps
head(Temps)
mean(Temps$Temp)
Temps$temp_centered = Temps$Temp - mean(Temps$Temp)
plot(Temps$value~Temps$temp_centered,xlab="Change from mean temp (19.027)", ylab="Shannon diversity")+abline(1.504,0.02275)

m = glm(value~ temp_centered, data=Temps)
summary(m)
m2 = lmer(value~temp_centered+(1|Year),data=Temps)
m1 = lmer(value~1+(1|Year),data=Temps)

anova(m1,m2,test='Chisq')
confint(m2)
summary(m2)

m4 = lmer(value~temp_centered*MoonPhase+(1|Year),data=Temps)
m5 = lmer(value~temp_centered+MoonPhase+(1|Year),data=Temps)
summary(m2)



AIC(m1,m3,m4,m5)

m1=glm(value~temp_centered+Year,data=Temps)

plot(residuals(m1)~Temp,data=Temps)
plot(resid(m1) ~ predict(m1))


m2 = lmer(value~DPFS+(1|Year),data=ShannonRichness,family)
m1 = lmer(value~1+(1|Year),data=ShannonRichness)
m1<-glm(value~temp_centered+percillum,data=Temps)
anova(m1,m2,test='Chisq')
confint(m2)
summary(m2)
head(Temps)

Percillum<-subset(ShannonRichness, percillum!= "#VALUE!")

Percillum$ni
plot(Ninverts100~percillum,data=Percillum)
Percillum$percillum<-as.numeric(Percillum$percillum)
m2<-(glm(value~percillum+Q,data=Percillum))
m2
plot(residuals(m2)~percillum,data=Percillum)

confint(m2)
hist(Percillum$percillum)
Percillum$percillum
#Total abundance percent illum test
m4 = glmer(Ninverts100~percillum+(1|Year),data=Percillum,family="poisson")
m5 = glmer(Ninverts100~1+(1|Year),data=Percillum,family="poisson")
anova(m4,m5,test='chisq')
confint(m4)
exp(.01)
m5 = glmer(value~temp_centered+MoonPhase+(1|Year),data=Temps,family="poisson")
summary(m4)
confint(m4)
summary(glm(Ninverts100~percillum,data=Percillum))
    



#By discharge

ShannonRichnessQ = mutate(ShannonRichness, quantile_rank = ntile(ShannonRichness$Q,4))
head(ShannonRichnessQ)
Trtdata <- ddply(ShannonRichnessQ, c("quantile_rank"), summarise,
                 N    = length(value),
                 meanShannon = mean(value),
                 sd   = sd(value),
                 se   = sd / sqrt(N)
)
#Trtdata
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
ggplot(Trtdata, aes(x=quantile_rank,y=meanShannon))+geom_bar(aes(),colour="black", stat="identity")+xlab("Shannon Quartile")+ylab("Shannon (SEM)")+
  geom_errorbar(aes(ymin=meanShannon-se,ymax=meanShannon+se))

hist(ShannonRichnessQ$value)

####################
##Total Abundance larval fish
####################
#Sturgeon
ggplot(ShannonRichness,aes(x=CTUSturgeon,y=Nsturgeon))+geom_point()#+geom_bar(aes(fill=DPFS),stat = "identity")

SubsetCTU<-subset(ShannonRichness, CTUSturgeon < 500&CTUSturgeon>250)
sum(SubsetCTU$Nsturgeon)/sum(ShannonRichness$Nsturgeon)*100
max(ShannonRichness$CTUSturgeon)
SubsetDrift<-subset(ShannonRichness,Nsturgeon>0)
SubsetDrift$Nsturgeon
min(SubsetDrift$CTUSturgeon)
Trtdata <- ddply(ShannonRichness, c("DPFS","Year"), summarise,
                 N    = length(Nsturgeon),
                 meanSturgeon = mean(Nsturgeon),
                 sd   = sd(Nsturgeon),
                 se   = sd / sqrt(N)
)
Trtdata
#Subset<-subset(Trtdata,Year=="2018"|Year=="2017"|Year=="2016")

SturgeonByDPFS<-ggplot(Trtdata, aes(x=DPFS,y=meanSturgeon))+geom_bar(colour="black", stat="identity")+xlab("Days Post First Spawning")+ylab("Larval Sturgeon Abundance")+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)

dev.off()
tiff("Figures/SturgeonByDPFSYear.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
SturgeonByDPFS
dev.off()



#SturgeonByMoonPhase
SturgeonLarvaeSubset<-subset(ShannonRichness,Nsturgeon>0)
Trtdata <- ddply(SturgeonLarvaeSubset, c("DayOfYear","Year","MoonPhase"), summarise,
                 N    = length(Nsturgeon),
                 meanSturgeon = mean(Nsturgeon),
                 sd   = sd(Nsturgeon),
                 se   = sd / sqrt(N)
)
Trtdata

ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon,color=MoonPhase))+geom_point()+xlab("DoY")+ylab("Larval Sturgeon Abundance (SEM)")+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+facet_grid(~Year)#+scale_fill_manual(values=cbPalette)

ggplot(ShannonRichness, aes(x=CTUSturgeon,y=SturgeonConc))+geom_point()+xlab("DoY")+ylab("Larval Sturgeon Abundance (SEM)")#+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_boxplot()#+facet_grid(~Year)#+scale_fill_manual(values=cbPalette)


  
ggplot(ShannonRichness, aes(x=CTUSturgeon,y=Nsuckers100))+geom_point()
ggplot(ShannonRichness, aes(x=Temp,y=Nsturgeon))+geom_point()
ggplot(ShannonRichness, aes(x=Temp,y=Ninverts100))+geom_point()

head(ShannonRichness)
gam_y <- gam(log(Nsuckers100+10e-5) ~ s(CTUSturgeon)+s(percillum)+s(Q)+s(Year,bs="re"), method = "REML",data=ShannonRichness)
gam_x <- gam(log(Nsuckers100+10e-5) ~ s(CTUSturgeon)+s(Year,bs="re")+s(percillum)+s(DischargeSampledByNight), method = "REML",data=ShannonRichness)
#gam.check(gam_y)
summary(gam_y)
gam.check(gam_y)

acf(resid(gam_y), lag.max = 36, main = "ACF")



library(mgcViz)
plot(gam_y)


dat <- gamSim(6,n=200,scale=.2,dist="poisson")
b2 <- gamm(y~s(x0)+s(x1)+s(x2),family=poisson,
           data=dat,random=list(fac=~1))


fac <- dat$fac
ShannonCatGamSubset<-subset(ShannonRichness, CTUSturgeon!="NA")
head(ShannonCatGamSubset)
fac<-ShannonRichness$Year


ShannonRichness$Nsuckers100

SuckerGAM<-ggplot(ShannonRichness, aes(CTUSturgeon, (log(Nsuckers100)) )) + geom_point(color="grey") +stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE,color="black")+
  stat_smooth(method = "gam", formula = (y) ~ s(x),geom="ribbon",linetype="dashed",fill=NA,color="black",level = 0.95)+ylab("log(Larval Catostomidae Abundance)")+xlab("Cumulative Temperature Units")
SuckerGAM$layers



ShannonRichness$CTUSturgeon
SuckerAbuGamma <- gamm(((Nsuckers100+10e-5))~ s(CTUSturgeon), data = ShannonRichness, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))
summary(SuckerAbuGamma$gam)
gam.check(SuckerAbuGamma$gam)

fac<-ShannonRichness$Year
hist(log(ShannonRichness$Nsuckers100))
ShannonRichness$Year<-as.factor(ShannonRichness$Year)
SuckerAbuLog <- gam(log(Nsuckers100+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonRichness, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
summary(SuckerAbuLog)
gam.check(SuckerAbuLog)#default K doing ok no k significance, model residuals mostly ok around zero but a couple outliers around -10 to -15 (<5)
appraise(SuckerAbuLog)
response1 <- predict(SuckerAbuLog, type="response", se.fit=T)
head(response1)

plot(0,type="n",xlim=c(150,800),ylim=c(-10,500))
lines((smooth.spline(SuckerAbuLog$model$CTUSturgeon,response1$fit)),col="red")

data<-smooth.spline(SuckerAbuLog$model$CTUSturgeon,response1$fit)
head(data)
GAMConfints<-confint(SuckerAbuLog,parm="CTUSturgeon",type="confidence",nsim=1000)
head(GAMConfints)
PredictedDataFrame<-data.frame(exp(GAMConfints$est),GAMMConfints$CTUSturgeon,exp(GAMMConfints$upper),exp(GAMMConfints$lower))
colnames(PredictedDataFrame)<-c("Nsuckers100","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)




plot(SuckerAbuLog)






observed_fitted_plot(
  SuckerAbuLog,
  ylab = NULL,
  xlab = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL
)

library(mgcViz)

Random<-ranef(SuckerAbuLog$lme)
Random$fac

newCTUDataata<-seq(min(ShannonRichness$CTUSturgeon, max(ShannonRichness$CTUSturgeon,1000)))
predict(SuckerAbuLog,newCTUDataata)
plot(SuckerAbuLog$gam)


AICctab(SuckerAbuLog$lme,SuckerAbuGamma$lme)
summary(SuckerAbuLog$gam)

new.CTU = seq(min(ShannonRichness$CTUSturgeon),max(ShannonRichness$CTUSturgeon),length= 10000)
SuckerAbuLog <- gam(log(Nsuckers100+10e-5)~ s(CTUSturgeon), data = ShannonRichness, method = "REML")

predict(SuckerAbuLog,newdata=new.CTU,type="response",se.fit=T)


newCTUData<-seq(min(ShannonRichness$CTUSturgeon, max(ShannonRichness$CTUSturgeon,length=1000)))
min(ShannonRichness$CTUSturgeon)
head(newCTUData)
fac<-ShannonRichness$Year
library(mgcViz)


library(gratia)
SuckerAbuLogViz <- gammV(log(Nsuckers100+10e-5)~ s(CTUSturgeon),random=list(fac=~1), data = ShannonRichness, method = "REML")
summary(SuckerAbuLogViz)

SuckerAbuLog <- gamm((log(Nsuckers100+10e-5))~ s(CTUSturgeon), data = ShannonRichness, method = "REML")

GAMMConfints<-confint(SuckerAbuGamma,parm="CTUSturgeon",type="confidence",nsim=1000)
head(GAMMConfints)

PredictedDataFrame<-data.frame(exp(GAMMConfints$est-10e-5),GAMMConfints$CTUSturgeon,exp(GAMMConfints$upper),exp(GAMMConfints$lower))
colnames(PredictedDataFrame)<-c("Nsuckers100","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)

GAMMSucker<-ggplot(ShannonRichness,aes(CTUSturgeon,Nsuckers100))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)#+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab(expression(Discharge~(m^3/sec)~Sampled))+ylab("Shannon Diversity")#+ theme(axis.title.y = element_text(size = 9))
GAMMSucker

ShannonSuckerSubset<-subset(ShannonRichness, Nsuckers100< 150000)
SuckerAbuGamma <- gam((Nsuckers100+10e-5)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonRichness, method = "REML",correlation=corAR1(form = ~DayOfYear|Year),family= Gamma(link = "log"))

appraise(SuckerAbuGamma)




SuckerAbuGamma$residuals
SuckerAbuGamma$aic
hist(ShannonRichness$Nsuckers100)
###########
#Total InvertN By Moon Phase Normalized 
############

Model<-aov(log(DriftInvertConc)~MoonPhase,data=ShannonRichness)
summary(Model)
hist(resid(Model))

kruskal.test(DriftInvertConc~MoonPhase, data=ShannonRichness)
Tukey<-TukeyHSD(Model,"MoonPhase")

Tukey
difference<-Tukey$MoonPhase[,"p adj"]
Letters<-multcompLetters(difference)

Letters
  
# 
# 
#  compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")
# 
#  Means=compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")
# 
#  Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
#  difference<-Means$p.adj
#  names(difference)<-Hyphenated
#  Letters<-multcompLetters(difference)
#  Letters
#  Letters$Letters
#  #manually renamed due to some weirdness with plotting where new moon donsn't start with a

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

TotalInvertAbuMoonPhase<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(fill=MoonPhase),stat="identity")+xlab("Moon Phase")+ylab(expression(Invertebrates~Per~100~m^3~Drift~(SE)))+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 10))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+1,label=LettersRearranged))+theme(legend.position = "none")+geom_text(aes(x=5,y=25,label= "ANOVA, F = 5.99, P < 0.001"),size=4)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
TotalInvertAbuMoonPhase
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TotalInvertebrateAbuByPhaseNormalized.tiff", width = 74, height = 74, units = 'mm', res = 1200)
TotalInvertAbuMoonPhase
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

compare_means(DriftBiomassConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Means=compare_means(DriftBiomassConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
Letters<-multcompLetters(difference)
Letters
Trtdata
vector<-c("bc","bc","abc","a","ab","b","ab","c")
DriftTotalInvertBiomass<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(fill=MoonPhase),stat="identity")+xlab("Moon Phase")+ylab(expression(Invertebrate~Biomass~(g)~Per~100~m^3~Drift~(SE)))+#Invertebrate Biomass (g) per 100 m3 drift (SEM)
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 8))+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+.08,label=vector))+geom_text(aes(x=5,y=1.5,label= "KW, chi-squared = 49.3, P < 0.001"),size=3)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
DriftTotalInvertBiomass
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/BiomassByPhase100m3.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
DriftTotalInvertBiomass
dev.off()




ShannonRichness


###############
#Relative abundance invert families by moon phase
######
RelTaxa3Per  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
RelTaxa3Per = filter_taxa(RelTaxa3Per, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 0.1%

df <- psmelt(RelTaxa3Per)
df$Abundance<-df$Abundance*100
head(df)
Trtdata <- ddply(df, c("MoonPhase","Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Trtdata
SubsetTotalAbuFamilies<-subset(Trtdata,Family=="Chironomidae"|Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")
SubsetTotalAbuFamilies
write.csv(SubsetTotalAbuFamilies,file="FamilyLevelTotalAbundanceMoonPhase.csv") #Write in other category sum to 100
Trtdata<-read.csv("FamilyLevelTotalAbundanceMoonPhase.csv", header=T)
head(Trtdata)
Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))
#Trtdata

compare_means(Abundance ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")
Means<-compare_means(Abundance ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="wilcox.test")

SigList<-length(unique(Trtdata$Family))
for (i in levels(Means$Family)){
  Tax<-i
  TaxAbundance<-subset(Means,Family==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  #print(Letters)
  SigList[i]<-Letters
  
}
vec<-unlist(SigList)
vec<-vec[-1]

(vec)
#For unadjusted p values
# vector<-c("","","","","","","","",
#           "a","ab","ab","b","b","ab","a","a",
#           "","","","","","","","",
#           "abc","ab","abc","a","abc","c","bc","bc",
#           "","","","","","","","",
#           "a","ab","cd","c","bcd","d","abd","a",
#           "","","","","","","","")
# length(vec)
# length(vector)

TrtdataSorted<-Trtdata[order(Trtdata$Family),]
head(TrtdataSorted)

KruskalLabel<- c("Kruskal-Wallis,\n P-adj = 0.45","KW, P-adj = 0.04","KW, P-adj = 1", "KW, P-adj = 0.016","KW, P-adj = 0.5","     P-adj <0.001","KW, P-adj = 0.62")
FamilyRelativeAbu=ggplot(TrtdataSorted, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab("Relative Invertebrate Abundance (%, SEM)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+ annotate("text", label = KruskalLabel, size = 1.75, x = 4, y = 50)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
FamilyRelativeAbu
theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FamiyLevelRelAbu.tiff", width = 84, height = 84, units = 'mm', res = 1200)
FamilyRelativeAbu
dev.off()



############
#Total Abundance Family by Moon phase
############
DriftDataCombined<-read.csv("DataClean\\AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)
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
physeqSubset<-subset_taxa(physeq,Family=="Chironomidae"|Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")
physeqSubset
df <- psmelt(physeqSubset)
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
 "","","","","","","","",
"a","a","ab","b","ab","a","a","a",
"a","a","a","a","a","a","a","a",
"a","ab","cd","d","bcd","bc","abc","abc",
"","","","","","","","")


length(vector)
length(TrtdataSorted$N)
KruskalLabel<- c("Kruskal-Wallis,\n P-adj = 0.19","KW, P-adj = 0.053","KW, P-adj = 0.89", "KW, P-adj < 0.001","KW, P-adj = 0.01","     P-adj < 0.001","KW, P-adj = 0.17")

dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.19","KW, P-adj = 0.053","KW, P-adj = 0.89", "KW, P-adj < 0.001","KW, P-adj = 0.01","     P-adj < 0.001","KW, P-adj = 0.17"),
  Family   = c("Chironomidae","Crayfish","Ephemerillidae","Heptageniidae","Hydropsychidae","Isonychiidae","Leptoceridae")
)
TrtdataSorted

TopFamilyAbuPer100=ggplot(TrtdataSorted, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Invertebrates~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+geom_text(aes(x=MoonPhase,y= mean+se+1.5),label=vector,size=1.6)
TopFamilyAbuPer100
TopFamilyAbuPer100<-TopFamilyAbuPer100+ geom_text(data=dat_text,size=2,mapping = aes(x = 4, y = 15, label = label))
TopFamilyAbuPer100

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FamilyLevelAbuPer100m3.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TopFamilyAbuPer100
dev.off()

##############
#Biomass By Family
##############

BiomassAverages<-read.table("DataClean\\BiomassFamilyAveragesR.txt",header=T)
head(BiomassAverages)
length(BiomassAverages)
DriftDataCombined<-read.csv("DataClean\\AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)
head(DriftDataCombined)
head(metadata)
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

Trtdata <- ddply(dfSubset, c("Family"), summarise,
                 N    = length(BiomassPer100),
                 mean = mean(BiomassPer100))
Trtdata            
TrtdataSorted<-Trtdata[order(-Trtdata$mean),]
head(TrtdataSorted)         


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
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))#+ annotate("text", label = KruskalLabel, size = 2, x = 4, y = .7)+scale_fill_manual(values=cbPalette)
AllFamilyBiomass
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/AllFamilyBiomass.tiff", width = 174, height = 174, units = 'mm', res = 1200)
AllFamilyBiomass
dev.off()


#Percentage biomass
AllFamilyBiomassPercent=ggplot(Trtdata, aes(x=MoonPhase,y=meanPercent))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab("Relative Biomass (%, SE)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=meanPercent-sePercent,ymax=meanPercent+sePercent))+
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))#+ annotate("text", label = KruskalLabel, size = 2, x = 4, y = .7)+scale_fill_manual(values=cbPalette)
AllFamilyBiomassPercent
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/AllFamilyBiomassPercent.tiff", width = 174, height = 174, units = 'mm', res = 1200)
AllFamilyBiomassPercent
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
#kruskal.test(Heptageniidae~MoonPhase,data=DriftDataCombined) Use adj-p values
Means<-compare_means(Heptageniidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHep<-multcompLetters(difference)
LettersHep

#Isonychiidae
#kruskal.test(Isonychiidae~MoonPhase,data=DriftDataCombined) Use adj-p values
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
vector<-c("","","","","","","","",
          "","","","","","","","",
          "","","","","","","","",
          "a","a","ab","b","ab","a","a","a",
          "a","ab","cd","d","bcd","bc","abc","abc",
          "","","","","","","","")


dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.2","KW, P-adj = 0.89","KW, P-adj = 0.55", "KW, P-adj < 0.001",  "      P-adj < 0.001","KW, P-adj = 0.89"),
  Family   = c("Crayfish","Ephemerillidae","Gomphidae","Heptageniidae","Isonychiidae","Lepidostomatidae")
)
TopFamilyBiomass=ggplot(TrtdataSubset, aes(x=MoonPhase,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Biomass~(g)~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+ geom_text(data=dat_text,size=2,mapping = aes(x = 4, y = 0.9, label = label))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase,y= mean+se+.05),label=vector,size=1.6)
TopFamilyBiomass
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TopFamiliesBiomass.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TopFamilyBiomass
dev.off()

cor(ShannonSubset$DPFS,ShannonSubset$percillum)






###########
#Beta diversity/Envfit
##############

physeq
ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="MoonPhase")+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplotTax<-ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 0, aes(fill = MoonPhase))+
  theme(legend.text = element_text(size = 8))+ theme(legend.background=element_blank())+geom_point(size=2.5)
ordplotTax
dev.off()
tiff("Figures/PCoATax.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ordplotTax
dev.off()


#Envfit
#

plot(value~DischargeSampledByNight,data=ShannonSubset)
physeqSubset<- subset_samples(physeq, percillum!="NA")
physeqSubset
physeqSubset<- subset_samples(physeqSubset, DischargeSampledByNight !="NA") #Remove samples without discharge info

physeqSubset<- subset_samples(physeqSubset, DischargeSampledByNight < 2) #Remove discharge outliers

physeqSubset
GPdist=phyloseq::distance(physeqSubset, "bray")

vare.mds= ordinate(physeqSubset, "NMDS",GPdist)
#vare.mds <- metaMDS(VeganDist, trace = FALSE)
#vare.mds
metadataEnvfitSubset<-subset(metadata,DischargeSampledByNight !="NA"&DischargeSampledByNight < 2)
EnvFitMeta=data.frame(metadataEnvfitSubset$percillum,metadataEnvfitSubset$DayOfYear,metadataEnvfitSubset$DischargeSampledByNight,metadataEnvfitSubset$Q)

head(EnvFitMeta)
colnames(EnvFitMeta)<-c("Percent \n Illumination","Day of Year","Discharge\n Sampled","Q")
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
plot(ef, p.max = 0.05,cex=0.7)

#mtext("NMDS2", side = 2, line = 1, cex = 1)

envplot=plot(ef, p.max = 0.05,cex=0.7)
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

#mtext("NMDS2", side = 2, line = 1, cex = 1)

envplot=plot(ef, p.max = 0.05,cex=0.7)



GPdist=phyloseq::distance(physeq, "jaccard")
adonis(GPdist ~ Year, as(sample_data(physeq), "data.frame"))


ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="MoonPhase")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = MoonPhase))+ #theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~Year)
ordplot








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

(ShannonSubset$Year)

InvertsByDischargeSampledDL0 = lme(log(DriftInvertConc)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL1 = lme(log(DriftInvertConc)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL2 = lme(log(DriftInvertConc)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL3 = lme(log(DriftInvertConc)~DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

InvertsByDischargeSampledDL4 = lme(log(DriftInvertConc)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL5 = lme(log(DriftInvertConc)~percillum+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
InvertsByDischargeSampledDL6 = lme(log(DriftInvertConc)~temp_centered+DoYCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))

InvertsByDischargeSampledDL7 = lme(log(DriftInvertConc)~percillum+DoYCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DayOfYear|Year))
plot(intervals(InvertsByDischargeSampledDL3))
AICctab(InvertsByDischargeSampledDL0,InvertsByDischargeSampledDL1,InvertsByDischargeSampledDL2,InvertsByDischargeSampledDL3,InvertsByDischargeSampledDL4,InvertsByDischargeSampledDL5,InvertsByDischargeSampledDL6,InvertsByDischargeSampledDL7, weights=TRUE)
summary(InvertsByDischargeSampledDL3) #AIC weight = 0.68




intervals(InvertsByDischargeSampledDL3)



summary(InvertsByDischargeSampledDL2) #AIC weight = 0.188
summary(InvertsByDischargeSampledDL6) #phi1=0.597
summary(InvertsByDischargeSampledDL0) #phi1=0.727
summary(InvertsByDischargeSampledDL5) #phi1=0.592
summary(InvertsByDischargeSampledDL4) #phi1=0.658
summary(InvertsByDischargeSampledDL1) #phi1=0.723
summary(InvertsByDischargeSampledDL7) #phi1=0.594



###############
#Invert abundance Effects of discharge
################

ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
DischargeSubset<-subset(ShannonSubset, Q!="NA")
DischargeSubset$QCentered = DischargeSubset$Q - mean(DischargeSubset$Q)
head(DischargeSubset)

DischargeModel = lme(log(Ninverts100)~QCentered,random=~1|Year,data=DischargeSubset,correlation=corAR1(form = ~DPFS|Year))
summary(DischargeModel)

intervals(DischargeModel,which=c("fixed"))





DischargeModel = lme(log(Ninverts100)~QCentered+DayOfYear,random=~1|Year,data=DischargeSubset,correlation=corAR1(form = ~DPFS|Year))
summary(DischargeModel)

DischargeSampledSubset<-subset(ShannonSubset, DischargeSampledByNight!="NA")
DischargeSampledSubset$DSPerNightCentered = DischargeSampledSubset$DischargeSampledByNight - mean(DischargeSampledSubset$DischargeSampledByNight)


DischargeModel2 = lme(log(Ninverts100)~DSPerNightCentered,random=~1|Year,data=DischargeSampledSubset,correlation=corAR1(form = ~DPFS|Year))
summary(DischargeModel2)



intervals(DischargeModel2,which=c("fixed"))




###############
#Best LMM Invertebrate Abundance Model
#########

plot(resid(InvertsByDischargeSampledDL3) ~ temp_centered, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$temp_centered), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ percillum, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$percillum), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ DayOfYear, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$DayOfYear), col=2)


hist(resid(InvertsByDischargeSampledDL3))
plot(resid(InvertsByDischargeSampledDL3) ~ DischargeSampledByNight, data=ShannonSubset)
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ ShannonSubset$DischargeSampledByNight), col=2)

plot(resid(InvertsByDischargeSampledDL3) ~ predict(InvertsByDischargeSampledDL12))
lines(lowess(resid(InvertsByDischargeSampledDL3) ~ predict(InvertsByDischargeSampledDL12)), col=2)

summary(InvertsByDischargeSampledDL3)
confint(InvertsByDischargeSampledDL3)
fixef(InvertsByDischargeSampledDL3)
ranef(InvertsByDischargeSampledDL3)
coef(InvertsByDischargeSampledDL3)

#Days Post First Spawning
summary(InvertsByDischargeSampledDL3)
#FixedEffect DPFSCentered = -0.0325416
#Intercept   2.2887841

InvertsByDischargeSampledDL3 = lme(log(DriftInvertConc)~DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
intervals(InvertsByDischargeSampledDL3,which=c("fixed"))


#Lower est DPFS -0.04604928
#Upper est DPFS -0.01903387


InvertsByDischargeSampledDL3NoCor = lmer(log(DriftInvertConc)~DayOfYear+(1|Year),data=ShannonSubset)
confint(InvertsByDischargeSampledDL3NoCor)

samps<-sim(InvertsByDischargeSampledDL3NoCor,n.sim=1000)

(quantile(samps@fixef[,'DayOfYear'],c(0.025,0.975)))

head(samps@fixef)

m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"DayOfYear"]*X
new.DPFS = seq(min(ShannonSubset$DayOfYear),max(ShannonSubset$DayOfYear),length= 10000)
m.fixef.out<-sapply(new.DPFS,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)

head(m.pred)
PredictedDataFrame<-data.frame(exp(m.pred),new.DPFS,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("DriftInvertConc","DayOfYear","UpperCI","LowerCI")
#PredictedDataFrame$DPFS<-PredictedDataFrame$DayOfYear+mean(ShannonSubset$DayOfYear) #Remove centering on DPFS for plotting
head(PredictedDataFrame)

#ggplot
DPFSBestInvertModel<-ggplot(ShannonSubset,aes(DayOfYear,DriftInvertConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Calendar Date")+ylab(expression(Invertebrates~Per~100~m^3~Drift))
DPFSBestInvertModel

dev.off()
tiff("Figures/LMM_DPFS_BestNInvertModel.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DPFSBestInvertModel
dev.off()
###############
#Parameter estimates DPFS model invert abundance
###################
#Top model
summary(InvertsByDischargeSampledDL3)
intervals(InvertsByDischargeSampledDL3,which=c("fixed"))

(exp(-0.03254157)-1)*100 #-3.20178% estimate DPFS
(exp(-0.04604928)-1)*100 #-4.501% CI
(exp(-0.01903387)-1)*100 #-1.885% CI

#Second best model
summary(InvertsByDischargeSampledDL2)
intervals(InvertsByDischargeSampledDL2,which=c("fixed"))

(exp(-0.1212155)-1)*100 #-11.415% estimate temp centered
(exp(-0.1926145)-1)*100 #-17.52% CI
(exp(-0.04981647)-1)*100 #-4.8595% CI



######
#InvertsByDischargePlotWith CI bootstrapping
######

InvertsByDischargeSampledDL3 = lme(log(DriftInvertConc)~DoYCentered,random=~1|Year,data=ShannonSubset,)#correlation=corAR1(form = ~DayOfYear|Year)

xvals <-  with(ShannonSubset,seq(min(DoYCentered),max(DoYCentered),length.out=100))
nresamp<-1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked <- mvrnorm(nresamp, mu = fixef(InvertsByDischargeSampledDL3), Sigma = vcov(InvertsByDischargeSampledDL3))
head(pars.picked)
#pars.picked[,1]
## predicted values: useful below
pframe <- with(ShannonSubset,data.frame(DoYCentered=xvals))
#pframe
pframe$DriftInvertConc <- predict(InvertsByDischargeSampledDL3,newdata=pframe,level=0)

## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}
head(pars.picked)
head(xvals)
set.seed(101)




head(yvals)
yvals <- apply(pars.picked,1,
               function(x) { SSasymp(xvals,x[1], x[2], x[3]) }
)

c1 <- get_CI(yvals)

################3
#Invert Abu Lmm temp plot
################

InvertsByDischargeSampledDL3NoCor = lmer(log(DriftInvertConc)~temp_centered+(1|Year),data=ShannonSubset)
confint(InvertsByDischargeSampledDL3NoCor)

samps<-sim(InvertsByDischargeSampledDL3NoCor,n.sim=1000)

(quantile(samps@fixef[,'temp_centered'],c(0.025,0.975)))

head(samps@fixef)

m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"temp_centered"]*X
new.Temp = seq(min(ShannonSubset$temp_centered),max(ShannonSubset$temp_centered),length= 10000)
m.fixef.out<-sapply(new.Temp,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)

head(m.pred)
PredictedDataFrame<-data.frame(exp(m.pred),new.Temp,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("DriftInvertConc","Temp","UpperCI","LowerCI")
PredictedDataFrame$Temp<-PredictedDataFrame$Temp+mean(ShannonSubset$Temp) #Remove centering on Temp for plotting
head(PredictedDataFrame)

#ggplot
TempBestInvertModel<-ggplot(ShannonSubset,aes(Temp,DriftInvertConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Temperature (°C)")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift))
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

TempBestInvertModel

dev.off()
tiff("Figures/LMM_Temp_BestNInvertModel.tiff", width = 84, height = 84, units = 'mm', res = 1000)
TempBestInvertModel
dev.off()

##########
#Join together Figure 2
##########

dev.off()
tiff("Figures/Figure2InvertAbundance.tiff", width = 174, height = 174, units = 'mm', res = 1200)
ggarrange(DPFSBestInvertModel,TempBestInvertModel,TotalInvertAbuMoonPhase,HourlyInvertAbu,
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


BiomassLog0 = lme(log(DriftBiomassConc)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

BiomassLog1 = lme(log(DriftBiomassConc)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
BiomassLog2 = lme(log(DriftBiomassConc)~DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
BiomassLog3 = lme(log(DriftBiomassConc)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

BiomassLog4 = lme(log(DriftBiomassConc)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
BiomassLog5 = lme(log(DriftBiomassConc)~percillum+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
BiomassLog6 = lme(log(DriftBiomassConc)~temp_centered+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

BiomassLog7 = lme(log(DriftBiomassConc)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))


AICctab(BiomassLog0,BiomassLog1,BiomassLog2,BiomassLog3,BiomassLog4,BiomassLog5,BiomassLog6,BiomassLog7, weights=TRUE) 
summary(BiomassLog2)
summary(BiomassLog3)



intervals(BiomassLog2,which=c("fixed"))
(exp(-0.0454973)-1)*100 #-4.45% estimate DPFS
(exp(-0.06147474)-1)*100 #-5.96 CI
(exp(-0.02951985)-1)*100 #-2.91% CI



summary(BiomassLog6)

intervals(BiomassLog6,which=c("fixed"))
(exp(-0.03705362)-1)*100 #-3.63755 estimate DPFS
(exp(-0.05652343)-1)*100 #-5.4955% CI
(exp(-0.01758381)-1)*100 #-1.7739% CI

intervals(BiomassLog6,which=c("fixed"))
(exp(-0.06976298)-1)*100 #-6.7385 estimate temp
(exp(-0.16277468)-1)*100 #-15.02% CI
(exp(0.02324872)-1)*100 #2.352% CI

summary(BiomassLog3) #t=-3.40
intervals(BiomassLog3,which=c("fixed"))
(exp(-0.1497092)-1)*100 #-13.904% estimate temp only
(exp(-0.2365605)-1)*100 #-21.07 CI
(exp(-0.06285792)-1)*100 #-6.09% CI


cor(ShannonSubset$DriftInvertConc,ShannonSubset$DriftBiomassConc) #Drift biomass and total inverts are highly correlated 0.9647
cor(ShannonSubset$temp_centered,ShannonSubset$DPFSCentered) #Temp and discharge are relatively highlty correlated (0.4933)
######################
#Biomass LMM plot
#####################
BiomassLog2NoCor = lmer(log(DriftBiomassConc)~DPFSCentered+(1|Year),data=ShannonSubset)

samps <-sim(BiomassLog2NoCor,n.sims=10000)
(quantile(samps@fixef[,'DPFSCentered'],c(0.025,0.975)))

head(samps@fixef)

m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"DPFSCentered"]*X
new.DPFS = seq(min(ShannonSubset$DPFSCentered),max(ShannonSubset$DPFSCentered),length= 10000)
m.fixef.out<-sapply(new.DPFS,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)


PredictedDataFrame<-data.frame(exp(m.pred),new.DPFS,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("DriftBiomassConc","DPFS","UpperCI","LowerCI")
head(PredictedDataFrame)

PredictedDataFrame$DPFS<- PredictedDataFrame$DPFS+mean(ShannonSubset$DPFS)#Remove centering of DPFS
head(PredictedDataFrame)

#ggplot
DPFSBestBiomassModel<-ggplot(ShannonSubset,aes(DPFS,DriftBiomassConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Days Post First Spawning")+ylab(expression(Invertebrate~Biomass~(g)~Per~100~m^3~Drift))+ theme(axis.title.y = element_text(size = 9))
DPFSBestBiomassModel
dev.off()
tiff("Figures/LMM_DPFS_BestBiomassModel.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DPFSBestBiomassModel
dev.off()

###############
#Discharge effects on invert biomass
###############
DischargeSubset<-subset(ShannonSubset, Q!="NA")
DischargeSubset$QCentered = DischargeSubset$Q - mean(DischargeSubset$Q)

DischargeModel = lme(log(InvertBiomass100)~QCentered,random=~1|Year,data=DischargeSubset,correlation=corAR1(form = ~DPFS|Year))
summary(DischargeModel)



intervals(DischargeModel,which=c("fixed"))
(exp(0.04696164)-1)*100 #4.81 estimate River discharge
(exp(-0.02652691)-1)*100 #-2.617% CI
(exp(0.1204502)-1)*100 #12.80% CI


DischargeModel = lme(log(InvertBiomass100)~DischargeCentered,random=~1|Year,data=DischargeSubset,correlation=corAR1(form = ~DPFS|Year))
summary(DischargeModel)











############
#Best LMM model biomass diagnostic plots
###########



plot(resid(BiomassLog2) ~ DPFS, data=ShannonSubset)
lines(lowess(resid(BiomassLog2) ~ ShannonSubset$DPFS), col=2)


plot(resid(BiomassLog2) ~ temp_centered, data=ShannonSubset)
lines(lowess(resid(BiomassLog2) ~ ShannonSubset$temp_centered), col=2)


hist(resid(BiomassLog2))
plot(resid(BiomassLog2) ~ DischargeSampledByNight, data=ShannonSubset)
lines(lowess(resid(BiomassLog2) ~ ShannonSubset$DischargeSampledByNight), col=2)

plot(resid(BiomassLog2) ~ predict(BiomassLog2))
lines(lowess(resid(BiomassLog2) ~ predict(BiomassLog2)), col=2)

summary(BiomassLog2)
confint(BiomassLog2)
fixef(BiomassLog2)
ranef(BiomassLog2)


############
#Join Together Figure 4
############

DPFSBestBiomassModel
DriftTotalInvertBiomass
DPFSBestRichnessModel
library(gridGraphics)

dev.off()
tiff("Figures/Figure4Biomass.tiff", width = 6.85, height = 6.85, units = 'in', res = 1000)
ggarrange(DPFSBestBiomassModel,DriftTotalInvertBiomass,DPFSBestRichnessModel,p1,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()



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
CrayfishLog0 = lme(log(CrayfishPer100+10e-5)~1,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))

CrayfishLog1 = lme(log(CrayfishPer100+10e-5)~percillum,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))

CrayfishLog2 = lme(log(CrayfishPer100+10e-5)~DPFSCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
# CrayfishLog2Poly = lme((CrayfishPer100)~poly(DPFS,2),random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
# summary(CrayfishLog2Poly)



CrayfishLog3 = lme(log(CrayfishPer100+10e-5)~temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))

CrayfishLog4 = lme(log(CrayfishPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
CrayfishLog5 = lme(log(CrayfishPer100+10e-5)~percillum+DPFSCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
CrayfishLog6 = lme(log(CrayfishPer100+10e-5)~temp_centered+DPFSCentered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))

CrayfishLog7 = lme(log(CrayfishPer100+10e-5)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
AICctab(CrayfishLog0,CrayfishLog1,CrayfishLog2,CrayfishLog3,CrayfishLog4,CrayfishLog5,CrayfishLog6,CrayfishLog7, weights=TRUE)

plot(resid(CrayfishLog2) ~ DPFS, data=CrayfishAbu)
lines(lowess(resid(CrayfishLog2) ~ CrayfishAbu$DPFS), col=2)

plot(resid(CrayfishLog2) ~ predict(CrayfishLog2Poly))
lines(lowess(resid(CrayfishLog2) ~ predict(CrayfishLog2)), col=2)


CrayfishGAM<-ggplot(CrayfishAbu, aes(DPFS, CrayfishPer100) ) + geom_point(color="grey") +stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE,color="black")+
  stat_smooth(method = "gam", formula = y ~ s(x),geom="ribbon",linetype="dashed",fill=NA,color="black",level = 0.95)+ylab(expression(Crayfish~Per~100~m^3~Drift))

CrayfishGAM

dev.off()
tiff("Figures/LMM_DPFSCrayfishGAM.tiff", width = 84, height = 84, units = 'mm', res = 600)
CrayfishGAM
dev.off()

head(CrayfishAbu)

library(tidymv)
library(mgcv)
gam_y <- gam(log(CrayfishPer100+10e-5) ~ s(DPFS)+s(Year,bs="re")+s(DischargeCentered)+s(percillum), method = "REML",data=CrayfishAbu)
#gam.check(gam_y)
summary(gam_y)

library(mgcViz)
plot(gam_y)

gam_yViz <- getViz(gam_y)
o <- plot( sm(gam_yViz, 1) )
o+ l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

qq(gam_yViz, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))

m2 = glmer(Abundance~DPFSCentered+(1|Year),data=CrayfishAbu,family = poisson)

model_p <- predict_gam(gam_y)
model_p

model_p %>%
  ggplot(aes(DPFS, fit))

check(gam_yViz,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))




CrayfishAbu$PresenceSuckers <- 0
CrayfishAbu$PresenceSturgeon <- 0


CrayfishAbu$PresenceSuckers[CrayfishAbu$Nsuckers100 > 0] <- 1
CrayfishAbu$PresenceSturgeon[CrayfishAbu$NSturgeon > 0] <- 1
CrayfishSucker = lme(log(CrayfishPer100+10e-5)~SuckerPer100,random=~1|Year,data=CrayfishAbu,correlation=corAR1(form = ~DPFS|Year))
summary(CrayfishSucker)
#
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
head(FamilyLMMData)
IsoAbu<-subset(FamilyLMMData,Family=="Isonychiidae")

IsoAbu$Abundance<-IsoAbu$Abundance*20 #Account for 5% subsample

hist(IsoAbu$Abundance)
hist(log(IsoAbu$Abundance))
IsoAbu$IsoPer100<-((IsoAbu$Abundance*100)/(60*4*60*IsoAbu$AreaSampled.m2.*IsoAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
hist(IsoAbu$IsoPer100)
hist(log(IsoAbu$IsoPer100))


IsoLog0 = lme(log(IsoPer100+10e-5)~1,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))

IsoLog1 = lme(log(IsoPer100+10e-5)~percillum,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))
IsoLog2 = lme(log(IsoPer100+10e-5)~DPFSCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))
IsoLog3 = lme(log(IsoPer100+10e-5)~temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))

IsoLog4 = lme(log(IsoPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))
IsoLog5 = lme(log(IsoPer100+10e-5)~percillum+DPFSCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))
IsoLog6 = lme(log(IsoPer100+10e-5)~temp_centered+DPFSCentered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))

IsoLog7 = lme(log(IsoPer100+10e-5)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=IsoAbu,correlation=corAR1(form = ~DPFS|Year))


AICctab(IsoLog0,IsoLog1,IsoLog2,IsoLog3,IsoLog4,IsoLog5,IsoLog6,IsoLog7, weights=TRUE)


summary(IsoLog2)
summary(IsoLog1)



# summary(L6)
# summary(n6)

plot(resid(IsoLog2) ~ DPFS, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$DPFS), col=2)



plot(resid(IsoLog2) ~ temp_centered, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$temp_centered), col=2)

plot(resid(IsoLog2) ~ percillum, data=IsoAbu)
lines(lowess(resid(IsoLog2) ~ IsoAbu$percillum), col=2)

intervals(IsoLog2,which=c("fixed"))
(exp(-0.08181913)-1)*100 #-7.85614% estimate  DPFS

(exp(-0.1138543)-1)*100 #-10.7612% CI
(exp(-0.04978396)-1)*100 #-4.856505% CI 
#DPFS

IsoLog2NoCor = lmer(log(IsoPer100+10e-5)~DPFSCentered+(1|Year),data=IsoAbu)

samps <-sim(IsoLog2NoCor,n.sims=10000)

m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"DPFSCentered"]*X
new.DPFS = seq(min(IsoAbu$DPFSCentered),max(IsoAbu$DPFSCentered),length= 1000)
m.fixef.out<-sapply(new.DPFS,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)

head(m.pred)
PredictedDataFrame<-data.frame(exp(m.pred),new.DPFS,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("IsoPer100","DPFS","UpperCI","LowerCI")

PredictedDataFrame$DPFS<-PredictedDataFrame$DPFS+mean(IsoAbu$DPFS) #Remove centering for plotting
#Remove the +1 from the effect for plotting
# PredictedDataFrame$IsoPer100<-PredictedDataFrame$IsoPer100-1
# PredictedDataFrame$UpperCI<-PredictedDataFrame$UpperCI-1
# PredictedDataFrame$LowerCI<-PredictedDataFrame$LowerCI-1
head(PredictedDataFrame)
head(IsoAbu)
#ggplot
DPFSBestIsoModel<-ggplot(IsoAbu,aes(DPFS,IsoPer100))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Days post first spawning")+ylab(expression(Isonychiidae~Per~100~m^3~Drift))
DPFSBestIsoModel
dev.off()
tiff("Figures/LMM_DPFSIsoModel.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DPFSBestIsoModel
dev.off()





##############
#Heptageniidae model
############
HepAbu<-subset(FamilyLMMData,Family=="Heptageniidae")
head(HepAbu)
hist(HepAbu$Abundance)
HepAbu$HepPer100<-((HepAbu$Abundance*100)/(60*4*60*HepAbu$AreaSampled.m2.*HepAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
hist(log(HepAbu$Abundance))


HepLog0 = lme(log(HepPer100+10e-5)~1,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))

HepLog1 = lme(log(HepPer100+10e-5)~percillum,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))
HepLog2 = lme(log(HepPer100+10e-5)~DPFSCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))
HepLog3 = lme(log(HepPer100+10e-5)~temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))

HepLog4 = lme(log(HepPer100+10e-5)~percillum+temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))
HepLog5 = lme(log(HepPer100+10e-5)~percillum+DPFSCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))
HepLog6 = lme(log(HepPer100+10e-5)~temp_centered+DPFSCentered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))

HepLog7 = lme(log(HepPer100+10e-5)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=HepAbu,correlation=corAR1(form = ~DPFS|Year))



AICctab(HepLog0,HepLog1,HepLog2,HepLog3,HepLog4,HepLog5,HepLog6,HepLog7, weights=TRUE)


summary(HepLog2)


plot(resid(HepLog2) ~ DPFS, data=HepAbu)
lines(lowess(resid(HepLog2) ~ HepAbu$DPFS), col=2)


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
(exp(-0.05397137)-1)*100 #-5.254% estimate  DPFS

(exp(-0.07527147)-1)*100 #-7.2508% CI
(exp(-0.03267128)-1)*100 #-3.2143% CI 



############
#LMM models Shannon
#############
head(ShannonSubset)

ShannonLog0 = lmer(log(value)~1+(1|Year),data=ShannonSubset,REML=F)
ShannonLog1 = lmer(log(value)~percillum+(1|Year),data=ShannonSubset,REML=F)
ShannonLog2 = lmer(log(value)~temp_centered+(1|Year),data=ShannonSubset,REML=F)
ShannonLog3 = lmer(log(value)~DPFSCentered+(1|Year),data=ShannonSubset,REML=F)

ShannonLog4 = lmer(log(value)~percillum+temp_centered+(1|Year),data=ShannonSubset,REML=F)
ShannonLog5 = lmer(log(value)~percillum+DPFSCentered+(1|Year),data=ShannonSubset,REML=F)
ShannonLog6 = lmer(log(value)~temp_centered+DPFSCentered+(1|Year),data=ShannonSubset,REML=F)

ShannonLog7 = lmer(log(value)~percillum+DPFSCentered+temp_centered+(1|Year),data=ShannonSubset,REML=F)

AICctab(ShannonLog0,ShannonLog1,ShannonLog2,ShannonLog3,ShannonLog4,ShannonLog5,ShannonLog6,ShannonLog7, weights=TRUE)
summary(ShannonLog2) #Top shannon model not including autocorrelation AIC =71.9


ShannonLog0A = lme(log(value)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

ShannonLog1A = lme(log(value)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog2A = lme(log(value)~DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog3A = lme(log(value)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

ShannonLog4A = lme(log(value)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog5A = lme(log(value)~percillum+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog6A = lme(log(value)~temp_centered+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

ShannonLog7A = lme(log(value)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog8A = lme(log(value)~DischargeCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog9A = lme(log(value)~DischargeCentered+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog10A = lme(log(value)~DischargeCentered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog11A = lme(log(value)~DischargeCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
ShannonLog12A = lme(log(value)~DischargeCentered+temp_centered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))


AICctab(ShannonLog0A,ShannonLog1A,ShannonLog2A,ShannonLog3A,ShannonLog4A,ShannonLog5A,ShannonLog6A,ShannonLog7A,ShannonLog8A,ShannonLog9A,ShannonLog10A,ShannonLog11A,ShannonLog12A, weights=TRUE)


summary(ShannonLog8A)
summary(ShannonLog0A)


intervals(ShannonLog8A,which=c("fixed"))
(exp(0.4876861)-1)*100 #62.85% estimate 
(exp(0.15201898)-1)*100 #16.418% CI 
(exp(0.8233531)-1)*100 #127.81% CI

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
#Shannon LMM Plot
#############
ShannonLog8NoCorr = lmer(log(value)~DischargeCentered+ (1|Year),data=ShannonSubset)

summary(ShannonLog8NoCorr)
samps <-sim(ShannonLog8NoCorr,n.sims=10000)
#Temperature graph
m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"DischargeCentered"]*X
new.discharge = seq(min(ShannonSubset$DischargeCentered),max(ShannonSubset$DischargeCentered),length= 10000)
m.fixef.out<-sapply(new.discharge,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)

PredictedDataFrame<-data.frame(exp(m.pred),new.discharge,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("value","DischargeSampledByNight","UpperCI","LowerCI")
head(PredictedDataFrame)
PredictedDataFrame$DischargeSampledByNight<-PredictedDataFrame$DischargeSampledByNight+mean(ShannonSubset$DischargeSampledByNight)



#ggplot
DischargeShannonLMM<-ggplot(ShannonSubset,aes(DischargeSampledByNight,value))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab(expression(Discharge~(m^3/sec)~Sampled))+ylab("Shannon Diversity")#+ theme(axis.title.y = element_text(size = 9))
DischargeShannonLMM

dev.off()
tiff("Figures/LMM_ShannonDischarge.tiff", width = 84, height = 84, units = 'mm', res = 1000)
DischargeShannonLMM
dev.off()

####################
#Invertebrate family richness
####################
RichnessLog0A = lme(log(Nfamilies)~1,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

RichnessLog1A = lme(log(Nfamilies)~percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog2A = lme(log(Nfamilies)~DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog3A = lme(log(Nfamilies)~temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

RichnessLog4A = lme(log(Nfamilies)~percillum+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog5A = lme(log(Nfamilies)~percillum+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog6A = lme(log(Nfamilies)~temp_centered+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

RichnessLog7A = lme(log(Nfamilies)~percillum+DPFSCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog8A = lme(log(Nfamilies)~DischargeCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog9A = lme(log(Nfamilies)~DischargeCentered+DPFSCentered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog10A = lme(log(Nfamilies)~DischargeCentered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog11A = lme(log(Nfamilies)~DischargeCentered+temp_centered,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))
RichnessLog12A = lme(log(Nfamilies)~DischargeCentered+temp_centered+percillum,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))


AICctab(RichnessLog0A,RichnessLog1A,RichnessLog2A,RichnessLog3A,RichnessLog4A,RichnessLog5A,RichnessLog6A,RichnessLog7A,RichnessLog8A,RichnessLog9A,RichnessLog10A,RichnessLog11A,RichnessLog12A, weights=TRUE)


summary(RichnessLog2A)

intervals(RichnessLog2A,which=c("fixed"))

#Percentages
(exp(-0.01022522)-1)*100 #-1.0277% estimate 
(exp(-0.01509176)-1)*100 #-1.5206% estimate 
(exp(-0.005358687)-1)*100 #-0.5373% estimate 



#Day 15
min(ShannonSubset$DPFS)
mean(ShannonSubset$DPFS) #-21.98624 to account for centering 
exp(2.16784519+(-21.98624*-0.01022522)) #N families Estimate at 15 days 10.94
exp(2.16784519+(-21.98624*-0.01509176)) #Upper CI at 15 days 12.18
exp(2.16784519+(-21.98624*-0.005358687)) #Lower CI at 15 days 9.83


#Day 68
max(ShannonSubset$DPFS) #68
60-36.98624
exp(2.16784519+(23.01376*-0.01022522)) #N families Estimate at 60 days 6.906
exp(2.16784519+(23.01376*-0.01509176)) #Upper CI at 60 days 6.175
exp(2.16784519+(23.01376*-0.005358687)) #Lower CI at 60 days 7.725

###############
#LMM Richness plot DPFS
#################

RichnessLog2NoCorr = lmer(log(Nfamilies)~DPFSCentered+(1|Year),data=ShannonSubset)

samps <-sim(RichnessLog2NoCorr,n.sims=10000)

m.fixef <-function(X) samps@fixef[,'(Intercept)']+samps@fixef[,"DPFSCentered"]*X
new.DPFS = seq(min(ShannonSubset$DPFSCentered),max(ShannonSubset$DPFSCentered),length= 10000)
m.fixef.out<-sapply(new.DPFS,m.fixef)

m.pred<-colMeans(m.fixef.out)
m.975<- apply(m.fixef.out,2,quantile, 0.975)
m.025 <-apply(m.fixef.out, 2 ,quantile, 0.025)


PredictedDataFrame<-data.frame(exp(m.pred),new.DPFS,exp(m.975),exp(m.025))
colnames(PredictedDataFrame)<-c("Nfamilies","DPFS","UpperCI","LowerCI")
PredictedDataFrame$DPFS<-PredictedDataFrame$DPFS+mean(ShannonSubset$DPFS) #Remove centering for plotting


head(PredictedDataFrame)
head(ShannonSubset)
#ggplot
DPFSBestRichnessModel<-ggplot(ShannonSubset,aes(DPFS,Nfamilies))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Days Post First Spawning")+ylab("Family Richness")
DPFSBestRichnessModel

dev.off()
tiff("Figures/LMM_NfamiliesDPFS.tiff", width = 84, height = 84, units = 'mm', res = 800)
DPFSBestRichnessModel
dev.off()



##########
#Join together figure 5
##########

dev.off()
tiff("Figures/Figure5.tiff", width = 84, height = 174, units = 'mm', res = 1200)
ggarrange(DischargeShannonLMM,TemperatureShannonLMM,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
dev.off()





#############
#Autocorrelation
############
Data2011<-subset(ShannonSubset,Year==2011)
head(Data2011)
model<-lmer(log(DriftInvertConc)~DPFSCentered+(1|Year),data=Data2011)

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2011$DPFSCentered)

Data2012<-subset(ShannonSubset,Year==2012)
head(Data2012)
model<-glm(DriftInvertConc~DPFSCentered,data=Data2012,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2012$DPFSCentered)

Data2013<-subset(ShannonSubset,Year==2013)
head(Data2013)
model<-glm(DriftInvertConc~DPFSCentered,data=Data2013,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2013$DPFSCentered)

Data2014<-subset(ShannonSubset,Year==2014)
head(Data2014)
model<-glm(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale,data=Data2014,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2014$DPFS)

Data2015<-subset(ShannonSubset,Year==2015)
head(Data2015)
model<-glm(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale,data=Data2015,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2015$DPFS)

Data2016<-subset(ShannonSubset,Year==2016)
head(Data2016)
model<-glm(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale,data=Data2016,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2016$DPFS)


Data2017<-subset(ShannonSubset,Year==2017)
head(Data2017)
model<-glm(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale,data=Data2017,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2017$DPFS)



Data2018<-subset(ShannonSubset,Year==2018)
head(Data2018)
model<-glm(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale,data=Data2018,family = Gamma(link=log))

res<-simulateResiduals(model)
testTemporalAutocorrelation(res, time =  Data2018$DPFS)


acf(resid(BiomassLog13), main="")
acf(resid(n13), plot=FALSE)$acf[2]




library(nlme)

BiomassLog13 = lmer(log(DriftInvertConc)~percillum+DPFSCentered+temp_centered+(1|Year),data=ShannonSubset,REML=F)

BiomassLog3Auto = lme(log(DriftInvertConc)~DPFS,random=~1|Year,data=ShannonSubset,correlation=corAR1(form = ~DPFS|Year))

summary(BiomassLog3Auto)#
AICtab(BiomassLog13,BiomassLog3Auto)
hist(resid(BiomassLog13Auto))



acf(resid(BiomassLog13))
pacf(resid(BiomassLog13Auto))

model<-glmmTMB(DriftInvertConc~DischargeByNightScale+percillumScale+DPFSScale+(1|Year),data=simdat,family = Gamma(link=log))
acf(resid(BiomassLog13), main="acf(resid(m1))")



################
#Hourly results
###############
Hourly<-read.csv("DataClean\\SturgeonDriftDataHourly12.21.19.csv",header=T)
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
Model<-aov(log(DriftInvertConc)~Time,data=StatSubset)

summary(Model)
plot(Model)
Tukey<-TukeyHSD(Model,"Time")
plot(Tukey)
Tukey
(exp(1.510104612)-1)*100 #352.7%Difference between 10 pm collection and midnight
(exp(1.1965172)-1)*100 #230.8574% lower CI
(exp(1.8236920)-1)*100 #519.4687% upper CI

difference<-Tukey$Time[,"p adj"]
Letters<-multcompLetters(difference)
Letters                         
#rearrange letters so a is at ten pm, c=a,a=b,c=c
LetterRearranged<-c("a","b","c","c","c")


hist(resid(Model))


Trtdata <- ddply(StatSubset, c("Time"), summarise,
                 N    = length(DriftInvertConc),
                 meanInverts = mean(DriftInvertConc),
                 sd   = sd(DriftInvertConc),
                 se   = sd / sqrt(N)
)
Trtdata

HourlyInvertAbu<-ggplot(Trtdata,aes(x=Time,y=meanInverts))+geom_point()+xlab("Collection Time")+
  geom_errorbar(aes(ymin=meanInverts-se,ymax=meanInverts+se))+ylab(expression(Invertebrates~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanInverts+se+0.5,label=LetterRearranged))+
  annotate("text", label = "ANOVA, F = 63.7,\n P < 0.001", size = 3.5, x = 1.5, y = 15)+ theme(axis.title.y = element_text(size = 10))
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
  geom_errorbar(aes(ymin=meanSucker-se,ymax=meanSucker+se))+ylab(expression(Catastomidae~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanSucker+se+5,label=LettersSucker$Letters))+
  annotate("text", label = "Kruskal-Wallis,\n chi-sq = 126.1,\n P < 0.001", size = 3.5, x = 4.5, y = 90)+ theme(axis.title.y = element_text(size = 10))
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

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersSturgeon<-multcompLetters(difference)
LettersSturgeon
TrtdataSturgeon
LettersSturgeon<-c("a","b","c","bc","bc")
HourlySturgeonAbu<-ggplot(TrtdataSturgeon,aes(x=Time,y=meanSturgeon))+geom_point()+xlab("Collection Time")+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ylab(expression(Sturgeon~larvae~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanSturgeon+se+0.5,label=LettersSturgeon))+
  annotate("text", label = "Kruskal-Wallis,\n chi-sq = 107.7,\n P < 0.001", size = 3.5, x = 2, y = 7.5)+ theme(axis.title.y = element_text(size = 10))
HourlySturgeonAbu


theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/Hourly_SturgeonAbundance.tiff", width = 84, height = 84, units = 'mm', res = 1200)
HourlySturgeonAbu
dev.off()


###########
#Join together hourly figure
###########
dev.off()
tiff("Figures/Hourly_Abundance.tiff", width = 174, height = 174, units = 'mm', res = 1200)

ggarrange(HourlyInvertAbu,HourlySuckerAbu,HourlySturgeonAbu,
          labels = c("a", "b","c"),
          ncol = 2, nrow = 2)

dev.off()

##################
#Sturgeon/Sucker Biomass
###################
ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)

head(ShannonRichness)
#Sturgeon average dry weight = 0.005 g
#Sucker average dry weight = 0.00204
ShannonRichness$SturgeonBiomass<-ShannonRichness$SturgeonConc*0.005

ShannonRichness$SuckerBiomass<-ShannonRichness$SuckerConc*0.00204

ShannonRichness$CombinedBiomassConc<-ShannonRichness$SturgeonBiomass+ShannonRichness$SuckerBiomass+ShannonRichness$DriftBiomassConc

#Sturgeon Biomass
mean(ShannonRichness$SturgeonBiomass,na.rm=T) #0.0167 g per 100 m3 drift
min(ShannonRichness$SturgeonBiomass,na.rm=T) #0
max(ShannonRichness$SturgeonBiomass,na.rm=T) #0.453


ShannonSubsetSampled<-subset(ShannonRichness,SturgeonBiomass!="NA")
length(ShannonSubsetSampled$SampleID)
sd(ShannonSubsetSampled$SturgeonBiomass)/sqrt(229) #0.003


#Sucker Biomass
mean(ShannonRichness$SuckerBiomass,na.rm=T) #0.111 g
min(ShannonRichness$SuckerBiomass,na.rm=T) #0
max(ShannonRichness$SuckerBiomass,na.rm=T) #3.09 g


ShannonSubsetSampled<-subset(ShannonRichness,SuckerBiomass!="NA")
length(ShannonSubsetSampled$SampleID)
sd(ShannonSubsetSampled$SuckerBiomass)/sqrt(229) #0.022



#InvertBiomass
mean(ShannonRichness$DriftBiomassConc,na.rm=T) #0.746 g
min(ShannonRichness$DriftBiomassConc,na.rm=T) #0.005297
max(ShannonRichness$DriftBiomassConc,na.rm=T) #6.781


ShannonSubsetSampled<-subset(ShannonRichness,DriftBiomassConc!="NA")
length(ShannonSubsetSampled$SampleID)
sd(ShannonSubsetSampled$DriftBiomassConc)/sqrt(229) #0.056




#Total Biomass Conc
mean(ShannonRichness$CombinedBiomassConc,na.rm=T) #0.874 g per 100 m3 drift
min(ShannonRichness$CombinedBiomassConc,na.rm=T) #0.0105
max(ShannonRichness$CombinedBiomassConc,na.rm=T) #7.688


#







#Extrapolating out to total river discharge

ShannonRichness$SuckerBiomassPerNight<-ShannonRichness$Nsuckers100*0.00204
ShannonRichness$SturgeonBiomassPerNight<-ShannonRichness$Nsturgeon*0.005

ShannonRichness$TotalSampledBiomassNight<-ShannonRichness$SuckerBiomassPerNight+ShannonRichness$SturgeonBiomassPerNight+ShannonRichness$InvertBiomass100
ShannonRichness$TotalRiverBiomassNightly<-ShannonRichness$TotalSampledBiomassNight/(ShannonRichness$PercentRiverDischargeSampled/100)
#Sampled Nightly Biomass
mean(ShannonRichness$TotalSampledBiomassNight,na.rm=T) #112.43 g per night sampled
min(ShannonRichness$TotalSampledBiomassNight,na.rm=T) #1.246
max(ShannonRichness$TotalSampledBiomassNight,na.rm=T) #940.43



ShannonSubsetSampled<-subset(ShannonRichness,TotalSampledBiomassNight!="NA")
length(ShannonSubsetSampled$SampleID)
sd(ShannonSubsetSampled$TotalSampledBiomassNight)/sqrt(240) #7.78




#River Biomass Nightly
mean(ShannonRichness$TotalRiverBiomassNightly,na.rm=T) #919.93 g per night sampled
min(ShannonRichness$TotalRiverBiomassNightly,na.rm=T) #13.31
max(ShannonRichness$TotalRiverBiomassNightly,na.rm=T) #7885.5

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
###


















SturgeonDataFrame<-data.frame("Sturgeon",ShannonRichness$DayOfYear, ShannonRichness$CTUMay2,ShannonRichness$Year,ShannonRichness$SturgeonBiomass)
colnames(SturgeonDataFrame)<-c("Taxa","DPFS","CTUMay2","Year","BiomassConc")

SuckerDataFrame<-data.frame("Catostomidae",ShannonRichness$DPFS, ShannonRichness$CTUMay2,ShannonRichness$Year,ShannonRichness$SuckerBiomass)
colnames(SuckerDataFrame)<-c("Taxa","DPFS","CTUMay2","Year","BiomassConc")

InvertDataFrame<-data.frame("Invertebrate",ShannonRichness$DPFS,ShannonRichness$CTUMay2,ShannonRichness$Year,ShannonRichness$DriftBiomassConc)
colnames(InvertDataFrame)<-c("Taxa","DPFS","CTUMay2","Year","BiomassConc")

PlottingDataframe<-rbind(SturgeonDataFrame,SuckerDataFrame,InvertDataFrame)
PlottingDataframe[is.na(PlottingDataframe)] <- 0 #Change NA values to 0 as no taxa were observed on those dates but sampling occured

head(PlottingDataframe)
PlottingDataframe$Taxa = factor(PlottingDataframe$Taxa, levels = c("Invertebrate","Sturgeon","Catostomidae"))

theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

BiomassAllPlot<-ggplot(PlottingDataframe, aes(x=DPFS, y= BiomassConc,color=Taxa,shape=Taxa,group=Taxa))+facet_wrap(~Year)+geom_point(size=1)+
  theme(legend.justification=c(1,0), legend.position=c(1,-0.01))+ylab(expression(Biomass~(g)~Per~100~m^3~Drift))+xlab("Days Post First Spawning")+
  scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+theme(legend.text = element_text(size = 7),legend.title = element_blank())+ 
  theme(legend.background=element_blank())+ guides(shape = guide_legend(override.aes = list(size=2)))#+geom_line()
BiomassAllPlot

BiomassAllPlot2<-ggplot(PlottingDataframe, aes(x=CTUMay2, y= BiomassConc,color=Taxa,shape=Taxa,group=Taxa))+facet_wrap(~Year)+geom_point(size=1)+
  theme(legend.justification=c(1,0), legend.position=c(1,-0.01),axis.text.x = element_text(angle = 45, hjust = 1))+ylab(expression(Biomass~(g)~Per~100~m^3~Drift))+xlab("Cumulative Thermal Units")+
  scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+theme(legend.text = element_text(size = 7),legend.title = element_blank())+ 
  theme(legend.background=element_blank())+ guides(shape = guide_legend(override.aes = list(size=2)))+ scale_x_continuous( limits=c(100, 1200))#+geom_line()
BiomassAllPlot2




dev.off()
tiff("Figures/BiomassAllPlot.tiff", width = 84, height = 84, units = 'mm', res = 1200)
BiomassAllPlot
dev.off()

dev.off()
tiff("Figures/BiomassAllPlot2.tiff", width = 84, height = 84, units = 'mm', res = 1200)
BiomassAllPlot2
dev.off()


############
#Combine fig 6
################
TopFamilyAbuPer100 #A
CrayfishGAM #B
TopFamilyBiomass #C
BiomassAllPlot #D

dev.off()
tiff("Figures/Figure6.tiff", width = 174, height = 174, units = 'mm', res = 1000)
ggarrange(TopFamilyAbuPer100,CrayfishGAM,TopFamilyBiomass,BiomassAllPlot,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()


#############
#Beta diversity exploration
##############


#Homogeneity of multivariate dispersions

GPdist=phyloseq::distance(physeq, "jaccard")
beta=betadisper(GPdist, sample_data(physeq)$Year)
permutest(beta)
boxplot(beta)


adonis(GPdist ~ Year, as(sample_data(physeq), "data.frame"))

ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="MoonPhase")+geom_point(size=4)+
  stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = MoonPhase))+facet_wrap(~Year)
ordplot


###############
#Nightly Total River Biomass
###############
head(metadata)

#############
#Sturgeon Abu Gam
################
#Model residuals not great (maybe zero inflated model would work?), not included
SturgeonAbuLog <- gam(log(Nsturgeon+10e-5)~ s(CTUSturgeon)+s(Year,bs="re")+s(DischargeSampled), data = ShannonRichness, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
appraise(SturgeonAbuLog)

summary(SturgeonAbuLog)
hist(resid(SturgeonAbuLog))

plot<-plot_smooth(SturgeonAbuLog, view="CTUSturgeon",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("Nsturgeon","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=Nsturgeon))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")


ggplot(ShannonRichness,aes(CTUSturgeon,Nsturgeon))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression("Drifting Larval Sturgeon Abundance"))






##############
#Scratch
##############
SuckerAbuLog <- gam(log(SuckerConc+10e-5)~ s(percillum)+s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

appraise(SuckerAbuLog)
plot<-plot_smooth(SuckerAbuLog, view="percillum",rm.ranef=F,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$percillum,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","percillum","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=percillum,y=SuckerConc))+geom_point()+geom_line(aes(percillum, UpperCI),color="red")+geom_line(aes(percillum,LowerCI),color="red")


SuckerGAMPlot<-ggplot(ShannonSubset,aes(percillum,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Larval~Catostomidae~Per~100~m^3~Drift))+ylim(NA, 500)
SuckerGAMPlot

dev.off()
tiff("Figures/SuckerGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SuckerGAMPlot
dev.off()

