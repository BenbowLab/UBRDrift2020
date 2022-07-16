#UBR R Code Hydrobiologia revisions (Last updated Mar 2022) 
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

library(itsadug)
library(gratia)

set.seed(45682)
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

metadata=read.csv("DataClean/SturgeonDriftMetadata10.28.2020.csv",header=TRUE)
metadata$InvertsByRiverDischarge<-metadata$Ninverts100/metadata$QManual
metadata$BiomassByRiverDischarge<- metadata$InvertBiomass100/metadata$QManual
metadata$DriftInvertConc<-((metadata$Ninverts100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
metadata$DriftBiomassConc<-((metadata$InvertBiomass100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water
metadata$SturgeonConc<-((metadata$Nsturgeon*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate sturgeon larvae/ 100 m3 water
metadata$SuckerConc<-((metadata$Nsuckers100*100)/(60*4*60*metadata$AreaSampled.m2.*metadata$AverageNetFlowByNight)) #Calculate sturgeon larvae/ 100 m3 water
metadata$DischargeSampledByNight<-metadata$AverageNetFlowByNight*2.145 #Flow x total area of nets (5 x 0.429 m2)
  
  
metadata$InvertsByDischargeSampled<-metadata$Ninverts100/metadata$DischargeSampledByNight
metadata$PercentRiverDischargeSampled<-metadata$DischargeSampledByNight/metadata$QManual*100
  
metadata$BiomassByDischargeSampled<-metadata$InvertBiomass100/metadata$DischargeSampledByNight
ShannonRichness<-metadata
# 
HistDischarge<-ggplot(metadata,aes(x=DischargeSampledByNight))+geom_histogram(fill="grey",color="black")+xlab(expression(Discharge~(m^3/sec)~Sampled))+ylab("Frequency")#+facet_wrap(~Year)
HistDischarge

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

min(metadata$QManual,na.rm=T)
max(metadata$QManual,na.rm=T)

mean(metadata$QManual,na.rm=T)
SD<-sd(metadata$QManual,na.rm=T)

N<-length(na.omit(metadata$QManual))
SD / sqrt(N)

min(metadata$DischargeSampledByNight,na.rm=T)
max(metadata$DischargeSampledByNight,na.rm=T)

mean(metadata$DischargeSampledByNight,na.rm=T)
SD<-sd(metadata$DischargeSampledByNight,na.rm=T)

N<-length(na.omit(metadata$DischargeSampledByNight))
SD / sqrt(N)
# Temperature
min(metadata$Temp,na.rm=T)
max(metadata$Temp,na.rm=T)

mean(metadata$Temp,na.rm=T)
SD<-sd(metadata$Temp,na.rm=T)

N<-length(na.omit(metadata$Temp))
SD / sqrt(N)



dev.off()
tiff("Figures/Discharge_Sampled.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
HistDischarge
dev.off()


HistQ<-ggplot(metadata,aes(x=QManual))+geom_histogram(fill="grey",color="black")+xlab(expression(River~Discharge~(m^3/sec)))+ylab("Frequency")#+facet_wrap(~Year)
HistQ


ggplot(metadata,aes(x=Q,y=DischargeSampledByNight))+geom_point()+xlab("River discharge (m3/sec)")+ylab("Average sampled discharge (m3/sec)")
# theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
# 
# dev.off()
# tiff("Figures/Discharge_Q.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
# HistQ
# dev.off()


ggplot(metadata,aes(x=Q,y=QManual))+geom_point()+xlab("River Discharge (HOBO)")+ylab("River Discharge (Manual)")+annotate("text",x=14,y=10,label="Pearson Correlation = 0.014")


ggplot(metadata,aes(x=DischargeSampledByNight,y=Qmanual))+geom_point()#+xlab("River Discharge (HOBO)")+ylab("River Discharge (Manual)")#+annotate("text",x=14,y=10,label="Pearson Correlation = 0.014")

DischargeSubset<-subset(metadata, DischargeSampledByNight!="NA" &Qmanual !="NA")

cor(DischargeSubset$DischargeSampledByNight,DischargeSubset$Qmanual)



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#sets the plotting theme

# ShannonRichness<-read.csv("SturgeonMetadataWDiversity.csv",header=T)
ShannonRichness<-metadata

AllData<-metadata
AllData$Date2<-as.Date(AllData$Date,format= "%d-%B")
#AllData

AllDataMay1toJuly15





#######
# CorrelationPlot
######
library(corrplot)

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




cor(ShannonSubset$DriftInvertConc,ShannonSubset$DriftBiomassConc)

CorrPlotDataframe<-data.frame(Temperature=ShannonSubset$Temp,DayOfYear=ShannonSubset$DayOfYear,CTU=ShannonSubset$CTUSturgeon,RiverDischarge=ShannonSubset$QManual,LunarIllumination=ShannonSubset$percillum)
CorrPlotDataframe<-subset(CorrPlotDataframe, RiverDischarge!="NA")
col_order <- c("Temperature", "DayOfYear", "RiverDischarge","LunarIllumination",
                "CTU")
CorrPlotDataframe2 <- CorrPlotDataframe[, col_order]


head(CorrPlotDataframe2)
M<-cor(CorrPlotDataframe2)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
)



plot(x=metadata$Q,y=metadata$DischargeSampledByNight, xlab="Mean River Discharge (m3/sec)",ylab="Mean Sampled Discharge (m3/sec)")


dev.off()
tiff("Figures/CorrPlot.tiff", width = 174, height = 174, units = 'mm', res = 1200)

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
)
dev.off()






####################
#Figure 1
########################

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

############
#Fig 1A Invert By day and year
############
Trtdata <- ddply(ShannonRichness, c("Date2","Year"), summarise,
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

Trtdata <- ddply(ShannonRichness, c("Date2","Year"), summarise,
                 N    = length(Nsuckers100),
                 meanSuckers = mean(Nsuckers100)
)

SuckersByDayOfYear<-ggplot(Trtdata, aes(x=Date2,y=meanSuckers))+geom_bar(colour="black", stat="identity")+xlab("Date")+ylab("Larval Sucker Abundance")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+facet_grid(Year~.)#+scale_fill_manual(values=cbPalette)
SuckersByDayOfYear

##########
# Figure 1C sturgeon larvae abu by day and year 
##########
Trtdata <- ddply(ShannonRichness, c("Date2","Year"), summarise,
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



###########
#Fig 1D Densities by CTU
##########
ShannonRichness$Date2<-as.Date(ShannonRichness$Date,format= "%d-%B")

DischargePlot<-ggplot(ShannonRichness,aes(x=Date2,y=QManual))+geom_bar(stat="identity",color="black")+geom_point(aes(y=Temp/1.75,col=Temp),size=0.5,data=ShannonRichness,color="red")+scale_y_continuous(sec.axis = sec_axis(~.*1.75,name="Water Temperature (\u00B0C)"))+xlab("Date")+facet_grid(Year~.)+ylab(expression(River~Discharge~(m^3/sec)))
DischargePlot


############
#Join together Figure 1
###########

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/Fig1Mar2022.tiff", width = 174, height = 190, units = 'mm', res = 600)
ggarrange(InvertByDayOfYear,SuckersByDayOfYear,SturgeonByDayOfYear,DischargePlot,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)
dev.off()






##############
# Sucker GAM
##############
#citation(gratia)


ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, DischargeSampledByNight!= "NA")
ShannonSubset$Year<-as.factor(ShannonSubset$Year)

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")
ShannonSubset<-subset(ShannonSubset, QManual!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)

ShannonSubset$Qcentered = ShannonSubset$QManual - mean(ShannonSubset$QManual)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

head(ShannonSubset)
# Null model
SuckerAbuLog0 <- gam(log(SuckerConc+0.001)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
SuckerAbuLog1 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog2 <- gam(log(SuckerConc+0.001)~ s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog3 <- gam(log(SuckerConc+0.001)~ (MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog4 <- gam(log(SuckerConc+0.001)~ s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
SuckerAbuLog5 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog6 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog7 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


SuckerAbuLog8 <- gam(log(SuckerConc+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog9 <- gam(log(SuckerConc+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SuckerAbuLog10 <- gam(log(SuckerConc+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
SuckerAbuLog11 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SuckerAbuLog12 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SuckerAbuLog13 <- gam(log(SuckerConc+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
SuckerAbuLog14 <- gam(log(SuckerConc+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))




AICctab(SuckerAbuLog0,SuckerAbuLog1,SuckerAbuLog2,SuckerAbuLog3,SuckerAbuLog4,SuckerAbuLog5,SuckerAbuLog6,SuckerAbuLog7,SuckerAbuLog8,
        SuckerAbuLog9,SuckerAbuLog10,SuckerAbuLog11,SuckerAbuLog12,SuckerAbuLog13,SuckerAbuLog14,
        weights=T)

concurvity(SuckerAbuLog12)
summary(SuckerAbuLog12)
summary(SuckerAbuLog7)

SuckerAbuLog12$aic

ShannonRichness$SturgeonBiomassNightly<-ShannonRichness$Nsturgeon*0.005 #Average sturgeon Nightly biomass (early and late spawning combined)
kruskal.test(SuckerConc~MoonPhase,data=ShannonRichness)

kruskal.test(SturgeonConc~MoonPhase,data=ShannonRichness)
compare_means(SturgeonConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

kruskal.test(DriftInvertConc~MoonPhase,data=ShannonRichness)
kruskal.test(DriftBiomassConc~MoonPhase,data=ShannonRichness)


ggplot(ShannonRichness,aes(x=MoonPhase,y=Nsturgeon))+geom_boxplot()+ylim(NA,50)
length(ShannonSubset$SampleID)
ggplot(ShannonRichness,aes(x=MoonPhase,y=DriftInvertConc))+geom_boxplot()+ylim(NA,50)

##########
#Top sucker model
##########

dev.off()
tiff("Figures/SuckerGAMAppraisal.tiff", width = 174, height = 174, units = 'mm', res = 1200)
appraise(SuckerAbuLog12)
dev.off()



summary(SuckerAbuLog12)
plot<-plot_smooth(SuckerAbuLog12, view="CTUSturgeon",rm.ranef=T,sim.ci = T)
plot_smooth(SuckerAbuLog12, view="CTUSturgeon",rm.ranef=T,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=SuckerConc))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")

head(ShannonSubset)
SuckerGAMPlot<-ggplot(ShannonSubset,aes(CTUSturgeon,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Larval~Suckers~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))+ylim(0, 1050)#+ylim(NA, 1000)
SuckerGAMPlot
subset(ShannonSubset, SuckerConc>1050)
# One samples was omitted from the plot for clarity (but included in all models), 1,515 Suckers per 100m3 at 169 CTU.

PredictedDataFramelog<-data.frame((plot$fv$fit),plot$fv$CTUSturgeon,(plot$fv$ul),(plot$fv$ll))
colnames(PredictedDataFramelog)<-c("SuckerConc","CTUSturgeon","UpperCI","LowerCI")

ggplot(ShannonSubset,aes(CTUSturgeon,log(SuckerConc)))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5,aes(y=log(SuckerConc)))+  
  geom_line(data = PredictedDataFrame, aes(y = log(LowerCI)), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=log(UpperCI)),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Log(Larval~Sucker~Per~100~m^3~Drift)))



dev.off()
tiff("Figures/SuckerGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SuckerGAMPlot
dev.off()

#Discharge
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")

plot<-plot_smooth(SuckerAbuLog12, view="QManual",rm.ranef=T,sim.ci = T,cond=list(CTUSturgeon = 300))
plot2<-plot_smooth(SuckerAbuLog12, view="CTUSturgeon",rm.ranef=F,sim.ci = T,plot_all = "Year",transform = exp, col=cbPalette,shade=TRUE)

tiff("Figures/SuckerRiverDischarge.tiff", width = 150, height = 150, units = 'mm', res = 1200)
plot_smooth(SuckerAbuLog12, view="QManual",rm.ranef=T,sim.ci = T,transform=exp,xlab = expression(River~Discharge~(m^3/sec)),ylab=expression(Suckers~Per~100~m^3~Drift,lwd=2))
dev.off()



plot_smooth(SuckerAbuLog12, view="CTUSturgeon",rm.ranef=F,sim.ci = T,plot_all = "Year",transform = exp, col=cbPalette,)

#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$QManual,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","QManual","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=QManual,y=SuckerConc))+geom_point()+geom_line(aes(QManual, UpperCI),color="red")+geom_line(aes(QManual,LowerCI),color="red")


SuckerGAMPlot2<-ggplot(ShannonSubset,aes(QManual,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("River Discharge (m^3/sec)")+ylab(expression(Larval~Suckers~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))+ylim(NA, 100)
SuckerGAMPlot2



dev.off()
tiff("Figures/SuckerGAM2.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SuckerGAMPlot2
dev.off()



plot<-plot_smooth(SuckerAbuLog12, view="Temp",rm.ranef=T,sim.ci = T,cond=list(CTUSturgeon = 300),
                  ylab = expression(log(Larval~Suckers~Per~100~m^3~Drift)),xlab="Water Temperature (\u00B0C)")

tiff("Figures/SuckerGAMTempPremade.tiff", width = 100, height = 100, units = 'mm', res = 600)
plot_smooth(SuckerAbuLog12, view="Temp",rm.ranef=T,sim.ci = T,cond=list(CTUSturgeon = 300),
            ylab = expression(Larval~Suckers~Per~100~m^3~Drift),xlab="Water Temperature (\u00B0C)",transform = exp)
dev.off()

#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$Temp,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SuckerConc","Temp","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=Temp,y=SuckerConc))+geom_point()+geom_line(aes(Temp, UpperCI),color="red")+geom_line(aes(Temp,LowerCI),color="red")


SuckerGAMPlot3<-ggplot(ShannonSubset,aes(Temp,SuckerConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Water Temperature (\u00B0C)")+ylab(expression(Larval~Suckers~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))+ylim(NA, 100)
SuckerGAMPlot3




# Combine Premade Sucker GAM

tiff("Figures/SuckerRiverGAMCombined.tiff", width = 160, height = 84, units = 'mm', res = 1200)
par(mfrow = c(1,2))
plot_smooth(SuckerAbuLog12, view="QManual",rm.ranef=T,sim.ci = T,transform=exp,xlab = "River Discharge (m^3/sec)",ylab=expression(Larval~Suckers~Per~100~m^3~Drift,lwd=2))

plot_smooth(SuckerAbuLog12, view="Temp",rm.ranef=T,sim.ci = T,cond=list(CTUSturgeon = 300),
            ylab = expression(Larval~Suckers~Per~100~m^3~Drift),xlab="Water Temperature (\u00B0C)",transform = exp)

dev.off()
par(mfrow = c(1,1))


par(mfrow = c(1,2))
plot_smooth(SuckerAbuLog12, view="QManual",rm.ranef=T,sim.ci = T,transform=exp,xlab = "River Discharge (m^3/sec)",ylab=expression(Suckers~Per~100~m^3~Drift,lwd=2))

plot_smooth(SuckerAbuLog12, view="Temp",rm.ranef=T,sim.ci = T,cond=list(CTUSturgeon = 300),
            ylab = expression(Larval~Suckers~Per~100~m^3~Drift),xlab="Water Temperature (\u00B0C)",transform = exp)





##############
#Sturgeon GAM
##############


ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset$Year<-as.factor(ShannonSubset$Year)
levels(ShannonSubset$Year)

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")
ShannonSubset<-subset(ShannonSubset, QManual!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)

ShannonSubset$Qcentered = ShannonSubset$QManual - mean(ShannonSubset$QManual)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)

# Null model
SturgeonAbuLog0 <- gam(log(SturgeonConc+0.001)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
SturgeonAbuLog1 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog2 <- gam(log(SturgeonConc+0.001)~ s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog3 <- gam(log(SturgeonConc+0.001)~ (MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog4 <- gam(log(SturgeonConc+0.001)~ s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
SturgeonAbuLog5 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog6 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog7 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


SturgeonAbuLog8 <- gam(log(SturgeonConc+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog9 <- gam(log(SturgeonConc+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SturgeonAbuLog10 <- gam(log(SturgeonConc+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
SturgeonAbuLog11 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
SturgeonAbuLog12 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

SturgeonAbuLog13 <- gam(log(SturgeonConc+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
SturgeonAbuLog14 <- gam(log(SturgeonConc+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))



AICctab(SturgeonAbuLog0,SturgeonAbuLog1,SturgeonAbuLog2,SturgeonAbuLog3,SturgeonAbuLog4,SturgeonAbuLog5,SturgeonAbuLog6,SturgeonAbuLog7,SturgeonAbuLog8,
        SturgeonAbuLog9,SturgeonAbuLog10,SturgeonAbuLog11,SturgeonAbuLog12,SturgeonAbuLog13,SturgeonAbuLog14,
        weights=T)

SturgeonAbuLog7$aic
summary(SturgeonAbuLog7)
summary(SturgeonAbuLog12)
cor(metadata$Temp,metadata$CTUSturgeon)


concurvity(SturgeonAbuLog12)
##########
#Top Sturgeon model
#########

summary(SturgeonAbuLog7)
gam.check(SturgeonAbuLog7)
plot(SturgeonAbuLog7)
appraise(SturgeonAbuLog7) #Model ok at best, not great




plot<-plot_smooth(SturgeonAbuLog7, view="CTUSturgeon",rm.ranef=T,sim.ci = T)
#plot_smooth(SturgeonAbuLog7, view="DischargeSampledByNight",rm.ranef=T,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SturgeonConc","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=SturgeonConc))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+
  geom_line(aes(CTUSturgeon,LowerCI),color="red")

ggplot(ShannonSubset,aes(CTUSturgeon,SturgeonConc))+geom_point()

SturgeonGAMPlot<-ggplot(ShannonSubset,aes(CTUSturgeon,SturgeonConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Larval~Sturgeon~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))#+ylim(NA, 1000)
SturgeonGAMPlot
# mean(ShannonSubset$Temp)
# PredictedDataFramelog<-data.frame((plot$fv$fit),plot$fv$temp_centered+19.207,(plot$fv$ul),(plot$fv$ll))
# colnames(PredictedDataFramelog)<-c("SturgeonConc","temp_centered","UpperCI","LowerCI")
# 
# #ggplot(ShannonSubset,aes(CTUSturgeon,log(SturgeonConc)))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5,aes(y=log(SturgeonConc)))+  
#   geom_line(data = PredictedDataFrame, aes(y = log(LowerCI)), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=log(UpperCI)),size=0.75,linetype="dashed")+
#   xlab("Temperature Centered")+ylab(expression(Log(Larval~Catostomidae~Per~100~m^3~Drift)))
# 


dev.off()
tiff("Figures/SturgeonGAM.tiff", width = 84, height = 84, units = 'mm', res = 1200)
SturgeonGAMPlot
dev.off()


dev.off()
tiff("Figures/SturgeonGAMAppraisal.tiff", width = 174, height = 174, units = 'mm', res = 1200)
appraise(SturgeonAbuLog7)
dev.off()

concurvity(SturgeonAbuLog7)

#####
#  Best Sturgeon GAM Temp Plot
#####

plot<-plot_smooth(SturgeonAbuLog7, view="Temp",rm.ranef=T,sim.ci = T)
#plot_smooth(SturgeonAbuLog7, view="DischargeSampledByNight",rm.ranef=T,sim.ci = T)
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$Temp,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("SturgeonConc","Temp","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=Temp,y=SturgeonConc))+geom_point()+geom_line(aes(Temp, UpperCI),color="red")+
  geom_line(aes(Temp,LowerCI),color="red")

ggplot(ShannonSubset,aes(Temp,SturgeonConc))+geom_point()

SturgeonGAMTempPlot<-ggplot(ShannonSubset,aes(Temp,SturgeonConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Water Temperature (\u00B0C)")+ylab(expression(Larval~Sturgeon~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))+ylim(NA, 15)
SturgeonGAMTempPlot


dev.off()
tiff("Figures/SturgeonGAMTempPlot.tiff", width = 174, height = 174, units = 'mm', res = 1200)
SturgeonGAMTempPlot
dev.off()




#######
# Invert GAM Models
#######


ShannonSubset<-ShannonRichness
ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, DischargeSampledByNight!= "NA")
ShannonSubset$Year<-as.factor(ShannonSubset$Year)
levels(ShannonSubset$Year)

ShannonSubset<-subset(ShannonSubset, AverageNetFlowByNight!= "NA")
ShannonSubset<-subset(ShannonSubset, QManual!= "NA")

ShannonSubset<-subset(ShannonSubset,ShannonSubset$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
ShannonSubset$temp_centered = ShannonSubset$Temp - mean(ShannonSubset$Temp)

ShannonSubset$Qcentered = ShannonSubset$QManual - mean(ShannonSubset$QManual)
ShannonSubset$AverageFlowCentered = ShannonSubset$AverageNetFlowByNight - mean(ShannonSubset$AverageNetFlowByNight)
ShannonSubset$DischargeCentered = ShannonSubset$DischargeSampledByNight - mean(ShannonSubset$DischargeSampledByNight)
ShannonSubset$DPFSCentered = ShannonSubset$DPFS - mean(ShannonSubset$DPFS)
ShannonSubset$DoYCentered = ShannonSubset$DayOfYear - mean(ShannonSubset$DayOfYear)
ShannonSubset$InvertConc<-ShannonSubset$DriftInvertConc # Changed to same name format as sturgeon/sucker so find/replace can be  used


head(ShannonSubset)
InvertAbuLog0 <- gam(log(InvertConc+0.001)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
InvertAbuLog1 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog2 <- gam(log(InvertConc+0.001)~ s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog3 <- gam(log(InvertConc+0.001)~ (MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog4 <- gam(log(InvertConc+0.001)~ s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
InvertAbuLog5 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog6 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog7 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


InvertAbuLog8 <- gam(log(InvertConc+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog9 <- gam(log(InvertConc+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

InvertAbuLog10 <- gam(log(InvertConc+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
InvertAbuLog11 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
InvertAbuLog12 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

InvertAbuLog13 <- gam(log(InvertConc+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
InvertAbuLog14 <- gam(log(InvertConc+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))



AICctab(InvertAbuLog0,InvertAbuLog1,InvertAbuLog2,InvertAbuLog3,InvertAbuLog4,InvertAbuLog5,InvertAbuLog6,InvertAbuLog7,InvertAbuLog8,
        InvertAbuLog9,InvertAbuLog10,InvertAbuLog11,InvertAbuLog12,InvertAbuLog13,InvertAbuLog14,
        weights=T)


summary(InvertAbuLog14)
appraise(InvertAbuLog14)
concurvity(InvertAbuLog14)
############
# Figure 2C Invert GAM plot
###########
summary(InvertAbuLog14)
concurvity(InvertAbuLog14)
plot<-plot_smooth(InvertAbuLog14, view="CTUSturgeon",rm.ranef=T,sim.ci = T)
InvertAbuLog14$coefficients
plot2<-plot_smooth(InvertAbuLog14, view="Temp",rm.ranef=T,sim.ci = F,plot_all = "MoonPhase",col=cbPalette,transform = exp,shade = F,ylab="Macroinvertebrate Concentration",xlab="Water Temperature (C)") 


levels(plot2$fv$Year)

plot_parametric(InvertAbuLog14, pred=list(MoonPhase=c("Full Moon", "Last Quarter","New Moon","Waning Crescent","Waning Gibbous","Waxing Crescent","Waxing Gibbous")))
coefs <- get_coefs(InvertAbuLog14)
CoefDF<-as.data.frame(coefs)
CoefDF$Estimate<-exp(CoefDF$Estimate)
CoefDF$`Std. Error`<-exp(CoefDF$`Std. Error`)
CoefDF$Group<-rownames(CoefDF)

CoefDF<-subset(CoefDF,Group!="(Intercept)")
ggplot(CoefDF,aes(x=Group,y=Estimate))+geom_point()+geom_errorbar(aes(ymin=Estimate-`Std. Error`,ymax=Estimate+`Std. Error`))


#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("DriftInvertConc","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


InvertGAMPlot<-ggplot(ShannonSubset,aes(CTUSturgeon,DriftInvertConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))#+ylim(NA, 1000)
InvertGAMPlot


dev.off()
tiff("Figures/InvertGAM.tiff", width = 174, height = 174, units = 'mm', res = 600)
InvertGAMPlot
dev.off()




appraise(InvertAbuLog14)

dev.off()
tiff("Figures/InvertGAMAppraisal.tiff", width = 174, height = 174, units = 'mm', res = 600)
appraise(InvertAbuLog14)
dev.off()


plot<-plot_smooth(InvertAbuLog14, view="Temp",rm.ranef=T,sim.ci = T)
plot_smooth(InvertAbuLog14, view="Temp",rm.ranef=F,sim.ci = T) 
#plot_sm
#plot$fv
FittedValues<-exp(plot$fv$fit)
head(FittedValues)

exp(7.5)
PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$Temp,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("DriftInvertConc","Temp","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=Temp,y=DriftInvertConc))+geom_point()+geom_line(aes(Temp, UpperCI),color="red")+
  geom_line(aes(Temp,LowerCI),color="red")

ggplot(ShannonSubset,aes(Temp,DriftInvertConc))+geom_point()

InvertGAMPlotTemp<-ggplot(ShannonSubset,aes(Temp,DriftInvertConc))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Water Temperature (\u00B0C)")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift))+theme(axis.title.x=element_text(size = 8))#+ylim(NA, 1000)
InvertGAMPlotTemp



dev.off()
tiff("Figures/InvertGAMTemp.tiff", width = 174, height = 174, units = 'mm', res = 600)
InvertGAMPlotTemp 
dev.off()










#####################
#Combine Abundance Model Figure
#####################
InvertGAMPlot
SturgeonGAMPlot
SuckerGAMPlot

dev.off()
tiff("Figures/GAMAbundanceModels.tiff", width = 174, height = 84, units = 'mm', res = 1200)

ggarrange(SturgeonGAMPlot,SuckerGAMPlot,InvertGAMPlot,
          labels = c("a", "b","c"),
          ncol = 3, nrow = 1)

dev.off()




###########
#Total InvertN By Moon Phase Normalized 
############


kruskal.test(DriftInvertConc~MoonPhase, data=ShannonRichness)

compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Means=compare_means(DriftInvertConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
Letters<-multcompLetters(difference)
Letters


DriftPlotNAsRemoved<-subset(ShannonRichness, DriftInvertConc!="NA")
DriftPlotNAsRemoved
Trtdata <- ddply(DriftPlotNAsRemoved, c("MoonPhase"), summarise,
                 N    = length(DriftInvertConc),
                 meanSturgeon = mean(DriftInvertConc),
                 sd   = sd(DriftInvertConc),
                 se   = sd / sqrt(N),na.rm =T
                 
)
Trtdata

#manually renamed due to plotting where new moon doesn't start with a
Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))
Letters$Letters# b->a, a->b, c=c

Letters2<-c("bc","bc","abc","a","ab","ab","ab","c")

Trtdata

TotalInvertAbuMoonPhase<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(),stat="identity")+xlab("Moon Phase")+ylab(expression(Macroinvertebrates~Per~100~m^3~Drift~(SE)))+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 10))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+1,label=Letters2))+theme(legend.position = "none")+geom_text(aes(x=5,y=25,label= "KW, Chi-sq = 39.98, P < 0.001"),size=4)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
TotalInvertAbuMoonPhase
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TotalInvertebrateAbuByPhaseNormalized.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TotalInvertAbuMoonPhase
dev.off()


###############
#Hourly Results
###############


Hourly<-read.csv("DataClean/SturgeonDriftDataHourly12.21.19.csv",header=T)
head(Hourly)
Hourly$Time = factor(Hourly$Time, levels = c("Ten","Eleven","Twelve","One","Two"))
Hourly$DriftInvertConc<-((Hourly$InvertSample100*100)/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
Hourly$DriftSuckerConc<-(Hourly$Suckers100*100/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
Hourly$DriftSturgeonConc<-(Hourly$LiveSturgeon*100/(60*60*Hourly$AreaSampled*Hourly$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)




StatSubset<-subset(Hourly, DriftInvertConc!= "NA") # Remove NAs for hours where no discharge was recorded

hist(log(Hourly$InvertSample100))
plot(DriftInvertConc~Time,data= StatSubset)


kruskal.test((DriftInvertConc)~Time,data=StatSubset)

compare_means(DriftInvertConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")


Means=compare_means(DriftInvertConc ~ Time, data = StatSubset, p.adjust.method = "fdr",method="wilcox.test")
Means
Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersInvert<-multcompLetters(difference)
LettersInvert<-LettersInvert$Letters
LettersInvert

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
  annotate("text", label = "Kruskal-Wallis,\n Chi-sq = 188.4, \nP < 0.001", size = 3.5, x = 1.5, y = 15)+ theme(axis.title.y = element_text(size = 10))
HourlyInvertAbu

dev.off()
tiff("Figures/Hourly_InvertAbundance.tiff", width = 84, height = 84, units = 'mm', res = 1200)
HourlyInvertAbu
dev.off()


#Sucker by collection time
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
  geom_errorbar(aes(ymin=meanSucker-se,ymax=meanSucker+se))+ylab(expression(Sucker~Larvae~Per~100~m^3~Drift~(SE)))+
  scale_x_discrete(labels=c("22:00","23:00","0:00","1:00","2:00"))+geom_text(aes(x=Time, y=meanSucker+se+5,label=LettersSucker$Letters))+
  annotate("text", label = "Kruskal-Wallis,\n Chi-sq = 126.1,\n P < 0.001", size = 3.5, x = 4.5, y = 90)+ theme(axis.title.y = element_text(size = 10))
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


#######
# SturgeonAndSuckerByMoonPhase
#######
kruskal.test(SuckerConc~MoonPhase, data=ShannonRichness)
kruskal.test(SturgeonConc~MoonPhase, data=ShannonRichness)

compare_means(SturgeonConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")

Means=compare_means(SturgeonConc ~ MoonPhase, data = ShannonRichness, p.adjust.method = "fdr",method="wilcox.test")
subset(Means, p.adj<0.1)

DriftPlotNAsRemoved<-subset(ShannonRichness, SturgeonConc!="NA")
DriftPlotNAsRemoved
Trtdata <- ddply(DriftPlotNAsRemoved, c("MoonPhase"), summarise,
                 N    = length(SturgeonConc),
                 meanSturgeon = mean(SturgeonConc),
                 sd   = sd(SturgeonConc),
                 se   = sd / sqrt(N),na.rm =T
                 
)
Trtdata

Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

SturgeonMoonPhase<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(),stat="identity")+xlab("Moon Phase")+ylab(expression(Larval~Sturgeon~Per~100~m^3~Drift~(SE)))+
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 10))+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+geom_text(aes(x=5,y=8,label= "KW, Chi-sq = 15.17, P = 0.03"),size=4)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
SturgeonMoonPhase


dev.off()
tiff("Figures/SturgeonMoonPhase.tiff", width = 84, height = 84, units = 'mm', res = 600)
SturgeonMoonPhase
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
DriftTotalInvertBiomass<-ggplot(Trtdata, aes(x=MoonPhase,y=meanSturgeon))+geom_bar(aes(),stat="identity")+xlab("Moon Phase")+ylab(expression(Macroinvertebrate~Biomass~(g)~Per~100~m^3~Drift~(SE)))+#Invertebrate Biomass (g) per 100 m3 drift (SEM)
  geom_errorbar(aes(ymin=meanSturgeon-se,ymax=meanSturgeon+se))+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.y = element_text(size = 8))+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  geom_text(aes(x=MoonPhase, y=meanSturgeon+se+.08,label=vector))+geom_text(aes(x=5,y=1.5,label= "KW, chi-squared = 49.3, P < 0.001"),size=3)+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))
DriftTotalInvertBiomass
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/BiomassByPhase100m3.tiff", width = 3.3, height = 3.3, units = 'in', res = 800)
DriftTotalInvertBiomass
dev.off()



###############
#Join together fig 3
##############

TotalInvertAbuMoonPhase
HourlyInvertAbu
HourlySuckerAbu
HourlySturgeonAbu
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/Fig3Mar2022.tiff", width = 174, height = 174, units = 'mm', res = 1200)

ggarrange(DriftTotalInvertBiomass,HourlyInvertAbu,HourlySturgeonAbu,HourlySuckerAbu,
          labels = c("a", "b","c","d"),
          ncol = 2, nrow = 2)

dev.off()



######
# Macroinvertebrate biomass GAM
######

head(ShannonSubset)
BiomassAbuLog0 <- gam(log(InvertBiomass100+0.001)~ 1+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
BiomassAbuLog1 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog2 <- gam(log(InvertBiomass100+0.001)~ s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog3 <- gam(log(InvertBiomass100+0.001)~ (MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog4 <- gam(log(InvertBiomass100+0.001)~ s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
BiomassAbuLog5 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog6 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog7 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


BiomassAbuLog8 <- gam(log(InvertBiomass100+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog9 <- gam(log(InvertBiomass100+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

BiomassAbuLog10 <- gam(log(InvertBiomass100+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
BiomassAbuLog11 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
BiomassAbuLog12 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

BiomassAbuLog13 <- gam(log(InvertBiomass100+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
BiomassAbuLog14 <- gam(log(InvertBiomass100+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = ShannonSubset, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


AICctab(BiomassAbuLog0,BiomassAbuLog1,BiomassAbuLog2,BiomassAbuLog3,BiomassAbuLog4,BiomassAbuLog5,BiomassAbuLog6,BiomassAbuLog7,BiomassAbuLog8,
        BiomassAbuLog9,BiomassAbuLog10,BiomassAbuLog11,BiomassAbuLog12,BiomassAbuLog13,BiomassAbuLog14,
        weights=T)

summary(BiomassAbuLog14)



############
#Total Abundance Family by Moon phase
############
otufull=read.table("DataClean/SturgeonDriftInvertAbundanceMatrix8.15.19.txt",header=TRUE)
taxmatrixfull=as.matrix(read.table("DataClean/SturgeonDriftInvertTaxNames8.15.19.txt"))
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





DriftDataCombined<-read.csv("DataClean/AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)
head(DriftDataCombined)
colnames(DriftDataCombined)


DriftDataCombined<-subset(DriftDataCombined,DriftDataCombined$DischargeSampledByNight < 2) #Remove discharge sampled outliers 

#Note: CTUSturgeon in this file still using base 0!!!



physeqSubset<-subset_taxa(physeq,Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")
physeqSubset

#df <- psmelt(physeqSubset)
df<-psmelt(physeq)
df$Abundance<-df$Abundance*20 #To get total number based on 5% subsample ID
head(df)
df$AbuPer100<-((df$Abundance*100)/(60*4*60*df$AreaSampled.m2.*df$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
df<-subset(df,AbuPer100!="NA")
df<-subset(df,DischargeSampledByNight<2) # Remove discharge sampled outliers

head(df)
levels(df$Family)
FamiliesPer100m3<-compare_means(AbuPer100 ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")


subset(FamiliesPer100m3,p.adj < 0.1)
subset(FamiliesPer100m3,Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")


HepAbuPer100<- subset(df,Family=="Heptageniidae")
IsoAbuPer100<- subset(df,Family=="Isonychiidae")
Means<-compare_means(AbuPer100 ~ MoonPhase, data = HepAbuPer100, group.by = "Family", p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHep<-multcompLetters(difference)
LettersHep

Means<-compare_means(AbuPer100 ~ MoonPhase, data = IsoAbuPer100, group.by = "Family", p.adjust.method = "fdr",method="wilcox.test")

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersIso<-multcompLetters(difference)
LettersIso


Trtdata <- ddply(df, c("MoonPhase","Family"), summarise,
                 N    = length(AbuPer100),
                 mean = mean(AbuPer100),
                 sd   = sd(AbuPer100),
                 se   = sd / sqrt(N)
)

Trtdata
Trtdata<-subset(Trtdata,Family=="Crayfish"|Family=="Ephemerillidae"|Family=="Heptageniidae"|Family=="Hydropsychidae"|Family=="Isonychiidae"|Family=="Leptoceridae")

Trtdata$MoonPhase = factor(Trtdata$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))


LettersHep
LettersIso

TrtdataSorted<-Trtdata[order(Trtdata$Family),]
head(TrtdataSorted)
vector<-c("","","","","","","","",
          "","","","","","","","",
          "ab"," abc","cd","d","acd ","b","abc","abc",
          "","","","","","","","",
          "a"," ab","cd","c","bcd ","bd","abd","ab",
          "","","","","","","","")


length(vector)
length(TrtdataSorted$N)
KruskalLabel<- c("Kruskal-Wallis,\n P-adj = 0.31","KW, P-adj = 0.92", "KW, P-adj < 0.001","KW, P-adj = 0.054","     P-adj < 0.001","KW, P-adj = 0.51")

dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.31","KW, P-adj = 0.92", "  KW, P-adj < 0.001","  KW, P-adj = 0.054","     P-adj < 0.001","KW, P-adj = 0.51"),
  Family   = c("Crayfish","Ephemerillidae","Heptageniidae","Hydropsychidae","Isonychiidae","Leptoceridae")
)
#TrtdataSorted

TopFamilyAbuPer100=ggplot(TrtdataSorted, aes(x=MoonPhase,y=mean))+geom_bar(aes(),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Macroinvertebrates~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+scale_fill_manual(values=cbPalette)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+geom_text(aes(x=MoonPhase,y= mean+se+1),label=vector,size=2)
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

BiomassAverages<-read.table("DataClean/BiomassFamilyAveragesR.txt",header=T)
head(BiomassAverages)
DriftDataCombined<-read.csv("DataClean/AllDriftDataCombined2011-2018FamilyRichness.csv",header=T)
otufull=read.table("DataClean/SturgeonDriftInvertAbundanceMatrix8.15.19.txt",header=TRUE)
taxmatrixfull=as.matrix(read.table("DataClean/SturgeonDriftInvertTaxNames8.15.19.txt"))

otubiomass<-otufull

row.names(otubiomass)<-taxmatrixfull

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
#Trtdata
#options(scipen = 999)
#Trtdata



AllFamilyBiomass=ggplot(Trtdata, aes(x=MoonPhase,y=mean))+geom_bar(aes(),colour="black", stat="identity")+xlab("Moon Phase")+
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

df <- psmelt(physeqBiomass)
df$Abundance<-df$Abundance*20# To account for 5% subsample used for ID
df$BiomassPer100<-((df$Abundance*100)/(60*4*60*df$AreaSampled.m2.*df$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
df$PercentageTotalBiomass<-df$BiomassPer100/df$DriftBiomassConc*100 #Percentage of total biomass for each family

dfSubset<-subset(df,BiomassPer100!="NA")
dfSubset<-subset(dfSubset,DischargeSampledByNight< 2)
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
sum(Top6BiomassFam$meanPercent) #85% of total biomass comes from top six families
Top6BiomassFam
BiomassPer100m3<-compare_means(BiomassPer100 ~ MoonPhase, data = df, group.by = "Family", p.adjust.method = "fdr",method="kruskal.test")
#BiomassPer100m3
subset(BiomassPer100m3,Family %in% Top6BiomassFam$Family)


Trtdata <- ddply(dfSubset, c("Family","MoonPhase"), summarise,
                 N    = length(BiomassPer100),
                 mean = mean(BiomassPer100),
                 sd   = sd(BiomassPer100),
                 se   = sd / sqrt(N))
TrtdataSubset<-subset(Trtdata,Family %in% Top6BiomassFam$Family)
TrtdataSubset


Hep<-subset(dfSubset,Family=="Heptageniidae")
Iso<-subset(dfSubset,Family=="Isonychiidae")




#Heptageniidae

kruskal.test(BiomassPer100~MoonPhase,data=Hep)
Means<-compare_means(Heptageniidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersHep<-multcompLetters(difference)
LettersHep
# 
# #Isonychiidae
kruskal.test(BiomassPer100~MoonPhase,data=Iso)
Means<-compare_means(Isonychiidae ~ MoonPhase, data = DriftDataCombined, p.adjust.method = "fdr")
Means

Hyphenated<-as.character(paste0(Means$group1,"-",Means$group2))
difference<-Means$p.adj
names(difference)<-Hyphenated
LettersIso<-multcompLetters(difference)
LettersIso

LettersHep
LettersIso
vector<-c("","","","","","","","",
          "","","","","","","","",
          "","","","","","","","",
          "abc"," abc"," bcd","d","cd","a","abc","ab",
          "a"," abc","de","e","cde ","bcd ","abcd","ab",
          "","","","","","","","")

# Isonychiidae d=a,a=b,b=c,c=d,e=e


dat_text <- data.frame(
  label = c("Kruskal-Wallis,\n P-adj = 0.2","KW, P-adj = 0.89","KW, P-adj = 0.55", "    KW, P-adj < 0.001",  "      P-adj < 0.001","KW, P-adj = 0.89"),
  Family   = c("Crayfish","Ephemerillidae","Gomphidae","Heptageniidae","Isonychiidae","Lepidostomatidae")
)
TopFamilyBiomass=ggplot(TrtdataSubset, aes(x=MoonPhase,y=mean))+geom_bar(aes(),colour="black", stat="identity")+xlab("Moon Phase")+
  ylab(expression(Biomass~(g)~Per~100~m^3~Drift~(SE))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  facet_wrap(Family~.)+theme(legend.position = "none")+
  scale_x_discrete(labels=c("New","WXC","FQ","WXG","Full","WAG","LQ","WNC"))+ geom_text(data=dat_text,size=2.4,mapping = aes(x = 4, y = 0.92, label = label))+scale_fill_manual(values=cbPalette)+
  geom_text(aes(x=MoonPhase,y= mean+se+.07),label=vector,size=2)
TopFamilyBiomass


theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/TopFamiliesFamiliesBiomass.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TopFamilyBiomass
dev.off()



#########
#Crayfish model
#########
library(itsadug)
library(gratia)

FamilyGAMData<-psmelt(physeq)
head(FamilyGAMData)

FamilyGAMData<-subset(FamilyGAMData, Temp!= "NA")
FamilyGAMData<-subset(FamilyGAMData, AverageNetFlowByNight!= "NA")
FamilyGAMData<-subset(FamilyGAMData, QManual!= "NA")

FamilyGAMData<-subset(FamilyGAMData,FamilyGAMData$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
FamilyGAMData$temp_centered = FamilyGAMData$Temp - mean(FamilyGAMData$Temp)
FamilyGAMData$Qcentered = FamilyGAMData$QManual - mean(FamilyGAMData$QManual)
FamilyGAMData$AverageFlowCentered = FamilyGAMData$AverageNetFlowByNight - mean(FamilyGAMData$AverageNetFlowByNight)
FamilyGAMData$DischargeCentered = FamilyGAMData$DischargeSampledByNight - mean(FamilyGAMData$DischargeSampledByNight)
FamilyGAMData$DPFSCentered = FamilyGAMData$DPFS - mean(FamilyGAMData$DPFS)
FamilyGAMData$DOYCentered = FamilyGAMData$DayOfYear - mean(FamilyGAMData$DayOfYear)

head(FamilyGAMData)


CrayfishAbu<-subset(FamilyGAMData,Family=="Crayfish")


head(CrayfishAbu)
CrayfishAbu$Abundance<-CrayfishAbu$Abundance*20 #Account for 5% sampling
CrayfishAbu$CrayfishPer100<-((CrayfishAbu$Abundance*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
CrayfishAbu$SturgeonPer100<-((CrayfishAbu$Nsturgeon*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
CrayfishAbu$SuckerPer100<-((CrayfishAbu$Nsuckers100*100)/(60*4*60*CrayfishAbu$AreaSampled.m2.*CrayfishAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
kruskal.test(CrayfishPer100~MoonPhase, data=CrayfishAbu)


head(CrayfishAbu)
CrayfishAbuLog0 <- gam(log(CrayfishPer100+0.001)~ 1+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
CrayfishAbuLog1 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog2 <- gam(log(CrayfishPer100+0.001)~ s(QManual)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog3 <- gam(log(CrayfishPer100+0.001)~ (MoonPhase)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog4 <- gam(log(CrayfishPer100+0.001)~ s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
CrayfishAbuLog5 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog6 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog7 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


CrayfishAbuLog8 <- gam(log(CrayfishPer100+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog9 <- gam(log(CrayfishPer100+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

CrayfishAbuLog10 <- gam(log(CrayfishPer100+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
CrayfishAbuLog11 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
CrayfishAbuLog12 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

CrayfishAbuLog13 <- gam(log(CrayfishPer100+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
CrayfishAbuLog14 <- gam(log(CrayfishPer100+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = CrayfishAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


AICctab(CrayfishAbuLog0,CrayfishAbuLog1,CrayfishAbuLog2,CrayfishAbuLog3,CrayfishAbuLog4,CrayfishAbuLog5,CrayfishAbuLog6,CrayfishAbuLog7,CrayfishAbuLog8,
CrayfishAbuLog9,CrayfishAbuLog10,CrayfishAbuLog11,CrayfishAbuLog12,CrayfishAbuLog13,CrayfishAbuLog14,
        weights=T)

summary(CrayfishAbuLog1)

summary(CrayfishAbuLog6)

appraise(CrayfishAbuLog1)

##########
#Crayfish model
#########
#Note: Both Gamma and log based GAM models performed poorly, but models show a very strong effect of CTU on crayfish density
plot<-plot_smooth(CrayfishAbuLog1, view="CTUSturgeon",rm.ranef=T,sim.ci = T)
#plot_sm
#plot$fv

FittedValues<-exp(plot$fv$fit)
head(FittedValues)
min(FittedValues)

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

FamilyGAMData<-psmelt(physeq)
head(FamilyGAMData)

FamilyGAMData<-subset(FamilyGAMData, Temp!= "NA")
FamilyGAMData<-subset(FamilyGAMData, AverageNetFlowByNight!= "NA")
FamilyGAMData<-subset(FamilyGAMData, QManual!= "NA")

FamilyGAMData<-subset(FamilyGAMData,FamilyGAMData$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
FamilyGAMData$temp_centered = FamilyGAMData$Temp - mean(FamilyGAMData$Temp)
FamilyGAMData$Qcentered = FamilyGAMData$QManual - mean(FamilyGAMData$QManual)
FamilyGAMData$AverageFlowCentered = FamilyGAMData$AverageNetFlowByNight - mean(FamilyGAMData$AverageNetFlowByNight)
FamilyGAMData$DischargeCentered = FamilyGAMData$DischargeSampledByNight - mean(FamilyGAMData$DischargeSampledByNight)
FamilyGAMData$DPFSCentered = FamilyGAMData$DPFS - mean(FamilyGAMData$DPFS)
FamilyGAMData$DOYCentered = FamilyGAMData$DayOfYear - mean(FamilyGAMData$DayOfYear)

head(FamilyGAMData)


IsonychiidaeAbu<-subset(FamilyGAMData,Family=="Isonychiidae")

head(IsonychiidaeAbu)
IsonychiidaeAbu$Abundance<-IsonychiidaeAbu$Abundance*20 #Account for 5% sampling
IsonychiidaeAbu$IsonychiidaePer100<-((IsonychiidaeAbu$Abundance*100)/(60*4*60*IsonychiidaeAbu$AreaSampled.m2.*IsonychiidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
IsonychiidaeAbu$SturgeonPer100<-((IsonychiidaeAbu$Nsturgeon*100)/(60*4*60*IsonychiidaeAbu$AreaSampled.m2.*IsonychiidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
IsonychiidaeAbu$SuckerPer100<-((IsonychiidaeAbu$Nsuckers100*100)/(60*4*60*IsonychiidaeAbu$AreaSampled.m2.*IsonychiidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)

kruskal.test(IsonychiidaePer100~MoonPhase, data=IsonychiidaeAbu)


head(IsonychiidaeAbu)
IsonychiidaeAbuLog0 <- gam(log(IsonychiidaePer100+0.001)~ 1+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
IsonychiidaeAbuLog1 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog2 <- gam(log(IsonychiidaePer100+0.001)~ s(QManual)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog3 <- gam(log(IsonychiidaePer100+0.001)~ (MoonPhase)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog4 <- gam(log(IsonychiidaePer100+0.001)~ s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
IsonychiidaeAbuLog5 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog6 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog7 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


IsonychiidaeAbuLog8 <- gam(log(IsonychiidaePer100+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog9 <- gam(log(IsonychiidaePer100+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

IsonychiidaeAbuLog10 <- gam(log(IsonychiidaePer100+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
IsonychiidaeAbuLog11 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
IsonychiidaeAbuLog12 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

IsonychiidaeAbuLog13 <- gam(log(IsonychiidaePer100+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
IsonychiidaeAbuLog14 <- gam(log(IsonychiidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = IsonychiidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


AICctab(IsonychiidaeAbuLog0,IsonychiidaeAbuLog1,IsonychiidaeAbuLog2,IsonychiidaeAbuLog3,IsonychiidaeAbuLog4,IsonychiidaeAbuLog5,IsonychiidaeAbuLog6,IsonychiidaeAbuLog7,IsonychiidaeAbuLog8,
        IsonychiidaeAbuLog9,IsonychiidaeAbuLog10,IsonychiidaeAbuLog11,IsonychiidaeAbuLog12,IsonychiidaeAbuLog13,IsonychiidaeAbuLog14,
        weights=T)


summary(IsonychiidaeAbuLog5)
appraise(IsonychiidaeAbuLog5)



plot<-plot_smooth(HeptageniidaeAbuLog5, view="CTUSturgeon",rm.ranef=T,sim.ci = T,cond = list(MoonPhase = "First Quarter"))

plot_smooth(IsonychiidaeAbuLog5, view="CTUSturgeon",rm.ranef=T,sim.ci = T,transform=exp)
FittedValues<-exp(plot$fv$fit)

PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("IsonychiidaePer100","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=IsonychiidaePer100))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")


IsonychiidaeGAMPlot<-ggplot(IsonychiidaeAbu,aes(CTUSturgeon,IsonychiidaePer100))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Isonychiidae~Per~100~m^3~Drift))+ylim(NA,50)
IsonychiidaeGAMPlot





##############
#Heptageniidae model
############
FamilyGAMData<-psmelt(physeq)
head(FamilyGAMData)

FamilyGAMData<-subset(FamilyGAMData, Temp!= "NA")
FamilyGAMData<-subset(FamilyGAMData, AverageNetFlowByNight!= "NA")
FamilyGAMData<-subset(FamilyGAMData, QManual!= "NA")

FamilyGAMData<-subset(FamilyGAMData,FamilyGAMData$DischargeSampledByNight < 2) #Remove discharge sampled outliers 
FamilyGAMData$temp_centered = FamilyGAMData$Temp - mean(FamilyGAMData$Temp)
FamilyGAMData$Qcentered = FamilyGAMData$QManual - mean(FamilyGAMData$QManual)
FamilyGAMData$AverageFlowCentered = FamilyGAMData$AverageNetFlowByNight - mean(FamilyGAMData$AverageNetFlowByNight)
FamilyGAMData$DischargeCentered = FamilyGAMData$DischargeSampledByNight - mean(FamilyGAMData$DischargeSampledByNight)
FamilyGAMData$DPFSCentered = FamilyGAMData$DPFS - mean(FamilyGAMData$DPFS)
FamilyGAMData$DOYCentered = FamilyGAMData$DayOfYear - mean(FamilyGAMData$DayOfYear)

head(FamilyGAMData)


HeptageniidaeAbu<-subset(FamilyGAMData,Family=="Heptageniidae")


head(HeptageniidaeAbu)
HeptageniidaeAbu$Abundance<-HeptageniidaeAbu$Abundance*20 #Account for 5% sampling
HeptageniidaeAbu$HeptageniidaePer100<-((HeptageniidaeAbu$Abundance*100)/(60*4*60*HeptageniidaeAbu$AreaSampled.m2.*HeptageniidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
HeptageniidaeAbu$SturgeonPer100<-((HeptageniidaeAbu$Nsturgeon*100)/(60*4*60*HeptageniidaeAbu$AreaSampled.m2.*HeptageniidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)
HeptageniidaeAbu$SuckerPer100<-((HeptageniidaeAbu$Nsuckers100*100)/(60*4*60*HeptageniidaeAbu$AreaSampled.m2.*HeptageniidaeAbu$AverageNetFlowByNight)) #Calculate inverts/ 100 m3 water (N*100(final vol )/time in sec*flow*area sampled)



head(HeptageniidaeAbu)
HeptageniidaeAbuLog0 <- gam(log(HeptageniidaePer100+0.001)~ 1+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
# Individual terms
# CTUSturgeon, Q,MoonPhase, Temp
HeptageniidaeAbuLog1 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog2 <- gam(log(HeptageniidaePer100+0.001)~ s(QManual)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog3 <- gam(log(HeptageniidaePer100+0.001)~ (MoonPhase)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog4 <- gam(log(HeptageniidaePer100+0.001)~ s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


# Two term models
HeptageniidaeAbuLog5 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+(MoonPhase)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog6 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog7 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


HeptageniidaeAbuLog8 <- gam(log(HeptageniidaePer100+0.001)~ s(QManual)+(MoonPhase)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog9 <- gam(log(HeptageniidaePer100+0.001)~ s(QManual)+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

HeptageniidaeAbuLog10 <- gam(log(HeptageniidaePer100+0.001)~ MoonPhase+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 3 term models
# CTUSturgeon, Q,MoonPhase,Temp
HeptageniidaeAbuLog11 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+(MoonPhase)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))
HeptageniidaeAbuLog12 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

HeptageniidaeAbuLog13 <- gam(log(HeptageniidaePer100+0.001)~ s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))

# 4 term models
# CTUSturgeon, Q,MoonPhase,Temp
HeptageniidaeAbuLog14 <- gam(log(HeptageniidaePer100+0.001)~ s(CTUSturgeon)+s(QManual)+MoonPhase+s(Temp)+s(Year,bs="re"), data = HeptageniidaeAbu, method = "REML",correlation=corAR1(form = ~DayOfYear|Year))


AICctab(HeptageniidaeAbuLog0,HeptageniidaeAbuLog1,HeptageniidaeAbuLog2,HeptageniidaeAbuLog3,HeptageniidaeAbuLog4,HeptageniidaeAbuLog5,HeptageniidaeAbuLog6,HeptageniidaeAbuLog7,HeptageniidaeAbuLog8,
        HeptageniidaeAbuLog9,HeptageniidaeAbuLog10,HeptageniidaeAbuLog11,HeptageniidaeAbuLog12,HeptageniidaeAbuLog13,HeptageniidaeAbuLog14,
        weights=T)


summary(HeptageniidaeAbuLog5)
appraise(HeptageniidaeAbuLog5)
plot_smooth(HeptageniidaeAbuLog5, view="CTUSturgeon",rm.ranef=T,sim.ci = T,cond = list(MoonPhase = "First Quarter"))

plot<-plot_smooth(HeptageniidaeAbuLog5, view="CTUSturgeon",rm.ranef=T,sim.ci = T,cond = list(MoonPhase = "First Quarter"))


FittedValues<-exp(plot$fv$fit)

PredictedDataFrame<-data.frame(exp(plot$fv$fit),plot$fv$CTUSturgeon,exp(plot$fv$ul),exp(plot$fv$ll))

colnames(PredictedDataFrame)<-c("HeptageniidaePer100","CTUSturgeon","UpperCI","LowerCI")
head(PredictedDataFrame)


ggplot(PredictedDataFrame, aes(x=CTUSturgeon,y=HeptageniidaePer100))+geom_point()+geom_line(aes(CTUSturgeon, UpperCI),color="red")+geom_line(aes(CTUSturgeon,LowerCI),color="red")
ggplot(HeptageniidaeAbu, aes(x=CTUSturgeon,y=HeptageniidaePer100))+geom_point(color="black")+ylab(expression(Heptageniidae~Per~100~m^3~Drift))+xlab("Day of Year")#+facet_wrap(~Year)


HeptageniidaeGAMPlot<-ggplot(HeptageniidaeAbu,aes(CTUSturgeon,HeptageniidaePer100))+geom_point(color="darkgrey")+geom_line(data=PredictedDataFrame,size=1.5)+  
  geom_line(data = PredictedDataFrame, aes(y = LowerCI), size = .75,linetype="dashed")+geom_line(data = PredictedDataFrame,aes(y=UpperCI),size=0.75,linetype="dashed")+
  xlab("Cumulative Temperature Units")+ylab(expression(Heptageniidae~Per~100~m^3~Drift))
HeptageniidaeGAMPlot









#######
# Combine Hep and Iso plot
######

dev.off()
tiff("Figures/HepIsoPlot.tiff", width = 174, height = 84, units = 'mm', res = 1200)

ggarrange(IsonychiidaeGAMPlot,HeptageniidaeGAMPlot,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

dev.off()



################
#Percentage Nightly Biomass (Fig 6 in latest revision)
################
#Percentage Nightly Biomass

ShannonRichness$SturgeonBiomassNightly<-ShannonRichness$Nsturgeon*0.005 #Average sturgeon Nightly biomass (early and late spawning combined)
#Sucker Biomass already calculated individually for early and late season sucker species 
#dry weight 1.53+/-0.1 mg white sucker
#dry weight 2.56 +/- 0.1 silver redhorse
ShannonRichness$SuckerBiomassNightly<-ShannonRichness$TotalSuckerBiomass.g
ShannonRichness$CombinedNightlyBiomass<-ShannonRichness$SturgeonBiomassNightly+ShannonRichness$SuckerBiomassNightly+ShannonRichness$InvertBiomass100

#Nightly combined biomass
mean(ShannonRichness$CombinedNightlyBiomass) #112.34 g
min(ShannonRichness$CombinedNightlyBiomass) #1.36 g
max(ShannonRichness$CombinedNightlyBiomass) #845.9

#Sucker Biomass
mean(ShannonRichness$SuckerBiomassNightly) #13.67 g
min(ShannonRichness$SuckerBiomassNightly) #0 g
max(ShannonRichness$SuckerBiomassNightly) #312.1 g

#Sturgeon Biomass
mean(ShannonRichness$SturgeonBiomassNightly) #2.27 g
min(ShannonRichness$SturgeonBiomassNightly) #0 g
max(ShannonRichness$SturgeonBiomassNightly) #84.5 g


#Invert Biomass
mean(ShannonRichness$InvertBiomass100) #96.4 g
min(ShannonRichness$InvertBiomass100) #0.59 g
max(ShannonRichness$InvertBiomass100) #815.5 g

#a 23% decrease in A. fulvescens predation for each additional 52 g of aquatic invertebrates present in the drift 
PercentDecreasePredation<-23/52 #Percent decrease in predation for each g of inverts
mean(ShannonRichness$InvertBiomass100)*PercentDecreasePredation #42.6
min(ShannonRichness$InvertBiomass100)*PercentDecreasePredation #0.26
max(ShannonRichness$InvertBiomass100)*PercentDecreasePredation #360.7
42.6-(77.12*PercentDecreasePredation)


#Mean Invert biomass during a new moon
NewMoon<-subset(ShannonRichness,MoonPhase=="New Moon")
#head(NewMoon)
mean(NewMoon$InvertBiomass100) # 171.3g
min(NewMoon$InvertBiomass100) #48.59 g
max(NewMoon$InvertBiomass100) #327.9 g

#Mean Invert biomass during a full moon
FullMoon<-subset(ShannonRichness,MoonPhase=="Full Moon")
#head(fullMoon)
mean(FullMoon$InvertBiomass100) # 43.198g
min(FullMoon$InvertBiomass100) #11.91 g
max(FullMoon$InvertBiomass100) #99.86 g


#Difference between new and full moon % probability predation
#Mean Invert biomass new moon 171.3 g
#Mean Invert biomass full moon 43.2 g
#23% decrease for 52 g Inverts
(171.3-43.2)*PercentDecreasePredation #Diff in g Inverts* %Reduction in predation per g =56.66%


#28% reduction in A. fulvescens predation for each 10% of drift biomass made up by larval catostomids
#24.3 g of A. fulvescens biomass present in the drift there was an estimated 82% increase incidences of larval A. fulvescens predation

#Transform nightly biomass to biomass per 100 m3
ShannonRichness$DriftBiomassConc<-((ShannonRichness$InvertBiomass100*100)/(60*4*60*ShannonRichness$AreaSampled.m2.*ShannonRichness$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water
ShannonRichness$DriftSuckerBiomassConc<-((ShannonRichness$SuckerBiomassNightly*100)/(60*4*60*ShannonRichness$AreaSampled.m2.*ShannonRichness$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water
ShannonRichness$DriftSturgeonBiomassConc<-((ShannonRichness$SturgeonBiomassNightly*100)/(60*4*60*ShannonRichness$AreaSampled.m2.*ShannonRichness$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water
ShannonRichness$DriftCombinedBiomassConc<-((ShannonRichness$CombinedNightlyBiomass*100)/(60*4*60*ShannonRichness$AreaSampled.m2.*ShannonRichness$AverageNetFlowByNight)) #Calculate biomass/ 100 m3 water


kruskal.test(SuckerConc~MoonPhase,data=ShannonRichness)


#Nightly combined biomass per 100 m3
ShannonSubsetSampled<-subset(ShannonRichness,DriftCombinedBiomassConc!="NA")

mean(ShannonSubsetSampled$DriftCombinedBiomassConc) #0.87 g
min(ShannonSubsetSampled$DriftCombinedBiomassConc) #0.0111 g
max(ShannonSubsetSampled$DriftCombinedBiomassConc) #6.915

#Nightly  invert biomass per 100 m3
ShannonSubsetSampled<-subset(ShannonRichness,DriftBiomassConc!="NA")

mean(ShannonSubsetSampled$DriftBiomassConc) #0.746 g
min(ShannonSubsetSampled$DriftBiomassConc) #0.005 g
max(ShannonSubsetSampled$DriftBiomassConc) #6.78

#Nightly sucker biomass per 100 m3
ShannonSubsetSampled<-subset(ShannonRichness,DriftSuckerBiomassConc!="NA")

mean(ShannonSubsetSampled$DriftSuckerBiomassConc) #.11 g
min(ShannonSubsetSampled$DriftSuckerBiomassConc) #0 g
max(ShannonSubsetSampled$DriftSuckerBiomassConc) #2.41
length(ShannonSubsetSampled$DriftSuckerBiomassConc) #229

sd(ShannonSubsetSampled$DriftSuckerBiomassConc)/sqrt(229)
head(ShannonSubsetSampled)

#Nightly Sturgeon biomass per 100 m3
ShannonSubsetSampled<-subset(ShannonRichness,DriftSturgeonBiomassConc!="NA")

mean(ShannonSubsetSampled$DriftSturgeonBiomassConc) #0.0167
min(ShannonSubsetSampled$DriftSturgeonBiomassConc) #0 g
max(ShannonSubsetSampled$DriftSturgeonBiomassConc) #0.453
length(ShannonSubsetSampled$DriftSturgeonBiomassConc) #229

sd(ShannonSubsetSampled$DriftSturgeonBiomassConc)/sqrt(229) #0.00318
head(ShannonSubsetSampled)




ShannonRichness$PercentNightlyBiomassSturgeon<- ShannonRichness$SturgeonBiomassNightly/ShannonRichness$CombinedNightlyBiomass*100
ShannonRichness$PercentNightlyBiomassSucker<- ShannonRichness$SuckerBiomassNightly/ShannonRichness$CombinedNightlyBiomass*100
ShannonRichness$PercentNightlyBiomassInvert<- ShannonRichness$InvertBiomass100/ShannonRichness$CombinedNightlyBiomass*100


#Count number of nights where invert biomass is greater than sucker and sturgeon combined (i.e. nightly biomass > 50% invertebrates)
length(ShannonRichness$PercentNightlyBiomassInvert) #240
Count<-subset(ShannonRichness,ShannonRichness$PercentNightlyBiomassInvert> 50)
length(Count$PercentNightlyBiomassInvert) #228
228/240 #On 228/240 nights invertebrate biomass was greater than sucker and sturgeon combined


Count<-subset(ShannonRichness,ShannonRichness$PercentNightlyBiomassSturgeon> 50)
Count[,1:10]
#Sturgeon biomass comprised >50% of total biomass on 25 and 26 May 2018




SturgeonDataFrame<-data.frame("Sturgeon",ShannonRichness$DPFS, ShannonRichness$Date2,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassSturgeon)
colnames(SturgeonDataFrame)<-c("Taxa","DPFS","Date2","Year","PercentNightly")

SuckerDataFrame<-data.frame("Sucker",ShannonRichness$DPFS, ShannonRichness$Date2,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassSucker)
colnames(SuckerDataFrame)<-c("Taxa","DPFS","Date2","Year","PercentNightly")

InvertDataFrame<-data.frame("Invertebrate",ShannonRichness$DPFS,ShannonRichness$Date2,ShannonRichness$Year,ShannonRichness$PercentNightlyBiomassInvert)
colnames(InvertDataFrame)<-c("Taxa","DPFS","Date2","Year","PercentNightly")

PlottingDataframe<-rbind(SturgeonDataFrame,SuckerDataFrame,InvertDataFrame)
PlottingDataframe[is.na(PlottingDataframe)] <- 0 #Change NA values to 0 as no taxa were observed on those dates but sampling occured

head(PlottingDataframe)
PlottingDataframe$Taxa = factor(PlottingDataframe$Taxa, levels = c("Invertebrate","Sturgeon","Sucker"))

cbPalette2 <- c("#E69F00", "#000000", "#0072B2")
cbrewer<-c("#66c2a5","#fc8d62","#8da0cb")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")

PercentBiomassPlot<-ggplot(PlottingDataframe, aes(x=Date2, y= PercentNightly,fill=Taxa))+facet_wrap(~Year)+geom_bar(stat="identity",lwd=0.1)+
  theme(legend.justification=c(1,0), legend.position=c(1,-0.0))+ylab("Nightly Biomass (%)")+xlab("Date")+
  scale_fill_manual(values=cbPalette2)+scale_color_manual(values=cbPalette2)+theme(legend.text = element_text(size = 19),legend.title = element_blank())+ 
  theme(legend.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))#+ guides(shape = guide_legend(override.aes = list(size=2)))#+geom_line()
PercentBiomassPlot

theme_set(theme_bw(base_size = 16)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/PercentBiomassPlot.tiff", width = 174, height = 174, units = 'mm', res = 1200)
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

# dev.off()
# tiff("Figures/Fig5July2020.tiff", width = 174, height = 190, units = 'mm', res = 1200)
# ggarrange(TopFamilyAbuPer100,CrayfishGAMPlot,PercentBiomassPlot,HepIsoBiomass,
#           labels = c("a", "b","c","d"),
#           ncol = 2, nrow = 2)
# dev.off()

dev.off()
tiff("Figures/Fig4Mar2022.tiff", width = 174, height = 100, units = 'mm', res = 1200)
ggarrange(TopFamilyAbuPer100,CrayfishGAMPlot,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)
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



##############
#Sucker COI
#############
SuckerCOI<-read.csv("DataClean//SuckerCOI2011-2018.csv",header=T)
head(SuckerCOI)
length(SuckerCOI$UniqueID)
SuckerCOI$Date2<-as.Date(SuckerCOI$Collection_Date,format= "%d-%B") #Date2 thinks all the dates are 2020 but 

Trtdata <- ddply(SuckerCOI, c("Date2","Year","Species_of_top_BLAST_hit"), summarise,
                 N    = length(Species_of_top_BLAST_hit)
)
unique(Trtdata$Species_of_top_BLAST_hit)


#Trtdata
COIAll<-ggplot(Trtdata,aes(x=Date2,y=N, fill=Species_of_top_BLAST_hit))+geom_bar(stat="identity",lwd=0.01,color="black")+facet_wrap(~Year)+ 
  scale_fill_discrete(name="Species",
                      labels = c(expression(italic("Acipenser fulvescens")),expression(italic("Catostomus commersonii")),"FAIL",expression(italic("Micropterus dolomieu")),expression(italic("Moxostoma anisurum")),expression(italic("Moxostoma valenciennesi")),expression(italic("Rhinichthys cataractae")),expression(italic("Semotilus atromaculatus"))))+
  xlab("Date")+ylab("Number of Individuals")+
  theme(legend.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme( legend.position=c(0.85,0.12))+theme(legend.text.align = 0,legend.title = element_blank())+scale_fill_viridis_d()

COIAll

dev.off()
tiff("Figures/COIAll.tiff", width = 174, height = 174, units = 'mm', res = 1200)
COIAll
dev.off()

#With Other category
Trtdata2 <- ddply(SuckerCOI, c("Date2","Year","Species2"), summarise,
                  N    = length(Species2)
)
unique(Trtdata2$Species2)
Trtdata2<-subset(Trtdata2,Species2!="FAIL")

#Trtdata
COIOther<-ggplot(Trtdata2,aes(x=Date2,y=N, fill=Species2))+geom_bar(stat="identity",lwd=0.01,color="black")+facet_wrap(~Year)+
  xlab("Date")+ylab("Number of Individuals")+scale_fill_manual(labels = c(expression(italic("Catostomus commersonii")),expression(italic("Moxostoma anisurum")),"Other"),values=c("#440154FF","#FDE725FF","#21908CFF"))+
  theme(legend.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme( legend.position=c(0.85,0.12))+theme(legend.text.align = 0,legend.title = element_blank())

COIOther
COIOther$theme
levels(Trtdata2)

g <- ggplot_build(COIOther)
unique(g$data[[1]]["fill"])



dev.off()
tiff("Figures/COIOther.tiff", width = 174, height = 174, units = 'mm', res = 1200)
COIOther
dev.off()

#Not faceted

Trtdata3 <- ddply(SuckerCOI, c("Date2","Species2"), summarise,
                  N    = length(Species2)
)
Trtdata3<-subset(Trtdata3,Species2!="FAIL")


COIOther2<-ggplot(Trtdata3,aes(x=Date2,y=N, fill=Species2))+geom_bar(stat="identity",lwd=0.01,color="black")+
  xlab("Date")+ylab("Number of Individuals")+scale_fill_manual(labels = c(expression(italic("Catostomus commersonii")),expression(italic("Moxostoma anisurum")),"Other"),values=c("#440154FF","#FDE725FF","#21908CFF"))+
  theme(legend.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme( legend.position=c(0.65,0.85))+theme(legend.text.align = 0,legend.title = element_blank())

COIOther2

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/COIOtherNoFacet.tiff", width = 84, height = 84, units = 'mm', res = 1200)
COIOther2
dev.off()






###Sucker CTU plot

labels<-c("Catostomus commersonii"=expression(italic("C. commersonii")),"Moxostoma anisurum"=expression(italic("M. anisurum")))
labels2<-c("Catostomus commersonii"="C. commersonii","Moxostoma anisurum"="M. anisurum")


head(SuckerCOI)
SuckerCOINOFAIL<-subset(SuckerCOI,Species2!="FAIL"&Species2!="Other")



SuckerCOICTU<-ggplot(SuckerCOINOFAIL, aes(x=SturgeonCTU,y=Species2,color=Species2,shape=Species2))+geom_jitter()+facet_grid(rows = vars(Species2),scales="free_y",switch="y",labeller=labeller(Species2=labels2))+
  theme(axis.text.y=element_blank())+theme(legend.position = "none")+
  scale_discrete_manual(aesthetics=c("color"),name="A",values=c("#440154FF","#21908CFF","#FDE725FF"),labels = c(expression(italic("C. commersonii")),expression(italic("M. anisurum"))))+theme( legend.position=c(0.75,0.82),legend.title = element_blank(),legend.text.align = 0)+
  scale_shape_manual(name="A",labels = c(expression(italic("C. commersonii")),expression(italic("M. anisurum"))),values = c(19,17))+ylab("Species")+xlab("Cumulative Temperature Units (CTU)")+theme(legend.position = "none")
SuckerCOICTU
theme_set(theme_bw(base_size = 14)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



Trtdata3 <- ddply(SuckerCOI, c("SturgeonCTU","Species2"), summarise,
                  N    = length(Species2)
)
Trtdata3<-subset(Trtdata3,Species2!="FAIL")

SuckerCOI<-subset(SuckerCOI,Species2!="FAIL")
COIOther3<-ggplot(SuckerCOI,aes(x=SturgeonCTU, fill=Species2))+geom_histogram(color="black")+
  xlab("Cumulative Temperature Units")+ylab("Number of Individuals")+scale_fill_manual(labels = c(expression(italic("C. commersonii")),expression(italic("M. anisurum")),"Other"),values=c("#440154FF","#FDE725FF","#21908CFF"))+
  theme(legend.background=element_blank())+
  theme( legend.position=c(0.65,0.85))+theme(legend.text.align = 0,legend.title = element_blank(),legend.text=element_text(size=10))
COIOther3


dev.off()
tiff("Figures/COISuckerCTU.tiff", width = 84, height = 84, units = 'mm', res = 1200)
COIOther3
dev.off()