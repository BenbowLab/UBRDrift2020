##Import Step
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
#library(mctoolsr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(lubridate)
library(ggpubr)
#detach("package:mctoolsr", unload=TRUE)
#DRIFT FILES
file.choose()
#read.csv2("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014All_AD_12.09.2017.csv",header=TRUE)
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllOtherOtu.txt",header=TRUE)
head(otufull)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllMeta.txt",header=TRUE)
taxmatrixfull=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllOtherTax.txt"))
head(taxmatrixfull)
#all taxa
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllOtu.txt",header=TRUE)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllMeta.txt",header=TRUE)
taxmatrixfull=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllTax.txt"))
metadataAll=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\DriftAllOther.txt",header=TRUE)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Mozzie\\map_mozBin2.txt",header=TRUE)#map_moz.txt

metadata$Date=as.Date(parse_date_time(metadata$Date, orders="dmy"))# fixes parsing for date function
metadata$CommonDate <- as.Date(paste0("2000-",format(metadata$Sampling_Date, "%j")), "%Y-%j")
metadata$CommonDate
head(metadata)
head(taxmatrixfull)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 16))
#sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
#head(OTU) #should be 6 taxa and 40 samples
TAX=tax_table(taxmatrixfull)
colnames(TAX)=("Family")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$NAME
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)#joins together OTU,TAX, and metadata into a 4D  object
levels(sample_data(physeq)$MoonPhase)
sample_data(physeq)$MoonPhase = factor(sample_data(physeq)$MoonPhase, levels = c("NewMoon","WaxingCrescent","FirstQuarter","WaxingGibbous","FullMoon","WaningGibbous","LastQuarter","WaningCrescent"))
levels(sample_data(physeq)$CODE)=c("NM","WXC","FQ","WXG","FM","WAG","LQ","WNC")
levels(sample_data(physeq)$NsturQ)=c("<8","9-40","41-131",">131")

levels(sample_data(physeq)$SDDQ)=c("<0.25","0.25-0.51","0.52-0.77",">0.78")
physeq#will 
sample_variables(physeq)

Richness=plot_richness(physeq, x="SDDQ",color="SDDQ", measures=c("Simpson","Shannon"))+geom_boxplot(aes(x=SDDQ, y=value, color=SDDQ), alpha=0.05)
Richness+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("SDD Quantile")
##Total Abundance
plot_bar(physeq, "MoonPhase","Abundance", "Family",facet_grid="Year~.")+xlab("Moon Phase")+ylab("Total Invertebrate Abundance")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,facet_grid="MoonPhase~."facet_grid will place each group of the specified category in a different plot



GPdist=phyloseq::distance(physeq, "bray")
adonis(GPdist ~ Date*MoonPhase, as(sample_data(physeq), "data.frame"))

##Ordellipse
ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="X72hrQuant",shape="X72hrQuant")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = X72hrQuant))+ #theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ordplot
#Relative abundance
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
RelTaxa1Per = filter_taxa(GPr, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%
RelTaxa1Per  = transform_sample_counts(RelTaxa1Per, function(x) x / sum(x) ) #transform samples based on relative abundance

plot_bar(RelTaxa1Per, "MoonPhase","Abundance", "Family")+scale_colour_manual(values=cbPalette)#+scale_fill_manual(values=cbPalette),facet_grid="MoonPhase~."facet_grid will place each group of the specified category in a different plot

#Rel abu plot
VelocityMerge=merge_samples(physeq,"MoonPhase")
sample_names(VelocityMerge)
sample_data(VelocityMerge)$VelocityBin=sample_names(VelocityMerge)
VelocityMerge=transform_sample_counts(VelocityMerge, function(x) 100*x/sum(x)) #merging samples #(averaging)
VelocityMerge
sample_data(VelocityMerge)$VelocityBin = factor(sample_data(VelocityMerge)$VelocityBin, levels = c("NewMoon","WaxingCrescent","FirstQuarter","WaxingGibbous","FullMoon","WaningGibbous","LastQuarter","WaningCrescent"))

Velocityplot=plot_bar(VelocityMerge, "VelocityBin","Abundance", fill='Family')+ylab("Total Abundance")+xlab("Moon Phase")+facet_wrap(~Family)# +scale_fill_manual(values=cbPalette) 
Velocityplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))


##Env fit
vare.mds= ordinate(physeq, "NMDS",GPdist)
#vare.mds <- metaMDS(VeganDist, trace = FALSE)
vare.mds
EnvFitMeta=metadata$SDD
ef =envfit(vare.mds, EnvFitMeta, na.rm=TRUE, permu=999)
#tax=envfit(vare.mds, EnvFitMeta,na.rm=TRUE)
ef
tax
plot(vare.mds,display="sites")
envplot=plot(ef, p.max = 0.05)


##Random Forest
ForestData=physeq#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$MoonPhase)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important Taxa for classifying  samples\n by treatment")#\n in a string tells it to start a new line
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification


levels(metadataAll$MoonPhase)= c("NewMoon","WaxingCrescent","FirstQuarter","WaxingGibbous","FullMoon","WaningGibbous","LastQuarter","WaningCrescent")
levels(metadata$SDDQ)=c("<0.25","0.25-0.51","0.52-0.77",">0.77")
levels(metadata$X72hrQuant)=c("<-1.09","-1.08_-0.19","-0.18-0.42",">0.42")
my_comparisons <- list( c("<0.25", ">0.77"), c("0.25-0.51","0.52-0.77"),c("0.25-0.51","<0.25"))# 
my_comparisons <- list( c("2015","2016"), c("2015","2017"),c("2016","2017"))# 
my_comparisons <- list(c("<-1.09","-1.08_-0.19"),c("<-1.09","-0.18-0.42"),c("-0.18-0.42",">0.42"),C("-1.08_-0.19","-0.18-0.42"))
Hep=metadataAll$Isonychiidae
Moon=metadataAll$MoonPhase


ggboxplot(metadataAll, x = "MoonPhase", y = "Heptaginiidae",
          color = "MoonPhase", palette = "jco",xlab= "Moon Phase" ,legend = "none")+ 
 geom_point()+  stat_compare_means(comparisons = my_comparisons)+# Add pairwise comparisons p-value
  stat_compare_means(label.y = 130)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Total Isonychiidae Abundance")
#ANOVA+TUKEY
SDDQ=metadata$SDDQ
Simpsons=metadata$Simpsons
Hep=metadataAll$Heptaginiidae
MoonPhase=metadata$MoonPhase
Ninverts=metadata$Ninverts
av=aov( Hep~ SDDQ)
summary(av)
posthoc <- TukeyHSD(av, 'MoonPhase', conf.level=0.95)
posthoc
#other anova option
fit=aov(Hep~MoonPhase,data=metadata)
summary(fit)
theme_set(theme_bw())
Tukey=TukeyHSD(fit,which='MoonPhase',ordered='TRUE')
plot(Tukey,las=1)
library(multcompView)

hsd <- TukeyHSD(fit) 
multcompLetters(extract_p(hsd$MoonPhase))
boxplot(Hep~MoonPhase,data=metadata, main="", xlab="Velocity", ylab="Shannon Diversity", col=MoonPhase)

##Heatmap
plot_heatmap(physeq, "NMDS","bray", "MoonPhase","Family")
#jaccard hclust
jaccCLC <- hclust(distance(physeq, "jaccard"))
plot(as.phylo(jaccCLC), show.tip.label = TRUE, tip.color = "white")
colorScale <- rainbow(length(levels(sample_data(physeq)$Year)))
cols <- colorScale[sample_data(physeq)$Year]
GP.tip.labels <- as(sample_data(physeq)$Year, "character")
tiplabels(GP.tip.labels, col = cols, frame = "none", adj = -0.05, 
          cex = 0.7)

##Family wise testing
Fam.p.table= mt(physeq,"Treatment",test="f")
print(head(Fam.p.table, 10))




#######################################################################################
##Transect Files
#####################################################################################
VeganDist=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\veganDist.txt",header=TRUE)
EnvFitMeta=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\EnvFitMetaFroude.txt",header=TRUE)
head(EnvFitMeta)
#metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014All_AD_12.09.2017.txt",header=TRUE)
head(metadata)
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014OTUOther.txt",header=TRUE)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014Meta.txt",header=TRUE)
#metadataAll=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\veganDist.txt",header=TRUE)

taxmatrixfull=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014TaxOther.txt"))
head(metadata)
head(taxmatrixfull)
head(otufull)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw())
#sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
#head(OTU) #should be 6 taxa and 40 samples
TAX=tax_table(taxmatrixfull)
colnames(TAX)=("Family")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)#joins together OTU,TAX, and metadata into a 4D  object
physeq
levels(sample_data(physeq)$IndQuant)
sample_data(physeq)$TransectGroup = factor(sample_data(physeq)$TransectGroup, levels = c("One", "Two", "Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven","Twelve","Thirteen","Fourteen","Fifteen","Sixteen","Seventeen","Eighteen","Nineteen","Twenty","Twenty-one","Twenty-two")) #fixes x-axis labels
sample_data(physeq)$SubQuart =factor(sample_data(physeq)$SubQuart, levels= c("0_0.03","0.03_0.039","0.039_0.045","0.045_0.06"))
sample_variables(physeq)



#####################
##
###################
heatmap=plot_heatmap(physeq,"NMDS","jaccard", sample.label= "med_RS","Family",sample.order = "med_RS")
heatmap
heatmap+facet_wrap(~FroudeQuart,scales = "free_x")
###################
#Network
#################
plot_net(physeq, color="FroudeQuart", shape="FroudeQuart",title="Network using Jaccard Distance",maxdist=1)
ig <- make_network(physeq, max.dist=0.3)
plot_network(ig, physeq, color= "ForudeQuart")
plot_ordination(postTreatment, ord,"samples", color="Treatment")+geom_point(size=4)+stat_ellipse(type="t",geom="polygon",alpha=0.5,aes(fill=Treatment))



sample.variables(physeq)
taxa_names(physeq)
plot_richness(physeq, x="MedQuart", measures=c("Simpson", "Shannon"))+geom_boxplot(aes(x=MedQuart, y=value, color=MedQuart), alpha=0.05)
##Total Abundance
plot_bar(physeq, "CoQuart","Abundance", "Family")+ylab("Relative Abundance (%)")+xlab("Velocity")#,facet_grid="MoonPhase~."facet_grid will place each group of the specified category in a different plot
head(sample_data(physeq))
GPdist=phyloseq::distance(physeq, "bray")
#vare.mds= ordinate(physeq, "NMDS",GPdist)
#sample_variables(physeq)
adonis(GPdist ~ MedQuart, as(sample_data(physeq), "data.frame"))
physeq=subset_samples(physeq, SampleID != "TS192")
physeq=subset_samples(physeq, SampleID != "TS16")
physeq=subset_samples(physeq, SampleID != "TS57")
physeq=subset_samples(physeq, SampleID != "TS58")
sample_variables(physeq)
##Ordellipse
ord=ordinate(physeq,"PCoA", "bray")
ordplot=plot_ordination(physeq, ord,"samples", color="MedQuart",shape="MedQuart")+geom_point(size=3)#+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = MedQuart))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




#Relative abundance
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
RelTaxa1Per = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE) #filter out any taxa lower tha 0.1%
RelTaxa1Per  = transform_sample_counts(RelTaxa1Per, function(x) x / sum(x) ) #transform samples based on relative abundance
#RelTaxa1Per
plot_bar(RelTaxa1Per, "VelocityQuart","Abundance", "Family")+scale_colour_manual(values=cbPalette)#+scale_fill_manual(values=cbPalette),facet_grid="MoonPhase~."facet_grid will place each group of the specified category in a different plot

plot_richness(physeq, x="IndQuant", measures=c("Simpson", "Shannon"))+geom_boxplot(aes(x=IndQuant, y=value, color=IndQuant), alpha=0.05)



VelocityMerge=merge_samples(physeq,"CoQuart")
sample_names(VelocityMerge)
sample_data(VelocityMerge)$VelocityBin=sample_names(VelocityMerge)
VelocityMerge=transform_sample_counts(VelocityMerge, function(x) 100*x/sum(x)) #merging samples #(averaging)
VelocityMerge
#sample_data(VelocityMerge)$VelocityBin =factor(sample_data(VelocityMerge)$VelocityBin, levels= c(">0.029","0.030-0.039","0.040-0.045", "<0.045"))
#sample_data(VelocityMerge)$VelocityBin =factor(sample_data(VelocityMerge)$VelocityBin, levels= c(">0.029","0.030-0.039","0.040-0.045", "<0.045"))

Velocityplot=plot_bar(VelocityMerge, "VelocityBin","Abundance", fill='Family')+ylab("Relative  Abundance (%)")+xlab("Substrate Index Quantile")+facet_wrap(~Family)# +scale_fill_manual(values=cbPalette) 
Velocityplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))

CoverMerge=merge_samples(physeq, "Cover_Class")
sample_names(CoverMerge)
sample_data(CoverMerge)$CoverBin=sample_names(CoverMerge)
CoverMerge=transform_sample_counts(CoverMerge, function(x) 100*x/sum(x)) #merging samples #(averaging)
sample_data(CoverMerge)$CoverBin
sample_data(CoverMerge)$CoverBin = factor(sample_data(CoverMerge)$CoverBin, levels = c("Zero", "One", "Two","Four","Five","Six","Seven","Eight","Eleven","Twelve","None"))

Coverplot=plot_bar(CoverMerge, "CoverBin","Abundance", fill='Family')+ylab("Relative  Abundance (%)")+xlab("Cover Type")+facet_wrap(~Family)# +scale_fill_manual(values=cbPalette) 
Coverplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))

DepthMerge=merge_samples(physeq,"Depth_Bin")
sample_names(DepthMerge)
sample_data(DepthMerge)$DepthBin=sample_names(DepthMerge)
DepthMerge=transform_sample_counts(DepthMerge, function(x) 100*x/sum(x)) #merging samples #(averaging)
DepthMerge
Depthplot=plot_bar(DepthMerge, "DepthBin","Abundance", fill='Family')+ylab("Relative  Abundance (%)")+xlab("Depth (m)")+facet_wrap(~Family)# +scale_fill_manual(values=cbPalette) 
Depthplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))

TransectMerge=merge_samples(physeq,"TransectGroup")
sample_data(TransectMerge)
sample_data(TransectMerge)$TransectBin = c("One", "Two", "Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven","Twelve","Thirteen","Fourteen","Fifteen","Sixteen","Seventeen","Eighteen","Nineteen","Twenty","Twenty-one","Twenty-two") #fixes x-axis labels

TransectMerge=transform_sample_counts(TransectMerge, function(x) 100*x/sum(x)) #merging samples #(averaging)
TransectMerge
sample_data(TransectMerge)$TransectBin = factor(sample_data(TransectMerge)$TransectBin, levels = c("One", "Two", "Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven","Twelve","Thirteen","Fourteen","Fifteen","Sixteen","Seventeen","Eighteen","Nineteen","Twenty","Twenty-one","Twenty-two")) #fixes x-axis labels

Transectplot=plot_bar(TransectMerge, "TransectBin","Abundance", fill='Family')+ylab("Relative  Abundance (%)")+xlab("Transect Group")+facet_wrap(~Family)# +scale_fill_manual(values=cbPalette) 
Transectplot+ theme(axis.text.x = element_text(angle = 90, hjust = 1))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))

####EnvFit 
vare.mds= ordinate(physeq, "NMDS",GPdist)
#vare.mds <- metaMDS(VeganDist, trace = FALSE)
vare.mds
EnvFitMeta=sample_data(physeq)$SubInd
ef =envfit(vare.mds, EnvFitMeta, na.rm=TRUE, permu=999)
#tax=envfit(vare.mds, EnvFitMeta,na.rm=TRUE)
ef
#tax
plot(vare.mds,display="sites")
envplot=plot(ef, p.max = 0.05)

#plot(tax, p.max =0.05)


plot(vare.mds,display="sites",type="p")

###########multiple taxa comparisons

Fam.p.table= mt(RelTaxa1Per,"SubQuart",test="f")
print(head(Fam.p.table, 10))

#########
# convert your processed phyloseq object into a dataframe
df$FamilySubInd= paste0(df$Family,df$IndQuant) #creates DateTreat variable
ggplot(df, aes(IndQuant, Abundance)) + facet_wrap(~Family)+geom_bar(stat="summary",fun.y = "mean")+
   stat_compare_means(comparisons = my_comparisons,label.y=6,label.x=5)+stat_compare_means(method="kruskal.test",label.y=10,label = "p.adj")#+theme(axis.text.x=element_blank())
#+ scale_x_continuous(limits = c(2, 6))
#stat_compare_means(label.y = 130)


  geom_bar(aes(FamilySubInd, Abundance, fill = as.factor(groups)), position = "dodge", stat = "summary", fun.y = "mean") #geom_bar(position='dodge') +
# For Sig Substrate Quart #forsiggraph=subset_taxa(physeq, Family!="Baetidae" & Family!="Brachycentridae"&Family!="Chironomidae"&Family!="Elmidae"&Family!="Empididae"&Family!="Other"&Family!="Perlidae")
forsiggraph=subset_taxa(physeq,Family =="Baetidae"|Family=="Hydropsychidae"|Family=="Lepidostomatidae"| Family=="Brachycentridae"|Family=="Gomphidae"|Family=="Chironomidae")

sample_variables(physeq)

physeq=subset_samples(physeq, SubQuart!= "NA")
physeq
df <- psmelt(RelTaxa1Per)
#head(df)
df$Abundance=df$Abundance
length(df$Abundance)
levels(df$MedQuart)=c("0_5.17","5.18_7.83","7.84_10.20","_10.21")
levels(df$CoQuart)=c("0_0.52","0.53_0.72","0.73_0.93","_0.93")
levels(df$SubQuart)=c("0_0.03","0.03_0.039","0.039_0.045","0.045_0.06")
levels(df$TransectGroup) = c("One", "Two", "Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven","Twelve","Thirteen","Fourteen","Fifteen","Sixteen","Seventeen","Eighteen","Nineteen","Twenty","Twenty-one","Twenty-two") #fixes x-axis labels

p <- ggbarplot(df, x = "MedQuart", y = "Abundance",add = c("mean_se"),
              color = "black", palette = "cbPalette", 
              line.color = "gray", line.size = 0.4,
              facet.by = "Family", short.panel.labs = FALSE,label.y=7, p.adjust.method = "bonferroni")
p+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Median Rock Size")+ylab("Relative Abundance (SE)")#+ylim = c(0, 20)
my_comparisons <- list(c("0_0.03","0.03_0.039"))
my_comparisons <- list(c("0_0.03","0.039_0.045"),c("0.03_0.039","0.039_0.045"), c("0_0.03","0.045_0.06"),c("0.039_0.045","0.045_0.06"))
my_comparisons <- list(c("0_0.03","0.045_0.06"))
my_comparisons <- list(c("0_0.03","0.03_0.039"))
p=p + stat_compare_means(comparisons = my_comparisons)#+ stat_compare_means(aes(group = SubQuart,label = ..p.adj..))#
p #+ annotate("text", label = "lab" , size = 4, x = 1, y = )

data=compare_means(Abundance~SubQuart, df, method = "anova", paired = FALSE,
              group.by = "Family", ref.group = NULL, symnum.args = list(),
              p.adjust.method = "bonferroni")
data
cdata$Family
lab=data$p.adj
AdjP=c("Adj P < 0.001","Adj P = 0.020", "Adj P = 0.019","Adj P < 0.001","Adj P < 0.001","Adj P < 0.001"," Adj P = 0.002")
sample_data(MozziePhyloseq)$DateTreat= paste0(sample_data(MozziePhyloseq)$Treatment,sample_data(MozziePhyloseq)$Sampling_Date) #creates DateTreat variable


cdata <- ddply(df, c("Family", "SubQuart"), summarise,
               N    = length(Abundance),
               mean = mean(Abundance),
               sd   = sd(Abundance),
               se   = sd / sqrt(N)
)
df
head(cdata)
#cdata$mean
# plot bar graph with standard deviation as error bars
ggplot(cdata, aes(SubQuart, mean)) + facet_wrap(~Family)+geom_point()+ #geom_bar(position='dodge') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))+ annotate("text", label = AdjP, size = 4, x = 2.5, y = 6)+
  geom_text(aes(label=N), vjust=-1.75)

head(cdata)

# barplots! (replace SampleType with treatment)
Plot=ggplot(cdata, aes(x=SubQuart, y=mean,lab=SubQuart))+facet_wrap(~Family)
Plot+geom_bar(stat="identity")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se),col="black")+  geom_text(aes(label=N), vjust=-1.75)
  theme(legend.justification=c(1,0), legend.position=c(1,0))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Relative Abundance (SE)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Substrate Index")+stat_compare_means(comparisons = my_comparisons)

my_comparisons <- list(c(">0.029","0.030-0.039"),c("0.040-0.045","<0.045")) #,c("-0.18-0.42",">0.42"),C("-1.08_-0.19","-0.18-0.42"))


ggboxplot(df, x = "SubQuart", y = "Abundance",color = "SubQuart", palette = "jco",xlab= "SubQuart")+ facet_wrap(~Family)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  geom_point()+  stat_compare_means(comparisons = my_comparisons)+# Add pairwise comparisons p-value
  stat_compare_means(label.y = 130)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Total Isonychiidae Abundance")


Plot=  ggplot(df, aes(x=SubQuart, y=Abundance,lab=SubQuart))+facet_wrap(~Family) +
    geom_signif(comparisons = split(t(combn(levels(df$SubQuart), 2)), seq(nrow(t(combn(levels(df$SubQuart), 2))))), 
              map_signif_level = TRUE,step_increase=-0.1)
Plot+geom_bar( position = "dodge", stat = "summary", fun.y = "mean")
newdata <- subset(df, Family=="Ephemerellidae" ,
                  select=OTU:Med_Rock_Size)

head(newdata)
levels(df$SubQuart)=c("0_0.03","0.03_0.039","0.039_0.045","0.045_0.06")
library(multcompView)
fit=aov(Abundance~SubQuart,data=newdata)
hsd <- TukeyHSD(fit) 
hsd$SubQuart
multcompLetters(extract_p(hsd$SubQuart))


library(ggsignif)
################################################################
#distance methods

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods <- dist_methods[-(1:3)]
dist_methods["designdist"]
dist_methods=dist_methods[-which(dist_methods=="ANY")]
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(physeq, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeq, "PCoA", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeq, iMDS, color="VelocityBin")
  # Add title to each plot
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=VelocityBin))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("Distance metrics for Benthos data ")
p

print(plist[["bray"]])

#####################################
#Alpha diversity
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Sturgeon\\Transect2014All.txt",header=TRUE)

library(gridExtra)
library(ggsignif)
head(metadata)
SDI=metadata$SDI
Ephemerellidae=metadata$Ephemerellidae
FroudeQuart=metadata$FroudeQuart
VBin=metadata$VelocityBin

av=aov( Ephemerellidae~ FroudeQuart)
summary(av)
posthoc <- TukeyHSD(av, 'FroudeQuart', conf.level=0.95)
posthoc

fit=aov(SDI~VelocityBin,data=metadata)
theme_set(theme_bw())
Tukey=TukeyHSD(fit,which='VBin',ordered='TRUE')
plot(Tukey,las=1)
library(multcompView)
fit=aov(SDI~VelocityBin,data=metadata)
hsd <- TukeyHSD(fit) 
multcompLetters(extract_p(hsd$VelocityBin))
boxplot(SDI~VelocityBi,data=metadata, main="", xlab="Velocity", ylab="Shannon Diversity", col=VBin)
  
  plot=(ggplot(metadata, aes(MoonPhase,Shannon))+stat_summary(fun.y=mean,geom="point", size=2)
        +stat_summary(fun.data=mean_se,geom="errorbar")+xlab("Froude Quartile")+ylab("Abundance Ephemerellidae"))# +geom_signif(comparisons = list(c("0.51-0.75", "0-0.25")), map_signif_level=TRUE) 
 # +geom_signif(comparisons = list(c("0.75-1.0", "1.0-1.25")), map_signif_level=TRUE)
  #+geom_signif(comparisons = list(c("0.51-0.75", "0.75-1.0")), map_signif_level=TRUE))
 plot 
##############

list=compare_means(Simpsons ~ MoonPhase, data = metadata)
list$p.signif
my_comparisons <- list( c("FM", "FQ"), c("LQ","FQ"),c("WXC","FQ"), c("WNC", "WAG") )


ggboxplot(metadata, x = "CODE", y = "Simpsons",
          color = "CODE", palette = "jco",xlab= "Moon Phase" ,legend = "none")+geom_point()+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.41)+geom_point()+ylab("Simpsons Diversity")+
    stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons p-value
   # Add global p-value




ForestData=physeq#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$FroudeQuart)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying  samples\n by treatment")#\n in a string tells it to start a new line
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification






################################################################################33
################### Example Code from other projects
###################################################################################
##filtering
sample_variables(physeq)
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
MozziePhyloseq = filter_taxa(GPr, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%
MozziePhyloseq  = transform_sample_counts(MozziePhyloseq, function(x) x / sum(x) )#transform samples so they are based on relative abundance
MozziePhyloseq <- subset_taxa(MozziePhyloseq, Family != "mitochondria" & Class != "Chloroplast")

MozziePhyloseq # you will see a difference in the number of taxa compared to the physeq option above
plot_bar(MozziePhyloseq)# will plot the bar graph based on abundances 
#Note: for groups with different numbers of samples the overall abundance will look different)
plot_bar(MozziePhyloseq, "Sampling_Date","Abundance", "Phylum",facet_grid="Instar~.")#facet_grid will place each group of the specified category in a different plot

Instar=merge_samples(MozziePhyloseq,"Instar")
sample_names(Instar)
sample_data(Instar)$Instar=sample_names(Instar)
#sample_data(Instar)
Instar=transform_sample_counts(Instar, function(x) 100*x/sum(x)) #merging samples #(averaging)
Instar
Instarplot=plot_bar(Instar, "Instar","Abundance", fill='Family') +scale_fill_manual(values=cbPalette) +ylab("Relative Bacterial Abundance (> 1%)")
Instarplot+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))
sample_data(MozziePhyloseq)$Sampling_Date




plot_bar(MozziePhyloseq, "Phylum","Abundance")#, "Sampling_Date",facet_grid="Treatment~.")
Chlorop=subset_taxa(MozziePhyloseq,Class=="Chloroplast")#subset out taxa that have the class Chloroplast (can replace chloroplast with whatever taxa you want)
Chlorop=subset_samples(Chlorop, Sampling_Date=="4/18/15")
Chlorop=subset_samples(Chlorop, Treatment!="BTI")
Chlorop
GPdist=phyloseq::distance(Chlorop, "wunifrac")
MONMDS= ordinate(Chlorop, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(Chlorop), "data.frame"))


PostTreatment=subset_samples(MozziePhyloseq, Sampling_Date=="4/18/15")
PostTreatment
plot_bar(Chlorop, "Treatment", "Abundance", "Order", facet_grid="Sampling_Date~.",title= "Chloroplast Abundance")

## NMDS
set.seed(1)
GP.ord=ordinate(MozziePhyloseq, "NMDS", "wunifrac") #multidimentional scaling

p1=plot_ordination(MozziePhyloseq,GP.ord, type="samples", color="Instar",shape="Instar")#plot of ordination by taxa,color="Treatment",shape="Treatment"
p1 #+ facet_wrap(~Phylum,3)+ geom_point(size=3)
print(p1)
p2=plot_ordination(MozziePhyloseq,GP.ord,type="samples",shape="Sampling_Date", color="Sampling_Date",facet_grid="Sampling_Date~.")
p2 + geom_point(size=5) + ggtitle("NMDS by Sampling Date"))

#PCoA
ordu = ordinate(MozziePhyloseq, "PCoA", "")
PCoA1=plot_ordination(MozziePhyloseq, ordu, color="Instar",shape="Treatment", title="PCoA")
PCoA1
PCoA1 + geom_point(size=5) + ggtitle("PCoA") + facet_wrap(~Sampling_Date, 3)
plot_ordination(MozziePhyloseq, ordu,color="Sampling_Date")


Instar2=subset_samples(MozziePhyloseq,Instar=="2nd")
Instar3=subset_samples(MozziePhyloseq,Instar=="3rd")
Instar4=subset_samples(MozziePhyloseq,Instar=="4th")
Instar34=subset_samples(MozziePhyloseq, Instar !="2nd")
Instar34
Q1v3=subset_samples(MozziePhyloseq, DensityQ !=2)
Q1v3
################################
##to vegan adonis (different package used for PERMANOVA)
data=ControlBtipost
GPdist=phyloseq::distance(data, "wunifrac")
MONMDS= ordinate(data, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(data), "data.frame"))
adonis(GPdist ~Instar*DensityBin*Treatment, as(sample_data(data), "data.frame"))

head(sample_data(MozziePhyloseq))
adonis(GPdist ~Sampling_Date+Instar+DensityBin+Treatment, as(sample_data(data), "data.frame"))
adonis(GPdist ~Treatment+DensityBin+Instar+Sampling_Date, as(sample_data(data),"data.frame"))
adonis(GPdist ~Instar*DensityBin*Treatment, as(sample_data(data), "data.frame"))

Instar2.3=subset_samples(MozziePhyloseq, Instar!=4)
Instar2.3
##Random forest (used to determine what bacterial taxa differ between different groups)
ControlMethpostglom=tax_glom(ControlMethpost, "Order")

ControlMethpost
ForestData=physeq#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$Instar)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying  samples\n by treatment")#\n in a string tells it to start a new line
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
g#heatmap
plot_heatmap(physeq,"NMDS","wunifrac", sample.label= "Instar",sample.order="Instar")
plot_bar(TreatmentSumPhylo, fill="TAXA","TAXA",facet_grid=Sampling_Date~Treatment)
###################
otu=read.table(file.choose(),header=TRUE)#L6otu.txt
metadata=read.table(file.choose(),header=TRUE)#map_moz.txt
taxmatrix=as.matrix(read.table(file.choose(),header=TRUE))#TAXAtableL6.txt

##Taxa summary import

otu=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\SumByTreatmentOTU.txt",header=TRUE)#SumbyTreatmentOTU.txt
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Tree divot project\\TreatmentSumMeta.txt",header=TRUE)#treatmentsummeta
taxmatrix=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\SumByTreatmentTAXA.txt",header=TRUE))#sumbytreatm
#head(otu)
#head(taxmatrix)
OTU=otu_table(otu, taxa_are_rows=TRUE)
TAX1=tax_table(taxmatrix)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$ID
sample_names(sampdat)=(sample_names(OTU))
taxa_names(TAX1)= c("Acidobacteria","Actinobacteria","Bacteroidetes",	"Cyanobacteria",	"Firmicutes",	"Alphaproteobacteria",	"Betaproteobacteria",	"Deltaproteobacteria",	"Gammaproteobacteria",	"Other")
row.names(OTU)=taxa_names(TAX1)
#head(TAX1)
#row.names(OTU)
#OTU
TreatmentSumPhylo=phyloseq(TAX1,OTU, sampdat)
TreatmentSumPhylo
TreatmentSumPhylo  = transform_sample_counts(TreatmentSumPhylo, function(x) x / sum(x) )
plot_bar(TreatmentSumPhylo, "Sampling_Date","Abundance","TAXA",facet_grid="Treatment~.")+ xlab("Sampling Date")+ylab("Mean Relative Abundance")
Control=subset_samples(MozziePhyloseq,Treatment=="Control")
Control
GPdist=phyloseq::distance(MozziePhyloseq, "bray")
MONMDS= ordinate(MozziePhyloseq, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(MozziePhyloseq), "data.frame"))
adonis(GPdist ~ Sampling_Date*Treatment, as(sample_data(MozziePhyloseq), "data.frame"),permutations=999)
adonis(GPdist ~ Sampling_Date*Density*HeadCapsule*Treatment, as(sample_data(MozziePhyloseq), "data.frame"))



#control vs bti
ControlBti=subset_samples(MozziePhyloseq,Treatment!="Methoprene")
GPdist=phyloseq::distance(ControlBti, "bray")
MONMDS= ordinate(ControlBti, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(ControlBti), "data.frame"))
adonis(GPdist ~Sampling_Date*Treatment, as(sample_data(ControlBti), "data.frame"))
#Control vs methoprene
ControlMeth=subset_samples(MozziePhyloseq,Treatment!="BTI")
GPdist=phyloseq::distance(ControlMeth, "bray")
MONMDS= ordinate(ControlMeth, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(ControlMeth), "data.frame"))
adonis(GPdist ~ Sampling_Date*Treatment, as(sample_data(ControlMeth), "data.frame"))
#Bti vs methoprene
BtiMeth=subset_samples(MozziePhyloseq,Treatment!="Control")
GPdist=phyloseq::distance(BtiMeth, "bray")
MONMDS= ordinate(BtiMeth, "NMDS",GPdist)
adonis(GPdist ~ Treatment, as(sample_data(BtiMeth), "data.frame"))
adonis(GPdist ~ Sampling_Date*Treatment, as(sample_data(BtiMeth), "data.frame"))
##vegan ordellipse
Samp4.18=subset_samples(MozziePhyloseq,Sampling_Date=="4/18/2015")
Samp4.18
Samp5.2=subset_samples(MozziePhyloseq,Sampling_Date=="5/2/2015")
Samp5.8=subset_samples(MozziePhyloseq,Sampling_Date=="5/8/2015")
names(metadata)
Samp4.18.NMDS=phyloseq::distance(Samp4.18, "bray")
ef=envfit(GPdist, metadata,permutations = 999)
ef
sample_variables(MozziePhyloseq)
sample_data(MozziePhyloseq)$"Date"=sample_data(MozziePhyloseq)$"Sampling_Date"

##ordellipse phyloseq
ord=ordinate(MozziePhyloseq,"PCoA", "wunifrac")
ordplot=plot_ordination(MozziePhyloseq, ord,"samples", color="Instar", shape="Instar")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/2, aes(fill = Instar))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ord4.18=ordinate(Samp4.18, "NMDS","bray")
ord4.18plot= plot_ordination(Samp4.18,ord4.18, "samples",color="Treatment")

##plot network
plot_net(MozziePhyloseq, color="Instar", shape="Instar",title="Network using Jaccard Distance",maxdist=1)
ig <- make_network(MozziePhyloseq, max.dist=0.3)
plot_network(ig, MozziePhyloseq, color= "Sampling_Date")
plot_ordination(postTreatment, ord,"samples", color="Treatment")+geom_point(size=4)+stat_ellipse(type="t",geom="polygon",alpha=0.5,aes(fill=Treatment))



##various distance methods
dist_methods <- unlist(distanceMethodList)
dist_methods <- dist_methods[-(1:3)]
dist_methods = dist_methods[-which(dist_methods=="ANY")]
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(physeq, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeq, "PCoA", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeq, iMDS, color="Treatment", shape="Sampling_Date")
  # Add title to each plot
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Treatment, shape=Sampling_Date))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p
p
print(plist[["e"]])
##try unifrac
betadisper(GPdist, metadata, type = c("median"), bias.adjust = FALSE,
           sqrt.dist = FALSE, add = FALSE)
##betadisper
sampledd <- data.frame(sample_data(sampdat))

beta <- betadisper(GPdist, sampledd$)
permutest(beta)
##CAP Ordination
cap_ord <- ordinate(
  physeq = MozziePhyloseq, 
  method = "CAP",
  distance = GPdist,
  formula = ~ Treatment)
cap_plot <- plot_ordination(
  physeq = MozziePhyloseq, 
  ordination = cap_ord, 
  color = "Sampling_Date") + 
  aes(shape = Treatment)
anova(cap_ord)
GPdisc=phyloseq::distance(NoChloro,method="bray")
adonis(GPdisc ~ Sampling_Date*Density*Treatment, as(sample_data(MozziePhyloseq), "data.frame"),permutations=999)
adonis(GPdist ~ Treatment, as(sample_data(MozziePhyloseq), "data.frame"),permutations=999)
ChloroSummary=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\ChloroplastSummary.txt",header=TRUE)

ChloroSummary
ChloroSummary$Chloroplast=ChloroSummary$Chloroplast*100
#ggplot(ChloroSummary, aes(x=Date, y=Mean,color=Treatment, shape=Treatment,group=Treatment))+geom_errorbar(aes(ymin=Mean-Sem/2,ymax=Mean+Sem/2),position = dodge,color="black",width=1)+geom_line(size=1)+geom_point(size=3)+ylab("Rel Abundance Algae (%)")


cdata <- ddply(ChloroSummary, c("Treatment", "Date"), summarise,
               N    = length(Chloroplast),
               mean = mean(Chloroplast),
               sd   = sd(Chloroplast),
               se   = sd / sqrt(N)
)
cdata
cdata$Date=c("1D Post","7D Post","14D Post","1D Post","7D Post","14D Post","1D Post","7D Post","14D Post")
# Set the "anchoring point" of the legend (bottom-left is 0,0; top-right is 1,1)
# Put bottom-left corner of legend box in bottom-left corner of graph
cdata$Date = factor(cdata$Date, levels = c("1D Post", "7D Post", "14D Post")) #fixes x-axis labels

cdataplot=ggplot(cdata, aes(x=Date, y=mean,color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=1.5,linetype="dashed")+geom_point(size=6)+ylab("Rel Abundance Algal Taxa (%)")+xlab("")+scale_colour_manual(values=cbPalette)
cdataplot+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Enterobac across dates
Enterobac=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\EnteroMeanSem.txt",header=TRUE)
Enterobac
dodge <- position_dodge(width=0.05)#ggplot position dodge
ggplot(Enterobac, aes(x=Date, y=mean,color=Treatment, shape=Treatment,group=Treatment))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),position = dodge,color="black",width=1)+geom_line(size=1)+geom_point(size=3)+ylab("Rel Abundance Enterobacteriaceae")



##Alpha Richness
AlphaRichness=file.choose()
AlphaRichness
AlphaRichness=read.csv("C:\\Users\\Joe Receveur\\Downloads\\alphadiversity (2).csv",header=TRUE)
head(AlphaRichness)
Simpson= ddply(AlphaRichness, c("POST2"),summarise,N    = length(AlphaRichness),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)
#plot_richness(physeq, x="Instar",color="Treatment",shape="Treatment", measures=c("Chao1", "Shannon","Simpson"),)+geom_boxplot(aes(x=Instar, y=value, color=Treatment), alpha=0.05)
ggplot(AlphaRichness, aes(x=Treatment, y=Average,color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=Average-SE,ymax=Average+SE),color="black",width=1)+geom_line(size=1)+geom_point(size=3)+ylab("Relative Abundance")+xlab("Sampling Date")
Simpson
ggplot(Simpson, aes(x=Treatment, y=mean,color=Treatment, shape=Treatment))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=3)+geom_point(size=3)+facet_wrap(~Date,1)

GammaSum=ddply(SUMbyTreatment, c("Treatment", "Sampling.Date"), summarize,N=length(SUMbyTreatment), meanG=mean(Gammaproteobacteria), sdG=sd(Gammaproteobacteria),seG=sdG/sqrt(N),meanF=mean(Firmicutes), sdF=sd(Firmicutes),seF=sdF/sqrt(N))

GammaSum#Sampling_Date	Treatment	Acidobacteria	Actinobacteria	Bacteroidetes	Cyanobacteria	Firmicutes	Alphaproteobacteria	Betaproteobacteria	Deltaproteobacteria	Gammaproteobacteria	Other
SUMbyTreatment=read.table(file.choose(),header=TRUE)
head(SUMbyTreatment)
ggplot(SUMbyTreatment, aes(x=Sampling.Date, y=mean,color=Taxa, shape=Taxa, group=Taxa))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=1)+geom_point(size=3)+facet_wrap(~Treatment,1)+ylab("Relative Abundance")+xlab("Sampling Date")
ggplot(SUMbyTreatment, aes(x=Sampling.Date, y=mean,color=Taxa, shape=Taxa, group=Taxa))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=1)+geom_point(size=3)+facet_wrap(~Treatment,1)+ylab("Relative Abundance")+xlab("Sampling Date")

MozDensity=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\MozDensity.txt",header=TRUE)
MozDensity
MozDensity$Date=c("7D Pre", "7D Pre", "7D Pre","1D Post","1D Post","1D Post","7D Post","7D Post","7D Post","14D Post","14D Post","14D Post")
MozDensity$Date = factor(MozDensity$Date, levels = c("7D Pre","1D Post", "7D Post", "14D Post")) #fixes x-axis labels
dodge=0.01
Density=ggplot(MozDensity, aes(x=Date, y=Density,color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=Density-SE_mean,ymax=Density+SE_mean),color="black",width=1,position = dodge)+geom_line(size=1.5,linetype="dashed")+geom_point(size=5)+ylab("Larval Density (per 1 L)")+xlab("")+scale_colour_manual(values=cbPalette)
Density+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")
MozInstar=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\HeadCapsule.txt",header=TRUE)
MozInstar
MozInstar$Date=c("1D Post","1D Post","1D Post","7D Post","7D Post","7D Post","14D Post","14D Post","14D Post")
MozInstar$Date = factor(MozInstar$Date, levels = c("1D Post", "7D Post", "14D Post")) #fixes x-axis labels

InstarPlot=ggplot(MozInstar, aes(x=Treatment, y=Instar,color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=Instar-HCSEM,ymax=Instar+HCSEM),color="black",width=1,position = dodge)+geom_line(size=1,linetype="dashed")+geom_point(size=5)+ylab("Mean Instar (SE)")+xlab("")+facet_wrap(~Date,1)+ theme(legend.justification=c(1,0), legend.position=c(1,0))+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_colour_manual(values=cbPalette)
InstarPlot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
VolTemp=read.csv("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\TemppHVolR.csv",header=TRUE)
VolTemp
VolTemp$Date = factor(VolTemp$Date, levels = c("7D Pre","1D Post", "7D Post", "14D Post")) #fixes x-axis labels

VolPlot=ggplot(VolTemp, aes(x=Date, y=Mean, shape=Var, group=Var))+geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),color="black",width=1,position = dodge)+geom_line(size=1.5,linetype="dashed")+geom_point(size=4)+ylab("")+xlab("Sampling Date")+facet_wrap(~Var,1,scales="free") + theme(legend.position="none")
VolPlot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
##faith's distance
faith=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\faithDwithMeta.txt",header=TRUE)
faithMeans= ddply(faith, c("Treatment","Date"),summarise,N    = length(faith),
                  mean = mean(PD_whole_tree),
                  sd   = sd(PD_whole_tree),
                  se   = sd / sqrt(N))
faithMeans
ggplot(faithMeans, aes(x=Date,y=mean,color=Treatment,shape=Treatment,group=Treatment))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=1)+geom_point(size=3)
####labdsv
library(labdsv)
library(vegan)
head(otufull)
dim(otufull)
colnames(taxmatrixfull) <- c("Domain", "Phylum", "Class", "Order", "Family","Genus")
veg=t(otufull)
head(taxmatrixfull)

site=metadata
site=subset(site, Treatment=="Methoprene")#how to subset samples
dim(veg)
## presence 
#cor(Density,HeadCapsule) correlation of two numeric vars
abuocc(veg)


#print(taxmatrixfull[374,])
#attach(site)
vegtab(veg,pltlbl=(site$Sampling_Date))

const(veg,site$Sampling_Date,min=0.1)#looks only at presence absence
importance(veg,site$Sampling_Date,min=1)#uses rel abundance

speciesDiscrimination <- importance(veg,site$Sampling_Date)
spcdisc(speciesDiscrimination,sort=TRUE)
detach(site)
#bray nmds using labdsv

dis.bray <- vegdist(veg,method="bray") 
#disana(dis.bray)
dis.bin <- dist(veg,"binary")
bin.pco <- pco(dis.bin,k=10)
#plot(bin.pco)
plot(bin.pco,title="Sampling_Date")
surf(bin.pco,site$Sampling_Date)
bc.nmds=nmds(dis.bray)
plot(bc.nmds)
bc4d.nmds <- nmds(dis.bray,4)
ordcomp(bc4d.nmds,dis.bray,dim=4)
ordcomp(bc.nmds,dis.bray,dim=2)
plot(bc.nmds)
surf(bc.nmds,Temp)
##indicator species analysis 
library(indicspecies)
Samp4.18=subset_samples(MozziePhyloseq,Sampling_Date=="4/18/2015")
#Samp4.18
Samp5.2=subset_samples(MozziePhyloseq,Sampling_Date=="5/2/2015")
Samp5.8=subset_samples(MozziePhyloseq,Sampling_Date=="5/8/2015")
#Samp5.8

#print(otu_table(Samp5.2))
#print(otu_table(Samp5.8))
groups4.18=c(3,1,1,3,2,2,2,3)#group 1=Bti 2=Control 3=Meth
groups4.18
ISA4.18=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\4.18ISA.txt",header=TRUE)
#(ISA4.18)
indval=multipatt(ISA4.18,groups4.18,func="r.g",control=how(nperm=999))
summary(indval,alpha=0.05)
ISA5.2=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\5.2ISA.txt",header=TRUE)
#head(ISA5.2)
groups5.2=c(1,2,3,3,3,2,3,1,3,3,3,3,1,1,2,2,1,2)
#groups5.2
indval2=multipatt(ISA5.2,groups5.2,func="r.g",control=how(nperm=999))
summary(indval2, alpha=0.05)
#options(max.print=1000)
ISA5.8=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\5.8ISA.txt",header=TRUE)
#head(ISA5.8)
groups5.8=c(1,3,2,3,3,3,3,1,2,1,3,1,1,2)
indval3=multipatt(ISA5.8,groups5.8,func="r.g",control=how(nperm=999))
summary(indval3)
print(taxa_names(MozziePhyloseq))
r <- rownames(tax_table(MozziePhyloseq))
kable(tax_table(MozziePhyloseq)[r, ])
##

#########Family level
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\L5otu.txt",header=TRUE)#L6otu.txt
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Tree divot project\\map_moz.txt",header=TRUE)#map_moz.txt
taxmatrixfull=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\L5TAXA.txt"))#TA.txt
#head(metadata)
##alternative import
#biome=import_biom(file.choose())
#physeq2=phyloseq(biome,sampdat)
##formatting
theme_set(theme_bw()) #sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
#head(OTU) #should be 6 taxa and 40 samples
TAX=tax_table(taxmatrixfull)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
#head(TAX) #will say six taxa by six ranks, theres more in the file it is just counting the header
#rownames(OTU)
#(TAX)
colnames(TAX) <- c("Domain", "Phylum", "Class", "Order", "Family")#assigns names to the taxa levels
#taxa_names(TAX)
#row.names(OTU)
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)#joins together OTU,TAX, and metadata into a 4D  object
physeq#will return an overview of the physeq object, check that it has 40 samples, 441 taxa by six taxanomic rank, if not you need to troubleshoot the import steps
#plot_richness(physeq)#plot richness before filtering to remove singletons or low samples
#plot_richness(physeq, x="Sampling_Date",color="Treatment",shape="Treatment", measures=c("Chao1", "Shannon","Simpson"),title="Richness by Date")
##filtering
sample_variables(physeq)
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
MozziePhyloseq = filter_taxa(GPr, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower than
MozziePhyloseq  = transform_sample_counts(MozziePhyloseq, function(x) x / sum(x) )#transform samples so they are based on relative abundance
MozziePhyloseq <- subset_taxa(MozziePhyloseq, Family != "mitochondria" & Class != "Chloroplast")


#Family Level Graphs
sample_data(MozziePhyloseq)$DateTreat= paste0(sample_data(MozziePhyloseq)$Treatment,sample_data(MozziePhyloseq)$Sampling_Date) #creates DateTreat variable
head(sample_data(MozziePhyloseq))
DateTreatMerge=merge_samples(MozziePhyloseq, "DateTreat") #creates 4d object with DateTreat column
#sample_variables(DateTreatMerge)
sample_data(DateTreatMerge)
DateTreatMerge=transform_sample_counts(DateTreatMerge, function(x) 100 *x/sum(x)) #merging samples #(averaging)
(sample_data(DateTreatMerge))
sample_data(DateTreatMerge)$DateTreat=sample_names(DateTreatMerge)
sample_data(DateTreatMerge)$Treatment=c("BTI","BTI","BTI","Control","Control","Control","Methoprene","Methoprene","Methoprene")
sample_data(DateTreatMerge)$Date=c("1D Post", "7D Post", "14D Post","1D Post", "7D Post", "14D Post","1D Post", "7D Post", "14D Post" )
sample_data(DateTreatMerge)$Date = factor(sample_data(DateTreatMerge)$Date, levels = c("1D Post","7D Post", "14D Post")) #fixes x-axis labels
DateTreatMerge
DateTreatMergeglom=tax_glom(DateTreatMerge, "Genus")

sample_variables(DateTreatMergeglom)
#DateTreatMergeglom=transform_sample_counts(DateTreatMergeglom, function(x) 100 *x/sum(x)) #merging samples #(averaging)
DateTreatMergeglom=subset_samples(DateTreatMergeglom, Treatment!="BTI")
DateTreatMergeglom=subset_samples(DateTreatMergeglom, Date!="1D Post")
Genusplot=plot_bar(DateTreatMerge, "Treatment","Abundance", "Genus")+xlab("")+ylab("Relative Bacterial Abundance (> 1%)") +facet_grid(Date~.)+ scale_fill_manual(values=cbPalette)
Genusplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p$data
p=plot_bar(DateTreatMergeglom, "Treatment","Abundance", "Genus",facet_grid="Date")+xlab("Treatment")+ylab("Relative Abundance (> 1%)")+ scale_fill_manual(values=cbPalette) 
p
p+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))
#Top OTUs
####Chloroplast Graph
sample_data(Chlorop)$DateTreat= paste0(sample_data(MozziePhyloseq)$Treatment,sample_data(MozziePhyloseq)$Sampling_Date) #creates DateTreat variable
DateTreatChlorop=merge_samples(Chlorop, "DateTreat") #creates 4d object with DateTreat column
sample_data(DateTreatChlorop)$DateTreat=sample_names(DateTreatChlorop)
sample_data(DateTreatChlorop)
DateTreatChlorop=transform_sample_counts(DateTreatChlorop, function(x) 100 *x/sum(x)) #merging samples #(averaging)
sample_data(DateTreatChlorop)$Treatment=c("BTI","BTI","BTI","Control","Control","Control","Methoprene","Methoprene","Methoprene")
sample_data(DateTreatChlorop)$Date=c("1D Post", "7D Post", "14D Post","1D Post", "7D Post", "14D Post","1D Post", "7D Post", "14D Post" )
sample_data(DateTreatChlorop)$Date = factor(sample_data(DateTreatChlorop)$Date, levels = c("1D Post","7D Post", "14D Post")) #fixes x-axis labels

DateTreatChloropglom=tax_glom(DateTreatChlorop, "Order")
p=plot_bar(DateTreatChloropglom, "Treatment","Abundance","Order",facet_grid="Date~.")+xlab("")+ylab("Rel Abundance Algal Taxa (> 1%)") +scale_fill_manual(values=cbPalette) 
p+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(ggplot2)

#Indicator Species analysis
library(indicspecies) #don't load this and randomforest at same time (indicator function masked)
Samp4.18=subset_samples(MozziePhyloseq,Sampling_Date=="4/18/2015")
#Samp4.18
Samp5.2=subset_samples(MozziePhyloseq,Sampling_Date=="5/2/2015")
Samp5.8=subset_samples(MozziePhyloseq,Sampling_Date=="5/8/2015")
#options(max.print=1000000)
#print(otu_table(Samp4.18))
#print(otu_table(Samp5.2))
#print(otu_table(Samp5.8))
groups4.18=c(3,1,1,3,2,2,2,3)#group 1=Bti 2=Control 3=Meth
groups4.18
ISA4.18=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\4.18FamilyISA.txt",header=TRUE)
#(ISA4.18)
indval=multipatt(ISA4.18,groups4.18,func="r.g",control=how(nperm=999))
summary(indval,alpha=0.1)
ISA5.2=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\5.2FamilyISA.txt",header=TRUE)
#head(ISA5.2)
groups5.2=c(1,2,3,3,3,2,3,1,3,3,3,3,1,1,2,2,1,2)
#groups5.2
indval2=multipatt(ISA5.2,groups5.2,func="r.g",control=how(nperm=999))
summary(indval2, alpha=0.05)
#options(max.print=1000)
ISA5.8=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\5.8FamilyISA.txt",header=TRUE)
#head(ISA5.8)
groups5.8=c(1,3,2,3,3,3,3,1,2,1,3,1,1,2)
indval3=multipatt(ISA5.8,groups5.8,func="r.g",control=how(nperm=999))
summary(indval3)
#print(taxa_names(MozziePhyloseq))
r <- rownames(tax_table(MozziePhyloseq))
kable(tax_table(MozziePhyloseq)[r, ])
##


library(mctoolsr)
####mctoolsr
taxaTable="C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Mozzie\\merged_moz.biom"
editedTaxa=("C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Mozzie\\otu_table_exportp__removed.txt")
metadata="C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Mozzie\\map_mozBin3.txt"
input=load_taxa_table(editedTaxa,metadata)
input$map_loaded
Post=filter_data(input,filter_cat='Sampling_Date',keep_vals= c("7 D Post","14 D Post"))
PostMeth=filter_data(Post,filter_cat='Treatment',keep_vals= c('Control','Methoprene'))
Post2week=filter_data(input,filter_cat='Sampling_Date',keep_vals= c('5/8/2015'))

taxa_sum_phyla=summarize_taxonomy(PostMeth,level=6, report_higher_tax = FALSE)
taxa_summary_by_sample_type(taxa_sum_phyla, Post$map_loaded, 
                            type_header = 'Treatment', filter_level = 0.05, test_type = 'KW')
#taxa_sum_phyla
plot_ts_heatmap(taxa_sum_phyla, input$map_loaded,0.001, 'Treatment')
#export_taxa_table(input, "C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Mozzie\\otu_table_export.txt")
input_proteobact = filter_taxa_from_input(input, taxa_to_keep = 'Proteobacteria')
input_chloro = filter_taxa_from_input(input, taxa_to_keep = 'Chloroplast')
taxa_sum_chloro=summarize_taxonomy(input_chloro,level=4, report_higher_tax = FALSE)
taxa_sum_chloro
plot_ts_heatmap(taxa_sum_chloro, input$map_loaded,0.001, 'Treatment')
taxa_summary_by_sample_type(taxa_sum_chloro, input$map_loaded, 
                            type_header = 'Treatment', test_type = 'KW')
ChloroBTI=filter_data(input_chloro,filter_cat='Treatment',keep_vals= c('Control','Methoprene'))
taxa_sum_chloro=summarize_taxonomy(ChloroBTI,level=4, report_higher_tax = FALSE)
taxa_summary_by_sample_type(taxa_sum_chloro, input$map_loaded, 
                            type_header = 'Treatment', test_type = 'KW')

# Pairwise PERMANOVAS
dm = calc_dm(input$data_loaded)
ord = calc_ordination(dm, 'nmds')
calc_pairwise_permanovas(dm, metadata_map=input$map_loaded, 'Treatment')
calc_pairwise_permanovas(dm, metadata_map=input$map_loaded, 'Sampling_Date')

adonis(dm ~ Treatment*Sampling_Date, input$map_loaded)
adonis(dm ~ Sample_Type, input$map_loaded)
adonis(dm ~ Density, input$map_loaded)

#Faith D 
faith=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\faithDsum.txt",header=TRUE)
#faith
ggplot(faith, aes(x=Date, y=mean, color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem))+geom_point(size=3)+ ylab("Faith's D (SEM)")+geom_line(size=1)


cdata
ggplot(cdata, aes(x=Date, y=mean,color=Treatment, shape=Treatment, group=Treatment))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+geom_line(size=1)+geom_point(size=3)+ylab("Rel Abundance Photosynthetic Algae (%)")


mydata=read.table(file.choose(),header=TRUE)
(mydata)

myData.mean <- aggregate(mydata$Volume,
                         by = list(mydata$Treatment, mydata$Date),
                         FUN = 'mean')
head(myData.mean)
colnames(myData.mean) <- c("Treatment","Date", "Volume")
stress.aov <- with(myData.mean,aov(Volume ~ Treatment* Date))
stress.aov
summary(stress.aov)


##Dist methods
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

dist_methods["designdist"]
dist_methods=dist_methods[-which(dist_methods=="ANY")]
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(MozziePhyloseq, method=i)
  # Calculate ordination
  iMDS  <- ordinate(MozziePhyloseq, "PCoA", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(MozziePhyloseq, iMDS, color="Date", shape="Treatment")
  # Add title to each plot
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Sampling_Date, shape=Treatment))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("Distance metrics for MozziePhyloseq ")
p

print(plist[["e"]])
?distance

#PCoA from microbiome analyst
PCoAScore=read.table("C:\\Users\\Joe Receveur\\Documents\\Wallace\\Rdata\\pcoa_score.txt",header=TRUE)
ggplot(PCoAScore, aes(x=Axis.1, y=Axis.2))

sample_data(ControlMethpost)$DateTreat= paste0(sample_data(ControlMethpost)$Treatment,sample_data(ControlMethpost)$Sampling_Date) #creates DateTreat variable
ControlMethpost
head(sample_data(ControlMethpost))
DateTreatMerge=merge_samples(ControlMethpost, "DateTreat") #creates 4d object with DateTreat column
#sample_variables(DateTreatMerge)
sample_data(DateTreatMerge)
DateTreatMerge=transform_sample_counts(DateTreatMerge, function(x) 100 *x/sum(x)) #merging samples #(averaging)
(sample_data(DateTreatMerge))

sample_data(DateTreatMerge)$DateTreat=sample_names(DateTreatMerge)
sample_data(DateTreatMerge)$Treatment=c("Control", "Control", "Methoprene", "Methoprene")
sample_data(DateTreatMerge)$Date=c("7D Post", "14D Post", "7D Post", "14D Post" )
DateTreatChloropglom=tax_glom(DateTreatMerge, "Genus")

plot_bar(DateTreatMerge, "Treatment","Abundance","Order",facet_grid="Date~.")+xlab("Treatment")+ylab("Relative Microbial Abundance  (> 1%)") scale_fill_manual(values=cbPalette) 
TopNOTUs <- names(sort(taxa_sums(MozziePhyloseq), TRUE)[1:10])
ent10   <- prune_taxa(TopNOTUs, MozziePhyloseq)
plot_bar(SubsetControlMeth, "Treatment","Abundance", "Order", facet_grid="Date~.")




##SUBSET CONTROL METH
SubsetControlMeth=ControlMethpost
sample_data(SubsetControlMeth)
sample_data(SubsetControlMeth)$DateTreat= paste0(sample_data(SubsetControlMeth)$Treatment,sample_data(SubsetControlMeth)$Sampling_Date) #creates DateTreat variable
subset_samples(SubsetControlMeth, Date!= "4/18/2015")

DateTreatMerge=merge_samples(SubsetControlMeth, "DateTreat") #creates 4d object with DateTreat column
DateTreatMerge=transform_sample_counts(DateTreatMerge, function(x) 100 *x/sum(x)) #merging samples #(averaging)

DateTreatMerge
sample_data(DateTreatMerge)
(sample_data(DateTreatMerge))
ToKeep=c("Bacillus","Corynebacterium","Legionella","Lysinibacillus","Staphylococcus")
DateTreatMerge=subset_taxa(DateTreatMerge,Genus %in% ToKeep)
sample_data(DateTreatMerge)$DateTreat=sample_names(DateTreatMerge)
sample_data(DateTreatMerge)$Treatment=c("Control", "Control", "Methoprene", "Methoprene")
sample_data(DateTreatMerge)$Date=c("7D Post", "14D Post", "7D Post", "14D Post" )
sample_data(DateTreatMerge)$Date = factor(sample_data(DateTreatMerge)$Date, levels = c("7D Post", "14D Post")) #fixes x-axis labels
DateTreatMergeglom=tax_glom(DateTreatMerge, "Genus")

p=plot_bar(DateTreatMergeglom, "Treatment","Abundance","Genus",facet_grid="~Date")+xlab("")+ylab("Relative Bacterial Abundance  (> 1%)")+ scale_fill_manual(values=cbPalette)
p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))


#Late Instar by treatment
sample_data(MozziePhyloseq)
DateTreatMerge=merge_samples(MozziePhyloseq, "InstarTreat") #creates 4d object with DateTreat column
#sample_variables(DateTreatMerge)
sample_data(DateTreatMerge)
DateTreatMerge=transform_sample_counts(DateTreatMerge, function(x) 100 *x/sum(x)) #merging samples #(averaging)
(sample_data(DateTreatMerge))
sample_data(DateTreatMerge)$InstarTreat=sample_names(DateTreatMerge)
sample_data(DateTreatMerge)$Treatment=c("BTI","Control","Methoprene","BTI","Control","Methoprene","Control","Methoprene")
sample_data(DateTreatMerge)$Instar=c("2nd","2nd","2nd","3rd","3rd","3rd","4th","4th" )
DateTreatMerge
sample_data(DateTreatMerge)
DateTreatMergeglom=tax_glom(DateTreatMerge, "Genus")

sample_variables(DateTreatMergeglom)
#DateTreatMergeglom=transform_sample_counts(DateTreatMergeglom, function(x) 100 *x/sum(x)) #merging samples #(averaging)
DateTreatMergeglom=subset_samples(DateTreatMergeglom, Treatment!="BTI")
DateTreatMergeglom=subset_samples(DateTreatMergeglom, Date!="1D Post")
Genusplot=plot_bar(DateTreatMerge, "Treatment","Abundance", "Genus")+xlab("")+ylab("Relative Bacterial Abundance (> 1%)") +facet_grid(Instar~.)+ scale_fill_manual(values=cbPalette)
Genusplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


glm()