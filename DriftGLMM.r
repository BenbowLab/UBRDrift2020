#Sturgeon Drift GLMM TESTING
#install.packages('TMB',type='source')
library(glmmTMB)
library(bbmle)
# 6. Fit Poisson GLMMs
ShannonSubset<-ShannonRichness
head(ShannonSubset)
# With count data, a good place to start can be the Poisson distribution. It's a very
# standard distribution for count data that can take integers from 0 to + infinity.
# Because our observations are grouped together within watersheds, streams, and years,
# we should use random effects to represent that variation. Let's also say that during the
# course of the study you noticed that the riparian habitat in different streams was of
# different quality. But you weren't able to assess riparian habitat quality. You just
# measured amount of riparian habitat. Your observation that riparian habitat quality 
# varies among streams has led you to think that the effect of the amount of riparian
# habitat on fish numbers might vary among streams (because of the unmeasured variation in
# habitat quality).

# Build a set of models with all possible combinations of fixed effects with your two
# predictors: riparianc and invasivec. The models should all have the same random effects
# structure: random intercepts for your 3 grouping factors, and add a random slope to
# represent variation in the effect of riparian habitat from stream to stream.
# Use the glmer function from lme4, and remember to tell glmer you want a Poisson dist.

#Centering Flow
hist(ShannonSubset$percillum)
ShannonSubset$QCentered = ShannonSubset$Q - mean(ShannonSubset$Q)
head(ShannonSubset)

ShannonSubset$MoonPhase = factor(ShannonSubset$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

ShannonSubset<-subset(ShannonSubset, percillum!= "NA")

ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
#ShannonSubset<-subset(ShannonSubset, Q!= "NA")
ShannonSubset<-subset(ShannonSubset, InvertByFlow!= "NA")
ShannonSubset$InvertByFlowByMeter<- as.numeric(ShannonSubset$InvertByFlowByMeter)
ShannonSubset$InvertByFlow<- as.numeric(ShannonSubset$InvertByFlow)
ShannonSubset$BiomassByFlow<- as.numeric(ShannonSubset$BiomassByFlow)
ShannonSubset$InvertByFlowCount<- as.numeric(ShannonSubset$InvertByFlowCount)
ShannonSubset#ShannonSubset<-subset(ShannonSubset, Ninverts100 < 4000)

length(ShannonSubset$SampleID)
p0 = glmer(Ninverts100~1+(1|Year)+(1|DPFS),data=ShannonSubset,family = poisson)

p1 = glmer(Ninverts100~percillum+(1|Year),data=ShannonSubset,family = poisson)
p2 = glmer(Ninverts100~percillum+Q+(1|Year)+(1|DPFS),data=ShannonSubset,family = poisson)
summary(p2)
plot(ShannonSubset$percillum~ShannonSubset$Ninverts100)
# The first thing to note (if you ran the same models I did) is the warning:
# "boundary (singular) fit: see ?isSingular"
# That warning is telling us that one of the parameters hsa been estimated to be on the
# boundary of what it is possible for it to be. In most cases, this means that one of the
# random effects variances we tried to estimate was estimated to be close to zero. Zero is
# the lower boundary for a variance parameter, because variances can't be negative.

# Use AIC to find the best model and examine that one to see what's going on:
AICctab(p0, p1, p2, weights=TRUE)

# Now examine the model with the lowest AIC
summary(p2)
#Check deviance vs df. resid
#deviance much higher than df.resid == overdispersed
#Try neg binomial

# 7. Fit negative binomial models in maximum likelihood

# One very useful way to model count data with more variance than the Poisson is with the
# negative binomial distribution. Elise introduced the negative binomial near the end of
# her section. The best way for us to think about the negative binomial is as a Poisson
# distribution with an extra parameter to model the extra variance in the counts. That
# extra parameter goes by a few different names, including theta or k or size (in R). I'll
# call it theta. When theta is small, it (counterintuitively) means that there is a lot of
# variance. When theta is large, there is less variance. When theta is greater than ~10,
# the Poisson and negative binomial look very similar. With ecological count data, theta
# is often around 1 or < 1, which means that the variance is much higher than the mean,
# and the data have much more variance than is assumed by the Poisson.

# Note: the negative binomial isn't the only way to model count data with extra variance, 
# but many find it to do a good job with ecological data, and it's generally easier than 
# the other options.

# Unfortunately, the negative binomial is slightly harder to fit than the Poisson because
# of that theta parameter. This means our friend glmer() cannot fit models with a negative
# binomial distribution. There is a glmer.nb() for negative binomials, but it is still 
# experimental as far as I know. And it fits the model in a stepwise fashion that means
# you can't really make inference on the theta parameter. Instead, let's use the glmmTMB()
# function from the glmmTMB package. The syntax is exactly the same as the lme4 syntax
# (some of the same people made both packages). But glmmTMB() can fit much more
# complicated models, including negative binomial models (and beta, gamma, etc., and
# zero-inflated models, hurdle models, etc.).

# OK, now use glmmTMB to fit the same set of models above (without the random effect 
# estimated to be zero) but use a negative binomial distribution to account for the high 
# variance in our fish count data. The syntax is exactly the same as our Poisson glmers
# except we need to specify family='nbinom2'. (nbinom2 because there are different ways
# to set up the negative binomial, but the standard way in ecology is the second way.)
head(ShannonSubset)
n0 = glmmTMB(InvertByFlowCount~1+(1|Year),data=ShannonSubset,family = 'nbinom2')

n1 = glmmTMB(InvertByFlowCount~percillum+(1|Year),data=ShannonSubset,family = 'nbinom2')
n2 = glmmTMB(InvertByFlowCount~percillum+Temp+(1|Year),data=ShannonSubset,family = 'nbinom2')
n3 = glmmTMB(Ninverts100~Q+Temp+(1|Year),data=ShannonSubset,family = 'nbinom2')
n4 = glmmTMB(Ninverts100~percillum+Temp+(1|Year),data=ShannonSubset,family = 'nbinom2')
summary(n1)
# Compare models with AICc
AICctab(n0,n1,n2, weights=TRUE)

# Examine best model
summary(n2)

# QUESTION What is the best model telling you about the population-level (fixed) effects
# of riparian habitat amount and number of invasive species?
# Strong positive effect of riparian habitat on fish counts. Negative (but weaker) effect
# of invasive predator species richness on fish counts. Perhaps a very week negative
# interaction between riparian habitat and invasive richness, such that increasing
# riparian habitat slightly lessens the negative effects of invasive predators on fish
# counts. But 95% CI of the estimate for the interaction overlaps zero (-0.032, 0.0003).

# QUESTION What are the estimated random effects variances telling you about how things
# vary across groups?
# Large variance from stream to stream in both fish counts at the mean of the predictors
# and in the relationship between riparian habitat amount and fish counts. Also
# considerable variance from year to year.

# Compare the best negative binomial model to the best Poisson model using AICc
AICctab(n1, p2, weights=TRUE)
# The negative binomial model is doing a MUCH better job job. Massive difference in AICc.
summary(n1)
plot(resid(n1) ~ percillum, data=ShannonSubset)
lines(lowess(resid(n1) ~ ShannonSubset$percillum), col=2)
exp(fixef(n1[1]))

plot(exp(predict(n1))~ShannonSubset$Ninverts100)+abline(0,1)

hist(ShannonSubset$percillum)
newillum = seq(0,1,length=500)

# Second, on the transformed scale (if you transformed), plot the data
plot(Ninverts100 ~ percillum, data=ShannonSubset)

# Fourth, add line for predicted mean relationship across all genotypes
exp(-0.5827)
predy = 7.5567 + exp(-0.5827) * newillum
lines(predy ~ newillum, lwd=3, col=2)

# Now do the same thing again but plotting on the arithmetic scale, back-transforming the
# predictions, if needed.  There are two blanks to fill.
plot(defense ~ cats, data=d)
for(i in 1:120) {
  newy = coef(m1)$genoID[i,1] + coef(m1)$genoID[i,2]*newcats
  lines(exp(newy) ~ newcats, col=grey(0.5, alpha=0.75), lwd=0.5)
}
lines(exp(predy) ~ newcats, lwd=3, col=2)

#ShannonSubset<-subset(ShannonSubset, Ninverts100 < 4000)

GLMILLum<-glm(Ninverts100~percillum,data=ShannonSubset)
summary(GLMILLum)
plot(resid(GLMILLum) ~ percillum, data=ShannonSubset)
lines(lowess(resid(GLMILLum) ~ ShannonSubset$percillum), col=2)

confint(GLMILLum)
hist(resid(GLMILLum))
newcats = seq(0,1, length=1000)

l1.pred = predict(GLMILLum, newdata=data.frame(percillum=newcats), se.fit=TRUE)
l1.pred.mean = l1.pred$fit
l1.pred.975 = l1.pred$fit + 1.96*l1.pred$se.fit
l1.pred.025 = l1.pred$fit - 1.96*l1.pred$se.fit

plot(InvertByFlow ~ percillum, data=ShannonSubset)
lines(l1.pred.mean ~ newcats, col=2)
lines(l1.pred.975 ~ newcats, col=2, lty=2)
lines(l1.pred.025 ~ newcats, col=2, lty=2)



plot(Biomass_g100 ~ Q, data=ShannonSubset)
ggplot(ShannonSubset,aes(x=LunarDay,y=InvertByFlow))+geom_smooth()
ggplot(ShannonSubset,aes(x=percillum,y=Biomass_g100))+geom_smooth()+geom_point()
GlmInvertByFlow<-glm(Ninverts100~percillum,data=ShannonSubset)
summary(GlmInvertByFlow)

variable = seq(min(ShannonSubset$percillum),max(ShannonSubset$percillum), length=1000)

l1.pred = predict(GlmInvertByFlow, newdata=data.frame(percillum=variable), se.fit=TRUE)
l1.pred.mean = l1.pred$fit
l1.pred.975 = l1.pred$fit + 1.96*l1.pred$se.fit
l1.pred.025 = l1.pred$fit - 1.96*l1.pred$se.fit

plot(Ninverts100 ~ percillum, data=ShannonSubset,ylab="Nightly Invertebrate Abundance",xlab= "Percent illumination")
#text(x = 20,y=6000,label="Abundance= -208.66 * Temp+5454.9") #add text to plot
lines(l1.pred.mean ~ variable, col=2)
lines(l1.pred.975 ~ variable, col=2, lty=2)
lines(l1.pred.025 ~ variable, col=2, lty=2)



plot(Temp~DPFS,data=ShannonRichness)




plot(resid(GlmInvertByFlow) ~ percillum, data=ShannonSubset)
lines(lowess(resid(GlmInvertByFlow) ~ ShannonSubset$percillum), col=2)
hist(resid(GlmInvertByFlow))



ggplot(ShannonSubset,aes(x=DPFS,y=Temp))+geom_smooth()+geom_point()+xlab("Days Post First Spawning")+ylab("Temperature (C)")
ggplot(ShannonSubset,aes(x=DPFS,y=Q))+geom_smooth()+geom_point()+xlab("Days Post First Spawning")+ylab("Discharge (m3/s)")

ggplot(ShannonSubset,aes(x=DPFS,y=Ninverts100))+geom_smooth()

cor(ShannonSubset$percillum,ShannonSubset$DPFS)
cor(ShannonSubset$LunarDay,ShannonSubset$DPFS)

mean(ShannonSubset$percillum)
hist(ShannonSubset$percillum)


ShannonSubset<-ShannonRichness
ShannonSubset$MoonPhase = factor(ShannonSubset$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

ShannonSubset<-subset(ShannonSubset, percillum!= "NA")

#ShannonSubset<-subset(ShannonSubset, Temp!= "NA")
ShannonSubset<-subset(ShannonSubset, Q!= "NA")
ShannonSubset<-subset(ShannonSubset, InvertByFlow!= "NA")
ShannonSubset$InvertByFlowByMeter<- as.numeric(ShannonSubset$InvertByFlowByMeter)
ShannonSubset$InvertByFlow<- as.numeric(ShannonSubset$InvertByFlow)



ggplot(ShannonSubset,aes(x=MoonPhase,y=InvertByFlow))+geom_boxplot()+geom_point()

Trtdata <- ddply(ShannonSubset, c("MoonPhase"), summarise,
                 N    = length(InvertByFlow),
                 meanInvert = mean(InvertByFlow),
                 sd   = sd(InvertByFlow),
                 se   = sd / sqrt(N)
)
Trtdata
ggplot(Trtdata, aes(x=MoonPhase,y=meanInvert))+geom_bar(aes(),colour="black", stat="identity")+xlab("Moon Phase")+ylab("Inverts by flow")+
  geom_errorbar(aes(ymin=meanInvert-se,ymax=meanInvert+se))


ShannonRichness$MoonPhase = factor(ShannonRichness$MoonPhase, levels = c("New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"))

ggplot(ShannonRichness,aes(x=MoonPhase, y=percillum))+geom_boxplot()+geom_point()
ggplot(ShannonSubset,aes(x=percillum, y=InvertByFlow))+geom_point()+geom_smooth()


#ShannonSubset<-subset(ShannonSubset, Q!= "NA")
ShannonRichness
Trtdata <- ddply(ShannonRichness, c("MoonPhase"), summarise,
                 N    = length(Q),
                 meanInvert = mean(Q),
                 sd   = sd(Q),
                 se   = sd / sqrt(N)
)
Trtdata
ggplot(Trtdata, aes(x=Q,y=meanInvert))+geom_bar(aes(),colour="black", stat="identity")+xlab("Moon Phase")+ylab("Inverts by flow")+
  geom_errorbar(aes(ymin=meanInvert-se,ymax=meanInvert+se))
