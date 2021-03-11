#Ch1-Analysis-final

#Lindsey Kemmerling
#Nov 12, 2020
#updated January 29, 2021.

#Analysis for pollinator data across all treatments in the GLBRC in 2018.

rm(list=ls())

setwd("/Users/lindseykemmerling/Dropbox/Research2018/GLBRC/data entry/biofuels")

#####Load needed libraries#####
library(lme4) #For constructing the models.
library(car) #For testing significance of predictors.
library(reshape2)
library(ggplot2)
library(Rmisc)
library(multcomp)
library(glmmADMB)
library(iNEXT)
library(outliers)

poll<-read.csv("data-entry_pollinators_GLBRC-group_analysis.csv", header=TRUE)
site<-read.csv("GLBRC_site_biofuel.csv", header=TRUE)
site.poll<-read.csv("GLBRC_site_biofuel_poll.csv", header=TRUE)
glbrc<-read.csv("GLBRC_yield_2018.csv", header=TRUE)
corn<-read.csv("GLBRC_yield_2018_plus-corn.csv", header=TRUE)
flowers<-read.csv("data-entry_flowers_GLBRC_20Dec19.csv", header=TRUE)

site$block=as.factor(site$block)
site$sampling_round=as.factor(site$sampling_round)

poll$treatment[poll$treatment=="native grass"] <- "native grasses"
site$treatment[site$treatment=="native grass"] <- "native grasses"
site.poll$treatment[site.poll$treatment=="native grass"] <- "native grasses"
glbrc$treatment[glbrc$treatment=="native grass"] <- "native grasses"
flowers$treatment[flowers$treatment=="native grass"] <- "native grasses"


##################################################
###### FLORAL RESOURCES VS TREATMENT #############
##################################################

#reshape data
floral <- flowers[which(!flowers$treatment %in% c("poplar")),]
floral$floral_resources=as.numeric(floral$floral_resources)
floral <- aggregate(floral_resources ~ block + treatment + sampling_round, data=floral, FUN=sum)

#add zeroes in floral where currently NA
floral.site<-merge(floral,site,all=TRUE) #merge with site info
floral.site$floral_resources[is.na(floral.site$floral_resources)]=0 #create 0s where currently NA

floral.site$block=as.factor(floral.site$block)
floral.site$sampling_round=as.factor(floral.site$sampling_round)
floral.site$floral_resources=as.numeric(floral.site$floral_resources)
floral.site$treatment=as.factor(floral.site$treatment)

hist(floral.site$floral_resources)

#test for dispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#find outlier/z score
z.scores.floral <- scores(floral.site$floral_resources, type = c("z"))
floral.outlier <- which(z.scores.floral > 3)
z.scores.floral
floral.outlier
#point 12, 75 is an outlier

################
model.flowers.negbi<-glmmadmb(floral_resources ~ treatment + sampling_round + (1|block), 
                               data=floral.site, family="nbinom")



overdisp_fun(model.flowers.negbi)
#looks good

Anova(model.flowers.negbi, type="III")
summary(model.flowers.negbi)

contrast.flowers.negbi <-  glht(model.flowers.negbi, linfct=mcp(treatment="Tukey"))

contrast.flowers.negbi
cld(contrast.flowers.negbi, Letters=Letters)

#plot model
abundance.floral.sum <- summarySE(floral, measurevar="floral_resources", groupvars=c("treatment"))
abundance.floral.sum

abundance.floral.plot = ggplot(abundance.floral.sum, aes(x=treatment, y=floral_resources), geom="bar", stat="identity",position="dodge")+
  geom_bar(stat="identity", width=0.7, fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73"), colour = "black") +
  ylim(0,24000)+
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Biofuel crop")+
  ylab("Flower abundance")+
  ggtitle("Flower abundance vs. biofuel crop")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(plot.title = element_blank())

abundance.floral.plot=abundance.floral.plot + 
  geom_errorbar(aes(ymax=floral_resources+se,ymin=floral_resources-se), position=position_dodge(0.6), width=.2, data=abundance.floral.sum)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

abundance.floral.plot


ggsave("figure2a.jpg", dpi=300, height=5, width=6, units="in")


#####plant abundance############
#########Not used in paper!!####
################################
#reshape data
plants <- flowers[which(!flowers$treatment %in% c("poplar")),]
plants$number_plants=as.numeric(plants$number_plants)
plants <- aggregate(number_plants ~ block + treatment + sampling_round, data=plants, FUN=sum)

#add zeroes in plants where currently NA
plants.site<-merge(plants,site,all=TRUE) #merge with site info
plants.site$number_plants[is.na(plants.site$number_plants)]=0 #create 0s where currently NA

plants.site$block=as.factor(plants.site$block)
plants.site$sampling_round=as.factor(plants.site$sampling_round)
plants.site$number_plants=as.numeric(plants.site$number_plants)
plants.site$treatment=as.factor(plants.site$treatment)

hist(plants.site$number_plants)
#negative binomial, use glmmadmb for model

#find outlier/z score
z.scores.plants <- scores(plants.site$number_plants, type = c("z"))
plants.outlier <- which(z.scores.plants > 3)
z.scores.plants
plants.outlier

#model.plants<-glmmadmb(number_plants ~ treatment + sampling_round + (1|block),
#                   data=plants.site, family = "nbinom")

#overdisp_fun(model.plants)
#looks good

#Anova(model.plants, type="III")
#summary(model.plants)

#contrast.plants <-  glht(model.plants, linfct=mcp(treatment="Tukey"))

#contrast.plants
#cld(contrast.plants, Letters=Letters)

#plot model
#plant.plants.sum <- summarySE(plants, measurevar="number_plants", groupvars=c("treatment"))
#plant.plants.sum

#plant.plants.plot = ggplot(plant.plants.sum, aes(x=treatment, y=number_plants), geom="bar", stat="identity",position="dodge")+
#  geom_bar(stat="identity", width=0.7, fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73"), colour = "black") +
#  ylim(0,675)+
#  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
#  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
#  xlab("Biofuel crop")+
#  ylab("Plant abundance")+
#  ggtitle("Plant abundance vs. biofuel crop")+
#  theme(axis.title.x = element_text(size=16, vjust=1)) + theme(axis.title.y = element_text(margin=margin(r=10), size=16)) +
#  theme(plot.title = element_text(margin=margin(b=15), size=16))+theme(plot.title = element_text(hjust = 0.5))

#plant.plants.plot=plant.plants.plot + 
#  geom_errorbar(aes(ymax=number_plants+se,ymin=number_plants-se), position=position_dodge(0.6), width=.2, data=plant.plants.sum)+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#plant.plants.plot
####################





###################################
##RAREFY PLANT DATA BY ABUNDANCE###
###################################

#create plant/site datafram
plants <- flowers[which(!flowers$treatment %in% c("poplar")),]
plant.rich <- aggregate(flower ~ block + treatment + sampling_round + plot, data=plants, FUN=function(x) length(unique(x)))
plant.rich <- rename(plant.rich, c("flower"="plant.richness"))

#add zeroes in poll where currently NA
plant.rich.site<-merge(plant.rich,site,all=TRUE) #merge with site info
plant.rich.site$plant.richness[is.na(plant.rich.site$plant.richness)]=0 #create 0s where currently NA

plant.rich.site$block=as.factor(plant.rich.site$block)
plant.rich.site$sampling_round=as.factor(plant.rich.site$sampling_round)

#reshape data
plants <- flowers[which(!flowers$treatment %in% c("poplar")),]
plant.rich <- aggregate(flower ~ block + treatment + sampling_round + plot, data=plants, FUN=function(x) length(unique(x)))
plant.rich <- rename(plant.rich, c("flower"="plant.richness"))

Castdata <- dcast(plants, flower ~ plot, value.var = "number_plants", 
                  fun.aggregate = sum)
Castdata <- Castdata[-c(1), ] #remove first row of just zeroes, not sure why that is there
Castdata2 <- Castdata[,-1] #to change first column to rownames, instead of a value, remove the row then re-add it as a rowname 
rownames(Castdata2) <- Castdata[,1]

####################
#Rarefy by abundance
####################

output.abundance <- iNEXT(Castdata2, q=0, datatype="abundance")
output.abundance
list2env(output.abundance,.GlobalEnv) #convert iNext output into independent dataframes
AsyEst <- AsyEst[which(AsyEst$Diversity %in% c("Species richness")),] #use only species richness estimates
plant.rarefied.est <- AsyEst[ -c(2, 3, 6, 7) ] #remove rows to merge with plant data
plant.rarefied.est <- rename(plant.rarefied.est, c("Site"="plot"))#make column names match
plant.rarefied<-merge(plant.rich.site,plant.rarefied.est,all=TRUE) #merge with plant data

hist(plant.rarefied$Estimator)
#hist(flowers.pollinators$poll.abundance)
#normal dist, use lmer for model

#model with rarefied plant richness
#model
model.plant.rich.rare<-lmer(Estimator ~ treatment + (1|block), 
                            data=plant.rarefied)

Anova(model.plant.rich.rare, type="III")
summary(model.plant.rich.rare)

contrast.plant.rich.rare <-  glht(model.plant.rich.rare, linfct=mcp(treatment="Tukey"))

contrast.plant.rich.rare
cld(contrast.plant.rich.rare, Letters=Letters)

#plot model
rich.plant.rare.sum <- summarySE(plant.rarefied, measurevar="Estimator", groupvars=c("treatment"))
rich.plant.rare.sum

rich.plant.rare.plot = ggplot(rich.plant.rare.sum, aes(x=treatment, y=Estimator), geom="bar", stat="identity",position="dodge")+
  geom_bar(stat="identity", width=0.7, fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73"), colour = "black") +
  ylim(0,20)+
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Biofuel crop")+
  ylab("Flowering plant richness")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())

rich.plant.rare.plot=rich.plant.rare.plot + 
  geom_errorbar(aes(ymax=Estimator+se,ymin=Estimator-se), position=position_dodge(0.6), width=.2, data=rich.plant.rare.sum)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

rich.plant.rare.plot

ggsave("figure2b.jpg", dpi=300, height=5, width=6, units="in")



#########################################
##POLLINATOR ABUNDANCE VS TREATMENT #####
#########################################

#reshape pollinator data
pollinators <- poll[which(!poll$treatment %in% c("poplar")),]
pollinators <- pollinators[which(!pollinators$group_analysis %in% c("0", "no")),]
pollinators <- aggregate( unique_ID ~ block + treatment + sampling_round + plot, data=pollinators, FUN=function(x) length(unique(x)))
pollinators <- rename(pollinators, c("unique_ID"="poll.abundance"))

#add zeroes in poll where currently NA
pollinators.site<-merge(pollinators,site,all=TRUE) #merge with site info
pollinators.site$poll.abundance[is.na(pollinators.site$poll.abundance)]=0 #create 0s where currently NA

pollinators.site$block=as.factor(pollinators.site$block)
pollinators.site$sampling_round=as.factor(pollinators.site$sampling_round)

#merge pollinator and rarefied floral richness data
poll.flowers <- merge(plant.rarefied, pollinators.site, all=TRUE)
poll.flowers <- rename(poll.flowers, c("Estimator"="flow.rare.rich"))

#merge in floral abundance data
flowers.pollinators <- merge(poll.flowers, plants.site, all=TRUE)
flowers.pollinators$treatment=as.factor(flowers.pollinators$treatment)

hist(flowers.pollinators$poll.abundance)
#negative binomial, use glmmadmb for model

#remove outlier
#flowers.pollinators <- flowers.pollinators[-75,]

#model pollinator abundance by treatment
model.poll.flow.abun<-glmmadmb(poll.abundance ~ sampling_round + treatment + (1|block), 
                               data=flowers.pollinators, family="nbinom")

overdisp_fun(model.poll.flow.abun)
#looks good

Anova(model.poll.flow.abun, type="III")
summary(model.poll.flow.abun)

contrast.poll.flow.abun <-  glht(model.poll.flow.abun, linfct=mcp(treatment="Tukey"))

contrast.poll.flow.abun
cld(contrast.poll.flow.abun, Letters=Letters)

#plot treatment by poll.abun
treat.poll.sum <- summarySE(flowers.pollinators, measurevar="poll.abundance", groupvars=c("treatment"))
treat.poll.sum

treat.poll.plot = ggplot(treat.poll.sum, aes(x=treatment, y=poll.abundance), geom="bar", stat="identity",position="dodge")+
  geom_bar(stat="identity", width=0.7, fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73"), colour="black") +
  ylim(0,85)+
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Biofuel crop")+
  ylab("Pollinator abundance")+
  ggtitle("Pollinator abundance vs. biofuel crop")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())

treat.poll.plot=treat.poll.plot + 
  geom_errorbar(aes(ymax=poll.abundance+se,ymin=poll.abundance-se), position=position_dodge(0.6), width=.2, data=treat.poll.sum)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

treat.poll.plot

ggsave("figure2c.jpg", dpi=300, height=5, width=6, units="in")



#######################################
###POLLINATOR DIVERSITY BY TREATMENT###
#######################################

library(MuMIn)

#reshape data
poll.rich <- poll[which(!poll$treatment %in% c("poplar")),]
poll.rich <- aggregate(group_analysis ~ block + treatment + sampling_round + plot, data=poll.rich, FUN=function(x) length(unique(x)))
poll.rich <- rename(poll.rich, c("group_analysis"="richness"))

#add zeroes in poll where currently NA
rich.site<-merge(poll.rich,site,all=TRUE) #merge with site info
rich.site$richness[is.na(rich.site$richness)]=0 #create 0s where currently NA

rich.site$block=as.factor(rich.site$block)
rich.site$treatment=as.factor(rich.site$treatment)
rich.site$sampling_round=as.factor(rich.site$sampling_round)

hist(rich.site$richness, breaks = seq(0,12,1))

#model
model.rich<-lmer(richness ~ treatment + sampling_round + (1|block), 
                     data=rich.site)
AICc(model.rich)

overdisp_fun(model.rich)

Anova(model.rich, type="III")
summary(model.rich)

contrast.rich <-  glht(model.rich, linfct=mcp(treatment="Tukey"))

contrast.rich
cld(contrast.rich, Letters=Letters)



#plot model treatment by richness
rich.poll.sum <- summarySE(poll.rich, measurevar="richness", groupvars=c("treatment"))
rich.poll.sum

rich.poll.plot = ggplot(rich.poll.sum, aes(x=treatment, y=richness), geom="bar", stat="identity",position="dodge")+
  geom_bar(stat="identity", width=0.7, fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73"), colour="black") +
  ylim(0,8)+
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Biofuel crop")+
  ylab("Pollinator group richness")+
  ggtitle("Pollinator group richness vs. biofuel crop")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())

rich.poll.plot=rich.poll.plot + 
  geom_errorbar(aes(ymax=richness+se,ymin=richness-se), position=position_dodge(0.6), width=.2, data=rich.poll.sum)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

rich.poll.plot


ggsave("figure2d.jpg", dpi=300, height=5, width=6, units="in")


#test correlation
#plot <- ggplot(flowers.pollinators, aes(flow.rare.rich, number_plants)) + 
 # geom_point() +
 # theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
 # theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
 # geom_smooth(method='lm', colour="black")

#plot



#########################################
#########################################
##POLLINATORS EXPLAINED BY FLOWERS ######
#########################################
#########################################

#reshape pollinator data
pollinators <- poll[which(!poll$treatment %in% c("poplar")),]
pollinators <- pollinators[which(!pollinators$group_analysis %in% c("0", "no")),]
pollinators <- aggregate( unique_ID ~ block + treatment + sampling_round + plot, data=pollinators, FUN=function(x) length(unique(x)))
pollinators <- rename(pollinators, c("unique_ID"="poll.abundance"))

#add zeroes in poll where currently NA
pollinators.site<-merge(pollinators,site,all=TRUE) #merge with site info
pollinators.site$poll.abundance[is.na(pollinators.site$poll.abundance)]=0 #create 0s where currently NA

pollinators.site$block=as.factor(pollinators.site$block)
pollinators.site$sampling_round=as.factor(pollinators.site$sampling_round)

#merge pollinator and rarefied floral richness data
poll.flowers <- merge(plant.rarefied, pollinators.site, all=TRUE)
poll.flowers <- rename(poll.flowers, c("Estimator"="flow.rare.rich"))
poll.flowers <- merge (rich.site, poll.flowers, all=TRUE)
poll.flowers <- rename(poll.flowers, c("richness"="poll.div"))

#merge in floral abundance data
flowers.pollinators <- merge(poll.flowers, plants.site, all=TRUE)

#find outlier/z score
z.scores.poll <- scores(flowers.pollinators$poll.abundance, type = c("z"))
poll.outlier <- which(z.scores.poll > 3)
#point 75 in an outlier

#remove outlier
flowers.pollinators <- flowers.pollinators[-75,]

#model pollinator abundance by flower abundance
poll.flow.abun<-lm(poll.abundance ~ number_plants, 
                           data=flowers.pollinators)
summary(poll.flow.abun)

plot <- ggplot(flowers.pollinators, aes(number_plants, poll.abundance)) + 
  geom_point() +
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Flowering plant abundance")+
  ylab("Pollinator abundance")+
  ggtitle("Pollinator abundance vs. flowering plant abundance")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method='lm', colour="black")

plot

#resid.plot <- resid(poll.flow.abun)
#plot(flowers.pollinators$floral_resources, resid.plot) 
#abline(0, 0)

ggsave("figure3a.jpg", dpi=300, height=5, width=6, units="in")


###############
poll.flow.div<-lm(poll.div ~ flow.rare.rich, data = flowers.pollinators)
summary(poll.flow.div)

plot <- ggplot(flowers.pollinators, aes(flow.rare.rich, poll.div)) + 
  geom_point() +
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Flower richness")+
  ylab("Pollinator group richness")+
  ggtitle("Pollinator group richness vs. flower richness")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method='lm', colour="black")

plot

ggsave("figure3d.jpg", dpi=300, height=5, width=6, units="in")


#model pollinator abundance by flower abundance
poll.flow.abun<-lm(poll.div ~ number_plants, 
                   data=flowers.pollinators)
summary(poll.flow.abun)

plot <- ggplot(flowers.pollinators, aes(number_plants, poll.div)) + 
  geom_point() +
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Flowering plant abundance")+
  ylab("Pollinator group richness")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method='lm', colour="black")

plot

ggsave("figure3c.jpg", dpi=300, height=5, width=6, units="in")

###############
poll.flow.div<-lm(poll.abundance ~ flow.rare.rich, data = flowers.pollinators)
summary(poll.flow.div)

plot <- ggplot(flowers.pollinators, aes(flow.rare.rich, poll.abundance)) + 
  geom_point() +
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Flower richness")+
  ylab("Pollinator abundance")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method='lm', colour="black")

plot

ggsave("figure3b.jpg", dpi=300, height=5, width=6, units="in")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#########################################
#########################################
##POLLINATORS EXPLAINED BY FLOWERS ######
#########################################
#########################################








###############################
######YIELD VS TREATMENT#######
###############################

#Reshape data
yield <- glbrc[which(glbrc$main_or_microplot %in% c("main")),]
#yield <- aggregate(dry_matter_yield_Mg_ha ~ replicate + treatment, data=glbrc, FUN=sum)
yield <- rename(yield, c("dry_matter_yield_Mg_ha"="Total_yield"))


hist(yield$Total_yield)
#normal dist, lmer

model.yield<-lmer(Total_yield ~ treatment + (1|replicate), 
                  data=yield)

Anova(model.yield, type="III")
summary(model.yield)

contrast.yield <-  glht(model.yield, linfct=mcp(treatment="Tukey"))

contrast.yield
cld(contrast.yield, Letters=Letters)


#plot model
yield.sum <- summarySE(yield, measurevar="Total_yield", groupvars=c("treatment"))
yield.sum

yield.plot = ggplot(yield.sum, aes(x=treatment, y=Total_yield), geom="bar", stat="identity",position="dodge")+
  geom_bar(stat="identity", width=0.7, colour="black", fill= c("native grass" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73")) +
  ylim(0,8)+
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  xlab("Biofuel crop")+
  ylab("Crop yield (Mg/ha)")+
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.title.y = element_text(size=18)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(plot.title = element_blank())

yield.plot=yield.plot + 
  geom_errorbar(aes(ymax=Total_yield+se,ymin=Total_yield-se), position=position_dodge(0.6), width=.2, data=yield.sum)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

yield.plot

ggsave("yield.plot.300.jpg", dpi=300, height=6, width=7, units="in")

###########
#average corn yield in 2018 at KBS
###########
#Reshape data
corn <- corn[which(corn$main_or_microplot %in% c("main")),]
corn <- corn[which(corn$treatment %in% c("Continuous Corn")),]
corn <- aggregate(dry_matter_yield_Mg_ha ~ replicate + treatment, data=corn, FUN=sum)
corn <- rename(corn, c("dry_matter_yield_Mg_ha"="Total_yield"))

avg.corn.yield <- mean(corn[,3])

##########################################################
######YIELD VS POLLINATOR DIVERSITY PARETO FRONTIER#######
##########################################################

#merge yield and pollinator richness data
yield.2 <- glbrc[which(glbrc$main_or_microplot %in% c("main")),]
yield.2 <- aggregate(dry_matter_yield_Mg_ha ~ replicate + treatment + plot, data=yield.2, FUN=sum)
yield.2 <- rename(yield.2, c("dry_matter_yield_Mg_ha"="Total_yield"))
yield.3 <- aggregate(Total_yield ~ treatment, data=yield.2, FUN=mean)

poll.rich.2 <- poll[which(!poll$treatment %in% c("poplar")),]
poll.rich.2 <- aggregate(group_analysis ~ block + treatment + plot, data=poll.rich.2, FUN=function(x) length(unique(x)))
poll.rich.2 <- rename(poll.rich.2, c("group_analysis"="richness"))
poll.rich.3 <- aggregate(richness ~ treatment, data=poll.rich.2, FUN=mean)

yield.rich<-merge(poll.rich.2,yield.2,all=TRUE) #merge yield and richness data
yield.rich.means <- merge(poll.rich.3, yield.3, all=TRUE)

yield.rich.model<-lm(Total_yield ~ richness, data = yield.rich)
summary(yield.rich.model)


plot <- ggplot(yield.rich, aes(richness, Total_yield, color=treatment, group=1)) + 
  geom_point(alpha = .8) +
  scale_color_manual(values = c("native grasses" = "#E69F00", "prairie" = "#D55E00", "successional" = "#0072B2", "switchgrass" = "#009E73")) +
  geom_point(data = yield.rich.means, size=4) +
  theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0)))+
  theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))+
  theme(legend.key=element_blank())+
  theme(legend.title = element_blank())+
  xlab("Pollinator group richness")+
  ylab("Crop yield")+
  ggtitle("Optimization of yield and pollinator group richness")+
  theme(title = element_text(size=12)) +
  theme(axis.title.x = element_text(size=15)) + 
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.text.x=element_text(size=13)) +
  theme(axis.text.y=element_text(size=13)) +
  theme(legend.text=element_text(size=13))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
  
plot


library(rPref) 
library(igraph)
library(dplyr)
# Calculate Skyline 
sky1 <- psel(yield.rich.means, high(richness) * high(Total_yield)) 
sky1

plot <- plot + geom_point(data = sky1, size = 3) + geom_line(data=sky1, colour="black")
plot

ggsave("figure5.jpg", dpi=300, height=4, width=5, units="in")

