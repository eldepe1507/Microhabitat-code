# Analysis of effect of Crevice size on Lobophora recruit numbers - Tile experiment
# Last edited 29/04/2018


# Libraries
library(ggplot2)
library(lme4)
library(Rmisc)
library(dplyr)
library (MuMIn)
library(effects)
library(tidyverse) # for data manipulation
#devtools::install_github('bbolker/glmmadmb')
library(glmmADMB) # better if there are convergence issues
library(glmmTMB) # better if there are convergence issues



# Import and set-up
remove(list=ls()) # clear everything
getwd() # check working directory
recruit.data1 <- read.csv("data/Lobophora-recruitment-crevicesize.csv") # import data set



# Let's check data and prepare it
str(recruit.data1)
recruit.data1$Tile.number <- as.factor(recruit.data1$Tile.number) # number doesn't conceil information so it will be a factor


# =======================================
# Graph preparation 

 
# To make a graph with Recruit number per Area, I will make new variable for each location
# First only for Lobophora recruits and only without taking sides into account for area

# target crevice:
recruit.data1$LOB.target.percm2.nosides <- recruit.data1$LOB.targetcrevice/recruit.data1$Area.targetcrevice.cm2.nosides
head(recruit.data1)

recruit.data1$LOB.crown.percm2 <- recruit.data1$LOB.crown/recruit.data1$Area.crown.cm2
head(recruit.data1)

recruit.data1$LOB.other.percm2.nosides <- recruit.data1$LOB.othercrevice/recruit.data1$Area.othercrevice.cm2.nosides
head(recruit.data1)

recruit.data1$LOB.open.percm2.nosides <- recruit.data1$LOB.opencrevice/recruit.data1$Area.opencrevice.cm2.nosides
head(recruit.data1)

recruit.data1$LOB.outside.percm2.nosides <- recruit.data1$LOB.outside/recruit.data1$Area.outside.nosides.cm2
head(recruit.data1)

recruit.data1$LOB.overall.percm2.nosides <- recruit.data1$LOB.overall/recruit.data1$Area.overall.cm2.nosides
head(recruit.data1)


# Let's do the same including sides

recruit.data1$LOB.target.percm2.inclsides <- recruit.data1$LOB.targetcrevice/recruit.data1$Area.targetcrevice.cm2.inclsides
head(recruit.data1)

recruit.data1$LOB.other.percm2.inclsides <- recruit.data1$LOB.othercrevice/recruit.data1$Area.othercrevice.cm2.inclsides
head(recruit.data1)

recruit.data1$LOB.open.percm2.inclsides <- recruit.data1$LOB.opencrevice/recruit.data1$Area.opencrevice.cm2.inclsides
head(recruit.data1)

recruit.data1$LOB.outside.percm2.inclsides <- recruit.data1$LOB.outside/recruit.data1$Area.outside.inclsides.cm2
head(recruit.data1)

recruit.data1$LOB.overall.percm2.inclsides <- recruit.data1$LOB.overall/recruit.data1$Area.overall.cm2.inclsides
head(recruit.data1)

#========================================
# GRAPHS

#----------------------------------------
# Let's have a first look at the number of recruits in different crevices without including sides

ggplot(recruit.data1, aes(x=crevice.kind, y=LOB.target.percm2.nosides))+
  geom_point()+
  facet_wrap(~Treatment)+
  theme_classic()


# Let's look at mean and errorbars
LOB.target.sum1 <- summarySE(recruit.data1, measurevar="LOB.target.percm2.nosides", groupvars=c("crevice.kind", "Treatment"))

# And in a plot
ggplot(LOB.target.sum1, aes(x=crevice.kind, y=LOB.target.percm2.nosides))+
  geom_col()+
  geom_errorbar(aes(ymin=LOB.target.percm2.nosides-se, ymax=LOB.target.percm2.nosides+se), width=0.1)+
  facet_wrap(~Treatment)+
  theme_classic()+
  xlab("Crevice kind")+
  ylab("Lobophora recruits in crevice (per cm2)")



# And now the same graph but including sides in area
ggplot(recruit.data1, aes(x=crevice.kind, y=LOB.target.percm2.inclsides))+
  geom_point()+
  facet_wrap(~Treatment)+
  theme_classic()


# Let's look at mean and errorbars
LOB.target.sum2 <- summarySE(recruit.data1, measurevar="LOB.target.percm2.inclsides", groupvars=c("crevice.kind", "Treatment"))

# And in a plot
ggplot(LOB.target.sum2, aes(x=crevice.kind, y=LOB.target.percm2.inclsides))+
  geom_col()+
  geom_errorbar(aes(ymin=LOB.target.percm2.inclsides-se, ymax=LOB.target.percm2.inclsides+se), width=0.1)+
  facet_wrap(~Treatment)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  xlab("Crevice kind")+
  ylab("Lobophora recruits in crevice (per cm2)")
ggsave('graphs/LOBrecruits-targetcrevice-treatments.eps')

# --> I will have to include the sides for area!


# Make a proportion of target:overal to account for influence of overall recruit number on results
recruit.data1$propLOB <- recruit.data1$LOB.target.percm2.inclsides/recruit.data1$LOB.overall.percm2.inclsides
str(recruit.data1)
LOB.prop.sum1 <- summarySE(recruit.data1, measurevar="propLOB", groupvars=c("crevice.kind", "Treatment"))
recruit.data1$propLOB
write.csv(recruit.data1, file="data/Recruitment.data.extendedR.csv")

ggplot(LOB.prop.sum1, aes(x=crevice.kind, y=propLOB))+
  geom_col()+
  geom_errorbar(aes(ymin=propLOB-se, ymax=propLOB+se), width=0.1)+
  facet_wrap(~Treatment)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  xlab("Crevice kind")+
  ylab("Proportion Lobophora recruits (cm2) in target crevice:overall")
ggsave("graphs/proportion_LOBtarget:LOBoverall.jpg")

# proportion of recruits in crevice to crown
recruit.data1$propLOBcrown <- recruit.data1$LOB.target.percm2.inclsides/recruit.data1$LOB.crown.percm2
str(recruit.data1)
LOB.prop.sum1 <- summarySE(recruit.data1, measurevar="propLOBcrown", groupvars=c("crevice.kind", "Treatment"))
recruit.data1$propLOBcrown
write.csv(recruit.data1, file="data/Recruitment.data.extendedR.csv")

ggplot(LOB.prop.sum1, aes(x=crevice.kind, y=propLOBcrown))+
  geom_col()+
  geom_errorbar(aes(ymin=propLOBcrown-se, ymax=propLOBcrown+se), width=0.1)+
  facet_wrap(~Treatment)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  xlab("Crevice kind")+
  ylab("Proportion Lobophora recruits (cm2) in target crevice:crown")
ggsave("graphs/proportion_LOBtarget:LOBcrown.jpg")
# doesn't work because there's too many zeros on the crown of partial cages!


#====================================
# Data exploration - graphs

# Relationship between Lobophora recruits in target crevice and overall number of Lobophora recruits
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x=LOB.overall.percm2.inclsides))+
  geom_point()+
  geom_smooth(method=lm)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  ylab("Lobophora recruits in target crevice (per cm2)")+
  xlab("Lobophora recruits on tile (per cm2)")
ggsave("graphs/Relationship-recruits-targetcrevice-overall.jpg")


# Relationship between Lobophora recruits in target crevice and overall number of Lobophora recruits
# split up 
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x=LOB.overall.percm2.inclsides, fill= crevice.kind))+
  geom_point()+
  geom_smooth(method=lm)+
  facet_wrap(~Treatment)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  ylab("Lobophora recruits in target crevice (per cm2)")+
  xlab("Lobophora recruits on tile (per cm2)")
ggsave("graphs/Relationship-recruits-targetcrevice-overall-splitup.jpg")


# Relationship between Lobophora recruits in target crevice and on crown
# split up 
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x=LOB.crown.percm2, fill= crevice.kind))+
  geom_point()+
  geom_smooth(method=lm)+
  facet_wrap(~Treatment)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  ylab("Lobophora recruits in target crevice (per cm2)")+
  xlab("Lobophora recruits on crown (per cm2)")
ggsave("graphs/Relationship-recruits-targetcrevice-crown-splitup.jpg")

# average recruits on crown
LOB.crown.sum <- summarySE(recruit.data1, measurevar="LOB.crown.percm2", groupvars= c("Treatment","crevice.kind"))
str(LOB.crown.sum)
ggplot(LOB.crown.sum, aes(y=LOB.crown.percm2, x= crevice.kind))+
  geom_col()+
  theme_bw() +
  facet_wrap(~Treatment)+
  geom_errorbar(aes(ymin=LOB.crown.percm2-se, ymax=LOB.crown.percm2+se), width=0.1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  ylab("Lobophora recruits on crown (per cm2)")+
  xlab("Tile kind")


# Average Lobophora recruits in target crevice per treatment
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x=Treatment))+
  geom_point()


# Average Lobophora recruits in target crevice per crevice
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x= crevice.kind))+
  geom_point()

# Average overall Lobophora recruits per treatment
ggplot(recruit.data1, aes(y=LOB.overall.percm2.inclsides, x=Treatment))+
  geom_point()

# Average overall Lobophora recruits per treatment
LOB.overall.sum <- summarySE(recruit.data1, measurevar="LOB.overall.percm2.inclsides", groupvars= c("Treatment", "crevice.kind"))
LOB.overall.sum
ggplot(LOB.overall.sum, aes(y=LOB.overall.percm2.inclsides, x=crevice.kind))+
  geom_point()+
  facet_wrap(~Treatment)+
  geom_errorbar(aes(ymin=LOB.overall.percm2.inclsides-se, ymax=LOB.overall.percm2.inclsides+se), width=0.1)+
  ylab("Lobophora recruits per cm2 on tile")+
  xlab("Tile with crevice kind")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave("graphs/Lobrecruits-overall-pertreatment-and-crevicekind.jpg")


# Average Lobophora recruits on crowns per treatment - points
ggplot(recruit.data1, aes(y=LOB.crown.percm2, x=Treatment))+
  geom_point()


# Average Lobophora recruits on crowns per treatment - column
LOB.crown.sum <- summarySE(recruit.data1, measurevar="LOB.crown.percm2", groupvars= "Treatment")

ggplot(LOB.crown.sum, aes(y=LOB.crown.percm2, x=Treatment))+
  geom_col()+
  geom_errorbar(aes(ymin=LOB.crown.percm2-se, ymax=LOB.crown.percm2+se), width=0.1)+
  ylab("Lobophora recruits on crown (per cm2)")+
  xlab("Treatment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave("graphs/Lobrecruits-oncrown-pertreatment.jpg")


# Look at the recruit numbers in target crevice for each set-up
ggplot(recruit.data1, aes(y=LOB.target.percm2.inclsides, x=crevice.kind))+
  geom_point()+
  facet_wrap(~Assigned.location.code)

# Linegraph average Lobophora recruits in target crevice over crevice kind, coloured by treatment
ggplot(LOB.target.sum2, aes(y=LOB.target.percm2.inclsides, x=crevice.kind, colour=Treatment))+
  geom_point()+
  geom_line(aes(y=LOB.target.percm2.inclsides, x=crevice.kind, group=Treatment))+
  geom_errorbar(aes(ymin=LOB.target.percm2.inclsides-se, ymax=LOB.target.percm2.inclsides+se), width=0.1)+
  ylab("Lobophora recruits (per cm2)")+
  xlab("Crevice")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
  ggsave("graphs/Lobrecruits-targetcrevice-treatment-linegraph.jpg")
#------------------------------------------
# DATA EXPLORATION - statistics

# Question: Can I pool partial and open treatments as open?

# First step: Data manipulation
# 1) Select columns of interest
str(recruit.data1) # I need LOB.target.percm2.inclsides, Treatment, crevice.kind, maybe assigned location code and tile.number
data.pairwise1 <- dplyr::select(recruit.data1, LOB.target.percm2.inclsides, Treatment, crevice.kind, Assigned.location.code, Tile.number)

# 2) turn into wide format
str(data.pairwise1)
data.pairwise2 <- spread(data.pairwise1, key=crevice.kind, value=LOB.target.percm2.inclsides)
head(data.pairwise2)
# other way of doing it:
data.pairwise2 <- data.pairwise1 %>% spread(key=crevice.kind, value=LOB.target.percm2.inclsides)

# 3) Filter by treatment to only keep 2
data.pairwise.cagedpartial <- filter(data.pairwise2, Treatment %in% c('Caged', 'Partial')) %>% droplevels
str(data.pairwise.cagedpartial)
data.pairwise.cagedopen <- filter(data.pairwise2, Treatment %in% c('Caged', 'Open')) %>% droplevels
str(data.pairwise.cagedopen)
data.pairwise.openpartial <- filter(data.pairwise2, Treatment %in% c('Open', 'Partial')) %>% droplevels
str(data.pairwise.openpartial)

# Are large creivces significantly different between caged and partial?
?t.test()
t.test(large ~ Treatment, data.pairwise.cagedpartial, na.action=na.omit) # no

# Are medium crevices significantly different between caged and partial?
t.test(medium ~ Treatment, data.pairwise.cagedpartial, na.action=na.omit) # no

# Are small crevices significantly different between caged and partial?
t.test(small ~ Treatment, data.pairwise.cagedpartial, na.action=na.omit) # no

# Are large crevices significantly different between caged and open?
t.test(large ~ Treatment, data.pairwise.cagedopen, na.action=na.omit) # no

# Are medium crevices significantly different between caged and open?
t.test(medium ~ Treatment, data.pairwise.cagedopen, na.action=na.omit) # no

# Are small crevices significantly different between caged and open?
t.test(small ~ Treatment, data.pairwise.cagedopen, na.action=na.omit) # yes

# Are large crevices significantly different between partial and open?
t.test(large ~ Treatment, data.pairwise.openpartial, na.action=na.omit) # no

# Are medium crevices significantly different between partial and open?
t.test(medium ~ Treatment, data.pairwise.openpartial, na.action=na.omit) # yes

# Are small crevices significantly different between partial and open?
t.test(small ~ Treatment, data.pairwise.openpartial, na.action=na.omit) # no


# Maybe it's better to do glmer but only with two treatments
# We will need to use raw data, not per cm2
# we don't need to spread it to wide format but can go straight to filtering treatments
# Filter by treatment to only keep 2
data.pairwise.cagedpartial <- filter(recruit.data1, Treatment %in% c('Caged', 'Partial')) %>% droplevels
str(data.pairwise.cagedpartial)
data.pairwise.cagedopen <- filter(recruit.data1, Treatment %in% c('Caged', 'Open')) %>% droplevels
str(data.pairwise.cagedopen)
data.pairwise.openpartial <- filter(recruit.data1, Treatment %in% c('Open', 'Partial')) %>% droplevels
str(data.pairwise.openpartial)

# glmer
# Caged vs partial
pairwise.cp.glmer1 <- glmer(data=data.pairwise.cagedpartial, LOB.targetcrevice ~ Treatment*crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) +(1|Assigned.location.code/Tile.number), family=poisson)
summary(pairwise.cp.glmer1) 
# --> There is a significant interaction term, they are therefore different

# Caged vs open
pairwise.co.glmer1 <- glmer(data=data.pairwise.cagedopen, LOB.targetcrevice ~ Treatment*crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) +(1|Assigned.location.code/Tile.number),  family=poisson)
summary(pairwise.co.glmer1) 
# --> Both interaction terms are significant, the treatments are therefore different
dev <- deviance(pairwise.co.glmer1) # glmer we need to call deviance like this, different to glm
rdf <- nrow(data.pairwise.cagedopen) - length(fixef(pairwise.co.glmer1)) -2
dev/rdf # no overdispersion - poisson is correct


# Partial vs Open
pairwise.po.glmer1 <- glmer(data=data.pairwise.openpartial, LOB.targetcrevice ~ Treatment*crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), family=poisson)
summary(pairwise.po.glmer1) # Doesn't converge


pairwise.po.glmmTMB <- glmmTMB(data=data.pairwise.openpartial, LOB.targetcrevice ~ Treatment*crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), family=poisson)
summary(pairwise.po.glmmTMB)
# Significant interaction terms

# --> I cannot combine partial and open treatments!



#==========================================
# MODELS

#------------------------------------------
# glmer and glmmTMB with interaction

# First Poisson because it's count data
# Just as a start let's see what the end model would look like
head(recruit.data1)
LOBrecruit.glmer1 <- glmer(LOB.targetcrevice ~ Treatment * crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.glmer1)
# an error was due to not logging Area before - since it is a numerical variable and a poisson model I need to log it
# there are still converging warnings - maybe scaling Area will help
# Results sometimes make sense

# Let's try glmmTMB which has less convergence issues

LOBrecruit.glmmTMB1 <- glmmTMB(LOB.targetcrevice ~ Treatment * crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.glmmTMB1)
# no more convergence issues


# Model assumptions - glmer:
#---------------------------------------------------
# glmer
# Let's check overdispersion
dev <- deviance(LOBrecruit.glmer1)
ddf <- nrow(LOBrecruit.glmer1@frame) - length(fixef(LOBrecruit.glmer1))- length(VarCorr(LOBrecruit.glmer1)) # if it's a poisson, you don't need to substract an additional -1 at the end!!
adf <- df.residual(LOBrecruit.glmer1) # it can now calculate the degrees of freedom like that
dev/ddf # no overdispersion!
dev/adf # same result as with ddf

# Residual plot:
plot(LOBrecruit.glmer1) # hmm, at least one point far off

# R squared
r.squaredGLMM(LOBrecruit.glmer1) # 10.1 % explained by fixed factors


# All looks alright, but just to check model fit, I will also run a negative binomial model
LOBrecruit.mod.nb <- glmer.nb(LOB.targetcrevice ~ Treatment * crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), data=recruit.data1)
summary(LOBrecruit.mod.nb)
# Results are looking very similar to poisson model

# Let's check AIC:
AICc(LOBrecruit.glmer1, LOBrecruit.mod.nb) # poisson is better AND simpler!

# The Poisson model is the best one if it converges!
# I will go with the glmmTMB because it always converges!

#----------------------------------------

# Interpretation Model output:
# - More variance between location code than between tiles - good
# Significant interaction for three terms: TreatmentOpen:crevicekindmedium, TreatmentOpen:crevicekindsmall and TreatmentPartial:crevice.kindsmall
# Significant interaction for TreatmentOpen:Crevicemedium means:
  # the number of Lobophora recruits in the medium crevice of the open treatment
  # is log(1.4) higher than would be expected if I summed the effect of Open treatments if it was averaged over all crevices
  # and for medium crevices if it was averaged over all treatments
# Significant interaction of TreatmentOpen:crevicesmall and TreatmentPartial:crevice small are similarly interpreted


# I will have to back transform my Estimates to get actual useful numbers
# My reference is large crevice and caged treatment

# Difference with Open treatment and medium crevice:
exp(1.4) # = 4.0552
# --> There are 4.1 times more recruits in the medium crevice of the open treatment than you would expect
# if you summed the baseline effect of Open Treatment vs Caged treatment (for reference crevice: large) 
# and medium crevice against large crevice (for reference treatment: caged)


# Difference with Open treatment and small crevices:
exp(2.06852) # 7.9131
# --> There are 7.9 times more recruits in the small crevice of the open treatment than you would expect 
# if you summed the baseline effect of Open Treatment vs Caged treatment (for reference crevice: large) 
# and small crevice against large crevice (for reference treatment: caged)

# Difference with Partial Treatment and small crevice:
exp(1.19743) # 3.3115
# --> There are 3.3 times more recruits in the small crevice of the partial treatment than you would expect 
# if you summed the baseline effect of Partial Treatment vs Caged treatment (for reference crevice large)
# and small crevice against large crevice (for reference treatment caged)


# Difference with Partial treatment and medium crevices (insignificant):
exp(-0.099) # 0.9057 
# --> there would be 0.91 times more recruits (so actually less) in the medium crevice of the partial treatment than you would expect
# if you summed the baseline effect of Partial Treatment  vs Caged treatment (for reference crevice: large)
# and medium crevice vs large crevice (for reference treatment: caged)

# Both small and medium crevices have more recruits than large crevices if they are caged or partially caged. 
# These are significant except for the medium crevice in a partial treatment
# Within the open treatment the effect is stronger in the small crevice than the medium crevice
# Within the partial treatment the effect is stronger for the smaller crevice than the medium crevice
# BUT small crevice in partial treatment is less strong than medium crevice in open treatment
# Partial treatment is somewhere between caged and open


## For my question, nesting my fixed factors may actually get me closer to what I want to know:
#----------------------------------------
# NESTED glmmTMB model

LOBrecruit.glmmTMB2 <- glmmTMB(LOB.targetcrevice ~ Treatment/crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.glmmTMB2)

# Backtransformation and interpretation:
# baseline is large crevice and caged treatment

# medium crevice in caged treatment:
exp(-1.33) # 0.2645
# There are 0.27 times more recruits (so actually less) in the medium crevice of the caged treatment 
# than in large crevices in the caged treatment

# medium crevice in open treatment:
exp(0.070) # 1.0725
# There are 1.07 times more recruits in the medium crevice of the open treatment 
# than in large crevices in the open treatment

# medium crevice in partial treatment:
exp(-1.43) # 0.2393
# There are 0.24 times more recruits (so actually less) in the medium crevice of the partial treatment 
# than in the large crevice in the partial treatment

# small crevice in caged treatment
exp(-1.98) # 0.1380
# There are 0.14 times more recruits (so actually less) in the small crevice of the caged treatment 
# than in the large crevice of the caged treatment

# small crevice in open treatment
exp(0.09) # 1.0942
# There are 1.09 times more recruits in the small crevice of the open treatment
# than in the large crevice of the open treatment

# small crevice in partial treatment
exp(-0.78) # 0.4584
# There are 0.46 times more recruits (so actually less) in the small crevice of the partial treatment
# than in the large crevice of the partial treatment


# Model taking overall Lobophora numbers into account
# Interaction
LOBrecruit.glmmTMB3 <- glmmTMB(LOB.targetcrevice ~ Treatment * crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + offset(log(LOB.overall.percm2.inclsides))+ (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.glmmTMB3)

# Nested
LOBrecruit.glmmTMB4 <- glmmTMB(LOB.targetcrevice ~ Treatment / crevice.kind + offset(log(Area.targetcrevice.cm2.inclsides)) + offset(log(LOB.overall.percm2.inclsides))+ (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.glmmTMB4)

#-----------------------------------------
# Something I tried to make converging problem go away, but didn't work:
# let's log it first and then center it - can't log scaled values (0s)
recruit.data1$Area.targetcrevice.cm2.inclsides.log<-log(recruit.data1$Area.targetcrevice.cm2.inclsides)
recruit.data1 <- recruit.data1 %>% mutate(cArea.targetcrevice.cm2.inclsides.log = Area.targetcrevice.cm2.inclsides.log-mean(Area.targetcrevice.cm2.inclsides.log))
LOBrecruit.mod2 <- glmer(LOB.targetcrevice ~ Treatment * crevice.kind + offset(cArea.targetcrevice.cm2.inclsides.log) + (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
summary(LOBrecruit.mod2)
# works but converging errors
# Results don't make sense, it clearly didn't work
LOBrecruit.mod3 <- glmer(LOB.targetcrevice ~ Treatment + crevice.kind + offset(cArea.targetcrevice.cm2.inclsides.log) + (1|Assigned.location.code/Tile.number), data=recruit.data1, family=poisson)
# without interaction it converges, could we use this?
AIC(LOBrecruit.mod1, LOBrecruit.mod3) # No, the model including the interaction is much better
#-----------------------------------------

#==========================================
