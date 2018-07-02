# Project: Analyse influence of different fish groups on Lobophora through microhabitat accessibility
# Survey data
# Date: last edited 02 July 2018

rm(list=ls()) 
options(repos="https://cran.curtin.edu.au/")
install.packages('lme4')
library(lme4)
install.packages("MuMIn")
library(MuMIn)
library(nlme) 
install.packages("lmerTest")
library(lmerTest) 
install.packages('reshape')
library(reshape) 
install.packages("piecewiseSEM")
library(piecewiseSEM)
install.packages("caret")
library(caret)
library(dplyr)
require(MASS)
require(car)
install.packages('devtools')
library(devtools)
install.packages('ggpmisc')
library(ggpmisc)
#install.packages('grid')
library(grid)
library(ggplot2)
install.packages('tidyverse')
library(tidyverse)
install.packages('ggfortify')
library(ggfortify)
install.packages('effects')
library(effects)
library(broom)
install.packages('multcomp')
library(multcomp)


## Set working directory
getwd()
setwd("/Users/Laura/Documents/PhD/Data/Chapter 3-Microhabitat/Microhabitat utilisation")

micro.raw<-read.csv("Microhabitat_raw.csv", header=T)
micro.biom.data <- read.csv("MH_avgbiomass.csv", header = T)
micro.raw.na <- read.csv('Microhabitat_raw_openingNA.csv')

summary(micro.raw)
str(micro.raw)

## drop the Lighthouse Reef site
micro.data<-droplevels(subset(micro.raw, site=="East Sheltered"))
str(micro.data)
micro.data$quadrat <- as.factor(micro.data$quadrat) 
micro.data$subquadrat <- as.factor(micro.data$subquadrat)
micro.data.na <- droplevels(subset(micro.raw.na, site=='East Sheltered'))
str(micro.data.na)
micro.data.na$quadrat <- as.factor(micro.data.na$quadrat) 
micro.data.na$subquadrat <- as.factor(micro.data.na$subquadrat)
micro.data.na$minimum_opening_mm <- as.numeric(micro.data.na$minimum_opening_mm)
micro.data.na$crevice_ratio_opening <- as.numeric(micro.data.na$crevice_ratio_opening)

# Standardise data for the fixed factors to get around warning message
micro.data$sdt_GrazingPressure_biomass <- scale(micro.data$GrazingPressure_biomass)
micro.data$sdt_crevice_width_mm <- scale(micro.data$crevice_width_mm)
micro.data$sdt_crevice_length_mm <- scale(micro.data$crevice_length_mm)
micro.data$sdt_crevice_depth_mm <- scale(micro.data$crevice_depth_mm)
micro.data$sdt_minimum_opening_mm <- scale(micro.data$minimum_opening_mm)
micro.data$sdt_crevice_opening_mm2 <- scale(micro.data$crevice_opening_mm2)

micro.data.na$sdt_GrazingPressure_biomass <- scale(micro.data.na$GrazingPressure_biomass)
micro.data.na$sdt_crevice_width_mm <- scale(micro.data.na$crevice_width_mm)
micro.data.na$sdt_crevice_length_mm <- scale(micro.data.na$crevice_length_mm)
micro.data.na$sdt_crevice_depth_mm <- scale(micro.data.na$crevice_depth_mm)
micro.data.na$sdt_minimum_opening_mm <- scale(micro.data.na$minimum_opening_mm)
micro.data.na$sdt_crevice_opening_mm2 <- scale(micro.data.na$crevice_opening_mm2)
## summary of micros with lobophora in each quadrat
with(micro.data, tapply(lobophora, quadrat, mean))



#---------------------------------------------------------------------------------------------------------------------------------------------------
# STATISTICS


# ABIOTIC PARAMETERS
# FIND BEST CREVICE MEASUREMENTS TO REPRESENT LOBOPHORA LIKELIHOOD

# let's check for correlation between grazing pressure and lobophora likelihood

# One way to test for correlation is the Pearson 's product-moment correlation
# Pearson analysis linear correlation # 1 is completely correlated, 0 is not at all

cor.graz.depth <-cor.test(micro.data$GrazingPressure_biomass, micro.data$crevice_depth_mm)
cor.graz.depth # there is a significant negative correlation between depth and grazing pressure , but only 0.24

cor.graz.length <- cor.test(micro.data$GrazingPressure_biomass, micro.data$crevice_length_mm)
cor.graz.length # there is a significant positive correlation between depth and grazing pressure of 0.47

cor.graz.width <- cor.test(micro.data$GrazingPressure_biomass, micro.data$crevice_width_mm)
cor.graz.width # there is a significant positive correlation between depth and grazing pressure of 0.49

cor.graz.min.opening <- cor.test(micro.data$GrazingPressure_biomass, micro.data$minimum_opening_mm)
cor.graz.min.opening # there is a significant positive correlation between depth and grazing pressure of 0.20

cor.graz.opening <- cor.test(micro.data$GrazingPressure_biomass, micro.data$crevice_opening_mm2)
cor.graz.opening # there is a significant positive correlation between depth and grazing pressure of 0.49



# Which crevice dimensions can explain Lobophora likelihood?

# relationship of lobophora with crevice depth, width and length
abiotic.m1 <- glmer(lobophora ~ sdt_crevice_depth_mm + sdt_crevice_width_mm + sdt_crevice_length_mm + (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m1) # AIC 482.1

# relationship of lobophora with crevice depth and minimum opening
abiotic.m2 <- glmer(lobophora ~ crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m2) # AIC 478.3 best model fit, incresing width decreases Lobophora likelihood significantly, depth is marginally insignificant
exp(-0.75) # 0.47 - one increase in width means 0.47 times the previous likelihood of Lobophora

# relationship of lobophora with crevice depth
abiotic.m3 <- glmer(lobophora ~ crevice_depth_mm + (1|quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m3) # non significant, AIC = 502.3

# relationship of lobophora with minimum opening
abiotic.m4 <- glmer(lobophora ~ sdt_minimum_opening_mm + (1|quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m4) # negative, significant, AIC = 480.7

# relationship of lobophora with crevice opening (length * width)
abiotic.m5 <- glmer(lobophora ~ sdt_crevice_opening_mm2 + (1|quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m5) #significant AIC 481.5

# relationship of lobophora with crevice ratio (depth/minimum opening)
abiotic.m6 <- glmer(lobophora ~ crevice_ratio_min + (1|quadrat/subquadrat), family = binomial, data = micro.data)
summary(abiotic.m6) # significant, actual effect is larger than just minimum opening, but higher AIC 495.5


# relationship of lobophora with crevice orientation
abiotic.m9 <- glmer(lobophora ~ surface_orientation + (1|quadrat/subquadrat), family= binomial, data = micro.data)
summary(abiotic.m9) # AIC 497.9

# relationship of lobophora with surface shape
abiotic.m10 <- glmer(lobophora ~ surface_shape + (1|quadrat/subquadrat), family= binomial, data= micro.data)
summary(abiotic.m10) # AIC 485.8

# relationsip of lobophora with crevice minimum opening and surface shape
abiotic.m11 <- glmer(lobophora ~ surface_shape + sdt_minimum_opening_mm + (1|quadrat/subquadrat), family= binomial, data= micro.data)
summary(abiotic.m11) # AIC 483.2

# relationship of lobophora with openness
abiotic.m12 <- glmer(lobophora ~ openness + (1|quadrat/subquadrat), family= binomial, data=micro.data)
summary(abiotic.m12) # AIC 487.1

# relationship of lobophora with openness + crevice depth and width
abiotic.m13 <- glmer(lobophora ~ openness + sdt_minimum_opening_mm + sdt_crevice_depth_mm + (1|quadrat/subquadrat), family= binomial, data=micro.data)
summary(abiotic.m13) # AIC 479.9 - would be similar as crevice depth and width (deltaAIC , 2),but it's more complicated, so we are choosing the easier model


# relationship of lobophora with crevice ratio (depth/crevice opening(=length*width))
#WARNING: model takes very long to calculate
#abiotic.m7 <- glmer(lobophora ~ crevice_ratio_opening + (1|quadrat/subquadrat), family = binomial, data = micro.data)
#summary(abiotic.m7) # model failed to converge majorly, non significant, AIC 608.4

# Based on AIC, best predictor for Lobophora likelihood in microhabitats seems to be minimum opening and crevice depth
# Let's see if we can visualise it

ggplot(micro.data, aes(y= lobophora, x=minimum_opening_mm))+
  geom_point()+
  geom_smooth() # could be linear

ggplot(micro.data, aes(y= lobophora, x=crevice_depth_mm))+
  geom_point()+
  geom_smooth() # not linear - shows the weird increase and then decrease 


# Rerun model without random effects for plotting
abiotic.m8 <- glm(lobophora ~ minimum_opening_mm, data = micro.data)
summary(abiotic.m8)

range(micro.data$minimum_opening_mm) ## gives range of GrazingPressure_abundance
xgrazing <- seq(3,300,0.1) ## creates a sequence of values to produce fitted values to
ygrazing <- predict(abiotic.m8, list(minimum_opening_mm=xgrazing), type="response") ## creates model for all values xgrazing
plot(micro.data$minimum_opening_mm, micro.data$lobophora, pch=16, col="aquamarine3", xlab="Crevice width", ylab="Lobophora likelihood") ## plots data points
lines(xgrazing, ygrazing, col="aquamarine4") ## adds trendline
# -> Probability of Lobophora decreases with increasing size of the crevice


# same for depth
abiotic.m9 <- glm(lobophora ~ crevice_depth_mm, data = micro.data)
summary(abiotic.m9)

range(micro.data$crevice_depth_mm) ## gives range of GrazingPressure_abundance
xgrazing <- seq(3,300,0.1) ## creates a sequence of values to produce fitted values to
ygrazing <- predict(abiotic.m9, list(crevice_depth_mm=xgrazing), type="response") ## creates model for all values xgrazing
plot(micro.data$crevice_depth_mm, micro.data$lobophora, pch=16, col="aquamarine3", xlab="Crevice depth", ylab="Lobophora likelihood") ## plots data points
lines(xgrazing, ygrazing, col="aquamarine4") ## adds trendline
# -> Probability of Lobophora decreases with increasing size of the crevice

dataless300 <- subset(micro.data, !minimum_opening_mm==300)
dataless300


# Let's plot just the data
# including all observations, i.e. 300 for the ones that are completely open

ggplot(micro.data, aes(y= lobophora, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Lobophora likelihood") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of Lobophora and crevize size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/lobo-mincrevice-all.eps')

# And the same graph excluding the ones marked 300, as they are technically open

ggplot(dataless300, aes(y= lobophora, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Lobophora likelihood") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of Lobophora and crevize size") +
  theme_classic()+
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))

ggsave('graphs/lobo-mincrevice-no300.eps')


# And the same kind of graph for the second parameter added to the model: crevice depth
ggplot(micro.data, aes(y= lobophora, x = crevice_depth_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Lobophora likelihood") +
  xlab("Crevice depth (mm)") +
  ggtitle("Relationship of Lobophora and crevize depth") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/lobo-crevice.depth.png')
# interesting. Likelihood of Lobophora first increases with depth (makes sense), but seems to decrease once the crevie gets deeper than ~ 20 mm


# Density plot with two predictors - crevice depth and crevice width
#install.packages('MASS')
require(MASS)
#install.packages('ks')
library(ks)

# kernel density plot with all crevices
crevice.kde <- kde2d(x= micro.data$minimum_opening_mm, y = micro.data$crevice_depth_mm,  n=100, lims=c(0, 300, 0, 150))
crevice.kde

# graph 1 - rainbow - Density of all microhabitats
par(mar=c(7,8,5,5))
filled.contour(crevice.kde, 
               color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')), 
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# graph 2 - red yellow
jpeg('graphs/crevice.density.all.jpg')
image(crevice.kde, xlab= 'minimum crevice width(mm)', ylab='crevice depth(mm)')
contour(crevice.kde, add=TRUE)
dev.off()


crevice.lobo <- subset(micro.data, lobophora==1)
crevice.nolobo <- subset(micro.data, lobophora==0)

# kernel density plot with lobophora
crevice.lobo.kde <- kde2d(x=crevice.lobo$minimum_opening_mm, y= crevice.lobo$crevice_depth_mm, n=100, lims=c(0, 300, 0, 150))
crevice.lobo.kde

# red and yellow plot
jpeg('graphs/crevice.density.lobo.jpg')
image(crevice.lobo.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.lobo.kde,  add=TRUE)
dev.off()

# rainbow plot - crevices with Lobophora

par(mar=c(7,8,5,5))
filled.contour(crevice.lobo.kde, 
               color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')), 
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })
#mtext(side=3, text='Density of microhabitats with Lobophora', cex=3, line=1, adj=0)



# Difference in densities
difference.data <- crevice.lobo.kde
difference.data$z <-  crevice.lobo.kde$z -crevice.kde$z

par(mar=c(7,8,5,5))
filled.contour(difference.data, 
               color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')), 
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })
#mtext(side=3, text='Difference between all microhabitats and crevices with Lobophora', cex=1.5, line=1, adj=-0.5)



#==========================================================================================

# RELATIONSHIP BETWEEN GRAZING AND ABIOTIC PARAMETERS 

#-------------------------------------------------------------------------------------------
# surface shape and orientation first


# let's check for CORRELATION between GRAZING PRESSURE and SURFACE SHAPE
# following a description online (https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable)

ggplot(micro.data, aes(y= GrazingPressure_biomass, x = surface_shape))+
  geom_boxplot() +
  scale_y_continuous(limits=c(0,150)) # definitely not normally distributed

model.cor.grazbio.shape <- lm(GrazingPressure_biomass ~ surface_shape, data = micro.data)
summary(model.cor.grazbio.shape) # estimate gives estimated grazing pressure within these
# coefficient of determination R^2 = 0.2137, correlation usually above 0.3 so should be fine
# this means 21 % of variance is explained by our model
# by square-rooting, we can find the multiple correlation coefficient R
rsq <- summary(model.cor.grazbio.shape)$r.squared
rsq # this is the correlation between the observed durations and the ones predicted by our model

# alternative is to check with anova
aov.cor.grazbio.shape <- aov(GrazingPressure_biomass ~ surface_shape, data = micro.data)
summary(aov.cor.grazbio.shape) # significant -- probably correlated
# equivalent to R^2 is eta squared in anova:

install.packages("heplots")
library(heplots)
etasq(aov.cor.grazbio.shape, partial = TRUE) # 0.2137 -- the same as with correlation above

# probably not the best idea to put these into the same model


# let's check for CORRELATION between GRAZING PRESSURE and surface orientation
# following a description online (https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable)

micro.data1 <- micro.data[!is.na(micro.data$surface_orientation), ]
ggplot(micro.data1, aes(y= GrazingPressure_biomass, x = surface_orientation))+
  geom_boxplot() +
  scale_y_continuous(limits=c(0,150)) 

model.cor.grazbio.orientation <- lm(GrazingPressure_biomass ~ surface_orientation, data = micro.data1)
summary(model.cor.grazbio.orientation) # estimate gives estimated grazing pressure within these
# coefficient of determination R^2 = 0.035
# this means 3.5 % of variance is explained by our model
# by square-rooting, we can find the multiple correlation coefficient R
rsq1 <- summary(model.cor.grazbio.orientation)$r.squared
rsq1 # 0.039 this is the correlation between the observed durations and the ones predicted by our model

# alternative is to check with anova
aov.cor.grazbio.orientation <- aov(GrazingPressure_biomass ~ surface_orientation, data = micro.data1)
summary(aov.cor.grazbio.orientation) # significant -- probably correlated
# equivalent to R^2 is eta squared in anova:

#install.packages("heplots")
#library(heplots)
etasq(aov.cor.grazbio.orientation, partial = TRUE) # 0.039 -- the same as with correlation above


#-------------------------------------------------------------------------------------------

# RELATIONSHIP BETWEEN CREVICE ATTRIBUTES AND ACCESS BY HERBIVORE SPECIES

# 1) NASO 15 CM 

#removing the 'open' microhabitats, i.e. 300mm
nasoaccess15.crevice.m1 <- glmer(naso_lituratus_15_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(nasoaccess15.crevice.m1) # significant negative impact of minimum opening, AIC 227.8

capture.output(summary(nasoaccess15.crevice.m1), file = 'stats/Routput/nasoaccess15.crevice.doc')

# including all microhabitats
nasoaccess15.crevice.m2 <- glmer(naso_lituratus_15_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(nasoaccess15.crevice.m2) # significant negative impact of depth and positive impact of increasing minimum opening, AIC 219.0

capture.output(summary(nasoaccess15.crevice.m2), file = 'stats/Routput/nasoaccess15.allcrevice.crevice.doc')


# make separate data sets of no acces by Naso lituratus 15 cm for kernel density plots

crevice.nonaso15 <- subset(micro.data, naso_lituratus_15_cm==0)

# kernel density plot without naso 15 cm access

crevice.nonaso15.kde <- kde2d(x=crevice.nonaso15$minimum_opening_mm, y= crevice.nonaso15$crevice_depth_mm, n=200, lims=c(0, 300, 0, 150))

# a new way to make a graph
filled.contour(crevice.nonaso15.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# normal kernel density plot
jpeg('graphs/crevice.density.nonaso15.jpg')
image(crevice.nonaso15.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nonaso15.kde,  add=TRUE)
dev.off()

ggplot(crevice.nonaso15, aes(minimum_opening_mm, crevice_depth_mm))+
  geom_density2d(aes(colour='rainbow'))+
  scale_x_continuous(limits= c(0,300))+
  scale_y_continuous(limits=c(0,140))+
  xlab('Microhabitat width (mm)')+
  ylab('Microhabitat depth (mm)')

?geom_density2d()
# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open

ggplot(micro.data, aes(y= naso_lituratus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (15 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (15 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso15-mincrevice-all.jpg')

# excluding open microhabitats, i.e. > 300 mm
ggplot(micro.data.na, aes(y= naso_lituratus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (15 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (15 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso15-mincrevice-no300.jpg')


#depth was also significant, therefore:
ggplot(micro.data, aes(y= naso_lituratus_15_cm, x = crevice_depth_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (15 cm) access')) +
  xlab("crevice depth (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (15 cm) and crevice depth')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso15-crevicedepth.jpg')



# 2) NASO 20 CM
#removing the 'open' microhabitats, i.e. 300mm
nasoaccess20.crevice.m1 <- glmer(naso_lituratus_20_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(nasoaccess20.crevice.m1) # marginally insignificant (0.055) negative impact of minimum opening, AIC 259.1

capture.output(summary(nasoaccess20.crevice.m1), file='stats/Routput/nasoaccess20.crevice.doc')

#including all microhabitats
nasoaccess20.crevice.m2 <- glmer(naso_lituratus_20_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(nasoaccess20.crevice.m2) # significant positive impact of incresing minimum opening, AIC 280

capture.output(summary(nasoaccess20.crevice.m2), file='stats/Routput/nasoaccess20.allcrevice.crevice.doc')


# make separate data sets of no access by Naso lituratus 20 cm for kernel density plots

crevice.nonaso20 <- subset(micro.data, naso_lituratus_20_cm==0)

# kernel density plot without naso 20 cm access

crevice.nonaso20.kde <- kde2d(x=crevice.nonaso20$minimum_opening_mm, y= crevice.nonaso20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


jpeg('graphs/crevice.density.nonaso20.jpg')
image(crevice.nonaso20.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nonaso20.kde,  add=TRUE)
dev.off()

# a new way to make a graph
filled.contour(crevice.nonaso20.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= naso_lituratus_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (20 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (20 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso20-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= naso_lituratus_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (20 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (20 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso20-mincrevice-no300.jpg')


# 3) NASO 25 CM
#excluding 'open' microhabitats, i.e. 300mm
nasoaccess25.crevice.m1 <- glmer(naso_lituratus_25_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(nasoaccess25.crevice.m1) # significant negative impact of minimum opening and crevice depth AIC 272.8

capture.output(summary(nasoaccess25.crevice.m1), file='stats/Routput/nasoaccess25.crevice.doc')

#including all microhabitats
nasoaccess25.crevice.m2 <- glmer(naso_lituratus_25_cm ~ crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(nasoaccess25.crevice.m2) # significant negative impact of crevice depth and positive impact of increasing minimum opening AIC 322.4

capture.output(summary(nasoaccess25.crevice.m2), file='stats/Routput/nasoaccess25.allcrevice.crevice.doc')


# make separate data sets of no access by Naso lituratus 25 cm for kernel density plots

crevice.nonaso25 <- subset(micro.data, naso_lituratus_25_cm==0)

# kernel density plot without naso 25 cm access

crevice.nonaso25.kde <- kde2d(x=crevice.nonaso25$minimum_opening_mm, y= crevice.nonaso25$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nonaso25.jpg')
image(crevice.nonaso25.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nonaso25.kde,  add=TRUE)
dev.off()

# a new way to make a graph
filled.contour(crevice.nonaso25.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= naso_lituratus_25_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (25 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (25 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso25-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= naso_lituratus_25_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (25 cm) access')) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression('Relationship of access by '*italic(N.~lituratus)*' (25 cm) and crevice size')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso25-mincrevice-no300.jpg')


#depth was also significant, therefore:
ggplot(micro.data, aes(y= naso_lituratus_25_cm, x = crevice_depth_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression('Likelihood of '*italic(N.~lituratus)*' (25 cm) access')) +
  xlab("crevice depth (mm)") +
  ggtitle(expression("Relationship of access by "*italic(N.~lituratus)*" (25 cm) and crevice depth")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/naso25-crevicedepth.jpg')

# 4) ZEBRASOMA SCOPAS 5 CM

#excluding 'open' microhabitats, i.e. > 300 mm
zebraaccess5.crevice.m1 <- glmer(zebrasoma_scopas_5_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(zebraaccess5.crevice.m1) # not significant, AIC 161.4

capture.output(summary(zebraaccess5.crevice.m1), file='stats/Routput/Zebrasomaaccess5.crevice.doc')

#including all microhabitats
zebraaccess5.crevice.m2 <- glmer(zebrasoma_scopas_5_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(zebraaccess5.crevice.m2) # not significant, AIC 300.8

capture.output(summary(zebraaccess5.crevice.m2), file='stats/Routput/Zebrasomaaccess5.allcrevices.crevice.doc')


# make separate data sets of no access by Zebrasoma scopas 5 cm for kernel density plots

crevice.nozebra5 <- subset(micro.data, zebrasoma_scopas_5_cm==0)

# kernel density plot without zebrasoma 5 cm access

crevice.nozebra5.kde <- kde2d(x=crevice.nozebra5$minimum_opening_mm, y= crevice.nozebra5$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nozebra5.jpg')
image(crevice.nozebra5.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nozebra5.kde,  add=TRUE)
dev.off()




# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= zebrasoma_scopas_5_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Z.~scopas)*" (5 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Z.~scopas)*" (5cm) and crevize size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/zebrasoma5-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= zebrasoma_scopas_5_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Z.~scopas)*" (5 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Z.~scopas)*" (5cm) and crevize size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/zebrasoma5-mincrevice-no300.jpg')


# 5) ZEBRASOMA SCOPAS 10 CM

#excluding 'open' microhabitats, i.e. > 300mm
zebraaccess10.crevice.m1 <- glmer(zebrasoma_scopas_10_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(zebraaccess10.crevice.m1) # not significant, AIC 174.8

capture.output(summary(zebraaccess10.crevice.m1), file='stats/Routput/Zebrasomaaccess10.crevice.doc')

#including all microhabitats
zebraaccess10.crevice.m2 <- glmer(zebrasoma_scopas_10_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(zebraaccess10.crevice.m2) # not significant, AIC 317.2

capture.output(summary(zebraaccess10.crevice.m2), file='stats/Routput/Zebrasomaaccess10.allcrevices.crevice.doc')

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= zebrasoma_scopas_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Z.~scopas)*" (10 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Z.~scopas)*" (10 cm) and crevize size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/zebrasoma10-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= zebrasoma_scopas_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Z.~scopas)*" (10 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Z.~scopas)*" (10 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/zebrasoma10-mincrevice-no300.jpg')



# make separate data sets of no access by Zebrasoma scopas 10 cm for kernel density plots

crevice.nozebra10 <- subset(micro.data, zebrasoma_scopas_10_cm==0)

# kernel density plot without zebrasoma scopas 10 cm access

crevice.nozebra10.kde <- kde2d(x=crevice.nozebra10$minimum_opening_mm, y= crevice.nozebra10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nozebra10.jpg')
image(crevice.nozebra10.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nozebra10.kde,  add=TRUE)
dev.off()




# 6) SIGANUS VULPINUS 10 CM

#excluding 'open' microhabitats, i.e. > 300mm
siganusaccess10.crevice.m1 <- glmer(siganus_vulpinus_10_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(siganusaccess10.crevice.m1) # not significant, AIC 45.0

capture.output(summary(siganusaccess10.crevice.m1), file='stats/Routput/siganusaccess10.crevice.doc')

#including all microhabitats
siganusaccess10.crevice.m2 <- glmer(siganus_vulpinus_10_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(siganusaccess10.crevice.m2) # not significant, AIC 46.6

capture.output(summary(siganusaccess10.crevice.m2), file='stats/Routput/siganusaccess10.allcrevices.crevice.doc')

# make separate data sets of no access by Siganus vulpinus 10 cm for kernel density plots

crevice.nosiganus10 <- subset(micro.data, siganus_vulpinus_10_cm==0)

# kernel density plot without Siganus vulpinus 10 cm access

crevice.nosiganus10.kde <- kde2d(x=crevice.nosiganus10$minimum_opening_mm, y= crevice.nosiganus10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nosiganus10.jpg')
image(crevice.nosiganus10.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nosiganus10.kde,  add=TRUE)
dev.off()

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= siganus_vulpinus_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(S.~vulpinus)*" (10 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(S.~vulpinus)*" (10 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/siganus10-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= siganus_vulpinus_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(S.~vulpinus)*" (10 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(S.~vulpinus)*" (10 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/siganus10-mincrevice-no300.jpg')


# 7) SIGANUS VULPINUS 15 CM

#excluding 'open' microhabitats, i.e. > 300mm
siganusaccess15.crevice.m1 <- glmer(siganus_vulpinus_15_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(siganusaccess15.crevice.m1) # not significant, AIC 59.9

capture.output(summary(siganusaccess15.crevice.m1), file='stats/Routput/siganusaccess15.crevice.doc')

#including all microhabitats
siganusaccess15.crevice.m2 <- glmer(siganus_vulpinus_15_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(siganusaccess15.crevice.m2) # not significant, AIC 61.2

capture.output(summary(siganusaccess15.crevice.m2), file='stats/Routput/siganusaccess15.allcrevices.crevice.doc')


# make separate data sets of no access by Siganus vulpinus 15 cm for kernel density plots

crevice.nosiganus15 <- subset(micro.data, siganus_vulpinus_15_cm==0)

# kernel density plot without Siganus vulpinus 15 cm access

crevice.nosiganus15.kde <- kde2d(x=crevice.nosiganus15$minimum_opening_mm, y= crevice.nosiganus15$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nosiganus15.jpg')
image(crevice.nosiganus15.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nosiganus15.kde,  add=TRUE)
dev.off()


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= siganus_vulpinus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(S.~vulpinus)*" (15 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(S.~vulpinus)*" (15 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/siganus15-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= siganus_vulpinus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(S.~vulpinus)*" (15 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(S.~vulpinus)*" (15 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/siganus15-mincrevice-no300.jpg')

# 8) CTENOCHAETUS STRIATUS 10 CM

#excluding 'open' microhabitats, i.e. > 300mm
ctenoaccess10.crevice.m1 <- glmer(ctenochaetus_striatus_10_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(ctenoaccess10.crevice.m1) # not significant, AIC 243.7

capture.output(summary(ctenoaccess10.crevice.m1), file='stats/Routput/ctenoaccess10.crevice.doc')

#including all microhabitats
ctenoaccess10.crevice.m2 <- glmer(ctenochaetus_striatus_10_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(ctenoaccess10.crevice.m2) # significant positive impact of increasing minimum opening, AIC 388.1

capture.output(summary(ctenoaccess10.crevice.m2), file='stats/Routput/ctenoaccess10.allcrevices.crevice.doc')

# make separate data sets of no access by Ctenochaetus striatus 10 cm for kernel density plots

crevice.nocteno10 <- subset(micro.data, ctenochaetus_striatus_10_cm==0)

# kernel density plot without Ctenochaetus striatus 10 cm access

crevice.nocteno10.kde <- kde2d(x=crevice.nocteno10$minimum_opening_mm, y= crevice.nocteno10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

par(mar=c(7,8,5,5))
filled.contour(crevice.nocteno10.kde, 
               color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')), 
               plot.axes={axis(1,cex.axis=2)
                 axis(2,cex.axis=2)
               },
               plot.title={title(xlab='Microhabitat width(mm)', cex.lab=2, line=4)
                 title(ylab='Microhabitat depth (mm)', cex.lab=2, line=5)
                 
               })
mtext(side=3, text=substitute(paste('Microhabitats without access by ' ,italic('Ctenochaetus sp.'),  ' (10cm)')), cex=2, line=1, adj=-0.5)


jpeg('graphs/crevice.density.nocteno10.jpg')
image(crevice.nocteno10.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nocteno10.kde,  add=TRUE)
dev.off()

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= ctenochaetus_striatus_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (10 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (10 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno10-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= ctenochaetus_striatus_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of Ctenochaetus striatus (10 cm) access") +
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (10 cm) access")) +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (10 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno10-mincrevice-no300.jpg')


# 9) CTENOCHAETUS STRIATUS 15 CM

#excluding 'open' microhabitats, i.e. > 300 mm
ctenoaccess15.crevice.m1 <- glmer(ctenochaetus_striatus_15_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(ctenoaccess15.crevice.m1) # not significant, AIC 284.8

capture.output(summary(ctenoaccess15.crevice.m1), file='stats/Routput/ctenoaccess15.crevice.doc')


#including all microhabitats
ctenoaccess15.crevice.m2 <- glmer(ctenochaetus_striatus_15_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(ctenoaccess15.crevice.m2) # significant positive impact of increasing minimum opening, AIC 430.5

capture.output(summary(ctenoaccess15.crevice.m2), file='stats/Routput/ctenoaccess15.allcrevices.crevice.doc')


# make separate data sets of no access by Ctenochaetus striatus 15 cm for kernel density plots

crevice.nocteno15 <- subset(micro.data, ctenochaetus_striatus_15_cm==0)

# kernel density plot without Ctenochaetus striatus 15 cm access

crevice.nocteno15.kde <- kde2d(x=crevice.nocteno15$minimum_opening_mm, y= crevice.nocteno15$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nocteno15.jpg')
image(crevice.nocteno15.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nocteno15.kde,  add=TRUE)
dev.off()


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= ctenochaetus_striatus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (15 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (15 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno15-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= ctenochaetus_striatus_15_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (15 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (15 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno15-mincrevice-no300.jpg')





# 10) CTENOCHAETUS STRIATUS 20 CM

#excluding 'open' microhabitats, i.e. > 300 mm
ctenoaccess20.crevice.m1 <- glmer(ctenochaetus_striatus_20_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(ctenoaccess20.crevice.m1) # not significant, AIC 301.0

capture.output(summary(ctenoaccess20.crevice.m1), file='stats/Routput/ctenoaccess20.crevice.doc')

#including all microhabitats
ctenoaccess20.crevice.m2 <- glmer(ctenochaetus_striatus_20_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(ctenoaccess20.crevice.m2) # significant positive impact of minimum opening, AIC 469.6

capture.output(summary(ctenoaccess20.crevice.m2), file='stats/Routput/ctenoaccess20.allcrevices.crevice.doc')

# make separate data sets of no access by Ctenochaetus striatus 20 cm for kernel density plots

crevice.nocteno20 <- subset(micro.data, ctenochaetus_striatus_20_cm==0)

# kernel density plot without Ctenochaetus striatus 20 cm access

crevice.nocteno20.kde <- kde2d(x=crevice.nocteno20$minimum_opening_mm, y= crevice.nocteno20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nocteno20.jpg')
image(crevice.nocteno20.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nocteno20.kde,  add=TRUE)
dev.off()


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= ctenochaetus_striatus_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (20 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (20 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno12-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= ctenochaetus_striatus_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(C.~striatus)*" (20 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(C.~striatus)*" (20 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/cteno20-mincrevice-no300.jpg')





# 11) PARROTFISH 10 CM

#excluding 'open' microhabitats, i.e. > 300 mm
genparrotaccess10.crevice.m1 <- glmer(general_parrotfish_10_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(genparrotaccess10.crevice.m1) # not significant, AIC 276.3

capture.output(summary(genparrotaccess10.crevice.m1), file='stats/Routput/genparrotaccess10.crevice.doc')

#including all microhabitats
genparrotaccess10.crevice.m2 <- glmer(general_parrotfish_10_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(genparrotaccess10.crevice.m2) # significant positive impact of increasing opening AIC 283.2

capture.output(summary(genparrotaccess10.crevice.m2), file='stats/Routput/genparrotaccess10.allcrevices.crevice.doc')


# make separate data sets of no access by parrotfish 10 cm for kernel density plots

crevice.nogenparrot10 <- subset(micro.data, general_parrotfish_10_cm==0)

# kernel density plot without parrotfish 10 cm access

crevice.nogenparrot10.kde <- kde2d(x=crevice.nogenparrot10$minimum_opening_mm, y= crevice.nogenparrot10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

par(mar=c(7,8,5,5))
filled.contour(crevice.nogenparrot10.kde, 
               color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')), 
               plot.axes={axis(1,cex.axis=2)
                 axis(2,cex.axis=2)
               },
               plot.title={title(xlab='Microhabitat width(mm)', cex.lab=2, line=4)
                 title(ylab='Microhabitat depth (mm)', cex.lab=2, line=5)
                 
               })
mtext(side=3, text=substitute(paste('Microhabitats without access by parrotfish (10cm)')), cex=2.5, line=1, adj=-0.5)

jpeg('graphs/crevice.density.nogenparrot10.jpg')
image(crevice.nogenparrot10.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nogenparrot10.kde,  add=TRUE)
dev.off()

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= general_parrotfish_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (10 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (10 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot10-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= general_parrotfish_10_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (10 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (10 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot10-mincrevice-no300.jpg')




# 12) PARROTFISH 20 CM

#excluding 'open' microhabitats, i.e. > 300mm
genparrotaccess20.crevice.m1 <- glmer(general_parrotfish_20_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(genparrotaccess20.crevice.m1) # not significant, AIC 293.8

capture.output(summary(genparrotaccess20.crevice.m1), file='stats/Routput/genparrotaccess20.crevice.doc')

#including all microhabitats
genparrotaccess20.crevice.m2 <- glmer(general_parrotfish_20_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(genparrotaccess20.crevice.m2) # significant positive impact of minimum opening, AIC 386.6

capture.output(summary(genparrotaccess20.crevice.m2), file='stats/Routput/genparrotaccess20.allcrevices.crevice.doc')


# make separate data sets of no access by parrotfish 20 cm for kernel density plots

crevice.nogenparrot20 <- subset(micro.data, general_parrotfish_20_cm==0)

# kernel density plot without parrotfish 20 cm access

crevice.nogenparrot20.kde <- kde2d(x=crevice.nogenparrot20$minimum_opening_mm, y= crevice.nogenparrot20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nogenparrot20.jpg')
image(crevice.nogenparrot20.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nogenparrot20.kde,  add=TRUE)
dev.off()

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= general_parrotfish_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (20 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (20 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot20-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= general_parrotfish_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (20 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (20 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot20-mincrevice-no300.jpg')



# 13) PARROTFISH 30 CM

#excluding 'open' microhabitats, ie. > 300mm
genparrotaccess30.crevice.m1 <- glmer(general_parrotfish_30_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(genparrotaccess30.crevice.m1) # not significant, AIC 273.5

capture.output(summary(genparrotaccess30.crevice.m1), file='stats/Routput/genparrotaccess30.crevice.doc')

#including all microhabitats
genparrotaccess30.crevice.m2 <- glmer(general_parrotfish_30_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(genparrotaccess30.crevice.m2) # significant positive impact of minimum opening, AIC 382.2

capture.output(summary(genparrotaccess30.crevice.m2), file='stats/Routput/genparrotaccess30.allcrevices.crevice.doc')


# make separate data sets of no access by parrotfish 30 cm for kernel density plots

crevice.nogenparrot30 <- subset(micro.data, general_parrotfish_30_cm==0)

# kernel density plot without parrotfish 30 cm access

crevice.nogenparrot30.kde <- kde2d(x=crevice.nogenparrot30$minimum_opening_mm, y= crevice.nogenparrot30$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nogenparrot30.jpg')
image(crevice.nogenparrot30.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nogenparrot30.kde,  add=TRUE)
dev.off()

# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= general_parrotfish_30_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (30 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (30 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot30-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= general_parrotfish_30_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab("Likelihood of parrotfish (30 cm) access") +
  xlab("Minimum crevice opening (mm)") +
  ggtitle("Relationship of access by parrotfish (30 cm) and crevice size") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/genparrot30-mincrevice-no300.jpg')



# 14) MICRORHINUS 20 CM

#excluding 'open' microhabitats, i.e. > 300 mm
microaccess20.crevice.m1 <- glmer(chlorurus_microrhinos_20_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(microaccess20.crevice.m1) # not significant, AIC 292.3

capture.output(summary(microaccess20.crevice.m1), file='stats/Routput/microaccess20.crevice.doc')

#including all microhabitats
microaccess20.crevice.m2 <- glmer(chlorurus_microrhinos_20_cm ~ crevice_depth_mm + minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(microaccess20.crevice.m2) # significant positive impact of minimum opening AIC 326.7

capture.output(summary(microaccess20.crevice.m2), file='stats/Routput/microaccess20.allcrevices.crevice.doc')


# make separate data sets of no access by Chlorurus microrhinos 20 cm for kernel density plots

crevice.nomicro20 <- subset(micro.data, chlorurus_microrhinos_20_cm==0)

# kernel density plot without Chlorurus microrhinos 20 cm access

crevice.nomicro20.kde <- kde2d(x=crevice.nomicro20$minimum_opening_mm, y= crevice.nomicro20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nomicro20.jpg')
image(crevice.nomicro20.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nomicro20.kde,  add=TRUE)
dev.off()


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= chlorurus_microrhinos_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Ch.~microrhinus)*" (20 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Ch.~microrhinus)*" (20 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/Chlmicro20-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= chlorurus_microrhinos_20_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Ch.~microrhinus)*" (20 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Ch.~microrhinus)*" (20 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/Chlmicro20-mincrevice-no300.jpg')


# 15) MICRORHINUS 30 CM

#excluding 'open' microhabitats, i.e. > 300 mm
microaccess30.crevice.m1 <- glmer(chlorurus_microrhinos_30_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data.na)
summary(microaccess30.crevice.m1) # not significant, AIC 235.8

capture.output(summary(microaccess30.crevice.m1), file='stats/Routput/microaccess30.crevice.doc')

#including all microhabitats
microaccess30.crevice.m2 <- glmer(chlorurus_microrhinos_30_cm ~ sdt_crevice_depth_mm + sdt_minimum_opening_mm +  (1| quadrat/subquadrat), family = binomial, data = micro.data)
summary(microaccess30.crevice.m2) # significant positive impact of minimum opening, marginally insignificant negative impact of depth, AIC = 394.5

capture.output(summary(microaccess30.crevice.m2), file='stats/Routput/microaccess30.allcrevices.crevice.doc')

# make separate data sets of no access by Chlorurus microrhinos 30 cm for kernel density plots

crevice.nomicro30 <- subset(micro.data, chlorurus_microrhinos_30_cm==0)

# kernel density plot without Chlorurus microrhinos 30 cm access

crevice.nomicro30.kde <- kde2d(x=crevice.nomicro30$minimum_opening_mm, y= crevice.nomicro30$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))

jpeg('graphs/crevice.density.nomicro30.jpg')
image(crevice.nomicro30.kde,  xlab='minimum crevice width (mm)', ylab= 'crevice depth (mm)')
contour(crevice.nomicro30.kde,  add=TRUE)
dev.off()


# Plots of individual crevice attributes
# including all observations, i.e. 300 for the ones that are completely open
ggplot(micro.data, aes(y= chlorurus_microrhinos_30_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Ch.~microrhinus)*" (30 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Ch.~microrhinus)*" (30 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/Chlmicro30-mincrevice-all.jpg')

# excluding open microhabitats, i.e. >300mm
ggplot(micro.data.na, aes(y= chlorurus_microrhinos_30_cm, x = minimum_opening_mm))+
  geom_point() +
  geom_smooth(se=FALSE)+
  ylab(expression("Likelihood of "*italic(Ch.~microrhinus)*" (30 cm) access")) +
  xlab("Minimum crevice opening (mm)") +
  ggtitle(expression("Relationship of access by "*italic(Ch.~microrhinus)*" (30 cm) and crevice size")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/Chlmicro30-mincrevice-no300.jpg')








# Rainbow density plots fish exclusion ----------------------------------


# Naso lituratus 15 cm
# make data excluding Naso lituratus 15 cm for kdplot
crevice.nonaso15 <- subset(micro.data, naso_lituratus_15_cm==0)

crevice.nonaso15.kde <- kde2d(x=crevice.nonaso15$minimum_opening_mm, y= crevice.nonaso15$crevice_depth_mm, n=200, lims=c(0, 300, 0, 150))

# Naso lituratus 20 cm
# make data excluding Naso lituratus 20 cm for kdplot
crevice.nonaso20 <- subset(micro.data, naso_lituratus_20_cm==0)

crevice.nonaso20.kde <- kde2d(x=crevice.nonaso20$minimum_opening_mm, y= crevice.nonaso20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))




# Naso lituratus 25 cm
# make data excluding Naso lituratus 25 cm for kdplot

crevice.nonaso25 <- subset(micro.data, naso_lituratus_25_cm==0)

crevice.nonaso25.kde <- kde2d(x=crevice.nonaso25$minimum_opening_mm, y= crevice.nonaso25$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))



# Zebrasoma scopas 5 cm
# make data excluding Zebrasoma scopas 5 cm for kdplot
crevice.nozebra5 <- subset(micro.data, zebrasoma_scopas_5_cm==0)

crevice.nozebra5.kde <- kde2d(x=crevice.nozebra5$minimum_opening_mm, y= crevice.nozebra5$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Zebrasoma scopas 10 cm
# make data excluding Zebrasoma scopas 10 cm for kdplot
crevice.nozebra10 <- subset(micro.data, zebrasoma_scopas_10_cm==0)

crevice.nozebra10.kde <- kde2d(x=crevice.nozebra10$minimum_opening_mm, y= crevice.nozebra10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Siganus vulpinus 10 cm
# make data excluding Siganus vulpinus 10 cm for kdplot
crevice.nosiganus10 <- subset(micro.data, siganus_vulpinus_10_cm==0)

# kernel density plot without Siganus vulpinus 10 cm access

crevice.nosiganus10.kde <- kde2d(x=crevice.nosiganus10$minimum_opening_mm, y= crevice.nosiganus10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))




# Siganus vulpinus 15 cm
# make data excluding Siganus vulpinus 15 cm for kdplot
crevice.nosiganus15 <- subset(micro.data, siganus_vulpinus_15_cm==0)

crevice.nosiganus15.kde <- kde2d(x=crevice.nosiganus15$minimum_opening_mm, y= crevice.nosiganus15$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))



# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 10 cm
# make data excluding surgeonfish 10 cm for kdplot
crevice.nocteno10 <- subset(micro.data, ctenochaetus_striatus_10_cm==0)

crevice.nocteno10.kde <- kde2d(x=crevice.nocteno10$minimum_opening_mm, y= crevice.nocteno10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 15 cm
# make data excluding surgeonfish 15 cm for kdplot
crevice.nocteno15 <- subset(micro.data, ctenochaetus_striatus_15_cm==0)

crevice.nocteno15.kde <- kde2d(x=crevice.nocteno15$minimum_opening_mm, y= crevice.nocteno15$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 20 cm
# make data excluding surgeonfish 20 cm for kdplot

crevice.nocteno20 <- subset(micro.data, ctenochaetus_striatus_20_cm==0)

crevice.nocteno20.kde <- kde2d(x=crevice.nocteno20$minimum_opening_mm, y= crevice.nocteno20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))



# Parrotfish 10 cm
# make data excluding general parrotfish 10 cm for kdplot
crevice.nogenparrot10 <- subset(micro.data, general_parrotfish_10_cm==0)

crevice.nogenparrot10.kde <- kde2d(x=crevice.nogenparrot10$minimum_opening_mm, y= crevice.nogenparrot10$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))



# Parrotfish 20 cm
# make data excluding general parrotfish 20 cm for kdplot
crevice.nogenparrot20 <- subset(micro.data, general_parrotfish_20_cm==0)

crevice.nogenparrot20.kde <- kde2d(x=crevice.nogenparrot20$minimum_opening_mm, y= crevice.nogenparrot20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Parrotfish 30 cm
# make data excluding general parrotfish 30 cm for kdplot
crevice.nogenparrot30 <- subset(micro.data, general_parrotfish_30_cm==0)

crevice.nogenparrot30.kde <- kde2d(x=crevice.nogenparrot30$minimum_opening_mm, y= crevice.nogenparrot30$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Chlorurus microrhinos 20 cm 
# make data excluding Chl. microrhinos 20 cm for kdplot
crevice.nomicro20 <- subset(micro.data, chlorurus_microrhinos_20_cm==0)

crevice.nomicro20.kde <- kde2d(x=crevice.nomicro20$minimum_opening_mm, y= crevice.nomicro20$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))



# Chlorurus microrhinos 30 cm
# make data excluding Chl. microrhinos 30 cm for kd plot
crevice.nomicro30 <- subset(micro.data, chlorurus_microrhinos_30_cm==0)

crevice.nomicro30.kde <- kde2d(x=crevice.nomicro30$minimum_opening_mm, y= crevice.nomicro30$crevice_depth_mm, n=50, lims=c(0, 300, 0, 150))


# Make levels to scale all the graphs to the same density measures
lvls <- pretty(range(crevice.nomicro30.kde$z, crevice.nomicro20.kde$z, crevice.nocteno10.kde$z, crevice.nocteno15.kde$z, crevice.nocteno20.kde$z, crevice.nogenparrot10.kde$z, crevice.nogenparrot20.kde$z, crevice.nogenparrot30.kde$z, crevice.nozebra10.kde$z, crevice.nozebra5.kde$z, crevice.nosiganus10.kde$z, crevice.nosiganus15.kde$z, crevice.nonaso15.kde$z, crevice.nonaso20.kde$z, crevice.nonaso25.kde$z), 20)

# Naso lituratus 15 cm
# rainbow density graph
filled.contour(crevice.nonaso15.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Naso lituratus 20 cm
filled.contour(crevice.nonaso20.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })




# Naso lituratus 25 cm
filled.contour(crevice.nonaso25.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# Zebrasoma scopas 5 cm
filled.contour(crevice.nozebra5.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Zebrasoma scopas 10 cm
filled.contour(crevice.nozebra10.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Siganus vulpinus 10 cm
filled.contour(crevice.nosiganus10.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })




# Siganus vulpinus 15 cm
filled.contour(crevice.nosiganus15.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 10 cm
filled.contour(crevice.nocteno10.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 15 cm
filled.contour(crevice.nocteno15.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Surgeonfish: Ctenocheatus striatus / Acanthurus nigrofuscus 20 cm
filled.contour(crevice.nocteno20.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Parrotfish 10 cm
filled.contour(crevice.nogenparrot10.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# Parrotfish 20 cm
filled.contour(crevice.nogenparrot20.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })



# Parrotfish 30 cm
filled.contour(crevice.nogenparrot30.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Chlorurus microrhinos 20 cm 
filled.contour(crevice.nomicro20.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })


# Chlorurus microrhinos 30 cm
# rainbow density graph
filled.contour(crevice.nomicro30.kde, color.palette=colorRampPalette(c('white','blue','yellow', 'red', 'darkred')),
               levels=lvls,
               plot.axes={axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
               },
               plot.title={title(xlab='Microhabitat width (mm)', cex.lab=1.5, line=3)
                 title(ylab='Microhabitat depth (mm)', cex.lab=1.5, line=4)
                 
               })




#===========================================================================================================
#-----------------------------------------------------------------------------------------------------------
# IMPACT OF FISH GRAZING ON LOBOPHORA
# GRAZING PRESSURE MODELS

## model to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)
grazbio.m1<-glmer(lobophora ~ GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro.data)
summary(grazbio.m1) # AIC 481.5
exp(fixef(grazbio.m1)) # 0.987 - per unit increase in grazing pressure, it is 0.99 of the original value 
plogis(fixef(grazbio.m1)) 
linkinv=binomial()$linkinv
linkinv(-0.013) # 0.49675 - not sure which one is the correct back calculation - above seems to make more sense


grazbio.m2 <- glm(lobophora ~  GrazingPressure_biomass, family = binomial, data=micro.data)
summary(grazbio.m2)

#Assumptions
# residuals
plot(grazbio.m1)


# To compare AIC with models not converging in glmer:
#install.packages('glmmTMB')
library(glmmTMB)

grazbio.glmmTMB<-glmmTMB(lobophora ~ GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro.data)
summary(grazbio.glmmTMB) # AIC 481.5


##plot data (without mixed-effects)

jpeg('graphs/Lobo~GrazPressure-combined.jpg')
range(micro.data$GrazingPressure_biomass) ## gives range of GrazingPressure_abundance
xgrazing <- seq(0,110,0.1) ## creates a sequence of values to produce fitted values to
ygrazing <- predict(grazbio.m2, list(GrazingPressure_biomass=xgrazing), type="response") ## creates model for all values xgrazing
plot(micro.data$GrazingPressure_biomass, micro.data$lobophora, pch=16, col="aquamarine3", xlab="Grazing Pressure", ylab="lobophora likelihood") ## plots data points
lines(xgrazing, ygrazing, col="aquamarine4") ## adds trendline
dev.off()
# -> Probability of Lobophora decreases with increasing Grazing pressure as measured by biomass


# Let's do a proper plot
str(micro.data)
newdata <- with(micro.data, expand.grid(GrazingPressure_biomass=seq(min(GrazingPressure_biomass),
                                                                    max(GrazingPressure_biomass), len=100))) # creates x values for GrazingPressure

Xmat <- model.matrix(~GrazingPressure_biomass, newdata) # creates a model matrix


coefs <- fixef(grazbio.m1) # extracts and stores coefficients
fit <- as.vector(coefs %*% t(Xmat)) # multiplies coefficients with new model matrix - creating actual y data
se <- sqrt(diag(Xmat %*% vcov(grazbio.m1) %*% t(Xmat))) # calculate standard error
q <- qt(0.975, df=df.residual(grazbio.m1)) # converges at 1.96 if large data set
q
linkinv=binomial()$linkinv # undoes the link function

newdata <- cbind(newdata, fit=linkinv(fit),
                 lower=linkinv(fit-q*se),
                 upper=linkinv(fit+q*se))

ggplot(newdata, aes(y=fit, x=GrazingPressure_biomass)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3)+
  geom_line() +
  geom_point(data=micro.data, aes(y=lobophora, x=GrazingPressure_biomass))+
  theme_classic()+
  ylab('Lobophora likelihood')+
  xlab('Grazing Pressure')+
  theme(axis.title = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20)) +
  theme(plot.title = element_text(size = 25)) +
  theme(strip.text = element_text(size=20))
ggsave('graphs/Lobophoralikelihood-GrazingPressure-7stepgraph.jpg')

###########

## NO SURGEONFISH
# Import and prepare data
micro_nosurgeon.raw <- read.csv("Microhabitat-data-woAcanthurids.csv", header=T)
micro_nosurgeon.data <-droplevels(subset(micro_nosurgeon.raw, site=='East Sheltered'))
str(micro_nosurgeon.data)
micro_nosurgeon.data$sdt_GrazingPressure_biomass <- scale(micro_nosurgeon.data$GrazingPressure_biomass)

## model without surgeonfish grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)

grazbiom.nosurg.m1<-glmer(lobophora ~ sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosurgeon.data)
summary(grazbiom.nosurg.m1) # doesn't converge
plogis(fixef(grazbiom.nosurg.m1))
## has lower AIC than with all species (481.5 vs 479.4/481.2)

# the model didn't converge. There is more variability between subquadrats within quadrats
# than among quadrats, so I will remove quadrats as a random effet
grazbiom.nosurg.m2<-glmer(lobophora ~ sdt_GrazingPressure_biomass + (1|subquadrat), family=binomial, data=micro_nosurgeon.data)
# the model converges now
summary(grazbiom.nosurg.m2) # AIC is 480.7, compared to 483.1 of all fish included (when quadrat not in model) 
plogis(fixef(grazbiom.nosurg.m2))


# Let's try another type of glmer that usually converges better'

grazbiom.nosurgeon.glmmTMB <- glmmTMB(lobophora ~ sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosurgeon.data)
summary(grazbiom.nosurgeon.glmmTMB) # converges, AIC 481.9


AICsurgeonbio <- exp((481.1-481.9)/2)
AICsurgeonbio ## 0.67 ##  model excluding surgeonfish is 0.39 times as likely as the model including allspecies to minimize information loss
## -> surgeonfish biomass helps explain Lobophora occurrence, but AIC difference to small to make this call
## It doesn't help explain Lobophora occurrence


#########

# NO PARROTFISH
# Import and prepare data
micro_noparrot.raw <- read.csv("Microhabitat-data-woparrot.csv", header=T)
micro_noparrot.data <-droplevels(subset(micro_noparrot.raw, site=="East Sheltered"))

## model without parrotfish grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)

grazbiom.noparrot.m1<-glmer(lobophora ~ GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_noparrot.data)
summary(grazbiom.noparrot.m1)
plogis(fixef(grazbiom.noparrot.m1))
## higher AIC than all species (497.2 vs 481.5)
AICparrotbio <- exp((481.1-497.2)/2)
AICparrotbio ## 0.00031 ## model excluding parrot is 0.00031 times as likely as the model including all species to minimize information loss
## -> parrotfish biomass helps explain Lobophora occurrence


# to compare with models failing to converge using glmmTMB
grazbiom.noparrot.glmmTMB<-glmmTMB(lobophora ~ GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_noparrot.data)
summary(grazbiom.noparrot.glmmTMB) # AIC 496.5




#########

# NO NASOS
# Import and prepare data
micro_nonaso.raw <- read.csv("Microhabitat-data-woNaso.csv", header=T)
micro_nonaso.data <-droplevels(subset(micro_nonaso.raw, site=="East Sheltered"))
micro_nonaso.data$sdt_GrazingPressure_biomass <- scale(micro_nonaso.data$GrazingPressure_biomass)
## model without naso grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)

grazbiom.nonaso.m1<-glmer(lobophora ~  sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nonaso.data)
summary(grazbiom.nonaso.m1)
plogis(fixef(grazbiom.nonaso.m1))
## lower AIC than all species (479.6 vs 481.5)
AICnasobio <- exp((479.5-481.5)/2)
AICnasobio ## 0.37 ## model excluding nasos is 0.37 times as likely as the model including all species to minimize information loss
## -> naso biomass doesn't help explain Lobophora occurrence

# to compare with models failing to converge using glmmTMB
grazbiom.nonaso.glmmTMB<-glmmTMB(lobophora ~  GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nonaso.data)
summary(grazbiom.nonaso.glmmTMB) # AIC 481.5




#########

# NO ZEBRASOMA
# Import and prepare data
micro_nozebrasoma.raw <- read.csv("Microhabitat-data-wozebrasoma.csv", header=T)
micro_nozebrasoma.data <-droplevels(subset(micro_nozebrasoma.raw, site=="East Sheltered"))


## model without zebrasoma grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)
# model has scaling issues, so let's rescale GrazingPressure
micro_nozebrasoma.data$sdt_GrazingPressure_biomass <- scale(micro_nozebrasoma.data$GrazingPressure_biomass)
grazbiom.nozebrasoma.m1<-glmer(lobophora ~  sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nozebrasoma.data)
summary(grazbiom.nozebrasoma.m1)
plogis(fixef(grazbiom.nozebrasoma.m1))
## lower AIC than all species (479.7 vs 481.5)
AICzebrasomabio <- exp((479.7-481.5)/2)
AICzebrasomabio ## 0.41 ## model excluding zebrasomas is 0.41 times as likely as the model including all species to minimize information loss
## -> zebrasoma biomass doesn't help explain Lobophora occurrence

#  model failing to converge
grazbiom.nozebrasoma.glmmTMB<-glmmTMB(lobophora ~  sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nozebrasoma.data)
summary(grazbiom.nozebrasoma.glmmTMB) # 481.5


#########

# NO SIGANUS
# Import and prepare data
micro_nosiganus.raw <- read.csv("Microhabitat-data-wosiganus.csv", header=T)
micro_nosiganus.data <-droplevels(subset(micro_nosiganus.raw, site=="East Sheltered"))

## model without siganus grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)
micro_nosiganus.data$sdt_GrazingPressure_biomass <- scale(micro_nosiganus.data$GrazingPressure_biomass)
grazbiom.nosiganus.m1<-glmer(lobophora ~   sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosiganus.data)
summary(grazbiom.nosiganus.m1)
plogis(fixef(grazbiom.nosiganus.m1))
## lower AIC than all species (479.7 vs 481.5)
AICsiganusbio <- exp((479.7-481.5)/2)
AICsiganusbio ## 0.41 ## model excluding siganuss is 0.41 times as likely as the model including all species to minimize information loss
## -> siganus biomass doesn't help explain Lobophora occurrence

# to compare to models not converging:
grazbiom.nosiganus.glmmTMB<-glmmTMB(lobophora ~   sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosiganus.data)
summary(grazbiom.nosiganus.glmmTMB) # 481.5


#########

# NO ZEBRASOMA NOR SIGANUS
# Import and prepare data
micro_nosigzeb.raw <- read.csv("Microhabitat-data-wosigzeb.csv", header=T)
micro_nosigzeb.data <-droplevels(subset(micro_nosigzeb.raw, site=="East Sheltered"))

## model without siganus and zebrasoma grazing to calculate mean probabilities of Lobophora occurence by grazing pressure(biomass)
#model had scaling issues, let's rescale
micro_nosigzeb.data$sdt_GrazingPressure_biomass <- scale(micro_nosigzeb.data$GrazingPressure_biomass)
grazbiom.nosigzeb.m1<-glmer(lobophora ~  sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosigzeb.data)
summary(grazbiom.nosigzeb.m1)
plogis(fixef(grazbiom.nosigzeb.m1))
## lower AIC than all species (479.7 vs 481.5)
AICsigzebbio <- exp((479.7-481.5)/2)
AICsigzebbio ## 0.41 ## model excluding siganus and zebrasoma is 0.41 times as likely as the model including all species to minimize information loss
## -> siganus and zebrasoma biomass combined doesn't help explain Lobophora occurrence

# because model was not converging:
grazbiom.nosigzeb.glmmTMB<-glmmTMB(lobophora ~  sdt_GrazingPressure_biomass + (1|quadrat/subquadrat), family=binomial, data=micro_nosigzeb.data)
summary(grazbiom.nosigzeb.glmmTMB) # 481.4





#=================================================================================

# CHISQUARE TO ANALYSE DEVIANCE OF LIKELIHOOD OF LOBOPHORA FROM EXPECTED PER CREVICE

#---------------------------------------------------------------------------------

# 1st Step: Make crevice categories

str(micro.data)

# show me the values minimum_opening can have
width_unique <- unique(micro.data$minimum_opening_mm)
width_unique <- as.numeric(width_unique)
sort(width_unique)
# [1]   3   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
#[19]  22  23  24  25  26  27  28  29  30  32  33  35  36  40  43  44  45  53
#[37]  54  55  64  65  71  73  79 300

# categories of 20 mm increments
# 1 - 20; 21 - 40; 41 - 60; 61 - 80; 300 --> 5 categories



# show me the values depth can have
depth_unique <- unique(micro.data$crevice_depth_mm)
depth_unique <- as.numeric(depth_unique)
sort(depth_unique)
#[1]   0   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
#[19]  20  21  22  23  24  25  26  27  28  29  30  31  32  33  35  36  38  39
#[37]  40  42  43  50  51  56  59  73  74  75  80  87 100 125

# categories of 20 mm increments
# 1 - 20; 21 - 40; 41 - 60; 61 - 80; 81 - 100; 121 - 140; 0 --> 7 categories

# number of overall categories: 5 * 7 = 35
# 1) width 1 - 20 & depth 1 - 20
# 2) width 1 - 20 & depth 21 - 40
# 3) width 1 - 20 & depth 41 - 60
# 4) width 1 - 20 & depth 61 - 80
# 5) width 1 - 20 & depth 81 - 100
# 6) width 1 - 20 & depth 121 - 140
# 7) width 1 - 20 & depth 0
# 8) width 21 - 40 & depth 1 - 20
# 9) width 21 - 40 & depth 21 - 40
# 10) width 21 - 40 & depth 41 - 60
# 11) width 21 - 40 & depth 61 - 80
# 12) width 21 - 40 & depth 81 - 100
# 13) width 21 - 40 & depth 121 - 140
# 14) width 21 - 40 & depth 0
# 15) width 41 - 60 & depth 1 - 20
# 16) width 41 - 60 & depth 21 - 40
# 17) width 41 - 60 & depth 41 - 60
# 18) width 41 - 60 & depth 61 - 80
# 19) width 41 - 60 & depth 81 - 100
# 20) width 41 - 60 & depth 121 - 140
# 21) width 41 - 60 & depth 0
# 22) width 61 - 80 & depth 1 - 20
# 23) width 61 - 80 & depth 21 - 40
# 24) width 61 - 80 & depth 41 - 60
# 25) width 61 - 80 & depth 61 - 80
# 26) width 61 - 80 & depth 81 - 100
# 27) width 61 - 80 & depth 121 - 140
# 28) width 61 - 80 & depth 0
# 29) width 300 & depth 1 - 20
# 30) width 300 & depth 21 - 40
# 31) width 300 & depth 41 - 60
# 32) width 300 & depth 61 - 80
# 33) width 300 & depth 81 - 100
# 34) width 300 & depth 121 - 140
# 35) width 300 & depth 0


# Let's try to print this as categories

micro.data$category <- ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '1',
                              ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '2',
                                     ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '3', 
                                            ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '4',
                                                   ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '5',
                                                          ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '6',
                                                                 ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm == 0, '7',
                                                                        ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '8',
                                                                               ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '9',
                                                                                      ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '10', 
                                                                                             ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '11',
                                                                                                    ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '12',
                                                                                                           ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '13',
                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm == 0, '14',
                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '15',
                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '16',
                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '17', 
                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '18',
                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '19',
                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '20',
                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm == 0, '21',
                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '22',
                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '23',
                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '24', 
                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '25',
                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '26',
                                                                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '27',
                                                                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm == 0, '28',
                                                                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '29',
                                                                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '30',
                                                                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '31', 
                                                                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '32',
                                                                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '33',
                                                                                                                                                                                                                                                              ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '34',
                                                                                                                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm == 0, '35', NA)))))))))))))))))))))))))))))))))))
# Let's have a look at it
micro.data$category
str(micro.data)
# what categories do realistically exist?
category_unique <- unique(micro.data$category)
category_unique <- as.numeric(category_unique)
sort(category_unique) # only 22 of 35 possible categories do exist

#-----------------------------------------------------------------------------------------
# 2nd Step: Let's get the proportion of crevices in each category

# there are 450 crevices in total

# let's count the categories
num_crevice_cat <- micro.data %>% 
  group_by(category) %>%
  tally()
sort(num_crevice_cat)
num_cat <-data.frame(num_crevice_cat)
num_cat

# data was entered into excel sheet crevice_categories.csv



#-------------------------------------------------------------------------------------------
# 3rd Step: Proportion of crevices with Lobophora growing in them

# overall number of crevices still 450
micro.data$lobophora
num_crevice_lobo <- micro.data %>% 
  group_by(lobophora) %>%
  tally()
data.frame(num_crevice_lobo)

#View(num_crevice_lobo) # only works in RStudio
# 341 crevices don't have Lobophora
# 109 crevices do have Lobophora

str(num_crevice_lobo)
num_crevice_lobo$lobophora <- as.factor(num_crevice_lobo$lobophora)
prop_lobo <- num_crevice_lobo$n/450
prop_lobo
# 0.76 don't have Lobophora
# 0.24 do have Lobophora


#-------------------------------------------------------------------------------------------
# 4th Step: Proportion of each crevice with Lobophora
str(micro.data)
micro.data$category <- as.factor(micro.data$category)

# overall number of crevices 450
# How many crevices of each category have Lobophora?
num_crevice_cat_lobo <- micro.data %>% 
  group_by(category,lobophora) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_lobo) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet



#-------------------------------------------------------------------------------------------
# 5th Step: Calculate Chi-Square with ChiSquare = (observed number - expected number)^2/ expected number
# Has to be done using full numbers, not proportions
# Calculated in excel sheet 

# import this excel sheet
crevice.categories <- read.csv('crevice_categories.csv')



#-------------------------------------------------------------------------------------------
# NOW THE SAME THING FOR FISH ACCESS BY EACH SHAPE
# Step1 and 2 are the same

#-------------------------------------------------------------------------------------------
# 3rd Step: Proportion of crevices with access by fish species

# NASO 15 CM
# overall number of crevices still 450
micro.data$naso_lituratus_15_cm
num_crevice_naso15 <- micro.data %>% 
  group_by(naso_lituratus_15_cm) %>%
  tally()
data.frame(num_crevice_naso15)

# Naso 15 doesn't have access to 50 crevices
# Naso 15 does have access to 400

prop_naso15 <- num_crevice_naso15$n/450
prop_naso15
# 0.11 aren't accessible
# 0.89 are accessible



# NASO 20 CM
# overall number of crevices still 450
micro.data$naso_lituratus_20_cm
num_crevice_naso20 <- micro.data %>% 
  group_by(naso_lituratus_20_cm) %>%
  tally()
data.frame(num_crevice_naso20)

# Naso 20 doesn't have access to 63 crevices
# Naso 20 does have access to 387

prop_naso20 <- num_crevice_naso20$n/450
prop_naso20
# 0.14 aren't accessible
# 0.86 are accessible


# NASO 25 CM
# overall number of crevices still 450
micro.data$naso_lituratus_25_cm
num_crevice_naso25 <- micro.data %>% 
  group_by(naso_lituratus_25_cm) %>%
  tally()
data.frame(num_crevice_naso25)

# Naso 25 doesn't have access to 79 crevices
# Naso 25 does have access to 371

prop_naso25 <- num_crevice_naso25$n/450
prop_naso25
# 0.18 aren't accessible
# 0.82 are accessible




# ZEBRASOMA SCOPAS 5 CM
# overall number of crevices still 450
micro.data$zebrasoma_scopas_5_cm
num_crevice_zebra5 <- micro.data %>% 
  group_by(zebrasoma_scopas_5_cm) %>%
  tally()
data.frame(num_crevice_zebra5)

# Zebrasoma 5 doesn't have access to 46 crevices
# Zebrasoma 5 does have access to 404

prop_zebra5 <- num_crevice_zebra5$n/450
prop_zebra5
# 0.10 aren't accessible
# 0.90 are accessible



# ZEBRASOMA SCOPAS 10 CM
# overall number of crevices still 450
micro.data$zebrasoma_scopas_10_cm
num_crevice_zebra10 <- micro.data %>% 
  group_by(zebrasoma_scopas_10_cm) %>%
  tally()
data.frame(num_crevice_zebra10)

# Zebrasoma 10 doesn't have access to 49 crevices
# Zebrasoma 10 does have access to 401

prop_zebra10 <- num_crevice_zebra10$n/450
prop_zebra10
# 0.11 aren't accessible
# 0.89 are accessible




# SIGANUS VULPINS 10 CM
# overall number of crevices still 450
micro.data$siganus_vulpinus_10_cm
num_crevice_siganus10 <- micro.data %>% 
  group_by(siganus_vulpinus_10_cm) %>%
  tally()
data.frame(num_crevice_siganus10)

# Siganus vulpinus 10 doesn't have access to 5 crevices
# Siganus vulpinus 10 does have access to 445

prop_siganus10 <- num_crevice_siganus10$n/450
prop_siganus10
# 0.01 aren't accessible
# 0.99 are accessible




# SIGANUS VULPINS 15 CM
# overall number of crevices still 450
micro.data$siganus_vulpinus_15_cm
num_crevice_siganus15 <- micro.data %>% 
  group_by(siganus_vulpinus_15_cm) %>%
  tally()
data.frame(num_crevice_siganus15)

# Siganus vulpinus 15 doesn't have access to 7 crevices
# Siganus vulpinus 15 does have access to 443

prop_siganus15 <- num_crevice_siganus15$n/450
prop_siganus15
# 0.02 aren't accessible
# 0.98 are accessible



# CTENOCHAETUS STRIATUS 10 CM
# overall number of crevices still 450
micro.data$ctenochaetus_striatus_10_cm
num_crevice_cteno10 <- micro.data %>% 
  group_by(ctenochaetus_striatus_10_cm) %>%
  tally()
data.frame(num_crevice_cteno10)

# Ctenochaetus striatus 10 doesn't have access to 75 crevices
# Ctenochaetus striatus 10 does have access to 375

prop_cteno10 <- num_crevice_cteno10$n/450
prop_cteno10
# 0.17 aren't accessible
# 0.83 are accessible


# CTENOCHAETUS STRIATUS 15 CM
# overall number of crevices still 450
micro.data$ctenochaetus_striatus_15_cm
num_crevice_cteno15 <- micro.data %>% 
  group_by(ctenochaetus_striatus_15_cm) %>%
  tally()
data.frame(num_crevice_cteno15)

# Ctenochaetus striatus 15 doesn't have access to 104 crevices
# Ctenochaetus striatus 15 does have access to 346

prop_cteno15 <- num_crevice_cteno15$n/450
prop_cteno15
# 0.23 aren't accessible
# 0.76 are accessible




# CTENOCHAETUS STRIATUS 20 CM
# overall number of crevices still 450
micro.data$ctenochaetus_striatus_20_cm
num_crevice_cteno20 <- micro.data %>% 
  group_by(ctenochaetus_striatus_20_cm) %>%
  tally()
data.frame(num_crevice_cteno20)

# Ctenochaetus striatus 20 doesn't have access to 145 crevices
# Ctenochaetus striatus 20 does have access to 305

prop_cteno20 <- num_crevice_cteno20$n/450
prop_cteno20
# 0.32 aren't accessible
# 0.68 are accessible



# GENERAL PARROTFISH 10 CM
# overall number of crevices still 450
micro.data$general_parrotfish_10_cm
num_crevice_parrot10 <- micro.data %>% 
  group_by(general_parrotfish_10_cm) %>%
  tally()
data.frame(num_crevice_parrot10)

# General parrotfish 10 cm doesn't have access to 77 crevices
# General parrotfish 10 cm does have access to 373

prop_parrot10 <- num_crevice_parrot10$n/450
prop_parrot10
# 0.17 aren't accessible
# 0.83 are accessible



# GENERAL PARROTFISH 20 CM
# overall number of crevices still 450
micro.data$general_parrotfish_20_cm
num_crevice_parrot20 <- micro.data %>% 
  group_by(general_parrotfish_20_cm) %>%
  tally()
data.frame(num_crevice_parrot20)

# General parrotfish 20 cm doesn't have access to 131 crevices
# General parrotfish 20 cm does have access to 319

prop_parrot20 <- num_crevice_parrot20$n/450
prop_parrot20
# 0.29 aren't accessible
# 0.71 are accessible



# GENERAL PARROTFISH 30 CM
# overall number of crevices still 450
micro.data$general_parrotfish_30_cm
num_crevice_parrot30 <- micro.data %>% 
  group_by(general_parrotfish_30_cm) %>%
  tally()
data.frame(num_crevice_parrot30)

# General parrotfish 30 cm doesn't have access to 163 crevices
# General parrotfish 30 cm does have access to 287

prop_parrot30 <- num_crevice_parrot30$n/450
prop_parrot30
# 0.36 aren't accessible
# 0.64 are accessible



# CHLORURUS MICRORHINUS 20 CM
# overall number of crevices still 450
micro.data$chlorurus_microrhinos_20_cm
num_crevice_chlorurus20 <- micro.data %>% 
  group_by(chlorurus_microrhinos_20_cm) %>%
  tally()
data.frame(num_crevice_chlorurus20)

# Chlorurus microrhinos 20 cm doesn't have access to 108 crevices
# Chlorurus microrhinos 20 cm does have access to 342

prop_chlorurus20 <- num_crevice_chlorurus20$n/450
prop_chlorurus20
# 0.24 aren't accessible
# 0.76 are accessible



# CHLORURUS MICRORHINUS 30 CM
# overall number of crevices still 450
micro.data$chlorurus_microrhinos_30_cm
num_crevice_chlorurus30 <- micro.data %>% 
  group_by(chlorurus_microrhinos_30_cm) %>%
  tally()
data.frame(num_crevice_chlorurus30)

# Chlorurus microrhinos 30 cm doesn't have access to 189 crevices
# Chlorurus microrhinos 30 cm does have access to 261

prop_chlorurus30 <- num_crevice_chlorurus30$n/450
prop_chlorurus30
# 0.42 aren't accessible
# 0.58 are accessible



#-------------------------------------------------------------------------------------------
# 4th Step: Proportion of each crevice with access by each fish group (i.e. surgeon, parrot, naso, zebrasoma, siganus)
str(micro.data)

# NASOS ACCESS
# overall number of crevices 450
# How many crevices of each category are not accessible by Nasos?
num_crevice_cat_naso <- micro.data %>% 
  group_by(category,naso_access) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_naso) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet


# ZEBRASOMA ACCESS
# overall number of crevices 450
# How many crevices of each category are not accessible by Nasos?
num_crevice_cat_zebra <- micro.data %>% 
  group_by(category,zebrasoma_access) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_zebra) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet



# SIGANUS ACCESS
# overall number of crevices 450
# How many crevices of each category are not accessible by Nasos?
num_crevice_cat_siganus <- micro.data %>% 
  group_by(category,siganus_access) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_siganus) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet



# SURGEON ACCESS
# overall number of crevices 450
# How many crevices of each category are not accessible by Nasos?
num_crevice_cat_surgeon <- micro.data %>% 
  group_by(category,surgeon_access) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_surgeon) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet




# PARROT ACCESS
# overall number of crevices 450
# How many crevices of each category are not accessible by Nasos?
num_crevice_cat_parrot <- micro.data %>% 
  group_by(category,parrot_access) %>%
  tally()
#View(num_crevice_cat_lobo) # only works in RStudio
data.frame(num_crevice_cat_parrot) # as an alternative to print it into the console

# proportions have been entered from this count into an excel sheet



#=======================================================================================

# Let's look at the relationship between lobophora Chi and fish group Chi

crevice.cat.noNA <- na.exclude(crevice.categories)
histogram(crevice.cat.noNA$ChiSquare.lobo, breaks=100)
mean(crevice.cat.noNA$ChiSquare.lobo)
var(crevice.cat.noNA$ChiSquare.lobo) # overdispersed data
crevice.cat.noNA$log.ChiSquare.lobo <- log(crevice.cat.noNA$ChiSquare.lobo)
histogram(crevice.cat.noNA$log.ChiSquare.lobo, breaks=10) # this is getting somewhat closer to normal distribution

# Naso

formula <- y~x
ggplot(crevice.cat.noNA, aes(y=ChiSquare.lobo, x=ChiSquare.nonaso))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  stat_poly_eq(aes(label=paste(..eq.label.., ..rr.label.., sep='~~~')),
               label.x.npc='right', label.y.npc=0.15,
               formula= formula, parse= TRUE, size=6)+
  ylab("ChiSquare Lobophora") +
  xlab("ChiSquare Naso sp.") +
  ggtitle("Correlation of unexpected crevice behaviour") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_relationship_Naso.jpg')

# gaussian model
Chirel.nonaso.glm1 <- glm(ChiSquare.lobo ~ChiSquare.nonaso, data=crevice.cat.noNA, family=gaussian)
autoplot(Chirel.nonaso.glm1, which=1:6, ncol=2, label.size=3) # not good!

# log-gaussian model
Chirel.nonaso.glm2 <- glm(log.ChiSquare.lobo ~ChiSquare.nonaso, data=crevice.cat.noNA, family=gaussian)
autoplot(Chirel.nonaso.glm2, which=1:6, ncol=2, label.size=3) # looks better but still not good

# gamma distribution
Chirel.nonaso.glm3 <- glm(ChiSquare.lobo ~ChiSquare.nonaso, data=crevice.cat.noNA, family=Gamma(link=inverse))
autoplot(Chirel.nonaso.glm3, which=1:6, ncol=2, label.size=3) # one cooks D is ridiculously high

# tweedie distribution
library(statmod)

Chirel.nonaso.glm4 <- glm(ChiSquare.lobo ~ChiSquare.nonaso, data=crevice.cat.noNA, family=tweedie)
autoplot(Chirel.nonaso.glm4, which=1:6, ncol=2, label.size=3) 

# negative binomial model
Chirel.nonaso.nb <- glm.nb(ChiSquare.lobo ~ ChiSquare.nonaso, data = crevice.cat.noNA) # tried many different distributions, neg binomial is the best fit even though still bad
summary(Chirel.nonaso.nb)

# compare model fits
AIC(Chirel.nonaso, Chirel.nonaso.glm1,Chirel.nonaso.glm2,Chirel.nonaso.glm3,Chirel.nonaso.glm4)


# best model seems to be binomial because of assumptions, gamma is a better fit but Cook's D is much much higher


plot(Chirel.nonaso.nb, which=1)
qqnorm(resid(Chirel.nonaso.nb))
qqline(resid(Chirel.nonaso.nb))


# Zebrasoma
ggplot(crevice.cat.noNA, aes(y=ChiSquare.lobo, x=ChiSquare.nozebra))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  stat_poly_eq(aes(label=paste(..eq.label.., ..rr.label.., sep='~~~')),
               label.x.npc='right', label.y.npc=0.15,
               formula= formula, parse= TRUE, size=6)+
  ylab("ChiSquare Lobophora") +
  xlab("ChiSquare Zebrasoma sp.") +
  ggtitle("Correlation of unexpected crevice behaviour") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_relationship_Zebra.jpg')

Chirel.nozebra <- glm.nb(ChiSquare.lobo ~ ChiSquare.nozebra, data = crevice.cat.noNA)
summary(Chirel.nozebra)

plot(Chirel.nozebra, which=1)
qqnorm(resid(Chirel.nozebra))
qqline(resid(Chirel.nozebra))

# Siganus
ggplot(crevice.cat.noNA, aes(y=ChiSquare.lobo, x=ChiSquare.nosiganus))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  stat_poly_eq(aes(label=paste(..eq.label.., ..rr.label.., sep='~~~')),
               label.x.npc='right', label.y.npc=0.15,
               formula= formula, parse= TRUE, size=6)+
  ylab("ChiSquare Lobophora") +
  xlab("ChiSquare Siganus sp.") +
  ggtitle("Correlation of unexpected crevice behaviour") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_relationship_Siganus.jpg')



Chirel.nosiganus <- glm.nb(ChiSquare.lobo ~ ChiSquare.nosiganus, data = crevice.cat.noNA)
summary(Chirel.nosiganus)

plot(Chirel.nosiganus, which=1)
qqnorm(resid(Chirel.nosiganus))
qqline(resid(Chirel.nosiganus))

# Surgeon
ggplot(crevice.cat.noNA, aes(y=ChiSquare.lobo, x=ChiSquare.nosurgeon))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  stat_poly_eq(aes(label=paste(..eq.label.., ..rr.label.., sep='~~~')),
               label.x.npc='right', label.y.npc=0.15,
               formula= formula, parse= TRUE, size=6)+
  ylab("ChiSquare Lobophora") +
  xlab("ChiSquare Acanthurus sp. & Ctenochaetus sp.") +
  ggtitle("Correlation of unexpected crevice behaviour") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_relationship_Surgeon.jpg')

Chirel.nosurgeon <- glm.nb(ChiSquare.lobo ~ ChiSquare.nosurgeon, data = crevice.cat.noNA)
summary(Chirel.nosurgeon)

plot(Chirel.nosurgeon, which=1)
qqnorm(resid(Chirel.nosurgeon))
qqline(resid(Chirel.nosurgeon))

# Parrot
ggplot(crevice.cat.noNA, aes(y=ChiSquare.lobo, x=ChiSquare.noparrot))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  stat_poly_eq(aes(label=paste(..eq.label.., ..rr.label.., sep='~~~')),
               label.x.npc='right', label.y.npc=0.15,
               formula= formula, parse= TRUE, size=6)+
  ylab("ChiSquare Lobophora") +
  xlab("ChiSquare Scarus sp. and Chlorurus sp.") +
  ggtitle("Correlation of unexpected crevice behaviour") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_relationship_parrot.jpg')

Chirel.noparrot <- glm.nb(ChiSquare.lobo ~ ChiSquare.noparrot, data = crevice.cat.noNA)
summary(Chirel.noparrot)


#======================================================================================
# Looking at data to see whether other crevice categories are better suited

# Excluding crevice categories with fewer than 5 categories in them
crevice.cat.larger5 <- subset(crevice.cat.noNA, num.category.crevices >=5)
head(crevice.cat.larger5)
crevice.cat.larger5$num.category.crevices
# only 7 crevice categories if all below 5 crevices are excluded


#--------------------------------------------------------------------------------------------
# Make crevice categories, split width into 1-10 and 11-20mm, everything else stays the same
micro.data$category2 <- ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '1',
                               ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '2',
                                      ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '3', 
                                             ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '4',
                                                    ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '5',
                                                           ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '6',
                                                                  ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm == 0, '7',
                                                                         ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '8',
                                                                                ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '9',
                                                                                       ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '10', 
                                                                                              ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '11',
                                                                                                     ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '12',
                                                                                                            ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '13',
                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm == 0, '14',
                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '15',
                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '16',
                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '17', 
                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '18',
                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '19',
                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '20',
                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm == 0, '21',
                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '22',
                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '23',
                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '24', 
                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '25',
                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '26',
                                                                                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '27',
                                                                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm == 0, '28',
                                                                                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '29',
                                                                                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '30',
                                                                                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '31', 
                                                                                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '32',
                                                                                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '33',
                                                                                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '34',
                                                                                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm == 0, '35',
                                                                                                                                                                                                                                                                             ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '36',
                                                                                                                                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '37',
                                                                                                                                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '38', 
                                                                                                                                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '39',
                                                                                                                                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100,'40',
                                                                                                                                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '41',
                                                                                                                                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm == 0, '42', NA))))))))))))))))))))))))))))))))))))))))))

# Let's have a look at it
micro.data$category2
str(micro.data)
# what categories do realistically exist?
category_unique <- unique(micro.data$category)
category_unique <- as.numeric(category_unique)
sort(category_unique) # only 22 of 35 possible categories do exist

#-----------------------------------------------------------------------------------------
# 2nd Step: Let's get the proportion of crevices in each category

# there are 450 crevices in total

# let's count the categories
num_crevice_cat2 <- micro.data %>% 
  group_by(category2) %>%
  tally()
sort(num_crevice_cat2)
num_cat2 <-data.frame(num_crevice_cat2)
num_cat2
test <- subset(num_crevice_cat2, n >=5)
test


micro.data$category2 <- ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '1',
                               ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '2',
                                      ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '3', 
                                             ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '4',
                                                    ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '5',
                                                           ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '6',
                                                                  ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm == 0, '7',
                                                                         ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '8',
                                                                                ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '9',
                                                                                       ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '10', 
                                                                                              ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '11',
                                                                                                     ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '12',
                                                                                                            ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '13',
                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm == 0, '14',
                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '15',
                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '16',
                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '17', 
                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '18',
                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '19',
                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '20',
                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm == 0, '21',
                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '22',
                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '23',
                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '24', 
                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '25',
                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '26',
                                                                                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '27',
                                                                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm == 0, '28',
                                                                                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '29',
                                                                                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '30',
                                                                                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '31', 
                                                                                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '32',
                                                                                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100, '33',
                                                                                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '34',
                                                                                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm == 0, '35',
                                                                                                                                                                                                                                                                             ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 20, '36',
                                                                                                                                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 40, '37',
                                                                                                                                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 60, '38', 
                                                                                                                                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 80, '39',
                                                                                                                                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 100,'40',
                                                                                                                                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 140, '41',
                                                                                                                                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm == 0, '42', NA))))))))))))))))))))))))))))))))))))))))))

# Let's have a look at it
micro.data$category2
str(micro.data)
# what categories do realistically exist?
category_unique <- unique(micro.data$category)
category_unique <- as.numeric(category_unique)
sort(category_unique) # only 22 of 35 possible categories do exist

# let's count the categories
num_crevice_cat2 <- micro.data %>% 
  group_by(category2) %>%
  tally()
sort(num_crevice_cat2)
num_cat2 <-data.frame(num_crevice_cat2)
num_cat2
#which ones have more than 5 crevices per category?
test <- subset(num_crevice_cat2, n >=5)
test # 9 crevice categories




#------------------------------------------------------------------------------------
# Let's try splitting everything up into smaller categories of 1-10, 11-20 etc.

micro.data$category2 <- ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '1',
                               ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '2',
                                      ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '3',
                                             ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '4',
                                                    ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '5', 
                                                           ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '6', 
                                                                  ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '7',
                                                                         ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '8',
                                                                                ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '9',
                                                                                       ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '10',
                                                                                              ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '11',
                                                                                                     ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '12',
                                                                                                            ifelse(micro.data$minimum_opening_mm >= 1 & micro.data$minimum_opening_mm <= 10 & micro.data$crevice_depth_mm == 0, '13',
                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '14',
                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '15',
                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '16',
                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '17',
                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '18',
                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '19',
                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '20',
                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '21',
                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '22',
                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '23',
                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '24',
                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '25',
                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 11 & micro.data$minimum_opening_mm <= 20 & micro.data$crevice_depth_mm == 0, '26',
                                                                                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '27',
                                                                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '28',
                                                                                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '29',
                                                                                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '30',
                                                                                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '31', 
                                                                                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '32', 
                                                                                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '33',
                                                                                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '34',
                                                                                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '35', NA)))))))))))))))))))))))))))))))))))
micro.data$category3 <-  ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '36',
                                ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '37',
                                       ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '38',
                                              ifelse(micro.data$minimum_opening_mm >= 21 & micro.data$minimum_opening_mm <= 30 & micro.data$crevice_depth_mm == 0, '39',
                                                     ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '40',
                                                            ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '41',
                                                                   ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '42',
                                                                          ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '43',
                                                                                 ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '44', 
                                                                                        ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '45', 
                                                                                               ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '46',
                                                                                                      ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '47',
                                                                                                             ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '48',
                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '49',
                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '50',
                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '51',
                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 31 & micro.data$minimum_opening_mm <= 40 & micro.data$crevice_depth_mm == 0, '52',
                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '53',
                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '54',
                                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '55',
                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '56',
                                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '57', 
                                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '58', 
                                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '59',
                                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '60',
                                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '61',
                                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '62',
                                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '63',
                                                                                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '64',
                                                                                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 41 & micro.data$minimum_opening_mm <= 50 & micro.data$crevice_depth_mm == 0, '65',
                                                                                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '66',
                                                                                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '67',
                                                                                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '68',
                                                                                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '69',
                                                                                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '70', NA)))))))))))))))))))))))))))))))))))
micro.data$category4 <-  ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '71', 
                                ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '72',
                                       ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '73',
                                              ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '74',
                                                     ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '75',
                                                            ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '76',
                                                                   ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '77',
                                                                          ifelse(micro.data$minimum_opening_mm >= 51 & micro.data$minimum_opening_mm <= 60 & micro.data$crevice_depth_mm == 0, '78',
                                                                                 ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '79',
                                                                                        ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '80',
                                                                                               ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '81',
                                                                                                      ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '82',
                                                                                                             ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '83', 
                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '84', 
                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '85',
                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '86',
                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '87',
                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '88',
                                                                                                                                                       ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '89',
                                                                                                                                                              ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '90',
                                                                                                                                                                     ifelse(micro.data$minimum_opening_mm >= 61 & micro.data$minimum_opening_mm <= 70 & micro.data$crevice_depth_mm == 0, '91',
                                                                                                                                                                            ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '92',
                                                                                                                                                                                   ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '93',
                                                                                                                                                                                          ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '94',
                                                                                                                                                                                                 ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '95',
                                                                                                                                                                                                        ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50, '96', 
                                                                                                                                                                                                               ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '97', 
                                                                                                                                                                                                                      ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '98',
                                                                                                                                                                                                                             ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '99',
                                                                                                                                                                                                                                    ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '100',
                                                                                                                                                                                                                                           ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '101',
                                                                                                                                                                                                                                                  ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130, '102',
                                                                                                                                                                                                                                                         ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '103',
                                                                                                                                                                                                                                                                ifelse(micro.data$minimum_opening_mm >= 71 & micro.data$minimum_opening_mm <= 80 & micro.data$crevice_depth_mm == 0, '104',
                                                                                                                                                                                                                                                                       ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 1 & micro.data$crevice_depth_mm <= 10, '105', NA)))))))))))))))))))))))))))))))))))
micro.data$category5 <-   ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 11 & micro.data$crevice_depth_mm <= 20, '106',
                                 ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 21 & micro.data$crevice_depth_mm <= 30, '107', 
                                        ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 31 & micro.data$crevice_depth_mm <= 40, '108',
                                               ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 41 & micro.data$crevice_depth_mm <= 50,'109',
                                                      ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 51 & micro.data$crevice_depth_mm <= 60, '110', 
                                                             ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 61 & micro.data$crevice_depth_mm <= 70, '111',
                                                                    ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 71 & micro.data$crevice_depth_mm <= 80, '112',
                                                                           ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 81 & micro.data$crevice_depth_mm <= 90, '113', 
                                                                                  ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 91 & micro.data$crevice_depth_mm <= 100, '114',
                                                                                         ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 121 & micro.data$crevice_depth_mm <= 130,'115',
                                                                                                ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm >= 131 & micro.data$crevice_depth_mm <= 140, '116',
                                                                                                       ifelse(micro.data$minimum_opening_mm == 300 & micro.data$crevice_depth_mm == 0, '117', NA))))))))))))
# for some reason I can't run it as one, so ran it as 4 individual ifelse loops and combined hereafter
str(micro.data)
# Let's have a look at it
micro.data$category2
micro.data$category3
micro.data$category4
micro.data$category5
micro.data$catcomb <- micro.data$category2
micro.data$catcomb[!is.na(micro.data$category3)] = micro.data$category3[!is.na(micro.data$category3)]
micro.data$catcomb[!is.na(micro.data$category4)] = micro.data$category4[!is.na(micro.data$category4)]
micro.data$catcomb[!is.na(micro.data$category5)] = micro.data$category5[!is.na(micro.data$category5)]
micro.data$catcomb


# what categories do realistically exist?
catcomb_unique <- unique(micro.data$catcomb)
catcomb_unique <- as.numeric(catcomb_unique)
sort(catcomb_unique) # only 46 out of 117 actually exist

# let's count the categories
num_crevice_catcomb <- micro.data %>% 
  group_by(catcomb) %>%
  tally()
sort(num_crevice_catcomb)
num_catcomb <-data.frame(num_crevice_catcomb)
num_catcomb
#which ones have more than 5 crevices per category?
catcomb.larger5 <- subset(num_crevice_catcomb, n >=5)
catcomb.larger5 # 11 crevice categories


#---------------------------------------------------------------------------

# Let's plot ChiSquare over crevice size category

plot1 <-ggplot(data=crevice.cat.noNA, aes(x=Category, y = ChiSquare.lobo))+
  geom_col()+
  ylab("ChiSquare Lobophora") +
  xlab("Category") +
  ggtitle("Unexpectedly high or low Lobophora likelihood per crevice") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))
ggsave('graphs/ChiSquare_lobo_Category_column.jpg')

# I could plot them above each other using the following method and assigning 'plot1' to the plot above
plot2 <- ggplot(data=crevice.cat.noNA, aes(x=Category, y = proportion.with.lobo))+
  geom_col()+
  ylab("Proportion of crevices with Lobophora") +
  xlab("Category") +
  ggtitle("Proportion of crevices with Lobophora") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))

# binding them requires libraries 'ggplot2', 'grid' and 'dplyr'

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = 'last'))


# Let's see if I can make a graph that shows 'surprise' Chi and whether it is higher or lower than expected
# this is not possible with ggplot

jpeg('graphs/ChiLobo-with-proportionsfound-allcrev.jpg')
par(mar = c(5,5,2,5))
with(crevice.cat.noNA, plot(Category, ChiSquare.lobo, type = 'h', col='black', lwd = 15,
                            ylab='ChiSquare Lobophora',
                            ylim = c(0,20)))
par(new=T)
with(crevice.cat.noNA, plot(Category, proportion.with.lobo, type = 'p',  pch= 16, col= 'red3', bg='red3',  axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = c(0,1)))
axis(side=4)
mtext(side = 4, line = 3, 'Proportion of Lobophora in crevice')

par(new=T)
with(crevice.cat.noNA, plot(Category, overall.exp.prop.lobo, type= 'l', axes=F, col='red3', lwd=5, xlab=NA, ylab=NA, cex=1.2, ylim = c(0,1)))
axis(side=4)
mtext(side = 4, col='red', line = 3, 'Proportion of Lobophora in crevice')
legend('topleft', 
       legend=c('ChiSquare Lobophora', 'proportion found', 'proportion expected'), lty=c(1,0), pch=c(NA, 16), col=c('black', 'red3', 'red3'))
dev.off()

# do the same thing only with categories that  have more than 5 crevices in each of them

jpeg('graphs/ChiLobo-with-proportionsfound-crevlarger5.jpg')
par(mar = c(5,5,2,5))
with(crevice.cat.larger5, plot(Category, ChiSquare.lobo, type = 'h', col='black', lwd = 15,
                               ylab='ChiSquare Lobophora',
                               ylim = c(0,20)))
par(new=T)
with(crevice.cat.larger5, plot(Category, proportion.with.lobo, type = 'p',  pch= 16, col= 'red3', bg='red3',  axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = c(0,1)))
axis(side=4)
mtext(side = 4, line = 3, 'Proportion of Lobophora in crevice')

par(new=T)
with(crevice.cat.larger5, plot(Category, overall.exp.prop.lobo, type= 'l', axes=F, col='red3', lwd=5, xlab=NA, ylab=NA, cex=1.2, ylim = c(0,1)))
axis(side=4)
mtext(side = 4, col='red', line = 3, 'Proportion of Lobophora in crevice')
legend('topleft', 
       legend=c('ChiSquare Lobophora', 'proportion found', 'proportion expected'), lty=c(1,0), pch=c(NA, 16), col=c('black', 'red3', 'red3'))
dev.off()
