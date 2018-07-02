# Check if there is a correlation between turf height and Lobophora recruit numbers
# Date edited: 04 May 2018

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
recruit.turf.data1 <- read.csv("data/Turf-Lobrecruit-correlation.csv") # import data set


# Let's check data and prepare it
str(recruit.turf.data1)
recruit.turf.data1$Tile.number <- as.factor(recruit.turf.data1$Tile.number) # number doesn't conceil information so it will be a factor

# make column with recruits per cm2 data
recruit.turf.data1$LOB.target.percm2.inclsides <- recruit.turf.data1$LOB.targetcrevice/recruit.turf.data1$Area.targetcrevice.cm2.inclsides
head(recruit.turf.data1)
recruit.turf.data1$LOB.overall.percm2.inclsides <- recruit.turf.data1$LOB.overall/recruit.turf.data1$Area.overall.cm2.inclsides
head(recruit.turf.data1)

# Reshape data and remove NAs
cor.data1 <- recruit.turf.data1 %>% na.omit()
cor.data.long <- cor.data1 %>% gather(key=Turf.measurement, value=Turf.height, turf.meas.1, turf.meas.2, turf.meas.3)
str(cor.data.long)

# maybe better for graph to just make average beforehand
cor.data1$Turf.height.avg <- rowMeans(cor.data1[,c('turf.meas.1', 'turf.meas.2', 'turf.meas.3')])

# Graph to show potential correlation between turf height and Lobophora recruit numbers
ggplot(cor.data1, aes(y= Turf.height.avg, x= LOB.target.percm2.inclsides, fill= crevice.kind))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~Treatment)
  
ggplot(cor.data1, aes(y=Turf.height.avg, x=crevice.kind))+
  geom_col()+
  facet_wrap(~Treatment)

