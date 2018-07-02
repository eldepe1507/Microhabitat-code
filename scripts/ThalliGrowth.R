# Analysis of adult Lobophora growth in manipulated tiles
# Last edited 14 June 2018

rm(list=ls())
library(dplyr) # data wrangling
library(ggplot2) # graphs
library(reshape2) # changing between wide and long format
library(data.table) # further data manipulation, including setnames()
library(MuMIn) # get r squared for mixed models
library(lubridate) # good for working with Dates
library(Rmisc) # allows you to use summarySE

# Import data
getwd() # check working directory
growth.data1 <- read.csv("data/ThalliMeasurements.csv") # import data set
str(growth.data1)
growth.data1$width.mm <- as.integer(growth.data1$width.mm)
growth.data1$Tile.Number <- as.factor(growth.data1$Tile.Number)
growth.data1$Date <- dmy(growth.data1$Date)
growth.data1$Crevice.number <- as.factor(growth.data1$Crevice.number)


# summarize data
growth.data2 <- growth.data1 %>%
  group_by(Date, Treatment, Crevice.tile) %>%
  filter(!is.na(length.mm)) %>%
  filter(!is.na(width.mm)) %>%
  summarise(avg.length = mean(length.mm),
            sd.length = sd(length.mm),
            avg.width = mean(width.mm),
            sd.width = sd(width.mm))

str(growth.data2)


growth.data3 <- growth.data1 %>%
  filter(!is.na(length.mm)) %>%
  filter(!is.na(width.mm))
summary.length.data <- summarySE(growth.data3, 
                          measurevar="length.mm", 
                          groupvars=c("Date", "Treatment", "Crevice.tile"))

head(summary.length.data)
summary.width.data <- summarySE(growth.data3,
                                measurevar="width.mm",
                                groupvars=c("Date", "Treatment", "Crevice.tile"))


head(summary.width.data)

# graph length
dp <- position_dodge(1.5) # errorbars overlap, so I'm moving them slightly apart
par(mfrow=c(2,1))
ggplot(summary.length.data, aes(y=length.mm, x= Date, col= Crevice.tile))+
  geom_point(position=dp)+
  geom_errorbar(aes(ymin=length.mm-se, ymax=length.mm+se), width=4, position=dp)+
  geom_line(position=dp)+
  facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1))+
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  scale_y_continuous(limits=c(0,10))+
  xlab("Date")+
  ylab("length (mm)")


# graph width
ggplot(summary.width.data, aes(y=width.mm, x= Date, col= Crevice.tile))+
  geom_point(position=dp)+
  geom_errorbar(aes(ymin=width.mm-se, ymax=width.mm+se), width=4, position=dp)+
  geom_line(position=dp)+
  facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1))+
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  scale_y_continuous(limits=c(0,10))+
  xlab("Date")+
  ylab("width (mm)")



