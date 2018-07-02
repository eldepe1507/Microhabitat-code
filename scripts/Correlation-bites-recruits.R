# Analyse correlation of fish bites with removal of Lobophora recruits
# Lots of data manipulation first
# Last edited 21/06/2018

# set-up
rm(list=ls())
library(dplyr) # data wrangling
library(ggplot2) # graphs
library(reshape2) # changing between wide and long format
library(data.table) # further data manipulation, including setnames()
library(MuMIn) # get r squared for mixed models
library(lme4) # glmer


# Including all data:
#==================================================================

# Import data
getwd() # check working directory
cor.data.raw <- read.csv("data/Correlation-bites-raw.csv") # import data set
str(cor.data.raw)
cor.data.raw$tile.number <- as.factor(cor.data.raw$tile.number)
cor.data.raw$Individual <- as.factor(cor.data.raw$Individual)
recruit.data <- read.csv("data/Correlation-recruitnumber.csv")
str(recruit.data)
recruit.data$tile.number <- as.factor(recruit.data$tile.number)



# Data manipulation
# Summarize by conditions: tile.number, Location and Species
bites.data <- cor.data.raw %>%
  group_by(tile.number, Location, Species) %>%
  summarise(sumbites = sum(number.bites))
str(bites.data)

# Check whether the manipulation worked - compare to manual sum
write.csv(bites.data, file="data/checkmanipulation.csv") # worked!


# Now I want to merge it with other data set
# binding by column Location
data.comb <- left_join(recruit.data, bites.data, by=c("tile.number", "Location"))
str(data.comb)
data.comb
# check data set again
write.csv(data.comb, file="data/checkmerging.csv") # looks ok
levels(data.comb$Species) # ignores NA, good

# Reshaping to wide format - making #bites by species x' the columns
data.comb.wide <- dcast(data.comb, colour.in.videos + Treatment + Crevice.kind +tile.number + Location + number.LOB.before + number.otherMA.before + number.unknown.before+ number.Lob.after + number.otherMA.after + number.unknown.after+ number.LOB.removed + number.otherMA.removed + number.unknown.removed ~ Species, value.var="sumbites")
data.comb.wide


# Now I need to remove NA column and turn NAs into 0s below fish species
# remove NA column
str(data.comb.wide) # NA is in last column, number changes, so we can't just say column number 21 etc.
data.comb.wide <- data.comb.wide[,!grepl("NA", names(data.comb.wide))] # NA column removed
str(data.comb.wide)


# Turn NA into 0
# We first need to get rid of space in column names
setnames(data.comb.wide, old=c("Acanthurus nigrofuscus", "Chlorurus sordidus", "Ctenochaetus binotatus", "Ctenochaetus striatus", "Naso lituratus", "Pomacentrus bankanensis", "Scarus niger", "unident wrasse", "Zebrasoma scopas"), new= c("Bites.Acanthurus.nigrofuscus", "Bites.Chlorurus.sordidus", "Bites.Ctenochaetus.binotatus", "Bites.Ctenochaetus.striatus", "Bites.Naso.lituratus", "Bites.Pomacentrus.bankanensis", "Bites.Scarus.niger", "Bites.unident.wrasse", "Bites.Zebrasoma.scopas"))
head(data.comb.wide)
str(data.comb.wide)

# in order to turn NAs into 0s I need to turn columns into numerical/double
data.comb.wide$Bites.Acanthurus.nigrofuscus <- as.double(data.comb.wide$Bites.Acanthurus.nigrofuscus)
data.comb.wide$Bites.Chlorurus.sordidus <- as.double(data.comb.wide$Bites.Chlorurus.sordidus)
data.comb.wide$Bites.Ctenochaetus.binotatus <- as.double(data.comb.wide$Bites.Ctenochaetus.binotatus)
data.comb.wide$Bites.Ctenochaetus.striatus <- as.double(data.comb.wide$Bites.Ctenochaetus.striatus)
data.comb.wide$Bites.Naso.lituratus <- as.double(data.comb.wide$Bites.Naso.lituratus)
data.comb.wide$Bites.Pomacentrus.bankanensis <- as.double(data.comb.wide$Bites.Pomacentrus.bankanensis)
data.comb.wide$Bites.Scarus.niger <- as.double(data.comb.wide$Bites.Scarus.niger)
data.comb.wide$Bites.unident.wrasse <- as.double(data.comb.wide$Bites.unident.wrasse)
data.comb.wide$Bites.Zebrasoma.scopas <- as.double(data.comb.wide$Bites.Zebrasoma.scopas)
str(data.comb.wide)


data.comb.wide <- data.comb.wide %>% 
  mutate(Bites.Acanthurus.nigrofuscus = if_else(is.na(Bites.Acanthurus.nigrofuscus), 0, Bites.Acanthurus.nigrofuscus),
         Bites.Chlorurus.sordidus = if_else(is.na(Bites.Chlorurus.sordidus), 0, Bites.Chlorurus.sordidus),
         Bites.Ctenochaetus.binotatus = if_else(is.na(Bites.Ctenochaetus.binotatus), 0, Bites.Ctenochaetus.binotatus),
         Bites.Ctenochaetus.striatus = if_else(is.na(Bites.Ctenochaetus.striatus), 0, Bites.Ctenochaetus.striatus),
         Bites.Naso.lituratus = if_else(is.na(Bites.Naso.lituratus), 0, Bites.Naso.lituratus),
         Bites.Pomacentrus.bankanensis = if_else(is.na(Bites.Pomacentrus.bankanensis), 0, Bites.Pomacentrus.bankanensis),
         Bites.Scarus.niger = if_else(is.na(Bites.Scarus.niger), 0, Bites.Scarus.niger),
         Bites.unident.wrasse = if_else(is.na(Bites.unident.wrasse), 0, Bites.unident.wrasse),
         Bites.Zebrasoma.scopas = if_else(is.na(Bites.Zebrasoma.scopas), 0, Bites.Zebrasoma.scopas)
  )
head(data.comb.wide)
str(data.comb.wide)


# Let's also make a long format data frame
data.comb.long <- melt(data.comb.wide, id.vars=c("colour.in.videos", "Treatment", "Crevice.kind", "tile.number", "Location", "number.LOB.before", "number.otherMA.before", "number.unknown.before", "number.Lob.after", "number.otherMA.after", "number.unknown.after", "number.LOB.removed", "number.otherMA.removed", "number.unknown.removed"), measure.vars=c("Bites.Acanthurus.nigrofuscus","Bites.Chlorurus.sordidus", "Bites.Ctenochaetus.binotatus", "Bites.Ctenochaetus.striatus", "Bites.Naso.lituratus", "Bites.Pomacentrus.bankanensis", "Bites.Scarus.niger", "Bites.unident.wrasse", "Bites.Zebrasoma.scopas"),variable.name='Species', value.name='Bites')
str(data.comb.long)

# Let's turn this into a excel table csv format
write.csv(data.comb.wide, file="data/Correlation-combined.csv")

# Remove all rows without Lobophora recruits before videos
# long data
str(data.comb.long)
data.comb.long2 <- data.comb.long[!data.comb.long$number.LOB.before==0,]
data.comb.long2
#wide data
str(data.comb.wide)
data.comb.wide2 <- data.comb.wide[!data.comb.wide$number.LOB.before==0,]


# Quick look at some plots
ggplot(data.comb.long2, aes(y= number.LOB.removed, x= Bites))+
  geom_point()+
  facet_wrap(~Species)+
  geom_smooth(method=lm)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  xlab("Bites")+
  ylab("Recruits removed")
ggsave("graphs/bitesrecruit-correlation.jpg")


ggplot(data.comb.long2, aes(y= number.LOB.removed, x= Bites))+
  geom_point()+
  geom_smooth(method=lm)


# Let's see if I can add the equation and R^2 to each plot
# so far only overall equation, not by facets

p <- ggplot(data.comb.long2, aes(y= number.LOB.removed, x= Bites))+
  geom_point()+
  facet_wrap(~Species)+
  geom_smooth(method=lm)
p

linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 3));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)  
  }
  as.character(as.expression(eq));              
}

fo = number.LOB.removed ~ Bites
p1 = p + annotate("text", x = 10, y = 5, label = linear(lm(fo, data.comb.long2)), colour="black", size = 5, parse=TRUE)
p1 # so far this is the same equation in all of the graphs

# What if we look at percentage removal of recruits?
data.comb.long2$percent.recruit.removal <- (data.comb.long2$number.LOB.removed/data.comb.long2$number.LOB.before)*100
str(data.comb.long2)

# let's plot the percentage removal of recruits
ggplot(data.comb.long2, aes(y= percent.recruit.removal, x= Bites))+
  geom_point()+
  facet_wrap(~Species)+
  geom_smooth(method=lm)

ggplot(data.comb.long2, aes(y= percent.recruit.removal, x= Bites))+
  geom_point()+
  geom_smooth(method=lm)


#------------------------------------------------

# Run simple regression models - glmer

# look at data
str(data.comb.wide2)

# Naso unicornis
Naso.cor.mod <- glmer(number.LOB.removed ~ Bites.Naso.lituratus + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Naso.cor.mod) # no significantl correlation
r.squaredGLMM(Naso.cor.mod) # doesn't seem to be reliable - all get same value...

# Ctenochaetus striatus
Ctstriatus.cor.mod <- glmer(number.LOB.removed ~ Bites.Ctenochaetus.striatus + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Ctstriatus.cor.mod) # significant correlation
r.squaredGLMM(Ctstriatus.cor.mod) # doesn't seem to be reliable

# Zebrasoma scopas
Zebrasoma.cor.mod <- glmer(number.LOB.removed ~ Bites.Zebrasoma.scopas + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Zebrasoma.cor.mod) # significant correlation
r.squaredGLMM(Ctstriatus.cor.mod) # doesn't seem to be reliable

# Chlorurus sordidus
Chlsordidus.cor.mod <- glmer(number.LOB.removed ~ Bites.Chlorurus.sordidus + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Chlsordidus.cor.mod) # no significantl correlation
r.squaredGLMM(Chlsordidus.cor.mod) # doesn't seem to be reliable

# Scarus niger
Scniger.cor.mod <- glmer(number.LOB.removed ~ Bites.Scarus.niger + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Scniger.cor.mod) # significant correlation
r.squaredGLMM(Scniger.cor.mod) # doesn't seem to be reliable


# Acanthurus nigrofuscus
Anigrofuscus.cor.mod <- glmer(number.LOB.removed ~ Bites.Acanthurus.nigrofuscus + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Anigrofuscus.cor.mod) # no significant correlation
r.squaredGLMM(Anigrofuscus.cor.mod) # doesn't seem to be reliable

# Ctenochaetus binotatus
Ctbinotatus.cor.mod <- glmer(number.LOB.removed ~ Bites.Ctenochaetus.binotatus + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Ctbinotatus.cor.mod) # converging issues, result looks odd, no significant correlation
r.squaredGLMM(Ctbinotatus.cor.mod) # doesn't seem to be reliable

# Pomacentrus bankanensis
Pombank.cor.mod <- glmer(number.LOB.removed ~ Bites.Pomacentrus.bankanensis + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(Pombank.cor.mod) # marginally insignificant correlation
r.squaredGLMM(Pombank.cor.mod) # doesn't seem to be reliable

# unidentified wrasse
wrasse.cor.mod <- glmer(number.LOB.removed ~ Bites.unident.wrasse + (1|colour.in.videos/tile.number), family = poisson, data = data.comb.wide2)
summary(wrasse.cor.mod) # doesn't converge, output odd, not significant
r.squaredGLMM(wrasse.cor.mod) # doesn't seem to be reliable





#------------------------------------------------
# Only using locations where one individual species has bitten:
#===================================================================

# import data set
# bites
cor.data.raw1 <- read.csv("data/Lobophora-recruitment-singlespecies-correlation.csv")
str(cor.data.raw1)
# recruits
cor.recruit.data.raw1 <- read.csv('data/Correlation-recruitnumber.csv')

# turn data I want into wide format
cor.data.wide1 <- dcast(cor.data.raw1, Unique.location ~ Species, value.var="number.bites", fun.aggregate=sum)
head(cor.data.wide1)
str(cor.data.wide1)
levels(cor.data.wide1$Unique.location) # looks all good

# drop column NA
cor.data.wide2 <- cor.data.wide1[, !names(cor.data.wide1) %in% "NA"]
str(cor.data.wide2)


# print new column that indicates how many species take bites on one location
cor.data.wide2$count <- rowSums(cor.data.wide2[-1] > 0)
head(cor.data.wide2) # worked!

# Make new dataset only using rows with count=1
cor.data.wide3 <- subset(cor.data.wide2, cor.data.wide2[,"count"] ==1)
cor.data.wide3 # worked, to check:
#write.csv(cor.data.wide3, file="data/check-selectingcolumnswithcount1.csv")
# data set is much smaller than before! only about 200 rows
str(cor.data.wide3)

# I can now remove the count column again
cor.data.wide4 <- cor.data.wide3[,!names(cor.data.wide3) %in% "count"]
str(cor.data.wide4)

# Now turn data set back into long format before merging with recruit data
cor.data.long1 <- melt(cor.data.wide4, id.vars = "Unique.location",
                       measur.vars=c("Acanthurus nigrofuscus", "Chlorurus sordidus", "Ctenochaetus binotatus", 
                                     "Ctenochaetus striatus", "Naso lituratus", "Pomacentrus bankanensis",
                                     "Scarus niger", "unident wrasse", "Zebrasoma scopas"),
                       variable.name= "Species", value.name="Bites")
head(cor.data.long1)

# Let's merge the data
str(cor.recruit.data.raw1)
cor.data.comb1 <- left_join(cor.data.long1, cor.recruit.data.raw1, by="Unique.location")
head(cor.data.comb1)
write.csv(cor.data.comb1, file="data/check-mergeofcorrelationdatawithonly1species.csv")




# Remove all rows without Lobophora recruits before videos
str(cor.data.comb1)
cor.data.comb2 <- cor.data.comb1[!cor.data.comb1$number.LOB.before==0,]
str(cor.data.comb2)
write.csv(cor.data.comb2, file="data/check-selectingrowswhenLOBbefore==0indataonly1species.csv")


# What's the overall number of bites per species?
sumofbites <- aggregate(cor.data.comb2$Bites, by=list(cor.data.comb2$Species), FUN=sum)
sumofbites

sumofbites <- 
  cor.data.comb2 %>% 
  group_by(Species) %>%
  summarise(., Bites= sum(Bites))
write.csv(sumofbites, file="data/sumbites-1speciesbiting.csv")


#######

# PLOTS
# show trend of recruit removal over fish bites - seperated by species
ggplot(cor.data.comb2, aes(y= number.LOB.removed, x= Bites))+
  geom_point()+
  facet_wrap(~Species)+
  geom_smooth(method=lm, se=FALSE)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15))+
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
  scale_y_continuous(breaks=c(0,1,2))+
  xlab("Bites")+
  ylab("# Recruits removed")
#ggsave("graphs/Correlation-bites-recruits_only1speciesbitingperlocation.jpg")
# Candidate species for removal are Acanthurus nigrofuscus, Zebrasoma scopas, Scarus niger
# No indication that Ctenochaetus striatus can remove Lobophora recruits!!
# No removal of recruits from C. sordidus may be due to overall low number of bites
# Overall number of bites should probably be included in this graph!

# Average bites per species
ggplot(cor.data.comb2, aes(y=Bites, x= Species))+
  geom_col()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 15, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size= 15))+
  xlab("Species")+
  ylab("Total # bites")+
  theme(plot.margin=unit(c(1,1,1.5,1.2), "cm"))
#ggsave("graphs/SumBitesSpecies-only1speciesbiting.jpg")


#-----------------------------------------------------------------

# Run simple regression models - glmer

# look at data
str(cor.data.comb2)
head(cor.data.comb2)

# Turn into wide format to run each one individually
cor.data.comb2.wide <- dcast(cor.data.comb2, Treatment + colour.in.videos + tile.number +  Location + Unique.location +  number.LOB.removed ~ Species, value.var="Bites", fun.aggregate=sum)
head(cor.data.comb2.wide)

# Turn column names into something usable in R
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Acanthurus nigrofuscus"] <- "Bites.Acanthurus.nigrofuscus"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Chlorurus sordidus"] <- "Bites.Chlorurus.sordidus"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Ctenochaetus binotatus"] <- "Bites.Ctenochaetus.binotatus"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Ctenochaetus striatus"] <- "Bites.Ctenochaetus.striatus"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Naso lituratus"] <- "Bites.Naso.lituratus"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Pomacentrus bankanensis"] <- "Bites.Pomacentrus.bankanensis"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Scarus niger"] <- "Bites.Scarus.niger"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "unident wrasse"] <- "Bites.unidentified.wrasse"
colnames(cor.data.comb2.wide)[colnames(cor.data.comb2.wide) == "Zebrasoma scopas"] <- "Bites.Zebrasoma.scopas"



# Let's run some models - GLMER
str(cor.data.comb2.wide)
cor.data.comb2.wide$tile.number <- as.factor(cor.data.comb2.wide$tile.number)

# Acanthurus nigrofuscus
Anigrofuscus.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Acanthurus.nigrofuscus + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Anigrofuscus.cor.mod2) # significant correlation

# Chlorurus sordidus
Chlsordidus.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Chlorurus.sordidus + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Chlsordidus.cor.mod2) # doesn't work....

# Ctenochaetus binotatus
Ctbinotatus.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Ctenochaetus.binotatus + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Ctbinotatus.cor.mod2) # doesn't work

# Ctenochaetus striatus
Ctstriatus.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Ctenochaetus.striatus + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Ctstriatus.cor.mod2) # doesn't work

# Naso unicornis
Naso.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Naso.lituratus + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Naso.cor.mod2) # doesn't work

# Pomacentrus bankanensis
Pombank.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Pomacentrus.bankanensis + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Pombank.cor.mod2) # rank deficient - doesn#t work

# Scarus niger
Scniger.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Scarus.niger + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Scniger.cor.mod2) # significant correlation

# unidentified wrasse
wrasse.cor.mod2 <- glmer(number.LOB.removed ~ Bites.unidentified.wrasse + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(wrasse.cor.mod2) # doesn't work


# Zebrasoma scopas
Zebrasoma.cor.mod2 <- glmer(number.LOB.removed ~ Bites.Zebrasoma.scopas + (1|colour.in.videos/tile.number/Location), family = poisson, data = cor.data.comb2.wide)
summary(Zebrasoma.cor.mod2) # not significant











  