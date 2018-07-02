# Analyse where fish bites on the tiles
# Analyse whether there is a difference depending on Lobophora being there or not
# Last edited 21.06.2018


# set-up
rm(list=ls())
#install.packages('Rmisc')
library(Rmisc) #summarySE
library(dplyr) # data wrangling
library(ggplot2) # graphs
library(reshape2) # changing between wide and long format
library(data.table) # further data manipulation, including setnames()
library(MuMIn) # get r squared for mixed models


# Import data
getwd() # check working directory
bites.data1 <- read.csv("data/Fishbites-MHtiles-Lob.vs.noLob.csv") # import data set
str(bites.data1)
bites.data1$Individual <- as.factor(bites.data1$Individual)



#----------------------------

# Create column for bites per cm^2 for each category
# crown
bites.data1$Bites.cm2.crown <- bites.data1$Bites.taken.on.crown/bites.data1$Area.crown.cm2
# Target crevice
bites.data1$Bites.cm2.targetcrevice <- bites.data1$Bites.taken.in.target.crevice/bites.data1$Area.target.cm2
# Other crevice
bites.data1$Bites.cm2.othercrevice <- bites.data1$Bites.on.other.crevice/bites.data1$Area.other.cm2
# Open crevice
bites.data1$Bites.cm2.opencrevice <- bites.data1$Bites.taken.open/bites.data1$Area.open.cm2
# Outside
bites.data1$Bites.cm2.outside <- bites.data1$Bites.on.outside.tile/bites.data1$Area.outside.cm2


# Create column for bites per cm^2 per minutes for each category
# Crown
bites.data1$Bites.cm2.min.crown <- bites.data1$Bites.cm2.crown/bites.data1$video.length.min
# Target crevice
bites.data1$Bites.cm2.min.targetcrevice <- bites.data1$Bites.cm2.targetcrevice/bites.data1$video.length.min
# Other crevice
bites.data1$Bites.cm2.min.othercrevice <- bites.data1$Bites.cm2.othercrevice/bites.data1$video.length.min
# Open crevice
bites.data1$Bites.cm2.min.opencrevice <- bites.data1$Bites.cm2.opencrevice/bites.data1$video.length.min
# Outside
bites.data1$Bites.cm2.min.outside <- bites.data1$Bites.cm2.outside/bites.data1$video.length.min
str(bites.data1)


bites.data1$crown <- bites.data1$Bites.cm2.min.crown 
bites.data1$targetcrevice <- bites.data1$Bites.cm2.min.targetcrevice 
bites.data1$othercrevice <- bites.data1$Bites.cm2.min.othercrevice
bites.data1$opencrevice <- bites.data1$Bites.cm2.min.opencrevice 
bites.data1$outside <- bites.data1$Bites.cm2.min.outside 
str(bites.data1)
  
# Before graphing we need to turn it into long format 
bites.data1.long <- melt(bites.data1, id.vars=c("ColourCode", "Treatment", "Tile", "Genus.species"),
                         measure.vars=c("crown", "targetcrevice", 
                                        "othercrevice", "opencrevice",
                                        "outside"),
                         variable.name="Location.of.bites", value.name="Bites.cm2.min")
str(bites.data1.long)
head(bites.data1.long)

# summarize data for graphs
bites.data2 <- bites.data1.long %>%
  group_by(ColourCode, Treatment, Tile, Genus.species, Location.of.bites) %>%
  summarise(sumbites.cm2.min = sum(Bites.cm2.min),
            meanbites.cm2.min = mean(Bites.cm2.min),
            sdbites.cm2.min = sd(Bites.cm2.min))
  
str(bites.data2)



# Graph it
# Bites per cm2 per minute split up by species, Treatment and crevice size
ggplot(bites.data2, aes(y=sumbites.cm2.min, x=Location.of.bites, fill=Genus.species))+
  geom_col(position="dodge")+
  facet_grid(Tile~ Treatment)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1))+
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")
ggsave("graphs/fishbitescrevicetile-LobvsnoLob.jpg")

# Bites per cm2 per minute split up by species and crevice size
ggplot(bites.data2, aes(y=sumbites.cm2.min, x=Location.of.bites, fill=Genus.species))+
  geom_col(position="dodge")+
  facet_grid(~ Tile)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")
ggsave("graphs/fishbitescrevicetile-comparecrevicekinds.jpg")


# Let's look at it from a different angle: split by fish species
# and coloured by crevice type
ggplot(bites.data2, aes(y=sumbites.cm2.min, x=Location.of.bites, fill=Tile))+
  geom_col(position="dodge")+
  facet_grid(~ Genus.species)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")


# Let's look at it from a different angle: split by fish species and Tile
# and coloured by crevice location
ggplot(bites.data2, aes(y=sumbites.cm2.min, x=Tile, fill=Location.of.bites))+
  geom_col(position="fill")+
  facet_grid(~ Genus.species)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Tile kind")+
  ylab("Proportion of bites on each tile")+
  theme(legend.text= element_text(size=15)) +
  theme(legend.title = element_text(size= 20))+
  scale_fill_discrete(name="Location of bites")
ggsave("graphs/comparisonspeciesbites-crevicekinds.jpg")


#----------------------------------------------------------

# Because I don't have any parrotfish bites on the tiles, I am going to include data 
# from other video set-up

# Start with importing and merging data


# Import data
cor.data.raw <- read.csv("data/Correlation-bites-raw.csv") # import data set
str(cor.data.raw)
cor.data.raw$tile.number <- as.factor(cor.data.raw$tile.number)
cor.data.raw$Individual <- as.factor(cor.data.raw$Individual)
# import translation of Location(a1, b2, etc.) to crevice location (target crevice, other crevice, crown etc.)
translation <- read.csv("data/translation-location.csv")
str(translation)
translation$tile.number <- as.factor(translation$tile.number)


# Merge the two data sets - I now know how many bites by which species
# in which crevice location (target etc.)
bites.comb1 <- full_join(cor.data.raw, translation, by=c("tile.number", "Location"))
head(bites.comb1)
write.csv(bites.comb1, file="data/check-translationbind.csv")
# Looks good

# to merge it with the other table, first turn it into wide format
bites.comb1.wide <- dcast(bites.comb1, Date + Reef + Video_code + colour + Species + Individual + fishlength.cm + tile.kind + tile.number + Location + length.video.min ~ Location.kind, value.var="number.bites", fun.aggregate=sum)
head(bites.comb1.wide)

# Next step is to merge it with the fish bite data from 1st set of videos
bites.data.all1 <- read.csv("data/Fishbites-MHtiles-Lob.vs.noLob.csv") # import data set again with different name
str(bites.data.all1)
bites.data.all1$Individual <- as.factor(bites.data.all1$Individual)
head(bites.data.all1)



#Modify datasets so they can be merged
# rename columns (Bites.taken.on... etc.) so they match between both datasets
colnames(bites.data.all1)[colnames(bites.data.all1) == "Bites.taken.on.crown"] <- "crown"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Bites.taken.in.target.crevice"] <- "targetcrevice"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Bites.on.other.crevice"] <- "othercrevice"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Bites.taken.open"] <- "opencrevice"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Bites.on.outside.tile"] <- "outside"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Location"] <- "Reef"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Tile"] <- "tile.kind"
colnames(bites.data.all1)[colnames(bites.data.all1) == "Genus.species"] <- "Genus.species"


colnames(bites.comb1.wide)[colnames(bites.comb1.wide) == "Species"] <- "Genus.species"

colnames(bites.comb1.wide)[colnames(bites.comb1.wide) == "fishlength.cm"] <- "Body.length"
colnames(bites.comb1.wide)[colnames(bites.comb1.wide) == "length.video.min"] <- "video.length.min"

bites.data.all1$tile.number <- as.character(bites.data.all1$tile.number)
# Let's merge

bites.comb2 <- full_join(bites.data.all1, bites.comb1.wide, by=c("Reef", "Genus.species", "Video_code",
                                                             "tile.number", "Individual", "tile.kind",
                                                             "video.length.min", "Body.length", "crown", "outside", 
                                                             "targetcrevice", "opencrevice", 
                                                             "othercrevice"))



head(bites.comb2)

#Let's check this data set
#write.csv(bites.comb2, file="data/checkjoining-bothvideosets-raw.csv") 

# What's the average body length for each species?
fish.length <- bites.comb2 %>%
  group_by(Genus.species) %>%
  summarise(avg.bodylength = mean(Body.length))
fish.length

# File looked relatively good but small adjustments had to be made,
# So let's reimport it and work with the new version

#----------------------------
# If I want to include parrotfish from other video set-up
# I can run everything from here on - because new dataset is created
bites.data.comb <- read.csv("data/checkjoining-bothvideosets.csv")




# Create column for bites per cm^2 for each category
# crown
bites.data.comb$Bites.cm2.crown <- bites.data.comb$crown/bites.data.comb$Area.crown.cm2
# Target crevice
bites.data.comb$Bites.cm2.targetcrevice <- bites.data.comb$targetcrevice/bites.data.comb$Area.target.cm2
# Other crevice
bites.data.comb$Bites.cm2.othercrevice <- bites.data.comb$othercrevice/bites.data.comb$Area.other.cm2
# Open crevice
bites.data.comb$Bites.cm2.opencrevice <- bites.data.comb$open/bites.data.comb$Area.open.cm2
# Outside
bites.data.comb$Bites.cm2.outside <- bites.data.comb$outside/bites.data.comb$Area.outside.cm2
head(bites.data.comb)

# Create column for bites per cm^2 per minutes for each category
# Crown
bites.data.comb$Bites.cm2.min.crown <- bites.data.comb$Bites.cm2.crown/bites.data.comb$video.length.min
# Target crevice
bites.data.comb$Bites.cm2.min.targetcrevice <- bites.data.comb$Bites.cm2.targetcrevice/bites.data.comb$video.length.min
# Other crevice
bites.data.comb$Bites.cm2.min.othercrevice <- bites.data.comb$Bites.cm2.othercrevice/bites.data.comb$video.length.min
# Open crevice
bites.data.comb$Bites.cm2.min.opencrevice <- bites.data.comb$Bites.cm2.opencrevice/bites.data.comb$video.length.min
# Outside
bites.data.comb$Bites.cm2.min.outside <- bites.data.comb$Bites.cm2.outside/bites.data.comb$video.length.min
str(bites.data.comb)
#write.csv(bites.data.comb, file="data/check-mergingofall.csv")


# Before graphing we need to turn it into long format 
bites.data.comb.long <- melt(bites.data.comb, id.vars=c("ColourCode", "Treatment", "tile.kind", "Genus.species"),
                         measure.vars=c("Bites.cm2.min.crown", "Bites.cm2.min.targetcrevice", 
                                        "Bites.cm2.min.othercrevice", "Bites.cm2.min.opencrevice",
                                        "Bites.cm2.min.outside"),
                         variable.name="Location.of.bites", value.name="Bites.cm2.min")
str(bites.data.comb.long)
head(bites.data.comb.long)


# What's the total # bites per cm2 and min?
bites.data.comb.long.summed <-
  bites.data.comb.long %>%
  dplyr::group_by(Genus.species) %>%
  dplyr::summarise(sumBites.cm2.min= sum(Bites.cm2.min))

bites.data.comb.long.summed
#write.csv(bites.data.comb.long.summed, file="data/summedBites-ontile-perspecies.csv")

# Graph it
# Bites per cm2 per minute split up by species, Treatment and crevice size
ggplot(bites.data.comb.long, aes(y=Bites.cm2.min, x=Location.of.bites, fill=Genus.species))+
  geom_col(position="dodge")+
  facet_grid(tile.kind~ Treatment)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1))+
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")
#ggsave("graphs/fishbitescrevicetile-LobvsnoLob-IncludingAllVideoSets.jpg")

# Bites per cm2 per minute split up by species and crevice size
ggplot(bites.data.comb.long, aes(y=Bites.cm2.min, x=Location.of.bites, fill=Genus.species))+
  geom_col(position="dodge")+
  facet_wrap(~ tile.kind)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")
#ggsave("graphs/fishbitescrevicetile-comparecrevicekinds-IncludingAllVideoSets.jpg")


# Let's look at it from a different angle: split by fish species
# and coloured by crevice type
ggplot(bites.data.comb.long, aes(y=Bites.cm2.min, x=Location.of.bites, fill=tile.kind))+
  geom_col(position="dodge")+
  facet_grid(~ Genus.species)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=15)) +
  xlab("Location of bites")+
  ylab("Bites taken (cm^-2/min-1)")


# Let's look at it from a different angle: split by fish species and Tile
# and coloured by crevice location
ggplot(bites.data.comb.long, aes(y=Bites.cm2.min, x=tile.kind, fill=Location.of.bites))+
  geom_col(position="fill")+
  facet_grid(~ Genus.species)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 45, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(angle=90, size=15))+
  xlab("Tile kind")+
  ylab("Proportion of bites on each tile")+
  theme(legend.text= element_text(size=15)) +
  theme(legend.title = element_text(size= 20))+
  scale_fill_discrete(name="Location of bites")
#ggsave("graphs/Proportion-bites-in-crevices-allvideosets.jpg")





#-----------------------------------------------------
# To make a different graph (showing medium, small and large crevices beside each other)
# I need to manipulate data further

str(bites.data.comb)
# grab important columns
newdata1 <- bites.data.comb[c("Genus.species", "tile.number", 
                              "tile.kind", "Bites.cm2.min.targetcrevice" )]
str(newdata1)
newdata1$tile.number <- as.factor(newdata1$tile.number)
head(newdata1)
# turn into wide format
newdata1.wide <- dcast(newdata1, tile.number + Genus.species ~ tile.kind, value.var="Bites.cm2.min.targetcrevice", fun.aggregate=sum)
head(newdata1.wide)

# grab important columns to merge with
newdata2 <- bites.data.comb[c("Genus.species", "tile.number", "Bites.cm2.min.crown")]
str(newdata2)
newdata2$tile.number <- as.factor(newdata2$tile.number)
head(newdata2)
# summ bites on crown by tile.number and Genus.species
newdata3 <- newdata2 %>%
  dplyr::group_by(tile.number, Genus.species) %>%
  dplyr::summarise(Bites.cm2.min.crown=sum(Bites.cm2.min.crown))
head(newdata3)            

# Now merge them
newdata4 <- full_join(newdata1.wide, newdata3, by=c("Genus.species", "tile.number"))
head(newdata4)
str(newdata4)

# rename Bites.cm2.min.crown to crown
colnames(newdata4)[colnames(newdata4) == "Bites.cm2.min.crown"] <- "crown"
str(newdata4)

# now turn back into long format for plotting
newdata4.long <- melt(newdata4, id.vars=c("tile.number", "Genus.species"),
                         measure.vars=c("small crevice", 
                                        "medium crevice", "large crevice",
                                        "crown"),
                         variable.name="Location.of.bites", value.name="Bites.cm2.min")
head(newdata4.long)

# need to sum the bites over all tiles 
newdata4.long.sum <- newdata4.long %>%
  dplyr::group_by(Location.of.bites, Genus.species) %>%
  dplyr::summarize(Bites.cm2.min.sum=sum(Bites.cm2.min))

head(newdata4.long.sum)

# Plot
ggplot(newdata4.long.sum, aes(y=Bites.cm2.min.sum, x= Location.of.bites, fill=Location.of.bites))+
  geom_col(position="dodge")+
  scale_fill_manual(values=c("#99CCFF","#6699FF","#0000CC","#000066"))+
  facet_wrap(~Genus.species, ncol=3)+
  theme_bw()+
  theme(panel.grid.major= element_blank(), panel.grid.minor= element_blank())+
  theme(axis.text.x = element_text(angle= 90, hjust=1)) +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 22)) +
  theme(strip.text = element_text(size=12, face="italic"))+
  xlab("Location of bites")+
  ylab(bquote("Bites (cm"^-2~min^-1*")"))+
  theme(legend.position="none")+
  coord_cartesian(ylim=c(0,2.5))
#ggsave("graphs/Bitesindifferentcrevices.jpg")
