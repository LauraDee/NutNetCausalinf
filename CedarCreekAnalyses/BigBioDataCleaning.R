# EWS: Read in BigBio (LTER E120) biomass data and calculate 
# mass and observed richness in each sample

#If you run the e120-biomass-data-12-04.R script – 2 csv files will be made.
# e120-biomass-data-output.csv is the one Eric made – we don’t want to use that because it has all the species.
# e120-plantedbiomass-data-output.csv is the one you want to use because it only includes planted species.

# Clear all existing data
rm(list=ls())
# Close graphics windows
graphics.off()

library(plyr)
library(here)

#######    READ IN DATA    #######
d<-read.csv(here::here("CedarCreekAnalyses", "data", "e120-biomass-data-2020-12-04.csv"), strip.white=T)
d$mass <- d$Biomass..g.m2.  
d$Species <- toupper(d$Species)
names(d)[] <- tolower(names(d))
summary(d)

bbspecies <- read.csv(here::here("CedarCreekAnalyses", "data", "CDRSPDat.csv"))
# Remove oak trees and seedlings
d <- d[grep("QUERCUS", d$species, invert=TRUE),]

# Flag unsorted, non-vascular, and dead masses
d$live <- TRUE
d$live[grep("LITTER", d$species)] <- FALSE
d$live[grep("NEST", d$species)] <- FALSE

d$sorted <- TRUE
d$sorted[grep("UNSORTED", d$species)] <- FALSE
d$sorted[grep("MISCELLANEOUS", d$species)] <- FALSE
d$sorted[grep("WEEDS", d$species)] <- FALSE
d$sorted[grep("GRASSES", d$species)] <- FALSE
d$sorted[grep("MATTER", d$species)] <- FALSE

d$vasc <- TRUE
d$vasc[grep("MOSS", d$species)] <- FALSE
d$vasc[grep("FUNGI", d$species)] <- FALSE

with(d[d$live,],sort(unique(species)))
with(d[d$sorted,],sort(unique(species)))
with(d[d$vasc,],sort(unique(species)))

# Sort out strips & substrips
with(d, table(year, strip))
with(d, table(year, is.na(strip)))
with(d, table(year, substrip))
with(d[d$year==2007 | d$year==2009,], table(substrip, strip, year))

# When was data sorted
with(d, table(year, sorted))
with(d, table(year, live))
with(d, table(year, vasc))
with(d, table(year,month)) 

# Samples not sorted before 2001 and in 2009
d <- d[d$year > 2000,]
d <- d[d$year != 2009,]
#in 2001, sampled in june and august - want august only data
d <- d[d$month != 6,]

# We need strip designation
d <- d[!is.na(d$strip),]

# Sum live mass and richness by strip

# For richness look at live, sorted, vascular plants
# Make sure there is only one measurement per species per strip
d$rich <- 1
# rich.strip.sp <- ddply(d[d$sorted & d$vasc & d$live,], .(year, plot, strip, numsp, spnum, species), colwise(max, .(rich)))
# # Summary species in each strip
# rich.strip <- ddply(rich.strip.sp, .(year, plot, strip, numsp, spnum), colwise(sum, .(rich)))
# summary(rich.strip)

# For mass look at live, vascular plants
# Make sure there is only one measurement per species per strip
d$mass.live <- d$mass
# mass.strip.sp <- ddply(d[d$vasc & d$live,], .(year, plot, strip, species), colwise(mean, .(mass.live)))
# # Sum up mass in each strip
# mass.strip <- ddply(mass.strip.sp, .(year, plot, strip), colwise(sum, .(mass.live)))
# summary(mass.strip)
# 
# # Merge mass and richness data
# comb <- merge(rich.strip, mass.strip, by=c("year", "plot", "strip"))
# summary(comb)
# 
# # Output data
# write.csv(comb, file="e120-biomass-data-output.csv")


### Planted species only ###
library(tidyr)

planted.sp <- as.data.frame(names(d[,c(17:34)])) # pull planted species from 
names(planted.sp) <- "X5Lspecid"
bbspecies[,8] <- tolower(bbspecies[,8])
bbspecies <- bbspecies[,c(1,8)]
planted.sp <- merge(planted.sp, bbspecies)
names(planted.sp)[2] <- "planted.sp"


newdf <- d %>% pivot_longer(c(17:34),names_to = 'X5Lspecid', values_to = 'planted')
newdf <- merge(newdf, planted.sp, all.x = TRUE) 
### amopet and monsol not matching over with any species
### assuming amopet = Amorpha canescens and monsol = Monarda fistulosa 
newdf$planted.sp[newdf$X5Lspecid == "amopet"] <- "Amorpha canescens"
newdf$planted.sp[newdf$X5Lspecid == "monsol"] <- "Monarda fistulosa"
newdf$planted.sp <- toupper(newdf$planted.sp)


newdf <- newdf[newdf$planted.sp == newdf$species & newdf$planted == 1,]

rich.strip.sp1 <- ddply(newdf, .(year, plot, strip, numsp, spnum, species), colwise(max, .(rich)))
# Summary species in each strip
rich.strip1 <- ddply(rich.strip.sp1, .(year, plot, strip, numsp, spnum), colwise(sum, .(rich)))
summary(rich.strip1)

mass.strip.sp1 <- ddply(newdf, .(year, plot, strip, species), colwise(mean, .(mass.live)))
# Sum up mass in each strip
mass.strip1 <- ddply(mass.strip.sp1, .(year, plot, strip), colwise(sum, .(mass.live)))
summary(mass.strip1)

# Merge mass and richness data
comb1 <- merge(rich.strip1, mass.strip1, by=c("year", "plot", "strip"))
summary(comb1)

comb1$spnum <- NULL

# Output data
write.csv(comb1, file=here::here("CedarCreekAnalyses", "data", "e120-plantedbiomass-data-output.csv"))
