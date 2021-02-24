### Cleaning percent cover and biomass for final analyses
library(tidyr) 
library(here)
library(plyr)
expdes <- read.csv(here::here("CedarCreekAnalyses","data", "BioCONExpDes.csv"))
tbiomass <- read.delim(here::here("CedarCreekAnalyses","data","AbvGrndBiomass_2017.txt"))

tbiomass$Date <- as.Date(tbiomass$Date,"%m/%d/%Y")
tbiomass$Month <- as.numeric(format(tbiomass$Date, format = "%m"))
tbiomass$Year <- as.numeric(format(tbiomass$Date, format = "%Y"))
tbiomass <- tbiomass[tbiomass$Month == 8,]

tbiomass$Species<- gsub(pattern = "Achillea millefolium ", replacement = "Achillea millefolium", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Bouteloua gracilis ", replacement = "Bouteloua gracilis", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Asclepias tuberosa ", replacement = "Asclepias tuberosa", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Schizachyrium scoparium ", replacement = "Schizachyrium scoparium", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Amorpha canescens ", replacement = "Amorpha canescens", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Bromus inermis ", replacement = "Bromus inermis", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Agropyron repens ", replacement = "Agropyron repens", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Lespedeza capitata ", replacement = "Lespedeza capitata", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Petalostemum villosum ", replacement = "Petalostemum villosum", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Poa pratensis ", replacement = "Poa pratensis", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "poa pratensis", replacement = "Poa pratensis", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Solidago rigida ", replacement = "Solidago rigida", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Koeleria cristata ", replacement = "Koeleria cristata", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Lupinus perennis ", replacement = "Lupinus perennis", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Andropogon gerardi ", replacement = "Andropogon gerardi", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Sorghastrum nutans ", replacement = "Sorghastrum nutans", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "Anemone cylindrica ", replacement = "Anemone cylindrica", tbiomass$Species)
tbiomass$Species<- gsub(pattern = "bromus inermis", replacement = "Bromus inermis", tbiomass$Species)
tbiomass$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "Amorpha  canescens", replacement = "Amorpha canescens", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "Solidago sigida", replacement = "Solidago rigida", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "AnemoneCylindrica", replacement = "Anemone cylindrica", tbiomass$monospecies)
tbiomass$monospecies<- gsub(pattern = "bromusInermis", replacement = "Bromus inermis", tbiomass$monospecies)
# Change Cenriched to Cenrich
tbiomass$CO2.Treatment <- gsub(pattern = "Cenriched", replacement = "Cenrich", tbiomass$CO2.Treatment)
tbiomass$CO2.Treatment <- gsub(pattern = "Cenrich ", replacement = "Cenrich", tbiomass$CO2.Treatment)


newdf <- expdes %>% pivot_longer(c(3:18),names_to = 'Species', values_to = 'planted')
newdf$Species <- gsub("[.]", " ", newdf$Species)
newdf <- merge(tbiomass[tbiomass$CO2.Treatment == "Camb" & tbiomass$Nitrogen.Treatment == "Namb",], newdf, by = c("Plot", "Ring", "Species"), all.x = TRUE)
newdf <- newdf[which(newdf$planted == 1),]
newdf <- newdf[-which(newdf$Aboveground.Biomass..g.m.2. == 0 & newdf$CountOfSpecies > 1),]
newdf$planted[newdf$planted == 1 & newdf$Aboveground.Biomass..g.m.2. ==0] <- 0

rich.sp <- ddply(newdf, .(Year, Plot, CountOfSpecies, Species), colwise(max, .(planted)))
rich<- ddply(rich.sp, .(Year, Plot, CountOfSpecies), colwise(sum, .(planted)))
summary(rich)
names(rich)[c(3,4)] <- c("TreatmentSR", "ObservedSR")

mass.sp <- ddply(newdf, .(Year, Plot, Species), colwise(mean, .(Aboveground.Biomass..g.m.2.)))
# Sum up mass in each strip
mass<- ddply(mass.sp, .(Year, Plot), colwise(sum, .(Aboveground.Biomass..g.m.2.)))
names(mass)[3] <- "live.mass"
summary(mass)
comb <- merge(rich, mass, by=c("Year", "Plot"))
summary(comb)

#save data 
write.csv(tbiomass, here::here("CedarCreekAnalyses","data","biocon_for_ld.csv"))
write.csv(comb, here::here("CedarCreekAnalyses","data", "biocon_plantedbiomass_output.csv"))

          