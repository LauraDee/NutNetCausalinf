###############################################################
## Create SM Data figures for Dee et al ###################
###############################################################
#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)

setwd("~/Dropbox/IV in ecology/NutNet")
##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-28-Apr-2017.csv'
comb <- fread("NutNetControlPlotDataToUseApril2018.csv",na.strings='NA')
merged.data <- fread("merged_NutNet_Oct2019.csv", na.strings = 'NA')

####################################################################################################
## Filter cover data to only control plots used  ###############################
##################################################################################################
comb = comb[is.control==T & has.5.yrs.data==T,]

##########################################################################################
### Process comb data to prep to use it in the models & merge #############################
############################################################################################
## Filter data to records with non-NA live_mass and non-NA richness.
## Models can generally handle NAs by implicitly dropping, but clustered SE
# # function currently does not.
 comb = comb[!is.na(live_mass) & !is.na(rich)]

# # make year a character, to be a dummy variable: 
comb$year <- as.character(comb$year)
cover$year <- as.character(cover$year)
# same with plot
comb$plot <- as.character(comb$plot)

###############################################################
## Figure S1  Plot-level (a &b) and site-level (c & d) biodiversity (as species richness) and productivity (as above-ground live mass) from between 2007-2017. 
###############################################################
par(mfrow=c(2,2))
plot(live_mass ~ rich, data = comb, main = "a)", xlab = "species richness", ylab = "live mass")
plot(log(live_mass) ~ log(rich), data = comb, main = "b)", xlab = "log(species richness)", ylab = "log(live mass)")
## also plot site-level data 
plot(comb$site_richness, comb$ave_site_live_mass, main = "c)", xlab = "species richness", ylab = "average site live mass")
plot(log(comb$site_richness), log(comb$ave_site_live_mass), main = "d)", ylab = "log(average site live mass)", xlab = "log(site richness)") 
mtext("Plot-level", side=3, outer=TRUE, line=-3)
mtext("Site-level", side=3, outer=TRUE, line=-17)
###############################################################
# Figure S2 ###################################################
###############################################################
FigS2 <- ggplot(data = comb, aes(x = changerich)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigS2
############################################################### 
## Figure S3. Change in species evenness per plot (overall data)
###############################################################
FigS3 <- ggplot(data = comb, aes(x = changeEvenness)) + geom_histogram() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in evenness") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigS3








##########################################################################
## Plot the # of invasive species and ovreall SR: #####################
####################################################################################
## SM Section 7 plot distribution of metrics 
## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI], xlab = "Dominance indicator (DI)", main = "Dominance indicator metric")
# summary(cover[live == 1,DI])

plot( mech.data$rich, mech.data$sr_INT, xlab = "overall species richness", ylab = "INT species richness", main = "Overall vs INT species richness")

##########################################################################c
## Plot rare SR  vs overall SR ##########################################

#**do I need to recreate these variables?

plot( mech.data$rich, mech.data$sr_nat_rare, xlab = "Overall species richness", ylab = "native rare species richness", main = "Overall vs native rare species richness")
plot( mech.data$rich, mech.data$sr_rarespp, xlab = "overall species richness", ylab = "all rare species richness", main = "Overall vs all rare species richness")




#hist of overall change in native, rare species  
hist(merged.data$change_rare.native, main = "Plot-level change in rare species richness (native)" ,
     xlab = "Change in number of species")
#print site level #s
table(merged.data$site, merged.data$change_rare.native)

#hist of overall change in native, non rare species  
hist(merged.data$change_nonrare.native, main = "Plot-level change in non-rare species richness (native)" , 
     xlab = "Change in number of species")

# # function currently does not.
merged.data = merged.data[!is.na(change_rare.native)]

## plot change in rare species by site 
FigSX <- ggplot(data = merged.data, aes(x = change_rare.native)) + geom_histogram()+ facet_wrap(~site) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

#** cut 
plot( mech.data$rich, mech.data$sr_rarespp, xlab = "overall species richness", ylab = "all rare species richness", main = "Overall vs all rare species richness")


############################################################################################################################
#### Figure Extra. Plot Correlations between all of the SR grouping variables ###################################################################
#############################################################################################################################
richnessvars.fig4c <- c("sr_non.rare_nat", "sr_non.rare_non.nat", "sr_non.nat_rare", "sr_nat_rare")
richnessvars.to.plot <- mech.data[,richnessvars.fig4c, with=F] #with = F will then use the data.frame conventions, using "" for col names
pairs(richnessvars.to.plot)

sr.metrics <- mech.data[, .(sr_non.nat_rare, sr_nat_rare, non_rare_spp,
                            sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom,
                            sr_non.nat_dom, sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare )]
cor(sr.metrics)
corrplot(cor(sr.metrics), method = "square", tl.cex = .5)


# compute change in richness in each group 
cover[order(year), change_sr_domspp := sr_domspp -shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), change_sr_rarespp := sr_rarespp -shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), change_sr_subordspp := sr_subordspp -shift(sr_subordspp), by =.(plot, site_code)]
cover[order(year), change_sr_non_rare_spp := sr_non_rare_spp -shift(sr_non_rare_spp), by =.(plot, site_code)]



summary(cover$change_sr_domspp) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.0000  0.0000  0.0000  0.0012  0.0000  1.0000    3005 
summary(cover$change_sr_subordspp) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0327   0.0000  10.0000     3005 
summary(cover$change_sr_rare) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0199   0.0000   7.0000     3005 


