##########################################
### Process NutNet Cover Data ############
## Laura Dee  Sept 14, 2018  ###############
# some updates on April 23, 2019 to visualize different abundance metrics #
# updates on Aug 28 2019 to include other ways of processing the rare, dom, etc variables (based on frequency and relative cover)
# also other groups: rare & non-rare, and rare-non-native, rare-native, non-native/non-rare
#Updates on Dec 28 2019 to include other cut offs for rare groups
##########################################

#June 22 2020 - Fix the DIgroup2, prior versions didnt have the cut off 2 for the DIgroup that matched the others. 
# check w/: DI group 1 vs 2
# table(cover$DIgroup2)
# 
# > table(cover$DIgroup)
# 
# Dominant        Rare Subordinate 
# 432        4521       16993 

#Close graphics and clear local memory
graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(foreign)
library(rmarkdown)

# setwd("~/Documents/Research")  # kaitlin
setwd("~/Dropbox/IV in ecology/NutNet/")  # laura
cover <- fread('full-cover-09-April-2018.csv',na.strings= 'NA')

## need to make max_cover NOT a character
cover$max_cover <- as.numeric(cover$max_cover)

#########################################################################################
## Compare Taxon in Live (live==1) or non-live (live == 0) ################################
##########################################################################################
## determine how many Taxon are listed as live==0 ##
table(cover$live==1)
a <- table(cover$Taxon, cover$live==0)
head(a)

cover = cover[which(cover$live == 1),]
######################################################################################################
### Native vs Non-Native Variables ###################################################################
#######################################################################################################
View(cover[which(cover$local_provenance == "Naturalised"),]) 
# one species at one site is categorized as "Naturalised" so combined this with "INT"

# covert Naturalised to INT (This use of Naturalised as a category is from one site)
cover[local_provenance=="Naturalised", local_provenance:="INT"]
# convert blank an NA entries to NULL
cover[local_provenance=="", local_provenance:="NULL"]
cover[local_provenance=="NA", local_provenance:="NULL"]
#convert NULL to UNK to combine 
cover[local_provenance=="NULL", local_provenance:="UNK"]

# all are all of the unknowns are a single site or across sites? Check with this line:
View(cover[which(cover$local_provenance == "UNK"),]) # they are across sites.

# Compute native and non-native richness by plot, site, year
cover[, sr_INT := length(unique(Taxon[local_provenance=="INT"])), by = .(plot, site_code, year)]
cover[, sr_NAT := length(unique(Taxon[local_provenance=="NAT"])), by = .(plot, site_code, year)]
cover[, sr_UNK := length(unique(Taxon[local_provenance=="UNK"])), by = .(plot, site_code, year)]

## do first differences (for the marginal efcomb[order(year), changeGround_PAR := Ground_PAR -shift(Ground_PAR), by =.(plot, site_code)]
cover[order(year), change_sr_INT := sr_INT-shift(sr_INT), by =.(plot, site_code)]
cover[order(year), change_sr_NAT := sr_NAT-shift(sr_NAT), by =.(plot, site_code)]
cover[order(year), change_sr_UNK := sr_NULL-shift(sr_UNK), by =.(plot, site_code)]

##########################################################################################################
##### Compute TOTAL & TOTAL LIVE RELATIVE COVER PER PLOT MEASURES ########################################
##########################################################################################################
# make a total cover in a plot, site, year
cover[,totplotcover.yr := sum(max_cover, na.rm=T), by=.(plot, site_code, year)]

# same for live species only, which is what we end up using to compute the other variables
cover[,totplotcover.yr.live := sum(max_cover[live=="1"], na.rm=T), by=.(plot, site_code, year)]

#Make a relative cover for each species in each plot, site, year
# based on TOTAL cover (including dead cover)
cover[,relative_sp_cover.yr := max_cover/totplotcover.yr]

# same for live species only
cover[,relative_sp_cover.yr.live := (max_cover*(live=="1"))/totplotcover.yr.live]
# in some cases, no live species in plot and year, so getting NA since totplotcover.yr.live is zero. Set these to zero.
cover[is.na(relative_sp_cover.yr.live),relative_sp_cover.yr.live:=0]

####################################################################################################################
### Dominance Variables using Dominance Indicator (DI) ############################################################################################
####################################################################################################################
## We use the Dominance Indicator (DI) metric from Avolio et al (2019), which separates dominance indication from its impact.

### Dominance indicator calculation ###
  # DI = (average relative abundance + relative frequency)/2
    #Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
    #Relative frequency = number of sampling units a species occurred / total number of sampling units
# Note: DI ranges from 0-1; Relative abundance can be any measure of abundance (here is it based on relative cover). Does not incorporate a measure of impact.
# There is not a cutoff for "which range from 0-1 = dominant species versus subordinate or rare, in the Avolio et al paper # 
# so if we want to group species in these groups, we will need to make one (then also should test robustness to that decision)

## Compute Relative Frequency per spp AT THE SITE (defining dominance in space, not time/across years)
# "Relative frequency = number of sampling units a species occurred / total number of sampling units" 
# this should be # of plots within a site that the species occurred in / total # of plots within a site, for pre-treatment year 

#total # of plots within a site, for pre-treatment year  & to filter to do just on the live species:
# again we use the pre-treatment year because we calculate the metrics at the site level and want to avoid classifying species post treatment
cover[, tot.num.plots := length(unique(plot[year_trt == 0])), by =.(site_code)]

#number of plots within a site, in the pre-treatment year, a species occurred in:
# & to filter to do just on the live species:
cover[, tot.num.plots.with.spp := length(unique(plot[year_trt == 0 & live==1])), by =.(site_code, Taxon)]
  # test to see if it works
    abisko.test = cover[site_code == "abisko.se" , ]

##Compute Relative Frequency: Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots, by = .( site_code)]

### Compute relative abundance in terms of live cover ###
#Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
# we compute the above per plot and then take the average for each species at the site and year:
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, site_code, year)]

year0.sp <- cover[cover$year_trt == 0,c(3,10,34)] # gets the species present in year 0
names(year0.sp)[3] <- "ave_rel_abundance_year0" # renaming last column for merge
year0.sp <- unique(year0.sp) # getting unique values - so one value per species per site
year0.sp$present_year0<-1
require(dplyr)
cover <- dplyr::left_join(cover,year0.sp) # joining the two dataframes together. NA values now where species was not present year 0

# if the species isnt present at a site in year_trt == 0, give the species a relative abundance of 0 in that year:
cover$ave_rel_abundance_year0[is.na(cover$ave_rel_abundance_year0)] = 0
cover$rel_freq.space[is.na(cover$rel_freq.space)] = 0

#check that it worked
summary(cover$ave_rel_abundance_year0)
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

# see and plot the summary of these metrics 
summary(cover$rel_freq.space)
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq.space, xlab = "Frequency at the site in pre-treatment year", main = "Frequency of occurrence")
table(cover$rel_freq.space)
summary(cover$ave_rel_abundance_year0 )
hist(cover$ave_rel_abundance_year0, xlab = "Average relative abundance at a site", main ="")

#plot correlation between relative abundance and frequency metrics
plot(cover$rel_freq.space,cover$rel_abundance_year0)

#######################################################################################################################
## Compute the DI per species per sie defined as:  DI = (average relative abundance + relative frequency)/2 #########
## considering live cover and pre-treatment year  ################################################################
##############################################################################################################################
cover[, DI := (ave_rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, site_code)]

## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI], xlab = "Dominance indicator (DI)", main = "Dominance indicator metric")

###########################################################################################################
####### Make Categorical Variables of Dominance too  - based on the DI metric #############################
###########################################################################################################
# then given each species a ranking into different categories,
# to look at changes in the number of those types of species and the impact on productivity. 
cover[, DIgroup:=cut(DI, breaks=c(-0.001,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#Calculate Species Richness in each group 
cover[, sr_domspp := length(unique(Taxon[DIgroup == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp := length(unique(Taxon[DIgroup == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp := length(unique(Taxon[DIgroup == "Subordinate"])), by = .(plot, site_code, year)]

### create a category that groups together all the non-rare
cover[, sr_non_rare_spp := length(unique(Taxon[DIgroup %in% c("Subordinate", "Dominant")])), by = .(plot, site_code, year)]
#create a rare grouping (including native and invasive)
cover[, rareSR_spp.DI1 := length(unique(Taxon[DIgroup %in% c("Rare")])), by = .(plot, site_code, year)]

summary(cover$sr_domspp) # max is 2 species that are dominant 
summary(cover$sr_subordspp) 
summary(cover$sr_rarespp) 
summary(cover$sr_non_rare_spp)

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

#########################################################################################################################
### Create variables that are combined groups of: #######################################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
##############################################################################################################################
# create combinations of all - as factor in a column
cover[, status.NN.RareDom := paste(DIgroup,local_provenance, sep = "_")]
table(cover$status.NN.RareDom)

# create a non-rare variable
cover[, non_rare_spp := DIgroup %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

###########################################
### SR variables by combined groupings ####
###########################################
##### for the unknown species are will the data processing three ways. #################
# We will run each scenario as sensitivity analysis.  #################################

## 1. Excluding them (as above) #### 
#do SR for non-native, rare:
cover[, sr_non.nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover[, sr_nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

## 2. Including the unknown spp origin all as native: ####
cover[, sr_nat_unk_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

# 3.Including them all as non-native: 
cover[, sr_non.nat_unk_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## look a the data 
  hist(cover$sr_non.nat_rare)
  table(cover$sr_non.nat_rare)
  hist(cover$sr_nat_rare)
  table(cover$sr_nat_rare)
  hist(cover$sr_non.nat_unk_rare)
  table(cover$sr_non.nat_unk_rare)
  hist(cover$sr_nat_unk_rare)
  table(cover$sr_nat_unk_rare)
  
## do the same for the non-rare variables: 
# 1. Create SR non-rare native and non-native excluding unknown species origin species
cover[, sr_non.rare_non.nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

## 2. Include the unknown spp origin all as native: ####
cover[, sr_non.rare_nat_unk := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## 3. Include the unknown spp origin all as nonnative: ####
cover[, sr_non.rare_non.nat_unk := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

### Make extra for sensitivity analyses:
#do SR for native and non-native for dom
cover[, sr_nat_dom := length(unique(Taxon[DIgroup == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom := length(unique(Taxon[DIgroup == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub := length(unique(Taxon[DIgroup == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub := length(unique(Taxon[DIgroup == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do for native, non-rare:
#cover[, sr_nat_non.rare := length(unique(Taxon[DIgroup %in% c("Subordinate", "Dominant") & local_provenance == "NAT"])), by = .(plot, site_code, year)]
# do for non-native, non-rare:
# cover[, sr_non.nat_non.rare := length(unique(Taxon[DIgroup %in%  c("Subordinate", "Dominant") & local_provenance == "INT"])), by = .(plot, site_code, year)]


######################################################################################################################
## Dominance and Rarity Variables based on relative abundance and frequency groups separately (versus DI above) #######
######################################################################################################################

######################################################################################################################
## Compute for Frequency: sub-ordinate and rare - Cut Off 1 based on breaks=c(0.0,0.2,0.8,1.0)   #####################
######################################################################################################################
cover[, Freq_group:=cut(rel_freq.space, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group  #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, freq_sr_domspp := length(unique(Taxon[Freq_group == "Dominant"])), by = .(plot, site_code, year)]
cover[, freq_sr_rarespp := length(unique(Taxon[Freq_group == "Rare"])), by = .(plot, site_code, year)]
cover[, freq_sr_subordspp := length(unique(Taxon[Freq_group == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$freq_sr_domspp)  #max 23
summary(cover$freq_sr_rarespp) #max 13
summary(cover$freq_sr_subordspp) #max 24

# create a rare variable for frequency
cover[, sr_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "FALSE"])), by = .(plot, site_code, year)]
cover[, sr_rare_non.nat.Freq:= length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.Freq := length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a non-rare variable for frequency
cover[, non_rare_spp.Freq := Freq_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.Freq  := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

######################################################################################################### #############
## Compute for Average Relative Abundance: sub-ordinate and rare - Cut off 1
######################################################################################################### #############
#summary(cover$rel_abundance_year0)
hist(cover$rel_abundance_year0)

cover[, RelAbund_group:=cut(rel_abundance_year0, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, relabund_sr_domspp := length(unique(Taxon[RelAbund_group == "Dominant"])), by = .(plot, site_code, year)]
cover[, relabund_sr_rarespp := length(unique(Taxon[RelAbund_group == "Rare"])), by = .(plot, site_code, year)]
cover[, relabund_sr_subordspp := length(unique(Taxon[RelAbund_group == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$relabund_sr_domspp) # max 2 
summary(cover$relabund_sr_rarespp) # max 31 (mean ~9)
summary(cover$relabund_sr_subordspp) #max 3, min 0 ?? mean 2

# create a non-rare variable for relative abundance
cover[, non_rare_spp.RelA := RelAbund_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.RelA := length(unique(Taxon[non_rare_spp.RelA == "TRUE"])), by = .(plot, site_code, year)]

#rare species measures
cover[, sr_rare_spp.RelA := length(unique(Taxon[non_rare_spp.RelA == "FALSE"])), by = .(plot, site_code, year)]
# compute the rare native vs non group!
cover[, sr_rare_non.nat.RelA:= length(unique(Taxon[non_rare_spp.RelA == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.RelA := length(unique(Taxon[non_rare_spp.RelA == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.RelA:= length(unique(Taxon[non_rare_spp.RelA == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.RelA := length(unique(Taxon[non_rare_spp.RelA == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.RelA := length(unique(Taxon[RelAbund_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.RelA := length(unique(Taxon[RelAbund_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.RelA := length(unique(Taxon[RelAbund_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.RelA := length(unique(Taxon[RelAbund_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

###########################################################################################################################
#### Compute Groups based on DI Metric and Cut Off 2 [] Rare, Dom, Subord Groups]  #################################################################################################
############################################################################################################################

## Do a different cut off for dominance, subordinate and rare:
cover[, DIgroup2:=cut(DI, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]
#richness in each group with different cutoff
cover[, sr_domspp2 := length(unique(Taxon[DIgroup2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp2 := length(unique(Taxon[DIgroup2 == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp2 := length(unique(Taxon[DIgroup2 == "Subordinate"])), by = .(plot, site_code, year)]
summary(cover$sr_domspp2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SR variables by combined groupings ####### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#do SR for non-native, rare:
cover[, sr_non.nat_rare2 := length(unique(Taxon[DIgroup2 == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover[, sr_nat_rare2 := length(unique(Taxon[DIgroup2 == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a non-rare variable  for DIgroup 2
cover[, non_rare_spp2 :=  length(unique(Taxon[DIgroup2 %in% c("Subordinate", "Dominant")])), by = .(plot, site_code, year)]
#create a rare grouping (including native and invasive) for DIgroup 2
cover[, rare_spp.DI2 := length(unique(Taxon[DIgroup %in% c("Rare")])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom2 := length(unique(Taxon[DIgroup2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom2 := length(unique(Taxon[DIgroup2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub2 := length(unique(Taxon[DIgroup2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub2 := length(unique(Taxon[DIgroup2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

# # compute change in richness in each group for this diff cutoff (though only the dominance one should change)
# cover[order(year), change_sr_domspp2 := sr_domspp2 -shift(sr_domspp2), by =.(plot, site_code)]
# cover[order(year), change_sr_rarespp2 := sr_rarespp2 -shift(sr_rarespp2), by =.(plot, site_code)]
# cover[order(year), change_sr_subordspp2 := sr_subordspp2 -shift(sr_subordspp2), by =.(plot, site_code)]

######################################################################################################################
## Compute for Frequency: sub-ordinate and rare - Cut Off 2 based on breaks=c(0.0,0.4,0.8,1.0).
######################################################################################################################
cover[, Freq_group2:=cut(rel_freq.space, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group  #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, freq2_sr_domspp := length(unique(Taxon[Freq_group2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, freq2_sr_rarespp := length(unique(Taxon[Freq_group2 == "Rare"])), by = .(plot, site_code, year)]
cover[, freq2_sr_subordspp := length(unique(Taxon[Freq_group2 == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$freq2_sr_domspp)  #previous was max 23; now: 23
summary(cover$freq2_sr_rarespp) #max 13; 19
summary(cover$freq2_sr_subordspp) #max 24; 18

# create a non-rare variable for frequency
cover[, non_rare_spp.Freq2 := Freq_group2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.Freq2  := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a rare variable for frequency
cover[, sr_rare_spp.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "FALSE"])), by = .(plot, site_code, year)]
cover[, sr_rare_non.nat.Freq2:= length(unique(Taxon[non_rare_spp.Freq2 == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.Freq2 := length(unique(Taxon[Freq_group2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.Freq2 := length(unique(Taxon[Freq_group2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.Freq2 := length(unique(Taxon[Freq_group2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.Freq2 := length(unique(Taxon[Freq_group2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

######################################################################################################### #############
## Compute for Average Relative Abundance: sub-ordinate and rare -- CUT OFF 2: breaks=c(0.0,0.4,0.8,1.0)
######################################################################################################### #############
cover[, RelAbund_group2 :=cut(rel_abundance_year0, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group 
cover[, relabund_sr_domspp2 := length(unique(Taxon[RelAbund_group2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, relabund_sr_rarespp2 := length(unique(Taxon[RelAbund_group2 == "Rare"])), by = .(plot, site_code, year)]
cover[, relabund_sr_subordspp2 := length(unique(Taxon[RelAbund_group2 == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$relabund_sr_domspp2) # max 2 
summary(cover$relabund_sr_rarespp2) # max 43
summary(cover$relabund_sr_subordspp2) #max 3, min 0  mean 2 

# create a non-rare variable for relative abundance
cover[, non_rare_spp.RelA2 := RelAbund_group2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE"])), by = .(plot, site_code, year)]

#rare species measures
cover[, sr_rare_spp.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE"])), by = .(plot, site_code, year)]
# compute the rare native vs non group!
cover[, sr_rare_non.nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.RelA2 := length(unique(Taxon[RelAbund_group2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.RelA2 := length(unique(Taxon[RelAbund_group2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.RelA2 := length(unique(Taxon[RelAbund_group2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.RelA2 := length(unique(Taxon[RelAbund_group2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]


######################################################################################################### #############
## Compute for variables for DI = NA species 
######################################################################################################### #############

table(cover$DI == "NA")
table(cover$DIgroup == "NA", cover$local_provenance)

table(cover$DI, cover$local_provenance)

cover[, NA_DI_NAT := length(unique(Taxon[DI == "NA" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, NA_DI_INT := length(unique(Taxon[DI == "NA" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, NA_DI_UNK := length(unique(Taxon[DI == "NA" & local_provenance == "UNK"])), by = .(plot, site_code, year)]


summary(cover$NA_DI_INT)
table(cover$NA_DI_INT)
table(cover$NA_DI_UNK)

table(cover$NA_DI_NAT)

NAspp = table(cover$DI, cover$Taxon)
write.csv(NAspp, "NAspp.csv")
na = cover[DI == "NA",]
NAspp2 = table(na$DI, na$Taxon)
write.csv(NAspp2, "NAspp2.csv")

########################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above . 
#**** NOTE: If using the full dataset, likely need to check for others since there could be more in the full dataset!!!! 
### but here is OK because this will be merge with the comb data for control plots that I fixed 
cover = cover[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

# write as csv datafile to use for R
# write.csv(cover, "NutNetCoverData_ProcessedAug2019.csv")

#############################################################################
### References ############################################################
############################################################################

# Avolio, M.L., et al. (2019). Demystifying dominant species. New Phytol., 223, 1106–1126.