##########################################
### Process NutNet Cover Data ############
## Laura Dee  # May 2021      ############
##########################################
#Code written by Laura Dee, Revised and Checked by Kaitlin Kimmel, Checked by Steve Miller and Carlos A.

#Close graphics and clear local memory
graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(foreign)
library(rmarkdown)

#setwd("~/Documents/Research")  # kaitlin
#setwd("~/Dropbox/IV in ecology/NutNet/")  # laura
setwd("~/Documents/GitHub/NutNetCausalinf/data/")  

#read in file as a data table
cover <- fread('full-cover-09-April-2018.csv',na.strings= 'NA')

## Need to make max_cover NOT a character
#Visual cover estimates are made to the nearest 1% for every species contained within (or over-hanging) the subplot. 
cover$max_cover <- as.numeric(cover$max_cover)

#Note that Several sites sample multiple subplots within plots: elkh.us , ucsc.us, anti.ec -- and not every plot at a site.
#Elkhorn, Antisana, and  University of California - Santa Cruz. However, 
# these sites do not meet our criteria for inclusion in the analysis, which is at least 5 years of data and a pre-treatment year 
# as of the April-2018 data version. 
# IF these end up being included in other analyses using this code, there needs to be an agreed-upon way of aggregating max cover to the plot level from the subplot
# level that is consistent. Examples might include 1) using a representative subplot that is in all years of data, or 
# 2) averaging max cover across subplots within a plot since subplots are of equal size.

########################################################################################################
## Compare Taxon in Live (live==1) or non-live (live == 0) #############################################
#######################################################################################################
#** Cut in btwn: # determine how many Taxon are listed as live==0 ##
table(cover$live==1)

#to print or view what cover is dead 
a <- table(cover$Taxon, cover$live==0) 
head(a)
#*****

#filter data table to live cover
cover = cover[live==1,]

############################################################################################################
### Native vs Non-Native Variables #########################################################################
#############################################################################################################
# Running this next line, we see that one species at *one* site is categorized as "Naturalised," so we combined this category with "INT"
 # View(cover[which(cover$local_provenance == "Naturalised"),]) 

# covert Naturalised to INT (This use of Naturalised as a category is only from one site, marc.ar)
cover[local_provenance =="Naturalised", local_provenance:="INT"]

# convert blank an NA entries to NULL
cover[local_provenance=="", local_provenance:="NULL"] 
cover[local_provenance=="NA", local_provenance:="NULL"] #these NAs are stored as a string in the data.table so can run this.
#convert NULL to UNK to combine  #these NULLs are stored as a string in the data.table so can run this. Checked with:
#cover[is.null(local_provenance),]
cover[local_provenance=="NULL", local_provenance:="UNK"] 

# all are all of the unknowns are a single site or across sites? Check with this line:
 # View(cover[which(cover$local_provenance == "UNK"),]) # the species of unknown origins are across different sites.

############################################################################################################
### Prep data and species list  #########################################################################
#############################################################################################################
## First, only the species PRESENT in a plot are recorded in the cover data, so species that are present at a site but not a plot 
# are *not* listed (i.e.  with max_cover = 0). We need to fix that before computing average relative abundance at a site.

# First, we compute list of all species at a site over the time period in the data:
sp.at.site = unique(cover[,.(site_code, Taxon)])
# Want to create and merge one record per (Taxon, site_code, plot, year), using merge to flag records that weren't present 
# in original data.
site.plot.year.combos = unique(cover[,.(site_code, plot, year)])
expanded.spp.recs = merge(site.plot.year.combos, sp.at.site, by=c("site_code"), allow.cartesian = T)
cover = merge(cover, expanded.spp.recs, by=c("site_code", "plot", "year", "Taxon"), all.y=T)

# need to update yr_trt, trt, live, local_provenance for the records that got added
cover[,live:=1]
cover[,year_trt:=min(year_trt[!is.na(year_trt)], na.rm=T), by=.(plot, site_code, year)]
cover[,trt:=min(trt[!is.na(trt)], na.rm=T), by=.(plot, site_code, year)]
cover[,local_provenance:=min(local_provenance[!is.na(local_provenance)]),by=.(site_code, Taxon)]

# now compute #s we care about, including max_cover,sr_INT, sr_NAT, sr_UNK, totplotcover.yr.live, relative_sp_cover.yr.live, 
# tot.num.plots, tot.num.plots.with.spp, rel_freq.space
cover[is.na(max_cover), max_cover:=0] # if NA, species wasn't there in that plot and year, so cover should be zero

# Compute native, non-native, and unknown origin species richness by plot, site, year. Note filter to max_cover > 0 to 
# consider only species that were actually present.
cover[, sr_INT := length(unique(Taxon[local_provenance=="INT" & max_cover>0])), by = .(plot, site_code, year)]
cover[, sr_NAT := length(unique(Taxon[local_provenance=="NAT" & max_cover>0])), by = .(plot, site_code, year)]
cover[, sr_UNK := length(unique(Taxon[local_provenance=="UNK" & max_cover>0])), by = .(plot, site_code, year)]

#####
##### Compute a convenience column that says whether a species was present in a plot in a site in the pre-treatment year
cover[,present_year0:=max_cover[year_trt==0]>0, by=.(Taxon,site_code,plot)]

#compute the SR of species that were not present in year 0 for each plot and year but are present in a given year in that plot
cover[, sr_NA := length(unique(Taxon[present_year0 == FALSE & max_cover>0])), by = .(plot, site_code, year)]


#*** ultimately cut btwn these lines to clean code*****
printNA <- table(cover$sr_NA, cover$site_code)
write.csv(printNA, "printNAbysite.csv")

cover.NA.unique = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt,  sr_NA)])
#plot 
Nasp <- ggplot(data =cover.NA.unique, aes(x = sr_NA)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Count of Number of Species that would be NA per plot and year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
Nasp
#*** ultimlately cut btwn these lines to clean code*****
#*

##########################################################################################################
##### Compute TOTAL & TOTAL LIVE RELATIVE COVER PER PLOT MEASURES ########################################
##########################################################################################################
# make a total cover in a plot, site, year. This includes live cover only.
cover[,totplotcover.yr.live := sum(max_cover, na.rm= T), by=.(plot, site_code, year)]

#Make anrelative cover for each species in each plot and year
# based on TOTAL cover (including only live cover as we already filtered to the data table to live above).
cover[,relative_sp_cover.yr.live := max_cover/totplotcover.yr.live]
#***to checK, run: sum(is.na(cover$relative_sp_cover.yr.live))

# make a site-level relative abundance for each species and year. This requires
#first summing the total cover per site and year and the total cover of each species in all plots at site 
cover[, totsitecover.yr := sum(totplotcover.yr.live, na.rm= T), by=.(site_code, year)]
cover[, tot_maxcover_site.yr  := sum(max_cover, na.rm= T), by=.(Taxon, site_code, year)]

#**need to set to 0 if cover$tot_maxcover_site.yr or totsitecover.yris NA....
cover[is.na(tot_maxcover_site.yr),tot_maxcover_site.yr := 0]
cover[is.na(totsitecover.yr),totsitecover.yr  := 0] #not necessary

#then divide these two numbers: 
cover[, relative_abundance_spp_site.yr  :=  tot_maxcover_site.yr/totsitecover.yr, by=.(Taxon, site_code, year)]

#create a variable for just year 0 
cover[, relative_abundance_spp_site.yr0 := min(relative_abundance_spp_site.yr[year_trt==0]), by=.(Taxon, site_code)]
cover[is.infinite(relative_abundance_spp_site.yr0),relative_abundance_spp_site.yr0 := NA]

# if the species isn't present at a site in year_trt == 0, give the species a relative abundance of 0 in that year:


### Next step --create a relative frequency in year 0 variable #####
#total # of plots within a site, for pre-treatment year:
# again we use the pre-treatment year because we calculate the metrics at the site level and want to avoid classifying species post treatment
cover[, tot.num.plots := length(unique(plot[year_trt == 0])), by =.(site_code)]  #this will work because no records of max_cover = 0.

#number of plots within a site, in the pre-treatment year, that a species occurred in:
cover[, tot.num.plots.with.spp := length(unique(plot[year_trt == 0 & max_cover>0])), by =.(site_code, Taxon)]
#*** test to see if it works
azi.cn.test = cover[site_code == "azi.cn" , ]
#****azi.cn.test[Taxon=="AGROSTIS STOLONIFERA",]


##Compute Relative Frequency in year 0.
## Relative frequency = number of plots at a site in year 0 a species occurred / total number of plots at a site in year 0" 
# If a site has no records for plots in a pre-treatment year (year_trt==0), rel_freq.space will be NA.
# That's fine -- these sites will be filtered out later
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots]
cover[is.na(rel_freq.space),rel_freq.space  := 0] #32520 NAs


#################################################################################################################################
## Run this Code (and all following code)  if you want to Filter data to only species present in year 0 and save that dataset  #####
#################################################################################################################################
cover_present_year0 = cover[present_year0 == TRUE,]
# write.csv(cover_present_year0, "cover_present_year0May142021.csv")
# # cover_present_year0[, DI := (ave_rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, site_code)]
# cover = cover_present_year0 
# 
## **test that this worked 
# cover[, sr_NA.test := length(unique(Taxon[present_year0 == FALSE & max_cover>0])), by = .(plot, site_code, year)]
# summary(cover$sr_NA.test)

###################################################################################################################################################
####### Make Categorical Variables to Label Spp as Dominant, Subordinant, and Rare   - based on the relative abundance Quantiles per Site #########################################
#############################################################################################################################################

#**Note to self ****
#0 quartile = 0 quantile = 0 percentile
# 1 quartile = 0.25 quantile = 25 percentile
# 2 quartile = .5 quantile = 50 percentile (median)
# 3 quartile = .75 quantile = 75 percentile
# 4 quartile = 1 quantile = 100 percentile

unique.ras = unique(cover_present_year0[, .(site_code, Taxon, relative_abundance_spp_site.yr0)])

unique.ras[,RAquant0.6:=quantile(relative_abundance_spp_site.yr0, probs=0.6), by=site_code]
unique.ras[,RAquant0.95:=quantile(relative_abundance_spp_site.yr0, probs=0.95), by=site_code]
unique.ras[,RAsite_group:=ifelse(relative_abundance_spp_site.yr0<RAquant0.6,"Rare",
                                     ifelse(relative_abundance_spp_site.yr0<RAquant0.95, "Subordinate","Dominant"))]
#re-merge the quantiles and classifications into the cover_present_year0 dataset
unique.ras[,relative_abundance_spp_site.yr0:=NULL] # drop before re-merge
cover_present_year0 = merge(cover_present_year0, unique.ras, by=c("site_code", "Taxon"))

#sanity check
cdcr = cover_present_year0[site_code =="cdcr.us",]
table(cdcr$Taxon, cdcr$RAsite_group.y)
konz  =  cover_present_year0[site_code =="konz.us",]
table(konz$Taxon, konz$RAsite_group.y)

#whats the breakdown of species classified in each group overall 
table(unique.ras$RAsite_group)
#whats the breakdown of species classified in each group by site
table(unique.ras$site_code, unique.ras$RAsite_group)

##Where did the saline.us site go?

########################################################################################################################
### Create variables that are combined groups of: #######################################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
##############################################################################################################################
# create combinations of all - as factor in a column
cover_present_year0[, status.NN.RareDom := paste(RAsite_group,local_provenance, sep = "_")]
table(cover_present_year0$status.NN.RareDom)

# 3.Including them all as non-native: 
# create a non-rare variable
cover_present_year0[, non_rare_spp := RAsite_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]


# From here on out, we don't need/want the extra records that were created in order
# to properly calculate rarity, relative abundance, etc if a species showed up
# in a plot after the first year or existed in some plots within a site.
#
# Filter back down to (original) records with max_cover > 0 
cover_present_year0 = cover_present_year0[max_cover>0,]
cover = cover[max_cover>0,]

## to make a plot of number of introduced species over time --left off here
cover[, sr_INT:= length(unique(Taxon[local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_INT.site:= length(unique(Taxon[local_provenance == "INT"])), by = .(site_code, year)]

#*** CUT IN BTWN *****
cover.int.unique = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt,  sr_INT.site, sr_INT)])

INTsp <- ggplot(data = cover.int.unique , aes(x = year, y = sr_INT.site) + geom_point() # + facet_wrap(~site_code) + theme_bw() +
  + labs(x = "Count of Number of Species that would be NA per plot and year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
INTsp

invasionmod <- lm(sr_INT.site ~ year_trt, data = cover.int.unique)
summary(invasionmod)
#*** CUT IN BTWN *****


#################################################################################################################################
### SR variables by combined groupings based on site relative abundance ###############################################
#################################################################################################################################

##### for the unknown species are will the data processing three ways. #################
# We will run each scenario as sensitivity analysis.  #################################

## 1. Excluding them (as above) #### 
#do SR for non-native, rare:
cover_present_year0[, sr_non.nat_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover_present_year0[, sr_nat_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

## 2. Including the unknown spp origin all as native: ####
cover_present_year0[, sr_nat_unk_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

# 3.Including them all as non-native: 
cover_present_year0[, sr_non.nat_unk_rare := length(unique(Taxon[RAsite_group== "Rare" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## do the same for the non-rare variables: 
# 1. Create SR non-rare native and non-native excluding unknown species origin species
cover_present_year0[, sr_non.rare_non.nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.rare_nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
#**check cover_present_year0[site_code=="yarra.au" & plot==8,.(plot, subplot, Taxon, year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]

## 2. Include the unknown spp origin all as native: ####
cover_present_year0[, sr_non.rare_nat_unk := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## 3. Include the unknown spp origin all as nonnative: ####
cover_present_year0[, sr_non.rare_non.nat_unk := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

### Make extra variables for sensitivity analyses:
#do SR for native and non-native for dom
cover_present_year0[, sr_nat_dom := length(unique(Taxon[RAsite_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.nat_dom := length(unique(Taxon[RAsite_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover_present_year0[, sr_nat_sub := length(unique(Taxon[RAsite_group== "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.nat_sub := length(unique(Taxon[RAsite_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

### Create a variable for the unknown SR
cover_present_year0[, sr_unk_rare := length(unique(Taxon[RAsite_group== "Rare" &  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## look a the data 
#***ultimately cut between these lines
cover_present_year0[site_code=="yarra.au" & plot==8,.(plot, year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare)]
hist(cover_present_year0$sr_non.nat_rare)
table(cover_present_year0$sr_non.nat_rare, cover_present_year0$year)
hist(cover_present_year0$sr_nat_rare)
table(cover_present_year0$sr_nat_rare)
hist(cover_present_year0$sr_non.nat_unk_rare)
table(cover_present_year0$sr_non.nat_unk_rare)
hist(cover_present_year0$sr_nat_unk_rare)
table(cover_present_year0$sr_nat_unk_rare)
#***
#*



#################################################################################################################################
### SR variables by combined groupings based on site relative Frequency ###############################################
#################################################################################################################################

# first need to create the cut offs for groups based on the site-level distributions of relative frequency:
unique.freq = unique(cover_present_year0[, .(site_code, Taxon, rel_freq.space)])

unique.freq[,Fquant0.6:=quantile(rel_freq.space, probs=0.6), by=site_code]
unique.freq[,Fquant0.95:=quantile(rel_freq.space, probs=0.95), by=site_code]
unique.freq[,Freq_group :=ifelse(rel_freq.space <Fquant0.6,"Rare",
                                 ifelse( rel_freq.space <Fquant0.95, "Subordinate","Dominant"))]

#re-merge the quantiles and classifications into the cover_present_year0 dataset
unique.freq[, rel_freq.space :=NULL] # drop before re-merge
cover_present_year0 = merge(cover_present_year0, unique.freq, by=c("site_code", "Taxon"))

#sanity check -- ** KAILIN LOOK HERE ***
cdcr = cover_present_year0[site_code =="cdcr.us",]
table(cdcr$Taxon, cdcr$Freq_group)
konz  =  cover_present_year0[site_code =="konz.us",]
table(konz$Taxon, konz$Freq_group)

#whats the breakdown of species classified in each group overall 
table(unique.freq$Freq_group)
#whats the breakdown of species classified in each group by site
table(unique.freq$site_code, unique.freq$Freq_group)  ##Where did the saline.us site go?


##################################################################
### Compute SR variables based on Frequency Groups above :#########
#####################################################################
# create a non-rare variable for frequency
cover_present_year0[, non_rare_spp.Freq := Freq_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover_present_year0[, sr_non_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE"])), by = .(plot, site_code, year)]

# create a rare variable for frequency
cover_present_year0[, sr_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "FALSE"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_rare_non.nat.Freq:= length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_rare_nat.Freq := length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover_present_year0[, sr_non.rare_non.nat.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.rare_nat.Freq  := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover_present_year0[, sr_nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover_present_year0[, sr_nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]


##### Now create variables that deal with species of unknown origin. ##### 
cover_present_year0[, status.NN.FreqGroup := paste(Freq_group,local_provenance, sep = "_")]
table(cover_present_year0$status.NN.FreqGroup)

## 2. ## group all of the unknown species origins as native for both rare and non-rare groups:
#rare
cover_present_year0[, sr_rare_unk_nat.Freq := length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]
#non rare
cover_present_year0[, sr_non.rare_nat_unk.Freq  := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## 3. ## group all of the unknown species origins as non-native for both rare and non-rare groups: 
# rare non native
cover_present_year0[, sr_rare_non.nat_unk.Freq:= length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]
#non-rare non native
cover_present_year0[, sr_non.rare_non.nat_unk.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]


#####################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above . 
#**** NOTE FOR DATA-REUSE: If using the full dataset, likely need to check for others since there could be more in the full dataset!!!! 
### but here is OK because this will be merge with the comb data for control plots that I fixed****** 
cover_present_year0 = cover_present_year0[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

#create a version of the data running the above code with cover = cover_present_year0  
#to process all variables of SR counts and groupings on only species present in year 0 then printing as:
write.csv(cover_present_year0, "NutNetCoverData_ProcessedAug2019-PresentYear0only_RelAbundanceOnly.csv")











#################################################################################################################################################
### Dominance Variables using Dominance Indicator (DI) ############################################################################################
#################################################################################################################################################
## We use a Dominance Indicator, the DC_i metric, from Avolio et al (2019), which separates dominance indication from its impact.

### Dominance indicator calculation ###
# DI = (average relative abundance + relative frequency)/2
#Relative abundance = abundance of a species a in sampling unit (here a plot) / total abundance of all species in a sampling unit (plot)
# so we do an average plot level Relative abundance at each site per species.

#Relative frequency = number of sampling units a species occurred / total number of sampling units
# Note: DI ranges from 0-1; Relative abundance can be any measure of abundance (here is it based on relative cover). Does not incorporate a measure of impact.
# There is not a cutoff for "which range from 0-1 = dominant species versus subordinate or rare, in the Avolio et al paper # 
# so if we want to group species in these groups, we will need to make one (then also should test robustness to that decision)

#####
#### Compute Relative Frequency per spp AT THE SITE (defining dominance in space, not time/across years) ###
# "Relative frequency = number of sampling units a species occurred / total number of sampling units"   
# this should be # of plots within a site that the species occurred in / total # of plots within a site, for pre-treatment year 

#total # of plots within a site, for pre-treatment year:
# again we use the pre-treatment year because we calculate the metrics at the site level and want to avoid classifying species post treatment
cover[, tot.num.plots := length(unique(plot[year_trt == 0])), by =.(site_code)]  #this will work because no records of max_cover = 0.

#number of plots within a site, in the pre-treatment year, that a species occurred in:
cover[, tot.num.plots.with.spp := length(unique(plot[year_trt == 0 & max_cover>0])), by =.(site_code, Taxon)]
  #*** test to see if it works
  azi.cn.test = cover[site_code == "azi.cn" , ]
  #****azi.cn.test[Taxon=="AGROSTIS STOLONIFERA",]

##Compute Relative Frequency in year 0.
## Relative frequency = number of plots at a site in year 0 a species occurred / total number of plots at a site in year 0" 
# If a site has no records for plots in a pre-treatment year (year_trt==0), rel_freq.space will be NA.
# That's fine -- these sites will be filtered out later
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots]

#################################################################################################################################
### Compute average plot-level relative abundance in terms of live cover per the Avolio et al paper: #################################

#Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
# we compute the above per plot and then take the average for each species at the site and year:
cover[, ave_rel_abundance_over_time.live := mean(relative_sp_cover.yr.live), by = .(Taxon, site_code, year)]
# compute average site-level relative cover in year 0 per species 
cover[, ave_rel_abundance_year0 := ave_rel_abundance_over_time.live[year_trt == 0], by = .(Taxon, site_code, year)]

# if the species isn't present at a site in year_trt == 0, give the species a relative abundance of 0 in that year:
cover[is.na(ave_rel_abundance_year0),ave_rel_abundance_year0 := 0]

#*******check that it worked ******
summary(cover$ave_rel_abundance_year0)
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 
hist(cover$ave_rel_abundance_over_time.live) 

# see and plot the summary of these metrics 
summary(cover$rel_freq.space)

#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq.space, xlab = "Frequency at the site in pre-treatment year", main = "Frequency of occurrence")
table(cover$rel_freq.space)
summary(cover$ave_rel_abundance_year0 )
hist(cover$ave_rel_abundance_year0, xlab = "Average relative abundance at a site in pre-treatment year", main ="Average relative abundance")

#plot correlation between relative abundance and frequency metrics
plot(cover$rel_freq.space, cover$rel_abundance_year0)

#######################################################################################################################
## Compute the DI per species per sie defined as:  DI = (average relative abundance + relative frequency)/2 ###########
## considering live cover and pre-treatment year  #####################################################################
#######################################################################################################################
cover[, DI := (ave_rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, site_code)]

## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI], xlab = "Dominance indicator (DI)", main = "Dominance indicator metric")

#plot the DI values by site
DIsp <- ggplot(data = cover, aes(x = DI)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "DI") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
DIsp


#######################################################################################################################
####### Make Categorical Variables of Dominance too  - based on the DI metric #########################################
#######################################################################################################################
# Then given each species a ranking into different categories,
# to look at changes in the number of those types of species and the impact on productivity. 
cover[, DIgroup:=cut(DI, breaks=c(-0.001,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant")), by = site_code]

#* take a site level distribution and then 25% of sample 
#* 
#* 
#check the data
table(cover$DIgroup)
table(cover$DIgroup, cover$site_code)

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

#########################################################################################################################
### Create variables that are combined groups of: #######################################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
##############################################################################################################################
# create combinations of all - as factor in a column
cover[, status.NN.RareDom := paste(DIgroup,local_provenance, sep = "_")]
table(cover$status.NN.RareDom)

# 3.Including them all as non-native: 
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

### Create a variable for the unknown SR
cover[, sr_unk_rare := length(unique(Taxon[DIgroup == "Rare" &  local_provenance == "UNK"])), by = .(plot, site_code, year)]

cover[, rare_known := sum(sr_nat_rare + sr_non.nat_rare), by = site_code ]
cover[, frac_unk_rare := sr_unk_rare/rare_known, by = site_code ]

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

# create a non-rare variable for frequency
cover[, non_rare_spp.Freq := Freq_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE"])), by = .(plot, site_code, year)]

# create a rare variable for frequency
cover[, sr_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "FALSE"])), by = .(plot, site_code, year)]
cover[, sr_rare_non.nat.Freq:= length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.Freq := length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.Freq  := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

######################################################################################################################
## Compute for Average Relative Abundance: sub-ordinate and rare - Cut off 1  #######################################
######################################################################################################################
summary(cover$rel_abundance_year0)
hist(cover$ave_rel_abundance_year0)

cover[, RelAbund_group:=cut(ave_rel_abundance_year0, breaks=c(-0.00001,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

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

#################################################################################################################################
#### Compute Groups based on DI Metric and Cut Off 2 for Rare, Dom, Subord Groups  ###############################################
#################################################################################################################################

## Do a different cut off for dominance, subordinate and rare:
cover[, DIgroup2:=cut(DI, breaks=c(-0.00001,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]
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
## Compute for Frequency: sub-ordinate and rare - Cut Off 2 based on breaks=c(0.0,0.4,0.8,1.0). #####################
######################################################################################################################
cover[, Freq_group2:=cut(rel_freq.space, breaks=c(-0.00001,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

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
cover[, RelAbund_group2 :=cut(ave_rel_abundance_year0, breaks=c(-0.00001,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

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

########################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above . 
#**** NOTE FOR DATA-REUSE: If using the full dataset, likely need to check for others since there could be more in the full dataset!!!! 
### but here is OK because this will be merge with the comb data for control plots that I fixed****** 
cover = cover[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

# write as csv datafile to use for R
#with ALL species including ones not present in year 0, which are then considered 0s and thus rare
# in all of the subsequent count and variable calculations
# write.csv(cover, "NutNetCoverData_ProcessedAug2019-AllSpecies.csv")

#create a version of the data running the above code with cover = cover_present_year0  
#to process all variables of SR counts and groupings on only species present in year 0 then printing as:
 write.csv(cover, "NutNetCoverData_ProcessedAug2019-PresentYear0only.csv")
 
 # write.csv(cover, "NutNetCoverData_ProcessedAug2019-2.csv")
 
 
#############################################################################
### References ############################################################
############################################################################

# Avolio, M.L., et al. (2019). Demystifying dominant species. New Phytol., 223, 1106–1126.