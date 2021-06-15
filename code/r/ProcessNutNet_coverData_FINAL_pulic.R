##########################################
### Process NutNet Cover Data ############
## Laura Dee  # June 14 2021      ############
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

## to make a plot of number of introduced species over time at a site 
cover[, sr_INT.site := length(unique(Taxon[local_provenance == "INT" & max_cover>0])), by = .(site_code, year)]


#*** CUT IN BTWN *****
cover.int.unique = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt,  sr_INT.site, sr_INT)])

INTsp <- ggplot(data = cover.int.unique , aes(x = year, y = sr_INT.site) + geom_point() +  # + facet_wrap(~site_code) + theme_bw() +
                  labs(x = "Count of Non-Native Species") +  theme_bw() +
                  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
                  theme(axis.text.y = element_text(size = 14)) + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                  theme(axis.text.x = element_text(size=14)) 
                INTsp
                
invasionmod <- lm(sr_INT.site ~ year_trt, data = cover.int.unique)
                summary(invasionmod)
                
plot(cover.int.unique$year_trt, cover.int.unique$sr_INT.site )
abline(invasionmod)
 #*** CUT IN BTWN / Do for the final dataset*****

                
#####
##### Compute a convenience column that says whether a species was present in a plot in a site in the pre-treatment year
cover[,present_year0:=max_cover[year_trt==0]>0, by=.(Taxon,site_code,plot)]

#compute the SR of species that were not present in year 0 for each plot and year but are present in a given year in that plot
cover[, sr_NA := length(unique(Taxon[present_year0 == FALSE & max_cover>0])), by = .(plot, site_code, year)]

#*** ultimately cut btwn these lines to clean code*****
printNA <- table(cover$sr_NA, cover$site_code)
# summary(cover$sr_NA.test)
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
#*** ultimately cut btwn these lines to clean code*****
#*

##########################################################################################################
##### Compute TOTAL & TOTAL LIVE RELATIVE COVER PER PLOT MEASURES ########################################
##########################################################################################################
# make a total cover in a plot, site, year. This includes live cover only.
cover[,totplotcover.yr.live := sum(max_cover, na.rm= T), by=.(plot, site_code, year)]

#Make anrelative cover for each species in each plot and year
# based on TOTAL cover (including only live cover as we already filtered to the data table to live above).
cover[,relative_sp_cover.yr.live := max_cover/totplotcover.yr.live]
#***to check, run: sum(is.na(cover$relative_sp_cover.yr.live))

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
   # azi.cn.test = cover[site_code == "azi.cn" , ]
# ****azi.cn.test[Taxon=="AGROSTIS STOLONIFERA",]

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
unique.ras[,RAsite_group := ifelse(relative_abundance_spp_site.yr0<RAquant0.6,"Rare",
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
#check for NAs 
table(unique.ras$site_code, unique.ras$RAsite_group, useNA = "ifany")
table(cover_present_year0$site_code, cover_present_year0$RAsite_group, useNA = "ifany")

## saline.us site does not have a pre-treatment year so its dropped

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

#################################################################################################################################
### SR variables by combined groupings based on site relative abundance ###############################################
#################################################################################################################################

##### for the unknown species are will the data processing three ways. #################
# We will run each scenario as sensitivity analysis.  #################################

## 1. Excluding them (as above) #### 
#do SR for non-native, rare:
cover_present_year0[, sr_non.nat_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
 # table(cover_present_year0$site_code, cover_present_year0$sr_non.nat_rare, useNA = "ifany")

#do SR for native, rare:
cover_present_year0[, sr_nat_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
# table(cover_present_year0$site_code, cover_present_year0$sr_nat_rare, useNA = "ifany")

## 2. Including the unknown spp origin all as native: ####
cover_present_year0[, sr_nat_unk_rare := length(unique(Taxon[RAsite_group == "Rare" & local_provenance == "NAT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

# 3.Including them all as non-native: 
cover_present_year0[, sr_non.nat_unk_rare := length(unique(Taxon[RAsite_group== "Rare" & local_provenance == "INT" |  local_provenance == "UNK"])), by = .(plot, site_code, year)]

## do the same for the non-rare variables: 
# 1. Create SR non-rare native and non-native excluding unknown species origin species
cover_present_year0[, sr_non.rare_non.nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
# table(cover_present_year0$site_code, cover_present_year0$sr_non.rare_non.nat, useNA = "ifany")

cover_present_year0[, sr_non.rare_nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
#**check cover_present_year0[site_code=="yarra.au" & plot==8,.(plot, subplot, Taxon, year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]
# table(cover_present_year0$site_code, cover_present_year0$sr_non.rare_nat, useNA = "ifany")

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
### SR variables by combined groupings based on site Relative Frequency ###############################################
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


#####################################################################
### Compute SR variables based on Frequency Groups above : #########
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

######################################################################################################################
## Make a 2nd Cut-off for Relative Abundance Groupings For Sensitivity Analyses 
###################################################################################################################################################
## Make Categorical Variables to Label Spp as Dominant, Subordinant, and Rare   - based on the relative abundance Quantiles per Site ############

#**Note to self ****
#0 quartile = 0 quantile = 0 percentile
# 1 quartile = 0.25 quantile = 25 percentile
# 2 quartile = .5 quantile = 50 percentile (median)
# 3 quartile = .75 quantile = 75 percentile
# 4 quartile = 1 quantile = 100 percentile

unique.ras = unique(cover_present_year0[, .(site_code, Taxon, relative_abundance_spp_site.yr0)])

unique.ras[,RAquant0.7:=quantile(relative_abundance_spp_site.yr0, probs=0.7), by=site_code]
unique.ras[,RAquant0.95:=quantile(relative_abundance_spp_site.yr0, probs=0.95), by=site_code]
unique.ras[,RAsite_group2 :=ifelse(relative_abundance_spp_site.yr0<RAquant0.7,"Rare",
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
table(unique.ras$RAsite_group2)
#whats the breakdown of species classified in each group by site
table(unique.ras$site_code, unique.ras$RAsite_group2)

##########################################################################################################
### Create variables (factor) that are combined groups of:  ###################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
cover_present_year0[, status.NN.RareDom2 := paste(RAsite_group2,local_provenance, sep = "_")]
table(cover_present_year0$status.NN.RareDom2)

# create a non-rare variable
cover_present_year0[, non_rare_spp2 := RAsite_group2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

# From here on out, we don't need/want the extra records that were created in order
# to properly calculate rarity, relative abundance, etc if a species showed up
# in a plot after the first year or existed in some plots within a site.
# Filter back down to (original) records with max_cover > 0 
cover_present_year0 = cover_present_year0[max_cover>0,]
cover = cover[max_cover>0,]

## 1. Excluding the species of unknown origin as in the main analysis #### 
#do SR for non-native, rare:
cover_present_year0[, sr_non.nat_rare2 := length(unique(Taxon[RAsite_group2 == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover_present_year0[, sr_nat_rare2 := length(unique(Taxon[RAsite_group2 == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# Create SR non-rare native and non-native excluding unknown species origin species
cover_present_year0[, sr_non.rare_non.nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.rare_nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
#**check cover_present_year0[site_code=="yarra.au" & plot==8,.(plot, subplot, Taxon, year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]

# include in the final code for data processing "finalprocess_coverdata.R"
  # sr_non.rare_nat2, sr_non.rare_non.nat2 , sr_nat_rare2, sr_non.nat_rare2


######################################################################################################################
## Make a 3rd Cut-off for Relative Abundance Groupings For Sensitivity Analyses 
###################################################################################################################################################
## Make Categorical Variables to Label Spp as Dominant, Subordinant, and Rare   - based on the relative abundance Quantiles per Site ############
unique.ras = unique(cover_present_year0[, .(site_code, Taxon, relative_abundance_spp_site.yr0)])

unique.ras[,RAquant0.5:=quantile(relative_abundance_spp_site.yr0, probs=0.5), by=site_code]
unique.ras[,RAquant0.95:=quantile(relative_abundance_spp_site.yr0, probs=0.95), by=site_code]
unique.ras[,RAsite_group3 :=ifelse(relative_abundance_spp_site.yr0<RAquant0.5,"Rare",
                                   ifelse(relative_abundance_spp_site.yr0<RAquant0.95, "Subordinate","Dominant"))]
#re-merge the quantiles and classifications into the cover_present_year0 dataset
unique.ras[,relative_abundance_spp_site.yr0:=NULL] # drop before re-merge
cover_present_year0 = merge(cover_present_year0, unique.ras, by=c("site_code", "Taxon"))

#sanity check
cdcr = cover_present_year0[site_code =="cdcr.us",]
table(cdcr$Taxon, cdcr$RAsite_group3)
konz  =  cover_present_year0[site_code =="konz.us",]
table(konz$Taxon, konz$RAsite_group3)

#whats the breakdown of species classified in each group overall 
table(unique.ras$RAsite_group3)
#whats the breakdown of species classified in each group by site
table(unique.ras$site_code, unique.ras$RAsite_group3)

##########################################################################################################
### Create variables (factor) that are combined groups of:  ###################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
cover_present_year0[, status.NN.RareDom3 := paste(RAsite_group3,local_provenance, sep = "_")]
table(cover_present_year0$status.NN.RareDom3)

# create a non-rare variable
cover_present_year0[, non_rare_spp3 := RAsite_group3 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

# From here on out, we don't need/want the extra records that were created in order
# to properly calculate rarity, relative abundance, etc if a species showed up
# in a plot after the first year or existed in some plots within a site.
# Filter back down to (original) records with max_cover > 0 
cover_present_year0 = cover_present_year0[max_cover>0,]
cover = cover[max_cover>0,]

## 1. Excluding the species of unknown origin as in the main analysis #### 
#do SR for non-native, rare:
cover_present_year0[, sr_non.nat_rare3 := length(unique(Taxon[RAsite_group3 == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover_present_year0[, sr_nat_rare3 := length(unique(Taxon[RAsite_group3 == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# Create SR non-rare native and non-native excluding unknown species origin species
cover_present_year0[, sr_non.rare_non.nat3 := length(unique(Taxon[non_rare_spp3 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover_present_year0[, sr_non.rare_nat3 := length(unique(Taxon[non_rare_spp3 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
#**check cover_present_year0[site_code=="yarra.au" & plot==8,.(plot, subplot, Taxon, year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]

# include in the final code for data processing in "finalprocess_coverdata.R"
# sr_non.rare_nat3, sr_non.rare_non.nat3 , sr_nat_rare3, sr_non.nat_rare3

#####################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above . 
#**** NOTE FOR DATA-REUSE: If using the full dataset, likely need to check for others since there could be more in the full dataset!!!! 
### but here is OK because this will be merge with the comb data for control plots that I fixed****** 
cover_present_year0 = cover_present_year0[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

#create a version of the data running the above code with cover = cover_present_year0  
#to process all variables of SR counts and groupings on only species present in year 0 then printing as:
write.csv(cover_present_year0, "NutNetCoverData_ProcessedFinal.csv")


################################################################################################################################
## Examine missing cover data for sier.us, kiny.au, and mcla.us for Control plots by years #######################################
##################################################################################################################################
# to then compare with coverage of observations of plots/years with comb data
sier.cover = cover_present_year0[site_code =="sier.us" & trt == "Control", ]
table(sier.cover$plot, sier.cover$year)

kiny.cover = cover_present_year0[site_code == "kiny.au" & trt == "Control", ]
table(kiny.cover$plot, kiny.cover$year)

mcla.cover = cover_present_year0[site_code == "mcla.us" & trt == "Control", ]
table(mcla.cover$plot, mcla.cover$year)

