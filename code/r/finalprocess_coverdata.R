####################################################################################################
## Filter cover data to only control plots & prep for merge and models ###############################
##################################################################################################
# The data use policy of the Nutrient Network is to publish the minimum dataset required to produce the analyses, from the derived data. 
# This code is used to filter to the data to meet their requirements, though we also provide the full code to show all steps from our data processing from the raw data for reproducibility.
# The raw data for unmanipulated plots that were not included in the analyses, because they did not meet the inclusion criteria, 
# are available under restricted access for which permission can be obtained by contacting the Nutrient Network at https://nutnet.org. 

# Updated by Laura Dee Feb 9 2023

#comb = comb[is.control==T & has.5.yrs.data==T,]  #this should be done already in the saved file from "NutNetControlPlotData_derived.csv"

### Confirm that only control plots are in the data
table(comb$trt)
# ## Confirm the # of years 
table(comb$year)

cover = cover[trt == "Control",]

### Confirm that only control plots are in the data
table(cover$trt)
### Confirm the # of years 
table(cover$year)

###############################################################################
#### Filter Data by years   ################################################
################################################################################ 

# in cover - filter to control plots with 5 years of data to match the main analyses
#Get maximum number of treatment years per site
cover[,max.trt.yr:=max(year_trt), by=.(site_code)]
#Find min value of treatment years; 0 = pre-treatment data
cover[,min.trt.yr:=min(year_trt), by=.(site_code)]

#filter to those sites and years
cover[, has.5.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 4)|(min.trt.yr==1 & max.trt.yr>=5)] #this keeps control plots without a pre-treatment survey
cover = cover[has.5.yrs.data==T,]

## the filter missed two sites that only have 4 years of data - but over the span of 5 years 
# Azi is missing live_mass data in 2012. Barta.us is missing data for all variables 2010, so only has 4 years of data.
# remove these 2 sites from the data with the filter that AT LEAST 5 years of data is needed for the site to be included:
cover = cover[site_code != "azi.cn",]
cover = cover[site_code != "barta.us",]

##########################################################################################
### Process comb data to prep to use it in the models & merge #############################
############################################################################################
## make year a character in both datasets, to be a dummy variable: 
comb$year <- as.character(comb$year)
cover$year <- as.character(cover$year)

# same with plot
comb$plot <- as.character(comb$plot)
cover$plot <- as.character(cover$plot)
cover$V1 = NULL
comb$V1 = NULL

#################################################################################################
## Merge comb with Processed Cover Data #######################################################
################################################################################################
# Merge the cover data with the comb data. Processed the April 2018 versions of the datasets.
# merge to keep the data that is in the comb file (i.e., sites with at least 5 years of data for control plots)
mech.data = merge(comb, cover, by=c("site_code","plot","year"), all.x=T)

#what I had been doing: mech.data = merge(comb, cover, by=c("site_code","plot","year"), all.x=T)
nrow(mech.data) #should be 1231
dim(mech.data)
nrow(comb) #1231
nrow(cover)

##########################################################################################
### Process merge data to prep to use it in the models  #################################
##########################################################################################
# make year a character, to be a dummy variable: 
mech.data$year <- as.character(mech.data$year)
## *Check to make sure that site and year are factors/characters* important*
class(mech.data$year)
class(mech.data$site_code)
class(mech.data$newplotid)

# Filter data to records with non-NA live_mass and non-NA richness.
mech.data = mech.data[!is.na(live_mass) & !is.na(rich)]

## Confirm that only control plots are in the data
table(mech.data$trt.x)

## Confirm the # of years 
table(mech.data$year)

## Confirm its the same site list
list(unique(mech.data$site_code))

#make a factor that is site by year
mech.data[, site.by.yeardummy := paste(site_code, year, sep = "_")]

# write out for STATA
write.csv(mech.data, "ProcessedMechanismAnalysisData.csv")

########################################################################################################################
##### Get some summary numbers on this final dataset for the counts of different groups of species ####################
########################################################################################################################
#the analyses in "analysis_fig5_smsection7.R" are dropping 56 observations. 24 are from saline.us as this site does not have a pre-treatment year
# and thus precludes use of our grouping procedure. What are th sources of the other 32 dropped obs?

# see which sites have NAs for these variables
table(mech.data$site_code, mech.data$sr_non.nat_rare, useNA = "ifany") 
table(mech.data$site_code, mech.data$sr_nat_rare, useNA = "ifany") 
table(mech.data$site_code, mech.data$sr_non.rare_non.nat, useNA = "ifany") 
table(mech.data$site_code, mech.data$ sr_non.rare_nat, useNA = "ifany") 

#check  the sites other than saline.us which have NAs and why: 
#sier.us ecploration of Nas
mech.data[site_code =="sier.us",.(plot,  year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]

sier.comb = comb[site_code =="sier.us" & trt == "Control", ]
table(sier.comb$plot, sier.comb$year)

#kiny.au exploration of NAs
mech.data[site_code =="kiny.au",.(plot,  year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]
#compare what is available for overall richness and biomass versus cover of particular Taxon:
kiny.comb = comb[site_code =="kiny.au" & trt == "Control", ]
table(kiny.comb$plot , kiny.comb$year)
kiny.cover = cover[site_code =="kiny.au" & trt == "Control", ]
table(kiny.cover$plot, kiny.cover$year)


#mcla.us exploration of NAs
mech.data[site_code =="mcla.us",.(plot,  year, sr_nat_rare, sr_non.nat_rare, sr_nat_unk_rare, sr_non.nat_unk_rare, sr_non.rare_non.nat, sr_non.rare_nat)]
mcla.comb = comb[site_code =="mcla.us" & trt == "Control", ]
table(mcla.comb$plot , mcla.comb$year)

mcla.cover = cover[site_code == "mcla.us" & trt == "Control", ]
table(mcla.cover$plot, mcla.cover$year)

#####################################################
###### Plot Data for Figure S11 #####################
#####################################################
par(mfrow=c(1,2))
#A
plot( mech.data$rich, mech.data$sr_nat_rare, pch=19, xlab = "Species richness (plot and year)", ylab = "Native rare species richness (plot and year)", main = "(A) Higher SR is positively associated with native rare SR")
#B
plot( mech.data$rich, mech.data$sr_INT, pch=19, xlab = "Species richness (plot and year)", ylab = "Non-native species richness (plot and year)", main = "(B) Higher SR is positively associated with non-native SR")

