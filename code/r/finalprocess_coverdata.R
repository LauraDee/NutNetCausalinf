####################################################################################################
## Filter cover data to only control plots & prep for merge and models ###############################
##################################################################################################
#comb = comb[is.control==T & has.5.yrs.data==T,]  #this should be done already in the saved file from "NutNetControlPlotDataToUseApril2018.csv"

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
cover$V1 = NULL

#################################################################################################
## Grab the summary columns from cover to merge  #######################################################
################################################################################################

#to subset columns and also remove duplicate rows from the cover file so that there is one observation per plot and year 
# and the data isn't artificially replicated 
# coversummaries = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt, totplotcover.yr.live, sr_NA,
#                                  sr_INT, sr_NAT, sr_domspp, sr_rarespp, sr_subordspp, sr_non_rare_spp, 
#                                   sr_non.nat_rare,  sr_nat_rare, sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom, sr_non.nat_dom, 
#                                 #  sr_non_rare_spp.RelA, 
#                                   sr_non_rare_spp.Freq, sr_non.rare_nat.Freq, sr_non.rare_non.nat.Freq,
#                                   sr_rare_non.nat.Freq, sr_rare_nat.Freq,
#                                   sr_domspp2, sr_rarespp2 , sr_subordspp2, non_rare_spp.DI2,
#                                   sr_nat_sub.Freq, sr_non.nat_sub.Freq, 
#                                   sr_non.rare_nat.RelA, sr_non.rare_non.nat.RelA, sr_rare_non.nat.RelA, sr_rare_nat.RelA,
#                                   sr_non.rare_nat.Freq2, sr_non.rare_non.nat.Freq2, sr_rare_non.nat.Freq2, sr_rare_nat.Freq2,
#                                   sr_non.rare_nat.RelA2, sr_non.rare_non.nat.RelA2, sr_rare_non.nat.RelA2, sr_rare_nat.RelA2, # non_rare_spp2,
#                                   sr_non.rare_nat2, sr_non.rare_non.nat2, sr_non.nat_rare2, sr_nat_rare2, 
#                                   # sr_rarespp2, 
#                                   #rare_spp.DI2,
#                                   freq_sr_domspp, freq_sr_rarespp, freq_sr_subordspp, 
#                                   sr_nat_unk_rare, ## 2. Including the unknown spp origin all as native: ####
#                                   sr_non.nat_unk_rare, # 3.Including them all as non-native: 
#                                   sr_non.rare_nat_unk, ## 2. Include the unknown spp origin all as native: ####
#                                   sr_non.rare_non.nat_unk # 3.Including them all as non-native: 
# )])
# #make sure number of rows isn't inflated 
# nrow(coversummaries)


#do for the new analyses with relative abundance in "NutNetCoverData_ProcessedFinal.csv"
coversummaries = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt, totplotcover.yr.live, sr_NA,
                                  sr_INT, sr_NAT, sr_UNK, sr_INT.site, 
                                  sr_non.nat_rare,  sr_nat_rare, sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom, sr_non.nat_dom, 
                                  sr_nat_unk_rare, ## 2. Including the unknown spp origin all as native: ####
                                  sr_non.nat_unk_rare, # 3.Including them all as non-native: 
                                  sr_non.rare_nat_unk, ## 2. Include the unknown spp origin all as native: ####
                                  sr_non.rare_non.nat_unk , # 3.Including them all as non-native
                                  sr_non_rare_spp.Freq, sr_non.rare_nat.Freq, sr_non.rare_non.nat.Freq,
                                  sr_rare_non.nat.Freq, sr_rare_nat.Freq,
                                  sr_rare_unk_nat.Freq , sr_non.rare_nat_unk.Freq, ## 2. Including the unknown spp origin all as native: ####
                                  sr_non.rare_non.nat_unk.Freq, sr_rare_non.nat_unk.Freq,   # 3.Including them all as non-native
                                  sr_non.rare_nat2, sr_non.rare_non.nat2 , sr_nat_rare2, sr_non.nat_rare2 , #include variables for cut-off 2 as sensitivity test for main model for Figure 5
                                  sr_non.rare_nat3, sr_non.rare_non.nat3 , sr_nat_rare3, sr_non.nat_rare3 #include variables for cut-off 2 as sensitvity test for main model for Figure 5
                                  )])
nrow(coversummaries)


#################################################################################################
## Merge comb with Processed Cover Data #######################################################
################################################################################################
# Merge the cover data with the comb data. Processed the April 2018 versions of the datasets.
# merge to keep the data that is in the comb file (i.e., sites with at least 5 years of data for control plots)
mech.data = merge(comb, coversummaries, by=c("site_code","plot","year"), all.x=T)

#what I had been doing: mech.data = merge(comb, cover, by=c("site_code","plot","year"), all.x=T)
nrow(mech.data) #should be 1231
dim(mech.data)
nrow(comb)
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


########################################################################################################################
##### Get some summary numbers on this final dataset for the counts of different groups of species ####################
########################################################################################################################

#whats the breakdown of species classified in each group overall 
table(unique.ras$RAsite_group2)
#whats the breakdown of species classified in each group by site
table(unique.ras$site_code, unique.ras$RAsite_group2)




