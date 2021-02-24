# Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-09-Apr-2018.csv'

length(comb$live_mass) #74 NAs for live mass
summary(comb$live_mass)
length(comb$rich)
summary(comb$rich) # 2 NAs for richness

#################################################################
### Process data to prep to use it in the models ###############
###############################################################
comb$site <- comb$site_code

# Filter data to records with non-NA live_mass and non-NA richness.
# Models can generally handle NAs by implicitly dropping, but clustered SE
# function currently does not.
comb = comb[!is.na(live_mass) & !is.na(rich)]

#count for rich and live_mass
length(comb$rich)
# 1252
length(comb$live_mass)
# 1252
# Determine the number of Obs.
nrow(comb) #1252

## Confirm that only control plots are in the data
table(comb$trt)
## Confirm the # of years 
table(comb$year)

## Check to make sure that site and year are factors/characters
class(comb$year)
class(comb$site_code) # 43
class(comb$newplotid) #172

## 
length(unique(comb$newplotid)) #172

#see if there are some plots with only one year - one control observation:
comb[,.N, by=c("newplotid")][N==1,]  # there are 21 that need to be removed.

#Create a new version of comb filtering out singletons- plots with only one year that is one control observation:
# flag the singleton
comb[,singleton:=(.N==1), by=c("newplotid")]
comb[singleton==T, .(newplotid)]  # these are the singletons

#remove the singletons, select rows that are NOT singletons 
comb = comb[singleton == F, ]

# Determine the number of Obs. removing the singletons and obs with NA or rich or live_mass
nrow(comb) #1231 
length(unique(comb$newplotid)) #151

##############################################################################
### Make dummy variables for the panel regression  analyses #################
##############################################################################

# make year a character, to be a dummy variable: 
comb$year <- as.character(comb$year)

#make a factor that is site by year
comb[, site.by.yeardummy := paste(site_code, year, sep = "_")]


#########################################################################################################################################
###Print Table of Data ##########################################################################################################
##############################################################################################################################################

# Write out what is covered in this dataset for these models with dataset version 1:
# for a SI Table: #print which sites and years are in this version of the filtered dataset.
comb.descript.v1 =  table(comb$site_name, comb$year)
# to view: 
# table(comb$site_code, comb$year)
write.csv(comb.descript.v1, "./output/Table_S1.csv")
# length(unique(comb$site_code))