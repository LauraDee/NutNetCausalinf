##############################################################################
# Process NutNet Data for biodiversity-productivity analyses on control plots #
# LDee April 27 2018 , May 11 2020 ############################################################
#############################################################################
# updated May 11 2020 to add:  #an IV of average neighbor TREATED richness at the site
     #variables identifying if richness change is  positive, negative, or no change

# Meta-data files are attached, but here are some key variables used:
      # Field	column type	description
      # site_name	varchar(45)	Expanded site name
      # site_code	varchar(45)	code in the form of "shorthand.country" e.g., cdcr.us, kiny.au, frue.ch
      # continent	varchar(45)	site continent
      # region	varchar(45)	region 
      # precipitation	double	Mean annual precip in mm
      # elevation	double	mean elevation (m)
      # latitude	double	latitude from -90 (S) to +90 (N) in decimal degrees
      # longitude	double	longitude from -180 (W) to 180 (E) in decimal degrees
      # plot	int(11)	plot number
      # N	tinyint(1)	1=N added, 0 = control
      # P	tinyint(1)	1=P added, 0 = control
      # K	tinyint(1)	1=K added, 0 = control
      # Exclose	tinyint(1)	1=fenced, 0=control
      # year_trt	decimal(6,0)	0 = pretreatment observation; 1=first observations after treatments; 2=2nd year after treatment; etc
      # rich	decimal(24,4)	COUNT of unique taxa in the plot
      # site_year_rich	bigint(21)	COUNT of unique taxa observed across all plots at the site in that year
      # site_richness	bigint(21)	COUNT of unique taxa ever observed across all plots in all years at that site
      # dead_mass	double	SUM of dead plant biomass/litter per plot (g)
      # live_mass	double	SUM of live plant biomass per plot (g)
      # total_mass	double	SUM of biomass per plot (g)
      # proportion_par	double	Ground_PAR / Ambient_PAR

##  Notes on variables:
  #  Ground_PAR varies at plot-yr level (but 1/3 of obs missing)
  #  ppm_P and ppm_K vary at plot level, but not across time

#Close graphics and clear local memory
graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(sandwich)
library(foreign)

# Load plot level plant, soil, and climate data
setwd("~/Dropbox/IV in ecology/NutNet")
## load latest data file
comb <- fread('comb-by-plot-clim-soil-diversity-09-Apr-2018.csv',na.strings= 'NA')
comb$site <- comb$site_code
comb$plot <- as.factor(comb$plot)

##############################################################################################
### Process Data To be Used for Analyses #####################################################
###############################################################################################
#Get maximum number of treatment years per site
comb[,max.trt.yr:=max(year_trt), by=.(site)]
#Find min value of treatment years; 0 = pre-treatment data
comb[,min.trt.yr:=min(year_trt), by=.(site)]

## Make columns that flag different observations (# of years, etc) to filter to run models on different data subsets
## has AT LEAST 5 years of data (including pre-treatment year = 0)
# another way to do this:
#kMinNumYears = 5
#comb <- comb[comb$min.trt.yr == 0 & comb$max.trt.yr > kMinNumYears-1,]

#comb[, has.5.yrs.data:=(max.trt.yr - min.trt.yr >= 4)]  # a short cut, if we didnt care about if min.trt year = 0 or 1.
comb[, has.5.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 4)|(min.trt.yr==1 & max.trt.yr>=5)] #this keeps control plots without a pre-treatment survey
# comb[, has.5.yrs.data:=(min.trt.yr == 0 & max.trt.yr > 4)] #without the min.trt.year==1 
comb[, has.4.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 3)|(min.trt.yr==1 & max.trt.yr>=4)] #this keeps control plots without a pre-treatment survey]
comb[, has.3.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 2)|(min.trt.yr==1 & max.trt.yr>=3)]
comb[, has.6.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 5)|(min.trt.yr==1 & max.trt.yr>=6)]

##############################################################################################################################
### Filter to sites with at least 5 years of data ############################################################################
##############################################################################################################################
comb = comb[has.5.yrs.data==T,]

#############################################################################################################
#### Identify THE CONTROL PLOTS ############################################################################
##########################################################################################################

# Mark which records are for control plots <-- check to make sure I have this right bc looks like there are more treatmnets
# trt groups are: Control, Fence, K, N, NK, NP, NPK, NPK+Fence, P, PK 
comb[ ,is.control:=(N==0 & P==0 & K==0 & Exclose==0),]

###################################################################################################
#### Use to filter pre-treatment years or treated years ############################################
#################################################################################################### 
# Make a variable that identifies if it is a pre-treatment year or not (should be 2951 obs with year = 0 that are TRUE; confirmed this)
comb[,is.PretreatmentYr := (year_trt == 0),]

#make this also a variable that is a 0 or 1 to indicate 1 (TRUE) if it is a treated year or not for processing the full dataset
# to only sites that also have pre-treatment years. 
comb[,is.TreatedYear := (is.PretreatmentYr != TRUE),]

##############################################################################################################
## Identify if fenced or not for mechanisms analysis ###########################################################
##########################################################################################################
# Make a variable that identifies if a plot is fenced or not.
# Oh, that already existed: Exclose:	1=fenced, 0=control
comb[,is.Fenced := Exclose != 0, ] 

##################################################################################################
### Making a unique plot id and year as factor ####################################################
##################################################################################################
#in general: my.dt[,newplotid:=as.factor(paste(site.id.col, plot.id.col, sep="_"))]
comb[,newplotid:=as.factor(paste(site_code, plot, sep="_"))]

#################################################################################################
### Making the block id as factor & with a unique id #############################################
##################################################################################################
comb$block = as.factor(comb$block)
# this makes it a unique block-site combination
comb[,newblockid:=as.factor(paste(site_code, block, sep="_"))]

################################################################################################
### Create instrument: average neighbor richness of the site #################################################
#####################################################################################################
# To calculate average of plot-level richness of all other plots in the same site and year, we'll group 
# by site and year, add up the richness across all rows in the group, remove own richness, and divide 
# by the number of items in the group minus one
comb[, avg_neighbor_rich :=(sum(rich, na.rm=T) - rich)/(.N-1) , by=.(site_code, year)]
plot(comb[,.(avg_neighbor_rich, rich)])

  # spot check this with an example
  comb[site_code=="bnch.us" & year==2007, .(site_code, year, plot, rich, avg_neighbor_rich),]
  # mean(c(9,12,9,10,8))
  # mean(c(11,12,9,10,8))

##################################################################################################################
### Create a lagged average neighborhood richness of the site  ###################################################
##################################################################################################################
comb[order(year), lagged_avg_neighbor_rich := shift(avg_neighbor_rich), by =.(plot, site_code)]

#####################################################################################################################
### Create average neighbor TREATED richness within the block IV #########################################################
######################################################################################################################
# do some blocks have multiple control plots? yes:
# whats max number of controls per block x year
max(comb[,sum(trt=="Control"), by=.(newblockid, year)][,3, with=F])
#gives counts of controls per block x year for those with multiple controls per block x year
comb[,sum(trt=="Control"), by=.(newblockid, year)][V1>1,]

## take the average neighbor richness of all the plots within block i, within site j, that are NOT controls: 
comb[,avg.trt.neigh.rich.within.block:=mean(rich[trt!="Control"]),by=.(newblockid, year)]

#create a lag
comb[order(year), lagged_avg.trt.neigh.rich.within.block := shift(avg.trt.neigh.rich.within.block), by =.(plot, site_code)]

#####################################################################################################################
### Create IV: average neighbor TREATED richness AT THE SITE LEVEL IV #################################################
######################################################################################################################
## take the average neighbor richness of all the plots within a site j that are NOT control plots (i.e. the treated richness): 
comb[,avg.trt.neigh.rich.within.site:= mean(rich[trt!="Control"], na.rm = T),by=.(site_code, year)]

#create a lag
comb[order(year), lagged_avg.trt.neigh.rich.within.site := shift(avg.trt.neigh.rich.within.site), by =.(plot, site_code)]

###################################################################################################################
### Make a lagged, plot-level richness & evenness term of own plot ###############################################
####################################################################################################################
# creating a lag with data table: https://rdrr.io/cran/data.table/man/shift.html
# shift(x, n=1L, fill=NA, type=c("lag", "lead"), give.names=FALSE)
# be sure to be sure to sort by year before doing shift with "setorder" fct
#sytax:DT[order(year), (cols) := shift(.SD, 1, type="lag"), .SDcols=cols]
comb[order(year), laggedrich := shift(rich), by =.(plot, site_code)]

############################################################################################
### Make a lagged, plot-level productivity term ###############################################
##############################################################################################
comb[order(year), laggedlive_mass := shift(live_mass), by =.(plot, site_code)]

################################################################################################
### Create a first-difference ########################################################################
#################################################################################################
#as for changes you can use the shift function in data.table; but things need to be sorted first
# if you want changes, rows need to be sorted by time. Sort your data table by whatever the time column is. then use the shift function.
#if you already computed laggedrich, you don't need the by part or the order part. 
# they are only needed for constructing the lags   #comb[, changerich := rich-laggedrich,]
comb[order(year), changerich := rich-shift(rich), by =.(plot, site_code)]
comb[order(year), changelive_mass := live_mass-shift(live_mass), by =.(plot, site_code)]

## Make a lag of the differences ##
comb[order(year), lagged_changelive_mass := shift(changelive_mass), by =.(plot, site_code)]
comb[order(year), lagged_changerich:= shift(changerich), by =.(plot, site_code)]

#summary(comb$changerich)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
  #-16.0000  -2.0000   0.0000  -0.2395   2.0000  20.0000      250 

################################################################################################
### Create variables that identify if the richness change is increase, decrease or no change  ##############
#################################################################################################
##### Whether richness chane is an increase or decrease #######
comb[, rich_increase := changerich>0]
comb[, rich_decrease := changerich<0]
comb[, rich_nochange := changerich == 0] 

#####################################################################################################################
## Also create a differenced variable for the PAR control variables ########################################
#####################################################################################################################
comb[order(year), changeGround_PAR := Ground_PAR -shift(Ground_PAR), by =.(plot, site_code)]
comb[order(year), changeAmbient_PAR := Ambient_PAR-shift(Ambient_PAR), by =.(plot, site_code)]

###############################################################################################################
### Make a lagged, plot-level evenness term & for other BD terms  ###############################################
################################################################################################################
# creating a lag with data table: https://rdrr.io/cran/data.table/man/shift.html
# shift(x, n=1L, fill=NA, type=c("lag", "lead"), give.names=FALSE)
# be sure to be sure to sort by year before doing shift with "setorder" fct
#sytax:DT[order(year), (cols) := shift(.SD, 1, type="lag"), .SDcols=cols]

# do for evenness
comb[order(year), laggedeven := shift(even), by =.(plot, site_code)]
# do for shannon's
comb[order(year), laggedshan := shift(shan), by =.(plot, site_code)]
#do for simpson's
comb[order(year), laggedsimpson := shift(simpson), by =.(plot, site_code)]

#####################################################################################################################
## also create a differenced variable for evenness & other BD variables in comb version april 2018 ###############################################
#####################################################################################################################
comb[order(year), changeEvenness := even-shift(even), by =.(plot, site_code)]
# do for shannon's
comb[order(year), changeShan := shan-shift(shan), by =.(plot, site_code)]
#do for simpson's
comb[order(year), changeSimpson := simpson-shift(simpson), by =.(plot, site_code)]

############################################################################################################
## Compute site-level productivity per yr and average (added June 29, 2018) ###############################
##########################################################################################################

# sum of live mass per year, summed over plots by year
comb[ ,site_live_mass.yr:= sum(live_mass, na.rm=T), by=.(site_code, year)]

#average site-level productivity, averaged across all years
comb[, ave_site_live_mass := mean(site_live_mass.yr), by = site_code]

############################################################################################################
## Compute site-level richness variables  (added May 11 2020 ) ###############################
##########################################################################################################

#alt way to do it:  comb[ , initial_site_rich := site_year_rich[is.PretreatmentYr == TRUE], by=.(site_code)]
comb[ , initial_site_rich := site_year_rich[year_trt == 0], by=.(site_code)]

site_rich_trend <- lm(site_year_rich ~ year_trt:site_code, data = comb)
summary(site_rich_trend)

#######################################################################################################
#### Filter Data by years and control or not  ########################################################
#####################################################################################################

## filter to control plots with 5 years of data
comb = comb[is.control==T & has.5.yrs.data==T,]

#confirm this worked and only control plots remain
table(comb$trt)

## the filter missed two sites that only have 4 years of data - but over the span of 5 years 
# Azi is missing live_mass data in 2012. Barta.us is missing data for all variables 2010, so only has 4 years of data.
# remove these 2 sites from the data with the filter that AT LEAST 5 years of data is needed for the site to be included:
comb = comb[site_code != "azi.cn",]
comb = comb[site_code != "barta.us",]

# check to make sure these are removed
list(unique(comb$site_code))  # list of remaining site code names
length(unique(comb$site_code))   # with this filter,  we have 43 sites. 

## Found an error in the dataset for in site "comp.pt", plots 5, 19, and 34, years 2013-2016 are repeated. 
#They are identical except for their indicators of treatment (year_trt). 
# Those years have a duplicate that says year_trt = 0, which isn't true since looks like the year 0 was 2012. 
# Need to remove those plots that are mistakes -
comb = comb[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

#check for duplicates
comb[,.N, by=c("site_code", "plot", "year")][N>1,]
# Empty data.table (0 rows and 4 cols): site_code,plot,year,N
table(comb[,.N, by=c("site_code", "plot", "year")][,N])

comb[,.N, by=c("site_code", "newplotid", "year")][N>1,]
# Empty data.table (0 rows and 4 cols): site_code,newplotid,year,N

comb[,.N, by=c("site_code", "newplotid", "year_trt")][N>1,]
# Empty data.table (0 rows and 4 cols): site_code,newplotid,year_trt,N

comb[,.N, by=c("site_code", "plot", "year_trt")][N>1,]
#Empty data.table (0 rows and 4 cols): site_code,plot,year_trt,N

# compare to Paul's STATA file to make sure duplicates are gone

#pf.comb2 = fread("~/Dropbox/Nutnet causality paper draft/Ferraro/NutNetControlPlotDataToUseApril2018.csv")
# pf.comb2[,.N, by=c("site_code", "plot", "year_trt")][N>1,]
# Empty data.table (0 rows and 4 cols): site_code,plot,year_trt,N

######################################################################
## Write Out files as .csv and .dta for STATA and R #################
####################################################################

## see which sites and years & write out the site and year list as a table: 
tab =  table(comb$site_name, comb$year)
write.csv(tab, "DatasetDescript-ControlPlotsSiteYearList.csv")

# write as csv datafile to use for R
write.csv(comb, "NutNetControlPlotDataToUseApril2018.csv")


