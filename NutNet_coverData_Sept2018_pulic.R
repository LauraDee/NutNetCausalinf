##########################################
### Process Nut Net Cover Data ############
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
library(AER)
library(sandwich)
library(foreign)
library(rmarkdown)

setwd("~/Dropbox/IV in ecology/NutNet")
cover <- fread('full-cover-09-April-2018.csv',na.strings= 'NA')

## need to make max_cover NOT a character
cover$max_cover <- as.numeric(cover$max_cover)

#** CUT ***
# data.table cheat sheets
# cover[, N_fixer_cover.yr := sum(relative_sp_cover.yr, by )]
# comb[,sum(trt=="Control"), by=.(newblockid, year)][V1>1,]

#*** dont filter to live cover only bc we need to know absolute .... or filter to live only????
##*** CUT ****

###### If want to filter cover dataset to only the live cover ################# 
# cover <- cover[live == 1,]
# table(cover$live)

#########################################################################################
## Compare Taxon in Live (live==1) or non-live (live == 0) ################################
##########################################################################################
## determine how many Taxon are listed as live==0 ##
table(cover$live==1)
a <- table(cover$Taxon, cover$live==0)
head(a)
# write.csv(a, "live_v_dead_spplist_NutNet.csv")

# It looks like these categories are primarily the "live == 0 " ones;
# seems like I should NOT include them inthe SR counts?
        # OTHER ANIMAL DIGGING                            0   47
        # OTHER ANIMAL DIGGINGS                           0   68
        # OTHER ANIMAL DROPPINGS                          0   75
        # OTHER ANT NEST                                  0    1
        # OTHER CRUST                                     0   13
        # OTHER DISTURBED SOIL                            0   16
        # OTHER LIGNOTUBER                                0   11
        # OTHER LITTER                                    0 1542
        # OTHER ROCK                                      0  175
        # OTHER SOIL BIOCRUST                             0    5
        # OTHER STANDING DEAD                             0    2
        # OTHER TRIODIA BASEDOWII (DEAD)                  0    5
        # OTHER UNKNOWN SOIL_CRUST                        0     17
        # OTHER WOOD                                      0     27
        #       FUNGI                                       	0   	2
        #       FUNGI SP.                                   	0	    2
        #       GROUND                                      	0	   1039
        #       OTHER WOODY OVERSTORY                       	0	    4
        

##### TOTAL & TOTAL LIVE COVER MEASURES ########################################

# make a total cover in a plot, site, year
cover[,totplotcover.yr := sum(max_cover, na.rm=T), by=.(plot, site_code, year)]
# same for live species only
cover[,totplotcover.yr.live := sum(max_cover[live=="1"], na.rm=T), by=.(plot, site_code, year)]

#Make a relative cover for each species in each plot, site, year
# based on TOTAL cover (including dead cover)
cover[,relative_sp_cover.yr := max_cover/totplotcover.yr]
# same for live species only
cover[,relative_sp_cover.yr.live := (max_cover*(live=="1"))/totplotcover.yr.live]
# in some cases, no live species in plot and year, so getting NA since totplotcover.yr.live is zero. Set these to zero.
cover[is.na(relative_sp_cover.yr.live),relative_sp_cover.yr.live:=0]

## FILTER TO LIVE TO COMPUTE SR MEASURES???? ********

######################################################################################################
### Native vs Non-Native Variables ###################################################################
#######################################################################################################

# covert Naturalised to INT (This use of Naturalised as a category is from one site)
cover[local_provenance=="Naturalised", local_provenance:="INT"]
# convert blank entries to NULL
cover[local_provenance=="", local_provenance:="NULL"]

# Compute native and non-native richness by plot, site, year
cover[, sr_INT := length(unique(Taxon[local_provenance=="INT"])), by = .(plot, site_code, year)]
cover[, sr_NAT := length(unique(Taxon[local_provenance=="NAT"])), by = .(plot, site_code, year)]
cover[, sr_NULL := length(unique(Taxon[local_provenance=="NULL"])), by = .(plot, site_code, year)]
cover[, sr_UNK := length(unique(Taxon[local_provenance=="UNK"])), by = .(plot, site_code, year)]

## do first differences (for the marginal efcomb[order(year), changeGround_PAR := Ground_PAR -shift(Ground_PAR), by =.(plot, site_code)]
cover[order(year), change_sr_INT := sr_INT-shift(sr_INT), by =.(plot, site_code)]
cover[order(year), change_sr_NAT := sr_NAT-shift(sr_NAT), by =.(plot, site_code)]

cover[order(year), change_sr_NULL := sr_NULL-shift(sr_NULL), by =.(plot, site_code)]
cover[order(year), change_sr_UNK := sr_NULL-shift(sr_UNK), by =.(plot, site_code)]

# plot them
hist(cover$change_sr_INT)
hist(cover$change_sr_NAT)
hist(cover$change_sr_UNK)
hist(cover$change_sr_NULL)

## Percent cover of introduced species (non-native to the site) per plot and year  ## 
## compute total cover of non-natives # local_provenance == INT
cover[, NonNative_cover.yr := sum(relative_sp_cover.yr[local_provenance=="INT"]), by = .(plot, site_code, year)]
head(cover)

  #spot check with example
    arch.nonnativecover = cover[site_code == "arch.us", .(year, plot, NonNative_cover.yr),]
  # sites we know have a lot of invasives    
    cover[site_code == "ping.au", .(year, plot, NonNative_cover.yr),]
    cover[site_code == "pinj.au", .(year, plot, NonNative_cover.yr),]
    cover[site_code == "smith.us", .(year, plot, NonNative_cover.yr),]

## Compute Native species cover 
cover[, Native_cover.yr := sum(relative_sp_cover.yr[local_provenance=="NAT"]), by = .(plot, site_code, year)]

##########################################################################################################
### Annual vs Perennial Variables ########################################################################
###########################################################################################################

### Now do for an annual and perennial cover
# Compute richness in annuals and perennials
cover[, sr_annual := length(unique(Taxon[local_lifespan=="ANNUAL"])), by = .(plot, site_code, year)]
cover[, sr_peren := length(unique(Taxon[local_lifespan=="PERENNIAL"])), by = .(plot, site_code, year)]
cover[, sr_null.lspan := length(unique(Taxon[local_lifespan=="NULL"])), by = .(plot, site_code, year)]  # this includes a lot non-live. 
cover[, sr_biennial := length(unique(Taxon[local_lifespan=="BIENNIAL"])), by = .(plot, site_code, year)]
cover[, sr_indeter := length(unique(Taxon[local_lifespan=="INDETERMINATE"])), by = .(plot, site_code, year)]

# plot 
hist(cover$sr_annual)
summary(cover$sr_annual)
summary(cover$sr_annual)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    1.00    3.13    4.00   23.00 
hist(cover$sr_peren)
summary(cover$sr_peren)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00    4.00    8.00   10.21   14.00   41.00 

# Compute first differences 
cover[order(year), change_sr_annual := sr_annual-shift(sr_annual), by =.(plot, site_code)]
cover[order(year), change_sr_peren := sr_peren-shift(sr_peren), by =.(plot, site_code)]

# Plot first differences
hist(cover$change_sr_annual)
summary(cover$change_sr_annual)
hist(cover$change_sr_peren)
summary(cover$change_sr_peren)

# compute cover by annuals per plot and year
cover[, AnnualPercentcover.yr := sum(relative_sp_cover.yr[local_lifespan=="ANNUAL"]), by = .(plot, site_code, year)]
# compute cover by perennial per plot and year
cover[, PerenPercentcover.yr := sum(relative_sp_cover.yr[local_lifespan=="PERENNIAL"]), by = .(plot, site_code, year)]

hist(cover$PerenPercentcover.yr)
hist(cover$AnnualPercentcover.yr)
#######################################################################################################
### Functional Group variables (e.g., Grass vs Forbs) ####################################################
#######################################################################################################

## Combine Graminoid and Grass into a single category 
# Graminoid + Grass = Graminoids
cover[functional_group == "GRASS" , functional_group:= "GRAMINOID"]

### Richness of grasses, woody, forbs  from the "functional_group" column
# Grasses & Graminoids
cover[, sr_graminoid := length(unique(Taxon[functional_group=="GRAMINOID"])), by = .(plot, site_code, year)]
# forbs
cover[, sr_forbs := length(unique(Taxon[functional_group=="FORB"])), by = .(plot, site_code, year)]
# woody
cover[, sr_woody := length(unique(Taxon[functional_group=="WOODY"])), by = .(plot, site_code, year)]
# legumes
cover[, sr_legume := length(unique(Taxon[functional_group=="LEGUME"])), by = .(plot, site_code, year)]
# Bryophyte
cover[, sr_bryophyte := length(unique(Taxon[functional_group=="BRYOPHYTE"])), by = .(plot, site_code, year)]
# Cactus
cover[, sr_cactus := length(unique(Taxon[functional_group=="CACTUS"])), by = .(plot, site_code, year)]
# Null
cover[, sr_null.fctgroup := length(unique(Taxon[functional_group=="NULL"])), by = .(plot, site_code, year)]
# Non-live
cover[, sr_non.live := length(unique(Taxon[functional_group=="NON-LIVE"])), by = .(plot, site_code, year)]

## Compute first differences ###
cover[order(year), change_sr_graminoid := sr_graminoid - shift(sr_graminoid), by =.(plot, site_code)]
cover[order(year), change_sr_woody := sr_woody - shift(sr_woody), by =.(plot, site_code)]
#summary(cover$change_sr_woody)
cover[order(year), change_sr_forbs := sr_forbs -shift(sr_forbs), by =.(plot, site_code)]
#summary(cover$change_sr_forbs)
cover[order(year), change_sr_legume := sr_legume -shift(sr_legume), by =.(plot, site_code)]
cover[order(year), change_sr_bryophyte := sr_bryophyte -shift(sr_bryophyte), by =.(plot, site_code)]
cover[order(year), change_sr_cactus := sr_cactus  -shift(sr_cactus), by =.(plot, site_code)]
cover[order(year), change_sr_null.fctgroup := sr_null.fctgroup - shift(sr_null.fctgroup), by =.(plot, site_code)]
cover[order(year), change_sr_non.live := sr_non.live - shift(sr_non.live), by =.(plot, site_code)]

### Compute cover by Functional groups - grasses, woody, forbs, legumes, etc -- per plot and year
cover[, GrassPercentcover.yr := sum(relative_sp_cover.yr[functional_group=="GRASS"]), by = .(plot, site_code, year)]
cover[, ForbPercentcover.yr := sum(relative_sp_cover.yr[functional_group=="FORB"]), by = .(plot, site_code, year)]
cover[, WoodyPercentcover.yr := sum(relative_sp_cover.yr[functional_group=="WOODY"]), by = .(plot, site_code, year)]
cover[, LegumePercentcover.yr := sum(relative_sp_cover.yr[functional_group=="LEGUME"]), by = .(plot, site_code, year)]

hist(cover$LegumePercentcover.yr)
plot(cover$LegumePercentcover.yr ~ cover$year)
lm(cover$LegumePercentcover.yr ~ as.numeric(cover$year))
# Coefficients:
#   (Intercept)  as.numeric(cover$year)  
# 2.107210               -0.001027 
abline(lm(cover$LegumePercentcover.yr ~ as.numeric(cover$year)))

# Lagged Legume cover per year
cover[order(year), lagged_LegumePercentCover.yr := shift(LegumePercentcover.yr), by =.(plot, site_code)]

##############################################################################################################################
### N-Fixing Species Variables ##############################################################################################
##############################################################################################################################

### Number of N-fixer in a plot over time
cover[, sr_Nfixer := length(unique(Taxon[N_fixer==1])), by = .(plot, site_code, year)]
cover[, sr_non.Nfixer := length(unique(Taxon[N_fixer==0])), by = .(plot, site_code, year)]
# what are the plots where there are 0 species that are NOT n-fixers?
  hist(cover$sr_non.Nfixer)
  table(cover$sr_non.Nfixer)
  
# change in N-fixers in a plot
cover[order(year), change_sr_Nfixer := sr_Nfixer-shift(sr_Nfixer), by =.(plot, site_code)]
hist(cover$change_sr_Nfixer)
summary(cover$change_sr_Nfixer)

# lagged N-fixer richness
cover[order(year), lagged_sr_Nfixer := shift(sr_Nfixer), by =.(plot, site_code)]

### Percent cover of N-fixers ###
cover[, N_fixer_cover.yr := sum(relative_sp_cover.yr[N_fixer==1]), by = .(plot, site_code, year)]
# change in N-fixer cover
cover[order(year), change_N_fixer_cover := N_fixer_cover.yr - shift(N_fixer_cover.yr), by = .(plot, site_code, year)]

# Lagged N_fixer cover
cover[order(year), lagged_N_fixer_cover.yr := shift(N_fixer_cover.yr), by =.(plot, site_code)]

plot(cover$lagged_N_fixer_cover.yr ~ cover$year)
abline(lm(cover$lagged_N_fixer_cover.yr ~ as.numeric(cover$year)))
lm(cover$lagged_N_fixer_cover.yr ~ as.numeric(cover$year))
#   Coefficients:
#     (Intercept)  as.numeric(cover$year)  
#   2.612394               -0.001279

summary(cover$change_N_fixer_cover)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0       0       0       0       0       0   14822 

#############################################################################################
### Compute site-level richness variables  ###################################################
#################################################################################################

# # average site richness over time
# cover[, ave_siterich := mean(site_richness, na.rm = TRUE), by = site_code]
# # average for the change year to year
# cover[, ave_siterich_change := mean(changerich, na.rm = TRUE), by = site_code]
# # initial site richness
# cover[, init_siterich := site_richness[year_trt == 0] ,  by = site_code]

####################################################################################################################
### Dominance Variables using DI ############################################################################################
####################################################################################################################

### Do relative abundance of live cover only. 
# Do SR of live only, throughout?
### do this for live cover only (discussed w Kim in June)

### Dominance indicator calculation ###
# Use the dominance indicator metric from Avoilio et al (seperates dominance indication from impact)
  # DI = (average relative abundance + relative frequency)/2
    #Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
    #Relative frequency = number of sampling units a species occurred / total number of sampling units
# note: ranges from 0-1; Relative abundance can be any measure of abundance. Does not incorporate a measure of impact.
# There is not a cutoff for "which range from 0-1 = dominant species versus subordinate or rare, in the Avolio et al paper # 
# so if we want to group species in these groups, we will need to make one (then also should test robustness to that decision)

#for LIVE cover
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

# relative abundance per species and plot in pre-treatment year (year_trt == 0)
cover[, rel_abundance_year0 := relative_sp_cover.yr.live[year_trt == 0], by = .(Taxon, plot, site_code)]
hist(cover$rel_abundance_year0, xlab = "Average relative abundance at a site", main ="")
summary(cover$rel_abundance_year0)

## Compute Relative Frequency per spp AT THE SITE (defining dominance in space, not time/across years)
# "Relative frequency = number of sampling units a species occurred / total number of sampling units" 
# this should be # of plots within a site that the species occured in / total # of plots within a site, for pre-treatment year 

#total # of plots within a site, for pre-treatment year 
# & to filter to do just on the live species:
cover[, tot.num.plots := length(unique(plot[year_trt == 0])), by =.(site_code)]

#number of plots within a site, in the pre-treatment year, a species occurred in:
# & to filter to do just on the live species:
cover[, tot.num.plots.with.spp := length(unique(plot[year_trt == 0 & live==1])), by =.(site_code, Taxon)]

# test to see if it works
  abisko.test = cover[site_code == "abisko.se" , ]

#Compute Relative Frequency
# " Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots, by = .(plot, site_code)]
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq.space, xlab = "Frequency at the site in pre-treatment year", main = "Frequency of occurrence")
summary(cover$rel_freq.space)

## Compute the DI per species defined as:  DI = (average relative abundance + relative frequency)/2

#FOR LIVE COVER -- the dead cover will be 0.
cover[, DI := (rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, plot, site_code)]
 # summary(cover$DI)

## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI], xlab = "Dominance indicator (DI)", main = "Dominance indicator metric")
# summary(cover[live == 1,DI])

####### ####### ####### ####### ####### ####### ####### ####### 
####### Make Categorical Variables of Dominance too ########
####### ####### ####### ####### ####### ####### ####### ####### 

# then given them a ranking into different categories? 
# to look at changes in those types of species and the impact on productivity?
## sub-ordinate and rare -- how can I compute this? ###########
cover[, DIgroup:=cut(DI, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group 
#*MAKE SURE ONLY TO DO FOR LIVE 
cover[, sr_domspp := length(unique(Taxon[DIgroup == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp := length(unique(Taxon[DIgroup == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp := length(unique(Taxon[DIgroup == "Subordinate"])), by = .(plot, site_code, year)]

### create a category that groups together all the non-rare
cover[, sr_non_rare_spp := length(unique(Taxon[DIgroup %in% c("Subordinate", "Dominant")])), by = .(plot, site_code, year)]
#create a rare grouping (including native and invasive)
cover[, rareSR_spp.DI1 := length(unique(Taxon[DIgroup %in% c("Rare")])), by = .(plot, site_code, year)]


summary(cover$sr_domspp) # max is 2 species that are dominant 
summary(cover$sr_subordspp) 
summary(cover$sr_rare) 
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
# > summary(cover$change_sr_subordspp) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0327   0.0000  10.0000     3005 

summary(cover$change_sr_rare) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0199   0.0000   7.0000     3005 

# lagged effects of sr in these groups
cover[order(year), lagged_sr_rarespp := shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), lagged_sr_domspp := shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), lagged_sr_subordspp := shift(sr_subordspp), by =.(plot, site_code)]

plot(ave_rel_abundance_over_time.live ~ DI, data= cover)
abline(lm(ave_rel_abundance_over_time.live ~ DI, data= cover), col = "yellow")
plot(relative_sp_cover.yr.live ~ DI, data = cover)

## Make a variable that captures changes in *cover* of species classified as dominance, rare, subordinate each year 
# (so we can look at variations from year to year)
cover[, Dom_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Dominant"]), by = .(plot, site_code, year)]
cover[, Subord_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Subordinate"]), by = .(plot, site_code, year)]
cover[, Rare_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Rare"]), by = .(plot, site_code, year)]

hist(cover$Dom_cover.yr)
hist(cover$Rare_cover.yr)
hist(cover$Subord_cover.yr)
plot(cover$Rare_cover.yr, cover$sr_rarespp)
plot(cover$Dom_cover.yr, cover$sr_rarespp)

# make a change in the relative_sp_cover.yr.live and then look at relationship with DI?
cover[order(year), change_rel_abundance_time.live := relative_sp_cover.yr.live - shift(relative_sp_cover.yr.live), by =.(plot, site_code)]

hist(cover$change_rel_abundance_time.live, na.omit = TRUE)
plot(change_rel_abundance_time.live ~ DI, na.omit = TRUE, data= cover)

#plot this change_rel_abundance_time.live by each DI grouping category
# need to label each spp as rare, dominate or subordinate to do this, and then put in facet r

#maybe should do this as average change per spp?
ggplot(data = cover, aes(x = change_rel_abundance_time.live)) + geom_histogram()+ facet_wrap(~DIgroup) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in relative abundance yr-to-yr") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 

ggplot(data = cover, aes(y = relative_sp_cover.yr.live, x = year)) + geom_point() + facet_wrap(~DIgroup) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in relative abundance yr-to-yr") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 

trends_by_group <- ggplot(data = cover, aes(y = change_rel_abundance_time.live[na.omit= TRUE], x = year)) + geom_point() + facet_wrap(~DIgroup) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in relative abundance yr-to-yr") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 

trends_by_group + geom_smooth(method = "glm", 
              formula = change_rel_abundance_time.live[na.omit = TRUE] ~ year,
              family = gaussian(link = 'log'))

# look at N-fixer cover over time
# change_N_fixer_cover -- not much at all.

# ggplot(data = cover, aes(x = change_N_fixer_cover)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
#   geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
#   labs(x = "Plot-level change in relative abundance yr-to-yr") +  theme_bw() +
#   theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
#   theme(axis.text.y = element_text(size = 14)) + 
#   #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text.x = element_text(size=14)) 
# 
# ggplot(data = cover, aes(x = N_fixer_cover.yr)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
#   +   geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
#   +   labs(x = "Plot-level change in relative abundance yr-to-yr") +  theme_bw() +
#   +   theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
#   +   theme(axis.text.y = element_text(size = 14)) + 
#   +   #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   +   theme(axis.text.x = element_text(size=14))

#########################################################################################################################
### Create variables that are combined groups of: #######################################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
##############################################################################################################################
# create combinations of all - as factor in a column
cover[, status.NN.RareDom := paste(DIgroup,local_provenance, sep = "_")]
# table(cover$status.NN.RareDom)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SR variables by combined groupings ####### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#do SR for non-native, rare:
cover[, sr_non.nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover[, sr_nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

## look a the data 
  hist(cover$sr_non.nat_rare)
  hist(cover$sr_nat_rare)

# create a non-rare variable
cover[, non_rare_spp := DIgroup %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

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

### cover variables by each group and the combined groupings ###
#using relative_sp_cover.yr.live based on year 0 groups

#Compute overall sum of relative cover by groups over time variables
cover[, cover_tot_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" ]), by = .(plot, site_code, year)]
cover[, cover_tot_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" ]), by = .(plot, site_code, year)]
cover[, cover_tot_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" ]), by = .(plot, site_code, year)]
cover[, cover_tot_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate"]), by = .(plot, site_code, year)]
cover[, cover_tot_NAT := sum(relative_sp_cover.yr.live[local_provenance == "NAT"]), by = .(plot, site_code, year)]
cover[, cover_tot_INT := sum(relative_sp_cover.yr.live[local_provenance == "INT"]), by = .(plot, site_code, year)]


#for NON rare (grouping dominant and subordinate)
cover[, cover_non.nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

#for rare
cover[, cover_nat_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]
cover[, cover_non.nat_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" & local_provenance == "INT"]), by = .(plot, site_code, year)]

# for dominant
cover[, cover_non.nat_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

#native dominant cover is declining through time on average
trend1 <- lm(cover_nat_dom  ~ site_code + as.numeric(year_trt) , data = cover)
summary(trend1)

# for subordinate
cover[, cover_non.nat_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate"& local_provenance == "NAT"]), by = .(plot, site_code, year)]

#cover[, cover_nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

# plot cover
par(mfrow = c(2,2), pty = "s")
hist(cover$cover_non.nat_non.rare, main = "cover non-native non-rare", xlab = "cover non-native non-rare")
hist(cover$cover_nat_non.rare, main = "cover native non-rare", xlab = "cover native non-rare")
hist(cover$cover_nat_rare, main = "cover native rare")
hist(cover$cover_non.nat_rare, main = "cover non-native rare")

#plot changes 
plot(cover$cover_non.nat_non.rare ~ cover$year)
abline(lm(cover$cover_non.nat_non.rare ~ length(unique(cover$year))), col = "red")

trend1 <- lm(cover_nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)
summary(trend1)
trend2 <- lm(cover_nat_dom  ~ site_code:as.numeric(year_trt) + as.numeric(year_trt) , data = cover)
summary(trend2)

# non-native rare trends - increaisng
trend1 <- lm(cover_non.nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)

#non native dom cover also decreasing
trend1 <- lm(cover_non.nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)

#EXAMPLE: explo[, pos_dev_thresh := dev_from_thresh[dev_from_thresh>1], by = .(site, year)]
#cover[, Fenced := trt[trt == "Fenced"], by = .(plot, site_code, year)]

#non native dom cover also decreasing
trend1 <- lm(cover_non.nat_rare  ~ site_code +  as.numeric(year_trt), data = cover)
summary(trend1)


############################################################################################################## #############
######################################################################################################### #############
## Dominance and Rarity Variables based on relative abundance and frequency groups seperately (versus DI above) #############
######################################################################################################### #############

######################################################################################################### #############
## Compute for Frequency: sub-ordinate and rare - Cut Off 1 based on breaks=c(0.0,0.2,0.8,1.0)
######################################################################################################### #############

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
#### DI Metric: Rare, Dom Subord Group Cut Off 2 #################################################################################################
############################################################################################################################
# Dec 29 2019 - fixed cutoff 2 on June 22 2020

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
#summary(cover$rel_abundance_year0)
hist(cover$rel_abundance_year0)

cover[, RelAbund_group2 :=cut(rel_abundance_year0, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, relabund_sr_domspp2 := length(unique(Taxon[RelAbund_group2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, relabund_sr_rarespp2 := length(unique(Taxon[RelAbund_group2 == "Rare"])), by = .(plot, site_code, year)]
cover[, relabund_sr_subordspp2 := length(unique(Taxon[RelAbund_group2 == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$relabund_sr_domspp2) # max 2 
summary(cover$relabund_sr_rarespp2) # max 43
summary(cover$relabund_sr_subordspp2) #max 3, min 0 ?? mean 2 

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


###########################################################################################################################
#### Compute variables for site-level richness to assess: ########################################################
# How do results depend on average site-level rich (averaged over time for the site)? #################################################################################################
############################################################################################################################
# April 24, 2019 - moved over from mechanisms analysis Aug 29 2019.

#compute the site richness from the cover data
cover[, new_siterich := uniqueN(Taxon), by = .(site_code, year)]
cover[, new_initsiterich := new_siterich[year_trt == 0][1], by = .(site_code)]
plot(new_siterich ~ new_initsiterich, data = cover)

cover[, avesiterich  := mean(new_siterich, na.rm = TRUE), by = site_code]
plot(avesiterich ~ new_initsiterich, data = cover)
cor(cover$avesiterich, cover$new_initsiterich, use='complete.obs')
# 0.9470813

#compute site level change 
# *to do* mech.data[, ave_siterich_change := mean(changerich, na.rm = TRUE), by = site_code]
lag.data=cover[,.(new_siterich=new_siterich[1]),by=.(site_code, year)][,lag_new_siterich:=shift(new_siterich, n=1, type=), by=site_code]
cover = merge(cover, lag.data[,.(site_code,year,lag_new_siterich)], by=c("site_code","year"), all.x=T)

#check to make sure this worked
cover[site_code=='bldr.us', .(mean(new_siterich), mean(lag_new_siterich)), by=year ]

cover[, change_siterich := new_siterich - lag_new_siterich, by=year ]
hist(cover$change_siterich)




########################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above 
#****** & likely need to check for others since there could be more in the full dataset!!!!
cover = cover[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]



# write as csv datafile to use for R
write.csv(cover, "NutNetCoverData_ProcessedAug2019.csv")





#####################################################################################################################
### Create Variable to determine, year-to-year, if species are gained or lost or zero change, and if so who. 
# instead do this in the comb file? ##################
#####################################################################################################################

# df1[, z := y > 0]
# or DT[, flag := 6500<=to & 6500>=from]
# table(DT$flag)
# FALSE  TRUE 
# 5567  4433 

comb[, rich_increase := changerich>0]
comb[, rich_decrease := changerich<0]

barplot(prop.table(table(comb$rich_decrease)), main = "richness decrease")
barplot(prop.table(table(comb$rich_increase)), main = "richness increase")

### DO FOR RARE SPECIES and NON-NATIVEs!!!!!




##### This computes dominance over time, which is not what we want: 

## compute average relative abundance (averaged across years.)
# Relative abundance = abundance of a species in a sampling unit / total abundance of all species in a sampling unit
# here we need to average across the years.

#I did the following based on plot-years as sampling unit, but we decided it should be defined spatially.

# for TOTAL cover
# cover[, ave_rel_abundance_over_time := ave(relative_sp_cover.yr), by = .(Taxon, plot, site_code)]
# hist(cover$ave_rel_abundance_over_time)
# summary(cover$ave_rel_abundance_over_time) 

#for LIVE cover
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

## Compute Relative Frequency per spp per plot.
# "Relative frequency = number of sampling units a species occurred / total number of sampling units" 
# first calculate the total number of years for a given plot in a site "total # of sampling units"
cover[, tot.yr := length(unique(year)), by =.(plot, site_code)]
# to determine which is the most recent year: 
  # cover[, max.yr := max(year), by =.(plot, site_code)]

# spot check
cover[site_code == "abisko.se", .(year, plot, tot.yr),]
cover[site_code == "barta.us", .(year, plot, tot.yr),]
cover[site_code == "smith.us", .(year, plot, tot.yr),]

#then calculate number of sampling units (years) a species occurred
# with a list of taxa per plot and year, count the number of records per (species, plot) pair. 
# so for ex., if the datatable is called dt  
# dt[,.N, by=.(plot, species)] 
#will count the records per plot x species and if they only record a species in a plot and year if it's there
# cover[, .N , by =.(plot, site_code, Taxon)]

# to store this count as part of the cover datatable:
cover[, Num_year_spp_occur := length(unique(year)), by = .(plot, site_code, Taxon)]
cover[, Num_recs_spp_plot_site := .N, by = .(plot, site_code, Taxon)]
#***send this to Ashley to find duplicates
  # cover[Num_year_spp_occur!=Num_recs_spp_plot_site, ] #812 records with duplicates? 
   #585 from site_name=="Antisana" & year==2013:
      # cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013,]
# table(cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013, .(site_code, plot, Taxon)])
# cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013 & Taxon=="WERNERIA PUMILA", ]
# looks like this is repeated bc multiple subplots
#how many sites have multiple subplots and taxa recorded by it???? anything that would be recorded at the subplot level?
# the max cover for cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013 & Taxon=="WERNERIA PUMILA", ] 
# is the same
cover[,Num_subplots:=length(unique(subplot)),.(site_code,plot, Taxon)]
cover[Num_subplots>1,length(unique(site_code))]
unique(cover[Num_subplots>1,.(site_code, year)])
#  site_code year
#    anti.ec 2013
#    elkh.us 2007
#    sgs.us  2007
#    ucsc.us 2007

#compute relative frequency
# " Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.time :=  Num_year_spp_occur/tot.yr, by = .(plot, site_code)]
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq)
summary(cover$rel_freq)

## Compute the DI per species <-- at what scale does this apply? the plot over all years? 
# DI = (average relative abundance + relative frequency)/2
#FOR LIVE COVER
cover[, DI.time := (ave_rel_abundance_over_time.live + rel_freq)/2 , by =.(Taxon, plot, site_code)]
hist(cover$DI.time)
summary(cover$DI.time)

# # then given them a ranking into different categories? 
# # to look at changes in those types of species and the impact on productivity?
# ## sub-ordinate and rare -- how can I compute this? ###########
# cover[, DIgroup:=cut(DI, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]
# 
# #richness in each group
# cover[, sr_domspp := length(unique(Taxon[DIgroup == "Dominant"])), by = .(plot, site_code, year)]
# cover[, sr_rarespp := length(unique(Taxon[DIgroup == "Rare"])), by = .(plot, site_code, year)]
# cover[, sr_subordspp := length(unique(Taxon[DIgroup == "Subordinate"])), by = .(plot, site_code, year)]
# 
# summary(cover$sr_domspp)
# summary(cover$sr_subordspp)
# summary(cover$sr_rarespp)

# compute change in richness in each group 
cover[order(year), change_sr_domspp := sr_domspp -shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), change_sr_rarespp := sr_rarespp -shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), change_sr_subordspp := sr_subordspp -shift(sr_subordspp), by =.(plot, site_code)]

#?? compute cover in each group, or is that circular???? compute *CHANGE* in cover from year-to-year 



## [ Do this for Tim too, plot by treatment and site over time? ###
# see how -- for rare species -- if their DI score is more driven by the relative abundance or the low frequency #
# those might be different types of rare spp so it seems like we should not obscure them into one ## ] 



#**** NEED SOMETHING TO DO WITH THE CHANGES FROM YEAR TO YEAR AND HOW TO CHARACTERIZE THEM BY % of DIFFERENT GROUPS?**** 
# WHATS A 1% richness change and how does that map onto the proportion of spp in these different categories? esp bc 
# the mechanisms are not mutually exclusive.



