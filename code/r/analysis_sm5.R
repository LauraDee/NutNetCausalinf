#######################################################################################################################################
### Supplemental Analyses: Test robustness to moderators: SM sections S5c and S5d ################################################################################################
#######################################################################################################################################

##### Analyses for Section S5c: Moderating Effect of site-level species richness ####
## See Table S4: code to reproduce Table S4 is below.

ModrMod_AveSiteRich     <- feols(log(live_mass) ~ log(rich) + log(rich):site_richness | newplotid + site.by.yeardummy, comb) 
ModrMod_AveSiteIntRich  <- feols(log(live_mass) ~ log(rich) + log(rich):site_introduced_richness | newplotid + site.by.yeardummy, comb) 
ModrMod_YrSiteRich      <- feols(log(live_mass) ~ log(rich) + log(rich):site_year_rich | newplotid + site.by.yeardummy, comb) 
ModrMod_AveSiteNatRich  <- feols(log(live_mass) ~ log(rich) + log(rich):site_native_richness | newplotid + site.by.yeardummy, comb) 

esttex(ModrMod_AveSiteRich, ModrMod_AveSiteIntRich, ModrMod_YrSiteRich, ModrMod_AveSiteNatRich,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S4_R_se.tex")

esttex(ModrMod_AveSiteRich, ModrMod_AveSiteIntRich, ModrMod_YrSiteRich, ModrMod_AveSiteNatRich,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S4_R_ci.tex")

#######################################################################################################################################
### Supplemental Analyses:  S5d Moderating Effect of site-level productivity  ################################################################################################
#######################################################################################################################################

## See Table S5 & 6: code to reproduce results

###########################################################
### Compute average site-level productivity (live mass) ###
###########################################################
comb[, ave_site_livemass := ave(live_mass, na.rm = T), by = .(site_code)]
#hist(comb$ave_site_livemass)

comb[, ave_site_livemass.peryr := ave(live_mass, na.rm = T), by = .(site_code, year)]
#hist(comb$ave_site_livemass.peryr)

ModrMod_AveSiteProd <- feols(log(live_mass) ~ log(rich) + log(rich):ave_site_livemass  | newplotid + site.by.yeardummy, comb)
ModrMod_YrSiteProd  <- feols(log(live_mass) ~ log(rich) + log(rich):ave_site_livemass.peryr  | newplotid + site.by.yeardummy, comb)

esttex(ModrMod_AveSiteProd, ModrMod_YrSiteProd,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S5_R_se.tex")

esttex(ModrMod_AveSiteProd, ModrMod_YrSiteProd,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S5_R_ci.tex")



#### Analyses for Table S6 using cut-offs akin to Wang et al. 2019 Nature Communications

## Using the productivity groups from Wang et al 2019 Nature Communications
# Yongfans response:
# The 151 grids in HerbDivNet data were divided into three equal groups,
# depending on their mean productivity: low, medium, and high productivity (with 50 to 51 grids each).
# The cutoff values of primary productivity of each group:
#   
# Low (51 grids): 30.18-238.73 (g/m^2)
# Medium (50 grids): 239.67-409.69 (g/m^2)
# High (50 grids): 414.29-1382.42 (g/m^2)


# do the groups based on the overall average since Wang et al was a cross-section so that is more comparable
summary(comb$ave_site_livemass.peryr)
summary(comb$ave_site_livemass)
#check that it worked length(unique(comb$ave_site_livemass))  == 43 

## Do a different cut off Low, Medium, High groups of site livemass, based on Wang et al cut-offs
comb[, ProdGroup_WangCutoffs:=cut(ave_site_livemass, breaks=c(30.18,239.67,414.20,1609), labels=c("Low","Medium","High")), by = .(site_code) ]
head(comb)
table(comb$ProdGroup_WangCutoffs)
summary(comb$ProdGroup_WangCutoffs)
plot(comb$ProdGroup_WangCutoffs, main = "Productivity Groups (Average Across Years)")

## the Wang et al. paper uses a single year of data. Thus, we also need to create these cut-offs for each site and year in our data:
comb[, ProdGroup := cut(ave_site_livemass.peryr, breaks=c(30.18,239.67,414.20,1609), labels=c("Low","Medium","High"))]
head(comb)
table(comb$ProdGroup)
plot(comb$ProdGroup, main = "Productivity Groups Per Year")

## Table S6

ModrMod_ProdGroup <- feols(log(live_mass) ~ log(rich) + log(rich):ProdGroup | newplotid + site.by.yeardummy, comb)
ModrMod_ProdGroupWang  <- feols(log(live_mass) ~ log(rich) + log(rich):ProdGroup_WangCutoffs  | newplotid + site.by.yeardummy, comb)

esttex(ModrMod_ProdGroup, ModrMod_ProdGroupWang,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S6_R_se.tex")

esttex(ModrMod_ProdGroup, ModrMod_ProdGroupWang,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S6_R_ci.tex")