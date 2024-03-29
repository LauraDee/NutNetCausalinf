################################################################################################################
# Make some figures to get at importance of panel data methods################################################ 
# and sources of variation in NutNet #############################################################################
#############################################################################################################
graphics.off()
rm(list=ls())
#load libraries
require(ggplot2)
library(data.table)
library(viridis)

##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-28-Apr-2017.csv'
setwd("~/Documents/GitHub/NutNetCausalInf/data/processed/")
comb <- fread("NutNetControlPlotData_v201804.csv",na.strings='NA')

comb$site <- comb$site_code
comb$plot = as.factor(comb$plot)

# Goals of plots
# Baseline levels differ by site
# Variation through site after removing important differences.
# changerich and changelive_mass are first differences in rich
# and live_mass at the plot level. Preferred models use plot FE and
# site x year effects, so take the changes and demean by plots within
# same site and year

comb[,dm.changerich:=changerich-mean(changerich, na.rm=T), by=.(site,year)]
comb[,dm.changelive_mass:=changelive_mass-mean(changelive_mass, na.rm=T), by=.(site,year)]

# We'll also do (double) demeaning of logged values to mimic 
# what's done in log-log model
comb[,`:=`(log.rich=log(rich), log.live_mass=log(live_mass))]
comb[order(year), change.log.rich := log(rich)-shift(log(rich)), by =.(plot, site_code)]
comb[order(year), change.log.live_mass := log(live_mass)-shift(log(live_mass)), by =.(plot, site_code)]
comb[,dm.change.log.rich:=change.log.rich-mean(change.log.rich, na.rm=T), by=.(site,year)]
comb[,dm.change.log.live_mass:=change.log.live_mass-mean(change.log.live_mass, na.rm=T), by=.(site,year)]

# Plot for raw data relationship
ggplot(comb[!is.na(rich) & !is.na(live_mass),], 
       aes(x=rich, 
           y=live_mass)) + 
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point()

# Raw changes in richness vs changes in live mass 
ggplot(comb[!is.na(changerich) & !is.na(changelive_mass),], 
       aes(x=changerich, 
           y=changelive_mass)) + 
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point()


### LOG BIVARIATE PLOTS #### *** TO ADD TO TUTORIAL ******
# Plots for log-log models - punchline

# Cross-sectional analog - ignoring fixed effects and any other covariates
ggplot(comb[!is.na(log.rich) & !is.na(log.live_mass),], 
       aes(x=log.rich, 
           y=log.live_mass)) + 
  geom_smooth(method="lm", se=T) + labs(y = "log(Live biomass)") + labs(x = "log(Richness)") + 
  theme_bw() +
  geom_point()

# Plot FE only
ggplot(comb[!is.na(change.log.rich) & !is.na(change.log.live_mass),], 
       aes(x=change.log.rich, 
           y=change.log.live_mass)) + labs(y = "log(Live biomass)") + labs(x = "log(Richness)") + 
  geom_smooth(method="lm", se=F) + 
  theme_bw() +
  geom_point()

# Panel analog (plot FE and site-year effects)
ggplot(comb[!is.na(dm.change.log.rich) & !is.na(dm.change.log.live_mass),], 
       aes(x=dm.change.log.rich, 
           y=dm.change.log.live_mass)) + labs(y = "log(Live biomass)") + labs(x = "log(Richness)") + 
  geom_smooth(method="lm", se=F) +
  theme_bw() +
  geom_point()



######
# Decompose variation from one plot --- *****THESE FIGURES FOR THE TUTORIAL! *****
######
comb[,singledm.log.live_mass:=log.live_mass-mean(log.live_mass, na.rm=T), by=.(site, plot)]
comb[,doubledm.log.live_mass:=singledm.log.live_mass-mean(singledm.log.live_mass, na.rm=T), by=.(site, year)]

# plot raw variation - no fixed effects
ggplot(comb[(site=="sedg.us" & plot %in% c("1","17")) | (site=="sevi.us" & plot %in% c("8","12")),],
       aes(x=year, y=log.live_mass, group=plot, linetype=plot, color = site)) + 
  geom_line() +   scale_color_manual(values=c('#999999','#E69F00')) +
  ggtitle("Raw Variation in log(live biomass)") + 
  theme_bw() + 
  ylim(c(-1,7)) + labs(y = "log(Live biomass)") +  labs(x = "Year") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 


## why does this look different? Plot FEs
ggplot(comb[site=="sedg.us" & plot %in% c("1","17") | (site=="sevi.us" & plot %in% c("8","12")), ],
       aes(x=year, y=singledm.log.live_mass, group=plot, linetype=plot, color = site)) + 
  geom_line() +    scale_color_manual(values=c('#999999','#E69F00')) +
  ggtitle("Variation in log(live biomass) after removing plot fixed effects") + 
  theme_bw() + 
  ylim(c(-5,7)) + labs(y = "log(Live biomass)") +  labs(x = "Year") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 

+  geom_segment(
    x = 1, y = 1,
    xend = 4, yend = 7,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 2, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "#EC7014") # Also accepts "red", "blue' etc
  


### Variation after accounting for site by year effects  ### 
sitebyyear = ggplot(comb[site=="sedg.us" & plot %in% c("1","17") | (site=="sevi.us" & plot %in% c("8","12")),],
       aes(x=year, y=doubledm.log.live_mass, group=plot, linetype=plot, color = site)) + 
  geom_line() +   scale_color_manual(values=c('#999999','#E69F00')) +
  ggtitle("Variation in log(live biomass) after removing plot fixed and site by year effects") +
  theme_bw() +  
  ylim(c(-5,7)) 

sitebyyear +
  labs(y = "log(Live biomass)") +  labs(x = "Year") + 
   theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 



#Look at all plots in sevi to see what happened in 2009
ggplot(comb[site=="sevi.us",],
       aes(x=year, y=log.live_mass, group=plot, linetype=plot)) + 
  geom_line() + 
  ggtitle("Raw variation in log biomass -- all plots in sevi.us") + 
  theme_bw() + 
  ylim(c(-1,7)) 




# Plot histogram of average richness by site
variation.histogram.plot = ggplot(comb[,mean(rich, na.rm=T), by=plot],
       aes(x=V1-11, y=..ncount..)) + 
  geom_histogram(binwidth=1, fill='gray') + 
  geom_histogram(data=comb, fill="black",
                 aes(x=changerich, y= -..ncount..),
                 binwidth=1) + 
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=18)) + 
  geom_text(aes(label="Cross-sectional differences",x=10, y=0.75), size=6) + 
  geom_text(aes(label="Temporal differences",x=10, y=-0.75), size=6) + 
  xlab("") + 
  ylab("Normalized count")

pdf("cross_vs_time_variation_histograms.pdf", width=8, height=8)
variation.histogram.plot
dev.off()



#plot trends in residuals - removing variables like in Grace et al - site_live_mass, shading (ground_PAR), soil fertility
# ** THIS IS INCOMPLETE AND LIKELY NOT NEEDED FOR TUTORIAL

plot(comb$live_mass, comb$site_live_mass)
comb[, site_live_mass := sum(live_mass, na.rm = T), by = site]
plot_prod.Grace = lm(log(live_mass) ~ log(site_live_mass) + log(rich)  , data = comb)
summary(plot_prod.Grace)

plot_prod.Grace.soilfert = lm(log(live_mass) ~ log(site_live_mass) + log(rich) + log(ppm_P) + PercentSilt + soilts_pH , data = comb)
summary(plot_prod.Grace.soilfert)
#soil fertility: ppm_P, PercentSilt, soilts_ppm_P, soilts_pH
plot_prod.Grace.soilfert = lm(log(live_mass) ~ log(site_live_mass) + log(rich) + log(soilts_ppm_P) + PercentSilt + soilts_pH , data = comb)
summary(plot_prod.Grace.soilfert)

# make a soil fert variable a la Grace :
# # make the Grace composite vaeiable - SoilFertility    = +0.144*LogP - 0.123*PH +0.458*Silt    
plot_prod.Grace.soilfertcomposite = lm(log(live_mass) ~ log(site_live_mass) + I(0.144*log(ppm_P) - 0.123*soilts_pH + 0.458*PercentSilt) + log(rich) , data = comb)
summary(plot_prod.Grace.soilfertcomposite)
screenreg(plot_prod.Grace.soilfertcomposite)

#rich mod:   Ground_PAR, site_rich, soil suitability
plot_rich.Grace = lm(log(rich) ~ log(site_year_rich) +  Ground_PAR   , data = comb)
summary(plot_rich.Grace)



