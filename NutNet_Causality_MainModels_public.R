########################################################################
## NutNet Main Text Models ###########################################
## Laura Dee - in LFE ###############################################
##########################################################################
# updated Dec 27, 2019
#updated May 11, 2020 to include Simpsons div models 
#cleaned up to focus on main results for Figures 2 & 3 on Feb 5 2021 

#notes on using lfe versus previous implementation of models:
#it does, especially if you're in a "FE nested within clusters" setting
#which the code I'd given you before did not account for
# since i was unaware of it at the time
# basically when doing degrees of freedom corrections, lfe won't count fixed effects nested
# in clusters b/c that would sort of be double counting

#plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

#plotting resurces: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# theme(legend.position="top") 
#http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software

#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(AER)
library(sandwich)
library(foreign)
library(car)
library(lfe)
library(texreg)
library(broom)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

# When log(0) need to use inverse hyperbolic sine transformation (REF NEEDED)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
}
# v.s. log(x+1) <- is defined for a negative x. 

##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-09-Apr-2018.csv'
setwd("~/Dropbox/IV in ecology/NutNet")
comb <- fread("NutNetControlPlotDataToUseApril2018.csv",na.strings='NA')

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
comP = comb[site_code == "comp.pt ",]

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
write.csv(comb.descript.v1, "DatasetDescript-ControlPlots-July2020.csv")
# length(unique(comb$site_code))

############################################################################################################################
#### Panel FE Plus: Plot-FE and site*year effects  ###########################################################################
##### # RUN AND PLOT:  ## Preferred Models - Main Text ###########################################################################
############################################################################################################################

######################################################################################################################
### Main Models (Models 1). Log-Log and fixed effects/dummies only. ################################################################
##################################################################################################################

##### Models with productivity as the Y variable as log live mass ##########
#A.  Log-log and fixed effects/dummies only.
ModPFE <- felm(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy | 0 | newplotid, data = comb)
summary(ModPFE, robust = TRUE)

## TESTS OF THE TIDY FUNCTON HERE****
cof <- tidy(ModPFE, robust = TRUE,  conf.int = T)
cof
# coefs.test <-tidy(ModPFE, cluster = T, robust = T)
cof <- tidy(ModPFE, robust = T)
cof

#with evenness:
ModPFE.2 <- felm(log(live_mass) ~ log(rich) + ihs(even)  | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFE.2 , robust = TRUE, cluster = TRUE)

#with lagged richness:
ModPFE.3 <- felm(log(live_mass) ~ log(rich) + log(laggedrich) | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.3 , robust = TRUE, cluster = TRUE)

#with lagged rich and evenness
ModPFE.4 <- felm(log(live_mass) ~ log(rich) + log(laggedrich) + ihs(even) | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.4 , robust = TRUE, cluster = TRUE)

#print log-log results
screenreg(list(ModPFE, ModPFE.2, ModPFE.3, ModPFE.4),     # object with results 
          custom.model.names= c("Main Model 1" , "Incld Evenness", "Incld Lagged Richness", "Incld Both"))
      

#print log-log results with clustered robust SEs and corresponding p-values 
screenreg(list(ModPFE, ModPFE.2, ModPFE.3, ModPFE.4),     # object with results 
          custom.model.names= c("Main Model 1" , "Incld Evenness", "Incld Lagged Richness", "Incld Both"),
          override.se=list(summary(ModPFE,)$coef[,2],
                           summary(ModPFE.2,)$coef[,2],
                           summary(ModPFE.3,)$coef[,2],
                           summary(ModPFE.4,)$coef[,2]),
          override.pval=list(summary(ModPFE,)$coef[,4],
                                 summary(ModPFE.2,)$coef[,4],
                                 summary(ModPFE.3,)$coef[,4],
                                 summary(ModPFE.4,)$coef[,4]), 
          )

## The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

coefs_ModPFE <- tidy(ModPFE, conf.int = T, robust = T)
coefs_ModPFE.2 <- tidy(ModPFE.2, conf.int = T, robust = T)
coefs_ModPFE.3 <- tidy(ModPFE.3, conf.int = T,  robust = T)
coefs_ModPFE.4 <- tidy(ModPFE.4, conf.int = T,  robust = T)

panelFE.main <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Model 1"),
  coefs_ModPFE.2  %>% mutate(reg = "Model 1 with evenness"),
  coefs_ModPFE.3  %>% mutate(reg = "Model 1 with lagged richness"),
  coefs_ModPFE.4  %>% mutate(reg = "Model 1 with evenness & lagged richness"),
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange() + theme_classic() +
  scale_colour_discrete(name="Model") +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) +
  facet_wrap(~reg)
  theme(axis.title.x = element_blank())

panelFE.main + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

panelFE.main + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

# try to put all models on one line but group them
panelFE.main.2 <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Model 1"),
  coefs_ModPFE.2  %>% mutate(reg = "Model 1 with evenness"),
  coefs_ModPFE.3  %>% mutate(reg = "Model 1 with lagged richness"),
  coefs_ModPFE.4  %>% mutate(reg = "Model 1 with evenness & lagged richness"),
  coefs_Mod.R1 %>% mutate(reg = "Model 1 with total live cover")
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
 theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) 
 # + facet_wrap(~reg)
# + theme(axis.title.x = element_blank())

panelFE.main.2 + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate")
 # labs(fill = "reg")

############################################################################################
# Simpson's Diversity Models - Panel FE #######################################################
##########################################################################################
##### Models with productivity as the Y variable as log live mass ##########
#D.  Log-log and fixed effects/dummies only.
ModPFEsimpsonD <- felm(log(live_mass) ~ log(simpson)  | newplotid + site.by.yeardummy |0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFEsimpsonD, robust = TRUE, cluster = TRUE)

#print log-log results
#print results
screenreg(list(ModPFE, ModPFE.2, ModPFEsimpsonD),     # object with results 
          custom.model.names= c("Main Model with Richness" , "Main Model with Richness and Evenness", "Model with Simpson's D"))

#print results - simple table for main results 
screenreg(list(ModPFE, ModPFE.2, ModPFE.3, ModPFE.4, ModPFEsimpsonD),     # object with results 
          custom.model.names= c("Main Model 1" , "Incld Evenness", "Incld Lagged Richness", "Incld Both", "Simpson's D"))

#######################################################################################################################################
### Make Figure 2 for Main text  #################################################################################################################
#######################################################################################################################################
coefs_ModPFE <- tidy(ModPFE, conf.int = T, robust = T)
coefs_ModPFE.2 <- tidy(ModPFE.2, conf.int = T,  robust = T)
coefs_ModPFEsimpsonD <- tidy(ModPFEsimpsonD, conf.int = T,  robust = T)

# coefficients from the traditional ecological model design (Ferraro estimated in STATA)
coefs.ferraro = tibble(term="log(rich)",
                       estimate=0.3792555,
                       std.error=0.1909606,
                       statistic=1.99,
                       p.value=0.047,
                       conf.low=0.0049797,
                       conf.high=0.7535313)

plot.data.main <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Model with Species Richness"),
  coefs_ModPFE.2  %>% mutate(reg = "Model  with Richness and Evenness"),
  coefs_ModPFEsimpsonD  %>% mutate(reg = "Model with Simpson's Diversity"),
  coefs.ferraro %>% mutate(reg = "Traditional Ecological Approach")
 # coefs_ModPFE.4  %>% mutate(reg = "Model 1 with evenness & lagged richness"),
) 

panelFE.main <-
  ggplot(plot.data.main,
      #   aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
      aes(x=term, y=estimate, ymin= estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange() + theme_classic() +
  scale_colour_discrete(name="Model") +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) +
  facet_wrap(~reg) +
theme(axis.title.x = element_blank())

panelFE.main + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

panelFE.main + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

##############################################################################
### Figure 2 - Main Text #########################################################
#################################################################################
coefs_ModPFE <- tidy(ModPFE, conf.int = T, robust = T)
coefs_ModPFE.2 <- tidy(ModPFE.2, conf.int = T, robust = T)
coefs_ModPFEsimpsonD <- tidy(ModPFEsimpsonD, conf.int = T, robust = T)

coefs.ferraro = tibble(term="log(rich)",
                       estimate=0.3792555,
                       std.error=0.1909606,
                       statistic=1.99,
                       p.value=0.047,
                       conf.low=0.0049797,
                       conf.high=0.7535313)

# try to put all models on one line but group them
panelFE.main.2.data <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Richness Model"),
  coefs_ModPFE.2  %>% mutate(reg = "Richness Model controlling for Evenness"),
  coefs_ModPFEsimpsonD  %>% mutate(reg = "Simpson's Diversity Model"),
  coefs.ferraro %>% mutate(reg = "Traditional Ecological Approach")
) 
panelFE.main.2.data$term = factor(panelFE.main.2.data$term,
                                  levels=c("log(rich)", 
                                           "log(simpson)",
                                           "ihs(even)"))
panelFE.main.2 <-  ggplot(panelFE.main.2.data,
                          aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), size = 1, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
 # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) 
# + facet_wrap(~reg)
# + theme(axis.title.x = element_blank())

panelFE.main.2 + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate")
# labs(fill = "reg")

################################################################################################################################
###**** Final Figure****  Figure 2 - Main Text - without plotting evenness estimate #########################################################
############################################################################################################################################

#* to do: change the order of the models:
# Change the order of items - with limits fct p + scale_x_discrete(name ="Dose (mg)", limits=c("2","1","0.5"))
coefs_ModPFE <- tidy(ModPFE, conf.int = T, robust = T)
coefs_ModPFE.2 <- tidy(ModPFE.2, conf.int = T,  robust = T)
coefs_ModPFEsimpsonD <- tidy(ModPFEsimpsonD, conf.int = T,  robust = T)

#to pull out just the richness term for plotting with other models -- to avoid plotting evenness 
coefs_ModPFE.2  <- filter(coefs_ModPFE.2 , term == "log(rich)") 

#traditional ecological modeling approach estimate from Ferraro - in STATA
coefs.ferraro = tibble(term="log(rich)",
                       estimate=0.3792555,
                       std.error=0.1909606,
                       statistic=1.99,
                       p.value=0.047,
                       conf.low=0.0049797,
                       conf.high=0.7535313)

# try to put all models on one line but group them
panelFE.main.2.data <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Richness Model"),
 coefs_ModPFE.2  %>% mutate(reg = "Richness Model controlling for Evenness"),
 coefs_ModPFEsimpsonD  %>% mutate(reg = "Simpson's Diversity Model"),
 coefs.ferraro %>% mutate(reg = "Traditional Ecological Approach")
) 
panelFE.main.2.data$term = factor(panelFE.main.2.data$term,
                                  levels=c("log(rich)", 
                                           "log(simpson)",
                                           "ihs(even)"))

panelFE.main.2 <-  ggplot(panelFE.main.2.data,
          #                aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  aes(x=term, y=estimate, ymin= estimate - (1.96*std.error), ymax=estimate + (1.96*std.error) , colour = term)) +
  geom_pointrange(aes(col = reg), size = 1, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8) +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = "", 
  ) 
# + facet_wrap(~reg)
# + theme(axis.title.x = element_blank())

#Alternative y-axis label:
panelFE.main.2 + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")

p <- panelFE.main.2 + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="Model") + scale_color_manual(values=cbPalette[c(7,9,4,8)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
 theme(legend.background = element_rect(# fill="lightblue", 
                                        size=0.5, linetype="solid",
                                             colour ="black"))
pp
# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
              axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))

# print final figure:
ppp + scale_x_discrete(labels = c('log(Species Richness)','log(Simpsons Diversity)')) + theme(
  axis.title.x = element_text(size = 20),
   axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16))

###### End ######################

###########################################################################################
### Figure 2 A FINAL #################################################################
#################################################################

#### Plot Result from our Design separately:
coefs_ModPFE <- tidy(ModPFE, conf.int = T , robust = T)
coefs_ModPFE.2 <- tidy(ModPFE.2, conf.int = T, robust = T)
coefs_ModPFEsimpsonD <- tidy(ModPFEsimpsonD, conf.int = T , robust = T)

#to pull out just the richness term for plotting with other models -- to avoid plotting evenness 
coefs_ModPFE.2  <- filter(coefs_ModPFE.2 , term == "log(rich)") 

# try to put all models on one line but group them
panelFE.main.2.data <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Species Richness"),
  coefs_ModPFE.2  %>% mutate(reg = "Species Richness controlling for Evenness"),
  coefs_ModPFEsimpsonD  %>% mutate(reg = "Simpson's Diversity"),
 # coefs.ferraro %>% mutate(reg = "Traditional Ecological Approach")
) 
panelFE.main.2.data$term = factor(panelFE.main.2.data$term,
                                  levels=c("log(rich)", 
                                           "log(simpson)",
                                           "ihs(even)"))

panelFE.main.2 <-  ggplot(panelFE.main.2.data,
                       #   aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  aes(x=term, y=estimate, ymin= estimate - (1.96*std.error), ymax=estimate + (1.96*std.error) , colour = term)) +
  geom_pointrange(aes(col = reg), size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8) +
  scale_y_continuous(name =  "Estimated effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", 
) 
# + facet_wrap(~reg)
# + theme(axis.title.x = element_blank())
panelFE.main.2

#Alternative y-axis label:
panelFE.main.2 + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")
#panelFE.main.2  + scale_color_manual(values=cbPalette[c(7,3,4,8)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

p <- panelFE.main.2 + theme(legend.position = c(0.6, 0.77)) + scale_colour_discrete(name="Model") + scale_color_manual(values=cbPalette[c(7,9,4,8)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
pp
# # adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
                   axis.title=element_text(size=20,face="bold")) +
                 theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))

# print final figure:
ppp <- ppp + scale_x_discrete(labels = c('log(Species Richness)','log(Simpsons Diversity)')) + theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16))

# a title to point out these results are from "Our study Design"

ppp <- ppp + scale_x_discrete(labels = c("Species Richness", "Simpson's Diversity")) 
Fig2A <- ppp + labs(title="Our Main Study Design") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 
#print final
Fig2A 


##### Fig 2 B ############################################################################################
## Plot traditional ecol. model results on their own plot: #####################################
################################################################################################

#traditional ecological modeling approach estimate from Ferraro - in STATA
coefs.ferraro = tibble(term="log(rich)",
                       estimate=0.3792555,
                       std.error=0.1909606,
                       statistic=1.99,
                       p.value=0.047,
                       conf.low=0.0049797,
                       conf.high=0.7535313)

##  Include other STATA output:
# Is the first estimate here for Simpson's the one from xtreg (the one to report,I gather, based on your comment on the evenness estimate?
# coef estimate; std error ; t val ; p- val, 95% conf interval [low, high]:
#           l_simpson |   .0744933   .0891051     0.84   0.403    -.1001496    .2491361
#           l_simpson |   .0746331   .0804479     0.93   0.354    -.0830419    .2323081

coefs.ferraro.simpsons = tibble(term="log(simpsons)",
                                estimate=.0744933,
                                std.error=.0891051 ,
                                statistic= 0.84,
                                p.value=  0.403,
                                conf.low= -.1001496 ,
                                conf.high= .2491361)

#STATA output from Paul 
# richness model with EVENNESS: Then I replicated the mixed-model with the addition of evenness (using xtreg and xtmixed – until I figure out whether xtmixed is estimating SEs correctly, use the first estimate from xtreg with random effects estimator).
#  # coef estimate; std error ; t val ; p- val, 95% conf interval [low, high]:
#              l_rich |   .3907834   .1838099     2.13   0.034     .0305227    .7510441
#              l_rich |   .3870898   .0956131     4.05   0.000     .1996916     .574488
# 
coefs.ferraro.even  = tibble(term="log(rich)",
                             estimate= .3907834,
                             std.error= .1838099,
                             statistic= 2.13,
                             p.value= 0.034,
                             conf.low= .0305227,
                             conf.high= .7510441)

# try to put all models on one line but group them
Fig2B.data <-  bind_rows(
   coefs.ferraro %>% mutate(reg = "Species Richness"), # "Richness Model"
   coefs.ferraro.even  %>% mutate(reg = "Species Richness controlling for Evenness"),  # "Richness Model controlling for Evenness"
   coefs.ferraro.simpsons  %>% mutate(reg = "Simpson's Diversity") # "Simpson's Diversity Model"
) 

Fig2B.data$term = factor(Fig2B.data$term,
                                  levels=c("log(rich)", 
                                           "log(simpson)",
                                           "ihs(even)"))
Fig.2B <-  ggplot(Fig2B.data,
                          aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), size = 1.5, position = position_dodge(width = 0.5)) +
  # scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
 # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8) +
  scale_y_continuous(name =  "Estimated effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig.2B 

#Alternative y-axis label:
Fig.2B + labs(
  title = "Traditional Ecological Design",
  # caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")
  caption = "", x = "Variable", y = "Estimated effect size") 
Fig.2B

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1", "red")

Fig.2B <- Fig.2B + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="Model") + scale_color_manual(values=cbPalette[c(8, 3, 10)])   +  labs(
  title = "Traditional Ecological Design",
 # caption = "", x = "Variable", y = "Estimate for log(species richness) effect size") + 
 caption = "", x = "Variable", y = "Estimated effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=14)) 
Fig.2B

# adjusting the legend 
Fig.2B <- Fig.2B +  theme(legend.title = element_blank()) +
  theme(legend.position = c(.5, 0.3)) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
Fig.2B

# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
Fig.2B <- Fig.2B + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
Fig.2B

# simplified variable names
Fig.2B <- Fig.2B + scale_x_discrete(labels = c('Species Richness', "Simpson's Diversity")) + theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
Fig.2B

###############################################################################################################################
### Put Fig 2 plots in same plot #####################################################################################################
#####################################################################################################
plot_grid(Fig2A , Fig.2B )

common.ylab = ylab("Estimated effect size")  #Estimated % Change in Productivity from a 1% Change in Diversity
plot_grid(Fig2A  + common.ylab,
          Fig.2B + common.ylab)

###################################################################################################
###### For information on the following, see the file NutNutAnalyses_SMSection5.R              ####
###### SM section S5b: FUNCTIONAL FORM ASSUMPTION CHECK Table S3                                ####
###### SM section 5c Analyses for Table S4: Moderating Effect of site-level species richness    ####
###### SM  SM section 5d Analyses for Table S5 &v S6: Moderating Effect of site-level productivity #
####################################################################################################

####################################################################################
##** Models for Figure 3 in the Main Text ##########################################
####################################################################################

#############################################################################################################################################
## Robustness Models: LDV models and IV models, R2 for Oster Test #####################################################################################
#################################################################################################################################

###############################################################################################
#### Print Info for Oster analysis ###########################################################
################################################################################################
# what is R^2 of richness equation?
#Log-log and fixed effects/dummies only - for Oster analysis.
#Richness Selection Model
richFE.forOster = lm(log(rich) ~ newplotid + site_code:year, data = comb) 
# Rsquared 0.9081 
summary(richFE.forOster)

#######################################################################################################################################
### Lagged-dependent models ############################################################################################################
#########################################################################################################################################

## there are NAs for lagged live mass so we need to subset this data
summary(comb$laggedlive_mass)
comb.laggedmod.dat = comb[!is.na(laggedlive_mass)]

# write out what this dataset includes
comb.lagged.descript =  table(comb.laggedmod.dat$site_name, comb.laggedmod.dat$year)
#write.csv(comb.descript.v1, "DatasetDescript-ControlPlots_laggedanal.csv") 
# Determine the number of Obs.
nrow(comb.laggedmod.dat) #confirm it's 1075 observations  # ??? now says 1075? 

#A.  Log-log and fixed effects/dummies only.
ModLD <- felm(log(live_mass) ~ log(rich)  + log(laggedlive_mass) | site.by.yeardummy | 0 | newplotid, data = comb.laggedmod.dat, exactDOF='rM')
summary(ModLD, robust = TRUE, cluster = TRUE)

#with evenness:
ModLD.2 <- felm(log(live_mass) ~ log(rich)  + log(laggedlive_mass) + ihs(even) | site.by.yeardummy | 0 | newplotid, data = comb.laggedmod.dat, exactDOF='rM')
summary(ModLD.2, robust = TRUE, cluster = TRUE)

#with laggged rich
ModLD.3 <- felm(log(live_mass) ~ log(rich)  + log(laggedlive_mass) + log(laggedrich) | site.by.yeardummy | 0 | newplotid, data = comb.laggedmod.dat, exactDOF='rM')
summary(ModLD.3, robust = TRUE, cluster = TRUE)

#with laggged rich and even
ModLD.4 <- felm(log(live_mass) ~ log(rich)  + log(laggedlive_mass) + log(laggedrich) + ihs(even) | site.by.yeardummy | 0 | newplotid, data = comb.laggedmod.dat, exactDOF='rM')
summary(ModLD.4, robust = TRUE, cluster = TRUE)

# print lagged dependent model results into a single table
screenreg(list(ModLD, ModLD.2, ModLD.3, ModLD.4),    
          custom.model.names=c("LDV", "LDV with Evenness", "LDV with lagged richness","LDV with evenness &lagged richnes"),
          omit.coef=c("(site_code)|(newplotid)"))  # object from estimation (unclustered) for BIC

### Plot Results #plotting coefficient estimates from felm objects:
coefs_ModLD <- tidy(ModLD, conf.int = T, robust = T)
coefs_ModLD.2 <- tidy(ModLD.2, conf.int = T, robust = T)
coefs_ModLD.3 <- tidy(ModLD.3, conf.int = T, robust = T)
coefs_ModLD.4 <- tidy(ModLD.4, conf.int = T, robust = T)

# try to put all models on one line but group them
LDM.main <-  bind_rows(
  coefs_ModLD  %>% mutate(reg = "Model LD"),
  coefs_ModLD.2  %>% mutate(reg = "Model LD with evenness"),
  coefs_ModLD.3  %>% mutate(reg = "Model LD with lagged richness"),
  coefs_ModLD.4  %>% mutate(reg = "Model LD with evenness & lagged richness"),
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
 # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.5, .5) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) 
# + facet_wrap(~reg)
# + theme(axis.title.x = element_blank())

LDM.main + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate")
# labs(fill = "reg")


# to plot all of the coefficients within one model
ggplot(data = coefs_ModLD.2, 
       mapping = aes(x = term , y = estimate, ymin = conf.low, ymax = conf.high, colour = term)) +
  geom_pointrange(fatten = 4) + theme_classic() + geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  labs(
    xlab = "Model X", 
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = "" ) +
  theme(axis.title.x = element_blank())


#to pull out just the richness term for plotting with other models
coefs_ModLD <- filter(coefs_ModLD, term == "log(rich)") 
coefs_ModLD.2 <- filter(coefs_ModLD.2, term == "log(rich)") 

bind_rows(
  coefs_ModLD %>% mutate(reg = "LD Model 1"),
  coefs_ModLD.2  %>% mutate(reg = "LD Model (with evenness)"),
) %>%
  ggplot(aes(x=reg, y=estimate, ymin=conf.low, ymax=conf.high, colour = reg)) +
  geom_pointrange() + theme_classic() +
  scale_colour_discrete(name="Model") +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-0.5, .3) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) +
  theme(axis.title.x = element_blank()) +  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
 # labs(title = "Log-Log Model Results") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(face="bold", size = 14))

#######################################################################################################################################
### IV models ############################################################################################################
#########################################################################################################################################
## IV here is average neighbor richness of the site (for each control plot and year) ##
#no NA's # summary(comb$avg_neighbor_rich)
modiv1 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
                 (log(rich) ~  log(avg_neighbor_rich)) | newplotid , # | newplotid + site.by.yeardummy), # first stage   Note the surrounding parentheses
               data = comb)
summary(modiv1, cluster = TRUE, robust = TRUE)

# check weak instrument test
summary(modiv1$stage1, cluster = TRUE, robust = TRUE)

# modiv1 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
#                  (log(rich) ~  ihs(avg_neighbor_rich)) | newplotid , # | newplotid + site.by.yeardummy), # first stage   Note the surrounding parentheses
#                data = comb)
# summary(modiv1, cluster = TRUE, robust = TRUE)


## IV 2 - in a block, an average neighbor's treated richness IV 
# log
## IV 2 - in a block, an average neighbor's treated richness IV 
modiv2 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
                 (log(rich) ~  log(avg.trt.neigh.rich.within.block)) | newplotid , # | newplotid + site.by.yeardummy), # first stage   Note the surrounding parentheses
               data = comb)
summary(modiv2, cluster = TRUE, robust = TRUE)
# check weak instrument test
summary(modiv2$stage1, cluster = TRUE, robust = TRUE)

#should this include | newplotid + site.by.yeardummy in the firsr stag?
modiv2 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
                 (log(rich) ~  log(avg.trt.neigh.rich.within.block)) | newplotid  + site.by.yeardummy, # first stage   Note the surrounding parentheses
               data = comb)
summary(modiv2, cluster = TRUE, robust = TRUE)
# check weak instrument test
summary(modiv2$stage1, cluster = TRUE, robust = TRUE)


#level for IV
# modiv2 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
#                  (log(rich) ~   avg.trt.neigh.rich.wi
 # thin.block) | newplotid , # | newplotid + site.by.yeardummy), # first stage   Note the surrounding parentheses
#                data = comb)
# summary(modiv2, cluster = TRUE, robust = TRUE)

# IV 2 - BOTH IVs. in a block, an average neighbor's treated richness IV 
modiv3 <- felm(log(live_mass) ~ 0 | newplotid + site.by.yeardummy |
                 (log(rich) ~ log(avg_neighbor_rich) + log(avg.trt.neigh.rich.within.block)) | newplotid , # | newplotid + site.by.yeardummy), # first stage   Note the surrounding parentheses
               data = comb)
summary(modiv3, cluster = TRUE, robust = TRUE)
# check weak instrument test
summary(modiv3$stage1, cluster = TRUE, robust = TRUE)

### Plot Results #plotting coefficient estimates from felm objects:
coefs_modiv1 <- tidy(modiv1, conf.int = T, robust = T)
coefs_modiv2 <- tidy(modiv2, conf.int = T, robust = T)
coefs_modiv3  <- tidy(modiv3, conf.int = T, robust = T)

bind_rows(
  coefs_modiv1 %>% mutate(reg = "IV 1"),
  coefs_modiv2  %>% mutate(reg = "IV 2: treated neighbor rich"),
  coefs_modiv3  %>% mutate(reg = "Model 3 both IVs")
) %>%
  ggplot(aes(x=reg, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-1.1, .6) +
  labs(
    title = "Effect size of Log Species Richness on Log Productivitiy",
    caption = ""
  ) +
  theme(axis.title.x = element_blank())

## IV 2 - in a block, an average neighbor's treated richness IV 
modivT.levels<- ivreg(log(live_mass) ~ log(rich) + newplotid + site_code:year 
                             | log(avg.trt.neigh.rich.within.block) + newplotid +  site_code:year, data = comb)
summary(modivT.levels, diagnostics=T)  # to see relevance tests


#IV 3 in a site an average neighbor's treated richness IV 
modivSite<- ivreg(log(live_mass) ~ log(rich) + newplotid + site_code:year 
                  | log(avg.trt.neigh.rich.within.site) + newplotid +  site_code:year, data = comb)
summary(modivSite, diagnostics=T)  # to see relevance tests


#levels
modivTlev <- ivreg(live_mass ~ rich + newplotid + site_code:year 
             | avg.trt.neigh.rich.within.block + newplotid +  site_code:year, data = comb)
summary(modivTlev, diagnostics=T)  # to see relevance tests


#output IV results into a single table
screenreg(list(clus.res.modiv1lev$cl.res, clus.res.modivTreatedlev$cl.res, clus.res.modiv1lev.groundpar$cl.res, clus.res.modivTlev.groundpar$cl.res),       # object with results from clx
          custom.model.names=c("IV 1 Levels", "IV 2 Levels", "IV 1 Levels'","IV 2 Levels'"),
          omit.coef=c("(site_code)|(newplotid)"))  # object from estimation (unclustered) for BIC

###########################################################################################################
#### Plotting main results for SR all on one plot  ###########################################################
############################################################################################################


###############################################################
### Main Text Figure 3 Plot - Robustness Results ################
##############################################################
#colors # scale_fill_brewer(palette="Greys")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Oster Results from Paul Ferraro from STATA
coef.oster = tibble(term="log(rich)",
                     estimate= -0.19962, 
                    std.error = 0,
                    conf.low = -0.19962,
                    conf.high= -0.19962)

#mechanism blocking analysis to test and control for reverse causality that flows through shading driven 
#by productivity, from Paul Ferraro's STATA analysis
#PUNCHLINE: NO CHANGE IN ESTIMATE AND THUS NO INDICATION OF REVERSE CAUSALITY.
#(d)	Might we not instead just use the regression of P on R and add ground_par?  
#If the estimated negative coefficient on R were coming from reverse causality that is mediated by ground_par, 
#then sticking ground_par into the model would eliminate the estimated negative effect of R on P, no?  
#In other words, the same logic we used for the regression of R on P should also hold for our original regression of 
# P on R, no?  This number can go in Fig 3 in the manuscript as another robustness: 

##* From STATA:   l_rich |  Beta = -.254443   SE = .0917572  
##*    t-stat = -2.77  pval = 0.006  Conf_lower =   -.435808  conf_upper =  -.073078
##*    
coef.mechanism = tibble(term="log(rich)",
                      estimate=  -.254443 , 
                      std.error = .0917572 ,
                      conf.low = -.435808,
                      conf.high= -.073078)

### Plot Results #plotting coefficient estimates from felm objects:
coefs_ModLD <- tidy(ModLD, conf.int = T, robust = T)
#coefs_modiv1 <- tidy(modiv1, conf.int = T, robust = T)
coefs_ModPFE <- tidy(ModPFE, conf.int = T, robust = T)
coefs_modiv2 <- tidy(modiv2, conf.int = T, robust = T) #ave treated neighbor richness within the block 

#to pull out just the richness term for plotting with other models
coefs_ModLD <- filter(coefs_ModLD, term == "log(rich)") 

# try to put all models on one line but group them
Figure3.data <-  bind_rows(
  coefs_ModPFE %>% mutate(reg = "Richness Model"),
  coefs_ModLD %>% mutate(reg = "Lagged-Dependent Variable Model"),
   coef.mechanism %>% mutate(reg = "Mechanism Blocking to test for reverse causality"),
  coefs_modiv2   %>% mutate(reg = "Instrumental Variables Regression"),
  coef.oster %>% mutate(reg = "Sensitivity Test")
) 

Figure3.data$term = factor(Figure3.data$term,
                             levels=c("log(rich)", 
                                     "log(simpson)",
                                     "ihs(even)"))
                    
Figure3 <- ## Just the models with richness in it:
  bind_rows(
    coefs_ModPFE %>% mutate(reg = "1. Our Main Study Design"),  #updated to call Our main design vs Panel FE Model or 2-way fixed effect 
    coef.oster %>% mutate(reg = "3. Sensitivity Test"),
     coefs_ModLD  %>% mutate(reg = "2. Dynamic Panel Design"),  # Lagged-Dependent Variable Model
    coef.mechanism %>% mutate(reg = "4. Mechanism Blocking Design"), # 4. Mechanism Blocking: a test for reverse causality"),
   coefs_modiv2 %>% mutate(reg = "5. Instrumental Variable Design"),
    # coefs_modiv1  %>% mutate(reg = "IV 2:  neighbor rich"),
    # coefs_modiv3  %>% mutate(reg = "Model 3 both IVs")
  ) %>%
  # mutate(name = fct_relevel(reg, "Lagged-Dependent Model", "Panel FE Model")) %>%

 # ggplot(aes(x=reg, y=estimate, ymin=conf.low, ymax=conf.high, colour = reg)) +
  ggplot(aes(x=reg, y=estimate, ymin= estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = reg)) +
  geom_pointrange(size =1.5) + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black" ) + # , size = 1) +
  #geom_hline(yintercept = .2, col = "grey", linetype = "dotdash", size = 1.1) +
  ylim(-1.1, 1.1) +
   scale_colour_discrete(name="Model") +
  labs(
   # title = "Effect size of ln Species Richness on ln Productivity",
    caption = "", x = "Design", y = "Estimated effect size")  +   theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(face="bold", size = 18))

# Estimated % Change in Productivity from a 1% Change in Diversity
Figure3

# Figure3 + scale_fill_brewer(palette="Greys")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Figure3 <- Figure3 + scale_color_manual(values=cbPalette[c(7,6,2, 3, 1)])  +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + 
  theme(axis.title.y= element_text(size=20)) + theme(axis.title.x= element_text(size=20))

# adjusting the legend 
Figure3 <- Figure3 +  theme(legend.position = c(0.4, 0.8)) + theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))

# adjusting the text size for the axis and the main title 
# https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
Figure3 <- Figure3 + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))

#remove the model names on x-axis (too busy with long names)
Fig3 <- Figure3 + theme(axis.text.x = element_blank()) + theme(legend.title = element_blank()) 

# print final figure:
Fig3
