##############################################
## Models for Figure 2A and Table S2 ETX XXX #
##############################################

MainMod_Rich     <- feols(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy, comb) 
MainMod_RichEven <- feols(log(live_mass) ~ log(rich) + ihs(even) | newplotid + site.by.yeardummy, comb) 
MainMod_Simpson  <- feols(log(live_mass) ~ log(simpson) | newplotid + site.by.yeardummy, comb) 
MainMod_RichLag  <- feols(log(live_mass) ~ log(rich) + log(laggedrich) | newplotid + site.by.yeardummy, comb) 
MainMod_RichEvenLag <- feols(log(live_mass) ~ log(rich) + log(laggedrich) + ihs(even) | newplotid + site.by.yeardummy, comb)

##############################################
## Models for Figure 2B and Table S2 ETX XXX #
##############################################

varnames <- "live_mass rich even simpson newplotid site_code country habitat year elevation managed burned grazed anthropogenic TEMP_VAR_v2 MIN_TEMP_v2 MAX_TEMP_v2 TEMP_WET_Q_v2 TEMP_DRY_Q_v2 TEMP_WARM_Q_v2 TEMP_COLD_Q_v2 pct_C pct_N ppm_P ppm_K ppm_Na ppm_Mg ppm_S ppm_Na ppm_Zn ppm_Mn ppm_Fe ppm_Cu ppm_B pH PercentSand PercentSilt PercentClay"

vv <- unlist(strsplit(varnames," "))
for (i in 1:38) {
  print(paste(i, grep(vv[i],names(ml_comb)), sep = " , "))
}

ml_comb <- as.data.frame(comb)
ml_comb <- ml_comb[,vv]
ml_comb <- ml_comb[complete.cases(ml_comb),]

## Run mixed models


MixedMod_Rich <- lmer(log(live_mass) ~ log(rich) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                        elevation + managed + burned + grazed + anthropogenic + 
                        TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                        pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                        pH + PercentSand + PercentSilt + PercentClay + 
                        (1|newplotid) + (1|site_code), ml_comb, REML = F)

MixedMod_RichEven <- lmer(log(live_mass) ~ log(rich) + ihs(even) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                            elevation + managed + burned + grazed + anthropogenic + 
                            TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                            pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                            pH + PercentSand + PercentSilt + PercentClay + 
                            (1|newplotid) + (1|site_code), ml_comb, REML = F)

MixedMod_Simpson <- lmer(log(live_mass) ~ log(simpson) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                           elevation + managed + burned + grazed + anthropogenic + 
                           TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                           pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                           pH + PercentSand + PercentSilt + PercentClay + 
                           (1|newplotid) + (1|site_code), ml_comb, REML = F)


## Because I couldn't find heteroskedasticity and autocorrelation robust variances with lme4, let's
## read in these more conservative CIs from stata to adjust

Fig2B.1 <- tidy(MixedMod_Rich) %>%
  filter(term == "log(rich)")
Fig2B.2 <- tidy(MixedMod_RichEven) %>%
  filter(term == "log(rich)")
Fig2B.3 <- tidy(MixedMod_Simpson) %>%
  filter(term == "log(simpson)")

Fig2B.1$conf.low <- .0049797 
Fig2B.1$conf.high <- .7535313
Fig2B.2$conf.low <- -.0117337
Fig2B.2$conf.high <- .8204197
Fig2B.3$conf.low <- -.0855739
Fig2B.3$conf.high <- .2348401


################################################
## Table S2
#######################################
etable(MainMod_Rich, MainMod_RichEven, MainMod_RichLag, MainMod_RichEvenLag, MainMod_Simpson,
       coefstat = "se")

etable(MainMod_Rich, MainMod_RichEven, MainMod_RichLag, MainMod_RichEvenLag, MainMod_Simpson,
       coefstat = "confint")

#esttex(MainMod_Rich, MainMod_RichEven, MainMod_RichLag, MainMod_RichEvenLag, MainMod_Simpson,
#      coefstat = c("se", "confint"))

################################################
## Table ?? FOR TRADITIOL RESULTS
#######################################




######################################
## PLOTTING YEAH! ####
###############################

## Master Palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1", "red")

################################################
## Plot Figure 2A

# Prep Data
Fig2A.1 <- tidy(MainMod_Rich)
Fig2A.2 <- tidy(MainMod_RichEven) %>%
  filter(term == "log(rich)")
Fig2A.3 <- tidy(MainMod_Simpson) 

Fig2A.data <-  bind_rows(
  Fig2A.1 %>% mutate(reg = "Species Richness"),
  Fig2A.2 %>% mutate(reg = "Species Richness controlling for Evenness"),
  Fig2A.3 %>% mutate(reg = "Simpson's Diversity") )

Fig2A.data$term = factor(Fig2A.data$term,
                         levels=c("log(rich)", 
                                  "log(simpson)",
                                  "ihs(even)"))

# Plot

Fig2A.plot <- Fig2A.data %>%
  ggplot(aes(x=term, y=estimate, ymin = conf.low, ymax = conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), size = 1.5, position = position_dodge(width = 0.5)) +
  scale_colour_discrete() +
  scale_color_manual(values=cbPalette[c(7,9,4,8)]) +
  theme_classic() +
  theme(legend.position = c(0.6, 0.77),
        legend.title = element_blank(), 
        legend.text  = element_text(size=14),
        legend.background = element_rect(size=0.5, 
                                         linetype="solid",
                                         colour ="black" ),
        axis.text=element_text(size=22),
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5) ) + 
  geom_hline(yintercept = 0, col = "black") +
  ylim(-.7, .8) +
  scale_x_discrete(labels = c("Species Richness", "Simpson's Diversity")) + 
  scale_y_continuous(limits=c(-.8, .85), 
                     breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8)
                     ) %>%
  labs(title = "Main Study Design",
       caption = "", x = "Variable", y = "Estimated effect size"
       )

Fig2A.plot

################################################
## Plot Figure 2B

Fig2B.data <-  bind_rows(
  Fig2B.1 %>% mutate(reg = "Species Richness"),
  Fig2B.2 %>% mutate(reg = "Species Richness controlling for Evenness"),
  Fig2B.3 %>% mutate(reg = "Simpson's Diversity") )

Fig2B.data$term = factor(Fig2B.data$term,
                         levels=c("log(rich)", 
                                  "log(simpson)",
                                  "ihs(even)"))

# Plot

Fig2B.plot <- Fig2B.data %>%
  ggplot(aes(x=term, y=estimate, ymin = conf.low, ymax = conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), size = 1.5, position = position_dodge(width = 0.5)) +
  scale_colour_discrete() +
  scale_color_manual(values=cbPalette[c(8,3,10)]) +
  theme_classic() +
  theme(legend.position = c(.5, 0.3),
        legend.title = element_blank(), 
        legend.text  = element_text(size=14),
        legend.background = element_rect(size=0.5, 
                                         linetype="solid",
                                         colour ="black" ),
        axis.text=element_text(size=22),
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5) ) + 
  geom_hline(yintercept = 0, col = "black") +
  ylim(-.7, .85) +
  scale_x_discrete(labels = c("Species Richness", "Simpson's Diversity")) + 
  scale_y_continuous(limits=c(-.8, .85), 
                     breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8)
  ) %>%
  labs(title = "Traditional Ecological Design",
       caption = "", x = "Variable", y = "Estimated effect size"
  )

Fig2B.plot

#################################
######## COMBINE 

plot_grid(Fig2A.plot , Fig.2B.plot )

common.ylab = ylab("Estimated effect size")  #Estimated % Change in Productivity from a 1% Change in Diversity
plot_grid(Fig2A.plot  + common.ylab,
          Fig2B.plot + common.ylab)


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
