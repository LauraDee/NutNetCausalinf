
######################################################################################################################
### Native vs Non-Native Richness (in SM only) ################################################################
##################################################################################################################

SpecMod_NatvNonNat     <- feols(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) | newplotid + site.by.yeardummy, mech.data) 
SpecMod_NatvNonNatEven <- feols(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + ihs(even) | newplotid + site.by.yeardummy, mech.data) 
SpecMod_4Cats          <- feols(log(live_mass) ~ ihs(sr_INT) + ihs(sr_nat_dom) +  ihs(sr_nat_sub) + ihs(sr_nat_rare) | newplotid + site.by.yeardummy, mech.data) 
SpecMod_RICHNatvNonNat <- feols(log(rich) ~ ihs(sr_INT) + ihs(sr_NAT) | newplotid + site.by.yeardummy, mech.data) 

linearHypothesis(SpecMod_NatvNonNat, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)", 
                 test = "F", singular.ok = T,
                 vcov = vcov(SpecMod_NatvNonNat, cluster = "newplotid"))

linearHypothesis(SpecMod_NatvNonNatEven, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)", 
                 test = "F", singular.ok = T,
                 vcov = vcov(SpecMod_NatvNonNatEven, cluster = "newplotid"))


linearHypothesis(SpecMod_4Cats, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_nat_dom)", 
                 test = "F", singular.ok = T,
                 vcov = vcov(SpecMod_4Cats, cluster = "newplotid"))

esttex(SpecMod_NatvNonNat,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S10_R_se.tex")

esttex(SpecMod_NatvNonNat,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S10_R_ci.tex")








################################################
## Plot Figure Native vs NonNative

invcols <- c("gray69", "saddlebrown")

# Prep Data
FigNatvNonNat.data <-  tidy(SpecMod_NatvNonNat) %>%
  mutate(reg = "Native vs Non-Native")

# Plot

FigNatvNonNat.plot <- FigNatvNonNat.data %>%
  ggplot(aes(x=term, y=estimate, ymin = conf.low, ymax = conf.high, colour = term)) +
  geom_pointrange(aes(col = reg), size = 1.5, position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="term") +
  scale_color_manual(values=invcols[c(1,2)]) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text=element_text(size=22),
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5) ) + 
  geom_hline(yintercept = 0, col = "black") +
  ylim(-.4, .2) +
  scale_x_discrete(labels = c('Non-native','Native'))  %>%
  labs(title = "A. Effect size of Native vs Invasive Species Richness on Productivity",
       caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size"
  )

FigNatvNonNat.plot
ggsave("./output/FigNatvNonNat.pdf", FigNatvNonNat.plot)





##################################################################################
## Model 1 Robustness Analyses for SI #############################################
###################################################################################

#control for total cover and evenness
Mod1A.2 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + ihs(even) + total_cover  | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1A.2, robust = TRUE, cluster = TRUE)

#control for total cover only
Mod1A.3 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + total_cover  | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1A.3, robust = TRUE, cluster = TRUE)

#B. Log-Level
Mod1B <- felm(log(live_mass) ~ sr_INT + sr_NAT | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1B, robust = TRUE, cluster = TRUE)

Mod1Bb <- felm(log(live_mass) ~ sr_INT + sr_NAT + ihs(even) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1Bb, robust = TRUE, cluster = TRUE)

#C. Levels #B. Log-Level
Mod1C <- felm(live_mass ~ sr_INT + sr_NAT | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1C, robust = TRUE, cluster = TRUE)

Mod1Cb <- felm(live_mass ~ sr_INT + sr_NAT + ihs(even) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1Cb, robust = TRUE, cluster = TRUE)

#print results

screenreg(list(Mod1A, Mod1A.1, Mod1B, Mod1Bb, Mod1C, Mod1Cb),     # object with results 
          #custom.model.names= "Native vs Non-Native",
          omit.coef=c("(site_code)|(newplotid)")) 

#output models into a single table
# screenreg(list(Mod1A, Mod1B ,
#                clus.res.PlotFEs.SiteYear.Log.groundPAR, clust.res.PlotFEs.SiteYear.Levels.groundPAR,
#                clus.res.PlotFEs.SiteYear.LevelLog.groundPAR),       # object with results from clx
#           custom.model.names=c("Log-Log'", "Levels'", "Level-Log'","Log-Log''", "Levels''", "Level-Log''" ),
#           omit.coef=c("(site_code)|(newplotid)"))  # object from estimation (unclustered) for BIC

# linearHypothesis(PlotFEs.SiteYear.Log.Log.Orgin, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)",
#                  test = "F", 
#                  vcov. = clust.res.PlotFEs.SiteYear.Log.Log.Orgin,
#                  singular.ok = T, # ignore fact that some variables are getting dropped implicitly
#                  white.adjust = "hc1") # heteroskedasticity robust F-Test 
# #read about the corrected SEs in this here - https://www.rdocumentation.org/packages/car/versions/3.0-2/topics/hccm
# # https://www.econometrics-with-r.org/7-3-joint-hypothesis-testing-using-the-f-statistic.html

### Plot Results
#plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering
coefs_Mod1A <- tidy(Mod1A, conf.int = T, robust = T)
coefs_Mod1A.1 <- tidy(Mod1A.1, conf.int = T, robust = T)
coefs_Mod1A.2 <- tidy(Mod1A.2, conf.int = T, robust = T)
coefs_Mod1A.3 <- tidy(Mod1A.3, conf.int = T, robust = T)

nvnn <- bind_rows(
  coefs_Mod1A %>% mutate(reg = "Native vs Non-Native"),
  coefs_Mod1A.1 %>% mutate(reg = "Native vs Non-Native Incld Evenness"),
  coefs_Mod1A.2 %>% mutate(reg = "Native vs Non-Native Incld Evenness & Total Cover"),
  coefs_Mod1A.3 %>% mutate(reg = "Native vs Non-Native Incld Total Cover")
) %>%
  #ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(aes(col = reg), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  # labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.4, .2) +
  labs(
    title = "Effect size of Native vs Invasive Species Richness on Productivitiy",
    caption = ""
  ) 

nvnn + labs(
  # title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate")

######################################################################################################################
### 2. Rare vs Subordinate vs. Dominant (in SM only) ########################################################################
##################################################################################################################
# Do analyses based on DI then based on relative abundance and frequency groups seperately (versus DI in A.)

###########
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of 1
##########

##### Models with productivity as the Y variable ##########
#log-log
Mod2A.1 <- felm(log(live_mass) ~ ihs(sr_domspp) + ihs(sr_rarespp) + ihs(sr_subordspp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod2A.1, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_subordspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_subordspp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)

################################################################################################
### Plot Figure 4 B ########################################################################
################################################################################################
coefs_Mod2A.1<- tidy(Mod2A.1 , conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4b.data <-  bind_rows(
  coefs_Mod2A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4b.data$term = factor(panelFE.Fig4b.data$term,
                                 levels=c( "ihs(sr_rarespp)", 
                                           "ihs(sr_domspp)" ,
                                           "ihs(sr_subordspp)" ))

Fig4B <-  #ggplot(panelFE.Fig4b.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4b.data, aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", 
) 
Fig4B 

#Alternative y-axis label:
Fig4B < + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")
#panelFE.main.2  + scale_color_manual(values=cbPalette[c(7,3,4,8)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

p <- Fig4B  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="Model") + scale_color_manual(values=cbPalette[c(7,9,4,8)])   +  labs(
  # title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
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
ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('ihs(rare species)','ihs(dominant species)', 'ihs(subordinate species)')) 

#alt simplified variable labels on x-axis
Fig4b <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) +
  scale_x_discrete(labels = c('Rare','Dominant', 'Subordinate')) 
Fig4b 
#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))


# do model with cut-off 2
Mod2A.DI2 <-felm(log(live_mass) ~ ihs( sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod2A.DI2, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(sr_domspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(sr_subordspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_subordspp2) = ihs(sr_domspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)

# print results for cut off 1 and 2:
screenreg(list( Mod2A.1, Mod2A.DI2),      # object with results 
          custom.model.names= c("Cutoff 1", "Cutoff 2")) #name the models

#log-level
Mod2B.2 <- felm(log(live_mass) ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.2, robust = TRUE, cluster = TRUE)

#levels
Mod2B.3 <- felm(live_mass ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.3, robust = TRUE, cluster = TRUE)

screenreg(list(Mod2B.1, Mod2B.2, Mod2B.3),      # object with results 
          #custom.model.names= "Native vs Non-Native"
) 

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Models with richness as the Y variable - how do these vars affect overall richness? ##########
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#log-log
RichMod2A.1 <- felm(log(rich) ~ ihs(sr_domspp) + ihs(sr_rarespp) + ihs(sr_subordspp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.1, robust = TRUE, cluster = TRUE)

#log-log for cut off 2
RichMod2A.DI2 <-  felm(log(rich) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.DI2, robust = TRUE, cluster = TRUE)

# print results for cut off 1 and 2:
screenreg(list(RichMod2A.DI2, RichMod2A.1 ),      # object with results 
          custom.model.names= c("Richness as Response: Cutoff 2", "Richness as Response: Cutoff 1"))

#log-levels
RichMod2A.2 <- felm(log(rich) ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.2, robust = TRUE, cluster = TRUE)

#levels
RichMod2A.3 <- felm(rich ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.3, robust = TRUE, cluster = TRUE)


###########
### AA. Grouped based on the Dominance Indicator (DI) and cutoffs of 2
##########
## do for Cut off 2 # sr_domspp2 sr_rarespp2 sr_subordspp2

##### Models with productivity as the Y variable ##########
#log-log
Mod2AA.1 <- felm(log(live_mass) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2AA.1, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2AA.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_subordspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_subordspp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Models with richness as the Y variable - how do these vars affect overall richness? ##########
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#log-log
RichMod2AA.1 <- felm(log(rich) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2AA.1, robust = TRUE, cluster = TRUE)

########################################################################################
### B. Grouped based on Relative Abundance in year 0 and cutoffs 1 of:  #################
#######################################################################################
## run models with groups defined based on relative abundance
# results are consistent with the results using the DI metrics for groupings of species

Mod2B.1  = felm(log(live_mass) ~ ihs(relabund_sr_domspp) + ihs(relabund_sr_rarespp) + ihs(relabund_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.1, robust = TRUE, cluster = TRUE)

#how does rare spp richness after SR?
RichMod2B.1 = felm(log(rich) ~ ihs(relabund_sr_domspp) + ihs(relabund_sr_rarespp) + ihs(relabund_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2B.1, robust = TRUE, cluster = TRUE)

# do linear hypothesis tests
linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_domspp) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_subordspp) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_subordspp) = ihs(relabund_sr_domspp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

##############################################################################################################
### C. Grouped based on Relative Frequency in year 0 and cutoff 1 of: 0, .2, .8, 1
################################################################################################################
## run models with groups defined based on frequency
# results are consistent with the results using the DI metrics for groupings of species
Mod2C.1  = felm(log(live_mass) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2C.1 , robust = TRUE, cluster = TRUE)

# Now linear hypotheses tests to assess whether the species' group effect on live biomass is equal
# Do linear hypothesis test - are the effects of rare and dominant species the same?
#linearHypothesis function in car package

linearHypothesis(Mod2C.1 , hypothesis.matrix = "ihs(freq_sr_domspp) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2C.1, hypothesis.matrix = "ihs(freq_sr_subordspp) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. = Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2C.1, hypothesis.matrix = "ihs(freq_sr_domspp) = ihs(freq_sr_subordspp)",
                 test = "F", 
                 vcov. = Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

## now assess how rare spp affects SR (for mechanisms tests)
#how does rare spp richness after SR?
RichMod2C.1  = felm(log(rich) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) | newplotid + site.by.yeardummy,  data = mech.data)
summary(RichMod2C.1 , robust = TRUE, cluster = TRUE)

## putting in other controls
# control for N-fixer cover 
Mod2C.2 = felm(log(live_mass) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) + ihs(N_fixer_cover.yr) + ihs(sr_Nfixer) | newplotid + site.by.yeardummy,  data = mech.data)
summary(Mod2C.2 , robust = TRUE, cluster = TRUE)

# are SR rare and SR N-fixing highly correlated? No. r = 0.187
cor(mech.data$sr_Nfixer, mech.data$freq_sr_rarespp)
# 0.1878382
plot(mech.data$sr_Nfixer, mech.data$freq_sr_rarespp, xlab = "SR of N-fixing species", ylab = "SR of rare species")

cor(mech.data$N_fixer_cover.yr, mech.data$freq_sr_rarespp)
plot(mech.data$N_fixer_cover.yr, mech.data$freq_sr_rarespp)

#what about with lagged N-fixers? 
plot(mech.data$lagged_sr_Nfixer, mech.data$freq_sr_rarespp)
cor(mech.data$lagged_sr_Nfixer, mech.data$freq_sr_rarespp)

######################################################################################################################
### 3. Rare vs Non-Rare (Aggregated) - in SM  #############################################################
####################################################################################################

########################################################################################
### A. Grouped based on the Dominance Indicator (DI) and cutoff 1 of:  0, .2, .8, 1
#######################################################################################

##### Models with productivity as the Y variable ##########
#log-log
Mod3A.1 <- felm(log(live_mass) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3A.1, robust = TRUE)

#hypothesis test 
linearHypothesis(Mod3A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_non_rare_spp)", 
                 test = "F", vcov = Mod3A.1$fevcov,  singular.ok = T)

#log-level
Mod3A.2 <- felm(log(live_mass) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy|  0 | newplotid, data = mech.data)
summary(Mod3A.1, robust = TRUE)

##### Models with richness as the Y variable- how do these vars affect overall richness? ##########
#log-log
RichMod3A.1 <- felm(log(rich) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(RichMod3A.1, robust = TRUE)

#print results from richness model
screenreg(list(RichMod3A.1),      # object with results 
          custom.model.names= c("Richness model as outcome, Cut off 1"))

#print results 
screenreg(list(Mod3A.1, Mod3AA.1 ),      # object with results 
          custom.model.names= c("Cut off 1", "Cut off 2"))

###################################################################################################################################
#### Plot Fig. 4D ##################################################################################################################
####################################################################################################################################
coefs_Mod3A.1<- tidy(Mod3A.1 , conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4d.data <-  bind_rows(
  coefs_Mod3A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4d.data$term = factor(panelFE.Fig4d.data$term,
                                 levels=c( "ihs(sr_rarespp)", 
                                           "ihs(sr_non_rare_spp)" ))

Fig4d <- # ggplot(panelFE.Fig4d.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4d.data , aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig4d

#Alternative y-axis label:
Fig4d <- Fig4d + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")
#panelFE.main.2  + scale_color_manual(values=cbPalette[c(7,3,4,8)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

#ggplot2 colors to pick from: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
palette <- c("darkslateblue", "green4")

p <- Fig4d  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="term") + scale_color_manual(values=palette[c(1,2)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
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

#simplified variable labels on x-axis
Fig4d <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('Rare species','Non-rare species')) 

Fig4d <-Fig4d + labs(title="B. Rare vs Non-Rare Richness")
Fig4d 
#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))

########################################################################################
### AA. Grouped based on the Dominance Indicator (DI) and cutoffs 2
#######################################################################################

##### Models with productivity as the Y variable ##########
#log-log
Mod3AA.1 <- felm(log(live_mass) ~ ihs(non_rare_spp.DI2) + ihs(sr_rarespp2)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3AA.1, robust = TRUE)

#hypothesis test 
linearHypothesis(Mod3AA.1, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(non_rare_spp.DI2)", 
                 test = "F", vcov = Mod3AA.1$fevcov,  singular.ok = T)

#print results 
screenreg(list(Mod3A.1, Mod3AA.1 ),      # object with results 
          custom.model.names= c("Cut off 1", "Cut off 2"))

###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:
##########
## run models with groups defined based on relative abundance
# results are consistent with the results using the DI metrics for groupings of species
Mod3B.1  = felm(log(live_mass) ~ ihs(sr_non_rare_spp.RelA) + ihs(relabund_sr_rarespp)  |newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3B.1, robust = TRUE)

#how does rare spp richness after SR?
RichMod3B.1 = felm(log(rich) ~ ihs(sr_non_rare_spp.RelA) + ihs(relabund_sr_rarespp)  |newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod3B.1, robust = TRUE, cluster = TRUE)

# do linear hypothesis tests
linearHypothesis(Mod3B.1 , hypothesis.matrix = "ihs(sr_non_rare_spp.RelA) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod3B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

###########
### C. Grouped based on Relative Frequency in year 0 and cutoff 1 of:
##########
## Run models with groups defined based on frequency
# non_rare_spp.Freq
# results are consistent with the results using the DI metrics for groupings of species
Mod3C.1  = felm(log(live_mass) ~ ihs(sr_non_rare_spp.Freq) + ihs(freq_sr_rarespp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod3C.1 , robust = TRUE)

# Now linear hypotheses tests to assess whether the species' group effect on live biomass is equal
# Do linear hypothesis test - are the effects of rare and dominant species the same?
linearHypothesis(Mod3C.1 , hypothesis.matrix = "ihs(sr_non_rare_spp.Freq) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod3C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

######################################################################################################################
### 4. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################

###########
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# these are the variable names/groups:
# sr_non.rare_nat  + sr_non.rare_non.nat + sr_non.nat_rare +  sr_nat_rare
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4A.1, robust = TRUE)

## Hypothesis tests
#Note: fevcov returns a square matrix with the bias corrected covariances. An attribute 'bias' contains
# the biases. The bias corrections have been subtracted from the bias estimates. I.e. vc = vc’ - b,
# where vc’ is the biased variance and b is the bias.

# not rare: native vs non-native 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

# Native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native vs native rare 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#hypothesis tests using lefm objective and a wald test:
# https://www.rdocumentation.org/packages/lfe/versions/2.8-3/topics/waldtest
#waldtest(object, R, r, type = c("default", "iid", "robust", "cluster"),
#lhs = NULL, df1, df2)
# hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)"
# R = hypothesis.matrix 
# waldtest(Mod4A.1, ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat),  type = c("robust", "cluster"))
#    , lhs = NULL, df1, df2)

#### Richness as Outcome Variable 
ModRich4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(ModRich4A.1, robust = TRUE)


###################################################################################################################################
### Plot Figure 4C ######################################################################################################################
#####################################################################################################################################
### Plot Results - Plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

coefs_Mod4A.1 <- tidy(Mod4A.1, conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4C.data <-  bind_rows(
  coefs_Mod4A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4C.data$term = factor(panelFE.Fig4C.data$term,
                                 levels=c( "ihs(sr_nat_rare)", 
                                           "ihs(sr_non.rare_nat)",
                                           "ihs(sr_non.rare_non.nat)",
                                           "ihs(sr_non.nat_rare)"))

Fig4C <-  # ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="term") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.7, .7), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig4C

#Alternative y-axis label:
Fig4C <- Fig4C + labs(
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size")
Fig4C

#ggplot2 colors to pick from: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
palette <- c("darkslateblue", "green4", "grey69", "maroon4" )

p <- Fig4C  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="term") + scale_color_manual(values=palette[c(1,2, 3,4)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
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
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) 
ppp

# print final figure:
ppp <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('ihs(Rare native)','ihs(Non-rare native)', 'ihs(Non-rare Non-native', 'ihs(Rare Non-native)')) 
ppp

#alt simplified variable labels on x-axis
Fig4C <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('Rare & Native','Non-rare & Native', 'Non-rare & Non-native', 'Rare & Non-native')) 
Fig4C

Fig4C <- Fig4C + labs(title="C. Rare Native vs. Non-rare Native vs. Non-rare non-native vs Rare non-native") +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0)) 
#print final
Fig4C

#simplify title
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 
#print final
Fig4C
#title to the left
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = -0.1)) 
#print final
Fig4C


#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))

################################################################################################################################################################################################################################
### Put Fig 4 plots in same plot #####################################################################################################
#####################################################################################################
# first make the titles A, B, C, D and left align them

#A invasive vs native
#simplified plot labels
Fig4A <- inv + labs( title =  "A. Native vs Invasive Species Richness",
                     caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") 
+ theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

Fig4A  <- Fig4A  + labs (title =  "A.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#rare vs non rare - B
Fig4d <- Fig4d  + labs(title="B.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#  C. dom, sub, rare 
Fig4b <- Fig4b + labs(title="C.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#all groups -D
Fig4C <- Fig4C + labs(title="D.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 

##print final
plot_grid(Fig4A, Fig4d, Fig4b, Fig4C)


#########################################
# plot only three of the plots ###########
############################################
Fig4d <- Fig4d  + labs(title="B.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 

#print just the Figure 4 C
Fig4C <- Fig4C + labs(title="") +   theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 
Fig4C

#print final
plot_grid(Fig4A, Fig4d,  Fig4C)

#Messing around with plot grid to make the plots look better together:
Fig4Apg  <- Fig4A  + labs (title =  "") 
Fig4dpg <- Fig4d  + labs(title="")
Fig4Cpg <- Fig4C + labs(title="")


### FINAL:


## put the three on the same and make
common.ylim = ylim(-0.4, 0.4)
kLabelSize = 25 #Add labels with plotgrid 
common.ylab = ylab("Estimated effect of species richness")  #Estimated coefficient of species richness
plot_grid(
  plot_grid(Fig4Apg + common.ylim + common.ylab, 
            Fig4dpg + common.ylim + common.ylab, 
            labels=c("A","B"), label_size=kLabelSize),
  plot_grid(NULL, Fig4Cpg + common.ylim + common.ylab, NULL, labels = c("", "C", ""), label_size=kLabelSize, rel_widths = c(0.2, 0.6, 0.2), nrow=1), 
  nrow = 2)


#tutorial - and way to add joint plot title https://wilkelab.org/cowplot/articles/plot_grid.html

##########################################################################################
#### Robustness analyses for fig 4C analyses ###################################################
##########################################################################################

#need to run some robustness checks with cover variables too
#****cover_tot_dom and cover_tot_sub didn't work in data processing
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) +  cover_tot_INT + cover_tot_NAT + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod4A.1, robust = TRUE, cluster = TRUE)

##### Models with richness as the Y variable- how do these vars affect overall richness? ##########
RichMod4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) +  cover_tot_INT + cover_tot_NAT + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                    +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod4A.1, robust = TRUE, cluster = TRUE)

RichMod4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                    +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod4A.1, robust = TRUE, cluster = TRUE)

########## breakng non-rare native into SR of native dom and native subordinate spp.
Mod4A.2 <- felm(log(live_mass) ~ ihs(sr_nat_dom) + ihs(sr_nat_sub) +  ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod4A.2, robust = TRUE, cluster = TRUE)


###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# sr_rare_non.nat.RelA sr_rare_nat.RelA sr_rare_non.nat.Freq sr_rare_nat.Freq

Mod4B.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.RelA) + ihs(sr_non.rare_non.nat.RelA)  + ihs(sr_rare_non.nat.RelA) 
                +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4B.1, robust = TRUE)

linearHypothesis(Mod4B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA) = ihs(sr_non.rare_non.nat.RelA)", 
                 test = "F", vcov = Mod4B.1$fevcov,  singular.ok = T)

linearHypothesis(Mod4B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA) = ihs(sr_rare_nat.RelA)", 
                 test = "F", vcov = Mod4B.1$fevcov,  singular.ok = T)

## Richness as the response
RichMod4B.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.RelA) + ihs(sr_non.rare_non.nat.RelA)  + ihs(sr_rare_non.nat.RelA) 
                    +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod4B.1, robust = TRUE, cluster = TRUE)

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0)
##########
Mod4C.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.Freq) + ihs(sr_non.rare_non.nat.Freq)  + ihs(sr_rare_non.nat.Freq) +  ihs(sr_rare_nat.Freq) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4C.1, robust = TRUE)

#make table of results
screenreg(Mod4C.1, custom.model.names= "Rel. Freq Cut off 1")

# hypothesis tests

#native vs non-native, non-rare
linearHypothesis(Mod4C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq) = ihs(sr_non.rare_non.nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

#non rare native vs rare native 
linearHypothesis(Mod4C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq) = ihs(sr_rare_nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

# rare native vs rare non-native
linearHypothesis(Mod4C.1, hypothesis.matrix = " ihs(sr_rare_non.nat.Freq)  = ihs(sr_rare_nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

# non-native rare vs non-rare  
linearHypothesis(Mod4C.1, hypothesis.matrix = " ihs(sr_rare_non.nat.Freq)  = ihs(sr_non.rare_non.nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

#richness as response
RichMod4C.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.Freq) + ihs(sr_non.rare_non.nat.Freq)  + ihs(sr_rare_non.nat.Freq) 
                    +  ihs(sr_rare_nat.Freq) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod4C.1, robust = TRUE, cluster = TRUE)

### Plot Results Using the Different ways of defining rarity 

coefs_Mod4A.1 <- tidy(Mod4A.1, conf.int = T, robust = T)
coefs_Mod4B.1 <- tidy(Mod4B.1, conf.int = T, robust = T)
coefs_Mod4C.1 <- tidy(Mod4C.1, conf.int = T, robust = T)

panelMech.all <-  bind_rows(
  coefs_Mod4A.1 %>% mutate(reg = "Model Using Dominance Indicator"),
  coefs_Mod4B.1 %>% mutate(reg = "Model Using Relative Abundance"),
  coefs_Mod4C.1 %>% mutate(reg = "Model Using Frequency")
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange() + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  scale_colour_discrete(name="Grouping Definition") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness By Group on Log Productivitiy",
    caption = "" ) + facet_wrap(~reg) 
#theme(axis.title.x = element_blank())

panelMech.all + labs(
  # title = "Effect size of Log Group Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

#####################################################################################################################
#### 5.Run Models from 4 with different grouping cutoffs of breaks=c(0.0,0.4,0.8,1.0).
# Comparing Rare vs Non-rare, Native or Invasive ###################################################################################
#####################################################################################################################

#####################
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of: breaks=c(0.0,0.4,0.8,1.0),
####################
#*for some reason this is now showing up so that the first two variables are 0..
Mod5A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat2) +   ihs(sr_non.rare_non.nat2)  + ihs(sr_non.nat_rare2)  +  ihs(sr_nat_rare2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5A.1, robust = TRUE)

#make table of results
screenreg(Mod5A.1, custom.model.names= "Dominance Indicator: Cut off 2")
#make table of results
screenreg(Mod4A.1, custom.model.names= "Dominance Indicator: cut off 1")
#put both in one table
screenreg(list(Mod4A.1, Mod5A.1), custom.model.names= c("Dominance Indicator: cut off 1","Dominance Indicator: Cut off 2"))

#hypothesis tests
#non-rare native vs non-native 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.rare_nat2) = ihs(sr_non.rare_non.nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#native rare vs not rare 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_nat_rare2) = ihs(sr_non.rare_nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#non-native: rare vs non rare 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.nat_rare2) = ihs(sr_non.rare_non.nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#rare native vs rare non native
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.nat_rare2) = ihs(sr_nat_rare2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0),
##########
Mod5B.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5B.1, robust = TRUE)

linearHypothesis(Mod5B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA2) = ihs(sr_non.rare_non.nat.RelA2)", 
                 test = "F", vcov = Mod5B.1$fevcov,  singular.ok = T)

linearHypothesis(Mod5B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA2) = ihs(sr_rare_nat.RelA2)", 
                 test = "F", vcov = Mod5B.1$fevcov,  singular.ok = T)

# also control for Nfixers
Mod5B.2 <- felm(log(live_mass) ~ ihs(sr_Nfixer) + ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod5B.2, robust = TRUE, cluster = TRUE)

#Richness as the response
RichMod5B.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                    +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod5B.1, robust = TRUE, cluster = TRUE)

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0)
##########
Mod5C.1 <- felm(log(live_mass) ~  ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5C.1, robust = TRUE)

# hypothesis tests
linearHypothesis(Mod5C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq2) = ihs(sr_non.rare_non.nat.Freq2)", 
                 test = "F", vcov = Mod5C.1$fevcov,  singular.ok = T)

linearHypothesis(Mod5C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq2) = ihs(sr_rare_nat.Freq2)", 
                 test = "F", vcov = Mod5C.1$fevcov,  singular.ok = T)

# also control for Nfixers
Mod5C.2 <- felm(log(live_mass) ~ ihs(sr_Nfixer) + ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod5C.2, robust = TRUE, cluster = TRUE)

#richness as response
RichMod5C.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                    +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod5C.1, robust = TRUE, cluster = TRUE)

### Plot Results Using the Different ways of defining rarity 
coefs_Mod5A.1 <- tidy(Mod5A.1, conf.int = T) 
#%>% as_tibble(rownames = "ihs(sr_non.rare_nat2)") -> Non_Rare.Native
coefs_Mod5B.1 <- tidy(Mod5B.1, conf.int = T)
coefs_Mod5C.1 <- tidy(Mod5C.1, conf.int = T)

panelMech.all_cutoff2 <-  bind_rows(
  coefs_Mod5A.1 %>% mutate(reg = "Model Using Dominance Indicator"),
  coefs_Mod5B.1 %>% mutate(reg = "Model Using Relative Abundance"),
  coefs_Mod5C.1 %>% mutate(reg = "Model Using Frequency")
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange() + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  scale_colour_discrete(name="Grouping Definition") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness By Group on Log Productivitiy",
    caption = "" ) + facet_wrap(~reg) 
#theme(axis.title.x = element_blank())

panelMech.all_cutoff2 + labs(
  # title = "Effect size of Log Group Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

# save plots my_plot <- stops_facet_plot 
ggsave("Fig4C", plot = panelMech.all_cutoff2, width=NA, height=NA)

#make table of results
screenreg(list(Mod5A.1, Mod5B.1, Mod5C.1),     # object with results 
          custom.model.names= c("Dominance Indicator",  "Relative abundance", "Frequency"))
# omit.coef=c("(site_code)|(newplotid)")) 

# Print Results Comparing The Cut-Offs 
#DI
screenreg(list(Mod4A.1, Mod5A.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

# Print Results Comparing The Cut-Offs -- for relative abundance. 
screenreg(list(Mod4B.1, Mod5B.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

# Print Results Comparing The Cut-Offs 
screenreg(list(Mod4C.1, Mod5C.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

############################################################################################################################
#### Plot Correlations between all of the SR groupings ###################################################################
#############################################################################################################################
richnessvars.fig4c <- c("sr_non.rare_nat", "sr_non.rare_non.nat", "sr_non.nat_rare", "sr_nat_rare")
richnessvars.to.plot <- mech.data[,richnessvars.fig4c, with=F] #with = F will then use the data.frame conventions, using "" for col names
pairs(richnessvars.to.plot)

######################################################################################################################
#### Visualize Correlations between SR vars ##########################################################################################################
######################################################################################################################
library(corrplot)
sr.metrics <- mech.data[, .(sr_non.nat_rare, sr_nat_rare, non_rare_spp,
                            sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom,
                            sr_non.nat_dom, sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare )]
cor(sr.metrics)
corrplot(cor(sr.metrics), method = "square", tl.cex = .5)