
####################################################################################################################################
#### Models for the role of rare and non-native species  ##########################################################################################################
####################################################################################################################################

# This code includes analysis for reproducing Figure 5 in the main text and SM analyses for Section 7.

######################################################################################################################
###  Figure 5 - Main text. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################
###########
### Figure 5 - Main text. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# We first present the analyses shown in the main text Figure 5, which include 4 groups of species, which are: 
#1) rare, native: sr_nat_rare
#2) rare non-native: sr_non.nat_rare
#3) non-rare, native: sr_non.rare_nat
#4) non-rare, non-native: 
# the analyses presented in the main text Figure 5 use classify rare versus non-rare groups based 
# on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0), where non-rare are
#both subordinate and dominant species (e.g., species with DI greater than the .2 cutoff). 

# 1. Create SR non-rare native and non-native excluding unknown species origin species
MechMod_All <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) + ihs(sr_NA)
                    | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")

vcov_MechMod <- vcov(MechMod_All, cluster = "newplotid")


MechMod_All_noNAs <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) 
                          | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechMod_noNAs <- vcov(MechMod_All_noNAs, cluster = "newplotid")


# Hypotheses Tests for the groups: are their effects on productivity significantly differently?
#Note: fevcov returns a square matrix with the bias corrected covariances. An attribute 'bias' contains
# the biases. The bias corrections have been subtracted from the bias estimates. I.e. vc = vc’ - b,
# where vc’ is the biased variance and b is the bias.

#testing if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechMod_All, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)",
                                       "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)


linearHypothesis(MechMod_All_noNAs, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)",
                                       "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod_noNAs,  singular.ok = T)


# not rare: native vs non-native - are their effects on productivity significantly differently?
linearHypothesis(MechMod_All, hypothesis.matrix = "ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

# native rare vs non-rare:
linearHypothesis(MechMod_All, hypothesis.matrix = "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)", 
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

# non-native rare vs non-rare
linearHypothesis(MechMod_All, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

# non-native vs native rare 
linearHypothesis(MechMod_All, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)", 
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

###################
## Export Table ## 

esttex(MechMod_All, 
       coefstat = "se", replace = TRUE,
       file = "./output/Table_forFig5_R_se-May202021.tex")

esttex(MechMod_All,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_forFig5_R_ci-May202021.tex")


esttex(MechMod_All, MechMod_All_noNAs,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_forFig5_R_se-no-NAs-May202021.tex")


esttex(MechMod_All, MechMod_All_noNAs,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_forFig5_R_ci-no-NAs-May202021.tex")

###################################################################################################################################
### Plot Figure 5 ######################################################################################################################
#####################################################################################################################################

Fig5.data <- tidy(MechMod_All)

Fig5.data <-  bind_rows(
  Fig5.data %>% mutate(reg = "Richness Model"),
) 

Fig5.data$term = factor(Fig5.data$term,
                        levels=c( "ihs(sr_nat_rare)", 
                                  "ihs(sr_non.rare_nat)",
                                  "ihs(sr_non.rare_non.nat)",
                                  "ihs(sr_non.nat_rare)"))

# Plot                    
palette <- c("darkslateblue", "green4", "grey69", "maroon4" )

Fig5.plot <- Fig5.data %>%
  ggplot(aes(x=term, y=estimate, ymin = conf.low, ymax = conf.high, colour = term)) +
  geom_pointrange(size = 1.5) +
  scale_colour_discrete(name="term") +
  scale_color_manual(values=palette[c(1,2,3,4)]) +
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
  scale_x_discrete(labels = c('Rare & Native','Non-rare & Native', 'Non-rare & Non-native', 'Rare & Non-native')) +
  scale_y_continuous(limits=c(-0.7, 0.7)) +
  ylim(-0.7,0.8) %>%
  labs(title = element_blank(),
       caption = "", x = "Species Type", y = "Estimate for log(species richness) effect size"
  )

Fig5.plot
ggsave("./output/Fig5.pdf", Fig5.plot)



###################################################################################################################################################
### SM Analyses for Section S7  ######################################################################################################################
#####################################################################################################################################################

#####################################################################################################################
#### Run Models Comparing with DI metrics for Rare vs Non-rare with different grouping cutoffs of using breaks=c(0.0,0.4,0.8,1.0).
#####################################################################################################################

#####################
### Table S16 Models.  Grouped based on the Dominance Indicator (DI) and cutoffs of: breaks=c(0.0,0.4,0.8,1.0), versus cut-offs from Figure 5
####################

MechMod_All2 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat2) +   ihs(sr_non.rare_non.nat2)  + ihs(sr_non.nat_rare2)  +  ihs(sr_nat_rare2)
                     | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")

vcov_MechMod2 <- vcov(MechMod_All2, cluster = "newplotid")

#testing if rare native has different effect than non-rare native
linearHypothesis(MechMod_All2, 
                 hypothesis.matrix = "ihs(sr_nat_rare2) = ihs(sr_non.rare_nat2)", 
                 test = "F", vcov = vcov_MechMod2,  singular.ok = T)

#test if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechMod_All2, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat2) = ihs(sr_non.rare_non.nat2)", "ihs(sr_nat_rare2) = ihs(sr_non.rare_nat2)",
                                       "ihs(sr_non.nat_rare2) = ihs(sr_nat_rare)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod2,  singular.ok = T)

################################################
## Table S9
#######################################
esttex(MechMod_All, MechMod_All2, 
       coefstat = "se", replace = TRUE,
       file = "./output/TableS16_R_se.tex")

esttex(MechMod_All, MechMod_All2, 
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS16_R_CI.tex")

#########################################################################################################################################################################
## Table S10 and S11 Models. Compare estimates using other metrics for defining rare vs non-rare,  ##################################
# based only on relative abundance and relative frequency for both cut-offs                       ##################################
#########################################################################################################################################################################

###########
### B. Use a metric of rare or non-rare based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0)  (same cut-offs as in Figure 5)
##########

MechRelA1 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat.RelA) + ihs(sr_non.rare_non.nat.RelA)  + ihs(sr_rare_non.nat.RelA) 
                  +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechRelA1 <- vcov(MechRelA1, cluster = "newplotid")

#testing if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechRelA1, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat.RelA) = ihs(sr_non.rare_non.nat.RelA)", "ihs(sr_non.rare_nat.RelA) = ihs(sr_rare_nat.RelA)", 
                                       " ihs(sr_rare_non.nat.RelA)  = ihs(sr_rare_nat.RelA)"), # by transitivity last group is dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov =vcov_MechRelA1,  singular.ok = T)

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0)
##########
MechFreq1 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat.Freq) + ihs(sr_non.rare_non.nat.Freq)  + ihs(sr_rare_non.nat.Freq) +  ihs(sr_rare_nat.Freq) 
               | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechFreq1 <- vcov(MechFreq1, cluster = "newplotid")

################################################
## Table S10
#######################################
esttex(MechRelA1, MechRelA2,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS10_R_se.tex")

esttex(MechRelA1, MechRelA2,
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS10_R_CI.tex")

###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0),
##########

MechRelA2 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                  +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechRelA2 <- vcov(MechRelA2, cluster = "newplotid")

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0)
##########
MechFreq2 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                  +  ihs(sr_rare_nat.Freq2)| newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechFreq2 <- vcov(MechFreq2, cluster = "newplotid")

################################################
## Table S11  #################################
##############################################
esttex(MechFreq1, MechFreq2,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS11_R_se.tex")

esttex(MechFreq1, MechFreq2,
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS11_R_CI.tex")


################################################
## Compare all metrics for each cut-off #########
################################################
#Cut off 1 
esttex(MechMod_All, MechRelA1, MechFreq1,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS9through11_ComparisonCutoff1_R_se.tex")

esttex(MechMod_All, MechRelA1, MechFreq1,
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS9through11_ComparisonCutoff1_R_CI.tex")


#Cut off 2
esttex(MechMod_All2, MechRelA2, MechFreq2,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS9through11_ComparisonCutoff2_R_se.tex")

esttex(MechMod_All2, MechRelA2, MechFreq2,
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS9through11_ComparisonCutoff2_R_CI.tex")


######################################################################################################################################################
#### Sensitivity Analyses: Run Models that categorize species coming into plots after year 0 as native or non-native in different ways #################
######################################################################################################################################################
# We test the sensitivity of our results to data processing decisions, with respect to species that have unknown origins (e.g., site coordinators did not know if the species
# was native or introduced).
## 1. Excluding them, as done in MechMod_All analyses above for results in Figure 5 and in Table S9 (Cut-off1) #### 

## 2. Including the unknown spp origin all as *native* : ####
## to do so, we replace: ihs(sr_non.rare_nat)   with  ihs(sr_non.rare_nat_unk)
## and replace ihs(sr_nat_rare)  with ihs(sr_nat_unk_rare)
MechMod_S2 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat_unk) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_unk_rare) + ihs(sr_NA)
                   | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechModS2 <- vcov(MechMod_S2, cluster = "newplotid")
summary(MechMod_S2)

#**need to update 
#testing if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechMod_S2, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)",
                                       "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)


#without controlling for NA species:
MechMod_S2.noNA <-feols(log(live_mass) ~ ihs(sr_non.rare_nat_unk) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_unk_rare)
                   | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechModS2.noNA <- vcov(MechMod_S2, cluster = "newplotid")


# 3. Including the species with unknown origin all as non-native, using these variables:
# sr_non.nat_unk_rare 
# sr_non.rare_non.nat_unk 

## replace: ihs(sr_non.rare_non.nat)     with  ihs(sr_non.rare_non.nat_unk)
## replace:  ihs(sr_non.nat_rare)    with    ihs(sr_non.nat_unk_rare )
MechMod_S3 <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat_unk) +  ihs(sr_non.nat_unk_rare) +  ihs(sr_nat_rare) + ihs(sr_NA)
                   | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechModS3 <- vcov(MechMod_S3, cluster = "newplotid")
summary(MechMod_S3)

#without controlling for counts of NA species
MechMod_S3.noNA <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat_unk) +  ihs(sr_non.nat_unk_rare) +  ihs(sr_nat_rare) 
| newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechModS3.noNA <- vcov(MechMod_S3, cluster = "newplotid")
summary(MechMod_S3.noNA)

#**update 
#testing if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechMod_S3, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)",
                                       "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

################################################
## Table S12  ##################################
################################################

esttex(MechMod_All, MechMod_S2, MechMod_S3,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS12_SensitivityAnal_R_se.tex")

esttex(MechMod_All, MechMod_S2, MechMod_S3,
       coefstat = "confint", replace = TRUE,
       file = "./output/TableS12_SensitivityAnal_R_CI.tex")

esttex(MechMod_All, MechMod_S2, MechMod_S2.noNA,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS12_SensitivityAnal_R_seMay202021.tex")

esttex(MechMod_All, MechMod_S3, MechMod_S3.noNA,
       coefstat = "se", replace = TRUE,
       file = "./output/TableS13_SensitivityAnal_R_seMay202021.tex")


############################################################################################################################
#### Figure Extra. Plot Correlations between all of the SR grouping variables ###################################################################
#############################################################################################################################
richnessvars.fig4c <- c("sr_non.rare_nat", "sr_non.rare_non.nat", "sr_non.nat_rare", "sr_nat_rare")
richnessvars.to.plot <- mech.data[,richnessvars.fig4c, with=F] #with = F will then use the data.frame conventions, using "" for col names
pairs(richnessvars.to.plot)

sr.metrics <- mech.data[, .(sr_non.nat_rare, sr_nat_rare, non_rare_spp,
                            sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom,
                            sr_non.nat_dom, sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare )]
cor(sr.metrics)
corrplot(cor(sr.metrics), method = "square", tl.cex = .5)





####################################################################################################################
###  Supplemental Figures of changes in species groups by site  ##########################################################
####################################################################################################################
mech.data[order(year), change_sr_nonrare.nonative_spp := sr_non.rare_non.nat -shift(sr_non.rare_non.nat), by =.(plot, site_code)]
mech.data[order(year), change_sr_nonrare.native_spp := sr_non.rare_nat -shift(sr_non.rare_nat), by =.(plot, site_code)]
mech.data[order(year), change_sr_rare.native_spp := sr_nat_rare -shift(sr_nat_rare), by =.(plot, site_code)]
mech.data[order(year), change_sr_rare.nonnative_spp := sr_non.nat_rare -shift(sr_non.nat_rare), by =.(plot, site_code)]


# change nonrare nonnative 
FigSX <- ggplot(data = mech.data, aes(x = change_sr_nonrare.nonative_spp)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in non-native non-rare species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

#native nonrare
FigSX <- ggplot(data = mech.data, aes(x =  change_sr_nonrare.native_spp)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in native non-rare species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

#native rare
FigSX <- ggplot(data = mech.data, aes(x =  change_sr_rare.native_spp)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in native rare species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

#non native rare - change_sr_rare.nonnative_spp 
FigSX <- ggplot(data = mech.data, aes(x =  change_sr_rare.nonnative_spp )) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in non-native rare species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX


## do first differences of the above species richness variables: 

#** this is wrong and shpuld be changed to the mech.data with the unique entry per plot and site and year
setorder(mech.data, year)
mech.data[, change_sr_INT := sr_INT-shift(sr_INT), by =.(plot, site_code)]
mech.data[, change_sr_NAT := sr_NAT-shift(sr_NAT), by =.(plot, site_code)]
mech.data[, change_sr_UNK := sr_UNK-shift(sr_UNK), by =.(plot, site_code)]

# compute change in richness in each group 
cover[order(year), change_sr_domspp := sr_domspp -shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), change_sr_rarespp := sr_rarespp -shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), change_sr_subordspp := sr_subordspp -shift(sr_subordspp), by =.(plot, site_code)]
cover[order(year), change_sr_non_rare_spp := sr_non_rare_spp -shift(sr_non_rare_spp), by =.(plot, site_code)]

summary(cover$change_sr_domspp) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.0000  0.0000  0.0000  0.0012  0.0000  1.0000    3005 
summary(cover$change_sr_subordspp) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0327   0.0000  10.0000     3005 
summary(cover$change_sr_rare) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0199   0.0000   7.0000     3005 


ggplot(data = mech.data, aes(x =  year, y = site_introduced_richness )) + geom_point()+ facet_wrap(~site_code) + theme_bw() +
  +     geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  +     labs(x = "site introduced species richness") +  theme_bw() +
  +     theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  +     theme(axis.text.y = element_text(size = 14)) + 
  +     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  +     theme(axis.text.x = element_text(size=14)) 
