
####################################################################################################################################
#### Models for the role of rare and non-native species  ##########################################################################################################
####################################################################################################################################

# This section is organized as follows, and accompanies section S7 in the SM. 
# Analysis for reproducing Figure 5 in the main text. 
# SM analyses for Section 7, where we do separate analyses for:
#....... *** to fill in ****


######################################################################################################################
###  Figure 5 - Main text. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################
###########
### Figure 5 - Main text. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
#. We first present the analyses shown in the main text figure 5, which include 4 groups of species, which are: 
# 1) rare, native: sr_nat_rare
#2) rare non-native: sr_non.nat_rare
#3) non-rare, native: sr_non.rare_nat
#4) non-rare, non-native: 
# the analyses presented in the main text Figure 5 use classify rare versus non-rare groups based 
# on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0), where non-rare are
#both subordinate and dominant species (e.g., species with DI greater than the .2 cutoff). 


MechMod_All <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) 
                    | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")

vcov_MechMod <- vcov(MechMod_All, cluster = "newplotid")

# Hypotheses Tests for the groups: are their effects on productivity significantly differently?
#Note: fevcov returns a square matrix with the bias corrected covariances. An attribute 'bias' contains
# the biases. The bias corrections have been subtracted from the bias estimates. I.e. vc = vc’ - b,
# where vc’ is the biased variance and b is the bias.

#testing if the groups are not all the same: rejecting the null that they are the same
linearHypothesis(MechMod_All, 
                 hypothesis.matrix = c("ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)",
                                       "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)


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


##################
## Export Table

esttex(MechMod_All, 
       coefstat = "se", replace = TRUE,
       file = "./output/Table_forFig5_R_se.tex")

esttex(MechMod_All,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_forFig5_R_ci.tex")

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
       caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size"
  )

Fig5.plot
ggsave("./output/Fig5.pdf", Fig5.plot)

###################################################################################################################################
### SM Analyses for Section S7  ######################################################################################################################
#####################################################################################################################################

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
                  +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")
vcov_MechFreq1 <- vcov(MechFreq1, cluster = "newplotid")

# hypothesis tests


#####################################################################################################################
#### Run Models Comparing Rare vs Non-rare, Native or Invasive with different grouping cutoffs of breaks=c(0.0,0.4,0.8,1.0).
#####################################################################################################################

#####################
### Table SX.  Grouped based on the Dominance Indicator (DI) and cutoffs of: breaks=c(0.0,0.4,0.8,1.0), versus cut-offs from Figure 5
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
                                       "ihs(sr_non.nat_rare2) = ihs(sr_non.rare_non.nat2)"), # by transitivity this is included but needs to be dropped: "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)"),  
                 test = "F", vcov = vcov_MechMod2,  singular.ok = T)





#########################################################################################################################################################################
## compare results with other ways of defining rare vs non-rare, based only on relative abundance and relative frequency for both cut-offs ##################################
#########################################################################################################################################################################

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
sr.metrics <- mech.data[, .(sr_non.nat_rare, sr_nat_rare, non_rare_spp,
                            sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom,
                            sr_non.nat_dom, sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare )]
cor(sr.metrics)
corrplot(cor(sr.metrics), method = "square", tl.cex = .5)

