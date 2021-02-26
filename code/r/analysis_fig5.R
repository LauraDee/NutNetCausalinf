
####################################################################################################################################
#### Models for the role of rare and non-native species  ##########################################################################################################
####################################################################################################################################

# This section is organized as follows, and accompanies section S7 in the SM. 
# Analysis for reproducing Figure 5 in the main text. 
# SM analyses for Section 7, where we do seperate analyses for:
#....... *** to fill in ****


######################################################################################################################
###  Figure 5 - Main text. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################
###########
### Figure 5 - Main text. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
#. We first present the analyses shown in the main text figure 5, which include 4 groups of species: 
# 1) rare, native: sr_nat_rare
#2) rare non-native: sr_non.nat_rare
#3) non-rare, native: sr_non.rare_nat
#4), non-rare, non-native: 
# the analyses presented in the main text Figure 5 use classify rare versus non-rare groups based 
# on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0), where non-rare are
#both subordinate and dominant species (e.g., species with DI greater than the .2 cutoff). 


MechMod_All <-feols(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) 
                    | newplotid + site.by.yeardummy, mech.data, cluster = "newplotid")

vcov_MechMod <- vcov(MechMod_All, cluster = "newplotid")

# not rare: native vs non-native 
linearHypothesis(MechMod_All, hypothesis.matrix = "ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = vcov_MechMod,  singular.ok = T)

# native rare vs non-rare
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
       coefstat = "se",
       file = "./output/Table_forFig5_R_se.tex")

esttex(MechMod_All,
       coefstat = "confint",
       file = "./output/Table_forFig5_R_ci.tex")

###################################################################################################################################
### Plot Figure 5 ######################################################################################################################
#####################################################################################################################################

Fig5.data <- tidy(MechMod_All)

# try to put all models on one line but group them
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
