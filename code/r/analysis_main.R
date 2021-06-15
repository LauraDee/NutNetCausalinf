##############################################
##############################################
## Run Models #
##############################################
##############################################

##############################################
## Models for Figure 2A and Table S2 


MainMod_Rich     <- feols(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy, comb) 
MainMod_Rich_sitecluster  <- feols(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy, comb, cluster = "site_code")
MainMod_RichEven <- feols(log(live_mass) ~ log(rich) + ihs(even) | newplotid + site.by.yeardummy, comb) 
MainMod_Simpson  <- feols(log(live_mass) ~ log(simpson) | newplotid + site.by.yeardummy, comb) 
MainMod_RichLag  <- feols(log(live_mass) ~ log(rich) + log(laggedrich) | newplotid + site.by.yeardummy, comb) 
MainMod_RichEvenLag <- feols(log(live_mass) ~ log(rich) + log(laggedrich) + ihs(even) | newplotid + site.by.yeardummy, comb)

Fig2A.1 <- tidy(MainMod_Rich) %>%
  filter(term == "log(rich)")
Fig2A.2 <- tidy(MainMod_RichEven) %>%
  filter(term == "log(rich)")
Fig2A.3 <- tidy(MainMod_Simpson) %>%
  filter(term == "log(simpson)")

##############################################
## Models for Figure 2B and Table S2, and Bivariate and Common Multivariate Analyses

MixedMod_Rich <- lmer(log(live_mass) ~ log(rich) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                        elevation + managed + burned + grazed + anthropogenic + 
                        TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                        pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                        pH + PercentSand + PercentSilt + PercentClay + 
                        (1|newplotid) + (1|site_code), comb, REML = F)

MixedMod_RichEven <- lmer(log(live_mass) ~ log(rich) + ihs(even) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                            elevation + managed + burned + grazed + anthropogenic + 
                            TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                            pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                            pH + PercentSand + PercentSilt + PercentClay + 
                            (1|newplotid) + (1|site_code), comb, REML = F)

MixedMod_Simpson <- lmer(log(live_mass) ~ log(simpson) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                           elevation + managed + burned + grazed + anthropogenic + 
                           TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                           pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                           pH + PercentSand + PercentSilt + PercentClay + 
                           (1|newplotid) + (1|site_code), comb, REML = F)

#  Bivariate and Common Multivariate Analyses Models

SimpleRE_TwoVar <- lmer(log(live_mass) ~ log(rich) + as.factor(year) + (1|newplotid), comb, REML = F)

SimpleRE_MultiVar <- lmer(log(live_mass) ~ log(rich) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                        elevation + managed + burned + grazed + anthropogenic + 
                        TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                        pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                        pH + PercentSand + PercentSilt + PercentClay + 
                        (1|newplotid), comb, REML = F)


## CS: Because I couldn't find heteroskedasticity and autocorrelation robust variances with lme4, let's
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


##############################################
## Models for Figure 3

MainMod_Rich #=Done above
MainMod_LagLiveMass  <- feols(log(live_mass) ~ log(rich) + log(laggedlive_mass) | site.by.yeardummy, comb, cluster = "newplotid")
MainMod_Oster <- "Done in Stata" 
MainMod_MechBlocking <- feols(log(live_mass) ~ log(rich) + ihs(even) + log(proportion_par) | newplotid + site.by.yeardummy, comb)


MainMod_IVRevCaus <- feols(log(live_mass) ~ 1 | newplotid + site.by.yeardummy | log(rich) ~ log(avg.trt.neigh.rich.within.block), comb)

Fig3.1 <- tidy(MainMod_Rich) %>%
  filter(term == "log(rich)")
Fig3.2 <- tidy(MainMod_LagLiveMass) %>%
  filter(term == "log(rich)")
Fig3.3 <- tibble(term="log(rich)",
                 estimate= -0.19962, 
                 std.error = 0,
                 conf.low = -0.19962,
                 conf.high= -0.19962)
Fig3.4 <- tidy(MainMod_MechBlocking) %>%
  filter(term == "log(rich)")
Fig3.5 <- tidy(MainMod_IVRevCaus) %>%
  filter(term == "fit_log(rich)")

##############################################
## Models for Table S3 (functional form)

MainMod_LogLog  <- MainMod_Rich
MainMod_LogLev  <- feols(log(live_mass) ~ rich  | newplotid + site.by.yeardummy, comb) 
MainMod_LevLev  <- feols(live_mass ~ rich | newplotid + site.by.yeardummy, comb) 
MainMod_LevQuad <- feols(live_mass ~ rich + I(rich^2) | newplotid + site.by.yeardummy, comb) 


##############################################
##############################################
## Tables #
##############################################
##############################################

################################################
## Table S2
#######################################

esttex(MainMod_Rich, MainMod_Rich_sitecluster, MainMod_RichEven, MainMod_RichLag, MainMod_RichEvenLag, MainMod_Simpson,
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S2_R_se.tex")

esttex(MainMod_Rich, MainMod_Rich_sitecluster,MainMod_RichEven, MainMod_RichLag, MainMod_RichEvenLag, MainMod_Simpson,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S2_R_ci.tex")

################################################
## Table S3
#######################################

esttex(MainMod_LogLog, MainMod_LogLev, MainMod_LevLev, MainMod_LevQuad, 
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S3_R_se.tex")

esttex(MainMod_LogLog, MainMod_LogLev, MainMod_LevLev, MainMod_LevQuad, 
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_S3_R_ci.tex")

################################################
## Table for Fig 3
#######################################

esttex(MainMod_Rich, MainMod_LagLiveMass, MainMod_Oster, MainMod_MechBlocking, MainMod_IVRevCaus, 
       coefstat = "se", replace = TRUE,
       file = "./output/Table_forFig3_R_se.tex")

esttex(MainMod_Rich, MainMod_LagLiveMass, MainMod_Oster, MainMod_MechBlocking, MainMod_IVRevCaus,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_forFig3_R_ci.tex")

################################################
## Table S9: Dynamic Panel Results 
#######################################

esttex( MainMod_LagLiveMass, 
       coefstat = "se", replace = TRUE,
       file = "./output/Table_S9_R_se.tex")

esttex(MainMod_LagLiveMass,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_Table_S9_R_ci.tex")

# Print results for mechanism blocking analysis, reported in the SM:
esttex(MainMod_MechBlocking,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_MechBlocking_R_ci.tex")

esttex(MainMod_MechBlocking,
       coefstat = "confint", replace = TRUE,
       file = "./output/Table_MechBlocking_R_ci.tex")

##############################################
##############################################
## Figures #
##############################################
##############################################

## Master Palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1", "red")

################################################
## Plot Figure 2A
################################################

# Prep Data
Fig2A.data <-  bind_rows(
  Fig2A.1 %>% mutate(reg = "Species Richness"),
  Fig2A.2 %>% mutate(reg = "Species Richness controlling for Evenness"))

Fig2A.data$term = factor(Fig2A.data$term,
                         levels=c("log(rich)",
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
  scale_x_discrete(labels = c("Species Richness" , "Simpson's Diversity")) + 
  scale_y_continuous(limits=c(-.8, .85), 
                     breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8)
  ) %>%
  labs(title = "Main Study Design",
       caption = "", x = "Variable", y = "Estimated effect size"
  )

Fig2A.plot
ggsave("./output/Fig2A.pdf", Fig2A.plot)

################################################
## Plot Figure 2B
################################################
Fig2B.data <-  bind_rows(
  Fig2B.1 %>% mutate(reg = "Species Richness"),
  Fig2B.2 %>% mutate(reg = "Species Richness controlling for Evenness"))

Fig2B.data$term = factor(Fig2B.data$term,
                         levels=c("log(rich)", 
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
  labs(title = "Common Design in Ecology",
       caption = "", x = "Variable", y = "Estimated effect size"
  )

Fig2B.plot
ggsave("./output/Fig2B.pdf", Fig2B.plot)


#################################
######## Combine 2A and 2B

common.ylab = ylab("Estimated effect size")  #Estimated % Change in Productivity from a 1% Change in Diversity
Fig2.both <- plot_grid(Fig2A.plot  + common.ylab,
                       Fig2B.plot + common.ylab)
Fig2.both
ggsave("./output/Fig2.pdf", Fig2.both, width=13, height=6)

#################################
######## Figure 3

Fig3.data <-  bind_rows(
  Fig3.1 %>% mutate(reg = "1. Main Study Design"),
  Fig3.2 %>% mutate(reg = "2. Dynamic Panel Design"),
  Fig3.3 %>% mutate(reg = "3. Sensitivity Test"),
  Fig3.4 %>% mutate(reg = "4. Mechanism Blocking Design"),
  Fig3.5 %>% mutate(reg = "5. Instrumental Variable Design") )

# Plot                    
Fig3.plot <- Fig3.data %>%
  ggplot(aes(x=reg, y=estimate, ymin = conf.low, ymax = conf.high, colour = reg)) +
  geom_pointrange(size = 1.5) +
  scale_colour_discrete(name="Model") +
  scale_color_manual(values=cbPalette[c(7,6,2, 3, 1)]) +
  theme_classic() +
  theme(legend.position = c(0.4, 0.8),
        legend.title = element_blank(), 
        legend.text  = element_text(size=16),
        legend.background = element_rect(size=0.5, 
                                         linetype="solid",
                                         colour ="black" ),
        axis.text=element_text(size=22),
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5) ) + 
  geom_hline(yintercept = 0, col = "black") +
  scale_y_continuous(limits=c(-1.1, 1.1)) +
  ylim(-1.1, 1.1) %>%
  labs(title = element_blank(),
       caption = "", x = "Design", y = "Estimated effect size"
  )

Fig3.plot
ggsave("./output/Fig3.pdf", Fig3.plot)


##################################################################
# Figure S4: Including Simpson's Diversity Estimate 
##################################################################
# for Main Design Figure S4A:
# Prep Data
FigS4.data <-  bind_rows(
  Fig2A.1 %>% mutate(reg = "Species Richness"),
  Fig2A.2 %>% mutate(reg = "Species Richness controlling for Evenness"),
  Fig2A.3 %>% mutate(reg = "Simpson's Diversity") )

FigS4.data$term = factor(FigS4.data$term,
                         levels=c("log(rich)", 
                                  "log(simpson)",
                                  "ihs(even)"))
# Plot
FigS4.plot <- FigS4.data %>%
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
  scale_x_discrete(labels = c("Species Richness" , "Simpson's Diversity")) + 
  scale_y_continuous(limits=c(-.8, .85), 
                     breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8)
  ) %>%
  labs(title = "Main Study Design",
       caption = "", x = "", y = "Estimated effect size"
  )

FigS4.plot
ggsave("./output/FigS4.pdf", FigS4.plot)

# plot traditional ecological model with simpsons for the SM Figure S4B
#prep model output data
FigS4B.data <-  bind_rows(
  Fig2B.1 %>% mutate(reg = "Species Richness"),
  Fig2B.2 %>% mutate(reg = "Species Richness controlling for Evenness"),
  Fig2B.3 %>% mutate(reg = "Simpson's Diversity") )

FigS4B.data$term = factor(FigS4B.data$term,
                         levels=c("log(rich)", 
                                  "log(simpson)",
                                  "ihs(even)"))
# Plot
FigS4B.plot <- FigS4B.data %>%
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
  labs(title = "Common Design in Ecology",
       caption = "", x = "", y = "Estimated effect size"
  )

FigS4B.plot
ggsave("./output/FigS4B.pdf", FigS4B.plot)

#################################
######## Combine Figure S4 A and B
common.ylab = ylab("Estimated effect size") #Estimated % Change in Productivity from a 1% Change in Diversity
common.xlab = xlab("") 
FigS4.both <- plot_grid(FigS4.plot  + common.ylab,
                       FigS4B.plot + common.ylab)
FigS4.both
ggsave("./output/FigS4.pdf", Fig2.both, width=13, height=6)