library(ggpubr)

#############################################
## Rank abundance curves - Figure S11 ####
############################################

sites <- unique(comb$site_code)
sites <- sites[-29]

rac.dat <- cover[which(cover$year_trt == 0 & cover$site_code %in% sites),]

rac.dat <- unique(rac.dat[,c("site_code","Taxon","relative_abundance_spp_site.yr0")])
# create lists for each site

newdf <- data.frame()
for (i in 1:length(sites)){
  temp <- rac.dat[which(rac.dat$site_code == sites[i]),]
  temp <- temp[order(temp$relative_abundance_spp_site.yr0, decreasing = TRUE)]
  temp$rank <- seq(1:nrow(temp))
  newdf <- rbind(temp, newdf)
  }
 
# Option 1
ggplot(aes(x = rank, y = relative_abundance_spp_site.yr0), data = newdf) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~site_code) + 
  theme_classic()

#option 2
dat.list <- list()
for (i in 1:length(sites)){
  dat.list[[i]] <- newdf[which(newdf$site_code == sites[i]),]
}

grfun = function(df){
  site <- unique(df$site_code)
  ggplot(aes(x = rank, y = relative_abundance_spp_site.yr0), data = df) + 
    geom_point() + 
    geom_line() +
    labs(x = "Species Rank", y = "Species Abundance", title  = site) +
    theme_classic()
}

gr.list <- lapply(dat.list, grfun)

ggarrange(plotlist = gr.list)
