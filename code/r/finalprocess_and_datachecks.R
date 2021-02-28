# Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-09-Apr-2018.csv'
length(comb$live_mass) #74 NAs for live mass
summary(comb$live_mass)
length(comb$rich)
summary(comb$rich) # 2 NAs for richness

# ---- process ----
comb$site <- comb$site_code
comb = comb[!is.na(live_mass) & !is.na(rich)]

# ---- ignore ----
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

# ---- classcheck ----
class(comb$year)
class(comb$site_code) 
class(comb$newplotid) 

# ---- ignore ----
## 
length(unique(comb$newplotid)) #172

# ---- removesingleton ----
comb[,.N, by=c("newplotid")][N==1,]  # there are 21 that need to be removed.

# flag the singleton and remove the singletons
comb[,singleton:=(.N==1), by=c("newplotid")]
comb = comb[singleton == F, ]

# ---- ignore ----

# Determine the number of Obs. removing the singletons and obs with NA or rich or live_mass
nrow(comb) #1231 
length(unique(comb$newplotid)) #151

# ---- morefactors ----
comb[, site.by.yeardummy := paste(site_code, year, sep = "_")]
comb$year <- as.character(comb$year)
 
# ---- tableS1 ----
  comb.descript.v1 =  table(comb$site_name, comb$year)
  print(comb.descript.v1)
  
  write.csv(comb.descript.v1, "./output/Table_S1.csv")
  