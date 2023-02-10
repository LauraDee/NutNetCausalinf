insheet using  "$datadir/NutNetControlPlotData_derived.csv", clear

** Check characteristics of data 
destring total_mass live_mass simpson shan rich ground_par ppm_* ph avg_neighbor_rich avgtrtneighrichwithinblock lagged_avg_neighbor_rich even initial_site_rich, i("NA") replace
count 
sum 	live_mass rich
tab 	year
rename  year year_actual
gen 	year = year_actual - 2006
bys plot site: egen yr_in_plot = count(total_mass)
tab 	yr_in_plot
drop if yr_in_plot<5
tab 	yr_in_plot

** Prep data and variables 
egen 	plst_id 	= group(site_code plot)
egen 	site_code2  = group(site_code)
egen 	styr_id 	= group(site_code year)
codebook site_code plot plst_id styr_id

** Remove duplicate rows
duplicates report 	plst_id year
duplicates tag 		plst_id year, generate(duplicates)
tab 	duplicates
list plst_id year year_trt rich live_mass if duplicates==1
drop if duplicates==1 & year_trt==0

** Set data as panel data
xtset plst_id year

** Create analysis variables
gen 	l_lmass 	= ln(live_mass)
gen 	l_rich 		= ln(rich)
gen 	l_nbrich 	= ln(avg_neighbor_rich)
gen		l_nbrichblock	= ln(avgtrtneighrichwithinblock)
gen 	l_simpson	= ln(simpson)
gen		l_even 		= ln(even)
gen 	l_ground_par=ln(ground_par)
gen		ihs_even 	= ln(even + sqrt(1+even^2))
gen 	rich2		= rich^2

** some soil attributes need to be destringed for use as covariates
destring pct_c pct_n percentsand percentsilt percentclay , gen(pct_c2 pct_n2 percentsand_2 percentsilt_2 percentclay_2 ) ignore("NA")

** country and habitat need to be destringed for use as covariates
encode country, gen(pais)
encode habitat, gen(tierra)

** SAVE OUTPUT

save 	"$datadir/NutNet_Prepped.dta", replace
outsheet using "$datadir/NutNet_Prepped.csv", comma replace
clear all

** PROCESSING NOTES FROM PAUL AND CHRIS

/* We started with 1340 observations. 2 don’t have richness values. 
** 74 don’t have biomass values.  
** 37 observations were in the plot for fewer than 5 years
** (specifically, all have only 1 year in plot). 
** 16 of these 37 obs have no biomass (all have richness). 
** We also have 3 plots that have duplicate year entries, creating 24 obs
** where there should only be 12. We drop the duplicates – drop 12 
** (none of these duplicates are part of the 37 with fewer than 5 yrs). 
** 1,340 -> Drop 12 duplicates -> 1,328 -> Drop 37 with <5yr -> 1,291 -> 
** Two don’t have richness -> 1,289 - > 58 don’t have biomass values 
** (the other 16 without biomass values were dropped in earlier steps) -> 
** 1,231 observations in the main analysis. */

** both live_mass and rich are right-skewed, thus a ln transform of both is
** preferred; yield more precise estimates & reduce the influence of outliers