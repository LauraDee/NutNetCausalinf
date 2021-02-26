////////////////////////////////////////////////////////////////////////////////
// Version Date: 02.16.2021
// Creation Date: 07.16.2020
// Paper: Dee et al. Species Richness and Productivity
//
// // This code builds on code created by C. Severen
// which is name FEandIV_20180826_control
//
////////////////////////////////////////////////////////////////////////////////


clear all
set more off 

//install psacalc and dataout
ssc install psacalc 
ssc install dataout

/*-----------------------------------------------*/
/* Global User Directory */
/*-----------------------------------------------*/

global userdir "C:\Users\pferrar5\Dropbox\Nutnet causality paper draft\Ferraro\Sensitivity"
global wmatch "$userdir/without matching"
global output "$userdir/output"

** log file for analysis
log using results_with_FE_and_IV_20180605_controlsonlySensitivity11, replace 

** Ferraro is commenting out the following three lines from Severin file to 
** show what files he used
** cd "/Users/yggdrasil/Dropbox/Research_Active/IV in ecology/NutNet/working_CS/stata_work/"
** insheet using NutNetControlPlotDataToUseApril2018.csv, c
*use nutnetcontrolplotsApril2018Data.dta

local date "20200716"

** Ferraro note: Original do file implements fixed-effects estimator with Stata
** user-developed command reghdfe for linear Models With Many Levels of 
** Fixed Effects. I keep that code but re-run with xtreg to ensure I get the
** the same result. Then run estimation with areg to run sensitivity test
** reghdfe requires the ftools package - if not installed,use next line command
** ssc install ftools
** ivreghdfe requires the ivreg2 and ranktest packages
** ssc install ivreg2
** ssc install ranktest
** for comparison, install xtivreg2
** ssc install xtivreg2
** ssc install reghdfe
reghdfe, compile

 
** use the following data set **
insheet using  "$wmatch/NutNetControlPlotDataToUseApril2018.csv", clear

/* Check characteristics of data */
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

/* Prep data and variables */
egen 	plst_id 	= group(site_code plot)
egen 	site_code2  = group(site_code)
egen 	styr_id 	= group(site_code year)
codebook site_code plot plst_id styr_id

gen 	l_lmass 	= ln(live_mass)
gen 	l_rich 		= ln(rich)
gen 	l_nbrich 	= ln(avg_neighbor_rich)
gen		l_nbrichblock	= ln(avgtrtneighrichwithinblock)
gen 	l_simpson	= ln(simpson)
gen		l_even = ln(even)

duplicates report 	plst_id year
duplicates tag 		plst_id year, generate(duplicates)
tab 	duplicates
list plst_id year year_trt rich live_mass if duplicates==1
drop if duplicates==1 & year_trt==0

xtset plst_id year


** some more data cleaning
** some soil attributes need to be destringed for use as covariates
destring pct_c pct_n percentsand percentsilt percentclay , gen(pct_c2 pct_n2 percentsand_2 percentsilt_2 percentclay_2 ) ignore("NA")

** country and habitat need to be destringed for use as covariates
encode country, gen(pais)
encode habitat, gen(tierra)

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

set matsize 1000



////////////////////////////////////////////////////////////////////////////////
//
** Bivariate and Traditional Multi-variate Approaches
//
////////////////////////////////////////////////////////////////////////////////

** Here are specifications that match as closely as we think one could do
** how traditional ecological analyses would do the analysis
** recognizing that they typically do not have panel data

** Simple bivariate correlation
** Ecologists would look for correlation between richness and biomass in a mixed 
** (or multi-level)modeling approach, which would be like a random effects
** model in economics
** They would estimate SEs clustered at the plot level and add year indicator
** A very simple model
xtreg l_lmass l_rich i.year, re cluster(plst_id)

** Now move on to multivariate approach that tries to control for 
** observable confounders
** They would condition on site management history, year and a plethora of 
** site-level weather and soil attributes; in other words, observable attributes 
** to which they had access (not necessarily the "best" observables to condition
** on. note that usually people would not add this many variables because the 
** cross-sectional sample sizes they use are too small, particularly after units 
** are dropped because of missing covariate values
** Whether country and habitat would be added is uncertain. 
** Are they potential confounders or words that describe the overall system that 
** includes the richness and biomass. But we can add them too. Does not change
** estimated coefficient by much.
xtreg l_lmass l_rich i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 , re cluster(plst_id)

** Note that this last regression includes 59 covariates and explains 57% of the 
** overall (spatial and temporal) variation in biomass. Explains 95% of between 
** plot variation; controls include 7 country variables, 11 habitat variables, 
** 11 year variables, 4 historical management variables, 7 weather variables, 
** 17 soil variables, and elevation

** Let's do the same regression specifically as a multi-level modeling with 
** random intercept exercise. need to confirm that SEs are estimated correctly 
** in this setup. but i believe this should be "analogous" to xtreg, re 
** according to my read of the Stata manual help xtmixed.  
** We specify that plots are nested within sites.

xtmixed l_lmass l_rich i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 || site_code2: || plst_id:

** same but with cluster robust SEs, which might be closer in spirit to what 
** we're doing in our main estimation, although we acknowledge that the number 
** of clusters in this reduced sample of obs is small
xtmixed l_lmass l_rich i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 || site_code2: || plst_id:, vce(robust)

** given Ferraro knows xtreg better than he knows xtmixed, we use
** results from xtreg in the main text of the paper

** Using Simpson's Index of Diversity instead of richness for multi-level models
xtreg l_lmass l_simpson i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 , re cluster(plst_id)
xtmixed l_lmass l_simpson i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 || site_code2: || plst_id:

** for multi-level models, Let's add ln evenness to model to show we're picking 
** up effect of richness on productivity not evenness 
xtreg l_lmass l_rich l_even i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 , re cluster(plst_id)
xtmixed l_lmass l_rich l_even i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 || site_code2: || plst_id:


////////////////////////////////////////////////////////////////////////////////
//
** Main Design
//
////////////////////////////////////////////////////////////////////////////////

** This is the Headline Result
reghdfe l_lmass l_rich, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** Run the same headline model but with newplotid simply to confirm that Chris' 
** plot code and Laura's plot code are the same
reghdfe l_lmass l_rich, keepsingletons  absorb(newplotid styr_id) cluster(newplotid)

** Confirm similar results when clustering SEs at the site level, 
** acknowledging that we do not have a lot of sites 
** worth reporting somewhere in manuscript
reghdfe l_lmass l_rich, keepsingletons  absorb(plst_id styr_id) cluster(site_code2)

** Using same specification as above but not keeping singletons
** Note that the seven singletons are seven years of observation of a single 
** plot at the pape.de site, so when we specify site-by-year FE in the reghdfe 
** command this plot becomes a singleton because no other plots from site are in
** the data set
reghdfe l_lmass l_rich, absorb(plst_id styr_id) cluster(plst_id)
** Dropping 7 singletons decreases SE at fourth decimal place and thus 
** for ease of exposition and to be conservative, we keep singletons
** Note that even if we cluster SEs at the site level, p=0.01 and CIs get only slightly bigger
reghdfe l_lmass l_rich, absorb(plst_id styr_id) cluster( site_code2 )

** confirm that we get the same answer from xtreg for outcome as we got with 
** reghdfe
xtreg l_lmass l_rich i.styr_id, fe cluster(plst_id)

** for comparison with the random effects estimator with covariates above we 
** re-run FE estimation but only for subsample that has values for covariates
** the estimated effect is still negative, indicating that the positive value 
** we obtained above is not a function of difference in sample composition after
** dropping obs, but rather a function of the identification strategy
reghdfe l_lmass l_rich if elevation~=. & managed~=. & burned~=. & grazed~=. &  anthropogenic~=. &  temp_var_v2~=. &   min_temp_v2~=. &  max_temp_v2~=. &  temp_wet_q_v2~=. &  temp_dry_q_v2~=. &  temp_warm_q_v2~=. &  temp_cold_q_v2~=. &  pct_c2~=. &  pct_n2~=. &  ppm_p~=. &  ppm_k~=. &  ppm_ca~=. &  ppm_mg~=. &  ppm_s~=. &  ppm_na~=. &  ppm_zn~=. &  ppm_mn~=. &  ppm_fe~=. &  ppm_cu~=. &  ppm_b~=. &  ph~=. & percentsand_2~=. &  percentsilt_2~=. &  percentclay_2~=., keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** Using Simpson's Index of Diversity instead of richness
reghdfe l_lmass l_simpson, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
reghdfe l_lmass l_simpson, absorb(plst_id styr_id) cluster(plst_id)

** Let's add ln evenness to model to show we're picking up effect of richness 
** on productivity not evenness
reghdfe l_lmass l_rich l_even, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** one way to see that the plot-level FEs have little effect is to just do 
** fixed effects estimation without site-by-year effects
** the site-by-year effects are the driver of the sign switching
reghdfe l_lmass l_rich, keepsingletons  absorb(plst_id) cluster(plst_id)
xtreg l_lmass l_rich, fe cluster(plst_id)

** that last result suggests that plot-level confounding is not the culprit
** it is the time-by-site fixed effects that is causing the issue
** they have a positive effect on richness and a positive effect on biomass

** There are 2 observations without richness values and 58 observations without 
** biomass. the former is too few to matter. the latter is relatively few 
** but potentially concerning if correlated with richness
** so let's create some variables to look more closely at attrition
gen attrition=1 if live_mass==.
recode attrition (.=0)
ttest rich, by(attrition)

** although the difference is not huge, the attriting plots have higher richness
** so let's stress test the results by imputing biomass for the 58 plots
** let's let their values be equal to the 95th percentile of biomass values

sum     live_mass, detail
** 95th percentile is 827.5
gen live_massimpute= live_mass
recode live_massimpute (.=827.5)
gen l_lmassimpute= ln(live_massimpute)

** confirm similar results to Headline Result when addressing attrition 
** through bounded imputation
reghdfe l_lmassimpute l_rich, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)


////////////////////////////////////////////////////////////////////////////////
//
** Lagged Dependent Variable Analysis
//
////////////////////////////////////////////////////////////////////////////////

** a dynamic model which also includes site-by-effects
** a.k.a. the lagged dependent variable model
reghdfe l_lmass l_rich L.l_lmass, keepsingletons absorb(styr_id) cluster(plst_id)
reghdfe l_lmass l_rich L.l_lmass, a(styr_id) cluster(plst_id)
** use areg just to confirm that the reghdfe output makes sense
areg l_lmass l_rich L.l_lmass , absorb(styr_id) cluster(plst_id)


////////////////////////////////////////////////////////////////////////////////
//
** "Control" for Shading Mechanism to assess reverse causality
//
////////////////////////////////////////////////////////////////////////////////

** shut down shade mechanism to see if it affects estimated effect of 
** richness on productivity
** not sure if I should put ground_par in level or log form, so do both
gen l_ground_par=ln(ground_par)
reghdfe l_lmass l_rich l_even ground_par , keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
reghdfe l_lmass l_rich l_even l_ground_par , keepsingletons  absorb(plst_id styr_id) cluster(plst_id)



////////////////////////////////////////////////////////////////////////////////
//
** IV analysis
//
////////////////////////////////////////////////////////////////////////////////

** the best IV would seem to be the average neighbor richness in the same block 
** - i.e., the interference is assumed to be local
** even better would the the average nitrogen enriched plot richness in the same
** block - i.e., interference is local and coming from plots for whom richness 
** is being exogenously manipulated through nitrogen additions
** we must assume that the nitrogen addition to a neighbor's plot in the block  
** has no effect on own productivity except through its effect on own richness

** Ferraro cannot get Chris' code to run
** His code: reghdfe l_lmass 	(l_rich = l_nbrich), a(plst_id styr_id) cluster(plst_id) ff
** It appears that the IV aspect of this command has migrated to ivreghdfe and 
** so I install that command and run

ivreghdfe l_lmass (l_rich = l_nbrichblock),a(plst_id styr_id) cluster(plst_id) first

** keepsingletons option is not allowed with ivreghdfe
** Note that we're assuming excludability of the IV conditional on site-by-year effects

** Laura and Chris also used the richness of the entire site, 
** but is this appropriate?
** ivreghdfe l_lmass (l_rich = l_nbrich), a(plst_id styr_id) cluster(plst_id)

** Running xtivreg instead. Should not change estimated effect. 
** estimated SE differences should be trivial.
** Note sample size changes a bit from ivreghdfe b'c "keepsingletons" 
** is not allowed
xtivreg l_lmass i.styr_id (l_rich = l_nbrichblock),fe vce(cluster plst_id) first


////////////////////////////////////////////////////////////////////////////////
//
//  Oster Sensitivity Analysis
//
////////////////////////////////////////////////////////////////////////////////

** Prep work for the sensitivity analysis that will use areg for two regressions
** one for the treatment "selection" equation (richness as dependent variable)
** one for the outcome equation (live mass as the dependent variable)

** areg on outcome equation log-log specification **
areg l_lmass l_rich i.styr_id, absorb(plst_id) cluster(plst_id)
estimates store m1, title(Areg1)

** the following specification does the same thing as the one above: 
** areg l_lmass l_rich i.year#i.sitenumber, absorb(plst_id) cluster(plst_id)

** areg on selection equation **
areg l_rich i.styr_id, absorb(plst_id) cluster(plst_id)
estimates store m2, title(Areg2)

esttab * using "$output/TableSens_a.tex", cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label collabels(, none) varlabels(_cons Constant) posthead("") prefoot("") postfoot("") varwidth(16) modelwidth(12) title(Sensitivity Test: Hidden Bias) replace

** confirm that you get the same answer from xtreg for selection as you get 
** below with areg
xtreg l_rich i.styr_id, fe cluster(plst_id)

** linear regression should give you same answer as areg
reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)

***************** post-estimation sensitivity analysis: psacalc ****************

** I cannot get the bootstrapped SEs using bs r(beta), rep(100): prior to
** launching psacalc not sure why. receive error message that there are not
** enough observations. so going to just estimate the beta
** The R2 in the outcome equation above is 0.8675 and the R2 in the selection
** equation is 0.9059
** so the unobserved confounder can only be about 1/10 as strong (0.0941) as 
** the observables in explaining variation in richness 
** so delta, the measure of proportional selection, is roughly -0.10. 
** In our published SM, delta in Oster's set-up is designated pi to avoid 
** confusion with the delta that we use to designate plot level fixed effects

psacalc beta l_rich, model(areg l_lmass l_rich i.styr_id, absorb(plst_id) cluster(plst_id)) delta(-0.1) rmax(1) 

** could also write in the following way
psacalc beta l_rich, model(reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)) delta(-0.1) rmax(1)

** for completeness let's run the sensitivity test assuming positive value of delta
psacalc beta l_rich, model(reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)) delta(0.1) rmax(1)

** for completeness let's run the sensitivity test using the lagged dependent variable model

psacalc beta l_rich, model(areg l_lmass l_rich L.l_lmass , absorb(styr_id) cluster(plst_id)) delta(-0.1) rmax(1)

** the site-year fixed effects keep me from being able to bootstrap the SEs on beta. I believe it's because there are too few observations for boostrapping with this many variables in the regression
** bs r(beta), rep(100): psacalc beta l_rich, model(areg.....) delta(-0.1) rmax(1) 


////////////////////////////////////////////////
//
//  Dominant Rare Native Non-native species
/////////////////////////////////////

** Ferraro does not have the breakdown by species type and cannot run this 
** analysis in Stata


////////////////////////////////////////////////
//
//  SM analyses
/////////////////////////////////////

** for the SM
** let's try the ADL AR(1) model
reghdfe l_lmass l_rich L.l_rich L.l_lmass L2.l_lmass, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
** use areg just to confirm that the reghdfe output makes sense
areg l_lmass l_rich L.l_rich L.l_lmass L2.l_lmass , absorb(styr_id) cluster(plst_id)

** let's try a model that blocks effect of lagged richness on current richness
reghdfe l_lmass l_rich L.l_rich, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)


** For the SM
** in log-linear and linear forms
reghdfe l_lmass rich, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
reghdfe live_mass rich, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** one argument for doing ln-ln model is non-linearity. we could skip log 
** transform and add higher order terms
gen rich2 = rich*rich
reghdfe live_mass rich rich2, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
** productivity hits min at richness = 20, which is about the 92nd percentile, 
** implying the negative relationship holds for most of the support
sum rich

** cubic 
gen rich3 = rich*rich*rich
reghdfe live_mass rich rich2 rich3, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)




**************
log close

translate results_with_FE_and_IV_20180605_controlsonlySensitivity11.smcl results_with_FE_and_IV_20180605_controlsonlySensitivity11.log, replace
