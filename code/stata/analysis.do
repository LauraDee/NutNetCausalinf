use 	"$datadir/NutNet_Prepped.dta"

////////////////////////////////////////////////////////////////////////////////
//
** Figure 2A and Table S2 (Main Design)
//
////////////////////////////////////////////////////////////////////////////////


** Headline results: Figure 2A and Table S2
** Notes: We keep singletons below slightly increases standard errors (at the 
** 		fourth decimal place). We keep singeltons to be conservative and ease 
** 		exposition. There are singletons, which correspond to seven years of 
**		data for a single plot at the pape.de site; when we include site-by-year 
** 		FEs this plot becomes a singleton because no other plots from site are 
**  	in the data.
** Including l_even controls for log evenness, l_rich barely changes
** l_simpson is Simpson's diversity index


** Figure 2A shows the following:
reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
reghdfe l_lmass l_rich ihs_even, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** Table S2 (Columns 1-4) includes the following 

est clear

eststo: reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
eststo: reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(site_code2)
eststo: reghdfe l_lmass l_rich ihs_even, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
eststo: reghdfe l_lmass l_simpson, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
eststo: reghdfe l_lmass l_rich l_laggedrich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
eststo: reghdfe l_lmass l_rich ihs_even l_laggedrich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

esttab * using "$userdir/output/Table_S2.tex", fragment cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label replace
est clear

// ADDITIONAL SUPPORT //

** Supporting: Drop singletons
** Inference is not affected by dropping singletons
reghdfe l_lmass l_rich, absorb(plst_id styr_id) cluster(plst_id)

** Supporting: Cluster at site
** Confirms that inference is robust to clustering SEs at the site level, 
** 		acknowledging that we do not have a lot of sites. Worth reporting.
reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(site_code2)


////////////////////////////////////////////////////////////////////////////////
//
** Figure 2B: Traditional Multi-variate Approaches
//
////////////////////////////////////////////////////////////////////////////////

local control_vars i.pais i.tierra i.year elevation managed burned grazed anthropogenic ///
					temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 ///
					temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ///
					ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 

xtmixed l_lmass l_rich `control_vars' || site_code2: || plst_id:, vce(robust)
xtmixed l_lmass l_rich ihs_even `control_vars' || site_code2: || plst_id:, vce(robust)

// ADDITIONAL SUPPORT //

** NOTE: to compare with the FEs and random effects methods we 
**   re-run FE estimation but only for subsample that has values for covariates
**   the estimated effect is still negative, indicating that the positive value 
**   we obtained above is not a function of difference in sample composition after
**   dropping obs, but rather a function of the identification strategy
reghdfe l_lmass l_rich if elevation~=. & managed~=. & burned~=. & grazed~=. & ///
							anthropogenic~=. &  temp_var_v2~=. &   min_temp_v2~=. & ///
							max_temp_v2~=. &  temp_wet_q_v2~=. &  temp_dry_q_v2~=. & ///
							temp_warm_q_v2~=. &  temp_cold_q_v2~=. &  pct_c2~=. &  ///
							pct_n2~=. &  ppm_p~=. &  ppm_k~=. &  ppm_ca~=. & ///
							ppm_mg~=. &  ppm_s~=. &  ppm_na~=. &  ppm_zn~=. & ///
							ppm_mn~=. &  ppm_fe~=. &  ppm_cu~=. &  ppm_b~=. & ///
							ph~=. & percentsand_2~=. &  percentsilt_2~=. &  percentclay_2~=., ///
						keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** "Simple" bivariate analysis done with random effects
						
** From Paul F:
/* Ecologists would look for correlation between richness and biomass in a mixed 
	(or multi-level)modeling approach, which would be like a random effects
	model in economics
	They would estimate SEs clustered at the plot level and add year indicator
*/
** A very simple model
xtreg l_lmass l_rich i.year, mle 

** From Paul F:
/* Now move on to multivariate approach that tries to control for observable confounders
	They would condition on site management history, year and a plethora of 
	site-level weather and soil attributes; in other words, observable attributes 
	to which they had access (not necessarily the "best" observables to condition
	on. note that usually people would not add this many variables because the 
	cross-sectional sample sizes they use are too small, particularly after units 
	are dropped because of missing covariate values
	Whether country and habitat would be added is uncertain. 
	Are they potential confounders or words that describe the overall system that 
	includes the richness and biomass. But we can add them too. Does not change
	estimated coefficient by much.
*/
xtreg l_lmass l_rich i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2, mle

xtreg l_lmass l_rich i.pais i.tierra i.year elevation managed burned grazed anthropogenic temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2, re

/* Note that this last regression includes 59 covariates and explains 57% of the 
	overall (spatial and temporal) variation in biomass. Explains 95% of between 
	plot variation; controls include 7 country variables, 11 habitat variables, 
	11 year variables, 4 historical management variables, 7 weather variables, 
	17 soil variables, and elevation	
*/
						
////////////////////////////////////////////////////////////////////////////////
//
** Figure 3: Robustness across various designs
//
////////////////////////////////////////////////////////////////////////////////


** Main Study Design
reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Dynamic Panel Design Lagged Dependent Variable
reghdfe l_lmass l_rich l_laggedlive_mass, keepsingletons absorb(styr_id) cluster(plst_id)

** Sensitivity Test (note: no standard errors)
psacalc beta l_rich, model(reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)) delta(-0.1) rmax(1)

** Mechanism Blocking Design
reghdfe l_lmass l_rich ihs_even l_ground_par , keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

** Instrumental Variable Design
ivreghdfe l_lmass (l_rich = l_nbrichblock), a(plst_id styr_id) cluster(plst_id) first

// ADDITIONAL SUPPORT //

** Sensitivity Test: How low delta makes effect 0? About -1.69:
psacalc beta l_rich, model(reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)) delta(-1.69) rmax(1)


////////////////////////////////////////////////////////////////////////////////
//
** Table S3 (Functional Form Tests)
//
////////////////////////////////////////////////////////////////////////////////

est clear

** Log-Log
eststo: reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Log-Level
eststo: reghdfe l_lmass rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Level-Level
eststo: reghdfe live_mass rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Level-Quadratic
eststo: reghdfe live_mass rich rich2, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

esttab * using "$userdir/output/Table_S3.tex", fragment cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label replace
est clear


////////////////////////////////////////////////////////////////////////////////
//
** Table S7 (IV)
//
////////////////////////////////////////////////////////////////////////////////

est clear
eststo: ivreghdfe l_lmass (l_rich = l_nbrichblock), a(plst_id styr_id) cluster(plst_id) first

gen 	includeinfirst = e(sample)
eststo: reghdfe l_rich l_nbrichblock if includeinfirst==1, a(plst_id styr_id) cluster(plst_id) 

esttab * using "$userdir/output/Table_S7.tex", fragment cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label replace
est clear
drop includeinfirst



////////////////////////////////////////////////////////////////////////////////
//
** Table S8: Lagged Dependent Varaibles (also replicates a result in Table 3)
//
////////////////////////////////////////////////////////////////////////////////

est clear

eststo: reghdfe l_lmass l_rich l_laggedlive_mass, keepsingletons a(styr_id) cluster(plst_id)

esttab * using "$userdir/output/Table_S8.tex", fragment cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label replace
est clear

// ADDITIONAL SUPPORT //
/* These are additional models to test dynamics; all return a negative 
	coefficient on l_rich */

/* SM Equation S3 ADL AR(1) model in the SM */
reghdfe l_lmass l_rich l_laggedrich l_laggedlive_mass L.l_laggedlive_mass, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

/* SM Equation S3 with more lags in the SM */
reghdfe l_lmass l_rich l_laggedrich L.l_laggedrich l_laggedlive_mass L.l_laggedlive_mass L2.l_laggedlive_mass, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

/* Hybrid Analysis in SM for Reviewer 5 that combines Lagged Dependent Variable Design with Main Design */
reghdfe l_lmass l_rich l_laggedlive_mass, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)

////////////////////////////////////////////////////////////////////////////////
//
** Table S9 (Sensitivity Analysis a la Oster)
//
////////////////////////////////////////////////////////////////////////////////

** areg on outcome equation log-log specification **
est clear
eststo: areg l_lmass l_rich i.styr_id, absorb(plst_id) cluster(plst_id)

** areg on selection equation **
eststo: areg l_rich i.styr_id, absorb(plst_id) cluster(plst_id)

esttab * using "$userdir/output/Table_S7.tex", fragment cells(b(fmt(%9.3f)) se(par) ci(par)) stats(r2 N, fmt(%9.3f %9.0g) label(R-squared)) legend label replace
est clear

psacalc beta l_rich, model(reg l_lmass l_rich i.plst_id i.styr_id, cluster(plst_id)) delta(-0.1) rmax(1)

