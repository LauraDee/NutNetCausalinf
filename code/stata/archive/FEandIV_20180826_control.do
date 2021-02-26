clear
cd "C:\Users\RNCNS02\Dropbox\Research_Active\IV in ecology\NutNet\working_CS\stata_work"
*cd "/Users/yggdrasil/Dropbox/Research_Active/IV in ecology/NutNet/working_CS/stata_work/"
insheet using NutNetControlPlotDataToUseApril2018.csv, c

/* Check characteristics of data */
destring total_mass live_mass rich ground_par ppm_* ph avg_neighbor_rich, i("NA") replace
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

duplicates report 	plst_id year
duplicates tag 		plst_id year, generate(duplicates)
tab 	duplicates
drop if duplicates==1 & year_trt==0

xtset plst_id year

/* Notes on variables
	ground_par varies at plot-yr level (but 1/3 of obs missing)
	ppm_P and ppm_K vary at plot level, but not across time
*/

gen 	res_sampleA = 1
replace res_sampleA = 0 if ground_par==. | ppm_p==. | ppm_k==.

gen 	res_sampleB = 1
replace res_sampleB = 0 if ground_par==.

****************************************
/* Cross sectional year-by-year regressions, no FE */

tab year
tab year_actual

** YR by YR **
* Levels
reg 	live_mass 	rich 							if year==4, robust
eststo t1a_1
reg 	live_mass 	rich 							if year==5, robust
eststo t1a_2
reg 	live_mass 	rich 							if year==6, robust
eststo t1a_3
reg 	live_mass 	rich 							if year==7, robust
eststo t1a_4

* Log-Log
reg 	l_lmass 	l_rich	 						if year==4, robust
eststo t1a_5
reg 	l_lmass 	l_rich	 						if year==5, robust
eststo t1a_6
reg 	l_lmass 	l_rich	 						if year==6, robust
eststo t1a_7
reg 	l_lmass 	l_rich	 						if year==7, robust
eststo t1a_8

estout t1a_* using table1a_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)


****************************************
/* Cross sectional year-by-year regressions, site(yr) FE */

** YR by YR **
* Levels
reghdfe live_mass 	rich 							if year==4, a(site_code2) vce(robust)
eststo t1b_1
reghdfe live_mass 	rich 							if year==5, a(site_code2) vce(robust)
eststo t1b_2
reghdfe live_mass 	rich 							if year==6, a(site_code2) vce(robust)
eststo t1b_3
reghdfe live_mass 	rich 							if year==7, a(site_code2) vce(robust)
eststo t1b_4

* Log-Log
reghdfe l_lmass 	l_rich	 						if year==4, a(site_code2) vce(robust)
eststo t1b_5
reghdfe l_lmass 	l_rich	 						if year==5, a(site_code2) vce(robust)
eststo t1b_6
reghdfe l_lmass 	l_rich	 						if year==6, a(site_code2) vce(robust)
eststo t1b_7
reghdfe l_lmass 	l_rich	 						if year==7, a(site_code2) vce(robust)
eststo t1b_8

estout t1b_* using table1b_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)
	

*********************************************
/* Cross sectional regressions, site-yr FE */

* Levels
reghdfe live_mass 	rich,			 				a(styr_id) cluster(plst_id)
eststo t2c_1
reghdfe live_mass 	rich			 				if res_sampleA==1, a(styr_id) cluster(plst_id)
eststo t2c_2
reghdfe live_mass 	rich ground_par ppm_p ppm_k,	a(styr_id) cluster(plst_id)
eststo t2c_3

* Log-log
reghdfe l_lmass 	l_rich,			 				a(styr_id) cluster(plst_id)
eststo t2c_4
reghdfe l_lmass 	l_rich			 				if res_sampleA==1, a(styr_id) cluster(plst_id)
eststo t2c_5
reghdfe l_lmass 	l_rich ground_par ppm_p ppm_k,	a(styr_id) cluster(plst_id)
eststo t2c_6

estout t2c_* using table2c_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)

** HOW MUCH IS EXPLAINED BY WEATHER **

/* Are there not time varying site levels measures of weather? */


*********************************************
/* Panel FE regressions, plot & site-yr FE */
/* exclude ppm_P and ppm_K bc time-invariant */

* Levels
reghdfe live_mass 	rich,			 	a(plst_id styr_id) cluster(plst_id)
eststo t3b_1
reghdfe live_mass 	rich			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id)
eststo t3b_2
reghdfe live_mass 	rich ground_par,	a(plst_id styr_id) cluster(plst_id)
eststo t3b_3

* Log-log
reghdfe l_lmass 	l_rich,			 	a(plst_id styr_id) cluster(plst_id)
eststo t3b_4
reghdfe l_lmass 	l_rich			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id)
eststo t3b_5
reghdfe l_lmass 	l_rich ground_par,	a(plst_id styr_id) cluster(plst_id)
eststo t3b_6

estout t3b_* using table3b_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)
	
reghdfe l_lmass if l_rich!=.,			 	a(plst_id styr_id) cluster(plst_id)		
reghdfe l_lmass 	l_rich,			 	a(plst_id styr_id) cluster(plst_id)	
reghdfe l_lmass 	l_rich ground_par,	a(plst_id styr_id) cluster(plst_id)

reghdfe l_rich,			 	a(plst_id) cluster(plst_id)	
reghdfe l_rich,			 	a(styr_id) cluster(plst_id)	
reghdfe l_rich,			 	a(plst_id styr_id) cluster(plst_id)	


reg  l_lmass l_rich if l_rich!=.,			 	 cluster(plst_id)	
reghdfe l_lmass l_rich if l_rich!=.,			 	a(plst_id) cluster(plst_id)	
reghdfe l_lmass l_rich if l_rich!=.,			 	a(styr_id) cluster(plst_id)		
reghdfe l_lmass l_rich if l_rich!=.,			 	a(plst_id styr_id) cluster(plst_id)		


*********************************************
/* Dynamic checks, include lagged dependent variable instead of FE, with site-yr FE */

* Levels
reghdfe live_mass 	rich L.live_mass,			 	a(styr_id) cluster(plst_id)
eststo t4b_1
reghdfe live_mass 	rich L.live_mass			 	if res_sampleB==1, a(styr_id) cluster(plst_id)
eststo t4b_2
reghdfe live_mass 	rich L.live_mass ground_par,	a(styr_id) cluster(plst_id)
eststo t4b_3

* Log-log
reghdfe l_lmass 	l_rich L.l_lmass,			 	a(styr_id) cluster(plst_id)
eststo t4b_4
reghdfe l_lmass 	l_rich L.l_lmass			 	if res_sampleB==1, a(styr_id) cluster(plst_id)
eststo t4b_5
reghdfe l_lmass 	l_rich L.l_lmass ground_par,	a(styr_id) cluster(plst_id)
eststo t4b_6

estout t4b_* using table4b_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)


*********************************************
/* Dynamics check, use lagged independent variable, plot & site-yr FE */

* Levels
reghdfe live_mass 	l.rich,			 	a(plst_id styr_id) cluster(plst_id)
eststo t5b_1
reghdfe live_mass 	l.rich			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id)
eststo t5b_2
reghdfe live_mass 	l.rich ground_par,	a(plst_id styr_id) cluster(plst_id)
eststo t5b_3

* Log-log
reghdfe l_lmass 	l.l_rich,			 	a(plst_id styr_id) cluster(plst_id)
eststo t5b_4
reghdfe l_lmass 	l.l_rich			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id)
eststo t5b_5
reghdfe l_lmass 	l.l_rich ground_par,	a(plst_id styr_id) cluster(plst_id)
eststo t5b_6

estout t5b_* using table5b_`date'.tex, style(tex) replace ///
	starlevels(+ 0.10 * 0.05 ** 0.01) cells(b(star fmt(3)) ///
	se(par fmt(3))) stats(r2 N)


*********************************************
/* Panel IV regressions, plot & site-yr FE */
/* using current neighbor richness as IV */

* Levels
reghdfe live_mass 	(rich = avg_neighbor_rich),				a(plst_id styr_id) cluster(plst_id) ff
reghdfe live_mass 	(rich = avg_neighbor_rich)			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id) ff
reghdfe live_mass 	ground_par (rich = avg_neighbor_rich),	a(plst_id styr_id) cluster(plst_id) ff

* Log-log
reghdfe l_lmass 	(l_rich = l_nbrich),				a(plst_id styr_id) cluster(plst_id) ff
reghdfe l_lmass 	(l_rich = l_nbrich)			 		if res_sampleB==1, a(plst_id styr_id) cluster(plst_id) ff
reghdfe l_lmass 	ground_par (l_rich = l_nbrich),		a(plst_id styr_id) cluster(plst_id) ff


*********************************************
/* Panel IV regressions, plot & site-yr FE */
/* using lagged neighbor richness as IV */

* Levels
reghdfe live_mass 	(rich = l.avg_neighbor_rich),				a(plst_id styr_id) cluster(plst_id) ff
reghdfe live_mass 	(rich = l.avg_neighbor_rich)			 	if res_sampleB==1, a(plst_id styr_id) cluster(plst_id) ff
reghdfe live_mass 	ground_par (rich = l.avg_neighbor_rich),	a(plst_id styr_id) cluster(plst_id) ff

* Log-log
reghdfe l_lmass 	(l_rich = l.l_nbrich),				a(plst_id styr_id) cluster(plst_id) ff
reghdfe l_lmass 	(l_rich = l.l_nbrich)				if res_sampleB==1, a(plst_id styr_id) cluster(plst_id) ff
reghdfe l_lmass 	ground_par (l_rich = l.l_nbrich),	a(plst_id styr_id) cluster(plst_id) ff


*********************************************
/* Investigating/decomposing Fixed effects */

* Run main specification, but save plot FE values and residuals
reghdfe l_lmass 	l_rich ground_par,	a(plotfe=plst_id styr_id) cluster(plst_id)
reghdfe l_lmass 	l_rich ground_par,	a(plst_id styr_id) cluster(plst_id) res(plotyrres)
bys plst_id: egen plotmeanres=mean(plotyrres)
sum plotmeanres, d /* Should be 0 with some computational error */

* Do some cell copying to fill out missing years 
bys plst_id: egen plotfe_fill = mean(plotfe)

* Similar action with soil values doesn't improve situation


*Year 4 5 6 7 have the most non-zeros /* tab sitecode2 year */

reghdfe plotfe_fill ppm_* ph if year==4, a(site_code2) vce(robust)
reghdfe plotfe_fill ppm_* ph if year==5, a(site_code2) vce(robust)
reghdfe plotfe_fill ppm_* ph if year==6, a(site_code2) vce(robust)
reghdfe plotfe_fill ppm_* ph if year==7, a(site_code2) vce(robust)

reghdfe plotfe_fill ppm_* ph, a(styr_id) cluster(plst_id)
reg plotfe_fill ppm_* ph, cluster(plst_id)


**************
log close

