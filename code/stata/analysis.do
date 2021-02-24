use 	"$datadir/processed/NutNet_Prepped.dta"

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

reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)
reghdfe l_lmass l_rich ihs_even, keepsingletons  absorb(plst_id styr_id) cluster(plst_id)
reghdfe l_lmass l_simpson, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

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
** Figure 2B: Bivariate and Traditional Multi-variate Approaches
//
////////////////////////////////////////////////////////////////////////////////

local control_vars i.pais i.tierra i.year elevation managed burned grazed anthropogenic ///
					temp_var_v2  min_temp_v2 max_temp_v2 temp_wet_q_v2 temp_dry_q_v2 temp_warm_q_v2 ///
					temp_cold_q_v2 pct_c2 pct_n2 ppm_p ppm_k ppm_ca ppm_mg ppm_s ppm_na ppm_zn ///
					ppm_mn ppm_fe ppm_cu ppm_b ph percentsand_2 percentsilt_2 percentclay_2 

xtmixed l_lmass l_rich `control_vars' || site_code2: || plst_id:, vce(robust)
xtmixed l_lmass l_rich ihs_even `control_vars' || site_code2: || plst_id:, vce(robust)
xtmixed l_lmass l_simpson `control_vars' || site_code2: || plst_id:, vce(robust)

** Possible to Replication with <lmer> in R? 

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


////////////////////////////////////////////////////////////////////////////////
//
** Figure 3: Robustness across various designs
//
////////////////////////////////////////////////////////////////////////////////


** Main Study Design
reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Dynamic Panel Design Lagged Dependent Variable
reghdfe l_lmass l_rich L.l_lmass, keepsingletons absorb(styr_id) cluster(plst_id)

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
** Table S4 (Functional Form Tests)
//
////////////////////////////////////////////////////////////////////////////////


** Log-Log
reghdfe l_lmass l_rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Log-Level
reghdfe l_lmass rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Level-Level
reghdfe live_mass rich, keepsingletons absorb(plst_id styr_id) cluster(plst_id)

** Level-Quadratic
reghdfe live_mass rich rich2, keepsingletons absorb(plst_id styr_id) cluster(plst_id)





////////////////////////////////////////////////////////////////////////////////
//
** Table S5 Not in SM
//
////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////
//
** Table S11: IV (replicates IV results of Table 3)
//
////////////////////////////////////////////////////////////////////////////////
ivreghdfe l_lmass (l_rich = l_nbrichblock), a(plst_id styr_id) cluster(plst_id) first





////////////////////////////////////////////////////////////////////////////////
//
** Table S12: Lagged Dependent Varaibles (also replicates a result in Table 3)
//
////////////////////////////////////////////////////////////////////////////////

reghdfe l_lmass l_rich L.l_lmass, keepsingletons a(styr_id) cluster(plst_id)
reghdfe l_lmass l_rich L.l_lmass ihs_even, keepsingletons a(styr_id) cluster(plst_id)

reghdfe l_lmass l_rich L.l_lmass L.l_rich, keepsingletons a(styr_id) cluster(plst_id)
reghdfe l_lmass l_rich L.l_lmass L.l_rich ihs_even, keepsingletons a(styr_id) cluster(plst_id)


