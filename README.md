# NutNetCausalinf
This repository contains a tutorial for the methods in and replication code for Increases in grassland biodiversity in 11 countries decrease average ecosystem productivity by Dee, Ferraro, Severen et al.

### Tutorial Instructions

To access the tutorial, see the folder "./tutorial in R/". The folder includes an html ("Ecology_FEs_for_Causality.html") and a Rmarkdown file that uses the data "cleaned_comb_data.csv" in the data.zip file for reviewers. The html can be accessed without using this data. 

### Replication Instructions

The analysis for this project was done in R and Stata, with the majority in R. Below, under *Provenance of Results*, is a list of all results in the paper and supplemental materials and where the authoritative version can be found. 

R Replication Instructions: Set cdir to your local repo location. For the review process, the data files are includes in the data.zip for reviewers as an attachments. Post acceptance: Note that only *comb* data is included; *cover* is larger and must be obtained separately, and requires settings its own path. The master.R runs through each step in the code and data processing for analysis; readers should start there.

Stata Replication Instructions: Set the global macros to local repo location.

#### Provenance of Results

Figure 2A

* All results in R 
* Figure produced in R
* Replicable in Stata 

Figure 2B

* All results in Stata 
* Figure produced in R 
* NOT replicated in R

Figure 3 

* Results 1, 2, and 4 in R 
* Result 3 and 5 in Stata
* Replicable in Stata 
Figure 5

* Results only in R
* Figure produced in R

Table S1

* Dataset Description 

Table S2

* Results and Table produced in R (analysis_main.R) 
* Replicable in Stata 

Table S3

* Results and Table produced in Stata and R (analysis_main.R)

Table S4

* Results and Table produced in R (analysis_sm5.R)

Table S5

* Results and Table produced in R (analysis_sm5.R)

Table S6 

* Results and Table produced in R (analysis_sm5.R)

Table S7

* Results and Table produced in Stata 

Table S8

* Results and Table produced in Stata and R (analysis_main.R)

Table S9

* Results and Table produced in R (analysis_main.R)
* Replicable in Stata 

Table S10

* Results and Table produced in R (analysis_fig5_smsection9.R)

Table S11

* Results and Table produced in R (analysis_fig5_smsection9.R)

Table S12

* Results and Table produced in R (analysis_fig5_smsection9.R)

Table S13

* Results and Table produced in R (analysis_fig5_smsection9.R)

Table S14

* Results and Table produced in R (analysis_fig5_smsection9.R)


SM Section 9: Cedar Creek Analysis 
* Results, Figures and Tables produced in R
* Table S15 -- produced in R (S10analyses_e120-planted_rich-mass.R)
* Table S16 -- produced in R (Biocon_rich_biomass_FINAL_NUTNETPAPER.R)

Supplemental Figures

* Figure S1 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S2 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S3 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S10 -- Produced in R (rank_abundance_curves.R)
* Figure S11 -- Produced in R  (finalprocess_coverdata.R)
* Figure S12 -- Produced in R  (analysis_fig5_smsection9.R)
* Figure S13 -- Produced in R  (analysis_fig5_smsection9.R)
* Figure S14 -- Produced in R (Biocon_rich_biomass_FINAL_NUTNETPAPER.R and S10analyses_e120-planted_rich-mass.R) 

