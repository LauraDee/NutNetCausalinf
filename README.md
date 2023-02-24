# NutNetCausalinf
This repository contains a tutorial for the methods in and replication code for "Clarifying the effect of biodiversity on productivity in natural ecosystems with longitudinal data and methods for causal inference" Nature Communicatins, by Dee, Ferraro, Severen et al.

### Tutorial Instructions

To access the tutorial, see the folder "./tutorial in R/". This folder includes an html ("RMarkdown_panel_designs_Dee_and_Severen.html") and a Rmarkdown file that uses the data "cleaned_comb_data.csv". The html can be accessed without using this data. 

### Replication Instructions

The analysis for this project was done in R and Stata, with the majority in R. Below, under *Provenance of Results*, is a list of all results in the paper and supplemental materials and where the authoritative version can be found. 

R Replication Instructions: Set cdir to your local repo location. For the review process, the data files are includes in the data.zip for reviewers as an attachments.The master.R runs through each step in the code and data processing for analysis; readers should start there.

Stata Replication Instructions: Set the global macros to local repo location.

Data availability statement: The processed data from the Nutrient Network generated in this study have been deposited in this repository.  The raw data for unmanipulated plots that were not included in the analyses, because they did not meet the inclusion criteria, are available under restricted access for which permission can be obtained by contacting the Nutrient Network at https://nutnet.org. The code for full data processing pipeline is available in this repository in the folder "./code/r/ProcessDataFromRaw/".

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

* Dataset Description (.csv produced from the code run in master.R)

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

Supplemental Figures

* Figure S1 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S2 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S3 -- Produced in R (FiguresS1_S2_S3.R) 
* Figure S10 -- Produced in R (rank_abundance_curves.R), code available.
* Figure S11 -- Produced in R  (finalprocess_coverdata.R)
* Figure S12 -- Produced in R  (analysis_fig5_smsection8.R)
* Figure S13 -- Produced in R  (analysis_fig5_smsection8.R)

