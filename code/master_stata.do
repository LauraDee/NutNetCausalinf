////////////////////////////////////////////////////////////////////////////////
// Version Date: 02.23.2021
// Paper: Dee et al. "Increases in Species Richness Decrease Average Eco..."
//
// Code Custody:
// C Severen original code named: FEandIV_20180826_control
// P Ferraro greatly revised code: FEandIV_20180826_control_REVISEDforSens11
// C Severen revisited for replication and cleaning: current
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// To replication results:
//
////////////////////////////////////////////////////////////////////////////////

/*-----------------------------------------------*/
/* Set Global User Directory */
/*-----------------------------------------------*/

global userdir "C:\GitHub\NutNetCausalinf"
clear all

////////////////////////////////////////////////////////////////////////////////
// Setup

** Set useful directories consistent with repo
global 	datadir "$userdir/data"
global 	codedir "$userdir/code"
global 	logdir  "$userdir/output/log"

** Check for required SSC packages and install if not present
do 		"$codedir/stata/config_stata.do" 

** Useful data macros
local date: display %td_CCYY_NN_DD date(c(current_date), "DMY")
local date_string = subinstr(trim("`date'"), " " , "", .)

** Initial log file
log using "$logdir/stata_results_d`date_string'", replace 

////////////////////////////////////////////////////////////////////////////////
// Main Block

do 		"$codedir/stata/final_processing.do"

do 		"$codedir/stata/analysis.do"

////////////////////////////////////////////////////////////////////////////////
// Clean up

log close
translate "$logdir/stata_results_d`date_string'.smcl" "$logdir/stata_results_d`date_string'.log", replace

clear all