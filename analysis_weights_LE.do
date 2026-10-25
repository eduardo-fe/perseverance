/*
The program "run_model_weights" is designed specifically for the current analysis and is not
intended as a general-purpose IPW routine. It relies on the structure of
our data: perseverance (`perserverance_2_sd') and cognitive skills (`pc1')
are the variables of interest with missing values, and the observation
indicator is defined using these measures and the outcome being analysed.
The outcome models include both measures as regressors, so Stata
automatically excludes observations with missing values in either measure.
The inverse-probability weights are based on a probit model for this
complete-case indicator. Any adaptation to other data or missingness
patterns requires reviewing the indicator, selection model and weights.

IMPORTANT: Do not use the if/in qualifiers when calling "run_model_weights".
The selection model is not designed to handle sample restrictions
passed through these qualifiers. To restrict the analysis sample,
use preserve, restrict the data (e.g. with keep if ...), run the
program, and then restore the original data. Ensure that the
restriction is appropriate for the intended target population.

*/

* Code for archiving


* The following is from "cinelli_extended.do"
/*
 * Heteroskedasticity-robust standard errors for the adapted
 * Cinelli-Hazlett sensitivity analysis.
 *
 * The adjusted standard errors are derived using a delta-method
 * approximation within the weighted-regression framework.
 *
 * Note: the published manuscript contains a typo in the definition
 * of the sandwich covariance matrix in the Appendix. The implementation
 * in this code is intended to use the correct expression.
 *
 * However, the variance estimates reported by this code and in
 * Appendix C may not fully account for the sampling variability
 * of the estimated sensitivity quantities (the tilde-R-squared
 * measures), which are not the usual regression R-squared measures.
 * In particular, their sampling variability and covariance with
 * the other estimated quantities may not be fully reflected in
 * the delta-method calculation. The resulting standard errors and
 * confidence intervals should therefore be treated with caution
 * until their coverage properties have been assessed.
 */

capture program drop senseweighter_hetero
program define senseweighter_hetero, rclass
    syntax varlist(min=2 max=2) [if] , ///
		KEYvars(string) ///
        Xvar(string) ///
        METHOD(string) ///
		PSVARS(string) ///
		[BENCHmark(varname)] ///
		[KD(real 1.0)] ///
		[KY(real 1.0)] ///
		 

	qui{ 
		tempvar Y D indicator ipw xb Ystar mX ///
            Y_perp_DX D_perp_X u_Y u_D

		tokenize `varlist'
		gen `Y' = `1'
		gen `D' = `2'
		
	    marksample touse
     
        * -----------------------
        * Selection / missingness model
        * -----------------------
            
		if "`psvars'"=="" {
			di as error "psvars() must be specified to estimate attrition weights"
			exit 198
		}
				
		gen `indicator' = 1
		
		replace `indicator' = 0 if missing(`Y') | missing(`D')
 
		foreach x of local keyvars {
				replace `indicator' = 0 if missing(`x')
		}
						
    
		qui probit `indicator' `psvars'
		predict `xb', pr
		qui gen `ipw' = `indicator'/`xb'
		qui sum `ipw' if `indicator' ==1 
		replace `ipw' = `ipw' / r(mean) 
		
		
		* -----------------------
        * Outcome or AIPW pseudo-outcome
        * -----------------------
        
        if "`method'"=="aipw" {
			 qui reg `Y' `xvar' [aw = `ipw'] if `indicator'==1
             predict double `mX'
             qui gen double `Ystar' = `mX' + `ipw' * (`Y' - `mX')
			 
        }
        else {
			gen `Ystar' = `Y'
        }
		
		* -----------------------
        * Point estimates (full sample)
        * -----------------------
		
		* Regression with robust SE
		qui regress `Ystar' `D' `xvar' if `indicator'==1, robust
		local tau_res = _b[`D']
		local se_tau = _se[`D']
		local dof = e(df_r)
		
		* Get residuals for outcome regression
		predict `u_Y', resid
		
		qui regress `D' `xvar'   if `indicator'==1, robust
		predict `D_perp_X', resid
		
		* Get residuals for treatment regression
		predict `u_D', resid
		
		* AUX
		qui reg `Ystar' `xvar' if `indicator' == 1
		local aux1 = e(r2) 
		qui reg `Ystar' `D' `xvar' if `indicator' == 1
		local aux2 = e(r2)
		local partYD = (`aux2' - `aux1')/(1-`aux1')
		
		qui regress `Ystar' `D' `xvar' if `indicator'==1
		predict `Y_perp_DX', resid

		* Compute weighted standard deviations
		su `Y_perp_DX' if `indicator'==1
		local sd_Y = sqrt(r(Var))
		
		su `D_perp_X'  if `indicator'==1
		local sd_D = sqrt(r(Var))
		 
		* Compute RV
		local kappa = abs(`tau_res') * `sd_D' / `sd_Y'
		local kappa2 = `kappa'^2
		local kappa4 = `kappa'^4
		local RV = (sqrt(`kappa4' + 4*`kappa2') - `kappa2')/2
	}
			
		
	if "`benchmark'" != ""  {
		
		local X_minus_j : list xvar - benchmark
		
	
		 * ---- Bounding R2_D~Z|X
		qui{ 
			reg `D' `X_minus_j'  if `indicator'==1
			local R2_D_Xminusj = e(r2)

			reg `D' `X_minus_j' `benchmark'  if `indicator'==1
			local R2_D_X = e(r2)

			local partialR2_D_Xj = ///
				(`R2_D_X' - `R2_D_Xminusj') / (1 - `R2_D_Xminusj')
			
			local f2_D_Xj = ///
				`partialR2_D_Xj' / (1 - `partialR2_D_Xj')
				
			local f_D_Xj = sqrt(`f2_D_Xj') 
				
			local R2_D_Z_condX = `kd' * `f2_D_Xj'
		}
		
		
		* ---- Bounding R2_Y~Z|X
		 
		* Part A: f2_Y~X_j|X_-j, D
		* Step 1: Regress Y on controls + D only
		qui{
			reg `Ystar' `D' `X_minus_j'  if `indicator'==1
			local R2_Y_Xminusj_D = e(r2)

			* Step 2: Regress Y on controls + D + benchmark
			reg `Ystar' `D' `X_minus_j' `benchmark'  if `indicator'==1
			local R2_Y_full = e(r2)

			* Step 3: Compute partial R^2 of benchmark
			local partialR2_Y_Xj_D = (`R2_Y_full' - `R2_Y_Xminusj_D') / (1 - `R2_Y_Xminusj_D')

			* Step 4: Convert to f^2
			local f2_Y_Xj_D = `partialR2_Y_Xj_D' / (1 - `partialR2_Y_Xj_D')
			
			* Part B: f_Kd
			
			local f2_Kd = `kd'*`partialR2_D_Xj' /(1 -`kd'*`partialR2_D_Xj')
			local f_Kd = sqrt(`f2_Kd')
			
			* Part C: eta
			
			local eta = (sqrt(`ky') + abs(`f_Kd'* `f_D_Xj') )/ (sqrt(1-`f2_Kd'*`f2_D_Xj'))
			
			* Part D: bound
			
			local R2_Y_Z_condDX = `eta'^2* `f2_Y_Xj_D'
		}
		
		local bias =(`sd_Y' / `sd_D')* sqrt( (`R2_Y_Z_condDX' * `R2_D_Z_condX')/(1- `R2_D_Z_condX'))
		local tau_adjusted = `tau_res' - sign(`tau_res') * abs(`bias')
		
		
		******************************************************************************
		* HETEROSKEDASTICITY-ROBUST ADJUSTED SE (Delta Method)
		******************************************************************************
		
		* Step 1: SE of treatment effect (already have from robust regression)
		local se_tau_robust = `se_tau'
		
		* Step 2: Compute SE of outcome residual SD
		* Var(sigma_Y^2) using HC-robust formula
		
		tempvar u_Y_sq u_D_sq
		gen `u_Y_sq' = `u_Y'^2
		gen `u_D_sq' = `u_D'^2
		
		qui sum `u_Y_sq' if `indicator'==1
		local mean_u_Y_sq = r(mean)
		local n_complete = r(N)
		
		* HC-robust variance of residual variance
		qui gen double temp_var_Y = (`u_Y_sq' - `mean_u_Y_sq')^2
		qui sum temp_var_Y if `indicator'==1
		local var_u_Y_sq = r(mean) / `n_complete'
		
		* SE of sigma_Y (delta method: d(sqrt(x))/dx = 1/(2*sqrt(x)))
		local se_sigma_Y = sqrt(`var_u_Y_sq') / (2 * `sd_Y')
		
		* Step 3: Compute SE of treatment residual SD
		qui sum `u_D_sq' if `indicator'==1
		local mean_u_D_sq = r(mean)
		
		qui gen double temp_var_D = (`u_D_sq' - `mean_u_D_sq')^2
		qui sum temp_var_D if `indicator'==1
		local var_u_D_sq = r(mean) / `n_complete'
		
		local se_sigma_D = sqrt(`var_u_D_sq') / (2 * `sd_D')
		
		* Step 4: Partial derivatives of bias bound w.r.t. sigma_Y and sigma_D
		* d(Bias)/d(sigma_Y) = (1/sigma_D) * sqrt(R2_Y * R2_D / (1 - R2_D))
		
		* Check for numerical issues: R2_D close to 1 or R2_Y/R2_D near 0
		if (`R2_D_Z_condX' >= 0.99) {
			local se_bias = .
			di as text "Warning: R²_D >= 0.99; SE(Bias) cannot be computed reliably"
		}
		else if (`R2_Y_Z_condDX' < 0.00001 | `R2_D_Z_condX' < 0.00001) {
			* If either R² is essentially 0, bias is ~0 with minimal uncertainty
			local se_bias = 0
		}
		else {
			local denom = (1 - `R2_D_Z_condX')
			
			* Compute the ratio carefully
			if (`denom' > 0.001) {
				local r2_ratio = sqrt((`R2_Y_Z_condDX' * `R2_D_Z_condX') / `denom')
				
				* Instead of computing partial derivatives separately,
				* use a simpler approximation based on the ratio of bias to SDs
				
				* The bias scales with the outcome and treatment SDs, so
				* SE(Bias) ≈ Bias * sqrt(SE(sigma_Y)^2/sigma_Y^2 + SE(sigma_D)^2/sigma_D^2)
				* This is a log-scale relative error formula
				
				local rel_se_sigma_Y = `se_sigma_Y' / `sd_Y'
				local rel_se_sigma_D = `se_sigma_D' / `sd_D'
				
				local rel_variance_bias = (`rel_se_sigma_Y'^2) + (`rel_se_sigma_D'^2)
				
				if (`rel_variance_bias' >= 0) {
					local se_bias = `bias' * sqrt(`rel_variance_bias')
				}
				else {
					local se_bias = .
					di as text "Warning: Numerical error in relative SE calculation"
				}
			}
			else {
				local se_bias = .
				di as text "Warning: 1 - R²_D too small; SE(Bias) is numerically unstable"
			}
		}
		
		* Step 6: Combined SE of adjusted estimate
		* SE(tau_adj) = sqrt(SE(tau)^2 + SE(Bias)^2)
		* Note: We ignore the covariance term as it's typically small
		
		if (missing(`se_bias')) {
			local se_adjusted = .
			local tau_adj_lo = .
			local tau_adj_hi = .
		}
		else {
			local se_adjusted = sqrt(`se_tau_robust'^2 + `se_bias'^2)
			
			* Step 7: Confidence intervals
			local tau_adj_lo = `tau_adjusted' - 1.96 * `se_adjusted'
			local tau_adj_hi = `tau_adjusted' + 1.96 * `se_adjusted'
		}
		
		******************************************************************************
		* Display Results
		******************************************************************************
		di as text _newline "Sensitivity Analysis Results (Heteroskedasticity-Robust)"
		di as text "{hline 70}"
		di as text "Method: " as result "`method'"
		di as text "Treatment effect: " as result %9.4f `tau_res'
		di as text "Standard error (robust): " as result %9.4f `se_tau_robust'
		di as text "Kappa: " as result %9.4f `kappa'
		di as text "Robustness Value (RV): " as result %9.4f `RV'
		di as text "R2yd.x: " as result %9.4f `partYD'
		di as text "{hline 70}"
	
		di as text _newline "Benchmark Bounds for Sensitivity Analysis"
		di as text "{hline 70}"
		di as text "Outcome: " as result "`1'"
		di as text "Treatment: " as result "`2'"
		di as text "Benchmark covariate: " as result "`benchmark'"

		di as text "{hline 70}"
		di as text _newline "Benchmark Covariate Explanatory Power:"
		di as text "  R̃²_{D ~ `benchmark' | X_{-j}}: " as result %9.4f `partialR2_D_Xj'
		di as text "  R̃²_{Y ~ `benchmark' | X_{-j}}: " as result %9.4f `partialR2_Y_Xj_D'
		di as text _newline "Confounder Strength Multipliers:"
		di as text "  k_D (treatment): " as result %9.2f `kd'
		di as text "  k_Y (outcome): " as result %9.2f `ky'
		di as text _newline "Implied Bounds on Confounder Z (assuming X ⊥ Z):"
		di as text "  R̃²_{D ~ Z | X}: " as result %9.4f `R2_D_Z_condX'
		di as text "  R̃²_{Y ~ Z | X}: " as result %9.4f `R2_Y_Z_condDX'
		di as text _newline "Treatment Effect:"
		di as text "  Estimated τ: " as result %9.4f `tau_res' as text " (SE: " as result %6.4f `se_tau_robust' as text ")"
		di as text "  |Bias| bound: " as result %9.4f `bias'
		di as text _newline "Bias bound SE components:"
		di as text "  SE(σ_Y): " as result %9.4f `se_sigma_Y'
		di as text "  SE(σ_D): " as result %9.4f `se_sigma_D'
		di as text "  SE(Bias): " as result %9.4f `se_bias'
		di as text "{hline 70}"
		
		di as text _newline "Adjusted Treatment Effect (Delta Method CIs):"
		di as text "  Adjusted τ: " as result %9.4f `tau_adjusted' 
		di as text "  Adjusted SE: " as result %9.4f `se_adjusted'
		di as text "  95% CI: [" as result %9.4f `tau_adj_lo' as text ", " as result %9.4f `tau_adj_hi' as text "]"
		di as text "{hline 70}"
		di as text "Interpretation: If an unobserved confounder Z explains at most"
		di as text "k_D=" as result `kd' as text " times as much variance in D and k_Y=" as result `ky' as text " times"
		di as text "as much variance in Y as `benchmark' (conditional on other X),"
		di as text "then the bias is bounded by ±" as result %6.4f `bias' as text "."
		di as text "{hline 70}"
		di as text _newline "{it:Note:} This SE computation uses the delta method applied to the bias bound formula,"
		di as text "accounting for heteroskedasticity from AIPW weighting. The covariance term between"
		di as text "tau and the bias bound is neglected as it is typically small."
		di as text "{hline 70}"
		
		 
	 
	******************************************************************************
	* Return Values
	******************************************************************************
	
	return scalar partialR2_D_Xj = `partialR2_D_Xj'
	return scalar partialR2_Y_Xj = `partialR2_Y_Xj_D'
	return scalar bound_R2_D_Z = `R2_D_Z_condX'
	return scalar bound_R2_Y_Z = `R2_Y_Z_condDX'
	return scalar tau_res = `tau_res'
	return scalar se_tau = `se_tau_robust'
	return scalar tau_adjusted = `tau_adjusted'
	return scalar se_adjusted = `se_adjusted'
	return scalar tau_lower = `tau_adj_lo'
	return scalar tau_upper = `tau_adj_hi'
	return scalar kd = `kd'
	return scalar ky = `ky' 
	return scalar bias = `bias'
	return scalar se_bias = `se_bias'
	return scalar RV = `RV'
	return scalar R2ydx = `partYD'
	
	drop temp_var_Y temp_var_D `u_Y_sq' `u_D_sq'
	}
	
end

* Taken from files "analysis_weights.do" and "cinelli_extended.do"

* The following is from "analysis_weights.do"
capture program drop run_models_weights
program define run_models_weights, rclass
    syntax [if] ,KEYvars(string) Binary(string) Continuous(string) COVariates(string) 
	qui {
			
		tempname mRES
		matrix `mRES' = J(100, 7, .)  // allow up to 100 outcomes; resize later
		local i = 1
		local rownames
		
		marksample touse

		// Loop over binary outcomes using xpologit
		foreach var of local binary {
			
			* Create Weights
			gen indicator = 0
			replace indicator = 1 if perserverance_1_sd!=. & pc1 != . & `var' !=.
			probit indicator $set1 $set2 $set3  
			predict xb, pr
			gen ipw = 1/xb

			
			logit `var' `keyvars' `covariates' [pweight=ipw] if `touse', cluster(bcsid)
				
			margins, dydx(`keyvars')

			matrix aux = r(table)
			
			matrix `mRES'[`i', 1] = aux[1,1]
			matrix `mRES'[`i', 2] = aux[2,1]
			matrix `mRES'[`i', 3] = aux[4,1]	
			matrix `mRES'[`i', 4] = aux[1,2]		 
			matrix `mRES'[`i', 5] = aux[2,2]	     
			matrix `mRES'[`i', 6] = aux[4,2]
			matrix `mRES'[`i', 7] = r(N)

	  
			local rownames `"`rownames' `var'"'
			local i = `i' + 1
			drop indicator xb ipw
		}

		// Loop over continuous outcomes using poregress
		foreach var of local continuous {
			
			* Create weights
			gen indicator = 0
			replace indicator = 1 if perserverance_1_sd!=. & pc1 != . & `var' !=.
			probit indicator $set1 $set2 $set3  
			predict xb
			gen ipw = 1/xb
			
			regress `var' `keyvars' `covariates' [pweight=ipw] if `var' > 0 &`touse', cluster(bcsid)

			margins, dydx(`keyvars')

			matrix aux = r(table)
			
			matrix `mRES'[`i', 1] = aux[1,1]
			matrix `mRES'[`i', 2] = aux[2,1]
			matrix `mRES'[`i', 3] = aux[4,1]	
			matrix `mRES'[`i', 4] = aux[1,2]		 
			matrix `mRES'[`i', 5] = aux[2,2]	     
			matrix `mRES'[`i', 6] = aux[4,2]
			matrix `mRES'[`i', 7] = r(N)

			local rownames `"`rownames' `var'"'
			local i = `i' + 1
			drop indicator xb ipw
		}

		local rowcount = `i' - 1
		matrix `mRES' = `mRES'[1..`rowcount', 1..7]
		matrix rownames `mRES' = `rownames'
		matrix colnames `mRES' = dydx1 se_z pval dydx2 se_z pval N
	}
    matrix list `mRES'
    
    // Return matrix
    return matrix results = `mRES'
end
 

 


* ------------------------------------------------------------------------------
* Load the covariates (run these files in that order 3, 1, 2)
* ------------------------------------------------------------------------------

qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)


* ------------------------------------------------------------------------------
* Who perseveres?
* ------------------------------------------------------------------------------

gen indicator = 0
replace indicator = 1 if perserverance_1_sd!=. & pc1 != . 
probit indicator $set1 $set2 $set3  c_delivery_age_mum_2
predict xb
gen ipw = 1/xb
		
regress perserverance_2_sd pc1 $set1 $set2  $set3  c_delivery_age_mum_2 [pweight=ipw] , cluster(bcsid)

estimates store model1

esttab model1 using "/Users/user/Dropbox/Econometrics/perserverance/regression_table.tex", replace ///
    booktabs se star(* 0.10 ** 0.05 *** 0.01) ///
    label title("Regression of Perseverance on Principal Components and Controls") ///
    mtitles("Model 1") scalars(r2 N) ///
    coeflabels(pc1 "First Principal Component" _cons "Constant") 


drop indicator xb ipw

* ------------------------------------------------------------------------------
* Estimate models 
* ------------------------------------------------------------------------------


run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary(ismanager employed employedAt46 employedAt51 unemployedAnytime ///
               manual_worker profORmanagerialClass degree noqualifications ///
               noLiquidSavingsAt46 noLiquidSavingsAt51 personalLoan46 personalLoan51 ///
               noDebtAt46 ismanagerAt46 luckvar_k080 luckvar_k088 luckvar_k094 luckvar_k075 ///
			   skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a skill_skill6a ///
			   skill_skill7a skill_skill8a skill_skill9a good_at_k036 good_at_k037 good_at_k038 ///
			   good_at_k039 good_at_k040 good_at_k041 good_at_k042 good_at_k043) ///
    continuous(loginc logWeeklyInc46 logTotalSavingsCapped51) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6  )

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/skills", mat(mResults) nobox format(%9.3f) replace

 

* ------------------------------------------------------------------------------
* Re-estimate models using the additional principal components.
* ------------------------------------------------------------------------------

run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary(ismanager employed employedAt46 employedAt51 unemployedAnytime ///
               manual_worker profORmanagerialClass degree noqualifications ///
               noLiquidSavingsAt46 noLiquidSavingsAt51 personalLoan46 personalLoan51 ///
               noDebtAt46 ismanagerAt46 luckvar_k080 luckvar_k088 luckvar_k094 luckvar_k075 ///
			   skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a skill_skill6a ///
			   skill_skill7a skill_skill8a skill_skill9a good_at_k036 good_at_k037 good_at_k038 ///
			   good_at_k039 good_at_k040 good_at_k041 good_at_k042 good_at_k043) ///
    continuous(loginc logWeeklyInc46 logTotalSavingsCapped51) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/skills_set7.tex", mat(mResults) nobox format(%9.3f) replace



* ------------------------------------------------------------------------------
* Labour market outcomes from the panel.
* ------------------------------------------------------------------------------




preserve

merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
gen logwage = log(net_pay_week_95)

run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary( employedFT unemployed isManagerSupervisor) ///
    continuous(hours net_pay_week_95) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 )

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour.tex", mat(mResults) nobox format(%9.3f) replace

run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary( employedFT unemployed isManagerSupervisor) ///
    continuous(hours logwage) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour_set7.tex", mat(mResults) nobox format(%9.3f) replace

restore


 
		
**# Sensitivity analysis


qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)

* Redefine the globals to avoid the "i." operator (fwl cannot handle it for some reason).

tab c_region_birth, gen(region_)
tab c_sexbirth, gen(sex_)
global set1 c_delivery_age_mum c_delivery_age_mum_2 region_1 region_2 region_3 region_4 region_5 region_6 region_7 region_8 region_9 region_10 c_mumsmoked c_abnormal_gest sex_2 c_birthweight 


global set4 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband1 c_incband2 c_incband3 c_incband4 c_incband5 c_incband6 c_incband7 c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10

 

********************************************************************************
* Loop over outcomes and store sensitivity results
********************************************************************************

* Define list of outcome variables
local outcomes luckvar_k075 luckvar_k080 luckvar_k094 ///
	skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a /// 
	skill_skill6a skill_skill7a  skill_skill9a ///
	degree noqualifications ///
	good_at_k036 good_at_k037 good_at_k038 good_at_k039 good_at_k040 



* Initialize empty matrix to store results
matrix results = J(1, 10, .)
matrix drop results  // Drop to start fresh
local row_names ""

* Loop over each outcome
local i = 1
foreach outcome of local outcomes {
    
    di as text _newline(2) "{hline 80}"
    di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
    di as text "{hline 80}"
    
    * Run sensitivity analysis
    capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
        xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
        psvars($set1 $set2 $set3) ///
        method("aipw") ///
        keyvars(pc1 perserverance_2_sd) ///
        benchmark(pc1)
    
    * Store results if successful
    if _rc == 0 {
        matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
                          r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
                          r(tau_lower), r(tau_upper), r(se_adjusted))
        
        if `i' == 1 {
            matrix results = row_`i'
        }
        else {
            matrix results = results \ row_`i'
        }
        
        local row_names "`row_names' `outcome'"
    }
    else {
        di as error "Failed to process `outcome' - skipping"
    }
    
    local i = `i' + 1
}

* Label matrix
matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
matrix rownames results = `row_names'

* Display summary table
di as text _newline(2) "{hline 100}"
di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
di as text "{hline 100}"
matlist results, format(%9.4f) border(all) rowtitle("Outcome")
outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_edu.tex", mat(results) nobox format(%9.3f) replace

**# Sensitivity Panel 


merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
gen logwage = log(net_pay_week_95)

local outcomes employedFT unemployed isManagerSupervisor hours logwage


* Initialize empty matrix to store results
matrix results = J(1, 10, .)
matrix drop results  // Drop to start fresh
local row_names ""

* Loop over each outcome
local i = 1
foreach outcome of local outcomes {
    
    di as text _newline(2) "{hline 80}"
    di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
    di as text "{hline 80}"
    
    * Run sensitivity analysis
    capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
        xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
        psvars($set1 $set2 $set3) ///
        method("aipw") ///
        keyvars(pc1 perserverance_2_sd) ///
        benchmark(pc1)
    
    * Store results if successful
    if _rc == 0 {
        matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
                          r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
                          r(tau_lower), r(tau_upper), r(se_adjusted))
        
        if `i' == 1 {
            matrix results = row_`i'
        }
        else {
            matrix results = results \ row_`i'
        }
        
        local row_names "`row_names' `outcome'"
    }
    else {
        di as error "Failed to process `outcome' - skipping"
    }
    
    local i = `i' + 1
}

* Label matrix
matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
matrix rownames results = `row_names'

* Display summary table
di as text _newline(2) "{hline 100}"
di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
di as text "{hline 100}"
matlist results, format(%9.4f) border(all) rowtitle("Outcome")
outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_lab.tex", mat(results) nobox format(%9.3f) replace



**# Measurement error: We repeat the sensitivity analysis, but we add the
* residuals from a regression of score in BAS-verbal on teachers' 
* assessment of vocabulary competence. 



qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)

* Redefine the globals to avoid the "i." operator (fwl cannot handle it for some reason).

tab c_region_birth, gen(region_)
tab c_sexbirth, gen(sex_)
global set1 c_delivery_age_mum c_delivery_age_mum_2 region_1 region_2 region_3 region_4 region_5 region_6 region_7 region_8 region_9 region_10 c_mumsmoked c_abnormal_gest sex_2 c_birthweight 


global set4 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband1 c_incband2 c_incband3 c_incband4 c_incband5 c_incband6 c_incband7 c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10

 
* Measurment error proxy

reg bas_def_total_sd j065 j066 j067 j068 j069 pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7 
predict measurement_error, res

********************************************************************************
* Loop over outcomes and store sensitivity results
********************************************************************************

* Define list of outcome variables
local outcomes luckvar_k075 luckvar_k080 luckvar_k094 ///
	skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a /// 
	skill_skill6a skill_skill7a  skill_skill9a ///
	degree noqualifications ///
	good_at_k036 good_at_k037 good_at_k038 good_at_k039 good_at_k040 



* Initialize empty matrix to store results
matrix results = J(1, 10, .)
matrix drop results  // Drop to start fresh
local row_names ""

* Loop over each outcome
local i = 1
foreach outcome of local outcomes {
    
    di as text _newline(2) "{hline 80}"
    di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
    di as text "{hline 80}"
    
    * Run sensitivity analysis
    capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
        xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7 measurement_error) ///
        psvars($set1 $set2 $set3) ///
        method("aipw") ///
        keyvars(pc1 perserverance_2_sd) ///
        benchmark(pc1)
    
    * Store results if successful
    if _rc == 0 {
        matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
                          r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
                          r(tau_lower), r(tau_upper), r(se_adjusted))
        
        if `i' == 1 {
            matrix results = row_`i'
        }
        else {
            matrix results = results \ row_`i'
        }
        
        local row_names "`row_names' `outcome'"
    }
    else {
        di as error "Failed to process `outcome' - skipping"
    }
    
    local i = `i' + 1
}

* Label matrix
matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
matrix rownames results = `row_names'

* Display summary table
di as text _newline(2) "{hline 100}"
di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
di as text "{hline 100}"
matlist results, format(%9.4f) border(all) rowtitle("Outcome")
outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_edu_MesErr.tex", mat(results) nobox format(%9.3f) replace


**# Measurement error panel


merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
gen logwage = log(net_pay_week_95)

local outcomes employedFT unemployed isManagerSupervisor hours logwage


* Initialize empty matrix to store results
matrix results = J(1, 10, .)
matrix drop results  // Drop to start fresh
local row_names ""

* Loop over each outcome
local i = 1
foreach outcome of local outcomes {
    
    di as text _newline(2) "{hline 80}"
    di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
    di as text "{hline 80}"
    
    * Run sensitivity analysis
    capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
        xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7 measurement_error) ///
        psvars($set1 $set2 $set3) ///
        method("aipw") ///
        keyvars(pc1 perserverance_2_sd) ///
        benchmark(pc1)
    
    * Store results if successful
    if _rc == 0 {
        matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
                          r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
                          r(tau_lower), r(tau_upper), r(se_adjusted))
        
        if `i' == 1 {
            matrix results = row_`i'
        }
        else {
            matrix results = results \ row_`i'
        }
        
        local row_names "`row_names' `outcome'"
    }
    else {
        di as error "Failed to process `outcome' - skipping"
    }
    
    local i = `i' + 1
}

* Label matrix
matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
matrix rownames results = `row_names'

* Display summary table
di as text _newline(2) "{hline 100}"
di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
di as text "{hline 100}"
matlist results, format(%9.4f) border(all) rowtitle("Outcome")
outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_lab_MesErr.tex", mat(results) nobox format(%9.3f) replace





**# Heterogeneity



qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)

* Redefine the globals to avoid the "i." operator (fwl cannot handle it for some reason).

tab c_region_birth, gen(region_)
tab c_sexbirth, gen(sex_)
global set1 c_delivery_age_mum c_delivery_age_mum_2 region_1 region_2 region_3 region_4 region_5 region_6 region_7 region_8 region_9 region_10 c_mumsmoked c_abnormal_gest sex_2 c_birthweight 


global set4 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband1 c_incband2 c_incband3 c_incband4 c_incband5 c_incband6 c_incband7 c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10

forvalues kappa = 0/1{
	
	
	preserve
	keep if shitlife== `kappa'

	* Define list of outcome variables
	local outcomes luckvar_k075 luckvar_k080 luckvar_k094 ///
		skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a /// 
		skill_skill6a skill_skill7a  skill_skill9a ///
		degree noqualifications ///
		good_at_k036 good_at_k037 good_at_k038 good_at_k039 good_at_k040 



		* Initialize empty matrix to store results
		matrix results = J(1, 10, .)
		matrix drop results  // Drop to start fresh
		local row_names ""

		* Loop over each outcome
		local i = 1
		foreach outcome of local outcomes {
			
			di as text _newline(2) "{hline 80}"
			di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
			di as text "{hline 80}"
			
			* Run sensitivity analysis
			capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
				xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
				psvars($set1 $set2 $set3) ///
				method("aipw") ///
				keyvars(pc1 perserverance_2_sd) ///
				benchmark(pc1)
			
			* Store results if successful
			if _rc == 0 {
				matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
								  r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
								  r(tau_lower), r(tau_upper), r(se_adjusted))
				
				if `i' == 1 {
					matrix results = row_`i'
				}
				else {
					matrix results = results \ row_`i'
				}
				
				local row_names "`row_names' `outcome'"
			}
			else {
				di as error "Failed to process `outcome' - skipping"
			}
			
			local i = `i' + 1
		}

		* Label matrix
		matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
		matrix rownames results = `row_names'

		* Display summary table
		di as text _newline(2) "{hline 100}"
		di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
		di as text "{hline 100}"
		matlist results, format(%9.4f) border(all) rowtitle("Outcome")
		outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_educ_Adv`kappa'.tex", mat(results) nobox format(%9.3f) replace
	
	restore

	** Sensitivity Panel 

	preserve
	
	merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
	gen logwage = log(net_pay_week_95)

	local outcomes employedFT unemployed isManagerSupervisor hours logwage

	keep if shitlife==`kappa'
	
		* Initialize empty matrix to store results
		matrix results = J(1, 10, .)
		matrix drop results  // Drop to start fresh
		local row_names ""

		* Loop over each outcome
		local i = 1
		foreach outcome of local outcomes {
			
			di as text _newline(2) "{hline 80}"
			di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
			di as text "{hline 80}"
			
			* Run sensitivity analysis
			capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
				xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
				psvars($set1 $set2 $set3) ///
				method("aipw") ///
				keyvars(pc1 perserverance_2_sd) ///
				benchmark(pc1)
			
			* Store results if successful
			if _rc == 0 {
				matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
								  r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
								  r(tau_lower), r(tau_upper), r(se_adjusted))
				
				if `i' == 1 {
					matrix results = row_`i'
				}
				else {
					matrix results = results \ row_`i'
				}
				
				local row_names "`row_names' `outcome'"
			}
			else {
				di as error "Failed to process `outcome' - skipping"
			}
			
			local i = `i' + 1
		}

		* Label matrix
		matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
		matrix rownames results = `row_names'

		* Display summary table
		di as text _newline(2) "{hline 100}"
		di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
		di as text "{hline 100}"
		matlist results, format(%9.4f) border(all) rowtitle("Outcome")
		outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_lab_Adv_`kappa'.tex", mat(results) nobox format(%9.3f) replace
	restore

	
}


** Gender 


qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)

* Redefine the globals to avoid the "i." operator (fwl cannot handle it for some reason).

tab c_region_birth, gen(region_)
tab c_sexbirth, gen(sex_)
global set1 c_delivery_age_mum c_delivery_age_mum_2 region_1 region_2 region_3 region_4 region_5 region_6 region_7 region_8 region_9 region_10 c_mumsmoked c_abnormal_gest sex_2 c_birthweight 


global set4 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband1 c_incband2 c_incband3 c_incband4 c_incband5 c_incband6 c_incband7 c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10

forvalues kappa = 1/2{
	
	
	preserve
	keep if c_sexbirth== `kappa'

	* Define list of outcome variables
	local outcomes luckvar_k075 luckvar_k080 luckvar_k094 ///
		skill_skill1a skill_skill2a skill_skill3a skill_skill4a skill_skill5a /// 
		skill_skill6a skill_skill7a  skill_skill9a ///
		degree noqualifications ///
		good_at_k036 good_at_k037 good_at_k038 good_at_k039 good_at_k040 



		* Initialize empty matrix to store results
		matrix results = J(1, 10, .)
		matrix drop results  // Drop to start fresh
		local row_names ""

		* Loop over each outcome
		local i = 1
		foreach outcome of local outcomes {
			
			di as text _newline(2) "{hline 80}"
			di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
			di as text "{hline 80}"
			
			* Run sensitivity analysis
			capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
				xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
				psvars($set1 $set2 $set3) ///
				method("aipw") ///
				keyvars(pc1 perserverance_2_sd) ///
				benchmark(pc1)
			
			* Store results if successful
			if _rc == 0 {
				matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
								  r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
								  r(tau_lower), r(tau_upper), r(se_adjusted))
				
				if `i' == 1 {
					matrix results = row_`i'
				}
				else {
					matrix results = results \ row_`i'
				}
				
				local row_names "`row_names' `outcome'"
			}
			else {
				di as error "Failed to process `outcome' - skipping"
			}
			
			local i = `i' + 1
		}

		* Label matrix
		matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
		matrix rownames results = `row_names'

		* Display summary table
		di as text _newline(2) "{hline 100}"
		di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
		di as text "{hline 100}"
		matlist results, format(%9.4f) border(all) rowtitle("Outcome")
		outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_educ_Gen_`kappa'.tex", mat(results) nobox format(%9.3f) replace
	
	restore

	** Sensitivity Panel 

	preserve
	
	merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
	gen logwage = log(net_pay_week_95)

	local outcomes employedFT unemployed isManagerSupervisor hours logwage

	keep if c_sexbirth==`kappa'
	
		* Initialize empty matrix to store results
		matrix results = J(1, 10, .)
		matrix drop results  // Drop to start fresh
		local row_names ""

		* Loop over each outcome
		local i = 1
		foreach outcome of local outcomes {
			
			di as text _newline(2) "{hline 80}"
			di as result "Processing: `outcome' (`i' of `: word count `outcomes'')"
			di as text "{hline 80}"
			
			* Run sensitivity analysis
			capture noisily senseweighter_hetero `outcome' perserverance_2_sd, ///
				xvar(pc1 $set1 $set2 $set3 $set4 $set5 $set6 $set7) ///
				psvars($set1 $set2 $set3) ///
				method("aipw") ///
				keyvars(pc1 perserverance_2_sd) ///
				benchmark(pc1)
			
			* Store results if successful
			if _rc == 0 {
				matrix row_`i' = (r(tau_res), r(se_tau), r(bias), r(R2ydx), r(RV), ///
								  r(bound_R2_D_Z), r(bound_R2_Y_Z), ///
								  r(tau_lower), r(tau_upper), r(se_adjusted))
				
				if `i' == 1 {
					matrix results = row_`i'
				}
				else {
					matrix results = results \ row_`i'
				}
				
				local row_names "`row_names' `outcome'"
			}
			else {
				di as error "Failed to process `outcome' - skipping"
			}
			
			local i = `i' + 1
		}

		* Label matrix
		matrix colnames results = Coef SE Bias R2yd_x RV R2dz_x R2yz_dx Lower_CI Upper_CI Adj_SE
		matrix rownames results = `row_names'

		* Display summary table
		di as text _newline(2) "{hline 100}"
		di as result "SUMMARY: SENSITIVITY ANALYSIS RESULTS"
		di as text "{hline 100}"
		matlist results, format(%9.4f) border(all) rowtitle("Outcome")
		outtable using "/Users/user/Dropbox/Econometrics/perserverance/sensitivity_lab_Gen_`kappa'.tex", mat(results) nobox format(%9.3f) replace
	restore

	
}

 
