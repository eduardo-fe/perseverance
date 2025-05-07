


program drop run_models_weights
program define run_models_weights, rclass
    syntax [if] ,KEYvars(string) Binary(string) Continuous(string) COVariates(string) 

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
		predict xb
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
run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary( employedFT unemployed) ///
    continuous(hours net_pay_week_95) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 )

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour.tex", mat(mResults) nobox format(%9.3f) replace

run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary( employedFT unemployed) ///
    continuous(hours net_pay_week_95) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)

matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour_set7.tex", mat(mResults) nobox format(%9.3f) replace

restore



* ------------------------------------------------------------------------------
* Next, we study work histories with data from the Activity file -to age 46.
* This includes number of job spells, duration of the spells, date of the
* spells, as well as economic activity -FT, PT, unemployed, etc.
* For manager, use the variable JEMPST.
* THese data are useful to understand the nubmer of managerial roles,
* the number of unemployed spells, the duration of those spells. 
* you cannot use it to see the likelihood of beig employed or unemployed. 
* ------------------------------------------------------------------------------


merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-6943-stata/stata/stata13/bcs70_activity_histories_eul.dta", keepusing( CMSEX NUMACTS ACTNUM LSTSWCMC JSTCMC JENDCMC JDUR JFTPT JACTIV JSTMTH JSTYR JENDMTH JENDYR HISTLEN JYLEFT JYLEFTX JEMPST JEMPST)

gen manager = JEMPST==5 | JEMPST==5
replace manager =. if JEMPST <0


gen unemployed_hist = JACTIV== 11  
replace unemployed_hist =. if JACTIV<0.

gen unemployed_dur = unemployed_hist * JDUR
replace unemployed_dur=. if JDUR <0.  


encode bcsid, gen(person_id)
xtset person_id ACTNUM
bysort person_id: egen min_spell = min(JDUR)
bysort person_id: egen max_spell = max(JDUR)


* ------------------------------------------------------------------------------
* Estimate models.
* ------------------------------------------------------------------------------


run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary(manager unemployed_hist) ///
    continuous(unemployed_dur min_spell max_spell) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6)


matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/history.tex", mat(mResults) nobox format(%9.3f) replace 

* ------------------------------------------------------------------------------
* Re-estimate models using the additional principal components.
* ------------------------------------------------------------------------------


run_models_weights , ///
	keyvars(perserverance_2_sd pc1) ///
    binary(manager unemployed_hist) ///
    continuous(unemployed_dur min_spell max_spell) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)


matrix mResults = r(results)
matlist mResults , format(%9.3f)

outtable using "/Users/user/Dropbox/Econometrics/perserverance/history_set7.tex", mat(mResults) nobox format(%9.3f) replace 
	


	
* ------------------------------------------------------------------------------
* Now estudy sub-sets (adverse early life)
* ------------------------------------------------------------------------------


	
forvalues kappa = 0/1{
	
	

	* ------------------------------------------------------------------------------
	* Load the covariates (run these files in that order 3, 1, 2)
	* ------------------------------------------------------------------------------

	qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
	qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
	qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"


	* ------------------------------------------------------------------------------
	* Estimate models 
	* ------------------------------------------------------------------------------


	run_models_weights if shitlife== `kappa' , ///
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

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/skills_sl`kappa'", mat(mResults) nobox format(%9.3f) replace

	 

	* ------------------------------------------------------------------------------
	* Re-estimate models using the additional principal components.
	* ------------------------------------------------------------------------------

	run_models_weights if shitlife== `kappa' , ///
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

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/skills_set7_sl`kappa'.tex", mat(mResults) nobox format(%9.3f) replace



	* ------------------------------------------------------------------------------
	* Labour market outcomes from the panel.
	* ------------------------------------------------------------------------------




	preserve

	merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta"
	run_models_weights if shitlife== `kappa' , ///
		keyvars(perserverance_2_sd pc1) ///
		binary( employedFT unemployed) ///
		continuous(hours net_pay_week_95) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6 )

	matrix mResults = r(results)
	matlist mResults , format(%9.3f)

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour_sl`kappa'.tex", mat(mResults) nobox format(%9.3f) replace

	run_models_weights if shitlife== `kappa' , ///
		keyvars(perserverance_2_sd pc1) ///
		binary( employedFT unemployed) ///
		continuous(hours net_pay_week_95) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)

	matrix mResults = r(results)
	matlist mResults , format(%9.3f)

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/labour_set7_sl`kappa'.tex", mat(mResults) nobox format(%9.3f) replace

	restore



	* ------------------------------------------------------------------------------
	* Next, we study work histories with data from the Activity file -to age 46.
	* This includes number of job spells, duration of the spells, date of the
	* spells, as well as economic activity -FT, PT, unemployed, etc.
	* For manager, use the variable JEMPST.
	* THese data are useful to understand the nubmer of managerial roles,
	* the number of unemployed spells, the duration of those spells. 
	* you cannot use it to see the likelihood of beig employed or unemployed. 
	* ------------------------------------------------------------------------------


	merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-6943-stata/stata/stata13/bcs70_activity_histories_eul.dta", keepusing( CMSEX NUMACTS ACTNUM LSTSWCMC JSTCMC JENDCMC JDUR JFTPT JACTIV JSTMTH JSTYR JENDMTH JENDYR HISTLEN JYLEFT JYLEFTX JEMPST JEMPST)

	gen manager = JEMPST==5 | JEMPST==5
	replace manager =. if JEMPST <0


	gen unemployed_hist = JACTIV== 11  
	replace unemployed_hist =. if JACTIV<0.

	gen unemployed_dur = unemployed_hist * JDUR
	replace unemployed_dur=. if JDUR <0.  


	encode bcsid, gen(person_id)
	xtset person_id ACTNUM
	bysort person_id: egen min_spell = min(JDUR)
	bysort person_id: egen max_spell = max(JDUR)


	* ------------------------------------------------------------------------------
	* Estimate models.
	* ------------------------------------------------------------------------------


	run_models_weights if shitlife== `kappa' , ///
		keyvars(perserverance_2_sd pc1) ///
		binary(manager unemployed_hist) ///
		continuous(unemployed_dur min_spell max_spell) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6)


	matrix mResults = r(results)
	matlist mResults , format(%9.3f)

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/history_sl`kappa'.tex", mat(mResults) nobox format(%9.3f) replace 

	* ------------------------------------------------------------------------------
	* Re-estimate models using the additional principal components.
	* ------------------------------------------------------------------------------


	run_models_weights if shitlife== `kappa' , ///
		keyvars(perserverance_2_sd pc1) ///
		binary(manager unemployed_hist) ///
		continuous(unemployed_dur min_spell max_spell) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6 $set7)


	matrix mResults = r(results)
	matlist mResults , format(%9.3f)

	outtable using "/Users/user/Dropbox/Econometrics/perserverance/history_set7_sl`kappa'.tex", mat(mResults) nobox format(%9.3f) replace 
		

}
		
	
	
	
	
	


