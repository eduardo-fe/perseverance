*-------------------------------------------------------------------------------
* Create mini panel dataset for employment, hours and income
*-------------------------------------------------------------------------------

*--------------------------------------
* Wave 1996, Study 3833
*--------------------------------------

use b960312 b960318 b960277 bcsid empstat using ///
    "/Users/user/Documents/datasets/britCohortStudy70/UKDA-3833-stata/stata/stata13/bcs96x.dta", clear

gen wave = 1996

* Rename variables
rename (b960312 b960318 b960277) (netpay_raw netpay_period hours)

* Clean negative or invalid values
replace netpay_raw     = . if netpay_raw < 0
replace netpay_period  = . if netpay_period < 0 | netpay_period >= 6
replace hours          = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw * hours       if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw * 5       if netpay_period == 2  // Daily
replace net_pay_week = netpay_raw           if netpay_period == 3  // Weekly
replace net_pay_week = (netpay_raw * 12)/52 if netpay_period == 4  // Monthly
replace net_pay_week = netpay_raw / 52      if netpay_period == 5  // Annual

* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)

gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed

* Harmonise employment variable

replace empstat = empstat +1 if empstat >6

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_1996.dta", replace


*--------------------------------------
* Wave 2000, Study 5558
*--------------------------------------

use cnetpay cnetprd chours1 bcsid empstat using ///
"/Users/user/Documents/datasets/britCohortStudy70/UKDA-5558-stata/stata/stata13_se/bcs2000.dta", clear

gen wave = 2000

* Rename variables
rename (cnetpay cnetprd chours1) (netpay_raw netpay_period hours)

* Clean negative or invalid values
replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 3  // Quarterly
replace net_pay_week = (netpay_raw *12)/52 if netpay_period == 4  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 52    if netpay_period == 5  // Annual


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2000.dta", replace


*--------------------------------------
* Wave 2004, Study 5585
*--------------------------------------

use bd7ecact bcsid b7cnetpd b7cnetpy b7chour* b7otimny using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-5585-stata/stata/stata13_se/bcs_2004_followup.dta", clear

gen wave = 2004

forvalues i=1/3{ 
	replace b7chour`i' =. if b7chour`i'<0
	}

egen hours = rowtotal (b7chour1 b7chour2 b7chour3), missing

rename (bd7ecact b7cnetpy b7cnetpd ) (empstat netpay_raw netpay_period)

replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 3  // Quarterly
replace net_pay_week = (netpay_raw*12)/52        if netpay_period == 4  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 52    if netpay_period == 5  // Annual


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2004.dta", replace


*--------------------------------------
* Wave 2012, Study 7473 
*--------------------------------------

use  bcsid B9NETA B9NETP B9OTIMNY B9CHOUR1 B9CHOUR2 B9CHOUR3 B9CHOUR4 using"/Users/user/Documents/datasets/britCohortStudy70/UKDA-7473-stata/stata/stata13/bcs70_2012_flatfile.dta", clear

gen wave = 2012

merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-7473-stata/stata/stata13/bcs70_2012_derived.dta", keepusing(BD9ECACT)

drop _merge



forvalues i=1/3{ 
	replace B9CHOUR`i' =. if B9CHOUR`i'<0
	}

egen hours = rowtotal (B9CHOUR1 B9CHOUR2 B9CHOUR3), missing

rename (BD9ECACT B9NETA B9NETP ) (empstat netpay_raw netpay_period)

replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 4  // Quarterly
replace net_pay_week = (netpay_raw*12)/52        if netpay_period == 5  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 3    if netpay_period == 3  // THree weeks


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2012.dta", replace


*--------------------------------------
* Wave 2016, Study 8547 
*--------------------------------------


use B10CHOUR1 B10NETW BD10ECACT bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-8547-stata/stata/stata13/bcs_age46_main.dta", clear

gen wave = 2016 

rename (B10CHOUR1 B10NETW BD10ECACT )(hours net_pay_week empstat)

replace net_pay_week =. if net_pay_week <0
* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed
* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2016.dta", replace



*--------------------------------------
* Wave 2021, Study 9347 
*--------------------------------------

use b11econact2 b11netw b11chours* bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/bcs11_age51_main.dta", clear

gen wave = 2021


foreach i of numlist 1 3 4 5{ 
	replace b11chours`i' =. if b11chours`i'<0
	}

egen hours = rowtotal (b11chours1 b11chours3), missing

rename (b11econact2 b11netw ) (empstat net_pay_week)
replace net_pay_week =. if net_pay_week <0
* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed
* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2021.dta", replace


append using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_1996.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2000.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2004.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2012.dta" "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2016.dta"


save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta", replace
