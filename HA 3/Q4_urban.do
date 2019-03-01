cd "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA 3"
clear all
set more off
	

* QUESTION 1 
use dataUGA.dta, clear
 
keep hh year wave lnc age age_sq ethnic female urban lninctotal_trans familysize
** Correct year
bysort year hh: gen dup = _N
replace year = 2010 if wave=="2010-2011" & year==2011 & dup==2
drop dup
bysort year hh: gen dup = _N // 218 repeated years
replace year = 2009 if wave=="2009-2010" & year==2010 & dup==2
drop dup wave 	

keep if urban==1
	 
** Residuals for consumption
reg lnc age age_sq familysize i.ethnic i.female i.year
predict res
rename res res_c
 
** Residuals for income
reg lninctotal_trans age age_sq familysize i.ethnic i.female i.year
predict res
rename res res_i
rename lninctotal_trans income
 
** Aggregate consumption
bysort year: egen ag_c = sum(lnc)
 
** Panel
keep res_c res_i ag_c hh year income
xtset hh year
reshape wide res_c res_i ag_c income, i(hh) j(year)
forvalues y = 10(1)14 {
	egen ag_c20`y'_t = mean(ag_c20`y')
	drop ag_c20`y'
	rename ag_c20`y'_t ag_c20`y'
}
egen ag_c2009_t = mean(ag_c2009)
drop ag_c2009
rename ag_c2009_t ag_c2009
reshape long res_c res_i ag_c income, i(hh)
rename _j year
 
** Interpolate the missing obs and drop hh with just one year
bysort hh: ipolate res_c year, generate(res_ci) epolate
bysort hh: ipolate res_i year, generate(res_ii) epolate  
bysort hh: ipolate income year, generate(income_i) epolate  
gen ones = 1
replace ones = 0 if res_ci ==.
egen numy = sum(ones), by(hh)
drop if numy <= 1
drop res_c res_i ones numy

sort hh year
egen id = group(hh) 
 
** Regressions with random coefficients for each hh
gen beta = .
gen phi = .
quietly forvalues i = 1(1)615 {
	reg d.res_ci d.res_ii d.ag_c if id==`i', nocons
	replace beta = _b[d.res_ii] if id==`i'
	replace phi = _b[d.ag_c] if id==`i'
}

* QUESTION 3: Regressions with average coefficients
reg d.res_ci d.res_ii d.ag_c, nocons
display _b[d.res_ii]
display _b[d.ag_c]
 
** Histogram of betas and phis
preserve
collapse beta phi, by(hh)
drop if beta > 2
drop if beta < -2
sum beta, detail
histogram beta, xtitle ("Beta") graphregion(color(white)) bcolor(blue)
graph export "hist_beta_urban.png", replace
restore

preserve
collapse beta phi, by(hh)
drop if phi > 0.00002
drop if phi < -0.00002
sum phi, detail
histogram phi, xtitle ("Phi") graphregion(color(white)) bcolor(green)
graph export "hist_phi_urban.png", replace
restore

* QUESTION 2
 
** Average hh income
gen ones = 1
replace ones = 0 if income_i ==.
egen numy = sum(ones), by(hh)
drop if numy <= 1
drop ones numy income
collapse (mean) income_i beta, by(hh)
 
** Mean and median of beta in five income groups
sort income_i
gen nobs = _N
gen nhh = _n 
gen inc_g = 0
replace inc_g = 1 if nhh<=122
replace inc_g = 2 if nhh>122 & nhh<=244
replace inc_g = 3 if nhh>244 & nhh<=365
replace inc_g = 4 if nhh>365 & nhh<=487
replace inc_g = 5 if nhh>487 & nhh<=609
 
** Mean and median betas
forvalues i = 1(1)5 {
	sum beta if inc_g==`i', detail
}
drop nhh
 
** Mean and median of income in five beta groups
sort beta
gen nhh = _n 
gen beta_g = 0
replace beta_g = 1 if nhh<=122
replace beta_g = 2 if nhh>122 & nhh<=244
replace beta_g = 3 if nhh>244 & nhh<=365
replace beta_g = 4 if nhh>365 & nhh<=487
replace beta_g = 5 if nhh>487 & nhh<=609
 
** Mean and median betas
forvalues i = 1(1)5 {
	sum income_i if beta_g==`i', detail
}
