cd "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1"
set more off
gl data "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data"
gl data2013 "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data/UGA_2013_UNPS_v01_M_STATA8"
gl output "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Output"

use "$output/CIW.dta", clear
* Convert to 2013 USD 
foreach var in C W I {
	replace `var' = `var'/3696.24
}
 
* Part 1
 
** Average CIW for rural and urban areas
 
mean C[pw=wgt_X] if urban==0 // C in rural 
mean C[pw=wgt_X] if urban==1 // C in urban  
mean I [pw=wgt_X] if urban==0 // I in rural  
mean I [pw=wgt_X] if urban==1 // I in urban   
mean W [pw=wgt_X] if urban==0 // W in rural 
mean W [pw=wgt_X] if urban==1 // W in urban 

* Part 2 
 
** Histogram
* Trim
_pctile C, nq(100)
drop if C >r(r99) 
_pctile W, nq(100)
drop if W >r(r99) 
_pctile I, nq(100)
drop if I >r(r99) 

twoway (histogram C if urban==0, fcolor(none) lcolor(red))(histogram C if urban==1, ///
fcolor(none) lcolor(blue)), legend(order(1 "Rural" 2 "Urban")) xtitle(C) graphregion(color(white)) 
graph export "$output/graph121.png", replace 
 
twoway (histogram I if urban==0, fcolor(none) lcolor(red))(histogram I if urban==1, ///
 fcolor(none) lcolor(blue)),legend(order(1 "Rural" 2 "Urban")) xtitle(I) graphregion(color(white)) 
graph export "$output/graph122.png", replace 

twoway (histogram W if urban==0, fcolor(none) lcolor(red))(histogram W if urban==1, ///
fcolor(none) lcolor(blue)),legend(order(1 "Rural" 2 "Urban")) xtitle(W) graphregion(color(white)) 
graph export "$output/graph123.png", replace 
 
** Variance of logs
* Create logs
foreach var in C W I {
	gen log_`var' = log(`var')
	gen log_`var'_mean = .
	gen v_`var' = .
}
* Variance
foreach var in C W I {
	sum log_`var' [w=wgt_X]
	replace log_`var'_mean = r(mean)
	replace v_`var'= (log_`var'-log_`var'_mean)^2
	mean v_`var' [pw=wgt_X] if urban==0
	mean v_`var' [pw=wgt_X] if urban==1
}
mean v_C [pw=wgt_X] if urban== 0 
mean v_C [pw=wgt_X] if urban== 1 
mean v_I [pw=wgt_X] if urban== 0 
mean v_I [pw=wgt_X] if urban== 1 
mean v_W [pw=wgt_X] if urban== 0 
mean v_W [pw=wgt_X] if urban== 1 
drop log_* v_*
 
* Part 3
correlate C I W 
correlate C I W if urban==0 
correlate C I W if urban==1 
 
* Part 4
** CIW level over the lifecycle
keep if age<70 
gen I_lc = .
gen W_lc = .
gen C_lc = .
foreach var in I W C {
	forvalues i=15(1)105 {
		sum `var' [w=wgt_X] if age==`i'
		replace `var'_lc = r(mean) if age==`i'
	}
}

preserve
collapse (mean) I_lc W_lc C_lc, by(age) 	 
graph twoway (line I_lc age, fcolor(none) lcolor(purple))(line W_lc age, fcolor(none) /// 
lcolor(blue))(line C_lc age, fcolor(none) lcolor(red)),legend(order(1 "I" 2 "W" 3 "C")) ///
xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ///
ylabel(, labsize(medlarge) nogrid) graphregion(color(white))
graph export "$output/graph141.png", replace 
restore

** CIW inequality over the lifecycle

foreach var in C W I {
	gen log_`var'=log(`var')
	gen log_`var'_mean=.
	gen v_`var'=.
	gen v_`var'_all =.
}
foreach var in C W I {
	forvalues i=15(1)70{
		sum log_`var' [w=wgt_X] if age==`i'
		replace log_`var'_mean = r(mean) if age==`i'
		replace v_`var' = (log_`var' - log_`var'_mean)^2 if age==`i'
		sum v_`var' [w=wgt_X] if age==`i'
		replace v_`var'_all = r(mean) if age==`i'
	}
}
foreach var in C W I {
	forvalues i=15(1)70 {
		sum v_`var' [w=wgt_X] if age==`i'
		replace v_`var'_all = r(mean) if age==`i'
	}
}
preserve
collapse (mean) v_*, by(age) 
graph twoway (line v_I_all age, fcolor(none) lcolor(purple))(line v_W_all age, fcolor(none) ///
lcolor(blue))(line v_C_all age, fcolor(none) lcolor(red)),legend(order(1 "Var log(I)" 2 "Var log(W)" 3 "Var log(C)")) ///
xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ///
ylabel(, labsize(medlarge) nogrid) graphregion(color(white))
graph export "$output/graph142.png", replace
restore 
 
* Part 5

use "$output/CIW.dta", clear
* Convert to 2013 USD 
foreach var in C W I {
	replace `var' = `var'/3696.24
}
sort I
* Consumption
preserve
sort I
gen cum_I_C = C[1]
replace cum_I_C = C[_n]+cum_I_C[_n-1] if _n>1
sum cum_I_C, d
scalar list
restore
* Wealth
preserve
sort I
gen cum_I_W = W[1]
replace cum_I_W = W[_n]+cum_I_W[_n-1] if _n>1
sum cum_I_W, d
scalar list
restore
 
