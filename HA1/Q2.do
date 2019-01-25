cd "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1"
set more off
gl data "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data"
gl data2013 "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data/UGA_2013_UNPS_v01_M_STATA8"
gl output "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Output"

* Labor Supply
use "$data2013/GSEC8_1.dta", clear
keep HHID h8q30a h8q30b h8q31a h8q31b h8q31c h8q36a h8q36b h8q36c h8q36d h8q36e ///
h8q36f h8q36g h8q43 h8q44 h8q44b h8q45a h8q45b h8q45c
foreach var of varlist h8q30a h8q30b h8q31a h8q31b h8q31c h8q36a h8q36b h8q36c ///
	h8q36d h8q36e h8q36f h8q36g h8q43 h8q44 h8q44b h8q45a h8q45b h8q45c {
	replace `var' = 0 if `var'==.
}

** Main job
gen hours_week = h8q36a+h8q36b+h8q36c+h8q36d+h8q36e+h8q36f+h8q36g // h per week
gen hours_year = hours_week*h8q30a*h8q30b // h per year

** Second job
gen hours_week2 = h8q43
gen hours_year2 = hours_week2*h8q44*h8q44b
collapse (sum) hours_week hours_week2 hours_year hours_year2, by(HHID)
 
replace HHID = subinstr(HHID, "H", "", .)
replace HHID = subinstr(HHID, "-", "", .)
destring HHID, gen(hh)
drop HHID
rename hh HHID
 
merge 1:1 HHID using "$output/CIW.dta"
drop _merge
 
** Intensive and extensive margin of labor supply
keep if age>15 & age<70
gen intensive = hours_year + hours_year2
bysort urban: egen emp=sum(wgt_X) if intensive>0
bysort urban: egen total=sum(wgt_X) 
gen extensive = emp/total 
 
save "$output/laborsupply.dta", replace

* Part 1  
 
** Average labor supply for rural and urban areas
mean intensive[pw=wgt_X]
mean intensive[pw=wgt_X] if urban==0  
mean intensive[pw=wgt_X] if urban==1  

mean extensive[pw=wgt_X]
mean extensive[pw=wgt_X] if urban==0  
mean extensive[pw=wgt_X] if urban==1  

** Histogram
preserve
twoway (histogram intensive if urban==0, fcolor(none) lcolor(red)) (histogram intensive if urban==1, ///
fcolor(none) lcolor(blue)),legend(order(1 "Rural" 2 "Urban")) xtitle(Intensive margin (hours worked)) /// 
graphregion(color(white)) 
graph export "$output/graph211.png", replace 
restore

** Variance of logs
*** Create logs
foreach var in intensive {
	gen log_`var'=log(`var')
	gen log_`var'_mean=.
	gen v_`var'=.
}
*** Variance
foreach var in intensive {
	sum log_`var' [w=wgt_X]
	replace log_`var'_mean = r(mean)
	replace v_`var'=(log_`var'-log_`var'_mean)^2
}
mean v_intensive [pw=wgt_X]
mean v_intensive [pw=wgt_X] if urban== 0 
mean v_intensive [pw=wgt_X] if urban== 1 
 
** CIW level, inequality, and covariances over the lifecycle
preserve	  
collapse (mean) intensive v_intensive, by(age)
graph twoway (line intensive age, fcolor(none) lcolor(blue)), xtitle("Age") xlabel(15(10)70, ///
labsize(medlarge) noticks grid angle(0)) ylabel(, labsize(medlarge) nogrid) ///
graphregion(color(white))
graph export "$output/graph212.png", replace
graph twoway (line v_intensive age, fcolor(none) lcolor(blue)), xtitle("Age") xlabel(15(10)70, ///
labsize(medlarge) noticks grid angle(0)) ylabel(, labsize(medlarge) nogrid) ///
graphregion(color(white))
graph export "$output/graph213.png", replace
restore
  
*Part 2
 
* Gender
use "$output/laborsupply.dta", clear

** Average labor supply for rural and urban areas
forvalues i=1/2 {
	mean intensive[pw=wgt_X] if gender==`i'
	mean intensive[pw=wgt_X] if urban==0 & gender==`i'
	mean intensive[pw=wgt_X] if urban==1 & gender==`i'
}	 

** Histogram
preserve
gen log_intensive = log(intensive) 
twoway (histogram log_intensive if urban==0 & gender==1, fcolor(none) lcolor(blue)) ///
(histogram log_intensive if urban==1 & gender==1, fcolor(none) lcolor(red)), ///
legend(order(1 "Rural" 2 "Urban")) xtitle("Int margin (hours worked)") graphregion(color(white)) 
graph export "$output/graph221.png", replace 
twoway (histogram log_intensive if urban==0 & gender==2, fcolor(none) lcolor(blue)) ///
(histogram log_intensive if urban==1 & gender==2, fcolor(none) lcolor(red)), ///
legend(order(1 "Rural" 2 "Urban")) xtitle("Int margin (hours worked)") graphregion(color(white)) 
graph export "$output/graph222.png", replace 
restore

** Variance of logs
*** Create logs
forvalues i = 1/2 {
	foreach var in intensive {
		gen log_`var'_`i'=log(`var') if gender == `i'
		gen log_`var'_mean_`i'=. if gender == `i'
		gen v_`var'_`i'=. if gender == `i'
	}
}
*** Variance
forvalues i = 1/2 {
	foreach var in intensive {
		sum log_`var'_`i' [w=wgt_X] if gender == `i'
		replace log_`var'_mean_`i' = r(mean) if gender == `i'
		replace v_`var'_`i'=(log_`var'_`i'-log_`var'_mean_`i')^2 if gender == `i'
	}
}
forvalues i = 1/2 {
	mean v_intensive_`i' [pw=wgt_X]
	mean v_intensive_`i' [pw=wgt_X] if urban== 0 
	mean v_intensive_`i' [pw=wgt_X] if urban== 1 
}

** CIW level, inequality, and covariances over the lifecycle
preserve
collapse (mean) intensive v_intensive_*, by(age gender) 	 
graph twoway (line intensive age if gender==1, fcolor(none) lcolor(blue))(line intensive age if gender==2, ///
fcolor(none) lcolor(red)), xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ///
ylabel(, labsize(medlarge) nogrid) ytitle("Int margin (hours worked)", size(medlarge)) ///
legend(order(1 "Men" 2 "Women")) graphregion(color(white))
graph export "$output/graph223.png", replace
graph twoway (line v_intensive_1 age, fcolor(none) lcolor(blue))(line v_intensive_2 age if gender==2, ///
fcolor(none) lcolor(red)), xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ///
ylabel(, labsize(medlarge) nogrid) ytitle("Int margin (var of log)", size(medlarge)) ///
legend(order(1 "Men" 2 "Women")) graphregion(color(white))
graph export "$output/graph224.png", replace
restore

* Education Groups
use "$output/laborsupply.dta", clear
drop if education==. | education==99
rename education educ
gen education =.
replace education = 1 if educ < 17 // less than P.7
replace education = 2 if educ>=17 & educ < 34 //Primary education but less than high school
replace education = 3 if educ>=34 // Secondary and more

** Average labor supply for rural and urban areas
forvalues i=1/3 {
	mean intensive[pw=wgt_X] if education==`i'
	mean intensive[pw=wgt_X] if urban==0 & education==`i'
	mean intensive[pw=wgt_X] if urban==1 & education==`i'  
}

** Histogram
preserve
gen log_intensive=log(intensive)
twoway (histogram log_intensive if urban==0 & education==1, fcolor(none) lcolor(blue)) ///
(histogram log_intensive if urban==1 & education==1, fcolor(none) lcolor(red)), ///
legend(order(1 "Rural" 2 "Urban")) xtitle("Int margin (hours worked)") graphregion(color(white)) 
graph export "$output/graph225.png", replace 
twoway (histogram log_intensive if urban==0 & education==2, fcolor(none) lcolor(blue)) ///
(histogram log_intensive if urban==1 & education==2, fcolor(none) lcolor(red)), ///
legend(order(1 "Rural" 2 "Urban")) xtitle("Int margin (hours worked)") graphregion(color(white)) 
graph export "$output/graph226.png", replace 
twoway (histogram log_intensive if urban==0 & education==3, fcolor(none) lcolor(blue)) ///
(histogram log_intensive if urban==1 & education==3, fcolor(none) lcolor(red)), ///
legend(order(1 "Rural" 2 "Urban")) xtitle("Int margin (hours worked)") graphregion(color(white)) 
graph export "$output/graph227.png", replace 
restore

** Variance of logs
*** Creates logs
forvalues i = 1/3 {
	foreach var in intensive {
		gen log_`var'_`i'=log(`var') if education == `i'
		gen log_`var'_mean_`i'=. if education == `i'
		gen v_`var'_`i'=. if education == `i'
	}
}
*** Variance
forvalues i = 1/3 {
	foreach var in intensive {
		sum log_`var'_`i' [w=wgt_X] if education == `i'
		replace log_`var'_mean_`i' = r(mean) if education == `i'
		replace v_`var'_`i'=(log_`var'_`i'-log_`var'_mean_`i')^2 if education == `i'
	}
}
forvalues i = 1/3 {
	mean v_intensive_`i' [pw=wgt_X] if education==`i'
	mean v_intensive_`i' [pw=wgt_X] if urban== 0 & education==`i'
	mean v_intensive_`i' [pw=wgt_X] if urban== 1 & education==`i'
}

** CIW level, inequality, and covariances over the lifecycle
preserve  
collapse (mean) intensive v_intensive_*, by(age education)	 	 
graph twoway (line intensive age if education==1, fcolor(none) lcolor(blue))(line intensive age if education==2, /// 
fcolor(none) lcolor(red))(line intensive age if education==3, fcolor(none) lcolor(purple)), ///
xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ylabel(, labsize(medlarge) nogrid) ///
ytitle("Int margin (hours worked)", size(medlarge))legend(order(1 "Less than Primary" 2 "Less than High School" 3 "High School or more")) ///
graphregion(color(white))
graph export "$output/graph228.png", replace	 
graph twoway (line v_intensive_1 age if education==1, fcolor(none) lcolor(blue))(line v_intensive_2 age if education==2, /// 
fcolor(none) lcolor(red))(line v_intensive_3 age if education==3, fcolor(none) lcolor(purple)), ///
xtitle("Age") xlabel(15(10)70, labsize(medlarge) noticks grid angle(0)) ylabel(, labsize(medlarge) nogrid) ///
ytitle("Int margin (hours worked)", size(medlarge))legend(order(1 "Less than Primary" 2 "Less than High School" 3 "High School or more")) ///
graphregion(color(white))
graph export "$output/graph229.png", replace	 
restore
