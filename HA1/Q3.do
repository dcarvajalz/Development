cd "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1"
set more off
gl data "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data"
gl data2013 "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data/UGA_2013_UNPS_v01_M_STATA8"
gl output "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Output"

use "$output/laborsupply.dta", clear
collapse (mean) intensive C I W, by (district_code)
foreach var in C I W {
	replace `var' = log(`var')
}
 
* Part 1 (Level)
 
preserve
scatter intensive C, ytitle("Mean hours worked", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("Log Mean C", size(medlarge)) graphregion(color(white))
graph export "$output/graph311.png", replace
restore
preserve
drop if district_code =="413" // outlier
scatter intensive I, ytitle("Mean hours worked", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("Log Mean I", size(medlarge)) graphregion(color(white))
graph export "$output/graph312.png", replace
restore
preserve
scatter intensive W, ytitle("Mean hours worked", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("LogMean W", size(medlarge)) graphregion(color(white))
graph export "$output/graph313.png", replace
restore
 
* Part 2 (Inequality)
 
use "$output/laborsupply.dta", clear

* Create logs
foreach var in intensive {
	gen log_`var' = log(`var')
	gen log_`var'_mean = .
	gen v_`var' = .
}

* Variance
foreach var in intensive {
	 sum log_`var' [w=wgt_X]
	 replace log_`var'_mean = r(mean)
	 replace v_`var'= (log_`var'-log_`var'_mean)^2
}
collapse (mean) v_intensive I W C, by(district_code)
foreach var in I W C {
	replace `var' = log(`var')
}
preserve
scatter v_intensive C, ytitle("Mean Var log", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("Log Mean C", size(medlarge)) graphregion(color(white))
graph export "$output/graph321.png", replace
restore
 
preserve
drop if district_code =="413" // drop the outlier
scatter v_intensive I, ytitle("Mean Var log", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("Log Mean I", size(medlarge)) graphregion(color(white))
graph export "$output/graph322.png", replace
restore
 
preserve
scatter v_intensive W, ytitle("Mean Var log", size(medlarge)) ylabel(, labsize(medlarge) noticks nogrid) /// 
xtitle("Log Mean W", size(medlarge)) graphregion(color(white))
graph export "$output/graph323.png", replace
restore
	
	
* Part 3 (Covariance)

use "$output/laborsupply.dta", clear

** Districts
levelsof district_code, local(district)
bysort district: correlate I intensive
egen corr_I = corr(I intensive), by(district)
egen corr_W = corr(W intensive), by(district)
egen corr_C = corr(C intensive), by(district)
collapse (mean) corr_I corr_W corr_C I W C, by(district_code)
foreach var in C I W {
	replace `var' = log(`var')
}
preserve
scatter corr_I I || lfit corr_I I, ytitle("Corr hours worked", size(medlarge)) ///
ylabel(, labsize(medlarge) noticks nogrid) xtitle("Log Mean I", size(medlarge)) ///
graphregion(color(white)) legend(off)
graph export "$output/graph331.png", replace
restore
 
preserve
drop if district_code =="413" // drop the outlier
scatter corr_W I || lfit corr_I I, ytitle("Corr hours worked", size(medlarge)) ///
ylabel(, labsize(medlarge) noticks nogrid) xtitle("Log Mean W", size(medlarge)) ///
graphregion(color(white)) legend(off)
graph export "$output/graph332.png", replace
restore

preserve
scatter corr_C I || lfit corr_I I, ytitle("Corr hours worked", size(medlarge)) ///
ylabel(, labsize(medlarge) noticks nogrid) xtitle("Log Mean C", size(medlarge)) ///
graphregion(color(white)) legend(off)
graph export "$output/graph333.png", replace
restore 
 
