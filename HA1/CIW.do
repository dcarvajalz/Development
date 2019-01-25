cd "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1"
set more off
gl data "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data"
gl data2013 "/Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Data/UGA_2013_UNPS_v01_M_STATA8"
gl output "//Users/danielcarvajalz/Dropbox/5th Term/Development/HA/HA1/Output"

****************************** Consumption ****************************

use "$data2013/UNPS 2013-14 Consumption Aggregate.dta", clear
gen C = cpexp30*12
rename HHID hh
sort hh
quietly by hh:  gen dup = cond(_N==1,0,_n) 
drop if dup > 0 
keep hh district_code hsize C urban ea region regurb wgt_X 
replace hh = subinstr(hh, "H", "", .)
replace hh = subinstr(hh, "-", "", .)
destring hh, gen(hh1)
drop hh
rename hh1 hh
save "$output/Consumption.dta", replace

****************************** Income ****************************

********************** Agricultural net production

* Agricultural net production
use "$data2013/AGSEC5A.dta", clear
drop if cropID==. & a5aq5_2==2 & a5aq6a==. & a5aq6d==. & a5aq7a==. & a5aq7d==. ///
& a5aq8==. & a5aq10==. & a5aq16==.
replace a5aq6d = a5aq7d if  a5aq6b==a5aq7b & a5aq6c==a5aq7c & a5aq6d!=a5aq7d // Same conversion weights for buying and selling 
keep HHID parcelID plotID cropID a5aq6a a5aq6b a5aq6c a5aq6d a5aq7a ///
a5aq7b a5aq7c a5aq7d a5aq8 a5aq10 a5aq16 a5aq5_2

** Harvested crop
gen Q1_t = .
replace Q1_t = a5aq6a*a5aq6d if a5aq6d!=. // Q * Conversion factor
replace Q1_t = a5aq6a if a5aq6c==1

** Harvested crop sold
gen Q2_t = .
replace Q2_t = 0 if a5aq7a==0
replace Q2_t = 0 if a5aq7a==.
replace Q2_t = a5aq7a*a5aq7d if a5aq7d!=.
replace Q2_t = a5aq7a if a5aq7c==1
gen dif = Q1_t-Q2_t
replace Q2_t = Q1_t if dif<0
bysort HHID cropID: egen Q2 = sum(Q2_t) // total Q sold by HH and crop
bysort HHID cropID: egen Q1 = sum(Q1_t) //total Q by HH and crop  
bysort HHID cropID: egen RQ_t = sum(a5aq8) //revenue by HH and crop
gen P_t = (RQ_t/Q2)
bysort cropID: egen P = mean(P_t)
 
** Value of retained output
gen Q1_2 = P*(Q1 - Q2)
replace RQ_t = 0 if RQ_t==. 
replace Q1_2 = 0 if Q1_2==. 

** Costs
*** Transportation costs
bysort HHID cropID: gen Q3 = a5aq10 
replace Q3 = 0 if Q3==.
collapse (mean) Q1_2 Q3 RQ_t, by(HHID cropID) 
collapse (sum) Q1_2 Q3 RQ_t, by(HHID) 
save "$output/ANP1.dta", replace

*** Rent-in land
use "$data2013/AGSEC2B.dta", clear
bysort HHID parcelID: gen Q4 = a2bq9
replace Q4 = 0 if Q4==. 
collapse (sum) Q4, by(HHID)
merge 1:1 HHID using "$output/ANP1.dta"
drop _merge
save "$output/ANP1.dta", replace

*** Hired labor
use "$data2013/AGSEC3A.dta", clear
bysort HHID parcelID plotID: gen Q5 = a3aq36
replace Q5 = 0 if Q5==. 

*** Pesticides and fertilizers
bysort HHID parcelID plotID: gen Q6 = a3aq8 
bysort HHID parcelID plotID: gen Q7 = a3aq18 
bysort HHID parcelID plotID: gen Q8 = a3aq27 
replace Q6 = 0 if Q6==. 
replace Q7 = 0 if Q7==. 
replace Q8 = 0 if Q8==. 
collapse (sum) Q*, by(HHID)
merge 1:1 HHID using "$output/ANP1.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/ANP1.dta", replace
 
*** Seeds
use "$data2013/AGSEC4A.dta", clear
bysort HHID parcelID plotID cropID: gen Q9 = a4aq15 
replace Q9 = 0 if Q9==. 
collapse (sum) Q*, by(HHID)
merge 1:1 HHID using "$output/ANP1.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/ANP1.dta", replace


gen Q = Q1_2+RQ_t-Q3-Q4-Q5-Q6-Q7-Q8-Q9
save "$output/ANP1.dta", replace


* Agricultural net production (second visit)

use "$data2013/AGSEC5B.dta", clear
drop if cropID==. & a5bq5_2==2 & a5bq6a==. & a5bq6d==. & a5bq7a==. & a5bq7d==. ///
& a5bq8==. & a5bq10==. & a5bq16==.
replace a5bq6d = a5bq7d if a5bq6b==a5bq7b & a5bq6c==a5bq7c & a5bq6d!=a5bq7d // Same conversion weights for buying and selling 
keep HHID parcelID plotID cropID a5bq6a a5bq6b a5bq6c a5bq6d a5bq16 a5bq7a ///
a5bq7b a5bq7c a5bq7d a5bq8 a5bq10 a5bq5_2

** Harvested crop
gen Q1_t2 = .
replace Q1_t2 = a5bq6a*a5bq6d if a5bq6d!=. // Q * Conversion factor
replace Q1_t2 = a5bq6a if a5bq6c==1
  
** Harvested crop sold 
gen Q2_t2 = .
replace Q2_t2 = a5bq7a*a5bq7d if a5bq7d!=.
replace Q2_t2 = 0 if a5bq7a==0 
replace Q2_t2 = 0 if a5bq7a==. 
replace Q2_t2 = a5bq7a if a5bq7c==1
gen dif = (Q1_t2 - Q2_t2)
replace Q2_t = Q1_t2 if dif<0
bysort HHID cropID: egen Q2_2 = sum(Q2_t2) // total Q sold by HH and crop
bysort HHID cropID: egen Q1_2 = sum(Q1_t2) // total Q by HH and crop
bysort HHID cropID: egen RQ_t2 = sum(a5bq8) // revenue by HH and crop
gen P_t2 = (RQ_t2/Q2_2)
bysort cropID: egen P_2 = mean(P_t2)
 
** Value of retained output
gen Q1_2_2 = P_2*(Q1_2 - Q2_2) 
replace Q1_2_2 = 0 if Q1_2_2 == . 
replace RQ_t2 = 0 if RQ_t2 == . 

** Costs
*** Transportation costs
bysort HHID cropID: gen Q3_2 = a5bq10
replace Q3_2 = 0 if Q3_2 == .
collapse (mean) Q1_2_2 Q3_2 RQ_t2, by(HHID cropID) 
collapse (sum) Q1_2_2 Q3_2 RQ_t2, by(HHID) 
save "$output/ANP2.dta", replace

*** Hired labor
use "$data2013/AGSEC3B.dta", clear
bysort HHID parcelID plotID: gen Q5_2 = a3bq36
replace Q5_2 = 0 if Q5_2 == . 
 
*** Pesticides and fertilizers
bysort HHID parcelID plotID: gen Q6_2 = a3bq8 
bysort HHID parcelID plotID: gen Q7_2 = a3bq18 
bysort HHID parcelID plotID: gen Q8_2 = a3bq27 
replace Q6_2 = 0 if Q6_2==. 
replace Q7_2 = 0 if Q7_2==. 
replace Q8_2 = 0 if Q8_2==. 
collapse (sum) Q*, by(HHID)
merge 1:1 HHID using "$output/ANP2.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/ANP2.dta", replace
 
*** Seeds
use "$data2013/AGSEC4B.dta", clear
bysort HHID parcelID plotID cropID: gen Q9_2 = a4bq15 
replace Q9_2 = 0 if Q9_2==. 
collapse (sum) Q*, by(HHID)
merge 1:1 HHID using "$output/ANP2.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}

gen Q_2 = Q1_2_2+RQ_t2-Q3_2-Q5_2-Q6_2-Q7_2-Q8_2-Q9_2
save "$output/ANP2.dta", replace

* Livestock

** Cattle
use "$data2013/AGSEC6A.dta", clear
keep if a6aq2!=2 & a6aq3a!=0 & a6aq3a!=. 
gen L1 = a6aq14a*a6aq14b if a6aq14a!=. & a6aq14a!=0  & a6aq14b!=. & a6aq14b!=0 //revenues
gen L2 = . 
replace L2 = a6aq5c if a6aq5c>0 & a6aq5c!=. //labor cost
replace L2 = 0 if L2 ==. 
collapse (sum) L1 (mean) L2, by(HHID)
save "$output/L.dta", replace
 
** Small animals
use "$data2013/AGSEC6B.dta", clear
keep if a6bq2 != 2 & a6bq3a != 0 & a6bq3a != . // we keep only those who own 
gen L3 = .
replace L3 = a6bq14a*a6bq14b if a6bq14a !=. & a6bq14a != 0  & a6bq14b !=. & a6bq14b !=0 //revenues = quantity * revenue by unit
replace L3 = 0 if L3==.
gen L4 = . 
replace L4 = a6bq5c if a6bq5c >0 & a6bq5c != . //cost labor
replace L4 = 0 if L4 ==.  
collapse (sum) L3 (mean) L4, by(HHID)
merge 1:1 HHID using "$output/L.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/L.dta", replace
 
** Rabbits
use "$data2013/AGSEC6C.dta", clear
keep if a6cq2 != 2 & a6cq3a != 0 & a6cq3a != . // we keep only those who own 
gen L5 = .
replace L5 = a6cq14a*a6cq14b if a6cq14a !=. & a6cq14a != 0  & a6cq14b !=. & a6cq14b !=0 //revenues = quantity * revenue by unit
replace L5 = 0 if L5==.
gen L6 = . 
replace L6 = a6cq5c if a6cq5c >0 & a6cq5c != . 
replace L6 = 0 if L6 ==. //cost labor
collapse (sum) L5 (mean) L6, by(HHID)
merge 1:1 HHID using "$output/L.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/L.dta", replace

** Other costs
use "$data2013/AGSEC7.dta", clear
keep if a7aq1 == 1 // we keep only those who own or raise cattle
bysort HHID: egen L7 = sum(a7bq2e) 
bysort HHID: egen L8 = sum(a7bq3f) 
bysort HHID: egen L9 = sum(a7bq5d) 
bysort HHID: egen L10 = sum(a7bq6c) 
bysort HHID: egen L11 = sum(a7bq7c) 
bysort HHID: egen L12 = sum(a7bq8c) 
gen L13 = L7 + L8 + L9 + L10 + L11 + L12
collapse (mean) L13, by(HHID)
merge 1:1 HHID using "$output/L.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}


gen L = L1 + L3 + L5 - L2 - L4 - L6 - L13
save "$output/L.dta", replace

* Livestock product

** Meat
use "$data2013/AGSEC8A.dta", clear
gen PM_t = .
replace PM_t = a8aq5/a8aq3 if a8aq1!=0 & a8aq5!=0 & a8aq5!=. & a8aq3!=0 & a8aq3!=.
bysort AGroup_ID: egen PM = mean(PM_t)
gen LM =. 
replace LM = PM*((a8aq1*a8aq2)-a8aq3)+a8aq5 if a8aq5!=. 
replace LM = PM*((a8aq1*a8aq2)-a8aq3) if a8aq5==.
replace LM = a8aq5 if ((a8aq1*a8aq2)-a8aq3)==0 & a8aq5!=.
replace LM = 0 if LM==.
collapse (sum) LM, by(HHID)
save "$output/LP.dta", replace

** Milk
use "$data2013/AGSEC8B.dta", clear
gen Milk_day = a8bq1*a8bq3 // animals milked * average milk production per day
replace Milk_day = 0 if Milk_day==.
replace a8bq5_1 = Milk_day if a8bq5_1>Milk_day & a8bq5_1!=0 & a8bq5_1!=. //Sales = production for those who report more sales than production
replace a8bq7 = 0 if a8bq6==0 | a8bq6==.
replace a8bq7 = a8bq6 if a8bq7>a8bq6 // sales = amount converted for those who report selling more than they convert per day
replace a8bq5_1 = Milk_day if Milk_day<a8bq5 & a8bq5_1!=0 & a8bq9!=0 & a8bq9!=. & a8bq6==0
replace a8bq5 = 0 if Milk_day<a8bq5 & a8bq5_1!=0 & a8bq9!=0 & a8bq9!=.
replace a8bq5_1 = a8bq5_1*30 *a8bq2
replace a8bq7 = a8bq7*30*a8bq2
replace a8bq7 = 0 if a8bq5_1==a8bq1*a8bq2*30*a8bq3 // day to year
gen PMilk_t = .
replace PMilk_t = a8bq9/(a8bq5_1+a8bq7) if a8bq1!=0 & a8bq5_1!=0 & a8bq5_1!=. ///
& a8bq9!=0 & a8bq9!=.| a8bq1!=0 & a8bq6!=0 & a8bq6!=. & a8bq7!=0 & ///
a8bq7!=. & a8bq9!=0 & a8bq9!=.
bysort AGroup_ID: egen PMilk = mean(PMilk_t)
replace a8bq2 = a8bq2*30 // year
gen Milk = a8bq1*a8bq2*a8bq3 //l of milk = quantity * days * avg day production
replace Milk = 0 if Milk==. 
replace a8bq7 = 0 if a8bq7 ==. 
replace a8bq5_1 = 0 if a8bq5_1==. 
gen dif = (Milk-(a8bq5_1+a8bq7)) //quantity milk production - quantity sold
gen LMilk = .
replace LMilk = PMilk*(Milk-(a8bq5_1+a8bq7))
replace LMilk = PMilk*(Milk-(a8bq5_1+a8bq7))+a8bq9 if a8bq9!=. 
replace LMilk = a8bq9 if Milk-(a8bq5_1+a8bq7)==0 &  a8bq9!=. 
collapse (sum) LMilk, by(HHID)
merge 1:1 HHID using "$output/LP.dta"
drop _merge
replace LMilk = 0 if LMilk==.
save "$output/LP.dta", replace

** Eggs
use "$data2013/AGSEC8C.dta", clear 
replace a8cq2 = a8cq2*4 //Q eggs per year
replace a8cq3 = a8cq3*4 //Q sold per year
replace a8cq5 = a8cq5*4 //revenue per year
replace a8cq3=a8cq2 if a8cq3>a8cq2
gen PE_t = a8cq5/a8cq3 if a8cq1!=0 & a8cq1!=0 & a8cq2!=. & a8cq2!=0
bysort AGroup_ID: egen PE = mean(PE_t)
gen LE = .
replace LE = PE*(a8cq2-a8cq3)
replace LE = PE*(a8cq2-a8cq3)+a8cq5 if a8cq5!=. 
replace LE = 0 if LE==. 
collapse (sum) LE, by(HHID)
merge 1:1 HHID using "$output/LP.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
save "$output/LP.dta", replace
 
** Dung
use "$data2013/AGSEC11.dta", clear
gen LD = a11q1c + a11q5 //revenues: dung + ploughing
collapse (sum) LD, by(HHID)
merge 1:1 HHID using "$output/LP.dta"
drop _merge
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
gen LP = LM+LMilk+LE+LD
save "$output/LP.dta", replace
 
* Renting in agricultural equipment and capital
use "$data2013/AGSEC10.dta", clear
rename a10q8 rent //value rent
collapse (sum) rent, by(HHID)
merge 1:1 HHID using "$output/ANP1.dta"
drop _merge
merge 1:1 HHID using "$output/ANP2.dta"
drop _merge 
merge 1:1 HHID using "$output/L.dta"
drop _merge
merge 1:1 HHID using "$output/LP.dta"
drop _merge 
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}

gen ANP_total = Q+Q_2+L+LP-rent
save "$output/ANP.dta", replace
 
********************** Labor market income
use "$data2013/GSEC8_1.dta", clear

** Main job
gen LMI1 = .
replace LMI1 = (h8q31a+h8q31b)*56 if h8q31c==1 // Assuming 8hrs every day of the week
replace LMI1 = (h8q31a+h8q31b)*4 if h8q31c==2 // Assuming 30 days per month
replace LMI1 = (h8q31a+h8q31b) if h8q31c==3 // Assuming 4 weeks per month
replace LMI1 = (h8q31a+h8q31b)/4 if h8q31c==4 // Assuming 4 weeks per month
replace LMI1 = LMI1*h8q30b*h8q30a // last year

** Secondary Job 
gen LMI2 = .
replace LMI2 = (h8q45a+h8q45b)*56 if h8q45c==1 // Assuming 8hrs every day of the week
replace LMI2 = (h8q45a+h8q45b)*4 if h8q45c==2 // Assuming 30 days per month
replace LMI2 = (h8q45a+h8q45b) if h8q45c==3 // Assuming 4 weeks per month
replace LMI2 = (h8q45a+h8q45b)/4 if h8q45c==4 // Assuming 4 weeks per month
replace LMI2 = LMI2*h8q44b*h8q44 // last year
 
gen LMI = LMI1+LMI2
replace LMI = 0 if LMI==. 
collapse (sum) LMI, by(HHID)
save "$output/LMI.dta", replace

********************** Business income

use "$data2013/gsec12.dta", clear
rename hhid HHID
gen BI1 = .
replace BI1 = h12q13 // Monthly gross revenue
replace BI1= 0 if BI1==. 
gen BI2 = .
replace BI2 = h12q15 // Monthly labor cost
replace BI2= 0 if BI2==. 
gen BI3 = .
replace BI3 = h12q16+h12q17 // Monthly expenditure raw materials + others
replace BI3= 0 if BI3==. 
gen BI = (BI1-BI2-BI3)*h12q12
replace BI = 0 if BI==. 
collapse (sum) BI, by(HHID)
save "$output/BI.dta", replace
 
********************** Other income sources

use "$data2013/GSEC11A.dta", clear
gen OIS = h11q5+h11q6
replace OIS = 0 if OIS==.
collapse (sum) OIS, by(HHID)
save "$output/OIS.dta", replace

********************** Transfers

use "$data2013/GSEC15B.dta", clear
gen Transfer = .
replace Transfer = h15bq10*h15bq11
replace Transfer = 0 if Transfer==.
collapse (sum) Transfer, by(HHID)

********************** Total Income

merge 1:1 HHID using "$output/LMI.dta"
drop _merge
merge 1:1 HHID using "$output/BI.dta"
drop _merge
merge 1:1 HHID using "$output/OIS.dta"
drop _merge 
replace HHID = subinstr(HHID, "H", "", .)
replace HHID = subinstr(HHID, "-", "", .)
destring HHID, gen(hh)
drop HHID
rename hh HHID
merge 1:1 HHID using "$output/ANP.dta"
drop _merge 
rename HHID hh
save "$output/income.dta", replace
foreach var of varlist _all {
	replace `var' = 0 if `var'==.
}
gen I = ANP_total+LMI+BI+OIS
save "$output/Income.dta", replace


****************************** Wealth ****************************

use "$data2013/GSEC14A.dta", clear

* House assets
gen  H_t = . 
replace H_t = h14q5 if h14q3==1
replace H_t = 0 if h14q3==2
replace H_t = 0 if H_t==.
bysort HHID: egen H = sum(H_t)
collapse (mean) H, by(HHID)
rename HHID hh
save "$output/Wealth.dta", replace 

* Agricultural Land Value
use "$data2013/AGSEC2B.dta", clear
keep if a2bq9!=.
gen P_r = a2bq9/a2bq5 // Rental price / acres
drop if P_r==.
collapse (mean) P_r // Rental price = 68673.928
 
** Ownership data 
use "$data2013/AGSEC2A.dta", clear
gen AGL = . 
replace AGL = 69778.398 * 10 * a2aq5 if a2aq5!=0 
collapse (sum) AGL, by(hh) 
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
replace AGL = 0 if AGL==. 
save "$output/Wealth.dta", replace 
 
* Agricultural equipment and structure capital
use "$data2013/AGSEC10.dta", clear
gen  AK_t = . 
replace AK_t = a10q2 if  a10q1>0
replace AK_t= 0 if AK_t==. 
bysort HHID: egen AK = sum(AK_t)
collapse (mean) AK, by(hh)
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
replace AK = 0 if AK == .
save "$output/Wealth.dta", replace 

* Livestock capital
 
** Cattle 
use "$data2013/AGSEC6A.dta", clear
bysort LiveStockID: egen P_buy = mean(a6aq13b) // Mean buying price
replace P_buy = 0 if P_buy==. 
bysort LiveStockID: egen P_sell = mean(a6aq14b) // Mean selling price
replace P_sell = 0 if P_sell==. 
gen P_cat = (P_buy + P_sell)/2 if P_sell!=0 & P_buy!=0
replace P_cat = P_sell if P_sell>0 & P_buy==0
replace P_cat = P_buy if P_sell==0 & P_buy>0
replace P_cat = 0 if P_sell==0 & P_buy==0
gen cat_t = . 
replace cat_t = P_cat*a6aq3a if a6aq3a>0 // only for those who currently own it
replace cat_t = 0 if cat_t==. // missing values do not have nor own
bysort HHID: egen cattle = sum(cat_t)
collapse (mean) cattle, by(hh)
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
replace cattle = 0 if cattle == .
save "$output/Wealth.dta", replace 

** Small animals 
use "$data2013/AGSEC6B.dta", clear
bysort ALiveStock_Small_ID: egen P_buy = mean(a6bq13b) // Mean buying price
replace P_buy = 0 if P_buy==. 
bysort ALiveStock_Small_ID: egen P_sell = mean(a6bq14b) // Mean selling price
replace P_sell = 0 if P_sell==. 
gen P_small = (P_buy + P_sell)/2 if P_sell!=0 & P_buy!=0
replace P_small = P_sell if P_sell>0 & P_buy==0
replace P_small = P_buy if P_sell==0 & P_buy>0
replace P_small = 0 if P_sell==0 & P_buy==0
gen  small_t = . 
replace small_t = a6bq3a*P_small if a6bq3a>0 // only for those who currently own it
replace small_t = 0 if small_t==. //missing values do not have nor own
bysort HHID: egen small = sum(small_t)
collapse (mean) small, by(hh)
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
replace small = 0 if small==.
save "$output/Wealth.dta", replace 

 ** Poultry 
use "$data2013/AGSEC6C.dta", clear
bysort APCode: egen P_buy = mean(a6cq13b) // Mean buying price
replace P_buy = 0 if P_buy==. 
bysort APCode: egen P_sell = mean(a6cq14b) // Mean selling price
replace P_sell = 0 if P_sell==. 
gen P_poultry = (P_buy + P_sell)/2 if P_sell!=0 & P_buy!=0
replace P_poultry = P_sell if P_sell>0 & P_buy==0
replace P_poultry = P_buy if P_sell==0 & P_buy>0
replace P_poultry = 0 if P_sell==0 & P_buy==0
gen  poultry_t = . 
replace poultry_t = a6cq3a*P_poultry  if a6cq3a>0 
replace poultry_t= 0 if poultry_t== .
bysort HHID: egen poultry = sum(poultry_t)
collapse (mean) poultry, by(hh)
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
replace poultry = 0 if poultry==.

* Total Wealth

gen W = H+AGL+AK+cattle+small+poultry
gen HHID = hh
keep hh W HHID
replace hh = subinstr(hh, "H", "", .)
replace hh = subinstr(hh, "-", "", .)
destring hh, gen(hh1)
drop hh
rename hh1 hh
save "$output/Wealth.dta", replace 



****************************** CIW ****************************



use "$output/Consumption.dta", clear
merge 1:1 hh using "$output/Income.dta"
drop _merge
merge 1:1 hh using "$output/Wealth.dta"
drop _merge
save "$output/CIW.dta", replace

* Merge with HH roster to get hh head, age, education
** Gender, Age
merge m:m HHID using "$data2013/GSEC2.dta"
drop _merge
keep if h2q4==1
rename h2q8 age 
rename h2q3 gender
keep  HHID  PID district_code urban ea region regurb C I W ///
wgt_X hsize  age gender h2q4

** Education
merge 1:1 HHID PID using "$data2013/GSEC4.dta"
drop _merge
keep if h2q4 == 1 
rename h4q7 education
keep HHID district_code urban ea region regurb C I W wgt_X ///
hsize age gender education
replace HHID = subinstr(HHID, "H", "", .)
replace HHID = subinstr(HHID, "-", "", .)
destring HHID, gen(hhid)
drop HHID
rename hhid HHID
drop if C ==.
drop if I == 0
drop if I <0
drop if C > W+I 
save "$output/CIW.dta", replace
 
