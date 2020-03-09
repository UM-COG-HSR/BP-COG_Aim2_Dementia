/* Based on Tutorial: Survival Estimation for Cox Regression
Models with Time-Varying Coecients Using SAS
and R 
*/

data Surv2;
input obs id age female time0 time1 death;
cards;
 1 1 20 0  0  1 1
 2 2 21 1  0  1 0
 8 2 21 1  1  4 0
 3 3 19 0  0  1 0
 9 3 19 0  1  7 1
 4 4 22 1  0  1 0
10 4 22 1  1  7 0
16 4 22 1  7 10 1
 5 5 20 0  0  1 0
11 5 20 0  1  7 0
17 5 20 0  7 10 0
23 5 20 0 10 12 0
 6 6 24 1  0  1 0
12 6 24 1  1  7 0
18 6 24 1  7 10 0
24 6 24 1 10 13 1
;
run;

proc format;
 value gender 0="F"
              1="M";
run;

/*-- analytical dataset with formatted female variable --*/
data surv2a;
 set surv2;
  lt_age = age*log(time1); /* lt_age is time dependent covariate */
  format female gender.;
run;

title "Surv2a";
proc print data=surv2a;
run;

data covs;
 age = 0;
 lt_age = 0;
 female = 0;
 format female gender.; /* !!!! */
run;

Title "surv2a";
proc print data=surv2a;
run;

Title "baseline statement";
proc phreg data = SURV2a;
class female (ref = "F");
model (time0, time1) * death(0) = female age lt_age;
baseline out = outset survival = survival covariates = covs lower=S_lower upper=S_upper
    / method = emp;
run;

Title "covs";
proc print data = covs;
run;

Title "outset";
proc print data = outset;
run;


© 2020 GitHub, Inc.
