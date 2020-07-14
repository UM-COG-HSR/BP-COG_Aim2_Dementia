/***********************************************
* Program: xxxxxx.SAS                          *
* Folder: x:\xxxx\xxxxxx\xxxxx                 *
* Author: xxxxxxx                              *
* Created: xx/xx/xx                            *
* Summary: xxxxxx xxx xxxxx xxxxxxxxxxx xxx    *
* Revisions: xxxxx xxx xxxx xxxxxxxxxxxx       *
***********************************************/
/* Please provide info about this repository */
%let repo_name = BP-COG_Aim2_Dementia; /* Repository name on GitHub */
%let repo_maintainer = Nick Tilton;
%let repo_description = Produce survival function for a set of baseline covariates and randomly selected test cases;

%put repo_name := &repo_name;  /* Github Repository name */

/*****---- SAS Setup starts here -----*****/ 

/* Define global macro variable names */
%global 
 SAS_work_dir
 computer_name
 sas_batchmode   /* Y/N */
 sas_progname    
 sas_fullname
;

/* Store path to SAS working directory in `SAS_work_dir` macro variable */
/* Based on `https://communities.sas.com/t5/SAS-Communities-Library/Find-current-directory-path/ta-p/485785` */
/* Note the location of the file*/
filename setupc  "C:/Users/Public/SAS_work_directory.sas";
%include setupc;
%let SAS_work_dir =  %SAS_work_directory;
%put SAS_work_dir := &SAS_work_dir;   
%put sysuserid    := &sysuserid;   /* User id */

/*--- Load repository assets ----*/ 
filename fx "&SAS_work_dir/_load_repo_assets.inc";
%include fx;

%computer_name;  /* Stores computer name in `computer_name` global macro variable */
%put computer_name := &computer_name;
%our_sas_session_info;
/***** SAS Setup ends here *****/


%let BPCOG_path = S:\Intmed_Rsrch2\GenMed\Restricted\BP COG;
%let BPCOG_DEMpath = &BPCOG_path.\Aim 2\Dementia Model;
%let OutputPath = &BPCOG_DEMpath.\Output;

libname storemac "&BPCOG_DEMpath.\SAS Compiled Macros";
libname fmts "&BPCOG_path.\Aim 1\Data Management\Data\formats";
libname dem "&BPCOG_DEMpath.\SAS Data";
libname sr "&BPCOG_DEMpath.\SAS Data\SAS Results";
options fmtsearch=(fmts) nofmterr mstored sasmstore=storemac;

%macro today_YYMMDD();
%let z=0;
%let y2=%sysfunc(today(),year2.);
%let m2=%sysfunc(today(),month2.);
%let d2=%sysfunc(today(),day2.);
%if %eval(&m2)<=9 %then %let m2 = &z&m2;
%if %eval(&d2)<=9 %then %let d2 = &z&d2;
%let ymd = &y2&m2&d2;
&ymd;
%mend;

%let ymd = %today_YYMMDD();

/* produce survival function for a set of given baseline covariates */

data blcovs;
input age female0 educ racebpcog gcp_bl gcp_slope;
format female0 yna. educ educ. racebpcog racebpcog.;
datalines;
20 1 5 2 50 -1
;
run; 

proc phreg data = dem.coxreg_currgcp; *outest = ESTS;
ods rtf file="&OutputPath.\Cox_&ymd";
class female0 (ref="Yes") racebpcog (ref="Non-Hispanic White") 
educ (ref="College graduate or more (Technical School Certificate, Associate Degree, Bachelor's Degree, Graduate or Professional School)");	   	
model (t_start,t_end)*deminc(0) = gcp_bl gcp_slope age female0 educ racebpcog;
baseline covariates=blcovs out=sr.bl_surfunc_&ymd survival=S lower=S_lower upper=S_upper;
ods output parameterestimates = sr.coxparams_&ymd;
run;

/* produce survival function for two randomly selected participants from each of ARIC, CHS and FOS */

proc sort data=dem.coxreg_currgcp out=currgcp; by newid t_end; run;

data currgcp;
set currgcp;
by newid;
if first.newid;
run;

  PROC SURVEYSELECT DATA=currgcp OUT=aricsamp METHOD=SRS
  SAMPSIZE=2 SEED=1234567;
  where studyname='aric';
  RUN;

  PROC SURVEYSELECT DATA=currgcp OUT=chssamp METHOD=SRS
  SAMPSIZE=2 SEED=1234567;
  where studyname='chs';
  RUN;

  PROC SURVEYSELECT DATA=currgcp OUT=fossamp METHOD=SRS
  SAMPSIZE=2 SEED=1234567;
  where studyname='fos';
  RUN;

data allsamp;
set aricsamp chssamp fossamp;
format female0 yna. educ educ. racebpcog racebpcog.;
keep newid age female0 educ racebpcog gcp_bl gcp_slope;
run;

	proc phreg data = dem.coxreg_currgcp;
	class female0 (ref="Yes") racebpcog (ref="Non-Hispanic White") 
	educ (ref="College graduate or more (Technical School Certificate, Associate Degree, Bachelor's Degree, Graduate or Professional School)");	   	
	model (t_start,t_end)*deminc(0) = gcp_bl gcp_slope age female0 educ racebpcog;
	baseline covariates=allsamp out=sr.bl_surfunc_testcases_&ymd survival=S lower=S_lower upper=S_upper cumhaz=cumhaz lowercumhaz=lower_haz uppercumhaz=upper_haz;
	assess var=(gcp_bl gcp_slope age female0 educ racebpcog);
	run;

data sr.bl_testcases_1yr_&ymd;
set sr.Bl_surfunc_testcases_&ymd;
if t_end > 1 and t_end < 1.07;
run;
