/***********************************************
* Program: xxxxxx.SAS                          *
* Folder: x:\xxxx\xxxxxx\xxxxx                 *
* Author: xxxxxxx                              *
* Created: xx/xx/xx                            *
* Summary: xxxxxx xxx xxxxx xxxxxxxxxxx xxx    *
* Revisions: xxxxx xxx xxxx xxxxxxxxxxxx       *
***********************************************/
/* Please provide info about this repository */
%let repo_name =SAS-project-template; /* Repository name on GitHub */
%let repo_maintainer = !!!Insert-name;
%let repo_description = !!!One line description. Do NOT use special characters;

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


libname fmts 'S:\Intmed_Rsrch2\GenMed\Restricted\BP COG\Aim 1\Data Management\Data\formats';
libname dem 'S:\Intmed_Rsrch2\GenMed\Restricted\BP COG\Aim 2\Dementia Model\SAS Data';
options fmtsearch=(fmts);


/* NOTE: For this program to run, any categorical variables must be stored as either numeric variables with associated formats or character variables*/

%let data_dir = dem;
%let analyticfile = coxreg4;
%let OutputPath = S:\Intmed_Rsrch2\GenMed\Restricted\BP COG\Aim 2\Dementia Model\Output;
%let filedesc = meangcp_educ;
%let varlist = mean_gcp_all age0 female0 educ;
%let intxlist =;
%let reflist = 1~5; 
	/*	Reference value for categorical variables in the order they appear in above list, separated by ~ 
		If you have a character variable, include the reference text string in the list (including any spaces or special characters except ~)
		e.g.: 
			%let varlist = mean_gcp_all age0 female0 charvar educ;
			%let reflist = 1~Reference Level~5; 
	*/

%macro VarLists;

%let nvars = %sysfunc(countw(&varlist));
data varord;			
length NAME $32;
do i = 1 to &nvars;
	NAME = scan("&varlist",i);
	vord = i;
	output;
end;
drop i;
run;

proc contents data=&data_dir..&analyticfile noprint out=contents (keep=NAME TYPE FORMAT); run;

data varfmts;
set contents;
if index(" &varlist ", ' ' || strip(NAME) || ' ');
NAME = strip(NAME);
run;

proc delete data=contents; run;
proc sort data=varfmts; by NAME; run;
proc sort data=varord; by NAME; run;
data varfmts;
merge varfmts varord;
by NAME;
run;

proc delete data=varord; run;
proc sort data=varfmts; by vord; run;

data lists;
set varfmts;
length catlist numlist typlist fmtlist $ 200;
retain catlist numlist typlist fmtlist;
if _N_ = 1 then do;
	catlist=' '; numlist = ' '; typlist = ' '; fmtlist = ' ';
end;
if TYPE=2 or (TYPE=1 and not missing(FORMAT)) then do;
	catlist = strip(catx(' ',catlist,NAME));
	typlist = strip(catx(' ',typlist,TYPE));
	if not missing(FORMAT) then fmtlist = strip(catx(' ',fmtlist,FORMAT)); 
	else fmtlist = strip(catx(' ',fmtlist,'null'));
end;
else numlist = strip(catx(' ',numlist,NAME));
ncats=countw(numlist); nnums=countw(catlist);
if _N_ = &nvars;
keep catlist numlist typlist fmtlist ncats nnums;
run;

proc delete data=varfmts; run;
%mend;


%macro BLCovs;
data _null_;
set lists;
call symputx ("catlist",catlist); call symputx ("ncats",ncats);
call symputx ("numlist",numlist); call symputx ("nnums",nnums);
call symputx ("typlist",typlist);
call symputx ("fmtlist",fmtlist);
run;

proc means data=&data_dir..&analyticfile median;
var &numlist;
where t_start=0;
ods output summary=medians;
run;

data _null_;
set medians;
%do i = 1 %to &nnums;
	%let tnum = %scan(&numlist,&i);
	call symput("&tnum._Med",&tnum._Median);
%end;
run;

proc delete data=medians; run;

data blcovs;
%do i = 1 %to &ncats;
	%let ttyp = %scan(&typlist,&i);
	%let tfmt = %scan(&fmtlist,&i);
	%let tcat = %scan(&catlist,&i);
	length &tcat._txt $ 200;
	&tcat = %scan(&reflist,&i,~);
	%if &ttyp = 1 %then %do; 
		format &tcat &tfmt..;
		&tcat._txt = strip(put(&tcat,&tfmt..));
	%end;
	%else %do; 
		length &tcat $ 200; 
		&tcat._txt = &tcat;
	%end;	
%end;
%do i = 1 %to &nnums;
	%let tnum = %scan(&numlist,&i);
	&tnum = &&&tnum._Med;
%end;
run;

data refchar; 
set blcovs;
drop &varlist;
run;

data blcovs;
set blcovs;
keep &varlist;
run;

%mend;


%macro tst1;

%Varlists;
%BLcovs;

data _null_;
set lists;
call symputx ("catlist",catlist); call symputx ("ncats",ncats);
run;

data _null_;
set refchar;
%do i = 1 %to &ncats;
	%let tcat = %scan(&catlist,&i);
	call symputx ("&tcat._ref",&tcat._txt);
%end;
run;

/* Get today's date in YYMMDD format for file names */
%let y2=%sysfunc(today(),year2.);
%if %sysevalf(%sysfunc(today(),month2.)) <= 9 %then %do; %let m2 = %sysfunc(compress(0%sysfunc(today(),month2.))); %end;
%else %do; %let m2 = %sysfunc(today(),month2.); %end;
%if %sysevalf(%sysfunc(today(),day2.)) <= 9 %then %do; %let d2 = %sysfunc(compress(0%sysfunc(today(),day2.))); %end;
%else %do; %let d2 = %sysfunc(today(),day2.); %end;

ods graphics on;
*ods rtf file="&OutputPath.\Cox_&y2.&m2.&d2._&filedesc";
proc phreg data=&data_dir..&analyticfile plots(overlay)=(survival);
%do i = 1 %to &ncats;
	%let tcat = %scan(&catlist,&i);
	class &tcat (ref="&&&tcat._ref");
%end;
model (t_start,t_end)*deminc(0) = &varlist &intxlist;
*baseline covariates=blcov out=sr.bl_surfunc_200225_meangcp_edu survival=S lower=S_lower upper=S_upper;
*ods output parameterestimates = sr.coxparams_200225_meangcp_edu;
run;
ods rtf close;
ods graphics off; 

%mend;

%tst1;



