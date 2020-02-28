
%macro our_sas_session_info / des = "Defines global macro variables ";
 /*---- Feb 2020  atg -----*/
 /* For Windows only. Can be adopted to UNIX */
 /* Get the Host Name and Port Number of the             */
 /* Stored Process Server that executed this request.    */

%global sas_batchmode sas_progname sas_fullname;

/* sas_batchmode */
%if %length(%sysfunc(getoption(sysin))) > 0 %then  %let sas_batchmode=Y;
    %else %let sas_batchmode=N;; 

/* sas_full program name */
%if &sas_batchmode=Y %then %do;
 %let sas_fullname = %sysfunc(getoption(sysin));;
 /* Count separators \ */
 %let sep_cnt=%sysfunc(count(&sas_fullname,%str(\)));
 %let sas_progname = %qscan(&sas_fullname, %eval(&sep_cnt+1),%str(\));
%end;

%if &sas_batchmode=N %then %do;
 %let sas_fullname = %sysget(SAS_EXECFILEPATH);
 %let sas_progname = %sysget(SAS_EXECFILENAME);
%end;

%put ==== Macro our_sas_session_info executed ======;
%put sas_batchmode           := &sas_batchmode;   /*--- Y/N */
%put sas_fullname            := &sas_fullname;    /*--- Full path to SAS program */
%put sas_progname            := &sas_progname;    /*--- For example: __sas_session_info ---*/
%put ==== Macro our_sas_session_info ended ======;
%mend our_sas_session_info;
