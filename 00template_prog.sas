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

/*****----  Our program starts here  *****/

/* Basic info on executed SAs program is printed */
ods html file = "_last_sas_session_info.html"
         path = "&SAS_work_dir"
         ;
  %sas_session_info_print; /* Print basic info on executed SAS session */ 
ods html close;
