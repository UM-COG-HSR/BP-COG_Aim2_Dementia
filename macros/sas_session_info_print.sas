%macro sas_session_info_print/ des = "Prints basic info about SAS program execution";
proc odstext; 
p "Report created by &sas_progname" / style=[color=red fontsize=15pt];
p "Executed on: &sysdate"; 
p "Repository name: &repo_name"; /* Repository name */
p "Repository description: &repo_description"; 
p "Repository maintainer: &repo_maintainer"; 
p "SAS work directory: &SAS_work_dir";

p '';
p "Userid:                 &sysuserid";
p "Computer name:          &computer_name";
p "SAS in batch mode?      &sas_batchmode";
p "SAS Version:            &SYSVLONG4";
p "OS:                     &SYSSCPL";
run;

%mend sas_session_info_print;
