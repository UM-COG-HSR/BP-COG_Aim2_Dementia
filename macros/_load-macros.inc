%let tmp=&SAS_work_dir\macros;

/* Assigns value to `computer_name` global macro variable */
filename ff "&tmp\computer_name.sas"; 
%include ff;

/* Returns SAS work directory */
/** DO NOT INCLUDE `SAS_work_directory.sas` */
/** Instead copy it to `C:\Users\Public` drectory */
filename ff "&tmp\SAS_work_directory.sas"; 
%***include ff;     /* Typically not needed */

/*- Assigns global macro variables -*/
filename ff "&tmp\our_sas_session_info.sas"; 
%include ff;

/* proc odstext used to print basic info about executed SAS session */
filename ff "&tmp\sas_session_info_print.sas"; 
%include ff;


