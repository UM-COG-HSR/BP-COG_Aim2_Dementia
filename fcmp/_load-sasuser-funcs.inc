%let tmp = &SAS_work_dir/fcmp;

/*--- Load files (one at a time) with user defined functions ---*/
filename ff "&tmp/ufuncs-utils.sas";
%include ff; 

filename ff "&tmp/ufuncs-conversions.sas";
%include ff; 

options cmplib=(sasuser.ufuncs);
