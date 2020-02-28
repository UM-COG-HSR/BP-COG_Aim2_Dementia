%macro computer_name / des = "Assigns value to `computer_name` global macro variable";
data _null_;
  format obparms $char1000.;
  obparms=getoption('objectserverparms');
  x1 = index(obparms,"port=");
  port = substr(obparms,x1+5);
  x1 = index(port," ");
  port = substr(port,1,x1);
  *** call symput('_PORT',trim(left(put(port,6.))));  
  host =sysget('COMPUTERNAME');    /* Use this statement on WINDOWS */
  /*  host =sysget('HOST'); */     /* Use this statement on UNIX    */
  call symput('computer_name',trim(left(put(host,$80.)))); 
run;
%mend computer_name;
