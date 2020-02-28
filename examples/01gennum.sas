/***---  Illustrates `gennum` data option. Feb 2020 ***/
/* Maximum number of generations is set */
/* Generation #1 of a dataset is created */
data clss (genmax=100
           label = "Data clss#001 created on &sysdate"         
          );
  set sashelp.class;
run;

/* Generation #2 of a dataset is created */
data clss (label = "Data clss#002 created on &sysdate"         
          );
  set clss(gennum =1);
  weight_kg = weight*0.454;
run;


/* Compare two generations of a dataset*/
proc compare data = clss(gennum=1) 
             compare = clss(gennum=2);
run;
