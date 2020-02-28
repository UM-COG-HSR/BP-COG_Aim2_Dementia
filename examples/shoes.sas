 /* link    */
 data shoes; 
   set sashelp.shoes;
 run;
  
 proc contents data=shoes position;
 run;
 
 /*  Integrity Constraints */
 proc datasets library= work nolist;
 modify shoes;
  ic create pkey = primary key (storenumber);
  ic create regprodsub = distinct (region product subsidiary)
   message = "Region, Product, Subsidiary combination must be unique";
  ic create storelimit = check(where=(stores < 50))
   message = "Limit of 50 stores";
  ic create returnsales = check(where=(returns+sales <= inventory))
   message = "Returns + Sales cannot exceed Inventory";
 run;
 quit;
