proc fcmp outlib=sasuser.ufuncs.utils;

 function BMI(height, weight);
/*@ BMI(height,weight)  Calculates BMI. Height in m, weight in kg */
/* height: height in meters */
/* weight: weight in kg */
 return((weight)/(height*height));
 endsub;
 
 function BSA(height,weight);

/* height: height in meters */
/* weight: weight in kg */
/*@ BSA(height,weight)  Calculates BSA using Haycock formula. Height in m, weight in kg*/
 hght = height*100;  /* height expressed in cm */
 BSA=0.024265*(HGHT **0.3964) *(WEIGHT **0.5378);  /* Calculate BSA using Haycock formula */
 return(BSA); /* BSA in m^2 */
endsub;
run;
 

