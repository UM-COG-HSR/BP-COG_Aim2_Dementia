proc fcmp outlib=sasuser.functions.conversions;

function lb2kg(lb);
 kg = lb/2.2;
 return (kg);
endsub;
run;

