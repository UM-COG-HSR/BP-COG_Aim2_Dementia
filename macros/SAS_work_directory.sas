/*******************************************/ 
/*                                         */
/*                                         */ 
/*  Store this file in                     */
/*     C:/Users/Public/ folder             */
/*                                         */ 
/*******************************************/ 

%macro SAS_work_directory / des = "Returns SAS work directory";

/*----------------------------------------------------------------------
Returns the current SAS directory physical name.
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
Originally developed by Tom Hoffman.
Posted in memory of Tom and Fan.
-----------------------------------------------------------------------
Usage:

%put %SAS_work_diectory is the current directory.;
------------------------------------------------------------------------
Notes:

-----------------------------------------------------------------------
History:

11MAR99 TRHoffman Creation - with help from Tom Abernathy.
06DEC00 TRHoffman Used . notation to reference current directory as
                  suggested by Fan Zhou.
----------------------------------------------------------------------*/
%local fr rc SAS_work_dir;

%let rc = %sysfunc(filename(fr,.));
%let SAS_work_dir = %sysfunc(pathname(&fr));
%let rc = %sysfunc(filename(fr));
&SAS_work_dir
%mend SAS_work_directory;


