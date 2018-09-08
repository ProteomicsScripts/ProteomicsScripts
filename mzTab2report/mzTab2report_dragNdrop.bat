@echo OFF
ECHO.
ECHO This batch file is used to invoke the mzTab2report script by drag'n'dropping of
ECHO a mzTab file onto this .bat file.
ECHO Within the same txt folder, a QC report in Html/PDF format will be created.
ECHO.
ECHO First time installation of this script [for admins only]:
ECHO The script expects a certain structure to be present.
ECHO This .bat file expects a subfolder named '_internal'
ECHO which holds:
ECHO      'R-3.1.0'             Holds a complete R installation 
ECHO                            (including all packages required to run the 
ECHO                             PTXQC package)
ECHO      'compute_QC_report.R' The R script that is invoked by this .bat


REM A usage message with multiline (i.e. keep the empty line!)
setlocal EnableDelayedExpansion
set USAGE= ^

Usage: ^

%0 ^<mzTab-file^>  ^

Please try again.
REM end of USAGE

REM check number of arguments (must be 1 [mzTab file])
set argC=0
for %%x in (%*) do Set /A argC+=1
if %argC%==0 (
  ECHO Not enough arguments^^! (No arguments provided^)
  ECHO !USAGE!
  goto end
)

REM Check upon first argument (the mzTab file)
REM use bare %1 not "%1", it works fine, even with spaces
if not exist %1 (
  goto not_found
) 

REM ~dp (drive,path); the final '\' is required!; manual quoting required as well
set MZTABFILE="%~dp1\"

ECHO mzTab file is at '%MZTABFILE%'

rem The script path consists of drive (d) and path (p) of the zeroth argument (0) i.e. the script itselfs.
SET SCRIPT_PATH=%~dp0

REM Execute 'mzTab2report.bat' script with the given mzTab file.
CALL %SCRIPT_PATH%mzTab2report.bat %MZTABFILE%


:not_found
echo Could not find the given (network) folder '%1'^^! Please contact your admin/bioinformatician^^!
goto end

:end
pause
