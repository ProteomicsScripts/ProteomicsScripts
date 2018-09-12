@echo OFF
REM This batch file is used to invoke the mzTab2report script by drag'n'dropping of an mzTab file onto this .bat file.

REM A usage message with multiline (i.e. keep the empty line!)
setlocal EnableDelayedExpansion
set USAGE= ^

Usage: ^

%0 ^<mzTab-file^>  ^

Please try again.
REM end of USAGE

REM Check number of arguments (must be 1, the mzTab file)
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
  ECHO The file %1 does not exist.
  goto :end
)

REM We assume mzTab2report to be in the PATH.
ECHO Generating PDF report for %1 by calling mzTab2report.
start mzTab2report %1

REM TO-DO: The above line opens a new terminal window. How to close it automatically?

:end
REM pause
exit