@ECHO OFF

REM Script for executing mzTab2PRTreport by simply drag-n-drop an mzTab file to this batch script.

REM Let us assume mzTab2PRTreport is in the $PATH.
mzTab2PRTreport "%1"

REM Let us assume mzTab2PRTreport is not in the $PATH.
REM D:\data\Scripts\ProteomicsScripts\mzTab2PRTreport\mzTab2PRTreport.bat "%1"
