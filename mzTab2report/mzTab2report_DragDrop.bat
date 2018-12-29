@ECHO OFF

REM Script for executing mzTab2report by simply drag-n-drop an mzTab file to this batch script.

REM Let us assume mzTab2report is in the $PATH.
mzTab2report "%1"

REM Let us assume mzTab2report is not in the $PATH.
REM D:\data\Scripts\ProteomicsScripts\mzTab2report\mzTab2report.bat "%1"
