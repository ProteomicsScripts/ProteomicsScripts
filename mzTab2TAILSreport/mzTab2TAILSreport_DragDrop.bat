@ECHO OFF

REM Script for executing mzTab2TAILSreport by simply drag-n-drop an mzTab file to this batch script.

REM Let us assume mzTab2TAILSreport is in the $PATH.
mzTab2TAILSreport "%1"

REM Let us assume mzTab2TAILSreport is not in the $PATH.
REM D:\data\Scripts\ProteomicsScripts\mzTab2TAILSreport\mzTab2TAILSreport.bat "%1"
