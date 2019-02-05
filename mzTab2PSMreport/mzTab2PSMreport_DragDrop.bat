@ECHO OFF

REM Script for executing mzTab2PSMreport by simply drag-n-drop an mzTab file to this batch script.

REM Let us assume mzTab2PSMreport is in the $PATH.
mzTab2PSMreport "%1"

REM Let us assume mzTab2PSMreport is not in the $PATH.
REM D:\data\Scripts\ProteomicsScripts\mzTab2PSMreport\mzTab2PSMreport.bat "%1"
