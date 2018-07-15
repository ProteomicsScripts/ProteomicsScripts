@ECHO OFF

rem Script to convert MaxQuant result files to (OpenMS-compatible) mzTab format.

rem script directory
rem The script path consists of drive (d) and path (p) of the zeroth argument (0) i.e. the script itselfs.
SET SCRIPT_PATH=%~dp0

IF "%1"=="" (
  ECHO Please specify a folder with MaxQuant results.
  rem Exit the script without closing the terminal
  EXIT /B
)

IF NOT EXIST %1 (
  ECHO Folder does not exist.
  EXIT /B
)

SET FOLDER=%1
rem The absolute path consists of drive (d), path (p), name (n) and extension (x) of the first argument (1).
SET FOLDER_ABSOLUTE=%~dpnx1

ECHO Generating report from MaxQuant input folder %FOLDER_ABSOLUTE%.

rem  Run the R code
Rscript %SCRIPT_PATH%/misc/MQ2mzTab.R "%FOLDER_ABSOLUTE"
