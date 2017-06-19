@ECHO OFF

rem Script to create a report from mzTab file

rem script directory
SET SCRIPT_PATH=%~dp0

IF "%1"=="" (
  ECHO Please specify a file.
  rem Exit the script without closing the terminal
  EXIT /B
)

SET file=%1
SET file_absolute=%input_directory%\%1
SET file_base=%~n1
SET input_directory=%CD%

rem FILE=$1; shift
rem FILE_ABSOLUTE=$(readlink -f -- $FILE)
rem FILE_PATH=$(dirname $FILE_ABSOLUTE)
rem FILE_BASE=$(basename $FILE_ABSOLUTE)
rem FILE_NAME=${FILE_BASE%.*}

IF NOT EXIST %1 (
  ECHO File does not exist.
  EXIT /B
)

ECHO Generating report from mzTab file %FILE_BASE%.

CD /d %SCRIPT_PATH%

rem copy mzTab
rem cp $FILE_ABSOLUTE analysis.mzTab

rem  replace dummy by file name
rem sed -e 's/FILE_NAME_DUMMY/'$FILE_NAME'/g' mzTab2PSMreport.Snw > mzTab2PSMreport_temp.Snw

rem  Run the R code
rem R -e "Sweave('mzTab2PSMreport_temp.Snw')"

rem  Run LaTeX code
rem pdflatex mzTab2PSMreport_temp.tex

rem  Copy final report to the input folder
rem mv mzTab2PSMreport_temp.pdf $FILE_PATH/$FILE_NAME.pdf

rem  clean-up
rem rm analysis*
rem rm plot*
rem rm mzTab2PSMreport_temp*

rem  Jump back to input folder
CD /d %input_directory%
