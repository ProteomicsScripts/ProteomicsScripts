@ECHO OFF

rem OpenMS TOPP tool directory
SET OpenMSHome=C:\Code\git\OpenMS\bin\Release

rem Script directory
SET SCRIPT_PATH=%~dp0

rem Check if a file is passed as a parameter
IF "%1"=="" (
ECHO Please specify a file.
rem Exit the script without closing the terminal
EXIT /B
)
rem Save the directory where consensusXML file is saved
SET input_directory=%CD%

rem Input file
SET file=%1

rem Set absolute file path
SET file_absolute=%input_directory%\%1
rem Get the file name without extension
SET file_base=%~n1

IF NOT EXIST %1 (
ECHO File does not exist.
EXIT /B
)

ECHO "Generating report from OpenMS consensusXML file"

rem Changing the directory to that of the script
CD /d %SCRIPT_PATH%

COPY %file_absolute% analysis.consensusXML
START %OpenMSHome%\%TextExporter -separator , -in analysis.consensusXML -out analysis.csv

rem Small 5 sec pause.
timeout /t 5 /NOBREAK

rem Run the R code
R -e "Sweave('cXML2report.Snw')"

rem Small 5 sec pause.
timeout /t 5 /NOBREAK

rem Run LaTeX code
pdflatex cXML2report.tex

rem Copy final report to the input folder
MOVE cXML2report.pdf %input_directory%\%file_base%.pdf

rem Clean-up
DEL cXML2report.tex
DEL cXML2report.aux
DEL cXML2report.ini
DEL cXML2report.log
DEL cXML2report.out
DEL density*
DEL Ratio*
DEL analysis*

rem Jump back to input folder
CD /d %input_directory%











