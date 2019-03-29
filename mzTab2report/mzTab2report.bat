@ECHO OFF

REM Script to create a report from mzTab file

REM script directory
REM The script path consists of drive (d) and path (p) of the zeroth argument (0) i.e. the script itselfs.
SET SCRIPT_PATH=%~dp0

IF "%1"=="" (
  ECHO Please specify a file.
  REM Exit the script without closing the terminal
  EXIT /B
)

IF NOT EXIST %1 (
  ECHO File does not exist.
  EXIT /B
)

SET CURRENT_PATH=%CD%
SET FILE=%1
REM The absolute path consists of drive (d), path (p), name (n) and extension (x) of the first argument (1).
SET FILE_ABSOLUTE=%~dpnx1
REM The file path consists of drive (d) and path (p) of the first argument (1).
SET FILE_PATH=%~dp1
REM The base name consists only of the name (n) of the first argument (1).
SET FILE_BASE=%~n1

ECHO Generating report from mzTab file %FILE_BASE%.

REM Unique directory to avoid name clashes in `analysis.mzTab` etc. when 
REM running multiple processes at once.
SET WORK_DIRECTORY=%SCRIPT_PATH%\%FILE_BASE%
MKDIR %WORK_DIRECTORY%
CD /d %WORK_DIRECTORY%

REM copy mzTab
COPY %FILE_ABSOLUTE% analysis.mzTab
COPY %SCRIPT_PATH%\Sweave.sty Sweave.sty

REM  replace dummy FILE_NAME_DUMMY by file name %FILE_BASE%
(for /f "delims=" %%i in (%SCRIPT_PATH%\mzTab2report.Snw) do (
    set "line=%%i"
    setlocal enabledelayedexpansion
    set "line=!line:FILE_NAME_DUMMY=%FILE_BASE%!"
    echo(!line!
    endlocal
))>"mzTab2report_temp.Snw"

REM  Run the R code
R -e "Sweave('mzTab2report_temp.Snw')"

REM If sweave fails, exit with error code.
if %errorlevel% neq 0 (
    ECHO "Sweave encountered an error! Exiting."
    exit \b %errorlevel%
)

REM Small 5 sec pause.
timeout /t 5 /NOBREAK

REM  Run LaTeX code
pdflatex mzTab2report_temp.tex

REM  Copy final report to the input folder
MOVE mzTab2report_temp.pdf %FILE_PATH%\%FILE_BASE%.pdf

REM  clean-up
CD /d %SCRIPT_PATH%
RMDIR /Q /S %WORK_DIRECTORY%

REM  Jump back to original folder
CD /d %CURRENT_PATH%
