# Script to create duplex report from .consensusXML results

#!/bin/sh

# OpenMS TOPP tool directory
OpenMSHome='/home/lars/Code/git/OpenMS-build/bin'

# script directory
SCRIPT_PATH=$(dirname -- "$(readlink -f -- "$0")")

if [[ $1 == "" ]]
then
echo "Please specify a file."
exit
fi

# input file
FILE=$1; shift
FILE_ABSOLUTE=$(readlink -f -- $FILE)
FILE_PATH=$(dirname $FILE_ABSOLUTE)
FILE_BASE=$(basename $FILE_ABSOLUTE)
FILE_NAME=${FILE_BASE%.*}

if ! [[ -f $FILE ]]
then
echo "File does not exist."
exit
fi

echo 'Generating report from OpenMS consensusXML file '$FILE_ABSOLUTE'.'
cd $SCRIPT_PATH

# export to csv
cp $FILE_ABSOLUTE analysis.consensusXML
$OpenMSHome/TextExporter -separator , -in analysis.consensusXML -out analysis.csv

# Run the R code
R -e "Sweave('cXML2report.Snw')"

# Run LaTeX code
pdflatex cXML2report.tex

# Copy final report to the input folder
mv cXML2report.pdf $FILE_PATH/$FILE_NAME.pdf

# clean-up
rm *.consensusXML
rm *.csv
rm *.pdf
rm *.aux
rm *.log
rm *.out
rm *.tex

# Jump back to input folder
cd $CURRENT_PATH
