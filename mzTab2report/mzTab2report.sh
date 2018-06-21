# Script to create a report from mzTab file

#!/bin/sh

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

WORK_DIRECTORY="$SCRIPT_PATH/$FILE_NAME"

if ! [[ -f $FILE ]]
then
echo "File does not exist."
exit
fi


echo 'Generating report from mzTab file '$FILE_ABSOLUTE'.'
mkdir "$WORK_DIRECTORY"
cd "$WORK_DIRECTORY"

# copy mzTab
cp $FILE_ABSOLUTE analysis.mzTab

# replace dummy by file name
sed -e 's/FILE_NAME_DUMMY/'$FILE_NAME'/g' "$SCRIPT_PATH/mzTab2report.Snw" > mzTab2report_temp.Snw

# Run the R code
R -e "Sweave('mzTab2report_temp.Snw')"

# Run LaTeX code
pdflatex mzTab2report_temp.tex

# Copy final report to the input folder
mv mzTab2report_temp.pdf $FILE_PATH/$FILE_NAME.pdf

# clean-up
rm analysis*
rm plot*
rm mzTab2report_temp*

# Jump back to input folder
cd $CURRENT_PATH
