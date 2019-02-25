clear;

# OpenMS executables
OpenMSHome='/home/lars/Code/git/OpenMS-build/bin'

# parameter file
ParameterFile='params'
ParameterFile_temp='params_temp'

# input file
file='/home/lars/Code/git/ProteomicsScripts/OpenQuant/BM4351.mzML'

# loop parameters
LoopStart=2
LoopEnd=15
LoopStep=0.1

# loop over single parameter
param=$LoopStart
while (( $(echo "$param <= $LoopEnd" |bc -l) ))
do
    echo "parameter = $param (start = $LoopStart, end = $LoopEnd, step = $LoopStep)"
    
    # replace dummy by file name
    sed -e 's/TOBEREPLACED/'$param'/g' $ParameterFile.ini > $ParameterFile_temp.ini
    
    # run peptide detection
    $OpenMSHome/FeatureFinderMultiplex -ini $ParameterFile_temp.ini -in $file.mzML -out $file.featureXML -out_multiplets $file.consensusXML

    # export to mzTab
    $OpenMSHome/MzTabExporter -in $file.consensusXML -out $file_$param.mzTab

    # generate report
    mzTab2report $file_$param.mzTab

    # param = param + LoopStep
    param=$(echo "scale=10; $param + $LoopStep" | bc)
done
