clear;

# OpenMS executables
OpenMSHome='/home/lars/Code/OpenMS-develop/OpenMS-build/bin'

# parameter file
ParameterFile='params'
ParameterFile_temp='params_temp'

# input file
#file='/home/lars/Code/git/ProteomicsScripts/OpenQuant/BM4351_2000_2200'
file='BM4351_2000_2200'

# loop parameters
LoopStart=5
LoopEnd=15
LoopStep=0.5

# loop over single parameter
param=$LoopStart
while (( $(echo "$param <= $LoopEnd" |bc -l) ))
do
    echo "parameter = $param (start = $LoopStart, end = $LoopEnd, step = $LoopStep)"
    
    # replace dummy by file name
    sed -e 's/TOBEREPLACED/'$param'/g' $ParameterFile.ini > $ParameterFile_temp.ini
        
    # run peptide detection
    $OpenMSHome/FeatureFinderMultiplex -ini $ParameterFile_temp.ini -in $file.mzML -out $file.featureXML -out_multiplets $file.consensusXML

    # count quants
    $OpenMSHome/MzTabExporter -in $file.consensusXML -out $file\_$param.mzTab
    
    # param = param + LoopStep
    param=$(echo "scale=10; $param + $LoopStep" | bc)
done
