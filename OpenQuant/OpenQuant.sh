clear;

# OpenMS executables
OpenMSHome='/home/lars/Code/git/OpenMS-build/bin'

# template parameter file
params='params'
params_temp='params_temp'

# input file
file='/media/lars/HDD/data/20181031_FixedRatios_BM3363-70/BM3366'

for i in {0..100}
do
  echo "loop = $i"
  shift=$(echo "scale=4; $i / 17" | bc)
  echo "shift = $shift"

  # replace dummy by file name
  sed -e 's/TOBEREPLACED/'$shift'/g' $params.ini > $params_temp.ini

  # run peptide detection
  #$OpenMSHome/FeatureFinderMultiplex -ini $params_temp.ini -in $file.mzML -out $file.featureXML -out_multiplets $file.consensusXML

  # count quants
  #$OpenMSHome/FileInfo -in out_multiplets $file.consensusXML -out $file_$i.txt

done


