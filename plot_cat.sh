#!/bin/bash


testn="test45"
cat='mtmt'

#get category number
catn=`python getCatn.py $cat`
echo "catn is $catn"

#make processes_special file
cp  processes_special_mmtt.yaml processes_special_${cat}.yaml

#make new category yaml file.
sed "s/mtmt/${cat}/" cat_mtmt_2017.yaml > tmp.txt
sed "s/\"cat\",\"==\",16/\"cat\",\"==\",$catn/" tmp.txt > cat_${cat}_2017.yaml
rm tmp.txt

#run MakeDistributions twice (separated to avoid MemoryError).
cp mcsamples_0_2017.csv bpgMCsamples_2017_v7.csv
python MakeDistributions.py -o $testn -ch $cat 
cp mcsamples_1_2017.csv bpgMCsamples_2017_v7.csv
python MakeDistributions.py -o $testn -ch $cat

#run MakeDistributions one more time to make the root file.
python MakeDistributions.py -o $testn -ch $cat --comb 1

echo "Making final plots."
#run MakePlots to make all the plots.
python MakePlots_bpg.py -o $testn -ch $cat

echo "Last command:"
echo "python MakePlots_bpg.py -o $testn -ch $cat"
