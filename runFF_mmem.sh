echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
#common
input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/
process=processes_special_mmem.yaml
csv=MCsamples_2016_v6_yaml.csv
#channel specific
output0=2016_full_dm_mmem
output1=2016_full_mmem
fo=2016_full_mmem
cat=cat_mmem_2016.yaml

python MakeDistributions_HAA_2016.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmem

echo "Copying over Fake Factors ..."
cp -r out$output0/pt_*.root FFhistos_$fo/.

echo "Applying Fake Factors ..."
python MakeDistributions_HAA_2016.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmem -s -ddZH -fi $fo

echo "Making Plots ..."
python MakePlots_skimmed.py -i skimmed_$output1.root -o $output1 -c $cat --ch mmem -ddZH

echo "Copying Plots Over ..."
cp -r FFhistos_$fo/*.png outplots_$output1/.
cp -r outplots_$output1/ /eos/home-s/shigginb/HAA_Plots/.
