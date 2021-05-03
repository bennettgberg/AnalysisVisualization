echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
#common
input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/
process=processes_special_mmmt.yaml
csv=MCsamples_2016_v7.csv
#channel specific
output0=2016_FF_dm_mmmt
output1=2016_FF_mmmt
fo=2016_FF_mmmt
cat=cat_mmmt_2016.yaml


#python MakeDistributions_HAA_2016.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmmt
python MakeDistributions.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmmt
#python MakeDistributions.py -c $cat  -csv $csv  -i $input -p $process -o $output0 -fo $fo -ch mmmt

#echo "Copying over Fake Factors ..."
#cp -r out$output0/pt_*.root FFhistos_$fo/.

# echo "Applying Fake Factors ..."
# python MakeDistributions.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmmt -s -ddZH -fi $fo

# echo "Making Plots ..."
# python MakePlots_skimmed.py -i skimmed_${output1}_mmmt_inclusive.root -o $output1 -c $cat --ch mmmt -ddZH -p $process
# python MakePlots_skimmed_sys.py -i skimmed_mmmt.root -o $output1 -c $cat --ch mmmt -p $process

# echo "Copying Plots Over ..."
# cp -r FFhistos_$fo/*.png outplots_$output1/.
# cp -r outplots_$output1/ /eos/home-s/shigginb/HAA_Plots/.
