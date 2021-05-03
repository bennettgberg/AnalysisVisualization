echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
#common
input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/
#input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/
process=processes_trigger_mmmt.yaml
csv=MCsamples_2016_v6.csv
#csv=MCsamples_2016_v6_yaml.csv
#csv=MCsamples_2016_v6_JetsInc.csv
#csv=MCsamples_2016_v6_NJets.csv
#channel specific
if [ -z $1 ]
then
    mainout='test'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi 
output0=2016_${mainout}_dm_mmmt
output1=2016_${mainout}_mmmt
fo=2016_${mainout}_dm_mmmt
cat=cat_mmmt_trigger_2016.yaml


#python MakeDistributions_v6.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmmt
python MakeDistributions_v6.py -c $cat  -csv $csv  -i $input -p $process -o $output0 -fo $fo -ch mmmt
#python MakePlots_skimmed_sys.py -i skimmed_${output0}.root -o $output0 -c $cat --ch mmmt -p $process

#echo "Copying over Fake Factors ..."
#cp -r out$output0/pt_*.root FFhistos_$fo/.

#echo "Applying Fake Factors ..."
#python MakeDistributions_HAA_2016.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmmt -s -ddZH -fi $fo
#python MakeDistributions_v6.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmmt -s -ddZH -fi $fo -fo $output1

echo "Making Plots ..."
#python MakePlots_skimmed.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmmt -ddZH -p $process
#python MakePlots_skimmed.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmmt -p $process
#
python MakePlots_skimmed_sys.py -i skimmed_${output0}.root -o $output0 -c $cat --ch mmmt -p $process -de
#python MakePlots_skimmed_sys.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmmt -p $process -de



echo "Copying Plots Over ..."
#cp -r FFhistos_$fo/*.png outplots_${output1}_Nominal/.
cp -r outplots_${output0}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
#cp -r outplots_${output1}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
#cp -r outplots_${output1}_cats/ /eos/home-s/shigginb/HAA_Plots/.
