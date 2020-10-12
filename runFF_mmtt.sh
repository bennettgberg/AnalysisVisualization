echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmtt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmtt.yaml -dmZH -o 2016_dm_mmtt -fo 2016_mmtt -ch mmtt

echo "Copying over Fake Factors ..."
cp -r out2016_dm_mmtt/pt_*.root FFhistos_2016_mmtt/.

echo "Applying Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmtt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmtt.yaml -ddZH -o 2016_mmtt -fi 2016_mmtt -ch mmtt

echo "Making Plots ..."
python MakePlots_histos.py -i 2016_mmtt -o 2016_mmtt -ddZH -mhs -c cat_mmtt_2016.yaml --ch mmtt


echo "Copying Plots Over ..."
cp FFhistos_2016_mmtt/*.png outplots_2016_mmtt_mmtt_inclusive/.

echo "Copying fake rate plots Over ..."
cp -r outplots_2016_mmtt_mmtt_inclusive /eos/home-s/shigginb/HAA_Plots/.
