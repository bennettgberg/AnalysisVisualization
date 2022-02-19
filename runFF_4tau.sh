echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
year=2017 #2017
catname='mmmt' #'mmtt' #'mmem' #'mmet' #'mmtt'
sysstr=''
systematics=1
if [ systematics ]
then
    sysstr=' -sys '
fi
#common
#input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/
#input=/eos/uscms/store/user/bgreenbe/haa_4tau_${year}/all/
input=/eos/uscms/store/user/bgreenbe/haa_4tau_${year}/all/
#input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/
process=processes_special_${catname}_${year}.yaml
#csv=MCsamples_2016_v6.csv
#CORRECT csv below!!!!!
csv=samples_${year}_v7.csv
#csv=samples_${year}_v7_test.csv
#echo "WARNING: using zh paper csv!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#csv=MCsamples_2016_v6_yaml.csv
#csv=MCsamples_2016_v6_JetsInc.csv
#csv=MCsamples_2016_v6_NJets.csv
#channel specific
if [ -z $1 ]
then
    mainout='systest0'    #'testan'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi 
#output0=2016_${mainout}_dm_mmtt
#output1=2016_${mainout}_mmtt
output0=${year}_${mainout}_dm_${catname}
output1=${year}_${mainout}_${catname}
#fo=2016_${mainout}_dm_mmtt
#cat=cat_mmtt_2016.yaml
fo=${year}_${mainout}_dm_${catname}
####noff: remove 'dm' from fo
if [ $catname == 'mmem' ]
then    
    fo=${year}_${mainout}_${catname}
fi
cat=cat_${catname}_${year}.yaml


#python MakeDistributions_v6.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmtt

#trying NO ff for cat mmem (noff: use below instead)
if [ $catname == 'mmem' ]
then
    python MakeDistributions_v7.py -c $cat  -csv $csv  -i $input -p $process -o $output1 -fo $fo -ch $catname -yr $year -pc 4 --extract $sysstr -dbg
#else
##UNCOMMENT BELOW LINE to measure fake rate!!!!!!
#    python MakeDistributions_v7.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch $catname -yr $year -pc 4 --extract $sysstr -dbg
fi

#python MakeDistributions_v6.py -c $cat  -csv $csv  -i $input -p $process -o $output0 -fo $fo -ch mmtt
#python MakePlots_skimmed_sys.py -i skimmed_${output0}.root -o $output0 -c $cat --ch mmtt -p $process

#echo "Copying over Fake Factors ..."
#cp -r out$output0/pt_*.root FFhistos_$fo/.

#echo "Applying Fake Factors ..."
#python MakeDistributions_HAA_2016.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmtt -s -ddZH -fi $fo
#python MakeDistributions_v6.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmtt -s -ddZH -fi $fo -fo $output1
#python MakeDistributions_v7.py -csv $csv -i $input -p $process -c $cat -o $output1 -ch $catname -s -ddZH -fi $fo -fo $output1 -yr $year -pc 4 -dbg
#UNCOMMENT BELOW LINE to apply fake rate!!! And remove -co argument to actually use the measurements from above.
#if [ $catname != 'mmem' ]
#then
#    python MakeDistributions_v7.py -csv $csv -i $input -p $process -c $cat -o $output1 -ch $catname -s -ddSM -fi $fo -fo $output1 -yr $year -pc 4 --extract $sysstr -dbg
#fi
#noff: comment out above line!
#echo "python MakeDistributions_v6.py -csv $csv -i $input -p $process -c $cat -o $output1 -ch $catname -s -ddZH -fi $fo -fo $output1 -co"

#echo "Making Plots ..."
#python MakePlots_skimmed.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmtt -ddZH -p $process
#python MakePlots_skimmed.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmtt -p $process
#
#python MakePlots_skimmed_sys.py -i skimmed_${output1}.root -o $output1 -c $cat --ch mmtt -p $process
#python MakePlots_skimmed_sys.py -i skimmed_${output1}.root -o $output1 -c $cat --ch $catname -p $process

#Uncomment to plot!!
#noff: comment out validation
for fakecat in "${catname}_inclusive" #"${catname}_FF_SS_validation" #"${catname}_FF_SS_1_loose" "${catname}_FF_SS_1_tight" "${catname}_FF_SS_2_loose" "${catname}_FF_SS_2_tight"   
do
    #year=1718
    #output1=${year}_${mainout}_${catname}
    python MakePlots_bpg.py -i skimmed_${output1}.root -o $output1 --ch $catname -p $process -fc $fakecat -yr $year -tn $mainout -nd $sysstr
done


#echo "Copying Plots Over ..."
#cp -r FFhistos_$fo/*.png outplots_$output1/.
#cp -r outplots_${output0}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
#cp -r outplots_${output1}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
#cp -r outplots_${output1}_cats/ /eos/home-s/shigginb/HAA_Plots/.

#would be great to get it to make a sound at the end so I know it's done...
#python beep.py
