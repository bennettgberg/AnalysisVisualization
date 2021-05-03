
#script to run the whole routine for skimming and making plots for a new category.
# run MakeDistributions 3 times to get skimmed root files, then run MakePlots.

import os, sys
import tauFun as tf

#test name
testn = "test56"

for cat in ['mmet', 'mmem']:  #'mtet' 'mtem'

    #get category number
    catn = tf.catToNumber(cat)
    print( "catn is %s"%(cat) )

    #make processes_special file
    os.system("cp  processes_special_mmtt.yaml processes_special_%s.yaml"%(cat))

    #make new category yaml file.
    os.system( "sed \"s/xxxx/%s/\" cat_xxxx_2017.yaml > tmp.txt"%(cat))
    #replace the category number
    os.system( "sed \"s/\\\"cat\\\",\\\"==\\\",16/\\\"cat\\\",\\\"==\\\",%d/\" tmp.txt > tmp2.txt"%(catn))
    os.system("rm tmp.txt")
    #replace the category cuts (for deepTau and isGlobal/isTracker cuts).
    ccline = "categorycuts: ["
    #for muons, add isGlobal/isTracker requirement; for taus, add deeptau requirement
    # also add iso requirement for muons, and electrons! (depends on the channel)
    for i,ptype in enumerate(cat):
        if i != 0:
            ccline += ","
        
        #need extra \\ in front of " because will be quote-ception.
        if ptype == 'm':
            ccline += "\[\[\\\"EQT\\\"\],\[\\\"isGlobal_%d\\\",\\\"isTracker_%d\\\"\],\\\"add\\\",\\\">\\\",0\]"%(i+1, i+1)
            #add iso
            if i == 1 or i == 3:
                #mm
                if cat[i-1] == 'm':
                    iso = 0.25
                #em
                elif cat[i-1] == 'e':
                    iso = 0.4
                else:
                    print("Error: unrecognized channel {}".format(cat[i-1:i+1]))
                    sys.exit()
            else:
                #mt
                if cat[i+1] == 't':
                    iso = 0.3
                #mm (for 2nd half of the leg)
                elif cat[i+1] == 'm':
                    iso = 0.25
                else:
                    print("Error: unrecognized channel {}".format(cat[i:i+2]))
                    sys.exit()
            ccline += ",\[\\\"iso_%d\\\",\\\"<\\\",%f\]"%(i+1, iso)
                    
        elif ptype == 't':
            ccline += "\[\\\"idDeepTau2017v2p1VSjet_%d\\\",\\\">=\\\",\*ffworkingpoint\]"%(i+1)

        elif ptype == 'e':
            #em
            if cat[i+1] == 'm':
                iso = 0.3
            #et
            elif cat[i+1] == 't':
                iso = 0.3

            ccline += "\[\\\"iso_%d\\\",\\\"<\\\",%f\]"%(i+1, iso)

        if i == 3:
            ccline += "\]"

    os.system("sed \"s/categorycuts: \[\]/{}/\" tmp2.txt > cat_{}_2017.yaml".format(ccline, cat))
    os.system("rm tmp2.txt")
               
    print("new category files written!")

    #run MakeDistributions twice (separated to avoid MemoryError).
    os.system("cp mcsamples_0_2017.csv bpgMCsamples_2017_v7.csv")
    os.system("python MakeDistributions.py -o %s -ch %s"%(testn, cat)) 
    os.system("cp mcsamples_1_2017.csv bpgMCsamples_2017_v7.csv")
    os.system("python MakeDistributions.py -o %s -ch %s"%(testn, cat))

    #run MakeDistributions one more time to make the root file.
    os.system("python MakeDistributions.py -o %s -ch %s --comb 1"%(testn, cat))

    print( "Making final plots.")
    #run MakePlots to make all the plots.
    os.system("python MakePlots_bpg.py -o %s -ch %s"%(testn, cat))

