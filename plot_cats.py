
#script to run the whole routine for skimming and making plots for a new category.
# run MakeDistributions 3 times to get skimmed root files, then run MakePlots.

import os, sys
import tauFun2 as tf

#test name
testn = "test90"



def makeCatFiles(cat):
    #get category number
    catn = tf.catToNumber(cat)
    print( "catn is %s"%(cat) )

    #make processes_special file
    if not os.path.exists("processes_special_%s.yaml"%cat):
        os.system("cp  processes_special_mmtt.yaml processes_special_%s.yaml"%(cat))
    #if it's not already there

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

            #add medium Id requirement IF not tt** category. ##jk, always add it.
            #if not (cat[0] == 't' and cat[1] == 't'):
            ccline += ",\[\\\"mediumId_%d\\\",\\\"==\\\",1\]"%(i+1)
                    
        elif ptype == 't':
            for vnum,vsboi in enumerate(["jet", "e", "mu"]):
                ccline += "\[\\\"idDeepTau2017v2p1VS%s_%d\\\",\\\">=\\\",\*ffworkingpoint\]"%(vsboi, i+1)
                if vnum != 2:
                    ccline += ","

        elif ptype == 'e':
            #em
            if cat[i+1] == 'm':
                iso = 0.15 # 0.3
            #et
            elif cat[i+1] == 't':
                iso = 0.15 #0.3

            ccline += "\[\\\"iso_%d\\\",\\\"<\\\",%f\]"%(i+1, iso)

            #loosen the electron requirement for tt** categories. ##jk don't do this.
            #if not (cat[0] == 't' and cat[1] == 't'): 
            ccline += ",\[\\\"Electron_mvaFall17V2noIso_WP90_%d\\\",\\\">\\\",%f\]"%(i+1, 0.5)

        if i == 3:
            ccline += "\]"

    os.system("sed \"s/categorycuts: \[\]/{}/\" tmp2.txt > cat_{}_2017.yaml".format(ccline, cat))
    os.system("rm tmp2.txt")
               
    print("new category files written!")

#run the full test
def runTest(cat, allBkgs=True):
    #if allBkgs is true then need to use all backgrounds; otherwise just use what's already in bpgMCsamples.
    #run MakeDistributions twice (separated to avoid MemoryError).
    #  or only once if just using the contents of bpgMCsamples.
    if allBkgs:
        os.system("cp mcsamples_0_2017.csv bpgMCsamples_2017_v7.csv")
        os.system("python MakeDistributions.py -o %s -ch %s"%(testn, cat)) 
        os.system("cp mcsamples_1_2017.csv bpgMCsamples_2017_v7.csv")
    os.system("python MakeDistributions.py -o %s -ch %s"%(testn, cat))

    #run MakeDistributions one more time to make the root file.
    os.system("python MakeDistributions.py -o %s -ch %s --comb 1"%(testn, cat))

#noData should be true if don't want to include data in the plots.
def makePlots(cat, noData=False, sigOnly=False):
    print( "Making final plots.")
    #run MakePlots to make all the plots.
    os.system("python MakePlots_bpg.py -o %s -ch %s %s %s"%(testn, cat, "-nd" if noData else "", "-so" if sigOnly else ""))


#list of all categories: 'ttmt', 'ttet', 'ttem', 'mtmt', 'mtet', 'mtem', 'etet', 'etem', 'emem', 'mmtt', 'mmmt', 'mmet', 'mmem'

for ct in ['mmtt','mmmt','mmet','mmem']: #'ttmt','ttet','ttem','mtmt','mtet','mtem','etet','etem','mmtt','mmmt','mmet','mmem']:
    makeCatFiles(ct)
    allBkgs = True
    noData = True
    sigOnly = True
    #runTest(ct, allBkgs)
    makePlots(ct, noData, sigOnly)
