#!/usr/bin/env python
import ROOT
import re
import math
from array import array
import numpy as np
#from collections import OrderedDict
#import varCfgPlotter
import argparse
import os
import io
import yaml
#from HttStyles import GetStyleHtt
#from HttStyles import MakeCanvas
import uproot
import root_numpy
import sys as system

def add_lumi(year):
    lowX=0.60
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    #get correct luminosity depending on the year.
    luminum = 35.9
    if year == "2017":
        print("year 2017")
        luminum = 41.8
    elif year == "2018":
        luminum = 59.74
        print("year 2018")
    elif year == "1718":
        luminum = 59.74 + 41.8
    else:
        print("year: %d"%year)
    if year == "1718":
        yrstr = "2017+2018"
    else:
        yrstr = year
    lumi.AddText(yrstr + " " + str(luminum) + " fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.17
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary(channel="mmmt"):
    lowX=0.30
    #lowY=0.690
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary "+channel)
    return lumi

def make_legend():
        output = ROOT.TLegend(0.165, 0.45, 0.350, 0.75, "", "brNDC")
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetFillColor(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

def make_legend_inset():
        output = ROOT.TLegend(0.6, 0.5, 0.98, 0.85, "", "brNDC")
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetFillColor(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

def add_categ(text):
       categ  = ROOT.TPaveText(0.6, 0.2+0.013, 0.89, 0.40+0.155, "NDC")
       categ.SetBorderSize(   0 )
       categ.SetFillStyle(    0 )
       categ.SetTextAlign(   12 )
       categ.SetTextSize ( 0.05 )
       categ.SetTextColor(    1 )
       categ.SetTextFont (   41 )
       categ.AddText(text)
       return categ

def normalize_hist(h):
    k=h.Clone()
    for j in range(1,h.GetSize()-1):
        k.SetBinContent(j,k.GetBinContent(j)/k.GetBinWidth(j))
        k.SetBinError(j,k.GetBinError(j)/k.GetBinWidth(j))
    return k

def assignColor(h):
    from utils.Processes import HAA_processes
    name = h.GetName()
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name:
        h.SetLineColor(ROOT.kRed)
        h.SetFillColor(0)
        h.SetLineWidth(4)
    if "ZZ" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["ZZ"].color[0]))
    if "ZTT" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
    if "TT" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["TTJets_1LTbar"].color[0]))
    if "Z" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
    if "W" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["WJetsToLNu"].color[0]))
    if "EWK" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["EWKWMinus2Jets"].color[0]))
    if "ST" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["ST_s"].color[0]))
    if "TTL" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["TTJets_1LTbar"].color[0]))
    if "ZTL" in name:
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
    return h

def makeHisto(h,hS,hB,hD):
    #k=h.Clone()
    name = h.GetName()
    #print "setting attributes for ",name
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name:
        h = assignColor(h)
        hS.Add(h)
    if "ZL" in name or "ZTT" in name or "TTT" in name or "TTL" in name or "ZTL" in name or "jetFakes" in name:
        h = assignColor(h)
        hB.Add(h)




    return h,hS,hB,hD








if __name__ == "__main__":

    #Change these to input category files ?? via parser?
    from utils.Parametrization import Category
    from utils.Parametrization import Process
    import argparse

    parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-i",  "--input", default="",  help="input root file name")
    parser.add_argument("-o",  "--output", default="test43",  help="postfix string")
    parser.add_argument("-ch",  "--channel", default="mtet",  help="postfix string")
  #  parser.add_argument("-c",  "--categories", default="cat_tttt_2017.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="bpgMCsamples_2017_v7_yaml.csv",  help="csv file")
    parser.add_argument("-p",  "--processes", default="processes_special_mtmt.yaml",  help="processes yaml file")
    #parser.add_argument("--dist", default="rareBkg",  help="single distribution to plot")
    parser.add_argument("--dist", default="Bkg",  help="single distribution to plot")
    parser.add_argument("-mc",  "--mc", default=True,action='store_true',  help="Use only mc skip data")
    parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
    parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
    parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
    parser.add_argument("-y",  "--maxY", default=-1,  help="Maximum y-value to show in plots (-1 to auto-set)")
    parser.add_argument("-nd",  "--noData", default=False, action='store_true',  help="Show data on the plot (False) or don't (True)")
    parser.add_argument("-so",  "--sigOnly", default=False, action='store_true',  help="True to include only signal; else False.")
    parser.add_argument("-fc",  "--fakeCat", default="",  help="Category for the fake factor plotting.")
    parser.add_argument("-yr",  "--year", default="2017",  help="Which year (2016, 2017, or 2018")
    parser.add_argument("-nZ",  "--noZZ", default=False,action='store_true',  help="Include this argument to exclude the ZZ bkg")
    parser.add_argument("-tn",  "--testname", default="",  help="Testname (esp. for FF, so don't need to include filename, or can include multiple files)")
    args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    yr = str(args.year) #"2017"
    if yr == "1718":
        yr = "2017"

    #only plot this signal mass.
    mass = 40

    #show ONLY signal
    sigOnly = args.sigOnly # False #True
    noSig = False #True #False

    #maximum value to set the y-axis to (-1 to auto-set)
    y_max = int(args.maxY) #10000 #250

    if sigOnly and noSig:
        print("Error: both sigOnly and noSig cannot be specified.")
        system.exit()

    noData = args.noData #False


    catfile = "cat_%s_%s.yaml"%(args.channel, yr)
    #importing analysis categories and conventions
    with io.open(catfile,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(args.processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print category
        #print categories[category]['name']
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat

    print "the categories   ",allcats
    newvars=[]
    variabledic={}

    #for newvar in  HAA_Inc_mmmt.newvariables.keys():
    for newvar in  allcats[args.channel+"_inclusive"].newvariables:
        newvars.append([newvar])
        variabledic[newvar]=[newvar,allcats[args.channel+"_inclusive"].newvariables[newvar][1],allcats[args.channel+"_inclusive"].newvariables[newvar][3],allcats[args.channel+"_inclusive"].newvariables[newvar][4]]
    variables = []
    for varHandle in allcats[args.channel+"_inclusive"].vars.keys():
        variables.append([varHandle])
        variabledic[varHandle]=allcats[args.channel+"_inclusive"].vars[varHandle]
    #variables = variables + newvars
    #print variables

    #The skimmed root file containing all the TTrees
    #file = ROOT.TFile(args.input,"read")
    if args.input == "":
        input_file = "skimmed_%s_%s.root"%(args.channel, args.output)
    else:
        input_file = args.input
    #year 1718 signifies user wants us to add together 2017 and 2018.
    if args.fakeCat != "" and args.year == "1718" and args.testname != "":
        infiles = ["skimmed_%d_%s_%s.root"%(yyyy, args.testname, args.channel) for yyyy in [2017, 2018]]        
        #print("infiles: {}".format(infiles))
        fim = uproot.open(infiles[0]) #just open the first one for now just to get the names and junk.
        fin = fim[args.fakeCat]
        
    elif args.fakeCat == "":
        fin = uproot.open(input_file)
        infiles = [input_file] 
    else:
        infiles = [input_file]
        fim = uproot.open(input_file)
        fin = fim[args.fakeCat]
        #dists = fin.keys()
        #print("dists: {}".format(dists))
        #system.exit()

    #print("fin: {}, infiles: {}".format(fin, infiles))

    histolist = {}
    finalhists = {}
    processes = []
    histodict = {}
    systematics = ["Nominal"] #,"scale_m_etalt1p2Up"]
    #get all the distributions
    dists = fin.keys()
    #get rid of the useless ;1 at the end of each dist
    for dnum,dist in enumerate(dists):
        dists[dnum] = dist.split(';')[0]
    #remove all the distributions we aren't dealing with (all the signals except the one mass)
    for sys in systematics:
        #list of masses that we aren't doing rn.
        bad_masses = [sys + "_a%d"%mss for mss in range(15, mass, 5)]
        for mss in range(mass+5, 61, 5): bad_masses.append(sys + "_a%d"%mss)
        if noSig: bad_masses.append(sys + "_a%d"%mass)
        #remove the ones we won't use from dists.
        for bm in bad_masses: 
      #      print("bm to remove: {}".format(bm))
            if bm in dists:
                dists.remove(bm)
        #if we're only tryna plot the signal then remove everything else.
        if sigOnly:
            i = 0
            while i < len(dists):
                if sigOnly and "a%d"%mass not in dists[i]:
                    dists.remove(dists[i])
                    i -= 1
                i += 1
        elif noData and (sys+"_data_obs") in dists:
        #    print("removing data observed!")
            dists.remove(sys+"_data_obs")

   ## temprorarily remove all bkg's except one to debug and junk.
   # nd = 0
   # while nd < len(dists):
   #     
   #     dist = dists[nd]
   #     if not dist == "Nominal_rareBkg": #"Nominal_Bkg":
   #         dists.remove(dist)
   #         print("dist to remove: " + dist)
   #         print("now dists = " + str(dists))
   #         nd -= 1
   #     nd += 1

    print("dists: " + str(dists))
        
    #print("dists: {}".format(dists))
    cats = [args.channel+"_inclusive"]
    vars = allcats[cats[0]].vars.keys()
    passvars = []

    for newvar in allcats[cats[0]].newvariables.keys():
        vars.append(newvar)

    #keep the integrals of each distribution so we can sort by integral later.
    integrals = {}
    #for ivar,var in enumerate(cat.vars.keys()):
    #for dist in fin.GetListOfKeys():
    for sys in systematics:
   
        histodict[sys]={}
        integrals[sys] = {}
        for distLong in dists:
            #if the dist is one we're not doing rn, continue.
        #    if distLong in bad_masses: continue
            #print("Now looping over dist: {}".format(distLong))
            dist = distLong.split(";")[0]
            #distribution = dist.ReadObj()  #The TTree
            histodict[sys][dist]={}
            integrals[sys][dist]={}
            for cat in cats:
                histodict[sys][dist][cat]={}
                integrals[sys][dist][cat]={}
                for variableHandle in vars:
                   # print variableHandle
                    if variableHandle in allcats[cats[0]].newvariables.keys():
                        var = variableHandle
                        bins = allcats[cat].newvariables[variableHandle][1]
                    else:
                        var = allcats[cat].vars[variableHandle][0]
                        bins = allcats[cat].vars[variableHandle][1]
                    title = "%s_%s_%s"%(str(dist), str(cat), str(variableHandle))
                    if type(bins[0])==list:
                       # print("Making new TH1D! title = {}".format(title))
                        histodict[sys][dist][cat][variableHandle] = ROOT.TH1D(title,title,bins[0][0],bins[0][1],bins[0][2])
                        integrals[sys][dist][cat][variableHandle] = 0 
                    else:
                       # print("ERROR: Integral not calculated in this branch!")
                        tmpbin = np.asarray(bins)
                        histodict[sys][dist][cat][variableHandle] = ROOT.TH1D(title,title,len(tmpbin)-1,tmpbin)
                        integrals[sys][dist][cat][variableHandle] = 0 

    #print("histodict before filling: {}".format(histodict))
    for ii,ff in enumerate(infiles):
        fim = uproot.open(ff)
        if args.fakeCat == "":
            fin = fim
        else:
            fin = fim[args.fakeCat]
        for sys in systematics:
            for distLong in dists:
                dist = distLong.split(";")[0]
                print("dist: " + dist)
                try:
                    tree = fin[dist]
                except:
                    continue
                masterArray = tree.arrays()
                for cat in cats:
                    for variableHandle in vars:
                        if variableHandle in allcats[cats[0]].newvariables.keys():
                            var = variableHandle
                            #bins = allcats[cat].newvariables[variableHandle][1]
                        else:
                            var = allcats[cat].vars[variableHandle][0]
                            #bins = allcats[cat].vars[variableHandle][1]
        
                        try:
                            val = masterArray[var]
                            #print("val: {}".format(val))
                            #print("masterArray[finalweight]: {}".format(masterArray["finalweight"]))
                            root_numpy.fill_hist(histodict[sys][dist][cat][variableHandle],val,masterArray["finalweight"])
                            integrals[sys][dist][cat][variableHandle] = histodict[sys][dist][cat][variableHandle].Integral() 
                            #print("just filled hist {}, {}, {}, {}: {}, int: {}".format(sys, dist, cat, variableHandle, histodict[sys][dist][cat][variableHandle], integrals[sys][dist][cat][variableHandle]))
                            #print("3rd bin: {}".format(histodict[sys][dist][cat][variableHandle][3]))
                            if ii == 0:
                                passvars.append(var)
                        except:
                            print "problem with variable so skipping ",var
                            continue

    #print("histodict after filling: {}".format(histodict))
    #for varCatDist in histodict.keys():

    #var is now the variable handle
    for sys in systematics:
        dirname = "outplots_"+args.output+"_"+sys+"%s%s%s"%("" if not sigOnly else "_sigOnly", "" if noData else "_data", "_noZZ" if args.noZZ else "")
        try:
            os.mkdir(dirname)
        except:
            print "dir prob exists"


        #dictionary to hold all the variables that were already done (to make sure none are repeated ->segfaults)
        donevars = {}

        for var in passvars:
            if var in donevars: 
                continue
            donevars[var] = True
            if args.mc:
                if var=="AMass":
                    fileout = open(dirname+"/"+str(allcats[cats[0]].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cats[0]].name+"\n")
            else:
                if var=="AMass_blinded":
                    fileout = open(dirname+"/"+str(allcats[cats[0]].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cats[0]].name+"\n")
                #if var in ["AMass","mll_m15","mll_m20","mll_m25","mll_m30","mll_m35","mll_m40","mll_m45","mll_m50","mll_m55","mll_m60"]:
                #    continue


            #make a separate plot for each category.
            for cat in cats:
                #signal
                hSignals={}

                hDataDict={}
                hMCDict={}
                hMCFake1Dict={}
                hMCFake2Dict={}

                #make new canvas for this plot
                titl = "%s_%s"%(cat, var)
                c=ROOT.TCanvas(titl,titl,0,0,600,600)
                H = 600
                W = 600

                H_ref = 600
                W_ref = 600


                T = 0.08*H_ref
                B = 0.12*H_ref
                L = 0.16*W_ref
                R = 0.04*W_ref


                B_ratio = 0.1*H_ref
                T_ratio = 0.03*H_ref

                B_ratio_label = 0.3*H_ref


                doRatio=False

                if not doRatio:
                    c.SetLeftMargin(L/W)
                    c.SetRightMargin(R/W)
                    c.SetTopMargin(T/H)
                    c.SetBottomMargin(B/H)
                c.cd()
                #print("cd'd to c")
                #print("made c.")
                #histogram stack to hold the stack of hists from each dist.
             #   print("sigOnly???????????????????????? " + str(sigOnly))
                if not sigOnly:
             #       print("Making hBkgTot!!!!!")
                    hBkgTot = ROOT.THStack()
                    hbkgerr = ROOT.TH1F() #hist to capture sum of all bkgs so can get MC error.
             #       print("Made hBkgTot!!!!!")
             #   else:
             #       print("sigOnly true????!!!!!!!!!!???????? " + str(sigOnly))
               # print("Original hBkgTot = {}".format(hBkgTot))
               # print("THStack x-axis: {}".format(hBkgTot.GetXaxis()))
                #set up legend.
                #should move depending on the variable.
                xR=0.65
                l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
                if var=="mll_fine":
                    l=ROOT.TLegend(0.40,0.55,0.40+0.28,0.9);
                #signal isn't stacked with the backgrounds.
                hSig = ROOT.TH1F()
                hData = ROOT.TH1F()
                #sort dists by Integral!
                sorted_dists = []
                try:
                    #the sig will be the last item, to try to fix this root bug. data will be the first.
                    signame = "" #will be filled in later if it's in dists (else won't be inserted).
                    for di in dists:
                        #if it's the sig, continue bc we'll put this at the front after this loop.
                        if di == (sys + "_a%d"%mass):  
                            signame = di
                            continue
                        if di == (sys + "_data_obs"):
                            #data_obs will be inserted at the front right after this loop.
                            continue
                        idx = 0
                        #find greatest integral to go on the bottom.
                        while idx < len(sorted_dists) and integrals[sys][di][cat][var] < integrals[sys][sorted_dists[idx]][cat][var]:
                            idx += 1
                        sorted_dists.insert(idx, di)
                    #now append the sig, if it's here.
                    if signame != "":
                        sorted_dists.append(signame)
                    if not noData:
                        sorted_dists.insert(0, sys + "_data_obs")
                except:
                    print("Error: cat {}, var {} failed.".format(cat, var))
                    continue
                    
                #print("Integrals: {} \nsorted_dists: {}".format(integrals[sys], sorted_dists))
                    

                #repeat once for each different distribution.
                #for dnum,dist in enumerate(dists):
                for dnum,dist in enumerate(sorted_dists):
                    #print "divising MC into categories "
                    #hirBackground = ROOT.TH1F()

                    #print("sys: {}, dist: {}, cat: {}, var: {}".format(sys, dist, cat, var))
                    #hirBackground = histodict[sys][sys+"_"+dist][cat][var].Clone()
                    try:
                        hirBackground = histodict[sys][dist][cat][var].Clone()
                        #initialize hbkgerr if it hasn't already been initialized.
                    except:
                    #    print("var not working: ", var)
                        continue


                    #data

                    #set different colors for each distribution number.
                    colstr = ""
                    bkgtit = ""
                    print("dist = {}, dnum = {}".format(dist, dnum))
                    #Bkg
                    is_sig = False
                    is_data = False
                    if dist == sys + "_Bkg":
                        colstr = "#CF8AC8"
                        if args.fakeCat == "":
                            bkgtit = "DY+Jets" 
                        else:
                            bkgtit = "Jets faking taus"
#################                    #irBkg
                   # else: continue
                    elif dist == sys + "_irBkg":
                        if args.noZZ:
                            continue
                        colstr = "#13E2FE"
                        bkgtit = "ZZ or ZH to 4l"
                    #TrialphaBkg
                    elif dist == sys + "_TrialphaBkg":
                        colstr = "#65E114"
                        bkgtit = "3 #alpha"
                    #rareBkg
                    elif dist == sys + "_rareBkg":
                        colstr = "#FF6600"
                        bkgtit = "TTXX, WZ to Inv, 4#alpha"
                    elif dist == (sys + "_a%d"%mass):  #in ["%s_a%d"%(sys, mss) for mss in range(15, 61, 5)]:
                        is_sig = True
                        print("is_sig!!")
                        colstr = "#0000FF" #blue
                        bkgtit = "Signal " + dist
                    elif dist == (sys + "_data_obs"):
                        is_data = True
                       # colstr = "#000000" #black (?)
                        bkgtit = "Data observed"
                    else:
                    #    print("Error: still need to implement dnum = {}, dist = {}".format(dnum, dist))
                        system.exit() 
                        #break

                    if not is_sig and not is_data:
                    #background
                        hirBackground.SetLineColor(1)
                        hirBackground.SetFillStyle(1001)
                        hirBackground.SetFillColor(ROOT.TColor.GetColor(colstr))
                        hBkgTot.Add(hirBackground)
                        if (not noData and dnum == 1) or (noData and dnum == 0):
                            hbkgerr.Sumw2() #will this get the stat errors right? ans: no :(
                            hbkgerr = hirBackground.Clone()
                            #hbkgerr.Sumw2() #will this get the stat errors right? ans: no :(
                            if var == "AMass":
                                print("started with bkg " + str(hirBackground))
                                print("Starting bin 1 content: " + str(hbkgerr.GetBinContent(1)))
                                print("cumulative error on hbkgerr bin 1: " + str(hbkgerr.GetBinError(1)))
                        else:
                            hbkgerr.Add(hirBackground) #add this background onto the total
                            if var == "AMass":
                                print("added bkg " + str(hirBackground))
                                print("New bin 1 content: " + str(hbkgerr.GetBinContent(1)))
                                print("cumulative error on hbkgerr bin 1: " + str(hbkgerr.GetBinError(1)))
                    #    print("added to hBkgTot: {}".format(str(hirBackground)))
                    elif is_sig:
                    #signal
                        hirBackground.SetLineColor(ROOT.TColor.GetColor(colstr))
                        hirBackground.SetFillStyle(0)
                        hSig = hirBackground.Clone()
                    else:
                    #data
                        hirBackground.SetMarkerColor(1)
                        hirBackground.SetMarkerStyle(20)
                        hData = hirBackground.Clone()

                    hirBackground.SetTitle(bkgtit)

                    if not is_data:
                        l.AddEntry(hirBackground)
                    else:
                        l.AddEntry(hData, "Data observed", "PE")
                        
                    #print("Added {} to hBkgTot. Now hBkgTot = {}".format(hirBackground, hBkgTot))
                    #add entry to the legend.
                    


        
                ptitle="pad1_%s_%s"%(cat, var)
                if(doRatio):
                    pad1 = ROOT.TPad(ptitle,ptitle,0.0016,0.291,1.0,1.0)
                    pad1.SetTicks(0,0)
                    pad1.SetBottomMargin(B_ratio/H)
                    pad1.SetFillColor(0)
                    pad1.SetBottomMargin(0)

                else:
                    pad1 = ROOT.TPad(ptitle,ptitle,0,0,1,1)
                    pad1.SetBottomMargin   (B/H)

                pad1.SetLeftMargin(L/W)
                pad1.SetRightMargin(R/W)
                pad1.SetTopMargin(T/H)



                #print("boutta draw pad1")
                pad1.Draw()
                #print("Drew pad1.")
                pad1.cd()
                lumi=add_lumi(args.year)
                cms=add_CMS()
                pre=add_Preliminary(args.channel)

                #print "signal entries   ",hSignal.GetEntries()
                if not sigOnly:
                    hBkgTot.Draw("HIST")
                   # hBkgTot.Draw("PE")
                #set maximum to what it should be (idk why root can't do this automatically)
                    realmax = hBkgTot.GetMaximum()
                    if not noSig:
                        realmax = max(realmax, hSig.GetMaximum())
                    if not noData:
                        realmax = max(realmax, hData.GetMaximum())
                    hBkgTot.SetMaximum(1.1*realmax) 
                    if y_max > -1:
                    #MUST use SetMaximum and draw a second time for ROOT's stupid ass to understand
                        hBkgTot.SetMaximum(y_max) 
                        hBkgTot.Draw("HIST")
                    #draw MC error right after drawing the THStack.
                    hbkgerr.SetFillStyle(3013)
                    hbkgerr.SetFillColor(ROOT.TColor.GetColor("#263238"))
                    hbkgerr.Draw("same E2")  #hopefully this will show MC errors
                    if not noSig:
                        hSig.Draw("same HIST")
                    if not noData:
                        hData.Draw("PE same") # same")
                else:
                    hSig.Draw("HIST")
                
                try:
                    stack_title = variabledic[var][3]+variabledic[var][2]
                except:
                    print "variable not wokring", var
                    continue
                #print("hBkgTot = {}".format(str(hBkgTot)))
                #print("THStack x-axis: {}".format(hBkgTot.GetXaxis()))
                if not sigOnly:
                    hBkgTot.GetXaxis().SetTitle(stack_title)
                    hBkgTot.GetYaxis().SetTitle("Events")
                    hBkgTot.SetTitle("")
                #elif not noData:
                #    hData.GetXaxis().SetTitle(stack_title)
                #    hData.GetYaxis().SetTitle("Events")
                #    hData.SetTitle("")
                else:
                    hSig.GetXaxis().SetTitle(stack_title)
                    hSig.GetYaxis().SetTitle("Events")
                    hSig.SetTitle("")

                #print("boutta draw lumi")
                #hirBackground.Draw("same")
                lumi.Draw()
                cms.Draw()
                pre.Draw(args.channel)
                l.Draw()

                #TPad 2 for ratio
                c.cd()
                #print "with cuts ",allcats[cati].cuts
                #print "data entries ",hData.GetEntries()
                #print "background entries ",hBackground.GetEntries()
                if args.fakeCat == "":
                    fname = dirname+"/"+var+"_"+str(cats[0])+".png"
                else:
                    fname = dirname + "/" + var + "_" + args.fakeCat + ".png"
                print("trying to save as: {}".format(fname))
                try:
                    c.SaveAs(fname)
                except:
                    print("Error: could not save {}".format(fname))

                if not sigOnly and hirBackground:
                    hirBackground.Delete()
