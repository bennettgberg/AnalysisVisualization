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

def add_lumi(year):
    lowX=0.65
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    lumi.AddText(str(year)+" 35.9 fb^{-1} (13 TeV)")
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
    lowX=0.45
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
    parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
    parser.add_argument("-o",  "--output", default="",  help="postfix string")
    parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
    parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
    parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
    parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
    parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
    args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)




    #importing analysis categories and conventions
    with io.open(args.categories,'r') as catsyam:
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

    for cat in allcats.keys():
        print "working on cat",cat
        for newvar in  allcats[cat].newvariables:
            newvars.append([newvar])
            variabledic[newvar]=[newvar,allcats[cat].newvariables[newvar][1],allcats[cat].newvariables[newvar][3],allcats[cat].newvariables[newvar][4]]
        variables = []
        for varHandle in allcats[cat].vars.keys():
            variables.append([varHandle])
            variabledic[varHandle]=allcats[cat].vars[varHandle]

        #fin = uproot.open(args.input)
        print "opening skimmed_"+args.input+"_"+cat+".root"
        fin = uproot.open("skimmed_"+args.input+"_"+str(cat)+".root")


        histolist = {}
        finalhists = {}
        processes = []
        histodict = {}
        dists = fin.keys()
        vars = allcats[cat].vars.keys()
        print "variables ",vars

        #for ivar,var in enumerate(cat.vars.keys()):
        #for dist in fin.GetListOfKeys():
        for distLong in dists:
            dist = distLong.split(";")[0]
            #distribution = dist.ReadObj()  #The TTree
            tree = fin[dist]
            masterArray = tree.arrays()
            histodict[dist]={}
            histodict[dist][cat]={}
            for variableHandle in vars:
                print variableHandle
                var = allcats[cat].vars[variableHandle][0]
                bins = allcats[cat].vars[variableHandle][1]
                if type(bins[0])==list:
                    histodict[dist][cat][variableHandle] = ROOT.TH1D(str(dist),str(dist),bins[0][0],bins[0][1],bins[0][2])
                    val = masterArray[var]
                    root_numpy.fill_hist(histodict[dist][cat][variableHandle],val,masterArray["finalweight"])
                else:
                    tmpbin = np.asarray(bins)
                    histodict[dist][cat][variableHandle] = ROOT.TH1D(str(dist),str(dist),len(tmpbin)-1,tmpbin)
                    val = masterArray[var]
                    root_numpy.fill_hist(histodict[dist][cat][variableHandle],val,masterArray["finalweight"])
                #histodict[variablehandle+":"+allcats[cat].name+":"+dist] = root.TH1D(str(process),str(process),bins[0][0],bins[0][1],bins[0][2])
                #root_numpy.fill_hist(histodict[variablehandle+":"+allcats[cat].name+":"+dist],val,masterArray["finalweight"])



        #for varCatDist in histodict.keys():
        try:
            os.mkdir("outplots_"+args.output)
        except:
            print "dir prob exists"

        #var is now the variable handle
        for var in vars:
            if args.mc:
                if var=="AMass":
                    fileout = open("outplots_"+args.output+"/"+str(allcats[cat].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cat].name+"\n")
            else:
                if var=="AMass_blinded":
                    fileout = open("outplots_"+args.output+"/"+str(allcats[cat].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cat].name+"\n")
                #if var in ["AMass","mll_m15","mll_m20","mll_m25","mll_m30","mll_m35","mll_m40","mll_m45","mll_m50","mll_m55","mll_m60"]:
                #    continue


            #signal
            hSignal = histodict["a40"][cat][var]
            hSignals={}

            hDataDict={}
            hMCDict={}
            hMCFake1Dict={}
            hMCFake2Dict={}

            #print "divising MC into categories "
            hBackground = ROOT.TH1F()
            hirBackground = ROOT.TH1F()
            hrareBackground = ROOT.TH1F()
            h3alphBackground = ROOT.TH1F()

            hBackground = histodict["Bkg"][cat][var].Clone()
            hirBackground = histodict["irBkg"][cat][var].Clone()
            h3alphBackground = histodict["TrialphaBkg"][cat][var].Clone()
            hrareBackground = histodict["rareBkg"][cat][var].Clone()


            #data

            histolist=[]

            hBackground.SetLineColor(1)
            hirBackground.SetLineColor(1)
            hrareBackground.SetLineColor(1)
            hBackground.SetFillStyle(1001)
            hirBackground.SetFillStyle(1001)
            hrareBackground.SetFillStyle(1001)
            hBackground.SetFillColor(ROOT.TColor.GetColor("#CF8AC8"))
            hirBackground.SetFillColor(ROOT.TColor.GetColor("#13E2FE"))
            #h3alphBackground.SetFillColor(ROOT.TColor.GetColor("#65E114"))
            hrareBackground.SetFillColor(ROOT.TColor.GetColor("#65E114"))
            #hrareBackground.SetFillColor(ROOT.TColor.GetColor("#FF6600"))
            h3alphBackground.SetFillColor(ROOT.TColor.GetColor("#FF6600"))

            #hSignal.SetTitle()
            hBackground.SetTitle("DY+Jets, T#bar{T}, ST , EWK")
            #hirBackground.SetTitle("irreducible")
            hirBackground.SetTitle("ZZ or ZH to 4l")
            h3alphBackground.SetTitle("3 #alpha")
            hrareBackground.SetTitle("TTXX, WZ to Inv, 4#alpha")

            hBkgTot = ROOT.THStack()
            hBkgTot.Add(hirBackground)
            hBkgTot.Add(hrareBackground)
            hBkgTot.Add(h3alphBackground)

            if args.datadrivenZH:
                hData = histodict["data_obs"][cat][var]
                #hPrompt = histodict["prompt"][cat][var]
                #hFake1 = histodict["fake1")
                #hFake2 = histodict["fake2")
                #hirBackground.Add(hFake1,-1)
                #hirBackground.Add(hFake2,-1)
                hFF = histodict["Bkg"][cat][var]
                hFF.SetTitle("Jet faking #tau")
                #hFF.Add(hFF2)
                #hFF.Add(hFF12,-1)
                #hFF.Add(hPrompt,-1)
                #hFF.Add(hFF12,-1)
                #hFF.Scale(0.5)
                hFF.SetLineColor(1)
                hFF.SetFillStyle(1001)
                hFF.SetFillColor(ROOT.TColor.GetColor("#CF8AC8"))

            if not args.mc:
                #hData = histodict["data_obs"][cat][var].Clone()
                dataMax = hData.GetMaximum()*1.35
                dataMin = hData.GetMinimum()*1.25
                if var in ["mll_fine"]:
                    hData.SetBinContent(7,-999)
                    hData.SetBinContent(8,-999)
                    hData.SetBinContent(9,-999)
                    hData.SetBinContent(10,-999)
                    hData.SetBinContent(11,-999)
                if var in ["mll"]:
                    hData.SetBinContent(3,-999)
                    hData.SetBinContent(4,-999)
                    hData.SetBinContent(5,-999)
                    hData.SetBinContent(6,-999)
                if var in ["AMass"]:
                    hData.SetBinContent(0,-999)
                    hData.SetBinContent(1,-999)
                    hData.SetBinContent(2,-999)
                    hData.SetBinContent(3,-999)
                    hData.SetBinContent(4,-999)
                    hData.SetBinContent(5,-999)
                    hData.SetBinContent(6,-999)
                    hData.SetBinContent(7,-999)



            if not (args.datadriven or args.datadrivenZH):
                hbkg = hBackground.Clone()
                hBkgTot.Add(hBackground)
            else:
                hbkg = hFF.Clone()
                hBkgTot.Add(hFF)


            hbkg.Add(hirBackground)
            hbkg.Add(hrareBackground)
            hbkg.Add(h3alphBackground)




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

            c=ROOT.TCanvas("canvas","",0,0,600,600)

            doRatio=True

            if not doRatio:
                c.SetLeftMargin(L/W)
                c.SetRightMargin(R/W)
                c.SetTopMargin(T/H)
                c.SetBottomMargin(B/H)
            c.cd()

            pad1 = ROOT.TPad("pad1","pad1",0.0016,0.291,1.0,1.0)
            #pad2 = ROOT.TPad("pad2","pad2",0.45,0.35,0.92,0.95)
            pad2 = ROOT.TPad("pad2","The lower pad",0,0,1.0,0.29)

            if(doRatio):
                pad1.SetTicks(0,0)
                pad1.SetLeftMargin(L/W)
                pad1.SetRightMargin(R/W)
                pad1.SetTopMargin(T/H)
                pad1.SetBottomMargin(B_ratio/H)
                pad1.SetFillColor(0)
                pad1.SetBottomMargin(0)

                pad2.SetLeftMargin(L/W)
                pad2.SetRightMargin(R/W)
                pad2.SetTopMargin(T_ratio/H)
                pad2.SetTopMargin(0.007)
                pad2.SetBottomMargin(B_ratio_label/H)
                pad2.SetGridy(1)
                pad2.SetFillColor(4000)

            else:
                pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
                pad1.SetLeftMargin     (L/W)
                pad1.SetRightMargin    (R/W)
                pad1.SetTopMargin      (T/H)
                pad1.SetBottomMargin   (B/H)



            pad1.Draw()
            pad1.cd()
            lumi=add_lumi("2016")
            cms=add_CMS()
            pre=add_Preliminary(args.channel)
            xR=0.65
            l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
            if var=="mll_fine":
                l=ROOT.TLegend(0.40,0.55,0.40+0.28,0.9);
            if not args.mc:
                l.AddEntry(hData)
            l.AddEntry(hSignal)
            if not (args.datadriven or args.datadrivenZH):
                l.AddEntry(hBackground)
            else:
                l.AddEntry(hFF)
            l.AddEntry(h3alphBackground)
            l.AddEntry(hirBackground)
            l.AddEntry(hrareBackground)

            #print "signal entries   ",hSignal.GetEntries()
            hSignal.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
            hSignal.GetYaxis().SetTitle("Events")

            if not args.mc:
                hData.SetTitle("")
                hData.SetMarkerStyle(8)
                hData.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
                hData.GetYaxis().SetTitle("Events")
                hData.GetYaxis().SetRangeUser(dataMin,dataMax)
                if var in ["mll_m40","mll","mll_fine","mll_fit0","mll_fit1","mll_fit2","mll_fit3"]:
                    if(hBkgTot.GetMaximum() > hSignal.GetMaximum()):
                        hBkgTot.Draw("HIST")
                        hBkgTot.SetTitle("")
                        hBkgTot.GetXaxis().SetTitle(str(variabledic[var][3])+str(variabledic[var][2]))
                        hBkgTot.GetYaxis().SetTitle("Events")
                        hSignal.Draw("HIST same axis")
                    else:
                        hSignal.Draw("HIST")
                        hSignal.SetTitle("")
                        hSignal.GetXaxis().SetTitle(str(variabledic[var][3])+str(variabledic[var][2]))
                        hSignal.GetYaxis().SetTitle("Events")
                        hBkgTot.Draw("HIST same")
                else:
                    hData.Draw("ep")
                    hBkgTot.Draw("HIST same")
                    hSignal.Draw("HIST same")
                    hData.Draw("ep same")
            else:
                if(hBkgTot.GetMaximum() > hSignal.GetMaximum()):
                    hBkgTot.Draw("HIST")
                    hBkgTot.GetYaxis().SetRangeUser(hBkgTot.GetMinimum()*1.25,hBkgTot.GetMaximum()*1.25)
                    hSignal.Draw("HIST same")
                else:
                    hSignal.SetTitle("")
                    hSignal.Draw("HIST")
                    hBkgTot.Draw("HIST same")
            #hBackground.Draw("same")
            #hirBackground.Draw("same")
            #for hist in histolist:
            #    hist.Draw("same")
            lumi.Draw()
            cms.Draw()
            pre.Draw(args.channel)
            l.Draw()

            #TPad 2 for ratio
            c.cd()
            if not args.mc:
                hData2=hData.Clone()
                pad2.Draw()
                pad2.cd()
                hData2.Divide(hbkg)
                hData2.SetMarkerStyle(20)
                hData2.SetTitleSize  (0.12,"Y")
                hData2.SetTitleOffset(0.40,"Y")
                hData2.SetTitleSize  (0.12,"X")
                hData2.SetLabelSize  (0.10,"X")
                hData2.SetLabelSize  (0.08,"Y")
                #hData2.GetYaxis().SetRangeUser(0.62,1.38)
                hData2.GetYaxis().SetRangeUser(0.3,2.2)
                hData2.GetYaxis().SetNdivisions(305)
                hData2.GetYaxis().SetTitle("Obs/Exp   ")
                if var in ["mll_m40","mll","mll_fine","mll_fit0","mll_fit1","mll_fit2","mll_fit3"]:
                    hSignal2=hSignal.Clone()
                    #hSignal2.Divide(hbkg)
                    #hSignal2.SetMarkerStyle(20)
                    hSignal2.SetTitleSize  (0.12,"Y")
                    hSignal2.SetTitleOffset(0.40,"Y")
                    hSignal2.SetTitleSize  (0.12,"X")
                    hSignal2.SetLabelSize  (0.10,"X")
                    hSignal2.SetLabelSize  (0.08,"Y")
                    #hSignal2.GetYaxis().SetRangeUser(0.62,1.38)
                    #hSignal2.GetYaxis().SetRangeUser(0.3,2.2)
                    hSignal2.GetYaxis().SetNdivisions(305)
                    hSignal2.GetYaxis().SetTitle("Signal   ")
                    hSignal2.Draw("P")
                else:
                    hData2.Draw("P");
                pad2.Draw()
            else:
                hSignal2=hSignal.Clone()
                pad2.Draw()
                pad2.cd()
                hSignal2.Divide(hbkg)
                hSignal2.SetMarkerStyle(20)
                hSignal2.SetTitleSize  (0.12,"Y")
                hSignal2.SetTitleOffset(0.40,"Y")
                hSignal2.SetTitleSize  (0.12,"X")
                hSignal2.SetLabelSize  (0.10,"X")
                hSignal2.SetLabelSize  (0.08,"Y")
                #hSignal2.GetYaxis().SetRangeUser(0.62,1.38)
                hSignal2.GetYaxis().SetRangeUser(0.3,2.2)
                hSignal2.GetYaxis().SetNdivisions(305)
                hSignal2.GetYaxis().SetTitle("Sig/Bkg   ")
                hSignal2.Draw("P");
                pad2.Draw()

            if args.mc:
                if var=="AMass":
                    for key in allcats[cat].cuts.keys():
                        for cut in allcats[cat].cuts[key]:
                            fileout.write(str(cut))
                        fileout.write("\n")
                    fileout.write("signal entries "+str(hSignal.Integral())+"\n")
                    fileout.write("background entries "+str(hBackground.Integral())+"\n")
            else:
                if var=="AMass_blinded":
                    for key in allcats[cat].cuts.keys():
                        for cut in allcats[cat].cuts[key]:
                            fileout.write(str(cut))
                        fileout.write("\n")
                    fileout.write("data entries "+str(hData.Integral())+"\n")
                    fileout.write("background entries "+str(hbkg.Integral())+"\n")

            if args.mhs:
                if var in ["mll_m40","mll","mll_fine","mll_fit0","mll_fit1","mll_fit2","mll_fit3"]:
                    histoout = ROOT.TFile.Open("final_"+str(cat)+"_"+str(var)+".root","recreate")
                    histoout.mkdir(str(cat))
                    histoout.cd(str(cat))
                    if args.mc:
                        hBackground.Write(hBackground.GetName(),ROOT.TObject.kOverwrite)
                    else:
                        hFF.Write(hFF.GetName(),ROOT.TObject.kOverwrite)
                    hirBackground.Write(hirBackground.GetName(),ROOT.TObject.kOverwrite)
                    h3alphBackground.Write(h3alphBackground.GetName(),ROOT.TObject.kOverwrite)
                    hrareBackground.Write(hrareBackground.GetName(),ROOT.TObject.kOverwrite)
                    hData.Write(hData.GetName(),ROOT.TObject.kOverwrite)
                    hSignal.Scale(args.signalScale)
                    hSignal.Write(hSignal.GetName(),ROOT.TObject.kOverwrite)
                    histoout.Close()

            c.SaveAs("outplots_"+args.output+"/"+var+"_"+str(cat)+".png")

            hBackground.Delete()
            hirBackground.Delete()
            h3alphBackground.Delete()
            if not args.mc:
                hData.Delete()
                hData2.Delete()
            else:
                hSignal2.Delete()
            hSignal.Delete()
            hBkgTot.Delete()
            hbkg.Delete()
