# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDistributions_HAA_2016.py

* Purpose : For visualizing datasets. Will create ultra skimmed root files or histograms for visualization and for limit setting later

* Usage:
python MakeDistributions_HAA_2016.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /path/to/input/root/files/ -p processes_special_mmtt.yaml -o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt -s

* Creation Date : june-20-2020

* Last Modified : feb-16-2021

'''
#########################
import sys
import os
import ROOT
import uproot
import pandas
import root_numpy
import numpy as np
from array import array
import itertools
import operator
import csv
import datetime
import time
import threading
import yaml
import glob
#user definitions
from utils.Parametrization import *



#Setting up operators for cut string iterator
ops = { "==": operator.eq, "!=": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}

class myThread (threading.Thread):
   def __init__(self, threadID, name, counter):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.name = name
      self.counter = counter
   def run(self):
      print("Starting " + self.name)
      # Get lock to synchronize threads
      threadLock.acquire()
      print_time(self.name, self.counter, 3)
      # Free lock to release next thread
      threadLock.release()

def print_time(threadName, delay, counter):
   while counter:
      time.sleep(delay)
      print("%s: %s" % (threadName, time.ctime(time.time())))
      counter -= 1

def subArrayLogic(evt,array):
    boo=ops[array[1]](returnValue(evt,array[0]),array[2])
    return boo


def returnArray(masterArray,variable):
    #if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
    #    val = masterArray[variable][:,0]
    #else:
    if masterArray[variable].shape == 1:   # flat array important for jagged arrays of input data
        val = masterArray[variable]
    else:
        val = masterArray[variable][:,0]
    return val

#create mask for all the events that do or don't pass cuts
def cutOnArray(masterArray,cuts):

    mask=np.full(len(masterArray["evt"]),True)

    for cut in cuts:
        #Every USE CASE MUST MULTIPLY BY MASK!
        if type(cut[0][0])==type([0]):
            mask = cutOnArray(masterArray,cut)

        if cut[0][0]=="EQT":
            for i,var in enumerate(cut[1]):
                if cut[2]=="mult":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask = tempmask * returnArray(masterArray,var)

                if cut[2]=="div":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask = tempmask / returnArray(masterArray,var)
                if cut[2]=="add":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask = tempmask + returnArray(masterArray,var)
                if cut[2]=="sub":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask = tempmask - returnArray(masterArray,var)

            tempmask = ops[cut[3]](tempmask,cut[4])

            mask *= tempmask.astype(bool)
            continue

        if cut[0][0]=="OR":
            tempmasks=[]
            tempmask=np.full(len(masterArray["evt"]),False)
            for oe in range(1,len(cut)):
                oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
                tempmask = ops["bor"](tempmask,oneor).astype(bool)
            mask *= tempmask.astype(bool)
            continue

        else:
            statement = cut[1]
            if statement=="absg": # lessthan OR greater than
                tempmask = ops["<"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                tempmask2= ops[">"](returnArray(masterArray,cut[0]),cut[2])
                tempmask = ops["bor"](tempmask,tempmask2)
            if statement=="absl": # greater than AND less tahn
                tempmask = ops[">"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                tempmask *= ops["<"](returnArray(masterArray,cut[0]),cut[2])
            else:
                tempmask = ops[statement](returnArray(masterArray,cut[0]),cut[2])
            mask *= tempmask.astype(bool)


    return mask

#Depreciated but still useful later
#Function that combines neighboring bins into a list
def returnBinList(flatbins):
    comboBins=[]
    for sl in range(0,len(flatbins)):
        comboBins.append([])
        for bin in range(0,len(flatbins[sl])-1):
            comboBins[sl].append([flatbins[sl][bin],flatbins[sl][bin+1]])
    return comboBins

#Depreciated but still useful later
def returnPerms(setlist):

    combinations=list(itertools.product(*setlist))

    return combinations


#Depreciated but still useful later
def createUnrolledHistogram(cat,numvar,filelist,process,variable,weight,filedict):

    filedict[variable[0]].cd()
    filedict[variable[0]].mkdir(cat.name)
    filedict[variable[0]].cd(cat.name)

    #creating primary bins
    primaryBins = np.asarray(cat.binning[numvar][0]) #always the first binning set
    #creating secondary bin list
    secondaryBins = []

    for bins in range(1,len(cat.binning[numvar])):
        secondaryBins.append(cat.binning[numvar][bins])

    binlist = returnBinList(secondaryBins)
    combos = returnPerms(binlist)
    procut = HAA_ProCuts.get(process)
    histList = []
    urnum = 0
    for bincuts in combos:
        htmp=ROOT.TH1D("htemp"+str(process)+str(urnum),"htemp"+str(process)+str(urnum),len(primaryBins)-1,primaryBins)
        unrollCut=""

        for extraVarnum in range(1,len(variable)):
            unrollCut += "&&"+str(variable[extraVarnum])+">"+str(bincuts[extraVarnum-1][0])+"&&"+str(variable[extraVarnum])+"<"+str(bincuts[extraVarnum-1][1])
        if procut:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+"&&"+str(HAA_ProCuts[process])+unrollCut+")"
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
        else:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+unrollCut+")"
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)

        histList.append(htmp)
        htmp.Write(htmp.GetName(),ROOT.TObject.kOverwrite)
        urnum+=1
        htmp=0
    masterbins = cat.binning[numvar][0]*len(combos)
    masterUnroll = ROOT.TH1D(str(process)+"_unrolled",str(process)+"_unrolled",len(masterbins),0,len(masterbins))
    masterBinNum=0
    for hist in histList:
        for bin in range(0,hist.GetNbinsX()+1):
            masterUnroll.SetBinContent(masterBinNum,hist.GetBinContent(bin))
            masterUnroll.SetBinError(masterBinNum,hist.GetBinError(bin))
            masterBinNum+=1
    masterUnroll.Write(masterUnroll.GetName(),ROOT.TObject.kOverwrite)
    return

#for fake factor return the bin value from event content
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    newArr = np.vectorize(pt.FindBin)(numpyArr)
    newArr = np.vectorize(pt.GetBinContent)(newArr)

    return newArr

#Main function that does the skimming of input uproot files
def skimAndWeight(functs,tree,HAA_Inc_mmmt,cat,HAA_processes,process,
                    nickname,histodict,weightstring,commonweight,datadrivenPackage,
                    verbose,test,weightHisto,weightHistoDict,
                    jetWeightMultiplicity,finalDistributions,plottedVars):

    masterArray = tree.arrays()
    procut = processObj.cuts[process]
    cuts=[]
    for cuttype in allcats[cat].cuts.keys():
        for cut in allcats[cat].cuts[cuttype]:
            cuts.append(cut)
    if procut!="":
        cuts.append(procut[0])

    masterArray["finalweight"]=np.full(len(masterArray["evt"]),1.0)
    skimArray={}

    newVarVals={}
    for var in allcats[cat].newvariables.keys():
        newVarVals[var]=0.0


    if process=="data_obs":
        for var in newVarVals.keys():
            arguments = allcats[cat].newvariables[var][2]
            tempvals=[]
            for ag in arguments:
               tempvals.append(returnArray(masterArray,ag))
            masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)

        mask = cutOnArray(masterArray,cuts)
        masterArray["mask"]=mask
        masterArray["finalweight"] *= mask.astype(int)
        weightfinal = 1.0   #don't weight the data!!


        print("before skim", len(masterArray["finalweight"]))
        for key in masterArray.keys():
            skimArray[key] = masterArray[key][mask]
        print("after skim", len(skimArray["mll"]))



        for key in skimArray.keys():
            if key not in plottedVars and key != "finalweight":
                del skimArray[key]
        return skimArray

    if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
        for var in newVarVals.keys():
            arguments = allcats[cat].newvariables[var][2]
            tempvals=[]
            for ag in arguments:
               tempvals.append(returnArray(masterArray,ag))
            masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
        mask = cutOnArray(masterArray,cuts)
        masterArray["mask"]=mask
        masterArray["finalweight"] *= mask.astype(int)
        weightDict = processObj.weights
        weightfinal = commonweight
        for scalefactor in weightDict.keys():
            if scalefactor == "kfactor":
                weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
            elif scalefactor in ["PU"]:
                masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
            elif scalefactor =="theoryXsec":
                weightfinal =  weightfinal * float(weightDict[scalefactor])

        if nickname =="DY1JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY1JetsToLL"]/HAA_processes["DY1JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY2JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY2JetsToLL"]/HAA_processes["DY2JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY3JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY3JetsToLL"]/HAA_processes["DY3JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY4JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY4JetsToLL"]/HAA_processes["DY4JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="W1JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W1JetsToLNu"]/HAA_processes["W1JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        elif nickname =="W2JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W2JetsToLNu"]/HAA_processes["W2JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        elif nickname =="W3JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W3JetsToLNu"]/HAA_processes["W3JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)


        else:
            # simple multiplation like weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ weightHisto.GetSumOfWeights()
            # for cases with EXT in the name then we need to combine them
            tempDenom=1.0
            sumOfWeights = 1.0
            for nic in weightHistoDict.keys():
                if nickname==nic:
                    sumOfWeights = weightHisto.GetSumOfWeights()
                if nickname==nic.strip("_ext")[0]:
                    sumOfWeights*=tempDenom*weightHistoDict[nic].GetSumOfWeights()
            weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights

        #multiply by scalar weight
        masterArray["finalweight"] *= weightfinal





    if(process=="FF" and datadrivenPackage["bool"]):
        for var in newVarVals.keys():
            arguments = allcats[cat].newvariables[var][2]
            tempvals=[]
            for ag in arguments:
               tempvals.append(returnArray(masterArray,ag))
            masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
        tempmask=np.full(len(masterArray["evt"]),1.0)

        #the actual events that pass the FF_1 criteria
        tempmask_1 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_1"])
        #print tempmask_1[:1000]
        tempmask_2 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_2"])
        #print tempmask_2[:1000]
        tempmask_12 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_12"])
        #print tempmask_12[:1000]


        #FF_1
        #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
        fitmask_1 = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
        fitmask_1 = fitmask_1.astype(int)
        ptarr_1 = masterArray["pt_3"]

        ffweight_1 = ptFun(datadrivenPackage["fakerate1"],ptarr_1)
        ffweight_1 = ffweight_1/(1.0000000001 - ffweight_1)



        #FF_2
        fitmask_2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
        fitmask_2 = fitmask_2.astype(int)
        ptarr_2 = masterArray["pt_4"]

        ffweight_2 = ptFun(datadrivenPackage["fakerate2"],ptarr_2)
        ffweight_2 = ffweight_2/(1.0000000001 - ffweight_2)


        ffweight = ffweight_1 * ffweight_2

        #replace 0s with constant fit value
        ffweight_1 *= fitmask_1
        ffweight_1[np.where(ffweight_1==0)] =  datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0))
        ffweight_2 *= fitmask_2
        ffweight_2[np.where(ffweight_2==0)] =  datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))

        #FF_12
        #replace 0s with constant fit value
        ffweight[np.where(ffweight==0)] =  (datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0)))

        fitmask_1 *= fitmask_2
        ffweight *= fitmask_1
        ffweight *= tempmask_12

        ffweight_1[~tempmask_1] = 0.0
        #print "ffweight_1 ",ffweight_1[:1000]
        ffweight_2[~tempmask_2] = 0.0
        finalWeight = ffweight_1 + ffweight_2
        intersection = np.all((tempmask_1,tempmask_2), axis=0)

        finalWeight -= ffweight


        #print "check on fake factor final weight ",finalWeight[:1000]


        masterArray["finalweight"] *= finalWeight
        print "summed final weight ",np.sum(finalWeight)

        keepEvents = ~np.where(finalWeight==0.0)[0]

        skimArray={}
        for key in masterArray.keys():
            skimArray[key] = masterArray[key][keepEvents]


        for key in skimArray.keys():
            if key not in plottedVars and key != "finalweight":
                del skimArray[key]
        return skimArray

    if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
        eventWeightDict = processObj.eventWeights
        if eventWeightDict:
            for scalefactor in eventWeightDict.keys():

                cutlist = eventWeightDict[scalefactor][0]
                weightMask=cutOnArray(masterArray,cutlist)
                weightMask=weightMask.astype(float)

                if type(eventWeightDict[scalefactor][1][0])==float:
                    weightMask *= eventWeightDict[scalefactor][1][0]

                if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                    arguments = eventWeightDict[scalefactor][1][1]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))
                    weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)

                if scalefactor!="fake" and scalefactor!="fake1" and scalefactor!="fake2":
                    weightMask[np.where(weightMask==0.0)]=1.0
                else:
                    print "subtracting fakes "
                    print weightMask[:100]

                masterArray["finalweight"] *= weightMask
        skipEvents = np.where(mask==0)[0]
        skimArray={}
        print("before skim", len(masterArray["finalweight"]))
        for key in masterArray.keys():
            skimArray[key] = masterArray[key][mask]
        print("after skim", len(skimArray["mll"]))

        for key in skimArray.keys():
            if key not in plottedVars and key != "finalweight":
                del skimArray[key]
        return skimArray

    return 0

def calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,cat,HAA_processes,processObj,
                    nickname,histodict,weightstring,commonweight,datadrivenPackage,
                    verbose,test,weightHisto,weightHistoDict,jetWeightMultiplicity):


    newVarVals={}
    for var in allcats[cat].newvariables.keys():
        newVarVals[var]=0.0

    #usefl print block to check the cuts and processes
    #print("categories   ",allcats.keys())
    #print("process obj  ",processObj.nickname)
    #
    #for cat in allcats.keys():
    #    cuts = []
    #    for cuttype in allcats[cat].cuts.keys():
    #        for cut in allcats[cat].cuts[cuttype]:
    #            cuts.append(cut)
    #    print("category ",cat)
    #    print("cat cuts",cuts)
    #    print("process obj cuts ",processObj.cuts)
    #print("all cuts ",cuts)

    masterArray = tree.arrays()




    for cat in allcats.keys():

        for process in processObj.cuts.keys():
            procut = processObj.cuts[process]
            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)
            if procut!="":
                cuts.append(procut[0])

            masterArray["finalweight"]=np.full(len(masterArray["evt"]),1.0)

            if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
                for var in newVarVals.keys():
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]
                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))
                    masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
                mask = cutOnArray(masterArray,cuts)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)

            if process=="data_obs":
                for var in newVarVals.keys():
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]
                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))
                    masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
                mask = cutOnArray(masterArray,cuts)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!


            #print("working on ",process)

            if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
                weightDict = processObj.weights
                weightfinal = commonweight
                for scalefactor in weightDict.keys():
                    if scalefactor == "kfactor":
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    elif scalefactor in ["PU"]:
                        masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                    elif scalefactor =="theoryXsec":
                        weightfinal =  weightfinal * float(weightDict[scalefactor])

                if nickname =="DY1JetsToLL":
                    norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["DY1JetsToLL"]/HAA_processes["DY1JetsToLL"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                if nickname =="DY2JetsToLL":
                    norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["DY2JetsToLL"]/HAA_processes["DY2JetsToLL"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                if nickname =="DY3JetsToLL":
                    norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["DY3JetsToLL"]/HAA_processes["DY3JetsToLL"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                if nickname =="DY4JetsToLL":
                    norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["DY4JetsToLL"]/HAA_processes["DY4JetsToLL"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                if nickname =="W1JetsToLNu":
                    norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["W1JetsToLNu"]/HAA_processes["W1JetsToLNu"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                elif nickname =="W2JetsToLNu":
                    norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["W2JetsToLNu"]/HAA_processes["W2JetsToLNu"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                elif nickname =="W3JetsToLNu":
                    norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                    norm2 = jetWeightMultiplicity["W3JetsToLNu"]/HAA_processes["W3JetsToLNu"].weights["xsec"]
                    weightfinal = weightfinal * 1/(norm1+norm2)
                else:
                    tempDenom=1.0
                    sumOfWeights = 1.0
                    for nic in weightHistoDict.keys():
                        if nickname==nic.strip("_ext")[0]:
                            sumOfWeights*=tempDenom*weightHistoDict[nic].GetSumOfWeights()
                    weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights

                #multiply by scalar weight
                masterArray["finalweight"] *= weightfinal

            if(process=="FF" and datadrivenPackage["bool"]):
                for var in newVarVals.keys():
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]
                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))
                    masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
                tempmask=np.full(len(masterArray["evt"]),1.0)

                #the actual events that pass the FF_1 criteria
                tempmask_1 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_1"])
                tempmask_1 = tempmask.astype(int)
                tempmask_2 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_2"])
                tempmask_2 = tempmask.astype(int)
                tempmask_12 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_12"])
                tempmask_12 = tempmask.astype(int)

                #FF_1
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                fitmask_1 = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
                fitmask_1 = fitmask_1.astype(int)
                ptarr_1 = masterArray["pt_3"]

                #get pt binned value per event...
                ffweight_1 = ptFun(datadrivenPackage["fakerate1"],ptarr_1)
                ffweight_1 = ffweight_1/(1.0000000001 - ffweight_1)

                #replace 0s with constant fit value
                ffweight_1 *= fitmask_1
                ffweight_1[np.where(ffweight_1==0)] =  datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0))
                ffweight_1 *= tempmask_1

                #FF_2
                fitmask_2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
                fitmask_2 = fitmask_2.astype(int)
                ptarr_2 = masterArray["pt_4"]

                #get pt binned value per event...
                ffweight_2 = ptFun(datadrivenPackage["fakerate2"],ptarr_2)
                ffweight_2 = ffweight_2/(1.0000000001 - ffweight_2)

                #replace 0s with constant fit value
                ffweight_2 *= fitmask_2
                ffweight_2[np.where(ffweight_2==0)] =  datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))
                ffweight_2 *= tempmask_2

                ffweight_12 = ffweight_1 * ffweight_2
                ffweight_12 *= -1.0 #subtract events that fail both

                #concatenating the ffweights
                ffweight_12 *=tempmask_12

                masterArray["finalweight"] *= ffweight_12
                masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0

                for var in newVarVals.keys():
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))

                    root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


                for variableHandle in allcats[cat].vars.keys():

                    val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                    filedict[variableHandle].cd()
                    filedict[variableHandle].cd(allcats[cat].name)
                    root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])

                continue


            if process!="data_obs" and len(tree)!=0:
                eventWeightDict = processObj.eventWeights
                if eventWeightDict:
                    for scalefactor in eventWeightDict.keys():

                        cutlist = eventWeightDict[scalefactor][0]

                        weightMask=cutOnArray(masterArray,cutlist)
                        weightMask=weightMask.astype(float)

                        if type(eventWeightDict[scalefactor][1][0])==float:
                            weightMask *= eventWeightDict[scalefactor][1][0]

                        if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                            arguments = eventWeightDict[scalefactor][1][1]
                            tempvals=[]

                            for ag in arguments:
                               tempvals.append(returnArray(masterArray,ag))
                            weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)

                        #preseving the events that don't need new weight
                        if scalefactor!="fake":
                            weightMask[np.where(weightMask==0.0)]=1.0
                        masterArray["finalweight"] *= weightMask

            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]

                for ag in arguments:
                   tempvals.append(returnArray(masterArray,ag))

                root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


            for variableHandle in allcats[cat].vars.keys():
                val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                filedict[variableHandle].cd()
                filedict[variableHandle].cd(allcats[cat].name)
                root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])





    return

if __name__ == "__main__":

    begin_time = datetime.datetime.now()

    #Structure for mapping between root files and processes
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights
    from utils.Weights import HAAWeights

    sys.path.insert(1, '/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src')
    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
    #sys.path.append('../SFs/')
    import utils.SFs.ScaleFactor as SF

    import io
    import os
    from shutil import copyfile

    import argparse

    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
    parser.add_argument("-ch",  "--channel", default="mmmt",  help="Please list the channel for fake factor histograms")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/",  help="Input files")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
    parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
    parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
    parser.add_argument("-s",  "--skim", default=False,action='store_true',  help="skim input files to make more TTrees")
    parser.add_argument("-mt",  "--mt", default=False,action='store_true',  help="Use Multithreading")
    parser.add_argument("-pt",  "--maxprint", default=False,action='store_true',  help="Print Info on cats and processes")
    args = parser.parse_args()


    #gather functions for computing variables in the event loop
    from utils.functions import functs
    from ROOT import gInterpreter

    #Structure for plotting variables
    from utils.Parametrization import Category
    from utils.Parametrization import Process

    #Gather the analysis datasets and info
    sampleDict = {}

    for line in open(args.csvfile,'r').readlines() :
            # structure for csv to sampleDict conversion
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            if len(line.split(',')[2].split("*"))>1:
                tempval=1.0
                for val in line.split(',')[2].split("*"):
                    tempval = tempval * float(val)
                row = line.split(',')
                sampleDict[row[0]] = [row[1],tempval,row[3],row[4],row[5],row[6]]
            else:
                row = line.split(',')
                sampleDict[row[0]] = [row[1],float(row[2]),row[3],row[4],row[5],row[6]]

    #importing analysis categories and conventions
    with io.open(args.categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(args.processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print(category)
        #print(categories[category]['name'])
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat


    #loading standard processes
    HAA_processes={}

    #MC for background
    if not (args.datadrivenZH):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={sampleDict[sample][0]:""}
            if "ggTo2mu2tau" in sample:
                temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} # worked before
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001* 5.0)} # SM Higgs xsec x BR Haa  for signal extraction data MC control plots(AN 17029)
            if "W" in sample and "Jets" in sample:
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
                temppro.cuts={"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]}
            if "TT" in sample and not "TTTT" in sample and not "TTHH" in sample:
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
                temppro.cuts={"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]}
            HAA_processes[temppro.nickname]=temppro

    if (args.datadrivenZH or args.datameasureZH):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue"}
            if args.channel=="mmtt":
                truetau = [
                            [["OR"],
                            ["gen_match_3","==",5],
                            ["gen_match_4","==",5]
                            ]]
            if args.channel=="mmem":
                truetau = [ [["OR"],
                            ["gen_match_3","==",15]
                            ],
                            [["OR"],
                            ["gen_match_4","==",15]
                            ] ]
            if args.channel=="mmet":
                truetau = [ [["OR"],
                            ["gen_match_3","==",15]
                            ],
                            [
                            ["gen_match_4","!=",-1],
                            ["gen_match_4","!=",0]
                            ] ]
            if args.channel=="mmmt":
                truetau = [[["OR"],
                            ["gen_match_3","==",15]
                            ],
                            [
                            ["gen_match_4","!=",-1],
                            ["gen_match_4","!=",0]
                            ] ]
            temppro.cuts={sampleDict[sample][0]:truetau}
            if "ggTo2mu2tau" in sample:
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)}
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001*5.0)} # SM Higgs xsec x BR Haax100 x 5 for DataMC control plots
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haax100 x 5 for DataMC control plots - correct control plots??
                temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec x BR Haa (17-029 paper)
            #adding extra process splitting
            if "W" in sample and "Jets" in sample and "EWK" not in sample:
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","kfactor":1.221}
                temppro.cuts={"W":truetau,"WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]],"fake1_W":[["gen_match_3","==",0]],"fake2_W":[["gen_match_4","==",0]]}
            if "DY" in sample and "Jets" in sample:
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","kfactor":1.1637}
                temppro.cuts={"DY":truetau}
            if "TT" in sample and not "TTTT" in sample and not "TTHH" in sample:
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue"}
                temppro.cuts={"TT":truetau,"TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]],"fake1_TT":[["gen_match_3","==",0]],"fake2_TT":[["gen_match_4","==",0]]}
            HAA_processes[temppro.nickname]=temppro

    #loading special processes ... fake factor and data!
    for process in processes_special:
        temppro = Process()
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro

    if args.maxprint:
        print("categories ")
        for cat in allcats:
            print(cat.name)
            print(cat.cuts)


        print("processes ")
        for pro in HAA_processes.keys():
            print(HAA_processes[pro].nickname)
            print(HAA_processes[pro].cuts)
            print(HAA_processes[pro].weights)
            print(HAA_processes[pro].file)




    dir = args.dir
    filelist = {}

    cmsswbase="$CMSSW_BASE/src"
    #can add sv fit if needed
    #gInterpreter.ProcessLine(".include .")
    #for baseName in ['../SVFit/MeasuredTauLepton','../SVFit/svFitAuxFunctions','../SVFit/FastMTT', cmsswbase+'/HTT-utilities/RecoilCorrections/src/MEtSys', cmsswbase+'/HTT-utilities/RecoilCorrections/src/RecoilCorrector']:
    #    if os.path.isfile("{0:s}_cc.so".format(baseName)):
    #        gInterpreter.ProcessLine(".L {0:s}_cc.so".format(baseName))
    #    else:
    #        gInterpreter.ProcessLine(".L {0:s}.cc++".format(baseName))

    wpp = 'Medium'

    #tauSFTool = TauIDSFTool(campaign[args.year],'DeepTau2017v2p1VSjet',wpp)
    #testool = TauESTool(campaign[args.year],'DeepTau2017v2p1VSjet', TESSF['dir'])
    #festool = TauESTool(campaign[args.year],'DeepTau2017v2p1VSjet')
    ##antiEleSFToolVVL = TauIDSFTool(campaign[args.year],'DeepTau2017v2p1VSe','VVLoose')
    ##antiMuSFToolVL  = TauIDSFTool(campaign[args.year],'DeepTau2017v2p1VSmu','VLoose')
    #antiEleSFToolT = TauIDSFTool(campaign[args.year],'DeepTau2017v2p1VSe','Tight')
    #antiMuSFToolT  = TauIDSFTool(campaign[args.year],'DeepTau2017v2p1VSmu','Tight')

    #Sam
    tauIDSF = TauIDSFTool('2016Legacy','DeepTau2017v2p1VSjet','Medium').getSFvsPT


    #recoilCorrector  = ROOT.RecoilCorrector("HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root").CorrectByMeanResolution

    # muon scale factor as function of pt
    sf_MuonTrig = SF.SFs()
    sf_MuonTrig.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/TriggerEffs/Muon/SingleMuon_Run2016_IsoMu24orIsoMu27.root")
    sf_EleTrig = SF.SFs()
    sf_EleTrig.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/TriggerEffs/Electron/SingleElectron_Run2016_Ele25orEle27.root")
    sf_MuonId = SF.SFs()
    sf_MuonId.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/LeptonEffs/Muon/Muon_Run2016_IdIso_0p2.root")
    sf_ElectronId = SF.SFs()
    sf_ElectronId.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/LeptonEffs/Electron/Electron_Run2016_IdIso_0p1.root")
    EventWeights={
        #"name":[[if statements],[weight to apply]]
        "3_mt_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[0.80]],
        "3_mj_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[1.21]],
        "3_mt_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[0.81]],
        "3_mj_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[1.11]],
        "3_mt_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[0.79]],
        "3_mj_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[1.2]],
        "3_mt_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[0.68]],
        "3_mj_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[1.16]],
        "3_mt_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[0.68]],
        "3_mj_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[2.25]],

        "4_mt_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
        "4_mj_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
        "4_mt_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
        "4_mj_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
        "4_mt_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
        "4_mj_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
        "4_mt_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
        "4_mj_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
        "4_mt_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
        "4_mj_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

        "3_et_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
        "3_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
        "3_et_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
        "3_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
        "3_et_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
        "3_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
        "3_et_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
        "3_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],

        "4_et_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[0.80]],
        "4_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[1.18]],
        "4_et_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.72]],
        "4_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.93]],
        "4_et_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.14]],
        "4_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.18]],
        "4_et_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[0.64]],
        "4_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[1.07]],

        "8_3_et_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
        "8_3_ej_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
        "8_3_et_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
        "8_3_ej_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
        "8_3_et_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
        "8_3_ej_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
        "8_3_et_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
        "8_3_ej_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],

        "8_4_mt_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
        "8_4_mj_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
        "8_4_mt_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
        "8_4_mj_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
        "8_4_mt_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
        "8_4_mj_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
        "8_4_mt_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
        "8_4_mj_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
        "8_4_mt_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
        "8_4_mj_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

        #bound method ... last input list is func parameters! so this uses the tauIDSF function
        "3_tauSF_8":[[["cat","==",8],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_8":[[["cat","==",8],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_7":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_7":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_6":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_6":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_5":[[["cat","==",5],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_5":[[["cat","==",5],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]]

        #recoil corrections? not affecting current fit variable
        #"recoilcorr2":[[["cat","==",6],["njets","==",2]],[recoilCorrector,["met_x","met_y","gen_match_3"]]],
    }




    for objkey in HAA_processes.keys():
        if not objkey=="data":
            HAA_processes[objkey].eventWeights = EventWeights

    #for testing only... one cat a15 ... just to test
    #HAA_processes=HAA_processes_test

    if args.datameasure:
        filelist["data"]=HAA_processes["data"].file
        for category in allcats[:]:
            if "_inclusive" in category.name:
                allcats.remove(category)
    if not args.datameasure:
        for proObj in HAA_processes.keys():
           filelist[proObj]=HAA_processes[proObj].file


    if (args.datameasureZH):
        try:
            os.mkdir("FFhistos_"+str(args.ffout))
        except:
            print("directory exists")
        for proObj in HAA_processes.keys():
            if proObj!="data" and proObj!="FF" and not args.skim:
                if args.channel=="mmmt" or args.channel=="mmet":
                    #Make sure that these gen_match criteria match the ones from earlier for consistent prompt and fake definitions
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","!=",15]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [[["OR"],["gen_match_4","==",0],["gen_match_4","==",-1]]]
                    if "ggTo2mu2tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0],["gen_match_4","!=",-1]]
                if args.channel=="mmtt":
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [[["OR"],["gen_match_3","==",0],["gen_match_3","==",-1]]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [[["OR"],["gen_match_4","==",0],["gen_match_4","==",-1]]]
                    if "ggTo2mu2tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
                if args.channel=="mmem":
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","!=",15]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","!=",15]]
                    if "ggTo2mu2tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",15]]
                else: # default case is close enough
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                    if "ggTo2mu2tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]

        for category in allcats.keys():
            print(allcats[category].name)
            print("_inclusive" in allcats[category].name)
            if "_inclusive" in allcats[category].name:
                del allcats[category]
    elif (args.datadrivenZH):
        for proObj in HAA_processes.keys():
            if proObj!="data" and proObj!="FF" and not args.skim:
                HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                if "ggTo2mu2tau" in proObj:
                    print "skipping signal prompt mc"
                    continue  #this hasn't been added yet ... this is needed to omit signal!!
                HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]

    else:
        for category in allcats.keys():
            if "_inclusive" not in allcats[category].name:
                print("removing category ",allcats[category].name)
                del allcats[category]

    treename = "Events"

    print("initializing histograms")

    try:
        os.mkdir("out"+str(args.outname))
    except:
        print("directory exists")
    filedict = {}
    histodict = {}
    newhistodict = {}

    numvar=0
    if not args.skim:
        for variableHandle in allcats[allcats.keys()[0]].vars.keys():

            filedict[variableHandle]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variableHandle)+".root","RECREATE")
            filedict[variableHandle].cd()
            for cat in allcats.keys():
                if variableHandle not in allcats[cat].vars:
                    continue
                filedict[variableHandle].mkdir(allcats[cat].name)
                filedict[variableHandle].cd(allcats[cat].name)

                for nickname in filelist.keys():

                    processObj = HAA_processes[nickname]

                    for process in processObj.cuts.keys():


                        if filedict[variableHandle].Get(allcats[cat].name).GetListOfKeys().Contains(str(process)):
                            continue
                        else:
                            bins = allcats[cat].vars[variableHandle][1]
                            if type(bins[0])==list:
                                histodict[variableHandle+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),bins[0][0],bins[0][1],bins[0][2])
                                histodict[variableHandle+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                            else:
                                tmpbin = np.asarray(bins)
                                histodict[variableHandle+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                                histodict[variableHandle+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
            numvar=numvar+1

        #gathering new varibles
        numvar=0
        for variable in allcats[allcats.keys()[0]].newvariables.keys():

            filedict[variable]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variable)+".root","RECREATE")
            filedict[variable].cd()
            for cat in allcats.keys():
                filedict[variable].mkdir(allcats[cat].name)
                filedict[variable].cd(allcats[cat].name)

                for nickname in filelist.keys():

                    processObj = HAA_processes[nickname]

                    for process in processObj.cuts.keys():

                        if filedict[variable].Get(allcats[cat].name).GetListOfKeys().Contains(str(process)):
                            continue
                        else:
                            bins = allcats[allcats.keys()[0]].newvariables[variable][1]
                            tmpbin = np.asarray(bins)
                            if  allcats[allcats.keys()[0]].newvariables[variable][-1]=="multiob":
                                for var,ivar in enumerate(temp):
                                    newhistodict[variable+"_"+str(ivar)+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin[ivar])-1,tmpbin[ivar])
                                    newhistodict[variable+"_"+str(ivar)+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                            else:
                                newhistodict[variable+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                                newhistodict[variable+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
            numvar=numvar+1


    #Weights for luminosity
    weight = CommonWeights["lumi"][0]
    weightstring = CommonWeights["string"]


    #precise weights for jet multiplicity
    jetWeightMultiplicity = {}
    DYJetsFile = ROOT.TFile.Open(dir+filelist["DYJetsToLLext1"],"read")
    jetWeightMultiplicity["DYJetsToLLext1"]=DYJetsFile.Get("hWeights").GetSumOfWeights()
    DY1JetsFile = ROOT.TFile.Open(dir+filelist["DY1JetsToLL"],"read")
    jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
    DY2JetsFile = ROOT.TFile.Open(dir+filelist["DY2JetsToLL"],"read")
    jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
    DY3JetsFile = ROOT.TFile.Open(dir+filelist["DY3JetsToLL"],"read")
    jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
    DY4JetsFile = ROOT.TFile.Open(dir+filelist["DY4JetsToLL"],"read")
    jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

    WJetsFile = ROOT.TFile.Open(dir+filelist["WJetsToLNu"],"read")
    jetWeightMultiplicity["WJetsToLNu"]=WJetsFile.Get("hWeights").GetSumOfWeights()
    W1JetsFile = ROOT.TFile.Open(dir+filelist["W1JetsToLNu"],"read")
    jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()
    W2JetsFile = ROOT.TFile.Open(dir+filelist["W2JetsToLNu"],"read")
    jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
    W3JetsFile = ROOT.TFile.Open(dir+filelist["W3JetsToLNu"],"read")
    jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
    W4JetsFile = ROOT.TFile.Open(dir+filelist["W4JetsToLNu"],"read")
    jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

    #Calculate the scale factors and fill the histograms for the fake factor methods?
    datadrivenPackage={}
    datadrivenPackage["bool"]=False


    if args.datadrivenZH or args.makeFakeHistos:
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadrivenZH
        ff_file_3 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_3_ff.root","read")
        ff_file_4 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_4_ff.root","read")
        ss_1_tight = ff_file_3.Get(args.channel+"_FF_SS_1_tight/data_obs")
        ss_1_loose = ff_file_3.Get(args.channel+"_FF_SS_1_loose/data_obs")
        ss_2_tight = ff_file_4.Get(args.channel+"_FF_SS_2_tight/data_obs")
        ss_2_loose = ff_file_4.Get(args.channel+"_FF_SS_2_loose/data_obs")

        ss_1_tight_prompt = ff_file_3.Get(args.channel+"_FF_SS_1_tight/prompt1")
        ss_1_loose_prompt = ff_file_3.Get(args.channel+"_FF_SS_1_loose/prompt1")
        ss_2_tight_prompt = ff_file_4.Get(args.channel+"_FF_SS_2_tight/prompt2")
        ss_2_loose_prompt = ff_file_4.Get(args.channel+"_FF_SS_2_loose/prompt2")

        #subtracting prompt MC execpt for low stat channel
        if args.channel!="mmem":
            ss_1_tight.Add(ss_1_tight_prompt,-1)
            ss_2_tight.Add(ss_2_tight_prompt,-1)
            ss_1_loose.Add(ss_1_loose_prompt,-1)
            ss_2_loose.Add(ss_2_loose_prompt,-1)


        f_1= ss_1_tight.Clone()
        f_1.Divide(ss_1_loose)
        f_2 = ss_2_tight.Clone()
        f_2.Divide(ss_2_loose)

        f_1.SetName(args.channel+" FakeRateLeg1")
        f_1.GetXaxis().SetTitle("p_T Leg1")
        f_1.GetYaxis().SetTitle("Fake Rate for Leg1")
        f_2.SetName(args.channel+" FakeRateLeg2")
        f_2.GetXaxis().SetTitle("p_T Leg2")
        f_2.GetYaxis().SetTitle("Fake Rate for Leg2")

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())
        tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax())

        fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/fakemeasure.root","RECREATE")
        fakemeasurefile.cd()
        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptFit()
        f_1.Draw()
        f_1.Fit("tf_1")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate1.png")
        f_2.Draw()
        f_2.Fit("tf_2")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate2.png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
        f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)
        datadrivenPackage["fakerate1"]=f_1
        datadrivenPackage["fakerate2"]=f_2
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fitrate2"]=tf_2
        fakemeasurefile.Close()

    if args.skim:
        Bkg = ["DY","W","TT","ST","EWK"]
        irBkg = ["ZZ","ZHToTauTau","vbf","WHTT"]
        TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ"]
        rareBkg = ["Other","rare","WZ"]
        finalDistributions = {}
        finalDistributions["Bkg"]=Bkg
        finalDistributions["data_obs"]=["data_obs"]
        finalDistributions["a15"]=["a15"]
        finalDistributions["a20"]=["a20"]
        finalDistributions["a25"]=["a25"]
        finalDistributions["a30"]=["a30"]
        finalDistributions["a35"]=["a35"]
        finalDistributions["a40"]=["a40"]
        finalDistributions["a45"]=["a45"]
        finalDistributions["a50"]=["a50"]
        finalDistributions["a55"]=["a55"]
        finalDistributions["a60"]=["a60"]
        finalDistributions["irBkg"]=irBkg
        finalDistributions["TrialphaBkg"]=TrialphaBkg
        finalDistributions["rareBkg"]=rareBkg
        if args.datadrivenZH:
            Bkg = ["FF"]
            finalDistributions["Bkg"]=Bkg

        allSkims = {}
        finalSkims ={}


        for cat in allcats.keys():
            finalSkims[cat]={}

    weightHistoDict = {}
    for nickname in filelist.keys():
        try:
            frooin = ROOT.TFile.Open(dir+filelist[nickname],"read")
            weightHistoDict[nickname]=frooin.Get("hWeights")
        except:
            print nickname," hWeights prob doesn't exist?"

    ###MAIN LOOP
    for nickname in filelist.keys():
        print("working on file ",filelist[nickname])
        frootin = ROOT.TFile.Open(dir+filelist[nickname],"read")
        fin = uproot.open(dir+filelist[nickname])
        tree = fin[treename]
        weightHisto = frootin.Get("hWeights")
        #newtree=tree.CloneTree()
        print("tree entries ",len(tree))
        if len(tree)==0:
            print("tree length 0 ... skipping")
            continue

        processObj = HAA_processes[nickname]
        if args.skim:
            print "all categories ",allcats.keys()
            for cat in allcats.keys():
                print "working on cat ",cat
                plottedVars = []
                processSkims = {}

                for var in allcats[cat].vars.keys():
                    plottedVars.append(allcats[cat].vars[var][0])

                for newvar in allcats[cat].newvariables.keys():
                    plottedVars.append(newvar) # new vars just label

                for process in processObj.cuts.keys():
                    print("process ",process)
                    if process not in processSkims:
                        processSkims[process] = skimAndWeight(functs, tree, allcats[allcats.keys()[0]], cat, HAA_processes, process,
                                        nickname, histodict, weightstring, weight, datadrivenPackage, args.verbose,
                                        args.test, weightHisto,weightHistoDict, jetWeightMultiplicity,finalDistributions,plottedVars)
                    else: # add up the trees for each process
                        temp = skimAndWeight(functs, tree, allcats[allcats.keys()[0]], cat, HAA_processes, process,
                                        nickname, histodict, weightstring, weight, datadrivenPackage, args.verbose,
                                        args.test, weightHisto,weightHistoDict, jetWeightMultiplicity,finalDistributions,plottedVars)
                        for key in processSkims[process].keys():
                                    processSkims[process][key]=np.concatenate((processSkims[process][key],temp[key]))
                #print "found the following distributions ... ",processSkims.keys()

                #now that we have all the process skims we need to combine them for the actual groups
                print("grouping processes",nickname)
                for catDist in finalDistributions.keys():
                    for process in finalDistributions[catDist]:
                        print "looking for process ",process
                    #for process in processSkims.keys():
                        if (process in processSkims.keys()) and (catDist not in finalSkims[cat]):
                            print "first ",process # this happens twice ... which is bad...
                            finalSkims[cat][catDist] = processSkims[process]
                            continue
                        elif (process in processSkims.keys()) and (catDist in finalSkims[cat]):
                            print "adding ",process
                            for key in finalSkims[cat][catDist].keys():
                                finalSkims[cat][catDist][key]=np.concatenate((finalSkims[cat][catDist][key],processSkims[process][key]))
                        else:
                            continue
        else: #uising this for data measure
            calculateHistos(functs, tree, allcats[allcats.keys()[0]], allcats, cat, HAA_processes, processObj,
                            nickname, histodict, weightstring, weight, datadrivenPackage, args.verbose,
                            args.test, weightHisto,weightHistoDict, jetWeightMultiplicity)

    if not args.skim:
        print("writing the histograms to the file")
        for varcatPro in histodict.keys():
            filedict[varcatPro.split(":")[0]].cd()
            filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])

            histodict[varcatPro].Write(histodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
        print("writing the new variable histograms to the file")
        for varcatPro in newhistodict.keys():
            filedict[varcatPro.split(":")[0]].cd()
            filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
            newhistodict[varcatPro].Write(newhistodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)

        if args.datameasureZH:
            for file in glob.glob("out"+str(args.outname)+"/pt*.root"):
                print("copying over ",file)
                copyfile(file, "FFhistos_"+str(args.ffout)+"/"+file.split("/")[1])
    else:
        for cat in allcats.keys():
            skimFile = ROOT.TFile("skimmed_"+args.outname+"_"+cat+".root","recreate")
            skimFile.cd()

            dataTypes =[[],[]]
            print finalSkims[cat].values()
            random_sample = finalSkims[cat].values()[0]
            for branch in random_sample.keys():
                dataTypes[0].append(branch)
                dataTypes[1].append(random_sample[branch].dtype)
            for catDist in finalSkims[cat].keys():
                #print "on the final dist ",catDist
                data = np.zeros(len(finalSkims[cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
                for branch in data.dtype.names:
                    #print "working on branch ",branch
                    #print "branch datatype ",type(finalSkims[cat][catDist][branch])
                    if len(finalSkims[cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
                        data[branch] = finalSkims[cat][catDist][branch]
                    else:
                        data[branch] = finalSkims[cat][catDist][branch][:,0]
                treeOut = root_numpy.array2tree(data, name=catDist)
                treeOut.Write()
            skimFile.Close()



    del sf_MuonId
    del sf_MuonTrig
    del sf_EleTrig
    del sf_ElectronId
    print("computation time")
    print(datetime.datetime.now() - begin_time)
