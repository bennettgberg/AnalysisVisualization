# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDataCards_HAA_dataframe.py

* Purpose : For Histogramming. Will create histograms for visualization and for limit setting later

* Creation Date : june-20-2020

* Last Modified :

'''
#########################
import sys
import os
import ROOT
import uproot
import pandas
import root_numpy
import numpy as np
#import argsparse
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

#GLOBAL VARIABLES AND FUNCTIONS!


#Setting up operators for cut string iterator ... very useful
ops = { "==": operator.eq, "!=": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}

class myThread (threading.Thread):
   def __init__(self, threadID, name, counter):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.name = name
      self.counter = counter
   def run(self):
      print "Starting " + self.name
      # Get lock to synchronize threads
      threadLock.acquire()
      print_time(self.name, self.counter, 3)
      # Free lock to release next thread
      threadLock.release()

def print_time(threadName, delay, counter):
   while counter:
      time.sleep(delay)
      print "%s: %s" % (threadName, time.ctime(time.time()))
      counter -= 1

def subArrayLogic(evt,array):
    boo=ops[array[1]](returnValue(evt,array[0]),array[2])
    return boo

#most used function... 
def returnValue(evt,variable):
    if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
        val = getattr(evt,variable,None)[0]
    else:
        val = getattr(evt,variable,None)
    return val

def cutStringBool(evt,cuts):
    survive=False
    #print "cuts",cuts
    for cut in cuts:

        #recursive case ... keep er going...  
        if type(cut[0][0])==type([0]):
            sruvive=cutStringBool(evt,cuts)

        #print "single cut",cut,"  len ",len(cut)
        #Here there are cut types ... everything will be product ands... but these are the subtypes 
        # "OR" so we can AND(OR)AND(OR) 
        if cut[0][0]=="OR":
            tempbool=False
            for oe in range(1,len(cut)):
                #print cut[oe]
                #if ops[cut[oe][1]](returnValue(evt,cut[oe][0]),cut[oe][2]):
                if ops[cut[oe][1]](returnValue(evt,cut[oe][0]),cut[oe][2]):
                    tempbool=True
                    break
                    #print "passed trigger "
                
            survive=tempbool    
            if not survive:
                break
            else: continue
        # "EQT" does a particular set of variables when applied in a function pass a cut?
        if cut[0][0]=="EQT":
            for i,var in enumerate(cut[1]):
                if cut[2]=="mult":
                    if i==0:
                        tempvar = 1.0
                    tempvar = tempvar * returnValue(evt,var)
                if cut[2]=="div":
                    if(i==0 and returnValue(evt,var)!=0):
                        tempvar = 1.0
                    tempvar = tempvar / returnValue(evt,var)
                if cut[2]=="add":
                    if(i==0):
                        tempvar = 0.0
                    tempvar = tempvar+returnValue(evt,var)
                if cut[2]=="sub":
                    if(i==0):
                        tempvar = 0.0
                    tempvar = tempvar - returnValue(evt,var)

            survive = ops[cut[3]](tempvar,cut[4])
            #survive=True
            if not survive:
                #print cut[0],"     ",cut[1][0],"  ",returnValue(evt,cut[1][0]),"   ",cut[1][1],"   ",returnValue(evt,cut[1][1]),"  multicharge  ",tempvar,"   bool   ",ops[cut[3]](tempvar,cut[4])
                break
            else: continue 
        # "IF" statments in order to make cuts after variables in an event pass a selection criteria ... "THEN" included only for asethetics 
        if cut[0][0]=="IF":
            tempbool=True
            for i,state in enumerate(cut[1]):
                #if not ops[state[1]](returnValue(evt,state[0]),state[2]):
                #print state,ops[state[1]](returnValue(evt,state[0]),state[2])
                if ops[state[1]](returnValue(evt,state[0]),state[2]):
                    tempbool=True
                else:
                    tempbool=False
                    break
            if tempbool==True:
                for i,state in enumerate(cut[3]):
                    #print state,ops[state[1]](returnValue(evt,state[0]),state[2]),returnValue(evt,state[0])
                    if ops[state[1]](returnValue(evt,state[0]),state[2]):
                        survive=True 
                    else:
                        survive=False 
                        break
            else:
                #print "no if statement soo no cuts"
                continue 

            if not survive:
                break
            else: continue 

        # Default CASE! does it pass the cut?
        # should  provide this as separate function??? because the absg does exist above... 
        else:
            statement = cut[1]
            #print statement
            if statement=="absg":
                if ops["<"](returnValue(evt,cut[0]),(-1 * cut[2])): 
                    survive=True
                if ops[">"](returnValue(evt,cut[0]),cut[2]): 
                    survive=True 
                else: 
                    survive=False
                    break
                    #print "failed ",cut
            if statement=="absl":
                if (ops[">"](returnValue(evt,cut[0]),(-1 * cut[2])) and ops["<"](returnValue(evt,cut[0]),cut[2])): 
                    survive=True
                else: 
                    survive=False
                    break
                    
                    #print "failed ",cut
            else:  
                if ops[statement](returnValue(evt,cut[0]),cut[2]):
                    survive=True
                    #print "passed cut",cut
                else:
                    survive=False
                    #return survive
                    break
                    #print "failed ",cut
        
    return survive 

def returnArray(masterArray,variable):
    if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
        val = masterArray[variable][:,0]
    else:
        val = masterArray[variable]
    return val

#create mask for all the events that do or don't pass cuts
def cutOnArray(masterArray,cuts):

    mask=np.full(len(masterArray["evt"]),True)

    for cut in cuts:
        #print "working on cut ",cut
        #Every USE CASE MUST MULTIPLY BY MASK! 
        if type(cut[0][0])==type([0]):
            mask = cutOnArray(masterArray,cut)
        #
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
            #recasting tempmask to boolian/integer array
            tempmask = ops[cut[3]](tempmask,cut[4])
            #print "EQT tempmask",tempmask

            #print " EQT mask before ",mask
            mask *= tempmask.astype(bool)
            continue

        if cut[0][0]=="OR":
            tempmasks=[]
            tempmask=np.full(len(masterArray["evt"]),False)
            for oe in range(1,len(cut)):
                #print cut[oe]
                oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
                #print oneor 
                tempmask = ops["bor"](tempmask,oneor).astype(bool)
            #print "tempmask ",tempmask 
            #print "mask before ",mask
            mask *= tempmask.astype(bool) 
            #print "mask after or  ",mask
            continue
            
        #default case done? FIX ABS!
        else:
            statement = cut[1]
            if statement=="absg": # lessthan OR greater than
                tempmask = ops["<"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                tempmask2= ops[">"](returnArray(masterArray,cut[0]),cut[2])
                tempmask = ops["bor"](tempmask,tempmask2)
            if statement=="absl": # greater than AND less tahn
                #print cut[0],returnArray(masterArray,cut[0])
                tempmask = ops[">"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                #print "tempmask absl 1",tempmask 
                tempmask *= ops["<"](returnArray(masterArray,cut[0]),cut[2])
                #print "tempmask absl 2",tempmask 
            else:
                #print statement 
                #print cut[0],returnArray(masterArray,cut[0])
                tempmask = ops[statement](returnArray(masterArray,cut[0]),cut[2]) 
                #print "tempmask ",tempmask 
            #print "mask before ",mask
            mask *= tempmask.astype(bool) 
        #print "cut ",cut,"  mask ",mask

    #masterArray["mask"] = mask

    return mask

#Function that combines neighboring bins into a list 
def returnBinList(flatbins):
    comboBins=[]
    for sl in range(0,len(flatbins)):
        comboBins.append([])
        for bin in range(0,len(flatbins[sl])-1):
            comboBins[sl].append([flatbins[sl][bin],flatbins[sl][bin+1]])
    return comboBins

def returnPerms(setlist):

    combinations=list(itertools.product(*setlist))

    return combinations


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
    #print "Primary Bins  ",primaryBins
    #print "Combinations for secondary bins  ",combos
    procut = HAA_ProCuts.get(process)
    histList = []
    urnum = 0 
    for bincuts in combos:
        htmp=ROOT.TH1D("htemp"+str(process)+str(urnum),"htemp"+str(process)+str(urnum),len(primaryBins)-1,primaryBins)
        print "working on variable   ",str(variable[0]),"for htemp"+str(process)+str(urnum)
        unrollCut=""

        for extraVarnum in range(1,len(variable)):
            unrollCut += "&&"+str(variable[extraVarnum])+">"+str(bincuts[extraVarnum-1][0])+"&&"+str(variable[extraVarnum])+"<"+str(bincuts[extraVarnum-1][1])
        if procut:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+"&&"+str(HAA_ProCuts[process])+unrollCut+")"
            #print "entries after cut  ",tree.GetEntries(allcuts)
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
        else:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+unrollCut+")"
            #print "entries after cut  ",tree.GetEntries(allcuts)
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
         
        histList.append(htmp)
        htmp.Write(htmp.GetName(),ROOT.TObject.kOverwrite)
        urnum+=1
        #print "htmp entries  ",htmp.GetEntries()
        #print allcuts
        htmp=0
    #next create a list of histograms and try to chain them together bin by bin. 
    masterbins = cat.binning[numvar][0]*len(combos)
    print masterbins
    masterUnroll = ROOT.TH1D(str(process)+"_unrolled",str(process)+"_unrolled",len(masterbins),0,len(masterbins))
    masterBinNum=0 
    for hist in histList:
        for bin in range(0,hist.GetNbinsX()+1):
            masterUnroll.SetBinContent(masterBinNum,hist.GetBinContent(bin))
            #print "bin    ",masterBinNum,"    bin content ",hist.GetBinContent(bin)
            masterUnroll.SetBinError(masterBinNum,hist.GetBinError(bin))
            masterBinNum+=1
           
    
    
    masterUnroll.Write(masterUnroll.GetName(),ROOT.TObject.kOverwrite)
    return

#for fake factor return the bin value from event content ... only for loop thus far 
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    newArr = np.vectorize(pt.FindBin)(numpyArr)
    newArr = np.vectorize(pt.GetBinContent)(newArr)
    #for val in range(len(numpyArr)):
    #    newArr[val] = pt.GetBinContent(pt.FindBin(numpyArr[val]))   
    return newArr


def calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,HAA_processes,processObj,nickname,histodict,weightstring,commonweight,datadrivenPackage,verbose,test,weightHisto,jetWeightMultiplicity):

    newVarVals={}
    for var in HAA_Inc_mmmt.newvariables.keys():
        newVarVals[var]=0.0

    #print "categories   ",allcats.keys() 
    #print "process obj  ",processObj.nickname
    #  
    #for cat in allcats.keys():
    #    cuts = []
    #    for cuttype in allcats[cat].cuts.keys():
    #        for cut in allcats[cat].cuts[cuttype]:
    #            cuts.append(cut)
    #    print "category ",cat 
    #    print "cat cuts",cuts
    #    print "process obj cuts ",processObj.cuts
    #print "all cuts ",cuts

    #to dataframe?
    masterArray = tree.arrays()

    #creating final weight 
    #scalar weight

    for cat in allcats.keys():

        for process in processObj.cuts.keys():
            procut = processObj.cuts[process]
            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)
            if procut!="":
                cuts.append(procut[0])

            #creat array after all the cuts 
            masterArray["finalweight"]=np.full(len(masterArray["evt"]),1.0)
            
            if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
                mask = cutOnArray(masterArray,cuts)
                #np.set_printoptions(threshold=sys.maxsize)
                #print "full mask before weights ",mask
                #np.set_printoptions(threshold=1000)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                #np.set_printoptions(threshold=sys.maxsize)
                #print "full weight after mask ",masterArray["finalweight"]
                #np.set_printoptions(threshold=1000)

            if process=="data_obs":
                mask = cutOnArray(masterArray,cuts)
                #np.set_printoptions(threshold=sys.maxsize)
                #print "full mask before weights ",mask
                #np.set_printoptions(threshold=1000)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!
                

            #print "working on ",process
            if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
                weightDict = processObj.weights
                weightfinal = commonweight
                for scalefactor in weightDict.keys():
                    #if scalefactor == "nevents":
                    #    weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    if scalefactor == "kfactor":
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    elif scalefactor in ["PU"]:
                        masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                    elif scalefactor =="theoryXsec":
                        weightfinal =  weightfinal * float(weightDict[scalefactor])
                    #else:
                    #    weightfinal =  weightfinal * float(weightDict[scalefactor])
            
                #if nickname in ["DYJetsToLLM10to50","DYJetsToLLext1","DYJetsToLLext2"]:
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
                    weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ weightHisto.GetSumOfWeights()
            
                #multiply by scalar weight
                masterArray["finalweight"] *= weightfinal
                masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0

                
            if(process=="FF_1" and datadrivenPackage["bool"]):
                tempmask=np.full(len(masterArray["evt"]),1.0)

                #the actual events that pass the FF_1 criteria 
                tempmask = cutOnArray(masterArray,processObj.cuts["FF_1"])
                tempmask = tempmask.astype(int)
                #np.set_printoptions(threshold=sys.maxsize)
                #print "ff1 tempmask ",tempmask
                   
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                fitmask = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
                fitmask = fitmask.astype(int) 
                ptarr = masterArray["pt_3"]
                #ptarr *= tempmask 
                #ptarr *= fitmask 
                np.set_printoptions(threshold=sys.maxsize)

                #get pt binned value per event... unfortunately need a loop here
                ffweight = ptFun(datadrivenPackage["fakerate1"],ptarr)
                ffweight = ffweight/(1.0000000001 - ffweight)
                print "bincontent ffweight ", ffweight
                #apply fit mask 
                ##print "full mask before weights ",mask
                #np.set_printoptions(threshold=1000)
                np.set_printoptions(threshold=sys.maxsize)
                
                #replace 0s with constant fit value 
                ffweight *= fitmask
                ffweight[np.where(ffweight==0)] =  datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0))  
                ffweight *= tempmask
                print "final ffweight ", ffweight

                #weightfinal = 0.025/(1-0.025) # from previous paper const
                #weightfinal = datadrivenPackage["pseudofit1"].Eval(evt.pt_3)/(1.0000001-datadrivenPackage["pseudofit1"].Eval(evt.pt_3)) #const for now
                masterArray["finalweight"] *= ffweight
                masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0

                for var in newVarVals.keys():
                    #print allcats[cat].newvariables[var][1]
                    #print var+":"+allcats[cat].name+":"+process
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag)) 

                    root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


                for variableHandle in allcats[cat].vars.keys():
                    #move this outside variables

                    #obtaining the right variable... WHAT about changed variables in the event!!?
                    #val = returnValue(evt,allcats[cat].vars[variableHandle][0])
                    val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                    filedict[variableHandle].cd()
                    filedict[variableHandle].cd(allcats[cat].name)
                    #Fix the weight problem!!!
                    #print "filling histogram ",variableHandle," summed entries ",np.sum(val*masterArray["finalweight"])
                    root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])

                continue 

            if(process=="FF_2" and datadrivenPackage["bool"]):
                #the actual events that pass the FF_1 criteria 
                tempmask=np.full(len(masterArray["evt"]),1.0)
                tempmask = cutOnArray(masterArray,processObj.cuts["FF_2"])
                tempmask = tempmask.astype(int)
                print "ff1 tempmask ",tempmask
                   
                   
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                fitmask = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
                fitmask = fitmask.astype(int) 
                ptarr = masterArray["pt_4"]
                #ptarr *= tempmask 

                #get pt binned value per event... unfortunately need a loop here
                ffweight = ptFun(datadrivenPackage["fakerate2"],ptarr)
                ffweight = ffweight/(1.0000000001 - ffweight)
                #apply fit mask 
                #ffweight *= fitmask
                
                #replace 0s with constant fit value 
                ffweight *= fitmask
                ffweight[np.where(ffweight==0)] =  datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))  
                ffweight *= tempmask

                #weightfinal = 0.025/(1-0.025) # from previous paper const
                #weightfinal = datadrivenPackage["pseudofit1"].Eval(evt.pt_3)/(1.0000001-datadrivenPackage["pseudofit1"].Eval(evt.pt_3)) #const for now
                masterArray["finalweight"] *= ffweight
                masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0


                for var in newVarVals.keys():
                    #print allcats[cat].newvariables[var][1]
                    #print var+":"+allcats[cat].name+":"+process
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag)) 

                    root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


                for variableHandle in allcats[cat].vars.keys():
                    #move this outside variables

                    #obtaining the right variable... WHAT about changed variables in the event!!?
                    #val = returnValue(evt,allcats[cat].vars[variableHandle][0])
                    val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                    filedict[variableHandle].cd()
                    filedict[variableHandle].cd(allcats[cat].name)
                    #Fix the weight problem!!!
                    #print "filling histogram ",variableHandle," summed entries ",np.sum(val*masterArray["finalweight"])
                    root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])


                continue 

            if(process=="FF_12" and datadrivenPackage["bool"]):
                #the actual events that pass the FF_1 criteria 
                tempmask=np.full(len(masterArray["evt"]),1.0)
                tempmask = cutOnArray(masterArray,processObj.cuts["FF_12"])
                tempmask = tempmask.astype(int)
                   
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                fitmask = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
                fitmask = fitmask.astype(int) 

                fitmask2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
                fitmask2 = fitmask2.astype(int) 

                ptarr = masterArray["pt_3"]
                ptarr2 = masterArray["pt_4"]
                #ptarr *= tempmask 
                #ptarr2 *= tempmask2 

                #get pt binned value per event... unfortunately need a loop here
                ffweight = ptFun(datadrivenPackage["fakerate1"],ptarr)
                ffweight2 = ptFun(datadrivenPackage["fakerate2"],ptarr2)

                ffweight = ffweight/(1.0000000001 - ffweight)
                ffweight2 = ffweight2/(1.0000000001 - ffweight2)
  
                ffweight = ffweight * ffweight2 
                #apply fit mask 
                fitmask *= fitmask2
                ffweight *= fitmask
                
                #replace 0s with constant fit value 
                ffweight[np.where(ffweight==0)] =  (datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0)))  
                ffweight *= tempmask

                #weightfinal = 0.025/(1-0.025) # from previous paper const
                #weightfinal = datadrivenPackage["pseudofit1"].Eval(evt.pt_3)/(1.0000001-datadrivenPackage["pseudofit1"].Eval(evt.pt_3)) #const for now
                masterArray["finalweight"] *= ffweight
                masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0
                
                for var in newVarVals.keys():
                    #print allcats[cat].newvariables[var][1]
                    #print var+":"+allcats[cat].name+":"+process
                    arguments = allcats[cat].newvariables[var][2]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag)) 

                    root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


                for variableHandle in allcats[cat].vars.keys():
                    #move this outside variables

                    #obtaining the right variable... WHAT about changed variables in the event!!?
                    #val = returnValue(evt,allcats[cat].vars[variableHandle][0])
                    val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                    filedict[variableHandle].cd()
                    filedict[variableHandle].cd(allcats[cat].name)
                    #Fix the weight problem!!!
                    #print "filling histogram ",variableHandle," summed entries ",np.sum(val*masterArray["finalweight"])
                    root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])

                continue 

            if process!="data_obs" and len(tree)!=0:
                eventWeightDict = processObj.eventWeights 
                if eventWeightDict:
                    for scalefactor in eventWeightDict.keys():

                        cutlist = eventWeightDict[scalefactor][0]

                        weightMask=cutOnArray(masterArray,cutlist)
                        weightMask=weightMask.astype(float)

                        #print "update to finalWeight ",masterArray["finalweight"]
                        if type(eventWeightDict[scalefactor][1][0])==float:
                            weightMask *= eventWeightDict[scalefactor][1][0] 

                        if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                            arguments = eventWeightDict[scalefactor][1][1]
                            tempvals=[]

                            for ag in arguments:
                               tempvals.append(returnArray(masterArray,ag)) 
                            weightMask*= eventWeightDict[scalefactor][1][0](*tempvals) 

                        #preseving the events that don't need new weight
                        weightMask[np.where(weightMask ==0.0)]=1.0
                        #print "weight mask ",weightMask
                        #updating the finalweight
                        masterArray["finalweight"] *= weightMask
                        masterArray["finalweight"][masterArray["finalweight"] < 0 ] = 0
            #np.set_printoptions(threshold=sys.maxsize)
            ##print "full mask before weights ",mask
            #print " weights ",masterArray["finalweight"]
            #np.set_printoptions(threshold=1000)

            for var in newVarVals.keys():
                #print allcats[cat].newvariables[var][1]
                #print var+":"+allcats[cat].name+":"+process
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]

                for ag in arguments:
                   tempvals.append(returnArray(masterArray,ag)) 

                root_numpy.fill_hist(newhistodict[var+":"+allcats[cat].name+":"+process],functs[allcats[cat].newvariables[var][0]](*tempvals),masterArray["finalweight"])


            for variableHandle in allcats[cat].vars.keys():
                #move this outside variables

                #obtaining the right variable... WHAT about changed variables in the event!!?
                #val = returnValue(evt,allcats[cat].vars[variableHandle][0])
                val = returnArray(masterArray,allcats[cat].vars[variableHandle][0])
                filedict[variableHandle].cd()
                filedict[variableHandle].cd(allcats[cat].name)
                #Fix the weight problem!!!
                #print "filling histogram ",variableHandle," summed entries ",np.sum(val*masterArray["finalweight"])
                root_numpy.fill_hist(histodict[variableHandle+":"+allcats[cat].name+":"+process],val,masterArray["finalweight"])
    
    



    return

if __name__ == "__main__":
    
    begin_time = datetime.datetime.now()




    #Structure for mapping between root files and processes  
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights
    from utils.Weights import HAAWeights

    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
    #sys.path.append('../SFs/')
    import utils.SFs.ScaleFactor as SF

    import io
    import os
    from shutil import copyfile

    import argparse

    #consider adding yaml file here ... :)
    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/",  help="Input files")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
    parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
    parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
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
 
    #with open("MCsamples_2016_v6_yaml.csv")  as csvfile:
    #    reader = csv.reader(csvfile, delimiter=',')
    #    for row in reader:
    #        sampleDict[row[0]] = [row[1],row[2],row[3],row[4],row[5],row[6]]

    for line in open(args.csvfile,'r').readlines() :
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
        #print category
        #print categories[category]['name']
        tempcat = Category() 
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat

    
    #loading standard processes
    HAA_processes={}
    #place fake rate here 
    if not (args.datadrivenZH):
        for sample in sampleDict.keys():
            #print processes[process]
            temppro = Process() 
            temppro.nickname=sample
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={sampleDict[sample][0]:""}
            if "ggTo2mu2tau" in sample:
                temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)}
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
            #print processes[process]
            temppro = Process() 
            temppro.nickname=sample
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue"}
            temppro.cuts={sampleDict[sample][0]:"","fake1_"+sampleDict[sample][0]:[["gen_match_3","==",0]],"fake2_"+sampleDict[sample][0]:[["gen_match_4","==",0]]}
            if "ggTo2mu2tau" in sample:
                temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)}
            if "W" in sample and "Jets" in sample and "EWK" not in sample: 
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","kfactor":1.221}
                temppro.cuts={"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]],"fake1_W":[["gen_match_3","==",0]],"fake2_W":[["gen_match_4","==",0]]}
            if "DY" in sample and "Jets" in sample: 
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","kfactor":1.1637}
            if "TT" in sample and not "TTTT" in sample and not "TTHH" in sample: 
                temppro.file=sample+"_2016.root"
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue"}
                temppro.cuts={"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]],"fake1_TT":[["gen_match_3","==",0]],"fake2_TT":[["gen_match_4","==",0]]}
            HAA_processes[temppro.nickname]=temppro

    #loading special processes ... fake factor and data
    for process in processes_special:
        temppro = Process() 
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro
        
    if args.maxprint:
        print "categories "
        for cat in allcats:
            print cat.name
            print cat.cuts
        

        print "processes "
        for pro in HAA_processes.keys():
            print HAA_processes[pro].nickname
            print HAA_processes[pro].cuts
            print HAA_processes[pro].weights
            print HAA_processes[pro].file
    
             
    
        
    #This is where the plotting takes place!
    #categories.append(HAA_Inc_mmmt) 
    #categories = allcats

    #Gather the Analysis Files
    #dir = "/eos/user/s/shigginb/HAA_ntuples/March_2020/"
    #dir = "/eos/home-s/shigginb/HAA_ntuples/SignalOnly/"
    #dir = "/eos/home-s/shigginb/HAA_ntuples/June2020/" #eos problems
    #dir = "/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/"
    dir = args.dir
    #filelist = HAA
    filelist = {}

    #Complex set of event weights
    # tau id scale factor from object
    cmsswbase="$CMSSW_BASE/src"
    #gInterpreter.ProcessLine(".include .")
    #for baseName in ['../SVFit/MeasuredTauLepton','../SVFit/svFitAuxFunctions','../SVFit/FastMTT', cmsswbase+'/HTT-utilities/RecoilCorrections/src/MEtSys', cmsswbase+'/HTT-utilities/RecoilCorrections/src/RecoilCorrector']:
    #    if os.path.isfile("{0:s}_cc.so".format(baseName)):
    #        gInterpreter.ProcessLine(".L {0:s}_cc.so".format(baseName))
    #    else:
    #        gInterpreter.ProcessLine(".L {0:s}.cc++".format(baseName))   
    tauIDSF = TauIDSFTool('2016Legacy','DeepTau2017v2p1VSjet','Medium').getSFvsPT
    #recoilCorrector  = ROOT.RecoilCorrector("HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root").CorrectByMeanResolution
    # muon scale factor as function of pt
    sf_MuonId = SF.SFs()
    sf_MuonId.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/LeptonEffs/Muon/Muon_Run2016_IdIso_0p2.root")
    EventWeights={
        #"name":[[if statements],[weight to apply]]
        "3_mt_lt0p4":[[["decayMode_3","==",0],["eta_3","<",0.4]],[0.80]],
        "3_mj_lt0p4":[[["decayMode_3","==",0],["eta_3","<",0.4]],[1.21]],
        "3_mt_0p4to0p8":[[["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[0.81]],
        "3_mj_0p4to0p8":[[["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[1.11]],
        "3_mt_0p8to1p2":[[["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[0.79]],
        "3_mj_0p8to1p2":[[["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[1.2]],
        "3_mt_1p2to1p7":[[["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[0.68]],
        "3_mj_1p2to1p7":[[["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[1.16]],
        "3_mt_1p7to2p3":[[["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[0.68]],
        "3_mj_1p7to2p3":[[["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[2.25]],

        "4_mt_lt0p4":[[["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
        "4_mj_lt0p4":[[["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
        "4_mt_0p4to0p8":[[["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
        "4_mj_0p4to0p8":[[["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
        "4_mt_0p8to1p2":[[["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
        "4_mj_0p8to1p2":[[["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
        "4_mt_1p2to1p7":[[["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
        "4_mj_1p2to1p7":[[["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
        "4_mt_1p7to2.3":[[["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
        "4_mj_1p7to2.3":[[["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

        #bound method ... last input list is func parameters
        "3_tauSF":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]]

        #recoil correction
        #"recoilcorr2":[[["cat","==",6],["njets","==",2]],[recoilCorrector,["met_x","met_y","gen_match_3"]]],
    }
    for objkey in HAA_processes.keys():
        if not objkey=="data":
            HAA_processes[objkey].eventWeights = EventWeights
            #HAA_processes[objkey].cuts = {sampleDict[HAA_processes[objkey].nickname][0]:""}
            #HAA_processes[objkey].cuts["prompt"] = [["gen_match_3","!=",0],["gen_match_4","!=",0]] # doing this in other script 
            print HAA_processes[objkey].nickname+"_2016"
    print "added SF weights"

    #for testing only... one cat a15
    #HAA_processes=HAA_processes_test
    if args.datameasure:
        filelist["data"]=HAA_processes["data"].file
        for category in allcats[:]:
            if category.name=="mmmt_inclusive":
                allcats.remove(category)
    if not args.datameasure:
        for proObj in HAA_processes.keys():
           filelist[proObj]=HAA_processes[proObj].file

    print filelist

    if (args.datameasureZH):
        try:
            os.mkdir("FFhistos_"+str(args.ffout))
        except:
            print "directory exists"
        for proObj in HAA_processes.keys():
            #if proObj!="data":  #isn't this right ... we want Prompt MC contribution 
            if proObj!="data" and proObj!="FF":  #isn't this right ... we want Prompt MC contribution 
                HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
                #HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","==",0]]
                #HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","==",0]]
                HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
        for category in allcats.keys():
            if allcats[category].name=="mmmt_inclusive":
                del allcats[category]
    elif (args.datadrivenZH):
        for proObj in HAA_processes.keys():
            if proObj!="data" and proObj!="FF":  #isn't this right ... we want Prompt MC contribution 
                HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
                #HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","==",0]]
                #HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","==",0]]
                HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]

        for category in allcats.keys():
            if allcats[category].name!="mmmt_inclusive":
                print "removing category ",allcats[category].name
                del allcats[category]
    else:
        for category in allcats.keys():
            if allcats[category].name!="mmmt_inclusive":
                print "removing category ",allcats[category].name
                del allcats[category]

    treename = "Events"

    print "initializing histograms"
    
    try:
        os.mkdir("out"+str(args.outname))
    except:
        print "directory exists"
    filedict = {}
    histodict = {}
    newhistodict = {}
    #cat=HAA_Inc_mmmt

    numvar=0
    for variableHandle in allcats[allcats.keys()[0]].vars.keys():
    
        #variable=fullvariable.split(":") # for the unrolled cases

        #filedict[variable[0]]=ROOT.TFile.Open(cat.name+"_"+str(variable[0])+".root","RECREATE")
        filedict[variableHandle]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variableHandle)+".root","RECREATE")
        #print "on variableHandle ",variableHandle
        filedict[variableHandle].cd()
        for cat in allcats.keys():
            #print "Working on category.",allcats[cat].name
            #print "With cuts ",allcats[cat].cuts
            filedict[variableHandle].mkdir(allcats[cat].name)
            filedict[variableHandle].cd(allcats[cat].name)

            for nickname in filelist.keys():
                #print "for sample ",nickname

                processObj = HAA_processes[nickname]

                for process in processObj.cuts.keys():


                    if filedict[variableHandle].Get(allcats[cat].name).GetListOfKeys().Contains(str(process)):
                        continue
                    else:
                        #bins = allcats[cat].binning[numvar]
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

        #filedict[variable]=ROOT.TFile.Open(cat.name+"_"+str(variable)+".root","RECREATE")
        filedict[variable]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variable)+".root","RECREATE")
        #print "on variable ",variable
        filedict[variable].cd()
        #for cat in allcats:
        for cat in allcats.keys():
            filedict[variable].mkdir(allcats[cat].name)
            filedict[variable].cd(allcats[cat].name)

            for nickname in filelist.keys():

                processObj = HAA_processes[nickname]

                for process in processObj.cuts.keys():
                    #print "process ",process,"  cuts  ",processObj.cuts[process]

                    if filedict[variable].Get(allcats[cat].name).GetListOfKeys().Contains(str(process)):
                        continue
                    else:
                        #tmpbin = np.asarray(allcats[cat].newvariablesbins[numvar])
                        bins = allcats[allcats.keys()[0]].newvariables[variable][1]
                        tmpbin = np.asarray(bins)
                        #print variable+":"+allcats[cat].name+":"+process
                        if  allcats[allcats.keys()[0]].newvariables[variable][-1]=="multiob":
                            for var,ivar in enumerate(temp):
                                newhistodict[variable+"_"+str(ivar)+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin[ivar])-1,tmpbin[ivar])
                                newhistodict[variable+"_"+str(ivar)+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                        else:
                            newhistodict[variable+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                            newhistodict[variable+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                        #newhistodict[variable+":"+allcats[cat].name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                        #newhistodict[variable+":"+allcats[cat].name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
        numvar=numvar+1

    #gather extra global variables or weights
    
    weight = CommonWeights["lumi"][0]
    weightstring = CommonWeights["string"]
   
    #weight = weight+"*"+CommonWeights["mcweight"][0]

    #precise weights
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

    #Calculate the scale factors and fill the histograms 
    #print sampleDict.keys()
    datadrivenPackage={}
    datadrivenPackage["bool"]=False

    if args.datadriven:
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadriven
        ff_file_3 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/jpt_1_ff.root")
        ff_file_4 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/jpt_2_ff.root")
        #sstight_3 = ff_file_3.Get("mmmt_inclusive_FF_tight_SS/data_obs")
        #sstight_4 = ff_file_4.Get("mmmt_inclusive_FF_tight_SS/data_obs")
        #ssloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_SS/data_obs")
        #ssloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_SS/data_obs")
        #osloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_OS/data_obs")
        #osloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_OS/data_obs")
        ss_1_tight = ff_file_3.Get("mmmt_inclusive_FF_SS_1_tight/data_obs")
        ss_1_loose = ff_file_3.Get("mmmt_inclusive_FF_SS_1_loose/data_obs")
        ss_2_tight = ff_file_4.Get("mmmt_inclusive_FF_SS_2_tight/data_obs")
        ss_2_loose = ff_file_4.Get("mmmt_inclusive_FF_SS_2_loose/data_obs")
        
        #f_1 = sstight_3.Clone()
        f_1= ss_1_tight.Clone()
        #f_1.Divide(ssloose_3)
        f_1.Divide(ss_1_loose)
        f_1.SetName("FakeRateMuLeg")
        f_1.GetXaxis().SetTitle("p_T #mu")
        f_1.GetYaxis().SetTitle("Fake Rate for #mu")
        f_2 = ss_2_tight.Clone()
        f_2.Divide(ss_2_loose)
        f_2.SetName("FakeRateTauLeg")
        f_2.GetXaxis().SetTitle("p_T #tau")
        f_2.GetYaxis().SetTitle("Fake Rate for #tau")

        #tf_1 = ROOT.TF1("tf_1","[0]*expo[1]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
        #tf_2 = ROOT.TF1("tf_2","[0]*expo[1]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 
        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
        tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 

        fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/fakemeasure.root","RECREATE")
        fakemeasurefile.cd()                       
        f_1.Fit("tf_1")
        f_2.Fit("tf_2")
        f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
        f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
        datadrivenPackage["fakerate1"]=f_1
        datadrivenPackage["fakerate2"]=f_2
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fitrate2"]=tf_2
        fakemeasurefile.Close()
        #applyf = ROOT.TFile.Open("FFhistos/ffmeas.root","READ")
        #applyf.cd()
        #pg1 = applyf.Get("mupoints")
        #pg2 = applyf.Get("taupoints")
        #pf1 = ROOT.TF1("pf1","[0]*exp([1]*x)+[2]",0,100) 
        #pf2 = ROOT.TF1("pf2","[0]*exp([1]*x)+[2]",0,100) 
        ##pg1.Fit("pf1")
        ##pg2.Fit("pf2")
        #datadrivenPackage["pseudofit1"]=pf1
        #datadrivenPackage["pseudofit2"]=pf2

    if args.datadrivenZH or args.makeFakeHistos:
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadrivenZH
        ff_file_3 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_3_ff.root","read")
        ff_file_4 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_4_ff.root","read")
        #sstight_3 = ff_file_3.Get("mmmt_inclusive_FF_tight_SS/data_obs")
        #sstight_4 = ff_file_4.Get("mmmt_inclusive_FF_tight_SS/data_obs")
        #ssloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_SS/data_obs")
        #ssloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_SS/data_obs")
        #osloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_OS/data_obs")
        #osloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_OS/data_obs")
        ss_1_tight = ff_file_3.Get("mmmt_inclusive_FF_SS_1_tight/data_obs")
        ss_1_loose = ff_file_3.Get("mmmt_inclusive_FF_SS_1_loose/data_obs")
        ss_2_tight = ff_file_4.Get("mmmt_inclusive_FF_SS_2_tight/data_obs")
        ss_2_loose = ff_file_4.Get("mmmt_inclusive_FF_SS_2_loose/data_obs")

        ss_1_tight_prompt = ff_file_3.Get("mmmt_inclusive_FF_SS_1_tight/prompt1")
        ss_1_loose_prompt = ff_file_3.Get("mmmt_inclusive_FF_SS_1_loose/prompt1")
        ss_2_tight_prompt = ff_file_4.Get("mmmt_inclusive_FF_SS_2_tight/prompt2")
        ss_2_loose_prompt = ff_file_4.Get("mmmt_inclusive_FF_SS_2_loose/prompt2")

        #subtracting prompt MC
        ss_1_tight.Add(ss_1_tight_prompt,-1)
        ss_2_tight.Add(ss_2_tight_prompt,-1)
        ss_1_loose.Add(ss_1_loose_prompt,-1)
        ss_2_loose.Add(ss_2_loose_prompt,-1)
        
        #f_1 = sstight_3.Clone()
        f_1= ss_1_tight.Clone()
        #f_1.Divide(ssloose_3)
        f_1.Divide(ss_1_loose)
        f_2 = ss_2_tight.Clone()
        f_2.Divide(ss_2_loose)

        f_1.SetName("FakeRateMuLeg")
        f_1.GetXaxis().SetTitle("p_T #mu")
        f_1.GetYaxis().SetTitle("Fake Rate for #mu")
        f_2.SetName("FakeRateTauLeg")
        f_2.GetXaxis().SetTitle("p_T #tau")
        f_2.GetYaxis().SetTitle("Fake Rate for #tau")

        #tf_1 = ROOT.TF1("tf_1","[0]*expo[1]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
        #tf_2 = ROOT.TF1("tf_2","[0]*expo[1]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 
        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
        tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 

        fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/fakemeasure.root","RECREATE")
        fakemeasurefile.cd()
        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptFit()
        f_1.Draw()
        f_1.Fit("tf_1")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/fakerate1.png")
        f_2.Draw()
        f_2.Fit("tf_2")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/fakerate2.png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
        f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)
        datadrivenPackage["fakerate1"]=f_1
        datadrivenPackage["fakerate2"]=f_2
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fitrate2"]=tf_2
        fakemeasurefile.Close()
        #applyf = ROOT.TFile.Open("FFhistos/ffmeas.root","READ")
        #applyf.cd()
        #pg1 = applyf.Get("mupoints")
        #pg2 = applyf.Get("taupoints")
        #pf1 = ROOT.TF1("pf1","[0]*exp([1]*x)+[2]",0,100) 
        #pf2 = ROOT.TF1("pf2","[0]*exp([1]*x)+[2]",0,100) 
        ##pg1.Fit("pf1")
        ##pg2.Fit("pf2")
        #datadrivenPackage["pseudofit1"]=pf1
        #datadrivenPackage["pseudofit2"]=pf2




    #Where the magic starts! 
    for nickname in filelist.keys():
        print "working on file ",filelist[nickname]
        frootin = ROOT.TFile.Open(dir+filelist[nickname],"read")
        fin = uproot.open(dir+filelist[nickname])
        #dataframe
        tree = fin[treename]
        weightHisto = frootin.Get("hWeights")
        #newtree=tree.CloneTree()
        print "tree entries ",len(tree)

        processObj = HAA_processes[nickname]

        calculateHistos(functs, tree, allcats[allcats.keys()[0]], allcats, HAA_processes, processObj, 
                        nickname, histodict, weightstring, weight, datadrivenPackage, args.verbose, 
                        args.test, weightHisto, jetWeightMultiplicity)

    print "writing the histograms to the file"
    for varcatPro in histodict.keys():
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
            
        histodict[varcatPro].Write(histodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
    print "writing the new variable histograms to the file"
    for varcatPro in newhistodict.keys():
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        newhistodict[varcatPro].Write(newhistodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)

    if args.datameasureZH:
        #copyfile("out"+str(args.outname)+"/pt*.root", "FFhistos_"+str(args.ffout)+"/.")
        for file in glob.glob("out"+str(args.outname)+"/pt*.root"):
            print "copying over ",file
            copyfile(file, "FFhistos_"+str(args.ffout)+"/"+file.split("/")[1])
    del sf_MuonId
    print "computation time"        
    print(datetime.datetime.now() - begin_time)

