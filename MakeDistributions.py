# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDistributions.py

* Purpose : For visualizing datasets. Will create ultra skimmed root files or histograms for visualization and for limit setting later

* Usage:
python MakeDistributions.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /path/to/input/root/files/ -p processes_special_mmtt.yaml -o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt -s

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

year = 2018 #set by args later; this is just so it's a global variable.
#directory where all the ntuple files are stored.
#direc = "/eos/home-s/shigginb/HAA_ntuples/2016_v7/"
#direc = "/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/"


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
    #print "working on branch ",variable
    #if len(masterArray[variable].shape) == 1:   # flat array important for jagged arrays of input data
    #    val = masterArray[variable]
    if "[" in variable:
        basevar = variable.split("[")[0]
        index = int(variable.split("[")[1].split("]")[0])
        val = masterArray[basevar][:,index]
    else:
        val = masterArray[variable]
    return val

#create mask for all the events that do or don't pass cuts
def cutOnArray(masterArray,cuts):

    mask=np.full(len(masterArray["evt"]),True)

    for cut in cuts:
        #Every USE CASE MUST MULTIPLY BY MASK!
        #print cut
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


#for fake factor return the bin value from event content
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    newArr = np.vectorize(pt.FindBin)(numpyArr)
    newArr = np.vectorize(pt.GetBinContent)(newArr)

    return newArr

def getEventWeightDicitonary():
    #import copyreg, copy, pickle # for picking the bound methods due to the multiplrocess tool

    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
    tauIDSF = TauIDSFTool('2016Legacy','DeepTau2017v2p1VSjet','Medium').getSFvsPT

    #def _pickle_method(method):
    #    func_name = method.im_func.__name__
    #    obj = method.im_self
    #    cls = method.im_class
    #return _unpickle_method, (func_name, obj, cls)
#
    #def _unpickle_method(func_name, obj, cls):
    #    for cls in cls.mro():
    #    try:
    #    func = cls.__dict__[func_name]
    #    except KeyError:
    #    pass
    #    else:
    #    break
    #return func.__get__(obj, cls)
    #copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


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

    return EventWeights

def initialize(channel, direc):
    import os
    from copy import copy
    import yaml

    #Structure for mapping between root files and processes
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
   # from utils.Weights import CommonWeights
    from utils.Weights import HAAWeights


    #sys.path.append('../SFs/')
    #import utils.SFs.ScaleFactor as SF

    import io
    import os
    from shutil import copyfile


    #gather functions for computing variables in the event loop
    from utils.functions import functs
    from ROOT import gInterpreter

    #Structure for plotting variables
    from utils.Parametrization import Category#, DylanCategory
    from utils.Parametrization import Process

    #Gather the analysis datasets and info
    sampleDict = {}
    #csvfile = "MCsamples_2016_v7_ZZ.csv"
    #csvfile = "MCsamples_2016_v7.csv"
    csvfile = "bpgMCsamples_%d_v7.csv"%(year)
    categories = "cat_"+channel+"_%d.yaml"%(year)
    processes = "processes_special_"+channel+".yaml"
    #categories = "cat_mmet_2016.yaml"
    #processes = "processes_special_mmet.yaml"



    weightHistoDict = {}
    for nickname in filelist.keys():
        try:
            #frooin = ROOT.TFile.Open(direc+filelist[nickname],"read")
            #weightHistoDict[nickname]=frooin.Get("hWeights")
            with ROOT.TFile.Open(direc+filelist[nickname],"read") as frooin:
               weightHistoDict[nickname]=frooin.Get("hWeights")
        except:
            print nickname," hWeights prob doesn't exist?"
    for line in open(csvfile,'r').readlines() :
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
    with io.open(categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print(category)
        #print(categories[category]['name'])
        if categories[category]['name']==channel+"_inclusive":
        #if categories[category]['name']=="mmet_inclusive":
            tempcat = Category()
            tempcat.name=categories[category]['name']
            tempcat.cuts=categories[category]['cuts']
            tempcat.newvariables=categories[category]['newvariables']
            tempcat.vars=categories[category]['vars']
            #tempcat.systematics=categories[category]['systematics']
            allcats[tempcat.name]=tempcat

    #allcats = {name: DylanCategory(cat) for name, cat in categories.items()}


    #for cat in allcats.keys():
        #print(allcats[cat].name)
        #print(allcats[cat].cuts)
        #print(allcats[cat].systematics)

    #loading standard processes
    HAA_processes={}

    #MC for background

    #HAA_process = {name: Process(direc, name, info) for name, info in sampleDict.items()}
    weightHistoDict = {}


    for sample in sampleDict.keys():
        temppro = Process()
        temppro.nickname=sample
        temppro.file=direc+sample+"_%d.root"%(year)

        #with ROOT.TFile.Open(temppro.file,"read") as frooin:
        frooin = ROOT.TFile.Open(temppro.file,"read")
        #frooin.ls()
        try:
            htemp=frooin.Get("hWeights")
            print "entries in hWeights ",htemp.GetEntries()
            weightHistoDict[sample]=copy(htemp)
        except:
            print "SKIPPING FILE NOT WORKING ",sample
            continue
        #frooin.Close()
        extension = "_ext1" #extension added to file name for extended samples
        smp = sample.strip(extension)
        kf = 1.0 #kfactor only for the special processes with many jets.
        if smp in ["W%dJetsToLNu"%(nj) for nj in range(1, 5)]:
            kf = 1.221
        elif smp in ["DY%dJetsToLL"%(nj) for nj in range(1, 5)]:
            kf = 1.1637
        temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"kfactor":kf,"PU":"weightPUtrue"}
        temppro.cuts={sampleDict[sample][0]:""}
        #if "ggTo2mu2tau" in sample:
        if "HToAATo4Tau" in sample:
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} # worked before
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(31.05*0.0001)} #.01% br; no multiplication for visibility.
            temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(31.05*0.75)} #75% br; corresponds to 2016 paper!
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001* 5.0)} # SM Higgs xsec x BR Haa  for signal extraction data MC control plots(AN 17029)
        HAA_processes[temppro.nickname]=temppro


    jetWeightMultiplicity = {}
    #print("WJets file: {}".format(HAA_processes["WJetsToLNu"].file))
    try:
        #get sum of weights for special processes (sum reg and _ext).
        DYJetsFile = ROOT.TFile.Open(HAA_processes["DYJetsToLL"].file,"read")
        jetWeightMultiplicity["DYJetsToLL"]=DYJetsFile.Get("hWeights").GetSumOfWeights()
        if year == 2017:
            DYJetsFilex = ROOT.TFile.Open(HAA_processes["DYJetsToLL_ext1"].file,"read")
            jetWeightMultiplicity["DYJetsToLL"] += DYJetsFilex.Get("hWeights").GetSumOfWeights()
       # DY1JetsFile = ROOT.TFile.Open(HAA_processes["DY1JetsToLL"].file,"read")
       # jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
       # DY2JetsFile = ROOT.TFile.Open(HAA_processes["DY2JetsToLL"].file,"read")
       # jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
       # DY3JetsFile = ROOT.TFile.Open(HAA_processes["DY3JetsToLL"].file,"read")
       # jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
       # DY4JetsFile = ROOT.TFile.Open(HAA_processes["DY4JetsToLL"].file,"read")
       # jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()
    ####bpg added, not sure if needed??#####
       # DYJetsM10to50File = ROOT.TFile.Open(HAA_processes["DYJetsToLLM10to50"].file,"read")
       # jetWeightMultiplicity["DYJetsToLLM10to50"]=DYJetsM10to50File.Get("hWeights").GetSumOfWeights()
    #
        WJetsFile = ROOT.TFile.Open(HAA_processes["WJetsToLNu"].file,"read")
      #  print("opened WJets file.")
        jetWeightMultiplicity["WJetsToLNu"]=WJetsFile.Get("hWeights").GetSumOfWeights()
        if year == 2017:
            WJetsFilex = ROOT.TFile.Open(HAA_processes["WJetsToLNu_ext1"].file,"read")
            jetWeightMultiplicity["WJetsToLNu"] += WJetsFilex.Get("hWeights").GetSumOfWeights()
      #  W1JetsFile = ROOT.TFile.Open(HAA_processes["W1JetsToLNu"].file,"read")
      #  jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()
      #  W2JetsFile = ROOT.TFile.Open(HAA_processes["W2JetsToLNu"].file,"read")
      #  jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
      #  W3JetsFile = ROOT.TFile.Open(HAA_processes["W3JetsToLNu"].file,"read")
      #  jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
      #  W4JetsFile = ROOT.TFile.Open(HAA_processes["W4JetsToLNu"].file,"read")
      #  jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()
    except:
        print "check jetWeightMultiplicity ... the files may not be loaded"


    Bkg = ["DY","W","TT","ST","EWK","Top"]
    irBkg = ["ZZ","ZHToTauTau","vbf","WHTT","Signal","ggZH"]
    TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW4F","HZJ"] #doesn't exist for 2018 rn.
    rareBkg = ["Other","rare","WZ","QCD"]
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

    allSkims = {}
    finalSkims ={}


    #EventWeights = getEventWeightDicitonary()

    return allcats, HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity





def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def fstar(args):
    print("arggggg matey!")
    return f(*args)

def slimskimstar(args):
    #return slimskim(*args)
    return slimskimoutput(*args)

def f(process, categories):
    info('function f')
    print('hello', process.nickname)
    for cat in categories.keys():
        print("category ",cat)
    return 1

def makeCutsOnTreeArray(process, masterArray,allcats,weightHistoDict,systematic):

    from utils.functions import functs
   # from utils.Weights import CommonWeights
    from ROOT import gInterpreter
    #2018 lumi
    commonweight = 59740 #CommonWeights["lumi"][0]
    if year == 2017:
        commonweight = 41800 #2017 lumi
    elif year == 2016:
        commonweight = 35900 #2016 lumi
    skimArrayPerCat = {}
    print "working on process obj ",process.nickname
    for cat in allcats.keys():
        #print(process)
        #print(process.cuts)
        #print(masterArray.keys())
        plottedVars = []

        #if process.nickname=="data_obs":
        if "data" in process.nickname:
            print("data obs!")
            newVarVals = {}
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
                    plottedVars.append(variableHandle)

        #    print("newvariables now")

            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


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

            skipEvents = np.where(mask==0)[0]
            skimArray={}
            print("before skim", len(masterArray["finalweight"]))
            for key in masterArray.keys():
                try:
                    skimArray[key] = masterArray[key][mask]
                except:
                    #print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                    #print "skipping branch ",key
                    continue
            print("after skim", len(skimArray["mll"]), process.file)
            if len(skimArray["mll"])==0:
                #continue
                print("Warning: length of skim array for data is 0.")

            for key in skimArray.keys():
                if key not in plottedVars and key != "finalweight":
                    del skimArray[key]
            print "working on category ",cat
            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray




        if(process.nickname=="FF" and datadrivenPackage["bool"]):
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
                    plottedVars.append(variableHandle)


            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)

            tempmask=np.full(len(masterArray["evt"]),1.0)

            #events that pass the FF_1 criteria
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
#            print "summed final weight ",np.sum(finalWeight)

            keepEvents = ~np.where(finalWeight==0.0)[0]

            skimArray={}
            for key in masterArray.keys():
                skimArray[key] = masterArray[key][keepEvents]

            #print("after skim", len(skimArray["mll"]), process.file)
            if len(skimArray["mll"])==0:
                continue

            for key in skimArray.keys():
                if key not in plottedVars and key != "finalweight":
                    del skimArray[key]
            #return skimArray
            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray




        print("process nickname: {}".format(process.nickname))
        if process.nickname not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
            EventWeights = getEventWeightDicitonary()

            newVarVals={}
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
           #     print("variableHandle: {}".format(variableHandle))
                variable = allcats[cat].vars[variableHandle][0]
           #     print("variable: {}".format(variable))
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
           #         print("will now add {} to plottedVars".format(variableHandle))
                    plottedVars.append(variableHandle)
           #         print("plottedVars: {}".format(plottedVars))


            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
                plottedVars.append(var)

            mask = cutOnArray(masterArray,cuts)
            masterArray["mask"]=mask
            masterArray["finalweight"] *= mask.astype(int)
            weightfinal = 1.0   #don't weight the data!!

            skipEvents = np.where(mask==0)[0]
            skimArray={}
            if len(masterArray["mll"])==0:
                continue

            weightDict = process.weights
            weightfinal = commonweight
            nickname = process.nickname
            for scalefactor in weightDict.keys():
                if scalefactor == "kfactor":
                    #kfactor will be dealt with below.
                    continue
                    #weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                elif scalefactor in ["PU"]:
                    #masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                        #need to multiply by Generator_weight, no?????
                    masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor])*returnArray(masterArray,"Generator_weight"))
                    #pass
                elif scalefactor =="theoryXsec":
                    weightfinal =  weightfinal * float(weightDict[scalefactor])

            #sum of weights (will include this file + its _ext (if it exists) or its non-_ext (if this is an _ext))
            sow = weightHistoDict[nickname].GetSumOfWeights()

            #print("sow, weightfinal before _ext business: {}, {}".format(sow, weightfinal))
            #extension added to the end of the filename for extended files
            extensions = ["_ext1", "_ext2"]

            #name without the extension at the end
            bare_name = nickname
            for ex,extension in enumerate(extensions):
                #name of the extended version of this process (if it exists)
                extname = nickname + extension
                #name of the non-extended version of this process (if it exists)
                nonextname = nickname.strip(extension)
                #name of the file with the alternate extension (eg ext2 instead of ext1)
                xtname = nonextname
                if ex == 0:
                    xtname += extensions[1]
                else:
                    xtname += extensions[0]

                #add the sow from this extended version
                if extname in weightHistoDict:
                    sow += weightHistoDict[extname].GetSumOfWeights()

                
                #add sow from the brother extended file 
                #but need to make sure we don't add weights from this file twice.
                if nonextname != nickname and xtname in weightHistoDict:
            #        print("using xtname block! nickname={}, xtname={}, sow before = {}".format(nickname, xtname, sow))
                    sow += weightHistoDict[xtname].GetSumOfWeights()
            #        print("sow after: {}".format(sow))
                

                #add sow from non-extended file.
                if nonextname != nickname and nonextname in weightHistoDict:
                    sow += weightHistoDict[nonextname].GetSumOfWeights() 
                    bare_name = nonextname

            #print("sow, weightfinal before special block: {}, {}".format(sow, weightfinal))
            #weightfinal array needed, because weights may be different for some events.
            wf_arr = np.full(len(masterArray["finalweight"]), 1.0)
            #special processes require special weighting.
            if bare_name in ["DY%dJetsToLL"%(nj) for nj in range(1, 5)] and "DYJetsToLL" in weightHistoDict:
            #    print("block 0!")
                #norm1 is for the inclusive process DYJetsToLL. (Should already have its _ext added into its jetWeightMultiplicity.)
                norm1 = jetWeightMultiplicity["DYJetsToLL"]/HAA_processes["DYJetsToLL"].weights["xsec"]
                #norm2 is specific to this process (sum of weights calculated just above this loop).
                norm2 = sow / HAA_processes[nickname].weights["xsec"]
         #       print("norm1: {}, norm2: {}".format(norm1, norm2))
                #norm2 also needs to be divided by the kfactor.
                norm2 *= 1.0 / float(weightDict["kfactor"])
         #       print("norm2 adjusted: {}".format(norm2))
                weightfinal = weightfinal * 1/(norm1+norm2) 
         #       print("weightfinal after block 0: {}".format(weightfinal))
                wf_arr *= weightfinal
            elif bare_name in ["W%dJetsToLNu"%(nj) for nj in range(1, 5)] and "WJetsToLNu" in weightHistoDict:
          #      print("block 1!")
                norm1 = jetWeightMultiplicity["WJetsToLNu"] / HAA_processes["WJetsToLNu"].weights["xsec"]
                norm2 = sow / HAA_processes[nickname].weights["xsec"]
         #       print("norm1: {}, norm2: {}".format(norm1, norm2))
                norm2 *= 1.0 / float(weightDict["kfactor"])
         #       print("norm2 adjusted: {}".format(norm2))
                weightfinal = weightfinal * 1/(norm1+norm2)
         #       print("weightfinal after block 1: {}".format(weightfinal))
                wf_arr *= weightfinal
            elif bare_name in ["DYJetsToLL", "WJetsToLNu"] and "DY1JetsToLL" in weightHistoDict:
          #      print("block 2!")
                #any events with 0 LHE_Njets should use normal weighting.
                try:
                    mask0j = masterArray["LHE_Njets"]==0
                except KeyError:
         #           print("Error: Branch LHE_Njets does not exist in file for process {}.".format(nickname))
                    raise KeyError #sys.exit()
                weight0j = weightfinal * HAA_processes[nickname].weights["xsec"] / sow
         #       print("weight0j: {}".format(weight0j))
                wf_arr = weight0j * mask0j.astype(int)
                #events with > 0 LHE_Njets should use same weighting as the corresponding exclusive file.
                for nj in range(1, 5):
                    masknj = masterArray["LHE_Njets"]==nj
                    if debug:
                        print("bare_name: {}".format(bare_name))
                        print("jetWeightMultiplicities: {}".format(jetWeightMultiplicity))
                    norm1 = jetWeightMultiplicity[bare_name]/HAA_processes[bare_name].weights["xsec"]
                    njname = "DY%dJetsToLL"%nj if bare_name == "DYJetsToLL" else "W%dJetsToLNu"%nj
                    norm2 = sow / HAA_processes[njname].weights["xsec"]
                    #norm2 also needs to be divided by the kfactor.
                    norm2 *= 1.0 / float(HAA_processes[njname].weights["kfactor"])
                    weightnj = weightfinal * 1.0 / (norm1 + norm2)
                    wf_arr += weightnj*masknj.astype(int)
                    
            #non-special processes get normal weighting lumi(already in weightfinal) * xsec / sum of weights
            else:
            #    print("block 3!")
                weightfinal *= HAA_processes[nickname].weights["xsec"] / sow
            #    print("weightfinal: {}".format(weightfinal))
                wf_arr *= weightfinal

         #   print("first 200 final weights: {}".format(wf_arr[:min(200,len(wf_arr))]))

            #multiply by weightfinal array.
            masterArray["finalweight"] *= wf_arr #weightfinal
         #   print("very final weights: {}".format(masterArray["finalweight"][:20])) #[:min(200,len(masterArray["finalweight"]))]))
        #    print "finalweight before per event scaling ",masterArray["finalweight"][:100]

         #   print("past very final weights.")
            #print the very high weights.
            bigweightmask = masterArray["finalweight"]>5
         #   print("bigweightmask!!!")
         #   print masterArray["evt"][bigweightmask]
         #   print(masterArray["finalweight"][bigweightmask])
         #   print("weightfinal: {}".format(weightfinal))
            
            #eventWeightDict = process.eventWeights
            eventWeightDict = EventWeights
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
                        #print "subtracting fakes "
                        #print weightMask[:100]
                        pass

                    masterArray["finalweight"] *= weightMask
           # print "finalweight after per event scaling ",masterArray["finalweight"][:100]


         #   print("before skim", len(masterArray["finalweight"]))
            for key,value in masterArray.iteritems():
                if ( key in plottedVars or key == "finalweight" ) \
                    and (len(mask)==len(value)):
                    skimArray[key] = value[mask]
            #for key in masterArray.keys():
                #try:
                    #skimArray[key] = masterArray[key][mask]
                #    skimArray[key] = value[mask]
                #except:
                    #print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                    #print "skipping branch ",key
                #    continue
         #   print("after skim", len(skimArray["mll"]), process.file)
            #for key in skimArray.keys():
            #    if key not in plottedVars and key != "finalweight":
            #        del skimArray[key]
            #print "working on category ",cat

            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray
    return skimArrayPerCat



#def slimskim(process,allcats,weightHistoDict,systematic):
#
#    skimArrayPerSysCats={}
#    #print "working on systematic ",systematic
#    work_dict = {}
#    #fin = uproot.open(process.file)
#    with uproot.open(process.file) as fin:
#
#        try:
#            tree = fin[systematic]
#        except:
#            return skimArrayPerSysCats
#
#        if systematic!="Events":
#            syst_names = set(fin[systematic].keys())
#            nom_names = set(fin["Events"].keys()) - syst_names
#            work_dict.update(fin[systematic].arrays(list(syst_names)))
#            work_dict.update(fin["Events"].arrays(list(nom_names)))
#            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))
#        else:
#            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal"))
#
#    del work_dict
#    print("deleted work_dict")
#    del tree
#    return skimArrayPerSysCats
#


#This is the one we're actually usuing!!!!!!!!!!
def slimskimoutput(process,allcats,weightHistoDict,systematic,massoutputdir):

    skimArrayPerSysCats={}
    #print "working on systematic ",systematic
    work_dict = {}
    #fin = uproot.open(process.file)
    with uproot.open(process.file) as fin:
        try:
           tree = fin[systematic]
        except:
           return

        if systematic!="Events":
            syst_names = set(fin[systematic].keys())
            nom_names = set(fin["Events"].keys()) - syst_names
            work_dict.update(fin[systematic].arrays(list(syst_names)))
            work_dict.update(fin["Events"].arrays(list(nom_names)))
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))
        else:
            #try to skrrrt around the MemoryError
            try:
                #print("About to call makeCutsOnTreeArray.")
                skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal"))
                #print("Now done with makeCutsOnTreeArray.")
            except MemoryError:
                print("Error! {} failed.".format(process.file))
                open("memErrFiles.txt", "a").write(process.file + "\n")


        createSlimOutput(skimArrayPerSysCats,massoutputdir)
    del skimArrayPerSysCats
    del work_dict
    del tree
    return

def createSlimOutput(skimArrayPerSysCats,outputdir):
   try:
      os.mkdir(outputdir)
   except:
      #print "dir exists"
      pass

   for key,dictionary in skimArrayPerSysCats.iteritems():
      dataTypes =[[],[]]
      for branch in dictionary.keys():
          dataTypes[0].append(branch)
          dataTypes[1].append(dictionary[branch].dtype)

      data = np.zeros(len(dictionary[branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})

      #filename = key.replace(":","_")
      filename = key.replace(":",".")
      fileout = ROOT.TFile.Open(outputdir+"/"+filename,"recreate")
      fileout.cd()
      for branch in data.dtype.names:
          if len(dictionary[branch].shape) == 1:   # flat array important for jagged arrays of input data
              data[branch] = dictionary[branch]
          else:
              data[branch] = dictionary[branch][:,0]

      #skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray
      root_numpy.array2tree(data,name=key.split(":")[0]+"_"+key.split(":")[-1])
      fileout.Write()
      fileout.Close()
      del data
      del fileout
          #treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist)
   return


def slimskimsys(process,allcats,weightHistoDict):
    import multiprocessing as mp

    skimArrayPerSysCats={}
    fin = uproot.open(process.file)
    try:
        tree = fin["Events"]
    except:
        return skimArrayPerSysCats
    skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal"))
    #print skimArrayPerSysCats.keys()
  #  try:
  #      print "final weight after nominal for Nominal:mmtt_inclusive:ZZTo4L:ZZ",skimArrayPerSysCats["Nominal:mmtt_inclusive:ZZTo4L:ZZ"]["finalweight"][:100]
  #  except:
  #      print "Nominal:mmtt_inclusive:ZZTo4L:ZZ doesn't exist for sample"
    #treatment of systematics

    systematics =[  ] #,   ##should Events be listed or nah??
                 #   "scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down","scale_m_eta1p2to2p1Up",
                 #  "scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down","scale_t_1prongUp","scale_t_1prongDown",
                 #  "scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown","scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp",
                 #  "scale_t_3prong1pizeroDown"]

    for systematic in systematics:
        #print "working on systematic ",systematic
        work_dict = dict()
        try:
            tree = fin[systematic]
        except:
            return skimArrayPerSysCats
        syst_names = set(fin[systematic].keys())
        nom_names = set(fin["Events"].keys()) - syst_names
        work_dict.update(fin[systematic].arrays(list(syst_names)))
        work_dict.update(fin["Events"].arrays(list(nom_names)))

        skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))



    #print skimArrayPerSysCats.keys()


    #pool  = mp.Pool(2)
    ##nums=pool.map(fstar,payloads)
#
    ##skims = slimskimstar(payloads[0])
    #sytematics = pool.map(slimskimstar,payloads)
#
    #pool.close()
    #pool.join()
    ##print(nums)





    return skimArrayPerSysCats












def createOutputSystematics(skimmedArraysSet,finalDistributions):
    import root_numpy

    finalSkims={}
    cats = []
    systematics = []
    #print skimmedArrays

    #print skimmedArrays.keys()

    for skimmedArrays in skimmedArraysSet:
        for sysCatNicPro in skimmedArrays.keys():
            #print "key ",sysCatNicPro
            sys = sysCatNicPro.split(":")[0]
            systematics.append(sys)
            channel = sysCatNicPro.split(":")[1].split("_")[0]
            #print "channel ",channel
            #print "systematic ",sys
            cat = sysCatNicPro.split(":")[1]
            cats.append(cat)
            if sys in finalSkims:
                if not cat in finalSkims[sys]:
                    finalSkims[sys][cat]={}
            else:
                finalSkims[sys]={}
                if not cat in finalSkims[sys]:
                    finalSkims[sys][cat]={}
            process= sysCatNicPro.split(":")[-1]
            #print "cat ",cat
            #print "process ",process

            for catDist in finalDistributions.keys():
                for processOut in finalDistributions[catDist]:
                    #if (catDist not in finalSkims[sys][cat]):
                    #print " is catDist ",catDist,"  inside keys ? ",finalSkims[sys][cat].keys()
                    if (processOut==process) and (catDist not in finalSkims[sys][cat].keys()):
                        #finalSkims[sys][cat][catDist] = processSkims[process]
                        #print "first output to finalskims ", sysCatNicPro,"  for process ",process," finalDist cat ",catDist
                        finalSkims[sys][cat][catDist] = skimmedArrays[sysCatNicPro]
                    #elif (catDist in finalSkims[sys][cat]):
                    elif (processOut==process) and (catDist in finalSkims[sys][cat].keys()):
                        #print "adding to finalskims ", sysCatNicPro,"  for process ",process," finalDist cat ",catDist
                        for branch in finalSkims[sys][cat][catDist].keys():
                            #finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],processSkims[process][branch]))
                            finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],skimmedArrays[sysCatNicPro][branch]))
                    else:
                        continue

    #print finalSkims.keys()
    #print finalSkims["Nominal"].keys()

    for cat in cats:
        localtest="allsys"
        skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
        skimFile.cd()
        for sys in finalSkims.keys():
            dataTypes =[[],[]]
            #print finalSkims[sys][cat].values()
            random_sample = finalSkims[sys][cat].values()[0]
            for branch in random_sample.keys():
                dataTypes[0].append(branch)
                dataTypes[1].append(random_sample[branch].dtype)
            for catDist in finalSkims[sys][cat].keys():
                #print "on the final dist ",catDist
                data = np.zeros(len(finalSkims[sys][cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
                for branch in data.dtype.names:
                    #print "working on branch ",branch
                    #print "branch datatype ",type(finalSkims[sys][cat][catDist][branch])
                    if len(finalSkims[sys][cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
                        data[branch] = finalSkims[sys][cat][catDist][branch]
                    else:
                        data[branch] = finalSkims[sys][cat][catDist][branch][:,0]
                treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist)
                treeOut.Write()
        skimFile.Close()

    return 1

#def slimskimorg(process,allcats,weightHistoDict):
#
#    from utils.functions import functs
#    from utils.Weights import CommonWeights
#    from ROOT import gInterpreter
#    commonweight = CommonWeights["lumi"][0]
#    fin = uproot.open(process.file)
#    tree = fin["Events"]
#    #print(fin)
#    #print(tree)
#
#    #for cat, catInfo in allcats.iteritems():
#    #    for variableHandle in catInfo.vars.keys():
#    for cat in allcats.keys():
#        masterArray = tree.arrays()
#        #print(process)
#        #print(process.cuts)
#        #print(masterArray.keys())
#        plottedVars = []
#
#        EventWeights = getEventWeightDicitonary()
#
#        newVarVals={}
#        skimArrayPerCat = {}
#        masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)
#
#        for variableHandle in allcats[cat].vars.keys():
#            variable = allcats[cat].vars[variableHandle][0]
#            if "[" in variable:
#                #print "adding jagged variable to array ",variable
#                basevar = variable.split("[")[0]
#                index = int(variable.split("[")[1].split("]")[0])
#                val = masterArray[basevar][:,index]
#                #masterArray[basevar+"_"+str(index)]=val
#                #plottedVars.append(basevar+"_"+str(index))
#                masterArray[variableHandle+"_"+str(index)]=val
#                plottedVars.append(variableHandle+"_"+str(index))
#            else:
#                plottedVars.append(variableHandle)
#
#
#        for var in allcats[cat].newvariables.keys():
#            newVarVals[var]=0.0
#
#        cuts=[]
#        for cuttype in allcats[cat].cuts.keys():
#            for cut in allcats[cat].cuts[cuttype]:
#                cuts.append(cut)
#
#
#        for var in newVarVals.keys():
#            arguments = allcats[cat].newvariables[var][2]
#            tempvals=[]
#            for ag in arguments:
#                tempvals.append(returnArray(masterArray,ag))
#            masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
#
#        mask = cutOnArray(masterArray,cuts)
#        masterArray["mask"]=mask
#        masterArray["finalweight"] *= mask.astype(int)
#        weightfinal = 1.0   #don't weight the data!!
#        #print "finalweight right after mask ",masterArray["finalweight"][:100]
#
#        weightDict = process.weights
#        #print " weight dicitonary ? ",weightDict
#        #print "sum of weights dicitonary ? ",weightHistoDict
#        weightfinal = commonweight
#        nickname = process.nickname
#        #print "working on the process nickname ",nickname
#        #print "sum of weights dicitonary keys ? ",weightHistoDict.keys()
#
#        for scalefactor in weightDict.keys():
#            if scalefactor == "kfactor":
#                weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
#            elif scalefactor in ["PU"]:
#                masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
#            elif scalefactor =="theoryXsec":
#                weightfinal =  weightfinal * float(weightDict[scalefactor])
#
#        #print "finalweight after PU and kFactor ",masterArray["finalweight"][:100]
#
#        if nickname =="DY1JetsToLL":
#            norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"]/HAA_processes["DYJetsToLL_ext1"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["DY1JetsToLL"]/HAA_processes["DY1JetsToLL"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        if nickname =="DY2JetsToLL":
#            norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"]/HAA_processes["DYJetsToLL_ext1"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["DY2JetsToLL"]/HAA_processes["DY2JetsToLL"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        if nickname =="DY3JetsToLL":
#            norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"]/HAA_processes["DYJetsToLL_ext1"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["DY3JetsToLL"]/HAA_processes["DY3JetsToLL"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        if nickname =="DY4JetsToLL":
#            norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"]/HAA_processes["DYJetsToLL_ext1"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["DY4JetsToLL"]/HAA_processes["DY4JetsToLL"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        if nickname =="W1JetsToLNu":
#            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["W1JetsToLNu"]/HAA_processes["W1JetsToLNu"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        elif nickname =="W2JetsToLNu":
#            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["W2JetsToLNu"]/HAA_processes["W2JetsToLNu"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        elif nickname =="W3JetsToLNu":
#            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
#            norm2 = jetWeightMultiplicity["W3JetsToLNu"]/HAA_processes["W3JetsToLNu"].weights["xsec"]
#            weightfinal = weightfinal * 1/(norm1+norm2)
#        else:
#            tempDenom=1.0
#            sumOfWeights = 1.0
#            for nic in weightHistoDict.keys():
#                if nickname==nic:
#                    sumOfWeights = weightHistoDict[nickname].GetSumOfWeights()
#                if nickname==nic.strip("_ext")[0]:
#                    sumOfWeights*=tempDenom*weightHistoDict[nic].GetSumOfWeights()
#            weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights
#        print "xsec/SoW ",HAA_processes[nickname].weights["xsec"]/ sumOfWeights
#        
#
#        #multiply by scalar weight
#        masterArray["finalweight"] *= weightfinal
#        #print "finalweight before per event scaling ",masterArray["finalweight"][:100]
#
#
#        eventWeightDict = process.eventWeights
#        if eventWeightDict:
#            for scalefactor in eventWeightDict.keys():
#
#                cutlist = eventWeightDict[scalefactor][0]
#                weightMask=cutOnArray(masterArray,cutlist)
#                weightMask=weightMask.astype(float)
#
#                if type(eventWeightDict[scalefactor][1][0])==float:
#                    weightMask *= eventWeightDict[scalefactor][1][0]
#
#                if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
#                    arguments = eventWeightDict[scalefactor][1][1]
#                    tempvals=[]
#
#                    for ag in arguments:
#                       tempvals.append(returnArray(masterArray,ag))
#                    weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)
#
#                if scalefactor!="fake" and scalefactor!="fake1" and scalefactor!="fake2":
#                    weightMask[np.where(weightMask==0.0)]=1.0
#                else:
#                    #print "subtracting fakes "
#                    #print weightMask[:100]
#                    pass
#
#                masterArray["finalweight"] *= weightMask
#        #print "finalweight after per event scaling ",masterArray["finalweight"][:100]
#        skipEvents = np.where(mask==0)[0]
#        skimArray={}
#        #print("before skim", len(masterArray["finalweight"]))
#        for key in masterArray.keys():
#            skimArray[key] = masterArray[key][mask]
#        #print("after skim", len(skimArray["mll"]))
#
#        for key in skimArray.keys():
#            if key not in plottedVars and key != "finalweight":
#                del skimArray[key]
#        skimArrayPerCat[cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray
#
#
#    return skimArrayPerCat

#def createOutput(skimmedArrays,finalDistributions):
#    import root_numpy
#
#    finalSkims={}
#    cats = []
#    for skimArrayPerCat in skimmedArrays:
#        #print skimArrayPerCat.keys()
#        for catNicPro in skimArrayPerCat.keys():
#            #print "key ",catNicPro
#            channel = catNicPro.split(":")[0].split("_")[0]
#            #print "channel ",channel
#            cat = catNicPro.split(":")[0]
#            cats.append(cat)
#            try:
#                finalSkims[cat]={}
#            except:
#                continue
#            process= catNicPro.split(":")[-1]
#            #print "cat ",cat
#            #print "process ",process
#
#            for catDist in finalDistributions.keys():
#                for processOut in finalDistributions[catDist]:
#                    #if (catDist not in finalSkims[cat]):
#                    if (processOut==process) and (catDist not in finalSkims[cat]):
#                        #finalSkims[cat][catDist] = processSkims[process]
#                        #print "first output to finalskims ", catNicPro,"  for process ",process," finalDist cat ",catDist
#                        finalSkims[cat][catDist] = skimArrayPerCat[catNicPro]
#                        continue
#                    #elif (catDist in finalSkims[cat]):
#                    elif (processOut==process) and (catDist in finalSkims[cat]):
#                        #print "adding to finalskims ", catNicPro,"  for process ",process," finalDist cat ",catDist
#                        for branch in finalSkims[cat][catDist].keys():
#                            #finalSkims[cat][catDist][branch]=np.concatenate((finalSkims[cat][catDist][branch],processSkims[process][branch]))
#                            finalSkims[cat][catDist][branch]=np.concatenate((finalSkims[cat][catDist][branch],skimArrayPerCat[catNicPro][branch]))
#                    else:
#                        continue
#
#
#    for cat in cats:
#        localtest="ZZ"
#        skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
#        skimFile.cd()
#
#        dataTypes =[[],[]]
#        #print finalSkims[cat].values()
#        random_sample = finalSkims[cat].values()[0]
#        for branch in random_sample.keys():
#            dataTypes[0].append(branch)
#            dataTypes[1].append(random_sample[branch].dtype)
#        for catDist in finalSkims[cat].keys():
#            #print "on the final dist ",catDist
#            data = np.zeros(len(finalSkims[cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
#            for branch in data.dtype.names:
#                #print "working on branch ",branch
#                #print "branch datatype ",type(finalSkims[cat][catDist][branch])
#                if len(finalSkims[cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
#                    data[branch] = finalSkims[cat][catDist][branch]
#                else:
#                    data[branch] = finalSkims[cat][catDist][branch][:,0]
#            treeOut = root_numpy.array2tree(data, name=catDist)
#            treeOut.Write()
#        skimFile.Close()
#
#    return 1

def combineRootFiles(systematics, finalDistributions, rootfiledir, channel):
   import os
   import glob
   import uproot
   rootFiles = {}
#   print "the final distributions ",finalDistributions
   #cat = "mmmt_inclusive"
   finalSkims={}

#   print("boutta start systematics loop!!")
   for sys in systematics:
    #  print("now running on systematic: {}".format(sys))
      if sys == "Events":
    #     print("rootfiledir: {}".format(rootfiledir))
         #rootFiles["Nominal"] = glob.glob(rootfiledir +"/Nominal_*")
         rootFiles["Nominal"] = glob.glob(rootfiledir +"/Nominal.*")
         finalSkims["Nominal"]={}
      else:
        rootFiles[sys] = glob.glob(rootfiledir+"/"+sys+".*")
        finalSkims[sys]={}


   mainTList = {}
   mainOutputTree={}
   #print("rootFiles iteritems: {}".format(str(rootFiles.iteritems())))
   for sys, globfiles in rootFiles.iteritems():
     # print("all globfiles for sys {}: {}".format(sys, globfiles))
      for globfile in globfiles:
#         process = "".join(globfile.split('/')[-1].split("_")[:-1])
         #process = globfile.split("_")[-1]
         process = globfile.split(".")[-1]
         print("file: {}, process: {}".format(globfile, process))
         with uproot.open(globfile) as fin:
            majorkey = sys + "_" + process
    #        if sys == "Nominal":
    #            majorkey = "Events"
            #print("**********MAJOR KEY**************: {}".format(majorkey))
            try:
                tree = fin[majorkey]
            except:
                print("file: {}".format(globfile))
                sys.exit("Error!!!!!!!!!")
            mainArrays = tree.arrays()
            #print("globfile {}".format(globfile))
            for catDist, final in finalDistributions.iteritems():
         #       print("catDist: {}, final: {}".format(catDist, final))
                for processOut in final:
         #           print("processOut: {}".format(processOut))
                    if (processOut==process) and (catDist not in finalSkims[sys]):
                    #    print "first output for process ",process," finalDist cat ",catDist
                        finalSkims[sys][catDist] = mainArrays
                        continue
                    elif (processOut==process) and (catDist in finalSkims[sys]):
                    #    print "adding to finalskims ", catDist,"  for process ",process," finalDist cat ",catDist
                        for branch in finalSkims[sys][catDist].keys():
                            finalSkims[sys][catDist][branch]=np.concatenate((finalSkims[sys][catDist][branch],mainArrays[branch]))
                    else:
                        continue
   print("after first loop of combineRootFiles, finalSkims[\"Nominal\"]:")
   print(finalSkims["Nominal"].keys())
   skimoutname = "skimmed_"+channel+ "_" + args.outname+".root"
   skimFile = ROOT.TFile(skimoutname,"recreate")
   skimFile.cd()
   #print("cd'd to file {}!!!!!".format(skimoutname))

   for sys in rootFiles.keys():
       # print "systematic ",sys
        dataTypes =[[],[]]
       # print("*******************finalSkims[{}] values:**************************".format(sys))
       # print finalSkims[sys].values()
        random_sample = finalSkims[sys].values()[0]
       # print "random sample branches? ",random_sample.keys()
        for branch, sample in random_sample.iteritems():
            dataTypes[0].append(branch)
            dataTypes[1].append(sample.dtype)
        for catDist, finalSkim in finalSkims[sys].iteritems():
            #print "on the final dist ",catDist
            data = np.zeros(len(finalSkim[branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
            for branch in data.dtype.names:
                #print "working on branch ",branch
                #print "branch datatype ",type(finalSkims[sys][catDist][branch])
                if branch in finalSkim and len(finalSkim[branch].shape) == 1:   # flat array important for jagged arrays of input data
                    data[branch] = finalSkim[branch]
                elif branch in finalSkim:
                    data[branch] = finalSkim[branch][:,0]
            treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist)
            treeOut.Write()
   skimFile.Close()

   return

if __name__ == "__main__":
    import datetime
    begin_time = datetime.datetime.now()
    from multiprocessing import *
    import multiprocessing as mp
    import logging
    import os

    import argparse


    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    parser.add_argument("-o",  "--outname", default="test43",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
    parser.add_argument("-ch",  "--channel", default="mtet",  help="Please list the channel for fake factor histograms")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2017_v7_yaml.csv",  help="categories yaml file")
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
    parser.add_argument("-comb",  "--combineFiles", default=False, action='store_true', help="Combine root files mode (1) or generate root files mode (0)")
    parser.add_argument("-dbg",  "--debug", default=False, action='store_true', help="Switch on to turn off multithreading (makes finding errors easier).")
    parser.add_argument("-yr",  "--year", default=2017, help="which year")

    args = parser.parse_args()

    debug = args.debug
    year = int(args.year)
    direc = "/eos/uscms/store/user/bgreenbe/haa_4tau_%d/all/"%(year)    #_mmonly/"

    if args.combineFiles:
        print("Running in combineFiles mode.")
    else:
        print("Running in generateFiles mode.")
    allcats={}
    HAA_processes={}
    finalDistributions={}
    filelist = {}
    weightHistoDict={}
    jetWeightMultiplicity={}
    EventWeights={}


    allcats,HAA_processes,finalDistributions,\
    weightHistoDict,jetWeightMultiplicity = initialize(args.channel, direc)

    print("Done with initialize.")

    #print(allcats)
    info('main line')
    #print(f)
    nums=[]
    skims=[]
    #skims={}
    payloadsdict={}
    payloads=[]

    systematics =[ "Events" ] #,
                #    "scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
                #   "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
                #   "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
                #   "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]

    tempdir = "massOutputDir_"+args.channel+"_"+args.outname

    if not args.combineFiles:
###############################step 0###############################################################

        for nickname, process in HAA_processes.items():
            for syst in systematics:
                payloads.append((process,allcats,weightHistoDict,syst,tempdir))

        #print(payloads)

        m = mp.Manager()
      #  logger_q = m.Queue()
        #parallelable_data = [(1, logger_q), (2, logger_q)]

        pool  = mp.Pool(4)   #12)

        #skims = pool.map(slimskimstar,payloads)

        if not debug:
            pool.map(slimskimstar,payloads)#this works for root output!

            pool.close()
            pool.join()
          #  while not logger_q.empty():
          #      print logger_q.get()
        else:
        #don't multithread
            for pl in payloads:
                slimskimstar(pl)



        print("ready to combine the output! (if no prior errors)")
        print("prior error files would show up here:")
        os.system("cat memErrFiles.txt")
        ##if createOutput(skims,finalDistributions) : print "successful output"
        ##else: print "ouch ..."
        #if createOutputSystematics(skims,finalDistributions) : print "successful output"
        #else: print "ouch ..."
#######################################################################################################
######################################step 1###########################################################
    else:
        combineRootFiles(systematics, finalDistributions, tempdir, args.channel)
        #combineRootFiles(systematics, finalDistributions, direc, args.channel)
#######################################################################################################



    print("computation time")
    print(datetime.datetime.now() - begin_time)
