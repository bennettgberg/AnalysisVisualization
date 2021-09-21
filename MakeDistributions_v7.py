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
ROOT.gErrorIgnoreLevel = ROOT.kFatal
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


#use constant value instead of fitting for FF or nah (can be set True with -const)
const = False

#Setting up operators for cut string iterator
ops = { "==": operator.eq, "!=": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}


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

# def returnArray(masterArray,variable):
#     if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
#         val = masterArray[variable][:,0]
#     else:
#         val = masterArray[variable]
#     return val

#create mask for all the events that do or don't pass cuts
def cutOnArray(masterArray,cuts):

    mask=np.full(len(masterArray["evt"]),True)

    i = 0
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
            #test

            mask *= tempmask.astype(bool)
            continue

        #if cut[0][0]=="OR":
        #    tempmasks=[]
        #    tempmask=np.full(len(masterArray["evt"]),False)
        #    for oe in range(1,len(cut)):
        #        oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
        #        tempmask = ops["bor"](tempmask,oneor).astype(bool)
        #    mask *= tempmask.astype(bool)
        #    continue
        if cut[0][0]=="OR":
            tempmasks=[]
            tempmask=np.full(len(masterArray["evt"]),False)
            for oe in range(1,len(cut)):
                oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
                #tempmask = ops["bor"](tempmask,oneor).astype(bool)
                tempmask = (tempmask + oneor).astype(bool)
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
        # if mask.all()==False:
        #     print "cut that made it false ",cuts[i]
        #     return mask
        # i=i+1


    return mask


#for fake factor return the bin value from event content
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    newArr = np.vectorize(pt.FindBin)(numpyArr)
    newArr = np.vectorize(pt.GetBinContent)(newArr)

    return newArr

def getEventWeightDicitonary(year):
    #import copyreg, copy, pickle # for picking the bound methods due to the multiplrocess tool

    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool

    if year=="2016":
        yearhandle = "2016Legacy"
    if year=="2017":
        yearhandle = "2017ReReco"
    if year=="2018":
        yearhandle = "2018ReReco"
    tauIDSF = TauIDSFTool(yearhandle,'DeepTau2017v2p1VSjet','Medium').getSFvsPT
    antiEleSFToolVL = TauIDSFTool(yearhandle,'DeepTau2017v2p1VSe','VLoose').getSFvsEta
    antiMuSFToolVL  = TauIDSFTool(yearhandle,'DeepTau2017v2p1VSmu','VLoose').getSFvsEta
    antiEleSFToolT = TauIDSFTool(yearhandle,'DeepTau2017v2p1VSe','Tight').getSFvsEta
    antiMuSFToolT  = TauIDSFTool(yearhandle,'DeepTau2017v2p1VSmu','Tight').getSFvsEta

    EventWeights={
            "3_tauSF_8":[[["cat","==",8],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            "4_tauSF_8":[[["cat","==",8],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            "3_tauSF_7":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            "4_tauSF_7":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            "3_tauSF_6":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            "4_tauSF_6":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            "3_tauSF_5":[[["cat","==",5],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            "4_tauSF_5":[[["cat","==",5],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],

            # "3_muonSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolVL),["eta_3","gen_match_3"]]],
            # "4_muonSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolVL),["eta_4","gen_match_4"]]],
            # "3_muonSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolVL),["eta_3","gen_match_3"]]],
            # "4_muonSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolVL),["eta_4","gen_match_4"]]],
            # "3_muonSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolVL),["eta_3","gen_match_3"]]],
            # "4_muonSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolVL),["eta_4","gen_match_4"]]],
            # "3_muonSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolVL),["eta_3","gen_match_3"]]],
            # "4_muonSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolVL),["eta_4","gen_match_4"]]],
            #
            # "3_electronSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolVL),["eta_3","gen_match_3"]]],
            # "4_electronSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolVL),["eta_4","gen_match_4"]]],
            # "3_electronSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolVL),["eta_3","gen_match_3"]]],
            # "4_electronSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolVL),["eta_4","gen_match_4"]]],
            # "3_electronSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolVL),["eta_3","gen_match_3"]]],
            # "4_electronSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolVL),["eta_4","gen_match_4"]]],
            # "3_electronSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolVL),["eta_3","gen_match_3"]]],
            # "4_electronSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolVL),["eta_4","gen_match_4"]]]

            "3_muonSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolT),["eta_3","gen_match_3"]]],
            "4_muonSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolT),["eta_4","gen_match_4"]]],
            "3_muonSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolT),["eta_3","gen_match_3"]]],
            "4_muonSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolT),["eta_4","gen_match_4"]]],
            "3_muonSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolT),["eta_3","gen_match_3"]]],
            "4_muonSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolT),["eta_4","gen_match_4"]]],
            "3_muonSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(antiMuSFToolT),["eta_3","gen_match_3"]]],
            "4_muonSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(antiMuSFToolT),["eta_4","gen_match_4"]]],

            "3_electronSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolT),["eta_3","gen_match_3"]]],
            "4_electronSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolT),["eta_4","gen_match_4"]]],
            "3_electronSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolT),["eta_3","gen_match_3"]]],
            "4_electronSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolT),["eta_4","gen_match_4"]]],
            "3_electronSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolT),["eta_3","gen_match_3"]]],
            "4_electronSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolT),["eta_4","gen_match_4"]]],
            "3_electronSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(antiEleSFToolT),["eta_3","gen_match_3"]]],
            "4_electronSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(antiEleSFToolT),["eta_4","gen_match_4"]]]
            }
    # if year == "2016":
    #     EventWeights={
            #"name":[[if statements],[weight to apply]]
            # "3_mt_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[0.80]],
            # "3_mj_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[1.21]],
            # "3_mt_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[0.81]],
            # "3_mj_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[1.11]],
            # "3_mt_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[0.79]],
            # "3_mj_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[1.2]],
            # "3_mt_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[0.68]],
            # "3_mj_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[1.16]],
            # "3_mt_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[0.68]],
            # "3_mj_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[2.25]],
            #
            # "4_mt_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
            # "4_mj_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
            # "4_mt_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
            # "4_mj_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
            # "4_mt_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
            # "4_mj_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
            # "4_mt_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
            # "4_mj_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
            # "4_mt_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
            # "4_mj_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],
            #
            # "3_et_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
            # "3_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
            # "3_et_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
            # "3_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
            # "3_et_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
            # "3_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
            # "3_et_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
            # "3_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],
            #
            # "4_et_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[0.80]],
            # "4_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[1.18]],
            # "4_et_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.72]],
            # "4_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.93]],
            # "4_et_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.14]],
            # "4_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.18]],
            # "4_et_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[0.64]],
            # "4_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[1.07]],
            #
            # "8_3_et_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
            # "8_3_ej_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
            # "8_3_et_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
            # "8_3_ej_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
            # "8_3_et_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
            # "8_3_ej_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
            # "8_3_et_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
            # "8_3_ej_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],
            #
            # "8_4_mt_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
            # "8_4_mj_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
            # "8_4_mt_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
            # "8_4_mj_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
            # "8_4_mt_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
            # "8_4_mj_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
            # "8_4_mt_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
            # "8_4_mj_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
            # "8_4_mt_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
            # "8_4_mj_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

            #bound method ... last input list is func parameters! so this uses the tauIDSF function
            # "3_tauSF_8":[[["cat","==",8],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            # "4_tauSF_8":[[["cat","==",8],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            # "3_tauSF_7":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            # "4_tauSF_7":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            # "3_tauSF_6":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            # "4_tauSF_6":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
            # "3_tauSF_5":[[["cat","==",5],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
            # "4_tauSF_5":[[["cat","==",5],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]]

            #recoil corrections? not affecting current fit variable
            #"recoilcorr2":[[["cat","==",6],["njets","==",2]],[recoilCorrector,["met_x","met_y","gen_match_3"]]],
        # }

    return EventWeights

def initialize(args):
    import os
    from copy import copy
    import yaml

    #Structure for mapping between root files and processes
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights
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
    #categories = "cat_"+channel+"_2016.yaml"
    #processes = "processes_special_"+channel+".yaml"
    #categories = "cat_mmet_2016.yaml"
    #processes = "processes_special_mmet.yaml"
    csvfile = args.csvfile
    categories = args.categories
    processes = args.processes

    #dir = "/eos/home-s/shigginb/HAA_ntuples/2016_v7/"
    #dir = "/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/"
    dir = args.dir



    for line in open(csvfile,'r').readlines() :
            # structure for csv to sampleDict conversion
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            if "#" in line[0]: continue
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
        #if categories[category]['name']==args.channel+"_inclusive":
        #if categories[category]['name']=="mmet_inclusive":
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        #tempcat.systematics=categories[category]['systematics']
        allcats[tempcat.name]=tempcat
        #print allcats[tempcat.name].cuts
        # if args.extract and (args.channel=="mmmt" or args.channel=="mmet":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",120.0])
        # if args.extract and args.channel=="mmtt":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",130.0])
        # if args.extract and args.channel=="mmem":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",110.0])

    #allcats = {name: DylanCategory(cat) for name, cat in categories.items()}


    #for cat in allcats.keys():
        #print(allcats[cat].name)
        #print(allcats[cat].cuts)
        #print(allcats[cat].systematics)

    #loading standard processes
    HAA_processes={}

    #MC for background

    #HAA_process = {name: Process(dir, name, info) for name, info in sampleDict.items()}
    #weightHistoDict = {}


    if not (args.datadrivenZH):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            temppro.file=dir+sample+"_"+args.year+".root"

            frooin = ROOT.TFile.Open(temppro.file,"read")
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={sampleDict[sample][0]:""}
            if "HToAATo4Tau" in sample:
                if args.extract:
                    temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots

                else:
                    temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} # worked before
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
                #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001* 5.0)} # SM Higgs xsec x BR Haa  for signal extraction data MC control plots(AN 17029)
            HAA_processes[temppro.nickname]=temppro

    if (args.datadrivenZH):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            #temppro.file=sample+"_2016.root"
            #temppro.file=dir+sample+"_2016.root"
            temppro.file=dir+sample+"_"+args.year+".root"
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
            if args.channel=="mmmt" or args.channel=="mmet":
                truetau = [[["OR"],
                            ["gen_match_3","==",15],
                            ["gen_match_4","==",5]
                            ]]
            temppro.cuts={sampleDict[sample][0]:truetau} #ONLY SELECT PRMOPT FOR MC!
            if "HToAATo4Tau" in sample:
                #for visualization!
                if args.extract:
                    temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots

                else:
                    temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} # worked before
            HAA_processes[temppro.nickname]=temppro

    if (args.datameasureZH):
        try:
            os.mkdir("FFhistos_"+str(args.ffout))
        except:
            print("directory exists")
        for proObj in HAA_processes.keys():
            if proObj!="data_obs" and proObj!="FF" and "HToAATo4Tau" not in proObj:
                if args.channel=="mmmt" or args.channel=="mmet":
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",5]]
                    HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","!=",5]]
                    HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","!=",5]]

                if args.channel=="mmtt":
                    if "HToAATo4Tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",5]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",5]]
                    HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","!=",5]]
                    HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","!=",5]]
                if args.channel=="mmem":
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","!=",15]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","!=",15]]
                    if "HToAATo4Tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",15]]
                else: # default case is close enough
                    HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                    HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                    if "HToAATo4Tau" in proObj:
                        continue
                    HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                    HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
    elif (args.datadrivenZH):
        for proObj in HAA_processes.keys():
            if proObj!="data_obs" and proObj!="FF" and not args.skim:
                HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                if "HToAATo4Tau" in proObj:
                    print "skipping signal prompt mc"
                    continue  #this hasn't been added yet ... this is needed to omit signal!!
                HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]


    #loading special processes ... fake factor and data!
    for process in processes_special:
        temppro = Process()
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=dir+processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro
    #print "==========================================================================================="
    #print " All the processes !!!!!  ",HAA_processes
    #print "==========================================================================================="


    #precise weights for jet multiplicity
    from utils.Weights import jet_inclusive_samples
    from utils.Weights import jetIncOnly
    inc_samps = jetIncOnly[args.year]
    inclusive_samples = jet_inclusive_samples[args.year]

    #file list for easy reading
    for proObj in HAA_processes.keys():
        filelist[proObj]=HAA_processes[proObj].file

    # adding kfactor to relevent samples
    for sample in sampleDict.keys():
        if sample in inclusive_samples and sample.startswith("DY"):
            HAA_processes[sample].weights.update({"kfactor":1.1637})
        if sample in inclusive_samples and sample.startswith("W"):
            HAA_processes[sample].weights.update({"kfactor":1.221})

    weightHistoDict = {}
    for nickname, filename in filelist.iteritems():
        #frooin = ROOT.TFile(dir+filename)
        frooin = ROOT.TFile(filename)
        hists = None
        try:
            if nickname.startswith("DYJetsToLL"):
                hists=[frooin.Get("hWeights").Clone(),
                       frooin.Get("DY1genWeights").Clone(),
                       frooin.Get("DY2genWeights").Clone(),
                       frooin.Get("DY3genWeights").Clone(),
                       frooin.Get("DY4genWeights").Clone(),
                       ]
            elif nickname.startswith("WJetsToLNu"):
                hists=[frooin.Get("hWeights").Clone(),
                      frooin.Get("W1genWeights").Clone(),
                      frooin.Get("W2genWeights").Clone(),
                      frooin.Get("W3genWeights").Clone(),
                      frooin.Get("W4genWeights").Clone(),
                      ]
            else:
                hists = frooin.Get("hWeights").Clone()
                hists.SetDirectory(0)
            if isinstance(hists, list):
                for hist in hists:
                    hist.SetDirectory(0)
            weightHistoDict[nickname] = hists
            frooin.Close()
        except:
            print nickname," hWeights prob doesn't exist?"

    # print weightHistoDict
    #print weightHistoDict["DYJetsToLLext1"]
    #print weightHistoDict["WJetsToLNu"]


    jetWeightMultiplicity = {}
    NJetWeights = {}
    #print inc_samps
    for inc_group in inc_samps:
        sumOfWeights = np.zeros(5)
        for sample in inc_group:
            weights = weightHistoDict[sample]
            sumOfWeights += [w.GetSumOfWeights() for w in weights]
        sumOfWeights[sumOfWeights==0] = 1e-10 # if value is 0 to stop infs
        #print(sumOfWeights)
        for sample in inc_group:
            jetWeightMultiplicity[sample] = sumOfWeights


    if args.year=="2018":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

    if args.year=="2017":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

    if args.year=="2016":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
        W2JetsFileext1 = ROOT.TFile.Open(filelist["W2JetsToLNu_ext1"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]+=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W2JetsToLNu_ext1"]=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W2JetsToLNu_ext1"]+=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
        W3JetsFileext1 = ROOT.TFile.Open(filelist["W3JetsToLNu_ext1"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]+=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W3JetsToLNu_ext1"]=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W3JetsToLNu_ext1"]+=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        W4JetsFileext1 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext1"],"read")
        W4JetsFileext2 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext2"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        jetWeightMultiplicity["W4JetsToLNu_ext1"]=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        jetWeightMultiplicity["W4JetsToLNu_ext2"]=W4JetsFileext2.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()


    Bkg = ["DY","W","TT","ST","EWK","Top"]
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
    if args.datameasureZH:
        finalDistributions["prompt1"]=["prompt1"]
        finalDistributions["prompt2"]=["prompt2"]
    if args.datadrivenZH:
        Bkg = ["FF"]
        finalDistributions["Bkg"]=Bkg

    allSkims = {}
    finalSkims ={}

    datadrivenPackage={}
    datadrivenPackage["bool"]=False


    if args.datadrivenZH or args.makeFakeHistos:
        print "MUST HAVE CREATED FF DISTRIBUTIONS BEFORE!"
        histodict = {}
        if not const:
            #inputFFile = ROOT.TFile.Open("FFskim_"+str(args.ffin)+".root","read")
            with uproot.open("skimmed_"+str(args.ffin)+".root") as inputFFile:
                histodict = createFakeFactorHistos(allcats, inputFFile)
            #inputFFile.Close()
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadrivenZH
        if not const:
            ff_file_3 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_3_ff.root","read")
            ff_file_4 = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/pt_4_ff.root","read")

            # ss_1_tight = ff_file_3.Get(args.channel+"_FF_SS_1_tight/data_obs")
            # ss_1_loose = ff_file_3.Get(args.channel+"_FF_SS_1_loose/data_obs")
            # ss_2_tight = ff_file_4.Get(args.channel+"_FF_SS_2_tight/data_obs")
            # ss_2_loose = ff_file_4.Get(args.channel+"_FF_SS_2_loose/data_obs")

            ss_1_tight = histodict[args.channel+"_FF_SS_1_tight"]["Nominal_data_obs"]["pt_3_ff"]
            ss_1_loose = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_3_ff"]
            ss_2_tight = histodict[args.channel+"_FF_SS_2_tight"]["Nominal_data_obs"]["pt_4_ff"]
            ss_2_loose = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_4_ff"]

            # ss_1_tight_prompt = ff_file_3.Get(args.channel+"_FF_SS_1_tight/prompt1")
            # ss_1_loose_prompt = ff_file_3.Get(args.channel+"_FF_SS_1_loose/prompt1")
            # ss_2_tight_prompt = ff_file_4.Get(args.channel+"_FF_SS_2_tight/prompt2")
            # ss_2_loose_prompt = ff_file_4.Get(args.channel+"_FF_SS_2_loose/prompt2")

            ss_1_tight_prompt = histodict[args.channel+"_FF_SS_1_tight"]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_prompt = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_prompt1"]["pt_3_ff"]
            ss_2_tight_prompt = histodict[args.channel+"_FF_SS_2_tight"]["Nominal_prompt2"]["pt_4_ff"]
            ss_2_loose_prompt = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_prompt2"]["pt_4_ff"]

            #subtracting prompt MC execpt for low stat channel
            # if args.channel!="mmem":
            ss_1_tight.Add(ss_1_tight_prompt,-1)
            ss_2_tight.Add(ss_2_tight_prompt,-1)
            ss_1_loose.Add(ss_1_loose_prompt,-1)
            ss_2_loose.Add(ss_2_loose_prompt,-1)


            f_1= ss_1_tight.Clone()
            f_1.Divide(ss_1_loose)
            f_2 = ss_2_tight.Clone()
            f_2.Divide(ss_2_loose)

            f_1.GetYaxis().SetTitleOffset(1.4)
            f_2.GetYaxis().SetTitleOffset(1.4)
            f_1.GetYaxis().SetMaxDigits(2)
            f_2.GetYaxis().SetMaxDigits(2)
            #ROOT.TGaxis().SetMaxDigits(2)

            f_1.SetName(args.channel+" FakeRateLeg1")
            f_1.SetTitle(args.channel+" Fake Rate Measurement Leg1")
            f_1.GetXaxis().SetTitle("p_T Leg1")
            f_1.GetYaxis().SetTitle("Fake Rate for Leg1")
            f_2.SetName(args.channel+" FakeRateLeg2")
            f_2.SetTitle(args.channel+" Fake Rate Measurement Leg2")
            f_2.GetXaxis().SetTitle("p_T Leg2")
            f_2.GetYaxis().SetTitle("Fake Rate for Leg2")

            tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())
            tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax())

            print("opening fakemeasurefile.")
            fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/fakemeasure.root","RECREATE")
            fakemeasurefile.cd()
            c=ROOT.TCanvas("canvas","",0,0,600,600)
            #ROOT.gStyle.SetOptFit()
            ROOT.gStyle.SetOptStat(0)
            #ROOT.gStyle.SetOptFit(1)
            if 0 == 1: #not "mmem" in args.ffin:
                f_1.Draw()
            else:
                f_1.Draw("PE")
            f_1.Fit("tf_1")
            c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate1.png")
            if 0 == 1: # not "mmem" in args.ffin:
                f_2.Draw()
            else:
                f_2.Draw("PE")
            f_2.Fit("tf_2")
            c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate2.png")
            tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
            tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
            f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
            f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)
            datadrivenPackage["fakerate1"]=f_1.Clone()
            datadrivenPackage["fakerate2"]=f_2.Clone()
            datadrivenPackage["fitrate1"]=tf_1
            datadrivenPackage["fitrate2"]=tf_2
            datadrivenPackage["fakemeasurefile"]=fakemeasurefile
            #fakemeasurefile.Close()

            ##print("WARNING: exiting before doing anything besides making the FFhistos files!!!!!")
            ##sys.exit()

    #EventWeights = getEventWeightDicitonary()

    #exit()
    # ROOT.fail
    return allcats, HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity,datadrivenPackage





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

def createFakeFactorHistos(allcats, inputFFile):
    newVarVals={}
    treetypes = ["Nominal_data_obs","Nominal_prompt1","Nominal_prompt2"]
    # for cat in allcats.keys():
    #     if not "_inclusive" in cat:
    #         histodict={}
    #         for treename in treetypes:
    #             histodict[cat]={}
    #             histodict[cat][treename]={}
    histodict = {c:{t: dict() for t in treetypes} for c in allcats.keys()}
    #print "histodict ",histodict


    #creating the Fake Factor histograms from pre-defined numpy arrays
    for cat in histodict.keys():
        #print cat
        if "inclusive" in cat: continue
        for treename in treetypes:
            tree = inputFFile[cat][treename]
            fakefactorArray = tree.arrays()
           # print "cats ",cat,"  tree name ",treename
           # print "entries ",len(fakefactorArray["finalweight"])
            #vars = fakefactorArray.keys()
            vars = allcats[cat].vars
            for variableHandle in vars:
                #print variableHandle
                var = allcats[cat].vars[variableHandle][0]
                #var = variableHandle
                #print var,"  ",variableHandle
             #   if "jpt" in variableHandle and "ff" in variableHandle:            #why this here??????????????
             #       firstnum = variableHandle.split("_")[1]
             #       secondnum = str(int(firstnum)-1)
             #       var = variableHandle+"_"+secondnum
             #       #print var
                bins = allcats[cat].vars[variableHandle][1]
                if type(bins[0])==list:
                    histodict[cat][treename][variableHandle] = ROOT.TH1D(str(variableHandle),str(variableHandle),bins[0][0],bins[0][1],bins[0][2])
                    val = fakefactorArray[var]
                    #val = returnArray(fakefactorArray,var)
                    root_numpy.fill_hist(histodict[cat][treename][variableHandle],val,fakefactorArray["finalweight"])
                else:
                    tmpbin = np.asarray(bins)
                    histodict[cat][treename][variableHandle] = ROOT.TH1D(str(variableHandle),str(variableHandle),len(tmpbin)-1,tmpbin)
                    val = fakefactorArray[var]
                    #val = returnArray(fakefactorArray,var)
                    root_numpy.fill_hist(histodict[cat][treename][variableHandle],val,fakefactorArray["finalweight"])
    #now I have the histograms ... so using standard methods
    return histodict





def makeCutsOnTreeArray(processObj, inputArray,allcats,weightHistoDict,systematic,args):

    from utils.functions import functs
    from utils.Weights import CommonWeights
    from utils.Weights import jet_inclusive_samples

    inclusive_samples = jet_inclusive_samples[args.year]

    from ROOT import gInterpreter
    import copy
    commonweight = CommonWeights["lumi"+args.year][0]
    skimArrayPerCat = {}
    #print "working on process obj ",processObj.nickname

    for process in processObj.cuts.keys():
        procut = processObj.cuts[process]
        print "working on process ",process
        for cat in allcats.keys():
            masterArray = inputArray.copy()
            cuts=[]
            #print "working on category ",cat
            #print "with vars ",allcats[cat].vars
            #print "starting length of dictionary ",len(masterArray["mll"])

            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)
            if procut!="":
                cuts.append(procut[0])
            if args.extract and args.channel=="mmmt":
                cuts.append(["AMass","<=",120.0])
                cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmtt":
                cuts.append(["AMass","<=",130.0])
                cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmet":
                cuts.append(["AMass","<=",120.0])
                cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmem":
                cuts.append(["AMass","<=",110.0])
                cuts.append(["mll-mtt",">",0.0])
            #print(process)
            #print(process.cuts)
            #print(masterArray.keys())
            plottedVars = []
            newVarVals={}
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
                    plottedVars.append(variable)
                else:
                    plottedVars.append(variableHandle)
                    plottedVars.append(variable)

            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0
                plottedVars.append(var) # new vars just label
            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)


            if process=="data_obs":
                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

                #Comments with Dylan!
                #move this out of the loop! :)
                #perhaps combine the category and process object
                #for the skimArray[] = skimArray return ... write function to get the key?



                mask = cutOnArray(masterArray,cuts)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!

                skipEvents = np.where(mask==0)[0]
                skimArray={}
                #print("before skim", len(masterArray["finalweight"]))
                for key in masterArray.keys():
                    try:
                        skimArray[key] = masterArray[key][mask]
                    except:
                        print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                        print "skipping branch ",key
                        continue
                #print("after skim", len(skimArray["mll"]), processObj.file)
                if len(skimArray["mll"])==0:
                    continue

                for key in skimArray.keys():
                    if key not in plottedVars and key != "finalweight":
                        del skimArray[key]
                skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray




            if(process=="FF" and datadrivenPackage["bool"]):
                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)


                tempmask=np.full(len(masterArray["evt"]),1.0)

                #the actual events that pass the FF_1 criteria
                cuts_1 = HAA_processes["FF"].cuts["FF_1"]
                cuts_2 = HAA_processes["FF"].cuts["FF_2"]
                cuts_12 = HAA_processes["FF"].cuts["FF_12"]
                if args.extract and (args.channel=="mmmt" or args.channel=="mmet"):
                    cuts_1.append(["AMass","<=",120.0])
                    cuts_1.append(["mll-mtt",">",0.0])
                    cuts_2.append(["AMass","<=",120.0])
                    cuts_2.append(["mll-mtt",">",0.0])
                    cuts_12.append(["AMass","<=",120.0])
                    cuts_12.append(["mll-mtt",">",0.0])
                if args.extract and args.channel=="mmtt":
                    cuts_1.append(["AMass","<=",130.0])
                    cuts_1.append(["mll-mtt",">",0.0])
                    cuts_2.append(["AMass","<=",130.0])
                    cuts_2.append(["mll-mtt",">",0.0])
                    cuts_12.append(["AMass","<=",130.0])
                    cuts_12.append(["mll-mtt",">",0.0])
                if args.extract and args.channel=="mmem":
                    cuts_1.append(["AMass","<=",110.0])
                    cuts_1.append(["mll-mtt",">",0.0])
                    cuts_2.append(["AMass","<=",110.0])
                    cuts_2.append(["mll-mtt",">",0.0])
                    cuts_12.append(["AMass","<=",110.0])
                    cuts_12.append(["mll-mtt",">",0.0])
                tempmask_1 = cutOnArray(masterArray,cuts_1)
                #print tempmask_1[:1000]
                tempmask_2 = cutOnArray(masterArray,cuts_2)
                #print tempmask_2[:1000]
                tempmask_12 = cutOnArray(masterArray,cuts_12)
                #print tempmask_12[:1000]


                #FF_1
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                #fitmask_1 = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
                #fitmask_1 = fitmask_1.astype(int)
                ptarr_1 = masterArray["pt_3"]

                if const:
                    #ffweight_1 = ptFun(datadrivenPackage["fakerate1"],ptarr_1)
                    ffweight_1 = np.ones(len(tempmask_1))*(0.153/(1-0.153))
                    #ffweight_1 = ffweight_1/(1.0000000001 - ffweight_1)
                else:
                    #ffweight_1 = np.ones(len(ffweight_1))*(0.153/(1-0.153))
                    ffweight_1 = np.ones(len(tempmask_1))*(datadrivenPackage["fitrate1"].GetParameter(0)/(1-datadrivenPackage["fitrate1"].GetParameter(0)))



                #FF_2
                #fitmask_2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
                #fitmask_2 = fitmask_2.astype(int)
                ptarr_2 = masterArray["pt_4"]

                if const:
                    #ffweight_2 = ptFun(datadrivenPackage["fakerate2"],ptarr_2)
                    #ffweight_2 = ffweight_2/(1.0000000001 - ffweight_2)
                    ffweight_2 = np.ones(len(tempmask_2))*(0.153/(1-0.153))
                else:
                    ffweight_2 = np.ones(len(tempmask_2))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1-datadrivenPackage["fitrate2"].GetParameter(0)))


                ffweight = -1.0 * ffweight_1 * ffweight_2

                #replace 0s with constant fit value
                #ffweight_1 *= fitmask_1
                #ffweight_1[np.where(ffweight_1==0)] =  datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0))
                #ffweight_2 *= fitmask_2
                #ffweight_2[np.where(ffweight_2==0)] =  datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))

                #FF_12
                #replace 0s with constant fit value
                #ffweight[np.where(ffweight==0)] =  (datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0)))

                #ffweight[np.where(ffweight==0)] = 0.153**2

                #fitmask_1 *= fitmask_2
                #ffweight *= fitmask_1
                ffweight *= tempmask_12

                #ffweight_1[~tempmask_1] = 0.0
                ffweight_1 *= tempmask_1
                #print "ffweight_1 ",ffweight_1[:1000]
                #ffweight_2[~tempmask_2] = 0.0
                ffweight_2 *= tempmask_2
                #finalWeight = ffweight_1 + ffweight_2
                finalWeight = ffweight_1 + ffweight_2 + ffweight
                #print("pair 1-2: ",  np.any(np.all((tempmask_1,tempmask_2), axis=0)))
                #print("pair 1-12: ", np.any(np.all((tempmask_1,tempmask_12), axis=0)))
                #print("pair 12-2: ", np.any(np.all((tempmask_12,tempmask_2), axis=0)))

                #finalWeight -= ffweight


                #print "check on fake factor final weight ",finalWeight[:1000]


                masterArray["finalweight"] *= finalWeight
                #print "summed final weight ",np.sum(finalWeight)

                keepEvents = ~np.where(finalWeight==0.0)[0]

                skimArray={}
                for key in masterArray.keys():
                    skimArray[key] = masterArray[key][keepEvents]

                #print("after skim", len(skimArray["mll"]), processObj.file)
                if len(skimArray["mll"])==0:
                    continue

                for key in skimArray.keys():
                    if key not in plottedVars and key != "finalweight":
                        del skimArray[key]
                #return skimArray
                skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray




            if process not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
                EventWeights = getEventWeightDicitonary(args.year)

                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)


                mask = cutOnArray(masterArray,cuts)
                masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!

                skipEvents = np.where(mask==0)[0]
                skimArray={}
                if len(masterArray["mll"])==0:
                    continue

                weightDict = processObj.weights
                weightfinal = commonweight
                nickname = processObj.nickname
                for scalefactor in weightDict.keys():
                    if scalefactor == "kfactor":
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    elif scalefactor in ["PU"]:
                        masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                    elif scalefactor =="theoryXsec":
                        weightfinal =  weightfinal * float(weightDict[scalefactor])

                #print "finalweight after PU and kFactor ",masterArray["finalweight"][:100]

                if nickname in inclusive_samples:
                    jetweights = 5*[0]
                    if nickname.startswith("DY"):
                        jet0Xsec = 4673.65
                        if args.year=="2016":
                            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLLext1"][0]
                            for nj in range(1, 5):
                                njname = "DY%dJetsToLL"%nj
                                norm1 = jetWeightMultiplicity["DYJetsToLLext1"][nj]
                                norm2 = jetWeightMultiplicity[njname]
                                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                        elif args.year=="2017":
                            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLL_ext1"][0]
                            for nj in range(1, 5):
                                njname = "DY%dJetsToLL"%nj
                                norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"][nj]
                                norm2 = jetWeightMultiplicity[njname]
                                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                        elif args.year=="2018":
                            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLL"][0]
                            for nj in range(1, 5):
                                njname = "DY%dJetsToLL"%nj
                                norm1 = jetWeightMultiplicity["DYJetsToLL"][nj]
                                norm2 = jetWeightMultiplicity[njname]
                                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                    else:
                        jet0Xsec = 49033.2
                        jetweights[0] = jet0Xsec/jetWeightMultiplicity["WJetsToLNu"][0]
                        for nj in range(1, 5):
                            njname = "W%dJetsToLNu"%nj
                            norm1 = jetWeightMultiplicity["WJetsToLNu"][nj]
                            norm2 = jetWeightMultiplicity[njname]
                            jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)

                    for i_jet,weight in enumerate(jetweights):
                        if args.year==2016:
                            njetmask = masterArray["LHE_Njets"]==i_jet
                        else:
                            njetmask = masterArray["njets"]==i_jet
                        masterArray["finalweight"] [njetmask] *= weight
                        #print "events that pass ",i_jet," jets ",np.count_nonzero(masterArray["finalweight"] [njetmask])

                    #print "jet weight array ",jetweights
                    #print  " sample name ",nickname,"xsec ",HAA_processes[nickname].weights["xsec"]," events that pass ", np.count_nonzero(masterArray["finalweight"])


                elif not type(weightHistoDict[nickname])==list:
                #else: #Will the run separate from the NJet cases?
                    sumOfWeights = 0.0
                    for nic,sowhist in weightHistoDict.iteritems():
                        if nickname in nic:
                            sumOfWeights += sowhist.GetSumOfWeights()
                    if sumOfWeights != 0.0:
                        weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights
                        #print "xsec/SoW ",HAA_processes[nickname].weights["xsec"]/ sumOfWeights

                masterArray["finalweight"] *= weightfinal
                #print  " sample name ",nickname,"xsec ",HAA_processes[nickname].weights["xsec"]," SoW ",sumOfWeights," events that pass ", np.count_nonzero(masterArray["finalweight"])

                #multiply by scalar weight
                #print "finalweight before per event scaling ",masterArray["finalweight"][:100]


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
                            pass
                            #print "subtracting fakes "
                            #print weightMask[:100]

                        masterArray["finalweight"] *= weightMask
                #print "finalweight after per event scaling ",masterArray["finalweight"][:100]
                #print "summed final weight ",np.sum(masterArray["finalweight"])
                #print " scalar weightfinal ",weightfinal

                #print("before skim", len(masterArray["finalweight"]))
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
                #print("after skim", len(skimArray["mll"]), processObj.file)
                #for key in skimArray.keys():
                #    if key not in plottedVars and key != "finalweight":
                #        del skimArray[key]

                #catching events with strange weights... like just lumi ?? strange.
                #skimArray["finalweight"][np.where(skimArray["finalweight"]>2.0)]=0.0
                skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray
    #exit()
    return skimArrayPerCat



def slimskim(process,allcats,weightHistoDict,systematic):

    skimArrayPerSysCats={}
    #print "working on systematic ",systematic
    work_dict = {}
    #fin = uproot.open(process.file)
    with uproot.open(process.file) as fin:

        try:
            tree = fin[systematic]
        except:
            return skimArrayPerSysCats

        if systematic!="Events":
            syst_names = set(fin[systematic].keys())
            nom_names = set(fin["Events"].keys()) - syst_names
            work_dict.update(fin[systematic].arrays(list(syst_names)))
            work_dict.update(fin["Events"].arrays(list(nom_names)))
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))
        else:
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal"))

    del work_dict
    del tree
    return skimArrayPerSysCats




def slimskimoutput(process,allcats,weightHistoDict,systematic,massoutputdir,datadrivenPackage,args):

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
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic,args))
        else:
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal",args))


        createSlimOutput(skimArrayPerSysCats,massoutputdir)
    del skimArrayPerSysCats
    del work_dict
    del tree
    return

def createSlimOutput(skimArrayPerSysCats,outputdir):
   try:
       os.mkdir(outputdir)
   except:
       print "dir exists " , outputdir

   for key,dictionary in skimArrayPerSysCats.iteritems():
      dataTypes =[[],[]]
      for branch in dictionary.keys():
          dataTypes[0].append(branch)
          dataTypes[1].append(dictionary[branch].dtype)

      data = np.zeros(len(dictionary[branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})

      filename = key.replace(":","_")
      fileout = ROOT.TFile.Open(outputdir+"/"+filename,"recreate")
      fileout.cd()
      for branch in data.dtype.names:
          if len(dictionary[branch].shape) == 1:   # flat array important for jagged arrays of input data
              data[branch] = dictionary[branch]
          else:
              data[branch] = dictionary[branch][:,0]

      #skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process] = skimArray
      #name = key.split(":")
      #root_numpy.array2tree(data,name=key.split(":")[0]+"_"+key.split(":")[-1])
      root_numpy.array2tree(data,name=(filename+"_t"))
      fileout.Write()
      fileout.Close()
      del data
      del fileout
   return




def combineRootFiles(systematics, allcats, processes,
                     finalDistributions, rootfiledir,
                     channel, outputstring):
   import os
   import glob
   import uproot
   rootFiles = {}
   #print "the final distributions ",finalDistributions
   #cat = "mmmt_inclusive"
   finalSkims={}
   import shutil
   #list comprehension example:
   systematics[systematics.index("Events")]="Nominal"
   finalSkims = {s:{c: dict() for c in allcats.keys()} for s in systematics}


   rootFiles={sys:glob.glob(rootfiledir+"/"+sys+"_*") for sys in systematics}


   mainTList = {}
   mainOutputTree={}
   for sys, globfiles in rootFiles.iteritems():
      for globfile in globfiles:
         #print "working on glob ",globfile
         #process = globfile.split("_")[-1] # problem with data_obs or fake_W
         #process = globfile.split(sys)[0].split(channel+"_inclusive_")[1]
         #added loop for better processing?
         for cat,catObj in allcats.iteritems():
             #if cat==args.channel+"_inclusive":
            for nickname, processObj in processes.iteritems():
                for process in processObj.cuts.keys():
                    #if nickname in str(globfile) and nickname in str(globfile):
                    if str(globfile).split("/")[1]==sys+"_"+cat+"_"+nickname+"_"+process:
                        #print "found match ",sys+"_"+cat+"_"+nickname+"_"+process
                        with uproot.open(globfile) as fin:
                            tree = fin[sys+"_"+cat+"_"+nickname+"_"+process + "_t"]
                            mainArrays = tree.arrays()
                            #print "tree ",sys+"_"+process," entries ",len(mainArrays["mll"])
                            for catDist, final in finalDistributions.iteritems():
                                for processOut in final:
                                    if (processOut==process) and (catDist not in finalSkims[sys][cat]):
                                        #print "first output for process ",process," finalDist cat ",catDist
                                        finalSkims[sys][cat][catDist] = mainArrays
                                        continue
                                    elif (processOut==process) and (catDist in finalSkims[sys][cat]):
                                        #print "adding to finalskims ", catDist,"  for process ",process," finalDist cat ",catDist
                                        for branch in finalSkims[sys][cat][catDist].keys():
                                            finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],mainArrays[branch]))
                                    else:
                                        continue

   #print "final skims example ", finalSkims["Nominal"]["mmtt_inclusive"]["irBkg"]
   skimFile = ROOT.TFile("skimmed_"+outputstring+".root","recreate")
   skimFile.cd()
   for cat, catObj in allcats.iteritems():
       #localtest="allsys"
       #skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
       skimFile.cd()
       skimFile.mkdir(cat)
       skimFile.cd(cat)
       for sys in finalSkims.keys():
           dataTypes =[[],[]]
           #print finalSkims[sys][cat].values()
           print "combining ",sys, " cat ",cat
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


if __name__ == "__main__":
    import datetime
    begin_time = datetime.datetime.now()
    from multiprocessing import *
    import multiprocessing as mp
    import logging
    import os
    import shutil


    #trial of global settings
    #import AnalysisSetup as AS
    from utils.Weights import jet_inclusive_samples


    import argparse


    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-yr",  "--year", default="2016",  help="Year")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
    parser.add_argument("-ch",  "--channel", default="mmmt",  help="Please list the channel for fake factor histograms")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/",  help="Input files")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dbg",  "--debug", default=False,action='store_true',  help="Disable the parallel processing mode")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ex",  "--extract", default=False,action='store_true',  help="Additional Cuts for Extraction")
    parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
    parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
    parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
    parser.add_argument("-s",  "--skim", default=False,action='store_true',  help="skim input files to make more TTrees")
    parser.add_argument("-mt",  "--mt", default=False,action='store_true',  help="Use Multithreading")
    parser.add_argument("-pt",  "--maxprint", default=False,action='store_true',  help="Print Info on cats and processes")
    parser.add_argument("-co",  "--constant", default=False,action='store_true',  help="True to use constant value instead of fitting.")

    args = parser.parse_args()
    allcats={}
    HAA_processes={}
    finalDistributions={}
    filelist = {}
    weightHistoDict={}
    jetWeightMultiplicity={}
    EventWeights={}
    datadrivenPackage={}

    const = args.constant

    allcats,HAA_processes,finalDistributions,\
    weightHistoDict,jetWeightMultiplicity,datadrivenPackage = initialize(args)

    print datadrivenPackage

    #print(allcats)
    info('main line')
    #print(f)
    nums=[]
    skims=[]
    #skims={}
    payloadsdict={}
    payloads=[]

    systematics =[ "Events" ] #,"scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
               #    "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
               #    "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
               #    "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]
    #systematics =[ "Events","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]
    #if args.datameasureZH or args.datadrivenZH:
    #    systematics =[ "Events"]
    #systematics =[ "Events"]

    inclusive_samples = jet_inclusive_samples[args.year]


    for nickname, process in HAA_processes.items():
        #if args.datameasureZH and process=="FF": continue
        #if nickname in ["data","ZZTo4L"]:
        #if nickname in ["data","WJetsToLNu","DYJetsToLLext1"]:
        # if nickname not in inclusive_samples:
        #     continue
        for sys in systematics:
            payloads.append(
                (process,allcats,
                weightHistoDict,sys,
                "massOutputDir_"+args.outname,
                 datadrivenPackage,args))

    #print(" PAYLOADS   ",payloads)

    if not args.debug:
        m = mp.Manager()
        #logger_q = m.Queue()
        #parallelable_data = [(1, logger_q), (2, logger_q)]
        #pool  = mp.Pool(12)
        #pool  = mp.Pool(4)
        pool  = mp.Pool(4)

        pool.map(slimskimstar,payloads)#this works for root output!
        pool.close()
        pool.join()
        #while not logger_q.empty():
        #    print logger_q.get()
    else:
        for payload in payloads:
            #print "working on payload ", payload
            slimskimstar(payload)


    print("root file generation computation time")
    print(datetime.datetime.now() - begin_time)


    #print("combining the output")
    #if createOutput(skims,finalDistributions) : print "successful output"
    #else: print "ouch ..."
    #if createOutputSystematics(skims,finalDistributions) : print "successful output"
    #else: print "ouch ..."

    if combineRootFiles(systematics, allcats, HAA_processes,
                     finalDistributions, "massOutputDir_"+args.outname,
                     args.channel, args.ffout):
        print "combination successful"

        #shutil.rmtree("massOutputDir_"+args.outname)

    if args.datadrivenZH or args.makeFakeHistos:
        if not const:
            print("closing fakemeasurefile.")
            datadrivenPackage["fakemeasurefile"].Close()

    print("computation time")
    print(datetime.datetime.now() - begin_time)
