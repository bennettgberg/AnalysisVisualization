#########################
#Author: Sam Higginbotham
'''

* File Name : Parametrization.py

* Purpose : Loads the variables that will be plotted against one another for further analysis

* Creation Date : 04-02-2020

* Last Modified : Development ongoing - plot the same variable in multiple ways ... "extravariabels"

'''
#########################
import numpy as np
import ROOT

class Category():
    def __init__(self):
        self.name = []
        self.variables = []
        self.vars={}
        self.newvariables = {}
        self.newvariablesbins = []
        self.cuts = {}
        self.newvarcuts = {}
        self.binning = []
        self.extraplots = {}
        self.systematics = {}

class Process():
    def __init__(self):
        #Nickname corresponds to CSV input file ... wil be the unique identifier
        self.nickname = ""
        self.plotname = ""
        #what are some important weights or scale factors NOT COMMON to all processes????
        #xsec, nevents,
        self.weights = {}
        #eventscale weights... a selection string and what weight to apply
        self.eventWeights = {}
        #cuts per process... overrides the plot name ... and we loop over collection of these in main
        self.cuts = {}
        self.file = ""
        self.classification = ""
        self.color = []


# class cutFlow():
#     def __init__(self):
#         #add containors for histograms
#         from ROOT import TH1, TFile,
#         self.hCutFlow_mmet = ROOT.TH1D()
#         self.hCutFlowWeighted_mmet = ROOT.TH1D()
#         self.hCutFlow_mmmt = ROOT.TH1D()
#         self.hCutFlowWeighted_mmmt = ROOT.TH1D()
#         self.hCutFlow_mmtt = ROOT.TH1D()
#         self.hCutFlowWeighted_mmtt = ROOT.TH1D()
#         self.hCutFlow_mmem = ROOT.TH1D()
#         self.hCutFlowWeighted_mmem = ROOT.TH1D()
#
#     def gatherHistos(self,infile):
#
#
#         return
#         #hWeights



class fakefactor() :
    #class to read the Standard Model Higgs Tau Tau, Tau fake rates
    def __init__(self):
        self.w1  = 0
        self.w2  = 0
        self.w12 = 0
        self.year = 0
        self.fakefile_mu  = ROOT.TFile()
        self.fakefile_tau = ROOT.TFile()
        self.fakefile_ele = ROOT.TFile()
        self.fakehisto_mu =  {}
        self.fakehisto_ele = {}
        self.fakehisto_tau = {}

    def Print(self):
        print "the gathered fake factor histograms"
        print self.fakehisto_mu
        print self.fakehisto_ele
        print self.fakehisto_tau

    def loadHistograms(self,filepath,year):
        print "reading file ",filepath+"FitHistograms_muFR_"+year+".root"
        self.fakefile_mu = ROOT.TFile.Open(filepath+"FitHistograms_muFR_"+year+".root","READ")
        print self.fakefile_mu.GetListOfKeys()
        self.fakefile_tau = ROOT.TFile.Open(filepath+"FitHistograms_tauFR_"+year+".root","READ")
        self.fakefile_ele = ROOT.TFile.Open(filepath+"FitHistograms_eleFR_"+year+".root","READ")

        for key in enumerate(self.fakefile_mu.GetListOfKeys()):
            #obj = key.ReadObj()
            #self.fakehisto_mu[obj.GetName()] = obj
            obj = key[1].ReadObj()
            self.fakehisto_mu[0] = obj
        for key in enumerate(self.fakefile_ele.GetListOfKeys()):
            #obj = key.ReadObj()
            #self.fakehisto_ele[obj.GetName()] = obj
            obj = key[1].ReadObj()
            self.fakehisto_ele[0] = obj
        for key in enumerate(self.fakefile_tau.GetListOfKeys()):
            #obj = key.ReadObj()
            #self.fakehisto_tau[obj.GetName()] = obj
            obj = key[1].ReadObj()
            self.fakehisto_tau[key[0]] = obj

        return 1

        #map the histo in to the channel and cuts ....
    def getHistoMask(self,channel,processObj,masterArray):
        # nomenclature of histograms
        # e file: FitHistograms_eleFR_2016.root
        # efr_numerator_efr_denominator
        # mu file: FitHistograms_muFR_2016.root
        # mufr_numerator_mufr_denominator
        # tau file: FitHistograms_tauFR_2016.root
        # hpt_dmall_Te_VLmu_1jet_deepveryloose_hpt_dmall_Te_VLmu_1jet_deepveryveryveryloose
        if channel=="mmem":
            index1 = np.full(len(masterArray["evt"]),0)
            index2 = np.full(len(masterArray["evt"]),0)
        '''
        2  hpt_dm0_deepveryveryloose_hpt_dm0_deepveryveryveryloose
        4  hpt_dm0_deepveryloose_hpt_dm0_deepveryveryveryloose
        8  hpt_dm0_deeploose_hpt_dm0_deepveryveryveryloose
        16 hpt_dm0_deepmedium_hpt_dm0_deepveryveryveryloose
        32 hpt_dm0_deeptight_hpt_dm0_deepveryveryveryloose
        64 hpt_dm0_deepverytight_hpt_dm0_deepveryveryveryloose
        if cat=='tt' : gName = "hpt_dm{0:d}_VLe_VLmu_deep{1:s}_hpt_dm{0:d}_VLe_VLmu_deepveryveryveryloose".format(DM,WPdict[WP])
        if cat=='mt' : gName = "hpt_dm{0:d}_VLe_Tmu_deep{1:s}_hpt_dm{0:d}_VLe_Tmu_deepveryveryveryloose".format(DM,WPdict[WP])
        if cat=='et' : gName = "hpt_dm{0:d}_Te_VLmu_deep{1:s}_hpt_dm{0:d}_Te_VLmu_deepveryveryveryloose".format(DM,WPdict[WP])
        (168, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeptight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a7ab40>)
        (169, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeptight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a7ad40>)
        (170, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeptight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a7af40>)
        (171, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepverytight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a7b140>)
        (172, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepverytight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4d7e0>)
        (173, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepverytight_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4d9f0>)
        (174, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepmedium_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4dc00>)
        (175, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepmedium_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4de80>)
        (176, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepmedium_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e090>)
        (177, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeploose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e2a0>)
        (178, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeploose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e450>)
        (179, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deeploose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e650>)
        (180, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e850>)
        (181, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4ea60>)
        (182, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4ec70>)
        (183, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a185c0>)
        (184, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a187a0>)
        (185, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepveryveryloose_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a189b0>)
        '''
        # masterArray["mll"][(masterArray["mll"]>0.0)*(masterArray["cat"]==6)]
        if channel=="mmmt":
            index1 = np.full(len(masterArray["evt"]),0)
            index2_0 = ( # medium
                     (masterArray["decayMode_4"]==0)*(masterArray["cat"]==6)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*176
                      #(176, <ROOT.TKey object ("hpt_dm0_VLe_Tmu_deepmedium_hpt_dm0_VLe_Tmu_deepveryveryveryloose") at 0x4a4e090>)
            index2_1 = (
                     (masterArray["decayMode_4"]==1)*(masterArray["cat"]==6)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*194
                      #(194, <ROOT.TKey object ("hpt_dm1_VLe_Tmu_deepmedium_hpt_dm1_VLe_Tmu_deepveryveryveryloose") at 0x4a19c80>)
            index2_2 = (
                     (masterArray["decayMode_4"]==10)*(masterArray["cat"]==6)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*212
                      #(212, <ROOT.TKey object ("hpt_dm10_VLe_Tmu_deepmedium_hpt_dm10_VLe_Tmu_deepveryveryveryloose") at 0x4a7c020>)
            index2_3 = (
                     (masterArray["decayMode_4"]==11)*(masterArray["cat"]==6)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*230
                      #(230, <ROOT.TKey object ("hpt_dm11_VLe_Tmu_deepmedium_hpt_dm11_VLe_Tmu_deepveryveryveryloose") at 0x4a4b980>)
            index2 = index2_0+index2_1+index2_2+index2_3

        if channel=="mmtt":
            index1_0 = ( # medium
                     (masterArray["decayMode_3"]==0)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_3"]>=15) \
                      ).astype(int)*320
                      #(320, <ROOT.TKey object ("hpt_dm0_VLe_VLmu_deepmedium_hpt_dm0_VLe_VLmu_deepveryveryveryloose") at 0x4a5da10>)
            index1_1 = (
                     (masterArray["decayMode_3"]==1)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_3"]>=15) \
                      ).astype(int)*338
                      #(338, <ROOT.TKey object ("hpt_dm1_VLe_VLmu_deepmedium_hpt_dm1_VLe_VLmu_deepveryveryveryloose") at 0x4a7fe50>)
            index1_2 = (
                     (masterArray["decayMode_3"]==10)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_3"]>=15) \
                      ).astype(int)*356
                      #(356, <ROOT.TKey object ("hpt_dm10_VLe_VLmu_deepmedium_hpt_dm10_VLe_VLmu_deepveryveryveryloose") at 0x4a82450>)
            index1_3 = (
                     (masterArray["decayMode_3"]==11)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_3"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_3"]>=15) \
                      ).astype(int)*374
                      #(374, <ROOT.TKey object ("hpt_dm11_VLe_VLmu_deepmedium_hpt_dm11_VLe_VLmu_deepveryveryveryloose") at 0x49e05d0>)
            index2_0 = ( # medium
                     (masterArray["decayMode_4"]==0)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*320
                      #(320, <ROOT.TKey object ("hpt_dm0_VLe_VLmu_deepmedium_hpt_dm0_VLe_VLmu_deepveryveryveryloose") at 0x4a5da10>)
            index2_1 = (
                     (masterArray["decayMode_4"]==1)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*338
                      #(338, <ROOT.TKey object ("hpt_dm1_VLe_VLmu_deepmedium_hpt_dm1_VLe_VLmu_deepveryveryveryloose") at 0x4a7fe50>)
            index2_2 = (
                     (masterArray["decayMode_4"]==10)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*356
                      #(356, <ROOT.TKey object ("hpt_dm10_VLe_VLmu_deepmedium_hpt_dm10_VLe_VLmu_deepveryveryveryloose") at 0x4a82450>)
            index2_3 = (
                     (masterArray["decayMode_4"]==11)*(masterArray["cat"]==7)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*374
                      #(374, <ROOT.TKey object ("hpt_dm11_VLe_VLmu_deepmedium_hpt_dm11_VLe_VLmu_deepveryveryveryloose") at 0x49e05d0>)

            index1 = index1_0+index1_1+index1_2+index1_3
            index2 = index2_0+index2_1+index2_2+index2_3

        if channel=="mmet":
            #leg1histo = self.fakehisto_ele["efr_numerator_efr_denominator"] #why???
            index1 = np.full(len(masterArray["evt"]),0)
            index2_0 = ( # medium
                     (masterArray["decayMode_4"]==0)*(masterArray["cat"]==5)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*248
                      #(248, <ROOT.TKey object ("hpt_dm0_Te_VLmu_deepmedium_hpt_dm0_Te_VLmu_deepveryveryveryloose") at 0x4a1b020>)
            index2_1 = (
                     (masterArray["decayMode_4"]==1)*(masterArray["cat"]==5)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*266
                      #(266, <ROOT.TKey object ("hpt_dm1_Te_VLmu_deepmedium_hpt_dm1_Te_VLmu_deepveryveryveryloose") at 0x4a08e50>)
            index2_2 = (
                     (masterArray["decayMode_4"]==10)*(masterArray["cat"]==5)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*284
                      #(284, <ROOT.TKey object ("hpt_dm10_Te_VLmu_deepmedium_hpt_dm10_Te_VLmu_deepveryveryveryloose") at 0x4a5afe0>)
            index2_3 = (
                     (masterArray["decayMode_4"]==11)*(masterArray["cat"]==5)*\
                     (masterArray["idDeepTau2017v2p1VSjet_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSmu_4"]>=15)* \
                     (masterArray["idDeepTau2017v2p1VSe_4"]>=15) \
                      ).astype(int)*302
                      #(302, <ROOT.TKey object ("hpt_dm11_Te_VLmu_deepmedium_hpt_dm11_Te_VLmu_deepveryveryveryloose") at 0x4a774f0>)

            index2 = index2_0+index2_1+index2_2+index2_3

        index1[np.where(index1==0)]==-1
        index2[np.where(index2==0)]==-1
        return index1,index2

    def getFakeWeight(self,channel,pt_1,pt_2,index1,index2):

        if channel=="mmem":
            ffweight_1 = np.full(len(pt_1),1.0)
            ffweight_1 = np.vectorize(self.fakehisto_ele[0].Eval)(pt_1)
            ffweight_2 = np.full(len(pt_2),1.0)
            ffweight_2 = np.vectorize(self.fakehisto_mu[0].Eval)(pt_2)
        if channel=="mmmt":
            ffweight_1 = np.full(len(pt_1),1.0)
            ffweight_1 = np.vectorize(self.fakehisto_mu[0].Eval)(pt_1)
            ffweight_2 = np.full(len(pt_2),1.0)
            #ffweight_2 = np.array([])
            ffweight_2 = self.gatherFF(index2,pt_2,ffweight_2)
        if channel=="mmtt":
            ffweight_1 = np.full(len(pt_1),1.0)
            #ffweight_1 = np.array([])
            ffweight_1 = self.gatherFF(index1,pt_1,ffweight_1)
            ffweight_2 = np.full(len(pt_2),1.0)
            #ffweight_2 = np.array([])
            ffweight_2 = self.gatherFF(index2,pt_2,ffweight_2)
        if channel=="mmet":
            ffweight_1 = np.full(len(pt_1),1.0)
            ffweight_1 = np.vectorize(self.fakehisto_ele[0].Eval)(pt_1)
            ffweight_2 = np.full(len(pt_2),1.0)
            #ffweight_2 = np.array([])
            ffweight_2 = self.gatherFF(index2,pt_2,ffweight_2)

        #print "index2 ",index2[:100]
        print "fake weight 2 in class ",ffweight_2[:100]
        return ffweight_1,ffweight_2

    def gatherFF(self,indicies,pt,ffweight):
        #np.append(ffweight,[fakehisto_tau[id].Eval(masterArray["pt_3"][i])],axis=0)
        for num,index in enumerate(indicies):
            if index==-1:
                np.append(ffweight,0,axis=0)
            else:
                #print "pt ",pt[num]," measured fake rate  ",self.fakehisto_tau[index].Eval(pt[num])
                #np.append(ffweight,[self.fakehisto_tau[index].Eval(pt[num])],axis=0)
                ffweight[num]=self.fakehisto_tau[index].Eval(pt[num])


        return ffweight


from ROOT import TFile, TTree
import sys

class dupeDetector() :

    def __init__(self):
        self.nCalls = 0
        self.runEventList = {}

    def checkEvent(self,entry) :
        self.nCalls += 1
        runEvent = "{0:d}:{1:d}:{2:d}".format(entry.run,entry.evt,entry.cat)
        try :
            return self.runEventList[runEvent]
        except KeyError :
            self.runEventList[runEvent] = True
            return False

    def printSummary(self) :
        print("Duplicate Event Summary: Calls={0:d} Unique Events={1:d}".format(self.nCalls,len(self.runEventList.keys())))
        return
