#########################
#Author: Sam Higginbotham
'''

* File Name : Parametrization.py

* Purpose : Loads the variables that will be plotted against one another for further analysis

* Creation Date : 04-02-2020

* Last Modified : Development ongoing - plot the same variable in multiple ways ... "extravariabels"

'''
#########################

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
