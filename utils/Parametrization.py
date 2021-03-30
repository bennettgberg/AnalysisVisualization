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
