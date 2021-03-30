#########################
#Author: Sam Higginbotham
'''

* File Name :functions.py

* Purpose : The Categories in the analysis... typically a set of extra cuts that define a phase space region of interest and stores the variables of interest in the fit and visualization

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
import numpy as np
from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
#load these functions into a dictionary so we can call on them in the event loop

functs={}
minus = lambda evt,x: getattr(evt,x[0],None) - getattr(evt,x[1],None)
functs["minus"]=minus
plus = lambda evt,x: getattr(evt,x[0],None) + getattr(evt,x[1],None)
functs["plus"]=plus
minusArr = lambda x,y: x-y
functs["minusArr"]=minusArr
plusArr = lambda x,y: x+y
functs["plusArr"]=plusArr
def multi(evt,x):
    product = 1.0
    for i,j in enumerate(x):
        product = product * getattr(evt,x[i],None)
    #if product>0:
    #    print evt.evt,x[0],x[1],"  multi   ",product
    return product
functs["multi"] = multi

def multiArr(x,y):
    return x*y
functs["multiArr"] = multiArr

def TLorentz(evt,x):
    vector = ROOT.TLorentzVector()
    if type(x[3])==float:
	    vector.SetPtEtaPhiM(getattr(evt,x[0],None),getattr(evt,x[1],None),getattr(evt,x[2],None),x[3])  #set the mass of the TLorentz Vector directly
    else:
	    vector.SetPtEtaPhiM(getattr(evt,x[0],None),getattr(evt,x[1],None),getattr(evt,x[2],None),getattr(evt,x[3],None))
    vectorList.append(vector[0])
    vectorList.append(vector[1])
    vectorList.append(vector[2])
    vectorList.append(vector[3])
    return vectorList
functs["TLorentz"] = TLorentz

def TLorentzMT(vector):
    return vector.Mt()
functs["TLorentzMT"] = TLorentzMT

def addTLorentz(allvectors):
    sumvector = ROOT.TLorentzVector()
    for vector in allvectors:
        sumvector = sumvector + vector
    return sumvector
functs["addTLorentz"] = addTLorentz

def invariantMass(pt1,pt2,eta1,eta2,phi1,phi2):
    return np.sqrt(2*pt1*pt2*(np.cosh(eta1-eta2)-np.cos(phi1-phi2)))
functs["invariantMass"] = invariantMass
