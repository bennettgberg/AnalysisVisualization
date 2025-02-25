#########################
#Author: Sam Higginbotham
'''

* File Name :Categories.py 

* Purpose : The Categories in the analysis... typically a set of extra cuts that define a phase space region of interest and stores the variables of interest in the fit and visualization

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
from Parametrization import Category
from functions import *
import copy

#fakefactor working points and cuts 
#15.0 is medium
#8.0 is loose
ffworkingpoint = 15.
#ffworkingpoint = 8.0
#ffworkingpoint = 2.0
ffbasecuts = [["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],["eta_4","absl",2.3],[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],["iso_1","<=",0.2],["iso_2","<=",0.2],["cat","==",6]]

#With HAA front end ... MUST select category -> channel ... 
#number = { 'eeet':1, 'eemt':2, 'eett':3, 'eeem':4, 'mmet':5, 'mmmt':6, 'mmtt':7, 'mmem':8, 'et':9, 'mt':10, 'tt':11 }
HAA_Inc_mmmt = Category()
HAA_Inc_mmmt.name = "mmmt_inclusive"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_1","q_2"],"mult","<",0],
#[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult",">",0],
#[["OR"],["q_3",">",0],["q_4",">",0][["EQT"],["q_1","q_2"],"mult",">",0],
["iso_1","<=",0.2],["iso_2","<=",0.2],
["nbtag","<",1.0]]   #

HAA_Inc_mmmt.cuts["categoryCuts"]= [
#["cat","==",6],["AMass",">=",120.],
["cat","==",6],[["EQT"],["q_3","q_4"],"mult","<",0],
["iso_3","<=",0.15],["mediumId_3",">=",1],
[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
#["idDeepTau2017v2p1VSjet_4",">=",2.]]
["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]] #this is cutting most events! 

#Trigger bit mapping... bits = [e.HLT_Ele27_eta2p1_WPTight_Gsf, e.HLT_Ele25_eta2p1_WPTight_Gsf, e.HLT_IsoMu24, e.HLT_IsoTkMu24, e.HLT_IsoMu27,e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ]
#2 -> 4 IsoMu24 ,  3 -> 8 IsoTkMu24  , 4 -> 16 IsoMu27  , 5 -> 32  e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ 
#HAA_Inc_mmmt.cuts["trigger"]= [[["triggerWord","band",4],"||",["triggerWord","band",8],"||",["triggerWord","band",16],"||",["triggerWord","band",32]]]   #
#if the statements are nested like this ... they are ors 
HAA_Inc_mmmt.cuts["trigger"]=[
#[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
[["OR"],["muonTriggerWord","band",4],["muonTriggerWord","band",8],["muonTriggerWord","band",16],["muonTriggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

#new variables that are defined on the fly. These are not in the Tree, but are functions of variables in the tree
#see function definitions above
HAA_Inc_mmmt.newvariables["mll-mtt"] = ["minus",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["mll","m_vis"],"[GeV]","mll - m_{vis}"]
HAA_Inc_mmmt.newvariables["charge_12"] = ["multi",[-2.0,-0.75,0.75,2.0],["q_1","q_2"],"","charge_1 #cross charge_2"]
HAA_Inc_mmmt.newvariables["charge_34"] = ["multi",[-2.0,-0.75,0.75,2.0],["q_3","q_4"],"","charge_3 #cross charge_4"]
#HAA_Inc_mmmt.newvariables["mll-mtt"] = tauSFTool.getSFvsPT(e.pt_3,e.gen_match_3)["mll","mtt"]

HAA_Inc_mmmt.varis = {
        #handle : ["root variable",[binning]or[[binning]],units,label]
        "m_vis":["m_vis",[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],"[GeV]","M_{vis}"],
        "AMass_blinded":["AMass",[100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot} Blind"],
        "AMass":["AMass",[0.0,20.0,40.0,60.0,80.0,100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot}"],
        "mll":["mll",[0.0,10.,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.,100.,110.,120.,130.,140.],"[GeV]","M_{#mu#mu}"],
        #"mll_m15":["mll",[12.0,13.0,13.125,13.25,13.375,13.5,13.625,13.75,13.875,14.0,14.125,14.25,14.375,14.5,14.625,14.75,14.875,15.0,15.125,15.25,15.375,15.5,15.625,15.75,15.875,16.0,16.125,16.25,16.375,16.5,16.625,16.75,16.875,17.0 ,18.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m20":["mll",[17.0,18.0,18.125,18.25,18.375,18.5,18.625,18.75,18.875,19.0,19.125,19.25,19.375,19.5,19.625,19.75,19.875,20.0,20.125,20.25,20.375,20.5,20.625,20.75,20.875,21.0,21.125,21.25,21.375,21.5,21.625,21.75,21.875,22.0 ,23.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m25":["mll",[22.0,23.0,23.125,23.25,23.375,23.5,23.625,23.75,23.875,24.0,24.125,24.25,24.375,24.5,24.625,24.75,24.875,25.0,25.125,25.25,25.375,25.5,25.625,25.75,25.875,26.0,26.125,26.25,26.375,26.5,26.625,26.75,26.875,27.0 ,28.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m30":["mll",[27.0,28.0,28.125,28.25,28.375,28.5,28.625,28.75,28.875,29.0,29.125,29.25,29.375,29.5,29.625,29.75,29.875,30.0,30.125,30.25,30.375,30.5,30.625,30.75,30.875,31.0,31.125,31.25,31.375,31.5,31.625,31.75,31.875,32.0,33.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m35":["mll",[32.0,33.0,33.125,33.25,33.375,33.5,33.625,33.75,33.875,34.0,34.125,34.25,34.375,34.5,34.625,34.75,34.875,35.0,35.125,35.25,35.375,35.5,35.625,35.75,35.875,36.0,36.125,36.25,36.375,36.5,36.625,36.75,36.875,37.0,38.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m40":["mll",[37.0,38.0,38.125,38.25,38.375,38.5,38.625,38.75,38.875,39.0,39.125,39.25,39.375,39.5,39.625,39.75,39.875,40.0,40.125,40.25,40.375,40.5,40.625,40.75,40.875,41.0,41.125,41.25,41.375,41.5,41.625,41.75,41.875,42.0 ,43.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m45":["mll",[42.0,43.0,43.125,43.25,43.375,43.5,43.625,43.75,43.875,44.0,44.125,44.25,44.375,44.5,44.625,44.75,44.875,45.0,45.125,45.25,45.375,45.5,45.625,45.75,45.875,46.0,46.125,46.25,46.375,46.5,46.625,46.75,46.875,47.0 ,48.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m50":["mll",[47.0,48.0,48.125,48.25,48.375,48.5,48.625,48.75,48.875,49.0,49.125,49.25,49.375,49.5,49.625,49.75,49.875,50.0,50.125,50.25,50.375,50.5,50.625,50.75,50.875,51.0,51.125,51.25,51.375,51.5,51.625,51.75,51.875,52.0 ,53.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m55":["mll",[52.0,53.0,53.125,53.25,53.375,53.5,53.625,53.75,53.875,54.0,54.125,54.25,54.375,54.5,54.625,54.75,54.875,55.0,55.125,55.25,55.375,55.5,55.625,55.75,55.875,56.0,56.125,56.25,56.375,56.5,56.625,56.75,56.875,57.0 ,58.0],"[GeV]","M_{#mu#mu}"],
        #"mll_m60":["mll",[57.0,58.0,58.125,58.25,58.375,58.5,58.625,58.75,58.875,59.0,59.125,59.25,59.375,59.5,59.625,59.75,59.875,60.0,60.125,60.25,60.375,60.5,60.625,60.75,60.875,61.0,61.125,61.25,61.375,61.5,61.625,61.75,61.875,62.0 ,63.0],"[GeV]","M_{#mu#mu}"],
        "mll_m15":["mll",[13.0,13.125,13.25,13.375,13.5,13.625,13.75,13.875,14.0,14.125,14.25,14.375,14.5,14.625,14.75,14.875,15.0,15.125,15.25,15.375,15.5,15.625,15.75,15.875,16.0,16.125,16.25,16.375,16.5,16.625,16.75,16.875,17.0],"[GeV]","M_{#mu#mu}"],
        "mll_m20":["mll",[18.0,18.125,18.25,18.375,18.5,18.625,18.75,18.875,19.0,19.125,19.25,19.375,19.5,19.625,19.75,19.875,20.0,20.125,20.25,20.375,20.5,20.625,20.75,20.875,21.0,21.125,21.25,21.375,21.5,21.625,21.75,21.875,22.0],"[GeV]","M_{#mu#mu}"],
        "mll_m25":["mll",[23.0,23.125,23.25,23.375,23.5,23.625,23.75,23.875,24.0,24.125,24.25,24.375,24.5,24.625,24.75,24.875,25.0,25.125,25.25,25.375,25.5,25.625,25.75,25.875,26.0,26.125,26.25,26.375,26.5,26.625,26.75,26.875,27.0],"[GeV]","M_{#mu#mu}"],
        "mll_m30":["mll",[28.0,28.125,28.25,28.375,28.5,28.625,28.75,28.875,29.0,29.125,29.25,29.375,29.5,29.625,29.75,29.875,30.0,30.125,30.25,30.375,30.5,30.625,30.75,30.875,31.0,31.125,31.25,31.375,31.5,31.625,31.75,31.875,32.0],"[GeV]","M_{#mu#mu}"],
        "mll_m35":["mll",[33.0,33.125,33.25,33.375,33.5,33.625,33.75,33.875,34.0,34.125,34.25,34.375,34.5,34.625,34.75,34.875,35.0,35.125,35.25,35.375,35.5,35.625,35.75,35.875,36.0,36.125,36.25,36.375,36.5,36.625,36.75,36.875,37.0],"[GeV]","M_{#mu#mu}"],
        "mll_m40":["mll",[38.0,38.125,38.25,38.375,38.5,38.625,38.75,38.875,39.0,39.125,39.25,39.375,39.5,39.625,39.75,39.875,40.0,40.125,40.25,40.375,40.5,40.625,40.75,40.875,41.0,41.125,41.25,41.375,41.5,41.625,41.75,41.875,42.0],"[GeV]","M_{#mu#mu}"],
        "mll_m45":["mll",[43.125,43.25,43.375,43.5,43.625,43.75,43.875,44.0,44.125,44.25,44.375,44.5,44.625,44.75,44.875,45.0,45.125,45.25,45.375,45.5,45.625,45.75,45.875,46.0,46.125,46.25,46.375,46.5,46.625,46.75,46.875],"[GeV]","M_{#mu#mu}"],
        "mll_m50":["mll",[48.0,48.125,48.25,48.375,48.5,48.625,48.75,48.875,49.0,49.125,49.25,49.375,49.5,49.625,49.75,49.875,50.0,50.125,50.25,50.375,50.5,50.625,50.75,50.875,51.0,51.125,51.25,51.375,51.5,51.625,51.75,51.875,52.0],"[GeV]","M_{#mu#mu}"],
        "mll_m55":["mll",[53.0,53.125,53.25,53.375,53.5,53.625,53.75,53.875,54.0,54.125,54.25,54.375,54.5,54.625,54.75,54.875,55.0,55.125,55.25,55.375,55.5,55.625,55.75,55.875,56.0,56.125,56.25,56.375,56.5,56.625,56.75,56.875,57.0],"[GeV]","M_{#mu#mu}"],
        "mll_m60":["mll",[58.0,58.125,58.25,58.375,58.5,58.625,58.75,58.875,59.0,59.125,59.25,59.375,59.5,59.625,59.75,59.875,60.0,60.125,60.25,60.375,60.5,60.625,60.75,60.875,61.0,61.125,61.25,61.375,61.5,61.625,61.75,61.875,62.0],"[GeV]","M_{#mu#mu}"],
        "pt_1":["pt_1",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_2":["pt_2",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_3":["pt_3",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_4":["pt_4",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#tau_{pt}"],
        "pt_1_fine":["pt_1",[[40,0,200]],"[Gev]","P_{T}(#tau_{1})"],
        "eta_1":["eta_1",[[60,-3,3]],"","#eta(l_{1})"],
        "phi_1":["phi_1",[[60,-3,3]],"","#phi(l_{1})"],
        "iso_1":["iso_1",[[20,0,1]],"","relIso(l_{1})"],
        "dZ_1":["dZ_1",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{1})"],
        "d0_1":["d0_1",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{1})"],
        "q_1":["q_1",[[3,-1.5,1.5]],"","charge(l_{1})"],
        "q_3":["q_3",[[3,-1.5,1.5]],"","charge(l_{3})"],
        "q_4":["q_4",[[3,-1.5,1.5]],"","charge(l_{4})"],
        "pt_2":["pt_2",[[40,0,200]],"[Gev]","P_{T}(l_{2})"],
        "eta_2":["eta_2",[[60,-3,3]],"","#eta(l_{2})"],
        "phi_2":["phi_2",[[60,-3,3]],"","#phi(l_{2})"],
        "iso_2":["iso_2",[[20,0,1]],"","relIso(l_{2})"],
        "dZ_2":["dZ_2",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{2})"],
        "d0_2":["d0_2",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{2})"],
        "q_2":["q_2",[[3,-1.5,1.5]],"","charge(l_{2})"],
	    "iso_3":["iso_3",[[20,0,1]],"","relIso(l_{3})"],
        "pt_3_ff":["pt_3",[[5,0,100]],"[Gev]","P_{T}(l_{3})"],
        "pt_3_ff_fine":["pt_3",[[40,0,100]],"[Gev]","P_{T}(l_{3})"],
        "eta_3":["eta_3",[[60,-3,3]],"","#eta(l_{3})"],
        "phi_3":["phi_3",[[60,-3,3]],"","#phi(l_{3})"],
        "dZ_3":["dZ_3",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{3})"],
        "d0_3":["d0_3",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{3})"],
        "iso_4":["iso_4",[[20,0,1]],"","relIso(l_{4})"],
        "pt_4_ff":["pt_4",[[5,0,100]],"[Gev]","P_{T}(l_{4})"],
        "pt_4_ff_fine":["pt_4",[[40,0,100]],"[Gev]","P_{T}(l_{4})"],
        "eta_4":["eta_4",[[60,-3,3]],"","#eta(l_{4})"],
        "phi_4":["phi_4",[[60,-3,3]],"","#phi(l_{4})"],
        "dZ_4":["dZ_4",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{4})"],
        "d0_4":["d0_4",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{4})"],

        "njets":["njets",[[10,-0.5,9.5]],"","nJets"],
        "jpt_1":["jpt_1",[[60,0,300]],"[GeV]","Jet^{1} P_{T}"], 
        "jeta_1":["jeta_1",[[60,-3,3]],"","Jet^{1} #eta"],
        "jpt_2":["jpt_2",[[60,0,300]],"[GeV]","Jet^{2} P_{T}"], 
        "jeta_2":["jeta_2",[[6,-3,3]],"","Jet^{2} #eta"],
        "bpt_1":["bpt_1",[[40,0,200]],"[GeV]","BJet^{1} P_{T}"], 
        "bpt_2":["bpt_2",[[40,0,200]],"[GeV]","BJet^{2} P_{T}"], 
        "nbtag":["nbtag",[[5,-0.5,4.5]],"","nBTag"],
        "beta_1":["beta_1",[[60,-3,3]],"","BJet^{1} #eta"],
        "beta_2":["beta_2",[[60,-3,3]],"","BJet^{2} #eta"],

        "met":["met",[[20,0,200]],"[GeV]","#it{p}_{T}^{miss}"], 
        "met_phi":["met_phi",[[60,-3,3]],"","#it{p}_{T}^{miss} #phi"], 
        "puppi_phi":["puppi_phi",[[60,-3,3]],"","PUPPI#it{p}_{T}^{miss} #phi"], 
        "puppimet":["puppimet",[[20,0,200]],"[GeV]","#it{p}_{T}^{miss}"], 
        "mll_fine":["mll",[[40,20,100]],"[Gev]","M_{#mu#mu}"],
        "m_vis_fine":["m_vis",[[30,50,200]],"[Gev]","m(#tau#tau)"],
        "pt_tt":["pt_tt",[[40,0,200]],"[GeV]","P_{T}(#tau#tau)"],
        "A2_DR":["H_DR",[[60,0,6]],"","#Delta R(#tau#tau)"],
        "A2_tot":["H_tot",[[30,0,200]],"[GeV]","m_{T}tot(#tau#tau)"],
        "A1_Pt":["Z_Pt",[[60,0,300]],"[Gev]","P_T(l_{1}l_{2})"],
        "A1_DR":["Z_DR",[[60,0,6]],"[Gev]","#Delta R(l_{1}l_{2})"],
        "inTimeMuon_1":["inTimeMuon_1",[[3,-1.5,1.5]],"","inTimeMuon_1"],
        "isGlobal_1":["isGlobal_1",[[3,-1.5,1.5]],"","isGlobal_1"],
        "isTracker_1":["isTracker_1",[[3,-1.5,1.5]],"","isTracker_1"],
        "looseId_1":["looseId_1",[[3,-1.5,1.5]],"","looseId_1"],
        "mediumId_1":["mediumId_1",[[3,-1.5,1.5]],"","mediumId_1"],
        "Electron_mvaFall17V2noIso_WP90_1":["Electron_mvaFall17V2noIso_WP90_1",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_1"],
        "gen_match_1":["gen_match_1",[[30,-0.5,29.5]],"","gen_match_1"],
        "inTimeMuon_2":["inTimeMuon_2",[[3,-1.5,1.5]],"","inTimeMuon_2"],
        "isGlobal_2":["isGlobal_2",[[3,-1.5,1.5]],"","isGlobal_2"],
        "isTracker_2":["isTracker_2",[[3,-1.5,1.5]],"","isTracker_2"],
        "looseId_2":["looseId_2",[[3,-1.5,1.5]],"","looseId_2"],
        "mediumId_2":["mediumId_2",[[3,-1.5,1.5]],"","mediumId_2"],
        "Electron_mvaFall17V2noIso_WP90_2":["Electron_mvaFall17V2noIso_WP90_2",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_2"],
        "gen_match_2":["gen_match_2",[[30,-0.5,29.5]],"","gen_match_2"],
        "inTimeMuon_3":["inTimeMuon_3",[[3,-1.5,1.5]],"","inTimeMuon_3"],
        "isGlobal_3":["isGlobal_3",[[3,-1.5,1.5]],"","isGlobal_3"],
        "isTracker_3":["isTracker_3",[[3,-1.5,1.5]],"","isTracker_3"],
        "looseId_3":["looseId_3",[[3,-1.5,1.5]],"","looseId_3"],
        "mediumId_3":["mediumId_3",[[3,-1.5,1.5]],"","mediumId_3"],
        "Electron_mvaFall17V2noIso_WP90_3":["Electron_mvaFall17V2noIso_WP90_3",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_3"],
        "gen_match_3":["gen_match_3",[[30,-0.5,29.5]],"","gen_match_3"],
        "inTimeMuon_4":["inTimeMuon_4",[[3,-1.5,1.5]],"","inTimeMuon_4"],
        "isGlobal_4":["isGlobal_4",[[3,-1.5,1.5]],"","isGlobal_4"],
        "isTracker_4":["isTracker_4",[[3,-1.5,1.5]],"","isTracker_4"],
        "looseId_4":["looseId_4",[[3,-1.5,1.5]],"","looseId_4"],
        "mediumId_4":["mediumId_4",[[3,-1.5,1.5]],"","mediumId_4"],
        "Electron_mvaFall17V2noIso_WP90_4":["Electron_mvaFall17V2noIso_WP90_4",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_4"],
        "gen_match_4":["gen_match_4",[[30,-0.5,29.5]],"","gen_match_4"],
        #"pt_3_ff":["pt_3",[0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,100.0,120.0],"[Gev]","P_{T}(l_{3})"],
        "pt_3_ff":["pt_3",[0.0,5.0,20.0,40.0,60.0,100.0,120.0],"[Gev]","P_{T}(l_{3})"],
        "jpt_1_ff":["jpt_1",[[5,0,100]],"[Gev]","Jet P_{T}(l_{1})"],
        "pt_3_ff_fine":["pt_3",[[40,0,100]],"[Gev]","P_{T}(l_{3})"],
        "pt_4":["pt_4",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#tau_{pt}"],
        #"pt_4_ff":["pt_4",[0.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,120.0],"[Gev]","P_{T}(l_{4})"],
        "pt_4_ff":["pt_4",[0.0,15.0,30.0,50.0,80.0,100.0,120.0],"[Gev]","P_{T}(l_{4})"],
        "jpt_2_ff":["jpt_2",[[5,0,100]],"[Gev]","Jet P_{T}(l_{2})"],
        "pt_4_ff_fine":["pt_4",[[40,0,100]],"[Gev]","P_{T}(l_{4})"],
}

#HAA_Inc_mmmt.varis = {
#        #handle : ["root variable",[binning]or[[binning]],units,label]
#        "m_vis":["m_vis",[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],"[GeV]","M_{vis}"],
#        "AMass_blinded":["AMass",[100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot} Blind"],
#        "AMass":["AMass",[0.0,20.0,40.0,60.0,80.0,100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot}"],
#        "mll":["mll",[0.0,10.,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.,100.,110.,120.,130.,140.],"[GeV]","M_{#mu#mu}"],
#        "mll_m15":["mll",[12.0,13.0,13.125,13.25,13.375,13.5,13.625,13.75,13.875,14.0,14.125,14.25,14.375,14.5,14.625,14.75,14.875,15.0,15.125,15.25,15.375,15.5,15.625,15.75,15.875,16.0,16.125,16.25,16.375,16.5,16.625,16.75,16.875,17.0 ,18.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m20":["mll",[17.0,18.0,18.125,18.25,18.375,18.5,18.625,18.75,18.875,19.0,19.125,19.25,19.375,19.5,19.625,19.75,19.875,20.0,20.125,20.25,20.375,20.5,20.625,20.75,20.875,21.0,21.125,21.25,21.375,21.5,21.625,21.75,21.875,22.0 ,23.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m25":["mll",[22.0,23.0,23.125,23.25,23.375,23.5,23.625,23.75,23.875,24.0,24.125,24.25,24.375,24.5,24.625,24.75,24.875,25.0,25.125,25.25,25.375,25.5,25.625,25.75,25.875,26.0,26.125,26.25,26.375,26.5,26.625,26.75,26.875,27.0 ,28.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m30":["mll",[27.0,28.0,28.125,28.25,28.375,28.5,28.625,28.75,28.875,29.0,29.125,29.25,29.375,29.5,29.625,29.75,29.875,30.0,30.125,30.25,30.375,30.5,30.625,30.75,30.875,31.0,31.125,31.25,31.375,31.5,31.625,31.75,31.875,32.0,33.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m35":["mll",[32.0,33.0,33.125,33.25,33.375,33.5,33.625,33.75,33.875,34.0,34.125,34.25,34.375,34.5,34.625,34.75,34.875,35.0,35.125,35.25,35.375,35.5,35.625,35.75,35.875,36.0,36.125,36.25,36.375,36.5,36.625,36.75,36.875,37.0,38.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m40":["mll",[37.0,38.0,38.125,38.25,38.375,38.5,38.625,38.75,38.875,39.0,39.125,39.25,39.375,39.5,39.625,39.75,39.875,40.0,40.125,40.25,40.375,40.5,40.625,40.75,40.875,41.0,41.125,41.25,41.375,41.5,41.625,41.75,41.875,42.0 ,43.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m45":["mll",[42.0,43.0,43.125,43.25,43.375,43.5,43.625,43.75,43.875,44.0,44.125,44.25,44.375,44.5,44.625,44.75,44.875,45.0,45.125,45.25,45.375,45.5,45.625,45.75,45.875,46.0,46.125,46.25,46.375,46.5,46.625,46.75,46.875,47.0 ,48.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m50":["mll",[47.0,48.0,48.125,48.25,48.375,48.5,48.625,48.75,48.875,49.0,49.125,49.25,49.375,49.5,49.625,49.75,49.875,50.0,50.125,50.25,50.375,50.5,50.625,50.75,50.875,51.0,51.125,51.25,51.375,51.5,51.625,51.75,51.875,52.0 ,53.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m55":["mll",[52.0,53.0,53.125,53.25,53.375,53.5,53.625,53.75,53.875,54.0,54.125,54.25,54.375,54.5,54.625,54.75,54.875,55.0,55.125,55.25,55.375,55.5,55.625,55.75,55.875,56.0,56.125,56.25,56.375,56.5,56.625,56.75,56.875,57.0 ,58.0],"[GeV]","M_{#mu#mu}"],
#        "mll_m60":["mll",[57.0,58.0,58.125,58.25,58.375,58.5,58.625,58.75,58.875,59.0,59.125,59.25,59.375,59.5,59.625,59.75,59.875,60.0,60.125,60.25,60.375,60.5,60.625,60.75,60.875,61.0,61.125,61.25,61.375,61.5,61.625,61.75,61.875,62.0 ,63.0],"[GeV]","M_{#mu#mu}"],
#        "pt_1":["pt_1",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
#        "pt_1_fine":["pt_1",[[40,0,200]],"[Gev]","P_{T}(#tau_{1})"],
#        "pt_2":["pt_2",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
#        "pt_3":["pt_3",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
#        "pt_3_ff":["pt_3",[0.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,120.0],"[Gev]","P_{T}(l_{3})"],
#        "jpt_1_ff":["jpt_1",[[5,0,100]],"[Gev]","Jet P_{T}(l_{1})"],
#        "pt_3_ff_fine":["pt_3",[[40,0,100]],"[Gev]","P_{T}(l_{3})"],
#        "pt_4":["pt_4",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#tau_{pt}"],
#        "pt_4_ff":["pt_4",[0.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,120.0],"[Gev]","P_{T}(l_{4})"],
#        "jpt_2_ff":["jpt_2",[[5,0,100]],"[Gev]","Jet P_{T}(l_{2})"],
#        "pt_4_ff_fine":["pt_4",[[40,0,100]],"[Gev]","P_{T}(l_{4})"],
#        "mll_fine":["mll",[[40,20,100]],"[Gev]","M_{#mu#mu}"],
#        "m_vis_fine":["m_vis",[[30,50,200]],"[Gev]","m(#tau#tau)"],
#        "gen_match_1":["gen_match_1",[[30,-0.5,29.5]],"","gen_match_1"],
#        "gen_match_2":["gen_match_2",[[30,-0.5,29.5]],"","gen_match_2"],
#        "gen_match_3":["gen_match_3",[[30,-0.5,29.5]],"","gen_match_3"],
#        "gen_match_4":["gen_match_4",[[30,-0.5,29.5]],"","gen_match_4"],
#}

HAA_Inc_mmmt_noIds = copy.deepcopy(HAA_Inc_mmmt)

HAA_Inc_mmmt_noIds.name = "mmmt_inclusive_noIds"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_noIds.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3]]

HAA_Inc_mmmt_noIds.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.]]


HAA_Inc_mmmt_charge = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_charge.name = "mmmt_inclusive_charge"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_charge.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0]]
HAA_Inc_mmmt_charge.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.]]


HAA_Inc_mmmt_iso = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_iso.name = "mmmt_inclusive_iso"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_iso.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_iso.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15]]
#["mediumId_3",">=",1],
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#,["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint] this is cutting most events! 
HAA_Inc_mmmt_iso.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #



HAA_Inc_mmmt_muon = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_muon.name = "mmmt_inclusive_muon"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_muon.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_muon.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],
[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#,["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint] this is cutting most events! 
HAA_Inc_mmmt_muon.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_tauid = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_tauid.name = "mmmt_inclusive_tauid"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_tauid.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_tauid.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDecayModeNewDMs_4",">",0.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]this is cutting most events! 
HAA_Inc_mmmt_tauid.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_medium = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_medium.name = "mmmt_inclusive_deeptauid_medium"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_medium.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_medium.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",16.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_medium.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_loose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_loose.name = "mmmt_inclusive_deeptauid_loose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_loose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_loose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",8.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_loose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #
allcats=[]

HAA_Inc_mmmt_deeptauid_vloose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_vloose.name = "mmmt_inclusive_deeptauid_vloose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_vloose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_vloose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",4.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_vloose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_vvloose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_vvloose.name = "mmmt_inclusive_deeptauid_vvloose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_vvloose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_vvloose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",2.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_vvloose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #


HAA_Inc_mmmt_FF_SS_1_tight = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_FF_SS_1_tight.name = "mmmt_inclusive_FF_SS_1_tight"
HAA_Inc_mmmt_FF_SS_1_tight.cuts["categoryCuts"]= [
["cat","==",6],[["EQT"],["q_3","q_4"],"mult",">",0],
["iso_3","<=",0.15],["nbtag","<",1.0],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]]

#invert muon criteria
HAA_Inc_mmmt_FF_SS_1_loose = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_1_loose.name = "mmmt_inclusive_FF_SS_1_loose"
HAA_Inc_mmmt_FF_SS_1_loose.cuts["categoryCuts"]= [
["cat","==",6],[["EQT"],["q_3","q_4"],"mult",">",0],["nbtag","<",1.0]]
#[["OR"],["iso_3",">",0.15],["mediumId_3","<",1],["isGlobal_3","<",1],["isTracker_3","<",1]],
#["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]]

HAA_Inc_mmmt_FF_SS_2_tight = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_2_tight.name = "mmmt_inclusive_FF_SS_2_tight"
HAA_Inc_mmmt_FF_SS_2_tight.cuts["categoryCuts"]= [
["cat","==",6],[["EQT"],["q_3","q_4"],"mult",">",0],
["iso_3","<=",0.15],["mediumId_3",">=",1],["nbtag","<",1.0],
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]]

#invert tau criteria
HAA_Inc_mmmt_FF_SS_2_loose = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_2_loose.name = "mmmt_inclusive_FF_SS_2_loose"
HAA_Inc_mmmt_FF_SS_2_loose.cuts["categoryCuts"]= [
["cat","==",6],[["EQT"],["q_3","q_4"],"mult",">",0],["nbtag","<",1.0],
["iso_3","<=",0.15],["mediumId_3",">=",1]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
#["idDeepTau2017v2p1VSjet_4","<",ffworkingpoint]]




allcats.append(HAA_Inc_mmmt)
#current ff here may 8th 2020
allcats.append(HAA_Inc_mmmt_FF_SS_1_tight)
allcats.append(HAA_Inc_mmmt_FF_SS_1_loose)
allcats.append(HAA_Inc_mmmt_FF_SS_2_tight)
allcats.append(HAA_Inc_mmmt_FF_SS_2_loose)

#allcats.append(HAA_Inc_mmmt_deeptauid_medium)
#allcats.append(HAA_Inc_mmmt_deeptauid_loose)
#allcats.append(HAA_Inc_mmmt_deeptauid_vloose)
#allcats.append(HAA_Inc_mmmt_deeptauid_vvloose)

#allcats.append(HAA_Inc_mmmt_noIds)
#allcats.append(HAA_Inc_mmmt_charge)
#allcats.append(HAA_Inc_mmmt_iso)
#allcats.append(HAA_Inc_mmmt_muon)
#allcats.append(HAA_Inc_mmmt_deeptauid)
