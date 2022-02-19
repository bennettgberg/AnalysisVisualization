#########################
#Author: Sam Higginbotham
'''

* File Name :Weights.py

* Purpose :Contains the mapping between weights and processes

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
CommonWeights={}
CommonWeights["lumi2018"]=[59740]
CommonWeights["lumi2016"]=[35900]#pb^-1
CommonWeights["lumi2017"]=[41800]#pb^-1
#CommonWeights["lumi"]=["1"]#pb^-1
#CommonWeights["mcweight"]=["__WEIGHT__*GENWEIGHT*puweight"]
#CommonWeights["mcweight"]=["LHEweight*weight"]
CommonWeights["string"]=[["Generator_weight"],["weight"]]

jet_inclusive_samples = {}

jet_inclusive_samples["2018"] = [
                        "DYJetsToLL", #bpg added ???
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "W1JetsToLNu" ,
                        "W3JetsToLNu",
                        ]

jet_inclusive_samples["2017"] = ["DYJetsToLL_ext1",
                        "DYJetsToLL", #bpg added ???
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "WJetsToLNu_ext1" ,
                        "W1JetsToLNu" ,
                        #"W2JetsToLNu_ext1" ,
                        "W3JetsToLNu",
                        #"W4JetsToLNu_ext1",
                        #"W4JetsToLNu_ext2"
                        #bpg adding all the below ----v
                        "DYJetsToLL_ext1", #bpg added ???
                        "DY1JetsToLL_ext1" ,
                        "DY2JetsToLL_ext1" ,
                        "DY3JetsToLL_ext1" ,
                        #"DY4JetsToLL_ext1" , #DNE fsr
                        ]

jet_inclusive_samples["2016"] = ["DYJetsToLLext1",
                        "DYJetsToLLext2" ,
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "WJetsToLNu_ext2" ,
                        "W1JetsToLNu" ,
                        "W2JetsToLNu_ext1" ,
                        "W3JetsToLNu",
                        "W3JetsToLNu_ext1",
                        "W4JetsToLNu",
                        "W4JetsToLNu_ext1",
                        "W4JetsToLNu_ext2"
                        ]

# jet_inclusive_samples["2016"] = ["DYJetsToLLext1",    #v6
#                         "DYJetsToLLext2" ,
#                         "DY1JetsToLL" ,
#                         "DY2JetsToLL" ,
#                         "DY3JetsToLL" ,
#                         "DY4JetsToLL" ,
#                         "WJetsToLNu" ,
#                         "WJetsToLNuext" ,
#                         "WJetsToLNu_ext2" ,
#                         "W1JetsToLNu" ,
#                         "W2JetsToLNuext1" ,
#                         "W3JetsToLNu",
#                         "W4JetsToLNu_ext1",
#                         "W4JetsToLNu_ext2"
#                         ]
jetIncOnly = {}
jetIncOnly["2018"] = [["DYJetsToLL"],["WJetsToLNu"]]
#jetIncOnly["2017"] = [["DYJetsToLL_ext1"],["WJetsToLNu","WJetsToLNu_ext1"]] ##wtf is this -bpg
jetIncOnly["2017"] = [["DYJetsToLL","DYJetsToLL_ext1"],["WJetsToLNu","WJetsToLNu_ext1"]]
#jetIncOnly["2016"] = [["DYJetsToLLext1","DYJetsToLLext2"],["WJetsToLNu","WJetsToLNuext","WJetsToLNu_ext2"]] #v6
jetIncOnly["2016"] = [["DYJetsToLLext1","DYJetsToLLext2"],["WJetsToLNu","WJetsToLNu_ext2"]]


#To Do here: combine the extended samples for everything before ... then feed into cut array
#STEPS
#combine all the extended samples first!

#SMHTTWeights={}
#SMHTTWeights["data_obs"]=[""]
#SMHTTWeights["TTT"]=CommonWeights["mcweight"]
#SMHTTWeights["TTL"]=CommonWeights["mcweight"]
#SMHTTWeights["ZTT"]=CommonWeights["mcweight"]
#SMHTTWeights["ZTL"]=CommonWeights["mcweight"]
#SMHTTWeights["ZL"]=CommonWeights["mcweight"]

#Where do I place energy scale shifts?
HAAWeights={}
#HAAWeights["data_obs"]=[""]
#HAAWeights["TTT"]=CommonWeights["mcweight"]
#HAAWeights["TTL"]=CommonWeights["mcweight"]
#HAAWeights["ZTT"]=CommonWeights["mcweight"]
#HAAWeights["ZTL"]=CommonWeights["mcweight"]
#HAAWeights["ZL"]=CommonWeights["mcweight"]
HAAWeights["a15"]="16*00005."
HAAWeights["a20"]="16*00005."
HAAWeights["a25"]="16*00005."
HAAWeights["a30"]="16*00005."
HAAWeights["a35"]="16*00005."
HAAWeights["a40"]="16*00005."
HAAWeights["a45"]="16*00005."
HAAWeights["a50"]="16*00005."
HAAWeights["a55"]="16*00005."
HAAWeights["a60"]="16*00005."

HAAWeights["TES"]=["(((pt_4*1.011)>40&&decayMode_4==0)||((pt_4*0.995)>40&&decayMode_4==1)||((pt_4*0.988)>40&&decayMode_4==10))&&(((pt_3*1.011)>50&&decayMode_3==0)||((pt_3*0.995)>50&&decayMode_3==1)||((pt_3*0.988)>50&&decayMode_3==10))"]
