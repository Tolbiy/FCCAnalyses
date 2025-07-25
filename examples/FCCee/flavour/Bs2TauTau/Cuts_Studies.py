import matplotlib.pyplot as plt
import numpy as np
import uproot
import ROOT as r
import pandas as pd
from termcolor import colored
import json

r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Cut_OnLists_float(ROOT::VecOps::RVec<float> List, float Threshold1, float Threshold2){
        ROOT::VecOps::RVec<float> newList;
        for (size_t i=0; i<List.size(); ++i){
            if (List.at(i) > Threshold1 && List.at(i) < Threshold2) {
                newList.push_back(List.at(i));
            }
        }
        return newList;
    }
''')


def Build_FilterArg(CutList):
    
    FullList = ""
    
    if len(CutList) == 1:
        FullList = CutList
    
    else:
        for Cut in CutList:
            if Cut == CutList[0]:
                FullList = Cut
            else:
                FullList += f" && {Cut}"

    return FullList


#-----------------------------------------------------------------------------------------------------------------------------------

################################# Warning: we may be reaching the statistical limit for conducting these studies ##################################################################

#Modes to be scanned
modes = ["p8_ee_Zbb_ecm91",
         "p8_ee_Zcc_ecm91",
         "p8_ee_Zss_ecm91",
         "p8_ee_Zud_ecm91",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"]

#Storage Values
Ntot = {}
NCut = {}

NTot = 0
NCuts = 0

#The Cut list to be studied
CutList = ["EVT_ThrustEmin_NTau23PiCand > 1", 
           "EVT_ThrustEmin_Eneutral < 10",
           "EVT_ThrustEmin_Nneutral <= 12", #not much gain since Eneutral is strongly correlated with Nneutral
           "Tau23PiCandidates_rho1mass.size() > 0",
           "Tau23PiCandidates_rho2mass.size() > 0",
           "Tau23PiCandidates_q.size() == 2", #Lost a lot of signal (25%->42%) while not much has changed for background
           "Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0"] 
           
           #Tau23PiCandidates_rho1mass.at(0) < 0.9


#prevents adding same lines in case calling multiple time the scipt 
Previous = {}
AlreadyDone = False
with open("Cuts_Studied.json","r") as ofile:
    Previous = json.load(ofile)
    
    for Cuts in Previous.keys():
        if Build_FilterArg(CutList) == Cuts:
            AlreadyDone = True
            print("Already studied")
            exit()



for mode in modes:

    print(f"_____________ {mode} _______________")

    #Load the number of cuts before any cuts
    Ntot[mode] = Previous["Before"][mode]
    NTot += Ntot[mode]

    #Load the each files into an rdf
    inf = r.TFile(f"/afs/cern.ch/work/t/tomonnar/public/{mode}.root", 'read')
    inTree = inf.Get("events")
    rdf = r.RDataFrame(inTree)
    
    #flatten the lists, cut and count the survivors
    rdf = rdf.Redefine("Tau23PiCandidates_rho1mass","Cut_OnLists_float(Tau23PiCandidates_rho1mass,0.55,1.0)")
    rdf = rdf.Redefine("Tau23PiCandidates_rho2mass","Cut_OnLists_float(Tau23PiCandidates_rho2mass,0.55,1.0)")
    rdf = rdf.Filter(Build_FilterArg(CutList))
    NCut[mode] = rdf.Count().GetValue()
    NCuts += NCut[mode]

    if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":
        color = "blue"
    else:
        color = "red"
    
    #Efficiencies of how many got rejected so red should be close to 100% and blue close to 0%
    print(colored(f"Cuts eff = {round(100*(Ntot[mode]-NCut[mode])/Ntot[mode],5)} %",color))

print("\n-----------------------------------------")
print(colored(f"=====> Total eff = {round(100*(NTot-NCuts)/NTot,5)} %","magenta"))
print("-----------------------------------------\n")


#Store these numbers
with open("Cuts_Studied.json","w") as ofile:
    Previous[Build_FilterArg(CutList)] = NCut
    json.dump(Previous, ofile)

     
            
