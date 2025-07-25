import matplotlib.pyplot as plt
import numpy as np
import uproot
import ROOT as r
import pandas as pd
from termcolor import colored
import json



#Modes to be scanned
modes = ["p8_ee_Zbb_ecm91",
         "p8_ee_Zcc_ecm91",
         "p8_ee_Zss_ecm91",
         "p8_ee_Zud_ecm91",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"]


Previous = {}

with open("Cuts_Studied.json","r") as ofile:
    Previous = json.load(ofile)


print(colored("\n---------Recaping-------------\n","magenta"))
for CutsListed in Previous.keys():
    
    print(colored(CutsListed.replace(" && ","\n"),"yellow"))
    
    NBkg_Before = 0
    NBkg_After = 0
    
    NSig_Before = 0
    NSig_After = 0
    for mode in modes:
        
        if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":
            NSig_Before += Previous["Before"][mode]
            NSig_After += Previous[CutsListed][mode]
        else:
            NBkg_Before += Previous["Before"][mode]
            NBkg_After += Previous[CutsListed][mode]
    print("Signal efficiency: ",colored(f"{round(100*(NSig_Before-NSig_After)/NSig_Before,5)} %","blue"))
    print("Background efficiency: ",colored(f"{round(100*(NBkg_Before-NBkg_After)/NBkg_Before,5)} %","red"))
    print("\n")
