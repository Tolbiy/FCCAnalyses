import matplotlib.pyplot as plt
import numpy as np
import uproot
import ROOT as r
import pandas as pd
from termcolor import colored
import json

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Function to recover a tau candidates vertices property from a full vertices property list
#Helps then to filter out any events whose tau candidates didn't match the requirements: Filter(newList.size()<Threshold) for instance 

r.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> Scan_Tau_Vertex_Prop(ROOT::VecOps::RVec<float> Obs_perParticle, ROOT::VecOps::RVec<int> Tauid){
        ROOT::VecOps::RVec<float> newList;
        for (size_t i=0; i<Tauid.size(); ++i){
            newList.push_back(Obs_perParticle.at(Tauid.at(i)));
        }
        return newList;
    }
''')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Some bookeeping for looping and the various decay mode and then drawing them

modes = ["p8_ee_Zbb_ecm91",
         "p8_ee_Zcc_ecm91",
         "p8_ee_Zss_ecm91",
         "p8_ee_Zud_ecm91",
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU"]

colors = {"p8_ee_Zbb_ecm91":                          [617,618,619],
         "p8_ee_Zcc_ecm91":                           [826,825,824],
         "p8_ee_Zss_ecm91":                           [401,402,403],
         "p8_ee_Zud_ecm91":                           [920,921,922],
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":  [600,601,602]}

labels = {"p8_ee_Zbb_ecm91":                          "#it{b#bar{b}}",  #it{Z}#rightarrow#it{b#bar{b}}                               
         "p8_ee_Zcc_ecm91":                           "#it{c#bar{c}}",  #it{Z}#rightarrow#it{c#bar{c}}
         "p8_ee_Zss_ecm91":                           "#it{s#bar{s}}",  #it{Z}#rightarrow#it{s#bar{s}}
         "p8_ee_Zud_ecm91":                           "#it{q#bar{q}}, #it{q}=#it{u},#it{d}", #it{Z}#rightarrow#it{q#bar{q}}, #it{q} = #it{u}, #it{d}
         "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":  "Signal"}

#The full branch list obtained from analysis_stage1.py should be put else where for readability purpose
branchList_Stock = [
    "MC_PDG","MC_M1","MC_M2","MC_n","MC_D1","MC_D2",

    "EVT_ThrustEmin_E",            "EVT_ThrustEmax_E",
    "EVT_ThrustEmin_Echarged",     "EVT_ThrustEmax_Echarged",
    "EVT_ThrustEmin_Eneutral",     "EVT_ThrustEmax_Eneutral",
    "EVT_ThrustEmin_N",            "EVT_ThrustEmax_N",
    "EVT_ThrustEmin_Ncharged",     "EVT_ThrustEmax_Ncharged",
    "EVT_ThrustEmin_Nneutral",     "EVT_ThrustEmax_Nneutral",
    "EVT_ThrustEmin_NDV",          "EVT_ThrustEmax_NDV",
    "EVT_ThrustEmin_NTau23PiCand", "EVT_ThrustEmax_NTau23PiCand",
    "EVT_Thrust_Mag",
    "EVT_Thrust_X",  "EVT_Thrust_XErr",
    "EVT_Thrust_Y",  "EVT_Thrust_YErr",
    "EVT_Thrust_Z",  "EVT_Thrust_ZErr",

    "EVT_NtracksPV", "EVT_NVertex", "EVT_NTau23Pi",

    "EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",

    "MC_Vertex_x", "MC_Vertex_y", "MC_Vertex_z",
    "MC_Vertex_ntrk", "MC_Vertex_n",

    "Vertex_x", "Vertex_y", "Vertex_z",
    "Vertex_xErr", "Vertex_yErr", "Vertex_zErr",
    "Vertex_isPV", "Vertex_ntrk", "Vertex_chi2", "Vertex_n",
    "Vertex_thrust_angle", "Vertex_thrusthemis_emin", "Vertex_thrusthemis_emax",

    "Vertex_d2PV", "Vertex_d2PVx", "Vertex_d2PVy", "Vertex_d2PVz",
    "Vertex_d2PVErr", "Vertex_d2PVxErr", "Vertex_d2PVyErr", "Vertex_d2PVzErr",
    "Vertex_mass",
    "DV_d0","DV_z0",

    "Tau23PiCandidates_mass", "Tau23PiCandidates_vertex", "Tau23PiCandidates_mcvertex", "Tau23PiCandidates_B",
    "Tau23PiCandidates_px", "Tau23PiCandidates_py", "Tau23PiCandidates_pz", "Tau23PiCandidates_p", "Tau23PiCandidates_q",
    "Tau23PiCandidates_d0",  "Tau23PiCandidates_z0","Tau23PiCandidates_anglethrust",

    "Tau23PiCandidates_rho1px", "Tau23PiCandidates_rho1py", "Tau23PiCandidates_rho1pz","Tau23PiCandidates_rho1mass",
    "Tau23PiCandidates_rho2px", "Tau23PiCandidates_rho2py", "Tau23PiCandidates_rho2pz","Tau23PiCandidates_rho2mass",

    "Tau23PiCandidates_pion1px", "Tau23PiCandidates_pion1py", "Tau23PiCandidates_pion1pz",
    "Tau23PiCandidates_pion1p", "Tau23PiCandidates_pion1q", "Tau23PiCandidates_pion1d0", "Tau23PiCandidates_pion1z0",
    "Tau23PiCandidates_pion2px", "Tau23PiCandidates_pion2py", "Tau23PiCandidates_pion2pz",
    "Tau23PiCandidates_pion2p", "Tau23PiCandidates_pion2q", "Tau23PiCandidates_pion2d0", "Tau23PiCandidates_pion2z0",
    "Tau23PiCandidates_pion3px", "Tau23PiCandidates_pion3py", "Tau23PiCandidates_pion3pz",
    "Tau23PiCandidates_pion3p", "Tau23PiCandidates_pion3q", "Tau23PiCandidates_pion3d0", "Tau23PiCandidates_pion3z0",

    "Bs_Mass", "Delta_M"
]

#Invalid branch: MC_Vertex_PDG, MC_Vertex_PDGmother, MC_Vertex_PDGgmother


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Some newly defined branches to test some cuts


def Scan_Tau_Bs_Prop(rdf):

    rdf = rdf.Filter("Tau23PiCandidates_q.size() == 2").Filter("Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0")
    #Careful with these cuts here because we only accept events with exactly 2 opposite charged taus (no guarantee that they are coming from Bs)
    #Seems that the exact 2 taus removes some amount but not that much to be checked in more details with the first set of histograms
    #To be applied in general not only to study Bs, ideally define a reconstruction method to obtain the columns needed
    #To sum it up this is too rough for now and need a correct way to get approximate Bs properties
    
    #Adding the properties of the pions linked to each tau
    rdf = rdf.Define("Tau23PiCandidates_firstTau_pion1p",  "Tau23PiCandidates_pion1p.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion2p",  "Tau23PiCandidates_pion2p.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion3p",  "Tau23PiCandidates_pion3p.at(0)")\
             \
             .Define("Tau23PiCandidates_secondTau_pion1p",  "Tau23PiCandidates_pion1p.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion2p",  "Tau23PiCandidates_pion2p.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion3p",  "Tau23PiCandidates_pion3p.at(1)")\
             \
             .Define("Tau23PiCandidates_firstTau_pion1E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion1p,2) + pow(0.13957,2))")\
             .Define("Tau23PiCandidates_firstTau_pion2E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion2p,2) + pow(0.13957,2))")\
             .Define("Tau23PiCandidates_firstTau_pion3E",  "sqrt(pow(Tau23PiCandidates_firstTau_pion3p,2) + pow(0.13957,2))")\
             \
             .Define("Tau23PiCandidates_secondTau_pion1E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion1p,2) + pow(0.13957,2))")\
             .Define("Tau23PiCandidates_secondTau_pion2E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion2p,2) + pow(0.13957,2))")\
             .Define("Tau23PiCandidates_secondTau_pion3E",  "sqrt(pow(Tau23PiCandidates_secondTau_pion3p,2) + pow(0.13957,2))")\
             \
             .Define("Tau23PiCandidates_firstTau_pion1px",  "Tau23PiCandidates_pion1px.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion1py",  "Tau23PiCandidates_pion1py.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion1pz",  "Tau23PiCandidates_pion1pz.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion2px",  "Tau23PiCandidates_pion2px.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion2py",  "Tau23PiCandidates_pion2py.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion2pz",  "Tau23PiCandidates_pion2pz.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion3px",  "Tau23PiCandidates_pion3px.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion3py",  "Tau23PiCandidates_pion3py.at(0)")\
             .Define("Tau23PiCandidates_firstTau_pion3pz",  "Tau23PiCandidates_pion3pz.at(0)")\
             \
             .Define("Tau23PiCandidates_secondTau_pion1px",  "Tau23PiCandidates_pion1px.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion1py",  "Tau23PiCandidates_pion1py.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion1pz",  "Tau23PiCandidates_pion1pz.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion2px",  "Tau23PiCandidates_pion2px.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion2py",  "Tau23PiCandidates_pion2py.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion2pz",  "Tau23PiCandidates_pion2pz.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion3px",  "Tau23PiCandidates_pion3px.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion3py",  "Tau23PiCandidates_pion3py.at(1)")\
             .Define("Tau23PiCandidates_secondTau_pion3pz",  "Tau23PiCandidates_pion3pz.at(1)")\
             \
             .Define("Bs_Mass","sqrt(pow(Tau23PiCandidates_firstTau_pion1E+Tau23PiCandidates_firstTau_pion2E+Tau23PiCandidates_firstTau_pion3E+Tau23PiCandidates_secondTau_pion1E+Tau23PiCandidates_secondTau_pion2E+Tau23PiCandidates_secondTau_pion3E,2) - "+\
                                    "pow(Tau23PiCandidates_firstTau_pion1px+Tau23PiCandidates_firstTau_pion2px+Tau23PiCandidates_firstTau_pion3px+Tau23PiCandidates_secondTau_pion1px+Tau23PiCandidates_secondTau_pion2px+Tau23PiCandidates_secondTau_pion3px,2) - "+\
                                    "pow(Tau23PiCandidates_firstTau_pion1py+Tau23PiCandidates_firstTau_pion2py+Tau23PiCandidates_firstTau_pion3py+Tau23PiCandidates_secondTau_pion1py+Tau23PiCandidates_secondTau_pion2py+Tau23PiCandidates_secondTau_pion3py,2) - "+\
                                    "pow(Tau23PiCandidates_firstTau_pion1pz+Tau23PiCandidates_firstTau_pion2pz+Tau23PiCandidates_firstTau_pion3pz+Tau23PiCandidates_secondTau_pion1pz+Tau23PiCandidates_secondTau_pion2pz+Tau23PiCandidates_secondTau_pion3pz,2))")\
             \
             .Define("Delta_M","Bs_Mass - Tau23PiCandidates_mass.at(0) - Tau23PiCandidates_mass.at(1)")\
             \
             .Define("TausDotted","(Tau23PiCandidates_firstTau_pion1E+Tau23PiCandidates_firstTau_pion2E+Tau23PiCandidates_firstTau_pion3E)*(Tau23PiCandidates_secondTau_pion1E*Tau23PiCandidates_secondTau_pion2E*Tau23PiCandidates_secondTau_pion3E) - "+\
                                  "(Tau23PiCandidates_firstTau_pion1px+Tau23PiCandidates_firstTau_pion2px+Tau23PiCandidates_firstTau_pion3px)*(Tau23PiCandidates_secondTau_pion1px*Tau23PiCandidates_secondTau_pion2px*Tau23PiCandidates_secondTau_pion3px) - "+\
                                  "(Tau23PiCandidates_firstTau_pion1py+Tau23PiCandidates_firstTau_pion2py+Tau23PiCandidates_firstTau_pion3py)*(Tau23PiCandidates_secondTau_pion1py*Tau23PiCandidates_secondTau_pion2py*Tau23PiCandidates_secondTau_pion3py) - "+\
                                  "(Tau23PiCandidates_firstTau_pion1pz+Tau23PiCandidates_firstTau_pion2pz+Tau23PiCandidates_firstTau_pion3pz)*(Tau23PiCandidates_secondTau_pion1pz*Tau23PiCandidates_secondTau_pion2pz*Tau23PiCandidates_secondTau_pion3pz)   ")\
             \
             \
             .Define("Tau23PiCandidates_firstTau_B","Tau23PiCandidates_B.at(0)").Define("Tau23PiCandidates_secondTau_B","Tau23PiCandidates_B.at(1)")\
             \
             .Define("Tau_px","Tau23PiCandidates_firstTau_pion1px+Tau23PiCandidates_firstTau_pion2px+Tau23PiCandidates_firstTau_pion3px")

#----------------------------------------------------------------------------------------------------------------------------------------------------

def Scan_Tau_Res_Prop(rdf):
    
    rdf = rdf.Define("Tau23PiCandidates_pion1E",  "sqrt(pow(Tau23PiCandidates_pion1p,2.0) + pow(0.13957,2.0))")\
             .Define("Tau23PiCandidates_pion2E",  "sqrt(pow(Tau23PiCandidates_pion2p,2.0) + pow(0.13957,2.0))")\
             .Define("Tau23PiCandidates_pion3E",  "sqrt(pow(Tau23PiCandidates_pion3p,2.0) + pow(0.13957,2.0))")\
             \
             .Define("a1_Mass","sqrt(pow(Tau23PiCandidates_pion1E  + Tau23PiCandidates_pion2E  + Tau23PiCandidates_pion3E ,2.0)"+\
                                  "- pow(Tau23PiCandidates_pion1px + Tau23PiCandidates_pion2px + Tau23PiCandidates_pion3px,2.0)"+\
                                  "- pow(Tau23PiCandidates_pion1py + Tau23PiCandidates_pion2py + Tau23PiCandidates_pion3py,2.0)"+\
                                  "- pow(Tau23PiCandidates_pion1pz + Tau23PiCandidates_pion2pz + Tau23PiCandidates_pion3pz,2.0))")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def Scan_Tau_Vertex(rdf):
               
    rdf = rdf.Filter("Tau23PiCandidates_q.size() == 2").Filter("Tau23PiCandidates_q.at(0) + Tau23PiCandidates_q.at(1) == 0")
    rdf = rdf.Define("Tau23PiCandidates_d2PVx","Scan_Tau_Vertex_Prop(Vertex_d2PVx,Tau23PiCandidates_vertex)")
    rdf = rdf.Define("Tau23PiCandidates_d2PVy","Scan_Tau_Vertex_Prop(Vertex_d2PVy,Tau23PiCandidates_vertex)")
    rdf = rdf.Define("Tau23PiCandidates_d2PVz","Scan_Tau_Vertex_Prop(Vertex_d2PVz,Tau23PiCandidates_vertex)")
    rdf = rdf.Define("Tau23PiCandidates_d2PVT","sqrt(pow(Tau23PiCandidates_d2PVx,2.0)+pow(Tau23PiCandidates_d2PVy,2.0))")
            
    rdf = rdf.Define("Tau23PiCandidates_Vx","Scan_Tau_Vertex_Prop(Vertex_x,Tau23PiCandidates_vertex)")
    rdf = rdf.Define("Tau23PiCandidates_Vy","Scan_Tau_Vertex_Prop(Vertex_y,Tau23PiCandidates_vertex)")
    rdf = rdf.Define("Tau23PiCandidates_Vz","Scan_Tau_Vertex_Prop(Vertex_z,Tau23PiCandidates_vertex)")

    rdf = rdf.Define("Tau23PiCandidates_d2TauVertices","sqrt(pow(Tau23PiCandidates_Vx.at(0)-Tau23PiCandidates_Vx.at(1),2.0)"+\
                                                            "pow(Tau23PiCandidates_Vx.at(0)-Tau23PiCandidates_Vx.at(1),2.0)"+\
                                                            "pow(Tau23PiCandidates_Vx.at(0)-Tau23PiCandidates_Vx.at(1),2.0))")
      
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#===============================================================================================================================================

#Select the branch or the list of branch to be scanned
#=====================================================================================
branchList = ["EVT_missingEnergy_x","EVT_missingEnergy_y","EVT_missingEnergy_z","EVT_missingEnergy_e"]
#=====================================================================================

#Load the histogram Min and Max found previously and tweaked
BranchLimits = {}
with open("BranchMinMax.json","r") as MFile:
    BranchLimits = json.load(MFile)

#Start the loop on the branches to be scanned
for Branch in branchList:

    print(f"______________ {Branch} ________________")

    Delta = BranchLimits[Branch][1]-BranchLimits[Branch][0]
    Min = BranchLimits[Branch][0] - 0.05*Delta
    Max = BranchLimits[Branch][1] + 0.05*Delta
    Height = 0

    #Some variables to correctly normalise HStack
    h_bkg_list = {}
    Norm_bkg = 0
    
    hs = r.THStack("hs"," ")
    SigShape = r.TH1F()

    #Start to loop on the decay modes
    for mode in modes:
        if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":
            inrep = mode
        else:
            inrep = "Bs2TauTau/"+mode
        inf = r.TFile(f"/afs/cern.ch/work/t/tomonnar/public/{inrep}.root", 'read')
        inTree = inf.Get("events") #Optimise by only loading the actual branch, careful for cuts if only loading the considered branch/ f"events:{Branch}"
        rdf = r.RDataFrame(inTree)
        rdf = rdf.Filter("EVT_ThrustEmin_NTau23PiCand>1")

        #Use the custom branches defined before
        if Branch == "Bs_Mass" or "Delta_M" or "Tau23PiCandidates_firstTau_pion1E" or "TausDotted" or "Tau23PiCandidates_firstTau_B" or "Tau23PiCandidates_secondTau_B" or "Tau_px":
            Scan_Tau_Bs_Prop(rdf)
        if Branch == "a1_Mass" or "Tau23PiCandidates_pion1E":
            Scan_Tau_Res_Prop(rdf)
        if Branch == "Tau23PiCandidates_d2PVz" or "Tau23PiCandidates_d2PVy" or "Tau23PiCandidates_d2PVx" or "Tau23PiCandidates_d2PVT" or "Tau23PiCandidates_d2TauVertices":
            Scan_Tau_Vertex(rdf)
            


        #Load the histograms and normalise them (sig and bkg separately) for drawing
        if mode != "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":

            h_bkg_list[mode] = rdf.Histo1D((f"{Branch}_{mode}"," ",100,Min,Max), Branch)
            Norm_bkg += h_bkg_list[mode].GetEntries()

            h_bkg_list[mode].SetFillColor(colors[mode][1])
            h_bkg_list[mode].SetLineColor(colors[mode][1])

        else:
            h_m = rdf.Histo1D((f"{Branch}_{mode}"," ",100,Min,Max), Branch)

            h_m.SetLineColor(colors[mode][1])
            h_m.SetLineWidth(3)
            h_m.Scale(1.0/h_m.GetEntries())
            Height = h_m.GetMaximum()
            r.gROOT.cd()
            SigShape = h_m.GetValue().Clone()

    for mode in h_bkg_list.keys():
        
        h_bkg_list[mode].Scale(1.0/Norm_bkg)
        r.gROOT.cd()
        hs.Add(h_bkg_list[mode].GetValue().Clone())

    if hs.GetMaximum() > Height:
        Height = hs.GetMaximum()

    #Drawing
    c = r.TCanvas(f"{Branch}", "c", 900, 700)
    c.SetGrid()

    hs.Draw("hist")
    hs.SetMaximum(1.1*Height)
    hs.GetXaxis().SetTitle(f"{Branch}")
    hs.GetYaxis().SetTitle("Normalised Count")
    SigShape.Draw("same hist")

    #Leg = r.gPad.BuildLegend(0.5,0.65,0.9,0.9)
    Leg = r.TLegend(0.75,0.68,0.89,0.89)

    Boxes = []
    for mode in modes:
        Box = r.TBox(0.25,0.25,0.75,0.75)
        Box.SetLineColor(colors[mode][1])
        Box.SetFillColor(colors[mode][1])
        #Box.Draw("same")
        Boxes.append(Box)
        r.gROOT.cd()
        if mode != "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":
            DrawOpt = "f"
        else:
            DrawOpt = "l"
        Leg.AddEntry(Box,labels[mode],DrawOpt)

    Leg.Draw("same")

    c.Draw()
    c.SaveAs(f"Scanned_2TauCut/{Branch}.pdf")





    #rdf = r.RDataFrame(inTree)

    #if mode == "p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU":
    #    print(colored(f"After filtering 2 tau candidates in the signal hemisphere and rescaling: {Bs_Production*BR_Bs2TauTau*(BR_Tau2HADNU)**2*rdf.Count().GetValue()}","blue"))
    #else:
    #    print(colored(f"After filtering 2 tau candidates in the signal hemisphere: {rdf.Count().GetValue()}","red"))