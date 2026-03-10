#analysis_stage1

# list of samples to process
processList_full = {
    'p8_ee_Zbb_ecm91':{'chunks':10,'fraction':0.01},
    'p8_ee_Zcc_ecm91':{'chunks':10,'fraction':0.01},
    'p8_ee_Zss_ecm91':{'chunks':10,'fraction':0.01},
    'p8_ee_Zud_ecm91':{'chunks':10,'fraction':0.01},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':10},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':20,'fraction':1},
}

#processList_full = {
#    'p8_ee_Zbb_ecm91':{'chunks':100},
#    'p8_ee_Zcc_ecm91':{'chunks':100},
#    'p8_ee_Zss_ecm91':{'chunks':100},
#    'p8_ee_Zud_ecm91':{'chunks':100},
#    'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':20},
#    'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':20},
#}

processList_test = {
    'p8_ee_Zbb_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':1, 'fraction':0.01},
    #'p8_ee_Zss_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zud_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zcc_ecm91':{'chunks':1, 'fraction':0.000002},
}

nCPUS       = 8
runBatch    = False
batchQueue  = "nextweek"
compGroup   = "group_u_FCC.local_gen"

processList  = processList_full
if not runBatch:
    processList  = processList_test

# tag for MC production campaign
prodTag     = "FCCee/winter2023/IDEA/"

# if runBatch = True, save output on eos
outputDirEos   = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_mumu_noFilter"

# if runBatch = False, save output locally
outputDir   = "../../../../../../../../../../../work/t/tomonnar/public/Bs2TauTau/mumu/bkg_studies/"
#"fccanalysis_output/" To be used with runBatch = True

includePaths = ["functions.h"]

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
               
               #-- Trick to quickly test any script --
               #/!\ if multiple thread the result might be inconsistent over multiple runs because of how the events are distributed
               .Define("nevent","rdfentry_")
               .Filter("nevent%1 == 0")


               #############################################
               ##          Aliases for # in python        ##
               #############################################
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")

               #############################################
               ##MC record to study the Z->bb events types##
               #############################################
               .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               .Define("MC_n",   "int(MC_PDG.size())")
               .Define("MC_M1",  "myUtils::get_MCMother1(Particle,Particle0)")
               .Define("MC_M2",  "myUtils::get_MCMother2(Particle,Particle0)")
               .Define("MC_D1",  "myUtils::get_MCDaughter1(Particle,Particle1)")
               .Define("MC_D2",  "myUtils::get_MCDaughter2(Particle,Particle1)")

               #############################################
               ##   gen b quark and Bs meson info         ##
               #############################################

               .Define("genBottom",   "FCCAnalyses::MCParticle::sel_pdgID(5, true)(Particle)")
               .Define("n_genBottoms",     "FCCAnalyses::MCParticle::get_n(genBottom)")
               .Define("genBottom_px",     "FCCAnalyses::MCParticle::get_px(genBottom)")
               .Define("genBottom_py",     "FCCAnalyses::MCParticle::get_py(genBottom)")
               .Define("genBottom_pz",     "FCCAnalyses::MCParticle::get_pz(genBottom)")
               .Define("genBottom_phi",    "FCCAnalyses::MCParticle::get_phi(genBottom)")
               .Define("genBottom_eta",    "FCCAnalyses::MCParticle::get_eta(genBottom)")
               .Define("genBottom_energy", "FCCAnalyses::MCParticle::get_e(genBottom)")
               .Define("genBottom_mass",   "FCCAnalyses::MCParticle::get_mass(genBottom)")
               .Define("genBottom_pdg",    "FCCAnalyses::MCParticle::get_pdg(genBottom)")

               .Define("genBs",   "FCCAnalyses::MCParticle::sel_pdgID(531, true)(Particle)")
               .Define("n_genBs",     "FCCAnalyses::MCParticle::get_n   (genBs)")
               .Define("genBs_px",    "FCCAnalyses::MCParticle::get_px  (genBs)")
               .Define("genBs_py",    "FCCAnalyses::MCParticle::get_py  (genBs)")
               .Define("genBs_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs)")
               .Define("genBs_phi",   "FCCAnalyses::MCParticle::get_phi (genBs)")
               .Define("genBs_eta",   "FCCAnalyses::MCParticle::get_eta (genBs)")
               .Define("genBs_energy","FCCAnalyses::MCParticle::get_e   (genBs)")
               .Define("genBs_mass",  "FCCAnalyses::MCParticle::get_mass(genBs)")
               .Define("genBs_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs)")


               #############################################
               ##               Build MC Vertex           ##
               #############################################
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")

                ############################################
                ##           single out Bs vertices       ##
                ############################################

               .Define("MC_Vertex_isBs",    "ROOT::VecOps::RVec<int> result; for (size_t i=0; i < MC_Vertex_PDGmother.size(); ++i) {int isBs=0; for (size_t j=0; j < MC_Vertex_PDGmother[i].size(); ++j) {if (abs(MC_Vertex_PDGmother[i][j])==531) isBs+=1;} result.push_back(isBs);} return result;")
               .Define("genBs_Vertex_x",   "MC_Vertex_x  [MC_Vertex_isBs>0]") #find Bs meson
               .Define("genBs_Vertex_y",   "MC_Vertex_y  [MC_Vertex_isBs>0]") #no implementation of abs() for vector
               .Define("genBs_Vertex_z",   "MC_Vertex_z  [MC_Vertex_isBs>0]")
            

               #############################################
               ##            Find genBs2TauTau            ##
               #############################################
               
               .Define("genBs2TauTau","FCCAnalyses::ZHfunctions::Find_genBs2TauTau(Particle,Particle1)")
               .Define("n_genBs2TauTau",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau)")
               .Define("genBs2TauTau_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau)")
               .Define("genBs2TauTau_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau)")
               .Define("genBs2TauTau_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau)")
               .Define("genBs2TauTau_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau)")
               .Define("genBs2TauTau_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau)")
               .Define("genBs2TauTau_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau)")
               .Define("genBs2TauTau_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau)")
               .Define("genBs2TauTau_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau)")
               .Define("genBs2TauTau_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau)")
               .Define("genBs2TauTau_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau)")
               .Define("genBs2TauTau_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau)")

               #############################################
               ##            Find genTau2MuNuNu           ##
               #############################################
               
               .Define("genBs2TauTau_Tauplus_muon","FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Muons(genBs2TauTau,Particle,Particle1,-1)")
               .Define("n_genBs2TauTau_Tauplus_muon",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauplus_muon)")
               .Define("genBs2TauTau_Tauplus_muon_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauplus_muon)")

               .Define("genBs2TauTau_Tauminus_muon","FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Muons(genBs2TauTau,Particle,Particle1,1)")
               .Define("n_genBs2TauTau_Tauminus_muon",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauminus_muon)")
               .Define("genBs2TauTau_Tauminus_muon_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauminus_muon)")
               
               #############################################
               ##      Select Bs->Tau(->mu)Tau(->mu)      ##
               #############################################
               #Only for signal to properly select the correct events
               #.Filter("n_genBs2TauTau_Tauminus_muon > 0 && n_genBs2TauTau_Tauplus_muon > 0")

               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")

               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",    "myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "float(myUtils::get_PV_ntracks(VertexObject))")
               .Define("EVT_NVertex",   "float(VertexObject.size())")
               .Filter("EVT_hasPV==1")

               #############################################
               ##         Full 3D missing energy          ##
               #############################################

               .Define("missingEnergy", "FCCAnalyses::ZHfunctions::missingEnergy(91.188, ReconstructedParticles)") ## 91.188 GeV total energy in pythia cards
               .Define("recoEmiss_px",  "missingEnergy[0].momentum.x")
               .Define("recoEmiss_py",  "missingEnergy[0].momentum.y")
               .Define("recoEmiss_pz",  "missingEnergy[0].momentum.z")
               .Define("recoEmiss_e",   "missingEnergy[0].energy")
               .Define("recoEmiss_p4",  "TLorentzVector(recoEmiss_px, recoEmiss_py, recoEmiss_pz, recoEmiss_e)")
               .Define("recoEmiss_m",   "recoEmiss_p4.M()")

               #############################################
               ##          Build RECO P with PID          ##
               #############################################
               .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")

               #############################################
               ##    Build RECO P with PID at vertex      ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##         Build vertex variables          ##
               #############################################
               .Define("Vertex_x",        "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")

               .Define("Vertex_d2PV",     "myUtils::get_Vertex_d2PV(VertexObject,-1)")
               .Define("Vertex_d2PVx",    "myUtils::get_Vertex_d2PV(VertexObject,0)")
               .Define("Vertex_d2PVy",    "myUtils::get_Vertex_d2PV(VertexObject,1)")
               .Define("Vertex_d2PVz",    "myUtils::get_Vertex_d2PV(VertexObject,2)")

               .Define("Vertex_d2PVErr",  "myUtils::get_Vertex_d2PVError(VertexObject,-1)")
               .Define("Vertex_d2PVxErr", "myUtils::get_Vertex_d2PVError(VertexObject,0)")
               .Define("Vertex_d2PVyErr", "myUtils::get_Vertex_d2PVError(VertexObject,1)")
               .Define("Vertex_d2PVzErr", "myUtils::get_Vertex_d2PVError(VertexObject,2)")

               .Define("Vertex_d2PVSig",  "Vertex_d2PV/Vertex_d2PVErr")
               .Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
               .Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
               .Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

               .Define("Vertex_d2MC",     "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,-1)")
               .Define("Vertex_d2MCx",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,0)")
               .Define("Vertex_d2MCy",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,1)")
               .Define("Vertex_d2MCz",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,2)")

               .Define("EVT_dPV2DVmin",   "myUtils::get_dPV2DV_min(Vertex_d2PV)")
               .Define("EVT_dPV2DVmax",   "myUtils::get_dPV2DV_max(Vertex_d2PV)")
               .Define("EVT_dPV2DVave",   "myUtils::get_dPV2DV_ave(Vertex_d2PV)")


               #############################################
               ##              Build the thrust           ##
               #############################################
               .Define("RP_e",          "ReconstructedParticle::get_e(RecoPartPIDAtVertex)")
               .Define("RP_px",         "ReconstructedParticle::get_px(RecoPartPIDAtVertex)")
               .Define("RP_py",         "ReconstructedParticle::get_py(RecoPartPIDAtVertex)")
               .Define("RP_pz",         "ReconstructedParticle::get_pz(RecoPartPIDAtVertex)")
               .Define("RP_charge",     "ReconstructedParticle::get_charge(RecoPartPIDAtVertex)")

               .Define("EVT_thrustNP",      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define("RP_thrustangleNP",  'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
               .Define("EVT_thrust",        'Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e, EVT_thrustNP)')
               .Define("RP_thrustangle",    'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')


               #############################################
               ##      thrust angle of gen b, Bs          ##
               #############################################

               .Define("genBottom_thrustangle",    'Algorithms::getAxisCosTheta(EVT_thrust, genBottom_px, genBottom_py, genBottom_pz)')
               .Define("genBs_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genBs_px, genBs_py, genBs_pz)')
               .Define("genBs_Vertex_thrustangle", 'Algorithms::getAxisCosTheta(EVT_thrust, genBs_Vertex_x, genBs_Vertex_y, genBs_Vertex_z)')
               .Define("recoEmiss_thrustangle",    'Algorithms::getAxisCosTheta(EVT_thrust, recoEmiss_px, recoEmiss_py, recoEmiss_pz)')

               #############################################
               ##        Get thrust related values        ##
               #############################################
               ##hemis0 == negative angle == max energy hemisphere if pointing
               ##hemis1 == positive angle == min energy hemisphere if pointing
               .Define("EVT_thrusthemis0_n",    "Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis1_n",    "Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis0_e",    "Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)")
               .Define("EVT_thrusthemis1_e",    "Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)")
               .Define("EVT_thrusthemis0_p",    "Algorithms::getAxisMomentum(0)(RP_thrustangle, RP_px, RP_py, RP_pz)")
               .Define("EVT_thrusthemis1_p",    "Algorithms::getAxisMomentum(1)(RP_thrustangle, RP_px, RP_py, RP_pz)")

               .Define("EVT_ThrustEmax_E",         "EVT_thrusthemis0_e.at(0)")
               .Define("EVT_ThrustEmax_Echarged",  "EVT_thrusthemis0_e.at(1)")
               .Define("EVT_ThrustEmax_Eneutral",  "EVT_thrusthemis0_e.at(2)")
               .Define("EVT_ThrustEmax_N",         "float(EVT_thrusthemis0_n.at(0))")
               .Define("EVT_ThrustEmax_Ncharged",  "float(EVT_thrusthemis0_n.at(1))")
               .Define("EVT_ThrustEmax_Nneutral",  "float(EVT_thrusthemis0_n.at(2))")
               .Define("EVT_ThrustEmax_px",        "EVT_thrusthemis0_p.at(0)")
               .Define("EVT_ThrustEmax_py",        "EVT_thrusthemis0_p.at(1)")
               .Define("EVT_ThrustEmax_pz",        "EVT_thrusthemis0_p.at(2)")

               .Define("EVT_ThrustEmin_E",         "EVT_thrusthemis1_e.at(0)")
               .Define("EVT_ThrustEmin_Echarged",  "EVT_thrusthemis1_e.at(1)")
               .Define("EVT_ThrustEmin_Eneutral",  "EVT_thrusthemis1_e.at(2)")
               .Define("EVT_ThrustEmin_N",         "float(EVT_thrusthemis1_n.at(0))")
               .Define("EVT_ThrustEmin_Ncharged",  "float(EVT_thrusthemis1_n.at(1))")
               .Define("EVT_ThrustEmin_Nneutral",  "float(EVT_thrusthemis1_n.at(2))")
               .Define("EVT_ThrustEmin_px",        "EVT_thrusthemis1_p.at(0)")
               .Define("EVT_ThrustEmin_py",        "EVT_thrusthemis1_p.at(1)")
               .Define("EVT_ThrustEmin_pz",        "EVT_thrusthemis1_p.at(2)")


               .Define("Vertex_thrust_angle",   "myUtils::get_Vertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               .Define("DVertex_thrust_angle",  "myUtils::get_DVertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               ###0 == negative angle==max energy , 1 == positive angle == min energy
               .Define("Vertex_thrusthemis_emin",    "myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 1)")
               .Define("Vertex_thrusthemis_emax",    "myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 0)")

               .Define("EVT_ThrustEmin_NDV", "float(myUtils::get_Npos(DVertex_thrust_angle))")
               .Define("EVT_ThrustEmax_NDV", "float(myUtils::get_Nneg(DVertex_thrust_angle))")

               .Define("EVT_Thrust_Mag",  "EVT_thrust.at(0)")
               .Define("EVT_Thrust_X",    "EVT_thrust.at(1)")
               .Define("EVT_Thrust_XErr", "EVT_thrust.at(2)")
               .Define("EVT_Thrust_Y",    "EVT_thrust.at(3)")
               .Define("EVT_Thrust_YErr", "EVT_thrust.at(4)")
               .Define("EVT_Thrust_Z",    "EVT_thrust.at(5)")
               .Define("EVT_Thrust_ZErr", "EVT_thrust.at(6)")


               .Define("DV_tracks", "myUtils::get_pseudotrack(VertexObject,RecoPartPIDAtVertex)")

               .Define("DV_d0",            "myUtils::get_trackd0(DV_tracks)")
               .Define("DV_z0",            "myUtils::get_trackz0(DV_tracks)")


               ##############################
               ##         Get muons        ##
               ##############################

               .Alias("Muon0","Muon#0.index")
               .Define("muons",            "ReconstructedParticle::get(Muon0, ReconstructedParticles)")
               .Define("n_muons",          "ReconstructedParticle::get_n(muons)")
               .Define("muon_px",          "ReconstructedParticle::get_px(muons)")
               .Define("muon_py",          "ReconstructedParticle::get_py(muons)")
               .Define("muon_pz",          "ReconstructedParticle::get_pz(muons)")
               .Define("muon_phi",         "ReconstructedParticle::get_phi(muons)")
               .Define("muon_eta",         "ReconstructedParticle::get_eta(muons)")
               .Define("muon_energy",      "ReconstructedParticle::get_e(muons)")
               .Define("muon_mass",        "ReconstructedParticle::get_mass(muons)")
               .Define("muon_charge",      "ReconstructedParticle::get_charge(muons)")
               .Define("muon_PDG",         "ReconstructedParticle::get_type(muons)")
               .Define("muon_thrustangles",'Algorithms::getAxisCosTheta(EVT_thrust, muon_px, muon_py, muon_pz)')
               .Define("muon_thrustEmin_n","int result (0); for (size_t i=0; i<muon_thrustangles.size(); ++i){if (muon_thrustangles[i] > 0.0) ++result;} return result;")
               .Define("muon_thrustEmax_n","int result (0); for (size_t i=0; i<muon_thrustangles.size(); ++i){if (muon_thrustangles[i] <= 0.0) ++result;} return result;")
               .Define("muon_OpeningAngle","FCCAnalyses::ZHfunctions::Muons_ComputeOpeningAngle(muon_px,muon_py,muon_pz)")
               .Define("dimuon_ind", "FCCAnalyses::ZHfunctions::Muons_ID_SmallestOA(muon_OpeningAngle,n_muons)")

               #.Define("nevent","rdfentry_")
               #.Define("Selected","if (n_muons > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dimuon > 0 && has_dimuon_SameSide > 0 && has_dimuon_OppositeCharges > 0 && has_dimuon_SigHemi) return 1; else return 0;")

               
               #####################################
               ##   muon Perform Truth-Matching   ##
               #####################################

               .Define("TM_MC_muplus_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_Mothers(MCRecoAssociations1,Particle,Particle0,15,531,-13)")
               .Define("TM_MC_muminus_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_Mothers(MCRecoAssociations1,Particle,Particle0,15,531,13)")

               .Define("TM_RECO_muplus_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO_Subset(MCRecoAssociations0,Muon0,TM_MC_muplus_ind,MCRecoAssociations1)")
               .Define("TM_RECO_muminus_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO_Subset(MCRecoAssociations0,Muon0,TM_MC_muminus_ind,MCRecoAssociations1)")

               .Define("n_TM_muplus","ReconstructedParticle::get_n(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("n_TM_muminus","ReconstructedParticle::get_n(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("n_TM_muons","n_TM_muplus+n_TM_muminus")

               .Define("TM_muplus_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(TM_RECO_muplus_ind,ReconstructedParticles))")
               .Define("TM_muplus_thrustangle","Algorithms::getAxisCosTheta(EVT_thrust, TM_muplus_px, TM_muplus_py, TM_muplus_pz)")
               
               .Define("TM_muminus_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(TM_RECO_muminus_ind,ReconstructedParticles))")
               .Define("TM_muminus_thrustangle","Algorithms::getAxisCosTheta(EVT_thrust, TM_muminus_px, TM_muminus_py, TM_muminus_pz)")

               .Define("TM_muons_OpeningAngle","FCCAnalyses::ZHfunctions::TM_ComputeOpeningAngle(TM_muplus_px,TM_muplus_py,TM_muplus_pz,TM_muminus_px,TM_muminus_py,TM_muminus_pz)")

               #.Define("RP_MC_index", "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles)")

               #.Define("TruthMatched_muplus", "FCCAnalyses::ZHfunctions::Muon_TruthMatching(ReconstructedParticles,MCRecoAssociations0,MCRecoAssociations1,Particle,Particle0,1,Muon0)")
               #.Define("n_TruthMatched_muplus",           "ReconstructedParticle::get_n(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_px",          "ReconstructedParticle::get_px(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_py",          "ReconstructedParticle::get_py(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_pz",          "ReconstructedParticle::get_pz(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_phi",         "ReconstructedParticle::get_phi(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_eta",         "ReconstructedParticle::get_eta(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_energy",      "ReconstructedParticle::get_e(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_mass",        "ReconstructedParticle::get_mass(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_charge",      "ReconstructedParticle::get_charge(TruthMatched_muplus)")
               #.Define("TruthMatched_muplus_PDG",         "ReconstructedParticle::get_type(TruthMatched_muplus)")

               #.Define("TruthMatched_muminus", "FCCAnalyses::ZHfunctions::Muon_TruthMatching(ReconstructedParticles,MCRecoAssociations0,MCRecoAssociations1,Particle,Particle0,-1,Muon0)")
               #.Define("n_TruthMatched_muminus",           "ReconstructedParticle::get_n(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_px",          "ReconstructedParticle::get_px(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_py",          "ReconstructedParticle::get_py(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_pz",          "ReconstructedParticle::get_pz(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_phi",         "ReconstructedParticle::get_phi(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_eta",         "ReconstructedParticle::get_eta(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_energy",      "ReconstructedParticle::get_e(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_mass",        "ReconstructedParticle::get_mass(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_charge",      "ReconstructedParticle::get_charge(TruthMatched_muminus)")
               #.Define("TruthMatched_muminus_PDG",         "ReconstructedParticle::get_type(TruthMatched_muminus)")

               #.Define("n_TruthMatched_muons","n_TruthMatched_muplus + n_TruthMatched_muminus")


               ###########################
               ##     Stage 1 cuts      ##
               ###########################
               
               .Define("has_dimuon","FCCAnalyses::ZHfunctions::Check_dimuon_Presence(muon_OpeningAngle)")
               .Define("has_dimuon_SameSide","FCCAnalyses::ZHfunctions::Check_dimuon_SameSide(muon_OpeningAngle,has_dimuon)")
               .Define("has_dimuon_SigHemi","FCCAnalyses::ZHfunctions::Check_dimuon_SigHemi(dimuon_ind, muon_thrustangles,has_dimuon)")
               .Define("has_dimuon_OppositeCharges","FCCAnalyses::ZHfunctions::Check_dimuon_Charges(dimuon_ind, muon_charge,has_dimuon)")
               .Define("has_dimuon_Vertex","FCCAnalyses::ZHfunctions::Check_dimuon_Vertices(dimuon_ind,Muon0,muon_charge,VertexObject,has_dimuon)")
               .Filter("n_muons > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dimuon > 0 && has_dimuon_SameSide > 0 && has_dimuon_OppositeCharges > 0 && has_dimuon_SigHemi > 0 && has_dimuon_Vertex == 0")

               ############################
               ##   Background studies   ##
               ############################

               .Define("MC_dimuon_ind","FCCAnalyses::ZHfunctions::Finding_MC_dimuon(dimuon_ind,Muon0,MCRecoAssociations0,MCRecoAssociations1)") #This function might not work as intended if the dimuon is made of muon not coming from Bs2TauTau (can happen by chance if 3 muon in sig hemi)
               .Define("n_MC_dimuon","MC_dimuon_ind.size()") #Check to see if the truth-matching in the RECO->MC direction can fail (== 2 -> ok else failed)

               .Define("MC_dimuon_CA_ind","FCCAnalyses::ZHfunctions::Find_MC_CommonAncestor(MC_dimuon_ind,Particle,Particle0)")
               .Define("MC_dimuon_CADaughters_ind","FCCAnalyses::ZHfunctions::Find_MC_CommonAncestor_Daughters(MC_dimuon_CA_ind,Particle,Particle1)") #Always check the indices value before looking into the collection (-1 for failed searches)

           )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                "MC_PDG","MC_M1","MC_M2","MC_n","MC_D1","MC_D2",
    
                "n_genBottoms",
                "genBottom_px", "genBottom_py", "genBottom_pz", "genBottom_eta", "genBottom_phi",
                "genBottom_energy", "genBottom_mass", "genBottom_pdg",
                "genBottom_thrustangle",

                "n_genBs",
                "genBs_px", "genBs_py", "genBs_pz", "genBs_eta", "genBs_phi",
                "genBs_energy", "genBs_mass", "genBs_pdg",
                "genBs_thrustangle","genBs_Vertex_x", "genBs_Vertex_y", "genBs_Vertex_z",
                "genBs_Vertex_thrustangle",

                "n_genBs2TauTau","genBs2TauTau_px","genBs2TauTau_py","genBs2TauTau_pz",
                "genBs2TauTau_phi","genBs2TauTau_eta","genBs2TauTau_energy","genBs2TauTau_mass","genBs2TauTau_pdg",
                "genBs2TauTau_Vertex_x","genBs2TauTau_Vertex_y","genBs2TauTau_Vertex_z",

                "n_genBs2TauTau_Tauplus_muon","genBs2TauTau_Tauplus_muon_px","genBs2TauTau_Tauplus_muon_py","genBs2TauTau_Tauplus_muon_pz",
                "genBs2TauTau_Tauplus_muon_phi","genBs2TauTau_Tauplus_muon_eta","genBs2TauTau_Tauplus_muon_energy","genBs2TauTau_Tauplus_muon_mass","genBs2TauTau_Tauplus_muon_pdg",
                "genBs2TauTau_Tauplus_muon_Vertex_x","genBs2TauTau_Tauplus_muon_Vertex_y","genBs2TauTau_Tauplus_muon_Vertex_z",

                "n_genBs2TauTau_Tauminus_muon","genBs2TauTau_Tauminus_muon_px","genBs2TauTau_Tauminus_muon_py","genBs2TauTau_Tauminus_muon_pz",
                "genBs2TauTau_Tauminus_muon_phi","genBs2TauTau_Tauminus_muon_eta","genBs2TauTau_Tauminus_muon_energy","genBs2TauTau_Tauminus_muon_mass","genBs2TauTau_Tauminus_muon_pdg",
                "genBs2TauTau_Tauminus_muon_Vertex_x","genBs2TauTau_Tauminus_muon_Vertex_y","genBs2TauTau_Tauminus_muon_Vertex_z",

                "EVT_ThrustEmin_E",            "EVT_ThrustEmax_E",
                "EVT_ThrustEmin_Echarged",     "EVT_ThrustEmax_Echarged",
                "EVT_ThrustEmin_Eneutral",     "EVT_ThrustEmax_Eneutral",
                "EVT_ThrustEmin_N",            "EVT_ThrustEmax_N",
                "EVT_ThrustEmin_Ncharged",     "EVT_ThrustEmax_Ncharged",
                "EVT_ThrustEmin_Nneutral",     "EVT_ThrustEmax_Nneutral",
                "EVT_ThrustEmin_px",           "EVT_ThrustEmax_px",
                "EVT_ThrustEmin_py",           "EVT_ThrustEmax_py",
                "EVT_ThrustEmin_pz",           "EVT_ThrustEmax_pz",
                "EVT_ThrustEmin_NDV",          "EVT_ThrustEmax_NDV",
                "EVT_Thrust_Mag",
                "EVT_Thrust_X",  "EVT_Thrust_XErr",
                "EVT_Thrust_Y",  "EVT_Thrust_YErr",
                "EVT_Thrust_Z",  "EVT_Thrust_ZErr",

                "EVT_NtracksPV", "EVT_NVertex",

                "EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",

                "MC_Vertex_x", "MC_Vertex_y", "MC_Vertex_z",
                "MC_Vertex_ntrk", "MC_Vertex_n",

                "MC_Vertex_PDG","MC_Vertex_PDGmother","MC_Vertex_PDGgmother",

                "Vertex_x", "Vertex_y", "Vertex_z",
                "Vertex_xErr", "Vertex_yErr", "Vertex_zErr",
                "Vertex_isPV", "Vertex_ntrk", "Vertex_chi2", "Vertex_n",
                "Vertex_thrust_angle", "Vertex_thrusthemis_emin", "Vertex_thrusthemis_emax",

                "Vertex_d2PV", "Vertex_d2PVx", "Vertex_d2PVy", "Vertex_d2PVz",
                "Vertex_d2PVErr", "Vertex_d2PVxErr", "Vertex_d2PVyErr", "Vertex_d2PVzErr",
                "Vertex_mass",
                "DV_d0","DV_z0",

                "recoEmiss_px", "recoEmiss_py", "recoEmiss_pz", "recoEmiss_e", "recoEmiss_m",
                "recoEmiss_thrustangle",

                "n_muons","muon_px","muon_py","muon_pz","muon_phi","muon_eta","muon_energy","muon_mass","muon_charge","muon_PDG","muon_thrustangles","muon_thrustEmin_n","muon_thrustEmax_n","muon_OpeningAngle","dimuon_ind",
                "has_dimuon","has_dimuon_SameSide","has_dimuon_SigHemi","has_dimuon_OppositeCharges","has_dimuon_Vertex",

                "TM_MC_muplus_ind","TM_MC_muminus_ind","TM_RECO_muplus_ind","TM_RECO_muminus_ind","n_TM_muplus","n_TM_muminus","n_TM_muons",
                "TM_muplus_energy","TM_muplus_px","TM_muplus_py","TM_muplus_pz","TM_muplus_phi","TM_muplus_eta","TM_muplus_mass","TM_muplus_charge","TM_muplus_PDG","TM_muplus_thrustangle",
                "TM_muminus_energy","TM_muminus_px","TM_muminus_py","TM_muminus_pz","TM_muminus_phi","TM_muminus_eta","TM_muminus_mass","TM_muminus_charge","TM_muminus_PDG","TM_muminus_thrustangle",
                "TM_muons_OpeningAngle",

                "MC_dimuon_ind","n_MC_dimuon","MC_dimuon_CA_ind","MC_dimuon_CADaughters_ind",

                #"Selected","nevent"
                #"n_TruthMatched_muplus","TruthMatched_muplus_px","TruthMatched_muplus_py","TruthMatched_muplus_pz","TruthMatched_muplus_phi","TruthMatched_muplus_eta","TruthMatched_muplus_energy","TruthMatched_muplus_mass","TruthMatched_muplus_charge","TruthMatched_muplus_PDG",
                #"n_TruthMatched_muminus","TruthMatched_muminus_px","TruthMatched_muminus_py","TruthMatched_muminus_pz","TruthMatched_muminus_phi","TruthMatched_muminus_eta","TruthMatched_muminus_energy","TruthMatched_muminus_mass","TruthMatched_muminus_charge","TruthMatched_muminus_PDG",
                #"n_TruthMatched_muons",             
                ]
        return branchList
