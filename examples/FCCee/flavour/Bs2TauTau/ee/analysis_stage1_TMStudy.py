#analysis_stage1

# list of samples to process
processList_full = {
    #'p8_ee_Zbb_ecm91':{'chunks':10,'fraction':0.01},
    #'p8_ee_Zcc_ecm91':{'chunks':10,'fraction':0.01},
    #'p8_ee_Zss_ecm91':{'chunks':10,'fraction':0.01},
    #'p8_ee_Zud_ecm91':{'chunks':10,'fraction':0.01},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':10},
    'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':20,'fraction':1},
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
    #'p8_ee_Zbb_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':1, 'fraction':0.000002},
    'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':1, 'fraction':0.01},
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
outputDirEos   = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_ee_noFilter"

# if runBatch = False, save output locally
outputDir   = "DummyRepo"
#"fccanalysis_output/" To be used with runBatch = True

includePaths = ["functions.h"]

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df

               #.Define("nevent","rdfentry_")
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
               
               .Define("genBs2TauTau_Tauplus_electron","FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Electrons(genBs2TauTau,Particle,Particle1,-1)")
               .Define("n_genBs2TauTau_Tauplus_electron",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauplus_electron)")
               .Define("genBs2TauTau_Tauplus_electron_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauplus_electron)")

               .Define("genBs2TauTau_Tauminus_electron","FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Electrons(genBs2TauTau,Particle,Particle1,1)")
               .Define("n_genBs2TauTau_Tauminus_electron",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauminus_electron)")
               .Define("genBs2TauTau_Tauminus_electron_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauminus_electron)")
               
               #############################################
               ##      Select Bs->Tau(->mu)Tau(->mu)      ##
               #############################################
               #Only for signal to properly select the correct events
               #.Define("Selected","if (n_genBs2TauTau_Tauminus_electron > 0 && n_genBs2TauTau_Tauplus_electron > 0) return 1; else return 0;")
               .Filter("n_genBs2TauTau_Tauminus_electron > 0 && n_genBs2TauTau_Tauplus_electron > 0") 

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
               ##       Get electrons      ##
               ##############################

               .Alias("Electron0","Electron#0.index")
               .Define("electrons",            "ReconstructedParticle::get(Electron0, ReconstructedParticles)")
               .Define("n_electrons",          "ReconstructedParticle::get_n(electrons)")
               .Define("electron_px",          "ReconstructedParticle::get_px(electrons)")
               .Define("electron_py",          "ReconstructedParticle::get_py(electrons)")
               .Define("electron_pz",          "ReconstructedParticle::get_pz(electrons)")
               .Define("electron_phi",         "ReconstructedParticle::get_phi(electrons)")
               .Define("electron_eta",         "ReconstructedParticle::get_eta(electrons)")
               .Define("electron_energy",      "ReconstructedParticle::get_e(electrons)")
               .Define("electron_mass",        "ReconstructedParticle::get_mass(electrons)")
               .Define("electron_charge",      "ReconstructedParticle::get_charge(electrons)")
               .Define("electron_PDG",         "ReconstructedParticle::get_type(electrons)")
               .Define("electron_thrustangles",'Algorithms::getAxisCosTheta(EVT_thrust, electron_px, electron_py, electron_pz)')
               .Define("electron_thrustEmin_n","int result (0); for (size_t i=0; i<electron_thrustangles.size(); ++i){if (electron_thrustangles[i] > 0.0) ++result;} return result;")
               .Define("electron_thrustEmax_n","int result (0); for (size_t i=0; i<electron_thrustangles.size(); ++i){if (electron_thrustangles[i] <= 0.0) ++result;} return result;")
               .Define("electron_OpeningAngle","FCCAnalyses::ZHfunctions::Muons_ComputeOpeningAngle(electron_px,electron_py,electron_pz)")
               .Define("has_dielectron","FCCAnalyses::ZHfunctions::Check_dimuon_Presence(electron_OpeningAngle)")
               
               #subslect the two muons with smallest opening angle (correctly selects the muon from Bs2TauTau) and check their properties
               .Define("dielectron_ind", "FCCAnalyses::ZHfunctions::Muons_ID_SmallestOA(electron_OpeningAngle,n_electrons)")
               .Define("has_dielectron_SameSide","FCCAnalyses::ZHfunctions::Check_dimuon_SameSide(electron_OpeningAngle,has_dielectron)")
               .Define("has_dielectron_SigHemi","FCCAnalyses::ZHfunctions::Check_dimuon_SigHemi(dielectron_ind, electron_thrustangles,has_dielectron)")
               .Define("has_dielectron_OppositeCharges","FCCAnalyses::ZHfunctions::Check_dimuon_Charges(dielectron_ind, electron_charge,has_dielectron)")
               .Define("has_dielectron_Vertex","FCCAnalyses::ZHfunctions::Check_dimuon_Vertices(dielectron_ind,Electron0,electron_charge,VertexObject,has_dielectron)")

               
               #####################################
               ## electron Perform Truth-Matching ##
               #####################################

               .Define("TM_MC_all_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_Mothers(MCRecoAssociations1,Particle,Particle0,15,531,0)")
               
               .Define("TM_MC_positron_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_Mothers(MCRecoAssociations1,Particle,Particle0,15,531,-11)")
               .Define("TM_MC_electron_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_Mothers(MCRecoAssociations1,Particle,Particle0,15,531,11)")

               .Define("TM_RECO_positron_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO_Subset(MCRecoAssociations0,Electron0,TM_MC_positron_ind,MCRecoAssociations1)")
               .Define("TM_RECO_electron_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO_Subset(MCRecoAssociations0,Electron0,TM_MC_electron_ind,MCRecoAssociations1)")

               .Define("n_TM_positron","ReconstructedParticle::get_n(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("n_TM_electron","ReconstructedParticle::get_n(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("n_TM_electrons","n_TM_positron+n_TM_electron")

               .Define("TM_positron_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(TM_RECO_positron_ind,ReconstructedParticles))")
               .Define("TM_positron_thrustangle","Algorithms::getAxisCosTheta(EVT_thrust, TM_positron_px, TM_positron_py, TM_positron_pz)")
               
               .Define("TM_electron_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(TM_RECO_electron_ind,ReconstructedParticles))")
               .Define("TM_electron_thrustangle","Algorithms::getAxisCosTheta(EVT_thrust, TM_electron_px, TM_electron_py, TM_electron_pz)")

               .Define("TM_electrons_OpeningAngle","FCCAnalyses::ZHfunctions::TM_ComputeOpeningAngle(TM_positron_px,TM_positron_py,TM_positron_pz,TM_electron_px,TM_electron_py,TM_electron_pz)")

               ############################
               ##   TM studies (dupli)   ##
               ############################

               #Filter before checking the MC decay tree
               .Filter("n_electrons > 1 && EVT_ThrustEmin_E < 38 && recoEmiss_e > 10 && EVT_ThrustEmin_Eneutral < 10 && has_dielectron > 0 && has_dielectron_SameSide > 0 && has_dielectron_OppositeCharges > 0 && has_dielectron_SigHemi > 0 && has_dielectron_Vertex == 0")
               .Define("MC_dielectron_ind","FCCAnalyses::ZHfunctions::Finding_MC_dimuon(dielectron_ind,Electron0,MCRecoAssociations0,MCRecoAssociations1)")
               .Define("n_MC_dielectron","MC_dielectron_ind.size()") #Check to see if the truth-matching in the RECO->MC direction can fail (== 2 -> ok else failed)

               .Define("MC_dielectron_CA_ind","FCCAnalyses::ZHfunctions::Find_MC_CommonAncestor(MC_dielectron_ind,Particle,Particle0)")
               .Define("MC_dielectron_CADaughters_ind","FCCAnalyses::ZHfunctions::Find_MC_CommonAncestor_Daughters(MC_dielectron_CA_ind,Particle,Particle1)") #Always check the indices value before looking into the collection (-1 for failed searches)

               .Define("MCRecoAsso","ROOT::VecOps::RVec<int> result; for(size_t i=0; i<MCRecoAssociations1.size(); ++i){result.push_back(MCRecoAssociations1.at(i));} return result;")
               .Define("MCRecoAssoo","ROOT::VecOps::RVec<int> result; for(size_t i=0; i<MCRecoAssociations0.size(); ++i){result.push_back(MCRecoAssociations0.at(i));} return result;")
               .Define("Electroo","ROOT::VecOps::RVec<int> result; for(size_t i=0; i<Electron0.size(); ++i){result.push_back(Electron0.at(i));} return result;")
               .Alias("Photon0","Photon#0.index")
               .Define("Photoo","ROOT::VecOps::RVec<int> result; for(size_t i=0; i<Photon0.size(); ++i){result.push_back(Photon0.at(i));} return result;")

               
               .Define("TM_RECO_electronwithDuplicates_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO(MCRecoAssociations0,TM_MC_electron_ind,MCRecoAssociations1)")
               .Define("TM_RECO_positronwithDuplicates_ind","FCCAnalyses::ZHfunctions::TruthMatch_MC_RECO(MCRecoAssociations0,TM_MC_positron_ind,MCRecoAssociations1)")
               
               .Define("TM_positronwithDuplicates_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_positronwithDuplicates_ind,ReconstructedParticles))")
               .Define("TM_electronwithDuplicates_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(TM_RECO_electronwithDuplicates_ind,ReconstructedParticles))")
               .Define("TM_positronwithDuplicates_charge",   "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_positronwithDuplicates_ind,ReconstructedParticles))")
               .Define("TM_electronwithDuplicates_charge",   "ReconstructedParticle::get_charge(ReconstructedParticle::get(TM_RECO_electronwithDuplicates_ind,ReconstructedParticles))")


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

                "n_genBs2TauTau_Tauplus_electron","genBs2TauTau_Tauplus_electron_px","genBs2TauTau_Tauplus_electron_py","genBs2TauTau_Tauplus_electron_pz",
                "genBs2TauTau_Tauplus_electron_phi","genBs2TauTau_Tauplus_electron_eta","genBs2TauTau_Tauplus_electron_energy","genBs2TauTau_Tauplus_electron_mass","genBs2TauTau_Tauplus_electron_pdg",
                "genBs2TauTau_Tauplus_electron_Vertex_x","genBs2TauTau_Tauplus_electron_Vertex_y","genBs2TauTau_Tauplus_electron_Vertex_z",

                "n_genBs2TauTau_Tauminus_electron","genBs2TauTau_Tauminus_electron_px","genBs2TauTau_Tauminus_electron_py","genBs2TauTau_Tauminus_electron_pz",
                "genBs2TauTau_Tauminus_electron_phi","genBs2TauTau_Tauminus_electron_eta","genBs2TauTau_Tauminus_electron_energy","genBs2TauTau_Tauminus_electron_mass","genBs2TauTau_Tauminus_electron_pdg",
                "genBs2TauTau_Tauminus_electron_Vertex_x","genBs2TauTau_Tauminus_electron_Vertex_y","genBs2TauTau_Tauminus_electron_Vertex_z",

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

                "n_electrons","electron_px","electron_py","electron_pz","electron_phi","electron_eta","electron_energy","electron_mass","electron_charge","electron_PDG","electron_thrustangles","electron_thrustEmin_n","electron_thrustEmax_n","electron_OpeningAngle","dielectron_ind",
                "has_dielectron","has_dielectron_SameSide","has_dielectron_SigHemi","has_dielectron_OppositeCharges","has_dielectron_Vertex",

                "TM_MC_positron_ind","TM_MC_electron_ind","TM_RECO_positron_ind","TM_RECO_electron_ind","n_TM_positron","n_TM_electron","n_TM_electrons",
                "TM_positron_energy","TM_positron_px","TM_positron_py","TM_positron_pz","TM_positron_phi","TM_positron_eta","TM_positron_mass","TM_positron_charge","TM_positron_PDG","TM_positron_thrustangle",
                "TM_electron_energy","TM_electron_px","TM_electron_py","TM_electron_pz","TM_electron_phi","TM_electron_eta","TM_electron_mass","TM_electron_charge","TM_electron_PDG","TM_electron_thrustangle",
                "TM_electrons_OpeningAngle","TM_MC_all_ind",

                "MC_dielectron_ind","n_MC_dielectron","MC_dielectron_CA_ind","MC_dielectron_CADaughters_ind",

                "MCRecoAsso","MCRecoAssoo","Electroo","Photoo",
                "TM_RECO_electronwithDuplicates_ind","TM_RECO_positronwithDuplicates_ind","TM_electronwithDuplicates_mass","TM_positronwithDuplicates_mass","TM_positronwithDuplicates_charge","TM_electronwithDuplicates_charge"

                #"nevent","Selected",
                #"n_TruthMatched_muplus","TruthMatched_muplus_px","TruthMatched_muplus_py","TruthMatched_muplus_pz","TruthMatched_muplus_phi","TruthMatched_muplus_eta","TruthMatched_muplus_energy","TruthMatched_muplus_mass","TruthMatched_muplus_charge","TruthMatched_muplus_PDG",
                #"n_TruthMatched_muminus","TruthMatched_muminus_px","TruthMatched_muminus_py","TruthMatched_muminus_pz","TruthMatched_muminus_phi","TruthMatched_muminus_eta","TruthMatched_muminus_energy","TruthMatched_muminus_mass","TruthMatched_muminus_charge","TruthMatched_muminus_PDG",
                #"n_TruthMatched_muons",             
                ]
        return branchList
