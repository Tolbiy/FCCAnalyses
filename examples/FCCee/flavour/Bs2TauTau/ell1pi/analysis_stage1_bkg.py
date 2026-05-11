#analysis_stage1

# list of samples to process
processList_full = {
    'p8_ee_Zbb_ecm91':{'chunks':10,'fraction':2e-3},
    'p8_ee_Zcc_ecm91':{'chunks':10,'fraction':2e-3},
    'p8_ee_Zss_ecm91':{'chunks':10,'fraction':2e-3},
    'p8_ee_Zud_ecm91':{'chunks':10,'fraction':2e-3},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTauTAUHADNU':{'chunks':10},
    #'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':10,'fraction':0.1},
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
    'p8_ee_Zbb_ecm91_EvtGen_Bs2TauTau':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zss_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zud_ecm91':{'chunks':1, 'fraction':0.000002},
    #'p8_ee_Zcc_ecm91':{'chunks':1, 'fraction':0.000002},
}

nCPUS       = 8
runBatch    = True
batchQueue  = "nextweek"
compGroup   = "group_u_FCC.local_gen"

processList  = processList_full
if not runBatch:
    processList  = processList_test

# tag for MC production campaign
prodTag     = "FCCee/winter2023/IDEA/"

# if runBatch = True, save output on eos
outputDirEos   = "/eos/experiment/fcc/ee/analyses_storage/flavor/Bs2TauTau/flatNtuples/winter2023/analysis_stage1_ell1pi"

# if runBatch = False, save output locally
outputDir   = "./" #"../DummyRepo"
#"fccanalysis_output/" To be used with runBatch = True

includePaths = ["../functions.h"]


#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
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

               .Define("genBc",   "FCCAnalyses::MCParticle::sel_pdgID(541, true)(Particle)")
               .Define("n_genBc",     "FCCAnalyses::MCParticle::get_n   (genBc)")
               .Define("genBc_px",    "FCCAnalyses::MCParticle::get_px  (genBc)")
               .Define("genBc_py",    "FCCAnalyses::MCParticle::get_py  (genBc)")
               .Define("genBc_pz",    "FCCAnalyses::MCParticle::get_pz  (genBc)")

               .Define("genBu",   "FCCAnalyses::MCParticle::sel_pdgID(521, true)(Particle)")
               .Define("n_genBu",     "FCCAnalyses::MCParticle::get_n   (genBu)")
               .Define("genBu_px",    "FCCAnalyses::MCParticle::get_px  (genBu)")
               .Define("genBu_py",    "FCCAnalyses::MCParticle::get_py  (genBu)")
               .Define("genBu_pz",    "FCCAnalyses::MCParticle::get_pz  (genBu)")

               .Define("genBd",   "FCCAnalyses::MCParticle::sel_pdgID(511, true)(Particle)")
               .Define("n_genBd",     "FCCAnalyses::MCParticle::get_n   (genBd)")
               .Define("genBd_px",    "FCCAnalyses::MCParticle::get_px  (genBd)")
               .Define("genBd_py",    "FCCAnalyses::MCParticle::get_py  (genBd)")
               .Define("genBd_pz",    "FCCAnalyses::MCParticle::get_pz  (genBd)")

               .Define("genLb",   "FCCAnalyses::MCParticle::sel_pdgID(5122, true)(Particle)")
               .Define("n_genLb",     "FCCAnalyses::MCParticle::get_n   (genLb)")
               .Define("genLb_px",    "FCCAnalyses::MCParticle::get_px  (genLb)")
               .Define("genLb_py",    "FCCAnalyses::MCParticle::get_py  (genLb)")
               .Define("genLb_pz",    "FCCAnalyses::MCParticle::get_pz  (genLb)")


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

                ############################################
                ##           Get PV MC truth              ##
                ############################################

               .Define("MC_PV_xyzt",      "FCCAnalyses::MCParticle::get_EventPrimaryVertexP4()(Particle)")

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
               ## Find genTau2Pi(Pi0)Nu & genTau2EllNuNu  ##
               #############################################
               .Define("genBs2TauTau_Tauplus_1pi",  "FCCAnalyses::ZHfunctions::Find_genBs2TauTau_1Pi(genBs2TauTau,Particle,Particle1,-1)")
               .Define("genBs2TauTau_Tauminus_1pi", "FCCAnalyses::ZHfunctions::Find_genBs2TauTau_1Pi(genBs2TauTau,Particle,Particle1,1)")
               .Define("genBs2TauTau_Tauplus_ell", "FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Leptons(genBs2TauTau,Particle,Particle1,-1)")
               .Define("genBs2TauTau_Tauminus_ell","FCCAnalyses::ZHfunctions::Find_genBs2TauTau_Leptons(genBs2TauTau,Particle,Particle1,1)")

               .Define("n_genBs2TauTau_Tauplus_1pi",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauplus_1pi)")
               .Define("genBs2TauTau_Tauplus_1pi_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauplus_1pi)")

               .Define("n_genBs2TauTau_Tauminus_1pi",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauminus_1pi)")
               .Define("genBs2TauTau_Tauminus_1pi_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauminus_1pi)")
               
               .Define("n_genBs2TauTau_Tauplus_ell",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauplus_ell)")
               .Define("genBs2TauTau_Tauplus_ell_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauplus_ell)")

               .Define("n_genBs2TauTau_Tauminus_ell",     "FCCAnalyses::MCParticle::get_n   (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_px",    "FCCAnalyses::MCParticle::get_px  (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_py",    "FCCAnalyses::MCParticle::get_py  (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_pz",    "FCCAnalyses::MCParticle::get_pz  (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_phi",   "FCCAnalyses::MCParticle::get_phi (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_eta",   "FCCAnalyses::MCParticle::get_eta (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_energy","FCCAnalyses::MCParticle::get_e   (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_mass",  "FCCAnalyses::MCParticle::get_mass(genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_pdg",   "FCCAnalyses::MCParticle::get_pdg (genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_Vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_Vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genBs2TauTau_Tauminus_ell)")
               .Define("genBs2TauTau_Tauminus_ell_Vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genBs2TauTau_Tauminus_ell)")
               
               #############################################
               ##      Select Bs->Tau(->ell)Tau(->1pi)    ##
               #############################################
               #IF SIGNAL => FILTER
               #.Filter("(n_genBs2TauTau_Tauminus_1pi > 0 && n_genBs2TauTau_Tauplus_ell > 0) || (n_genBs2TauTau_Tauplus_1pi > 0 && n_genBs2TauTau_Tauminus_ell > 0)") 
               
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
               .Define("Vertex_x_pre",        "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y_pre",        "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z_pre",        "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr_pre",     "myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr_pre",     "myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr_pre",     "myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2_pre",     "myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV_pre",     "myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk_pre",     "myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x_pre.size())")
               .Define("Vertex_mass_pre",     "myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")
               .Define("Vertex_p4",       "FCCAnalyses::ZHfunctions::get_Vertex_p4(VertexObject,RecoPartPIDAtVertex)")
               .Define("Vertex_px_pre",       "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_p4) {result.push_back(p.Px());} return result;")
               .Define("Vertex_py_pre",       "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_p4) {result.push_back(p.Py());} return result;")
               .Define("Vertex_pz_pre",       "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_p4) {result.push_back(p.Pz());} return result;")
               .Define("Vertex_e_pre",        "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_p4) {result.push_back(p.E ());} return result;")
               .Define("Vertex_vec",      "FCCAnalyses::ZHfunctions::build_p4(Vertex_x_pre, Vertex_y_pre, Vertex_z_pre, Vertex_mass_pre)")
               .Define("Vertex_phi",      "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_vec) {result.push_back(p.Phi());} return result;")
               .Define("Vertex_theta",    "ROOT::VecOps::RVec<float> result; for (auto & p: Vertex_vec) {result.push_back(p.Theta());} return result;")

               .Define("Vertex_d2PV_pre",     "myUtils::get_Vertex_d2PV(VertexObject,-1)")
               .Define("Vertex_d2PVx",    "myUtils::get_Vertex_d2PV(VertexObject,0)")
               .Define("Vertex_d2PVy",    "myUtils::get_Vertex_d2PV(VertexObject,1)")
               .Define("Vertex_d2PVz",    "myUtils::get_Vertex_d2PV(VertexObject,2)")

               .Define("Vertex_d2PVErr",  "myUtils::get_Vertex_d2PVError(VertexObject,-1)")
               .Define("Vertex_d2PVxErr", "myUtils::get_Vertex_d2PVError(VertexObject,0)")
               .Define("Vertex_d2PVyErr", "myUtils::get_Vertex_d2PVError(VertexObject,1)")
               .Define("Vertex_d2PVzErr", "myUtils::get_Vertex_d2PVError(VertexObject,2)")

               .Define("Vertex_d2PVSig_pre",  "Vertex_d2PV_pre/Vertex_d2PVErr")
               .Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
               .Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
               .Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

               .Define("Vertex_d2MC",     "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,-1)")
               .Define("Vertex_d2MCx",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,0)")
               .Define("Vertex_d2MCy",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,1)")
               .Define("Vertex_d2MCz",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,2)")

               .Define("EVT_dPV2DVmin",   "myUtils::get_dPV2DV_min(Vertex_d2PV_pre)")
               .Define("EVT_dPV2DVmax",   "myUtils::get_dPV2DV_max(Vertex_d2PV_pre)")
               .Define("EVT_dPV2DVave",   "myUtils::get_dPV2DV_ave(Vertex_d2PV_pre)")


               #############################################
               ##              Get RECO info              ##
               #############################################
               
               .Define("RP_e_pre",              "ReconstructedParticle::get_e(RecoPartPIDAtVertex)")
               .Define("RP_m_true_pre",         "ReconstructedParticle::get_mass(RecoPartPIDAtVertex)")
               .Define("RP_m_reco_pre",         "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_px_pre",             "ReconstructedParticle::get_px(RecoPartPIDAtVertex)")
               .Define("RP_py_pre",             "ReconstructedParticle::get_py(RecoPartPIDAtVertex)")
               .Define("RP_pz_pre",             "ReconstructedParticle::get_pz(RecoPartPIDAtVertex)")
               .Define("RP_eta",            "ReconstructedParticle::get_eta(RecoPartPIDAtVertex)")
               .Define("RP_phi",            "ReconstructedParticle::get_phi(RecoPartPIDAtVertex)")
               .Define("RP_theta",          "ReconstructedParticle::get_theta(RecoPartPIDAtVertex)")
               .Define("RP_charge_pre",         "ReconstructedParticle::get_charge(RecoPartPIDAtVertex)")
               .Define("RP_fromPV_pre",         "FCCAnalyses::ZHfunctions::get_RP_isfromPV(VertexObject,RecoPartPIDAtVertex)")
               .Define("RP_vert_ind_pre",       "FCCAnalyses::ZHfunctions::get_RP_Vert_Ind(VertexObject,RecoPartPIDAtVertex)")
               .Define("RP_vert_e_pre",         "ROOT::VecOps::RVec<float> result; for (auto & i: RP_vert_ind_pre) {if (i==-1) result.push_back(-1); else result.push_back(Vertex_e_pre.at(i));} return result;")
               .Define("RP_vert_mass_pre",      "ROOT::VecOps::RVec<float> result; for (auto & i: RP_vert_ind_pre) {if (i==-1) result.push_back(-1); else result.push_back(Vertex_mass_pre.at(i));} return result;")

               .Define("RP_trk_d0_pre",         "ReconstructedParticle2Track::getRP2TRK_D0       (RecoPartPIDAtVertex,EFlowTrack_1)")
               .Define("RP_trk_z0_pre",         "ReconstructedParticle2Track::getRP2TRK_Z0       (RecoPartPIDAtVertex,EFlowTrack_1)")
               .Define("RP_trk_phi_pre",        "ReconstructedParticle2Track::getRP2TRK_phi      (RecoPartPIDAtVertex,EFlowTrack_1)")
               .Define("RP_trk_omega_pre",      "ReconstructedParticle2Track::getRP2TRK_omega    (RecoPartPIDAtVertex,EFlowTrack_1)")
               .Define("RP_trk_tanLambda_pre",  "ReconstructedParticle2Track::getRP2TRK_tanLambda(RecoPartPIDAtVertex,EFlowTrack_1)")
               .Define("RP_dndx_pre",           "FCCAnalyses::ZHfunctions::get_RP_dndx(RecoPartPIDAtVertex, EFlowTrack_2, EFlowTrack)")
               .Define("RP_mtof_pre",           "FCCAnalyses::ZHfunctions::get_RP_mtof(RecoPartPIDAtVertex, EFlowTrack_L, EFlowTrack, TrackerHits, EFlowPhoton, EFlowNeutralHadron, CalorimeterHits, MC_PV_xyzt)")

               ###################################################
               ##    check MC mathc with certain decay chains   ##
               ###################################################

               .Define("RP_nMC_pre",        "FCCAnalyses::ZHfunctions::getRP2MC_nMC(MCRecoAssociations0,MCRecoAssociations1,RecoPartPIDAtVertex)")
               .Define("RP_MCidx_pre",      "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0,MCRecoAssociations1,RecoPartPIDAtVertex)")
               .Define("RP_fromBc_pre",     "FCCAnalyses::ZHfunctions::get_RP_isDescendant(541, true)(RP_MCidx_pre, Particle, Particle1)")
               .Define("RP_fromBs_pre",     "FCCAnalyses::ZHfunctions::get_RP_isDescendant(531, true)(RP_MCidx_pre, Particle, Particle1)")
               .Define("RP_fromBu_pre",     "FCCAnalyses::ZHfunctions::get_RP_isDescendant(521, true)(RP_MCidx_pre, Particle, Particle1)")
               .Define("RP_fromBd_pre",     "FCCAnalyses::ZHfunctions::get_RP_isDescendant(511, true)(RP_MCidx_pre, Particle, Particle1)")
               .Define("RP_fromLb_pre",     "FCCAnalyses::ZHfunctions::get_RP_isDescendant(5122, true)(RP_MCidx_pre, Particle, Particle1)")

               .Define("Vertex_fromBc_pre", "FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBc_pre)")
               .Define("Vertex_fromBs_pre", "FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBs_pre)")
               .Define("Vertex_fromBu_pre", "FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBu_pre)")
               .Define("Vertex_fromBd_pre", "FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBd_pre)")
               .Define("Vertex_fromLb_pre", "FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromLb_pre)")

               #############################################
               ##              Build the thrust           ##
               #############################################

               .Define("EVT_thrustNP",      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px_pre, RP_py_pre, RP_pz_pre)')
               .Define("RP_thrustangleNP",  'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px_pre, RP_py_pre, RP_pz_pre)')
               .Define("EVT_thrust",        'Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e_pre, EVT_thrustNP)')
               .Define("RP_thrustangle_pre",    'Algorithms::getAxisCosTheta(EVT_thrust, RP_px_pre, RP_py_pre, RP_pz_pre)')
               .Define("EVT_thrust_phi",    'FCCAnalyses::ZHfunctions::getAxisPhi(EVT_thrust)')
               .Define("EVT_thrust_theta",  'FCCAnalyses::ZHfunctions::getAxisTheta(EVT_thrust)')


               #############################################
               ##      thrust angle of gen b, Bs          ##
               #############################################

               .Define("genBs_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genBs_px, genBs_py, genBs_pz)')
               .Define("genBd_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genBd_px, genBd_py, genBd_pz)')
               .Define("genBu_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genBu_px, genBu_py, genBu_pz)')
               .Define("genBc_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genBc_px, genBc_py, genBc_pz)')
               .Define("genLb_thrustangle",        'Algorithms::getAxisCosTheta(EVT_thrust, genLb_px, genLb_py, genLb_pz)')
               .Define("recoEmiss_thrustangle",    'Algorithms::getAxisCosTheta(EVT_thrust, recoEmiss_px, recoEmiss_py, recoEmiss_pz)')

               #############################################
               ##        Get thrust related values        ##
               #############################################
               ##hemis0 == negative angle == max energy hemisphere if pointing
               ##hemis1 == positive angle == min energy hemisphere if pointing
               .Define("EVT_thrusthemis0_n",    "Algorithms::getAxisN(0)(RP_thrustangle_pre, RP_charge_pre)")
               .Define("EVT_thrusthemis1_n",    "Algorithms::getAxisN(1)(RP_thrustangle_pre, RP_charge_pre)")
               .Define("EVT_thrusthemis0_e",    "Algorithms::getAxisEnergy(0)(RP_thrustangle_pre, RP_charge_pre, RP_e_pre)")
               .Define("EVT_thrusthemis1_e",    "Algorithms::getAxisEnergy(1)(RP_thrustangle_pre, RP_charge_pre, RP_e_pre)")
               .Define("EVT_thrusthemis0_p",    "Algorithms::getAxisMomentum(0)(RP_thrustangle_pre, RP_px_pre, RP_py_pre, RP_pz_pre)")
               .Define("EVT_thrusthemis1_p",    "Algorithms::getAxisMomentum(1)(RP_thrustangle_pre, RP_px_pre, RP_py_pre, RP_pz_pre)")

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

               #.Filter("EVT_ThrustEmin_NTau23PiCand > 1 && recoEmiss_e > 3.5 && EVT_ThrustEmin_Eneutral < 10 && EVT_ThrustEmin_Nneutral <= 13 && EVT_ThrustEmin_E < 38.")

               ##############################
               ##      Get the leptons     ##
               ##############################

               #concat the two indices collections to treat all of them in one go
               .Alias("Electron0","Electron#0.index")
               .Alias("Muon0","Muon#0.index")
               .Define("Lepton0","ROOT::VecOps::RVec<int> result; for(size_t i=0; i<Electron0.size(); ++i){result.push_back(Electron0.at(i));} for(size_t j=0; j<Muon0.size(); ++j){result.push_back(Muon0.at(j));} return result;")

               .Define("n_ell",      "ReconstructedParticle::get_n(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_d0",     "ReconstructedParticle2Track::getRP2TRK_D0(ReconstructedParticle::get(Lepton0,ReconstructedParticles),EFlowTrack_1)")
               .Define("ell_z0",     "ReconstructedParticle2Track::getRP2TRK_Z0(ReconstructedParticle::get(Lepton0,ReconstructedParticles),EFlowTrack_1)")
               .Define("ell_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(Lepton0,ReconstructedParticles))")
               .Define("ell_thrustangles",'Algorithms::getAxisCosTheta(EVT_thrust, ell_px, ell_py, ell_pz)')
               .Define("ell_thrustEmin_n","int result (0); for (size_t i=0; i<ell_thrustangles.size(); ++i){if (ell_thrustangles[i] > 0.0) ++result;} return result;")
               .Define("ell_thrustEmax_n","int result (0); for (size_t i=0; i<ell_thrustangles.size(); ++i){if (ell_thrustangles[i] <= 0.0) ++result;} return result;")


               ##############################
               ##       Get the pions      ##
               ##############################

               #Use perfect PID from PIDatVertex
               .Define("pi_ind","ROOT::VecOps::RVec<int> result; for (size_t i=0; i<RP_m_true_pre.size(); ++i) if (std::abs(RP_m_true_pre.at(i)-0.13957039) < 1e-3) result.push_back(i); return result;")

               .Define("n_pi",      "ReconstructedParticle::get_n(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_px",     "ReconstructedParticle::get_px(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_py",     "ReconstructedParticle::get_py(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_pz",     "ReconstructedParticle::get_pz(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_phi",    "ReconstructedParticle::get_phi(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_eta",    "ReconstructedParticle::get_eta(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_d0",     "ReconstructedParticle2Track::getRP2TRK_D0(ReconstructedParticle::get(pi_ind,ReconstructedParticles),EFlowTrack_1)")
               .Define("pi_z0",     "ReconstructedParticle2Track::getRP2TRK_Z0(ReconstructedParticle::get(pi_ind,ReconstructedParticles),EFlowTrack_1)")
               .Define("pi_energy", "ReconstructedParticle::get_e(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_mass",   "ReconstructedParticle::get_mass(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_charge", "ReconstructedParticle::get_charge(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_PDG",    "ReconstructedParticle::get_type(ReconstructedParticle::get(pi_ind,ReconstructedParticles))")
               .Define("pi_thrustangles",'Algorithms::getAxisCosTheta(EVT_thrust, pi_px, pi_py, pi_pz)')
               .Define("pi_thrustEmin_n","int result (0); for (size_t i=0; i<pi_thrustangles.size(); ++i){if (pi_thrustangles[i] > 0.0) ++result;} return result;")
               .Define("pi_thrustEmax_n","int result (0); for (size_t i=0; i<pi_thrustangles.size(); ++i){if (pi_thrustangles[i] <= 0.0) ++result;} return result;")
               
               ###############################
               ##   Had tagger columns      ##
               ###############################
               
               .Define("RP_e",              "RP_e_pre     [RP_thrustangle_pre>0]")
               .Define("RP_m_true",         "RP_m_true_pre[RP_thrustangle_pre>0]")
               .Define("RP_m_reco",         "RP_m_reco_pre[RP_thrustangle_pre>0]")
               .Define("RP_px",             "RP_px_pre    [RP_thrustangle_pre>0]")
               .Define("RP_py",             "RP_py_pre    [RP_thrustangle_pre>0]")
               .Define("RP_pz",             "RP_pz_pre    [RP_thrustangle_pre>0]")
               .Define("RP_Dphi",           "RP_phi   [RP_thrustangle_pre>0] - EVT_thrust_phi")
               .Define("RP_Dtheta",         "RP_theta [RP_thrustangle_pre>0] - EVT_thrust_theta")
               .Define("RP_charge",         "RP_charge_pre[RP_thrustangle_pre>0]")
               .Define("RP_fromPV",         "RP_fromPV_pre[RP_thrustangle_pre>0]")
               .Define("RP_vert_ind",       "RP_vert_ind_pre [RP_thrustangle_pre>0]")
               .Define("RP_vert_e",         "RP_vert_e_pre   [RP_thrustangle_pre>0]")
               .Define("RP_vert_mass",      "RP_vert_mass_pre[RP_thrustangle_pre>0]")
               .Define("RP_trk_d0",         "RP_trk_d0_pre       [RP_thrustangle_pre>0]")
               .Define("RP_trk_z0",         "RP_trk_z0_pre       [RP_thrustangle_pre>0]")
               .Define("RP_trk_phi",        "RP_trk_phi_pre      [RP_thrustangle_pre>0]")
               .Define("RP_trk_omega",      "RP_trk_omega_pre    [RP_thrustangle_pre>0]")
               .Define("RP_trk_tanLambda",  "RP_trk_tanLambda_pre[RP_thrustangle_pre>0]")
               .Define("RP_dndx",           "RP_dndx_pre         [RP_thrustangle_pre>0]")
               .Define("RP_mtof",           "RP_mtof_pre         [RP_thrustangle_pre>0]")
               .Define("RP_nMC",            "RP_nMC_pre   [RP_thrustangle_pre>0]")
               .Define("RP_MCidx",          "RP_MCidx_pre [RP_thrustangle_pre>0]")
               .Define("RP_fromBc",         "RP_fromBc_pre[RP_thrustangle_pre>0]")
               .Define("RP_fromBs",         "RP_fromBs_pre[RP_thrustangle_pre>0]")
               .Define("RP_fromBu",         "RP_fromBu_pre[RP_thrustangle_pre>0]")
               .Define("RP_fromBd",         "RP_fromBd_pre[RP_thrustangle_pre>0]")
               .Define("RP_fromLb",         "RP_fromLb_pre[RP_thrustangle_pre>0]")
               .Define("RP_thrustangle",    "RP_thrustangle_pre[RP_thrustangle_pre>0]")
               .Define("Vertex_isPV",       "Vertex_isPV_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_mass",       "Vertex_mass_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_px",         "Vertex_px_pre       [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_py",         "Vertex_py_pre       [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_pz",         "Vertex_pz_pre       [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_e",          "Vertex_e_pre        [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_x",          "Vertex_x_pre        [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_y",          "Vertex_y_pre        [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_z",          "Vertex_z_pre        [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_xErr",       "Vertex_xErr_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_yErr",       "Vertex_yErr_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_zErr",       "Vertex_zErr_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_ntrk",       "Vertex_ntrk_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_chi2",       "Vertex_chi2_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_d2PV",       "Vertex_d2PV_pre     [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_d2PVSig",    "Vertex_d2PVSig_pre  [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_Dphi",       "Vertex_phi      [Vertex_thrust_angle>0 || Vertex_isPV_pre==1] - EVT_thrust_phi")
               .Define("Vertex_Dtheta",     "Vertex_theta    [Vertex_thrust_angle>0 || Vertex_isPV_pre==1] - EVT_thrust_theta")
               .Define("Vertex_thrustangle","Vertex_thrust_angle    [Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_fromBc",     "Vertex_fromBc_pre[Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_fromBs",     "Vertex_fromBs_pre[Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_fromBu",     "Vertex_fromBu_pre[Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_fromBd",     "Vertex_fromBd_pre[Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")
               .Define("Vertex_fromLb",     "Vertex_fromLb_pre[Vertex_thrust_angle>0 || Vertex_isPV_pre==1]")

               .Define("n_Bc",              "int(genBc_thrustangle[genBc_thrustangle>0].size())")
               .Define("n_Bs",              "int(genBs_thrustangle[genBs_thrustangle>0].size())")
               .Define("n_Bu",              "int(genBu_thrustangle[genBu_thrustangle>0].size())")
               .Define("n_Bd",              "int(genBd_thrustangle[genBd_thrustangle>0].size())")
               .Define("n_Lb",              "int(genLb_thrustangle[genLb_thrustangle>0].size())")
               .Define("label_Bc",          "int(n_Bc==1 && n_Bs==0 && n_Bu==0 && n_Bd==0 && n_Lb==0)" )
               .Define("label_Bs",          "int(n_Bc==0 && n_Bs==1 && n_Bu==0 && n_Bd==0 && n_Lb==0)" )
               .Define("label_Bu",          "int(n_Bc==0 && n_Bs==0 && n_Bu==1 && n_Bd==0 && n_Lb==0)" )
               .Define("label_Bd",          "int(n_Bc==0 && n_Bs==0 && n_Bu==0 && n_Bd==1 && n_Lb==0)" )
               .Define("label_Lb",          "int(n_Bc==0 && n_Bs==0 && n_Bu==0 && n_Bd==0 && n_Lb==1)" )
               .Define("label_light",       "int(n_Bc==0 && n_Bs==0 && (n_Bu>0 || n_Bd>0) && n_Lb==0)" )
               .Define("label_hasBc",       "int(n_Bc>0)" )
               .Define("label_has1Bc",      "int(n_Bc)" )


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

                "n_genBs",
                "genBs_px", "genBs_py", "genBs_pz", "genBs_eta", "genBs_phi",
                "genBs_energy", "genBs_mass", "genBs_pdg",
                "genBs_thrustangle",

                "genBs_Vertex_x", "genBs_Vertex_y", "genBs_Vertex_z",

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

                "Vertex_x_pre", "Vertex_y_pre", "Vertex_z_pre",
                "Vertex_xErr_pre", "Vertex_yErr_pre", "Vertex_zErr_pre",
                "Vertex_isPV_pre", "Vertex_ntrk_pre", "Vertex_chi2_pre", "Vertex_n",
                "Vertex_thrust_angle", "Vertex_thrusthemis_emin", "Vertex_thrusthemis_emax",

                "Vertex_d2PV_pre", "Vertex_d2PVx", "Vertex_d2PVy", "Vertex_d2PVz",
                "Vertex_d2PVErr", "Vertex_d2PVxErr", "Vertex_d2PVyErr", "Vertex_d2PVzErr",
                "Vertex_mass_pre",
                "DV_d0","DV_z0",

                "recoEmiss_px", "recoEmiss_py", "recoEmiss_pz", "recoEmiss_e", "recoEmiss_m",
                "recoEmiss_thrustangle",

                "n_ell","ell_px","ell_py","ell_pz","ell_phi","ell_eta","ell_d0","ell_z0","ell_energy","ell_mass","ell_charge","ell_PDG","ell_thrustangles","ell_thrustEmin_n","ell_thrustEmax_n",
                "n_pi","pi_px","pi_py","pi_pz","pi_phi","pi_eta","pi_d0","pi_z0","pi_energy","pi_mass","pi_charge","pi_PDG","pi_thrustangles","pi_thrustEmin_n","pi_thrustEmax_n",

                "RP_e","RP_m_true","RP_m_reco","RP_px","RP_py","RP_pz","RP_Dphi","RP_Dtheta","RP_charge","RP_fromPV",
                "RP_vert_ind","RP_vert_e","RP_vert_mass","RP_trk_d0","RP_trk_z0","RP_trk_phi","RP_trk_omega","RP_trk_tanLambda","RP_dndx","RP_mtof",
                "RP_nMC","RP_MCidx","RP_fromBc","RP_fromBs","RP_fromBu","RP_fromBd","RP_fromLb","RP_thrustangle",
                "Vertex_isPV","Vertex_mass", "Vertex_px","Vertex_py","Vertex_pz","Vertex_e","Vertex_x","Vertex_y","Vertex_z","Vertex_xErr","Vertex_yErr","Vertex_zErr","Vertex_ntrk","Vertex_chi2",
                "Vertex_d2PV","Vertex_d2PVSig","Vertex_Dphi","Vertex_Dtheta","Vertex_thrustangle","Vertex_fromBc","Vertex_fromBs","Vertex_fromBu","Vertex_fromBd","Vertex_fromLb",
                "n_Bc","n_Bs","n_Bu","n_Bd","n_Lb","label_Bc","label_Bs","label_Bu","label_Bd","label_Lb","label_light","label_hasBc","label_has1Bc",

                ]
        return branchList
