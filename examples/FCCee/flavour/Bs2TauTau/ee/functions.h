#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"


namespace FCCAnalyses { namespace ZHfunctions {


// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass_recoil {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    Vec_rp operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) ;
};

resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

Vec_rp resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) {

    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            rp reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        Vec_rp bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }

     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}    




struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    Vec_rp operator() (Vec_rp in, Vec_f iso);
  };

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
    Vec_rp result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}

 
// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    Vec_f operator() (Vec_rp in, Vec_rp rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
  
    Vec_f result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }
    
    for(size_t i = 0; i < in.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }

    
    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) {
    
        double sumNeutral = 0.0;
        double sumCharged = 0.0;
    
        // charged
        for (auto & lv_charged_ : lv_charged) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
            if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }
        
        // neutral
        for (auto & lv_neutral_ : lv_neutral) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
            if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }
        
        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}
 
 
 
// returns missing energy vector, based on reco particles
Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for(auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e += p.energy;
    }
    
    Vec_rp ret;
    rp res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm-e;
    ret.emplace_back(res);
    return ret;
}

// calculate the cosine(theta) of the missing energy vector
float get_cosTheta_miss(Vec_rp met){
    
    float costheta = 0.;
    if(met.size() > 0) {
        
        TLorentzVector lv_met;
        lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
        costheta = fabs(std::cos(lv_met.Theta()));
    }
    return costheta;
}

 
ROOT::VecOps::RVec<TLorentzVector> build_p4(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz, ROOT::VecOps::RVec<float> mass) {
    ROOT::VecOps::RVec<TLorentzVector> p4;
    for (size_t i = 0; i < px.size(); ++i) {
        TLorentzVector tlv;
        tlv.SetXYZM(px[i], py[i], pz[i], mass[i]);
        p4.push_back(tlv);
    }
    return p4;
} 

//std::vector<int> get_MCpdgMCVertex(std::vector<std::vector<int>> vertex_mother_PDG,
//                                   std::vector<std::vector<int>> vertex_daughter_PDG,
//				   int require_mother,
//				   int require_daughter){  //TODO: need to change daughter to a vector
//  std::vector<int> result;
//  for (size_t i=0; i < vertex_mother_PDG.size(); ++i){
//    int mo_yes=0;
//    int da_yes=0;
//    for (size_t j=0; j < )
//    std::vector<int> tmp;
//    for (size_t i = 0; i < p.mc_ind.size(); ++i) tmp.push_back(mc.at(p.mc_ind.at(i)).PDG);
//    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) tmp.push_back(mc.at(p.mc_indneutral.at(i)).PDG);
//    result.push_back(tmp);
//  }
//  return result;
//}


//A set of functions to compute the angle between two tau vertices/momenta and then select pairs of tau candidates based on that angle

float Compute_CosTheta(float axis_x, float axis_y, float axis_z, float x, float y, float z){
    float num = axis_x*x + axis_y*y + axis_z*z;
    float den = sqrt(axis_x*axis_x+axis_y*axis_y+axis_z*axis_z)*sqrt(x*x+y*y+z*z);
    return num/den;
}

//Compute the angle between each pair of Tau candidates
ROOT::VecOps::RVec<float> Compute_TauCand_MomentaAngles(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz){
    ROOT::VecOps::RVec<float> result; 
    for (size_t i=0;i<px.size();++i) {
        for (size_t j=i+1;j<px.size();++j) result.push_back(Compute_CosTheta(px[i],py[i],pz[i],px[j],py[j],pz[j]));
    } 
    return result;
}

//Selects event if at least one angle between two tau candidates if within a certain range (signal peaks at cos(angle) ~ 0.95 while background peaks at 1)
bool CloseEnough(ROOT::VecOps::RVec<float> in, float cut_low, float cut_high){
    int Counter = 0;
    for (size_t i=0; i < in.size(); ++i){
        if (in[i] > cut_low and in[i] < cut_high) Counter += 1;
    }
    if (Counter == 0) return false;
    else return true;
}


//----------------- Select specific tau decays to analyse exclusive tau decay channels --------------------------------------------------

//Get the gen Bs2TauTau + additional photons (used as a basis to then look for the exclusive tau decays)
ROOT::VecOps::RVec<edm4hep::MCParticleData> Find_genBs2TauTau(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> daughter){
    
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    
    for (size_t i=0; i < in.size(); ++i) {
        
        if (std::abs(in.at(i).PDG) == 531 && in.at(i).daughters_end != in.at(i).daughters_begin){   // && in.at(i).daughters_end-in.at(i).daughters_begin == 2
            int IndTauplus(-999);
            int IndTauminus(-999);
            int TauCounter(0);
            ROOT::VecOps::RVec<int> IndRem;
                
            for (size_t j=in.at(i).daughters_begin; j < in.at(i).daughters_end; ++j){
                if (in[daughter.at(j)].PDG == 15) {
                    TauCounter += 1;
                    IndTauminus = daughter.at(j);
                }
                else if (in[daughter.at(j)].PDG == -15) {
                    TauCounter += 1;
                    IndTauplus = daughter.at(j);
                }
                else IndRem.push_back(daughter.at(j));
            }
                
            if (IndTauplus != -999 && IndTauminus != -999 && TauCounter == 2){
                result.push_back(in.at(i));
                result.push_back(in.at(IndTauplus));
                result.push_back(in.at(IndTauminus));
                for (size_t k = 0; k < IndRem.size(); ++k){
                    result.push_back(in.at(IndRem[k]));
                }                
            }
        }
    }
    return result; //Structure (genBs, genTauPlus, genTauMinus, Photon0, Photon1, ...)
}

//Get the exclusive Tau2eNuNu Tau decay in the Bs2TauTau decays WARNING: PM = 1 GIVES TAU_MINUS WHILE -1 GIVES TAU_PLUS
ROOT::VecOps::RVec<edm4hep::MCParticleData> Find_genBs2TauTau_Electrons(ROOT::VecOps::RVec<edm4hep::MCParticleData> genBs2TauTau_list, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> daughter, int pm){
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    for (size_t i = 0; i < genBs2TauTau_list.size(); ++i){
        if (genBs2TauTau_list[i].PDG == pm*15){
            ROOT::VecOps::RVec<int> Electrons;
            ROOT::VecOps::RVec<int> Rem;
            for (size_t j = genBs2TauTau_list[i].daughters_begin; j < genBs2TauTau_list[i].daughters_end; ++j){
                if (std::abs(in[daughter.at(j)].PDG) == 11){
                    Electrons.push_back(daughter.at(j));
                }
                else {
                    Rem.push_back(daughter.at(j));
                }
            }
            if (Electrons.size() == 1){
                for (size_t m = 0; m < Electrons.size(); ++m){
                    result.push_back(in.at(Electrons[m]));
                }
                for (size_t n = 0; n < Rem.size(); ++n){
                    result.push_back(in.at(Rem[n]));
                }
            }
        }
    }
    return result; //Structure (electron, nu0, nu1, photon0, photon1, ...)
}


//Truth-match the muons, use a pm flag (+1 = mu+, -1 = mu-) to create two collection of reco muons truth-matched to the MC ones
//Perform the truth matching by going through the reco part that got truth-matched
//Checks that the MC is a muon whose mother is a Tau AND grand mother is a Bs AND the associated reco is in the muon collection
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> Muon_TruthMatching(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco_in, 
                                                                          ROOT::VecOps::RVec<int> reco_ind, 
                                                                          ROOT::VecOps::RVec<int> mc_ind, 
                                                                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc_in, 
                                                                          ROOT::VecOps::RVec<int> mc_parents,
                                                                          int pm,
                                                                          ROOT::VecOps::RVec<int> muon_ind) {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

    for (size_t i=0; i < reco_ind.size(); ++i) {
        if (mc_in.at(mc_ind.at(i)).PDG == pm*(-13)) { //Check if the corresponding mc is a mu
            for (size_t j=mc_in.at(mc_ind.at(i)).parents_begin; j < mc_in.at(mc_ind.at(i)).parents_end; ++j){
                if (std::abs(mc_in.at(mc_parents.at(j)).PDG) == 15){ //Check if the parent of that mu is a tau
                    for (size_t k=mc_in.at(mc_parents.at(j)).parents_begin; k < mc_in.at(mc_parents.at(j)).parents_end; ++k){
                        if (std::abs(mc_in.at(mc_parents.at(k)).PDG) == 531){ //Check if the grand-parent of that mu is a Bs
                            for (size_t l=0; l<muon_ind.size(); ++l){
                                if (muon_ind.at(l) == reco_ind.at(i)){ //Check that the id'd index is pointing to a particle that is also being pointed at by an index in the muon collection
                                    result.push_back(reco_in.at(reco_ind.at(i)));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}


//!!!FAILED ATTEMPT!!! Perform Truth-matching procedure with the muon collection only (not the full reco particles)
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> Muon_TruthMatching2(ROOT::VecOps::RVec<int> muons_ind,
                                                                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco_in,
                                                                          ROOT::VecOps::RVec<float> RP2MC,
                                                                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc_in, 
                                                                          ROOT::VecOps::RVec<int> mc_parents,
                                                                          int pm) {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

    //Check MC parents from truth-matched muons to established if coming from Bs->TauTau or else
    for (size_t i=0; i < muons_ind.size(); ++i) {
        edm4hep::MCParticleData MC_muon (mc_in.at(RP2MC.at(muons_ind.at(i)))); //From my understanding, muons_ind should give back an index in the reco_in collection and it seems that reco_in and reco_ind are one-to-one
        for (size_t j=MC_muon.parents_begin; j < MC_muon.parents_end; ++j){
            if (std::abs(mc_in.at(mc_parents.at(j)).PDG) == pm*(-15)){
                for (size_t k=mc_in.at(mc_parents.at(j)).parents_begin; k < mc_in.at(mc_parents.at(j)).parents_end; ++k){
                    if (std::abs(mc_in.at(mc_parents.at(k)).PDG) == 531){
                        result.push_back(reco_in.at(muons_ind.at(i)));
                    }
                }
            }    
        }
    }
    return result;
}



//Truth-matching is a bit different from MC identifiction since we only look at a subset of MC part that have been truth-matched to reco
//Select the subset of MC part that are truth-matched involved in the decay of interest, disregard the charge for mothers while one can specify which type of particle it wants to keep at the end with charge taken into account
ROOT::VecOps::RVec<int> TruthMatch_MC_Mothers(ROOT::VecOps::RVec<int> mc_ind, 
                                              ROOT::VecOps::RVec<edm4hep::MCParticleData> mc_in, 
                                              ROOT::VecOps::RVec<int> mc_parents,
                                              int Mother_PDG,
                                              int gMother_PDG,
                                              int daughter_PDG=0){

    ROOT::VecOps::RVec<int> all_daughters;

    //Checks that the truth-matched daughters are coming from a specific decays
    for (size_t i=0; i < mc_ind.size(); ++i){
        for (size_t j=mc_in.at(mc_ind.at(i)).parents_begin; j<mc_in.at(mc_ind.at(i)).parents_end; ++j){
            if (std::abs(mc_in.at(mc_parents.at(j)).PDG) == Mother_PDG){
                for (size_t k=mc_in.at(mc_parents.at(j)).parents_begin; k<mc_in.at(mc_parents.at(j)).parents_end; ++k){
                    if (std::abs(mc_in.at(mc_parents.at(k)).PDG) == gMother_PDG){
                        all_daughters.push_back(mc_ind.at(i));
                    }
                }
            }
        }
    }

    //By default returns all the daughters, one can specify which particles it wants
    if (daughter_PDG == 0) return all_daughters;
    else {
        ROOT::VecOps::RVec<int> results;
        for (size_t i=0; i < all_daughters.size(); ++i){
            if (mc_in.at(all_daughters.at(i)).PDG == daughter_PDG){
                results.push_back(all_daughters.at(i));
            }
        }
        return results;
    }
}

//Truth-matching procedure in the case where the MC and RECO part subsets of interest are known (for instance from MC decay and muon RECO) 
ROOT::VecOps::RVec<int> TruthMatch_MC_RECO_Subset(ROOT::VecOps::RVec<int> reco_ind,
                                                  ROOT::VecOps::RVec<int> reco_ind_subset,
                                                  ROOT::VecOps::RVec<int> mc_ind_subset,
                                                  ROOT::VecOps::RVec<int> mc_ind){

    ROOT::VecOps::RVec<int> results;
    ROOT::VecOps::RVec<int> Reco_list;
    ROOT::VecOps::RVec<int> Mc_list;

    //Find indices position in the full truth-matched RECO indices collection by matching them to the reco indices subcollection
    for (size_t i=0; i<reco_ind_subset.size(); ++i){
        for (size_t j=0; j<reco_ind.size(); ++j){
            if (reco_ind_subset.at(i) == reco_ind.at(j)){
                Reco_list.push_back(j);
            }
        }
    }

    //Same as previous loop but in the full truth-matched MC indices
    for (size_t i=0; i<mc_ind_subset.size(); ++i){
        for (size_t j=0; j<mc_ind.size(); ++j){
            if (mc_ind_subset.at(i) == mc_ind.at(j)){
                Mc_list.push_back(j);
            }
        }
    }

    //Compare the MC and RECO indices position, if they are on the same line, save the reco index
    for (size_t i=0; i<Reco_list.size(); ++i){
        for (size_t j=0; j<Mc_list.size(); ++j){
            if (Mc_list.at(j) == Reco_list.at(i)){
                results.push_back(reco_ind.at(Reco_list.at(i)));
            }
        }
    }

    return results;
}

//Modification of the previous function to study truth-matched RECO from ALL MC particles (TM MC ind duplications) 
ROOT::VecOps::RVec<int> TruthMatch_MC_RECO(ROOT::VecOps::RVec<int> reco_ind,
                                           ROOT::VecOps::RVec<int> mc_ind_subset,
                                           ROOT::VecOps::RVec<int> mc_ind){

    ROOT::VecOps::RVec<int> results;
    ROOT::VecOps::RVec<int> ReadMCList;
    ROOT::VecOps::RVec<int> Mc_list;

    //Get the  indices (dupplication shows up as subMCind1->[MCind1,MCind2], subMCind2->[MCind1,MCind2]=>[MCind1,MCind2,MCind1,MCind2])
    for (size_t i=0; i<mc_ind_subset.size(); ++i){
        for (size_t j=0; j<mc_ind.size(); ++j){
            if (mc_ind_subset.at(i) == mc_ind.at(j)){
                Mc_list.push_back(j);
            }
        }
    }

    //Check duplication and returns unique indices
    for (size_t j=0; j<Mc_list.size(); ++j){
        bool New (true);
        for (size_t k=0; k<ReadMCList.size(); ++k){
            if (Mc_list.at(j) == ReadMCList.at(k)){ 
                New = false;
                break;
            }
        }
        if (New){
                results.push_back(reco_ind.at(Mc_list.at(j)));
                ReadMCList.push_back(Mc_list.at(j));
        }
    }

    return results;
}

//Compute the opening angle between the TM muons
float TM_ComputeOpeningAngle(ROOT::VecOps::RVec<float> px1,
                             ROOT::VecOps::RVec<float> py1,
                             ROOT::VecOps::RVec<float> pz1,
                             ROOT::VecOps::RVec<float> px2,
                             ROOT::VecOps::RVec<float> py2,
                             ROOT::VecOps::RVec<float> pz2){
    
    if (px1.size() == 0 || px2.size() == 0){
        return -2.0;
    }
    else{
        return (px1[0]*px2[0] + py1[0]*py2[0] + pz1[0]*pz2[0])/sqrt(px1[0]*px1[0]+py1[0]*py1[0]+pz1[0]*pz1[0])/sqrt(px2[0]*px2[0]+py2[0]*py2[0]+pz2[0]*pz2[0]);
    }
}

ROOT::VecOps::RVec<float> Muons_ComputeOpeningAngle(ROOT::VecOps::RVec<float> px,ROOT::VecOps::RVec<float> py,ROOT::VecOps::RVec<float> pz){

    ROOT::VecOps::RVec<float> results;

    if (px.size() < 2){
        results.push_back(-2.0);
        return results;
    }
    else{
        
        for (size_t i=0; i<px.size(); ++i){
            for (size_t j=i+1; j<px.size(); ++j){
                results.push_back((px[i]*px[j]+py[i]*py[j]+pz[i]*pz[j])/sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])/sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]));
            }
        }

        return results;
    }
}


//Select the smallest opening angle in the list to ID the di-muon system of interest
float GetMin(ROOT::VecOps::RVec<float> list){
    if (list.size() == 0) return -2.0;
    else if (list.size() == 1) return list[0];
    else{
        float result (list[0]);
        for (size_t i=1; i<list.size();++i){
            if (result < list[i]) result = list[i];
        }
        return result;
    }
}

//ID the two muons with smallest opening angle with their places index in the muon collection 
ROOT::VecOps::RVec<int> Muons_ID_SmallestOA(ROOT::VecOps::RVec<float> OAs, int n_muons){
    
    ROOT::VecOps::RVec<int> ind_list;
    
    //Not enough muon to compute an OA (== only 1 element and equal to -2)
    if (std::abs(OAs.at(0) + 2.0) < 1e-4) ind_list.push_back(-1);
    
    //Only one pair of muon (obvious case)
    else if (OAs.size() == 1) {
        ind_list.push_back(0); 
        ind_list.push_back(1);
    }
    
    //Proceed to find min and return the corresponding ind that produced the min
    else{

        //Create the indices map
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> map;
        for (int i=0; i<n_muons; ++i){
            for (int j=i+1; j<n_muons; ++j){
                ROOT::VecOps::RVec<int> temp = {i,j};
                map.push_back(temp);
            }
        }
        
        //find OA min and select corresponding muon indices through the map created previously
        float min (OAs[0]);
        ind_list = map[0];

        for (size_t i=1; i<OAs.size();++i){
            if (min < OAs[i]) {
                min = OAs[i];
                ind_list = map[i];
            }
        }
    }
    return ind_list;
}


//Functions used to define the Stage 1 cut criteria (dimuon properties: Opening angle, hemisphere emission, charge)
int Check_dimuon_Presence(ROOT::VecOps::RVec<float> muon_OAs){
    if (muon_OAs.size() == 1 && std::abs(muon_OAs[0]+2.0) < 1e-4) return 0;
    else return 1;
}

int Check_dimuon_SameSide(ROOT::VecOps::RVec<float> muon_OAs, int has_dimuon){
    if (has_dimuon < 1) return 0;
    if (GetMin(muon_OAs) < 0.0) return 0;
    else return 1;
}

int Check_dimuon_SigHemi(ROOT::VecOps::RVec<int> dimuon_ind, ROOT::VecOps::RVec<float> muon_thrustangles, int has_dimuon){
    if (has_dimuon < 1) return 0;
    if (muon_thrustangles[dimuon_ind.at(0)] < 0.0 && muon_thrustangles[dimuon_ind.at(1)] < 0.0) return 0;
    else return 1;
}

int Check_dimuon_Charges(ROOT::VecOps::RVec<int> dimuon_ind, ROOT::VecOps::RVec<int> muon_charges, int has_dimuon){
    if (has_dimuon < 1) return 0;
    if (muon_charges[dimuon_ind.at(0)]*muon_charges[dimuon_ind.at(1)] > 0) return 0;
    else return 1;
}

//Check if muons have associated vertices (-1 == no dimuon, 0 == no vertex associated, 10 == only muminus has a vertex, 1 (01) == only muplus has a vertex, 11 == both have associated vertices)
int Check_dimuon_Vertices(ROOT::VecOps::RVec<int> dimuon_ind, ROOT::VecOps::RVec<int> muon_ind, ROOT::VecOps::RVec<int> muon_q, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> Vertices, int has_dimuon){
    if (has_dimuon < 1) return -1;
    
    int muplus_ind;
    int muminus_ind;

    if (muon_q.at(dimuon_ind.at(0)) < 0){
        muplus_ind = muon_ind.at(dimuon_ind.at(1));
        muminus_ind = muon_ind.at(dimuon_ind.at(0));
    }
    else {
        muplus_ind = muon_ind.at(dimuon_ind.at(0));
        muminus_ind = muon_ind.at(dimuon_ind.at(1));
    }

    int muplus_hasVertex (0);
    int muminus_hasVertex (0);

    //Proceed with the full search among all particle coming from vertices
    for (auto &p:Vertices){
        for (auto &r:p.reco_ind){
            if (r == muplus_ind) muplus_hasVertex = 1;
            if (r == muminus_ind) muminus_hasVertex = 1;   
        }
    }

    return 10*muminus_hasVertex + muplus_hasVertex; 
}

//Finding the MC muons truth-matched to dimuon
ROOT::VecOps::RVec<int> Finding_MC_dimuon(ROOT::VecOps::RVec<int> dimuon_ind,
                                          ROOT::VecOps::RVec<int> reco_muon_ind,
                                          ROOT::VecOps::RVec<int> reco_ind,
                                          ROOT::VecOps::RVec<int> mc_ind){

    ROOT::VecOps::RVec<int> results;

    //Extended usage to events without dimuon                                        
    if (dimuon_ind.size() < 2) {
        results.push_back(-1);
    }
    
    else {
        for (size_t i=0;i<dimuon_ind.size();++i){
            for (size_t j=0;j<reco_ind.size();++j){
                if (reco_muon_ind.at(dimuon_ind.at(i)) == reco_ind.at(j)){ //the dimuon indices point in the reco muon subcollection 
                    results.push_back(mc_ind.at(j));
                }
            }
        }
    }
    //Shouldn't happen in that direction but the truth-matching might fail
    //Guaranteed two reco muon because of the cuts but no guarantee on the MC side although pretty sure that it should work
    //The fast sim use a lot the MC info to produce RECO info thus not sure that it can produce ghost particles
    //Add a check by checking the size of the MC_dimuon list because if truth-matching fails then it doesn't appear in the MCRecoAssociation lists, and thus not push_back'ed
    return results;
}

//Explore the MC tree of the dimuon: Find the common ancestor of both muon (assume only one mother)
int Find_MC_CommonAncestor(ROOT::VecOps::RVec<int> MC_dimuon_ind,
                           ROOT::VecOps::RVec<edm4hep::MCParticleData> Particle,
                           ROOT::VecOps::RVec<int> Parents_ind){


    if (MC_dimuon_ind.size() < 2) return -1; //Extended usage to events without dimuon


    //Assume only one mother for all particle (actually checks it and use it as a termination condition (-1))
    ROOT::VecOps::RVec<int> Parents_Mu1;
    ROOT::VecOps::RVec<int> Parents_Mu2;

    if (Particle.at(MC_dimuon_ind[0]).parents_begin+1 == Particle.at(MC_dimuon_ind[0]).parents_end){
        Parents_Mu1.push_back(Parents_ind.at(Particle.at(MC_dimuon_ind[0]).parents_begin)); 
        bool SingleParent1 (true);
        do {
            if (Particle.at(Parents_Mu1.back()).parents_begin+1 == Particle.at(Parents_Mu1.back()).parents_end){
                Parents_Mu1.push_back(Parents_ind.at(Particle.at(Parents_Mu1.back()).parents_begin));
            }
            else {
                Parents_Mu1.push_back(-1);
                SingleParent1 = false;
            }
        } while(SingleParent1);
    }
    else Parents_Mu1.push_back(-1);
         
    if (Particle.at(MC_dimuon_ind[1]).parents_begin+1 == Particle.at(MC_dimuon_ind[1]).parents_end){
        Parents_Mu2.push_back(Parents_ind.at(Particle.at(MC_dimuon_ind[1]).parents_begin)); 
        bool SingleParent2 (true);
        do {
            if (Particle.at(Parents_Mu2.back()).parents_begin+1 == Particle.at(Parents_Mu2.back()).parents_end){
                Parents_Mu2.push_back(Parents_ind.at(Particle.at(Parents_Mu2.back()).parents_begin));
            }
            else {
                Parents_Mu2.push_back(-1);
                SingleParent2 = false;
            }
        } while(SingleParent2);
    }
    else Parents_Mu2.push_back(-1);
    

    //Check who's the common ancestor, -1 == No common ancestor
    int CommonAncestor (-1);
    bool Found (false);
    for (size_t i=0; i<Parents_Mu1.size(); ++i){
        for (size_t j=0; j<Parents_Mu2.size(); ++j){
            if (Parents_Mu1.at(i) == -1 or Parents_Mu2.at(j) == -1) continue;
            if (Parents_Mu1.at(i) == Parents_Mu2.at(j)){
                CommonAncestor = Parents_Mu1.at(i);
                Found = true;
                break;
            }
        }
        if (Found) break;
    }
    return CommonAncestor;
}


//Reconstruct the whole decay chain from the common ancestor
ROOT::VecOps::RVec<int> Find_MC_CommonAncestor_Daughters(int CA_ind,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> Particle,
                                                         ROOT::VecOps::RVec<int> Daughters_ind){

    ROOT::VecOps::RVec<int> results;
    
    //Check if CA exists
    if (CA_ind == -1){
        results.push_back(-1);
        return results;
    }                                                       

    else{
        //Init with the decay of the common ancestor
        results.push_back(CA_ind);
        for (size_t k=Particle.at(CA_ind).daughters_begin; k<Particle.at(CA_ind).daughters_end; ++k) {
            results.push_back(Daughters_ind.at(k));
        }
        results.push_back(-2); //Split subdecays

        //Since the indices will not be unique in the list, keep track of those already read
        ROOT::VecOps::RVec<int> bookkeeping;
        bookkeeping.push_back(CA_ind);

        int i (1); //Start after the common ancestor
        do {
            //First check if already scanned
            bool InBook (false);
            for (size_t m=0; m<bookkeeping.size(); ++m){
                if (results.at(i) == bookkeeping.at(m)){ 
                    InBook = true;
                    break;
                }
            }

            //Do the scan by adding first the mother, then the daughters, then the splitter
            if (results.at(i)>0 && not InBook){
                    results.push_back(results.at(i));
                    bookkeeping.push_back(results.at(i));
                    for (size_t j=Particle.at(results.at(i)).daughters_begin; j<Particle.at(results.at(i)).daughters_end; ++j){
                        results.push_back(Daughters_ind.at(j));
                    }
                    results.push_back(-2);
            }

            ++i;
        } while(i<results.size());
        
        return results;
    } //Output as Mother,Daughter1,Daughter2,-2,Mother1(Daughter1),Daughter11,Daughter12,-2,...
}

}}
#endif
