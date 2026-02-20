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

//Get the exclusive Tau2MuNuNu Tau decay in the Bs2TauTau decays WARNING: PM = 1 GIVES TAU_MINUS WHILE -1 GIVES TAU_PLUS
ROOT::VecOps::RVec<edm4hep::MCParticleData> Find_genBs2TauTau_Muons(ROOT::VecOps::RVec<edm4hep::MCParticleData> genBs2TauTau_list, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> daughter, int pm){
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    for (size_t i = 0; i < genBs2TauTau_list.size(); ++i){
        if (genBs2TauTau_list[i].PDG == pm*15){
            ROOT::VecOps::RVec<int> Muons;
            ROOT::VecOps::RVec<int> Rem;
            for (size_t j = genBs2TauTau_list[i].daughters_begin; j < genBs2TauTau_list[i].daughters_end; ++j){
                if (std::abs(in[daughter.at(j)].PDG) == 13){
                    Muons.push_back(daughter.at(j));
                }
                else {
                    Rem.push_back(daughter.at(j));
                }
            }
            if (Muons.size() == 1){
                for (size_t m = 0; m < Muons.size(); ++m){
                    result.push_back(in.at(Muons[m]));
                }
                for (size_t n = 0; n < Rem.size(); ++n){
                    result.push_back(in.at(Rem[n]));
                }
            }
        }
    }
    return result; //Structure (muon, nu0, nu1, photon0, photon1, ...)
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


}}
#endif
