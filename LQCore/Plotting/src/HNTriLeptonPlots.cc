#include "HNTriLeptonPlots.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

HNTriLeptonPlots::HNTriLeptonPlots(TString name): StdPlots(name){

  TH1::SetDefaultSumw2(true);
  map_sig["h_Nevents"] = new TH1D("h_Nevents_"+name, "", 1, 0., 1.);
  map_sig["h_Njets"]                  =     new TH1D("h_Njets_"             + name,"number of jets",10,0,10);
  map_sig["h_Nbjets"]                 =     new TH1D("h_Nbjets_"            + name,"number of b jets",5,0,5);
  map_sig["h_jets_pt"]                =     new TH1D("h_jets_pt_"           + name,"jet pt",60,0,300);
  map_sig["h_jets_eta"]               =     new TH1D("h_jets_eta_"          + name,"#eta distribution of the two jets",120,-3,3);

  map_sig["h_osllmass"]               =     new TH1D("h_osllmass_"          + name,"Invariant mass of the two leading os electrons",100,0,500);
  map_sig["h_llljjmass"]               =     new TH1D("h_llljjmass_"         + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l1jjmass"]               =     new TH1D("h_l1jjmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);
  map_sig["h_l2jjmass"]               =     new TH1D("h_l2jjmass_"          + name,"Invariant mass of the two leading jets and second electron",100,0,1000);
  map_sig["h_dijetmass"]               =     new TH1D("h_dijetmass_"          + name,"Invariant mass of the two leading jets" ,100,0,1000);
  map_sig["h_mlll"] = new TH1D("h_mlll_"+name, "", 1000, 0., 1000.);
  map_sig["h_mT"] = new TH1D("h_mT_"+name, "", 1000, 0., 1000.);

  map_sig["h_LeptonPt"]               =     new TH1D("h_LeptonPt_"          + name,"lepton pt",100,0,500);
  map_sig["h_leadingLepton_Pt"]        =     new TH1D("h_leadingLepton_Pt_"   + name,"leading lepton pt",100,0,500);
  map_sig["h_secondLepton_Pt"]         =     new TH1D("h_secondLepton_Pt_"    + name,"secondary lepton pt",60,0,300);
  map_sig["h_thirdLepton_Pt"]         =     new TH1D("h_thirdLepton_Pt_"    + name,"thirdary lepton pt",60,0,300);

  map_sig["h_LeptonEta"]              =     new TH1D("h_LeptonEta_"         + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_leadingLepton_Eta"]       =     new TH1D("h_leadingLepton_Eta_"  + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_secondLepton_Eta"]        =     new TH1D("h_secondLepton_Eta_"   + name,"second lepton eta",60,-3.,3.);
  map_sig["h_thirdLepton_Eta"]        =     new TH1D("h_thirdLepton_Eta_"   + name,"third lepton eta",60,-3.,3.);
  
  //==== id variables
  map_sig["h_LeptonRelIso"] = new TH1D("h_LeptonRelIso_"+name, "lepton LeptonRelIso04", 100, 0, 0.5);
  map_sig["h_leadingLepton_LeptonRelIso"] = new TH1D("h_leadingLepton_LeptonRelIso_"+name, "leading lepton LeptonRelIso04", 100, 0, 0.5);
  map_sig["h_secondLepton_LeptonRelIso"] = new TH1D("h_secondLepton_LeptonRelIso_"+name, "secondary lepton LeptonRelIso04", 100, 0, 0.5);
  map_sig["h_thirdLepton_LeptonRelIso"] = new TH1D("h_thirdLepton_LeptonRelIso_"+name, "thirdary lepton LeptonRelIso04", 100, 0, 0.5);

  map_sig["h_dXY"] = new TH1D("h_dXY_"+name, "lepton dXY", 1000, 0, 0.1);
  map_sig["h_leadingLepton_dXY"] = new TH1D("h_leadingLepton_dXY_"+name, "leading lepton dXY", 1000, 0, 0.1);
  map_sig["h_secondLepton_dXY"] = new TH1D("h_secondLepton_dXY_"+name, "secondary lepton dXY", 1000, 0, 0.1);
  map_sig["h_thirdLepton_dXY"] = new TH1D("h_thirdLepton_dXY_"+name, "thirdary lepton dXY", 1000, 0, 0.1);

  map_sig["h_dXYSig"] = new TH1D("h_dXYSig_"+name, "lepton dXYSig", 80, 0., 8.);
  map_sig["h_leadingLepton_dXYSig"] = new TH1D("h_leadingLepton_dXYSig_"+name, "leading lepton dXYSig", 80, 0., 8.);
  map_sig["h_secondLepton_dXYSig"] = new TH1D("h_secondLepton_dXYSig_"+name, "secondary lepton dXYSig", 80, 0., 8.);
  map_sig["h_thirdLepton_dXYSig"] = new TH1D("h_thirdLepton_dXYSig_"+name, "thirdary lepton dXYSig", 80, 0., 8.);

  map_sig["h_dZ"] = new TH1D("h_dZ_"+name, "lepton dZ", 5000, 0, 0.5);
  map_sig["h_leadingLepton_dZ"] = new TH1D("h_leadingLepton_dZ_"+name, "leading lepton dZ", 5000, 0, 0.5);
  map_sig["h_secondLepton_dZ"] = new TH1D("h_secondLepton_dZ_"+name, "secondary lepton dZ", 5000, 0, 0.5);
  map_sig["h_thirdLepton_dZ"] = new TH1D("h_thirdLepton_dZ_"+name, "thirdary lepton dZ", 5000, 0, 0.5);

  map_sig["h_GlobalChi2"] = new TH1D("h_GlobalChi2_"+name, "lepton GlobalChi2", 100, 0, 100);
  map_sig["h_leadingLepton_GlobalChi2"] = new TH1D("h_leadingLepton_GlobalChi2_"+name, "leading lepton GlobalChi2", 100, 0, 100);
  map_sig["h_secondLepton_GlobalChi2"] = new TH1D("h_secondLepton_GlobalChi2_"+name, "secondary lepton GlobalChi2", 100, 0, 100);
  map_sig["h_thirdLepton_GlobalChi2"] = new TH1D("h_thirdLepton_GlobalChi2_"+name, "thirdary_lepton GlobalChi2", 100, 0, 100);

  map_sig["h_MuonType"] = new TH1D("h_MuonType_"+name, "lepton MuonType", 100, 0, 100);
  map_sig["h_leadingLepton_MuonType"] = new TH1D("h_leadingLepton_MuonType_"+name, "leading lepton MuonType", 100, 0, 100);
  map_sig["h_secondLepton_MuonType"] = new TH1D("h_secondLepton_MuonType_"+name, "secondary lepton MuonType", 100, 0, 100);
  map_sig["h_thirdLepton_MuonType"] = new TH1D("h_thirdLepton_MuonType_"+name, "thirdary_lepton MuonType", 100, 0, 100);

  map_sig["h_PFMET"]                  =     new TH1D("h_PFMET_"               + name,"Missing Et",100,0.0,500.0);
  map_sig["h_PFMET_phi"]              =     new TH1D("h_PFMET_phi_"           + name,"Missing Et",100,-3.2,3.2);


  map_sig["h_nVertices"]              =     new TH1D("h_nVertices_"         + name,"number of even vertices",60,0.0,60.0);
  /// Charge plot
  map_sig["h_sumcharge"]              =     new TH1D("h_sumcharge_"         + name,"Charge of the electron pair",6,-3,3);

  /// Number of objects
  map_sig["h_Nmuons"]                 =     new TH1D("h_Nmuons_"           + name,"number of mu",5,0,5);
  map_sig["h_Nelectrons"]             =     new TH1D("h_Nelectrons_"           + name,"number of el",5,0,5);

  map_sig["h_HN_mass_class1"] = new TH1D("h_HN_mass_class1_"+name, "", 2000, 0., 2000.);
  map_sig["h_HN_mass_class2"] = new TH1D("h_HN_mass_class2_"+name, "", 2000, 0., 2000.);
  map_sig["h_HN_mass_class3"] = new TH1D("h_HN_mass_class3_"+name, "", 2000, 0., 2000.);
  map_sig["h_HN_mass_class4"] = new TH1D("h_HN_mass_class4_"+name, "", 2000, 0., 2000.);
  map_sig["h_W_pri_lowmass_mass"] = new TH1D("h_W_pri_lowmass_mass_"+name, "", 2000, 0., 2000.);
  map_sig["h_W_pri_highmass_mass"] = new TH1D("h_W_pri_highmass_mass_"+name, "", 2000, 0., 2000.);
  map_sig["h_W_sec_highmass_mass"] = new TH1D("h_W_sec_highmass_mass_"+name, "", 2000, 0., 2000.);
  map_sig["h_deltaR_OS_min"] = new TH1D("h_deltaR_OS_min_"+name, "", 50, 0., 5.);
  map_sig["h_deltaR_OS"] = new TH1D("h_deltaR_OS_"+name, "", 50, 0., 5.);
  map_sig["h_deltaR_ll"] = new TH1D("h_deltaR_ll_"+name, "", 50, 0., 5.);
  map_sig["h_gamma_star_mass"] = new TH1D("h_gamma_star_mass_"+name, "", 200, 0., 200.);
  map_sig["h_z_candidate_mass"] = new TH1D("h_z_candidate_mass_"+name, "", 200, 0., 200.);




  MET = -999.;
  METphi = -999.;
  NBjet = -999;

}




void HNTriLeptonPlots::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, Double_t weight) {

  Fill("h_Nevents", 0., weight);
  Fill("h_Nmuons" ,muons.size(), weight);
  Fill("h_Nelectrons" ,electrons.size(), weight);

  Fill("h_HN_mass_class1", HN[0].M(), weight);
  Fill("h_HN_mass_class2", HN[1].M(), weight);
  Fill("h_HN_mass_class3", HN[2].M(), weight);
  Fill("h_HN_mass_class4", HN[3].M(), weight);
  Fill("h_W_pri_lowmass_mass", W_pri_lowmass.M(), weight);
  Fill("h_W_pri_highmass_mass", W_pri_highmass.M(), weight);
  Fill("h_W_sec_highmass_mass", W_sec.M(), weight);

  double thismet(MET), thismetphi(METphi);
  if(MET    == -999.) thismet = ev.PFMET();
  if(METphi == -999.) thismetphi = ev.METPhi(snu::KEvent::pfmet);

  Fill("h_PFMET",thismet, weight);
  Fill("h_PFMET_phi",thismetphi, weight);
  Fill("h_nVertices", ev.nVertices(), weight);

  if(muons.size()==3){

    Fill("h_deltaR_ll", muons[0].DeltaR(muons[1]),weight);
    Fill("h_deltaR_ll", muons[0].DeltaR(muons[2]),weight);
    Fill("h_deltaR_ll", muons[1].DeltaR(muons[2]),weight);

    Fill("h_deltaR_OS", muons[OppSign].DeltaR(muons[SameSign[0]]),weight);
    Fill("h_deltaR_OS", muons[OppSign].DeltaR(muons[SameSign[1]]),weight);

    snu::KParticle gamma_star, z_candidate;
    double deltaR_OS_min;
    if( muons[OppSign].DeltaR(muons[SameSign[0]]) < muons[OppSign].DeltaR(muons[SameSign[1]]) ){
      deltaR_OS_min = muons[OppSign].DeltaR(muons[SameSign[0]]);
      gamma_star = muons[OppSign] + muons[SameSign[0]];
    }
    else{
      deltaR_OS_min = muons[OppSign].DeltaR(muons[SameSign[1]]);
      gamma_star = muons[OppSign] + muons[SameSign[1]];
    }

    if( fabs( (muons[OppSign] + muons[SameSign[0]]).M() - 91.1876 ) <
        fabs( (muons[OppSign] + muons[SameSign[1]]).M() - 91.1876 )   ){
      z_candidate = muons[OppSign] + muons[SameSign[0]];
    }
    else{
      z_candidate = muons[OppSign] + muons[SameSign[1]];
    }

    Fill("h_deltaR_OS_min", deltaR_OS_min, weight);
    Fill("h_gamma_star_mass", gamma_star.M(), weight);
    Fill("h_z_candidate_mass", z_candidate.M(), weight);

    snu::KParticle nu;
    nu.SetPxPyPzE(thismet*TMath::Cos(thismetphi), thismet*TMath::Sin(thismetphi), 0, thismet);

    if(muons[0].Charge() != muons[1].Charge())    Fill("h_osllmass", (muons[0]+muons[1]).M(),weight);
    if(muons[0].Charge() != muons[2].Charge())    Fill("h_osllmass", (muons[0]+muons[2]).M(),weight);
    if(muons[1].Charge() != muons[2].Charge())    Fill("h_osllmass", (muons[1]+muons[2]).M(),weight);
    Fill("h_mlll", (muons[0]+muons[1]+muons[2]).M(), weight);
    int imu=0;
    for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++, imu++){
      double LeptonRelIso = muit->RelIso04();
      Fill("h_LeptonPt", muit->Pt(),weight);
      Fill("h_LeptonEta",muit->Eta(),weight);
      Fill("h_LeptonRelIso", LeptonRelIso, weight);
      Fill("h_dXY", fabs(muit->dXY()), weight);
      Fill("h_dXYSig", fabs(muit->dXYSig()), weight);
      Fill("h_dZ", muit->dZ(), weight);
      Fill("h_GlobalChi2", muit->GlobalChi2(), weight);
      Fill("h_mT", MT(*muit, nu), weight);
      Fill("h_MuonType", muit->GetType(), weight);
      if(imu ==0) {
        Fill("h_leadingLepton_Pt", muit->Pt(),weight);
        Fill("h_leadingLepton_Eta",muit->Eta(),weight);
        Fill("h_leadingLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_leadingLepton_dXY", fabs(muit->dXY()), weight);
        Fill("h_leadingLepton_dXYSig", fabs(muit->dXYSig()), weight);
        Fill("h_leadingLepton_dZ", muit->dZ(), weight);
        Fill("h_leadingLepton_GlobalChi2", muit->GlobalChi2(), weight);
        Fill("h_leadingLepton_MuonType", muit->GetType(), weight);
      }
      if(imu ==1) {
        Fill("h_secondLepton_Pt", muit->Pt(),weight);
        Fill("h_secondLepton_Eta",muit->Eta(),weight);
        Fill("h_secondLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_secondLepton_dXY", fabs(muit->dXY()), weight);
        Fill("h_secondLepton_dXYSig", fabs(muit->dXYSig()), weight);
        Fill("h_secondLepton_dZ", muit->dZ(), weight);
        Fill("h_secondLepton_GlobalChi2", muit->GlobalChi2(), weight);
        Fill("h_secondLepton_MuonType", muit->GetType(), weight);
      }
      if(imu ==2) {
        Fill("h_thirdLepton_Pt", muit->Pt(),weight);
        Fill("h_thirdLepton_Eta",muit->Eta(),weight);
        Fill("h_thirdLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_thirdLepton_dXY", fabs(muit->dXY()), weight);
        Fill("h_thirdLepton_dXYSig", fabs(muit->dXYSig()), weight);
        Fill("h_thirdLepton_dZ", muit->dZ(), weight);
        Fill("h_thirdLepton_GlobalChi2", muit->GlobalChi2(), weight);
        Fill("h_thirdLepton_MuonType", muit->GetType(), weight);
      }
      
    }

  }

  if(muons.size() == 2 && electrons.size()==1){

    Fill("h_mlll", (muons[0]+muons[1]+electrons[0]).M(), weight);

    Fill("h_deltaR_ll", muons[0].DeltaR(electrons[0]),weight);
    Fill("h_deltaR_ll", muons[1].DeltaR(electrons[0]),weight);
    Fill("h_deltaR_ll", muons[0].DeltaR(muons[1]),weight);

    int imu=0;
    for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++, imu++){
      double LeptonRelIso = muit->RelIso04();
      Fill("h_LeptonPt", muit->Pt(),weight);
      Fill("h_LeptonEta",muit->Eta(),weight);
      Fill("h_LeptonRelIso", LeptonRelIso, weight);
      Fill("h_dXY", fabs(muit->dXY()), weight);
      Fill("h_dXYSig", fabs(muit->dXYSig()), weight);
      Fill("h_dZ", muit->dZ(), weight);
      Fill("h_GlobalChi2", muit->GlobalChi2(), weight);
      if(imu ==0) {
        Fill("h_leadingLepton_Pt", muit->Pt(),weight);
        Fill("h_leadingLepton_Eta",muit->Eta(),weight);
        Fill("h_leadingLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_leadingLepton_dXY", fabs(muit->dXY()), weight);
        Fill("h_leadingLepton_dXYSig", fabs(muit->dXYSig()), weight);
        Fill("h_leadingLepton_dZ", muit->dZ(), weight);
        Fill("h_leadingLepton_GlobalChi2", muit->GlobalChi2(), weight);
      }
      if(imu ==1) {
        Fill("h_secondLepton_Pt", muit->Pt(),weight);
        Fill("h_secondLepton_Eta",muit->Eta(),weight);
        Fill("h_secondLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_secondLepton_dXY", fabs(muit->dXY()), weight);
        Fill("h_secondLepton_dXYSig", fabs(muit->dXYSig()), weight);
        Fill("h_secondLepton_dZ", muit->dZ(), weight);
        Fill("h_secondLepton_GlobalChi2", muit->GlobalChi2(), weight);
      }
    }
    Fill("h_thirdLepton_Pt", electrons[0].Pt(),weight);
    Fill("h_thirdLepton_Eta",electrons[0].Eta(),weight);
    Fill("h_thirdLepton_LeptonRelIso", electrons[0].PFRelIso(0.3), weight);
    Fill("h_thirdLepton_dXY", fabs(electrons[0].dxy()), weight);
    Fill("h_thirdLepton_dXYSig", fabs(electrons[0].dxySig()), weight);
    Fill("h_thirdLepton_dZ", electrons[0].dz(), weight);

  }
  
  Fill("h_Njets",jets.size(), weight);
  int nbjet(NBjet);
  if(NBjet==-999){
    int nbjet=0;
    for(UInt_t j=0; j < jets.size(); j++){
      if(jets.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) nbjet++;
    }
  }
  Fill("h_Nbjets", nbjet, weight); 

  for(UInt_t j=0; j < jets.size(); j++){
    Fill("h_jets_pt", jets[j].Pt(),weight);
    Fill("h_jets_eta",jets[j].Eta(),weight);
  }


  return;
}/// End of Fill



void HNTriLeptonPlots::Write() {
 
  for(map<TString, TH1*>::iterator it = map_sig.begin(); it != map_sig.end(); it++){
    it->second->Write();
  }

  for(map<TString, TH2*>::iterator it = map_sig2.begin(); it != map_sig2.end(); it++){
    it->second->Write();
  }


  for(map<TString, TH3*>::iterator it = map_sig3.begin(); it != map_sig3.end(); it++){
    it->second->Write();
  }

}


HNTriLeptonPlots::HNTriLeptonPlots(): StdPlots(){
}


/**
 * Copy constructor.
 */
HNTriLeptonPlots::HNTriLeptonPlots(const HNTriLeptonPlots& sp): StdPlots(sp)
{
  for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
    std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
    std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
    std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
    mit->second = mit2->second;
  }
  
}


HNTriLeptonPlots& HNTriLeptonPlots::operator= (const HNTriLeptonPlots& sp)
{
  if (this != &sp) {

    for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
      std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
      mit->second = mit2->second;
    }
    
    for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
      std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
      mit->second = mit2->second;
    }

    for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
      std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
      mit->second = mit2->second;
    }
  }
  return *this;
}

HNTriLeptonPlots::~HNTriLeptonPlots() {
   for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
     delete mit->second ;
  }

   for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
     delete mit->second ;
   }

   for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
     delete mit->second ;
   }
   
}

void HNTriLeptonPlots::Fill(TString name, double value, double w){
  std::map<TString, TH1*>::iterator it = map_sig.find(name);
  if(it!= map_sig.end())   it->second->Fill(value, w);

  else cout << name << " not found in map_sig" << endl;
  return;
}
 
void HNTriLeptonPlots::Fill(TString name, double value1, double value2, double w){
   std::map<TString, TH2*>::iterator it = map_sig2.find(name);
   if(it!= map_sig2.end()) it->second->Fill(value1, value2, w);
   else cout << name << " not found in map_sig" << endl;
   return;
 }



void HNTriLeptonPlots::Fill(TString name, double value1, double value2, double value3, double w){
  std::map<TString, TH3*>::iterator it = map_sig3.find(name);
  if(it!= map_sig3.end()) it->second->Fill(value1, value2, value3, w);
  else cout << name << " not found in map_sig" << endl;
  return;
}

void HNTriLeptonPlots::SetMetInfo(double met, double metphi){
  //cout << "[HNTriLeptonPlots::SetMetInfo] MET set" << endl;
  MET = met;
  METphi = metphi;
}

void HNTriLeptonPlots::SetParticleInfo(snu::KParticle *hn, snu::KParticle w_pri_lowmass, snu::KParticle w_pri_highmass, snu::KParticle w_sec){

  HN = hn;
  W_pri_lowmass = w_pri_lowmass;
  W_pri_highmass = w_pri_highmass;
  W_sec = w_sec;

}

void HNTriLeptonPlots::SetChargeSign(double os, double ss0, double ss1){

  OppSign = os;
  SameSign[0] = ss0;
  SameSign[1] = ss1;

}

void HNTriLeptonPlots::SetBJet(int nbjet){

  NBjet = nbjet;

}


double HNTriLeptonPlots::MT(TLorentzVector a, TLorentzVector b){
  double dphi = a.DeltaPhi(b);
  return TMath::Sqrt( 2.*a.Pt()*b.Pt()*(1.- TMath::Cos(dphi) ) );

}




