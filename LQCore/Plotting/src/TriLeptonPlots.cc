#include "TriLeptonPlots.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

TriLeptonPlots::TriLeptonPlots(TString name) : isPlotsFilled(false)
{

  TH1::SetDefaultSumw2(true);
  map_sig["h_Njets"]                  =     new TH1F("h_Njets_"             + name,"number of jets",10,0,10);
  map_sig["h_Nbjets"]                 =     new TH1F("h_Nbjets_"            + name,"number of b jets",5,0,5);
  map_sig["h_jets_pt"]                =     new TH1F("h_jets_pt_"           + name,"jet pt",60,0,300);
  map_sig["h_jets_eta"]               =     new TH1F("h_jets_eta_"          + name,"#eta distribution of the two jets",120,-3,3);

  map_sig["h_osllmass"]               =     new TH1F("h_osllmass_"          + name,"Invariant mass of the two leading os electrons",100,0,500);
  map_sig["h_llljjmass"]               =     new TH1F("h_llljjmass_"         + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l1jjmass"]               =     new TH1F("h_l1jjmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);
  map_sig["h_l2jjmass"]               =     new TH1F("h_l2jjmass_"          + name,"Invariant mass of the two leading jets and second electron",100,0,1000);
  map_sig["h_dijetmass"]               =     new TH1F("h_dijetmass_"          + name,"Invariant mass of the two leading jets" ,100,0,1000);

  map_sig["h_Lepton_Eta"]              =     new TH1F("h_Lepton_Eta_"         + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_Lepton_Pt"]               =     new TH1F("h_Lepton_Pt_"          + name,"lepton pt",100,0,500);
  map_sig["h_leadingLepton_Pt"]        =     new TH1F("h_leadingLepton_Pt_"   + name,"leading lepton pt",100,0,500);
  map_sig["h_secondLepton_Pt"]         =     new TH1F("h_secondLepton_Pt_"    + name,"secondary lepton pt",60,0,300);
  map_sig["h_thirdLepton_Pt"]         =     new TH1F("h_thirdLepton_Pt_"    + name,"thirdary lepton pt",60,0,300);
  map_sig["h_leadingLepton_Eta"]       =     new TH1F("h_leadingLepton_Eta_"  + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_secondLepton_Eta"]        =     new TH1F("h_secondLepton_Eta_"   + name,"second lepton eta",60,-3.,3.);
  map_sig["h_thirdLepton_Eta"]        =     new TH1F("h_thirdLepton_Eta_"   + name,"third lepton eta",60,-3.,3.);
  

  map_sig["h_PFMET"]                  =     new TH1F("h_PFMET_"               + name,"Missing Et",100,0.0,500.0);
  map_sig["h_PFMET_phi"]              =     new TH1F("h_PFMET_phi_"           + name,"Missing Et",100,0.0,500.0);
  map_sig["h_HT"]                     =     new TH1F("h_HT_"           + name,"H_{T}",300,0.0,300.0);





  map_sig["h_nVertices"]              =     new TH1F("h_nVertices_"         + name,"number of even vertices",60,0.0,60.0);
  /// Charge plot
  map_sig["h_sumcharge"]              =     new TH1F("h_sumcharge_"         + name,"Charge of the electron pair",6,-3,3);

  /// Number of objects
  map_sig["h_Nmuons"]                 =     new TH1F("h_Nmuons_"           + name,"number of mu",5,0,5);
  map_sig["h_Nelectrons"]             =     new TH1F("h_Nelectrons_"           + name,"number of el",5,0,5);

  //==== id variables
  map_sig["h_LeptonRelIso"] = new TH1F("h_LeptonRelIso_"+name, "lepton LeptonRelIso03", 100, 0, 0.5);
  map_sig["h_leadingLepton_LeptonRelIso"] = new TH1F("h_leadingLepton_LeptonRelIso_"+name, "leading lepton LeptonRelIso03", 100, 0, 0.5);
  map_sig["h_secondLepton_LeptonRelIso"] = new TH1F("h_secondLepton_LeptonRelIso_"+name, "secondary lepton LeptonRelIso03", 100, 0, 0.5);
  map_sig["h_thirdLepton_LeptonRelIso"] = new TH1F("h_thirdLepton_LeptonRelIso_"+name, "thirdary lepton LeptonRelIso03", 100, 0, 0.5);
  map_sig["h_dXY"] = new TH1F("h_dXY_"+name, "lepton dXY", 1000, 0, 0.1);
  map_sig["h_leadingLepton_dXY"] = new TH1F("h_leadingLepton_dXY_"+name, "leading lepton dXY", 1000, 0, 0.1);
  map_sig["h_secondLepton_dXY"] = new TH1F("h_secondLepton_dXY_"+name, "secondary lepton dXY", 1000, 0, 0.1);
  map_sig["h_thirdLepton_dXY"] = new TH1F("h_thirdLepton_dXY_"+name, "thirdary lepton dXY", 1000, 0, 0.1);
  map_sig["h_dZ"] = new TH1F("h_dZ_"+name, "lepton dZ", 5000, 0, 0.5);
  map_sig["h_leadingLepton_dZ"] = new TH1F("h_leadingLepton_dZ_"+name, "leading lepton dZ", 5000, 0, 0.5);
  map_sig["h_secondLepton_dZ"] = new TH1F("h_secondLepton_dZ_"+name, "secondary lepton dZ", 5000, 0, 0.5);
  map_sig["h_thirdLepton_dZ"] = new TH1F("h_thirdLepton_dZ_"+name, "thirdary lepton dZ", 5000, 0, 0.5);
  map_sig["h_GlobalChi2"] = new TH1F("h_GlobalChi2_"+name, "lepton GlobalChi2", 100, 0, 100);
  map_sig["h_leadingLepton_GlobalChi2"] = new TH1F("h_leadingLepton_GlobalChi2_"+name, "leading lepton GlobalChi2", 100, 0, 100);
  map_sig["h_secondLepton_GlobalChi2"] = new TH1F("h_secondLepton_GlobalChi2_"+name, "secondary lepton GlobalChi2", 100, 0, 100);
  map_sig["h_thirdLepton_GlobalChi2"] = new TH1F("h_thirdLepton_GlobalChi2_"+name, "thirdary_lepton GlobalChi2", 100, 0, 100);

  //==== TTL 
  map_sig["h_TTL_Lepton_Pt"]               =     new TH1F("h_TTL_Lepton_Pt_"          + name,"TTL lepton pt",100,0,500);
  map_sig["h_TTL_leadingLepton_Pt"]        =     new TH1F("h_TTL_leadingLepton_Pt_"   + name,"TTL leading lepton pt",100,0,500);
  map_sig["h_TTL_secondLepton_Pt"]         =     new TH1F("h_TTL_secondLepton_Pt_"    + name,"TTL secondary lepton pt",60,0,300);
  map_sig["h_TTL_thirdLepton_Pt"]          =     new TH1F("h_TTL_thirdLepton_Pt_"     + name,"TTL thirdary lepton pt",60,0,300);
  map_sig["h_TTL_Lepton_Eta"]              =     new TH1F("h_TTL_Lepton_Eta_"         + name,"TTL leading lepton eta",60,-3.,3.);
  map_sig["h_TTL_leadingLepton_Eta"]       =     new TH1F("h_TTL_leadingLepton_Eta_"  + name,"TTL leading lepton eta",60,-3.,3.);
  map_sig["h_TTL_secondLepton_Eta"]        =     new TH1F("h_TTL_secondLepton_Eta_"   + name,"TTL second lepton eta",60,-3.,3.);
  map_sig["h_TTL_thirdLepton_Eta"]         =     new TH1F("h_TTL_thirdLepton_Eta_"    + name,"TTL third lepton eta",60,-3.,3.);
  //==== TLL
  map_sig["h_TLL_Lepton_Pt"]               =     new TH1F("h_TLL_Lepton_Pt_"          + name,"TLL lepton pt",100,0,500);
  map_sig["h_TLL_leadingLepton_Pt"]        =     new TH1F("h_TLL_leadingLepton_Pt_"   + name,"TLL leading lepton pt",100,0,500);
  map_sig["h_TLL_secondLepton_Pt"]         =     new TH1F("h_TLL_secondLepton_Pt_"    + name,"TLL secondary lepton pt",60,0,300);
  map_sig["h_TLL_thirdLepton_Pt"]          =     new TH1F("h_TLL_thirdLepton_Pt_"     + name,"TLL thirdary lepton pt",60,0,300);
  map_sig["h_TLL_Lepton_Eta"]              =     new TH1F("h_TLL_Lepton_Eta_"         + name,"TLL leading lepton eta",60,-3.,3.);
  map_sig["h_TLL_leadingLepton_Eta"]       =     new TH1F("h_TLL_leadingLepton_Eta_"  + name,"TLL leading lepton eta",60,-3.,3.);
  map_sig["h_TLL_secondLepton_Eta"]        =     new TH1F("h_TLL_secondLepton_Eta_"   + name,"TLL second lepton eta",60,-3.,3.);
  map_sig["h_TLL_thirdLepton_Eta"]         =     new TH1F("h_TLL_thirdLepton_Eta_"    + name,"TLL third lepton eta",60,-3.,3.);
  //==== LLL
  map_sig["h_LLL_Lepton_Pt"]               =     new TH1F("h_LLL_Lepton_Pt_"          + name,"LLL lepton pt",100,0,500);
  map_sig["h_LLL_leadingLepton_Pt"]        =     new TH1F("h_LLL_leadingLepton_Pt_"   + name,"LLL leading lepton pt",100,0,500);
  map_sig["h_LLL_secondLepton_Pt"]         =     new TH1F("h_LLL_secondLepton_Pt_"    + name,"LLL secondary lepton pt",60,0,300);
  map_sig["h_LLL_thirdLepton_Pt"]          =     new TH1F("h_LLL_thirdLepton_Pt_"     + name,"LLL thirdary lepton pt",60,0,300);
  map_sig["h_LLL_Lepton_Eta"]              =     new TH1F("h_LLL_Lepton_Eta_"         + name,"LLL leading lepton eta",60,-3.,3.);
  map_sig["h_LLL_leadingLepton_Eta"]       =     new TH1F("h_LLL_leadingLepton_Eta_"  + name,"LLL leading lepton eta",60,-3.,3.);
  map_sig["h_LLL_secondLepton_Eta"]        =     new TH1F("h_LLL_secondLepton_Eta_"   + name,"LLL second lepton eta",60,-3.,3.);
  map_sig["h_LLL_thirdLepton_Eta"]         =     new TH1F("h_LLL_thirdLepton_Eta_"    + name,"LLL third lepton eta",60,-3.,3.);

  map_sig["h_n_tight"]         =     new TH1F("h_n_tight_"    + name,"# of Tight Muons", 4, 0.,4.);

}




void TriLeptonPlots::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, Double_t weight) {

  isPlotsFilled = true;
 
  Fill("h_Nmuons" ,muons.size(), weight);
  Fill("h_Nelectrons" ,electrons.size(), weight);
  //if(! (muons.size()==3 || electrons.size()==3 || (electrons.size() == 2 && muons.size() ==1) || (electrons.size() == 1 && muons.size() ==2) )) return;
  //if(electrons.size() == 3 && muons.size() > 0) return;
  //if(electrons.size() > 0 && muons.size() ==3) return;

  Fill("h_Nmuons" ,muons.size(), weight);
  Fill("h_Nelectrons" ,electrons.size(), weight);


  if(electrons.size()==2){
    if(electrons[0].Charge() != electrons[1].Charge())   Fill("h_osllmass", (electrons[0]+electrons[1]).M(),weight);
  }
  if(muons.size()==2){
    if(muons[0].Charge() != muons[1].Charge())   Fill("h_osllmass", (muons[0]+muons[1]).M(),weight);
  }

  if(muons.size()==3){
    if(muons[0].Charge() != muons[1].Charge())    Fill("h_osllmass", (muons[0]+muons[1]).M(),weight);
    if(muons[0].Charge() != muons[2].Charge())    Fill("h_osllmass", (muons[0]+muons[2]).M(),weight);
    if(muons[1].Charge() != muons[2].Charge())    Fill("h_osllmass", (muons[1]+muons[2]).M(),weight);

    int imu=0, n_tight=0;

    for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++){
      float LeptonRelIso = (muit->SumIsoCHDR03() + std::max(0.0, muit->SumIsoNHDR03() + muit->SumIsoPHDR03() - 0.5* muit->SumPUIsoR03()))/muit->Pt() ;
      if(LeptonRelIso<0.1) n_tight++;
    }
    Fill("h_n_tight", n_tight, 1);

    for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++, imu++){
      float LeptonRelIso = (muit->SumIsoCHDR03() + std::max(0.0, muit->SumIsoNHDR03() + muit->SumIsoPHDR03() - 0.5* muit->SumPUIsoR03()))/muit->Pt() ;
      Fill("h_Lepton_Pt", muit->Pt(),weight);
      Fill("h_Lepton_Eta",muit->Eta(),weight);
      Fill("h_LeptonRelIso", LeptonRelIso, weight);
      Fill("h_dXY", muit->dXY(), weight);
      Fill("h_dZ", muit->dZ(), weight);
      Fill("h_GlobalChi2", muit->GlobalChi2(), weight);
  
      //==== TTL
      if(n_tight==2){
        Fill("h_TTL_Lepton_Pt", muit->Pt(),1);
        Fill("h_TTL_Lepton_Eta",muit->Eta(),1);
      }
      //==== TLL
      if(n_tight==1){
        Fill("h_TLL_Lepton_Pt", muit->Pt(),1);
        Fill("h_TLL_Lepton_Eta",muit->Eta(),1);
      }
      //==== LLL
      if(n_tight==0){
        Fill("h_LLL_Lepton_Pt", muit->Pt(),1);
        Fill("h_LLL_Lepton_Eta",muit->Eta(),1);
      }

      if(imu ==0) {
        Fill("h_leadingLepton_Pt", muit->Pt(),weight);
        Fill("h_leadingLepton_Eta",muit->Eta(),weight);
        Fill("h_leadingLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_leadingLepton_dXY", muit->dXY(), weight);
        Fill("h_leadingLepton_dZ", muit->dZ(), weight);
        Fill("h_leadingLepton_GlobalChi2", muit->GlobalChi2(), weight);

        //==== TTL
        if(n_tight==2){
          Fill("h_TTL_leadingLepton_Pt", muit->Pt(),1);
          Fill("h_TTL_leadingLepton_Eta",muit->Eta(),1);          
        }
        //==== TLL
        if(n_tight==1){
          Fill("h_TLL_leadingLepton_Pt", muit->Pt(),1);
          Fill("h_TLL_leadingLepton_Eta",muit->Eta(),1);
        }
        //==== LLL
        if(n_tight==0){
          Fill("h_LLL_leadingLepton_Pt", muit->Pt(),1);
          Fill("h_LLL_leadingLepton_Eta",muit->Eta(),1);
        }

      }
      if(imu ==1) {
        Fill("h_secondLepton_Pt", muit->Pt(),weight);
        Fill("h_secondLepton_Eta",muit->Eta(),weight);
        Fill("h_secondLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_secondLepton_dXY", muit->dXY(), weight);
        Fill("h_secondLepton_dZ", muit->dZ(), weight);
        Fill("h_secondLepton_GlobalChi2", muit->GlobalChi2(), weight);

        //==== TTL
        if(n_tight==2){
          Fill("h_TTL_secondLepton_Pt", muit->Pt(),1);
          Fill("h_TTL_secondLepton_Eta",muit->Eta(),1);
        }
        //==== TLL
        if(n_tight==1){
          Fill("h_TLL_secondLepton_Pt", muit->Pt(),1);
          Fill("h_TLL_secondLepton_Eta",muit->Eta(),1);
        }
        //==== LLL
        if(n_tight==0){
          Fill("h_LLL_secondLepton_Pt", muit->Pt(),1);
          Fill("h_LLL_secondLepton_Eta",muit->Eta(),1);
        }

      }
      if(imu ==2) {
        Fill("h_thirdLepton_Pt", muit->Pt(),weight);
        Fill("h_thirdLepton_Eta",muit->Eta(),weight);
        Fill("h_thirdLepton_LeptonRelIso", LeptonRelIso, weight);
        Fill("h_thirdLepton_dXY", muit->dXY(), weight);
        Fill("h_thirdLepton_dZ", muit->dZ(), weight);
        Fill("h_thirdLepton_GlobalChi2", muit->GlobalChi2(), weight);

        //==== TTL
        if(n_tight==2){
          Fill("h_TTL_thirdLepton_Pt", muit->Pt(),1);
          Fill("h_TTL_thirdLepton_Eta",muit->Eta(),1);
        }
        //==== TLL
        if(n_tight==1){
          Fill("h_TLL_thirdLepton_Pt", muit->Pt(),1);
          Fill("h_TLL_thirdLepton_Eta",muit->Eta(),1);
        }
        //==== LLL
        if(n_tight==0){
          Fill("h_LLL_thirdLepton_Pt", muit->Pt(),1);
          Fill("h_LLL_thirdLepton_Eta",muit->Eta(),1);
        }

      }
      
    }
    if(jets.size() >1){
      snu::KParticle jj = jets[0] + jets[1];
      snu::KParticle l1jj = muons[0] + jets[0] + jets[1];
      snu::KParticle l2jj = muons[1] + jets[0] + jets[1];
      snu::KParticle llljj = muons[0] + muons[1] + muons[2] + jets[0] + jets[1];
      Fill("h_l1jjmass", l1jj.M() ,weight);
      Fill("h_l2jjmass", l2jj.M() ,weight);
      Fill("h_llljjmass", llljj.M() ,weight);
      Fill("h_dijetmass", jj.M() ,weight);

    }
  }

  if(electrons.size()==3){
    if(electrons[0].Charge() != electrons[1].Charge())    Fill("h_osllmass", (electrons[0]+electrons[1]).M(),weight);
    if(electrons[0].Charge() != electrons[2].Charge())    Fill("h_osllmass", (electrons[0]+electrons[2]).M(),weight);
    if(electrons[1].Charge() != electrons[2].Charge())    Fill("h_osllmass", (electrons[1]+electrons[2]).M(),weight);
    int imu=0;
    for(std::vector<snu::KElectron>::iterator muit = electrons.begin(); muit != electrons.end(); muit++, imu++){
      Fill("h_Lepton_Pt", muit->Pt(),weight);
      Fill("h_Lepton_Eta",muit->Eta(),weight);
      if(imu ==0) {
        Fill("h_leadingLepton_Pt", muit->Pt(),weight);
        Fill("h_leadingLepton_Eta",muit->Eta(),weight);
      }
      if(imu ==1) {
        Fill("h_secondLepton_Pt", muit->Pt(),weight);
        Fill("h_secondLepton_Eta",muit->Eta(),weight);
      }
      if(imu ==2) {
        Fill("h_thirdLepton_Pt", muit->Pt(),weight);
        Fill("h_thirdLepton_Eta",muit->Eta(),weight);
      }

    }
    if(jets.size() >1){
      snu::KParticle jj = jets[0] + jets[1];
      snu::KParticle l1jj = electrons[0] + jets[0] + jets[1];
      snu::KParticle l2jj = electrons[1] + jets[0] + jets[1];
      snu::KParticle llljj = electrons[0] + electrons[1] + electrons[2] + jets[0] + jets[1];
      Fill("h_l1jjmass", l1jj.M() ,weight);
      Fill("h_l2jjmass", l2jj.M() ,weight);
      Fill("h_llljjmass", llljj.M() ,weight);
      Fill("h_dijetmass", jj.M() ,weight);

    }
 
  }
  
  int sum_charge=0;
  for(std::vector<snu::KElectron>::iterator muit = electrons.begin(); muit != electrons.end(); muit++){
    sum_charge+= muit->Charge();
  }
  for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++){
    sum_charge+= muit->Charge();
  }
  if(fabs(sum_charge)!=1){
    for(std::vector<snu::KElectron>::iterator muit = electrons.begin(); muit != electrons.end(); muit++){
      cout << "Electron charge = " <<  muit->Charge() << endl;
    }
    for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++){
      cout << "Muon Charge = " << muit->Charge() << endl;
    }
    cout << "Sum charge = " << sum_charge << endl;
  }
  
  Fill("h_sumcharge",sum_charge, weight);
  
 
  Fill("h_Njets",jets.size(), weight);
  int nbjet=0;
  for(UInt_t j=0; j < jets.size(); j++){
    if(jets.at(j).CombinedSecVertexBtag() > 0.679) nbjet++;
  }
  Fill("h_Nbjets", nbjet, weight); 

  Fill("h_PFMET",ev.PFMET(), weight);
  Fill("h_PFMET_phi",ev.PFMETphi(), weight);
  Fill("h_nVertices", ev.nVertices(), weight);

  float sumpt=0.;
  for(UInt_t j=0; j < jets.size(); j++){
    Fill("h_jets_pt", jets[j].Pt(),weight);
    Fill("h_jets_eta",jets[j].Eta(),weight);
    sumpt += jets[j].Pt();
  }
  Fill("h_HT", sumpt, weight);


  return;
}/// End of Fill



void TriLeptonPlots::Write() {

  if(!isPlotsFilled) return;
 
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


TriLeptonPlots::TriLeptonPlots(): StdPlots(){
}


/**
 * Copy constructor.
 */
TriLeptonPlots::TriLeptonPlots(const TriLeptonPlots& sp): StdPlots(sp)
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


TriLeptonPlots& TriLeptonPlots::operator= (const TriLeptonPlots& sp)
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

TriLeptonPlots::~TriLeptonPlots() {
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

void TriLeptonPlots::Fill(TString name, double value, double w){
  std::map<TString, TH1*>::iterator it = map_sig.find(name);
  if(it!= map_sig.end())   it->second->Fill(value, w);

  else cout << name << " not found in map_sig" << endl;
  return;
}
 
void TriLeptonPlots::Fill(TString name, double value1, double value2, double w){
   std::map<TString, TH2*>::iterator it = map_sig2.find(name);
   if(it!= map_sig2.end()) it->second->Fill(value1, value2, w);
   else cout << name << " not found in map_sig" << endl;
   return;
 }



void TriLeptonPlots::Fill(TString name, double value1, double value2, double value3, double w){
  std::map<TString, TH3*>::iterator it = map_sig3.find(name);
  if(it!= map_sig3.end()) it->second->Fill(value1, value2, value3, w);
  else cout << name << " not found in map_sig" << endl;
  return;
}




