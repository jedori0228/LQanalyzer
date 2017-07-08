// $Id: FakeRateCalculator_Mu_dxysig_DILEP.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_Mu_dxysig_DILEP Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_Mu_dxysig_DILEP.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_Mu_dxysig_DILEP);

FakeRateCalculator_Mu_dxysig_DILEP::FakeRateCalculator_Mu_dxysig_DILEP() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_Mu_dxysig_DILEP");
  
  Message("In FakeRateCalculator_Mu_dxysig_DILEP constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_Mu_dxysig_DILEP::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:
  //ResetLumiMask(snu::KEvent::gold);

  TDirectory* origDir = gDirectory;

  TString lqdir = getenv("LQANALYZER_DIR");
  TFile *file_tmp = new TFile( lqdir+"/JskimData/HalfSample/FR_sampleA.root" );

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }
  tempDir->cd();

  FR_sampleA = (TH2D*)file_tmp->Get("FR_sampleA")->Clone();

  file_tmp->Close();
  delete file_tmp;

  origDir->cd();

  return;
}


void FakeRateCalculator_Mu_dxysig_DILEP::ExecuteEvents()throw( LQError ){

  double thisdXYCut = 0.005;
  double TightISO = 0.07;

  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.MET();

  //============================================
  //==== Apply the gen weight (for NLO, +1,-1)
  //============================================

  if(!isData) weight*=MCweight;
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  //===================
  //==== [CUT] No Cut
  //===================

  FillCutFlow("NoCut", 1.);

  //======================
  //==== [CUT] METFilter
  //======================
 
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);

  //=======================
  //==== [CUT] Vertex cut
  //=======================

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  //==== Has Good Primary vertex:
  //==== if ( vtx.ndof() > 4 &&
  //====   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //====   ( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //====   !(vtx.isFake() ) )
  FillCutFlow("VertexCut", 1.);

  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_NOLEPTONVETO", 30., 5.);

  int n_jets = jetColl_hn.size();
  int n_bjets=0;

  int For_HLT_Mu3_PFJet40_v=0;
/*
  //==== HLT_Mu3_PFJet40_v PFJet Pt Check
  if( PassTrigger("HLT_Mu3_PFJet40_v") ){
    if(jetColl_hn.size()!=0){
      FillHist("HLT_Mu3_PFJet40_v_LeadingJetPt", jetColl_hn.at(0).Pt(), 1., 0., 200., 200);
    }
  }
  return;
*/
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jetColl_hn.at(j).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
  }

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight=(1.0);
  if(!k_isdata){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  //================================================
  //==== Apply corrections (PU reweight only here)
  //================================================

  if(!isData && !k_running_nonprompt){
    weight*=pileup_reweight;
    weight*=GetKFactor();
  }

  //==========================
  //==== Call Loosest Muon
  //==========================

  std::vector<snu::KMuon> muontriNodXYCutVLooseColl_raw = GetMuons("MUON_HN_NODXYCUT_VLOOSE_lowestPtCut", true);

  //==================
  //==== FR binnings
  //==================

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {0., 5., 10., 15., 20., 25., 30., 35., 45., 60., 80., 100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  //===========
  //==== HLTs
  //===========

  std::map< TString, std::vector<double> > HLT_ptrange;
  HLT_ptrange["HLT_Mu3_PFJet40_v"].push_back(5.);
  HLT_ptrange["HLT_Mu3_PFJet40_v"].push_back(10.);
  HLT_ptrange["HLT_Mu8_v"].push_back(10.);
  HLT_ptrange["HLT_Mu8_v"].push_back(20.);
  HLT_ptrange["HLT_Mu8_TrkIsoVVL_v"].push_back(10.);
  HLT_ptrange["HLT_Mu8_TrkIsoVVL_v"].push_back(20.);
  HLT_ptrange["HLT_Mu17_v"].push_back(20.);
  HLT_ptrange["HLT_Mu17_v"].push_back(25.);
  HLT_ptrange["HLT_Mu17_TrkIsoVVL_v"].push_back(20.);
  HLT_ptrange["HLT_Mu17_TrkIsoVVL_v"].push_back(25.);
  HLT_ptrange["HLT_Mu20_v"].push_back(25.);
  HLT_ptrange["HLT_Mu20_v"].push_back(9999.);
  HLT_ptrange["HLT_IsoMu24_v"].push_back(25.);
  HLT_ptrange["HLT_IsoMu24_v"].push_back(9999.);

  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }

  //==== 2) back-to-back dijet topology

  //==== tag jet collections
  //==== pt > 40 GeV
  //==== |eta| > 2.4
  //==== LeptonVeto
  std::vector<snu::KJet> jetColl_tag = GetJets("JET_HN", 40., 2.4);
  std::vector<snu::KMuon> TEST_tight = GetMuons("MUON_HN_TIGHT", false);
  std::vector<snu::KMuon> TEST_loose = GetMuons("MUON_HN_LOOSE", false);

  if( PassTrigger("HLT_Mu17_v") || PassTrigger("HLT_Mu8_v") ){

    //==== dijet event seletion
    if(jetColl_tag.size() != 0){
      if(TEST_loose.size() == 1){

        Double_t this_weight_Loose = weight*GetPrescale(TEST_loose, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        snu::KMuon muon = TEST_loose.at(0);

        float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;
        bool histfilled = false;
        for(unsigned int i=0; i<jetColl_tag.size(); i++){

          if(histfilled) break;
          dR = jetColl_tag.at(i).DeltaR(muon);
          dPhi = fabs(jetColl_tag.at(i).DeltaPhi(muon));
          ptmu_ptjet = muon.Pt() / jetColl_tag.at(i).Pt();
          FillHist("SingleMuonTrigger_Dijet_dR", dR, this_weight_Loose, 0., 10., 100);

          if( dR > 1.0 && ptmu_ptjet < 1. && dPhi > 2.5 ){
            FillDenAndNum("SingleMuonTrigger_Dijet", muon, this_weight_Loose, (TEST_tight.size() == 1) );
            if(TEST_tight.size() == 1){
              histfilled=true;
            }
          }

        } // tag jet for loop

      } // only one muon
    } // tag jet exist
  } // trigger

  //=========================================================
  //==== Large dXYSig Muon definitions for systematic study
  //=========================================================

  double dXYMin_central = 4.0;
  double RelIsoMax_central = 0.4;

  const int n_dXYMins = 3;
  double dXYMins[n_dXYMins] = {3.0, 4.0, 5.0};
  const int n_RelIsoMaxs = 6;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.2, 0.3, 0.4, 0.6, 0.8, 1.0};

/*
  //==== Central values only for test
  const int n_dXYMins = 1;
  double dXYMins[n_dXYMins] = {4.0};
  const int n_RelIsoMaxs = 1;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.4};
*/
  //=================
  //==== Start Loop
  //=================

  for(int aaa=0; aaa<n_dXYMins; aaa++){
    for(int bbb=0; bbb<n_RelIsoMaxs; bbb++){

      int dXY_Digit1 = int(dXYMins[aaa]);
      int dXY_Digit0p1 = 10*dXYMins[aaa]-10*dXY_Digit1;
      TString str_dXYCut = "dXYSigMin_"+TString::Itoa(dXY_Digit1,10)+"p"+TString::Itoa(dXY_Digit0p1,10);

      int iso_Digit1 = int(RelIsoMaxs[bbb]);
      int iso_Digit0p1 = 10*RelIsoMaxs[bbb]-10*iso_Digit1;
      TString str_iso = "LooseRelIsoMax_"+TString::Itoa(iso_Digit1,10)+"p"+TString::Itoa(iso_Digit0p1,10);

      str_dXYCut = str_dXYCut+"_"+str_iso;

      //===============================
      //==== Prepare Muon Collections
      //===============================

      //==== NodXYCut muons
      std::vector<snu::KMuon> muontriNodXYCutLooseColl_raw;
      std::vector<snu::KMuon> muontriNodXYCutTightColl_raw;
      std::vector<snu::KMuon> muontriNodXYCutLooseColl;
      std::vector<snu::KMuon> muontriNodXYCutTightColl;
      //==== Small dXYSig muons
      std::vector<snu::KMuon> muontriLooseColl_raw;
      std::vector<snu::KMuon> muontriTightColl_raw;
      std::vector<snu::KMuon> muontriLooseColl;
      std::vector<snu::KMuon> muontriTightColl;
      //==== Large dXYSig muons
      std::vector<snu::KMuon> muontriHighdXYLooseColl_raw;
      std::vector<snu::KMuon> muontriHighdXYTightColl_raw;
      std::vector<snu::KMuon> muontriHighdXYLooseColl;
      std::vector<snu::KMuon> muontriHighdXYTightColl;

      for(unsigned int j=0; j<muontriNodXYCutVLooseColl_raw.size(); j++){
        snu::KMuon this_muon = muontriNodXYCutVLooseColl_raw.at(j);
        //==== MaxRelIso for Loose
        if( this_muon.RelIso04() >= RelIsoMaxs[bbb] ) continue;

        //==== Fill NodXYCut
        muontriNodXYCutLooseColl_raw.push_back( this_muon ); 
        if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl_raw.push_back( this_muon );
        //==== Fill Small dXYSig
        if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
          muontriLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriTightColl_raw.push_back( this_muon );
        }
        //==== Fill Large dXYSig
        if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
          muontriHighdXYLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl_raw.push_back( this_muon );
        }

        //==== for data, or QCD
        //==== do the same for prompt collection
        if( k_isdata || k_sample_name.Contains("QCD") ){
          //==== Fill NodXYCut
          muontriNodXYCutLooseColl.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl.push_back( this_muon );
          //==== Fill Small dXYSig
          if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
            muontriLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriTightColl.push_back( this_muon );
          }
          //==== Fill Large dXYSig
          if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
            muontriHighdXYLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl.push_back( this_muon );
          }
        }

        //==== for MC,
        //==== fill only MCMatched() muons
        else{
          if( TruthMatched(this_muon) ){
            //==== Fill NodXYCut
            muontriNodXYCutLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl.push_back( this_muon );
            //==== Fill Small dXYSig
            if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
              muontriLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < TightISO ) muontriTightColl.push_back( this_muon );
            }
            //==== Fill Large dXYSig
            if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              muontriHighdXYLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl.push_back( this_muon );
            }
          }
        }


      } //==== Loosest muon loop

      //===================
      //==== dXYCut study
      //===================

      if( k_sample_name.Contains("DY") || k_sample_name.Contains("HN") ){

        for(unsigned int i=0; i<muontriNodXYCutLooseColl.size(); i++){
          snu::KMuon this_muon = muontriNodXYCutLooseColl.at(i);
          if(fabs(this_muon.dXY()) < 1.){
            FillHist(str_dXYCut+"_prompt_Loose_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_prompt_Loose_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            if(this_muon.RelIso04()<TightISO){
              FillHist(str_dXYCut+"_prompt_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_prompt_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }
      if( k_sample_name.Contains("TTJets_aMC") || k_sample_name.Contains("QCD") ){

        for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
          snu::KMuon this_muon = muontriNodXYCutLooseColl_raw.at(i);
          if( fabs(this_muon.dXY()) < 1. && !TruthMatched(this_muon) ){
            FillHist(str_dXYCut+"_fake_Loose_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_fake_Loose_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            if(this_muon.RelIso04()<TightISO){
              FillHist(str_dXYCut+"_fake_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_fake_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }

      //=========================
      //==== SingleMuon Trigger
      //=========================

      if( PassTriggerOR(AllHLTs) ){

        std::map<TString, double> this_weight_Loose, this_weight_HighdXYLoose, this_weight_NodXYCutLoose;
        for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
          this_weight_Loose[it->first]         = weight*GetTriggerWeightByPtRange(it->first, it->second, muontriLooseColl, For_HLT_Mu3_PFJet40_v);
          this_weight_HighdXYLoose[it->first]  = weight*GetTriggerWeightByPtRange(it->first, it->second, muontriHighdXYLooseColl, For_HLT_Mu3_PFJet40_v);
          this_weight_NodXYCutLoose[it->first] = weight*GetTriggerWeightByPtRange(it->first, it->second, muontriNodXYCutLooseColl, For_HLT_Mu3_PFJet40_v);
        }

        //Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_NodXYCutLoose = weight*GetPrescale(muontriNodXYCutLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        //==== 1) LooseMuon study

        if( muontriLooseColl.size() == 1 ){
          snu::KMuon muon = muontriLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.1, 100);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_Loose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dZ", fabs(muon.dZ()), this_weight_Loose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_onebin", 0., this_weight_Loose, 0., 1., 1);
        }

        //==== 2) SingleMuon HighdXY LooseMuon study

        if( muontriHighdXYLooseColl.size() == 1 ){
          snu::KMuon muon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_onebin", 0., this_weight_HighdXYLoose, 0., 1., 1);
        }

        //==== 2') SingleMuon NodXY Loose Muon study
        if( muontriNodXYCutLooseColl.size() == 1 ){
          snu::KMuon muon = muontriNodXYCutLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_RelIso", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_NodXYCutLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_NodXYCutLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_NodXYCutLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_onebin", 0., this_weight_NodXYCutLoose, 0., 1., 1);
        }


        //==== 3) large dXY muon method

        if(muontriHighdXYLooseColl.size()==1){
          snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = HighdXYmuon.RelIso04();

          //==== half sample test
          if(dXYMins[aaa]== dXYMin_central && RelIsoMaxs[bbb]==RelIsoMax_central){

            int EventNumber = eventbase->GetEvent().EventNumber();

            double pt_for_FR = HighdXYmuon.Pt();
            double eta_for_FR = fabs(HighdXYmuon.Eta());
            if(pt_for_FR>=60.) pt_for_FR=59.;
            if(eta_for_FR>=2.5) eta_for_FR=2.2;

            snu::KEvent Evt = eventbase->GetEvent();
            double MET = Evt.MET();

            //==== sample A
            //==== obtain FR_A(pt, eta)
            if(abs(EventNumber%2==0)){
              FillHist("TEST_HalfSample", 0., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleA";
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0", 
                                HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              }

            }
            //==== sample B
            //==== 1) obtain FR_measured(MET)
            //==== 2-1) obtain MET distribution := MET_F0
            //==== 2-2) obtain MET weighted by FR_A(pt_B, eta_B) := MET_F_predicted
            //==== 2-3) MET_F_predicted / MET_F0 := FR_predicted(MET)
            else{
              FillHist("TEST_HalfSample", 1., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleB";
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0",
                                HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F0", MET, this_weight_HighdXYLoose, 0., 500., 500);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F", MET, this_weight_HighdXYLoose, 0., 500., 500);
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              }

              int this_FR_bin = FR_sampleA->FindBin(pt_for_FR,eta_for_FR);
              double this_FR_sampleA = FR_sampleA->GetBinContent(this_FR_bin);

              std::map<TString, double> this_weight_HighdXYLoose_FRweighted;
              for(std::map< TString, double >::iterator it=this_weight_HighdXYLoose.begin(); it!=this_weight_HighdXYLoose.end(); it++){
                this_weight_HighdXYLoose_FRweighted[it->first] = (it->second)*this_FR_sampleA;
              }

              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_Predicted", MET, this_weight_HighdXYLoose_FRweighted, 0., 500., 500);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_Predicted", n_jets, this_weight_HighdXYLoose_FRweighted, 0., 10., 10);

            }

          }

          bool isThisTight = (LeptonRelIso < TightISO);

          //==== all jet
          FillDenAndNum(str_dXYCut+"_HighdXY_alljet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);

          //==== no jet
          if( n_jets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0jet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }

          //==== no b-jet
          if( n_bjets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0bjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with b-jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withbjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }


        } // END one HighdXYSig Muon 

        //==== 5) MC Truth
        //==== here, we pick FAKE muons using GEN info.
        if( !k_isdata ){

          for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){

            snu::KMuon muon = muontriNodXYCutLooseColl_raw.at(i);

            //==== if prompt, skip
            if( TruthMatched(muon) ) continue;

            double LeptonRelIso = muon.RelIso04();
            bool isThisTight = (LeptonRelIso < TightISO);

            //=================================================
            //==== 1) |dXY| < 1 cm, |dXY/err| > dXYMins[aaa].
            //=================================================

            if( fabs( muon.dXY() ) < 1. &&
                fabs( muon.dXYSig() ) > dXYMins[aaa] ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_alljet", muon, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0jet", muon, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withjet", muon, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0bjet", muon, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withbjet", muon, 1., isThisTight);
              }

            } // small dXYSig region


            //=========================================
            //==== 2) |dXY| < thisdXYCut cm, |dXY/err| < 3.
            //=========================================

            if( fabs( muon.dXY() ) < thisdXYCut &&
                fabs( muon.dXYSig() ) < 3. ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_alljet", muon, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0jet", muon, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withjet", muon, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0bjet", muon, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withbjet", muon, 1., isThisTight);
              }

            } // Large dXYSig Region


          } // END Muon loop
        } // END k_isdata


      } // SingleMuon trigger fired


      //==================
      //==== TagZ method
      //==================

      std::vector<TString> triggerlist;
      triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

      if( PassTriggerOR(triggerlist) ){
        float trigger_ps_weight = WeightByTrigger(triggerlist, TargetLumi);
        double this_weight = weight*trigger_ps_weight;

        int n_nodxycut_loose = muontriNodXYCutLooseColl.size();
        if( n_nodxycut_loose >= 3){

          //==== Find Z-pair
          double m_OS_Zclosest = (91.2+15.0);
          int z_pair_index[2];
          bool isZTagged = false;

          for(unsigned i=0; i<n_nodxycut_loose-1; i++){
            for(unsigned j=i+1; j<n_nodxycut_loose; j++){
              if( muontriNodXYCutLooseColl.at(i).Charge() == muontriNodXYCutLooseColl.at(j).Charge() ) continue;
          bool i_muon_tight =  fabs(muontriNodXYCutLooseColl.at(i).dXY())<thisdXYCut && fabs(muontriNodXYCutLooseColl.at(i).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(i).RelIso04()<TightISO;
          bool j_muon_tight =  fabs(muontriNodXYCutLooseColl.at(j).dXY())<thisdXYCut && fabs(muontriNodXYCutLooseColl.at(j).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(j).RelIso04()<TightISO;

              double m_thisOS = ( muontriNodXYCutLooseColl.at(i)+muontriNodXYCutLooseColl.at(j) ).M();
              if(  i_muon_tight &&  j_muon_tight && fabs(m_thisOS-91.2) < fabs(m_OS_Zclosest-91.2) ){
                isZTagged = true;
                m_OS_Zclosest = m_thisOS;
                z_pair_index[0] = i;
                z_pair_index[1] = j;
              }
            }
          }
          if(!isZTagged) continue;
          for(unsigned i=0; i<n_nodxycut_loose; i++){
            if(i==z_pair_index[0] || i==z_pair_index[1]) continue;

            snu::KMuon this_muon =  muontriNodXYCutLooseColl.at(i);
            double conept = MuonConePt(this_muon,TightISO);

            //==== Large dXYSig
            if( fabs(this_muon.dXY()) < 1.0 && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_vs_eta_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_cone_vs_eta_F0", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);


              if( this_muon.RelIso04() < TightISO ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_vs_eta_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_cone_vs_eta_F", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);

              }
            }

            //==== Small dXYSig
            if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_vs_eta_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_cone_vs_eta_F0", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              if( this_muon.RelIso04() < TightISO ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_vs_eta_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_cone_vs_eta_F", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              }
            }

          }

        } //== END at least three muons


      } //==== DiMuon Trigger



    } //==== RelIso loop
  } //==== dXYMin loop


   return;
}// End of execute event loop
  


void FakeRateCalculator_Mu_dxysig_DILEP::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_Mu_dxysig_DILEP::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

FakeRateCalculator_Mu_dxysig_DILEP::~FakeRateCalculator_Mu_dxysig_DILEP() {
  
  Message("In FakeRateCalculator_Mu_dxysig_DILEP Destructor" , INFO);
  
}


void FakeRateCalculator_Mu_dxysig_DILEP::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 3,0.,3.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");

  }
}


void FakeRateCalculator_Mu_dxysig_DILEP::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_Mu_dxysig_DILEP::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_Mu_dxysig_DILEPCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_Mu_dxysig_DILEP::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

float FakeRateCalculator_Mu_dxysig_DILEP::GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh){
  float prescale_trigger = 0.;

  if(muon.size() == 1){

    if(muon.at(0).Pt() >= 20.){
      if(passhigh){
        prescale_trigger = WeightByTrigger("HLT_Mu17_v", TargetLumi) ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = WeightByTrigger("HLT_Mu8_v", TargetLumi) ;
      }
      else prescale_trigger = 0.;
    }
  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}

double FakeRateCalculator_Mu_dxysig_DILEP::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, std::vector<snu::KMuon> muons, int npfjet50){

  double prescale_trigger = 0.;

  if(muons.size()==1){
    snu::KMuon muon = muons.at(0);
    double minpt = ptrange.at(0);
    double maxpt = ptrange.at(1);

    if(PassTrigger(hltname)){
      if(muon.Pt() >= minpt && muon.Pt() < maxpt){
        prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
      }

      if(hltname=="HLT_Mu3_PFJet40_v"){
        if(npfjet50==0) prescale_trigger = 0.;
      }
    }

  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}

void FakeRateCalculator_Mu_dxysig_DILEP::FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value, it->second, xmin, xmax, nbins);
  }

}
void FakeRateCalculator_Mu_dxysig_DILEP::FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value1, value2, it->second, x, nbinsx, y, nbinsy);
  }

}

void FakeRateCalculator_Mu_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {0., 5., 10., 15., 20., 25., 30., 35., 45., 60., 80., 100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  double TightISO = 0.07;
  double conept = MuonConePt(muon,TightISO);

  FillHist(prefix+"_eta_F0", muon.Eta(), thisweight, -3., 3., 30);
  FillHist(prefix+"_pt_F0", muon.Pt(), thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_F0", conept, thisweight, 0., 200., 200);
  FillHist(prefix+"_RelIso_F0", muon.RelIso04(), thisweight, 0., 1., 100);
  FillHist(prefix+"_Chi2_F0", muon.GlobalChi2(), thisweight, 0., 50., 50);
  FillHist(prefix+"_dXY_F0", fabs(muon.dXY()), thisweight, 0., 1., 1000);
  FillHist(prefix+"_dXYSig_F0", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
  FillHist(prefix+"_dZ_F0", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
  FillHist(prefix+"_onebin_F0", 0., thisweight, 0., 1., 1);
  FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_PFMET_F0", METauto, thisweight, 0., 1000., 1000);

  if( isTight ){
    FillHist(prefix+"_eta_F", muon.Eta(), thisweight, -3., 3., 30);
    FillHist(prefix+"_pt_F", muon.Pt(), thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_F", conept, thisweight, 0., 200., 200);
    FillHist(prefix+"_RelIso_F", muon.RelIso04(), thisweight, 0., 1., 100);
    FillHist(prefix+"_Chi2_F", muon.GlobalChi2(), thisweight, 0., 50., 50);
    FillHist(prefix+"_dXY_F", fabs(muon.dXY()), thisweight, 0., 1., 1000);
    FillHist(prefix+"_dXYSig_F", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
    FillHist(prefix+"_dZ_F", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
    FillHist(prefix+"_onebin_F", 0., thisweight, 0., 1., 1);
    FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);

  }


}

void FakeRateCalculator_Mu_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KMuon muon, std::map<TString, double> hltweight, bool isTight){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){

    FillDenAndNum(it->first+"_"+prefix, muon, it->second, isTight);

  }

}


















