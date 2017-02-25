// $Id: FakeRateCalculator_Mu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_Mu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_Mu.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_Mu);

FakeRateCalculator_Mu::FakeRateCalculator_Mu() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_Mu");
  
  Message("In FakeRateCalculator_Mu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_Mu::InitialiseAnalysis() throw( LQError ) {
  
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


  return;
}


void FakeRateCalculator_Mu::ExecuteEvents()throw( LQError ){

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

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);

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
  }

  //==========================
  //==== Call Loosest Muon
  //==========================

  std::vector<snu::KMuon> muontriNodXYCutVLooseColl_raw = GetMuons("MUON_HN_TRI_NODXYCUT_VLOOSE", true);

  //==================
  //==== FR binnings
  //==================

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  //==== 2) back-to-back dijet topology

  //==== tag jet collections
  //==== pt > 40 GeV
  //==== |eta| > 2.4
  //==== LeptonVeto
  std::vector<snu::KJet> jetColl_tag = GetJets("JET_HN", 40., 2.4);
  std::vector<snu::KMuon> TEST_tight = GetMuons("MUON_HN_TRI_TIGHT");
  std::vector<snu::KMuon> TEST_loose = GetMuons("MUON_HN_TRI_LOOSE");

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
            FillHist("SingleMuonTrigger_Dijet_eta_F0", muon.Eta(), this_weight_Loose, -3, 3, 30);
            FillHist("SingleMuonTrigger_Dijet_pt_F0", muon.Pt(), this_weight_Loose, 0., 200., 200);
            FillHist("SingleMuonTrigger_Dijet_RelIso_F0", muon.RelIso04(), this_weight_Loose, 0., 1., 100);
            FillHist("SingleMuonTrigger_Dijet_Chi2_F0", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
            FillHist("SingleMuonTrigger_Dijet_dXY_F0", fabs(muon.dXY()), this_weight_Loose, 0., 1., 1000);
            FillHist("SingleMuonTrigger_Dijet_dXYSig_F0", fabs(muon.dXYSig()), this_weight_Loose, 0., 15., 150);
            FillHist("SingleMuonTrigger_Dijet_dZ_F0", fabs(muon.dZ()), this_weight_Loose, 0., 0.5, 50);
            FillHist("SingleMuonTrigger_Dijet_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
            if(TEST_tight.size() == 1){
              histfilled=true;
              FillHist("SingleMuonTrigger_Dijet_eta_F", muon.Eta(), this_weight_Loose, -3, 3, 30);
              FillHist("SingleMuonTrigger_Dijet_pt_F", muon.Pt(), this_weight_Loose, 0., 200., 200);
              FillHist("SingleMuonTrigger_Dijet_RelIso_F", muon.RelIso04(), this_weight_Loose, 0., 1., 100);
              FillHist("SingleMuonTrigger_Dijet_Chi2_F", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
              FillHist("SingleMuonTrigger_Dijet_dXY_F", fabs(muon.dXY()), this_weight_Loose, 0., 1, 1000);
              FillHist("SingleMuonTrigger_Dijet_dXYSig_F", fabs(muon.dXYSig()), this_weight_Loose, 0., 15., 150);
              FillHist("SingleMuonTrigger_Dijet_dZ_F", fabs(muon.dZ()), this_weight_Loose, 0., 0.5, 50);
              FillHist("SingleMuonTrigger_Dijet_events_F", muon.Pt(), fabs(muon.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
            }
          }

        } // tag jet for loop

      } // only one muon
    } // tag jet exist
  } // trigger

  //=========================================================
  //==== Large dXYSig Muon definitions for systematic study
  //=========================================================

  const int n_dXYMins = 3;
  double dXYMins[n_dXYMins] = {3.0, 4.0, 5.0};
  const int n_RelIsoMaxs = 6;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.2, 0.3, 0.4, 0.6, 0.8, 1.0};

  //==== Central values only for test
  //const int n_dXYMins = 1;
  //double dXYMins[n_dXYMins] = {4.0};
  //const int n_RelIsoMaxs = 1;
  //double RelIsoMaxs[n_RelIsoMaxs] = {0.4};

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
        if( this_muon.RelIso04() < 0.1 ) muontriNodXYCutTightColl_raw.push_back( this_muon );
        //==== Fill Small dXYSig
        if( fabs(this_muon.dXY()) < 0.05 && fabs(this_muon.dXYSig()) < 3.0 ){
          muontriLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < 0.1 ) muontriTightColl_raw.push_back( this_muon );
        }
        //==== Fill Large dXYSig
        if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
          muontriHighdXYLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < 0.1 ) muontriHighdXYTightColl_raw.push_back( this_muon );
        }

        //==== for data, or QCD
        //==== do the same for prompt collection
        if( k_isdata || k_sample_name.Contains("QCD") ){
          //==== Fill NodXYCut
          muontriNodXYCutLooseColl.push_back( this_muon );
          if( this_muon.RelIso04() < 0.1 ) muontriNodXYCutTightColl.push_back( this_muon );
          //==== Fill Small dXYSig
          if( fabs(this_muon.dXY()) < 0.05 && fabs(this_muon.dXYSig()) < 3.0 ){
            muontriLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < 0.1 ) muontriTightColl.push_back( this_muon );
          }
          //==== Fill Large dXYSig
          if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
            muontriHighdXYLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < 0.1 ) muontriHighdXYTightColl.push_back( this_muon );
          }
        }

        //==== for MC,
        //==== fill only MCMatched() muons
        else{
          if( this_muon.MCMatched() ){
            //==== Fill NodXYCut
            muontriNodXYCutLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < 0.1 ) muontriNodXYCutTightColl.push_back( this_muon );
            //==== Fill Small dXYSig
            if( fabs(this_muon.dXY()) < 0.05 && fabs(this_muon.dXYSig()) < 3.0 ){
              muontriLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < 0.1 ) muontriTightColl.push_back( this_muon );
            }
            //==== Fill Large dXYSig
            if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              muontriHighdXYLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < 0.1 ) muontriHighdXYTightColl.push_back( this_muon );
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
            if(this_muon.RelIso04()<0.1){
              FillHist(str_dXYCut+"_prompt_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_prompt_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }
      if( k_sample_name.Contains("TTJets_aMC") || k_sample_name.Contains("QCD") ){

        for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
          snu::KMuon this_muon = muontriNodXYCutLooseColl_raw.at(i);
          if( fabs(this_muon.dXY()) < 1. && !this_muon.MCMatched() ){
            FillHist(str_dXYCut+"_fake_Loose_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_fake_Loose_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            if(this_muon.RelIso04()<0.1){
              FillHist(str_dXYCut+"_fake_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_fake_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }

      //=========================
      //==== SingleMuon Trigger
      //=========================

      bool highpt_trigger(false), lowpt_trigger(false);
      highpt_trigger = PassTrigger("HLT_Mu17_v");
      lowpt_trigger = PassTrigger("HLT_Mu8_v");

      if( highpt_trigger || lowpt_trigger ){

        Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        Double_t this_weight_NodXYCutLoose = weight*GetPrescale(muontriNodXYCutLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        //==== 1) LooseMuon study

        if( muontriLooseColl.size() == 1 ){
          snu::KMuon muon = muontriLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.1, 100);
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_Loose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_dZ", fabs(muon.dZ()), this_weight_Loose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_LooseMuon_onebin", 0., this_weight_Loose, 0., 1., 1);
        }

        //==== 2) SingleMuon HighdXY LooseMuon study

        if( muontriHighdXYLooseColl.size() == 1 ){
          snu::KMuon muon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_LooseMuon_onebin", 0., this_weight_HighdXYLoose, 0., 1., 1);
        }

        //==== 2') SingleMuon NodXY Loose Muon study
        if( muontriNodXYCutLooseColl.size() == 1 ){
          snu::KMuon muon = muontriNodXYCutLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_RelIso", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_NodXYCutLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_NodXYCutLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_NodXYCutLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_NodXY_LooseMuon_onebin", 0., this_weight_NodXYCutLoose, 0., 1., 1);
        }


        //==== 3) large dXY muon method

        if(muontriHighdXYLooseColl.size()==1){
          snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = HighdXYmuon.RelIso04();
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dXYSig_F0", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
          FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dXYSig_F", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
          }

          //==== check # of jet dependecy
          //==== no jet

          int n_jets = jetColl_hn.size();
          int n_bjets=0;
          for(int j=0; j<n_jets; j++){
            if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)) n_bjets++;
          }

          if( n_jets == 0){
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dXYSig_F0", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < 0.1 ){
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dXYSig_F", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_0jet_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
          }
          //==== with jet
          if( n_jets != 0){
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dXYSig_F0", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < 0.1 ){
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dXYSig_F", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
            //==== with b-jet
            if( n_bjets == 0){
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dXYSig_F0", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dXYSig_F", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_0bjet_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              }
            }
            //==== no b-jet
            if( n_bjets != 0){
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dXYSig_F0", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0., 50., 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dXYSig_F", fabs(HighdXYmuon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
                FillHist(str_dXYCut+"_SingleMuonTrigger_HighdXY_withjet_withbjet_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              }
            }

          } // END with jet

        }

        //==== 5) MC Truth
        //==== here, we pick FAKE muons using GEN info.
        if( !k_isdata ){

          for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){

            snu::KMuon muon = muontriNodXYCutLooseColl_raw.at(i);

            //==== if prompt, skip
            if( muon.MCMatched() ) continue;

            double LeptonRelIso = muon.RelIso04();

            //==== 1) |dXY| < 1 cm, |dXY/err| > dXYMins[aaa].
            if( fabs( muon.dXY() ) < 1. &&
                fabs( muon.dXYSig() ) > dXYMins[aaa] ){

              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_eta_F0", muon.Eta(), 1., -3, 3, 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_pt_F0", muon.Pt(), 1., 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_Chi2_F0", muon.GlobalChi2(), 1., 0, 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dXY_F0", fabs(muon.dXY()),1., 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dXYSig_F0", fabs(muon.dXYSig()),1., 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dZ_F0", fabs(muon.dZ()), 1., 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_onebin_F0", 0., 1., 0., 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), 1., ptarray, 9, etaarray, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_eta_F", muon.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_pt_F", muon.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_Chi2_F", muon.GlobalChi2(), 1., 0, 50., 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dXY_F", fabs(muon.dXY()), 1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dXYSig_F", fabs(muon.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_dZ_F", fabs(muon.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_onebin_F", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), 1., ptarray, 9, etaarray, 4);
              }

            }

            //==== 2) |dXY| < 0.05 cm, |dXY/err| < 3.
            if( fabs( muon.dXY() ) < 0.05 &&
                fabs( muon.dXYSig() ) < 3. ){

              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_eta_F0", muon.Eta(), 1., -3, 3, 30);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_pt_F0", muon.Pt(), 1., 0., 200., 200);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_Chi2_F0", muon.GlobalChi2(), 1., 0, 50., 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dXY_F0", fabs(muon.dXY()), 1., 0., 0.1, 100);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dXYSig_F0", fabs(muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dZ_F0", fabs(muon.dZ()), 1., 0, 0.5, 50);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_onebin_F0", 0., 1., 0, 1., 1);
              FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_events_F0", muon.Pt(), fabs(muon.Eta()), 1., ptarray, 9, etaarray, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_eta_F", muon.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_pt_F", muon.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_Chi2_F", muon.GlobalChi2(), 1., 0, 50., 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dXY_F", fabs(muon.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dXYSig_F", fabs(muon.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_dZ_F", fabs(muon.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_onebin_F", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleMuonTrigger_MCTruth_events_F", muon.Pt(), fabs(muon.Eta()), 1., ptarray, 9, etaarray, 4);
              }

            }


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
          bool i_muon_tight =  fabs(muontriNodXYCutLooseColl.at(i).dXY())<0.05 && fabs(muontriNodXYCutLooseColl.at(i).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(i).RelIso04()<0.1;
          bool j_muon_tight =  fabs(muontriNodXYCutLooseColl.at(j).dXY())<0.05 && fabs(muontriNodXYCutLooseColl.at(j).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(j).RelIso04()<0.1;

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

            //==== Large dXYSig
            if( fabs(this_muon.dXY()) < 1.0 && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);


              if( this_muon.RelIso04() < 0.1 ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);

              }
            }

            //==== Small dXYSig
            if( fabs(this_muon.dXY()) < 0.05 && fabs(this_muon.dXYSig()) < 3.0 ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              if( this_muon.RelIso04() < 0.1 ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              }
            }

          }

        } //== END at least three muons


      } //==== DiMuon Trigger



    } //==== RelIso loop
  } //==== dXYMin loop


  if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ){
    float trigger_ps_weight = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);
    double this_weight = weight*trigger_ps_weight;
    std::vector<snu::KMuon> PtEtaMuons = GetMuons("MUON_PTETA", false);
    for(unsigned int i=0; i<PtEtaMuons.size(); i++){
      snu::KMuon this_muon = PtEtaMuons.at(i);
      if(this_muon.MCMatched()){
        FillHist("DiMuonTrigger_PromptRate_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
        FillHist("DiMuonTrigger_PromptRate_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
        FillHist("DiMuonTrigger_PromptRate_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
        FillHist("DiMuonTrigger_PromptRate_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
        FillHist("DiMuonTrigger_PromptRate_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
        FillHist("DiMuonTrigger_PromptRate_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
        FillHist("DiMuonTrigger_PromptRate_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
        FillHist("DiMuonTrigger_PromptRate_onebin_F0", 0., this_weight, 0, 1., 1);
        FillHist("DiMuonTrigger_PromptRate_events_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
        if( this_muon.RelIso04() < 0.1 ){
          FillHist("DiMuonTrigger_PromptRate_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
          FillHist("DiMuonTrigger_PromptRate_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
          FillHist("DiMuonTrigger_PromptRate_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
          FillHist("DiMuonTrigger_PromptRate_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
          FillHist("DiMuonTrigger_PromptRate_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
          FillHist("DiMuonTrigger_PromptRate_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
          FillHist("DiMuonTrigger_PromptRate_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
          FillHist("DiMuonTrigger_PromptRate_onebin_F", 0., this_weight, 0, 1., 1);
          FillHist("DiMuonTrigger_PromptRate_events_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
        }
      }
    }
  }


   return;
}// End of execute event loop
  


void FakeRateCalculator_Mu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_Mu::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_Mu::~FakeRateCalculator_Mu() {
  
  Message("In FakeRateCalculator_Mu Destructor" , INFO);
  
}


void FakeRateCalculator_Mu::FillCutFlow(TString cut, float weight){

  
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


void FakeRateCalculator_Mu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_Mu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_MuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_Mu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

float FakeRateCalculator_Mu::GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh){
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

