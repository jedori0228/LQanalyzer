// $Id: FakeRateCalculator_Mu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_Mu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond     <jalmond@cern.ch>        - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_Mu.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                           
#include "BaseSelection.h"

ClassImp (FakeRateCalculator_Mu);
FakeRateCalculator_Mu::FakeRateCalculator_Mu() :  AnalyzerCore(), out_muons(0)  {
  SetLogName("FakeRateCalculator_Mu");
  InitialiseAnalysis();
}


void FakeRateCalculator_Mu::InitialiseAnalysis() throw( LQError ) {
  MakeHistograms();  

  return;
 }


void FakeRateCalculator_Mu::ExecuteEvents()throw( LQError ){

  //==== To list all available triggers 
  //ListTriggersAvailable();
  //return;

  FillCutFlow("NoCut", weight);

  //==== FIXME Do we need this?
  if(!PassBasicEventCuts()) return;
  FillCutFlow("EventCut", weight);
 
  //==== DoubleMuon DataSet, SingleMuon trigger  
  std::vector<TString> triggerlist_Mu8, triggerlist_Mu17;
  triggerlist_Mu8.push_back("HLT_Mu8_v");
  triggerlist_Mu17.push_back("HLT_Mu17_v");
  //==== DoubleMuon DataSet, DoubleMuon trigger
  std::vector<TString> triggerlist_Mu17TkMu8;
  triggerlist_Mu17TkMu8.push_back("HLT_Mu17_TkMu8_v");
  //==== SingleMuon DataSet, SingleMuon trigger
  std::vector<TString> triggerlist_Mu5, triggerlist_Mu12, triggerlist_Mu24;
  triggerlist_Mu5.push_back("HLT_Mu5_v");
  triggerlist_Mu12.push_back("HLT_Mu12_v");
  triggerlist_Mu24.push_back("HLT_Mu24_v");

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; 
  FillCutFlow("VertexCut", weight);
  
  numberVertices = eventbase->GetEvent().nVertices();  
  
  if (MC_pu&&!k_isdata) {
    weight = weight* reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
  }
  

  //=====================
  //==== Select objetcs
  //=====================

  //std::vector<snu::KMuon> muontriTightColl;
  //eventbase->GetMuonSel()->HNtriTightMuonSelection(muontriTightColl);
  //std::vector<snu::KMuon> muontriLooseColl;
  //eventbase->GetMuonSel()->HNtriLooseMuonSelection(muontriLooseColl);
  //std::vector<snu::KMuon> muontriHighdXYTightColl;
  //eventbase->GetMuonSel()->HNtriHighdXYTightMuonSelection(muontriHighdXYTightColl);
  //std::vector<snu::KMuon> muontriHighdXYLooseColl;
  //eventbase->GetMuonSel()->HNtriHighdXYLooseMuonSelection(muontriHighdXYLooseColl);

  //==== Muons
  std::vector<snu::KMuon> muontriTightColl_raw;
  eventbase->GetMuonSel()->HNtriTightMuonSelection(muontriTightColl_raw);
  std::vector<snu::KMuon> muontriLooseColl_raw;
  eventbase->GetMuonSel()->HNtriLooseMuonSelection(muontriLooseColl_raw);
  //==== QCD muons
  std::vector<snu::KMuon> muontriHighdXYTightColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYTightMuonSelection(muontriHighdXYTightColl_raw);
  std::vector<snu::KMuon> muontriHighdXYLooseColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYLooseMuonSelection(muontriHighdXYLooseColl_raw);
  //==== No dXY cut muons
  std::vector<snu::KMuon> muontriNodXYCutTightColl_raw;
  eventbase->GetMuonSel()->HNtriNodXYCutTightMuonSelection(muontriNodXYCutTightColl_raw);
  std::vector<snu::KMuon> muontriNodXYCutLooseColl_raw;
  eventbase->GetMuonSel()->HNtriNodXYCutLooseMuonSelection(muontriNodXYCutLooseColl_raw);

  //==== For MC, we only keep the Prompt Mouns,
  //==== and then subtract the histrograms from the data,
  //==== to get the fake-only values

  //std::vector<snu::KMuon> muontriTightColl = GetTruePrompt(muontriTightColl_raw, false);
  //std::vector<snu::KMuon> muontriLooseColl = GetTruePrompt(muontriLooseColl_raw, false);
  //std::vector<snu::KMuon> muontriHighdXYTightColl = GetTruePrompt(muontriHighdXYTightColl_raw, false);
  //std::vector<snu::KMuon> muontriHighdXYLooseColl = GetTruePrompt(muontriHighdXYLooseColl_raw, false);  

  std::vector<snu::KMuon> muontriTightColl;
  std::vector<snu::KMuon> muontriLooseColl;
  std::vector<snu::KMuon> muontriHighdXYTightColl;
  std::vector<snu::KMuon> muontriHighdXYLooseColl;
  std::vector<snu::KMuon> muontriNodXYCutTightColl;
  std::vector<snu::KMuon> muontriNodXYCutLooseColl;

  for(unsigned int i=0; i<muontriTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriTightColl_raw.at(i);
    if(isData) muontriTightColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriLooseColl_raw.at(i);
    if(isData) muontriLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriLooseColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYTightColl_raw.at(i);
    if(isData) muontriHighdXYTightColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriHighdXYTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYLooseColl_raw.at(i);
    if(isData) muontriHighdXYLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriHighdXYLooseColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
    if(isData) muontriNodXYCutTightColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriNodXYCutTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
    if(isData) muontriNodXYCutLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriNodXYCutLooseColl.push_back(thismuon);
    }
  }


  //==== FIXME can we use this scale factor for our muons?
  //if(!isData){ 
  //   for(std::vector<snu::KMuon>::iterator it = muontriTightColl.begin(); it!= muontriTightColl.end(); it++){
  //     weight *= MuonScaleFactor(it->Eta(), it->Pt(), true);
  //   }
  //}

  std::vector<snu::KJet> jetColl_lepveto;
  eventbase->GetJetSel()->SetEta(5.0);
  //==== should remove the loosest leptons in the analysis
  eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl_lepveto, AnalyzerCore::GetMuons("veto"), AnalyzerCore::GetElectrons(false,false, "veto") );


  if( !PassTrigger(triggerlist_Mu8,prescale) && !PassTrigger(triggerlist_Mu17,prescale) && !PassTrigger(triggerlist_Mu17TkMu8,prescale) ) return;

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};
  //float ptarray [] = {10., 11., 12., 13., 14., 15.};


  //=========================
  //==== SingleMuon Trigger
  //=========================

  if( PassTrigger(triggerlist_Mu8,prescale) || PassTrigger(triggerlist_Mu17,prescale) ){

    Double_t prescale_trigger = GetPrescale(muontriLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));
    Double_t this_weight = weight*prescale_trigger;

    //==== 1) LooseMuon study

    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      snu::KMuon muon = muontriLooseColl.at(i);
      double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
      FillHist("SingleMuon_LooseMuon_RelIso", LeptonRelIso, this_weight, 0., 1., 100);
      FillHist("SingleMuon_LooseMuon_Chi2", muon.GlobalChi2(), this_weight, 0, 10, 10);
    }

    //==== 2) back-to-back dijet topology

    //==== tag jet collections
    std::vector<snu::KJet> jetColl = GetJets("HNtriFRTagJet");
    //==== dijet event seletion
    if(jetColl.size() != 0){
      if(muontriLooseColl.size() == 1){

        //snu::KEvent Evt = eventbase->GetEvent();
        //double MET = Evt.PFMET();
        //if( MET < 40 ) return; // Let Wjets dominate

        snu::KMuon muon;
        muon = muontriLooseColl.at(0);
        muon.SetPxPyPzE(muon.Px(),muon.Py(),muon.Pz(),muon.E());

        float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;

        for(unsigned int i=0; i<jetColl.size(); i++){
          dR = jetColl.at(i).DeltaR(muon);
          dPhi = fabs(jetColl.at(i).DeltaPhi(muon));
          ptmu_ptjet = muon.Pt() / jetColl.at(i).Pt();
          FillHist("Dijet_dR", dR, this_weight, 0., 10., 100);

          if( dR > 1.0 ){
            FillHist("Dijet_ptmu_ptjet", ptmu_ptjet, this_weight, 0., 2., 200);
            FillHist("Dijet_dPhi", dPhi, this_weight, 0., 4., 40);
            if( ptmu_ptjet < 1. ){
              if( dPhi > 2.5 ){
                double HT_loose = AnalyzerCore::SumPt(jetColl_lepveto);
                double HT_tag = AnalyzerCore::SumPt(jetColl);
                double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
                FillHist("Dijet_eta_F0", muon.Eta(), this_weight, -3, 3, 30);
                FillHist("Dijet_pt_F0", muon.Pt(), this_weight, 0., 200., 200);
                FillHist("Dijet_RelIso_F0", LeptonRelIso, this_weight, 0., 1., 100);
                FillHist("Dijet_Chi2_F0", muon.GlobalChi2(), this_weight, 0., 10., 10);
                FillHist("Dijet_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
                FillHist("Dijet_HT_loose_F0", HT_loose, this_weight, 0, 300, 300);
                FillHist("Dijet_HT_tag_F0", HT_tag, this_weight, 0, 300, 300);
                if(muontriTightColl.size() == 1){
                  FillHist("Dijet_eta_F", muon.Eta(), this_weight, -3, 3, 30);
                  FillHist("Dijet_pt_F", muon.Pt(), this_weight, 0., 200., 200);
                  FillHist("Dijet_RelIso_F", LeptonRelIso, this_weight, 0., 1., 100);
                  FillHist("Dijet_Chi2_F", muon.GlobalChi2(), this_weight, 0., 10., 10);
                  FillHist("Dijet_events_F", muon.Pt(), fabs(muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
                  FillHist("Dijet_HT_loose_F", HT_loose, this_weight, 0, 300, 300);
                  FillHist("Dijet_HT_tag_F", HT_tag, this_weight, 0, 300, 300);
                }
                goto stop;
              }
            }
          }
        }
        stop: ;
  
      }
    }

    //==== 3) SingleMuon HighdXY LooseMuon study

    for(unsigned int i=0; i<muontriHighdXYLooseColl.size(); i++){
      snu::KMuon muon = muontriHighdXYLooseColl.at(i);
      double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
      FillHist("SingleMuon_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight, 0., 1., 100);
      FillHist("SingleMuon_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight, 0, 10, 10);
    }

    //==== 4) large dXY muon method
  
    if( muontriHighdXYLooseColl.size() == 1 ){
      snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0); 
      double LeptonRelIso = (HighdXYmuon.SumIsoCHDR03() + std::max(0.0, HighdXYmuon.SumIsoNHDR03() + HighdXYmuon.SumIsoPHDR03() - 0.5* HighdXYmuon.SumPUIsoR03()))/HighdXYmuon.Pt(); 
      FillHist("SingleMuon_HighdXY_eta_F0", HighdXYmuon.Eta(), this_weight, -3, 3, 30);
      FillHist("SingleMuon_HighdXY_pt_F0", HighdXYmuon.Pt(), this_weight, 0., 200., 200);
      FillHist("SingleMuon_HighdXY_RelIso_F0", LeptonRelIso, this_weight, 0., 1., 100);
      FillHist("SingleMuon_HighdXY_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight, 0, 10, 10);
      FillHist("SingleMuon_HighdXY_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight, ptarray, 9, etaarray, 4);
      if( eventbase->GetMuonSel()->HNIstriHighdXYTight(HighdXYmuon) ){
        FillHist("SingleMuon_HighdXY_eta_F", HighdXYmuon.Eta(), this_weight, -3, 3, 30);
        FillHist("SingleMuon_HighdXY_pt_F", HighdXYmuon.Pt(), this_weight, 0., 200., 200);
        FillHist("SingleMuon_HighdXY_RelIso_F", LeptonRelIso, this_weight, 0., 1., 100);
        FillHist("SingleMuon_HighdXY_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight, 0, 10, 10);
        FillHist("SingleMuon_HighdXY_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight, ptarray, 9, etaarray, 4);
      }
    }

  } // SingleMuon trigger fired



  //=====================
  //==== DiMuon Trigger
  //=====================

  if( PassTrigger(triggerlist_Mu17TkMu8,prescale) ){

    //==== 1) dXY cut optimization

    //==== for QCD
    for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutTightColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) continue;
      FillHist("TightIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("TightIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("TightIsoMuon_fake_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1 ){ // 1 cm = 10 mm
        FillHist("TightIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("TightIsoMuon_fake_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) continue;
      FillHist("LooseIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("LooseIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("LooseIsoMuon_fake_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1 ){ // 1 cm = 10 mm
        FillHist("LooseIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("LooseIsoMuon_fake_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    //==== for DY (prompt matched)
    for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutTightColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
      if( thismuon.GetType() != 0 && thismuon.GetType() != 7 ) continue;
      FillHist("TightIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("TightIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("TightIsoMuon_prompt_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1 ){ // 1 cm = 10 mm
        FillHist("TightIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("TightIsoMuon_prompt_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetType() != 0 && thismuon.GetType() != 7 ) continue;
      FillHist("LooseIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("LooseIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("LooseIsoMuon_prompt_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1 ){ // 1 cm = 10 mm
        FillHist("LooseIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("LooseIsoMuon_prompt_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }

    //==== 2) DiMuon HighdXY LooseMuon study

    for(unsigned int i=0; i<muontriHighdXYLooseColl.size(); i++){
      if(muontriHighdXYLooseColl.at(0).Pt() < 20.) break;
      snu::KMuon muon = muontriHighdXYLooseColl.at(i);
      double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
      FillHist("DiMuon_HighdXY_LooseMuon_RelIso", LeptonRelIso, weight, 0., 1., 100);
      FillHist("DiMuon_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), weight, 0, 10, 10);
    }

    //==== 3) Filling num and den for FR

    float etaarray_2 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
    float ptarray_2 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

    int n_jets = jetColl_lepveto.size();

    //==== Two muons
    if( muontriHighdXYLooseColl.size() == 2){
      if( muontriHighdXYLooseColl.at(0).Pt() > 20. ){
        double dimuon_mass = ( muontriHighdXYLooseColl.at(0) + muontriHighdXYLooseColl.at(1) ).M();
        FillHist("DiMuon_HighdXY_dRdimuon", muontriHighdXYLooseColl.at(0).DeltaR( muontriHighdXYLooseColl.at(1) ), weight, 0, 4, 40);
        FillHist("DiMuon_HighdXY_mdimuon", dimuon_mass, weight, 0, 200, 40);
        FillHist("DiMuon_HighdXY_n_jets", n_jets, weight, 0, 10, 10);

        if( (dimuon_mass > 15 && dimuon_mass < 60) || dimuon_mass > 120 ){
          FillHist("DiMuon_HighdXY_mdimuon_after_cut", dimuon_mass, weight, 0, 200, 40);
          for(unsigned int i=0; i<2; i++){

            snu::KMuon muon = muontriHighdXYLooseColl.at(i);
            double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
            FillHist("DiMuon_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuon_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuon_HighdXY_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuon_HighdXY_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("DiMuon_HighdXY_n_jets_F0", n_jets, weight, 0, 10, 10);
            FillHist("DiMuon_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            if( eventbase->GetMuonSel()->HNIstriHighdXYTight(muon, false) ){
              FillHist("DiMuon_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuon_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuon_HighdXY_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuon_HighdXY_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_n_jets_F", n_jets, weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            }

            if( n_jets == 0){
              FillHist("DiMuon_HighdXY_0jet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuon_HighdXY_0jet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuon_HighdXY_0jet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuon_HighdXY_0jet_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_0jet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_0jet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( eventbase->GetMuonSel()->HNIstriHighdXYTight(muon, false) ){
                FillHist("DiMuon_HighdXY_0jet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuon_HighdXY_0jet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuon_HighdXY_0jet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuon_HighdXY_0jet_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
                FillHist("DiMuon_HighdXY_0jet_n_jets_F", n_jets, weight, 0, 10, 10);
                FillHist("DiMuon_HighdXY_0jet_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              }
            }
            if( n_jets != 0){
              FillHist("DiMuon_HighdXY_withjet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuon_HighdXY_withjet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuon_HighdXY_withjet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuon_HighdXY_withjet_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_withjet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuon_HighdXY_withjet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( eventbase->GetMuonSel()->HNIstriHighdXYTight(muon, false) ){
                FillHist("DiMuon_HighdXY_withjet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuon_HighdXY_withjet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuon_HighdXY_withjet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuon_HighdXY_withjet_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
                FillHist("DiMuon_HighdXY_withjet_n_jets_F", n_jets, weight, 0, 10, 10);
                FillHist("DiMuon_HighdXY_withjet_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              }
            }
          }

        } // END 15 < m(mumu) < 60 || m(mumu) > 120


      } // END lead pT > 20 GeV
    } // END if(dimuon)


    //==== 4) MCTruth

    if( !k_isdata ){

      float etaarray_2 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
      float ptarray_2 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

      for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
        //==== lead pT > 20 GeV
        if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

        snu::KMuon muon = muontriNodXYCutLooseColl_raw.at(i);
        //==== if prompt, skip
        if( muon.GetType() == 0 || muon.GetType() == 7 ) continue;
        //==== |dXY| < 1 cm
        if( fabs( muon.dXY() ) > 1. ) continue;

        double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();

        //==== 1) no |dXY|/err cut
        FillHist("MCTruth_HighdXY_no_sigcut_eta_F0", muon.Eta(), weight, -3, 3, 30);
        FillHist("MCTruth_HighdXY_no_sigcut_pt_F0", muon.Pt(), weight, 0., 200., 200);
        FillHist("MCTruth_HighdXY_no_sigcut_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
        FillHist("MCTruth_HighdXY_no_sigcut_pt_F0", muon.GlobalChi2(), weight, 0, 10, 10);
        FillHist("MCTruth_HighdXY_no_sigcut_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
        if( LeptonRelIso < 0.1 ){
          FillHist("MCTruth_HighdXY_no_sigcut_eta_F", muon.Eta(), weight, -3, 3, 30);
          FillHist("MCTruth_HighdXY_no_sigcut_pt_F", muon.Pt(), weight, 0., 200., 200);
          FillHist("MCTruth_HighdXY_no_sigcut_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
          FillHist("MCTruth_HighdXY_no_sigcut_pt_F", muon.GlobalChi2(), weight, 0, 10, 10);
          FillHist("MCTruth_HighdXY_no_sigcut_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
        }

        //==== 2) |dXY|/err > 4
        if( fabs( muon.dXY() / muon.dXYErrPat() ) > 4 ){
          FillHist("MCTruth_HighdXY_large_sig_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("MCTruth_HighdXY_large_sig_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("MCTruth_HighdXY_large_sig_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("MCTruth_HighdXY_large_sig_pt_F0", muon.GlobalChi2(), weight, 0, 10, 10);
          FillHist("MCTruth_HighdXY_large_sig_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("MCTruth_HighdXY_large_sig_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("MCTruth_HighdXY_large_sig_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("MCTruth_HighdXY_large_sig_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("MCTruth_HighdXY_large_sig_pt_F", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("MCTruth_HighdXY_large_sig_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }
        }

        //==== 3) |dXY|/err < 3
        if( fabs( muon.dXY() / muon.dXYErrPat() ) < 3 ){
          FillHist("MCTruth_HighdXY_small_sig_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("MCTruth_HighdXY_small_sig_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("MCTruth_HighdXY_small_sig_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("MCTruth_HighdXY_small_sig_pt_F0", muon.GlobalChi2(), weight, 0, 10, 10);
          FillHist("MCTruth_HighdXY_small_sig_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("MCTruth_HighdXY_small_sig_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("MCTruth_HighdXY_small_sig_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("MCTruth_HighdXY_small_sig_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("MCTruth_HighdXY_small_sig_pt_F", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("MCTruth_HighdXY_small_sig_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }
        }

      }


    } // END !k_isdata


  } // DiMuon trigger fired


  return;
}// End of execute event loop


void FakeRateCalculator_Mu::EndCycle()throw( LQError ){
  Message("In EndCycle" , INFO);
}


void FakeRateCalculator_Mu::BeginCycle() throw( LQError ){
  Message("In begin Cycle", INFO);
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "MyDataPileupHistogram.root").c_str());
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");
  return;
}

FakeRateCalculator_Mu::~FakeRateCalculator_Mu() {
  Message("In FakeRateCalculator_Mu Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
}


void FakeRateCalculator_Mu::FillCutFlow(TString cut, float weight){
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 5,0.,5.);
  }
}


void FakeRateCalculator_Mu::BeginEvent( )throw( LQError ){
  Message("In BeginEvent() " , DEBUG);
  return;
}


void FakeRateCalculator_Mu::MakeHistograms(){
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
}


void FakeRateCalculator_Mu::ClearOutputVectors() throw(LQError) {
  out_muons.clear();
  out_electrons.clear();
}

float FakeRateCalculator_Mu::GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh){
  float prescale_trigger = 0.;
  if(muon.size() == 1){
    if(muon.at(0).Pt() >= 20.){
      if(passhigh){
        prescale_trigger = (25.014) / 19674 ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = (2.021276) / 19674 ;
      }
      else prescale_trigger = 0.;
    }
  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}
