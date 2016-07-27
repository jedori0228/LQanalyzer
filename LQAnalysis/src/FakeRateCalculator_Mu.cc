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
  eventbase->GetJetSel()->SetEta(2.4);
  //==== should remove the loosest leptons in the analysis
  eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl_lepveto, AnalyzerCore::GetMuons("veto"), AnalyzerCore::GetElectrons(false,false, "veto") );


  //if( !PassTrigger(triggerlist_Mu8,prescale) && !PassTrigger(triggerlist_Mu17,prescale) && !PassTrigger(triggerlist_Mu17TkMu8,prescale) ) return;

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  //=========================
  //==== SingleMuon Trigger
  //=========================

  bool highpt_trigger(false), lowpt_trigger(false);
  //==== DoubleMuon Dataset (HLT_Mu8_v, HLT_Mu17_v)
  if( k_jskim_flag_1 == "DoubleMuon" ){
    highpt_trigger = PassTrigger(triggerlist_Mu17,prescale);
    lowpt_trigger = PassTrigger(triggerlist_Mu8,prescale);    
  }
  //==== SingleMuon Dataset (HLT_Mu12_v, HLT_Mu24_v)
  if( k_jskim_flag_1 == "SingleMuon" ){
    highpt_trigger = PassTrigger(triggerlist_Mu24,prescale);
    lowpt_trigger = PassTrigger(triggerlist_Mu12,prescale);
  }

  if( highpt_trigger || lowpt_trigger ){

    //Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, lowpt_trigger, highpt_trigger);
    //Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, lowpt_trigger, highpt_trigger);

  //if( PassTrigger(triggerlist_Mu8,prescale) || PassTrigger(triggerlist_Mu17,prescale) ){

    Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));
    Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));

    //==== 1) LooseMuon study

    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      snu::KMuon muon = muontriLooseColl.at(i);
      double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
      FillHist("SingleMuonTrigger_LooseMuon_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
      FillHist("SingleMuonTrigger_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_Loose, 0, 10, 10);
      FillHist("SingleMuonTrigger_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 20);
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
          FillHist("SingleMuonTrigger_Dijet_dR", dR, this_weight_Loose, 0., 10., 100);

          if( dR > 1.0 ){
            FillHist("SingleMuonTrigger_Dijet_ptmu_ptjet", ptmu_ptjet, this_weight_Loose, 0., 2., 200);
            FillHist("SingleMuonTrigger_Dijet_dPhi", dPhi, this_weight_Loose, 0., 4., 40);
            if( ptmu_ptjet < 1. ){
              if( dPhi > 2.5 ){
                double HT_loose = AnalyzerCore::SumPt(jetColl_lepveto);
                double HT_tag = AnalyzerCore::SumPt(jetColl);
                double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
                FillHist("SingleMuonTrigger_Dijet_eta_F0", muon.Eta(), this_weight_Loose, -3, 3, 30);
                FillHist("SingleMuonTrigger_Dijet_pt_F0", muon.Pt(), this_weight_Loose, 0., 200., 200);
                FillHist("SingleMuonTrigger_Dijet_RelIso_F0", LeptonRelIso, this_weight_Loose, 0., 1., 100);
                FillHist("SingleMuonTrigger_Dijet_Chi2_F0", muon.GlobalChi2(), this_weight_Loose, 0., 10., 10);
                FillHist("SingleMuonTrigger_Dijet_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
                FillHist("SingleMuonTrigger_Dijet_HT_loose_F0", HT_loose, this_weight_Loose, 0, 300, 300);
                FillHist("SingleMuonTrigger_Dijet_HT_tag_F0", HT_tag, this_weight_Loose, 0, 300, 300);
                if(muontriTightColl.size() == 1){
                  FillHist("SingleMuonTrigger_Dijet_eta_F", muon.Eta(), this_weight_Loose, -3, 3, 30);
                  FillHist("SingleMuonTrigger_Dijet_pt_F", muon.Pt(), this_weight_Loose, 0., 200., 200);
                  FillHist("SingleMuonTrigger_Dijet_RelIso_F", LeptonRelIso, this_weight_Loose, 0., 1., 100);
                  FillHist("SingleMuonTrigger_Dijet_Chi2_F", muon.GlobalChi2(), this_weight_Loose, 0., 10., 10);
                  FillHist("SingleMuonTrigger_Dijet_events_F", muon.Pt(), fabs(muon.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
                  FillHist("SingleMuonTrigger_Dijet_HT_loose_F", HT_loose, this_weight_Loose, 0, 300, 300);
                  FillHist("SingleMuonTrigger_Dijet_HT_tag_F", HT_tag, this_weight_Loose, 0, 300, 300);
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
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_HighdXYLoose, 0, 10, 10);
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 20);
    }

    //==== 4) large dXY muon method
  
    if( muontriHighdXYLooseColl.size() == 1 ){
      snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0); 
      double LeptonRelIso = (HighdXYmuon.SumIsoCHDR03() + std::max(0.0, HighdXYmuon.SumIsoNHDR03() + HighdXYmuon.SumIsoPHDR03() - 0.5* HighdXYmuon.SumPUIsoR03()))/HighdXYmuon.Pt(); 
      FillHist("SingleMuonTrigger_HighdXY_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3, 3, 30);
      FillHist("SingleMuonTrigger_HighdXY_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
      FillHist("SingleMuonTrigger_HighdXY_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
      FillHist("SingleMuonTrigger_HighdXY_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0, 10, 10);
      FillHist("SingleMuonTrigger_HighdXY_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
      if( LeptonRelIso < 0.1 ){
        FillHist("SingleMuonTrigger_HighdXY_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3, 3, 30);
        FillHist("SingleMuonTrigger_HighdXY_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
        FillHist("SingleMuonTrigger_HighdXY_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
        FillHist("SingleMuonTrigger_HighdXY_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0, 10, 10);
        FillHist("SingleMuonTrigger_HighdXY_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
      }
    }

    //==== 5) MC Truth
    //==== here, we pick FAKE muons using GEN info.
    if( !k_isdata ){

      Double_t this_weight_NodXYCutLoose = weight*GetPrescale(muontriNodXYCutLooseColl_raw, lowpt_trigger, highpt_trigger);

      bool HistoFilled = false;

      for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){

        //==== allow histo filled only once
        if(HistoFilled) break;

        snu::KMuon muon = muontriNodXYCutLooseColl_raw.at(i);
        //==== if prompt, skip
        if( muon.GetType() == 0 || muon.GetType() == 7 ) continue;

        double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();

        //==== 1) |dXY| < 1 cm, |dXY/err| > 4.
        if( fabs( muon.dXY() ) < 1. &&
            fabs( muon.dXY() / muon.dXYErrPat() ) > 4. ){

          FillHist("SingleMuonTrigger_MCTruth_HighdXY_eta_F0", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_pt_F0", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_Chi2_F0", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 10, 10);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_eta_F", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_pt_F", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_Chi2_F", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 10, 10);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          }

          HistoFilled = true;

        }

        //==== 2) |dXY| < 0.01 cm, |dXY/err| < 3.
        if( fabs( muon.dXY() ) < 0.01 &&
            fabs( muon.dXY() / muon.dXYErrPat() ) < 3. ){

          FillHist("SingleMuonTrigger_MCTruth_eta_F0", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
          FillHist("SingleMuonTrigger_MCTruth_pt_F0", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
          FillHist("SingleMuonTrigger_MCTruth_RelIso_F0", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist("SingleMuonTrigger_MCTruth_Chi2_F0", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 10, 10);
          FillHist("SingleMuonTrigger_MCTruth_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("SingleMuonTrigger_MCTruth_eta_F", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
            FillHist("SingleMuonTrigger_MCTruth_pt_F", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
            FillHist("SingleMuonTrigger_MCTruth_RelIso_F", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
            FillHist("SingleMuonTrigger_MCTruth_Chi2_F", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 10, 10);
            FillHist("SingleMuonTrigger_MCTruth_events_F", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          }

          HistoFilled = true;

        }


      } // END Muon loop


    } // END k_isdata



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
      FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) continue;
      FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    //==== for DY (prompt matched)
    for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutTightColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
      if( thismuon.GetType() != 0 && thismuon.GetType() != 7 ) continue;
      FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetType() != 0 && thismuon.GetType() != 7 ) continue;
      FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
      FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_D0Err", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.dXYErrPat(), weight, 0., 10., 50);
        FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_D0Err_dXYcut_10mm", fabs( thismuon.dXY() )/thismuon.D0Err(), weight, 0., 10., 50);
      }
    }

    //==== 2) DiMuon HighdXY LooseMuon study

    for(unsigned int i=0; i<muontriHighdXYLooseColl.size(); i++){
      if(muontriHighdXYLooseColl.at(0).Pt() < 20.) break;
      snu::KMuon muon = muontriHighdXYLooseColl.at(i);
      double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_RelIso", LeptonRelIso, weight, 0., 1., 100);
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), weight, 0, 10, 10);
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), weight, 0., 0.2, 20);
    }

    //==== 3) Filling num and den for FR

    float etaarray_2 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
    float ptarray_2 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

    int n_jets = jetColl_lepveto.size();

    //==== Two muons
    if( muontriHighdXYLooseColl.size() == 2){
      if( muontriHighdXYLooseColl.at(0).Pt() > 20. ){

        double dimuon_mass = ( muontriHighdXYLooseColl.at(0) + muontriHighdXYLooseColl.at(1) ).M();
        bool isOS = ( muontriHighdXYLooseColl.at(0).Charge() != muontriHighdXYLooseColl.at(1).Charge() );
        bool isOSandZresonance = false;
        if(isOS){
          if(dimuon_mass > 60 && dimuon_mass < 120) isOSandZresonance = true;
        }

        FillHist("DiMuonTrigger_HighdXY_dRdimuon", muontriHighdXYLooseColl.at(0).DeltaR( muontriHighdXYLooseColl.at(1) ), weight, 0, 4, 40);
        FillHist("DiMuonTrigger_HighdXY_mdimuon", dimuon_mass, weight, 0, 200, 40);
        FillHist("DiMuonTrigger_HighdXY_n_jets", n_jets, weight, 0, 10, 10);

        if( dimuon_mass > 15 && !isOSandZresonance ){
          FillHist("DiMuonTrigger_HighdXY_mdimuon_after_cut", dimuon_mass, weight, 0, 200, 40);
          int n_tight=0;
          for(unsigned int i=0; i<2; i++){

            snu::KMuon muon = muontriHighdXYLooseColl.at(i);
            double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
            FillHist("DiMuonTrigger_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_HighdXY_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_HighdXY_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("DiMuonTrigger_HighdXY_n_jets_F0", n_jets, weight, 0, 10, 10);
            FillHist("DiMuonTrigger_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            if( LeptonRelIso < 0.1 ){
              n_tight++;
              FillHist("DiMuonTrigger_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_n_jets_F", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            }

            if( n_jets == 0){
              FillHist("DiMuonTrigger_HighdXY_0jet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_0jet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_0jet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_0jet_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_0jet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_0jet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist("DiMuonTrigger_HighdXY_0jet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuonTrigger_HighdXY_0jet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuonTrigger_HighdXY_0jet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuonTrigger_HighdXY_0jet_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
                FillHist("DiMuonTrigger_HighdXY_0jet_n_jets_F", n_jets, weight, 0, 10, 10);
                FillHist("DiMuonTrigger_HighdXY_0jet_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              }
            }
            if( n_jets != 0){
              FillHist("DiMuonTrigger_HighdXY_withjet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_withjet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_withjet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_withjet_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_withjet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_withjet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( LeptonRelIso < 0.1 ){
                FillHist("DiMuonTrigger_HighdXY_withjet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuonTrigger_HighdXY_withjet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuonTrigger_HighdXY_withjet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuonTrigger_HighdXY_withjet_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
                FillHist("DiMuonTrigger_HighdXY_withjet_n_jets_F", n_jets, weight, 0, 10, 10);
                FillHist("DiMuonTrigger_HighdXY_withjet_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              }
            }
          }
          FillHist("DiMuonTrigger_HighdXY_n_tight", n_tight, weight, 0., 3., 3);

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

        double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();

        //==== 1) |dXY| < 1 cm, |dXY/err| > 4.
        if( fabs( muon.dXY() ) < 1. &&
            fabs( muon.dXY() / muon.dXYErrPat() ) > 4. ){

          FillHist("DiMuonTrigger_MCTruth_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("DiMuonTrigger_MCTruth_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }

        }

        //==== 2) |dXY| < 0.01 cm, |dXY/err| < 3.
        if( fabs( muon.dXY() ) < 0.01 &&
            fabs( muon.dXY() / muon.dXYErrPat() ) < 3. ){

          FillHist("DiMuonTrigger_MCTruth_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("DiMuonTrigger_MCTruth_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("DiMuonTrigger_MCTruth_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("DiMuonTrigger_MCTruth_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
          FillHist("DiMuonTrigger_MCTruth_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.1 ){
            FillHist("DiMuonTrigger_MCTruth_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_MCTruth_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_MCTruth_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_MCTruth_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("DiMuonTrigger_MCTruth_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }

        }

      }


    } // END !k_isdata


    //==== 5) Tagged Z + additional muon

    std::vector<snu::KMuon> TagZ_muons;
    for(unsigned int i=0; i<muontriNodXYCutLooseColl.size(); i++){
      snu::KMuon muon = muontriNodXYCutLooseColl.at(i);
      if( fabs(muon.dXY()) < 0.2 && fabs( muon.dXY() / muon.dXYErrPat() ) < 3. ) TagZ_muons.push_back(muon); 
    }

    if( TagZ_muons.size() == 3 ){
     int sum_charge = TagZ_muons.at(0).Charge() + TagZ_muons.at(1).Charge() + TagZ_muons.at(2).Charge();
     if( abs(sum_charge) != 3 && TagZ_muons.at(0).Pt() > 20. ){

				double m_z = 91.1876;
        //==== find opp sign and SS
        snu::KMuon muon_OS, muon_SS[2];
        if     ( TagZ_muons.at(0).Charge()==TagZ_muons.at(1).Charge() ){
          muon_SS[0] = TagZ_muons.at(0);
          muon_SS[1] = TagZ_muons.at(1);
          muon_OS    = TagZ_muons.at(2);
        }
        else if( TagZ_muons.at(0).Charge()==TagZ_muons.at(2).Charge() ){
          muon_SS[0] = TagZ_muons.at(0);
          muon_SS[1] = TagZ_muons.at(2);
          muon_OS    = TagZ_muons.at(1);
        }
        else if( TagZ_muons.at(1).Charge()==TagZ_muons.at(2).Charge() ){
          muon_SS[0] = TagZ_muons.at(1);
          muon_SS[1] = TagZ_muons.at(2);
          muon_OS    = TagZ_muons.at(0);
        }
        else{
          Message("This should not happen.." , INFO);
        }

        double m_z_candidate[2];
        bool both_tight[2];
        for(int i=0; i<2; i++){
          m_z_candidate[i] = (muon_SS[i] + muon_OS).M();  
          both_tight[i] = ( fabs(muon_SS[i].dXY())<0.01 && muon_SS[i].LeptonRelIso()<0.1 ) && 
                          ( fabs(   muon_OS.dXY())<0.01 &&    muon_OS.LeptonRelIso()<0.1 );
        }

				//==== z region test
				double z_width[6] = {5., 10., 15., 20., 25., 30.};
				for(unsigned int i=0; i<6; i++){
          TString str_this_width = TString::Itoa(z_width[i], 10);
          int index_z_candidate = 0;
          int n_z_candidate = 0;
          double temp_m = 0.;
					for(int j=0; j<2; j++){
						if( fabs(m_z_candidate[j] - m_z) < z_width[i] && both_tight[i]  ){
              index_z_candidate = j;
              n_z_candidate++;
              temp_m = m_z_candidate[j];
						}
					}

          FillHist("DiMuonTrigger_TagZ_n_m_z_candidate_width_"+str_this_width, n_z_candidate, weight, 0, 4, 4);
          if(n_z_candidate==1){
            FillHist("DiMuonTrigger_TagZ_m_z_candidate_width_"+str_this_width, temp_m, weight, 0, 200., 200);

            float etaarray_3 [] = {0.0, 1.479, 2.5};
            float ptarray_3 [] = {10., 40., 50., 60., 100.};

            snu::KMuon muon = TagZ_muons.at(index_z_candidate);
            double LeptonRelIso = (muon.SumIsoCHDR03() + std::max(0.0, muon.SumIsoNHDR03() + muon.SumIsoPHDR03() - 0.5* muon.SumPUIsoR03()))/muon.Pt();
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_eta_F0", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_pt_F0", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_Chi2_F0", muon.GlobalChi2(), weight, 0, 10, 10);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_n_jets_F0", n_jets, weight, 0, 10, 10);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_3, 4, etaarray_3, 2);
            if( fabs( muon.dXY() ) < 0.01 && LeptonRelIso < 0.1 ){
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_eta_F", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_pt_F", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_Chi2_F", muon.GlobalChi2(), weight, 0, 10, 10);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_n_jets_F", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_3, 4, etaarray_3, 2);
            }
            
          }

				}



			}

    }

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
    
    //==== DoubleMuon
    //====   HLT_Mu17 -> 20 < pT <= 25
    //====   HLT_Mu8  -> 10 < pT <= 20
    //==== SingllueMuon
    //====   HLT_Mu24 -> 25 < pT
    //====   HLT_Mu12 -> 10 < pT <= 20

/*
    bool high_pt_pass(false), low_pt_pass(false);
    float high_prescale(0.), low_prescale(0.);

    if( k_jskim_flag_1 == "DoubleMuon" ){
      high_pt_pass = muon.at(0).Pt() > 20. && muon.at(0).Pt() <= 25.;
      high_prescale = (25.014) / 19674;
      low_pt_pass = muon.at(0).Pt() <= 20.;
      low_prescale = (2.021276) / 19674;
    }
    if( k_jskim_flag_1 == "SingleMuon" ){
      high_pt_pass = muon.at(0).Pt() > 25.;
      high_prescale = (93.985) / 19674;
      low_pt_pass = muon.at(0).Pt() <= 20.;
      low_prescale = (2.369835) / 19674;
    }

    if( high_pt_pass ){
      if( passhigh ){
        prescale_trigger = high_prescale;
      }
      else prescale_trigger = 0.;
    }
    else if( low_pt_pass ){
      if( passlow ){
        prescale_trigger = low_prescale;
      }
      else prescale_trigger = 0.;
    }
*/

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


  } // END ONE Muon loop

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}
