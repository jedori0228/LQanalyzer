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


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
FakeRateCalculator_Mu::FakeRateCalculator_Mu() :  AnalyzerCore(), out_muons(0) {
  
  rmcor = new rochcor2015();
  
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

  
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  

  /// Acts on data to remove bad reconstructed event 
  if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;
  

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment
   
  if(!PassBasicEventCuts()) return;     /// Initial event cuts : 

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  //==== DoubleMuon DataSet, SingleMuon trigger  
  std::vector<TString> triggerlist_Mu8, triggerlist_Mu17;
  triggerlist_Mu8.push_back("HLT_Mu8_v");
  triggerlist_Mu17.push_back("HLT_Mu17_v");
  //==== DoubleMuon DataSet, DoubleMuon trigger
  std::vector<TString> triggerlist_Mu17TkMu8;
  triggerlist_Mu17TkMu8.push_back("HLT_Mu17_TkMu8_v");
 
  //float trigger_ps_weight= ApplyPrescale("HLT_IsoMu20", TargetLumi,lumimask);
   
  // Trigger matching is done using KMuon::TriggerMatched(TString) which returns a bool

  /* // #### CAT::: trigger matching information is stored for muons and electrons for:
  ///HLT_IsoMu24_eta2p1_v
  ///HLT_Mu17_Mu8_DZ_v
  ///HLT_Mu17_TkMu8_DZ_v
  ///HLT_IsoMu20
  ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
  ///HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele12_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
  ///HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_
  ///HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v
  ///HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v
  */

  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;


  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
   

  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed

  //==== TriMuon signal muons
  std::vector<snu::KMuon> muontriTightColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_TIGHT);
  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_LOOSE);
  //==== QCD muons
  std::vector<snu::KMuon> muontriHighdXYTightColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_HIGHDXY_TIGHT);
  std::vector<snu::KMuon> muontriHighdXYLooseColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_HIGHDXY_LOOSE);
  //==== No dXY cut muons
  std::vector<snu::KMuon> muontriNodXYCutTightColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_NODXYCUT_TIGHT);
  std::vector<snu::KMuon> muontriNodXYCutLooseColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_NODXYCUT_LOOSE);
   
  CorrectMuonMomentum(muonTightColl);
  float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonTightColl,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl             = GetJets(BaseSelection::JET_NOLEPTONVETO); // All jets
  std::vector<snu::KJet> jetColl_loose       = GetJets(BaseSelection::JET_LOOSE); // pt > 10; eta < 5. ; PFlep veto
  std::vector<snu::KJet> jetColl_tight       = GetJets(BaseSelection::JET_TIGHT);// pt > 20 ; eta < 2.5; PFlep veto
  std::vector<snu::KJet> jetColl_hn          = GetJets(BaseSelection::JET_HN);// pt > 20 ; eta < 2.5; PFlep veto; pileup ID
   
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");
   
  numberVertices = eventbase->GetEvent().nVertices();   
   
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    pileup_reweight = eventbase->GetEvent().PileUpWeight(lumimask);

  }
   
  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);

  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    weight*=pileup_reweight;
    //weight*=weight_trigger_sf;
    //weight*=trigger_ps_weight;
  }


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
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriLooseColl_raw.at(i);
    if(isData) muontriLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriLooseColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYTightColl_raw.at(i);
    if(isData) muontriHighdXYTightColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriHighdXYTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYLooseColl_raw.at(i);
    if(isData) muontriHighdXYLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriHighdXYLooseColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
    if(isData) muontriNodXYCutTightColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriNodXYCutTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
    if(isData) muontriNodXYCutLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriNodXYCutLooseColl.push_back(thismuon);
    }
  }

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  //=========================
  //==== SingleMuon Trigger
  //=========================

  bool highpt_trigger(false), lowpt_trigger(false);
  highpt_trigger = PassTrigger(triggerlist_Mu17,prescale);
  lowpt_trigger = PassTrigger(triggerlist_Mu8,prescale);

  if( highpt_trigger || lowpt_trigger ){

    Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));
    Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));

    //==== 1) LooseMuon study

    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      snu::KMuon muon = muontriLooseColl.at(i);
      double LeptonRelIso = muon.RelIso04();
      FillHist("SingleMuonTrigger_LooseMuon_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
      FillHist("SingleMuonTrigger_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
      FillHist("SingleMuonTrigger_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 20);
      FillHist("SingleMuonTrigger_LooseMuon_dZ", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 20);
    }

    //==== 2) back-to-back dijet topology

    //==== tag jet collections
    //==== pt > 40 GeV
    //==== |eta| > 2.4
    //==== LeptonVeto
    std::vector<snu::KJet> jetColl;
    eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
    eventbase->GetJetSel()->SetPt(40.);
    eventbase->GetJetSel()->SetEta(2.4);
    eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl, GetMuons(BaseSelection::MUON_HN_VETO), GetElectrons(false,false, BaseSelection::ELECTRON_HN_VETO));

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
                double HT_loose = AnalyzerCore::SumPt(jetColl_loose);
                double HT_tag = AnalyzerCore::SumPt(jetColl);
                double LeptonRelIso = muon.RelIso04();
                FillHist("SingleMuonTrigger_Dijet_eta_F0", muon.Eta(), this_weight_Loose, -3, 3, 30);
                FillHist("SingleMuonTrigger_Dijet_pt_F0", muon.Pt(), this_weight_Loose, 0., 200., 200);
                FillHist("SingleMuonTrigger_Dijet_RelIso_F0", LeptonRelIso, this_weight_Loose, 0., 1., 100);
                FillHist("SingleMuonTrigger_Dijet_Chi2_F0", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
                FillHist("SingleMuonTrigger_Dijet_dXY_F0", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 40);
                FillHist("SingleMuonTrigger_Dijet_dZ_F0", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 40);
                FillHist("SingleMuonTrigger_Dijet_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
                FillHist("SingleMuonTrigger_Dijet_HT_loose_F0", HT_loose, this_weight_Loose, 0, 300, 300);
                FillHist("SingleMuonTrigger_Dijet_HT_tag_F0", HT_tag, this_weight_Loose, 0, 300, 300);
                if(muontriTightColl.size() == 1){
                  FillHist("SingleMuonTrigger_Dijet_eta_F", muon.Eta(), this_weight_Loose, -3, 3, 30);
                  FillHist("SingleMuonTrigger_Dijet_pt_F", muon.Pt(), this_weight_Loose, 0., 200., 200);
                  FillHist("SingleMuonTrigger_Dijet_RelIso_F", LeptonRelIso, this_weight_Loose, 0., 1., 100);
                  FillHist("SingleMuonTrigger_Dijet_Chi2_F", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
                  FillHist("SingleMuonTrigger_Dijet_dXY_F", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 40);
                  FillHist("SingleMuonTrigger_Dijet_dZ_F", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 40);
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
      double LeptonRelIso = muon.RelIso04();
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 20);
      FillHist("SingleMuonTrigger_HighdXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 20);
    }

    //==== 4) large dXY muon method

    if( muontriHighdXYLooseColl.size() == 1 ){
      snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0);
      double LeptonRelIso = HighdXYmuon.RelIso04();
      FillHist("SingleMuonTrigger_HighdXY_eta_F0", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3, 3, 30);
      FillHist("SingleMuonTrigger_HighdXY_pt_F0", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
      FillHist("SingleMuonTrigger_HighdXY_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
      FillHist("SingleMuonTrigger_HighdXY_Chi2_F0", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
      FillHist("SingleMuonTrigger_HighdXY_dXY_F0", fabs(HighdXYmuon.dXY()), this_weight_Loose, 0., 0.2, 40);
      FillHist("SingleMuonTrigger_HighdXY_dZ_F0", fabs(HighdXYmuon.dZ()), this_weight_Loose, 0., 0.2, 40);
      FillHist("SingleMuonTrigger_HighdXY_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
      if( LeptonRelIso < 0.05 ){
        FillHist("SingleMuonTrigger_HighdXY_eta_F", HighdXYmuon.Eta(), this_weight_HighdXYLoose, -3, 3, 30);
        FillHist("SingleMuonTrigger_HighdXY_pt_F", HighdXYmuon.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
        FillHist("SingleMuonTrigger_HighdXY_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
        FillHist("SingleMuonTrigger_HighdXY_Chi2_F", HighdXYmuon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
        FillHist("SingleMuonTrigger_HighdXY_dXY_F", fabs(HighdXYmuon.dXY()), this_weight_Loose, 0., 0.2, 40);
        FillHist("SingleMuonTrigger_HighdXY_dZ_F", fabs(HighdXYmuon.dZ()), this_weight_Loose, 0., 0.2, 40);
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
        if( muon.GetParticleType() == snu::KMuon::PROMPT ) continue;

        double LeptonRelIso = muon.RelIso04();

        //==== 1) |dXY| < 1 cm, |dXY/err| > 4.
        if( fabs( muon.dXY() ) < 1. &&
            fabs( muon.dXYSig() ) > 4. ){

          FillHist("SingleMuonTrigger_MCTruth_HighdXY_eta_F0", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_pt_F0", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_Chi2_F0", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_dXY_F0", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 40);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_dZ_F0", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 40);
          FillHist("SingleMuonTrigger_MCTruth_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < 0.05 ){
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_eta_F", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_pt_F", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_Chi2_F", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_dXY_F", fabs(muon.dXY()), this_weight_Loose, 0., 0.2, 40);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_dZ_F", fabs(muon.dZ()), this_weight_Loose, 0., 0.2, 40);
            FillHist("SingleMuonTrigger_MCTruth_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          }

          HistoFilled = true;

        }

        //==== 2) |dXY| < 0.005 cm, |dXY/err| < 3.
        if( fabs( muon.dXY() ) < 0.005 &&
            fabs( muon.dXYSig() ) < 3. ){

          FillHist("SingleMuonTrigger_MCTruth_eta_F0", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
          FillHist("SingleMuonTrigger_MCTruth_pt_F0", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
          FillHist("SingleMuonTrigger_MCTruth_RelIso_F0", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist("SingleMuonTrigger_MCTruth_Chi2_F0", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
          FillHist("SingleMuonTrigger_MCTruth_dXY_F0", fabs(muon.dXY()), this_weight_NodXYCutLoose, 0, 0.2, 40);
          FillHist("SingleMuonTrigger_MCTruth_dZ_F0", fabs(muon.dZ()), this_weight_NodXYCutLoose, 0, 0.2, 40);
          FillHist("SingleMuonTrigger_MCTruth_events_F0", muon.Pt(), fabs(muon.Eta()), this_weight_NodXYCutLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < 0.05 ){
            FillHist("SingleMuonTrigger_MCTruth_eta_F", muon.Eta(), this_weight_NodXYCutLoose, -3, 3, 30);
            FillHist("SingleMuonTrigger_MCTruth_pt_F", muon.Pt(), this_weight_NodXYCutLoose, 0., 200., 200);
            FillHist("SingleMuonTrigger_MCTruth_RelIso_F", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
            FillHist("SingleMuonTrigger_MCTruth_Chi2_F", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
            FillHist("SingleMuonTrigger_MCTruth_dXY_F", fabs(muon.dXY()), this_weight_NodXYCutLoose, 0, 0.2, 40);
            FillHist("SingleMuonTrigger_MCTruth_dZ_F", fabs(muon.dZ()), this_weight_NodXYCutLoose, 0, 0.2, 40);
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
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) continue;
      FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_TightIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) continue;
      FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_dXYErrPat", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_LooseIsoMuon_fake_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      }
    }
    //==== for DY (prompt matched)
    for(unsigned int i=0; i<muontriNodXYCutTightColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutTightColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutTightColl_raw.at(i);
      if( thismuon.GetParticleType() != snu::KMuon::PROMPT ) continue;
      FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_TightIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      }
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
      //==== lead pT > 20 GeV
      if(muontriNodXYCutLooseColl_raw.at(0).Pt() < 20.) break;

      snu::KMuon thismuon = muontriNodXYCutLooseColl_raw.at(i);
      if( thismuon.GetParticleType() != snu::KMuon::PROMPT ) continue;
      FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY", fabs( thismuon.dXY() ), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_dXYErrPat", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      if( fabs( thismuon.dXY() ) < 1. ){
        FillHist("DiMuonTrigger_LooseIsoMuon_prompt_dXY_over_dXYErrPat_dXYcut_10mm", fabs( thismuon.dXYSig() ), weight, 0., 10., 50);
      }
    }

    //==== 2) DiMuon HighdXY LooseMuon study

    for(unsigned int i=0; i<muontriHighdXYLooseColl.size(); i++){
      if(muontriHighdXYLooseColl.at(0).Pt() < 20.) break;
      snu::KMuon muon = muontriHighdXYLooseColl.at(i);
      double LeptonRelIso = muon.RelIso04();
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_RelIso", LeptonRelIso, weight, 0., 1., 100);
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), weight, 0, 50., 50);
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), weight, 0., 0.2, 20);
      FillHist("DiMuonTrigger_HighdXY_LooseMuon_dZ", fabs(muon.dZ()), weight, 0., 0.2, 20);
    }

    //==== 3) Filling num and den for FR

    float etaarray_2 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
    float ptarray_2 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

    int n_jets = jetColl_loose.size();

    //==== Two muons
    if( muontriHighdXYLooseColl.size() == 2){
      if( muontriHighdXYLooseColl.at(0).Pt() > 20. ){

        double dimuon_mass = ( muontriHighdXYLooseColl.at(0) + muontriHighdXYLooseColl.at(1) ).M();
        bool isOS = ( muontriHighdXYLooseColl.at(0).Charge() != muontriHighdXYLooseColl.at(1).Charge() );
        bool isOSandZresonance = false;
        if(isOS){
          if(dimuon_mass > 60 && dimuon_mass < 120) isOSandZresonance = true;
        }

/*
        if(dimuon_mass > 35 && dimuon_mass < 45){
          std::vector<snu::KTruth> truthColl;
          eventbase->GetTruthSel()->Selection(truthColl);

          //==== print truth info
          cout << "=========================================================" << endl;
          cout << "truth size = " << truthColl.size() << endl;
          cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
          for(int i=2; i<truthColl.size(); i++){
            cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << endl;
          } 
          return;
        }*/
        FillHist("DiMuonTrigger_HighdXY_dRdimuon", muontriHighdXYLooseColl.at(0).DeltaR( muontriHighdXYLooseColl.at(1) ), weight, 0, 4, 40);
        FillHist("DiMuonTrigger_HighdXY_mdimuon", dimuon_mass, weight, 0, 200, 40);
        if(isOS) FillHist("DiMuonTrigger_HighdXY_OSmdimuon", dimuon_mass, weight, 0, 200, 40);
        else     FillHist("DiMuonTrigger_HighdXY_SSmdimuon", dimuon_mass, weight, 0, 200, 40);
        FillHist("DiMuonTrigger_HighdXY_n_jets", n_jets, weight, 0, 10, 10);

        if( dimuon_mass > 15 && !isOSandZresonance ){
          FillHist("DiMuonTrigger_HighdXY_mdimuon_after_cut", dimuon_mass, weight, 0, 200, 40);
          int n_tight=0;
          for(unsigned int i=0; i<2; i++){

            snu::KMuon muon = muontriHighdXYLooseColl.at(i);
            double LeptonRelIso = muon.RelIso04();
            FillHist("DiMuonTrigger_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_HighdXY_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_HighdXY_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
            FillHist("DiMuonTrigger_HighdXY_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_HighdXY_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_HighdXY_n_jets_F0", n_jets, weight, 0, 10, 10);
            FillHist("DiMuonTrigger_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            if( LeptonRelIso < 0.05 ){
              n_tight++;
              FillHist("DiMuonTrigger_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
              FillHist("DiMuonTrigger_HighdXY_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_n_jets_F", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
            }
           if( n_jets == 0){
              FillHist("DiMuonTrigger_HighdXY_0jet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_0jet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_0jet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_0jet_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
              FillHist("DiMuonTrigger_HighdXY_0jet_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_0jet_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_0jet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_0jet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( LeptonRelIso < 0.05 ){
                FillHist("DiMuonTrigger_HighdXY_0jet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuonTrigger_HighdXY_0jet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuonTrigger_HighdXY_0jet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuonTrigger_HighdXY_0jet_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
                FillHist("DiMuonTrigger_HighdXY_0jet_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
                FillHist("DiMuonTrigger_HighdXY_0jet_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
                FillHist("DiMuonTrigger_HighdXY_0jet_n_jets_F", n_jets, weight, 0, 10, 10);
                FillHist("DiMuonTrigger_HighdXY_0jet_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              }
            }
            if( n_jets != 0){
              FillHist("DiMuonTrigger_HighdXY_withjet_eta_F0", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_HighdXY_withjet_pt_F0", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_HighdXY_withjet_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_HighdXY_withjet_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
              FillHist("DiMuonTrigger_HighdXY_withjet_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_withjet_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_HighdXY_withjet_n_jets_F0", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_HighdXY_withjet_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
              if( LeptonRelIso < 0.05 ){
                FillHist("DiMuonTrigger_HighdXY_withjet_eta_F", muon.Eta(), weight, -3, 3, 30);
                FillHist("DiMuonTrigger_HighdXY_withjet_pt_F", muon.Pt(), weight, 0., 200., 200);
                FillHist("DiMuonTrigger_HighdXY_withjet_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
                FillHist("DiMuonTrigger_HighdXY_withjet_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
                FillHist("DiMuonTrigger_HighdXY_withjet_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
                FillHist("DiMuonTrigger_HighdXY_withjet_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
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
        if( muon.GetParticleType() == snu::KMuon::PROMPT ) continue;

        double LeptonRelIso = muon.RelIso04();

        //==== 1) |dXY| < 1 cm, |dXY/err| > 4.
        if( fabs( muon.dXY() ) < 1. &&
            fabs( muon.dXYSig() ) > 4. ){

          FillHist("DiMuonTrigger_MCTruth_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
          FillHist("DiMuonTrigger_MCTruth_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.05 ){
            FillHist("DiMuonTrigger_MCTruth_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_MCTruth_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }

        }
        //==== 2) |dXY| < 0.005 cm, |dXY/err| < 3.
        if( fabs( muon.dXY() ) < 0.005 &&
            fabs( muon.dXYSig() ) < 3. ){

          FillHist("DiMuonTrigger_MCTruth_eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("DiMuonTrigger_MCTruth_pt_F0", muon.Pt(), weight, 0., 200., 200);
          FillHist("DiMuonTrigger_MCTruth_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
          FillHist("DiMuonTrigger_MCTruth_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
          FillHist("DiMuonTrigger_MCTruth_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
          FillHist("DiMuonTrigger_MCTruth_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
          FillHist("DiMuonTrigger_MCTruth_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          if( LeptonRelIso < 0.05 ){
            FillHist("DiMuonTrigger_MCTruth_eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_MCTruth_pt_F", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_MCTruth_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_MCTruth_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
            FillHist("DiMuonTrigger_MCTruth_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_MCTruth_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_MCTruth_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 4);
          }

        }

      }


    } // END !k_isdata

    //==== 5) Tagged Z + additional muon

    std::vector<snu::KMuon> TagZ_muons;
    for(unsigned int i=0; i<muontriNodXYCutLooseColl.size(); i++){
      snu::KMuon muon = muontriNodXYCutLooseColl.at(i);
      if( fabs(muon.dXY()) < 0.2 && fabs( muon.dXYSig() ) < 3. ) TagZ_muons.push_back(muon);
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
          both_tight[i] = ( fabs(muon_SS[i].dXY())<0.005 && muon_SS[i].RelIso04()<0.05 ) &&
                          ( fabs(   muon_OS.dXY())<0.005 &&    muon_OS.RelIso04()<0.05 );
        }
        //==== z region test
        double z_width[6] = {5., 10., 15., 20., 25., 30.};
        for(unsigned int i=0; i<6; i++){
          TString str_this_width = TString::Itoa(z_width[i], 10);
          int index_z_candidate = 0;
          int n_z_candidate = 0;
          double temp_m = 0.;
          for(int j=0; j<2; j++){
            if( fabs(m_z_candidate[j] - m_z) < z_width[i] && both_tight[j]  ){
              index_z_candidate = j;
              n_z_candidate++;
              temp_m = m_z_candidate[j];
            }
          }

          FillHist("DiMuonTrigger_TagZ_n_m_z_candidate_width_"+str_this_width, n_z_candidate, weight, 0, 4, 4);
          if(n_z_candidate==1){
            FillHist("DiMuonTrigger_TagZ_m_z_candidate_width_"+str_this_width, temp_m, weight, 0, 200., 200);

            //float etaarray_3 [] = {0.0, 1.479, 2.5};
            //float ptarray_3 [] = {10., 40., 50., 60., 100.};
           float etaarray_3 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
           float ptarray_3 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

            snu::KMuon muon;
            if     (index_z_candidate==0)  muon = muon_SS[1];
            else if(index_z_candidate==1)  muon = muon_SS[0];

            double LeptonRelIso = muon.RelIso04();
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_eta_F0", muon.Eta(), weight, -3, 3, 30);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_pt_F0", muon.Pt(), weight, 0., 200., 200);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_RelIso_F0", LeptonRelIso, weight, 0., 1., 100);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_Chi2_F0", muon.GlobalChi2(), weight, 0, 50., 50);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_dXY_F0", fabs(muon.dXY()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_dZ_F0", fabs(muon.dZ()), weight, 0, 0.2, 40);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_n_jets_F0", n_jets, weight, 0, 10, 10);
            FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_3, 9, etaarray_3, 4);
            if( fabs( muon.dXY() ) < 0.005 && LeptonRelIso < 0.05 ){
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_eta_F", muon.Eta(), weight, -3, 3, 30);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_pt_F", muon.Pt(), weight, 0., 200., 200);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_RelIso_F", LeptonRelIso, weight, 0., 1., 100);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_Chi2_F", muon.GlobalChi2(), weight, 0, 50., 50);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_dXY_F", fabs(muon.dXY()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_dZ_F", fabs(muon.dZ()), weight, 0, 0.2, 40);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_n_jets_F", n_jets, weight, 0, 10, 10);
              FillHist("DiMuonTrigger_TagZ_width_"+str_this_width+"_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_3, 9, etaarray_3, 4);
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
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

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
  if(!k_isdata)delete reweightPU;
  
}


void FakeRateCalculator_Mu::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 7,0.,7.);

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
        prescale_trigger = ApplyPrescale("HLT_Mu17", TargetLumi,lumimask) ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = ApplyPrescale("HLT_Mu8", TargetLumi,lumimask) ;
      }
      else prescale_trigger = 0.;
    }


  } // END ONE Muon loop

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}

