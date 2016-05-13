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
  //==== Large dXY muons
  std::vector<snu::KMuon> muontriHighdXYTightColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYTightMuonSelection(muontriHighdXYTightColl_raw);
  std::vector<snu::KMuon> muontriHighdXYLooseColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYLooseMuonSelection(muontriHighdXYLooseColl_raw);

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

    //==== 1) back-to-back dijet topology

    //==== tag jet collections
    std::vector<snu::KJet> jetColl = GetJets("HNtriFRTagJet");

    //if(jetColl.size() == 0) return;
    if(jetColl.size() != 0){
      //if(muontriLooseColl.size() != 1) return;
      if(muontriLooseColl.size() == 1){

        //snu::KEvent Evt = eventbase->GetEvent();
        //double MET = Evt.PFMET();
        //if( MET < 40 ) return; // Let Wjets dominate

        snu::KParticle muon;
        muon = muontriLooseColl.at(0);
        muon.SetPxPyPzE(muon.Px(),muon.Py(),muon.Pz(),muon.E());

        float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;

        for(unsigned int i=0; i<jetColl.size(); i++){
          dR = jetColl.at(i).DeltaR(muon);
          dPhi = fabs(jetColl.at(i).DeltaPhi(muon));
          ptmu_ptjet = muon.Pt() / jetColl.at(i).Pt();
          FillHist("dR", dR, this_weight, 0., 10., 100);

          if( dR > 1.0 ){
            FillHist("ptmu_ptjet", ptmu_ptjet, this_weight, 0., 2., 200);
            FillHist("dPhi", dPhi, this_weight, 0., 4., 40);
            if( ptmu_ptjet < 1. ){
              if( dPhi > 2.5 ){
                double HT_loose = AnalyzerCore::SumPt(jetColl_lepveto);
                double HT_tag = AnalyzerCore::SumPt(jetColl);
                FillHist("eta_F0", muon.Eta(), this_weight, -3, 3, 30);
                FillHist("pt_F0", muon.Pt(), this_weight, 0., 200., 200./1.);
                FillHist("events_F0", muon.Pt(), fabs(muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
                FillHist("HT_loose_F0", HT_loose, this_weight, 0, 300, 300);
                FillHist("HT_tag_F0", HT_tag, this_weight, 0, 300, 300);
                if(muontriTightColl.size() == 1){
                  FillHist("eta_F", muon.Eta(), this_weight, -3, 3, 30);
                  FillHist("pt_F", muon.Pt(), this_weight, 0., 200., 200/1.);
                  FillHist("events_F", muon.Pt(), fabs(muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
                  FillHist("HT_loose_F", HT_loose, this_weight, 0, 300, 300);
                  FillHist("HT_tag_F", HT_tag, this_weight, 0, 300, 300);
                }
                goto stop;
              }
            }
          }
        }
        stop: ;
  
      }
    }

    //==== 2) large dXY muon method
  
    if( muontriHighdXYLooseColl.size() == 1 ){
      snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0); 
    
      FillHist("HighdXY_eta_F0", HighdXYmuon.Eta(), this_weight, -3, 3, 30);
      FillHist("HighdXY_pt_F0", HighdXYmuon.Pt(), this_weight, 0., 200., 200./1.);
      FillHist("HighdXY_events_F0", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight, ptarray, 9, etaarray, 4);
      if(muontriHighdXYTightColl.size() == 1){
        FillHist("HighdXY_eta_F", HighdXYmuon.Eta(), this_weight, -3, 3, 30);
        FillHist("HighdXY_pt_F", HighdXYmuon.Pt(), this_weight, 0., 200., 200/1.);
        FillHist("HighdXY_events_F", HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight, ptarray, 9, etaarray, 4);
      }
    }

    std::vector<snu::KMuon> muontriNodXYCutTightColl;
    eventbase->GetMuonSel()->HNtriNodXYCutTightMuonSelection(muontriNodXYCutTightColl);
    std::vector<snu::KMuon> muontriNodXYCutLooseColl;
    eventbase->GetMuonSel()->HNtriNodXYCutLooseMuonSelection(muontriNodXYCutLooseColl);

    for(unsigned int i=0; i<muontriNodXYCutTightColl.size(); i++){
      FillHist("TightIsoMuon_dXY", muontriNodXYCutTightColl.at(i).dXY(), this_weight, 0., 0.2, 0.2/0.01);
    }
    for(unsigned int i=0; i<muontriNodXYCutLooseColl.size(); i++){
      FillHist("LooseIsoMuon_dXY", muontriNodXYCutLooseColl.at(i).dXY(), this_weight, 0., 0.2, 0.2/0.01);
    }



  } // SingleMuon trigger fired



  //=====================
  //==== DiMuon Trigger
  //=====================

  if( PassTrigger(triggerlist_Mu17TkMu8,prescale) ){

    //float etaarray_2 [] = {0.0, 0.8, 1.479, 2.0, 2.5};
    float etaarray_2 [] = {0.0, 0.8, 1.479, 2.5};
    float ptarray_2 [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};
    //float ptarray_2 [] = {10,15,20,25,30,35,40,45,50,55,60,80,100};

    if( muontriHighdXYLooseColl.size() == 2){
      FillHist("DiMuon_HighdXY_dRdimuon", muontriHighdXYLooseColl.at(0).DeltaR( muontriHighdXYLooseColl.at(1) ), weight, 0, 4, 4./0.1);
      FillHist("DiMuon_HighdXY_mdimuon", (muontriHighdXYLooseColl.at(0)+muontriHighdXYLooseColl.at(1)).M(), weight, 0, 200, 200./5.);
      for(unsigned int i=0; i<2; i++){
        snu::KMuon muon = muontriHighdXYLooseColl.at(i);
        FillHist("DiMuon_HighdXY_eta_F0", muon.Eta(), weight, -3, 3, 30);
        FillHist("DiMuon_HighdXY_pt_F0", muon.Pt(), weight, 0., 200., 200./1.);
        FillHist("DiMuon_HighdXY_events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 3);
        if( eventbase->GetMuonSel()->HNIstriHighdXYTight(muon, false) ){
          FillHist("DiMuon_HighdXY_eta_F", muon.Eta(), weight, -3, 3, 30);
          FillHist("DiMuon_HighdXY_pt_F", muon.Pt(), weight, 0., 200., 200/1.);
          FillHist("DiMuon_HighdXY_events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray_2, 9, etaarray_2, 3);
        }
      } // two muon iterator

    } // if(dimuon)


  } // DiMuon trigger fired



/*
  //////////////////////////////////////
  ///////////// use MC Turth ///////////
  //////////////////////////////////////

  int n_triTight_muons = muontriTightColl.size();
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_lepveto.size();

  //FillHist("n_loose_muons", n_triLoose_muons, 1, 0, 10, 10);
  //FillHist("n_tight_muons", n_triTight_muons, 1, 0, 10, 10);
  //FillHist("n_jets", n_jets, 1, 0, 10, 10);

  if( n_triLoose_muons != 1 ) return;
  snu::KMuon muon = muontriLooseColl.at(0);
  //if( muon.GetType() == 1 || muon.GetType() == 2 || muon.GetType() == 3 ) return;
  if( muon.GetType() == 0 ) return;

  FillHist("eta_den", muon.Eta(), weight, -3, 3, 30);
  FillHist("pt_den", muon.Pt(), weight, 0., 200., 200./1.);
  FillHist("events_den", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);
  FillHist("n_jets_den", n_jets, 1, 0, 15, 15);
  double HT = AnalyzerCore::SumPt(jetColl_lepveto);
  FillHist("HT_den", HT, weight, 0, 300, 300);
  if( n_triTight_muons == 1 ){
    FillHist("eta_num", muon.Eta(), weight, -3, 3, 30);
    FillHist("pt_num", muon.Pt(), weight, 0., 200., 200./1.);
    FillHist("events_num", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);
    FillHist("n_jets_num", n_jets, 1, 0, 15, 15);
    FillHist("HT_num", HT, weight, 0, 300, 300);
  }
*/

/*
  ////////////////////////////////////////////////
  ///////////// impact parameter study ///////////
  ////////////////////////////////////////////////

  if( !PassTrigger(triggerlist_Mu8,prescale) && !PassTrigger(triggerlist_Mu17,prescale) ) return;
  float prescale_trigger = GetPrescale(muontriLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));

  weight *= prescale_trigger;

  vector<snu::KMuon> selected_muons;
  FillHist("events", 0, weight, 0, 2, 2);
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    snu::KMuon thismuon = muontriLooseColl.at(i);
    float LeptonRelIso = (thismuon.SumIsoCHDR03() + std::max(0.0, thismuon.SumIsoNHDR03() + thismuon.SumIsoPHDR03() - 0.5* thismuon.SumPUIsoR03()))/thismuon.Pt() ;
    FillHist("reliso", LeptonRelIso, weight, 0, 0.6, 0.6/0.01);
    FillHist("dXY", fabs( thismuon.dXY() ), weight, 0, 0.2, 0.2/0.01);
    FillHist("dXYPat", fabs( thismuon.dXYPat() ), weight, 0, 0.2, 0.2/0.01);
    FillHist("dXYPat_over_dXYErrPat", fabs( thismuon.dXYPat()/thismuon.dXYErrPat() ), weight, 0, 5, 5./0.1);
    FillHist("D0", fabs( thismuon.D0() ), weight, 0, 0.2, 0.2/0.01);
    FillHist("D0_over_D0Err", fabs( thismuon.D0()/thismuon.D0Err() ), weight, 0, 5, 5./0.1);
    if(    0.02 < fabs( thismuon.D0() )
        &&        fabs( thismuon.D0() ) < 1 
        && 3.0 < fabs( thismuon.D0()/thismuon.D0Err() )
        &&       fabs( thismuon.D0()/thismuon.D0Err() ) < 5.0 ){

      selected_muons.push_back(thismuon);

    }
  }

  if(selected_muons.size() != 0){
    FillHist("events", 1, weight, 0, 2, 2);
  }


  if(selected_muons.size() == 1){
    snu::KParticle muon = selected_muons.at(0);
    FillHist("eta_F0", muon.Eta(), weight, -3, 3, 30);
    FillHist("pt_F0", muon.Pt(), weight, 0., 200., 200./1.);
    FillHist("events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);
    if(muontriTightColl.size() == 1){
      FillHist("eta_F", muon.Eta(), weight, -3, 3, 30);
      FillHist("pt_F", muon.Pt(), weight, 0., 200., 200/1.);
      FillHist("events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);
    }
  }
*/

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
