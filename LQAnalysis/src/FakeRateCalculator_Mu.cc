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

  FillCutFlow("NoCut", weight);

  if(!PassBasicEventCuts()) return; //FIXME Do we need this?
  FillCutFlow("EventCut", weight);
  
  std::vector<TString> triggerlist_Mu5, triggerlist_Mu8, triggerlist_Mu17;
  triggerlist_Mu5.push_back("HLT_Mu5_v");
  triggerlist_Mu8.push_back("HLT_Mu8_v");
  triggerlist_Mu17.push_back("HLT_Mu17_v");

  //if( !PassTrigger(triggerlist_Mu5,prescale) ) return;
  //FillHist("Mu5", 0, 1, 0, 1, 1);
    
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; 
  FillCutFlow("VertexCut", weight);
  
  numberVertices = eventbase->GetEvent().nVertices();  
  
  if (MC_pu&&!k_isdata) {
    weight = weight* reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
  }
  

  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////  

  std::vector<snu::KMuon> muontriTightColl;
  eventbase->GetMuonSel()->HNtriTightMuonSelection(muontriTightColl);
  std::vector<snu::KMuon> muontriLooseColl;
  eventbase->GetMuonSel()->HNtriLooseMuonSelection(muontriLooseColl);


  //FIXME can we use thie scale factor for our muons?
/*
  if(!isData){ 
     for(std::vector<snu::KMuon>::iterator it = muontriTightColl.begin(); it!= muontriTightColl.end(); it++){
       weight *= MuonScaleFactor(it->Eta(), it->Pt(), true);
     }
   }
*/

  std::vector<snu::KMuon> muonTightColl;
  eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl); 
  std::vector<snu::KMuon> muonHighPtColl;
  eventbase->GetMuonSel()->HNTightHighPtMuonSelection(muonHighPtColl); 
  CorrectMuonMomentum(muonTightColl);
  CorrectMuonMomentum(muonHighPtColl);

  std::vector<snu::KMuon> muonLooseColl;
  eventbase->GetMuonSel()->HNLooseMuonSelection(muonLooseColl);
  std::vector<snu::KElectron> electronTightColl;
  eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);

  std::vector<snu::KJet> jetColl = GetJets("fakerate");

  ////////////////////////////////////////////////////
  ///////////////// Analysis
  ////////////////////////////////////////////////////

  if(jetColl.size() == 0) return;

  if(muontriLooseColl.size() != 1) return;

  if( !PassTrigger(triggerlist_Mu8,prescale) && !PassTrigger(triggerlist_Mu17,prescale) ) return;
  //if( ! (PassTrigger(triggerlist_Mu8,prescale) && !PassTrigger(triggerlist_Mu17,prescale)) ) return;

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.PFMET();
  if( MET < 40 ) return; // Let Wjets dominate

  float prescale_trigger = GetPrescale(muontriLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));

  FillHist("prescale_trigger", prescale_trigger, 1, 0, 1, 10000);
  
  cout << prescale_trigger << endl;
  weight *= prescale_trigger;  


  snu::KParticle muon;
  muon = muontriLooseColl.at(0);
  muon.SetPxPyPzE(muon.Px(),muon.Py(),muon.Pz(),muon.E());

  float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;

  float etaarray [] = {0.0,0.8,1.479,2.0,2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  for(unsigned int i=0; i<jetColl.size(); i++){
    dR = jetColl.at(i).DeltaR(muon);
    dPhi = fabs(jetColl.at(i).DeltaPhi(muon));
    ptmu_ptjet = muon.Pt() / jetColl.at(i).Pt();
    FillHist("dR", dR, weight, 0., 10., 100);

    if( dR > 1.0 ){
      FillHist("ptmu_ptjet", ptmu_ptjet, weight, 0., 2., 200);
      FillHist("dPhi", dPhi, weight, 0., 4., 40);
      if( ptmu_ptjet < 1. ){
        if( dPhi > 2.5 ){
          FillHist("eta_F0", muon.Eta(), weight, -3, 3, 30);
          FillHist("pt_F0", muon.Pt(), weight, 0., 200., 200./1.);
          FillHist("events_F0", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);
          if(muontriTightColl.size() == 1){
            FillHist("eta_F", muon.Eta(), weight, -3, 3, 30);
            FillHist("pt_F", muon.Pt(), weight, 0., 200., 200/1.);
            FillHist("events_F", muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);          
          }
          goto stop;
        }
      }
    }
  }
  stop: ;

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
