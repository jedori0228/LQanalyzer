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
   if(!PassBasicEventCuts()) return;    
   
   std::vector<TString> triggerlist_Mu8, triggerlist_Mu17;
   triggerlist_Mu8.push_back("HLT_Mu8_v");
	 triggerlist_Mu17.push_back("HLT_Mu17_v");
      
   m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
   
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; 
   
   numberVertices = eventbase->GetEvent().nVertices();   
   
   if (MC_pu&&!k_isdata) {
     weight = weight* reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
   }
   

   //////////////////////////////////////////////////////
   //////////// Select objetcs
   //////////////////////////////////////////////////////   

   std::vector<snu::KMuon> muonTightColl;
   eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl);
   
   std::vector<snu::KMuon> muonHighPtColl;
   eventbase->GetMuonSel()->HNTightHighPtMuonSelection(muonHighPtColl);
   if(!isData){ 
  	 for(std::vector<snu::KMuon>::iterator it = muonTightColl.begin(); it!= muonTightColl.end(); it++){
	     weight *= MuonScaleFactor(it->Eta(), it->Pt(), true);
   	 }
	 }
   
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
//	if(jetColl.Pt()<40 || jetColl.Eta()<2.4) return;

	FillHist("flow_rate", 0, 1, 0.,2., 2);
  if(muonLooseColl.size() != 1) return;
	FillHist("flow_rate", 1, 1, 0.,2., 2);
	float prescale_trigger = GetPrescale(muonLooseColl, PassTrigger(triggerlist_Mu8,prescale), PassTrigger(triggerlist_Mu17,prescale));
	weight *= prescale_trigger;	


  snu::KParticle muon;
  muon = muonLooseColl.at(0);
  muon.SetPxPyPzE(muon.Px(),muon.Py(),muon.Pz(),muon.E());

  float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;

  float etaarray [] = {0.0,0.8,1.479,2.0,2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  for(unsigned int i=0; i<jetColl.size(); i++){
    dR=jetColl.at(i).DeltaR(muon);
    dPhi=fabs(jetColl.at(i).DeltaPhi(muon));
    ptmu_ptjet = muon.Pt() / jetColl.at(i).Pt();
    FillHist("dR", dR, weight, 0., 10., 100);

    if(dR>1){
      FillHist("ptmu_ptjet", ptmu_ptjet, weight, 0., 500., 500);
      FillHist("dPhi", dPhi, weight, 0., 10., 100);
      if(ptmu_ptjet<1){
        if(dPhi>2.5){
          FillHist("eta_F0", fabs(muon.Eta()), weight, 0., 10., 100);
          FillHist("pt_F0", muon.Pt(), weight, 0., 500., 500);
					FillHist("events_F0", fabs(muon.Eta()), muon.Pt(), weight, etaarray, 4, ptarray, 9);
          if(muonTightColl.size() == 1){
            FillHist("eta_F", fabs(muon.Eta()), weight, 0., 10., 100);
            FillHist("pt_F", muon.Pt(), weight, 0., 500., 500);
						FillHist("events_F", fabs(muon.Eta()), muon.Pt(), weight, etaarray, 4, ptarray, 9);					
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
  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  DeclareVariable(out_muons, "Signal_Muons");
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

////////////////////////////////////////////

float FakeRateCalculator_Mu::GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh){
  float prescale_trigger= 0.;
  if(muon.size() ==1){
    if(muon.at(0).Pt() >= 20.){
      if(passhigh){
        prescale_trigger = (16.95) / 19789 ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = (3.546650) / 19789 ;
      }
      else prescale_trigger = 0.;
    }
  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}
