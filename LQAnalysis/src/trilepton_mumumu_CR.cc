// $Id: trilepton_mumumu_CR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_CR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond     <jalmond@cern.ch>        - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_CR.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                  
#include "BaseSelection.h"

#define FR_n_pt_bin 7
#define FR_n_eta_bin 4

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_CR);


/**
 *  This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_mumumu_CR::trilepton_mumumu_CR() :  AnalyzerCore(), out_muons(0)
{

  
  // To have the correct name in the log:                                                                                   
  SetLogName("trilepton_mumumu_CR");

  Message("In trilepton_mumumu_CR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  
}


void trilepton_mumumu_CR::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)  output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  return;
 }


void trilepton_mumumu_CR::ExecuteEvents()throw( LQError ){

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  /// FillCutFlow(cut, weight) fills a basic TH1 called cutflow. It is used to check number of events passing different cuts
  /// The string cut must match a bin label in FillCutFlow function
  FillCutFlow("NoCut", weight);
  
  //cout << "JES UNC GetJECUncertainty(FlavorPureQuark, 2.2, 100.,true) = " << GetJECUncertainty("FlavorPureQuark", 2.2, 100.,true) << endl;
  //cout << "JES UNC GetJECUncertainty(FlavorPureQuark, 5.2, 1000.,true) = " << GetJECUncertainty("FlavorPureQuark", 5.2, 1000.,true) << endl;
  //cout << "JES UNC GetJECUncertainty(FlavorPureGluon, 0.1, 500.,true) = " << GetJECUncertainty("FlavorPureGluon", 0.1, 500.,true) << endl;

  //return;


  ///// Apply some general cuts on event to clean MET
  /// Taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
  /// These are applied in AnalyzerCore::PassBasicEventCuts
  if(!PassBasicEventCuts()) return;    /// Initial event cuts  
  FillCutFlow("EventCut", weight);

  
  /// Trigger List (specific to muons channel)
  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TkMu8_v");
    
  if(!PassTrigger(triggerslist, prescale)) return;
  
  //// if the trigger that fired the event is prescaled you can reweight the event accordingly using the variable prescale
  
  FillCutFlow("TriggerCut", weight);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  
  
  /// Check the event has a "Good" Primary vertex
  /// Good is taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingPFGJob:
  /// defined as : !isFake && ndof > 4 && |z| <= 24 cm && position.Rho <= 2cm (rho = radius of vertex)
  /// Cut is coded in SKTreeFiller and stored in KEvent class as HasGoodPrimaryVertex()
  /// More info on primary vertex can be found https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction (LQNtuples use offlinePrimaryVertices)
  // isFake is true if the vertex is based on the beam spot (as no reconstructed vertex is found

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  
  FillCutFlow("VertexCut", weight);
  
  /// Use the number of vertices in the event to check effect of pileup reweighting
  numberVertices = eventbase->GetEvent().nVertices();  
  
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    /// Here is an alternative method to Fill a histogram. 
    /// The histogram with name "h_nvtx_norw"/"h_nvtx_rw" were not declared in the MakeHistogram code. 
    /// To avoid adding this by hand we can just use FillHist() function with 3 additional inputs i.e., xmin, xmax and nbinsx       
    pileup_reweight = reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
  }
  
  //std::vector<snu::KMuon> muonTightColl;
  //eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl);
  //std::vector<snu::KMuon> muonHighPtColl;
  //eventbase->GetMuonSel()->HNTightHighPtMuonSelection(muonHighPtColl);

/*
  for(std::vector<snu::KMuon>::iterator it = muonTightColl.begin(); it!= muonTightColl.end(); it++){
    cout << "Weight = " << weight << endl;
    weight *= MuonScaleFactor(it->Eta(), it->Pt());
    cout << "Weight = " << weight << endl;
    cout << "Tight muon pt = " << it->Pt() << " " << it->Eta() << " " << it->Phi() << endl; 
  }
*/
  
  /// Correct the muon momentum with rochester corrections
  //CorrectMuonMomentum(muonTightColl);
  //CorrectMuonMomentum(muonHighPtColl);
  
  /// Example of how to get fake weight for dimuon channel
  //std::vector<snu::KMuon> muonLooseColl;
  //eventbase->GetMuonSel()->HNLooseMuonSelection(muonLooseColl);
  
  std::vector<snu::KMuon> muontriTightColl;
  eventbase->GetMuonSel()->HNtriTightMuonSelection(muontriTightColl);
  std::vector<snu::KMuon> muontriLooseColl;
  eventbase->GetMuonSel()->HNtriLooseMuonSelection(muontriLooseColl);
  
  std::vector<snu::KElectron> electronTightColl;
  eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);

  // (jit->Pt() >= 20.) && fabs(jit->Eta()) < 2.5   && PassUserID(PFJET_LOOSE, *jit) && jit->PileupJetIDLoose() //
  std::vector<snu::KJet> jetColl_lepveto;
  //eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muontriTightColl, electronTightColl); // HNSelection is too tight
  eventbase->GetJetSel()->SetEta(2.4);
  //eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl_lepveto, muontriTightColl, electronTightColl);
  eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl_lepveto, AnalyzerCore::GetMuons("veto"), AnalyzerCore::GetElectrons(false,false, "veto") );

  int n_triTight_muons = muontriTightColl.size();
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_lepveto.size();

  // CR related variables //
  FillHist("n_tight_muons_control_PU", n_triTight_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_loose_muons_control_PU", n_triLoose_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_jets_control_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  int n_bjets=0;
  for(UInt_t j=0; j < n_jets; j++){
    if(jetColl_lepveto.at(j).CombinedSecVertexBtag() > 0.679) n_bjets++;
  }
  FillHist("n_bjets_control_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);

  //================
  //==== define CR
  //================

  //==== 1) SS dimuon
  if( ! (n_triLoose_muons == 2) ) return;
  if( ! (n_triTight_muons == 2) ) return;
  if( ! (muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge()) ) return;
  //==== 2) jets...
  //if( ! n_jets == 0 ) return;
  if( ! n_bjets > 0) return;

  snu::KMuon lep[2];
  double LeptonRelIso[2]; 
  for(unsigned int i=0; i<2; i++){
    lep[i] = muontriLooseColl.at(i);
    LeptonRelIso[i] = (lep[i].SumIsoCHDR03() + std::max(0.0, lep[i].SumIsoNHDR03() + lep[i].SumIsoPHDR03() - 0.5* lep[i].SumPUIsoR03()))/lep[i].Pt();
  }

  //==== FR weighted plots
  FillHist("TT_n_events_control_PU", 0, weight*pileup_reweight, 0, 1, 1);
  FillHist("TT_mll_control_PU", (lep[0]+lep[1]).M() , weight*pileup_reweight, 0, 200, 200);
  FillHist("TT_leadingLepton_Pt_control_PU", lep[0].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("TT_secondLepton_Pt_control_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("TT_leadingLepton_Eta_control_PU", lep[0].Eta() , weight*pileup_reweight, -3, 3, 60);
  FillHist("TT_secondLepton_Eta_control_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
  FillHist("TT_leadingLepton_RelIso_control_PU", LeptonRelIso[0] , weight*pileup_reweight, 0, 1.0, 10);
  FillHist("TT_secondLepton_RelIso_control_PU", LeptonRelIso[1] , weight*pileup_reweight, 0, 1.0, 10);
  FillHist("TT_leadingLepton_Chi2_control_PU", lep[0].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
  FillHist("TT_secondLepton_Chi2_control_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
  //==== weight = 1 plots
  FillHist("TT_n_events_1_control_PU", 0, 1, 0, 1, 1);
  FillHist("TT_mll_1_control_PU", (lep[0]+lep[1]).M() , 1, 0, 200, 200);
  FillHist("TT_leadingLepton_Pt_1_control_PU", lep[0].Pt() , 1, 0, 200, 200);
  FillHist("TT_secondLepton_Pt_1_control_PU", lep[1].Pt() , 1, 0, 200, 200);
  FillHist("TT_leadingLepton_Eta_1_control_PU", lep[0].Eta() , 1, -3, 3, 60);
  FillHist("TT_secondLepton_Eta_1_control_PU", lep[1].Eta() , 1, -3, 3, 60);
  FillHist("TT_leadingLepton_RelIso_1_control_PU", LeptonRelIso[0] , 1, 0, 1.0, 10);
  FillHist("TT_secondLepton_RelIso_1_control_PU", LeptonRelIso[1] , 1, 0, 1.0, 10);
  FillHist("TT_leadingLepton_Chi2_1_control_PU", lep[0].GlobalChi2() , 1, 0, 10, 100);
  FillHist("TT_secondLepton_Chi2_1_control_PU", lep[1].GlobalChi2() , 1, 0, 10, 100);

  return;

}// End of execute event loop
  


void trilepton_mumumu_CR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_CR::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "MyDataPileupHistogram.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

trilepton_mumumu_CR::~trilepton_mumumu_CR() {
  
  Message("In trilepton_mumumu_CR Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
}


void trilepton_mumumu_CR::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
   GetHist("cutflow")->Fill(cut,weight);
  
  }
  else{
   AnalyzerCore::MakeHistograms("cutflow", 7,0.,7.);

   GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3muon");
   
  }
}


void trilepton_mumumu_CR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_CR::MakeHistograms(){
  //// Additional plots to make
   
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
  *  Remove//Overide this trilepton_mumumu_CRCore::MakeHistograms() to make new hists for your analysis
  **/
  
}


void trilepton_mumumu_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
