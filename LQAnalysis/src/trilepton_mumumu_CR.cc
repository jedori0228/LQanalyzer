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
  
  //==== Signal muons
  std::vector<snu::KMuon> muontriTightColl_raw;
  eventbase->GetMuonSel()->HNtriTightMuonSelection(muontriTightColl_raw);
  std::vector<snu::KMuon> muontriLooseColl_raw;
  eventbase->GetMuonSel()->HNtriLooseMuonSelection(muontriLooseColl_raw);
  //==== Large dXY muons
  std::vector<snu::KMuon> muontriHighdXYTightColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYTightMuonSelection(muontriHighdXYTightColl_raw);
  std::vector<snu::KMuon> muontriHighdXYLooseColl_raw;
  eventbase->GetMuonSel()->HNtriHighdXYLooseMuonSelection(muontriHighdXYLooseColl_raw);

  std::vector<snu::KMuon> muontriTightColl;
  std::vector<snu::KMuon> muontriLooseColl;
  std::vector<snu::KMuon> muontriHighdXYTightColl;
  std::vector<snu::KMuon> muontriHighdXYLooseColl;

  for(unsigned int i=0; i<muontriTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriTightColl_raw.at(i);
    FillHist("TriTightMuon_PartType", thismuon.GetType(), 1, 0, 15, 15);
    if(isData) muontriTightColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriLooseColl_raw.at(i);
    FillHist("TriLooseMuon_PartType", thismuon.GetType(), 1, 0, 15, 15);
    if(isData) muontriLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriLooseColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYTightColl_raw.at(i);
    FillHist("TriHighdXYTightMuon_PartType", thismuon.GetType(), 1, 0, 15, 15);
    if(isData) muontriHighdXYTightColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriHighdXYTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriHighdXYLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriHighdXYLooseColl_raw.at(i);
    FillHist("TriHighdXYLooseMuon_PartType", thismuon.GetType(), 1, 0, 15, 15);
    if(isData) muontriHighdXYLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetType() == 0 || thismuon.GetType() == 7 ) muontriHighdXYLooseColl.push_back(thismuon);
    }
  }
  
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
  int n_triHighdXYTight_muons = muontriHighdXYTightColl.size();
  int n_triHighdXYLoose_muons = muontriHighdXYLooseColl.size();
  int n_jets = jetColl_lepveto.size();

  // CR related variables //
  FillHist("n_tight_muons_control_PU", n_triTight_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_loose_muons_control_PU", n_triLoose_muons, weight*pileup_reweight, 0, 10, 10);
  if( n_triTight_muons == 2 && n_triLoose_muons == 2){
    int isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ? 1 : 0;
    FillHist("2Muons_OS0_SS1_control_PU", isSS, weight*pileup_reweight, 0, 2, 2);
  }
  FillHist("n_jets_control_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if(jetColl_lepveto.at(j).CombinedSecVertexBtag() > 0.679) n_bjets++;
  }
  FillHist("n_bjets_control_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);

  //================
  //==== define CR
  //================

  bool isCR = true;

  int n_muons(0);

  bool isTwoMuon   = n_triLoose_muons == 2 && n_triTight_muons == 2;
  bool isThreeMuon = n_triLoose_muons == 3 && n_triTight_muons == 3;
  //bool isHighdXYTwoMuon   = n_triHighdXYLoose_muons == 2 && n_triHighdXYTight_muons == 2;
  bool isHighdXYThreeMuon = n_triHighdXYLoose_muons == 3 && n_triHighdXYTight_muons == 3;

  if(k_jskim_flag_2 == "CR1"){
    if( !isTwoMuon ) return;
    n_muons = 2;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();
    isCR = isTwoMuon && isSS && (n_jets == 0) && (m_dimuon > 15.);
  }
  if(k_jskim_flag_2 == "CR2"){
    if( !isTwoMuon ) return;
    n_muons = 2;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();
    if( isTwoMuon && isSS && (n_bjets > 0) ) FillHist("TT_mll_low_control_PU", m_dimuon , weight*pileup_reweight, 0, 20, 200);
    isCR = isTwoMuon && isSS && (n_bjets > 0) && (m_dimuon > 15.);
  }
  if(k_jskim_flag_2 == "CR3"){
    if( !isTwoMuon ) return;
    n_muons = 2;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();
    isCR = isTwoMuon && !isSS && (n_bjets > 0) && (m_dimuon > 15.);
  }
  if(k_jskim_flag_2 == "CR4"){
    if( !isThreeMuon ) return;
    n_muons = 3;
    FillHist("3Muon", 0, 1, 0, 1, 1); 
    snu::KMuon lep[3];
    for(unsigned int i=0; i<3; i++) lep[i] = muontriLooseColl.at(i);
    if( lep[0].Charge() == lep[1].Charge() && lep[1].Charge() == lep[2].Charge() ) return;
    FillHist("3Muon1OS", 0, 1, 0, 1, 1);
    int SS_indices[2] = {0, 1}, OS_index(2);
    if( lep[0].Charge() == lep[2].Charge() ){
      SS_indices[0] = 0;
      SS_indices[1] = 2;
      OS_index = 1;
    }
    else if( lep[1].Charge() == lep[2].Charge() ){
      SS_indices[0] = 1;
      SS_indices[1] = 2;
      OS_index = 0;
    }
    double m_OS[2] = { ( lep[SS_indices[0]]+lep[OS_index] ).M(), ( lep[SS_indices[1]]+lep[OS_index] ).M() };
    double m_Z = 91.1876;
    double m_Z_candidate = m_OS[0];
    int W_muons_candidate_index(1);
    if( fabs(m_OS[0]-m_Z) > fabs(m_OS[1]-m_Z) ){
      m_Z_candidate = m_OS[1];
      W_muons_candidate_index = 0;
    }
    FillHist("m_Z_candidate", m_Z_candidate, weight*pileup_reweight, 50., 150., 150);
    if( fabs(m_Z_candidate-m_Z) > 10. ) return;
    snu::KEvent Evt = eventbase->GetEvent();
    double MET = Evt.PFMET(), METphi = Evt.PFMETphi();
    snu::KParticle nu;
    nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
    if( MET < 20. ) return;
    if( lep[SS_indices[W_muons_candidate_index]].Pt() < 20. ) return;
    double mT_W = ( lep[SS_indices[W_muons_candidate_index]] + nu ).Mt();
    FillHist("TT_mT_W_control_PU", mT_W, weight*pileup_reweight, 0., 300., 300);
	  FillHist("TT_mT_W_1_control_PU", mT_W, 1, 0., 300., 300);
  }
  if(k_jskim_flag_2 == "CR5"){
    n_muons = 3;
    isCR = isHighdXYThreeMuon;
  }

  if(!isCR) return;

  if(k_jskim_flag_2 == "CR5"){
    if( muontriHighdXYLooseColl.at(0).Pt() < 20. ) return;
  }
  else{
    if( muontriLooseColl.at(0).Pt() < 20. ) return;
  }

  snu::KMuon lep[n_muons];
  double LeptonRelIso[n_muons]; 
  for(int i=0; i<n_muons; i++){
    if(k_jskim_flag_2 == "CR5"){
      lep[i] = muontriHighdXYLooseColl.at(i);
    }
    else{
      lep[i] = muontriLooseColl.at(i);
    }
    LeptonRelIso[i] = (lep[i].SumIsoCHDR03() + std::max(0.0, lep[i].SumIsoNHDR03() + lep[i].SumIsoPHDR03() - 0.5* lep[i].SumPUIsoR03()))/lep[i].Pt();
  }

  //==== FR weighted plots
  FillHist("TT_n_events_control_PU", 0, weight*pileup_reweight, 0, 1, 1);
  FillHist("TT_n_jets_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  FillHist("TT_n_bjets_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);
  FillHist("TT_mll_control_PU", (lep[0]+lep[1]).M() , weight*pileup_reweight, 0, 500, 500);
  FillHist("TT_leadingLepton_Pt_control_PU", lep[0].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("TT_leadingLepton_Eta_control_PU", lep[0].Eta() , weight*pileup_reweight, -3, 3, 60);
  FillHist("TT_leadingLepton_RelIso_control_PU", LeptonRelIso[0] , weight*pileup_reweight, 0, 1.0, 10);
  FillHist("TT_leadingLepton_Chi2_control_PU", lep[0].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
  FillHist("TT_secondLepton_Pt_control_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("TT_secondLepton_Eta_control_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
  FillHist("TT_secondLepton_RelIso_control_PU", LeptonRelIso[1] , weight*pileup_reweight, 0, 1.0, 10);
  FillHist("TT_secondLepton_Chi2_control_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);

  //==== weight = 1 plots
  FillHist("TT_n_events_1_control_PU", 0, 1., 0, 1, 1);
  FillHist("TT_n_jets_PU", n_jets, 1., 0, 10, 10);
  FillHist("TT_n_bjets_PU", n_bjets, 1., 0, 10, 10);
  FillHist("TT_mll_1_control_PU", (lep[0]+lep[1]).M(), 1., 0, 500, 500);
  FillHist("TT_leadingLepton_Pt_1_control_PU", lep[0].Pt(), 1., 0, 200, 200);
  FillHist("TT_leadingLepton_Eta_1_control_PU", lep[0].Eta(), 1., -3, 3, 60);
  FillHist("TT_leadingLepton_RelIso_1_control_PU", LeptonRelIso[0], 1., 0, 1.0, 10);
  FillHist("TT_leadingLepton_Chi2_1_control_PU", lep[0].GlobalChi2(), 1., 0, 10, 100);
  FillHist("TT_secondLepton_Pt_1_control_PU", lep[1].Pt(), 1., 0, 200, 200);
  FillHist("TT_secondLepton_Eta_1_control_PU", lep[1].Eta(), 1., -3, 3, 60);
  FillHist("TT_secondLepton_RelIso_1_control_PU", LeptonRelIso[1], 1., 0, 1.0, 10);
  FillHist("TT_secondLepton_Chi2_1_control_PU", lep[1].GlobalChi2(), 1., 0, 10, 100);

  if(n_muons==3){
    FillHist("TT_thirdLepton_Pt_control_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
    FillHist("TT_thirdLepton_Eta_control_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
    FillHist("TT_thirdLepton_RelIso_control_PU", LeptonRelIso[1] , weight*pileup_reweight, 0, 1.0, 10);
    FillHist("TT_thirdLepton_Chi2_control_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);

    FillHist("TT_thirdLepton_Pt_1_control_PU", lep[2].Pt(), 1., 0, 200, 200);
    FillHist("TT_thirdLepton_Eta_1_control_PU", lep[2].Eta(), 1., -3, 3, 60);
    FillHist("TT_thirdLepton_RelIso_1_control_PU", LeptonRelIso[2], 1., 0, 1.0, 10);
    FillHist("TT_thirdLepton_Chi2_1_control_PU", lep[2].GlobalChi2(), 1., 0, 10, 100);
  }

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
