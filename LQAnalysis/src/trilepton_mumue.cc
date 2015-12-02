// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzerDiElectron Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumue.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumue);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_mumue::trilepton_mumue() :  AnalyzerCore() {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumue");

  Message("In trilepton_mumue constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
}


void trilepton_mumue::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist
   MakeCleverHistograms(sighist, "ElectronMuonEvents");

   
   return;
}


void trilepton_mumue::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  
  /// FillCutFlow(cut, weight) fills a basic TH1 called cutflow. It is used to check number of events passing different cuts
  /// The string cut must match a bin label in FillCutFlow function
  FillCutFlow("NoCut", weight);
  
  ///// Apply some general cuts on event to clean MET
  /// Taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
  /// These are applied in AnalyzerCore::PassBasicEventCuts
  if(!PassBasicEventCuts()) return;     /// Initial event cuts  
  FillCutFlow("EventCut", weight);
  
  /// Trigger List 
  std::vector<TString> triggerslist;
  /// This is the analysis electron trigger 
  /// No Scale Factors are yet applied to correct MC
  triggerslist.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  triggerslist.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoV");
  //if(!PassTrigger(triggerslist, prescale)) return;
  
  //// if the trigger that fired the event is prescaled you can reweight the event accordingly using the variable prescale
  
  FillCutFlow("TriggerCut", weight);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
    
  FillCutFlow("VertexCut", weight);
  
  /// Use the number of vertices in the event to check effect of pileup reweighting
  numberVertices = eventbase->GetEvent().nVertices();   
    /// Correct MC for pileup   
  
  float pileup_reweight (1.);
  if (MC_pu&&!k_isdata) {
    pileup_reweight = reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
		//pileup_reweight = MCweight; // for Wtollln
  }
  
    
  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////   
  
  /// We want to select events with 2 medium electrons (we will also remove events with a looser third muon to show how it is done)
  /// We will use 4 different object collections
  /// 1) Tight Electrons  ||  eventbase->GetElectronSel()->HNTightElectronSelection
  /// 2) Tight Muons      || eventbase->GetMuonSel()->HNTightMuonSelection
  /// 3) Jets(with lepton veto)
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  /// 1) TightElectrons                                                                                                                                                     
  ///////////////////////////////////////////////////////////////////////////////////////////

  std::vector<snu::KElectron> electronTightColl;
  eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);
  std::vector<snu::KElectron> electronLooseColl;
  eventbase->GetElectronSel()->HNLooseElectronSelection(electronLooseColl);
	std::vector<snu::KElectron> electronVetoColl; 
	eventbase->GetElectronSel()->HNVetoElectronSelection(electronVetoColl);

  std::vector<snu::KMuon> muonTightColl;
  eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl);
  std::vector<snu::KMuon> muonLooseColl;
  eventbase->GetMuonSel()->HNLooseMuonSelection(muonLooseColl);
 
  std::vector<snu::KJet> jetColl_lepveto;
  eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
  eventbase->GetJetSel()->SetPt(20.);
  eventbase->GetJetSel()->SetEta(2.5);
  eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonTightColl, electronTightColl); 
  //eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonLooseColl, electronLooseColl);

	std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

	double Event[100];
 	double Muon[100];
	double Electron[100];
  double Jet[100];
	//double Gen[100];

  int n_loose_muon = muonLooseColl.size(), n_loose_electron = electronLooseColl.size(), n_jet = jetColl_lepveto.size();
  if(n_loose_electron < 1 || n_loose_muon < 2) return;
  FillCutFlow("mumue", weight);
	//if(muonLooseColl.at(0).Charge() != muonLooseColl.at(1).Charge() || electronLooseColl.at(0).Charge() == muonLooseColl.at(0).Charge()) return;
  //FillCutFlow("same_mumu_opp_e", weight);
 
  int n_bjet = 0;
  for(int i=0; i<n_jet; i++){
    if(jetColl_lepveto.at(i).CombinedSecVertexBtag() > 0.679){
      n_bjet++;
    }
  }
 
  snu::KEvent Evt = eventbase->GetEvent();
  double METv = Evt.PFMET();
  double METphi = Evt.PFMETphi();
  Event[0] = Evt.PFMET();
  Event[1] = Evt.PFMETphi();
  Event[2] = n_bjet;
	Event[3] = weight; 
  Event[4] = pileup_reweight;
  Event[5] = MCweight;
 
  Muon[0] = n_loose_muon; // <- only n_loose_muon <= 3 will be saved in the ntuple
  snu::KMuon mu;
  double LeptonRelIso;
	int n_muonvar = 20;
  for(int i=0; i<n_loose_muon; i++){
    mu = muonLooseColl.at(i);
    LeptonRelIso = (mu.SumIsoCHDR03() + std::max(0.0, mu.SumIsoNHDR03() + mu.SumIsoPHDR03() - 0.5*mu.SumPUIsoR03()))/mu.Pt();
    Muon[1+n_muonvar*i] = mu.Px();
    Muon[2+n_muonvar*i] = mu.Py();
    Muon[3+n_muonvar*i] = mu.Pz();
    Muon[4+n_muonvar*i] = mu.E();
    Muon[5+n_muonvar*i] = LeptonRelIso;
    Muon[6+n_muonvar*i] = mu.Eta();
    Muon[7+n_muonvar*i] = mu.dXY();
    Muon[8+n_muonvar*i] = mu.dZ();
    Muon[9+n_muonvar*i] = mu.GlobalChi2();
    Muon[10+n_muonvar*i] = mu.IsoHcalVeto();
    Muon[11+n_muonvar*i] = mu.IsoEcalVeto();
    Muon[12+n_muonvar*i] = mu.IsGlobal();
    Muon[13+n_muonvar*i] = mu.validHits();
    Muon[14+n_muonvar*i] = mu.validPixHits();
    Muon[15+n_muonvar*i] = mu.validStations();
    Muon[16+n_muonvar*i] = mu.ActiveLayer();
    Muon[17+n_muonvar*i] = mu.IsPF();
    Muon[18+n_muonvar*i] = mu.IsTracker();
    Muon[19+n_muonvar*i] = mu.Charge();
    Muon[20+n_muonvar*i] = mu.Pt();
	}
  
  
  snu::KElectron el;
  double PHONH_03[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14};
  double rho = eventbase->GetEvent().JetRho();
  double ElectronIsoDR03, LeptonRelIsoDR03;
  int ifid, n_elvar = 13;
  float id_scalefactor;

  Electron[0] = n_loose_electron; // <- only n_loose_electron <= 2 will be saved in the ntuple
  for(int i=0; i<n_loose_electron; i++){
    el = electronLooseColl.at(i);
    if (fabs(el.SCEta()) < 1.0) ifid = 0;
    else if (fabs(el.SCEta()) < 1.479) ifid = 1;
    else if (fabs(el.SCEta()) < 2.0) ifid = 2;
    else if (fabs(el.SCEta()) < 2.2) ifid = 3;
    else if (fabs(el.SCEta()) < 2.3) ifid = 4;
    else if (fabs(el.SCEta()) < 2.4) ifid = 5;
    else ifid = 6;
    ElectronIsoDR03 =  el.PFChargedHadronIso03() + max( el.PFNeutralHadronIso03() + el.PFPhotonIso03() - rho * PHONH_03[ifid],  0.);
    LeptonRelIsoDR03 = ElectronIsoDR03/  el.Pt();
    Electron[1+i*n_elvar] = el.HasMatchedConvPhot();
    Electron[2+i*n_elvar] = el.MissingHits();
    Electron[3+i*n_elvar] = LeptonRelIsoDR03;
    Electron[4+i*n_elvar] = el.GsfCtfScPixChargeConsistency();
    Electron[5+i*n_elvar] = el.Eta();
    Electron[6+i*n_elvar] = el.Px();
    Electron[7+i*n_elvar] = el.Py();
    Electron[8+i*n_elvar] = el.Pz();
    Electron[9+i*n_elvar] = el.E();
    Electron[10+i*n_elvar] = el.dxy();
    id_scalefactor = 1.0;
    if(!isData)	id_scalefactor *=  ElectronScaleFactor(el.Eta(), el.Pt());
    Electron[11+i*n_elvar] = id_scalefactor;
    Electron[12+i*n_elvar] = el.Charge();
    Electron[13+i*n_elvar] = el.Pt();
  }

  snu::KJet jet;
  int n_jetvar = 9;
  Jet[0] = n_jet; // <- only n_loose_electron <= 4 will be saved in the ntuple
  for(int i=0; i<n_jet; i++){
    jet = jetColl_lepveto.at(i);
    Jet[1+i*n_jetvar] = jet.Px();
    Jet[2+i*n_jetvar] = jet.Py();
    Jet[3+i*n_jetvar] = jet.Pz();
    Jet[4+i*n_jetvar] = jet.E();
    Jet[5+i*n_jetvar] = jet.Eta();
    Jet[6+i*n_jetvar] = jet.Pt();
    Jet[7+i*n_jetvar] = jet.PFJetTrackCountingHighPurBTag();
    Jet[8+i*n_jetvar] = jet.BtagProb();
    Jet[9+i*n_jetvar] = jet.CombinedSecVertexBtag();
  }
  
  
  FillNtp("Event", Event);
  FillNtp("Muon", Muon);
  FillNtp("Electron", Electron);
	FillNtp("Jet", Jet);

 /* 
   ////////////////////////
   /// Generation Level ///
   ////////////////////////
   
   int truthmax = truthColl.size();
   int gen_l_1_index, gen_W_real_index, gen_HN_index, gen_W_virtual_index, gen_l_3_index;
   snu::KParticle gen_nu, gen_W_real, gen_HN, gen_W_virtual;
   snu::KParticle gen_l_1, gen_l_2, gen_l_3;
   bool isWV=false;
   
   
   
   for(int i=2;i<truthmax;i++){ // find W_real index
   if(abs(truthColl.at(i).PdgId()) == 24){
   gen_W_real_index = i;
   gen_W_real = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){ // find HN index
   if(truthColl.at(i).PdgId() == 80000002){
   gen_HN_index = i;
   gen_HN = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){ // check W_virtual exists
   if(fabs(truthColl.at(i).PdgId()) == 24 && truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 80000002){
   gen_W_virtual_index = i;
   gen_W_virtual = truthColl.at(i);
   isWV = true;
   break;
   }
   }
   
   if(!isWV){ // no W_virtual in truthColl index case
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){ // find l_1 at gen. level (muon)
   gen_l_1 = truthColl.at(i);
   gen_l_1_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_HN_index){ // find nu at gen. level (muon)
   gen_nu = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_index).PdgId() && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_3 at gen. level (muon)
   gen_l_3 = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 11 && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_2 at gen. level (electron)
   gen_l_2 = truthColl.at(i);
   break;
   }
   }
   }
   else{ // W_virtual in truthColl index case
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){ // find l_1 at gen. level (muon)
   gen_l_1 = truthColl.at(i);
   gen_l_1_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_virtual_index){ // find nu at gen. level (muon)
   gen_nu = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(fabs(truthColl.at(i).PdgId()) == fabs(truthColl.at(gen_l_1_index).PdgId()) && truthColl.at(i).IndexMother() == gen_W_virtual_index){ // find l_3 at gen. level (muon)
   gen_l_3 = truthColl.at(i);
   gen_l_3_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(fabs(truthColl.at(i).PdgId()) == 11 && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_2 at gen. level (electron)
   gen_l_2 = truthColl.at(i);
   }
   }
   if(truthColl.at(gen_l_3_index).PdgId() == -truthColl.at(gen_l_1_index).PdgId()){
   snu::KParticle gen_temp = gen_l_2;
   gen_l_2 = gen_l_3;
   gen_l_3 = gen_temp;
   }
   }
   
   Gen[0] = gen_l_1.Px();
   Gen[1] = gen_l_1.Py();
   Gen[2] = gen_l_1.Pz();
   Gen[3] = gen_l_1.E();
   Gen[4] = gen_l_2.Px();
   Gen[5] = gen_l_2.Py();
   Gen[6] = gen_l_2.Pz();
   Gen[7] = gen_l_2.E();
   Gen[8] = gen_l_3.Px();
   Gen[9] = gen_l_3.Py();
   Gen[10] = gen_l_3.Pz();
   Gen[11] = gen_l_3.E();
   Gen[12] = gen_nu.Px();
   Gen[13] = gen_nu.Py();
   Gen[14] = gen_nu.Pz();
   Gen[15] = gen_nu.E();
   Gen[16] = gen_W_real.Px();
   Gen[17] = gen_W_real.Py();
   Gen[18] = gen_W_real.Pz();
   Gen[19] = gen_W_real.E();
   Gen[20] = gen_HN.Px();
   Gen[21] = gen_HN.Py();
   Gen[22] = gen_HN.Pz();
   Gen[23] = gen_HN.E();

  
  FillNtp("Gen", Gen);
*/  
  return;
}// End of execute event loop


void trilepton_mumue::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumue::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "MyDataPileupHistogram.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree


  
  return;
  
}

trilepton_mumue::~trilepton_mumue() {
  
  Message("In trilepton_mumue Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_mumue::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 5,0.,5.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"mumue");
    
  }
}


void trilepton_mumue::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void trilepton_mumue::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);

	MakeNtp("Event", "MET:METphi:n_bjet:weight:PU_reweight:MCweight");
  MakeNtp("Muon", "n_LooseMuon:0_Px:0_Py:0_Pz:0_E:0_LeptonRelIso:0_Eta:0_dXY:0_dZ:0_GlobalChi2:0_IsoHcalVeto:0_IsoEcalVeto:0_IsGlobal:0_validHits:0_validPixHits:0_validStations:0_ActiveLayer:0_IsPf:0_IsTracker:0_Charge:0_Pt:1_Px:1_Py:1_Pz:1_E:1_LeptonRelIso:1_Eta:1_dXY:1_dZ:1_GlobalChi2:1_IsoHcalVeto:1_IsoEcalVeto:1_IsGlobal:1_validHits:1_validPixHits:1_validStations:1_ActiveLayer:1_IsPf:1_IsTracker:1_Charge:1_Pt:2_Px:2_Py:2_Pz:2_E:2_LeptonRelIso:2_Eta:2_dXY:2_dZ:2_GlobalChi2:2_IsoHcalVeto:2_IsoEcalVeto:2_IsGlobal:2_validHits:2_validPixHits:2_validStations:2_ActiveLayer:2_IsPf:2_IsTracker:2_Charge:2_Pt");
  MakeNtp("Electron", "n_LooseElectron:0_HasMatchedConvPhot:0_MissingHits:0_LeptonRelIsoDR03:0_GsfCtfScPixChargeConsistency:0_Eta:0_Px:0_Py:0_Pz:0_E:0_dxy:0_id_scalefactor:0_Charge:0_Pt:1_HasMatchedConvPhot:1_MissingHits:1_LeptonRelIsoDR03:1_GsfCtfScPixChargeConsistency:1_Eta:1_Px:1_Py:1_Pz:1_E:1_dxy:1_id_scalefactor:1_Charge:1_Pt");
  MakeNtp("Jet","n_Jet:0_Px:0_Py:0_Pz:0_E:0_Eta:0_Pt:0_PFJetTrackCountingHighPurBTag:0_BtagProb:0_CombinedSecVertexBtag:1_Px:1_Py:1_Pz:1_E:1_Eta:1_Pt:1_PFJetTrackCountingHighPurBTag:1_BtagProb:1_CombinedSecVertexBtag:2_Px:2_Py:2_Pz:2_E:2_Eta:2_Pt:2_PFJetTrackCountingHighPurBTag:2_BtagProb:2_CombinedSecVertexBtag:3_Px:3_Py:3_Pz:3_E:3_Eta:3_Pt:3_PFJetTrackCountingHighPurBTag:3_BtagProb:3_CombinedSecVertexBtag:4_Px:4_Py:4_Pz:4_E:4_Eta:4_Pt:4_PFJetTrackCountingHighPurBTag:4_BtagProb:4_CombinedSecVertexBtag:5_Px:5_Py:5_Pz:5_E:5_Eta:5_Pt:5_PFJetTrackCountingHighPurBTag:5_BtagProb:5_CombinedSecVertexBtag:6_Px:6_Py:6_Pz:6_E:6_Eta:6_Pt:6_PFJetTrackCountingHighPurBTag:6_BtagProb:6_CombinedSecVertexBtag:7_Px:7_Py:7_Pz:7_E:7_Eta:7_Pt:7_PFJetTrackCountingHighPurBTag:7_BtagProb:7_CombinedSecVertexBtag:8_Px:8_Py:8_Pz:8_E:8_Eta:8_Pt:8_PFJetTrackCountingHighPurBTag:8_BtagProb:8_CombinedSecVertexBtag:9_Px:9_Py:9_Pz:9_E:9_Eta:9_Pt:9_PFJetTrackCountingHighPurBTag:9_BtagProb:9_CombinedSecVertexBtag");
  /**
   *  Remove//Overide this trilepton_mumueCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_mumue::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //

}



