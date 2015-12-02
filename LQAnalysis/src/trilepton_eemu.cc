// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzerDiElectron Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_eemu.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_eemu);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_eemu::trilepton_eemu() :  AnalyzerCore() {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_eemu");

  Message("In trilepton_eemu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
}


void trilepton_eemu::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_eemu::ExecuteEvents()throw( LQError ){
  
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
    //pileup_reweight = reweightPU->GetWeight(int(eventbase->GetEvent().PileUpInteractionsTrue()))* MCweight;
		pileup_reweight = MCweight; // for signals and Wtollln
  }

  
    
  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////   
  

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

	std::vector<snu::KTruth> truthColl;
	eventbase->GetTruthSel()->Selection(truthColl);

	double Event[100];
 	double Muon[100];
	double Electron[100];

  if(electronLooseColl.size()<2 || muonLooseColl.size()<1)	return;
  FillCutFlow("eemu", weight);
	if(electronLooseColl.at(0).Charge() != electronLooseColl.at(1).Charge() || electronLooseColl.at(0).Charge() == muonLooseColl.at(0).Charge()) return;
  FillCutFlow("same_ee_opp_mu", weight);
  
  
  snu::KEvent Evt = eventbase->GetEvent();
  double METv = Evt.PFMET();
  double METphi = Evt.PFMETphi();
  Event[0] = Evt.PFMET();
  Event[1] = Evt.PFMETphi();
  Event[2] = jetColl_lepveto.size();
 	Event[3] = weight;
  Event[4] = pileup_reweight;
  Event[5] = MCweight; 

  Muon[0] = muonLooseColl.size();
  snu::KMuon mu;
  double LeptonRelIso;
	int n_muonvar = 18;
  for(int i=0; i<1; i++){
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
  }
  
  
  snu::KElectron el;
  double PHONH_03[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14};
  double rho = eventbase->GetEvent().JetRho();
  double ElectronIsoDR03, LeptonRelIsoDR03;
  int ifid, n_elvar = 11;
	float id_scalefactor;

  Electron[0] = electronLooseColl.size();
  for(int i=0; i<2; i++){
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
        if(!isData) id_scalefactor *=  ElectronScaleFactor(el.Eta(), el.Pt());
        Electron[11+i*n_elvar] = id_scalefactor;
  }

  
  FillNtp("Event", Event);
  FillNtp("Muon", Muon);
  FillNtp("Electron", Electron);
  
  return;
}// End of execute event loop


void trilepton_eemu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_eemu::BeginCycle() throw( LQError ){
  
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

trilepton_eemu::~trilepton_eemu() {
  
  Message("In trilepton_eemu Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_eemu::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 6,0.,6.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"eemu");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"same_ee_opp_mu");
   
    
  }
}


void trilepton_eemu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void trilepton_eemu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);

  MakeNtp("Event", "MET:METphi:n_jets:weight:PU_reweight:MCweight");
  MakeNtp("Muon", "n_LooseMuon:0_Px:0_Py:0_Pz:0_E:0_LeptonRelIso:0_Eta:0_dXY:0_dZ:0_GlobalChi2:0_IsoHcalVeto:0_IsoEcalVeto:0_IsGlobal:0_validHits:0_validPixHits:0_validStations:0_ActiveLayer:0_IsPf:0_IsTracker");
  MakeNtp("Electron", "n_LooseElectron:0_HasMatchedConvPhot:0_MissingHits:0_LeptonRelIsoDR03:0_GsfCtfScPixChargeConsistency:0_Eta:0_Px:0_Py:0_Pz:0_E:0_dxy:0_id_scalefactor:1_HasMatchedConvPhot:1_MissingHits:1_LeptonRelIsoDR03:1_GsfCtfScPixChargeConsistency:1_Eta:1_Px:1_Py:1_Pz:1_E:1_dxy:1_id_scalefactor");

  /**
   *  Remove//Overide this trilepton_eemuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_eemu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //

}



