// $Id: ExampleAnalyzerDiMuon.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzerDiMuon Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond		 <jalmond@cern.ch>			  - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzerDiMuon.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"																																									
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ExampleAnalyzerDiMuon);


/**
 *	This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
ExampleAnalyzerDiMuon::ExampleAnalyzerDiMuon() :  AnalyzerCore(), out_muons(0)  {

  
  // To have the correct name in the log:																																									 
  SetLogName("ExampleAnalyzerDiMuon");

  Message("In ExampleAnalyzerDiMuon constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
}


void ExampleAnalyzerDiMuon::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)	output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

	Message("Making clever hists for Z ->ll test code", INFO);

	return;
 }


void ExampleAnalyzerDiMuon::ExecuteEvents()throw( LQError ){

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
	
	/// FillCutFlow(cut, weight) fills a basic TH1 called cutflow. It is used to check number of events passing different cuts
	/// The string cut must match a bin label in FillCutFlow function
	FillCutFlow("NoCut", weight);
	
	///// Apply some general cuts on event to clean MET
	/// Taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
	/// These are applied in AnalyzerCore::PassBasicEventCuts
	if(!PassBasicEventCuts()) return;	  /// Initial event cuts  
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
	  pileup_reweight = reweightPU->GetWeight(eventbase->GetEvent().PileUpInteractionsTrue())* MCweight;
	}
	
	std::vector<snu::KMuon> muonTightColl;
	eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl);
	
	std::vector<snu::KMuon> muonHighPtColl;
	eventbase->GetMuonSel()->HNTightHighPtMuonSelection(muonHighPtColl);
	
	CorrectMuonMomentum(muonTightColl);
	CorrectMuonMomentum(muonHighPtColl);
	
	std::vector<snu::KMuon> muonVetoColl;
	eventbase->GetMuonSel()->HNVetoMuonSelection(muonVetoColl);


	std::vector<snu::KElectron> electronTightColl;
	eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);

	std::vector<snu::KJet> jetColl_lepveto;
	eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonTightColl, electronTightColl);

	std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl); 

/*
	cout
	<< "======================================================================================================" << endl
	<< "Index" << "\t" << "PdgId" << "\t" << "Status" << "\t" <<  "Mother" << "\t" << "Mother PdgId" << endl;
	for(int i=2; i<truthColl.size(); i++){
		cout << i << "\t" << truthColl.at(i).PdgId() << "\t" << truthColl.at(i).GenStatus() << "\t" << truthColl.at(i).IndexMother() << "\t" << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << "\t" << endl; 
	}
*/

	int Mother_of_HN, Mother_of_mu_1;
	for(int i=0; i<truthColl.size(); i++){
		if(truthColl.at(i).PdgId() == 80000002){
			Mother_of_HN = truthColl.at(i).IndexMother();
			break;
		}
	}
	for(int i=0; i<truthColl.size(); i++){
		if( abs(truthColl.at(i).PdgId()) == 13 ){
			Mother_of_mu_1 = truthColl.at(i).IndexMother();
			break;
		}
	}

	if(Mother_of_HN != Mother_of_mu_1){
  cout
  << "======================================================================================================" << endl
  << "Index" << "\t" << "PdgId" << "\t" << "Status" << "\t" <<  "Mother" << "\t" << "Mother PdgId" << endl;
		cout
		<< "Mother_of_HN = " << Mother_of_HN << endl
		<< "Mother of mu = " << Mother_of_mu_1 << endl;
		for(int i=2; i<truthColl.size(); i++){
			cout << i << "\t" << truthColl.at(i).PdgId() << "\t" << truthColl.at(i).GenStatus() << "\t" << truthColl.at(i).IndexMother() << "\t" << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << "\t" << endl;
		}
	} 
  
  return;
}// End of execute event loop
  


void ExampleAnalyzerDiMuon::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
}


void ExampleAnalyzerDiMuon::BeginCycle() throw( LQError ){
  
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

ExampleAnalyzerDiMuon::~ExampleAnalyzerDiMuon() {
  
  Message("In ExampleAnalyzerDiMuon Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void ExampleAnalyzerDiMuon::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
	 GetHist("cutflow")->Fill(cut,weight);
	
  }
  else{
	 AnalyzerCore::MakeHistograms("cutflow", 5,0.,5.);

	 GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
	 GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
	 GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
	 GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
	 GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"DiMu_tight");
	
	 
  }
}


void ExampleAnalyzerDiMuon::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ExampleAnalyzerDiMuon::MakeHistograms(){
  //// Additional plots to make
	 
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
	*  Remove//Overide this ExampleAnalyzerDiMuonCore::MakeHistograms() to make new hists for your analysis
	**/
  
}


void ExampleAnalyzerDiMuon::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



