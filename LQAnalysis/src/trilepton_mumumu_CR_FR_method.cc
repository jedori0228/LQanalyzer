// $Id: trilepton_mumumu_CR_FR_method.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_CR_FR_method Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond     <jalmond@cern.ch>        - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_CR_FR_method.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                  
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_CR_FR_method);


/**
 *  This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_mumumu_CR_FR_method::trilepton_mumumu_CR_FR_method() :  AnalyzerCore(), out_muons(0)
{

  
  // To have the correct name in the log:                                                                                   
  SetLogName("trilepton_mumumu_CR_FR_method");

  Message("In trilepton_mumumu_CR_FR_method constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  
}


void trilepton_mumumu_CR_FR_method::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)  output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  string lqdir = getenv("LQANALYZER_DIR");

  //==== dijet topology
  //TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_dijet_topology.root").c_str() );
  //hist_trimuon_FR = (TH2F*)file->Get("events_F")->Clone();
  //==== HighdXY muons
  //TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_HighdXY.root").c_str() );
  //hist_trimuon_FR = (TH2F*)file->Get("HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons
  TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuon_HighdXY.root").c_str() );
  hist_trimuon_FR = (TH2F*)file->Get("DiMuon_HighdXY_events_F")->Clone();

  TH1I* hist_bins = (TH1I*)file->Get("hist_bins");
  FR_n_pt_bin = hist_bins->GetBinContent(1);
  FR_n_eta_bin = hist_bins->GetBinContent(2);

  return;
 }


void trilepton_mumumu_CR_FR_method::ExecuteEvents()throw( LQError ){

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
  eventbase->GetJetSel()->SetEta(5.0);
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

  // define CR         //
  // SS dimuon + 0 jet //
  if( n_triLoose_muons != 2 ) return;
  if( n_triTight_muons == 2 ) return; // return for TT case
  if( n_jets != 0 ) return;
  if( muontriLooseColl.at(0).Charge() != muontriLooseColl.at(1).Charge() ) return;
  if( muontriLooseColl.at(0).Pt() < 15 ) return;

  //if(n_triTight_muons != 3) return;
  //FillHist("TTT", 0, weight*pileup_reweight, 0, 1, 1);
  //vector<snu::KMuon> muons_nofake = GetTruePrompt(muontriTightColl, false);
  //if(muons_nofake.size() != 3){
    //return;
    //FillHist("fake_included", 0, weight*pileup_reweight, 0, 1, 1);
  //}
  //return;

  snu::KParticle lep[2];
  vector<double> FR_muon;
  for(int i=0;i<2;i++){
    lep[i] = muontriLooseColl.at(i);
    if( !eventbase->GetMuonSel()->HNIstriTight(muontriLooseColl.at(i) , false) ){
      FR_muon.push_back( get_FR(lep[i]) );
    }
  }

  // fake method weighting
  for(unsigned int i=0; i<FR_muon.size(); i++){
    weight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
  }
  if( FR_muon.size() == 2 ) weight *= -1.; // minus sign for LL

/*
  if( n_triTight_muons == 0 ) weight *= FR_muon.at(0)*FR_muon.at(1)*FR_muon.at(2)/( (1.-FR_muon.at(0))*(1.-FR_muon.at(1))*(1.-FR_muon.at(2)) ); // LLL
  else if ( n_triTight_muons == 1) weight *= -FR_muon.at(0)*FR_muon.at(1)/( (1.-FR_muon.at(0))*(1.-FR_muon.at(1)) ); // TLL
  else if ( n_triTight_muons == 2) weight *= FR_muon.at(0)/(1.-FR_muon.at(0)); // TTL
  else return; // return TTT event
*/
/*
  double FR_muon[3];
  for(int i=0; i<3; i++) FR_muon[i] = get_FR(lep[i]);
  if( n_triTight_muons == 0 ) weight *= FR_muon[0]*FR_muon[1]*FR_muon[2]/( (1.-FR_muon[0])*(1.-FR_muon[1])*(1.-FR_muon[2]) ); // LLL
  else if ( n_triTight_muons == 1) weight *= -FR_muon[1]*FR_muon[2]/( (1.-FR_muon[1])*(1.-FR_muon[2]) ); // TLL
  else if ( n_triTight_muons == 2) weight *= FR_muon[2]/(1.-FR_muon[2]); // TTL
  else return; // return TTT event
*/


  snu::KParticle CR_muon[2];
  CR_muon[0] = muontriLooseColl.at(0);
  CR_muon[1] = muontriLooseColl.at(1);
  FillHist("n_events_control_PU", 0, weight*pileup_reweight, 0, 1, 1);
  FillHist("mll_control_PU", (CR_muon[0]+CR_muon[1]).M() , weight*pileup_reweight, 0, 200, 200);
  FillHist("leadingLepton_Pt_control_PU", CR_muon[0].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("secondLepton_Pt_control_PU", CR_muon[1].Pt() , weight*pileup_reweight, 0, 200, 200);
  FillHist("leadingLepton_Eta_control_PU", CR_muon[0].Eta() , weight*pileup_reweight, -3, 3, 6./0.1);
  FillHist("secondLepton_Eta_control_PU", CR_muon[1].Eta() , weight*pileup_reweight, -3, 3, 6./0.1);

  return;
}// End of execute event loop
  


void trilepton_mumumu_CR_FR_method::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_CR_FR_method::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_CR_FR_method::~trilepton_mumumu_CR_FR_method() {
  
  Message("In trilepton_mumumu_CR_FR_method Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
}


void trilepton_mumumu_CR_FR_method::FillCutFlow(TString cut, float weight){

  
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
   GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"mllsf4");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"2SS1OS");  
   
  }
}


void trilepton_mumumu_CR_FR_method::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_CR_FR_method::MakeHistograms(){
  //// Additional plots to make
   
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
  *  Remove//Overide this trilepton_mumumu_CR_FR_methodCore::MakeHistograms() to make new hists for your analysis
  **/
  
}


void trilepton_mumumu_CR_FR_method::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

void trilepton_mumumu_CR_FR_method::SetBinInfo(int cut){
  if(cut==0){ // no cut
    HN_x_min=0; HN_x_max=500; HN_dx=10;
    W_pri_lowmass_x_min=70; W_pri_lowmass_x_max=500; W_pri_lowmass_dx=10;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1;
  }
  if(cut==1){ // "dR_OS_min > 0.5"
    HN_x_min=0; HN_x_max=500; HN_dx=10;
    W_pri_lowmass_x_min=70; W_pri_lowmass_x_max=500; W_pri_lowmass_dx=10;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1;
  }
  if(cut==2){ // "dR_OS_min > 0.5", "W_pri < 100 GeV"
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_pri_lowmass_x_min=70; W_pri_lowmass_x_max=120; W_pri_lowmass_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }
  if(cut==3){
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_pri_lowmass_x_min=70; W_pri_lowmass_x_max=120; W_pri_lowmass_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }
  if(cut==4){
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_pri_lowmass_x_min=70; W_pri_lowmass_x_max=120; W_pri_lowmass_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
   z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }

}

double trilepton_mumumu_CR_FR_method::get_FR(snu::KParticle muon){

  double this_pt = muon.Pt();
  double this_eta = fabs( muon.Eta() );

  // FR_n_pt_bin = 7
  // array index      0    1    2    3    4    5    6    7
  // bin numbe          1    2     3    4    5   6     7
  // ptarray[7+1] = {10., 15., 20., 25., 30., 35., 45., 60.}; 

  double ptarray[FR_n_pt_bin+1], etaarray[FR_n_eta_bin+1];
  //cout << "FR_n_pt_bin = " << FR_n_pt_bin << endl;
  for(int i=0; i<FR_n_pt_bin; i++){
    ptarray[i] = hist_trimuon_FR->GetXaxis()->GetBinLowEdge(i+1);
    //cout << " " << ptarray[i] << endl;
    if(i==FR_n_pt_bin-1){
      ptarray[FR_n_pt_bin] = hist_trimuon_FR->GetXaxis()->GetBinUpEdge(i+1);
      //cout << " " << ptarray[FR_n_pt_bin] << endl;
    }
  }
  //cout << "FR_n_eta_bin = " << FR_n_eta_bin << endl;
  for(int i=0; i<FR_n_eta_bin; i++){
    etaarray[i] = hist_trimuon_FR->GetYaxis()->GetBinLowEdge(i+1);
    //cout << " " << etaarray[i] << endl;
    if(i==FR_n_eta_bin-1){
      etaarray[FR_n_eta_bin] = hist_trimuon_FR->GetYaxis()->GetBinUpEdge(i+1);
      //cout << " " << etaarray[FR_n_eta_bin] << endl;
    }
  }

  int this_pt_bin;
  if( this_pt >= ptarray[FR_n_pt_bin] ) this_pt_bin = FR_n_pt_bin;
  else{
    for(int i=0; i<FR_n_pt_bin; i++){
      if( ptarray[i] <= this_pt && this_pt < ptarray[i+1] ){
        this_pt_bin = i+1;
        break;
      }
    }
  }
  int this_eta_bin;
  if( this_eta >= etaarray[FR_n_eta_bin] ) this_eta_bin = FR_n_eta_bin;
  else{
    for(int i=0; i<FR_n_eta_bin; i++){
      if( etaarray[i] <= this_eta && this_eta < etaarray[i+1] ){
        this_eta_bin = i+1;
        break;
      }
    }
  }

  return hist_trimuon_FR->GetBinContent(this_pt_bin, this_eta_bin);

}
