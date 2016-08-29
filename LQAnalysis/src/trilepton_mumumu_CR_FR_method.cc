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
  TFile* file[5];

  //==== dijet topology
  file[0] = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_SingleMuonTrigger_Dijet.root").c_str() );
  hist_trimuon_FR[0] = (TH2F*)file[0]->Get("SingleMuonTrigger_Dijet_events_F")->Clone();
  //==== HighdXY muons
  file[1] = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_SingleMuonTrigger_HighdXY.root").c_str() );
  hist_trimuon_FR[1] = (TH2F*)file[1]->Get("SingleMuonTrigger_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons
  file[2] = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuonTrigger_HighdXY.root").c_str() );
  hist_trimuon_FR[2] = (TH2F*)file[2]->Get("DiMuonTrigger_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons + n_jet bins
  file[3] = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuonTrigger_HighdXY_0jet.root").c_str() );
  hist_trimuon_FR[3] = (TH2F*)file[3]->Get("DiMuonTrigger_HighdXY_0jet_events_F")->Clone();
  file[4] = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuonTrigger_HighdXY_withjet.root").c_str() );
  hist_trimuon_FR[4] = (TH2F*)file[4]->Get("DiMuonTrigger_HighdXY_withjet_events_F")->Clone();

  for(int i=0; i<5; i++){
    TH1I* hist_bins = (TH1I*)file[i]->Get("hist_bins");
    FR_n_pt_bin[i] = hist_bins->GetBinContent(1);
    FR_n_eta_bin[i] = hist_bins->GetBinContent(2);
    delete hist_bins;
    file[i]->Close();
    delete file[i];
  }

  TFile* file_FR_SF = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_SF_SingleMuonTrigger_QCD_single.root").c_str() );
  hist_trimuon_FR_SF = (TH1F*)file_FR_SF->Get("SingleMuonTrigger_MCTruth_pt_F");

  return;
}


void trilepton_mumumu_CR_FR_method::ExecuteEvents()throw( LQError ){

  weight = 1.0; //initializing

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
  std::vector<snu::KMuon> muontriTightColl = GetMuons("HNtriTight");
  std::vector<snu::KMuon> muontriLooseColl = GetMuons("HNtriLoose");
 
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

  //==== CR related variables
  FillHist("n_tight_muons_control_PU", n_triTight_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_loose_muons_control_PU", n_triLoose_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_jets_control_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  int n_bjets=0;
  for(int j=0; j < n_jets; j++){
    if(jetColl_lepveto.at(j).CombinedSecVertexBtag() > 0.679) n_bjets++;
  }
  FillHist("n_bjets_control_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);

  //================
  //==== define CR
  //================

  bool isCR = true;

  int n_muons(0);

  bool isTwoMuon   = n_triLoose_muons == 2 && n_triTight_muons != 2; // No TT case
  bool isThreeMuon = n_triLoose_muons == 3 && n_triTight_muons != 3; // No TTT case

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.PFMET(), METphi = Evt.PFMETphi();

  //==== CR with Two Muons
  if(isTwoMuon){
    snu::KMuon lep[2];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);

		vector<double> FR_muon;
		FR_muon.clear();
		for(int i=0;i<2;i++){
			snu::KMuon this_muon = muontriLooseColl.at(i);
			if( !eventbase->GetMuonSel()->HNtriTightMuonSelection( this_muon ) ){
				FR_muon.push_back( get_FR(this_muon, k_jskim_flag_1, n_jets) );
			}
		}
    double FR_reweight = 1.;
		for(unsigned int i=0; i<FR_muon.size(); i++){
			FR_reweight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
		}
		if( FR_muon.size() == 2 ) FR_reweight *= -1.; // minus sign for LL

    bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();

    std::map< TString, bool > map_whichCR_to_isCR;
    map_whichCR_to_isCR.clear();
    map_whichCR_to_isCR["SS_0jet_vetoLowRes"] = leadPt20 && isSS && (n_jets == 0) && (m_dimuon > 15.);
    map_whichCR_to_isCR["SS_AL1bjet_vetoLowRes"] = leadPt20 && isSS && (n_bjets > 0) && (m_dimuon > 15.);
    map_whichCR_to_isCR["SS_AL2bjet_vetoLowRes"] = leadPt20 && isSS && (n_bjets > 1) && (m_dimuon > 15.);

    for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
      TString this_suffix = it->first;
      if(it->second){
        FillHist("n_events_"+this_suffix+"_PU", 0, weight*pileup_reweight*FR_reweight, 0, 1, 1);
        FillHist("n_jets_"+this_suffix+"_PU", n_jets, weight*pileup_reweight*FR_reweight, 0, 10, 10);
        FillHist("n_bjets_"+this_suffix+"_PU", n_bjets, weight*pileup_reweight*FR_reweight, 0, 10, 10);
        FillHist("PFMET_"+this_suffix+"_PU", MET, weight*pileup_reweight*FR_reweight, 0, 500, 500);
        FillHist("mll_"+this_suffix+"_PU", m_dimuon , weight*pileup_reweight*FR_reweight, 0, 500, 500);
        FillHist("leadingLepton_Pt_"+this_suffix+"_PU", lep[0].Pt() , weight*pileup_reweight*FR_reweight, 0, 200, 200);
        FillHist("leadingLepton_Eta_"+this_suffix+"_PU", lep[0].Eta() , weight*pileup_reweight*FR_reweight, -3, 3, 60);
        FillHist("leadingLepton_RelIso_"+this_suffix+"_PU", lep[0].LeptonRelIso() , weight*pileup_reweight*FR_reweight, 0, 1.0, 100);
        FillHist("leadingLepton_Chi2_"+this_suffix+"_PU", lep[0].GlobalChi2() , weight*pileup_reweight*FR_reweight, 0, 10, 100);
        FillHist("secondLepton_Pt_"+this_suffix+"_PU", lep[1].Pt() , weight*pileup_reweight*FR_reweight, 0, 200, 200);
        FillHist("secondLepton_Eta_"+this_suffix+"_PU", lep[1].Eta() , weight*pileup_reweight*FR_reweight, -3, 3, 60);
        FillHist("secondLepton_RelIso_"+this_suffix+"_PU", lep[1].LeptonRelIso() , weight*pileup_reweight*FR_reweight, 0, 1.0, 100);
        FillHist("secondLepton_Chi2_"+this_suffix+"_PU", lep[1].GlobalChi2() , weight*pileup_reweight*FR_reweight, 0, 10, 100);
      }
    }

  } // isTwoMuon

  if(isThreeMuon){
    snu::KMuon lep[3];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);
    lep[2] = muontriLooseColl.at(2);

    vector<double> FR_muon;
    FR_muon.clear();
    for(int i=0;i<3;i++){
      snu::KMuon this_muon = muontriLooseColl.at(i);
      if( !eventbase->GetMuonSel()->HNtriTightMuonSelection( this_muon ) ){
        FR_muon.push_back( get_FR(this_muon, k_jskim_flag_1, n_jets) );
      }
    }
    double FR_reweight = 1.;
    for(unsigned int i=0; i<FR_muon.size(); i++){
      FR_reweight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
    }
    if( FR_muon.size() == 2 ) FR_reweight *= -1.; // minus sign for TLL

    bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
    bool AllSameCharge = ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) &&
                         ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(2).Charge() );

    if( leadPt20 && !AllSameCharge ){
      snu::KMuon OS, SS[2];
      if     ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ){
        SS[0] = muontriLooseColl.at(0);
        SS[1] = muontriLooseColl.at(1);
        OS    = muontriLooseColl.at(2);
      }
      else if( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(2).Charge() ){
        SS[0] = muontriLooseColl.at(0);
        SS[1] = muontriLooseColl.at(2);
        OS    = muontriLooseColl.at(1);
      }
      else if( muontriLooseColl.at(1).Charge() == muontriLooseColl.at(2).Charge() ){
        SS[0] = muontriLooseColl.at(1);
        SS[1] = muontriLooseColl.at(2);
        OS    = muontriLooseColl.at(0);
      }
      else Message("?", INFO);

      double m_dimuon[2], m_Z = 91.1876;
      m_dimuon[0] = ( SS[0] + OS ).M();
      m_dimuon[1] = ( SS[1] + OS ).M();

      snu::KParticle Z_candidate;
      snu::KMuon ExtraMuon;
      if( fabs(m_dimuon[0]-m_Z) < fabs(m_dimuon[1]-m_Z) ){
        Z_candidate = SS[0] + OS;
        ExtraMuon = SS[1];
      }
      else{
        Z_candidate = SS[1] + OS;
        ExtraMuon = SS[0];
      }

      bool isZresonance = fabs(Z_candidate.M()-m_Z) < 20.;
      bool PtCutOnExtraMuon = ExtraMuon.Pt() > 20.;
      bool METCut = MET > 20.;

      if( isZresonance && PtCutOnExtraMuon && METCut ){
        TString this_suffix = "WZ";
        snu::KParticle nu;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
        snu::KParticle W_candidate = nu+ExtraMuon;

        FillHist("n_events_"+this_suffix+"_PU", 0, weight*pileup_reweight*FR_reweight, 0, 1, 1);
        FillHist("n_jets_"+this_suffix+"_PU", n_jets, weight*pileup_reweight*FR_reweight, 0, 10, 10);
        FillHist("n_bjets_"+this_suffix+"_PU", n_bjets, weight*pileup_reweight*FR_reweight, 0, 10, 10);
        FillHist("PFMET_"+this_suffix+"_PU", MET, weight*pileup_reweight*FR_reweight, 0, 500, 500);
        FillHist("m_Z_candidate_"+this_suffix+"_PU", Z_candidate.M(), weight*pileup_reweight*FR_reweight, 0, 150, 150);
        FillHist("mt_W_candidate_"+this_suffix+"_PU", W_candidate.Mt(), weight*pileup_reweight*FR_reweight, 0, 300, 300);

        FillHist("leadingLepton_Pt_"+this_suffix+"_PU", lep[0].Pt() , weight*pileup_reweight*FR_reweight, 0, 200, 200);
        FillHist("leadingLepton_Eta_"+this_suffix+"_PU", lep[0].Eta() , weight*pileup_reweight*FR_reweight, -3, 3, 60);
        FillHist("leadingLepton_RelIso_"+this_suffix+"_PU", lep[0].LeptonRelIso() , weight*pileup_reweight*FR_reweight, 0, 1.0, 100);
        FillHist("leadingLepton_Chi2_"+this_suffix+"_PU", lep[0].GlobalChi2() , weight*pileup_reweight*FR_reweight, 0, 10, 100);
        FillHist("secondLepton_Pt_"+this_suffix+"_PU", lep[1].Pt() , weight*pileup_reweight*FR_reweight, 0, 200, 200);
        FillHist("secondLepton_Eta_"+this_suffix+"_PU", lep[1].Eta() , weight*pileup_reweight*FR_reweight, -3, 3, 60);
        FillHist("secondLepton_RelIso_"+this_suffix+"_PU", lep[1].LeptonRelIso() , weight*pileup_reweight*FR_reweight, 0, 1.0, 100);
        FillHist("secondLepton_Chi2_"+this_suffix+"_PU", lep[1].GlobalChi2() , weight*pileup_reweight*FR_reweight, 0, 10, 100);
        FillHist("thirdLepton_Pt_"+this_suffix+"_PU", lep[2].Pt() , weight*pileup_reweight*FR_reweight, 0, 200, 200);
        FillHist("thirdLepton_Eta_"+this_suffix+"_PU", lep[2].Eta() , weight*pileup_reweight*FR_reweight, -3, 3, 60);
        FillHist("thirdLepton_RelIso_"+this_suffix+"_PU", lep[2].LeptonRelIso() , weight*pileup_reweight*FR_reweight, 0, 1.0, 100);
        FillHist("thirdLepton_Chi2_"+this_suffix+"_PU", lep[2].GlobalChi2() , weight*pileup_reweight*FR_reweight, 0, 10, 100);

      } // Z Resonance


    } // Not All Same Charge

  } // isThreeMuon


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

double trilepton_mumumu_CR_FR_method::get_FR(snu::KParticle muon, TString whichFR, int n_jets){

  int FR_index = 0;

  if(whichFR=="dijet_topology") FR_index = 0;
  if(whichFR=="HighdXY")        FR_index = 1;
  if(whichFR=="DiMuon_HighdXY") FR_index = 2;
  if(whichFR=="DiMuon_HighdXY_n_jets"){
    FR_index = 3;
    if(n_jets>0) FR_index = 4;
  }

  double this_pt = muon.Pt();
  double this_eta = fabs( muon.Eta() );

  // FR_n_pt_bin = 7
  // array index      0    1    2    3    4    5    6    7
  // bin numbe          1    2     3    4    5   6     7
  // ptarray[7+1] = {10., 15., 20., 25., 30., 35., 45., 60.}; 

  double ptarray[FR_n_pt_bin[FR_index]+1], etaarray[FR_n_eta_bin[FR_index]+1];
  //cout << "FR_n_pt_bin = " << FR_n_pt_bin[FR_index] << endl;
  for(int i=0; i<FR_n_pt_bin[FR_index]; i++){
    ptarray[i] = hist_trimuon_FR[FR_index]->GetXaxis()->GetBinLowEdge(i+1);
    //cout << " " << ptarray[i] << endl;
    if(i==FR_n_pt_bin[FR_index]-1){
      ptarray[FR_n_pt_bin[FR_index]] = hist_trimuon_FR[FR_index]->GetXaxis()->GetBinUpEdge(i+1);
      // cout << " " << ptarray[FR_n_pt_bin[FR_index]] << endl;
    }
  }
  //cout << "FR_n_eta_bin = " << FR_n_eta_bin[FR_index] << endl;
  for(int i=0; i<FR_n_eta_bin[FR_index]; i++){
    etaarray[i] = hist_trimuon_FR[FR_index]->GetYaxis()->GetBinLowEdge(i+1);
    //cout << " " << etaarray[i] << endl;
    if(i==FR_n_eta_bin[FR_index]-1){
      etaarray[FR_n_eta_bin[FR_index]] = hist_trimuon_FR[FR_index]->GetYaxis()->GetBinUpEdge(i+1);
      //cout << " " << etaarray[FR_n_eta_bin[FR_index]] << endl;
    }
  }

  int this_pt_bin;
  if( this_pt >= ptarray[FR_n_pt_bin[FR_index]] ) this_pt_bin = FR_n_pt_bin[FR_index];
  else{
    for(int i=0; i<FR_n_pt_bin[FR_index]; i++){
      if( ptarray[i] <= this_pt && this_pt < ptarray[i+1] ){
        this_pt_bin = i+1;
        break;
      }
    }
  }
  //cout << "this pt bin = " << this_pt_bin << endl;
  int this_eta_bin;
  if( this_eta >= etaarray[FR_n_eta_bin[FR_index]] ) this_eta_bin = FR_n_eta_bin[FR_index];
  else{
    for(int i=0; i<FR_n_eta_bin[FR_index]; i++){
      if( etaarray[i] <= this_eta && this_eta < etaarray[i+1] ){
        this_eta_bin = i+1;
        break;
      }
    }
  }
  //cout << "this eta bin = " << this_eta_bin << endl;

  double this_FR = hist_trimuon_FR[FR_index]->GetBinContent(this_pt_bin, this_eta_bin);
  double FR_SF = hist_trimuon_FR_SF->GetBinContent(this_pt_bin);

  return this_FR;
  //return this_FR*FR_SF;

}


