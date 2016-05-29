// $Id: trilepton_mumumu_FR_method.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_FR_method Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond     <jalmond@cern.ch>        - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_FR_method.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                  
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_FR_method);


/**
 *  This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_mumumu_FR_method::trilepton_mumumu_FR_method() :  AnalyzerCore(), out_muons(0)
{

  
  // To have the correct name in the log:                                                                                   
  SetLogName("trilepton_mumumu_FR_method");

  Message("In trilepton_mumumu_FR_method constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  
}


void trilepton_mumumu_FR_method::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)  output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  string lqdir = getenv("LQANALYZER_DIR");

  //==== dijet topology
  //TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_Dijet.root").c_str() );
  //hist_trimuon_FR[0] = (TH2F*)file->Get("Dijet_events_F")->Clone();
  //==== HighdXY muons
  //TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_SingleMuon_HighdXY.root").c_str() );
  //hist_trimuon_FR[0] = (TH2F*)file->Get("SingleMuon_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons
  //TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuon_HighdXY.root").c_str() );
  //hist_trimuon_FR[0] = (TH2F*)file->Get("DiMuon_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons + n_jet bins
  TFile* file = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuon_HighdXY_0jet.root").c_str() );
  hist_trimuon_FR[0] = (TH2F*)file->Get("DiMuon_HighdXY_0jet_events_F")->Clone();
  TFile* file_withjet = new TFile( (lqdir+"/data/rootfiles/8TeV_trimuon_FR_DiMuon_HighdXY_withjet.root").c_str() );
  hist_trimuon_FR[1] = (TH2F*)file_withjet->Get("DiMuon_HighdXY_withjet_events_F")->Clone();

  TH1I* hist_bins = (TH1I*)file->Get("hist_bins");
  FR_n_pt_bin = hist_bins->GetBinContent(1);
  FR_n_eta_bin = hist_bins->GetBinContent(2);

  return;
}


void trilepton_mumumu_FR_method::ExecuteEvents()throw( LQError ){

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


  //if(n_triTight_muons != 3) return;
  //FillHist("TTT", 0, weight*pileup_reweight, 0, 1, 1);
  //vector<snu::KMuon> muons_nofake = GetTruePrompt(muontriTightColl, false);
  //if(muons_nofake.size() != 3){
    //return;
    //FillHist("fake_included", 0, weight*pileup_reweight, 0, 1, 1);
  //}
  //return;

  if( n_triLoose_muons != 3 ) return;

  if( muontriLooseColl.at(0).Pt() < 20. ) return;
  FillCutFlow("3muon", weight);

  snu::KParticle lep[3], HN[4];
  vector<double> FR_muon;
  for(int i=0;i<3;i++){
    lep[i] = muontriLooseColl.at(i);
    if( !eventbase->GetMuonSel()->HNIstriTight(muontriLooseColl.at(i) , false) ){
      FR_muon.push_back( get_FR(lep[i], n_jets) ); 
    }
  }


  // fake method weighting
  if( n_triTight_muons == 3 ) return; // return TTT case
  for(unsigned int i=0; i<FR_muon.size(); i++){
    weight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
  }
  if( FR_muon.size() == 2 ) weight *= -1.; // minus sign for TLL

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


  // MC samples has m(ll)_saveflavour > 4 GeV cut at gen level
  // MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
  // POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
  if( ! (lep[0]+lep[1]).M() > 4 || ! (lep[0]+lep[2]).M() > 4 || ! (lep[1]+lep[2]).M() > 4 ) return;
  FillCutFlow("mllsf4", weight);

  int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
  if(lep[0].Charge() * lep[1].Charge() > 0){ // Q(0) = Q(1)
    if(lep[1].Charge() * lep[2].Charge() < 0){ // Q(1) != Q(2)
      OppSign = 2;
      SameSign[0] = 0;
      SameSign[1] = 1;
    }
    else return; // veto Q(0) = Q(1) = Q(2)
  }
  else{ // Q(0) != Q(1)
    if(lep[0].Charge() * lep[2].Charge() > 0){ // Q(0) = Q(2)
      OppSign = 1;
      SameSign[0] = 0;
      SameSign[1] = 2;
    }
    else if(lep[1].Charge() * lep[2].Charge() > 0){ // Q(1) = Q(2)
      OppSign = 0;
      SameSign[0] = 1;
      SameSign[1] = 2;
    }
  } // Find l2 and assign l1&l3 in ptorder 
  FillCutFlow("2SS1OS", weight);

  ///////////////////////////////////////////
  ////////// m(HN) < 80 GeV region //////////
  ///////////////////////////////////////////

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.PFMET(), METphi = Evt.PFMETphi();
  snu::KParticle W_pri_lowmass, nu_lowmass, gamma_star, z_candidate;
  nu_lowmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
  double pz_sol_lowmass[2];
  pz_sol_lowmass[0] = solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"); // 0 = minus
  pz_sol_lowmass[1] = solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 1 = plus
  //PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
  //PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus

  int solution_selection_lowmass = 0;
  if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
    // take the one with smaller magnitude
    if( fabs( pz_sol_lowmass[0] ) > fabs( pz_sol_lowmass[1] ) ){
      solution_selection_lowmass = 1;
    }
  }
  
  // reconstruct HN and W_real 4-vec with selected Pz solution
  PutNuPz(&nu_lowmass, pz_sol_lowmass[solution_selection_lowmass] );
  // SameSign[0] : leading among SS
  // SameSign[1] : subleading among SS
  // [class1]
  // HN40, HN50 - SS_leading is primary
  // [class2]
  // HN60       - SS_subleading is primary

  HN[0] = lep[OppSign] + lep[SameSign[1]] + nu_lowmass; // [class1]
  HN[1] = lep[OppSign] + lep[SameSign[0]] + nu_lowmass; // [class2]
  W_pri_lowmass = lep[0] + lep[1] + lep[2] + nu_lowmass;
  

  double deltaR_OS_min;
  if( lep[OppSign].DeltaR(lep[SameSign[0]]) < lep[OppSign].DeltaR(lep[SameSign[1]]) ){
    deltaR_OS_min = lep[OppSign].DeltaR(lep[SameSign[0]]);
    gamma_star = lep[OppSign] + lep[SameSign[0]];
  }
  else{
    deltaR_OS_min = lep[OppSign].DeltaR(lep[SameSign[1]]);
    gamma_star = lep[OppSign] + lep[SameSign[1]];
  }

  if( fabs( (lep[OppSign] + lep[SameSign[0]]).M() - 91.1876 ) <
      fabs( (lep[OppSign] + lep[SameSign[1]]).M() - 91.1876 )   ){
    z_candidate = lep[OppSign] + lep[SameSign[0]];
  }
  else{
    z_candidate = lep[OppSign] + lep[SameSign[1]];
  }

  ///////////////////////////////////////////
  ////////// m(HN) > 80 GeV region //////////
  ///////////////////////////////////////////

  snu::KParticle W_pri_highmass, nu_highmass, W_sec;
  nu_highmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
  int l_3_index = find_mlmet_closest_to_W(lep, nu_highmass);
  double pz_sol_highmass[2]; 
  pz_sol_highmass[0] = solveqdeq(80.4, lep[l_3_index], MET, METphi, "m"); // 0 = minus
  pz_sol_highmass[1] = solveqdeq(80.4, lep[l_3_index], MET, METphi, "p"); // 1 = plus
  int solution_selection_highmass = 0;
  if( pz_sol_highmass[0] != pz_sol_highmass[1] ){ 
    // take the one with smaller magnitude
    if( fabs( pz_sol_highmass[0] ) > fabs( pz_sol_highmass[1] ) ){
      solution_selection_highmass = 1;
    }
  }
  PutNuPz( &nu_highmass, pz_sol_highmass[solution_selection_highmass] );

  W_pri_highmass = lep[0] + lep[1] + lep[2] + nu_highmass;

  // [class3]
  // m(HN) : 90 ~ 500 GeV - primary lepton has larger pT
  // [class4]
  // m(HN) : 700 ~ 1000 GeV - primary lepton has smaller pT

  W_sec = lep[l_3_index] + nu_highmass;

  if(l_3_index == OppSign){
     
    HN[2] = W_sec + lep[SameSign[1]]; // [class3]
    HN[3] = W_sec + lep[SameSign[0]]; // [class4]
 
  }
  else{
      HN[2] = W_sec + lep[OppSign]; // [class3]
      HN[3] = W_sec + lep[OppSign]; // [class4]
  }

  bool is_deltaR_OS_min_0p5 = deltaR_OS_min > 0.5;
  bool is_W_pri_lowmass_100 = W_pri_lowmass.M() < 100;

  //==== No PU
  FillHist("HN_mass_class1_cut0", HN[0].M(), weight, 0., 500., 50);
  FillHist("HN_mass_class2_cut0", HN[1].M(), weight, 0., 500., 50);
  FillHist("HN_mass_class3_cut0", HN[2].M(), weight, 0, 500, 50);
  FillHist("HN_mass_class4_cut0", HN[3].M(), weight, 0, 1000, 100);
  FillHist("W_pri_lowmass_mass_cut0", W_pri_lowmass.M(), weight, 70., 500., 43);
  FillHist("W_pri_highmass_mass_cut0", W_pri_highmass.M(), weight, 0, 1000, 1000);
  FillHist("deltaR_OS_min_cut0", deltaR_OS_min, weight, 0, 5, 50);
  FillHist("gamma_star_mass_cut0", gamma_star.M(), weight, 0., 120., 120);
  FillHist("z_candidate_mass_cut0", z_candidate.M(), weight, 0., 120., 120);
  FillHist("n_jets_cut0", n_jets, weight, 0, 10, 10);
  FillCLHist(trilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight);
  FillHist("n_events_cut0", 0, weight, 0, 1, 1);
  //==== PU
  FillHist("HN_mass_class1_cut0_PU", HN[0].M(), weight*pileup_reweight, 0., 500., 50);
  FillHist("HN_mass_class2_cut0_PU", HN[1].M(), weight*pileup_reweight, 0., 500., 50);
  FillHist("HN_mass_class3_cut0_PU", HN[2].M(), weight*pileup_reweight, 0, 500, 50);
  FillHist("HN_mass_class4_cut0_PU", HN[3].M(), weight*pileup_reweight, 0, 1000, 100);
  FillHist("W_pri_lowmass_mass_cut0_PU", W_pri_lowmass.M(), weight*pileup_reweight, 70., 500., 43);
  FillHist("W_pri_highmass_mass_cut0_PU", W_pri_highmass.M(), weight*pileup_reweight, 0, 1000, 1000);
  FillHist("deltaR_OS_min_cut0_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 50);
  FillHist("gamma_star_mass_cut0_PU", gamma_star.M(), weight*pileup_reweight, 0., 120., 120);
  FillHist("z_candidate_mass_cut0_PU", z_candidate.M(), weight*pileup_reweight, 0., 120., 120);
  FillHist("n_jets_cut0_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  FillCLHist(trilephist, "cut0_PU", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight*pileup_reweight);
  FillHist("n_events_cut0_PU", 0, weight*pileup_reweight, 0, 1, 1);
  if( is_deltaR_OS_min_0p5 ){
    //==== No PU
    FillHist("HN_mass_class1_cutdR", HN[0].M(), weight, 0., 500., 50);
    FillHist("HN_mass_class2_cutdR", HN[1].M(), weight, 0., 500., 50);
    FillHist("HN_mass_class3_cutdR", HN[2].M(), weight, 0, 500, 50);
    FillHist("HN_mass_class4_cutdR", HN[3].M(), weight, 0, 1000, 100);
    FillHist("W_pri_lowmass_mass_cutdR", W_pri_lowmass.M(), weight, 70., 500., 43);
    FillHist("W_pri_highmass_mass_cutdR", W_pri_highmass.M(), weight, 0, 1000, 1000);
    FillHist("deltaR_OS_min_cutdR", deltaR_OS_min, weight, 0, 5, 50);
    FillHist("gamma_star_mass_cutdR", gamma_star.M(), weight, 0., 120., 120);
    FillHist("z_candidate_mass_cutdR", z_candidate.M(), weight, 0., 120., 120);
    FillHist("n_jets_cutdR", n_jets, weight, 0, 10, 10);
    FillCLHist(trilephist, "cutdR", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight);
    FillHist("n_events_cutdR", 0, weight, 0, 1, 1);
    //==== PU
    FillHist("HN_mass_class1_cutdR_PU", HN[0].M(), weight*pileup_reweight, 0., 500., 50);
    FillHist("HN_mass_class2_cutdR_PU", HN[1].M(), weight*pileup_reweight, 0., 500., 50);
    FillHist("HN_mass_class3_cutdR_PU", HN[2].M(), weight*pileup_reweight, 0, 500, 50);
    FillHist("HN_mass_class4_cutdR_PU", HN[3].M(), weight*pileup_reweight, 0, 1000, 100);
    FillHist("W_pri_lowmass_mass_cutdR_PU", W_pri_lowmass.M(), weight*pileup_reweight, 70., 500., 43);
    FillHist("W_pri_highmass_mass_cutdR_PU", W_pri_highmass.M(), weight*pileup_reweight, 0, 1000, 1000);
    FillHist("deltaR_OS_min_cutdR_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 50);
    FillHist("gamma_star_mass_cutdR_PU", gamma_star.M(), weight*pileup_reweight, 0., 120., 120);
    FillHist("z_candidate_mass_cutdR_PU", z_candidate.M(), weight*pileup_reweight, 0., 120., 120);
    FillHist("n_jets_cutdR_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
    FillCLHist(trilephist, "cutdR_PU", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight*pileup_reweight);
    FillHist("n_events_cutdR_PU", 0, weight*pileup_reweight, 0, 1, 1);
    if( is_W_pri_lowmass_100 ){
      //==== No PU
      FillHist("HN_mass_class1_cutdR_cutW", HN[0].M(), weight, 0., 100., 10);
      FillHist("HN_mass_class2_cutdR_cutW", HN[1].M(), weight, 0., 100., 10);
      FillHist("W_pri_lowmass_mass_cutdR_cutW", W_pri_lowmass.M(), weight, 70., 120., 10);
      FillHist("W_pri_highmass_mass_cutdR_cutW", W_pri_highmass.M(), weight, 0, 1000, 1000);
      FillHist("deltaR_OS_min_cutdR_cutW", deltaR_OS_min, weight, 0, 5, 50);
      FillHist("gamma_star_mass_cutdR_cutW", gamma_star.M(), weight, 0., 120., 24);
      FillHist("z_candidate_mass_cutdR_cutW", z_candidate.M(), weight, 0., 120., 24);
      FillHist("n_jets_cutdR_cutW", n_jets, weight, 0, 10, 10);
      FillCLHist(trilephist, "cutdR_cutW", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight);
      FillHist("n_events_cutdR_cutW", 0, weight, 0, 1, 1);
      //==== PU
      FillHist("HN_mass_class1_cutdR_cutW_PU", HN[0].M(), weight*pileup_reweight, 0., 100., 10);
      FillHist("HN_mass_class2_cutdR_cutW_PU", HN[1].M(), weight*pileup_reweight, 0., 100., 10);
      FillHist("W_pri_lowmass_mass_cutdR_cutW_PU", W_pri_lowmass.M(), weight*pileup_reweight, 70., 120., 10);
      FillHist("W_pri_highmass_mass_cutdR_cutW_PU", W_pri_highmass.M(), weight*pileup_reweight, 0, 1000, 1000);
      FillHist("deltaR_OS_min_cutdR_cutW_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 50);
      FillHist("gamma_star_mass_cutdR_cutW_PU", gamma_star.M(), weight*pileup_reweight, 0., 120., 24);
      FillHist("z_candidate_mass_cutdR_cutW_PU", z_candidate.M(), weight*pileup_reweight, 0., 120., 24);
      FillHist("n_jets_cutdR_cutW_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
      FillCLHist(trilephist, "cutdR_cutW_PU", eventbase->GetEvent(), muontriLooseColl, electronTightColl, jetColl_lepveto, weight*pileup_reweight);
      FillHist("n_events_cutdR_cutW_PU", 0, weight*pileup_reweight, 0, 1, 1);
    }
  }

  return;
}// End of execute event loop
  


void trilepton_mumumu_FR_method::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_FR_method::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_FR_method::~trilepton_mumumu_FR_method() {
  
  Message("In trilepton_mumumu_FR_method Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
}


void trilepton_mumumu_FR_method::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_FR_method::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_FR_method::MakeHistograms(){
  //// Additional plots to make
   
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
  *  Remove//Overide this trilepton_mumumu_FR_methodCore::MakeHistograms() to make new hists for your analysis
  **/
  
}


void trilepton_mumumu_FR_method::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

double trilepton_mumumu_FR_method::get_FR(snu::KParticle muon, int n_jets){

  double this_pt = muon.Pt();
  double this_eta = fabs( muon.Eta() );

  if(n_jets > 0) n_jets = 1;
  if(!hist_trimuon_FR[1]) n_jets = 0;

  // FR_n_pt_bin = 7
  // array index      0    1    2    3    4    5    6    7
  // bin numbe          1    2     3    4    5   6     7
  // ptarray[7+1] = {10., 15., 20., 25., 30., 35., 45., 60.}; 

  double ptarray[FR_n_pt_bin+1], etaarray[FR_n_eta_bin+1];
  //cout << "FR_n_pt_bin = " << FR_n_pt_bin << endl;
  for(int i=0; i<FR_n_pt_bin; i++){
    ptarray[i] = hist_trimuon_FR[n_jets]->GetXaxis()->GetBinLowEdge(i+1);
    //cout << " " << ptarray[i] << endl;
    if(i==FR_n_pt_bin-1){
      ptarray[FR_n_pt_bin] = hist_trimuon_FR[n_jets]->GetXaxis()->GetBinUpEdge(i+1);
      //cout << " " << ptarray[FR_n_pt_bin] << endl;
    }
  }
  //cout << "FR_n_eta_bin = " << FR_n_eta_bin << endl;
  for(int i=0; i<FR_n_eta_bin; i++){
    etaarray[i] = hist_trimuon_FR[n_jets]->GetYaxis()->GetBinLowEdge(i+1);
    //cout << " " << etaarray[i] << endl;
    if(i==FR_n_eta_bin-1){
      etaarray[FR_n_eta_bin] = hist_trimuon_FR[n_jets]->GetYaxis()->GetBinUpEdge(i+1);
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

  return hist_trimuon_FR[n_jets]->GetBinContent(this_pt_bin, this_eta_bin);

}




