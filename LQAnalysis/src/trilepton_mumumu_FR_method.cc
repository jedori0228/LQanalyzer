// $Id: trilepton_mumumu_FR_method.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_FR_method Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
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
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
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
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:
  //ResetLumiMask(snu::KEvent::gold);

  string lqdir = getenv("LQANALYZER_DIR");
  TFile* file[5];

  //==== dijet topology
  //file[0] = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_SingleMuonTrigger_Dijet.root").c_str() );
  //hist_trimuon_FR[0] = (TH2F*)file[0]->Get("SingleMuonTrigger_Dijet_events_F")->Clone();
  //==== HighdXY muons
  file[1] = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_SingleMuonTrigger_HighdXY.root").c_str() );
  hist_trimuon_FR[1] = (TH2F*)file[1]->Get("SingleMuonTrigger_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons
  //file[2] = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_DiMuonTrigger_HighdXY.root").c_str() );
  //hist_trimuon_FR[2] = (TH2F*)file[2]->Get("DiMuonTrigger_HighdXY_events_F")->Clone();
  //==== DiMuonHighdXY muons + n_jet bins
  //file[3] = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_DiMuonTrigger_HighdXY_0jet.root").c_str() );
  //hist_trimuon_FR[3] = (TH2F*)file[3]->Get("DiMuonTrigger_HighdXY_0jet_events_F")->Clone();
  //file[4] = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_DiMuonTrigger_HighdXY_withjet.root").c_str() );
  //hist_trimuon_FR[4] = (TH2F*)file[4]->Get("DiMuonTrigger_HighdXY_withjet_events_F")->Clone();

  for(int i=1; i<2; i++){
    TH1I* hist_bins = (TH1I*)file[i]->Get("hist_bins");
    FR_n_pt_bin[i] = hist_bins->GetBinContent(1);
    FR_n_eta_bin[i] = hist_bins->GetBinContent(2);
    delete hist_bins;
    file[i]->Close();
    delete file[i];
  }

  TFile* file_FR_SF = new TFile( (lqdir+"/data/rootfiles/13TeV_trimuon_FR_SF_SingleMuonTrigger_QCD_mu.root").c_str() );
  hist_trimuon_FR_SF = (TH2F*)file_FR_SF->Get("SingleMuonTrigger_MCTruth_events_F");
  hist_trimuon_FR_SF_pt = (TH1F*)file_FR_SF->Get("SingleMuonTrigger_MCTruth_pt_F");

  return;


}


void trilepton_mumumu_FR_method::ExecuteEvents()throw( LQError ){

  weight = 1.0; //initializing 
 
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");


  if(!PassTrigger(triggerslist, prescale)) return;
  FillCutFlow("TriggerCut", weight);
  // Trigger matching is done using KMuon::TriggerMatched(TString) which returns a bool

  /* // #### CAT::: trigger matching information is stored for muons and electrons for:
  ///HLT_IsoMu24_eta2p1_v
  ///HLT_Mu17_Mu8_DZ_v
  ///HLT_Mu17_TkMu8_DZ_v
  ///HLT_IsoMu20
  ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
  ///HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele12_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
  ///HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_
  ///HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v
  ///HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v
  */

  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;


  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", weight);


  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  //std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  //std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  //std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  //std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed
  std::vector<snu::KMuon> muontriTightColl = GetMuons("MUON_HN_TRI_TIGHT"); 
  std::vector<snu::KMuon> muontriLooseColl = GetMuons("MUON_HN_TRI_LOOSE");
  
 // CorrectMuonMomentum(muonTightColl);
  //float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonTightColl,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN");// pt > 20 ; eta < 2.5; PFlep veto; pileup ID
  
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  numberVertices = eventbase->GetEvent().nVertices();

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    //pileup_reweight = eventbase->GetEvent().PileUpWeight();

  }
  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);


  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    //weight*=weight_trigger_sf;
  }

  int n_triTight_muons = muontriTightColl.size();
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_hn.size();

  FillHist("GenWeight_NJet" , n_jets*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  if( n_triLoose_muons != 3 ) return;

  if( muontriLooseColl.at(0).Pt() < 20. ) return;
  FillCutFlow("3muon", weight);

  snu::KParticle lep[3], HN[4];
  vector<double> FR_muon, FR_error_muon;
  FR_muon.clear();
  FR_error_muon.clear();
  //for(int i=0;i<k_flags.size();i++) cout << "k_flags = " << k_flags.at(i) << endl;
  for(int i=0;i<3;i++){
    lep[i] = muontriLooseColl.at(i);
    //==== find loose but not tight muon ( 0.1 < RelIso (< 0.6) )
    if( muontriLooseColl.at(i).RelIso04() > 0.1 ){
      FR_muon.push_back( get_FR(lep[i], k_flags.at(0), n_jets, false) );
      FR_error_muon.push_back( get_FR(lep[i], k_flags.at(0), n_jets, true) );
    }
  }


  //==== fake method weighting
  if( n_triTight_muons == 3 ) return; // return TTT case
  for(unsigned int i=0; i<FR_muon.size(); i++){
    weight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
  }
  if( FR_muon.size() == 2 ) weight *= -1.; // minus sign for TLL
  //==== weight error
  double weight_err(0.);
  if( FR_muon.size() == 1 ){
    double fr1 = FR_muon.at(0);
    double fr1_err = FR_error_muon.at(0);
    weight_err = fr1_err/pow(fr1-1,2);
  }
  else if( FR_muon.size() == 2 ){
    double fr1 = FR_muon.at(0);
    double fr1_err = FR_error_muon.at(0);
    double fr2 = FR_muon.at(1);
    double fr2_err = FR_error_muon.at(1);
    weight_err = sqrt( pow( fr1_err*fr2*(1-fr2),2) +
                       pow( fr2_err*fr1*(1-fr1),2)   ) / pow( (1-fr1)*(1-fr2), 2 );
  }
  else if( FR_muon.size() == 3 ){
    double fr1 = FR_muon.at(0);
    double fr1_err = FR_error_muon.at(0);
    double fr2 = FR_muon.at(1);
    double fr2_err = FR_error_muon.at(1);
    double fr3 = FR_muon.at(2);
    double fr3_err = FR_error_muon.at(2);
    weight_err = sqrt( pow( fr1_err*fr2*(1-fr2)*fr3*(1-fr3), 2) +
                       pow( fr2_err*fr3*(1-fr3)*fr1*(1-fr1), 2) +
                       pow( fr3_err*fr1*(1-fr1)*fr2*(1-fr2), 2)   ) / pow( (1-fr1)*(1-fr2)*(1-fr3) ,2 );
  }
  else{
    Message("?", INFO);
  }

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

  // MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
  // MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
  // POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
  if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
      (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) return;
  FillCutFlow("mllsf4", weight);

  ///////////////////////////////////////////
  ////////// m(HN) < 80 GeV region //////////
  ///////////////////////////////////////////

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
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

  FillUpDownHist("HN_mass_class1_cut0", HN[0].M(), weight, weight_err, 0., 500., 50);
  FillUpDownHist("HN_mass_class2_cut0", HN[1].M(), weight, weight_err, 0., 500., 50);
  FillUpDownHist("HN_mass_class3_cut0", HN[2].M(), weight, weight_err, 0, 500, 50);
  FillUpDownHist("HN_mass_class4_cut0", HN[3].M(), weight, weight_err, 0, 1000, 100);
  FillUpDownHist("W_pri_lowmass_mass_cut0", W_pri_lowmass.M(), weight, weight_err, 70., 500., 43);
  FillUpDownHist("W_pri_highmass_mass_cut0", W_pri_highmass.M(), weight, weight_err, 0, 1000, 1000);
  FillUpDownHist("deltaR_OS_min_cut0", deltaR_OS_min, weight, weight_err, 0, 5, 50);
  FillUpDownHist("gamma_star_mass_cut0", gamma_star.M(), weight, weight_err, 0., 120., 120);
  FillUpDownHist("z_candidate_mass_cut0", z_candidate.M(), weight, weight_err, 0., 120., 120);
  FillUpDownHist("n_jets_cut0", n_jets, weight, weight_err, 0, 10, 10);
  FillUpDownHist("n_events_cut0", 0, weight, weight_err, 0, 1, 1);
  FillCLHist(trilephist, "cut0_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
  FillCLHist(trilephist, "cut0_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
  FillCLHist(trilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
  if( is_deltaR_OS_min_0p5 ){
    FillUpDownHist("HN_mass_class1_cutdR", HN[0].M(), weight, weight_err, 0., 500., 50);
    FillUpDownHist("HN_mass_class2_cutdR", HN[1].M(), weight, weight_err, 0., 500., 50);
    FillUpDownHist("HN_mass_class3_cutdR", HN[2].M(), weight, weight_err, 0, 500, 50);
    FillUpDownHist("HN_mass_class4_cutdR", HN[3].M(), weight, weight_err, 0, 1000, 100);
    FillUpDownHist("W_pri_lowmass_mass_cutdR", W_pri_lowmass.M(), weight, weight_err, 70., 500., 43);
    FillUpDownHist("W_pri_highmass_mass_cutdR", W_pri_highmass.M(), weight, weight_err, 0, 1000, 1000);
    FillUpDownHist("deltaR_OS_min_cutdR", deltaR_OS_min, weight, weight_err, 0, 5, 50);
    FillUpDownHist("gamma_star_mass_cutdR", gamma_star.M(), weight, weight_err, 0., 120., 120);
    FillUpDownHist("z_candidate_mass_cutdR", z_candidate.M(), weight, weight_err, 0., 120., 120);
    FillUpDownHist("n_jets_cutdR", n_jets, weight, weight_err, 0, 10, 10);
    FillUpDownHist("n_events_cutdR", 0, weight, weight_err, 0, 1, 1);
    FillCLHist(trilephist, "cutdR_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
    FillCLHist(trilephist, "cutdR_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
    FillCLHist(trilephist, "cutdR", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
    if( is_W_pri_lowmass_100 ){
      FillUpDownHist("HN_mass_class1_cutdR_cutW", HN[0].M(), weight, weight_err, 0., 100., 10);
      FillUpDownHist("HN_mass_class2_cutdR_cutW", HN[1].M(), weight, weight_err, 0., 100., 10);
      FillUpDownHist("W_pri_lowmass_mass_cutdR_cutW", W_pri_lowmass.M(), weight, weight_err, 70., 120., 10);
      FillUpDownHist("W_pri_highmass_mass_cutdR_cutW", W_pri_highmass.M(), weight, weight_err, 0, 1000, 1000);
      FillUpDownHist("deltaR_OS_min_cutdR_cutW", deltaR_OS_min, weight, weight_err, 0, 5, 50);
      FillUpDownHist("gamma_star_mass_cutdR_cutW", gamma_star.M(), weight, weight_err, 0., 120., 24);
      FillUpDownHist("z_candidate_mass_cutdR_cutW", z_candidate.M(), weight, weight_err, 0., 120., 24);
      FillUpDownHist("n_jets_cutdR_cutW", n_jets, weight, weight_err, 0, 10, 10);
      FillUpDownHist("n_events_cutdR_cutW", 0, weight, weight_err, 0, 1, 1);
      FillCLHist(trilephist, "cutdR_cutW_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
      FillCLHist(trilephist, "cutdR_cutW_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
      FillCLHist(trilephist, "cutdR_cutW", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
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
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
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
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS"); 
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"mllsf4");
    
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

double trilepton_mumumu_FR_method::get_FR(snu::KParticle muon, TString whichFR, int n_jets, bool geterror){

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
      //cout << " " << ptarray[FR_n_pt_bin[FR_index]] << endl;
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
  //cout << "this_FR = " << this_FR << endl;
  double FR_SF = hist_trimuon_FR_SF->GetBinContent(this_pt_bin, this_eta_bin);
  double FR_SF_pt = hist_trimuon_FR_SF_pt->GetBinContent(this_pt_bin+1); // +1 : SF starts from [0,10] bin..
  //cout << "this_SF = " << FR_SF << endl;
  double this_FR_error = hist_trimuon_FR[FR_index]->GetBinError(this_pt_bin, this_eta_bin);

  if(geterror){
    if(k_flags.at(1) == "SF") return this_FR_error*FR_SF;
    else if(k_flags.at(1) == "SF_pt") return this_FR_error*FR_SF_pt;
    else return this_FR_error;
  }
  else{
    if(k_flags.at(1) == "SF") return this_FR*FR_SF;
    else if(k_flags.at(1) == "SF_pt") return this_FR*FR_SF_pt;
    else return this_FR;
  }

}



