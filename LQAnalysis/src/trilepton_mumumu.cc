// $Id: trilepton_mumumu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu::trilepton_mumumu() :  AnalyzerCore(),
n_gen_pass(0), sol_sel_chi2_best(0), sol_sel_chi2_plus(0), sol_sel_chi2_minus(0), sol_sel_chi2_smaller(0), sol_sel_chi2_larger(0),
out_muons(0)
{
  
  rmcor = new rochcor2015();
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu");
  
  Message("In trilepton_mumumu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(trilephist,"cut0");
  MakeCleverHistograms(trilephist,"cut0_PU");
  MakeCleverHistograms(trilephist,"cutdR");
  MakeCleverHistograms(trilephist,"cutdR_PU");
  MakeCleverHistograms(trilephist,"cutdR_cutW");
  MakeCleverHistograms(trilephist,"cutdR_cutW_PU");

}


void trilepton_mumumu::InitialiseAnalysis() throw( LQError ) {
  
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


  return;
}


void trilepton_mumumu::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  

  /// Acts on data to remove bad reconstructed event 
  if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;
  

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment

  if(!PassBasicEventCuts()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");

  float trigger_ps_weight= ApplyPrescale("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi,lumimask);

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
  std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed
  std::vector<snu::KMuon> muontriTightColl, muontriLooseColl;
  if( k_sample_name.Contains("HN") ){
    muontriTightColl = GetMuons(BaseSelection::MUON_HN_TRI_TIGHT);
    muontriLooseColl = GetMuons(BaseSelection::MUON_HN_TRI_LOOSE);
  }
  else{
    muontriTightColl = GetMuons(BaseSelection::MUON_HN_TRI_TIGHT, false);
    muontriLooseColl = GetMuons(BaseSelection::MUON_HN_TRI_LOOSE, false);
  }

  CorrectMuonMomentum(muonTightColl);
  float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muontriTightColl, 0); ///MUON_POG_TIGHT == MUON_HN_TIGHT
  muon_id_iso_sf *= MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muontriTightColl, 0);

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl             = GetJets(BaseSelection::JET_NOLEPTONVETO); // All jets
  std::vector<snu::KJet> jetColl_loose       = GetJets(BaseSelection::JET_LOOSE); // pt > 10; eta < 5. ; PFlep veto
  std::vector<snu::KJet> jetColl_tight       = GetJets(BaseSelection::JET_TIGHT);// pt > 20 ; eta < 2.5; PFlep veto
  std::vector<snu::KJet> jetColl_hn          = GetJets(BaseSelection::JET_HN);// pt > 20 ; eta < 2.5; PFlep veto; pileup ID

  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  numberVertices = eventbase->GetEvent().nVertices();

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    pileup_reweight = eventbase->GetEvent().PileUpWeight(lumimask);
    //pileup_reweight = eventbase->GetEvent().AltPileUpWeight(lumimask);

  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);


  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    //weight*=weight_trigger_sf;
    weight*=trigger_ps_weight;
  }

  int n_triTight_muons = muontriTightColl.size();
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_loose.size();

  FillHist("GenWeight_NJet" , n_jets*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  if( n_triLoose_muons != 3 ) return;
  if( n_triTight_muons != 3 ) return;

  if( muontriLooseColl.at(0).Pt() < 20. ) return;
  FillCutFlow("3muon", weight);

  snu::KParticle lep[3], HN[4];
  for(int i=0;i<3;i++){
    lep[i] = muontriLooseColl.at(i);
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

  if(k_sample_name.Contains("HN")) gen_matching();

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
  FillCLHist(trilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight);
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
  FillCLHist(trilephist, "cut0_PU", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight*pileup_reweight);
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
    FillCLHist(trilephist, "cutdR", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight);
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
    FillCLHist(trilephist, "cutdR_PU", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight*pileup_reweight);
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
      FillCLHist(trilephist, "cutdR_cutW", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight);
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
      FillCLHist(trilephist, "cutdR_cutW_PU", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_loose, weight*pileup_reweight);
      FillHist("n_events_cutdR_cutW_PU", 0, weight*pileup_reweight, 0, 1, 1);
    }
  }   




   return;
}// End of execute event loop
  


void trilepton_mumumu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  cout
  << "n_gen_pass = " << n_gen_pass << endl
  << "best = " << sol_sel_chi2_best/n_gen_pass << endl
  << "plus = " << sol_sel_chi2_plus/n_gen_pass << endl
  << "minus = " << sol_sel_chi2_minus/n_gen_pass << endl
  << "smaller = " << sol_sel_chi2_smaller/n_gen_pass << endl
  << "larger = " << sol_sel_chi2_larger/n_gen_pass << endl;

  TH1F* GEN_solution_selection_chi2 = new TH1F("GEN_solution_selection_chi2", "", 6, 0, 6);
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(1, "n_gen_pass");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(2, "best");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(3, "plus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(4, "minus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(5, "smaller");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(6, "larger");
  GEN_solution_selection_chi2->SetBinContent(1, n_gen_pass);
  GEN_solution_selection_chi2->SetBinContent(2, sol_sel_chi2_best);
  GEN_solution_selection_chi2->SetBinContent(3, sol_sel_chi2_plus);
  GEN_solution_selection_chi2->SetBinContent(4, sol_sel_chi2_minus);
  GEN_solution_selection_chi2->SetBinContent(5, sol_sel_chi2_smaller);
  GEN_solution_selection_chi2->SetBinContent(6, sol_sel_chi2_larger);

  m_outputFile->cd();
  GEN_solution_selection_chi2->Write();

}


void trilepton_mumumu::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu::~trilepton_mumumu() {
  
  Message("In trilepton_mumumu Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_mumumu::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_mumumu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

void trilepton_mumumu::gen_matching(){

  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  //==== print truth info
  //cout << "=========================================================" << endl;
  //cout << "truth size = " << truthColl.size() << endl;
  //cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
  //for(int i=2; i<truthColl.size(); i++){
  //  cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << endl;
  //}

  //==== get reco info here
  std::vector<snu::KMuon> muontriTightColl = GetMuons(BaseSelection::MUON_HN_TRI_TIGHT);
  snu::KParticle reco_lep[3];
  for(int i=0;i<3;i++){
    reco_lep[i] = muontriTightColl.at(i);
  }
  int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
  if(reco_lep[0].Charge() * reco_lep[1].Charge() > 0){ // Q(0) = Q(1)
    if(reco_lep[1].Charge() * reco_lep[2].Charge() < 0){ // Q(1) != Q(2)
      OppSign = 2;
      SameSign[0] = 0;
      SameSign[1] = 1;
    }
    else return; // veto Q(0) = Q(1) = Q(2)
  }
  else{ // Q(0) != Q(1)
    if(reco_lep[0].Charge() * reco_lep[2].Charge() > 0){ // Q(0) = Q(2)
      OppSign = 1;
      SameSign[0] = 0;
      SameSign[1] = 2;
    }
    else if(reco_lep[1].Charge() * reco_lep[2].Charge() > 0){ // Q(1) = Q(2)
      OppSign = 0;
      SameSign[0] = 1;
      SameSign[1] = 2;
    }
  } // Find l2 and assign l1&l3 in ptorder 

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
  snu::KParticle reco_MET;
  reco_MET.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);

  //==== W_pri : first W
  //==== l_1 : lepton from first W
  //==== l_2 : lepton from HN
  //==== W_sec : W from HN
  //==== l_3 : lepton from second W
  
  int truthmax = truthColl.size();
  vector<int> gen_HN_indices, gen_W_pri_indices, gen_W_sec_indices, gen_l_1_indices, gen_nu_indices, gen_l_3_indices, gen_l_2_indices;
  snu::KParticle gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;
  bool W_sec_in_truth=false, isLowMass = false;

  //==== check if this is low/high mass region
  if(k_sample_name.Contains("HN40_") || k_sample_name.Contains("HN50_") || k_sample_name.Contains("HN60_")) isLowMass = true;
 
  //==== find HN index
  for(int i=2;i<truthmax;i++){
    if(truthColl.at(i).PdgId() == 80000002){
      gen_HN_indices.push_back(i);
      find_decay(truthColl, i, gen_HN_indices);
      break;
    }
  }
  //print_all_indices("gen_HN", gen_HN_indices);
  if(gen_HN_indices.size() == 0) return;

  //======================
  //==== low mass region
  //======================
  if(isLowMass){
    //==== for low mass, W_pri is on-shell
    //==== so we can find them in gen particle collections
    //==== find W_pri index
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 24){
        gen_W_pri_indices.push_back(i);
        find_decay(truthColl, i, gen_W_pri_indices);
        break;
      }
    }
    //print_all_indices("gen_W_pri", gen_W_pri_indices);
    if(gen_W_pri_indices.size() == 0) return;

    //==== find l_1 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_pri_indices.size();j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = W_pri
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_pri_indices.at(j) ){
          gen_l_1_indices.push_back(i);
          find_decay(truthColl, i, gen_l_1_indices);
          break;
        }
      }
      if(gen_l_1_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_1", gen_l_1_indices);
    if(gen_l_1_indices.size() == 0) return;

    //==== As m(HN) goes closer to W mass, W_sec starting to appear in truth coll
    //==== if W_sec is virtual, it may not appear in the gen particle collections
    //==== check W_sec exists
    for(int i=2;i<truthmax;i++){
      //==== 1) |PID| = 24
      //==== 2) Mother's PID = 80000002
      if(fabs(truthColl.at(i).PdgId()) == 24 && truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 80000002){
        gen_W_sec_indices.push_back(i);
        find_decay(truthColl, i, gen_W_sec_indices);
        W_sec_in_truth = true;
        break;
      }
    }
    //print_all_indices("gen_W_sec", gen_W_sec_indices);

    //==== no W_sec in truthColl index case
    if(!W_sec_in_truth){
      //cout << "[No W_sec!]" << endl;

      //==== find nu at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = 14
          //==== 2) Mother = {HN}
          if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_nu_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_nu_indices.size() != 0) break;
      }
      //print_all_indices("gen_nu", gen_nu_indices);
      if(gen_nu_indices.size() == 0) return;

      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = PID of l_1 (SS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_l_3_indices);
            break;
          }
        }
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0) return;

      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = - (PID of l_1) (OS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
        if(gen_l_2_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0) return;
    }
    //==== W_sec in truthColl index case
    else{
      //cout << "[W_sec!]" << endl;
      //==== find nu at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
          //==== 1) |PID| = 14
          //==== 2) Mother = {W_sec}
          if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
            gen_nu_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_nu_indices.size() != 0) break;
      }
      //print_all_indices("gen_nu", gen_nu_indices);
      if(gen_nu_indices.size() == 0) return;

      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
          //==== 1) |PID| = 13
          //==== 2) Mother = {W_sec}
          if(fabs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_l_3_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0) return;

      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) |PID| = 13
          //==== 2) Mother = {HN}
          if(fabs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
        if(gen_l_2_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0) return;
      
    }

    gen_l_1 = truthColl.at( gen_l_1_indices.back() );
    gen_l_2 = truthColl.at( gen_l_2_indices.back() );
    gen_l_3 = truthColl.at( gen_l_3_indices.back() );
    gen_nu = truthColl.at( gen_nu_indices.back() );

    //==== solution selection
    double pz_sol_lowmass[2];
    pz_sol_lowmass[0] = solveqdeq(80.4, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "m"); // 0 = minus
    pz_sol_lowmass[1] = solveqdeq(80.4, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "p"); // 1 = plus
    if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
      FillHist("GEN_all_found", 0, 1, 0, 1, 1);
      n_gen_pass++;
      int best_sel = fabs(pz_sol_lowmass[0]-gen_nu.Pz()) < fabs(pz_sol_lowmass[1]-gen_nu.Pz()) ? 0 : 1;
      int smaller = fabs(pz_sol_lowmass[0]) < fabs(pz_sol_lowmass[1]) ? 0 : 1;
      int larger = smaller == 0 ? 1 : 0;
      //cout
      //<< "gen_nu.Pz() = " << gen_nu.Pz() << endl
      //<< "pz_sol_lowmass[0] = " << pz_sol_lowmass[0] << endl
      //<< "pz_sol_lowmass[1] = " << pz_sol_lowmass[1] << endl;
      FillHist("GEN_chi2_best", pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2), 1, 0, 10000, 10000);
      FillHist("GEN_chi2_plus", pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2), 1, 0, 10000, 10000);
      FillHist("GEN_chi2_minus", pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2), 1, 0, 10000, 10000);
      FillHist("GEN_chi2_smaller", pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2), 1, 0, 10000, 10000);
      FillHist("GEN_chi2_larger", pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2), 1, 0, 10000, 10000);

      if(best_sel == 0) FillHist("GEN_solsel_minus_0_plus_1", 0, 1, 0, 2, 2);
      else              FillHist("GEN_solsel_minus_0_plus_1", 1, 1, 0, 2, 2);
      if(best_sel == smaller) FillHist("GEN_solsel_smaller_0_larger_1", 0, 1, 0, 2, 2);
      else                    FillHist("GEN_solsel_smaller_0_larger_1", 1, 1, 0, 2, 2);

      sol_sel_chi2_best += pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2);
      sol_sel_chi2_plus += pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2);
      sol_sel_chi2_minus += pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2);
      sol_sel_chi2_smaller += pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2);
      sol_sel_chi2_larger += pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2);
    }


  }

  //=======================
  //==== high mass region
  //=======================
  else{
    //==== find l_1 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<1;j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = Mother of HN.at(0) : becase they are generated at the same time.
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == truthColl.at(gen_HN_indices.at(j)).IndexMother()){
          gen_l_1_indices.push_back(i);
          find_decay(truthColl, i, gen_l_1_indices);
          break;
        }
      }
      if(gen_l_1_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_1", gen_l_1_indices);
    if(gen_l_1_indices.size() == 0) return;

    //==== fine l_2 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_HN_indices.size();j++){
      //==== 1) |PID| = 13
      //==== 2) Mother = {HN}
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
          gen_l_2_indices.push_back(i);
          find_decay(truthColl, i, gen_l_2_indices);
          break;
        }
      }
      if(gen_l_2_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_2", gen_l_2_indices);    
    if(gen_l_2_indices.size() == 0) return;

    //==== find W_sec
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_HN_indices.size();j++){
        //==== 1) |PID| = 24
        //==== 2) Mother = {HN}
        if(abs(truthColl.at(i).PdgId()) == 24 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
          gen_W_sec_indices.push_back(i);
          find_decay(truthColl, i, gen_W_sec_indices);
          break;
        }
      }
      if(gen_W_sec_indices.size() != 0) break;
    }
    //print_all_indices("gen_W_sec", gen_W_sec_indices);
    if(gen_W_sec_indices.size() == 0) return;

    //==== find nu at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
        //==== 1) |PID| = 14
        //==== 2) Mother = {W_sec}
        if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
          gen_nu_indices.push_back(i);
          find_decay(truthColl, i, gen_nu_indices);
          break;
        }
      }
      if(gen_nu_indices.size() != 0) break;
    }
    //print_all_indices("gen_nu", gen_nu_indices);
    if(gen_nu_indices.size() == 0) return;

    //==== find l_3 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = {W_sec}
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
          gen_l_3_indices.push_back(i);
          find_decay(truthColl, i, gen_l_3_indices);
          break;
        }
      }
      if(gen_l_3_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_3", gen_l_3_indices);
    if(gen_l_3_indices.size() == 0) return;

    gen_l_1 = truthColl.at( gen_l_1_indices.back() );
    gen_l_2 = truthColl.at( gen_l_2_indices.back() );
    gen_l_3 = truthColl.at( gen_l_3_indices.back() );
    gen_nu = truthColl.at( gen_nu_indices.back() );

    snu::KParticle reco_lep_tlv[3];
    for(int i=0; i<3; i++) reco_lep_tlv[i] = reco_lep[i];
    int l_3_cand = find_mlmet_closest_to_W(reco_lep_tlv, reco_MET);

    FillHist("GEN_reco_MET", reco_MET.Pt(), 1, 0, 120, 120);
    FillHist("GEN_reco_lep_1_MET", (reco_lep[0] + reco_MET).M() - 80.4, 1, -60, 60, 120);
    FillHist("GEN_reco_lep_2_MET", (reco_lep[1] + reco_MET).M() - 80.4, 1, -60, 60, 120);
    FillHist("GEN_reco_lep_3_MET", (reco_lep[2] + reco_MET).M() - 80.4, 1, -60, 60, 120);

    FillHist("GEN_l_3_cand", l_3_cand, 1, 0, 3, 3);
    if( gen_l_3.DeltaR(reco_lep[l_3_cand]) < 0.15 ) FillHist("GEN_highmass_mlmet_Wmass_match_gen_l_3", 1, 1, 0, 2, 2);
    else FillHist("GEN_highmass_mlmet_Wmass_match_gen_l_3", 0, 1, 0, 2, 2);

    //==== 1) pt ordering firstly done = m1
    
    int l_1_cand_m1 = SameSign[0], l_SS_rem = SameSign[1], signal_class = 3;
    //==== pt ordering reversed for m(HN) >= 700 GeV
    if( k_sample_name.Contains("HN700_") || k_sample_name.Contains("HN1000_") ){
      signal_class = 4;
      l_1_cand_m1 = SameSign[1];
      l_SS_rem = SameSign[0];
    }
    if( reco_lep[l_1_cand_m1].DeltaR(gen_l_1) < 0.15 ){
      FillHist("GEN_pt_order_first", 1, 1, 0, 2, 2);

      int l_2_cand_m1, l_3_cand_m1;
      if( fabs( (reco_lep[OppSign]+reco_MET).M() - 80.4 ) < fabs( (reco_lep[l_SS_rem]+reco_MET).M() - 80.4 ) ){
        l_3_cand_m1 = OppSign;
        l_2_cand_m1 = l_SS_rem;
      }
      else{
        l_3_cand_m1 = l_SS_rem;
        l_2_cand_m1 = OppSign;
      }
      if( gen_l_2.DeltaR( reco_lep[l_2_cand_m1] ) < 0.15 && gen_l_3.DeltaR( reco_lep[l_3_cand_m1] ) < 0.15 ) FillHist("GEN_pt_order_first_mlmet_next", 1, 1, 0, 2, 2);
      else FillHist("GEN_pt_order_first_mlmet_next", 0, 1, 0, 2, 2);

    }
    else{
      FillHist("GEN_pt_order_first", 0, 1, 0, 2, 2);
    }

    //==== 2) mlmet first = m2
    
    if( gen_l_3.DeltaR(reco_lep[l_3_cand]) < 0.15 ){
      FillHist("GEN_mlmet_first", 1, 1, 0, 2, 2);      

      int l_1_cand_m2, l_2_cand_m2;
      if( l_3_cand == OppSign ){
        if( signal_class == 3){
          l_1_cand_m2 = SameSign[0];
          l_2_cand_m2 = SameSign[1];
        }
        else{
          l_1_cand_m2 = SameSign[1];
          l_2_cand_m2 = SameSign[0];
        }
      }
      else{
        l_2_cand_m2 = OppSign;
        if( l_3_cand == SameSign[0] ) l_1_cand_m2 = SameSign[1];
        else l_1_cand_m2 = SameSign[0];
      }

      if( gen_l_1.DeltaR( reco_lep[l_1_cand_m2] ) < 0.15 && gen_l_2.DeltaR( reco_lep[l_2_cand_m2] ) < 0.15 ) FillHist("GEN_mlmet_first_pt_order_next", 1, 1, 0, 2, 2);
      else FillHist("GEN_mlmet_first_pt_order_next", 0, 1, 0, 2, 2);

    }
    else{
      FillHist("GEN_mlmet_first", 0, 1, 0, 2, 2);
    }


  }

  //==== histograms
 
  FillHist("GEN_matching_validation_W_pri", (gen_l_1+gen_l_2+gen_l_3+gen_nu).M(), 1, 0, 1100, 110);
  FillHist("GEN_matching_validation_HN", (gen_l_2+gen_l_3+gen_nu).M(), 1, 0, 1100, 110);
  FillHist("GEN_matching_validation_W_sec", (gen_l_3+gen_nu).M(), 1, 0, 100, 100);

  FillHist("GEN_gen_l_1_Pt", gen_l_1.Pt(), 1, 0, 100, 100);
  FillHist("GEN_gen_l_2_Pt", gen_l_2.Pt(), 1, 0, 100, 100);
  FillHist("GEN_gen_l_3_Pt", gen_l_3.Pt(), 1, 0, 100, 100);
  FillHist("GEN_gen_nu_Pt", gen_nu.Pt(), 1, 0, 100, 100);
 
  snu::KParticle gen_l_SS;
  if( gen_l_1.Charge() == gen_l_2.Charge() ) gen_l_SS = gen_l_2;
  else gen_l_SS = gen_l_3;

  FillHist("GEN_gen_l_SS_Pt", gen_l_SS.Pt(), 1, 0, 100, 100);
  
  //==== check in gen level
  if( gen_l_1.Pt() > gen_l_SS.Pt() ) FillHist("GEN_gen_pri_lep_pt_greater", 1, 1, 0, 2, 2);
  else FillHist("GEN_gen_pri_lep_pt_greater", 0, 1, 0, 2, 2);

  if( reco_lep[SameSign[0]].DeltaR(gen_l_1) < 0.15 
      //&& fabs(reco_lep[SameSign[0]].Pt()-gen_l_1.Pt())/gen_l_1.Pt() < 0.05 
    ) FillHist("GEN_reco_leading_SS_match_gen_l_1", 1, 1, 0, 2, 2);
  else FillHist("GEN_reco_leading_SS_match_gen_l_1", 0, 1, 0, 2, 2); 

  if( reco_lep[SameSign[1]].DeltaR(gen_l_1) < 0.15
      //&& fabs(reco_lep[SameSign[0]].Pt()-gen_l_1.Pt())/gen_l_1.Pt() < 0.05 
    ) FillHist("GEN_reco_subleading_SS_match_gen_l_1", 1, 1, 0, 2, 2);
  else FillHist("GEN_reco_subleading_SS_match_gen_l_1", 0, 1, 0, 2, 2);


}


void trilepton_mumumu::find_decay(std::vector<snu::KTruth> truthcoll, int target_index, std::vector<int>& indices){

  for(unsigned int i=target_index+1; i<truthcoll.size(); i++){ 
    if( truthcoll.at(i).IndexMother() == target_index && truthcoll.at(i).PdgId() == truthcoll.at(target_index).PdgId() ){
      indices.push_back(i);
      find_decay(truthcoll, i, indices);
    }
  }

}

void trilepton_mumumu::print_all_indices(TString particle, std::vector<int> vec){

  cout << particle+" indices" << endl;
  for(unsigned int i=0; i<vec.size(); i++) cout << " " << vec.at(i) << endl;

}

