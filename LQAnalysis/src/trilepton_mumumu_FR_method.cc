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

  this_dXYSig = 4.0;
  this_RelIso = 0.4;

  return;

}


void trilepton_mumumu_FR_method::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", 1.);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  bool DoCutOp = std::find(k_flags.begin(), k_flags.end(), "cutop") != k_flags.end();

  std::vector<TString> triggerlist;
  //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  double MinLeadingMuonPt = 20;
  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  bool trigger_pass = false;
  for(unsigned int i=0; i<triggerlist.size(); i++){
    if( PassTrigger(triggerlist.at(i)) ){
      trigger_pass = true;
      break;
    }
  }

  if(!trigger_pass) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", 1.);


  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  //std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  //std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  //std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  //std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed
  std::vector<snu::KMuon> muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  
  //CorrectMuonMomentum(muontriLooseColl); //FIXME do this for v8-0-4
  //float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonTightColl,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", true, 30., 2.4);
  
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  int n_bjets=0;
  for(int j=0; j<njet; j++){
    if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)) n_bjets++;
  }

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

  int n_triTight_muons(0);
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    if(muontriLooseColl.at(i).RelIso04() < 0.1) n_triTight_muons++;
  }
  int n_triLoose_muons = muontriLooseColl.size();

  if( n_triLoose_muons != 3 ) return;

  if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) return;
  FillCutFlow("3muon", 1.);

  snu::KParticle lep[3], HN[4];
  //for(int i=0;i<k_flags.size();i++) cout << "k_flags = " << k_flags.at(i) << endl;
  for(int i=0;i<3;i++){
    lep[i] = muontriLooseColl.at(i);
  }

  //==== fake method weighting
  if( n_triTight_muons == 3 ) return; // return TTT case
  std::vector<snu::KElectron> empty_electron;
  empty_electron.clear();
  weight = Get_DataDrivenWeight(false, muontriLooseColl, empty_electron, 3, 0);
  double weight_err = Get_DataDrivenWeight(true, muontriLooseColl, empty_electron, 3, 0);

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
  FillCutFlow("2SS1OS", 1.);

  FillHist("lowosllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("lowosllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowssllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);

  //==== MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
  //==== MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
  //==== POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
  if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
      (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) return;
  FillCutFlow("mllsf4", 1.);

  ///////////////////////////////////////////
  ////////// m(HN) < 80 GeV region //////////
  ///////////////////////////////////////////

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
  snu::KParticle W_pri_lowmass, nu_lowmass, gamma_star, z_candidate;
  nu_lowmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
  double pz_sol_lowmass[2];
  pz_sol_lowmass[0] = solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "m"); // 0 = minus
  pz_sol_lowmass[1] = solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 1 = plus
  //PutNuPz(&selection_nu[0], solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
  //PutNuPz(&selection_nu[1], solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus

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
  pz_sol_highmass[0] = solveqdeq(80.385, lep[l_3_index], MET, METphi, "m"); // 0 = minus
  pz_sol_highmass[1] = solveqdeq(80.385, lep[l_3_index], MET, METphi, "p"); // 1 = plus
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
  // m(HN) : 90 ~ 700 GeV - primary lepton has larger pT
  // [class4]
  // m(HN) : 1000 GeV - primary lepton has smaller pT

  W_sec = lep[l_3_index] + nu_highmass;

  if(l_3_index == OppSign){
     
    HN[2] = W_sec + lep[SameSign[1]]; // [class3]
    HN[3] = W_sec + lep[SameSign[0]]; // [class4]
 
  }
  else{
      HN[2] = W_sec + lep[OppSign]; // [class3]
      HN[3] = W_sec + lep[OppSign]; // [class4]
  }

  bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
  if(!VetoZResonance) return;
  FillCutFlow("ZVeto", 1.);

  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);

  //==== preselection is done

  if( DoCutOp ){

    double cutop[100];
    cutop[0] = muontriLooseColl.at(0).Pt();
    cutop[1] = muontriLooseColl.at(1).Pt();
    cutop[2] = muontriLooseColl.at(2).Pt();
    cutop[3] = deltaR_OS_min;
    cutop[4] = HN[0].M();
    cutop[5] = HN[1].M();
    cutop[6] = HN[2].M();
    cutop[7] = HN[3].M();
    cutop[8] = W_pri_lowmass.M();
    cutop[9] = W_pri_highmass.M();
    cutop[10] = weight;
    cutop[11] = muontriLooseColl.at(0).dXY();
    cutop[12] = muontriLooseColl.at(1).dXY();
    cutop[13] = muontriLooseColl.at(2).dXY();
    cutop[14] = muontriLooseColl.at(0).dZ();
    cutop[15] = muontriLooseColl.at(1).dZ();
    cutop[16] = muontriLooseColl.at(2).dZ();
    cutop[17] = muontriLooseColl.at(0).RelIso04();
    cutop[18] = muontriLooseColl.at(1).RelIso04();
    cutop[19] = muontriLooseColl.at(2).RelIso04();
    cutop[20] = W_sec.M();
    cutop[21] = MET;
    cutop[22] = weight_err;

    FillNtp("cutop",cutop);
    return;
  }

  bool isLowMass = (W_pri_lowmass.M() < 150.);
  bool isHighMass = (MET > 20.);

  FillUpDownHist("HN_mass_class1_cut0", HN[0].M(), weight, weight_err, 0., 2000., 2000);
  FillUpDownHist("HN_mass_class2_cut0", HN[1].M(), weight, weight_err, 0., 2000., 2000);
  FillUpDownHist("HN_mass_class3_cut0", HN[2].M(), weight, weight_err, 0, 2000, 2000);
  FillUpDownHist("HN_mass_class4_cut0", HN[3].M(), weight, weight_err, 0, 2000, 2000);
  FillUpDownHist("W_pri_lowmass_mass_cut0", W_pri_lowmass.M(), weight, weight_err, 0., 2000., 2000);
  FillUpDownHist("W_pri_highmass_mass_cut0", W_pri_highmass.M(), weight, weight_err, 0, 2000, 2000);
  FillUpDownHist("W_sec_highmass_mass_cut0", W_sec.M(), weight, weight_err, 0., 2000., 2000);
  FillUpDownHist("deltaR_OS_min_cut0", deltaR_OS_min, weight, weight_err, 0, 5, 50);
  FillUpDownHist("gamma_star_mass_cut0", gamma_star.M(), weight, weight_err, 0., 200., 200);
  FillUpDownHist("z_candidate_mass_cut0", z_candidate.M(), weight, weight_err, 0., 200., 200);
  FillUpDownHist("n_events_cut0", 0, weight, weight_err, 0, 1, 1);
  FillCLHist(hntrilephist, "cut0_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
  FillCLHist(hntrilephist, "cut0_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
  FillCLHist(hntrilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);

  if( isLowMass ){
    FillUpDownHist("HN_mass_class1_cutWlow", HN[0].M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("HN_mass_class2_cutWlow", HN[1].M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("HN_mass_class3_cutWlow", HN[2].M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("HN_mass_class4_cutWlow", HN[3].M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("W_pri_lowmass_mass_cutWlow", W_pri_lowmass.M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("W_pri_highmass_mass_cutWlow", W_pri_highmass.M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("W_sec_highmass_mass_cutWlow", W_sec.M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("deltaR_OS_min_cutWlow", deltaR_OS_min, weight, weight_err, 0, 5, 50);
    FillUpDownHist("gamma_star_mass_cutWlow", gamma_star.M(), weight, weight_err, 0., 200., 200);
    FillUpDownHist("z_candidate_mass_cutWlow", z_candidate.M(), weight, weight_err, 0., 200., 200);
    FillUpDownHist("n_events_cutWlow", 0, weight, weight_err, 0, 1, 1);
    FillCLHist(hntrilephist, "cutWlow_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
    FillCLHist(hntrilephist, "cutWlow_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
    FillCLHist(hntrilephist, "cutWlow", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);

    FillCutFlow("LowMass", 1.);

  }

  if( isHighMass ){
    FillUpDownHist("HN_mass_class1_cutWhigh", HN[0].M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("HN_mass_class2_cutWhigh", HN[1].M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("HN_mass_class3_cutWhigh", HN[2].M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("HN_mass_class4_cutWhigh", HN[3].M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("W_pri_lowmass_mass_cutWhigh", W_pri_lowmass.M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("W_pri_highmass_mass_cutWhigh", W_pri_highmass.M(), weight, weight_err, 0, 2000, 2000);
    FillUpDownHist("W_sec_highmass_mass_cutWhigh", W_sec.M(), weight, weight_err, 0., 2000., 2000);
    FillUpDownHist("deltaR_OS_min_cutWhigh", deltaR_OS_min, weight, weight_err, 0, 5, 50);
    FillUpDownHist("gamma_star_mass_cutWhigh", gamma_star.M(), weight, weight_err, 0., 200., 200);
    FillUpDownHist("z_candidate_mass_cutWhigh", z_candidate.M(), weight, weight_err, 0., 200., 200);
    FillUpDownHist("n_events_cutWhigh", 0, weight, weight_err, 0, 1, 1);
    FillCLHist(hntrilephist, "cutWhigh_up", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight+weight_err);
    FillCLHist(hntrilephist, "cutWhigh_down", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight-weight_err);
    FillCLHist(hntrilephist, "cutWhigh", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);

    FillCutFlow("HighMass", 1.);

  }


  return;

}// End of execute event loop
  


void trilepton_mumumu_FR_method::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_FR_method::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
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
  
}


void trilepton_mumumu_FR_method::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 11,0.,11.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3muon");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS"); 
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"mllsf4");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"ZVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"bjetVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"LowMass");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"HighMass");
    
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

  MakeNtp("cutop", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:first_dXY:second_dXY:third_dXY:first_dZ:second_dZ:third_dZ:first_RelIso:second_RelIso:third_RelIso:W_sec_highmass_mass:PFMET:weight_err");

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

