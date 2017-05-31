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
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_FR_method);

trilepton_mumumu_FR_method::trilepton_mumumu_FR_method() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_FR_method");
  
  Message("In trilepton_mumumu_FR_method constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  MakeCleverHistograms(hntrilephist, "cut0");
  MakeCleverHistograms(hntrilephist, "cut0_up");
  MakeCleverHistograms(hntrilephist, "cut0_down");
  MakeCleverHistograms(hntrilephist, "cutWlow");
  MakeCleverHistograms(hntrilephist, "cutWlow_up");
  MakeCleverHistograms(hntrilephist, "cutWlow_down");
  MakeCleverHistograms(hntrilephist, "cutWhigh");
  MakeCleverHistograms(hntrilephist, "cutWhigh_up");
  MakeCleverHistograms(hntrilephist, "cutWhigh_down");

  int signal_masses[] = {5, 10, 20, 30, 40, 50, 60, 70, 90, 100, 150, 200, 300, 400, 500, 700, 1000};
  for(int i=0; i<17; i++){
    TString thiscut = "cutHN"+TString::Itoa(signal_masses[i],10);
    MakeCleverHistograms(hntrilephist, thiscut);
    MakeCleverHistograms(hntrilephist, thiscut+"_up");
    MakeCleverHistograms(hntrilephist, thiscut+"_down");
  }

  MakeCleverHistograms(hntrilephist, "MuMuE");
  MakeCleverHistograms(hntrilephist, "MuMuE_up");
  MakeCleverHistograms(hntrilephist, "MuMuE_down");


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

  //============================================
  //==== Apply the gen weight (for NLO, +1,-1)
  //============================================

  if(!isData) weight*=MCweight;
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  //===================
  //==== [CUT] No Cut
  //===================

  FillCutFlow("NoCut", 1.);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  //======================
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);

  //====================
  //==== [CUT] Trigger
  //====================

  std::vector<TString> triggerlist;
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  //==== If this is a Cut Optimization study,
  //==== let's not use trigger pass here.
  bool DoCutOp = std::find(k_flags.begin(), k_flags.end(), "cutop") != k_flags.end();
  if(!DoCutOp){
    if(!PassTriggerOR(triggerlist)) return;
    FillCutFlow("TriggerCut", 1.);
    m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  }

  //=======================
  //==== [CUT] Vertex cut
  //=======================

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  //==== Has Good Primary vertex:
  //==== if ( vtx.ndof() > 4 &&
  //====   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //====   ( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //====   !(vtx.isFake() ) )
  FillCutFlow("VertexCut", 1.);

  //=============================================
  //==== Prepare Muons and Rochestor Correction
  //=============================================

  double this_RelIso = 0.4;
  std::vector<snu::KMuon> muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);
  
  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);
  //std::vector<snu::KJet> jetColl_hn = GetJets("JET_NOLEPTONVETO", 25., 5.0);
  int n_jets = jetColl_hn.size();
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      n_bjets++;
      FillHist("bjet_pt", jetColl_hn.at(j).Pt(), 1., 0., 200., 200);
    }
  }

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrontriLooseColl = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE");

  //==================================
  //==== Number of Loose/Tight Muons
  //==================================

  int n_triLoose_muons = muontriLooseColl.size();
  int n_triTight_muons(0);
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    if(PassID(muontriLooseColl.at(i), "MUON_HN_TRI_TIGHT")) n_triTight_muons++;
  }

  int n_triLoose_electrons = electrontriLooseColl.size();
  int n_triTight_electrons(0);
  for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
    if(PassID(electrontriLooseColl.at(i), "ELECTRON_MVA_TIGHT")) n_triTight_electrons++;
  }

  int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
  int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

  //=============================
  //==== [CUT] Three Tight mons
  //=============================

  //==== Three Muons, and no fourth loose lepton
  bool isThreeMuon     = (n_triLoose_leptons == 3)
                         && (n_triLoose_muons == 3 && n_triTight_muons != 3);
  bool isTwoMuonOneElectron = (n_triLoose_leptons == 3)
                              && (n_triLoose_muons == 2)
                              && (n_triLoose_electrons == 1)
                              && (n_triTight_leptons != 3);

  if(!isThreeMuon && !isTwoMuonOneElectron) return;

  double MinLeadingMuonPt = 20;
  if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) return;
  FillCutFlow("3muon", 1.);

  KLepton lep[3];
  snu::KParticle HN[4];
  for(unsigned int i=0;i<n_triLoose_muons;i++){
    lep[i] = muontriLooseColl.at(i);
  }
  for(unsigned int i=0;i<n_triLoose_electrons;i++){
    lep[i+n_triLoose_muons] = electrontriLooseColl.at(i);
  }

  //==== fake method weighting
  bool UsePtCone = std::find(k_flags.begin(), k_flags.end(), "UsePtCone") != k_flags.end();
  m_datadriven_bkg->SetUsePtCone(UsePtCone);
  double this_weight =     m_datadriven_bkg->Get_DataDrivenWeight(false, muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
  double this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");

  int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
  if(isThreeMuon){
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

  }
  else if(isTwoMuonOneElectron){
    //==== We want SS muon.
    //==== If OS, return
    if(muontriLooseColl.at(0).Charge() != muontriLooseColl.at(1).Charge()) return;
    //==== Remaining electron should have opposite sign
    //==== If SS, return
    if(muontriLooseColl.at(0).Charge() == electrontriLooseColl.at(0).Charge()) return;

    OppSign = 2;
    SameSign[0] = 0;
    SameSign[1] = 1;

  }
  else return;

  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), this_weight, 0., 20., 200);
  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), this_weight, 0., 20., 200);
  FillHist("CutStudy_lowssllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), this_weight, 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), this_weight, 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), this_weight, 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), this_weight, 0., 20., 200);

  if(isThreeMuon){
    //==== MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
    //==== MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
    //==== POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
    if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
        (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) return;
    FillCutFlow("mllsf4", 1.);
  }

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET();
  double METphi = Evt.METPhi();
  CorrectedMETRochester(muontriLooseColl, MET, METphi);
  m_HNgenmatch->SetMETInfo(MET, METphi);

  ///////////////////////////////////////////
  ////////// m(HN) < 80 GeV region //////////
  ///////////////////////////////////////////

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
  
  //==== reconstruct HN and W_real 4-vec with selected Pz solution
  PutNuPz(&nu_lowmass, pz_sol_lowmass[solution_selection_lowmass] );
  //==== SameSign[0] : leading among SS
  //==== SameSign[1] : subleading among SS
  //==== [class1]
  //==== m(HN) : 5 ~ 50 GeV - SS_leading is primary
  //==== [class2]
  //==== m(HN) : 60, 70 GeV - SS_subleading is primary

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
  int l_3_index(0);
  if(isThreeMuon) l_3_index = find_mlmet_closest_to_W(lep, nu_highmass);
  else if(isTwoMuonOneElectron) l_3_index = 2;
  else return;
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

    if(isThreeMuon){

    if(l_3_index == OppSign){
      HN[2] = W_sec + lep[SameSign[1]]; // [class3]
      HN[3] = W_sec + lep[SameSign[0]]; // [class4]
    }
    else{
      HN[2] = W_sec + lep[OppSign]; // [class3]
      HN[3] = W_sec + lep[OppSign]; // [class4]
    }

  }
  else if(isTwoMuonOneElectron){
    HN[2] = W_sec + lep[SameSign[1]]; // [class3]
    HN[3] = W_sec + lep[SameSign[0]]; // [class4]
  }
  else return;

  if(isThreeMuon){
    FillHist("CutStudy_m_Z_candidate", z_candidate.M(), this_weight, 0., 1000., 1000);
    bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
    if(!VetoZResonance) return;
    FillCutFlow("ZVeto", 1.);

    FillHist("CutStudy_mlll", (lep[0] + lep[1] + lep[2]).M(), this_weight, 0., 1000., 1000);
    bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;
    if(!mllloffZ) return;
    FillCutFlow("mllloffZ", 1.);
  }

  FillHist("CutStudy_nbjet", n_bjets, this_weight, 0., 10., 10);

  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);

  //==== preselection is done

  if( DoCutOp ){

    //==== HLT_Mu50_v 36810.752
    //==== HLT_IsoTkMu24_v 36810.751
    //==== HLT_IsoMu24_v 36810.751

    std::vector<TString> triggerlist_1, triggerlist_2, triggerlist_3, triggerlist_123;

    triggerlist_1 = triggerlist;
    triggerlist_1.push_back("HLT_Mu50_v");

    triggerlist_2 = triggerlist;
    triggerlist_2.push_back("HLT_IsoTkMu24_v");

    triggerlist_3 = triggerlist;
    triggerlist_3.push_back("HLT_IsoMu24_v");

    triggerlist_123 = triggerlist;
    triggerlist_123.push_back("HLT_Mu50_v");
    triggerlist_123.push_back("HLT_IsoTkMu24_v");
    triggerlist_123.push_back("HLT_IsoMu24_v");

    if(PassTriggerOR(triggerlist)){
      FillHist("TriggerStudy_unweighted", 0., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 0., this_weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_1)){
      FillHist("TriggerStudy_unweighted", 1., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 1., this_weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_2)){
      FillHist("TriggerStudy_unweighted", 2., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 2., this_weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_3)){
      FillHist("TriggerStudy_unweighted", 3., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 3., this_weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_123)){
      FillHist("TriggerStudy_unweighted", 4., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 4., this_weight, 0., 5., 5);
    }

    return;

  }

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);
  SetPlotHNTriLepMetInfo(MET, METphi);
  SetPlotHNTriLepParticleInfo(HN, W_pri_lowmass, W_pri_highmass, W_sec);
  SetPlotHNTriLepChargeSign(OppSign, SameSign[0], SameSign[1]);
  SetPlotHNTriBJet(n_bjets);

  bool isLowMass = (W_pri_lowmass.M() < 150.);
  bool isHighMass = (MET > 20.);

  if(isThreeMuon){

    FillCLHist(hntrilephist, "cut0_up", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight+this_weight_err);
    FillCLHist(hntrilephist, "cut0_down", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight-this_weight_err);
    FillCLHist(hntrilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight);

    if( isLowMass ){
      FillCLHist(hntrilephist, "cutWlow_up", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight+this_weight_err);
      FillCLHist(hntrilephist, "cutWlow_down", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight-this_weight_err);
      FillCLHist(hntrilephist, "cutWlow", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight);
      FillCutFlow("LowMass", 1.);
    }

    if( isHighMass ){
      FillCLHist(hntrilephist, "cutWhigh_up", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight+this_weight_err);
      FillCLHist(hntrilephist, "cutWhigh_down", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight-this_weight_err);
      FillCLHist(hntrilephist, "cutWhigh", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight);
      FillCutFlow("HighMass", 1.);
    }

    int signal_masses[] = {5, 10, 20, 30, 40, 50, 60, 70, 90, 100, 150, 200, 300, 400, 500, 700, 1000};
    for(int i=0; i<17; i++){

      TString thiscut = "cutHN"+TString::Itoa(signal_masses[i],10);
      double this_W_pri_mass = W_pri_lowmass.M();
      if( signal_masses[i] > 80 ) this_W_pri_mass = W_pri_highmass.M();

      double hnmass = -999.;
      if(signal_masses[i] <= 50) hnmass = HN[0].M();
      else if(signal_masses[i] <= 80) hnmass = HN[1].M();
      else if(signal_masses[i] <= 1000) hnmass = HN[2].M();
      else hnmass = HN[3].M();

      bool pass_op = PassOptimizedCut(signal_masses[i],
        muontriLooseColl.at(0).Pt(), muontriLooseColl.at(1).Pt(), muontriLooseColl.at(2).Pt(),
        this_W_pri_mass, hnmass,
        deltaR_OS_min, gamma_star.M(),
        MET
      );

      if(pass_op){
        FillCLHist(hntrilephist, thiscut+"_up", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight+this_weight_err);
        FillCLHist(hntrilephist, thiscut+"_down", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight-this_weight_err);
        FillCLHist(hntrilephist, thiscut, eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight);
      }

    }

  }
  else if(isTwoMuonOneElectron){
    FillCLHist(hntrilephist, "MuMuE_up", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight+this_weight_err);
    FillCLHist(hntrilephist, "MuMuE_down", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight-this_weight_err);
    FillCLHist(hntrilephist, "MuMuE", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, this_weight);  

    if(DoCutOp){
      bool PassIsoMu24 = PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v");

      if(PassIsoMu24){
        FillHist("TEST_MuMuE_IsoMu24", 0., this_weight, 0., 1., 1);
      }
      if(PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_DiMu", 0., this_weight, 0., 1., 1);
      }
      if(PassIsoMu24 || PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_IsoMu24_OR_DiMu", 0., this_weight, 0., 1., 1);
      }
      if(PassIsoMu24 && !PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_IsoMu24_AND_NOT_DiMu", 0., this_weight, 0., 1., 1);

        int n_IsoMu24_but_notTight=0;
        for(unsigned int i=0; i<muontriLooseColl.size(); i++){
          snu::KMuon thismuon = muontriLooseColl.at(i);
          if(thismuon.TriggerMatched("HLT_IsoMu24_v") || thismuon.TriggerMatched("HLT_IsoTkMu24_v")){
            FillHist("TEST_Mu24FiredMuon_RelIso", thismuon.RelIso04(), 1., 0., 1.0, 100);
            if(thismuon.RelIso04()>0.1) n_IsoMu24_but_notTight++;
          }
        }
        FillHist("TEST_n_Mu24FiredMuon", n_IsoMu24_but_notTight, 1., 0., 5., 5);
      }
    }
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
    AnalyzerCore::MakeHistograms("cutflow", 12,0.,12.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3muon");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS"); 
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"mllsf4");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"ZVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"mllloffZ");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"bjetVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"LowMass");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"HighMass");
    
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
















