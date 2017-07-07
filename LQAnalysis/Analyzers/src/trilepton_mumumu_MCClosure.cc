// $Id: trilepton_mumumu_MCClosure.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_MCClosure Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_MCClosure.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_MCClosure);

trilepton_mumumu_MCClosure::trilepton_mumumu_MCClosure() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_MCClosure");
  
  Message("In trilepton_mumumu_MCClosure constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_MCClosure::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_mumumu_MCClosure::ExecuteEvents()throw( LQError ){

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

  if(!PassTriggerOR(triggerlist)) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);

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

  //======================
  //==== Prepare Leptons
  //======================

  double this_RelIso = 0.4;

  std::vector<snu::KMuon> muontriLooseColl= GetHNTriMuonsByLooseRelIso(this_RelIso, true);;
  std::vector<snu::KElectron>  electrontriLooseColl = GetElectrons(false, true, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);

  std::vector<snu::KMuon> muontriLooseColl_prompt = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
  std::vector<snu::KElectron> electrontriLooseColl_prompt = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);

/*
  //====================================
  //==== Return all prompt lepton case
  //====================================

  if( (muontriLooseColl.size()==muontriLooseColl_prompt.size()) &&
      (electrontriLooseColl.size()==electrontriLooseColl_prompt.size()) ) return;
*/

  //=========================== 
  //==== Get Muon Corrections
  //===========================

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muontriLooseColl);

  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);
  int n_jets = jetColl_hn.size();
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      n_bjets++;
    }
  }

  //===========================
  //==== Trigger Scale Factor
  //===========================

  double trigger_sf = 1.;
  if(!k_isdata){
    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, "ELECTRON_HN_TIGHTv4", muontriLooseColl, "MUON_HN_TRI_TIGHT", 0, 0, 0);
    double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, "ELECTRON_HN_TIGHTv4", muontriLooseColl, "MUON_HN_TRI_TIGHT", 0, 1, 0);
    trigger_sf = trigger_eff_Data/trigger_eff_MC;
  }

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight=(1.0);
  if(!k_isdata){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  //===============================
  //==== Get Electron Corrections
  //===============================

  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrontriLooseColl, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrontriLooseColl);

  //========================
  //==== Apply corrections
  //========================

  weight*=muon_id_iso_sf;
  weight*=trigger_sf;
  weight*=MuTrkEffSF;
  weight*=trigger_ps_weight;
  weight*=pileup_reweight;
  weight*=GetKFactor();
  weight*=electron_sf;
  weight*=electron_RecoSF;

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

  FillHist("n_loose_muon", n_triLoose_muons, 1., 0., 10., 10);
  FillHist("n_tight_muon", n_triTight_muons, 1., 0., 10., 10);

  //=============================
  //==== [CUT] Three Tight mons
  //=============================

  //==== Three Tight Muons, and no fourth loose lepton
  bool isThreeMuon_Measured          = (n_triLoose_leptons == 3)
                                       && (n_triLoose_muons == 3 && n_triTight_muons == 3);
  bool isTwoMuonOneElectron_Measured = (n_triLoose_leptons == 3)
                                       && (n_triLoose_muons == 2 && n_triTight_muons == 2)
                                       && (n_triLoose_electrons == 1 && n_triTight_electrons == 1);

  bool isThreeMuon_Predicted          = (n_triLoose_leptons == 3)
                                        && (n_triLoose_muons == 3 && n_triTight_muons != 3);
  bool isTwoMuonOneElectron_Predicted = (n_triLoose_leptons == 3)
                                        && (n_triLoose_muons == 2)
                                        && (n_triLoose_electrons == 1)
                                        && (n_triTight_leptons != 3);

  //==== MuMuMu or MuMuE
  bool isThreeMuon = isThreeMuon_Measured || isThreeMuon_Predicted;
  bool isTwoMuonOneElectron = isTwoMuonOneElectron_Measured || isTwoMuonOneElectron_Predicted;

  //==== For measured or for predicted
  bool IsForMeasured  = isThreeMuon_Measured || isTwoMuonOneElectron_Measured;
  bool IsForPredicted = isThreeMuon_Predicted || isTwoMuonOneElectron_Predicted;
  if(!IsForMeasured && !IsForPredicted ) return;

  double MinLeadingMuonPt = 20;
  if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) return;
  FillCutFlow("3muon", 1.);

  std::vector<KLepton> lep;
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    KLepton this_lep( muontriLooseColl.at(i) );
    lep.push_back( this_lep );
  }
  for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
    KLepton this_lep( electrontriLooseColl.at(i) );
    lep.push_back( this_lep );
  }

  snu::KParticle HN[4];

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

    //==== To check OS+e
    //if(muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge()) return;

    OppSign = 2;
    SameSign[0] = 0;
    SameSign[1] = 1;

    
  }
  else return;

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
  if(k_sample_name.Contains("HN") && m_HNgenmatch->allgenfound && isThreeMuon) m_HNgenmatch->solution_selection_study(muontriLooseColl);

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
  
  // reconstruct HN and W_real 4-vec with selected Pz solution
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
  // m(HN) : 90 ~ 1000 GeV - primary lepton has larger pT
  // [class4]
  // m(HN) > 1000 GeV - primary lepton has smaller pT

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

/*
  if(isThreeMuon){
    bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
    if(!VetoZResonance) return;
    FillCutFlow("ZVeto", 1.);

    bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;
    if(!mllloffZ) return;
    FillCutFlow("mllloffZ", 1.);
  }


  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);
*/

  //==== preselection is done


  m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true);
  bool UsePtCone = std::find(k_flags.begin(), k_flags.end(), "UsePtCone") != k_flags.end();
  m_datadriven_bkg->SetUsePtCone(UsePtCone);
  double this_weight =     m_datadriven_bkg->Get_DataDrivenWeight(false, muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
  double this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");

  this_weight *= weight;
  this_weight_err *= weight;
  
  if(isThreeMuon){

    bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
    bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;

    std::map<TString, bool> map_to_AnalysisRegion;
    map_to_AnalysisRegion.clear();
    map_to_AnalysisRegion["ThreeMuon"] = true;
    map_to_AnalysisRegion["ThreeMuon_ZVeto"] = map_to_AnalysisRegion["ThreeMuon"] && VetoZResonance;
    map_to_AnalysisRegion["ThreeMuon_ZVeto_mllloffZ"] = map_to_AnalysisRegion["ThreeMuon_ZVeto"] && mllloffZ;
    map_to_AnalysisRegion["ThreeMuon_ZVeto_mllloffZ_bjetveto"] = map_to_AnalysisRegion["ThreeMuon_ZVeto_mllloffZ"] && (n_bjets==0);

    FillMCClosurePlot(
    map_to_AnalysisRegion,
    IsForMeasured,
    IsForPredicted,
    lep,
    weight,
    this_weight,
    this_weight_err,
    MET, METphi
    );

  }

  if(isTwoMuonOneElectron){

    std::map<TString, bool> map_to_AnalysisRegion;
    map_to_AnalysisRegion.clear();
    map_to_AnalysisRegion["TwoMuonOneElectron"] = true;
    map_to_AnalysisRegion["TwoMuonOneElectron_bjetveto"] = map_to_AnalysisRegion["TwoMuonOneElectron"] && (n_bjets==0);

    FillMCClosurePlot(
    map_to_AnalysisRegion,
    IsForMeasured,
    IsForPredicted,
    lep,
    weight,
    this_weight,
    this_weight_err,
    MET, METphi
    );


  }

   return;
}// End of execute event loop
  


void trilepton_mumumu_MCClosure::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_MCClosure::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_MCClosure::~trilepton_mumumu_MCClosure() {
  
  Message("In trilepton_mumumu_MCClosure Destructor" , INFO);
  
}


void trilepton_mumumu_MCClosure::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_MCClosure::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_MCClosure::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_MCClosureCore::MakeHistograms() to make new hists for your analysis
   **/

}


void trilepton_mumumu_MCClosure::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

void trilepton_mumumu_MCClosure::FillMCClosurePlot(
    std::map<TString, bool> map_to_AnalysisRegion,
    bool IsForMeasured,
    bool IsForPredicted,
    std::vector<KLepton> lep,
    double lumi_weight,
    double this_weight,
    double this_weight_err,
    double MET, double METPhi
  ){

  for(std::map<TString, bool>::iterator it=map_to_AnalysisRegion.begin(); it!=map_to_AnalysisRegion.end(); it++){

    if(!(it->second)) continue;
    TString this_Region = it->first;

    if(IsForMeasured){
      FillHist("Measured_"+this_Region, 0., lumi_weight, 0., 1., 1);
      FillHist("Measured_"+this_Region+"_weight", lumi_weight, 1., -10., 10., 200);

      TString leporder[3] = {"leading", "second", "third"};
      for(unsigned int i=0; i<3; i++){
        FillHist("Measured_"+this_Region+"_"+leporder[i]+"Lepton_Pt", lep[i].Pt(), lumi_weight, 0., 1000., 1000);
        FillHist("Measured_"+this_Region+"_"+leporder[i]+"Lepton_Eta", lep[i].Eta(), lumi_weight, -3., 3., 30);
      }

      FillHist("Measured_"+this_Region+"_PFMET", MET, lumi_weight, 0., 500., 500);
      FillHist("Measured_"+this_Region+"_PFMET_phi", METPhi, lumi_weight, -3.2, 3.2, 100);

    }
    else if(IsForPredicted){

      FillHist("Predicted_"+this_Region, 0., this_weight, 0., 1., 1);
      FillHist("Predicted_"+this_Region+"_up", 0., this_weight+this_weight_err, 0., 1., 1);
      FillHist("Predicted_"+this_Region+"_down", 0., this_weight-this_weight_err, 0., 1., 1);

      TString leporder[3] = {"leading", "second", "third"};
      for(unsigned int i=0; i<3; i++){
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Pt", lep[i].Pt(), this_weight, 0., 1000., 1000);
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Eta", lep[i].Eta(), this_weight, -3., 3., 30);
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Pt_up", lep[i].Pt(), this_weight+this_weight_err, 0., 1000., 1000);
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Eta_up", lep[i].Eta(), this_weight+this_weight_err, -3., 3., 30);
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Pt_down", lep[i].Pt(), this_weight-this_weight_err, 0., 1000., 1000);
        FillHist("Predicted_"+this_Region+"_"+leporder[i]+"Lepton_Eta_down", lep[i].Eta(), this_weight-this_weight_err, -3., 3., 30);
      }

      FillHist("Predicted_"+this_Region+"_PFMET", MET, this_weight, 0., 500., 500);
      FillHist("Predicted_"+this_Region+"_PFMET_up", MET, this_weight+this_weight_err, 0., 500., 500);
      FillHist("Predicted_"+this_Region+"_PFMET_down", MET, this_weight-this_weight_err, 0., 500., 500);

      FillHist("Predicted_"+this_Region+"_PFMET_phi", METPhi, this_weight, -3.2, 3.2, 100);
      FillHist("Predicted_"+this_Region+"_PFMET_phi_up", METPhi, this_weight+this_weight_err, -3.2, 3.2, 100);
      FillHist("Predicted_"+this_Region+"_PFMET_phi_down", METPhi, this_weight-this_weight_err, -3.2, 3.2, 100);

    }

  }




}



















