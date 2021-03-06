// $Id: trilepton_mumumu_CR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_CR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_CR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_CR);

trilepton_mumumu_CR::trilepton_mumumu_CR() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_CR");
  
  Message("In trilepton_mumumu_CR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_CR::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_mumumu_CR::ExecuteEvents()throw( LQError ){

/*
  std::vector<snu::KMuon> testmuon         = GetHNTriMuonsByLooseRelIso(0.4, true);
  std::vector<snu::KElectron> testelectron = GetElectrons(false, true, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);

  for(unsigned int i=0; i<testelectron.size(); i++){
    FillHist("TEST_Electron_F0", testelectron.at(i).GetType(), 1., 0., 100., 100);
    if(PassID(testelectron.at(i), "ELECTRON_MVA_TIGHT")){
      FillHist("TEST_Electron_F", testelectron.at(i).GetType(), 1., 0., 100., 100);
    }
    if( testelectron.at(i).MCIsExternalConversion() && testelectron.at(i).MCMatched()){
      FillHist("Conv_FR", 0., 1., 0., 2., 2);
      if(PassID(testelectron.at(i), "ELECTRON_MVA_TIGHT")){
        FillHist("Conv_FR", 1., 1., 0., 2., 2);
      }
    }
  }
  return;

  if(testmuon.size()==2 && testelectron.size()==1){
    snu::KMuon test_mu[2];
    test_mu[0] = testmuon.at(0);
    test_mu[1] = testmuon.at(1);
    snu::KElectron test_el = testelectron.at(0);
    snu::KParticle Zcand = test_mu[0]+test_mu[1]+test_el;
    if( fabs( 91.1876-Zcand.M() ) > 10. ) return;
    if( testelectron.at(0).GetType() != 36 ) return;
    cout << "RunNumber = " << eventbase->GetEvent().RunNumber() << endl;
    cout << "EventNumber = " << eventbase->GetEvent().EventNumber() << endl;
    PrintTruth();
    cout << "#### Reco Muon ####" << endl;
    for(unsigned int i=0;i<testmuon.size();i++){
      cout << testmuon.at(i).Pt() << "\t" << testmuon.at(i).Eta() << "\t" << testmuon.at(i).GetType() << "\t" << TruthMatched(testmuon.at(i)) << endl;
    }
    cout << "#### Reco Electron ####" << endl;
    for(unsigned int i=0;i<testelectron.size();i++){
      cout << testelectron.at(i).Pt() << "\t" << testelectron.at(i).Eta() << "\t" << testelectron.at(i).GetType() << "\t" << TruthMatched(testelectron.at(i),false) << endl;
    }
  }
  return;
*/

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

  std::vector<TString> dimutrigger, dieltrigger;

  dimutrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  //dimutrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  //dimutrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  dimutrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  dieltrigger.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");

  float trigger_ps_weight= WeightByTrigger(dimutrigger, TargetLumi);

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


  //====================
  //==== Prepare Muons
  //====================

  std::vector<snu::KMuon> muontriLooseColl;
  double this_RelIso = 0.4;

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();
  bool DoKeepFake = std::find(k_flags.begin(), k_flags.end(), "DoKeepFake") != k_flags.end();

  //==== Gen Matching is not done correctly for my private samples..
  if( k_sample_name.Contains("HN") ){
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  }
  //==== For MC Closure test, we want to reject P(rompt)P event
  else if( DoMCClosure ){
    std::vector<snu::KMuon> muontriLooseColl_prompt = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
    if(muontriLooseColl_prompt.size()==2) return;

    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  }
  //==== non-prompt : keep fake
  else if( DoKeepFake ){
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  }
  //==== For Prompt MC, collect prompt muons
  //==== For data, it automarically collect all muons :)
  else{
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
  }

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);

  //=========================== 
  //==== Get Muon Corrections
  //===========================

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muontriLooseColl);

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrontriLooseColl;
  if(DoMCClosure){
    std::vector<snu::KElectron> electrontriLooseColl_prompt = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);
    if(electrontriLooseColl_prompt.size()==2) return;

    electrontriLooseColl = GetElectrons(false, true, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);
  }
  else{
    electrontriLooseColl = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);
  }

  //===============================
  //==== Get Electron Corrections
  //===============================

  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrontriLooseColl, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrontriLooseColl);

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

  std::vector<snu::KJet> jetColl_hn_nolepveto = GetJets("JET_NOLEPTONVETO", 25., 2.4);
  std::vector<snu::KJet> jetColl_hn_nearby;
  for(unsigned int i=0; i<jetColl_hn_nolepveto.size(); i++){
    bool isNearByJet = false;
    for(unsigned int j=0; j<muontriLooseColl.size(); j++){
      if(jetColl_hn_nolepveto.at(i).DeltaR( muontriLooseColl.at(j) ) < 0.4){
        isNearByJet = true;
        break;
      }
    }
    if(isNearByJet) jetColl_hn_nearby.push_back( jetColl_hn_nolepveto.at(i) );
  }
  int n_jets_nearby = jetColl_hn_nearby.size();
  int n_bjets_nearby=0;
  for(int j=0; j<n_jets_nearby; j++){
    if( IsBTagged(jetColl_hn_nearby.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      n_bjets_nearby++;
      FillHist("bjet_nearby_pt", jetColl_hn_nearby.at(j).Pt(), 1., 0., 200., 200);
    }
  }

  //===========================
  //==== Trigger Scale Factor
  //===========================

  double trigger_sf = 1.;
  if(!k_isdata){
    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, "ELECTRON_MVA_TIGHT", muontriLooseColl, "MUON_HN_TRI_TIGHT", 0, 0, 0);
    double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, "ELECTRON_MVA_TIGHT", muontriLooseColl, "MUON_HN_TRI_TIGHT", 0, 1, 0);
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

  //========================
  //==== Apply corrections
  //========================

  if(!isData && !k_running_nonprompt){

    if(DoMCClosure){
      if(k_sample_name.Contains("QCD")){
        weight = weight;
      }
      else{
        weight = 1.*MCweight;
      }
    }
    else{
      weight*=muon_id_iso_sf;
      weight*=trigger_sf;
      weight*=trigger_ps_weight;
      weight*=pileup_reweight;
      weight*=MuTrkEffSF;
      weight*=electron_sf;
      weight*=electron_RecoSF;
      weight*=GetKFactor();
    }

  }

  //============================================
  //==== Number of Loose/Tight Muons/Electrons
  //============================================

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

  //===========================
  //==== CR related variables
  //===========================

  FillHist("n_tight_muons_control", n_triTight_muons, weight, 0., 10., 10);
  FillHist("n_loose_muons_control", n_triLoose_muons, weight, 0., 10., 10);
  FillHist("n_tight_electrons_control", n_triTight_electrons, weight, 0., 10., 10);
  FillHist("n_loose_electrons_control", n_triLoose_electrons, weight, 0., 10., 10);
  if( n_triTight_muons == 2 && n_triLoose_muons == 2){
    int isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ? 1 : 0;
    FillHist("2Muons_OS0_SS1_control", isSS, weight, 0., 2., 2);
  }
  FillHist("n_jets_control", n_jets, weight, 0., 10., 10);
  FillHist("n_bjets_control", n_bjets, weight, 0., 10., 10);

  //================
  //==== define CR
  //================

  bool isTwoMuon     = (n_triLoose_leptons == 2)
                       && (n_triLoose_muons == 2 && n_triTight_muons == 2);
  bool isTwoLepton = (n_triLoose_leptons == 2) && (n_triTight_leptons == 2);
  bool isThreeLepton = (n_triLoose_leptons == 3) && (n_triTight_leptons == 3);
  bool isFourLepton  = (n_triLoose_leptons == 4)
                       && (
                         (n_triLoose_muons == 4 && n_triTight_muons == 4) ||
                         (n_triLoose_electrons == 4 && n_triTight_electrons == 4) ||
                         (n_triLoose_muons == 2 && n_triTight_muons == 2 && n_triLoose_electrons == 2 && n_triTight_electrons == 2)
                       );

  if(n_triLoose_muons == 2 && n_triTight_muons == 0) FillHist("LL_TL_TT", 0., 1., 0., 3., 3);
  if(n_triLoose_muons == 2 && n_triTight_muons == 1) FillHist("LL_TL_TT", 1., 1., 0., 3., 3);
  if(n_triLoose_muons == 2 && n_triTight_muons == 2) FillHist("LL_TL_TT", 2., 1., 0., 3., 3);

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET();
  double METphi = Evt.METPhi();
  JSCorrectedMETRochester(muontriLooseColl, MET, METphi);

  double m_Z = 91.1876;

  if(isTwoMuon) FillHist("isTwoMuon", 0., 1., 0., 1., 1);
  if(isTwoLepton&&(n_triLoose_muons==2)) FillHist("isTwoLepton_2muon", 0., 1., 0., 1., 1);

  //==== CR with Two Muons
  if(DoMCClosure && isTwoLepton){

    KLepton lep[2];

    //==== 2 Muon : TwoLeptonConfig = 0;
    //==== 2 Electron : TwoLeptonConfig = 1;

    TString lepconfig="";
    int TwoLeptonConfig = 0;
    double m_dimuon(0.);
    if(n_triLoose_muons==2){
      TwoLeptonConfig = 0;

      lep[0] = muontriLooseColl.at(0);
      lep[1] = muontriLooseColl.at(1);

      if(!PassTriggerOR(dimutrigger)) return;
      if(lep[0].Pt() <= 20) return;

      m_dimuon = (muontriLooseColl.at(0)+muontriLooseColl.at(1)).M();

      lepconfig = "DiMuon";
    }
    else if(n_triLoose_electrons==2){
      TwoLeptonConfig = 1;

      lep[0] = electrontriLooseColl.at(0);
      lep[1] = electrontriLooseColl.at(1);

      if(!PassTriggerOR(dieltrigger)) return;
      if(lep[0].Pt() <= 20) return;
      if(lep[1].Pt() <= 15) return;

      m_dimuon = (electrontriLooseColl.at(0)+electrontriLooseColl.at(1)).M();
      lepconfig = "DiElectron";
    }
    else{
      return;
    }

    bool isSS = lep[0].Charge() == lep[1].Charge();

    if(k_sample_name.Contains("DY") && !isSS) return;

    bool ZResonance = fabs(m_dimuon-m_Z) < 10.;

    if(lepconfig=="DiMuon" && !isTwoMuon){
      FillHist("n_triLoose_muons", n_triLoose_muons, 1., 0., 5., 5);
      FillHist("n_triTight_muons", n_triTight_muons, 1., 0., 5., 5);
      FillHist("n_triLoose_electrons", n_triLoose_electrons, 1., 0., 5., 5);
      FillHist("n_triTight_electrons", n_triTight_electrons, 1., 0., 5., 5);
      FillHist("n_triLoose_leptons", n_triLoose_leptons, 1., 0., 5., 5);
      FillHist("n_triTight_leptons", n_triTight_leptons, 1., 0., 5., 5);
    }

    std::map< TString, bool > map_whichCR_to_isCR;
    map_whichCR_to_isCR.clear();
    map_whichCR_to_isCR[lepconfig] = true;
    map_whichCR_to_isCR["SS"+lepconfig] = isSS;
    map_whichCR_to_isCR["OS"+lepconfig] = !isSS;
    map_whichCR_to_isCR["OS"+lepconfig+"_Z_10GeV"] = !isSS && ZResonance;

    for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
      TString this_suffix = it->first;
      if(it->second){
        FillHist("weight_"+this_suffix, weight, 1., -1., 1., 1000);
        FillHist("n_events_"+this_suffix, 0, weight, 0., 1., 1);
        FillHist("n_jets_"+this_suffix, n_jets, weight, 0., 10., 10);
        FillHist("n_bjets_"+this_suffix, n_bjets, weight, 0., 10., 10);
        FillHist("PFMET_"+this_suffix, MET, weight, 0., 500., 500);
        FillHist("PFMET_phi_"+this_suffix, METphi, weight, -3.2, 3.2, 64);
        FillHist("mll_"+this_suffix, m_dimuon , weight, 0., 500., 500);
        FillHist("n_vertices_"+this_suffix, numberVertices, weight, 0., 50., 50);
        FillHist("leadingLepton_Pt_"+this_suffix, lep[0].Pt() , weight, 0., 200., 200);
        FillHist("leadingLepton_Eta_"+this_suffix, lep[0].Eta() , weight, -3., 3., 60);
        FillHist("leadingLepton_RelIso_"+this_suffix, lep[0].RelIso() , weight, 0., 1.0, 100);
        //FillHist("leadingLepton_Chi2_"+this_suffix, lep[0].GlobalChi2() , weight, 0., 10., 100);
        FillHist("leadingLepton_dXY_"+this_suffix+"", fabs(lep[0].dXY()) , weight, 0., 0.1, 100);
        FillHist("leadingLepton_dXYSig_"+this_suffix+"", fabs(lep[0].dXYSig()) , weight, 0., 4., 40);
        FillHist("secondLepton_Pt_"+this_suffix, lep[1].Pt() , weight, 0., 200., 200);
        FillHist("secondLepton_Eta_"+this_suffix, lep[1].Eta() , weight, -3., 3., 60);
        FillHist("secondLepton_RelIso_"+this_suffix, lep[1].RelIso() , weight, 0., 1.0, 100);
        //FillHist("secondLepton_Chi2_"+this_suffix, lep[1].GlobalChi2() , weight, 0., 10., 100);
        FillHist("secondLepton_dXY_"+this_suffix+"", fabs(lep[1].dXY()) , weight, 0., 0.1, 100);
        FillHist("secondLepton_dXYSig_"+this_suffix+"", fabs(lep[1].dXYSig()) , weight, 0., 4., 40);
      }
    } 

    return;

  } // MCClosure

  //==== CR with Three Muons
  if(!DoMCClosure && isThreeLepton){

    if(!PassTriggerOR(dimutrigger)) return;
    FillCutFlow("TriggerCut", 1.);
    m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

    FillCutFlow("3muon", 1.);
    FillHist("TEST_trigger_sf", trigger_sf, 1., 0., 2., 2000);

    std::vector<KLepton> lep;
    TString lepOrder[3] = {"leading", "second", "third"};
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      KLepton this_lep( muontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      KLepton this_lep( electrontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }

    bool AllSameCharge(false);
    //==== 3 Muon : ThreeLeptonConfig = 0;
    //==== 2 Muon + 1 Electron : ThreeLeptonConfig = 1;
    //==== 1 Muon + 2 Electron ; ThreeLeptonConfig = 2;
    //==== 3 Electron : ThreeLeptonConfig = 3;
    int ThreeLeptonConfig = 0;

    if(muontriLooseColl.size()==3 || electrontriLooseColl.size()==3){

      if(muontriLooseColl.size()==3)     ThreeLeptonConfig = 0;
      if(electrontriLooseColl.size()==3) ThreeLeptonConfig = 3;

      if( ( lep.at(0).Charge() == lep.at(1).Charge() ) &&
          ( lep.at(0).Charge() == lep.at(2).Charge() ) ) AllSameCharge = true;
    }
    else if(muontriLooseColl.size()==2){

      ThreeLeptonConfig = 1;

      if( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) AllSameCharge = true;
    }
    else if(electrontriLooseColl.size()==2){

      ThreeLeptonConfig = 2;

      if( electrontriLooseColl.at(0).Charge() == electrontriLooseColl.at(1).Charge() ) AllSameCharge = true;
    }
    else{
      Message("?", INFO);
    }

    FillHist("ThreeLeptonConfig", ThreeLeptonConfig, 1., 0., 10., 10);

    if( !AllSameCharge ){

      snu::KParticle Z_candidate;
      KLepton ZLepton_leading, ZLepton_subleading, WLepton;
      double m_OSSF[2];

      if(ThreeLeptonConfig == 0 || ThreeLeptonConfig == 3){
        KLepton OS, SS[2];
        if     ( lep.at(0).Charge() == lep.at(1).Charge() ){
          SS[0] = lep.at(0);
          SS[1] = lep.at(1);
          OS    = lep.at(2);
        }
        else if( lep.at(0).Charge() == lep.at(2).Charge() ){
          SS[0] = lep.at(0);
          SS[1] = lep.at(2);
          OS    = lep.at(1);
        }
        else if( lep.at(1).Charge() == lep.at(2).Charge() ){
          SS[0] = lep.at(1);
          SS[1] = lep.at(2);
          OS    = lep.at(0);
        }
        else Message("?", INFO);

        m_OSSF[0] = ( SS[0] + OS ).M();
        m_OSSF[1] = ( SS[1] + OS ).M();

        KLepton SS_ZMuon;
        if( fabs(m_OSSF[0]-m_Z) < fabs(m_OSSF[1]-m_Z) ){
          Z_candidate = SS[0] + OS;
          SS_ZMuon = SS[0];
          WLepton = SS[1];
        }
        else{
          Z_candidate = SS[1] + OS;
          SS_ZMuon = SS[1];
          WLepton = SS[0];
        }

        if( SS_ZMuon.Pt() > OS.Pt() ){
          ZLepton_leading = SS_ZMuon;
          ZLepton_subleading = OS;
        }
        else{
          ZLepton_leading = OS;
          ZLepton_subleading = SS_ZMuon;
        }

      }
      else if(ThreeLeptonConfig == 1){
        Z_candidate = muontriLooseColl.at(0)+muontriLooseColl.at(1);
        ZLepton_leading = muontriLooseColl.at(0);
        ZLepton_subleading = muontriLooseColl.at(1);
        WLepton = electrontriLooseColl.at(0);

        m_OSSF[0] = Z_candidate.M();
        m_OSSF[1] = Z_candidate.M();

        FillHist("2Mu1El_m_Z", Z_candidate.M(), 1., 0., 150., 150);
        FillHist("2Mu1El_ZLepton_leading_Pt", ZLepton_leading.Pt(), 1., 0., 150., 150);
        FillHist("2Mu1El_WLepton_Pt", WLepton.Pt(), 1., 0., 150., 150);
      }
      else if(ThreeLeptonConfig == 2){
        Z_candidate = electrontriLooseColl.at(0)+electrontriLooseColl.at(1);
        ZLepton_leading = electrontriLooseColl.at(0);
        ZLepton_subleading = electrontriLooseColl.at(1);
        WLepton = muontriLooseColl.at(0);

        m_OSSF[0] = Z_candidate.M();
        m_OSSF[1] = Z_candidate.M();
      }

      double mlll = (lep.at(0)+lep.at(1)+lep.at(2)).M();
      bool ZLeptonPtCut = (ZLepton_leading.Pt() > 20.);
      bool WLeptonPtCut = (WLepton.Pt() > 20.);
      bool isZresonance = (fabs(Z_candidate.M()-m_Z) < 10.);
      bool METCut = (MET > 30.);
      bool mlllCut = (mlll > 100.);
      bool mll4 = (m_OSSF[0] < 4.) || (m_OSSF[1] < 4.);
      bool bjetveto = (n_bjets == 0);

      //==== If you don't want to veto b-jet, set it true
      //bjetveto = true;

      FillHist("m_Z_candidate_before_cut_WZ", Z_candidate.M(), weight, 0., 150., 150);
      FillHist("m_lll_before_cut_WZ", mlll, weight, 0., 500., 500);
      FillHist("PFMET_before_cut_WZ", MET, weight, 0., 500., 500);
      FillHist("n_bjets_before_cut_WZ", n_bjets, weight, 0., 10., 10);
      FillHist("n_vertices_before_cut_WZ", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

      //==== N-1 plots
      vector<bool> WZ_cuts;
      WZ_cuts.push_back( fabs(Z_candidate.M()-m_Z) < 10. );
      WZ_cuts.push_back( mlll > 100. );
      WZ_cuts.push_back( bjetveto );
      WZ_cuts.push_back( MET > 30. );
      if( ZLeptonPtCut && WLeptonPtCut && !mll4 ){
        FillHist("N1_preselection_WZ", 0, weight, 0., 1., 1);
        for(unsigned int i=0; i<WZ_cuts.size(); i++){
          bool this_bool = true;
          for(unsigned int j=0; j<WZ_cuts.size(); j++){
            if(j==i) continue;
            if(!WZ_cuts.at(j)) this_bool = false;
          }
          if(this_bool){
            if(i==0) FillHist("N1_Z_mass_WZ", fabs(Z_candidate.M()-m_Z), weight, 0., 60., 60);
            if(i==1) FillHist("N1_mlll_WZ",  mlll, weight, 0., 500., 25);
            if(i==2) FillHist("N1_n_bjets_WZ", n_bjets, weight, 0., 4., 4);
            if(i==3) FillHist("N1_PFMET_WZ", MET, weight, 0., 200., 20);
          }
        }
      }

      snu::KParticle nu;
      nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
      snu::KParticle W_candidate = nu+WLepton;

      std::map< TString, bool > map_whichCR_to_isCR;
      map_whichCR_to_isCR.clear();
      map_whichCR_to_isCR["WZ"]    = ZLeptonPtCut && isZresonance && WLeptonPtCut && METCut      && mlllCut   && !mll4 && bjetveto;
      map_whichCR_to_isCR["ZJets"] = ZLeptonPtCut && isZresonance                 && (MET < 20.) && mlllCut   && !mll4 && bjetveto && MT(nu, WLepton) < 30.;
      map_whichCR_to_isCR["ZLep"]  = ZLeptonPtCut && isZresonance                                && mlllCut   && !mll4 && bjetveto;
      map_whichCR_to_isCR["ZGamma"]= ZLeptonPtCut && (fabs(Z_candidate.M()-m_Z) > 15.) && (MET < 50.) && (fabs(mlll-m_Z) < 10.) && !mll4 && bjetveto;
      map_whichCR_to_isCR["LowMassResonance"]= ZLeptonPtCut && (fabs(Z_candidate.M()-m_Z) > 15.) && (MET < 50.) && !mll4 && bjetveto;

      map_whichCR_to_isCR["WZ_3mu0el"] = map_whichCR_to_isCR["WZ"] && (ThreeLeptonConfig==0);
      map_whichCR_to_isCR["WZ_2mu1el"] = map_whichCR_to_isCR["WZ"] && (ThreeLeptonConfig==1);
      map_whichCR_to_isCR["WZ_1mu2el"] = map_whichCR_to_isCR["WZ"] && (ThreeLeptonConfig==2);
      map_whichCR_to_isCR["WZ_3mu3el"] = map_whichCR_to_isCR["WZ"] && (ThreeLeptonConfig==3);
      map_whichCR_to_isCR["ZJets_3mu0el"] = map_whichCR_to_isCR["ZJets"] && (ThreeLeptonConfig==0);
      map_whichCR_to_isCR["ZJets_2mu1el"] = map_whichCR_to_isCR["ZJets"] && (ThreeLeptonConfig==1);
      map_whichCR_to_isCR["ZJets_1mu2el"] = map_whichCR_to_isCR["ZJets"] && (ThreeLeptonConfig==2);
      map_whichCR_to_isCR["ZJets_3mu3el"] = map_whichCR_to_isCR["ZJets"] && (ThreeLeptonConfig==3);
      map_whichCR_to_isCR["ZLep_3mu0el"] = map_whichCR_to_isCR["ZLep"] && (ThreeLeptonConfig==0);
      map_whichCR_to_isCR["ZLep_2mu1el"] = map_whichCR_to_isCR["ZLep"] && (ThreeLeptonConfig==1);
      map_whichCR_to_isCR["ZLep_1mu2el"] = map_whichCR_to_isCR["ZLep"] && (ThreeLeptonConfig==2);
      map_whichCR_to_isCR["ZLep_3mu3el"] = map_whichCR_to_isCR["ZLep"] && (ThreeLeptonConfig==3);
      map_whichCR_to_isCR["ZGamma_3mu0el"] = map_whichCR_to_isCR["ZGamma"] && (ThreeLeptonConfig==0);
      map_whichCR_to_isCR["ZGamma_2mu1el"] = map_whichCR_to_isCR["ZGamma"] && (ThreeLeptonConfig==1);
      map_whichCR_to_isCR["ZGamma_1mu2el"] = map_whichCR_to_isCR["ZGamma"] && (ThreeLeptonConfig==2);
      map_whichCR_to_isCR["ZGamma_3mu3el"] = map_whichCR_to_isCR["ZGamma"] && (ThreeLeptonConfig==3);
      map_whichCR_to_isCR["LowMassResonance_3mu0el"] = map_whichCR_to_isCR["LowMassResonance"] && (ThreeLeptonConfig==0);
      map_whichCR_to_isCR["LowMassResonance_2mu1el"] = map_whichCR_to_isCR["LowMassResonance"] && (ThreeLeptonConfig==1);
      map_whichCR_to_isCR["LowMassResonance_1mu2el"] = map_whichCR_to_isCR["LowMassResonance"] && (ThreeLeptonConfig==2);
      map_whichCR_to_isCR["LowMassResonance_3mu3el"] = map_whichCR_to_isCR["LowMassResonance"] && (ThreeLeptonConfig==3);

      if(map_whichCR_to_isCR["ZGamma_2mu1el"]){
        FillHist("TEST_ZGamma_2mu1el_Electron_GeyType", electrontriLooseColl.at(0).GetType(), 1., 0., 100., 100);
      }
      if(map_whichCR_to_isCR["ZGamma_3mu0el"]){
        FillHist("TEST_ZGamma_3mu0el_Muon_GeyType", muontriLooseColl.at(2).GetType(), 1., 0., 100., 100);
      }

      for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
        TString this_suffix = it->first;
        if(it->second){

          FillHist("n_events_"+this_suffix, 0, weight, 0., 1., 1);
          FillHist("n_vertices_"+this_suffix, eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
          FillHist("n_jets_"+this_suffix, n_jets, weight, 0., 10., 10);
          FillHist("n_bjets_"+this_suffix, n_bjets, weight, 0., 10., 10);
          FillHist("n_jets_nearby_"+this_suffix, n_jets_nearby, weight, 0., 10., 10);
          FillHist("n_bjets_nearby_"+this_suffix, n_bjets_nearby, weight, 0., 10., 10);
          FillHist("PFMET_"+this_suffix, MET, weight, 0., 500., 500);
          FillHist("PFMET_phi_"+this_suffix, METphi, weight, -3.2, 3.2, 64);
          FillHist("osllmass_"+this_suffix, m_OSSF[0], weight, 0., 500., 500);
          FillHist("osllmass_"+this_suffix, m_OSSF[1], weight, 0., 500., 500);
          FillHist("m_Z_candidate_"+this_suffix, Z_candidate.M(), weight, 0., 150., 150);
          FillHist("mt_W_candidate_"+this_suffix, MT(nu, WLepton), weight, 0., 300., 300);
          FillHist("m_lll_"+this_suffix, mlll, weight, 0., 500., 500);
          FillHist("Z_candidate_Pt_"+this_suffix, Z_candidate.Pt(), weight, 0., 400., 400);
          FillHist("W_candidate_Pt_"+this_suffix, W_candidate.Pt(), weight, 0., 400., 400);
          FillHist("dRZLeptonWLepton_"+this_suffix, ZLepton_leading.DeltaR(WLepton), weight, 0., 6., 60);
          FillHist("dRZLeptonWLepton_"+this_suffix, ZLepton_subleading.DeltaR(WLepton), weight, 0., 6., 60);
          FillHist("dRMETWLepton_"+this_suffix, nu.DeltaR(WLepton), weight, 0., 6., 60);
          for(unsigned int j=0; j<jetColl_hn_nearby.size(); j++){
            FillHist("dRNearByJetWLepton_"+this_suffix, jetColl_hn_nearby.at(j).DeltaR(WLepton), weight, 0., 6., 60);
            if(jetColl_hn_nearby.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)){
              FillHist("dRNearByBJetWLepton_"+this_suffix, jetColl_hn_nearby.at(j).DeltaR(WLepton), weight, 0., 6., 60);
            }
          }

          FillLeptonKinematicPlot(lep, this_suffix, weight);

          FillHist("ZLepton_leading_Pt_"+this_suffix, ZLepton_leading.Pt(), weight, 0., 200., 200);
          FillHist("ZLepton_leading_Eta_"+this_suffix, ZLepton_leading.Eta(), weight, -3., 3., 60);
          FillHist("ZLepton_leading_RelIso_"+this_suffix, ZLepton_leading.RelIso(), weight, 0., 1.0, 100);
          FillHist("ZLepton_leading_dXY_"+this_suffix, fabs(ZLepton_leading.dXY()), weight, 0., 0.01, 100);
          FillHist("ZLepton_leading_dXYSig_"+this_suffix, fabs(ZLepton_leading.dXYSig()), weight, 0., 8., 80);
          FillHist("ZLepton_leading_dZ_"+this_suffix, fabs(ZLepton_leading.dZ()), weight, 0., 0.5, 500);

          FillHist("ZLepton_subleading_Pt_"+this_suffix, ZLepton_subleading.Pt(), weight, 0., 200., 200);
          FillHist("ZLepton_subleading_Eta_"+this_suffix, ZLepton_subleading.Eta(), weight, -3., 3., 60);
          FillHist("ZLepton_subleading_RelIso_"+this_suffix, ZLepton_subleading.RelIso(), weight, 0., 1.0, 100);
          FillHist("ZLepton_subleading_dXY_"+this_suffix, fabs(ZLepton_subleading.dXY()), weight, 0., 0.01, 100);
          FillHist("ZLepton_subleading_dXYSig_"+this_suffix, fabs(ZLepton_subleading.dXYSig()), weight, 0., 8., 80);
          FillHist("ZLepton_subleading_dZ_"+this_suffix, fabs(ZLepton_subleading.dZ()), weight, 0., 0.5, 500);

          FillHist("WLepton_Pt_"+this_suffix, WLepton.Pt(), weight, 0., 200., 200);
          FillHist("WLepton_Eta_"+this_suffix, WLepton.Eta(), weight, -3., 3., 60);
          FillHist("WLepton_RelIso_"+this_suffix, WLepton.RelIso(), weight, 0., 1.0, 100);
          FillHist("WLepton_dXY_"+this_suffix, fabs(WLepton.dXY()), weight, 0., 0.01, 100);
          FillHist("WLepton_dXYSig_"+this_suffix, fabs(WLepton.dXYSig()), weight, 0., 8., 80);
          FillHist("WLepton_dZ_"+this_suffix, fabs(WLepton.dZ()), weight, 0., 0.5, 500);

        }
      }

    } // Not All Same Charge

  } // isThreeLepton

  if(isFourLepton){

   std::vector<KLepton> lep;
    TString lepOrder[4] = {"leading", "second", "third", "fourth"};
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      KLepton this_lep( muontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      KLepton this_lep( electrontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }

    //==== 4 Muon : FourLeptonConfig = 0;
    //==== 2 Muon + 2 Electron : FourLeptonConfig = 1;
    //==== 4 Electron : FourLeptonConfig = 2;
    int FourLeptonConfig = 0;

    if(muontriLooseColl.size()==4 || electrontriLooseColl.size()==4){

      if(muontriLooseColl.size()==4)     FourLeptonConfig = 0;
      if(electrontriLooseColl.size()==4) FourLeptonConfig = 2;

    }
    else if(muontriLooseColl.size()==2){

      FourLeptonConfig = 1;

    }
    else{
      Message("?", INFO);
    }

    std::vector<KLepton> LepPlus, LepMinus;
    for(unsigned int i=0; i<lep.size(); i++){
      if(lep.at(i).Charge() > 0) LepPlus.push_back(lep.at(i));
      else                       LepMinus.push_back(lep.at(i));
    }

    if( (LepPlus.size() == 2) && (LepMinus.size() == 2) ){

      snu::KParticle ll_case1_1 = LepPlus.at(0)+LepMinus.at(0);
      snu::KParticle ll_case1_2 = LepPlus.at(1)+LepMinus.at(1);
      bool TwoOnZ_case1 = ( fabs( ll_case1_1.M() - m_Z ) < 10. ) && (LepPlus.at(0).LeptonFlavour()==LepMinus.at(0).LeptonFlavour()) &&
                          ( fabs( ll_case1_2.M() - m_Z ) < 10. ) && (LepPlus.at(1).LeptonFlavour()==LepMinus.at(1).LeptonFlavour());

      snu::KParticle ll_case2_1 = LepPlus.at(0)+LepMinus.at(1);
      snu::KParticle ll_case2_2 = LepPlus.at(1)+LepMinus.at(0);
      bool TwoOnZ_case2 = ( fabs( ll_case2_1.M() - m_Z ) < 10. ) && (LepPlus.at(0).LeptonFlavour()==LepMinus.at(1).LeptonFlavour()) &&
                          ( fabs( ll_case2_2.M() - m_Z ) < 10. ) && (LepPlus.at(1).LeptonFlavour()==LepMinus.at(0).LeptonFlavour());

      if( (TwoOnZ_case1 && max(LepPlus.at(0).Pt(),LepMinus.at(0).Pt())>20. && max(LepPlus.at(1).Pt(),LepMinus.at(1).Pt())>20. ) ||
          (TwoOnZ_case2 && max(LepPlus.at(0).Pt(),LepMinus.at(1).Pt())>20. && max(LepPlus.at(1).Pt(),LepMinus.at(0).Pt())>20. )     ){

        std::map< TString, bool > map_whichCR_to_isCR;
        map_whichCR_to_isCR.clear();
        map_whichCR_to_isCR["ZZ"] = true;
        map_whichCR_to_isCR["ZZ_4mu0el"] = map_whichCR_to_isCR["ZZ"] && (FourLeptonConfig==0);
        map_whichCR_to_isCR["ZZ_2mu2el"] = map_whichCR_to_isCR["ZZ"] && (FourLeptonConfig==1);
        map_whichCR_to_isCR["ZZ_0mu4el"] = map_whichCR_to_isCR["ZZ"] && (FourLeptonConfig==2);

        for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
          TString this_suffix = it->first;
          if(it->second){

            FillHist("n_events_"+this_suffix, 0, weight, 0., 1., 1);
            FillHist("n_vertices_"+this_suffix, eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
            FillHist("n_jets_"+this_suffix, n_jets, weight, 0., 10., 10);
            FillHist("n_bjets_"+this_suffix, n_bjets, weight, 0., 10., 10);
            FillHist("n_jets_nearby_"+this_suffix, n_jets_nearby, weight, 0., 10., 10);
            FillHist("n_bjets_nearby_"+this_suffix, n_bjets_nearby, weight, 0., 10., 10);
            FillHist("PFMET_"+this_suffix, MET, weight, 0., 500., 500);
            FillHist("PFMET_phi_"+this_suffix, METphi, weight, -3.2, 3.2, 64);
            if(TwoOnZ_case1){
              FillHist("osllmass_"+this_suffix, ll_case1_1.M(), weight, 0., 500., 500);
              FillHist("osllmass_"+this_suffix, ll_case1_2.M(), weight, 0., 500., 500);
            }
            if(TwoOnZ_case2){
              FillHist("osllmass_"+this_suffix, ll_case2_1.M(), weight, 0., 500., 500);
              FillHist("osllmass_"+this_suffix, ll_case2_2.M(), weight, 0., 500., 500);
            }
            FillHist("m_llll_"+this_suffix, (ll_case1_1+ll_case1_2).M(), weight, 0., 1000., 1000);

            for(unsigned int j=0; j<4; j++){
              TString this_order = lepOrder[j];
              FillHist(this_order+"Lepton_Pt_"+this_suffix, lep[j].Pt(), weight, 0., 200., 200);
              FillHist(this_order+"Lepton_Eta_"+this_suffix, lep[j].Eta(), weight, -3., 3., 60);
              FillHist(this_order+"Lepton_RelIso_"+this_suffix, lep[j].RelIso(), weight, 0., 1.0, 100);
              FillHist(this_order+"Lepton_dXY_"+this_suffix, fabs(lep[j].dXY()), weight, 0., 0.01, 100);
              FillHist(this_order+"Lepton_dXYSig_"+this_suffix, fabs(lep[j].dXYSig()), weight, 0., 8., 80);
              FillHist(this_order+"Lepton_dZ_"+this_suffix, fabs(lep[j].dZ()), weight, 0., 0.5, 500);
            }

          }

        } // Bool Loop

      }

    } // 2OS

  } // isFourLepton
 


  return;

}// End of execute event loop
  


void trilepton_mumumu_CR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_CR::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_CR::~trilepton_mumumu_CR() {
  
  Message("In trilepton_mumumu_CR Destructor" , INFO);
  
}


void trilepton_mumumu_CR::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_CR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_CR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_CRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_mumumu_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}




