// $Id: DiLeptonAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDiLeptonAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "DiLeptonAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"
#include "TSystem.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DiLeptonAnalyzer);

DiLeptonAnalyzer::DiLeptonAnalyzer() :
AnalyzerCore(),
weight_cf(-999), weight_err_cf(-999),
weight_fr(-999), weight_err_fr(-999),
MET(-999), METphi(-999),
ST(-999), HT(-999), LT(-999), contramass(-999),
nbjets(-999), nbjets_fwd(-999), nbjets_nolepveto(-999), n_vtx(-999),
index_jjW_j1(-999), index_jjW_j2(-999),
index_lljjW_j1(-999), index_lljjW_j2(-999),
index_fjW(-999),
RunNtp(false)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("DiLeptonAnalyzer");
  
  Message("In DiLeptonAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void DiLeptonAnalyzer::InitialiseAnalysis() throw( LQError ) {
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

  TDirectory* origDir = gDirectory;

  TString lqdir = getenv("LQANALYZER_DIR");

  //==== Get Fake

  TFile *file_Muon_FR = new TFile( lqdir+"/data/Fake/80X/FakeRate13TeV_muon_2016_opt_all.root");
  TFile *file_Electron_FR = new TFile( lqdir+"/data/Fake/80X/FakeRate13TeV_2016_hnid.root");

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }
  tempDir->cd();

  hist_Electron_FR = (TH2D*)file_Electron_FR->Get("FakeRate_ptcorr_eta40")->Clone();
  hist_Muon_FR = (TH2D*)file_Muon_FR->Get("FakeRate_40_ptcorr_etaSNUTightdijet_0.07_0.005_3_0.04")->Clone();

  file_Muon_FR->Close();
  file_Electron_FR->Close();

  delete file_Muon_FR;
  delete file_Electron_FR;

  origDir->cd();

  return;
}


void DiLeptonAnalyzer::ExecuteEvents()throw( LQError ){

/*
  std::vector< snu::KElectron > testelectrons = GetElectrons(false, true, "ELECTRON_HN_FAKELOOSE");
  //==== FRTEST
  for(unsigned int i=0; i<testelectrons.size(); i++){
    bool fromtau = testelectrons.at(i).MCFromTau();
    FillHist("TEST_ELECTRON_FR_TYPE_F0", testelectrons.at(i).GetType(), 1., 0., 40., 40);
    if(fromtau){
      FillHist("TEST_ELECTRON_fromtau_FR_TYPE_F0", testelectrons.at(i).GetType(), 1., 0., 40., 40);
    }
    if(PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4")){
      FillHist("TEST_ELECTRON_FR_TYPE_F", testelectrons.at(i).GetType(), 1., 0., 40., 40);
      if(fromtau){
        FillHist("TEST_ELECTRON_fromtau_FR_TYPE_F", testelectrons.at(i).GetType(), 1., 0., 40., 40);
      }
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
  //==== CutFlow setup..
  //======================

  std::map< TString, double > w_cutflow;
  std::vector<TString> triggerlist_DiMuon, triggerlist_DiElectron, triggerlist_EMu;

  w_cutflow.clear();
  triggerlist_DiMuon.clear();
  triggerlist_DiElectron.clear();
  triggerlist_EMu.clear();

  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  //triggerlist_DiElectron.push_back("HLT_Ele27_WPTight_Gsf_v");

  triggerlist_EMu.push_back("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v");

  float pileup_reweight=(1.0);
  if(!isData){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  w_cutflow["DiMuon"] = weight*WeightByTrigger(triggerlist_DiMuon, TargetLumi)*pileup_reweight;
  w_cutflow["DiElectron"] = weight*WeightByTrigger(triggerlist_DiElectron, TargetLumi)*pileup_reweight;
  w_cutflow["EMu"] = weight*pileup_reweight; //FIXME

  FillCutFlowByName("DiMuon", "NoCut", w_cutflow["DiMuon"], isData);
  FillCutFlowByName("DiElectron", "NoCut", w_cutflow["DiElectron"], isData);
  FillCutFlowByName("EMu", "NoCut", w_cutflow["EMu"], isData);

  //======================
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);

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

  //==== Muons
  std::vector< snu::KMuon > muons = GetMuons("MUON_HN_LOOSE", false);
  if(muons.size()>=2){
    if( muons.at(0).Pt() > muons.at(1).Pt() ) FillHist("Test_Muon01_CorrectPtOrder", 1., 1., 0., 2., 2);
    else                                      FillHist("Test_Muon01_CorrectPtOrder", 0., 1., 0., 2., 2);
  }
  std::sort(muons.begin(), muons.end(), MuonPtComparing);
  if(muons.size()>=2){
    if( muons.at(0).Pt() > muons.at(1).Pt() ) FillHist("TestAfter_Muon01_CorrectPtOrder", 1., 1., 0., 2., 2);
    else                                      FillHist("TestAfter_Muon01_CorrectPtOrder", 0., 1., 0., 2., 2);
  }

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", muons, 0);
  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muons);

  //==== Electrons
  std::vector< snu::KElectron > electrons = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  double electron_sf = 1.; //FIXME

  bool DoPOGElSF = std::find(k_flags.begin(), k_flags.end(), "DoPOGElSF") != k_flags.end();
  if(DoPOGElSF) electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrons, 0); //FIXME
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);

  //==== If Chargeflip, then shift down electron, and replace electrons with it
  //==== To get correct CF rate, we save original (OS event) electron to electrons_before_shift
  bool RunningChargeFlipData = k_running_chargeflip && isData;
  std::vector< snu::KElectron > electrons_before_shift;
  electrons_before_shift.clear();
  if(RunningChargeFlipData){
    for(unsigned int i=0; i<electrons.size(); i++){
      electrons_before_shift.push_back( electrons.at(i) );
      snu::KElectron tmp_el = electrons.at(i);
      double shift_ = 1.-0.01;
      tmp_el.SetPtEtaPhiM(shift_*tmp_el.Pt(), tmp_el.Eta(), tmp_el.Phi(), 0.511e-3);
      //tmp_el.SetPxPyPzE(shift_*tmp_el.Px(), shift_*tmp_el.Py(), shift_*tmp_el.Pz(), shift_*tmp_el.E());
      electrons.at(i) = tmp_el;
    }
  }

/*
  //==== Quick Fake Histogram name check
  if(muons.size()==2){
  cout << "#### Muon Fake ####" << endl;
  m_datadriven_bkg->Get_DataDrivenWeight_MM(
    false, muons, PassID(muons[0],"MUON_HN_TIGHT"),  PassID(muons[1],"MUON_HN_TIGHT"),
    "Tight0.07_0.005_3_0.04",
    "Tight0.07_0.005_3_0.04",
    true, true, "ptcorr_eta", 0.07, 0.07,false,false);
  }
  if(electrons.size()==2){
  cout << "#### Electron Fake ####" << endl;
  m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electrons, "ELECTRON_HN_FAKELOOSE","ELECTRON_HN_TIGHTv4","40",true);
  }
  return;
*/

  //==== Jets
  std::vector<snu::KJet> jets = GetJets("JET_HN_eta5", 20., 2.5);
  std::vector<snu::KJet> jets_nolepveto = GetJets("JET_HN_eta5_nolepveto", 20., 2.5);
  std::vector<snu::KJet> jets_fwd;
  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");

  std::vector<snu::KJet> jets_eta5 = GetJets("JET_HN_eta5", 20., 5.);
  for(unsigned int i=0; i<jets_eta5.size(); i++){
    if( fabs(jets_eta5.at(i).Eta()) > 2.5 ) jets_fwd.push_back( jets_eta5.at(i) );
  }

  nbjets = 0;
  for(int j=0; j<jets.size(); j++){
    if( IsBTagged(jets.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      nbjets++;
    }
  }

  nbjets_nolepveto = 0;
  for(int j=0; j<jets_nolepveto.size(); j++){
    if( IsBTagged(jets_nolepveto.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      nbjets_nolepveto++;
    }
  }

  nbjets_fwd = 0;
  for(int j=0; j<jets_fwd.size(); j++){
    if( IsBTagged(jets_fwd.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      nbjets_fwd++;
    }
  }

  //==== Lepton Numbers
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO", false);
  std::vector<snu::KMuon> muons_tight; muons_tight.clear();
  std::vector<snu::KElectron> electrons_veto = GetElectrons(false, false, "ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons_tight; electrons_tight.clear();
  std::vector<bool> isT;

  int n_veto_muons = muons_veto.size();
  int n_triLoose_muons = muons.size();
  int n_triTight_muons(0);
  for(unsigned int i=0; i<muons.size(); i++){
    if(PassID(muons.at(i), "MUON_HN_TIGHT")){
      isT.push_back(true);
      muons_tight.push_back(muons.at(i));
      n_triTight_muons++;
    }
    else{
      isT.push_back(false);
    }
  }

  int n_veto_electrons = electrons_veto.size();
  int n_triLoose_electrons = electrons.size();
  int n_triTight_electrons(0);
  for(unsigned int i=0; i<electrons.size(); i++){
    if(PassID(electrons.at(i), "ELECTRON_HN_TIGHTv4")){
      isT.push_back(true);
      electrons_tight.push_back(electrons.at(i));
      n_triTight_electrons++;
    }
    else{
      isT.push_back(false);
    }
  }

  int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
  int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

  bool isTwoMuon_TT    = (n_triTight_muons == 2); // veto third later
  bool isTwoMuon_Loose = (n_triLoose_muons == 2 && n_triTight_muons != 2);

  bool isTwoElectron_TT    = (n_triTight_electrons == 2); // veto third later 
  bool isTwoElectron_Loose = (n_triLoose_electrons == 2 && n_triTight_electrons != 2);

  bool isEMu_TT    = (n_triTight_muons == 1 && n_triTight_electrons == 1);
  bool isEMu_Loose = (n_triLoose_muons == 1)     &&
                     (n_triLoose_electrons == 1) &&
                     (n_triTight_leptons != 2);

  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  RunNtp = std::find(k_flags.begin(), k_flags.end(), "RunNtp") != k_flags.end();

  if(!isData && !NonPromptRun){
    weight*=muon_id_iso_sf;
    weight*=MuTrkEffSF;
    weight*=pileup_reweight;
    weight*=GetKFactor();
    weight*=electron_sf;
    weight*=electron_RecoSF;
  }

  snu::KEvent Evt = eventbase->GetEvent();
  MET = Evt.MET();
  METphi = Evt.METPhi();
  CorrectedMETRochester(muons, MET, METphi);
  CorrectedMETElectron(0, electrons, MET, METphi);

  n_vtx = Evt.nVertices();

  //==== Define Analysis Region

  std::vector< TString > Suffixs;
  std::vector< std::vector<TString> > Triggers;
  std::vector< bool > isTTs, isLOOSEs, isNoExtra, isNoExtraOtherFlavour;

  bool RunningNonPromptData = NonPromptRun && isData;

/*
  Suffixs.push_back("DiMuon");
  Triggers.push_back(triggerlist_DiMuon);
  isTTs.push_back( isTwoMuon_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isTwoMuon_Loose && RunningNonPromptData );
  isNoExtra.push_back( n_veto_muons == 2 );
  isNoExtraOtherFlavour.push_back( n_veto_electrons == 0 );
*/

  Suffixs.push_back("DiElectron");
  Triggers.push_back(triggerlist_DiElectron);
  isTTs.push_back( isTwoElectron_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isTwoElectron_Loose && RunningNonPromptData );
  isNoExtra.push_back( n_veto_electrons == 2 );
  isNoExtraOtherFlavour.push_back( n_veto_muons == 0 );

/*
  Suffixs.push_back("EMu");
  Triggers.push_back(triggerlist_EMu);
  isTTs.push_back( isEMu_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isEMu_Loose && RunningNonPromptData );
  isNoExtra.push_back( n_veto_muons == 1 && n_veto_electrons == 1 );
  isNoExtraOtherFlavour.push_back( true  );
*/
  for(unsigned int i=0; i<Suffixs.size(); i++){

    TString Suffix = Suffixs.at(i);

    //==== Trigger pass
    if(!PassTriggerOR( Triggers.at(i) )) continue;
    FillCutFlowByName(Suffix, "MET_PV_Trig", w_cutflow[Suffix], isData);

    //==== Two leptons
    if(!isTTs.at(i) && !isLOOSEs.at(i)) continue;

    //==== That two lepton pass basic pt cuts
    if(Suffix=="DiMuon"){
      //==== to properly veto below tricky event
      //==== ### failing third muon veto ###
      //==== muons.size() = 3
      //==== muons.at(0).Pt() = 51.8417 => isTight = 0
      //==== muons.at(1).Pt() = 27.8285 => isTight = 1
      //==== muons.at(2).Pt() = 8.72782 => isTight = 1
      if(isTTs.at(i)){
        if((muons_tight.at(0).Pt() < 20.) || (muons_tight.at(1).Pt() < 10.)) continue;
      }
      else{
        if((muons.at(0).Pt() < 20.) || (muons.at(1).Pt() < 10.)) continue;
      }
    }
    if(Suffix=="DiElectron"){
      if(isTTs.at(i)){
        if(electrons_tight.at(0).Pt() < 25. || electrons_tight.at(1).Pt() < 15.) continue;
      }
      else{
        if(electrons.at(0).Pt() < 25. || electrons.at(1).Pt() < 15.) continue;
      }
    }
    if(Suffix=="EMu"){
      if(muons.at(0).Pt() < 10. || electrons.at(1).Pt() < 15.) continue;
    }
    FillCutFlowByName(Suffix, "TwoLeptons", w_cutflow[Suffix], isData);

    //==== No Extra lepton
    if(!isNoExtra.at(i)) continue;
    FillCutFlowByName(Suffix, "NoExtraLepton", w_cutflow[Suffix], isData);

    //==== No Extra different flavour lepton
    if(!isNoExtraOtherFlavour.at(i)) continue;
    FillCutFlowByName(Suffix, "NoExtraFlavourLepton", w_cutflow[Suffix], isData);

    //==== If Loose event, we only consider Data
    if(isLOOSEs.at(i) && !isData) continue;

    //==== DiMuon-DoubleMuon PD / ...
    if(isData && k_channel != "DoubleMuon_CF"){
      if(Suffix == "DiMuon"){
        if(k_channel != "DoubleMuon") continue;
      }
      if(Suffix == "DiElectron"){
        if(k_channel != "DoubleEG") continue;
      }
      if(Suffix == "EMu"){
        if(k_channel != "MuonEG") continue;
      }
    }

    double trigger_ps_weight = WeightByTrigger(Triggers.at(i), TargetLumi);
    double this_weight = weight*trigger_ps_weight;

    double trigger_sf = 1.;
    if(!isData && Suffix=="DiMuon"){
      double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 0, 0, 0);
      double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 0, 1, 0);
      trigger_sf = trigger_eff_Data/trigger_eff_MC;
    }

    this_weight *= trigger_sf;

    std::vector<KLepton> lep;
    for(unsigned int j=0; j<muons.size(); j++){
      KLepton this_lep( muons.at(j) );
      lep.push_back( this_lep );
    }
    for(unsigned int j=0; j<electrons.size(); j++){
      KLepton this_lep( electrons.at(j) );
      lep.push_back( this_lep );
    }

    bool isSS = lep.at(0).Charge() == lep.at(1).Charge();
    bool isSSForCF = isSS;
    if(RunningChargeFlipData) isSSForCF = !isSS;

    double m_Z = 91.1876;
    bool isOffZ = fabs( (lep.at(0)+lep.at(1)).M() - m_Z ) > 10.;

    //==== mll Cut Study
    FillHist("CutStudy_m_ll_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
    if(isSS) FillHist("CutStudy_m_ll_SS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
    else FillHist("CutStudy_m_ll_OS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
    bool mll10GeV = ( lep.at(0)+lep.at(1) ).M() < 10.;
    bool mll12GeV = ( lep.at(0)+lep.at(1) ).M() < 12.;
    //if(!isSS && !isOffZ) continue;
    if(mll10GeV) continue;
    FillCutFlowByName(Suffix, "LowDileptonMass", w_cutflow[Suffix], isData);

    double this_weight_err(0.);
    if( isLOOSEs.at(i) ){
      //this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      //this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      //m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");

      get_eventweight(muons, electrons, isT, 0);

      this_weight = weight_fr;
      this_weight_err = weight_err_fr;
      
    }
    if( RunningChargeFlipData && Suffix=="DiElectron" && !isSS ){
      GetCFWeight(electrons_before_shift.at(0), electrons_before_shift.at(1));

      this_weight = weight_cf;
      this_weight_err = weight_err_cf;
    }

    //==== Fill Histogram

    std::map< TString, bool > map_Region_to_Bool;
    map_Region_to_Bool.clear();

    //==== Default Variation
    //==== 1) OnZ/OffZ/All
    //==== 2) OS/SS/All

    //==== Two Lepton
    map_Region_to_Bool[Suffix] = true;

    //==== # of jets
    map_Region_to_Bool[Suffix+"_0jets"] = (jets.size()==0);
    map_Region_to_Bool[Suffix+"_1jets"] = (jets.size()==1);
    map_Region_to_Bool[Suffix+"_2jets"] = (jets.size()==2); // SR
    //map_Region_to_Bool[Suffix+"_3jets"] = (jets.size()==3);
    //map_Region_to_Bool[Suffix+"_4jets"] = (jets.size()==4);
    map_Region_to_Bool[Suffix+"_Inclusive2jets"] = (jets.size()>=2); // SR

    //==== # of (normal) b jet
    //map_Region_to_Bool[Suffix+"_0bjets"] = (nbjets==0);
    //map_Region_to_Bool[Suffix+"_Inclusive1bjets"] = (nbjets>=1);
    //map_Region_to_Bool[Suffix+"_Inclusive1bjets_METge50"] = (nbjets>=1) && (MET>50.); // control region
    //map_Region_to_Bool[Suffix+"_Inclusive2bjets"] = (nbjets>=2);
    //map_Region_to_Bool[Suffix+"_Inclusive3bjets"] = (nbjets>=3);

    //==== # of forward b jet
    map_Region_to_Bool[Suffix+"_0nlbjets"] = (nbjets_nolepveto==0);
    map_Region_to_Bool[Suffix+"_Inclusive1nlbjets"] = (nbjets_nolepveto>=1);
    map_Region_to_Bool[Suffix+"_Inclusive1nlbjets_METge50"] = (nbjets_nolepveto>=1) && (MET>50.); // control region
    map_Region_to_Bool[Suffix+"_Inclusive2nlbjets"] = (nbjets_nolepveto>=2);
    //map_Region_to_Bool[Suffix+"_Inclusive3nlbjets"] = (nbjets_nolepveto>=3);

    //==== fatjet
    map_Region_to_Bool[Suffix+"_0fatjets"] = (fatjets.size()==0);
    map_Region_to_Bool[Suffix+"_1fatjets"] = (fatjets.size()==1); // together with onZ, we can get a W(->jj)Z(->ll) control region
    map_Region_to_Bool[Suffix+"_Inclusive1fatjets"] = (fatjets.size()>=1);
    map_Region_to_Bool[Suffix+"_Inclusive1fatjets_METlt50"] = (fatjets.size()>=1) && (MET<50.);

    ST = lep.at(0).Pt() + lep.at(1).Pt();
    //float looseST = lep.at(0).Pt() + lep.at(1).Pt();
    HT = 0.;
    for(unsigned int ij=0; ij <jets_nolepveto.size(); ij++){
      ST += jets_nolepveto.at(ij).Pt();
    }
    for(unsigned int ij=0; ij <jets.size(); ij++){
      HT += jets.at(ij).Pt();
    }
    LT = lep.at(0).Pt() + lep.at(1).Pt();


    //==== At least two jets
    if( jets.size()>=2 ){
      index_jjW_j1 = 0;
      index_jjW_j2 = 1;
      double mjj = GetDijetMassClosest(jets, 80.4, index_jjW_j1, index_jjW_j2);
      index_lljjW_j1 = 0;
      index_lljjW_j2 = 1;
      double mlljj = GetDileptonDijetMassClosest(lep, jets, 80.4, index_lljjW_j1, index_lljjW_j2);

      //==== More CutFlows
      if( jets.size()>=2 ){
        FillCutFlowByName(Suffix, "InclusiveTwoJets", w_cutflow[Suffix], isData);
        if( MET < 50. ){
          FillCutFlowByName(Suffix, "MET50", w_cutflow[Suffix], isData);
          if( mjj < 200. ){
            FillCutFlowByName(Suffix, "Mjj200", w_cutflow[Suffix], isData);
            if( nbjets == 0 ){
              FillCutFlowByName(Suffix, "NoBJet", w_cutflow[Suffix], isData);
              if( nbjets_nolepveto == 0 ){
                FillCutFlowByName(Suffix, "NoBJet_nolepveto", w_cutflow[Suffix], isData);
              }
            }
          }
        }
      }

      float dPhi = fabs(TVector2::Phi_mpi_pi(jets.at(index_jjW_j1).Phi() - jets.at(index_jjW_j2).Phi()));
      contramass=2*jets.at(index_jjW_j1).Pt()*jets.at(index_jjW_j2).Pt()*(1+cos(dPhi));
      contramass=sqrt(contramass);

/*
      map_Region_to_Bool[Suffix+"_Low"] = map_Region_to_Bool[Suffix+"_Preselection"] &&
                                          ( (lep.at(0)+lep.at(1)).M() < 70. ) &&
                                          ( contramass < 100. ) &&
                                          ( (jets.at(0)+jets.at(1)).M() < 250. ) &&
                                          ( lep.at(0).DeltaR( lep.at(1) ) < 3.5 ) &&
                                          ( MET*MET/ST < 15. && MET < 80. ) &&
                                          ( ST < 450. );

      map_Region_to_Bool[Suffix+"_Medium"] = map_Region_to_Bool[Suffix+"_Preselection"] &&
                                             ( MET*MET/ST < 12.5 && MET < 80. );

      map_Region_to_Bool[Suffix+"_High"] = map_Region_to_Bool[Suffix+"_Preselection"] &&
                                           ( MET*MET/ST < 8. ) &&
                                           ( ST > 400. ) &&
                                           ( (lep.at(0)+lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M() > 200. ) &&
                                           ( (lep.at(0).Pt() > 50.) && (lep.at(1).Pt() > 25.) );
*/
    }
    if(fatjets.size() >= 1){
      GetFatjetMassClosest(fatjets, 80.4, index_fjW);
    }

    //==== If OS, we should veto m(OS) < 12 GeV..
    //if(!isSS && mll12GeV) continue; //FIXME just for now. I will use this

    map_Region_to_Bool[Suffix+"_Preselection"] = isOffZ &&
                                                 (nbjets_nolepveto==0) &&
                                                 ( ( lep.at(0)+lep.at(1) ).M() > 20. ) &&
                                                 (  (jets.size()>=2) || (fatjets.size()>=1) );
    map_Region_to_Bool[Suffix+"_Preselection_Inclusive2jets"] = isOffZ &&
                                                 (nbjets_nolepveto==0) &&
                                                 ( ( lep.at(0)+lep.at(1) ).M() > 20. ) &&
                                                 ( (jets.size()>=2) );

    if( map_Region_to_Bool[Suffix+"_Preselection_Inclusive2jets"] ){

      map_Region_to_Bool[Suffix+"_Preselection_Inclusive2jets_mjj50to110"] = isOffZ &&
                                                   (nbjets_nolepveto==0) &&
                                                   ( ( lep.at(0)+lep.at(1) ).M() > 20. ) &&
                                                   ( (jets.size()>=2) ) &&
                                                   ( (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M() > 50. ) &&
                                                   ( (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M() < 110. );
    }
    map_Region_to_Bool[Suffix+"_Preselection_Inclusive1fatjets"] = isOffZ &&
                                                 (nbjets_nolepveto==0) &&
                                                 ( ( lep.at(0)+lep.at(1) ).M() > 20. ) &&
                                                 ( (fatjets.size()>=1) );


    bool NtupleSkim = map_Region_to_Bool[Suffix+"_Preselection"] && isSSForCF;
    if(RunNtp && NtupleSkim){
      double cutop[100];
      cutop[0] = lep.at(0).Pt();
      cutop[1] = lep.at(1).Pt();
      cutop[2] = lep.at(0).DeltaR( lep.at(1) );
      cutop[3] = (lep.at(0)+lep.at(1)).M();
      if(!RunningChargeFlipData) cutop[4] = isSS ? 1 : 0;
      else                       cutop[4] = isSS ? 0 : 1;
      cutop[5] = isOffZ ? 1 : 0;

      cutop[6] = jets.size();
      cutop[7] = nbjets;
      cutop[8] = jets_nolepveto.size();
      cutop[9] = nbjets_nolepveto;
      cutop[10] = jets_fwd.size();
      cutop[11] = nbjets_fwd;

      //==== Two Jets
      if(jets.size() >= 2){

        //==== pt order
        cutop[12] = jets.at(0).Pt();
        cutop[13] = jets.at(1).Pt();
        cutop[14] = jets.at(0).DeltaR( jets.at(1) );
        cutop[15] = (jets.at(0)+jets.at(1)).M();
        cutop[16] = (lep.at(0)+jets.at(0)+jets.at(1)).M();
        cutop[17] = (lep.at(1)+jets.at(0)+jets.at(1)).M();
        cutop[18] = (lep.at(0)+lep.at(1)+jets.at(0)+jets.at(1)).M();

        //==== m(jj)~W (high mass)
        cutop[19] = jets.at(index_jjW_j1).Pt();
        cutop[20] = jets.at(index_jjW_j2).Pt();
        cutop[21] = (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
        cutop[22] = (lep.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
        cutop[23] = (lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
        cutop[24] = (lep.at(0)+lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
        cutop[25] = jets.at(index_jjW_j1).DeltaR(jets.at(index_jjW_j2));
        cutop[26] = lep.at(0).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) );
        cutop[27] = lep.at(1).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) );
        cutop[28] = lep.at(0).DeltaR( lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) );
        cutop[29] = lep.at(1).DeltaR( lep.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) );

        //==== m(lljj)~W (low mass)
        cutop[30] = jets.at(index_lljjW_j1).Pt();
        cutop[31] = jets.at(index_lljjW_j2).Pt();
        cutop[32] = (jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        cutop[33] = (lep.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        cutop[34] = (lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        cutop[35] = (lep.at(0)+lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        cutop[36] = jets.at(index_lljjW_j1).DeltaR(jets.at(index_lljjW_j2));
        cutop[37] = lep.at(0).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );
        cutop[38] = lep.at(1).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );
        cutop[39] = lep.at(0).DeltaR( lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );
        cutop[40] = lep.at(1).DeltaR( lep.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );


      }
      else{
        for(int i=12;i<=40;i++) cutop[i] = -999.;
      }

      //==== Two Forward Jets
      if(jets_fwd.size() >= 2){
        cutop[41] = jets_fwd.at(0).DeltaR( jets_fwd.at(1) );
      }
      else{
        cutop[41] = -999;
      }

      //==== Transverse Energies
      cutop[42] = MET;
      cutop[43] = ST;
      cutop[44] = HT;
      cutop[45] = LT;

      cutop[46] = this_weight;
      cutop[47] = this_weight_err;
      cutop[48] = fatjets.size();

      if(fatjets.size() >= 1){

        snu::KFatJet fatjet = fatjets.at(0);
        snu::KFatJet fatjet_W = fatjets.at(index_fjW);

        cutop[49] = fatjet.Pt();
        cutop[50] = fatjet.M();
        cutop[51] = fatjet.PrunedMass();
        cutop[52] = fatjet.SoftDropMass();
        cutop[53] = fatjet.Tau2()/fatjet.Tau1();
        cutop[54] = fatjet.Tau3()/fatjet.Tau2();
        cutop[55] = (lep.at(0)+fatjet).M();
        cutop[56] = (lep.at(1)+fatjet).M();
        cutop[57] = (lep.at(0)+lep.at(1)+fatjet).M();
        cutop[58] = lep.at(0).DeltaR( fatjet );
        cutop[59] = lep.at(1).DeltaR( fatjet );
        cutop[60] = lep.at(0).DeltaR( lep.at(1)+fatjet );
        cutop[61] = lep.at(1).DeltaR( lep.at(0)+fatjet );

        cutop[62] = fatjet_W.Pt();
        cutop[63] = fatjet_W.M();
        cutop[64] = fatjet_W.PrunedMass();
        cutop[65] = fatjet_W.SoftDropMass();
        cutop[66] = fatjet_W.Tau2()/fatjet_W.Tau1();
        cutop[67] = fatjet_W.Tau3()/fatjet_W.Tau2();
        cutop[68] = (lep.at(0)+fatjet_W).M();
        cutop[69] = (lep.at(1)+fatjet_W).M();
        cutop[70] = (lep.at(0)+lep.at(1)+fatjet_W).M();
        cutop[71] = lep.at(0).DeltaR( fatjet_W );
        cutop[72] = lep.at(1).DeltaR( fatjet_W );
        cutop[73] = lep.at(0).DeltaR( lep.at(1)+fatjet_W );
        cutop[74] = lep.at(1).DeltaR( lep.at(0)+fatjet_W );

      }
      else{
        for(int i=49;i<=74;i++) cutop[i] = -999.;
      }

      FillNtp("Ntp_"+Suffix+"_Preselection_SS",cutop);

    }
    //if(RunNtp) continue;

    for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
      TString this_suffix = it->first;
      if(it->second){

        //=====================
        //==== NO Charge Flip
        //=====================

        if(!RunningChargeFlipData){

          //=========================
          //==== Filling Histograms
          //=========================

          //==== SS+OS
          FillDiLeptonPlot(this_suffix+"_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
          //==== SS
          if(isSS){
            FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            FillHist("LeptonType_"+this_suffix+"_SS", lep.at(0).GetType(), 1., 0., 50., 50);
            FillHist("LeptonType_"+this_suffix+"_SS", lep.at(1).GetType(), 1., 0., 50., 50);
          }
          //==== OS
          else{
            FillDiLeptonPlot(this_suffix+"_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
          }

          //==== OffZ
          if(isOffZ){

            //==== SS+OS
            FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

            //==== SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }
            //==== OS
            else{
              FillDiLeptonPlot(this_suffix+"_OffZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }

          } // Off-Z
          else{

            //==== SS+OS
            FillDiLeptonPlot(this_suffix+"_OnZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

            //==== SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }
            //==== OS
            else{
              FillDiLeptonPlot(this_suffix+"_OnZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }

          }


        }

        //=================
        //==== ChargeFlip
        //=================

        //==== To NOT DRAWING OS plot from CF (file size..)
        //==== Separate CF from others

        else{
          //==== using OS event, weight CF and estimate SS
          if(Suffix=="DiElectron" && !isSS){

            //=========================
            //==== Filling Histograms
            //=========================

            FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            if(isOffZ){
              FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }
            else{
              FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }
            
          } // END fill chargeflip only for DiElectron OS
        }


      } // END passing this region
    } // END Search Region loop


  } // END Suffix (Channel; DiMuon, DiElectron, EMu) loop



  return;

} // End of execute event loop
  


void DiLeptonAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void DiLeptonAnalyzer::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  return;
  
}

DiLeptonAnalyzer::~DiLeptonAnalyzer() {
  
  Message("In DiLeptonAnalyzer Destructor" , INFO);
  
}


void DiLeptonAnalyzer::FillCutFlow(TString cut, float w){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,w);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12,0.,12.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    
  }
}

void DiLeptonAnalyzer::FillCutFlowByName(TString histname, TString cut, float w, bool IsDATA){

  TString this_histname = "Cutflow_"+histname;

  FillHist(this_histname+"_"+cut, 0., w, 0., 1., 1);

}


void DiLeptonAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void DiLeptonAnalyzer::MakeHistograms(){
  //// Additional plots to make

  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);

  MakeNtp("Ntp_DiElectron_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:Nfatjets:FatJet_Pt:m_fatjet:FatJet_PrunedMass:FatJet_SoftDropMass:FatJet_Tau21:FatJet_Tau32:m_Leadlfatjet:m_SubLeadlfatjet:m_llfatjet:DeltaRLeadl_fatjet:DeltaRSubLeadl_fatjet:DeltaRLeadl_SubLeadlfatjet:DeltaRSubLeadl_Leadlfatjet:FatJet_Wclosest_Pt:m_fatjet_Wclosest:FatJet_Wclosest_PrunedMass:FatJet_Wclosest_SoftDropMass:FatJet_Wclosest_Tau21:FatJet_Wclosest_Tau32:m_Leadlfatjet_Wclosest:m_SubLeadlfatjet_Wclosest:m_llfatjet_Wclosest:DeltaRLeadl_fatjet_Wclosest:DeltaRSubLeadl_fatjet_Wclosest:DeltaRLeadl_SubLeadlfatjet_Wclosest:DeltaRSubLeadl_Leadlfatjet_Wclosest");



  /**
   *  Remove//Overide this DiLeptonAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/

}


void DiLeptonAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //


}

void DiLeptonAnalyzer::FillDiLeptonPlot(
  TString histsuffix,
  std::vector< KLepton> leptons,
  std::vector< snu::KJet > jets,
  std::vector< snu::KJet > jets_fwd,
  std::vector< snu::KJet > jets_nolepveto,
  std::vector< snu::KFatJet > fatjets,
  double thisweight,
  double thieweighterr
  ){

  TString leporder[4] = {"leading", "second", "third", "fourth"};

  FillHist("Nevents_"+histsuffix, 0., thisweight, 0., 1., 1);
  FillHist("PFMET_"+histsuffix, MET, thisweight, 0., 1000., 1000);
  FillHist("PFMET_phi_"+histsuffix, METphi, thisweight, -3.2, 3.2, 100);
  FillHist("Njets_"+histsuffix, jets.size(), thisweight, 0., 10., 10);
  FillHist("Njets_nolepveto_"+histsuffix, jets_nolepveto.size(), thisweight, 0., 10., 10);
  FillHist("Nfwdjets_"+histsuffix, jets_fwd.size(), thisweight, 0., 10., 10);
  FillHist("Nbjets_"+histsuffix, nbjets, thisweight, 0., 10., 10);
  FillHist("Nbjets_nolepveto_"+histsuffix, nbjets_nolepveto, thisweight, 0., 10., 10);
  FillHist("Nbfwdjets_"+histsuffix, nbjets_fwd, thisweight, 0., 10., 10);
  FillHist("Nfatjets_"+histsuffix, fatjets.size(), thisweight, 0., 10., 10);
  FillHist("Nvtx_"+histsuffix, n_vtx, thisweight, 0., 100., 100);

  FillHist("HT_"+histsuffix, HT, thisweight, 0., 2000., 2000);
  FillHist("ST_"+histsuffix, ST, thisweight, 0., 2000., 2000);
  FillHist("LT_"+histsuffix, LT, thisweight, 0., 2000., 2000);
  FillHist("MCT_"+histsuffix, contramass, thisweight, 0., 2000., 2000);
  FillHist("MET2overST_"+histsuffix, MET*MET/ST, thisweight, 0., 2000., 2000);
  FillHist("DeltaRl1l2_"+histsuffix, leptons.at(0).DeltaR(leptons.at(1)), thisweight, 0., 10., 100);

  FillHist("m_ll_"+histsuffix, (leptons.at(0)+leptons.at(1)).M(), thisweight, 0., 2000., 2000);

  for(int i=0; i<leptons.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"Lepton_Pt_"+histsuffix,  leptons.at(i).Pt(), thisweight, 0., 2000., 2000);
    FillHist(leporder[i]+"Lepton_Eta_"+histsuffix, leptons.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(int i=0; i<jets.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"Jet_Pt_"+histsuffix,  jets.at(i).Pt(), thisweight, 0., 2000., 2000);
    FillHist(leporder[i]+"Jet_Eta_"+histsuffix, jets.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(int i=0; i<jets_fwd.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"ForwardJet_Pt_"+histsuffix,  jets_fwd.at(i).Pt(), thisweight, 0., 2000., 2000);
    FillHist(leporder[i]+"ForwardJet_Eta_"+histsuffix, jets_fwd.at(i).Eta(), thisweight, -5., 5., 100);
  }
  for(int i=0; i<jets_nolepveto.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"NoLepVetoJet_Pt_"+histsuffix, jets_nolepveto.at(i).Pt(), thisweight, 0., 2000., 2000);
    FillHist(leporder[i]+"NoLepVetoJet_Eta_"+histsuffix, jets_nolepveto.at(i).Eta(), thisweight, -3., 3., 60);
  }

  if(jets.size() >= 2){
    //==== m(jj) closeset to m(W) : high mass scenario
    FillHist("m_jj_jjWclosest_"+histsuffix, (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_Leadljj_jjWclosest_"+histsuffix, (leptons.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadljj_jjWclosest_"+histsuffix, (leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_lljj_jjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("DeltaRjjWclosest_"+histsuffix, jets.at(index_jjW_j1).DeltaR(jets.at(index_jjW_j2)), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_jjWclosest_"+histsuffix, leptons.at(0).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_jjWclosest_"+histsuffix, leptons.at(1).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_SubLeadljjWclosest_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_LeadljjWclosest_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);

    //==== m(lljj) cloeset to m(W) : low mass scenario
    FillHist("m_jj_lljjWclosest_"+histsuffix, (jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_Leadljj_lljjWclosest_"+histsuffix, (leptons.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadljj_lljjWclosest_"+histsuffix, (leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_lljj_lljjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    FillHist("DeltaRlljjWclosest_"+histsuffix, jets.at(index_lljjW_j1).DeltaR(jets.at(index_lljjW_j2)), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_lljjWclosest_"+histsuffix, leptons.at(0).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_lljjWclosest_"+histsuffix, leptons.at(1).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_SubLeadllljjWclosest_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_LeadllljjWclosest_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);

    FillHist("DeltaRjjptorder_"+histsuffix, jets.at(0).DeltaR( jets.at(1) ), thisweight, 0., 10., 100);
    FillHist("m_jjptorder_"+histsuffix, (jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_Leadljjptorder_"+histsuffix, (leptons.at(0)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadljjptorder_"+histsuffix, (leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_lljjptorder_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
  }

  if(fatjets.size() >= 1){
    snu::KFatJet fatjet = fatjets.at(0);
    snu::KFatJet fatjet_W = fatjets.at(index_fjW);

    //==== Leading FatJet
    FillHist("m_fatjet_"+histsuffix, fatjet.M(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Pt_"+histsuffix, fatjet.Pt(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Eta_"+histsuffix, fatjet.Eta(), thisweight, -3., 3., 60);
    FillHist("FatJet_PrunedMass_"+histsuffix, fatjet.PrunedMass(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_SoftDropMass_"+histsuffix, fatjet.SoftDropMass(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Tau21_"+histsuffix, fatjet.Tau2()/fatjet.Tau1(), thisweight, 0., 2., 200);
    FillHist("FatJet_Tau32_"+histsuffix, fatjet.Tau3()/fatjet.Tau2(), thisweight, 0., 2., 200);
    FillHist("m_Leadlfatjet_"+histsuffix, (leptons.at(0)+fatjet).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadlfatjet_"+histsuffix, (leptons.at(1)+fatjet).M(), thisweight, 0., 2000., 2000);
    FillHist("m_llfatjet_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjet).M(), thisweight, 0., 2000., 2000);
    FillHist("DeltaRLeadl_fatjet_"+histsuffix, leptons.at(0).DeltaR( fatjet ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_fatjet_"+histsuffix, leptons.at(1).DeltaR( fatjet ), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_SubLeadlfatjet_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+fatjet ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_Leadlfatjet_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+fatjet ), thisweight, 0., 10., 100);

    //==== W closest FatJet
    FillHist("m_fatjet_Wclosest_"+histsuffix, fatjet_W.M(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Wclosest_Pt_"+histsuffix, fatjet_W.Eta(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Wclosest_Eta_"+histsuffix, fatjet_W.Eta(), thisweight, -3., 3., 60);
    FillHist("FatJet_Wclosest_PrunedMass_"+histsuffix, fatjet_W.PrunedMass(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Wclosest_SoftDropMass_"+histsuffix, fatjet_W.SoftDropMass(), thisweight, 0., 2000., 2000);
    FillHist("FatJet_Wclosest_Tau21_"+histsuffix, fatjet_W.Tau2()/fatjet_W.Tau1(), thisweight, 0., 2., 200);
    FillHist("FatJet_Wclosest_Tau32_"+histsuffix, fatjet_W.Tau3()/fatjet_W.Tau2(), thisweight, 0., 2., 200);
    FillHist("m_Leadlfatjet_Wclosest_"+histsuffix, (leptons.at(0)+fatjet_W).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadlfatjet_Wclosest_"+histsuffix, (leptons.at(1)+fatjet_W).M(), thisweight, 0., 2000., 2000);
    FillHist("m_llfatjet_Wclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjet_W).M(), thisweight, 0., 2000., 2000);
    FillHist("DeltaRLeadl_fatjet_Wclosest_"+histsuffix, leptons.at(0).DeltaR( fatjet_W ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_fatjet_Wclosest_"+histsuffix, leptons.at(1).DeltaR( fatjet_W ), thisweight, 0., 10., 100);
    FillHist("DeltaRLeadl_SubLeadlfatjet_Wclosest_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+fatjet_W ), thisweight, 0., 10., 100);
    FillHist("DeltaRSubLeadl_Leadlfatjet_Wclosest_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+fatjet_W ), thisweight, 0., 10., 100);


  }

  if(thieweighterr!=0.){
    FillDiLeptonPlot(histsuffix+"_up",   leptons, jets, jets_fwd, jets_nolepveto, fatjets, thisweight + thieweighterr, 0.);
    FillDiLeptonPlot(histsuffix+"_down", leptons, jets, jets_fwd, jets_nolepveto, fatjets, thisweight - thieweighterr, 0.);
  }

}




void DiLeptonAnalyzer::GetCFWeight(KLepton lep1, KLepton lep2){

  if(lep1.Charge()==lep2.Charge()) return;

  //==== Okay, now lep1 and lep2 are OS

  double cf1 = GetCF(lep1, false);
  double cf2 = GetCF(lep2, false);
  double cf1_err = GetCF(lep1, true);
  double cf2_err = GetCF(lep2, true);

  weight_cf = cf1/(1.-cf1) + cf2/(1.-cf2);
  weight_err_cf = sqrt( cf1_err/( (1.-cf1)*(1.-cf1) ) + cf2_err/( (1.-cf2)*(1.-cf2) ) );

}

double DiLeptonAnalyzer::GetCF(KLepton lep, bool geterr){

  double el_eta = fabs(lep.GetElectronPtr()->SCEta());
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./lep.Pt();
  double a = 999., b= 999.;
  if(el_eta < 0.9){
    if(invPt< 0.023){
      a=(-0.00138635);
      b=(4.35054e-05);
    }
    else{
      a=(0.00114356);
      b=(-1.55941e-05);
    }
  }
  else if(el_eta < 1.4442){
    if(invPt < 0.016){
      a=(-0.0369937);
      b=(0.000797434);
    }
    else if(invPt < 0.024){
      a=(-0.0159017);
      b=(0.00046038);
    }
    else{
      a=(-0.00214657);
      b=(0.000147245);
    }
  }
  else{
    if(invPt< 0.012){
      a=(-0.4293);
      b=(0.00641511);
    }
    else if(invPt< 0.020){
      a=(-0.104796);
      b=(0.00256146);
    }
    else{
      a=(-0.0161499);
      b=(0.00076872);
    }
  }

  double sf(1.);
  if(el_eta < 1.4442) sf = 0.75362822;
  else sf = 0.821682654;

  double sys=0.;
  double rate = (a)*invPt + (b);
  if(rate < 0) rate = 0.;

  rate *= sf;

  if(!geterr) return rate;
  else return 0.;

}

double DiLeptonAnalyzer::GetDijetMassClosest(std::vector<snu::KJet> js, double mass, int& m, int& n){

  if(js.size()<2) return -999;

  double m_close(9999999999999999.);
  for(unsigned int i=0; i<js.size()-1; i++){
    for(unsigned int j=i+1; j<js.size(); j++){

      double m_tmp = (js.at(i)+js.at(j)).M();
      if( fabs(m_tmp-mass) < fabs(m_close-mass) ){
        m = i;
        n = j;
        m_close = m_tmp;
      }

    }
  }

  return m_close;

}

double DiLeptonAnalyzer::GetDileptonDijetMassClosest(std::vector<KLepton> leps, std::vector<snu::KJet> js, double mass, int& m, int& n){

  if(leps.size() != 2) return -999;
  if(js.size()<2) return -999;

  double m_close(9999999999999999.);
  for(unsigned int i=0; i<js.size()-1; i++){
    for(unsigned int j=i+1; j<js.size(); j++){

      double m_tmp = (leps.at(0)+leps.at(1)+js.at(i)+js.at(j)).M();
      if( fabs(m_tmp-mass) < fabs(m_close-mass) ){
        m = i;
        n = j;
        m_close = m_tmp;
      }

    }
  }

  return m_close;

}

double DiLeptonAnalyzer::GetFatjetMassClosest(std::vector<snu::KFatJet> fjs, double mass, int& m){

  if(fjs.size()==0) return -999;

  double m_close(9999999999999999.);
  for(unsigned int i=0; i<fjs.size(); i++){

    double m_tmp = fjs.at(i).M();
    if( fabs(m_tmp-mass) < fabs(m_close-mass) ){
      m = i;
      m_close = m_tmp;
    }

  }

  return m_close;

}

double DiLeptonAnalyzer::CorrPt(KLepton lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso()-T_iso)));
  return ptcorr;
}

double DiLeptonAnalyzer::CorrPt(snu::KMuon lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso04()-T_iso)));
  return ptcorr;
}

double DiLeptonAnalyzer::CorrPt(snu::KElectron lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.PFRelIso(0.3)-T_iso)));
  return ptcorr;
}



double DiLeptonAnalyzer::GetMuonFR(bool geterr, float pt,  float eta){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.4) eta = 2.3;

  int binx = hist_Muon_FR->FindBin(pt, abs(eta));
  //cout << "[DiLeptonAnalyzer::GetMuonFR] pt = " << pt << ", eta = " << eta << endl;
  //cout << "[DiLeptonAnalyzer::GetMuonFR] FR = " << hist_Muon_FR->GetBinContent(binx) << endl;
  //cout << "[DiLeptonAnalyzer::GetMuonFR] FR_err = " << hist_Muon_FR->GetBinError(binx) << endl;
  if(geterr) return hist_Muon_FR->GetBinError(binx);
  else return hist_Muon_FR->GetBinContent(binx);

}

double DiLeptonAnalyzer::GetMuonPR(bool geterr, float pt,  float eta){
/*
  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.5) eta = 2.4;

  int binx = hist_Muon_PR->FindBin(pt, abs(eta));
  if(geterr) return hist_Muon_PR->GetBinError(binx);
  else return hist_Muon_PR->GetBinContent(binx);
*/
  return 1.;
}

double DiLeptonAnalyzer::GetElectronFR(bool geterr, float pt,  float eta){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.5) eta = 2.4;

  int binx = hist_Electron_FR->FindBin(pt, abs(eta));
  if(geterr) return hist_Electron_FR->GetBinError(binx);
  else return hist_Electron_FR->GetBinContent(binx);

}

double DiLeptonAnalyzer::GetElectronPR(bool geterr, float pt,  float eta){
/*
  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.5) eta = 2.4;

  int binx = hist_Electron_PR->FindBin(pt, abs(eta));
  if(geterr) return hist_Electron_PR->GetBinError(binx);
  else return hist_Electron_PR->GetBinContent(binx);
*/
  return 1.;
}


void DiLeptonAnalyzer::get_eventweight(std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<bool> isT, int HalfSampleErrorDir){

  unsigned int n_leptons = isT.size();
  cout << "[DiLeptonAnalyzer::get_eventweight] muons.size() = " << muons.size() << ", electrons.size() = " << electrons.size() << endl;

  vector<float> lep_pt, lep_eta;
  vector<bool> ismuon;
  for(unsigned int i=0; i<muons.size(); i++){
    //lep_pt.push_back(muons.at(i).Pt());
    lep_pt.push_back( CorrPt(muons.at(i), 0.07) );
    lep_eta.push_back(muons.at(i).Eta());
    ismuon.push_back(true);
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    //lep_pt.push_back(electrons.at(i).Pt());
    lep_pt.push_back( CorrPt(electrons.at(i), 0.08) );
    lep_eta.push_back(electrons.at(i).Eta());
    ismuon.push_back(false);

  }

  vector<float> fr, pr, fr_err, pr_err;

  for(unsigned int i=0; i<n_leptons; i++){
    //==== Muon
    if(ismuon.at(i)){
      fr.push_back( GetMuonFR(false, lep_pt.at(i), lep_eta.at(i)));
      pr.push_back( GetMuonPR(false, lep_pt.at(i), lep_eta.at(i)));
      fr_err.push_back( GetMuonFR(true, lep_pt.at(i), lep_eta.at(i)));
      pr_err.push_back( GetMuonPR(true, lep_pt.at(i), lep_eta.at(i)));
    }
    //==== If not, it's an electron
    else{
      fr.push_back( GetElectronFR(0, lep_pt.at(i), lep_eta.at(i)));
      pr.push_back( GetElectronPR(0, lep_pt.at(i), lep_eta.at(i)));
      fr_err.push_back( GetElectronFR(1, lep_pt.at(i), lep_eta.at(i)));
      pr_err.push_back( GetElectronPR(1, lep_pt.at(i), lep_eta.at(i)));
    }
  }

  //==== if HalfSampleErrorDir!=0,
  //==== assign "5 %" HalfSampleTest Systematic Uncertainty

  if(HalfSampleErrorDir>0){
    for(unsigned int i=0; i<fr.size(); i++){
      fr.at(i)     = fr.at(i)    *(1.+0.05);
      fr_err.at(i) = fr_err.at(i)*(1.+0.05);
    }
  }
  else if(HalfSampleErrorDir<0){
    for(unsigned int i=0; i<fr.size(); i++){
      fr.at(i)     = fr.at(i)    *(1.-0.05);
      fr_err.at(i) = fr_err.at(i)*(1.-0.05);
    }
  }

  //==== let a == f/(1-f)

  vector<float> a, fr_onlyLoose;

  for(unsigned int i=0; i<n_leptons; i++) a.push_back( fr.at(i)/(1.-fr.at(i)) );
  for(unsigned int i=0; i<n_leptons; i++){
    if(!isT.at(i)){
      cout << "[DiLeptonAnalyzer::get_eventweight] "<<i<<" th lepton is Loose" << endl;
      fr_onlyLoose.push_back( a.at(i) );
    }
  }

  //==== Initialise weight
  float this_weight=-1.;

  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }
  cout << "[DiLeptonAnalyzer::get_eventweight] this_weight = " << this_weight << endl;

  //==== d(a)/a = d(f)/f(1-f)
  //==== so, if w = a1*a2,
  //==== d(w)/w = d(a1)/a1 + d(a2)/a2

  vector<float> da_over_a;
  for(unsigned int i=0; i<n_leptons; i++) da_over_a.push_back( fr_err.at(i) / ( fr.at(i)*(1.-fr.at(i)) ) );
  float this_weight_err = 0.;
  for(unsigned int i=0; i<n_leptons; i++){
    if(!isT.at(i)) this_weight_err += da_over_a.at(i)*da_over_a.at(i);
  }

  this_weight_err = sqrt(this_weight_err);
  this_weight_err = this_weight_err*fabs(this_weight);

  weight_fr = this_weight;
  weight_err_fr =  this_weight_err;


}

/*

    double cutop[100];
    cutop[] = lep.at(0).Pt();
    cutop[] = lep.at(1).Pt();
    cutop[] = lep.at(0).DeltaR( lep.at(1) );
    cutop[] = (lep.at(0)+lep.at(1)).M();
    cutop[] = isSS ? 0 : 1;
    cutop[] = isOffZ ? 0 : 1;

    cutop[] = jets.size();
    cutop[] = nbjets;
    cutop[] = jets_nolepveto.size();
    cutop[] = nbjets_nolepveto;
    cutop[] = jets_fwd.size();
    cutop[] = nbjets_fwd;

    //==== Two Jets
    if(jets.size() >= 2){

      cutop[] = jets.at(0).Pt();
      cutop[] = jets.at(1).Pt();
      cutop[] = jets.at(0).DeltaR( jets.at(1) );
      cutop[] = (jets.at(0)+jets.at(1)).M();
      cutop[] = (lep.at(0)+jets.at(0)+jets.at(1)).M();
      cutop[] = (lep.at(1)+jets.at(0)+jets.at(1)).M();
      cutop[] = (lep.at(0)+lep.at(1)+jets.at(0)+jets.at(1)).M();

      cutop[] = jets.at(index_jjW_j1).Pt();
      cutop[] = jets.at(index_jjW_j2).Pt();
      cutop[] = jets.at(index_jjW_j1).DeltaR( jets.at(index_jjW_j2) );
      cutop[] = (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
      cutop[] = (lep.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
      cutop[] = (lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
      cutop[] = (lep.at(0)+lep.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();

      cutop[] = jets.at(index_lljjW_j1).Pt();
      cutop[] = jets.at(index_lljjW_j2).Pt();
      cutop[] = jets.at(index_lljjW_j1).DeltaR( jets.at(index_lljjW_j2) );
      cutop[] = (jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
      cutop[] = (lep.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
      cutop[] = (lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
      cutop[] = (lep.at(0)+lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();

    }
    else{

    }

    //==== Two Forward Jets
    if(jets_fwd.size() >= 2){
      cutop[] = jets_fwd.at(0).DeltaR( jets_fwd.at(1) );
    }

    //==== Transverse Energies
    cutop[] = MET;
    cutop[] = ST;
    cutop[] = HT;
    cutop[] = LT;

*/


