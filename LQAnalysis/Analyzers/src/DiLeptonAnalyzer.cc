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
weight_fr(-999), weight_err_fr(-999), NTightLeptons(0),
MuFR_key(""), ElFR_key(""),
MET(-999), METphi(-999),
ST(-999), HT(-999), LT(-999), contramass(-999),
nbjets(-999), nbjets_fwd(-999), nbjets_nolepveto(-999), n_vtx(-999),
index_jjW_j1(-999), index_jjW_j2(-999),
index_lljjW_j1(-999), index_lljjW_j2(-999),
index_fjW(-999),
RunNtp(false),NowRunningCentral(true)
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

  //===============
  //==== Get Fake
  //===============

  //==== Muon

  TString MuonFRType_QCD = "v7_SIP3_";
  //TString MuonFRType_QCD = "HighdXY_Large_";

  //TString MuonFRType_Data = "v7_SIP3_";
  TString MuonFRType_Data = MuonFRType_QCD;

  MuonLooseID_loosest = "MUON_HN_LOOSEv7_SIP3_loosest";
  MuonVetoID_loosest = "MUON_HN_VETO_loosest";
  MuonTightID = "MUON_HN_TIGHT";

  bool UesTightv2 = false;
  if(MuonFRType_QCD.Contains("v9") || MuonFRType_Data.Contains("v9")){
    MuonLooseID_loosest = "MUON_HN_LOOSEv9_loosest";
    UesTightv2 = true;
  }
  if(MuonFRType_QCD.Contains("v10") || MuonFRType_Data.Contains("v10")){
    MuonLooseID_loosest = "MUON_HN_LOOSEv10_loosest";
    UesTightv2 = true;
  }
  if(MuonFRType_QCD.Contains("HighdXY") || MuonFRType_Data.Contains("HighdXY")){
    MuonLooseID_loosest = "MUON_HN_Loose_HighdXY_Small_loosest";
  }

  if(UesTightv2){
    MuonVetoID_loosest = "MUON_HN_VETOv2_loosest"; // chi2 loosen to 999
    MuonTightID = "MUON_HN_TIGHTv2";
  }

  //==== Electron

  TString ElectronFRType_QCD = "v7_";

  //TString ElectronFRType_Data = "v7_";
  TString ElectronFRType_Data = ElectronFRType_QCD;

  ElectronLooseID_loosest = "ELECTRON_HN_FAKELOOSEv7_loosest";
  ElectronVetoID_loosest = "ELECTRON_HN_VETO_loosest";
  ElectronTightID = "ELECTRON_HN_TIGHTv4";

  //==== Summary
  cout << "## Muon Fake ##" << endl;
  cout << "MuonFRType_QCD = " << MuonFRType_QCD << endl;
  cout << "MuonFRType_Data = " << MuonFRType_Data << endl;
  cout << "MuonVetoID_loosest = " << MuonVetoID_loosest << endl;
  cout << "MuonLooseID_loosest = " << MuonLooseID_loosest << endl;
  cout << "MuonTightID = " << MuonTightID << endl;

  //==== Read rootfiles

  TFile *file_Muon_FR         = new TFile( lqdir+"/JskimData/FR/Muon_Data_"+MuonFRType_Data+"FR.root");
  TFile *file_Muon_FR_QCD     = new TFile( lqdir+"/JskimData/FR/Muon_QCD_" +MuonFRType_QCD+"FR.root"); 
  TFile *file_Electron_FR     = new TFile( lqdir+"/JskimData/FR/Electron_Data_"+ElectronFRType_Data+"FR.root");
  TFile *file_Electron_FR_QCD = new TFile( lqdir+"/JskimData/FR/Electron_QCD_" +ElectronFRType_QCD+"FR.root");

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

  TString CantralAwayJet = "40";

  TString btag_wp = "_Loose";
  //TString btag_wp = "";

  hist_Muon_FR                 = (TH2D*)file_Muon_FR->Get("Muon_Data_"+MuonFRType_Data+"FR_Awayjet"+CantralAwayJet)->Clone();
  hist_Muon_FR_withbjet        = (TH2D*)file_Muon_FR->Get("Muon_Data_"+MuonFRType_Data+"FR_Awayjet"+CantralAwayJet+"_withbjet"+btag_wp)->Clone();
  hist_Muon_FR_withoutbjet     = (TH2D*)file_Muon_FR->Get("Muon_Data_"+MuonFRType_Data+"FR_Awayjet"+CantralAwayJet+"_withoutbjet"+btag_wp)->Clone();
  hist_Muon_FR_QCD             = (TH2D*)file_Muon_FR_QCD->Get("Muon_QCD_"+MuonFRType_QCD+"FR_Awayjet"+CantralAwayJet)->Clone();
  hist_Muon_FR_QCD_withbjet    = (TH2D*)file_Muon_FR_QCD->Get("Muon_QCD_"+MuonFRType_QCD+"FR_Awayjet"+CantralAwayJet+"_withbjet"+btag_wp)->Clone();
  hist_Muon_FR_QCD_withoutbjet = (TH2D*)file_Muon_FR_QCD->Get("Muon_QCD_"+MuonFRType_QCD+"FR_Awayjet"+CantralAwayJet+"_withoutbjet"+btag_wp)->Clone();

  hist_Electron_FR                 = (TH2D*)file_Electron_FR->Get("Electron_Data_"+ElectronFRType_Data+"FR_Awayjet"+CantralAwayJet)->Clone();
  hist_Electron_FR_withbjet        = (TH2D*)file_Electron_FR->Get("Electron_Data_"+ElectronFRType_Data+"FR_Awayjet"+CantralAwayJet+"_withbjet"+btag_wp)->Clone();
  hist_Electron_FR_withoutbjet     = (TH2D*)file_Electron_FR->Get("Electron_Data_"+ElectronFRType_Data+"FR_Awayjet"+CantralAwayJet+"_withoutbjet"+btag_wp)->Clone();
  hist_Electron_FR_QCD             = (TH2D*)file_Electron_FR_QCD->Get("Electron_QCD_"+ElectronFRType_QCD+"FR_Awayjet"+CantralAwayJet)->Clone();
  hist_Electron_FR_QCD_withbjet    = (TH2D*)file_Electron_FR_QCD->Get("Electron_QCD_"+ElectronFRType_QCD+"FR_Awayjet"+CantralAwayJet+"_withbjet"+btag_wp)->Clone();
  hist_Electron_FR_QCD_withoutbjet = (TH2D*)file_Electron_FR_QCD->Get("Electron_QCD_"+ElectronFRType_QCD+"FR_Awayjet"+CantralAwayJet+"_withoutbjet"+btag_wp)->Clone();

  //==== away pt
  TString awayjetpt[3] = {"20", "30", "60"};
  for(int i=0; i<3; i++){
    if(!MuonFRType_Data.Contains("HighdXY")) hist_Muon_FR_syst["Awatjet_"+awayjetpt[i]]     = (TH2D*)file_Muon_FR->Get("Muon_Data_"+MuonFRType_Data+"FR_Awayjet"+awayjetpt[i])->Clone();
    hist_Electron_FR_syst["Awatjet_"+awayjetpt[i]] = (TH2D*)file_Electron_FR->Get("Electron_Data_"+ElectronFRType_Data+"FR_Awayjet"+awayjetpt[i])->Clone();
  }

  file_Muon_FR->Close();
  file_Muon_FR_QCD->Close();
  file_Electron_FR->Close();
  file_Electron_FR_QCD->Close();

  delete file_Muon_FR;
  delete file_Muon_FR_QCD;
  delete file_Electron_FR;
  delete file_Electron_FR_QCD;

  origDir->cd();

  ForTree_PdfWeights = new vector<float>();
  ForTree_ScaleWeights = new vector<float>();

  return;
}


void DiLeptonAnalyzer::ExecuteEvents()throw( LQError ){

/*
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_NOCUT", true);
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  std::vector<snu::KJet> jets_nolepveto2 = GetJets("JET_HN_eta5_nolepveto", 20., 2.5);
  bool HasZeroIsolationMuon = false;
  FillHist("FakeMuon_Nvtx", eventbase->GetEvent().nVertices(), 1., 0., 50., 50);
  for(unsigned int i=0; i<testmuons.size(); i++){
    if(TruthMatched(testmuons.at(i))) continue;

    FillHist("FakeMuon_RelIso", testmuons.at(i).RelIso04(), 1., 0., 0.6, 60);
    FillHist("FakeMuon_dXY", testmuons.at(i).dXY(), 1., 0., 1.0, 1000);
    FillHist("FakeMuon_dZ", testmuons.at(i).dXY(), 1., 0., 1.0, 1000);


    if(testmuons.at(i).RelIso04() < 0.001){
      HasZeroIsolationMuon = true;

      FillHist("ZeroRelIso_dXY", testmuons.at(i).dXY(), 1., 0., 1.0, 1000);
      FillHist("ZeroRelIso_dZ", testmuons.at(i).dZ(), 1., 0., 1.0, 1000);

    }

  }

  if(HasZeroIsolationMuon){
    FillHist("ZeroRelIso_Nvtx", eventbase->GetEvent().nVertices(), 1., 0., 50., 50);
  }
  return;
*/
/*
   std::vector< snu::KElectron > testelectrons = GetElectrons(false, true, "ELECTRON_HN_FAKELOOSEv1");
   for(unsigned int i=0; i<testelectrons.size(); i++){
     if( TruthMatched(testelectrons.at(i),false) ) continue;

     FillHist("ElectronIsoTest_RelIso", testelectrons.at(i).PFRelIso(0.3), 1., 0., 0.6, 60);
     if(testelectrons.at(i).Pt() > 15) FillHist("ElectronIsoTest_Above15_RelIso", testelectrons.at(i).PFRelIso(0.3), 1., 0., 0.6, 60);
     else FillHist("ElectronIsoTest_Below15_RelIso", testelectrons.at(i).PFRelIso(0.3), 1., 0., 0.6, 60);

     if( (testelectrons.at(i).PFRelIso(0.3) < 0.01) ){
       FillHist("ElectronIsoTest_SmallISO_RelIso", testelectrons.at(i).PFRelIso(0.3), 1., 0., 0.01, 100);
       FillHist("ElectronIsoTest_SmallISO_Type", testelectrons.at(i).GetType(), 1., 0., 50., 50);
     }

   }


   return;
*/
/*
  std::vector< snu::KElectron > testelectrons = GetElectrons(false, true, "ELECTRON_HN_FAKELOOSEv2");
  //==== FRTEST
  for(unsigned int i=0; i<testelectrons.size(); i++){
    FillHist("TEST_ELECTRON_FR_TYPE_F0", testelectrons.at(i).GetType(), 1., 0., 50., 50);
    if(PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4")){
      FillHist("TEST_ELECTRON_FR_TYPE_F", testelectrons.at(i).GetType(), 1., 0., 50., 50);
    }

    TString EtaRegion = "InnerBarrel";
    if(fabs(testelectrons.at(i).SCEta()) > 1.479) EtaRegion = "EndCap";
    else if(fabs(testelectrons.at(i).SCEta()) > 0.8) EtaRegion = "OuterBarrel";
    else EtaRegion = "InnerBarrel";

    if(!TruthMatched(testelectrons.at(i), false)){
      FillHist("TEST_ELECTRON_FAKE_MVA_"+EtaRegion, testelectrons.at(i).MVA(), 1., -1., 1., 200);
    }
  }

  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_LOOSEv3", true);
  //std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_LOOSE", true);
  //==== FRTEST
  for(unsigned int i=0; i<testmuons.size(); i++){
    FillHist("TEST_MUON_FR_TYPE_F0", testmuons.at(i).GetType(), 1., 0., 50., 50);
    if(PassID(testmuons.at(i), "MUON_HN_TIGHT")){
      FillHist("TEST_MUON_FR_TYPE_F", testmuons.at(i).GetType(), 1., 0., 50., 50);
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
  std::vector<TString> triggerlist_DiMuon, triggerlist_DiElectron;

  std::vector<TString> triggerlist_EMu, triggerlist_EMu_PeriodBtoG, triggerlist_EMu_PeriodH, triggerlist_EMu_Mu8Ele23, triggerlist_EMu_Mu23Ele8;

  w_cutflow.clear();
  triggerlist_DiMuon.clear();
  triggerlist_DiElectron.clear();
  triggerlist_EMu.clear();

  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  //triggerlist_DiElectron.push_back("HLT_Ele27_WPTight_Gsf_v");

  triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  triggerlist_EMu_PeriodH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_EMu_PeriodH.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

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

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "DoMCClosure") != k_flags.end();
  bool AllLepton = std::find(k_flags.begin(), k_flags.end(), "AllLepton") != k_flags.end();
  bool KeepFakeLepton = false;
  if(DoMCClosure){
    KeepFakeLepton = true;
    MuFR_key = "QCD";
    ElFR_key = "QCD";
  }
  if(AllLepton) KeepFakeLepton = true;

  bool DoAwayJet = std::find(k_flags.begin(), k_flags.end(), "DoAwayJet") != k_flags.end();
  bool LooseSampleFakeJetPt = std::find(k_flags.begin(), k_flags.end(), "LooseSampleFakeJetPt") != k_flags.end();

  //==== Muons
  std::vector<snu::KMuon> muons_loosest = GetMuons(MuonLooseID_loosest, KeepFakeLepton);
  std::vector<snu::KMuon> muons_veto_loosest = GetMuons(MuonVetoID_loosest, true);

  //==== Electrons
  std::vector<snu::KElectron> electrons_loosest = GetElectrons(false, KeepFakeLepton, ElectronLooseID_loosest);
  std::vector<snu::KElectron> electrons_veto_loosest = GetElectrons(true, true, ElectronVetoID_loosest);

  bool RunningChargeFlipData = k_running_chargeflip && isData;

  //==== Jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto_loosest", 10., 5.);
  //std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJetsWFT("JET_HN_eta5_nolepveto_loosest", 10., 5.);

  //==== FatJets
  std::vector<snu::KFatJet> fatjets_loosest = GetFatJets("FATJET_HN_loosest");

/*
  //==== jet HadFlavour vs BTagging
  for(int j=0; j<jets_nolepveto.size(); j++){
    snu::KJet jet = jets_nolepveto.at(j);
    FillHist("TEST_JET_HadronFlavour", jet.HadronFlavour(), 1., 0., 10., 10);
    if( IsBTagged(jets_nolepveto.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      FillHist("TEST_JET_HadronFlavour_CSVv2_Medium", jet.HadronFlavour(), 1., 0., 10., 10);
    }
    if( IsBTagged(jets_nolepveto.at(j), snu::KJet::CSVv2, snu::KJet::Loose) ){
      FillHist("TEST_JET_HadronFlavour_CSVv2_Loose", jet.HadronFlavour(), 1., 0., 10., 10);
    }
  }
  return;
*/

  if(!isData){
    weight*=pileup_reweight;
    weight*=GetKFactor();
  }

  bool DoFRBJET = std::find(k_flags.begin(), k_flags.end(), "DoFRBJET") != k_flags.end();
  snu::KJet::WORKING_POINT btag_wp = snu::KJet::Loose;

  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  RunNtp = std::find(k_flags.begin(), k_flags.end(), "RunNtp") != k_flags.end();

  ForTree_PdfWeights->clear();
  ForTree_ScaleWeights->clear();
  snu::KEvent Evt = eventbase->GetEvent();
  for(unsigned int i=0; i<Evt.PdfWeights().size(); i++){
    ForTree_PdfWeights->push_back(Evt.PdfWeights().at(i));
  }
  for(unsigned int i=0; i<Evt.ScaleWeights().size(); i++){
    ForTree_ScaleWeights->push_back(Evt.ScaleWeights().at(i));
  }
  n_vtx = Evt.nVertices();

  //==== Define Analysis Region

  std::vector< TString > Suffixs;
  std::vector< std::vector<TString> > Triggers;

  bool RunningNonPromptData = NonPromptRun && (isData||DoMCClosure);

  //==== Make Suffix 

  Suffixs.push_back("DiMuon");
  Triggers.push_back(triggerlist_DiMuon);

  Suffixs.push_back("DiElectron");
  Triggers.push_back(triggerlist_DiElectron);

  Suffixs.push_back("EMu");
  Triggers.push_back(triggerlist_EMu);

  for(unsigned int i=0; i<Suffixs.size(); i++){

    TString OriginalSuffix = Suffixs.at(i);

    //==== Trigger pass
    if(!PassTriggerOR( Triggers.at(i) )) continue;

    double EMu_MCTriggerWeight = 0.;
    if(OriginalSuffix.Contains("EMu")){
      //==== Data
      if(isData){

        EMu_MCTriggerWeight = 1.;

        //==== Period BCDEFG
        if(1 <= GetDataPeriod() && GetDataPeriod() <= 6){
          if(!PassTriggerOR( triggerlist_EMu_PeriodBtoG )) continue;
        }
        //==== Period H
        else if(GetDataPeriod() == 7){
          if(!PassTriggerOR( triggerlist_EMu_PeriodH )) continue;
        }
        else{
          cout << "GetDataPeriod() Wrong... return" << endl;
          return;
        }
      }
      //==== MC
      else{
        if(PassTriggerOR( triggerlist_EMu_PeriodBtoG )) EMu_MCTriggerWeight += 27261.369;
        if(PassTriggerOR( triggerlist_EMu_PeriodH )) EMu_MCTriggerWeight += 8605.69;
      }
    }
    FillCutFlowByName(OriginalSuffix, "MET_PV_Trig", w_cutflow[OriginalSuffix], isData);

    //====================================
    //==== Systematic source loop starts
    //==== 1) Muon Energy Scale
    //==== 2) Jet Energy Scale
    //==== 3) Jet Energy Resolution
    //==== 4) Unclustered Energy
    //==== 5) Muon ID Scale Factor
    //==== 6) Pileup
    //==== 7) Trigger SF
    //==== 8) Electron ID Scale Factor
    //==== 9) Electron Energy Scale
    //==== 10) BTagSF Eff
    //==== 11) BTagSF Miss 
    //====================================

    int N_sys = 2*11+1;
    int it_sys_start = 0;
    if(!RunNtp){
      it_sys_start = 0;
      N_sys = it_sys_start+1;
    }

    for(int it_sys=it_sys_start; it_sys<N_sys; it_sys++){

      //==== MET
      //==== also set string for this systematic type
      MET = Evt.MET();
      METphi = Evt.METPhi();
      TString this_syst;
      if(it_sys==0){
        this_syst = ""; //Central
      }
      else if(it_sys==1){
        this_syst = "_MuonEn_up";
        //MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::up);
      }
      else if(it_sys==2){
        this_syst = "_MuonEn_down";
        //MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::down);
      }
      else if(it_sys==3){
        this_syst = "_JetEn_up";
        MET = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::up);
      }
      else if(it_sys==4){
        this_syst = "_JetEn_down";
        MET = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::down);
      }
      else if(it_sys==5){
        this_syst = "_JetRes_up";
        MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);
      }
      else if(it_sys==6){
        this_syst = "_JetRes_down";
        MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);
      }
      else if(it_sys==7){
        this_syst = "_Unclustered_up";
        MET = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
      }
      else if(it_sys==8){
        this_syst = "_Unclustered_down";
        MET = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
      }
      else if(it_sys==9){
        this_syst = "_MuonIDSF_up";
      }
      else if(it_sys==10){
        this_syst = "_MuonIDSF_down";
      }
      else if(it_sys==11){
        this_syst = "_PU_down";
      }
      else if(it_sys==12){
        this_syst = "_PU_up";
      }
      else if(it_sys==13){
        this_syst = "_TriggerSF_down";
      }
      else if(it_sys==14){
        this_syst = "_TriggerSF_up";
      }
      else if(it_sys==15){
        this_syst = "_ElectronIDSF_up";
      }
      else if(it_sys==16){
        this_syst = "_ElectronIDSF_down";
      }
      else if(it_sys==17){
        this_syst = "_ElectronEn_up";
      }
      else if(it_sys==18){
        this_syst = "_ElectronEn_down";
      }
      else if(it_sys==19){
        this_syst = "_BTagSFEff_up";
      }
      else if(it_sys==20){
        this_syst = "_BTagSFEff_down";
      }
      else if(it_sys==21){
        this_syst = "_BTagSFMiss_up";
      }
      else if(it_sys==22){
        this_syst = "_BTagSFMiss_down";
      }
      else{
        Message("it_sys out of range!" , INFO);
        return;
      }

      if(this_syst==""){
        NowRunningCentral = true;
      }
      else{
        NowRunningCentral = false;
      }

      TString Suffix = OriginalSuffix+this_syst;

      //================
      //==== Make Muon
      //================

      double this_MuonLooseRelIso = 0.4;
      double this_MuonVetoRelIso = 0.6;
      std::vector<snu::KMuon> muons, muons_veto;
      if(this_syst == "_MuonEn_up"){
        //==== Signal Leptons
        for(unsigned int j=0; j<muons_loosest.size(); j++){
          snu::KMuon this_muon = muons_loosest.at(j);
          this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
          double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedUp();
          this_muon.SetRelIso(0.4, new_RelIso);
          if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
          snu::KMuon this_muon = muons_veto_loosest.at(j);
          this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
          double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedUp();
          this_muon.SetRelIso(0.4, new_RelIso);
          if( this_muon.Pt() >= 5. && new_RelIso < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
        }
      }
      else if(this_syst == "_MuonEn_down"){
        //==== Signal Leptons
        for(unsigned int j=0; j<muons_loosest.size(); j++){
          snu::KMuon this_muon = muons_loosest.at(j);
          this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
          double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedDown();
          this_muon.SetRelIso(0.4, new_RelIso);
          if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
          snu::KMuon this_muon = muons_veto_loosest.at(j);
          this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
          double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedDown();
          this_muon.SetRelIso(0.4, new_RelIso);
          if( this_muon.Pt() >= 5. && new_RelIso < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
        }
      }
      //==== normal muons
      else{
        //==== Signal Leptons
        for(unsigned int j=0; j<muons_loosest.size(); j++){
          snu::KMuon this_muon = muons_loosest.at(j);
          if( this_muon.Pt() >= 10. && this_muon.RelIso04() < this_MuonLooseRelIso ) muons.push_back( this_muon );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
          snu::KMuon this_muon = muons_veto_loosest.at(j);
          if( this_muon.Pt() >= 5. && this_muon.RelIso04() < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
        }
      }
      std::sort(muons.begin(), muons.end(), MuonPtComparing);

      //====================
      //==== Make Electron
      //====================

      std::vector<snu::KElectron> electrons, electrons_veto;
      std::vector<snu::KElectron> electrons_notshifted;
      double this_ElectronLooseRelIso = 0.6;
      double this_ElectronVetoRelIso = 0.6;
      int ElEnDir = 0;
      if(this_syst == "_ElectronEn_up"){
        ElEnDir = 1;
        //==== Signal Leptons
        for(unsigned int j=0; j<electrons_loosest.size(); j++){
          snu::KElectron this_electron = electrons_loosest.at(j);
          this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
          double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
          this_electron.SetPFRelIso(0.3, new_RelIso);
          if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
          snu::KElectron this_electron = electrons_veto_loosest.at(j);
          this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
          double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
          this_electron.SetPFRelIso(0.3, new_RelIso);
          if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
        }
      }
      else if(this_syst == "_ElectronEn_down"){
        ElEnDir = -1;
        //==== Signal Leptons
        for(unsigned int j=0; j<electrons_loosest.size(); j++){
          snu::KElectron this_electron = electrons_loosest.at(j);
          this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
          double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
          this_electron.SetPFRelIso(0.3, new_RelIso);
          if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
          snu::KElectron this_electron = electrons_veto_loosest.at(j);
          this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
          double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
          this_electron.SetPFRelIso(0.3, new_RelIso);
          if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
        }
      }
      //==== normal electrons
      else{
        //==== Signal Leptons
        for(unsigned int j=0; j<electrons_loosest.size(); j++){
          snu::KElectron this_electron = electrons_loosest.at(j);
          if( this_electron.Pt() >= 10. && this_electron.PFRelIso(0.3) < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
        }
        //==== Veto Leptons
        for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
          snu::KElectron this_electron = electrons_veto_loosest.at(j);
          if( this_electron.Pt() >= 10. && this_electron.PFRelIso(0.3) < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
        }
      }

      //==== If Chargeflip, then shift down electron, and replace electrons with it
      //==== To get correct CF rate, we save original (OS event) electron to electrons_before_shift
      std::vector< snu::KElectron > electrons_before_shift;
      electrons_before_shift.clear();
      if(RunningChargeFlipData){
        for(unsigned int j=0; j<electrons.size(); j++){
          electrons_before_shift.push_back( electrons.at(j) );
          snu::KElectron tmp_el = electrons.at(j);
          double shift_ = 1.-0.015;
          tmp_el.SetPtEtaPhiM(shift_*tmp_el.Pt(), tmp_el.Eta(), tmp_el.Phi(), tmp_el.M());
          //tmp_el.SetPxPyPzE(shift_*tmp_el.Px(), shift_*tmp_el.Py(), shift_*tmp_el.Pz(), shift_*tmp_el.E());
          electrons.at(j) = tmp_el;
        }
      }

      //==== We don't have MiniAODPt for electorn..

      if(ElEnDir!=0){
        for(unsigned int j=0; j<electrons.size(); j++){
          snu::KElectron new_electron = electrons.at(j);
          if(ElEnDir>0) new_electron.SetPtEtaPhiM( new_electron.Pt()/new_electron.PtShiftedUp(), new_electron.Eta(), new_electron.Phi(), new_electron.M() );
          else new_electron.SetPtEtaPhiM( new_electron.Pt()/new_electron.PtShiftedDown(), new_electron.Eta(), new_electron.Phi(), new_electron.M() );
          electrons_notshifted.push_back(new_electron);
        }
      }


      JSCorrectedMETElectron(ElEnDir, electrons_notshifted, MET, METphi);
      std::sort(electrons.begin(), electrons.end(), ElectronPtComparing);

      //==== MET is calculated with No-Rochestor-Corrected Muons
      //==== In this step, muons are 
      //==== 1) Rochestor corrected & Up/Down
      //==== 2) Rochestor corrected
      //==== Both cases, we can correct MET (w.r.t. muon) using
      //==== AnalyzerCore::JSCorrectedMETRochester(std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi)
      JSCorrectedMETRochester(muons, MET, METphi);

      //==================
      //==== Make FatJet
      //==================

      vector<snu::KFatJet> fatjets;
      if( (this_syst == "_JetEn_up") || (this_syst == "_JetRes_up") ){
        for(unsigned int j=0; j<fatjets_loosest.size(); j++){
          snu::KFatJet this_jet = fatjets_loosest.at(j);

          double this_scaling = 1.;
          if(this_syst == "_JetEn_up") this_scaling = this_jet.ScaledUpEnergy();
          if(this_syst == "_JetRes_up") this_scaling = this_jet.SmearedResUp()/this_jet.SmearedRes();
          this_jet *= this_scaling;

          if(JSFatJetID(this_jet)) fatjets.push_back(this_jet);
        }
      }
      else if( (this_syst == "_JetEn_down") || (this_syst == "_JetRes_down") ){
        for(unsigned int j=0; j<fatjets_loosest.size(); j++){
          snu::KFatJet this_jet = fatjets_loosest.at(j);

          double this_scaling = 1.;
          if(this_syst == "_JetEn_down") this_scaling = this_jet.ScaledDownEnergy();
          if(this_syst == "_JetRes_down") this_scaling = this_jet.SmearedResDown()/this_jet.SmearedRes();;
          this_jet *= this_scaling;

          if(JSFatJetID(this_jet)) fatjets.push_back(this_jet);
        }
      }
      else{
        for(unsigned int j=0; j<fatjets_loosest.size(); j++){
          snu::KFatJet this_jet = fatjets_loosest.at(j);

          if(k_sample_name.Contains("HN")){ //FIXME
            FillHist(Suffix+"_"+"BeforeIDCut_FatJetPrunedMass", this_jet.PrunedMass(), 1., 0., 500., 500);
          }

          if(JSFatJetID(this_jet)) fatjets.push_back(this_jet);
        }
      }
      //cout << "fatjets.size() = " << fatjets.size() << endl;

      //===============
      //==== Make Jet
      //===============

      std::vector<snu::KJet> jets_eta5_nolepveto; // eta < 5, NO lepton-veto
      if( (this_syst == "_JetEn_up") || (this_syst == "_JetRes_up") ){
        for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
          snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

          double this_scaling = 1.;
          if(this_syst == "_JetEn_up") this_scaling = this_jet.ScaledUpEnergy();
          if(this_syst == "_JetRes_up") this_scaling = this_jet.SmearedResUp()/this_jet.SmearedRes();

          //double this_E = this_jet.E()*this_scaling;
          //double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
          //double this_3p_sf = this_3p/this_jet.P();
          //this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);

          this_jet *= this_scaling;

          if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
        }
      }
      else if( (this_syst == "_JetEn_down") || (this_syst == "_JetRes_down") ){
        for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
          snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

          double this_scaling = 1.;
          if(this_syst == "_JetEn_down") this_scaling = this_jet.ScaledDownEnergy();
          if(this_syst == "_JetRes_down") this_scaling = this_jet.SmearedResDown()/this_jet.SmearedRes();;

          //double this_E = this_jet.E()*this_scaling;
          //double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
          //double this_3p_sf = this_3p/this_jet.P();
          //this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);

          this_jet *= this_scaling;

          if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
        }
      }
      else{
        for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
          snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);
          if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
        }
      }

      std::vector<snu::KJet> jets; // eta < 2.5, lepton-veto, away from fatjets
      std::vector<snu::KJet> jets_nolepveto; // eta < 2.5, NO lepton-veto
      std::vector<snu::KJet> jets_fwd; // 2.5 < eta < 5, lepton-veto => to make forward

      for(unsigned int j=0; j<jets_eta5_nolepveto.size(); j++){

        snu::KJet this_jet = jets_eta5_nolepveto.at(j);
        bool IsNormalJet = fabs( this_jet.Eta() ) < 2.5;
        bool IsForwardJet = fabs( this_jet.Eta() ) >= 2.5;
        bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
        bool awayfromfatjet = IsAwayFromFatJet(this_jet, fatjets);

        if(IsNormalJet && !lepinside && awayfromfatjet) jets.push_back( this_jet );
        if(IsNormalJet) jets_nolepveto.push_back( this_jet );
        if(IsForwardJet && !lepinside) jets_fwd.push_back( this_jet );

      }

      int BTagSFDir = 0;
      if(this_syst == "_BTagSFEff_up") BTagSFDir = +1;
      if(this_syst == "_BTagSFEff_down") BTagSFDir = -1;
      if(this_syst == "_BTagSFMiss_up") BTagSFDir = +3;
      if(this_syst == "_BTagSFMiss_down") BTagSFDir = -3;

      nbjets = 0;
      for(int j=0; j<jets.size(); j++){
        if( IsBTagged(jets.at(j), snu::KJet::CSVv2, snu::KJet::Medium, -1, BTagSFDir) && fabs(jets.at(j).Eta())<2.4 ){
          nbjets++;
        }
      }

      nbjets_nolepveto = 0;
      for(int j=0; j<jets_nolepveto.size(); j++){
        if( IsBTagged(jets_nolepveto.at(j), snu::KJet::CSVv2, snu::KJet::Medium, -1, BTagSFDir) && fabs(jets_nolepveto.at(j).Eta())<2.4 ){
          nbjets_nolepveto++;
        }
      }

      nbjets_fwd = 0;
      for(int j=0; j<jets_fwd.size(); j++){
        if( IsBTagged(jets_fwd.at(j), snu::KJet::CSVv2, snu::KJet::Medium, -1, BTagSFDir) && fabs(jets_fwd.at(j).Eta())<2.4 ){
          nbjets_fwd++;
        }
      }

      //==== Lepton Numbers

      std::vector<snu::KMuon> muons_tight; muons_tight.clear();
      std::vector<snu::KElectron> electrons_tight; electrons_tight.clear();
      std::vector<bool> isT;
      std::vector<int> NearBjet;

      int NPromptTight(0);

      int n_veto_muons = muons_veto.size();
      int n_triLoose_muons = muons.size();
      int n_triTight_muons(0);
      for(unsigned int j=0; j<muons.size(); j++){

        if(DoFRBJET) NearBjet.push_back( HasCloseBjet(muons.at(j), jets_nolepveto, btag_wp) );
        else NearBjet.push_back( -1 );

        if(TruthMatched(muons.at(j))) NPromptTight++;
        if(PassID(muons.at(j), MuonTightID)){
          isT.push_back(true);
          muons_tight.push_back(muons.at(j));
          n_triTight_muons++;
        }
        else{
          isT.push_back(false);
        }
      }

      int n_veto_electrons = electrons_veto.size();
      int n_triLoose_electrons = electrons.size();
      int n_triTight_electrons(0);
      for(unsigned int j=0; j<electrons.size(); j++){

        if(DoFRBJET) NearBjet.push_back( HasCloseBjet(electrons.at(j), jets_nolepveto, btag_wp) );
        else NearBjet.push_back( -1 );

        if(TruthMatched(electrons.at(j), false)) NPromptTight++;
        if(PassID(electrons.at(j), ElectronTightID)){
          isT.push_back(true);
          electrons_tight.push_back(electrons.at(j));
          n_triTight_electrons++;
        }
        else{
          isT.push_back(false);
        }
      }

      int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
      int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;
      NTightLeptons = n_triTight_leptons;

      bool isTwoMuon_TT    = (n_triTight_muons == 2); // veto third later
      bool isTwoMuon_Loose = (n_triLoose_muons == 2 && n_triTight_muons != 2);

      bool isTwoElectron_TT    = (n_triTight_electrons == 2); // veto third later 
      bool isTwoElectron_Loose = (n_triLoose_electrons == 2 && n_triTight_electrons != 2);

      bool isEMu_TT    = (n_triTight_muons == 1 && n_triTight_electrons == 1);
      bool isEMu_Loose = (n_triLoose_muons == 1)     &&
                         (n_triLoose_electrons == 1) &&
                         (n_triTight_leptons != 2);

      bool isTT, isLOOSE;
      bool isNoExtra, isNoExtraOtherFlavour;

      if(Suffix.Contains("DiMuon")){
        isTT = isTwoMuon_TT && !RunningNonPromptData;
        isLOOSE = isTwoMuon_Loose && RunningNonPromptData;
        isNoExtra = n_veto_muons == 2;
        isNoExtraOtherFlavour = n_veto_electrons == 0;
      }
      if(Suffix.Contains("DiElectron")){
        isTT = isTwoElectron_TT && !RunningNonPromptData;
        isLOOSE = isTwoElectron_Loose && RunningNonPromptData;
        isNoExtra = n_veto_electrons == 2;
        isNoExtraOtherFlavour = n_veto_muons == 0;
      }
      if(Suffix.Contains("EMu")){
        isTT = isEMu_TT && !RunningNonPromptData;
        isLOOSE = isEMu_Loose && RunningNonPromptData;
        isNoExtra = (n_veto_muons == 1) && (n_veto_electrons == 1);
        isNoExtraOtherFlavour = true;
      }
      
      //==== Two leptons
      if(!isTT && !isLOOSE) continue;

      //==== That two lepton pass basic pt cuts
      if(Suffix.Contains("DiMuon")){
        //==== to properly veto below tricky event
        //==== ### failing third muon veto ###
        //==== muons.size() = 3
        //==== muons.at(0).Pt() = 51.8417 => isTight = 0
        //==== muons.at(1).Pt() = 27.8285 => isTight = 1
        //==== muons.at(2).Pt() = 8.72782 => isTight = 1
        if(isTT){
          if((muons_tight.at(0).Pt() < 20.) || (muons_tight.at(1).Pt() < 10.)) continue;
        }
        else{
          if((muons.at(0).Pt() < 20.) || (muons.at(1).Pt() < 10.)) continue;
        }
      }
      if(Suffix.Contains("DiElectron")){
        if(isTT){
          if(electrons_tight.at(0).Pt() < 25. || electrons_tight.at(1).Pt() < 15.) continue;
        }
        else{
          if(electrons.at(0).Pt() < 25. || electrons.at(1).Pt() < 15.) continue;
        }
      }
      if(Suffix.Contains("EMu")){
        double MuMinPt = 9999., ElMinPt = 9999.;

        bool PtOkay = false;

        if(PassTriggerOR(triggerlist_EMu_Mu8Ele23)){

          MuMinPt = 10.;
          ElMinPt = 25.;

          if(isTT){
            if( (muons_tight.at(0).Pt() > MuMinPt) && (electrons_tight.at(0).Pt() > ElMinPt) ) PtOkay = true;
          }
          else{
            if( (muons.at(0).Pt() > MuMinPt) && (electrons.at(0).Pt() > ElMinPt) ) PtOkay = true;
          }

        }
        if(PassTriggerOR(triggerlist_EMu_Mu23Ele8)){

          MuMinPt = 25.;
          ElMinPt = 10.;

          if(isTT){
            if( (muons_tight.at(0).Pt() > MuMinPt) && (electrons_tight.at(0).Pt() > ElMinPt) ) PtOkay = true;
          }
          else{ 
            if( (muons.at(0).Pt() > MuMinPt) && (electrons.at(0).Pt() > ElMinPt) ) PtOkay = true;
          }


        }

        if( !PtOkay ) continue;

      }

      FillCutFlowByName(Suffix, "TwoLeptons", w_cutflow[Suffix], isData);

      //==== No Extra lepton
      if(!isNoExtra) continue;
      FillCutFlowByName(Suffix, "NoExtraLepton", w_cutflow[Suffix], isData);

      //==== No Extra different flavour lepton
      if(!isNoExtraOtherFlavour) continue;
      FillCutFlowByName(Suffix, "NoExtraFlavourLepton", w_cutflow[Suffix], isData);

      //==== DiMuon-DoubleMuon PD / ...
      if(isData && k_channel != "DoubleMuon_CF"){
        if(Suffix.Contains("DiMuon")){
          if(k_channel != "DoubleMuon") continue;
        }
        if(Suffix.Contains("DiElectron")){
          if(k_channel != "DoubleEG") continue;
        }
        if(Suffix.Contains("EMu")){
          if(k_channel != "MuonEG") continue;
        }
      }

      double trigger_ps_weight(1.);
      if(Suffix.Contains("EMu")) trigger_ps_weight = EMu_MCTriggerWeight;
      else              trigger_ps_weight = WeightByTrigger(Triggers.at(i), TargetLumi);

      double this_weight = weight*trigger_ps_weight;

      //double muon_id_iso_sf = mcdata_correction->MuonScaleFactor(MuonTightID, muons, 0); //FIXME
      double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", muons, 0);
      double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muons);
      //muon_id_iso_sf = 1.; //FIXME
      this_weight *= muon_id_iso_sf*MuTrkEffSF;

      double electron_sf = mcdata_correction->ElectronScaleFactor(ElectronTightID, electrons, 0);
      //electron_sf = 1.; //FIXME
      double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);
      this_weight *= electron_sf*electron_RecoSF;

      double trigger_sf = 1.;
      int TriggerSFDir = 0;
      if(this_syst == "_TriggerSF_up"){
        TriggerSFDir = +1;
      }
      else if(this_syst == "_TriggerSF_up"){
        TriggerSFDir = -1;
      }
      else{
        TriggerSFDir = 0;
      }


      if(!isData && Suffix.Contains("DiMuon")){
        //double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, MuonTightID, 0, 0, TriggerSFDir); //FIXME
        //double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, MuonTightID, 0, 1, -1*TriggerSFDir);

        double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 0, 0, TriggerSFDir);
        double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 0, 1, -1*TriggerSFDir);
        trigger_sf = trigger_eff_Data/trigger_eff_MC;
      }
      if(!isData && Suffix.Contains("DiElectron")){

        double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, MuonTightID, 1, 0, TriggerSFDir);
        double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, MuonTightID, 1, 1, -1*TriggerSFDir);
        trigger_sf = trigger_eff_Data/trigger_eff_MC;

        //cout << "## Calculating DiElectron Trigger SF ##" << endl;
        //cout << "1) Data" << endl;
        //cout << "trigger_eff_Data = " << trigger_eff_Data << endl;
        //cout << "2) MC" << endl;
        //cout << "trigger_eff_MC = " << trigger_eff_MC << endl;
        //cout << "=> sf = " << trigger_sf << endl;

      }
      //trigger_sf = 1.;//FIXME

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

      if(DoMCClosure){

        bool HasStrangeLeptons=false;
        bool HasFakeLeptons=false;

        for(unsigned int j=0; j<lep.size(); j++){
  /*
          int this_type = lep.at(j).GetType();
          if(lep.at(j).LeptonFlavour()==KLepton::ELECTRON){
            //if(this_type==2) HasStrangeLeptons=true; // From Z/W, but matchedpdgid != pdgid (conversion) //TODO test this
            if(this_type==3) HasStrangeLeptons=true; // e>eee
            if(this_type==16) HasStrangeLeptons=true; // gammastar
            if(this_type==40) HasStrangeLeptons=true; // photon matched
          }
          if(lep.at(j).LeptonFlavour()==KLepton::MUON){
            if(this_type==8) HasStrangeLeptons=true; // ??
            if(this_type==10) HasStrangeLeptons=true; // gammastar
          }
          if(lep.at(j).MCIsCF()) HasStrangeLeptons=true;
  */

          if( lep.at(j).GetType() == 0 ) HasStrangeLeptons = true;

          //==== Check if fake exist
          if(lep.at(j).LeptonFlavour()==KLepton::MUON){
            const snu::KMuon *mptr = lep.at(j).GetMuonPtr();
            if(!TruthMatched(*mptr)) HasFakeLeptons = true;
          }
          if(lep.at(j).LeptonFlavour()==KLepton::ELECTRON){
            const snu::KElectron *eptr = lep.at(j).GetElectronPtr();
            if(!TruthMatched(*eptr,false)) HasFakeLeptons = true;
          }


        } // END lepton loop


        if(HasStrangeLeptons) continue;
        if(!HasFakeLeptons) continue;

      }

      bool isSS = lep.at(0).Charge() == lep.at(1).Charge();
      bool isSSForCF = isSS;
      if(RunningChargeFlipData) isSSForCF = !isSS;

      double m_Z = 91.1876;
      bool isOffZ = fabs( (lep.at(0)+lep.at(1)).M() - m_Z ) > 10.;
      bool isAboveZ = ( (lep.at(0)+lep.at(1)).M() - m_Z ) > 10.;
      bool isBelowZ = ( (lep.at(0)+lep.at(1)).M() - m_Z ) < -10.;

  /*
      //==== TT/DY Same-Sign Prompt check
      if(DoMCClosure && NPromptTight==2 && isSS){
        FillHist("SSType_"+Suffix, lep.at(0).GetType(), 1., 0., 40., 40);
        FillHist("SSType_"+Suffix, lep.at(1).GetType(), 1., 0., 40., 40);
        FillHist("SSIsCF_"+Suffix, lep.at(0).MCIsCF(), 1., 0., 2., 2);
        FillHist("SSIsCF_"+Suffix, lep.at(1).MCIsCF(), 1., 0., 2., 2);
        cout << "########### SS Two Prompt ###########" << endl;
        cout << "Channel = " << Suffix << endl;
        cout << "<Truth>" << endl;
        TruthPrintOut();
        cout << "<Reco>" << endl;
        cout << "pt\teta\tphi\tcharge\ttype\n" << endl;
        for(unsigned int j=0; j<lep.size(); j++){
          cout << lep.at(j).Pt() << "\t" << lep.at(j).Eta() << "\t" << lep.at(j).Phi() << "\t" << lep.at(j).Charge() << "\t" << lep.at(j).GetType() << endl;
        }
      }
      continue;
  */

      //==== mll Cut Study
      if(this_syst == ""){
        FillHist("CutStudy_m_ll_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
        if(isSS) FillHist("CutStudy_m_ll_SS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
        else FillHist("CutStudy_m_ll_OS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);

        if( lep.at(0).DeltaR( lep.at(1) ) < 0.1 ){

          FillHist("CutStudy_m_ll_dR0p1_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
          if(isSS) FillHist("CutStudy_m_ll_dR0p1_SS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
          else FillHist("CutStudy_m_ll_dR0p1_OS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);

        }

      }
      bool mll10GeV = ( lep.at(0)+lep.at(1) ).M() < 10.;
      if(mll10GeV) continue;
      FillCutFlowByName(Suffix, "LowDileptonMass", w_cutflow[Suffix], isData);

      double this_weight_err(0.);
      vector<double> FRweights;
      vector<TString> FRweights_name;
      if( isLOOSE ){

        get_eventweight(muons, electrons, isT, NearBjet, 0);
        this_weight *= weight_fr;
        this_weight_err = this_weight*weight_err_fr;

        if(DoAwayJet){

          FRweights.push_back(weight_fr);
          FRweights_name.push_back("FakeRate_ptcorr_eta40");

          ElFR_key = "Awatjet_20";
          MuFR_key = "Awatjet_20";
          get_eventweight(muons, electrons, isT, NearBjet, 0);
          FRweights.push_back(weight_fr);
          FRweights_name.push_back("FakeRate_ptcorr_eta20");

          ElFR_key = "Awatjet_30";
          MuFR_key = "Awatjet_30";
          get_eventweight(muons, electrons, isT, NearBjet, 0);
          FRweights.push_back(weight_fr);
          FRweights_name.push_back("FakeRate_ptcorr_eta30");

          ElFR_key = "Awatjet_60";
          MuFR_key = "Awatjet_60";
          get_eventweight(muons, electrons, isT, NearBjet, 0);
          FRweights.push_back(weight_fr);
          FRweights_name.push_back("FakeRate_ptcorr_eta60");


          //==== reset
          ElFR_key = "";
          MuFR_key = "";
        }


      }
      if( RunningChargeFlipData && Suffix.Contains("DiElectron") && !isSS ){
        GetCFWeight(electrons_before_shift.at(0), electrons_before_shift.at(1));

        this_weight *= weight_cf;
        this_weight_err = this_weight*weight_err_cf;
      }

      //==== Now,
      //==== Fill Histogram
      //====

      //==== We have to sort lepton after we get fake weight,
      //==== because isT = {muon, electron}
      //==== If MCClosure, keep ordering as Muon-Electron (to see type)
      if(!DoMCClosure) std::sort(lep.begin(), lep.end(), LeptonPtComparing);

      std::map< TString, bool > map_Region_to_Bool;
      map_Region_to_Bool.clear();

      //==== Default Variation
      //==== 1) OnZ/OffZ/All
      //==== 2) OS/SS/All

      //==== # of jets
      map_Region_to_Bool[Suffix+"_0jets"] = (jets.size()==0);
      map_Region_to_Bool[Suffix+"_1jets"] = (jets.size()==1);
      map_Region_to_Bool[Suffix+"_Inclusive2jets"] = (jets.size()>=2);
      map_Region_to_Bool[Suffix+"_1jets_0nlbjets"] = (jets.size()==1) && (nbjets_nolepveto==0);
      map_Region_to_Bool[Suffix+"_1jets_0nlbjets_METge50"] = (jets.size()==1) && (nbjets_nolepveto==0) && (MET>50.);
      map_Region_to_Bool[Suffix+"_1jets_0nlbjets_mllge110"] = (jets.size()==1) && (nbjets_nolepveto==0) && (( lep.at(0)+lep.at(1) ).M() >= 110.);
      map_Region_to_Bool[Suffix+"_Inclusive1fatjets"] = (fatjets.size()>=1);

      //==== # of forward b jet
      map_Region_to_Bool[Suffix+"_0nlbjets"] = (nbjets_nolepveto==0);
      map_Region_to_Bool[Suffix+"_Inclusive1nlbjets"] = (nbjets_nolepveto>=1);
      map_Region_to_Bool[Suffix+"_Inclusive1nlbjets_METge50"] = (nbjets_nolepveto>=1) && (MET>50.); // control region

      //==== ST = lepton + jet + MET
      ST = lep.at(0).Pt() + lep.at(1).Pt() + MET;
      for(unsigned int ij=0; ij <jets.size(); ij++){
        ST += jets.at(ij).Pt();
      }
      for(unsigned int ij=0; ij <fatjets.size(); ij++){
        ST += fatjets.at(ij).Pt();
      }

      //==== HT = jet
      HT = 0.;
      for(unsigned int ij=0; ij <jets.size(); ij++){
        HT += jets.at(ij).Pt();
      }
      for(unsigned int ij=0; ij <fatjets.size(); ij++){
        HT += fatjets.at(ij).Pt();
      }

      //==== LT = lepton
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

      }
      if( fatjets.size()>=1 ){
        double mfatjet = GetFatJetMassClosest(fatjets, 80.4, index_fjW);
      }


      //==== Preselection

      bool TwoJet_NoFatJet = (jets.size()>=2) && (fatjets.size()==0);
      bool OneJet_NoFatJet = (jets.size()==1) && (fatjets.size()==0) && ( (lep.at(0)+lep.at(1)).M() < 80 ); // has m(ll) < 80 GeV
      bool OneFatJet       =                     (fatjets.size()>=1);

      map_Region_to_Bool[Suffix+"_Preselection"] = TwoJet_NoFatJet || OneJet_NoFatJet || OneFatJet;
      map_Region_to_Bool[Suffix+"_Preselection_secondptge20"] = map_Region_to_Bool[Suffix+"_Preselection"] && (lep.at(1).Pt() >= 20.);
      //==== For DiElectron, remove Z peak (CF)
      if(Suffix.Contains("DiElectron")){
        map_Region_to_Bool[Suffix+"_Preselection"] = map_Region_to_Bool[Suffix+"_Preselection"] && isOffZ;
      }

      //==== If Preselection, then define Low/High

      if(map_Region_to_Bool[Suffix+"_Preselection"]){

        //==== Low Mass
        //==== 1) TwoJet_NoFatJet
        //==== 2) OneJet_NoFatjet (has mll<80GeV cut)

        double mlljj_low = -999.;
        if(TwoJet_NoFatJet) mlljj_low = (lep.at(0)+lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        else if(OneJet_NoFatJet) mlljj_low = (lep.at(0)+lep.at(1)+jets.at(0)).M();

        map_Region_to_Bool[Suffix+"_Low"] = (TwoJet_NoFatJet || OneJet_NoFatJet) &&
                                            (nbjets_nolepveto == 0) &&
                                            ( mlljj_low < 300.) &&
                                            (MET < 80.);
        map_Region_to_Bool[Suffix+"_LowCR"] = (TwoJet_NoFatJet || OneJet_NoFatJet) &&
                                              ( (nbjets_nolepveto >= 1) || (MET > 100.) ) &&
                                              ( mlljj_low < 300.);

        map_Region_to_Bool[Suffix+"_Low_TwoJet_NoFatJet"] = TwoJet_NoFatJet &&
                                                            (nbjets_nolepveto == 0) &&
                                                            ( mlljj_low < 300.) &&
                                                            (MET < 80.);
        map_Region_to_Bool[Suffix+"_LowCR_TwoJet_NoFatJet"] = TwoJet_NoFatJet &&
                                                              ( (nbjets_nolepveto >= 1) || (MET > 100.) ) &&
                                                              ( mlljj_low < 300.);

        map_Region_to_Bool[Suffix+"_Low_OneJet_NoFatJet"] = OneJet_NoFatJet &&
                                                   (nbjets_nolepveto == 0) &&
                                                   ( mlljj_low < 300.) &&
                                                   (MET < 80.);
        map_Region_to_Bool[Suffix+"_LowCR_OneJet_NoFatJet"] = OneJet_NoFatJet &&
                                                     ( (nbjets_nolepveto >= 1) || (MET > 100.) ) &&
                                                     ( mlljj_low < 300.);

        //==== High Mass
        //==== 1) TwoJet_NoFatJet
        //==== 2) OneFatJet

        double mjj_high = -999.;
        if(TwoJet_NoFatJet){
          mjj_high = (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M();
        }
        else if(OneFatJet){
          mjj_high = fatjets.at(index_fjW).PrunedMass();
        }
        else{
          
        }

        map_Region_to_Bool[Suffix+"_High"] = ( TwoJet_NoFatJet || OneFatJet ) &&
                                             (nbjets_nolepveto == 0) &&
                                             ( mjj_high < 150. ) &&
                                             ( MET*MET/ST < 15. );
        map_Region_to_Bool[Suffix+"_HighCR"] = ( TwoJet_NoFatJet || OneFatJet ) &&
                                               ( (nbjets_nolepveto >= 1) || (MET*MET/ST > 20.) ) &&
                                               ( mjj_high < 150. );

        map_Region_to_Bool[Suffix+"_High_TwoJet_NoFatJet"]   = map_Region_to_Bool[Suffix+"_High"] && TwoJet_NoFatJet;
        map_Region_to_Bool[Suffix+"_HighCR_TwoJet_NoFatJet"] = map_Region_to_Bool[Suffix+"_HighCR"] && TwoJet_NoFatJet;

        map_Region_to_Bool[Suffix+"_High_OneFatJet"]   = map_Region_to_Bool[Suffix+"_High"] && OneFatJet;
        map_Region_to_Bool[Suffix+"_HighCR_OneFatJet"] = map_Region_to_Bool[Suffix+"_HighCR"] && OneFatJet;


        if(Suffix.Contains("EMu")){
          if(lep.at(1).LeptonFlavour()==KLepton::ELECTRON){
            map_Region_to_Bool[Suffix+"_Preselection_ElectronSubLead"] = true;
          }
          else{
            map_Region_to_Bool[Suffix+"_Preselection_MuonSubLead"] = true;
          }
        }
      }

      //==== Test 1) Loose sample faking jet pt distribution
      if(LooseSampleFakeJetPt){

        if(map_Region_to_Bool[Suffix+"_Preselection"] && jets_nolepveto.size() > 0){

          for(unsigned int j=0; j<lep.size(); j++){

            if(lep.at(j).LeptonFlavour() == KLepton::MUON){
              if( PassID( *(lep.at(j).GetMuonPtr()), MuonTightID ) ) continue;
            }
            else if(lep.at(j).LeptonFlavour() == KLepton::ELECTRON){
              if( PassID( *(lep.at(j).GetElectronPtr()), ElectronTightID ) ) continue;
            }
            else continue;

            for(unsigned int k=0; k<jets_nolepveto.size(); k++){
              double this_dr = (lep.at(j).DeltaR( jets_nolepveto.at(k) ) );
              if(this_dr < 0.3){
                FillHist(Suffix+"_LooseSample_LeptonClosestJet_Pt", jets_nolepveto.at(k).Pt(), 1, 0., 200., 200);
              }
            }


          }

        }

        continue;
      }

      //==== Test 2) Awayjet pt variation on fake yield
      if(DoAwayJet){
        if(map_Region_to_Bool[Suffix+"_Preselection"] && isSSForCF){

          for(unsigned int j=0; j<FRweights.size(); j++){
            FillHist(Suffix+"_FRsyst_"+FRweights_name.at(j), 0., FRweights.at(j), 0., 1., 1);
          }

        }
        continue;
      }

      //==== Test 3) If MCClosure, only save Preselection
      if(DoMCClosure){
        bool tempbool_presel = map_Region_to_Bool[Suffix+"_Preselection"];
        bool tempbool_alljets = true;
        bool tempbool_0jets = map_Region_to_Bool[Suffix+"_0jets"];
        bool tempbool_1jets = map_Region_to_Bool[Suffix+"_1jets"];
        map_Region_to_Bool.clear();
        map_Region_to_Bool[Suffix+"_Preselection"] = tempbool_presel;
        map_Region_to_Bool[Suffix+"_alljets"] = tempbool_alljets;
        map_Region_to_Bool[Suffix+"_0jets"] = tempbool_0jets;
        map_Region_to_Bool[Suffix+"_1jets"] = tempbool_1jets;
        bool PositiveMCWeight = std::find(k_flags.begin(), k_flags.end(), "PositiveMCWeight") != k_flags.end();
        if(PositiveMCWeight) this_weight *= MCweight;

        if( map_Region_to_Bool[Suffix+"_Preselection"] && isSSForCF ){
          FillHist("NLooseNotTight_"+Suffix+"_Preselection_SS", n_triLoose_leptons-n_triTight_leptons, this_weight, 0., 3., 3);
          FillHist("NLooseNotTight_weight1_"+Suffix+"_Preselection_SS", n_triLoose_leptons-n_triTight_leptons, 1., 0., 3., 3);
        }
      }

      //==== Make Ntuple
      bool NtupleSkim = (map_Region_to_Bool[Suffix+"_Preselection"]) && isSSForCF;

      if(RunNtp && NtupleSkim){
        double cutop[100];
        cutop[0] = lep.at(0).Pt();
        cutop[1] = lep.at(1).Pt();
        cutop[2] = lep.at(0).DeltaR( lep.at(1) );
        cutop[3] = (lep.at(0)+lep.at(1)).M();
        cutop[4] = isSSForCF;
        cutop[5] = isOffZ ? 1 : 0;

        cutop[6] = jets.size();
        cutop[7] = nbjets;
        cutop[8] = jets_nolepveto.size();
        cutop[9] = nbjets_nolepveto;
        cutop[10] = jets_fwd.size();
        cutop[11] = nbjets_fwd;

        //==== Two Jets no FatJet
        if(TwoJet_NoFatJet){

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
        //==== FatJet
        else if(OneFatJet){
          for(int j=12;j<=40;j++) cutop[j] = -999.;

          //==== pt order
          cutop[12] = fatjets.at(0).Pt();
          cutop[13] = -999.;
          cutop[14] = -999.;
          cutop[15] = (fatjets.at(0)).PrunedMass();
          cutop[16] = (lep.at(0)+fatjets.at(0)).M();
          cutop[17] = (lep.at(1)+fatjets.at(0)).M();
          cutop[18] = (lep.at(0)+lep.at(1)+fatjets.at(0)).M();

          //==== m(fj)~W (high mass)
          cutop[19] = fatjets.at(index_fjW).Pt();
          cutop[20] = -999.;
          cutop[21] = (fatjets.at(index_fjW)).PrunedMass();
          cutop[22] = (lep.at(0)+fatjets.at(index_fjW)).M();
          cutop[23] = (lep.at(1)+fatjets.at(index_fjW)).M();
          cutop[24] = (lep.at(0)+lep.at(1)+fatjets.at(index_fjW)).M();
          cutop[25] = -999.;
          cutop[26] = lep.at(0).DeltaR( fatjets.at(index_fjW) );
          cutop[27] = lep.at(1).DeltaR( fatjets.at(index_fjW) );
          cutop[28] = lep.at(0).DeltaR( lep.at(1)+fatjets.at(index_fjW) );
          cutop[29] = lep.at(1).DeltaR( lep.at(0)+fatjets.at(index_fjW) );

        }
        else if(OneJet_NoFatJet){
          for(int j=12;j<=40;j++) cutop[j] = -999.;

          cutop[30] = jets.at(0).Pt();
          //cutop[31] = jets.at(index_lljjW_j2).Pt();
          cutop[32] = (jets.at(0)).M();
          cutop[33] = (lep.at(0)+jets.at(0)).M();
          cutop[34] = (lep.at(1)+jets.at(0)).M();
          cutop[35] = (lep.at(0)+lep.at(1)+jets.at(0)).M();
          //cutop[36] = jets.at(index_lljjW_j1).DeltaR(jets.at(index_lljjW_j2));
          cutop[37] = lep.at(0).DeltaR( jets.at(0) );
          cutop[38] = lep.at(1).DeltaR( jets.at(0) );
          //cutop[39] = lep.at(0).DeltaR( lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );
          //cutop[40] = lep.at(1).DeltaR( lep.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) );

        }
        else{
          for(int j=12;j<=40;j++) cutop[j] = -999.;
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

        cutop[48] = lep.at(0).Eta();
        cutop[49] = lep.at(1).Eta();

        cutop[50] = fatjets.size();


        FillNtp("Ntp_"+Suffix+"_Preselection_SS",cutop);

      }
      if(RunNtp) continue;

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

            //==== All m(ll) / SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              FillHist("LeptonType_"+this_suffix+"_SS", lep.at(0).GetType(), 1., 0., 50., 50);
              FillHist("LeptonType_"+this_suffix+"_SS", lep.at(1).GetType(), 1., 0., 50., 50);
            }

            if(isOffZ){
              //==== OffZ / SS
              if(isSS){
                FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
              if(isAboveZ){
                //==== AboveZ / SS
                if(isSS){
                  FillDiLeptonPlot(this_suffix+"_AboveZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
              }
            }
            else{
              //==== OnZ / SS
              if(isSS){
                FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
            }

/*
            //==== All m(ll) / SS+OS
            FillDiLeptonPlot(this_suffix+"_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            //==== All m(ll) / SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              FillHist("LeptonType_"+this_suffix+"_SS", lep.at(0).GetType(), 1., 0., 50., 50);
              FillHist("LeptonType_"+this_suffix+"_SS", lep.at(1).GetType(), 1., 0., 50., 50);
            }
            //==== All m(ll) / OS
            else{
              FillDiLeptonPlot(this_suffix+"_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }

            if(isOffZ){

              //==== OffZ / SS+OS
              FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

              //==== OffZ / SS
              if(isSS){
                FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
              //==== OffZ / OS
              else{
                FillDiLeptonPlot(this_suffix+"_OffZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }

              if(isAboveZ){

                //==== AboveZ / SS+OS
                FillDiLeptonPlot(this_suffix+"_AboveZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

                //==== AboveZ / SS
                if(isSS){
                  FillDiLeptonPlot(this_suffix+"_AboveZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
                //==== AboveZ / OS
                else{
                  FillDiLeptonPlot(this_suffix+"_AboveZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }

              }

              if(isBelowZ){

                //==== BelowZ / SS+OS
                FillDiLeptonPlot(this_suffix+"_BelowZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

                //==== BelowZ / SS
                if(isSS){
                  FillDiLeptonPlot(this_suffix+"_BelowZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
                //==== BelowZ / OS
                else{
                  FillDiLeptonPlot(this_suffix+"_BelowZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }

              }

            }
            else{

              //==== OnZ / SS+OS
              FillDiLeptonPlot(this_suffix+"_OnZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);

              //==== OnZ / SS
              if(isSS){
                FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
              //==== OnZ /OS
              else{
                FillDiLeptonPlot(this_suffix+"_OnZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }

            }
*/


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
              //==== OffZ
              if(isOffZ){
                FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
              //==== OnZ
              else{
                FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
              }
              
            } // END fill chargeflip only for DiElectron OS
          }


        } // END passing this region
      } // END Search Region loop

    } // END Systematic loop

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
  
  //==== Prompt Rate
  delete hist_Muon_PR;
  delete hist_Electron_PR;
  //==== Fake Rate
  delete hist_Muon_FR;
  delete hist_Electron_FR;
  delete hist_Muon_FR_withbjet;
  delete hist_Electron_FR_withbjet;
  delete hist_Muon_FR_withoutbjet;
  delete hist_Electron_FR_withoutbjet;
  delete hist_Muon_FR_QCD;
  delete hist_Electron_FR_QCD;
  delete hist_Muon_FR_QCD_withbjet;
  delete hist_Electron_FR_QCD_withbjet;
  delete hist_Muon_FR_QCD_withoutbjet;
  delete hist_Electron_FR_QCD_withoutbjet;

  hist_Muon_FR_syst.clear();
  hist_Electron_FR_syst.clear();


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

  int N_sys = 2*11+1;
  TString systs[] = {"", "_MuonEn_up", "_MuonEn_down", "_JetEn_up", "_JetEn_down", "_JetRes_up", "_JetRes_down", "_Unclustered_up", "_Unclustered_down", "_MuonIDSF_up", "_MuonIDSF_down", "_PU_down", "_PU_up", "_TriggerSF_down", "_TriggerSF_up", "_ElectronIDSF_up", "_ElectronIDSF_down", "_ElectronEn_up", "_ElectronEn_down", "_BTagSFEff_up", "_BTagSFEff_down", "_BTagSFMiss_up", "_BTagSFMiss_down"};

  for(int i=0; i<N_sys; i++){

  MakeNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");

  MakeNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");

  MakeNtp("Ntp_EMu"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");


    GetNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);
    GetNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);
    GetNtp("Ntp_EMu"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_EMu"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);

  }


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

  if(!NowRunningCentral) return;

  TString leporder[4] = {"leading", "second", "third", "fourth"};

  JSFillHist(histsuffix, "Nevents_"+histsuffix, 0., thisweight, 0., 1., 1);
  JSFillHist(histsuffix, "weight_fr_"+histsuffix, weight_fr, 1., -2., 2., 400);
  JSFillHist(histsuffix, "PFMET_"+histsuffix, MET, thisweight, 0., 1000., 1000);
  JSFillHist(histsuffix, "PFMET_phi_"+histsuffix, METphi, thisweight, -3.2, 3.2, 100);
  JSFillHist(histsuffix, "Njets_"+histsuffix, jets.size(), thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Njets_nolepveto_"+histsuffix, jets_nolepveto.size(), thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nfwdjets_"+histsuffix, jets_fwd.size(), thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nbjets_"+histsuffix, nbjets, thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nbjets_nolepveto_"+histsuffix, nbjets_nolepveto, thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nbfwdjets_"+histsuffix, nbjets_fwd, thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nfatjets_"+histsuffix, fatjets.size(), thisweight, 0., 10., 10);
  JSFillHist(histsuffix, "Nvtx_"+histsuffix, n_vtx, thisweight, 0., 100., 100);

  JSFillHist(histsuffix, "HT_"+histsuffix, HT, thisweight, 0., 2000., 2000);
  JSFillHist(histsuffix, "ST_"+histsuffix, ST, thisweight, 0., 2000., 2000);
  JSFillHist(histsuffix, "LT_"+histsuffix, LT, thisweight, 0., 2000., 2000);
  JSFillHist(histsuffix, "MCT_"+histsuffix, contramass, thisweight, 0., 2000., 2000);
  JSFillHist(histsuffix, "MET2overST_"+histsuffix, MET*MET/ST, thisweight, 0., 2000., 2000);
  JSFillHist(histsuffix, "DeltaRl1l2_"+histsuffix, leptons.at(0).DeltaR(leptons.at(1)), thisweight, 0., 10., 100);

  JSFillHist(histsuffix, "m_ll_"+histsuffix, (leptons.at(0)+leptons.at(1)).M(), thisweight, 0., 2000., 2000);

  JSFillHist(histsuffix, "NTightLeptons_weighted_"+histsuffix,   NTightLeptons, thisweight, 0., 5., 5);
  JSFillHist(histsuffix, "NTightLeptons_unweighted_"+histsuffix, NTightLeptons, 1., 0., 5., 5);
  for(unsigned int i=0; i<leptons.size(); i++){
    if(i==4) break;
    JSFillHist(histsuffix, leporder[i]+"Lepton_Pt_"+histsuffix,  leptons.at(i).Pt(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"Lepton_Eta_"+histsuffix, leptons.at(i).Eta(), thisweight, -3., 3., 60);
    JSFillHist(histsuffix, leporder[i]+"Lepton_Type_"+histsuffix, leptons.at(i).GetType(), thisweight, 0., 50., 50);
    JSFillHist(histsuffix, leporder[i]+"Lepton_RelIso_"+histsuffix, leptons.at(i).RelIso(), thisweight, 0., 1., 100);
    JSFillHist(histsuffix, leporder[i]+"Lepton_dXY_"+histsuffix, fabs(leptons.at(i).dXY()), thisweight, 0., 1., 1000);
    JSFillHist(histsuffix, leporder[i]+"Lepton_dXYSig_"+histsuffix, fabs(leptons.at(i).dXYSig()), thisweight, 0., 40., 400);

    double TightIso = 0.07;
    if(leptons.at(i).LeptonFlavour()==KLepton::ELECTRON){
      TightIso = 0.08;
      JSFillHist(histsuffix, leporder[i]+"Lepton_mva_"+histsuffix, leptons.at(i).GetElectronPtr()->MVA(), thisweight, -1., 1., 200);
    }
    else{
      JSFillHist(histsuffix, leporder[i]+"Lepton_Chi2_"+histsuffix, leptons.at(i).GetMuonPtr()->GlobalChi2(), thisweight, 0., 200., 200);
    }
    JSFillHist(histsuffix, leporder[i]+"Lepton_Pt_cone_"+histsuffix, CorrPt(leptons.at(i), TightIso), thisweight, 0., 2000., 2000);

  }
  for(unsigned int i=0; i<jets.size(); i++){
    if(i==4) break;
    JSFillHist(histsuffix, leporder[i]+"Jet_Pt_"+histsuffix,  jets.at(i).Pt(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"Jet_Eta_"+histsuffix, jets.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(unsigned int i=0; i<jets_fwd.size(); i++){
    if(i==4) break;
    JSFillHist(histsuffix, leporder[i]+"ForwardJet_Pt_"+histsuffix,  jets_fwd.at(i).Pt(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"ForwardJet_Eta_"+histsuffix, jets_fwd.at(i).Eta(), thisweight, -5., 5., 100);
  }
  for(unsigned int i=0; i<jets_nolepveto.size(); i++){
    if(i==4) break;
    JSFillHist(histsuffix, leporder[i]+"NoLepVetoJet_Pt_"+histsuffix, jets_nolepveto.at(i).Pt(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"NoLepVetoJet_Eta_"+histsuffix, jets_nolepveto.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(unsigned int i=0; i<fatjets.size(); i++){
    if(i==index_fjW){
      JSFillHist(histsuffix, "WClosest_FatJet_Pt_"+histsuffix, fatjets.at(i).Pt(), thisweight, 0., 2000., 2000);
      JSFillHist(histsuffix, "WClosest_FatJet_Eta_"+histsuffix, fatjets.at(i).Eta(), thisweight, -3., 3., 60);
      JSFillHist(histsuffix, "WClosest_FatJet_Mass_"+histsuffix, fatjets.at(i).M(), thisweight, 0., 2000., 2000);
      JSFillHist(histsuffix, "WClosest_FatJet_PrunedMass_"+histsuffix, fatjets.at(i).PrunedMass(), thisweight, 0., 2000., 2000);
      JSFillHist(histsuffix, "WClosest_FatJet_Tau21_"+histsuffix, fatjets.at(i).Tau2()/fatjets.at(i).Tau1(), thisweight, 0., 1., 100);
    }
    if(i==4) break;
    JSFillHist(histsuffix, leporder[i]+"FatJet_Pt_"+histsuffix, fatjets.at(i).Pt(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"FatJet_Eta_"+histsuffix, fatjets.at(i).Eta(), thisweight, -3., 3., 60);
    JSFillHist(histsuffix, leporder[i]+"FatJet_Mass_"+histsuffix, fatjets.at(i).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"FatJet_PrunedMass_"+histsuffix, fatjets.at(i).PrunedMass(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, leporder[i]+"FatJet_Tau21_"+histsuffix, fatjets.at(i).Tau2()/fatjets.at(i).Tau1(), thisweight, 0., 1., 100);
  }

  if(jets.size() == 1){
    JSFillHist(histsuffix, "m_Leadlj_"+histsuffix, (leptons.at(0)+jets.at(0)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadlj_"+histsuffix, (leptons.at(1)+jets.at(0)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_llj_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(0)).M(), thisweight, 0., 2000., 2000);
  }
  if(jets.size() >= 2){
    //==== m(jj) closeset to m(W) : high mass scenario
    JSFillHist(histsuffix, "m_jj_jjWclosest_"+histsuffix, (jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_Leadljj_jjWclosest_"+histsuffix, (leptons.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadljj_jjWclosest_"+histsuffix, (leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_lljj_jjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "DeltaRjjWclosest_"+histsuffix, jets.at(index_jjW_j1).DeltaR(jets.at(index_jjW_j2)), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRLeadl_jjWclosest_"+histsuffix, leptons.at(0).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRSubLeadl_jjWclosest_"+histsuffix, leptons.at(1).DeltaR( jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRLeadl_SubLeadljjWclosest_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRSubLeadl_LeadljjWclosest_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+jets.at(index_jjW_j1)+jets.at(index_jjW_j2) ), thisweight, 0., 10., 100);

    //==== m(lljj) cloeset to m(W) : low mass scenario
    JSFillHist(histsuffix, "m_jj_lljjWclosest_"+histsuffix, (jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_Leadljj_lljjWclosest_"+histsuffix, (leptons.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadljj_lljjWclosest_"+histsuffix, (leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_lljj_lljjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "DeltaRlljjWclosest_"+histsuffix, jets.at(index_lljjW_j1).DeltaR(jets.at(index_lljjW_j2)), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRLeadl_lljjWclosest_"+histsuffix, leptons.at(0).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRSubLeadl_lljjWclosest_"+histsuffix, leptons.at(1).DeltaR( jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRLeadl_SubLeadllljjWclosest_"+histsuffix, leptons.at(0).DeltaR( leptons.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "DeltaRSubLeadl_LeadllljjWclosest_"+histsuffix, leptons.at(1).DeltaR( leptons.at(0)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2) ), thisweight, 0., 10., 100);

    //==== ptorder
    JSFillHist(histsuffix, "DeltaRjjptorder_"+histsuffix, jets.at(0).DeltaR( jets.at(1) ), thisweight, 0., 10., 100);
    JSFillHist(histsuffix, "m_jjptorder_"+histsuffix, (jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_Leadljjptorder_"+histsuffix, (leptons.at(0)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadljjptorder_"+histsuffix, (leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_lljjptorder_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
  }

  if(fatjets.size() > 0){
    //==== Leading
    JSFillHist(histsuffix, "m_Leadlfj_ptorder_"+histsuffix, (leptons.at(0)+fatjets.at(0)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadlfj_ptorder_"+histsuffix, (leptons.at(1)+fatjets.at(0)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_llfj_ptorder_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjets.at(0)).M(), thisweight, 0., 2000., 2000);

    //==== fjWclosest
    JSFillHist(histsuffix, "m_Leadlfj_fjWclosest_"+histsuffix, (leptons.at(0)+fatjets.at(index_fjW)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_SubLeadlfj_fjWclosest_"+histsuffix, (leptons.at(1)+fatjets.at(index_fjW)).M(), thisweight, 0., 2000., 2000);
    JSFillHist(histsuffix, "m_llfj_fjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjets.at(index_fjW)).M(), thisweight, 0., 2000., 2000);
  }

  if(thieweighterr!=0.){
    FillDiLeptonPlot(histsuffix+"_up",   leptons, jets, jets_fwd, jets_nolepveto, fatjets, thisweight + thieweighterr, 0.);
    FillDiLeptonPlot(histsuffix+"_down", leptons, jets, jets_fwd, jets_nolepveto, fatjets, thisweight - thieweighterr, 0.);

    //==== Check Single/Double Fake

    //==== 1) LL : thisweight = -e^2
    if(NTightLeptons==0){
      FillDiLeptonPlot(histsuffix+"_SingleFake", leptons, jets, jets_fwd, jets_nolepveto, fatjets, 2.*thisweight, 0.);
      FillDiLeptonPlot(histsuffix+"_DoubleFake", leptons, jets, jets_fwd, jets_nolepveto, fatjets, -1.*thisweight, 0.);
    }
    //==== 2) TL : thisweight = e
    if(NTightLeptons==1){
      FillDiLeptonPlot(histsuffix+"_SingleFake", leptons, jets, jets_fwd, jets_nolepveto, fatjets, thisweight, 0.);
    }

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
    if(invPt< 0.023){a=(-0.001381); b=(4.334e-05);}
    else{a=(0.001010); b=(-1.146e-05);}
  }
  else if(el_eta < 1.4442){
    if(invPt< 0.015){a=(-0.04296); b=(0.0008670);}
    else if(invPt< 0.023){a=(-0.01529); b=(0.0004522);}
    else{a=(-0.001546); b=(0.0001272);}
  }
  else{
    if(invPt< 0.012){a=(-0.4238); b=(0.006366);}
    else if(invPt< 0.020){a=(-0.1040); b=(0.002550);}
    else{a=(-0.01603); b=(0.0007672);}
  }

  double sf(1.);
  if(el_eta < 1.4442) sf = 0.693589;
  else sf = 0.684761;

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

double DiLeptonAnalyzer::GetFatJetMassClosest(std::vector<snu::KFatJet> fjs, double mass, int& m){

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



double DiLeptonAnalyzer::GetMuonFR(bool geterr, float pt,  float eta, int NearBjet){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.4) eta = 2.3;

  //cout << "[DiLeptonAnalyzer::GetMuonFR] pt = " << pt << endl;
  //cout << "[DiLeptonAnalyzer::GetMuonFR] eta = " << eta << endl;
  //cout << "[DiLeptonAnalyzer::GetMuonFR] MuFR_key = " << MuFR_key << endl;
  //cout << "[DiLeptonAnalyzer::GetMuonFR] NearBjet = " << NearBjet << endl;

  TH2D *THISFRHIST = NULL;

  if(MuFR_key==""){
    if(NearBjet==1)      THISFRHIST = hist_Muon_FR_withbjet;
    else if(NearBjet==0) THISFRHIST = hist_Muon_FR_withoutbjet;
    else                 THISFRHIST = hist_Muon_FR;
  }
  else if(MuFR_key=="QCD"){
    if(NearBjet==1)      THISFRHIST = hist_Muon_FR_QCD_withbjet;
    else if(NearBjet==0) THISFRHIST = hist_Muon_FR_QCD_withoutbjet;
    else                 THISFRHIST = hist_Muon_FR_QCD;
  }
  else{
    THISFRHIST = hist_Muon_FR_syst[MuFR_key];
  }

  int binx = THISFRHIST->FindBin(pt, abs(eta));
  //cout << "[DiLeptonAnalyzer::GetMuonFR] => FR = " << THISFRHIST->GetBinContent(binx) << endl;
  

  if(geterr) return THISFRHIST->GetBinError(binx);
  else return THISFRHIST->GetBinContent(binx);

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

double DiLeptonAnalyzer::GetElectronFR(bool geterr, float pt,  float eta, int NearBjet){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.5) eta = 2.4;

  TH2D *THISFRHIST = NULL;

  if(ElFR_key==""){
    if(NearBjet==1)      THISFRHIST = hist_Electron_FR_withbjet;
    else if(NearBjet==0) THISFRHIST = hist_Electron_FR_withoutbjet;
    else                 THISFRHIST = hist_Electron_FR;
  }
  else if(ElFR_key=="QCD"){
    if(NearBjet==1)      THISFRHIST = hist_Electron_FR_QCD_withbjet;
    else if(NearBjet==0) THISFRHIST = hist_Electron_FR_QCD_withoutbjet;
    else                 THISFRHIST = hist_Electron_FR_QCD;
  }
  else{
    THISFRHIST = hist_Electron_FR_syst[ElFR_key];
  }

  int binx = THISFRHIST->FindBin(pt, abs(eta));
  if(geterr) return THISFRHIST->GetBinError(binx);
  else return THISFRHIST->GetBinContent(binx);

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


void DiLeptonAnalyzer::get_eventweight(std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<bool> isT, std::vector<int> NearBjet, int HalfSampleErrorDir){

  unsigned int n_leptons = isT.size();
  //cout << "[DiLeptonAnalyzer::get_eventweight] muons.size() = " << muons.size() << ", electrons.size() = " << electrons.size() << endl;

  vector<float> lep_pt, lep_eta;
  vector<bool> ismuon;
  for(unsigned int i=0; i<muons.size(); i++){
    lep_pt.push_back( CorrPt(muons.at(i), 0.07) );
    lep_eta.push_back(muons.at(i).Eta());
    ismuon.push_back(true);
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    lep_pt.push_back( CorrPt(electrons.at(i), 0.08) );
    lep_eta.push_back(electrons.at(i).Eta());
    ismuon.push_back(false);

  }

  vector<float> fr, pr, fr_err, pr_err;

  for(unsigned int i=0; i<n_leptons; i++){
    //==== Muon
    if(ismuon.at(i)){
      fr.push_back( GetMuonFR(false, lep_pt.at(i), lep_eta.at(i), NearBjet.at(i)) );
      pr.push_back( GetMuonPR(false, lep_pt.at(i), lep_eta.at(i)) );
      fr_err.push_back( GetMuonFR(true, lep_pt.at(i), lep_eta.at(i), NearBjet.at(i)) );
      pr_err.push_back( GetMuonPR(true, lep_pt.at(i), lep_eta.at(i)) );
    }
    //==== If not, it's an electron
    else{
      fr.push_back( GetElectronFR(0, lep_pt.at(i), lep_eta.at(i), NearBjet.at(i)) );
      pr.push_back( GetElectronPR(0, lep_pt.at(i), lep_eta.at(i)));
      fr_err.push_back( GetElectronFR(1, lep_pt.at(i), lep_eta.at(i), NearBjet.at(i)) );
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
      //cout << "[DiLeptonAnalyzer::get_eventweight] "<<i<<" th lepton is Loose" << endl;
      fr_onlyLoose.push_back( a.at(i) );
    }
  }

  //==== Initialise weight
  float this_weight=-1.;

  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }
  //cout << "[DiLeptonAnalyzer::get_eventweight] this_weight = " << this_weight << endl;

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

bool DiLeptonAnalyzer::JSFatJetID(snu::KFatJet fatjet){

  if( !(fatjet.PrunedMass() > 40) ) return false;
  if( !(fatjet.PrunedMass() < 130) ) return false;
  if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;
  if( !( fatjet.Pt() > 100. ) ) return false;

  return true;

}

bool DiLeptonAnalyzer::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 0.8 ) return false;
  }

  return true;


}



















