// $Id: PairNAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQPairNAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "PairNAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"
#include "TSystem.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (PairNAnalyzer);

PairNAnalyzer::PairNAnalyzer() :
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
RunNtp(false),NowRunningCentral(true),
AUTO_N_syst(0), AUTO_it_syst(0), AUTO_syst_type("")
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("PairNAnalyzer");
  
  Message("In PairNAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void PairNAnalyzer::InitialiseAnalysis() throw( LQError ) {
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

  TString ElectronFRType_QCD = "v7_5_";

  //TString ElectronFRType_Data = "v7_4_";
  TString ElectronFRType_Data = ElectronFRType_QCD;

  ElectronLooseID_loosest = "ELECTRON_HN_FAKELOOSEv7_5_loosest";
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
  TFile *file_Muon_FR_suoh    = new TFile( lqdir+"/JskimData/FR/Muon_Data_fake_Rate_syst.root");
  TFile *file_Muon_FR_QCD     = new TFile( lqdir+"/JskimData/FR/Muon_QCD_" +MuonFRType_QCD+"FR.root"); 

  TFile *file_Electron_FR     = new TFile( lqdir+"/JskimData/FR/Electron_Data_"+ElectronFRType_Data+"FR.root");
  TFile *file_Electron_FR_suoh= new TFile( lqdir+"/JskimData/FR/Electron_Data_fake_Rate_syst.root");
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

  TString btag_wp = "_Medium";
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
    cout << awayjetpt[i] << endl;
    hist_Muon_FR_syst["Awayjet_"+awayjetpt[i]]     = (TH2D*)file_Muon_FR_suoh->Get("FR_ptcone_awayjet_pt_"+awayjetpt[i])->Clone();
    hist_Electron_FR_syst["Awayjet_"+awayjetpt[i]] = (TH2D*)file_Electron_FR_suoh->Get("FR_ptcone_awayjet_pt_"+awayjetpt[i])->Clone();
  }
  cout << "40" << endl;
  hist_Muon_FR_syst["Awayjet_40"]     = (TH2D*)file_Muon_FR_suoh->Get("FR_ptcone_central")->Clone();
  hist_Electron_FR_syst["Awayjet_40"] = (TH2D*)file_Electron_FR_suoh->Get("FR_ptcone_central_EC50syst")->Clone(); //Updated in 18/02/12, EC prompt syst increased to 50%


  TString Muon_FRsystsources[5] = {"dphi_1p5", "dphi_1", "dphi_3", "pj_over_pl_1", "pj_over_pl_2"};
  for(int i=0; i<5; i++){
    cout << Muon_FRsystsources[i] << endl;
    hist_Muon_FR_syst[Muon_FRsystsources[i]] = (TH2D*)file_Muon_FR_suoh->Get("FR_ptcone_"+Muon_FRsystsources[i])->Clone();
  }
  TString Electron_FRsystsources[5] = {"dphi_1p5", "dphi_2", "dphi_3", "pj_over_pl_0p8", "pj_over_pl_1p2"};
  for(int i=0; i<5; i++){
    cout << Electron_FRsystsources[i] << endl;
    hist_Electron_FR_syst[Electron_FRsystsources[i]] = (TH2D*)file_Electron_FR_suoh->Get("FR_ptcone_"+Electron_FRsystsources[i])->Clone();
  }


  file_Muon_FR->Close();
  file_Muon_FR_QCD->Close();
  file_Muon_FR_suoh->Close();
  file_Electron_FR->Close();
  file_Electron_FR_QCD->Close();
  file_Electron_FR_suoh->Close();

  delete file_Muon_FR;
  delete file_Muon_FR_QCD;
  delete file_Muon_FR_suoh;
  delete file_Electron_FR;
  delete file_Electron_FR_QCD;
  delete file_Electron_FR_suoh;

  origDir->cd();

  ForTree_PdfWeights = new vector<float>();
  ForTree_ScaleWeights = new vector<float>();

  //cout << "Muon, 17,0.4 = " << hist_Muon_FR->GetBinContent(hist_Muon_FR->FindBin(17,0.4)) << endl;
  //cout << "Electron, 17,0.4 = " << hist_Electron_FR->GetBinContent(hist_Electron_FR->FindBin(17,0.4)) << endl;

  return;
}


void PairNAnalyzer::ExecuteEvents()throw( LQError ){

  std::vector< snu::KMuon > testmuons_all = GetMuons("MUON_NOCUT", true, 10., 2.4);
  std::vector< snu::KMuon > testmuons;
  std::vector<snu::KMuon> SUSYmuons;

  std::vector< snu::KElectron > testelectrons_all = GetElectrons(true, true, "ELECTRON_NOCUT", 10, 2.5);
  std::vector< snu::KElectron > testelectrons;

  bool isHNPairSignal = k_sample_name.Contains("HNpair");
  for(unsigned int i=0; i<testmuons_all.size(); i++){
    if(isHNPairSignal)  testmuons.push_back( testmuons_all.at(i) );
    else{
      if( !TruthMatched( testmuons_all.at(i) ) ) testmuons.push_back( testmuons_all.at(i) );
    }
  }
  for(unsigned int i=0; i<testelectrons_all.size(); i++){
    if(isHNPairSignal)  testelectrons.push_back( testelectrons_all.at(i) );
    else{
      if( !TruthMatched( testelectrons_all.at(i),true ) ) testelectrons.push_back( testelectrons_all.at(i) );
    }
  }

  std::vector<snu::KJet> test_jets = GetJets("JET_HN_eta5_nolepveto");
  std::vector<snu::KFatJet> test_fatjets_chs = GetFatJets("FATJET_HN_nolepveto");
  std::vector<snu::KFatJet> test_fatjets; // Puppi + SD

  cout << "####### Making Puppi #######" << endl;
  cout << "Pt\tEta\tPhi\tM\tSD" << endl;
  for(unsigned int i=0; i < test_fatjets_chs.size(); i++){
    snu::KFatJet tmpj = test_fatjets_chs[i];

    cout << i << "th fatjet" << endl; 
    cout << "  " << tmpj.Pt() << "\t" << tmpj.Eta() << "\t" << tmpj.Phi() << "\t" << tmpj.M() << "\t" << tmpj.SoftDropMass() << endl;

	  tmpj.SetPtEtaPhiE(tmpj.PuppiPt(), tmpj.PuppiEta(),tmpj.PuppiPhi(),tmpj.E());
    cout << "->" << tmpj.Pt() << "\t" << tmpj.Eta() << "\t" << tmpj.Phi() << "\t" << tmpj.M() << "\t" << tmpj.SoftDropMass() << endl;

    test_fatjets.push_back(tmpj);
  }
  std::sort(test_fatjets.begin(), test_fatjets.end(), FatjetSoftDropMassComparing);

  vector<TString> test_trig_dimuon;
  test_trig_dimuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  test_trig_dimuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  test_trig_dimuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  test_trig_dimuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  vector<TString> test_trig_IsoMu24;
  test_trig_IsoMu24.push_back("HLT_IsoMu24_v");
  test_trig_IsoMu24.push_back("HLT_IsoTkMu24_v");

  vector<TString> test_trig_Mu50;
  test_trig_IsoMu24.push_back("HLT_Mu50_v");

  vector< vector<TString> > test_triggerlists;
  test_triggerlists.push_back(test_trig_dimuon);
  test_triggerlists.push_back(test_trig_IsoMu24);
  test_triggerlists.push_back(test_trig_Mu50);

  vector< TString > test_triggerlists_histnames;
  test_triggerlists_histnames.push_back("Mu17Mu8");
  test_triggerlists_histnames.push_back("IsoMu24");
  test_triggerlists_histnames.push_back("Mu50");

  FillHist("Muon_Size", testmuons.size(), 1., 0., 10., 10);
  for(unsigned int i=0; i<testmuons.size(); i++){

    snu::KMuon muon = testmuons.at(i);
    FillHist("Muon_"+TString::Itoa(i,10)+"_Pt", muon.Pt(), 1., 0., 5000., 5000);
    FillHist("Muon_"+TString::Itoa(i,10)+"_Eta", muon.Eta(), 1., -3., 3., 60);
    FillHist("Muon_"+TString::Itoa(i,10)+"_dXY", fabs(muon.dXY()), 1., 0., 0.5, 500);
    FillHist("Muon_"+TString::Itoa(i,10)+"_dXYSig", fabs(muon.dXYSig()), 1., 0., 10, 100);
    FillHist("Muon_"+TString::Itoa(i,10)+"_dZ", fabs(muon.dZ()), 1., 0., 0.5, 500);
    FillHist("Muon_"+TString::Itoa(i,10)+"_RelIso03", muon.RelIso03(), 1., 0., 1.0, 100);
    FillHist("Muon_"+TString::Itoa(i,10)+"_RelIso04", muon.RelIso04(), 1., 0., 1.0, 100);
    FillHist("Muon_"+TString::Itoa(i,10)+"_RelMiniIso", muon.RelMiniIso(), 1., 0., 1.0, 100);

    //==== j    0      1      2               18       19
    //==== pt 10-20, 20-30, 30-40, ... ,    190-200, 200-
    for(int j=0; j<20; j++){

      double this_pt = muon.Pt();
      double ptmin = 10+j*10;
      double ptmax = ptmin+10;
      if(j==19) ptmax=9999;

      if( ptmin < this_pt && this_pt < ptmax ){
        FillHist("Muon_"+TString::Itoa(i,10)+"_RelIso04_PtBin_"+TString::Itoa(j,10), muon.RelIso04(), 1., 0., 1.0, 100);
        FillHist("Muon_"+TString::Itoa(i,10)+"_RelMiniIso_PtBin_"+TString::Itoa(j,10), muon.RelMiniIso(), 1., 0., 1.0, 100);
      }

    }

    bool isFound = false;
    double dR = 0.4;
    double ptratio = 1.;
    double ptrel = 0.;
    for(unsigned int j=0; j<test_jets.size(); j++){
      double this_dR = muon.DeltaR( test_jets.at(j) );
      if(this_dR < dR){
        isFound = true;
        dR = this_dR;

        snu::KJet corjet = GetCorrectedJetCloseToLepton(muon, test_jets.at(j));

        ptratio = muon.Pt()/corjet.Pt();

        TVector3 el3 =  muon.Vect();
        TVector3 jet3 = corjet.Vect();
        TVector3 lepjetrel = jet3-el3;
        ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();

      }
    }
    if(isFound){
      FillHist("Muon_"+TString::Itoa(i,10)+"_PtRatio", ptratio, 1., 0., 3.0, 300);
      FillHist("Muon_"+TString::Itoa(i,10)+"_PtRel", ptrel, 1., 0., 100., 100);
      FillHist("Muon_"+TString::Itoa(i,10)+"_PtRel_vs_PtRatio", ptrel, ptratio, 1., 0., 30., 30, 0., 1.5, 150);

      FillHist("Muon_PtRatio", ptratio, 1., 0., 3.0, 300);
      FillHist("Muon_PtRel", ptrel, 1., 0., 100., 100);
      FillHist("Muon_PtRel_vs_PtRatio", ptrel, ptratio, 1., 0., 30., 30, 0., 1.5, 150);

    }
    double reliso04 = muon.RelIso04();
    double minireliso = muon.RelMiniIso();
    double WhichIsoPass = 0;

    //==== None
      FillHist("Muon_IsolationPass", 0, 1., 0., 20., 20);
    //==== RelIso
    if( reliso04 < 0.15 ){
      FillHist("Muon_IsolationPass", 1, 1., 0., 20., 20);
    }
    //==== MiniIso
    if( minireliso < 0.15 ){
      FillHist("Muon_IsolationPass", 2, 1., 0., 20., 20);
    }
    //==== MultiIso loose, Only MiniIso
    if( minireliso < 0.4 ){
      FillHist("Muon_IsolationPass", 3, 1., 0., 20., 20);
    }
    if( PassMultiIso("Muon_loose", minireliso, ptratio, ptrel) ){
      FillHist("Muon_IsolationPass", 4, 1., 0., 20., 20);
    }
    //==== MultiIso tight, Only MiniIso
    if( minireliso < 0.16 ){
      FillHist("Muon_IsolationPass", 5, 1., 0., 20., 20);
    }
    if( PassMultiIso("Muon_tight", minireliso, ptratio, ptrel) ){
      FillHist("Muon_IsolationPass", 6, 1., 0., 20., 20);

      if(PassID(muon, "MUON_SUSY_TIGHT")){
        SUSYmuons.push_back( muon );
      }

    }

    FillHist("Muon_"+TString::Itoa(i,10)+"_RelIso04_vs_Pt", muon.RelIso04(), muon.Pt(),     1., 0., 1.0, 100, 0, 500, 50);
    FillHist("Muon_"+TString::Itoa(i,10)+"_RelMiniIso_vs_Pt", muon.RelMiniIso(), muon.Pt(), 1., 0., 1.0, 100, 0, 500, 50);

    FillHist("Muon_"+TString::Itoa(i,10)+"_IsPF", muon.IsPF(), 1., 0., 2., 2);
    FillHist("Muon_"+TString::Itoa(i,10)+"_IsGlobal", muon.IsGlobal(), 1., 0., 2., 2);
    FillHist("Muon_"+TString::Itoa(i,10)+"_IsTracker", muon.IsTracker(), 1., 0., 2., 2);
    FillHist("Muon_"+TString::Itoa(i,10)+"_validHits", muon.validHits(), 1., 0., 100., 100);
    FillHist("Muon_"+TString::Itoa(i,10)+"_validPixHits", muon.validPixHits(), 1., 0., 20., 20);
    FillHist("Muon_"+TString::Itoa(i,10)+"_validStations", muon.validStations(), 1., 0., 20., 20);
    FillHist("Muon_"+TString::Itoa(i,10)+"_ActiveLayer", muon.ActiveLayer(), 1., 0., 30., 30);
    FillHist("Muon_"+TString::Itoa(i,10)+"_GlobalChi2", muon.GlobalChi2(), 1., 0., 30., 30);


    FillHist("Muon_Pt", muon.Pt(), 1., 0., 5000., 5000);
    FillHist("Muon_Eta", muon.Eta(), 1., -3., 3., 60);
    FillHist("Muon_dXY", fabs(muon.dXY()), 1., 0., 0.5, 500);
    FillHist("Muon_dXYSig", fabs(muon.dXYSig()), 1., 0., 10, 100);
    FillHist("Muon_dZ", fabs(muon.dZ()), 1., 0., 0.5, 500);
    FillHist("Muon_RelIso03", muon.RelIso03(), 1., 0., 1.0, 100);
    FillHist("Muon_RelIso04", muon.RelIso04(), 1., 0., 1.0, 100);
    FillHist("Muon_RelMiniIso", muon.RelMiniIso(), 1., 0., 1.0, 100);
    FillHist("Muon_IsPF", muon.IsPF(), 1., 0., 2., 2);
    FillHist("Muon_IsGlobal", muon.IsGlobal(), 1., 0., 2., 2);
    FillHist("Muon_IsTracker", muon.IsTracker(), 1., 0., 2., 2);
    FillHist("Muon_validHits", muon.validHits(), 1., 0., 100., 100);
    FillHist("Muon_validPixHits", muon.validPixHits(), 1., 0., 20., 20);
    FillHist("Muon_validStations", muon.validStations(), 1., 0., 20., 20);
    FillHist("Muon_ActiveLayer", muon.ActiveLayer(), 1., 0., 30., 30);
    FillHist("Muon_GlobalChi2", muon.GlobalChi2(), 1., 0., 30., 30);

  }

  if(testmuons.size()>=2){
    for(unsigned int i=0; i<test_triggerlists.size(); i++){

      vector<TString> this_triggerlist = test_triggerlists.at(i);
      TString this_triggerlist_name = test_triggerlists_histnames.at(i);

      FillHist("TriggerPass_"+this_triggerlist_name, PassTriggerOR(this_triggerlist), 1., 0., 2., 2);

    }
  }

  FillHist("Electron_Size", testelectrons.size(), 1., 0., 10., 10);
  for(unsigned int i=0; i<testelectrons.size(); i++){

    snu::KElectron electron = testelectrons.at(i);
    FillHist("Electron_"+TString::Itoa(i,10)+"_Pt", electron.Pt(), 1., 0., 5000., 5000);
    FillHist("Electron_"+TString::Itoa(i,10)+"_Eta", electron.Eta(), 1., -3., 3., 60);
    FillHist("Electron_"+TString::Itoa(i,10)+"_dXY", fabs(electron.dXY()), 1., 0., 0.5, 500);
    FillHist("Electron_"+TString::Itoa(i,10)+"_dXYSig", fabs(electron.dXYSig()), 1., 0., 10, 100);
    FillHist("Electron_"+TString::Itoa(i,10)+"_dZ", fabs(electron.dZ()), 1., 0., 0.5, 500);
    FillHist("Electron_"+TString::Itoa(i,10)+"_RelIso03", electron.PFRelIso(0.3), 1., 0., 1.0, 100);
    FillHist("Electron_"+TString::Itoa(i,10)+"_RelMiniIso", electron.PFRelMiniIso(), 1., 0., 1.0, 100);

    //==== j    0      1      2               18       19
    //==== pt 10-20, 20-30, 30-40, ... ,    190-200, 200-
    for(int j=0; j<20; j++){

      double this_pt = electron.Pt();
      double ptmin = 10+j*10;
      double ptmax = ptmin+10;
      if(j==19) ptmax=9999;

      if( ptmin < this_pt && this_pt < ptmax ){
        FillHist("Electron_"+TString::Itoa(i,10)+"_RelIso03_PtBin_"+TString::Itoa(j,10), electron.PFRelIso(0.3), 1., 0., 1.0, 100);
        FillHist("Electron_"+TString::Itoa(i,10)+"_RelMiniIso_PtBin_"+TString::Itoa(j,10), electron.PFRelMiniIso(), 1., 0., 1.0, 100);
      }

    }

    bool isFound = false;
    double dR = 0.4;
    double ptratio = 1.;
    double ptrel = 0.;
    for(unsigned int j=0; j<test_jets.size(); j++){
      double this_dR = electron.DeltaR( test_jets.at(j) );
      if(this_dR < dR){
        isFound = true;
        dR = this_dR;

        snu::KJet corjet = GetCorrectedJetCloseToLepton(electron, test_jets.at(j));

        ptratio = electron.Pt()/corjet.Pt();

        TVector3 el3 =  electron.Vect();
        TVector3 jet3 = corjet.Vect();
        TVector3 lepjetrel = jet3-el3;
        ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();

      }
    }
    if(isFound){
      FillHist("Electron_"+TString::Itoa(i,10)+"_PtRatio", ptratio, 1., 0., 3.0, 300);
      FillHist("Electron_"+TString::Itoa(i,10)+"_PtRel", ptrel, 1., 0., 100., 100);

      FillHist("Electron_PtRatio", ptratio, 1., 0., 3.0, 300);
      FillHist("Electron_PtRel", ptrel, 1., 0., 100., 100);
    }
    double RelIso03 = electron.PFRelIso(0.3);
    double minireliso = electron.PFRelMiniIso();
    double WhichIsoPass = 0;

    //==== None
      FillHist("Electron_IsolationPass", 0, 1., 0., 20., 20);
    //==== RelIso
    if( RelIso03 < 0.15 ){
      FillHist("Electron_IsolationPass", 1, 1., 0., 20., 20);
    }
    //==== MiniIso
    if( minireliso < 0.15 ){
      FillHist("Electron_IsolationPass", 2, 1., 0., 20., 20);
    }
    //==== MultiIso loose, Only MiniIso
    if( minireliso < 0.4 ){
      FillHist("Electron_IsolationPass", 3, 1., 0., 20., 20);
    }
    if( PassMultiIso("Electron_loose", minireliso, ptratio, ptrel) ){
      FillHist("Electron_IsolationPass", 4, 1., 0., 20., 20);
    }
    //==== MultiIso tight, Only MiniIso
    if( minireliso < 0.16 ){
      FillHist("Electron_IsolationPass", 5, 1., 0., 20., 20);
    }
    if( PassMultiIso("Electron_tight", minireliso, ptratio, ptrel) ){
      FillHist("Electron_IsolationPass", 6, 1., 0., 20., 20);
    }

    TString EtaRegion = "InnerBarrel";
    if(fabs(electron.SCEta()) > 1.479) EtaRegion = "EndCap";
    else if(fabs(electron.SCEta()) > 0.8) EtaRegion = "OuterBarrel";
    else EtaRegion = "InnerBarrel";


    FillHist("Electron_"+TString::Itoa(i,10)+"_RelIso03_vs_Pt", electron.PFRelIso(0.3), electron.Pt(),     1., 0., 1.0, 100, 0, 500, 50);
    FillHist("Electron_"+TString::Itoa(i,10)+"_RelMiniIso_vs_Pt", electron.PFRelMiniIso(), electron.Pt(), 1., 0., 1.0, 100, 0, 500, 50);

    FillHist("Electron_"+TString::Itoa(i,10)+"_GsfCtfScPixChargeConsistency", electron.GsfCtfScPixChargeConsistency(), 1., 0., 2., 2);
    FillHist("Electron_"+TString::Itoa(i,10)+"_PassesConvVeto", electron.PassesConvVeto(), 1., 0., 2., 2);
    FillHist("Electron_"+TString::Itoa(i,10)+"_MVA_"+EtaRegion, electron.MVA(), 1., -1., 1., 200);



    FillHist("Electron_Pt", electron.Pt(), 1., 0., 5000., 5000);
    FillHist("Electron_Eta", electron.Eta(), 1., -3., 3., 60);
    FillHist("Electron_dXY", fabs(electron.dXY()), 1., 0., 0.5, 500);
    FillHist("Electron_dXYSig", fabs(electron.dXYSig()), 1., 0., 10, 100);
    FillHist("Electron_dZ", fabs(electron.dZ()), 1., 0., 0.5, 500);
    FillHist("Electron_RelIso03", electron.PFRelIso(0.3), 1., 0., 1.0, 100);
    FillHist("Electron_RelMiniIso", electron.PFRelMiniIso(), 1., 0., 1.0, 100);
    FillHist("Electron_GsfCtfScPixChargeConsistency", electron.GsfCtfScPixChargeConsistency(), 1., 0., 2., 2);
    FillHist("Electron_PassesConvVeto", electron.PassesConvVeto(), 1., 0., 2., 2);
    FillHist("Electron_MVA_"+EtaRegion, electron.MVA(), 1., -1., 1., 200);


  }

  vector<bool> test_fatjets_IsLeptonInside;
  FillHist("FatJet_Size", test_fatjets.size(), 1., 0., 10., 10);
  for(unsigned int i=0; i<test_fatjets.size(); i++){

    snu::KFatJet jet = test_fatjets.at(i);

    bool this_leptoninside = false;
    for(unsigned int j=0; j<testmuons.size(); j++){
      if( jet.DeltaR( testmuons.at(j) )<0.4 ){
        this_leptoninside = true;
         break;
      }
    }
    for(unsigned int j=0; j<testelectrons.size(); j++){
      if( jet.DeltaR( testelectrons.at(j) )<0.4 ){
        this_leptoninside = true;
         break;
      }
    }
    test_fatjets_IsLeptonInside.push_back(this_leptoninside);

    FillHist("FatJet_"+TString::Itoa(i,10)+"_Pt", jet.Pt(), 1., 0., 5000., 5000);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Eta", jet.Eta(), 1., -3., 3., 60);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PrunedMass", jet.PrunedMass(), 1., 0., 5000., 5000);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau1", jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau2", jet.Tau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau3", jet.Tau3(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau21", jet.Tau2()/jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau31", jet.Tau3()/jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_Tau32", jet.Tau3()/jet.Tau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_SoftDropMass", jet.SoftDropMass(), 1., 0., 5000., 5000);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau1", jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau2", jet.PuppiTau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau3", jet.PuppiTau3(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau21", jet.PuppiTau2()/jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau31", jet.PuppiTau3()/jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_"+TString::Itoa(i,10)+"_PuppiTau32", jet.PuppiTau3()/jet.PuppiTau2(), 1., 0., 1.5, 150);

    FillHist("FatJet_Pt", jet.Pt(), 1., 0., 5000., 5000);
    FillHist("FatJet_Eta", jet.Eta(), 1., -3., 3., 60);
    FillHist("FatJet_PrunedMass", jet.PrunedMass(), 1., 0., 5000., 5000);
    FillHist("FatJet_Tau1", jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_Tau2", jet.Tau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_Tau3", jet.Tau3(), 1., 0., 1.5, 150);
    FillHist("FatJet_Tau21", jet.Tau2()/jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_Tau31", jet.Tau3()/jet.Tau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_Tau32", jet.Tau3()/jet.Tau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_SoftDropMass", jet.SoftDropMass(), 1., 0., 5000., 5000);
    FillHist("FatJet_PuppiTau1", jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_PuppiTau2", jet.PuppiTau2(), 1., 0., 1.5, 150);
    FillHist("FatJet_PuppiTau3", jet.PuppiTau3(), 1., 0., 1.5, 150);
    FillHist("FatJet_PuppiTau21", jet.PuppiTau2()/jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_PuppiTau31", jet.PuppiTau3()/jet.PuppiTau1(), 1., 0., 1.5, 150);
    FillHist("FatJet_PuppiTau32", jet.PuppiTau3()/jet.PuppiTau2(), 1., 0., 1.5, 150);

    FillHist("FatJet_Mass_vs_SoftDropMass", jet.M(), jet.SoftDropMass(), 1., 0., 5000., 5000, 0., 5000., 5000);

  }

  vector<bool> test_jets_IsLeptonInside, test_jets_PassPUMVA, test_jets_AwayFromFatJet;
  FillHist("Jet_Size", test_jets.size(), 1., 0., 10., 10);
  cout << "####### AK4 Jet #######" << endl;
  cout << "Pt\tEta\tPhi\tM" << endl;
  for(unsigned int i=0; i<test_jets.size(); i++){

    snu::KJet jet = test_jets.at(i);

    bool this_leptoninside = false;
    for(unsigned int j=0; j<testmuons.size(); j++){
      if( jet.DeltaR( testmuons.at(j) )<0.4 ){
        this_leptoninside = true;
         break;
      }
    }
    for(unsigned int j=0; j<testelectrons.size(); j++){
      if( jet.DeltaR( testelectrons.at(j) )<0.4 ){
        this_leptoninside = true;
         break;
      }
    }

    test_jets_IsLeptonInside.push_back(this_leptoninside);
    test_jets_PassPUMVA.push_back( jet.PassPileUpMVA("Loose") );
    test_jets_AwayFromFatJet.push_back( IsAwayFromFatJet(jet, test_fatjets) );

    cout << jet.Pt() << "\t" << jet.Eta() << "\t" << jet.Phi() << "\t" << jet.M() << "\t" << endl;
 

    FillHist("Jet_"+TString::Itoa(i,10)+"_Pt", jet.Pt(), 1., 0., 5000., 5000);
    FillHist("Jet_"+TString::Itoa(i,10)+"_Eta", jet.Eta(), 1., -3., 3., 60);

    FillHist("Jet_Pt", jet.Pt(), 1., 0., 5000., 5000);
    FillHist("Jet_Eta", jet.Eta(), 1., -3., 3., 60);

  }
  cout << endl;

  //==== Event Selection test
  FillCutFlow("NoCut", 1.);
  if(SUSYmuons.size()==2){

    FillCutFlow("DiMuon",1);

    snu::KMuon mu[2] = {SUSYmuons.at(0), SUSYmuons.at(1)};
    snu::KParticle N[2] = {mu[0], mu[1]};
    int n_Merged_FatJet[2] = {0, 0};
    double dR_Max_FatJet = 3.0;

    for(unsigned int i=0; i<test_fatjets.size(); i++){

      snu::KFatJet fatjet = test_fatjets.at(i);
      int CloserIndex = -1;
      double dR = 999;

      if(mu[0].DeltaR(fatjet) < mu[1].DeltaR(fatjet)){
        CloserIndex = 0;
        dR = mu[0].DeltaR(fatjet);
      }
      else{
        CloserIndex = 1;
        dR = mu[1].DeltaR(fatjet);
      }

      cout << "## Printing Current N["<<CloserIndex<<"] ##" << endl;
      cout << "Pt\tEta\tPhi\tM" << endl;
      cout << N[CloserIndex].Pt() << "\t" << N[CloserIndex].Eta() << "\t" << N[CloserIndex].Phi() << "\t" << N[CloserIndex].M() << endl;
      cout << "  Adding fatjet as follows" << endl;
      cout << "  Pt\tEta\tPhi\tM\tSD" << endl;
      cout << "  " << fatjet.Pt() << "\t" << fatjet.Eta() << "\t" << fatjet.Phi() << "\t" << fatjet.M() << "\t" << fatjet.SoftDropMass() << endl;

      N[CloserIndex] = N[CloserIndex]+fatjet;
      n_Merged_FatJet[CloserIndex]++;

      cout << "--> Resulting in" << endl;
      cout << "Pt\tEta\tPhi\tM" << endl;
      cout << N[CloserIndex].Pt() << "\t" << N[CloserIndex].Eta() << "\t" << N[CloserIndex].Phi() << "\t" << N[CloserIndex].M() << endl;
      cout << endl;

/*
      //==== Lepton Inside FatJet
      snu::KMuon CloseMuon = mu[CloserIndex];
      if(CloseMuon.DeltaR(fatjet) < 0.8){
        N[CloserIndex] = fatjet;
        n_Merged_FatJet[CloserIndex]++;
      }
      else if(dR < dR_Max_FatJet){
        N[CloserIndex] = N[CloserIndex]+fatjet;
        n_Merged_FatJet[CloserIndex]++;
      }
*/

    }

    double dR_Max_Jet = 4.0;
    for(unsigned int i=0; i<test_jets.size(); i++){

      snu::KJet jet = test_jets.at(i);

      if(!test_jets_AwayFromFatJet.at(i)) continue;
      if(test_jets_PassPUMVA.at(i)) continue;

      int CloserIndex = -1;
      double dR = 999;

      if(mu[0].DeltaR(jet) < mu[1].DeltaR(jet)){
        CloserIndex = 0;
        dR = mu[0].DeltaR(jet);
      }
      else{
        CloserIndex = 1;
        dR = mu[1].DeltaR(jet);
      }

      //==== Lepton Inside Jet
      snu::KMuon CloseMuon = mu[CloserIndex];

      //==== TEST FIXME

      cout << "## Printing Current N["<<CloserIndex<<"] ##" << endl;
      cout << "Pt\tEta\tPhi\tM" << endl;
      cout << N[CloserIndex].Pt() << "\t" << N[CloserIndex].Eta() << "\t" << N[CloserIndex].Phi() << "\t" << N[CloserIndex].M() << endl;
      cout << "  Adding jet as follows" << endl;
      cout << "  Pt\tEta\tPhi\tM" << endl;
      cout << "  " << jet.Pt() << "\t" << jet.Eta() << "\t" << jet.Phi() << "\t" << jet.M() << "\t" << endl;

      N[CloserIndex] = N[CloserIndex]+jet; 
      n_Merged_FatJet[CloserIndex]++;

      cout << "--> Resulting in" << endl;
      cout << "Pt\tEta\tPhi\tM" << endl;
      cout << N[CloserIndex].Pt() << "\t" << N[CloserIndex].Eta() << "\t" << N[CloserIndex].Phi() << "\t" << N[CloserIndex].M() << endl;

/*
      if(CloseMuon.DeltaR(jet) < 0.4){
        N[CloserIndex] = jet;
        n_Merged_FatJet[CloserIndex]++;
      }
      else if(dR < dR_Max_Jet){
        N[CloserIndex] = N[CloserIndex]+jet;
        n_Merged_FatJet[CloserIndex]++;
      }
*/
    }

    if(N[0].M()<10. || N[1].M()<10.){
      TruthPrintOut();
    }

    FillHist("TEST_N", N[0].M(), 1., 0., 3000., 3000);
    FillHist("TEST_N", N[1].M(), 1., 0., 3000., 3000);

    FillHist("TEST_n_Merged_FatJet", n_Merged_FatJet[0], 1., 0., 20., 20);
    FillHist("TEST_n_Merged_FatJet", n_Merged_FatJet[1], 1., 0., 20., 20);

  }

  return;



/*
  //==== Check Rochestor Syst
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HNGENT_TIGHT", true);
  for(unsigned int i=0; i<testmuons.size(); i++){
    double this_width = mcdata_correction->GetRochesterMomentumWidth(testmuons.at(i));
    cout << testmuons.at(i).Pt() << " -> " << this_width << endl;
    FillHist("RochDiff", this_width, 1., 0., 0.5, 500);
  }
  return;
*/
/*
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_TIGHT", true, 0, 2.4);
  std::vector< snu::KElectron > testelectrons = GetElectrons(true, true, "ELECTRON_HN_TIGHTv4", 0, 2.5);
  if(testmuons.size()==2){
    FillHist("DiMuon_lep1_pt", testmuons.at(0).Pt(), 1., 0., 200., 200);
    FillHist("DiMuon_lep2_pt", testmuons.at(1).Pt(), 1., 0., 200., 200);
  }
  return;
*/
/*
  //==== ZCR
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_TIGHT", true);
  std::vector< snu::KMuon > testmuons_veto = GetMuons("MUON_HN_VETO", true);
  std::vector< snu::KElectron > testelectrons = GetElectrons(true, true, "ELECTRON_HN_TIGHTv4");
  std::vector< snu::KElectron > testelectrons_veto = GetElectrons(true, true, "ELECTRON_HN_VETO");

  TString LepChannel = "";
  KLepton lep1, lep2;
  double IDSF = 1.;
  double TriggerSF = 1.;
  vector<TString> triggerlists;

  double test_weight = weight;
  if(!isData){
    test_weight *= mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    test_weight *= MCweight;
  }

  //==== DiMuon
  if(testmuons.size()==2 && testmuons_veto.size()==2 && testelectrons_veto.size()==0){

    LepChannel = "DiMuon";

    lep1 = testmuons.at(0);
    lep2 = testmuons.at(1);

    IDSF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", testmuons, 0);
    test_weight *= mcdata_correction->MuonTrackingEffScaleFactor(testmuons);

    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(testelectrons, "", testmuons, "MUON_HN_TIGHT", 0, 0, 0);
    double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(testelectrons, "", testmuons, "MUON_HN_TIGHT", 0, 1, 0);
    TriggerSF = trigger_eff_Data/trigger_eff_MC;


    triggerlists.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    triggerlists.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

    if( !( lep1.Pt()>20. && lep2.Pt() >10. ) ) return;

  }
  else if(testelectrons.size()==2 && testelectrons_veto.size()==2 && testmuons_veto.size()==0){

    LepChannel = "DiElectron";

    lep1 = testelectrons.at(0);
    lep2 = testelectrons.at(1);

    IDSF = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHTv4", testelectrons, 0);
    test_weight *= mcdata_correction->ElectronRecoScaleFactor(testelectrons);

    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(testelectrons, "", testmuons, "MUON_HN_TIGHT", 1, 0, 0);
    double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(testelectrons, "", testmuons, "MUON_HN_TIGHT", 1, 1, 0);
    TriggerSF = trigger_eff_Data/trigger_eff_MC;

    triggerlists.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    if( !( lep1.Pt()>25. && lep2.Pt() >15. ) ) return;
  }
  else{
    return;
  }
  test_weight *= WeightByTrigger(triggerlists, TargetLumi);
  double m_ll = (lep1+lep2).M();

  double weights[3] = {1., IDSF, IDSF*TriggerSF};
  TString str_weights[3] = {"NoSF", "IDSF", "IDandTrigSF"};

  cout << WeightByTrigger(triggerlists, TargetLumi) << endl;
  if(lep1.Charge()!=lep2.Charge() && fabs(m_ll-91.1876) < 20. && PassTriggerOR(triggerlists)){
    for(int i=0; i<3; i++){
      FillHist(LepChannel+"_"+str_weights[i]+"_Nevents", 0., test_weight*weights[i], 0., 1., 1);
      FillHist(LepChannel+"_"+str_weights[i]+"_m_ll", m_ll, test_weight*weights[i], 60., 120., 60);
      FillHist(LepChannel+"_"+str_weights[i]+"_lep1_pt", lep1.Pt(), test_weight*weights[i], 0., 200., 200);
      FillHist(LepChannel+"_"+str_weights[i]+"_lep1_eta", lep1.Eta(), test_weight*weights[i], -3., 3., 60);
      FillHist(LepChannel+"_"+str_weights[i]+"_lep2_pt", lep2.Pt(), test_weight*weights[i], 0., 200., 200);
      FillHist(LepChannel+"_"+str_weights[i]+"_lep2_eta", lep2.Eta(), test_weight*weights[i], -3., 3., 60);
    }

  }
  return;
*/
/*
  //==== Ghent Test
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HNGENT_TIGHT", true);
  std::vector< snu::KElectron > testelectrons = GetElectrons(true, true, "ELECTRON_GENT_TIGHT");
  std::vector< snu::KJet> testjets = GetJets("JET_GHENT");

  int testnbjets = 0;
  for(int j=0; j<testjets.size(); j++){
    if( IsBTagged(testjets.at(j), snu::KJet::CSVv2, snu::KJet::Loose, -1, 0) && fabs(testjets.at(j).Eta())<2.4 ){
      testnbjets++;
    }
  }
  if(testnbjets!=0) return;

  int n_muons = testmuons.size();
  int n_electrons = testelectrons.size();

  if(n_muons+n_electrons!=3) return; // lll
  if(n_muons<2) return; // mm+(m,e)
  if(testmuons.at(0).Pt() < 15) return;

  if(n_muons==2&&n_electrons==1){
    if(testmuons.at(1).Pt()<10){
      if(testmuons.at(0).Pt()<23) return;
    }
  }

  std::vector<KLepton> testlep;
  for(unsigned int j=0; j<testmuons.size(); j++){
    KLepton this_lep( testmuons.at(j) );
    testlep.push_back( this_lep );
  }
  for(unsigned int j=0; j<testelectrons.size(); j++){
    KLepton this_lep( testelectrons.at(j) );
    testlep.push_back( this_lep );
  }
  std::sort(testlep.begin(), testlep.end(), LeptonPtComparing);

  double testMET = eventbase->GetEvent().MET();
  bool UseEvent = (testMET<75) && ((testlep.at(0)+testlep.at(1)+testlep.at(2)).M()<80) && (testlep.at(0).Pt()<55);

  if(UseEvent){
    FillHist("test_Pass", 0., 1., 0., 1., 1);
  }
  bool HighUseEvent = ( !(testMET<75) || !((testlep.at(0)+testlep.at(1)+testlep.at(2)).M()<80) ) && (testlep.at(0).Pt()<55);
  if(HighUseEvent){
    if(testlep.at(1).Pt()>15&&testlep.at(2).Pt()>10){
      FillHist("test_Pass", 0., 1., 0., 1., 1);
    }
  }
  return;
*/




/*
  std::vector<snu::KFatJet> testfatjets = GetFatJets("FATJET_HN_loosest");
  for(unsigned int i=0; i<testfatjets.size(); i++){
    FillHist("test", testfatjets.at(i).L1JetCorr(), 1., 0., 2., 200);
  }
  return;
*/
/*
  std::vector< snu::KElectron > testelectrons = GetElectrons(false, true, "ELECTRON_HN_TIGHTv4");
  if(testelectrons.size()==2){
    FillHist("test", testelectrons.at(0).Pt(), 1., 0., 200., 200);
  }
  return;
*/
/*
  //==== Single Trigger OR

  vector<TString> trig_di, trig_si;
  trig_di.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  trig_di.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  trig_si.push_back("HLT_IsoMu24_v");

  bool Pass_trig_di = PassTriggerOR(trig_di);
  bool Pass_trig_si = PassTriggerOR(trig_si);

  if(Pass_trig_di){
    FillHist("TriggerPass", 0., 1., 0., 10., 10);
  }
  if(Pass_trig_si){
    FillHist("TriggerPass", 1., 1., 0., 10., 10);
  }
  if(Pass_trig_di || Pass_trig_si){
    FillHist("TriggerPass", 2., 1., 0., 10., 10);
  }
*/

/*
  //==== Isolation check for muon
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
   //==== Isolation check for electron
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
  //==== Conversion Rate
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  std::vector< snu::KElectron > testelectrons = GetElectrons(true, true, "ELECTRON_HN_FAKELOOSEv7_4");
  for(unsigned int i=0; i<testelectrons.size(); i++){
      
    snu::KElectron electron = testelectrons.at(i);
    
    int leptype =  GetLeptonType(electron, truthColl);
    bool IsTight = PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4");

    FillHist("Type_F0", leptype, 1., 20, -10., 10.);
    if(IsTight){
      FillHist("Type_F", leptype, 1., 20, -10., 10.);
    }

  };
  return;
*/
/*
  //==== Lepton Fake vs Type
  std::vector< snu::KElectron > testelectrons = GetElectrons(false, true, "ELECTRON_HN_FAKELOOSEv7_5");
  int testtypes[5] = {7, 22, 26, 27, 28};

  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  //==== FRTEST
  for(unsigned int i=0; i<testelectrons.size(); i++){

    snu::KElectron electron = testelectrons.at(i);

    int leptype =  GetLeptonType(electron, truthColl);

    //if(TruthMatched(electron, false)) continue;

    TString EtaRegion = "InnerBarrel";
    if(fabs(electron.SCEta()) > 1.479) EtaRegion = "EndCap";
    else if(fabs(electron.SCEta()) > 0.8) EtaRegion = "OuterBarrel";
    else EtaRegion = "InnerBarrel";

    double ptcone = CorrPt(testelectrons.at(i), 0.08);

    FillHist(EtaRegion+"_ELECTRON_FR_TYPE_F0", testelectrons.at(i).GetType(), 1., 0., 50., 50);
    FillHist(EtaRegion+"_ELECTRON_FR_MVA_F0", testelectrons.at(i).MVA(), 1., -1., 1., 200);
    FillHist(EtaRegion+"_ELECTRON_FR_nmisshit_F0", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
    FillHist(EtaRegion+"_ELECTRON_FR_JHTYPE_F0", leptype, 1., -10., 10., 20);

    if(PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4")){
      FillHist(EtaRegion+"_ELECTRON_FR_TYPE_F", testelectrons.at(i).GetType(), 1., 0., 50., 50);
      FillHist(EtaRegion+"_ELECTRON_FR_MVA_F", testelectrons.at(i).MVA(), 1., -1., 1., 200);
      FillHist(EtaRegion+"_ELECTRON_FR_nmisshit_F", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
      FillHist(EtaRegion+"_ELECTRON_FR_JHTYPE_F", leptype, 1., -10., 10., 20);
    }

    for(int j=0; j<5; j++){
      if(testelectrons.at(i).GetType()==testtypes[j]){
        FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_MVA_F0", testelectrons.at(i).MVA(), 1., -1., 1., 200);
        FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_dxy_F0", testelectrons.at(i).dxy(), 1., -0.2, 0.2, 400);
        FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_dz_F0", testelectrons.at(i).dz(), 1., -0.1, 0.1, 200);
        FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_nmisshit_F0", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
        if((PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4"))){
          FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_MVA_F", testelectrons.at(i).MVA(), 1., -1., 1., 200);
          FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_dxy_F", testelectrons.at(i).dxy(), 1., -0.2, 0.2, 400);
          FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_dz_F", testelectrons.at(i).dz(), 1., -0.1, 0.1, 200);
          FillHist(EtaRegion+"_"+TString::Itoa(testtypes[j],10)+"_ELECTRON_FR_nmisshit_F", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
        }
      }
    }

  }
  if(testelectrons.size()==2){
    if( testelectrons.at(0).Charge() == testelectrons.at(1).Charge() ){

      for(unsigned int i=0; i<testelectrons.size(); i++){

        snu::KElectron electron = testelectrons.at(i);

        int leptype =  GetLeptonType(electron, truthColl);

        //if(TruthMatched(electron, false)) continue;

        TString EtaRegion = "InnerBarrel";
        if(fabs(electron.SCEta()) > 1.479) EtaRegion = "EndCap";
        else if(fabs(electron.SCEta()) > 0.8) EtaRegion = "OuterBarrel";
        else EtaRegion = "InnerBarrel";

       EtaRegion = "SS_"+EtaRegion;

        double ptcone = CorrPt(testelectrons.at(i), 0.08);

        FillHist(EtaRegion+"_ELECTRON_FR_TYPE_F0", testelectrons.at(i).GetType(), 1., 0., 50., 50);
        FillHist(EtaRegion+"_ELECTRON_FR_MVA_F0", testelectrons.at(i).MVA(), 1., -1., 1., 200);
        FillHist(EtaRegion+"_ELECTRON_FR_nmisshit_F0", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
        FillHist(EtaRegion+"_ELECTRON_FR_JHTYPE_F0", leptype, 1., -10., 10., 20);

        if(PassID(testelectrons.at(i), "ELECTRON_HN_TIGHTv4")){
          FillHist(EtaRegion+"_ELECTRON_FR_TYPE_F", testelectrons.at(i).GetType(), 1., 0., 50., 50);
          FillHist(EtaRegion+"_ELECTRON_FR_MVA_F", testelectrons.at(i).MVA(), 1., -1., 1., 200);
          FillHist(EtaRegion+"_ELECTRON_FR_nmisshit_F", testelectrons.at(i).MissingHits(), 1., 0., 10., 10);
          FillHist(EtaRegion+"_ELECTRON_FR_JHTYPE_F", leptype, 1., -10., 10., 20);
        }

      }

    }
  }
*/
/*
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_LOOSEv7_SIP3", true);
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  //==== FRTEST
  for(unsigned int i=0; i<testmuons.size(); i++){
    int leptype =  GetLeptonType(testmuons.at(i), truthColl);
    FillHist("TEST_MUON_FR_TYPE_F0", testmuons.at(i).GetType(), 1., 0., 50., 50);
    FillHist("TEST_MUON_FR_LeptonType_F0", leptype, 1., -10., 10., 20);
    if(PassID(testmuons.at(i), "MUON_HN_TIGHT")){
      FillHist("TEST_MUON_FR_TYPE_F", testmuons.at(i).GetType(), 1., 0., 50., 50);
      FillHist("TEST_MUON_FR_LeptonType_F", leptype, 1., -10., 10., 20);
    }
  }
  return;
*/
/*
  //==== DZ filter

  std::vector< snu::KMuon > test_vetomuons = GetMuons("MUON_HN_VETO", true);
  std::vector< snu::KElectron > test_vetoelectrons = GetElectrons(true, true, "ELECTRON_HN_VETO");

  vector<TString> dimuontriggers_NonDZ, dimuontriggers_DZ;
  dimuontriggers_NonDZ.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  dimuontriggers_NonDZ.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  dimuontriggers_DZ.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  dimuontriggers_DZ.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  std::vector< snu::KMuon > testmuons = GetMuons("MUON_HN_TIGHT", true);

  float etabins_muon_dz[] = {0., 0.9, 1.2, 2.1, 2.4};
  if(testmuons.size() == 2 && test_vetomuons.size() == 2 && test_vetoelectrons.size() == 0){

    snu::KMuon mu1 = testmuons.at(0);
    snu::KMuon mu2 = testmuons.at(1);

    if( PassTriggerOR(dimuontriggers_NonDZ) ){

      bool UseEvent = false;
      UseEvent = ( fabs( (mu1+mu2).M() - 91.1876 ) < 15. ) && (mu1.Pt() > 20.) && (mu2.Pt() > 10.) && (mu1.Charge()!=mu2.Charge());

      if( UseEvent ){

        FillHist("DiMuonNonDZFired_onebin", 0., 1., 0., 1., 1);
        FillHist("DiMuonNonDZFired_mZ", (mu1+mu2).M(), 1., 50., 150., 100);
        FillHist("DiMuonNonDZFired_eta1_vs_eta2", fabs(mu1.Eta()), fabs(mu2.Eta()), 1., etabins_muon_dz, 4, etabins_muon_dz, 4);
        FillHist("DiMuonNonDZFired_dZ", fabs(mu1.dZ()-mu2.dZ()), 1., 0., 1., 1000);
        if( PassTriggerOR(dimuontriggers_DZ) ){
          FillHist("DiMuonDZFired_onebin", 0., 1., 0., 1., 1);
          FillHist("DiMuonDZFired_mZ", (mu1+mu2).M(), 1., 50., 150., 100);
          FillHist("DiMuonDZFired_eta1_vs_eta2", fabs(mu1.Eta()), fabs(mu2.Eta()), 1., etabins_muon_dz, 4, etabins_muon_dz, 4);
          FillHist("DiMuonDZFired_dZ", fabs(mu1.dZ()-mu2.dZ()), 1., 0., 1., 1000);
        }

      }

    }
  }

  vector<TString> dielectrontriggers_NonDZ, dielectrontriggers_DZ;
  dielectrontriggers_NonDZ.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  dielectrontriggers_DZ.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  std::vector< snu::KElectron > testelectrons = GetElectrons(true, true, "ELECTRON_HN_TIGHTv4");

  float etabins_electron_dz[] = {0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};
  if(testelectrons.size() == 2 && test_vetomuons.size() == 0 && test_vetoelectrons.size() == 2){

    snu::KElectron el1 = testelectrons.at(0);
    snu::KElectron el2 = testelectrons.at(1);

    if( PassTriggerOR(dielectrontriggers_NonDZ) ){

      bool UseEvent = false;
      UseEvent = ( fabs( (el1+el2).M() - 91.1876 ) < 15. ) && (el1.Pt() > 25.) && (el2.Pt() > 15.) && (el1.Charge()!=el2.Charge());

      if( UseEvent ){

        FillHist("DiElectronNonDZFired_onebin", 0., 1., 0., 1., 1);
        FillHist("DiElectronNonDZFired_mZ", (el1+el2).M(), 1., 50., 150., 100);
        FillHist("DiElectronNonDZFired_eta1_vs_eta2", fabs(el1.Eta()), fabs(el2.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
        FillHist("DiElectronNonDZFired_dZ", fabs(el1.dZ()-el2.dZ()), 1., 0., 1., 1000);
        if( PassTriggerOR(dielectrontriggers_DZ) ){
          FillHist("DiElectronDZFired_onebin", 0., 1., 0., 1., 1);
          FillHist("DiElectronDZFired_mZ", (el1+el2).M(), 1., 50., 150., 100);
          FillHist("DiElectronDZFired_eta1_vs_eta2", fabs(el1.Eta()), fabs(el2.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
          FillHist("DiElectronDZFired_dZ", fabs(el1.dZ()-el2.dZ()), 1., 0., 1., 1000);
        }

      }

    }
  }


  vector<TString> emutriggers_Mu8Ele23_NonDZ, emutriggers_Mu8Ele23_DZ;
  emutriggers_Mu8Ele23_NonDZ.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  emutriggers_Mu8Ele23_DZ.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(testmuons.size() == 1 && testelectrons.size() == 1 && test_vetomuons.size() == 1 && test_vetoelectrons.size() == 1){

    snu::KMuon mu = testmuons.at(0);
    snu::KElectron el = testelectrons.at(0);

    if( PassTriggerOR(emutriggers_Mu8Ele23_NonDZ) ){

      bool UseEvent = false;
      UseEvent = (mu.Pt() > 10.) && (el.Pt() > 25.) && (mu.Charge()!=el.Charge());

      if( UseEvent ){

        FillHist("Mu8Ele23_EMuNonDZFired_onebin", 0., 1., 0., 1., 1);
        FillHist("Mu8Ele23_EMuNonDZFired_mZ", (mu+el).M(), 1., 50., 150., 100);
        FillHist("Mu8Ele23_EMuNonDZFired_eta1_vs_eta2", fabs(mu.Eta()), fabs(el.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
        FillHist("Mu8Ele23_EMuNonDZFired_dZ", fabs(mu.dZ()-el.dZ()), 1., 0., 1., 1000);
        if( PassTriggerOR(emutriggers_Mu8Ele23_DZ) ){
          FillHist("Mu8Ele23_EMuDZFired_onebin", 0., 1., 0., 1., 1);
          FillHist("Mu8Ele23_EMuDZFired_mZ", (mu+el).M(), 1., 50., 150., 100);
          FillHist("Mu8Ele23_EMuDZFired_eta1_vs_eta2", fabs(mu.Eta()), fabs(el.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
          FillHist("Mu8Ele23_EMuDZFired_dZ", fabs(mu.dZ()-el.dZ()), 1., 0., 1., 1000);
        }
      }

    }

  }

  vector<TString> emutriggers_Mu23Ele8_NonDZ, emutriggers_Mu23Ele8_DZ;
  emutriggers_Mu23Ele8_NonDZ.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  emutriggers_Mu23Ele8_DZ.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(testmuons.size() == 1 && testelectrons.size() == 1 && test_vetomuons.size() == 1 && test_vetoelectrons.size() == 1){

    snu::KMuon mu = testmuons.at(0);
    snu::KElectron el = testelectrons.at(0);

    if( PassTriggerOR(emutriggers_Mu23Ele8_NonDZ) ){

      bool UseEvent = false;
      UseEvent = (mu.Pt() > 25.) && (el.Pt() > 10.) && (mu.Charge()!=el.Charge());

      if( UseEvent ){

        FillHist("Mu23Ele8_EMuNonDZFired_onebin", 0., 1., 0., 1., 1);
        FillHist("Mu23Ele8_EMuNonDZFired_mZ", (mu+el).M(), 1., 50., 150., 100);
        FillHist("Mu23Ele8_EMuNonDZFired_eta1_vs_eta2", fabs(mu.Eta()), fabs(el.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
        FillHist("Mu23Ele8_EMuNonDZFired_dZ", fabs(mu.dZ()-el.dZ()), 1., 0., 1., 1000);
        if( PassTriggerOR(emutriggers_Mu23Ele8_DZ) ){
          FillHist("Mu23Ele8_EMuDZFired_onebin", 0., 1., 0., 1., 1);
          FillHist("Mu23Ele8_EMuDZFired_mZ", (mu+el).M(), 1., 50., 150., 100);
          FillHist("Mu23Ele8_EMuDZFired_eta1_vs_eta2", fabs(mu.Eta()), fabs(el.Eta()), 1., etabins_electron_dz, 5, etabins_electron_dz, 5);
          FillHist("Mu23Ele8_EMuDZFired_dZ", fabs(mu.dZ()-el.dZ()), 1., 0., 1., 1000);
        }
      }

    }

  }



  return;
*/
/*
  //==== Colinear DiMuon test
  std::vector<snu::KMuon> testmuons = GetMuons("MUON_POG_TIGHT", true);
  std::vector<TString> testtriggers;
  testtriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  testtriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  testtriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  testtriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(PassTriggerOR(testtriggers)){
    if(testmuons.size()==2){
      snu::KMuon mu1 = testmuons.at(0);
      snu::KMuon mu2 = testmuons.at(1);
      bool temp_SS = false;
      if(mu1.Charge()==mu2.Charge()) temp_SS = true;

      if(temp_SS) return;

      //==== OS only..

      snu::KParticle dimu = mu1+mu2;
      double dRmm = mu1.DeltaR(mu2);   

      FillHist("m_dimuon", dimu.M(), 1., 0., 2000., 2000);
      FillHist("pt_dimuon", dimu.Pt(), 1., 0., 2000., 2000);
      FillHist("dR_dimuon", dRmm, 1., 0., 6., 600);
      FillHist("mu1pt", mu1.Pt(), 1., 0., 2000., 2000);
      FillHist("mu2pt", mu2.Pt(), 1., 0., 2000., 2000);

      if(dRmm<0.1){
        FillHist("Collinear_m_dimuon", dimu.M(), 1., 0., 50., 5000);
        FillHist("Collinear_pt_dimuon", dimu.Pt(), 1., 0., 2000., 2000);
        FillHist("Collinear_dR_dimuon", dRmm, 1., 0., 0.1, 100);
        FillHist("Collinear_mu1pt", mu1.Pt(), 1., 0., 2000., 2000);
        FillHist("Collinear_mu2pt", mu2.Pt(), 1., 0., 2000., 2000);
      }
  
      
    }
  }
  return;
*/

  //==== Signal PDF VERY FISRT

  float pileup_reweight=(1.0);
  if(!isData){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  ForTree_PdfWeights->clear();
  ForTree_ScaleWeights->clear();
  snu::KEvent Evt = eventbase->GetEvent();
  int N_RunNumber = Evt.RunNumber();
  int N_EventNumber = Evt.EventNumber();
  bool SignalPDFSyst = std::find(k_flags.begin(), k_flags.end(), "SignalPDFSyst") != k_flags.end();

  //==== PDF weight has MCweight applied
  //==== MCweight will be applied later
  double weight_for_pdf_den = weight*pileup_reweight*35863.3;
  for(unsigned int i=0; i<Evt.PdfWeights().size(); i++){
    ForTree_PdfWeights->push_back(Evt.PdfWeights().at(i));
    if(SignalPDFSyst){
      FillHist("ForTree_PdfWeights", i, weight_for_pdf_den*Evt.PdfWeights().at(i), 0., 1.*Evt.PdfWeights().size(), Evt.PdfWeights().size());
    }
  }
  for(unsigned int i=0; i<Evt.ScaleWeights().size(); i++){
    ForTree_ScaleWeights->push_back(Evt.ScaleWeights().at(i));
    if(SignalPDFSyst){
      FillHist("ForTree_ScaleWeights", i, weight_for_pdf_den*Evt.ScaleWeights().at(i), 0., 1.*Evt.ScaleWeights().size(), Evt.ScaleWeights().size());
    }
  }
  if(SignalPDFSyst){
    FillHist("ForTree_Central", 0, weight_for_pdf_den*MCweight, 0., 1., 1);
    return;
  }
  n_vtx = Evt.nVertices();

  //==== Initializing

  AUTO_N_syst = 0;
  AUTO_it_syst = 0;
  AUTO_syst_type = "";

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

  //==== DiMuon channel Trigger

  if(isData){
    if(k_channel.Contains("DoubleMuon")){
      triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    }
    if(k_channel.Contains("SingleMuon")){
      triggerlist_DiMuon.push_back("HLT_IsoMu24_v");
      triggerlist_DiMuon.push_back("HLT_IsoTkMu24_v");
    }
  }
  else{
    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    //triggerlist_DiMuon.push_back("HLT_IsoMu24_v");
    //triggerlist_DiMuon.push_back("HLT_IsoTkMu24_v");
  }

  vector<TString> triggerlist_DiMuon_Mu17Mu8, triggerlist_DiMuon_Mu24;
  triggerlist_DiMuon_Mu17Mu8.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  triggerlist_DiMuon_Mu17Mu8.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  triggerlist_DiMuon_Mu17Mu8.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon_Mu17Mu8.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon_Mu24.push_back("HLT_IsoMu24_v");
  triggerlist_DiMuon_Mu24.push_back("HLT_IsoTkMu24_v");

  vector<TString> triggerlist_DiMuon_PeriodH;
  triggerlist_DiMuon_PeriodH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon_PeriodH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  //==== DiElectron channel Trigger

  if(isData){
    if(k_channel.Contains("DoubleEG")){
      triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    }
    if(k_channel.Contains("SingleElectron")){
      triggerlist_DiElectron.push_back("HLT_Ele27_WPTight_Gsf_v");
    }
  }
  else{
    triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    //triggerlist_DiElectron.push_back("HLT_Ele27_WPTight_Gsf_v");
  }

  vector<TString> triggerlist_DiElectron_Ele23Ele12, triggerlist_DiElectron_Ele27;
  triggerlist_DiElectron_Ele23Ele12.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_DiElectron_Ele27.push_back("HLT_Ele27_WPTight_Gsf_v");

  //==== EMu channel Trigger

  if(isData){
    if(k_channel.Contains("MuonEG")){
      triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
      triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
      triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
    }
    if(k_channel.Contains("SingleMuon")){
      triggerlist_EMu.push_back("HLT_IsoMu24_v");
      triggerlist_EMu.push_back("HLT_IsoTkMu24_v");
    }
    if(k_channel.Contains("SingleElectron")){
      triggerlist_EMu.push_back("HLT_Ele27_WPTight_Gsf_v");
    }
  }
  else{
    triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
    //triggerlist_EMu.push_back("HLT_IsoMu24_v");
    //triggerlist_EMu.push_back("HLT_IsoTkMu24_v");
    //triggerlist_EMu.push_back("HLT_Ele27_WPTight_Gsf_v");
  }

  vector<TString> triggerlist_EMu_EMu;
  triggerlist_EMu_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_EMu_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  triggerlist_EMu_PeriodH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_EMu_PeriodH.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  w_cutflow["DiMuon"] = weight*WeightByTrigger(triggerlist_DiMuon, TargetLumi)*pileup_reweight;
  w_cutflow["DiElectron"] = weight*WeightByTrigger(triggerlist_DiElectron, TargetLumi)*pileup_reweight;
  w_cutflow["EMu"] = weight*35863.3*pileup_reweight; //FIXME

  FillCutFlowByName("DiMuon", "NoCut", w_cutflow["DiMuon"], isData);
  FillCutFlowByName("DiElectron", "NoCut", w_cutflow["DiElectron"], isData);
  FillCutFlowByName("EMu", "NoCut", w_cutflow["EMu"], isData);

  //======================
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);

  FillCutFlowByName("DiMuon", "MET", w_cutflow["DiMuon"], isData);
  FillCutFlowByName("DiElectron", "MET", w_cutflow["DiElectron"], isData);
  FillCutFlowByName("EMu", "MET", w_cutflow["EMu"], isData);

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

  FillCutFlowByName("DiMuon", "PV", w_cutflow["DiMuon"], isData);
  FillCutFlowByName("DiElectron", "PV", w_cutflow["DiElectron"], isData);
  FillCutFlowByName("EMu", "PV", w_cutflow["EMu"], isData);

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "DoMCClosure") != k_flags.end();
  bool AllLepton = std::find(k_flags.begin(), k_flags.end(), "AllLepton") != k_flags.end();
  bool DoNoSF = std::find(k_flags.begin(), k_flags.end(), "DoNoSF") != k_flags.end();
  bool DoIDOnly = std::find(k_flags.begin(), k_flags.end(), "DoIDOnly") != k_flags.end();
  bool DoConversion = std::find(k_flags.begin(), k_flags.end(), "DoConversion") != k_flags.end();
  bool MCFakeSubtract = std::find(k_flags.begin(), k_flags.end(), "MCFakeSubtract") != k_flags.end();
  bool MCCFSubtract = std::find(k_flags.begin(), k_flags.end(), "MCCFSubtract") != k_flags.end();

  if(DoConversion) weight *= 0.5;

  if(DoMCClosure) DoNoSF = true;

  bool KeepFakeLepton = false;
  if(DoMCClosure){
    KeepFakeLepton = true;
    MuFR_key = "QCD";
    ElFR_key = "QCD";
  }
  if(AllLepton) KeepFakeLepton = true;

  bool KeepCFLepton = false;
  if(MCCFSubtract){
    KeepCFLepton = true;
  }

  bool DoFakeSyst = std::find(k_flags.begin(), k_flags.end(), "DoFakeSyst") != k_flags.end();
  bool LooseSampleFakeJetPt = std::find(k_flags.begin(), k_flags.end(), "LooseSampleFakeJetPt") != k_flags.end();

  //==== Muons
  std::vector<snu::KMuon> muons_loosest = GetMuons(MuonLooseID_loosest, KeepFakeLepton);
  std::vector<snu::KMuon> muons_veto_loosest = GetMuons(MuonVetoID_loosest, true);

  //==== Electrons
  std::vector<snu::KElectron> electrons_loosest = GetElectrons(KeepCFLepton, KeepFakeLepton, ElectronLooseID_loosest);
  std::vector<snu::KElectron> electrons_veto_loosest = GetElectrons(true, true, ElectronVetoID_loosest);

/*
  //=====================
  //FIXME
  std::vector<snu::KElectron> temp_electrons;

  electrons_loosest = GetElectrons(false, KeepFakeLepton, "ELECTRON_HN_TIGHTv4");
  for(unsigned int i=0; i<electrons_loosest.size(); i++){
    if( fabs(electrons_loosest.at(i).SCEta()) < 1.444 ) temp_electrons.push_back( electrons_loosest.at(i) );
  }
  std::vector<snu::KElectron> temp_pogelectrons = GetElectrons(false, KeepFakeLepton, "ELECTRON_POG_TIGHT_Charge");
  for(unsigned int i=0; i<temp_pogelectrons.size(); i++){
    if( fabs(temp_pogelectrons.at(i).SCEta()) > 1.444 ) temp_electrons.push_back( temp_pogelectrons.at(i) );
  }
  electrons_loosest.clear();
  for(unsigned int i=0; i<temp_electrons.size(); i++){
    electrons_loosest.push_back( temp_electrons.at(i) );
  }
  //=====================
*/

  bool RunningChargeFlipData = k_running_chargeflip && (isData);

  //==== Jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto_loosest", 10., 5.);
  //std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJetsWFT("JET_HN_eta5_nolepveto_loosest", 10., 5.);

  //==== FatJets
  std::vector<snu::KFatJet> fatjets_loosest_UnSmeared = GetFatJets("FATJET_HN_loosest");
  std::vector<snu::KFatJet> fatjets_loosest = GetCorrectedFatJet(fatjets_loosest_UnSmeared); // Smear both energy and mass

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
    weight*=GetKFactor();
  }

  bool DoFRBJET = std::find(k_flags.begin(), k_flags.end(), "DoFRBJET") != k_flags.end();
  snu::KJet::WORKING_POINT btag_wp = snu::KJet::Medium;

  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  RunNtp = std::find(k_flags.begin(), k_flags.end(), "RunNtp") != k_flags.end();

  //==== Define Analysis Region

  std::vector< TString > Suffixs;
  std::vector< std::vector<TString> > Triggers;

  bool RunningNonPromptData = NonPromptRun && (isData||DoMCClosure||MCFakeSubtract);

  //==== Make Suffix 

  Suffixs.push_back("DiMuon");
  Triggers.push_back(triggerlist_DiMuon);

  Suffixs.push_back("DiElectron");
  Triggers.push_back(triggerlist_DiElectron);

  Suffixs.push_back("EMu");
  Triggers.push_back(triggerlist_EMu);

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
  //==== 12) Jet Mass Scale
  //==== 13) Jet Mass Res
  //==== 14) Tau21
  //====================================

  int N_sys = 2*14+1;
  int it_sys_start = 0;

  if( isData || DoMCClosure || MCFakeSubtract ){
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
    else if(it_sys==23){
      this_syst = "_JetMass_up";
    }
    else if(it_sys==24){
      this_syst = "_JetMass_down";
    }
    else if(it_sys==25){
      this_syst = "_JetMassRes_up";
    }
    else if(it_sys==26){
      this_syst = "_JetMassRes_down";
    }
    else if(it_sys==27){
      this_syst = "_Tau21_up";
    }
    else if(it_sys==28){
      this_syst = "_Tau21_down";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    AUTO_N_syst = N_sys;
    AUTO_it_syst = it_sys;
    AUTO_syst_type = this_syst;

    if(this_syst==""){
      NowRunningCentral = true;
    }
    else{
      NowRunningCentral = false;
    }

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
        double scale = 1.+mcdata_correction->GetRochesterMomentumWidth(this_muon);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*scale, this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/scale;
        this_muon.SetRelIso(0.4, new_RelIso);
        if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
        snu::KMuon this_muon = muons_veto_loosest.at(j);
        double scale = 1.+mcdata_correction->GetRochesterMomentumWidth(this_muon);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*scale, this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/scale;
        this_muon.SetRelIso(0.4, new_RelIso);
        if( this_muon.Pt() >= 5. && new_RelIso < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
      }
    }
    else if(this_syst == "_MuonEn_down"){
      //==== Signal Leptons
      for(unsigned int j=0; j<muons_loosest.size(); j++){
        snu::KMuon this_muon = muons_loosest.at(j);
        double scale = 1.-mcdata_correction->GetRochesterMomentumWidth(this_muon);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*scale, this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/scale;
        this_muon.SetRelIso(0.4, new_RelIso);
        if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
        snu::KMuon this_muon = muons_veto_loosest.at(j);
        double scale = 1.-mcdata_correction->GetRochesterMomentumWidth(this_muon);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*scale, this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/scale;
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
        //if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso && fabs(this_electron.SCEta()) < 1.444) electrons.push_back( this_electron );
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
        //if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso && fabs(this_electron.SCEta()) < 1.444) electrons.push_back( this_electron );
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
        //if( this_electron.Pt() >= 10. && this_electron.PFRelIso(0.3) < this_ElectronLooseRelIso && fabs(this_electron.SCEta()) < 1.444) electrons.push_back( this_electron );
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
    //==== 3) If non-prompt run, now pt is pt-cone. So this will be also corrected
    //==== Both cases, we can correct MET (w.r.t. muon) using
    //==== AnalyzerCore::JSCorrectedMETRochester(std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi)
    JSCorrectedMETRochester(muons, MET, METphi);

    //==================
    //==== Make FatJet
    //==================

    vector<snu::KFatJet> fatjets;
    if( (this_syst == "_JetEn_up") || (this_syst == "_JetRes_up") || (this_syst == "_JetMass_up") || (this_syst == "_JetMassRes_up") ){
      for(unsigned int j=0; j<fatjets_loosest.size(); j++){
        snu::KFatJet this_jet = fatjets_loosest.at(j);
        double this_scaling = 1.;
        double this_mass_scaling = 1.;
        if(this_syst == "_JetEn_up") this_scaling = this_jet.ScaledUpEnergy();
        if(this_syst == "_JetRes_up") this_scaling = this_jet.SmearedResUp()/this_jet.SmearedRes();
        if(this_syst == "_JetMass_up") this_mass_scaling = this_jet.ScaledMassUp();
        if(this_syst == "_JetMassRes_up") this_mass_scaling = this_jet.SmearedMassResUp()/this_jet.SmearedMassRes();
        this_jet *= this_scaling;
        this_jet.SetPrunedMass(this_jet.PrunedMass()*this_mass_scaling);

        if(JSFatJetID(this_jet) && this_jet.Pt() >= 200.) fatjets.push_back(this_jet);
      }
    }
    else if( (this_syst == "_JetEn_down") || (this_syst == "_JetRes_down") || (this_syst == "_JetMass_down") || (this_syst == "_JetMassRes_down") ){
      for(unsigned int j=0; j<fatjets_loosest.size(); j++){
        snu::KFatJet this_jet = fatjets_loosest.at(j);

        double this_scaling = 1.;
        double this_mass_scaling = 1.;
        if(this_syst == "_JetEn_down") this_scaling = this_jet.ScaledDownEnergy();
        if(this_syst == "_JetRes_down") this_scaling = this_jet.SmearedResDown()/this_jet.SmearedRes();
        if(this_syst == "_JetMass_down") this_mass_scaling = this_jet.ScaledMassDown();
        if(this_syst == "_JetMassRes_down") this_mass_scaling = this_jet.SmearedMassResDown()/this_jet.SmearedMassRes();
        this_jet *= this_scaling;
        this_jet.SetPrunedMass(this_jet.PrunedMass()*this_mass_scaling);

        if(JSFatJetID(this_jet) && this_jet.Pt() >= 200.) fatjets.push_back(this_jet);
      }
    }
    else{
      for(unsigned int j=0; j<fatjets_loosest.size(); j++){
        snu::KFatJet this_jet = fatjets_loosest.at(j);

        if(JSFatJetID(this_jet) && this_jet.Pt() >= 200.) fatjets.push_back(this_jet);
      }
    }
    //cout << "fatjets.size() = " << fatjets.size() << endl;
    double FatJetTau21_SF = 1.;

    int FatJetTau21_dir = 0;
    if(this_syst=="_Tau21_up") FatJetTau21_dir = 1;
    if(this_syst=="_Tau21_down") FatJetTau21_dir = -1;

    for(unsigned int j=0; j<fatjets.size(); j++){
      FatJetTau21_SF *= GetFatJetSF(fatjets.at(j), 0.6, FatJetTau21_dir);
    }

    JSCorrectedMETFatJet(fatjets, MET, METphi);

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

    std::vector<snu::KJet> jets_eta5; // eta < 5.0, lepton-veto, away from fatjets
    std::vector<snu::KJet> jets; // eta < 2.5, lepton-veto, away from fatjets
    std::vector<snu::KJet> jets_nonfatjetveto;
    std::vector<snu::KJet> jets_InSideFatJet; // If jets inside fatjet, remove it's smearing from MET. Because FatJet smearing is already propagted to MET
    std::vector<snu::KJet> jets_nolepveto; // eta < 2.5, NO lepton-veto
    std::vector<snu::KJet> jets_fwd; // 2.5 < eta < 5, lepton-veto => to make forward

    for(unsigned int j=0; j<jets_eta5_nolepveto.size(); j++){

      double NormalJetMaxEta = 2.7;

      snu::KJet this_jet = jets_eta5_nolepveto.at(j);
      bool IsNormalJet = fabs( this_jet.Eta() ) < NormalJetMaxEta;
      bool IsForwardJet = fabs( this_jet.Eta() ) >= NormalJetMaxEta;
      bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
      bool awayfromfatjet = IsAwayFromFatJet(this_jet, fatjets);
      bool PassPUID = this_jet.PassPileUpMVA("Loose");

      if(!lepinside && awayfromfatjet) jets_eta5.push_back( this_jet );
      if(IsNormalJet && !lepinside && PassPUID && awayfromfatjet) jets.push_back( this_jet );
      if(IsNormalJet && !lepinside && PassPUID) jets_nonfatjetveto.push_back( this_jet );
      if(IsNormalJet && !lepinside && !awayfromfatjet) jets_InSideFatJet.push_back( this_jet );
      if(IsNormalJet) jets_nolepveto.push_back( this_jet );
      if(IsForwardJet && !lepinside) jets_fwd.push_back( this_jet );

    }
    JSCorrectedMETJetInsideFatJet(jets_InSideFatJet, MET, METphi);

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
    std::vector<int> NearBjet, NearBjet_TRUE;

    int NPromptTight(0);

    int n_veto_muons = muons_veto.size();
    int n_triLoose_muons = muons.size();
    int n_triTight_muons(0);
    for(unsigned int j=0; j<muons.size(); j++){

      int hasclosebjet = HasCloseBjet(muons.at(j), jets_nolepveto, btag_wp);

      NearBjet_TRUE.push_back(hasclosebjet);

      if(DoFRBJET) NearBjet.push_back( hasclosebjet );
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
    bool HasLooseCF = false; // for cf-fake correction
    double LooseCFSF = 1.;
    for(unsigned int j=0; j<electrons.size(); j++){

      int hasclosebjet = HasCloseBjet(electrons.at(j), jets_nolepveto, btag_wp);

      NearBjet_TRUE.push_back(hasclosebjet);

      if(DoFRBJET) NearBjet.push_back( hasclosebjet );
      else NearBjet.push_back( -1 );

      if(TruthMatched(electrons.at(j), false)) NPromptTight++;
      if(PassID(electrons.at(j), ElectronTightID)){
        isT.push_back(true);
        electrons_tight.push_back(electrons.at(j));
        n_triTight_electrons++;
      }
      else{
        isT.push_back(false);

        if(MCIsCF(electrons.at(j))){
          HasLooseCF = true;


          if(fabs(electrons.at(j).SCEta())<1.4442) LooseCFSF = 0.7370;
          else LooseCFSF = 0.9190;

          //==== With Wrong CF type definition..
          //if(fabs(electrons.at(j).SCEta())<1.4442) LooseCFSF = 0.7138;
          //else LooseCFSF = 0.9151;

          //==== CF/(CF+40)
          //if(fabs(electrons.at(j).SCEta())<1.4442) LooseCFSF = 0.7233*0.8444;
          //else LooseCFSF = 0.9120*0.8968;

        }

      }
    }

    if(MCCFSubtract){
      if(!HasLooseCF) continue;
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

    //==== ID Scale Factros before pt->pt-cone / multiplied later

    int MuonIDDir = 0;
    if(this_syst=="_MuonIDSF_up") MuonIDDir = +1;
    else if(this_syst=="_MuonIDSF_down") MuonIDDir = -1;
    else MuonIDDir = 0;

    double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", muons, MuonIDDir);
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muons);
    if(DoNoSF) muon_id_iso_sf = 1.;

    int ElectronIDDir = 0;
    if(this_syst=="_ElectronIDSF_up") ElectronIDDir = +1;
    else if(this_syst=="_ElectronIDSF_down") ElectronIDDir = -1;
    else ElectronIDDir = 0;

    double electron_sf = mcdata_correction->ElectronScaleFactor(ElectronTightID, electrons, ElectronIDDir);
    if(DoNoSF) electron_sf = 1.;
    double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);

    //==== Non-prompt, change pt to pt-cone
    if(NonPromptRun){

      for(unsigned int j=0; j<muons.size(); j++){
        snu::KMuon this_muon = muons.at(j);
        this_muon.SetPtEtaPhiM( CorrPt(this_muon, 0.07), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        muons.at(j) = this_muon;
      }
      std::sort(muons.begin(), muons.end(), MuonPtComparing);

      for(unsigned int j=0; j<electrons.size(); j++){
        snu::KElectron this_electron = electrons.at(j);
        this_electron.SetPtEtaPhiM( CorrPt(this_electron, 0.08), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        electrons.at(j) = this_electron;
      }
      std::sort(electrons.begin(), electrons.end(), ElectronPtComparing);

    }

    for(unsigned int i=0; i<Suffixs.size(); i++){

      TString Suffix = Suffixs.at(i);

      FillCutFlowByName(Suffix, "MET_PV", w_cutflow[Suffix], isData);

      //==== Trigger pass
      if(!PassTriggerOR( Triggers.at(i) )) continue;

      double DiMuon_MCTriggerWeight = 0.;
      double EMu_MCTriggerWeight = 0.;
      if(Suffix.Contains("DiMuon")){
        //==== Data
        if(isData){

          DiMuon_MCTriggerWeight = 1.;

          //==== periodH : DZ must be fired. (nonDZ were off in this period)

          //==== Period H
          if(GetDataPeriod() == 7){
            if(!PassTriggerOR( triggerlist_DiMuon_PeriodH )) continue;
          }
        }
        else{

          //==== If single lepton triggers were fired,
          //==== Use full lumi.
          //==== But, if ony EMu triggers were fired,
          //==== nonDZ fired : add lumi of periodBtoG
          //==== DZ fired : add lumi of periodH
          if(PassTriggerOR( triggerlist_DiMuon )) DiMuon_MCTriggerWeight += 27257.617;
          if(PassTriggerOR( triggerlist_DiMuon_PeriodH )) DiMuon_MCTriggerWeight += 8605.69;

        }
      }
      if(Suffix.Contains("EMu")){
        //==== Data
        if(isData){

          EMu_MCTriggerWeight = 1.;

          //==== periodBtoG : nonDZ must be fired. (DZ were off in this period)
          //==== periodH : DZ must be fired. (nonDZ were off in this period)

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

          //==== nonDZ fired : add lumi of periodBtoG
          //==== DZ fired : add lumi of periodH

          if(PassTriggerOR( triggerlist_EMu_PeriodBtoG )) EMu_MCTriggerWeight += 27257.617;
          if(PassTriggerOR( triggerlist_EMu_PeriodH )) EMu_MCTriggerWeight += 8605.69;

        }
      }
      FillCutFlowByName(Suffix, "MET_PV_Trig", w_cutflow[Suffix], isData);

      if(Suffix.Contains("DiMuon")){
        isTT = isTwoMuon_TT;
        isLOOSE = isTwoMuon_Loose;
        isNoExtra = n_veto_muons == 2;
        isNoExtraOtherFlavour = n_veto_electrons == 0;
      }
      else if(Suffix.Contains("DiElectron")){
        isTT = isTwoElectron_TT;
        isLOOSE = isTwoElectron_Loose;
        isNoExtra = n_veto_electrons == 2;
        isNoExtraOtherFlavour = n_veto_muons == 0;
      }
      else if(Suffix.Contains("EMu")){
        isTT = isEMu_TT;
        isLOOSE = isEMu_Loose;
        isNoExtra = (n_veto_muons == 1) && (n_veto_electrons == 1);
        isNoExtraOtherFlavour = true;
      }
      else{
        cout << "Suffix Wrong" << endl;
        exit(EXIT_FAILURE);
      }

      isTT = isTT && !RunningNonPromptData;
      isLOOSE = isLOOSE && RunningNonPromptData;
      
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

        bool PtOkay = false;

        if(PassTriggerOR(triggerlist_DiMuon_Mu17Mu8)){
          if(isTT){
            if((muons_tight.at(0).Pt() > 20.) && (muons_tight.at(1).Pt() > 10.)) PtOkay = true;
          }
          else{
            if((muons.at(0).Pt() > 20.) && (muons.at(1).Pt() > 10.)) PtOkay = true;
          }
        }
/*
        if(PassTriggerOR(triggerlist_DiMuon_Mu24)){
          if(isTT){
            if( (muons_tight.at(0).Pt() > 26.) ) PtOkay = true;
          }
          else{
            if( (muons.at(0).Pt() > 26.) ) PtOkay = true;
          }
        }
*/

        if( !PtOkay ) continue;


      }
      if(Suffix.Contains("DiElectron")){

        bool PtOkay = false;

        if(PassTriggerOR(triggerlist_DiElectron_Ele23Ele12)){
          if(isTT){
            if((electrons_tight.at(0).Pt() > 25.) && (electrons_tight.at(1).Pt() > 15.)) PtOkay = true;
          }
          else{
            if((electrons.at(0).Pt() > 25.) && (electrons.at(1).Pt() > 15.)) PtOkay = true;
          }
        }
/*
        if(PassTriggerOR(triggerlist_DiElectron_Ele27)){
          if(isTT){
            if( (electrons_tight.at(0).Pt() > 30.) ) PtOkay = true;
          }
          else{
            if( (electrons.at(0).Pt() > 30.) ) PtOkay = true;
          }
        }
*/

        if( !PtOkay ) continue;

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
/*
        if(PassTriggerOR(triggerlist_DiMuon_Mu24)){
          
          if(isTT){
            if( (muons_tight.at(0).Pt() > 26.) ) PtOkay = true;
          }
          else{ 
            if( (muons.at(0).Pt() > 26.) ) PtOkay = true;
          }
        }
        if(PassTriggerOR(triggerlist_DiElectron_Ele27)){

          if(isTT){
            if( (electrons_tight.at(0).Pt() > 30.) ) PtOkay = true;
          }
          else{
            if( (electrons.at(0).Pt() > 30.) ) PtOkay = true;
          }
        }
*/

        if( !PtOkay ) continue;

      }

      //==== DiMuon-DoubleMuon PD / ...
      if(isData && !k_channel.Contains("DoubleMuon_CF")){

        //==== DoubleMuon PD
        if(k_channel.Contains("DoubleMuon")){
          if(!Suffix.Contains("DiMuon")) continue;
        }
        //==== DoubleEG PD
        if(k_channel.Contains("DoubleEG")){
          if(!Suffix.Contains("DiElectron")) continue;
        }
        //===== MuonEG PD
        if(k_channel.Contains("MuonEG")){
          if(!Suffix.Contains("EMu")) continue;
        }

        //==== SingleMuon PD
        if(k_channel.Contains("SingleMuon")){

          //==== For DiMuon channel
          //==== veto Mu17Mu8 fired event
          if(Suffix.Contains("DiMuon")){
            if(PassTriggerOR(triggerlist_DiMuon_Mu17Mu8)) continue;
          }

          //=== Do not fill DiElectron channel
          if(Suffix.Contains("DiElectron")) continue;

          //==== For EMu channel
          //==== veto Mu23Ele8 fired event
          if(Suffix.Contains("EMu")){
            if(PassTriggerOR(triggerlist_EMu_EMu)) continue;
          }

        }

        //==== SingleElectron PD
        if(k_channel.Contains("SingleElectron")){

          //==== Do not fill DiMuon channel
          if(Suffix.Contains("DiMuon")) continue;

          //==== For DiElectron channel
          //==== veto Ele23Ele12 fired event
          if(Suffix.Contains("DiElectron")){
            if(PassTriggerOR(triggerlist_DiElectron_Ele23Ele12)) continue;
          }

          //==== For EMu channel
          //==== veto Mu23Ele8 fired event
          //==== veto Mu24 fired event
          if(Suffix.Contains("EMu")){
            if(PassTriggerOR(triggerlist_EMu_EMu)||PassTriggerOR(triggerlist_DiMuon_Mu24)) continue;
          }
        }

      } // isData

      if(MCCFSubtract){
        if(!Suffix.Contains("DiElectron")) continue;
      }

      double trigger_ps_weight(1.);
      if(Suffix.Contains("DiMuon")) trigger_ps_weight = DiMuon_MCTriggerWeight;
      else if(Suffix.Contains("EMu")) trigger_ps_weight = EMu_MCTriggerWeight;
      else              trigger_ps_weight = WeightByTrigger(Triggers.at(i), TargetLumi);

      double this_weight = weight*trigger_ps_weight*FatJetTau21_SF;

      this_weight *= muon_id_iso_sf*MuTrkEffSF;
      this_weight *= electron_sf*electron_RecoSF;

      if(MCCFSubtract) this_weight *= LooseCFSF;

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
        //cout << "TriggerSFDir = " << TriggerSFDir << endl;
        //cout << "1) Data" << endl;
        //cout << "trigger_eff_Data = " << trigger_eff_Data << endl;
        //cout << "2) MC" << endl;
        //cout << "trigger_eff_MC = " << trigger_eff_MC << endl;
        //cout << "=> sf = " << trigger_sf << endl;

      }
      if(!isData && Suffix.Contains("EMu")){

        double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 2, 0, TriggerSFDir);
        double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muons, "MUON_HN_TIGHT", 2, 1, -1*TriggerSFDir);
        trigger_sf = trigger_eff_Data/trigger_eff_MC;
      }

      if(DoNoSF||DoIDOnly) trigger_sf = 1.;

      this_weight *= trigger_sf;

/*
      if(fabs(this_weight*pileup_reweight) < 1){
       //WTF
        cout << "FatJetTau21_SF = " << FatJetTau21_SF << endl;
        cout << "trigger_sf = " << trigger_sf << endl;
        cout << "electron_sf = " << electron_sf << endl;
        cout << "electron_RecoSF = " << electron_RecoSF << endl;
        cout << "muon_id_iso_sf = " << muon_id_iso_sf << endl;
        cout << "MuTrkEffSF = " << MuTrkEffSF << endl;
        cout << "trigger_ps_weight = " << trigger_ps_weight << endl;
        cout << "pileup_reweight = " << pileup_reweight << endl;
        cout << "n_vtx = " << n_vtx << endl;
        cout << "pu_up = " << mcdata_correction->CatPileupWeight(eventbase->GetEvent(),+1) << endl;
        cout << "pu_down = " << mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1) << endl;
        cout << "=> weight = " << this_weight << endl;
      }
*/

      //==== pileup reweight
      double purew = 1.;
      if(!isData){
        if(this_syst=="_PU_up"){
          purew = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),+1);
        }
        else if(this_syst=="_PU_down"){
          purew = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);
        }
        else{
          purew = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
        }
      }
      this_weight *= purew;

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

          //if( lep.at(j).GetType() == 0 ) HasStrangeLeptons = true;

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
      if(DoConversion) isSS = true;
      bool isSSForCF = isSS;
      if(RunningChargeFlipData) isSSForCF = !isSS;

      FillCutFlowByName(Suffix, "TwoLeptons", this_weight, isData);

      if(!DoConversion && !DoMCClosure && !isSSForCF) continue;

      //FIXME for v2
      if(k_sample_name.Contains("HeavyNeutrinoTo")){
        //==== SS:OS = 0.4736:0.5264
        if(isSSForCF) this_weight *= 0.5/0.4736;
        else          this_weight *= 0.5/0.5264;
      }

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

        if( lep.at(0).Pt()>=20 && lep.at(1).Pt() >= 15){
          FillHist("CutStudy_m_ll_pt2015"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
        }

        if( lep.at(0).DeltaR( lep.at(1) ) < 0.1 ){

          FillHist("CutStudy_m_ll_dR0p1_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
          if(isSS) FillHist("CutStudy_m_ll_dR0p1_SS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);
          else FillHist("CutStudy_m_ll_dR0p1_OS_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 40., 400);

        }

      }
      bool mll10GeV = ( lep.at(0)+lep.at(1) ).M() < 10.;
/*
      if(mll10GeV) continue;
      if(isSSForCF) FillCutFlowByName(Suffix, "LowDileptonMass", w_cutflow[Suffix], isData);
*/

      double this_weight_err(0.);
      vector<double> FRweights;
      vector<TString> FRweights_name;
      FRweights.clear();
      FRweights_name.clear();
      if( isLOOSE ){

        get_eventweight(muons, electrons, isT, NearBjet, 0);
        this_weight *= weight_fr;
        this_weight_err = this_weight*weight_err_fr;

        if(MCFakeSubtract){
          this_weight *= -1.;
          this_weight_err = 0.;
        }

        if(DoFakeSyst){

          ElFR_key = "";
          MuFR_key = "";
          get_eventweight(muons, electrons, isT, NearBjet, 0);
          FRweights.push_back(weight_fr);
          FRweights_name.push_back("Central");
          //cout << "Central : " << weight_fr << endl;

          FRweights.push_back(weight_fr+weight_fr*weight_err_fr);
          FRweights_name.push_back("StatUp");

          ElFR_key = "Awayjet_40";
          MuFR_key = "Awayjet_40";
          get_eventweight(muons, electrons, isT, NearBjet, 0);
          FRweights.push_back(weight_fr+weight_fr*weight_err_fr);
          FRweights_name.push_back("PromptUp");

          ElFR_key = "";
          MuFR_key = "";
          get_eventweight(muons, electrons, isT, NearBjet_TRUE, 0);
          //cout << "-> " << NearBjet_TRUE.at(0) << "\t" << NearBjet_TRUE.at(1) << " : " << weight_fr << endl;
          FRweights.push_back(weight_fr);
          FRweights_name.push_back("BJet");

          TString Muon_FRsystsources[8] = {"Awayjet_20", "Awayjet_30", "Awayjet_60", "dphi_1p5", "dphi_1", "dphi_3", "pj_over_pl_1", "pj_over_pl_2"};
          TString Electron_FRsystsources[8] = {"Awayjet_20", "Awayjet_30", "Awayjet_60", "dphi_1p5", "dphi_2", "dphi_3", "pj_over_pl_0p8", "pj_over_pl_1p2"};

          if(Suffix.Contains("DiMuon")){

            for(int aa=0; aa<8; aa++){
              ElFR_key = "";
              MuFR_key = Muon_FRsystsources[aa];
              get_eventweight(muons, electrons, isT, NearBjet, 0);
              FRweights.push_back(weight_fr);
              FRweights_name.push_back("Muon_"+Muon_FRsystsources[aa]);
            }

          }
          else if(Suffix.Contains("DiElectron")){
          
            for(int aa=0; aa<8; aa++){
              ElFR_key = Electron_FRsystsources[aa];
              MuFR_key = "";
              get_eventweight(muons, electrons, isT, NearBjet, 0);
              FRweights.push_back(weight_fr);
              FRweights_name.push_back("Electron_"+Electron_FRsystsources[aa]);
            }

          }
          else if(Suffix.Contains("EMu")){

            for(int aa=0; aa<8; aa++){
              ElFR_key = "";
              MuFR_key = Muon_FRsystsources[aa];
              get_eventweight(muons, electrons, isT, NearBjet, 0);
              FRweights.push_back(weight_fr);
              FRweights_name.push_back("Muon_"+Muon_FRsystsources[aa]);
            }

            for(int aa=0; aa<8; aa++){
              ElFR_key = Electron_FRsystsources[aa];
              MuFR_key = "";
              get_eventweight(muons, electrons, isT, NearBjet, 0);
              FRweights.push_back(weight_fr);
              FRweights_name.push_back("Electron_"+Electron_FRsystsources[aa]);
            }

          }

          //==== reset
          ElFR_key = "";
          MuFR_key = "";
        }


      }
      if( RunningChargeFlipData && !isSS ){
        if( Suffix.Contains("DiElectron") ){
          GetCFWeight(electrons_before_shift.at(0), electrons_before_shift.at(1));
        }
        if( Suffix.Contains("EMu") ){
          GetCFWeight(electrons_before_shift.at(0));
        }
        double tmp_weight = this_weight;
        this_weight     = tmp_weight*weight_cf;
        this_weight_err = tmp_weight*weight_err_cf;

      }

      //==== Now,
      //==== Fill Histogram
      //====

      //==== We have to sort lepton after we get fake weight,
      //==== because isT = {muon, electron}
      //==== If MCClosure, keep ordering as Muon-Electron (to see type)
      std::sort(lep.begin(), lep.end(), LeptonPtComparing);

      //==== Do CutFlow here..
      if(isSSForCF) FillCutFlowByName(Suffix, "SS", this_weight, isData);

      //==== No Extra lepton
      if(!isNoExtra) continue;
      if(isSSForCF) FillCutFlowByName(Suffix, "NoExtraLepton", this_weight, isData);

      //==== No Extra different flavour lepton
      if(!isNoExtraOtherFlavour) continue;
      if(isSSForCF) FillCutFlowByName(Suffix, "NoExtraFlavourLepton", this_weight, isData);

      if(mll10GeV) continue;
      if(isSSForCF) FillCutFlowByName(Suffix, "LowDileptonMass", this_weight, isData);



      std::map< TString, bool > map_Region_to_Bool;
      map_Region_to_Bool.clear();

      //==== Default Variation
      //==== 1) OnZ/OffZ/All
      //==== 2) OS/SS/All

      //==== SS-dilepton
      map_Region_to_Bool[""] = true;

/*
      //FIXME singleor
      if(Suffix.Contains("DiMuon")){
        map_Region_to_Bool["SingleButNotDouble"] = !PassTriggerOR(triggerlist_DiMuon_Mu17Mu8) && PassTriggerOR(triggerlist_DiMuon_Mu24);
      }
      else if(Suffix.Contains("DiElectron")){
        map_Region_to_Bool["SingleButNotDouble"] = !PassTriggerOR(triggerlist_DiElectron_Ele23Ele12) && PassTriggerOR(triggerlist_DiElectron_Ele27);
      }
      else if(Suffix.Contains("EMu")){
        map_Region_to_Bool["SingleButNotDouble"] = !PassTriggerOR(triggerlist_EMu_EMu) && (PassTriggerOR(triggerlist_DiMuon_Mu24) || PassTriggerOR(triggerlist_DiElectron_Ele27) );
      }

      if(this_syst=="" && map_Region_to_Bool["SingleButNotDouble"]){
        cout << "###########" << endl;
        cout << "## Event ##" << endl;
        cout << "###########" << endl << endl;

        cout << "leading :" << endl;
        cout << "  pt = " << lep.at(0).Pt() << endl;
        cout << "  eta = " << lep.at(0).Eta() << endl;
        cout << "  phi = " << lep.at(0).Phi() << endl;
        cout << "  reliso = " << lep.at(0).RelIso() << endl;
        cout << "  dXY = " << lep.at(0).dXY() << endl;
        cout << "  dZ = " << lep.at(0).dZ() << endl;
        TString lep0_trigmatch = (lep.at(0).GetMuonPtr())->TrigMatch();
        cout << "  HLT_Mu17_TrkIsoVVL_v = " << lep0_trigmatch.Contains("HLT_Mu17_TrkIsoVVL_v") << endl;
        cout << "  HLT_Mu17_v = " << lep0_trigmatch.Contains("HLT_Mu17_v") << endl;
        cout << "  HLT_Mu8_TrkIsoVVL_v = " << lep0_trigmatch.Contains("HLT_Mu8_TrkIsoVVL_v") << endl;
        cout << "  HLT_Mu8_v = " << lep0_trigmatch.Contains("HLT_Mu8_v") << endl;
        cout << "  HLT_IsoMu24_v = " << lep0_trigmatch.Contains("HLT_IsoMu24_v") << endl;


        cout << "subleading :" << endl;
        cout << "  pt = " << lep.at(1).Pt() << endl;
        cout << "  eta = " << lep.at(1).Eta() << endl;
        cout << "  phi = " << lep.at(1).Phi() << endl;
        cout << "  reliso = " << lep.at(1).RelIso() << endl;
        cout << "  dXY = " << lep.at(1).dXY() << endl;
        cout << "  dZ = " << lep.at(1).dZ() << endl;
        TString lep1_trigmatch = (lep.at(1).GetMuonPtr())->TrigMatch();
        cout << "  HLT_Mu17_TrkIsoVVL_v = " << lep1_trigmatch.Contains("HLT_Mu17_TrkIsoVVL_v") << endl;
        cout << "  HLT_Mu17_v = " << lep1_trigmatch.Contains("HLT_Mu17_v") << endl;
        cout << "  HLT_Mu8_TrkIsoVVL_v = " << lep1_trigmatch.Contains("HLT_Mu8_TrkIsoVVL_v") << endl;
        cout << "  HLT_Mu8_v = " << lep1_trigmatch.Contains("HLT_Mu8_v") << endl;
        cout << "  HLT_IsoMu24_v = " << lep1_trigmatch.Contains("HLT_IsoMu24_v") << endl;

      }
*/

      //==== # of jets
      map_Region_to_Bool["0jets"] = (jets.size()==0);
      map_Region_to_Bool["0jets_0nlbjets_dRllge2p5"] = (jets.size()==0) && (nbjets_nolepveto==0) && ( lep.at(0).DeltaR(lep.at(1)) > 2.5 );

      map_Region_to_Bool["1jets"] = (jets.size()==1);
      map_Region_to_Bool["1jets_0nlbjets"] = (jets.size()==1) && (nbjets_nolepveto==0);
      map_Region_to_Bool["1jets_0nlbjets_mllge110"] = (jets.size()==1) && (nbjets_nolepveto==0) && (( lep.at(0)+lep.at(1) ).M() >= 110.);

      //==== # of b jet
      map_Region_to_Bool["0nlbjets"] = (nbjets_nolepveto==0);
      map_Region_to_Bool["Inclusive1nlbjets"] = (nbjets_nolepveto>=1);

      if(k_sample_name.Contains("HN")){
        map_Region_to_Bool["LegacyTwoJets"] = (jets_nonfatjetveto.size() >= 2);
      }


      //==== Z-seleciton
      map_Region_to_Bool["Z_CR"] = (!isSSForCF) && ((lep.at(0)+lep.at(1)).M()>60.) && ((lep.at(0)+lep.at(1)).M()<12.);

      //==== W+W+ CR
      if(jets_eta5.size() >= 2){
        snu::KJet j1 = jets_eta5.at(0);
        snu::KJet j2 = jets_eta5.at(1);
        double dEtajj = fabs(j1.Eta()-j2.Eta());
        double MeanEtajj = (j1.Eta()+j2.Eta())/2.;

        map_Region_to_Bool["WpWp_CR"] = (j1.Pt() > 30.) && (j2.Pt() > 30.)
                                                && (( lep.at(0)+lep.at(1) ).M() >= 20.) && (MET>40.)
                                                && (nbjets_nolepveto==0)
                                                && ((j1+j2).M() > 500.)
                                                && (dEtajj>2.5)
                                                && ( (lep.at(0).Eta()-MeanEtajj)/dEtajj < 0.75 )
                                                && ( (lep.at(1).Eta()-MeanEtajj)/dEtajj < 0.75 );
        if(Suffix.Contains("DiElectron")){
          map_Region_to_Bool["WpWp_CR"] = map_Region_to_Bool["WpWp_CR"] && isOffZ;
        }

      }


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
      }
      if( fatjets.size()>=1 ){
        double mfatjet = GetFatJetMassClosest(fatjets, 80.4, index_fjW);
      }

      //==== Preselection
      bool TwoJet_NoFatJet = (jets.size()>=2) && (fatjets.size()==0);
      bool OneJet_NoFatJet = (jets.size()==1) && (fatjets.size()==0) && ( (lep.at(0)+lep.at(1)).M() < 80 ); // has m(ll) < 80 GeV
      bool OneFatJet       =                     (fatjets.size()>=1);

      map_Region_to_Bool["Preselection"] = TwoJet_NoFatJet || OneJet_NoFatJet || OneFatJet;
      map_Region_to_Bool["Preselection_secondptge20"] = map_Region_to_Bool["Preselection"] && (lep.at(1).Pt() >= 20.);
      //==== For DiElectron, remove Z peak (CF)
      if(Suffix.Contains("DiElectron")){
        map_Region_to_Bool["Preselection"] = map_Region_to_Bool["Preselection"] && isOffZ;
      }

      if(isSSForCF){
        if(TwoJet_NoFatJet || OneJet_NoFatJet || OneFatJet){
          FillCutFlowByName(Suffix, "JetRequirements", this_weight, isData);
          if(isOffZ){
            FillCutFlowByName(Suffix, "OffZ", this_weight, isData);
          }
        }
      }

      //==== If Preselection, then define Low/High

      if(map_Region_to_Bool["Preselection"]){

        //==== Low Mass
        //==== 1) TwoJet_NoFatJet
        //==== 2) OneJet_NoFatjet (has mll<80GeV cut)

        double mlljj_low = -999.;
        if(TwoJet_NoFatJet) mlljj_low = (lep.at(0)+lep.at(1)+jets.at(index_lljjW_j1)+jets.at(index_lljjW_j2)).M();
        else if(OneJet_NoFatJet) mlljj_low = (lep.at(0)+lep.at(1)+jets.at(0)).M();

        map_Region_to_Bool["Low"] = (TwoJet_NoFatJet || OneJet_NoFatJet) &&
                                            (nbjets_nolepveto == 0) &&
                                            ( mlljj_low < 300.) &&
                                            (MET < 80.);
        map_Region_to_Bool["LowCR"] = (TwoJet_NoFatJet || OneJet_NoFatJet) &&
                                              ( (nbjets_nolepveto >= 1) || (MET > 100.) ) &&
                                              ( mlljj_low < 300.);

        map_Region_to_Bool["Low_TwoJet_NoFatJet"] = TwoJet_NoFatJet &&
                                                            (nbjets_nolepveto == 0) &&
                                                            ( mlljj_low < 300.) &&
                                                            (MET < 80.);
        map_Region_to_Bool["LowCR_TwoJet_NoFatJet"] = TwoJet_NoFatJet &&
                                                              ( (nbjets_nolepveto >= 1) || (MET > 100.) ) &&
                                                              ( mlljj_low < 300.);

        map_Region_to_Bool["Low_OneJet_NoFatJet"] = OneJet_NoFatJet &&
                                                   (nbjets_nolepveto == 0) &&
                                                   ( mlljj_low < 300.) &&
                                                   (MET < 80.);
        map_Region_to_Bool["LowCR_OneJet_NoFatJet"] = OneJet_NoFatJet &&
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

        map_Region_to_Bool["High"] = ( TwoJet_NoFatJet || OneFatJet ) &&
                                             (nbjets_nolepveto == 0) &&
                                             ( mjj_high < 150. ) &&
                                             ( MET*MET/ST < 15. );
        map_Region_to_Bool["HighCR"] = ( TwoJet_NoFatJet || OneFatJet ) &&
                                               ( (nbjets_nolepveto >= 1) || (MET*MET/ST > 20.) ) &&
                                               ( mjj_high < 150. );

        map_Region_to_Bool["High_TwoJet_NoFatJet"]   = map_Region_to_Bool["High"] && TwoJet_NoFatJet;
        map_Region_to_Bool["HighCR_TwoJet_NoFatJet"] = map_Region_to_Bool["HighCR"] && TwoJet_NoFatJet;

        map_Region_to_Bool["High_OneFatJet"]   = map_Region_to_Bool["High"] && OneFatJet;
        map_Region_to_Bool["HighCR_OneFatJet"] = map_Region_to_Bool["HighCR"] && OneFatJet;


        if(Suffix.Contains("EMu")){
          if(lep.at(1).LeptonFlavour()==KLepton::ELECTRON){
            map_Region_to_Bool["Preselection_ElectronSubLead"] = true;
          }
          else{
            map_Region_to_Bool["Preselection_MuonSubLead"] = true;
          }
        }
      }

      //==== Test 1) Loose sample faking jet pt distribution
      if(LooseSampleFakeJetPt){

        if(map_Region_to_Bool["Preselection"] && jets_nolepveto.size() > 0){

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
      if(DoFakeSyst){

        if(isSSForCF){

          TString ToSaveRegions[] = {
            "Preselection",
            "Low", "High",
          };
          for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
            TString this_suffix = it->first;

            bool ToSave = false;
            for(int aaa=0;aaa<3;aaa++){
              if(this_suffix.Contains(ToSaveRegions[aaa])) ToSave = true;
            }

            if(!ToSave) it->second = false;

          }

          for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
            TString this_suffix = it->first;
            if(it->second){

              for(unsigned int j=0; j<FRweights.size(); j++){
                FillHist(this_suffix+"_SS_FRsyst_"+FRweights_name.at(j), 0., FRweights.at(j), 0., 1., 1);
              }

            }
          }

        } // If SS

        continue;
      }

      //==== Test 3) If MCClosure, only save Preselection
      if(DoMCClosure){

        for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
          TString this_suffix = it->first;

          bool ToSave = false;
          if(this_suffix.Contains("Preselection")) ToSave = true;
          if(this_suffix.Contains("Low")) ToSave = true;
          if(this_suffix.Contains("High")) ToSave = true;
          if(this_suffix.Contains("CR")) ToSave = false;

          if(!ToSave) it->second = false;
        }

        bool PositiveMCWeight = std::find(k_flags.begin(), k_flags.end(), "PositiveMCWeight") != k_flags.end();
        if(PositiveMCWeight) this_weight *= MCweight;

        if( map_Region_to_Bool["Preselection"] && isSSForCF ){
          FillHist("NLooseNotTight_"+Suffix+"_Preselection_SS", n_triLoose_leptons-n_triTight_leptons, this_weight, 0., 3., 3);
          FillHist("NLooseNotTight_weight1_"+Suffix+"_Preselection_SS", n_triLoose_leptons-n_triTight_leptons, 1., 0., 3., 3);
        }
      }

      //==== Make Ntuple
      bool NtupleSkim = (map_Region_to_Bool["Preselection"]) && isSSForCF;

      if(RunNtp && NtupleSkim){
        int temp_n_ntuplevar = 51;
        if(isData) temp_n_ntuplevar = 53;
        const int n_ntuplevar = temp_n_ntuplevar;

        double cutop[n_ntuplevar];
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

        if(isData){
          cutop[51] = N_RunNumber;
          cutop[52] = N_EventNumber;
        }

        FillNtp("Ntp_"+Suffix+this_syst+"_Preselection_SS",cutop);

      }
      if(RunNtp) continue;

      //==== Add Suffix
      std::map< TString, bool > ForPlot_map_Region_to_Bool;
      ForPlot_map_Region_to_Bool.clear();
      for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
        
        TString this_suffix = it->first;
        if(this_suffix==""){
          ForPlot_map_Region_to_Bool[Suffix] = it->second;
          ForPlot_map_Region_to_Bool["DiLepton"] = it->second;
        }
        else{
          ForPlot_map_Region_to_Bool[Suffix+"_"+this_suffix] = it->second;
          ForPlot_map_Region_to_Bool["DiLepton_"+this_suffix] = it->second;
        }

      }


      for(std::map< TString, bool >::iterator it = ForPlot_map_Region_to_Bool.begin(); it != ForPlot_map_Region_to_Bool.end(); it++){
        TString this_suffix = it->first;
        //cout << this_suffix << "\t" << it->second << endl;
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
            }
            else{
              FillDiLeptonPlot(this_suffix+"_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }

            //==== If Dielectron channel, fill OffZ and OnZ too
            if(Suffix.Contains("DiElectron")){

              if(isOffZ){
                //==== OffZ / SS
                if(isSS){
                  FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
                else{
                  FillDiLeptonPlot(this_suffix+"_OffZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
              }
              else{
                //==== OnZ / SS
                if(isSS){
                  FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
                else{
                  FillDiLeptonPlot(this_suffix+"_OnZ_OS", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
                }
              }

            }

            //==== For WJet MC Closure, fill OS and SS together
            if(DoMCClosure && k_sample_name.Contains("WJets")){
              FillDiLeptonPlot(this_suffix+"_AllCharge", lep, jets, jets_fwd, jets_nolepveto, fatjets, this_weight, this_weight_err);
            }


          }

          //=================
          //==== ChargeFlip
          //=================

          //==== To NOT DRAWING OS plot from CF (file size..)
          //==== Separate CF from others

          else{
            //==== using OS event, weight CF and estimate SS
            if( (Suffix.Contains("DiElectron")||Suffix.Contains("EMu")) && !isSS){

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
  


void PairNAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void PairNAnalyzer::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  return;
  
}

PairNAnalyzer::~PairNAnalyzer() {
  
  Message("In PairNAnalyzer Destructor" , INFO);
  
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


void PairNAnalyzer::FillCutFlow(TString cut, float w){

  
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

void PairNAnalyzer::FillCutFlowByName(TString histname, TString cut, float w, bool IsDATA){

  if(AUTO_syst_type!="") return;

  TString this_histname = "Cutflow_"+histname;

  FillHist(this_histname+"_"+cut, 0., w, 0., 1., 1);

}


void PairNAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void PairNAnalyzer::MakeHistograms(){
  //// Additional plots to make

  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);

  int N_sys = 2*14+1;
  TString systs[] = {"", "_MuonEn_up", "_MuonEn_down", "_JetEn_up", "_JetEn_down", "_JetRes_up", "_JetRes_down", "_Unclustered_up", "_Unclustered_down", "_MuonIDSF_up", "_MuonIDSF_down", "_PU_down", "_PU_up", "_TriggerSF_down", "_TriggerSF_up", "_ElectronIDSF_up", "_ElectronIDSF_down", "_ElectronEn_up", "_ElectronEn_down", "_BTagSFEff_up", "_BTagSFEff_down", "_BTagSFMiss_up", "_BTagSFMiss_down", "_JetMass_up", "_JetMass_down", "_JetMassRes_up", "_JetMassRes_down", "_Tau21_up", "_Tau21_down"};

  cout << "[JSKIM] isData = " << isData << endl;
  cout << "[JSKIM] k_isdata = " << k_isdata << endl;
  for(int i=0; i<N_sys; i++){

    if(isData){

  MakeNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets:N_RunNumber:N_EventNumber");

  MakeNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets:N_RunNumber:N_EventNumber");

  MakeNtp("Ntp_EMu"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets:N_RunNumber:N_EventNumber");

    }
    else{

  MakeNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");

  MakeNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");

  MakeNtp("Ntp_EMu"+systs[i]+"_Preselection_SS", "leadingLepton_Pt:secondLepton_Pt:DeltaRl1l2:m_ll:isSS:isOffZ:Njets:Nbjets:Njets_nolepveto:Nbjets_nolepveto:Nfwdjets:Nbfwdjets:leadingJet_Pt:secondJet_Pt:DeltaRjjptorder:m_jjptorder:m_Leadljjptorder:m_SubLeadljjptorder:m_lljjptorder:leadingJet_jjWclosest_pt:secondJet_jjWclosest_pt:m_jj_jjWclosest:m_Leadljj_jjWclosest:m_SubLeadljj_jjWclosest:m_lljj_jjWclosest:DeltaRjjWclosest:DeltaRLeadl_jjWclosest:DeltaRSubLeadl_jjWclosest:DeltaRLeadl_SubLeadljjWclosest:DeltaRSubLeadl_LeadljjWclosest:leadingJet_lljjWclosest_pt:secondJet_lljjWclosest_pt:m_jj_lljjWclosest:m_Leadljj_lljjWclosest:m_SubLeadljj_lljjWclosest:m_lljj_lljjWclosest:DeltaRlljjWclosest:DeltaRLeadl_lljjWclosest:DeltaRSubLeadl_lljjWclosest:DeltaRLeadl_SubLeadllljjWclosest:DeltaRSubLeadl_LeadllljjWclosest:fwd_dRjj:PFMET:ST:HT:LT:weight:weight_err:leadingLepton_Eta:secondLepton_Eta:Nfatjets");

    }

    GetNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_DiMuon"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);
    GetNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_DiElectron"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);
    GetNtp("Ntp_EMu"+systs[i]+"_Preselection_SS")->Branch("PdfWeights", "vector<float>",&ForTree_PdfWeights);
    GetNtp("Ntp_EMu"+systs[i]+"_Preselection_SS")->Branch("ScaleWeights", "vector<float>",&ForTree_ScaleWeights);

  }


  /**
   *  Remove//Overide this PairNAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/

}


void PairNAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //


}

void PairNAnalyzer::FillDiLeptonPlot(
  TString histsuffix,
  std::vector< KLepton> leptons,
  std::vector< snu::KJet > jets,
  std::vector< snu::KJet > jets_fwd,
  std::vector< snu::KJet > jets_nolepveto,
  std::vector< snu::KFatJet > fatjets,
  double thisweight,
  double thieweighterr
  ){

  JSFillHist("Systematics", "Yields_"+histsuffix, 1.*AUTO_it_syst, thisweight, 0., 1.*AUTO_N_syst, AUTO_N_syst);

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

  double MinDeltaR = 999.;
  for(unsigned int i=0; i<leptons.size(); i++){
    for(unsigned int j=0; j<jets.size(); j++){
      double this_dr = (leptons.at(i)).DeltaR( jets.at(j) );
      if(this_dr<MinDeltaR) MinDeltaR = this_dr;
    }
  }
  JSFillHist(histsuffix, "MinDeltaR_lj_"+histsuffix, MinDeltaR, thisweight, 0., 10., 100);

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

      TString EtaRegion = "InnerBarrel";
      if(fabs(leptons.at(i).GetElectronPtr()->SCEta()) > 1.479) EtaRegion = "EndCap";
      else if(fabs(leptons.at(i).GetElectronPtr()->SCEta()) > 0.8) EtaRegion = "OuterBarrel";
      else EtaRegion = "InnerBarrel";

      JSFillHist(histsuffix, EtaRegion+"Lepton_Pt_"+histsuffix,  leptons.at(i).Pt(), thisweight, 0., 2000., 2000);
      JSFillHist(histsuffix, EtaRegion+"Lepton_Eta_"+histsuffix, leptons.at(i).Eta(), thisweight, -3., 3., 60);

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
    JSFillHist(histsuffix, "m_Leadlfj_ptorder_"+histsuffix, (leptons.at(0)+fatjets.at(0)).M(), thisweight, 0., 4000., 4000);
    JSFillHist(histsuffix, "m_SubLeadlfj_ptorder_"+histsuffix, (leptons.at(1)+fatjets.at(0)).M(), thisweight, 0., 4000., 4000);
    JSFillHist(histsuffix, "m_llfj_ptorder_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjets.at(0)).M(), thisweight, 0., 4000., 4000);

    //==== fjWclosest
    JSFillHist(histsuffix, "m_Leadlfj_fjWclosest_"+histsuffix, (leptons.at(0)+fatjets.at(index_fjW)).M(), thisweight, 0., 4000., 4000);
    JSFillHist(histsuffix, "m_SubLeadlfj_fjWclosest_"+histsuffix, (leptons.at(1)+fatjets.at(index_fjW)).M(), thisweight, 0., 4000., 4000);
    JSFillHist(histsuffix, "m_llfj_fjWclosest_"+histsuffix, (leptons.at(0)+leptons.at(1)+fatjets.at(index_fjW)).M(), thisweight, 0., 4000., 4000);
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




void PairNAnalyzer::GetCFWeight(KLepton lep1, KLepton lep2){

  if(lep1.Charge()==lep2.Charge()) return;

  //==== Okay, now lep1 and lep2 are OS

  double cf1 = GetCF(lep1, false);
  double cf2 = GetCF(lep2, false);
  weight_cf = cf1/(1.-cf1) + cf2/(1.-cf2);

  double cf1_scale_up = 1.;
  if(1./lep1.Pt() < 0.013) cf1_scale_up = 2.;
  else cf1_scale_up = 1.2;

  double cf2_scale_up = 1.;
  if(1./lep2.Pt() < 0.013) cf2_scale_up = 2.;
  else cf2_scale_up = 1.2;

  cf1 *= cf1_scale_up;
  cf2 *= cf2_scale_up;
  weight_err_cf = cf1/(1.-cf1) + cf2/(1.-cf2);
  weight_err_cf = weight_err_cf-weight_cf;


}

void PairNAnalyzer::GetCFWeight(KLepton lep1){

  //==== Okay, now lep1 and lep2 are OS

  double cf1 = GetCF(lep1, false);
  double cf1_err = GetCF(lep1, true);

  weight_cf = cf1/(1.-cf1);
  weight_err_cf = sqrt( cf1_err/( (1.-cf1)*(1.-cf1) ) );
  

}

double PairNAnalyzer::GetCF(KLepton lep, bool geterr){

  double el_eta = fabs(lep.GetElectronPtr()->SCEta());
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./lep.Pt();
  double a = 999., b= 999.;
  if(el_eta < 0.8){
    if(invPt< 0.023){a=(-0.001008); b=(3.055e-05);}
    else{a=(-0.0002274); b=(1.208e-05);}
  }
  else if(el_eta < 1.4442){
    if(invPt< 0.015){a=(-0.03484); b=(0.0007158);}
    else if(invPt< 0.023){a=(-0.01439); b=(0.0004098);}
    else{a=(-0.002442); b=(0.0001343);}
  }
  else{
    if(invPt< 0.012){a=(-0.4065); b=(0.006359);}
    else if(invPt< 0.021){a=(-0.1111); b=(0.002914);}
    else{a=(-0.02092); b=(0.001051);}
  }

  double sf(1.);

  if(el_eta < 1.4442) sf = 0.8014;
  else sf = 0.8658;

  //if(el_eta < 1.4442) sf = 0.7146;
  //else sf = 0.8555;

  //==== Below using MC (CF)/(CF+type40)
  //if(el_eta < 1.4442) sf = 0.7052*0.8251;
  //else sf = 0.8548*0.8673;



/*
  if(el_eta < 0.8){
    if(invPt< 0.023){a=(-0.001196); b=(3.806e-05);}
    else{a=(0.0008751); b=(-9.844e-06);}
  }
  else if(el_eta < 1.4442){
    if(invPt< 0.015){a=(-0.03537); b=(0.0007231);}
    else if(invPt< 0.023){a=(-0.01381); b=(0.0004019);}
    else{a=(-0.0007848); b=(0.0000972);}
  }
  else{
    if(invPt< 0.011){a=(-0.5366); b=(0.007573);}
    else if(invPt< 0.020){a=(-0.09821); b=(0.002702);}
    else{a=(-0.02011); b=(0.00105);}
  }
  double sf(1.);
  if(el_eta < 1.4442) sf = 0.7112*0.8270;
  else sf = 0.850198*0.8673;
*/

  double sys=0.;
  double rate = (a)*invPt + (b);
  if(rate < 0) rate = 0.;

  rate *= sf;

  if(!geterr) return rate;
  else return 0.;

}

double PairNAnalyzer::GetDijetMassClosest(std::vector<snu::KJet> js, double mass, int& m, int& n){

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

double PairNAnalyzer::GetDileptonDijetMassClosest(std::vector<KLepton> leps, std::vector<snu::KJet> js, double mass, int& m, int& n){

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

double PairNAnalyzer::GetFatJetMassClosest(std::vector<snu::KFatJet> fjs, double mass, int& m){

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

double PairNAnalyzer::CorrPt(KLepton lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso()-T_iso)));
  return ptcorr;
}

double PairNAnalyzer::CorrPt(snu::KMuon lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso04()-T_iso)));
  return ptcorr;
}

double PairNAnalyzer::CorrPt(snu::KElectron lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.PFRelIso(0.3)-T_iso)));
  return ptcorr;
}



double PairNAnalyzer::GetMuonFR(bool geterr, float pt,  float eta, int NearBjet){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.4) eta = 2.3;

  //cout << "[PairNAnalyzer::GetMuonFR] pt = " << pt << endl;
  //cout << "[PairNAnalyzer::GetMuonFR] eta = " << eta << endl;
  //cout << "[PairNAnalyzer::GetMuonFR] MuFR_key = " << MuFR_key << endl;
  //cout << "[PairNAnalyzer::GetMuonFR] NearBjet = " << NearBjet << endl;

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
  //cout << "[PairNAnalyzer::GetMuonFR] => FR = " << THISFRHIST->GetBinContent(binx) << endl;
  

  if(geterr) return THISFRHIST->GetBinError(binx);
  else return THISFRHIST->GetBinContent(binx);

}

double PairNAnalyzer::GetMuonPR(bool geterr, float pt,  float eta){
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

double PairNAnalyzer::GetElectronFR(bool geterr, float pt,  float eta, int NearBjet){

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

double PairNAnalyzer::GetElectronPR(bool geterr, float pt,  float eta){
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


void PairNAnalyzer::get_eventweight(std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<bool> isT, std::vector<int> NearBjet, int HalfSampleErrorDir){

  unsigned int n_leptons = isT.size();
  //cout << "[PairNAnalyzer::get_eventweight] muons.size() = " << muons.size() << ", electrons.size() = " << electrons.size() << endl;

  vector<float> lep_pt, lep_eta;
  vector<bool> ismuon;
  for(unsigned int i=0; i<muons.size(); i++){
    //lep_pt.push_back( CorrPt(muons.at(i), 0.07) );
    lep_pt.push_back( muons.at(i).Pt() ); // now, fake lepton pt is replaced by pt-cone
    lep_eta.push_back(muons.at(i).Eta());
    ismuon.push_back(true);
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    //lep_pt.push_back( CorrPt(electrons.at(i), 0.08) );
    lep_pt.push_back( electrons.at(i).Pt() ); // now, fake lepton pt is replaced by pt-cone
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
      //cout << "[PairNAnalyzer::get_eventweight] "<<i<<" th lepton is Loose" << endl;
      fr_onlyLoose.push_back( a.at(i) );
    }
  }

  //==== Initialise weight
  float this_weight=-1.;

  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }
  //cout << "[PairNAnalyzer::get_eventweight] this_weight = " << this_weight << endl;

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

bool PairNAnalyzer::JSFatJetID(snu::KFatJet fatjet){

  if( !(fatjet.PrunedMass() > 40) ) return false;
  if( !(fatjet.PrunedMass() < 130) ) return false;
  if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;

  return true;

}

bool PairNAnalyzer::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 0.8 ) return false;
  }

  return true;


}



















