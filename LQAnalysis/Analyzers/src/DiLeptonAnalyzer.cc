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

DiLeptonAnalyzer::DiLeptonAnalyzer() :  AnalyzerCore(), weight_cf(-999), weight_err_cf(-999), weight_fr(-999), weight_err_fr(-999), MET(-999), METphi(-999), nbjets(-999), nbjets_fwd(-999), nbjets_nolepveto(-999), n_vtx(-999)
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

  string lqdir = getenv("LQANALYZER_DIR");

  //==== Get Fake

  TFile *file_Muon_FR = new TFile( (lqdir+"/data/Fake/80X/FakeRate13TeV_muon_2016_opt_all.root").c_str(), "");
  TFile *file_Electron_FR = new TFile( (lqdir+"/data/Fake/80X/FakeRate13TeV_2016_hnid.root").c_str(), "");

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
  hist_Muon_FR = (TH2D*)file_Muon_FR->Get("FakeRate_40_ptcorr_etaSNUTightisodijet_0.07_0.005_3_0.04")->Clone();

  file_Muon_FR->Close();
  file_Electron_FR->Close();

  delete file_Muon_FR;
  delete file_Electron_FR;

  origDir->cd();

  return;
}


void DiLeptonAnalyzer::ExecuteEvents()throw( LQError ){

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

  double muon_id_iso_sf = 1.; //FIXME
  //double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", muons, 0); //FIXME

  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muons);

  //==== Electrons
  std::vector< snu::KElectron > electrons = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  double electron_sf = 1.; //FIXME
  //double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrons, 0); //FIXME
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);

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
  std::vector<snu::KElectron> electrons_veto = GetElectrons(false, false, "ELECTRON_HN_VETO");
  std::vector<bool> isT;

  int n_veto_muons = muons_veto.size();
  int n_triLoose_muons = muons.size();
  int n_triTight_muons(0);
  for(unsigned int i=0; i<muons.size(); i++){
    if(PassID(muons.at(i), "MUON_HN_TIGHT")){
      isT.push_back(true);
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
      n_triTight_electrons++;
    }
    else{
      isT.push_back(false);
    }
  }

  int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
  int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

  bool isTwoMuon_TT    = (n_triLoose_muons == 2 && n_triTight_muons == 2); 
  bool isTwoMuon_Loose = (n_triLoose_muons == 2 && n_triTight_muons != 2);

  bool isTwoElectron_TT    = (n_triLoose_electrons == 2 && n_triTight_electrons == 2); 
  bool isTwoElectron_Loose = (n_triLoose_electrons == 2 && n_triTight_electrons != 2);

  bool isEMu_TT    = (n_triLoose_muons == 1     && n_triTight_muons == 1) &&
                     (n_triLoose_electrons == 1 && n_triTight_electrons == 1);
  bool isEMu_Loose = (n_triLoose_muons == 1)     &&
                     (n_triLoose_electrons == 1) &&
                     (n_triTight_leptons != 2);

/*
  bool isTwoMuon_TT = (n_triLoose_leptons == 2)
                      && (n_triLoose_muons == 2 && n_triTight_muons == 2);
  bool isTwoMuon_Loose = (n_triLoose_leptons == 2)
                         && (n_triLoose_muons == 2 && n_triTight_muons != 2);

  bool isTwoElectron_TT = (n_triLoose_leptons == 2)
                          && (n_triLoose_electrons == 2 && n_triTight_electrons == 2);
  bool isTwoElectron_Loose = (n_triLoose_leptons == 2)
                             && (n_triLoose_electrons == 2 && n_triTight_electrons != 2);

  bool isEMu_TT = (n_triLoose_leptons == 2)
                  && (n_triLoose_muons == 1 && n_triTight_muons == 1)
                  && (n_triLoose_electrons == 1 && n_triTight_electrons == 1);
  bool isEMu_Loose = (n_triLoose_leptons == 2)
                     && (n_triLoose_muons == 1)
                     && (n_triLoose_electrons == 1)
                     && (n_triTight_leptons != 2);
*/
  //if(n_triLoose_leptons < 2) return;
  //if(jets.size() <2) return;

  //bool NonPromptRun = k_running_nonprompt;
  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();

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

  n_vtx = Evt.nVertices();

  //==== Define Analysis Region

  std::vector< TString > Suffixs;
  std::vector< std::vector<TString> > Triggers;
  std::vector< bool > isTTs, isLOOSEs, isNoExtra, isNoExtraOtherFlavour;

  bool RunningNonPromptData = NonPromptRun && isData;
  bool RunningChargeFlipData = k_running_chargeflip && isData;

  Suffixs.push_back("DiMuon");
  Triggers.push_back(triggerlist_DiMuon);
  isTTs.push_back( isTwoMuon_TT && (!RunningNonPromptData||RunningChargeFlipData) );
  isLOOSEs.push_back( isTwoMuon_Loose && (RunningNonPromptData||RunningChargeFlipData) );
  isNoExtra.push_back( n_veto_muons == 2 );
  isNoExtraOtherFlavour.push_back( n_veto_electrons == 0 );

  Suffixs.push_back("DiElectron");
  Triggers.push_back(triggerlist_DiElectron);
  isTTs.push_back( isTwoElectron_TT && (!RunningNonPromptData||RunningChargeFlipData) );
  isLOOSEs.push_back( isTwoElectron_Loose && (RunningNonPromptData||RunningChargeFlipData) );
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
      if(muons.at(0).Pt() < 20. || muons.at(1).Pt() < 10.) continue;
    }
    if(Suffix=="DiElectron"){
/*
      bool PassDiElTrig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      bool PassSingleElTrig = PassTrigger("HLT_Ele27_WPTight_Gsf_v");
      if(PassDiElTrig && !PassSingleElTrig){
        if(electrons.at(0).Pt() < 25. || electrons.at(1).Pt() < 15.) continue;
      }
      else if(!PassDiElTrig && PassSingleElTrig){
        if(electrons.at(0).Pt() < 30.) continue;
      }
      else if(PassDiElTrig && PassSingleElTrig){
        if(electrons.at(0).Pt() < 25.) continue;
      }
*/
      if(electrons.at(0).Pt() < 25. || electrons.at(1).Pt() < 15.) continue;
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
/*
    if(!isData && Suffix=="DiMuon"){
      double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, muons, 0, 0, 0);
      double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, muons, 0, 1, 0);
      trigger_sf = trigger_eff_Data/trigger_eff_MC;
    } //FIXME NOT YET
*/
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
    double m_Z = 91.1876;
    bool isOffZ = fabs( (lep.at(0)+lep.at(1)).M() - m_Z ) > 10.;

    //==== mll Cut Study
    FillHist("CutStudy_m_ll_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 20., 200);
    bool mll4GeV = ( lep.at(0)+lep.at(1) ).M() < 12.;

    //if(!isSS && !isOffZ) continue;
    if(mll4GeV) continue;

    double this_weight_err(0.);
    if( isLOOSEs.at(i) ){
      //this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      //this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      //m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");

      get_eventweight(muons, electrons, isT, 0);

      this_weight = weight_fr;
      this_weight_err = weight_err_fr;
      
    }

    //==== Fill Histogram

    std::map< TString, bool > map_Region_to_Bool;
    map_Region_to_Bool.clear();
    map_Region_to_Bool[Suffix] = true;
    map_Region_to_Bool[Suffix+"_0jets"] = (jets.size()==0);
    map_Region_to_Bool[Suffix+"_1jets"] = (jets.size()==1);
    map_Region_to_Bool[Suffix+"_2jets"] = (jets.size()==2);
    map_Region_to_Bool[Suffix+"_3jets"] = (jets.size()==3);
    map_Region_to_Bool[Suffix+"_4jets"] = (jets.size()==4);
    map_Region_to_Bool[Suffix+"_Inclusive2jets"] = (jets.size()>=2);
    map_Region_to_Bool[Suffix+"_0bjets"] = (nbjets==0);
    map_Region_to_Bool[Suffix+"_Inclusive1bjets"] = (nbjets>=1);
    map_Region_to_Bool[Suffix+"_Inclusive2bjets"] = (nbjets>=2);
    map_Region_to_Bool[Suffix+"_Inclusive3bjets"] = (nbjets>=3);

    map_Region_to_Bool[Suffix+"_Preselection"] = false;
    //==== More CutFlows
    if( map_Region_to_Bool[Suffix+"_Inclusive2jets"] ){
      FillCutFlowByName(Suffix, "InclusiveTwoJets", w_cutflow[Suffix], isData);
      if( nbjets == 0 ){
        FillCutFlowByName(Suffix, "NoBJet", w_cutflow[Suffix], isData);
        if( nbjets_nolepveto == 0 ){
          FillCutFlowByName(Suffix, "NoBJet_nolepveto", w_cutflow[Suffix], isData);
          if( MET < 50. ){
            FillCutFlowByName(Suffix, "MET50", w_cutflow[Suffix], isData);
            if( GetDijetMassClosest(jets, 80.4) < 200. ){
              FillCutFlowByName(Suffix, "Mjj200", w_cutflow[Suffix], isData);
              map_Region_to_Bool[Suffix+"_Preselection"] = true;
            }
          }
        }
      }
    }

    for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
      TString this_suffix = it->first;
      if(it->second){

        if(!RunningChargeFlipData){
     
          //==== SS+OS
          FillDiLeptonPlot(this_suffix+"_AllCharge", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
          //==== SS
          if(isSS){
            FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
            FillHist("LeptonType_"+this_suffix+"_SS", lep.at(0).GetType(), 1., 0., 50., 50);
            FillHist("LeptonType_"+this_suffix+"_SS", lep.at(1).GetType(), 1., 0., 50., 50);

/*
            //==== DY SS event check
            bool OneOfThemType40 = false;
            if(Suffix=="DiElectron"){
              TruthPrintOut();
              cout << "## Reco ##" << endl;
              cout << "pt\teta\tphi\tcharge\ttype\tMCTruthIndex" << endl;
              cout << lep.at(0).Pt() << "\t" << lep.at(0).Eta() << "\t" << lep.at(0).Phi() << "\t" << lep.at(0).Charge() << "\t" << lep.at(0).GetType() << "\t" << lep.at(0).GetElectronPtr()->MCTruthIndex() << endl;
              cout << lep.at(1).Pt() << "\t" << lep.at(1).Eta() << "\t" << lep.at(1).Phi() << "\t" << lep.at(1).Charge() << "\t" << lep.at(1).GetType() << "\t" << lep.at(1).GetElectronPtr()->MCTruthIndex() << endl;
              if(lep.at(0).GetType()==40 || lep.at(1).GetType()==40) OneOfThemType40 = true;
            }
            FillHist("OneOfThemType40", OneOfThemType40, 1., 0., 2., 2);
*/
          }
          //==== OS
          else{
            FillDiLeptonPlot(this_suffix+"_OS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
          }

          //==== OffZ
          if(isOffZ){

            //==== SS+OS
            FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);

            //==== SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
            }
            //==== OS
            else{
              FillDiLeptonPlot(this_suffix+"_OffZ_OS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
            }

          } // Off-Z
          else{

            //==== SS+OS
            FillDiLeptonPlot(this_suffix+"_OnZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);

            //==== SS
            if(isSS){
              FillDiLeptonPlot(this_suffix+"_OnZ_SS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
            }
            //==== OS
            else{
              FillDiLeptonPlot(this_suffix+"_OnZ_OS", lep, jets, jets_fwd, jets_nolepveto, this_weight, this_weight_err);
            }

          }


        }
        //==== ChargeFlip
        else{
          //==== using OS event, weight CF and estimate SS
          if(Suffix=="DiElectron" && !isSS){
            GetCFWeight(lep.at(0), lep.at(1));
            //double weight_cf, weight_err_cf;
            std::vector<KLepton> shifted_lep;
            for(unsigned int i=0; i<lep.size(); i++){
              KLepton tmp_lep = lep.at(i);
              double shift_ = 1.-0.009; //FIXME
              tmp_lep.SetPxPyPzE(shift_*tmp_lep.Px(), shift_*tmp_lep.Py(), shift_*tmp_lep.Pz(), shift_*tmp_lep.E());
              shifted_lep.push_back(tmp_lep);
            }
            if( isLOOSEs.at(i) ){
              weight_cf     *= -1.;
              weight_err_cf *= -1.;
            }
            FillDiLeptonPlot(this_suffix+"_SS", shifted_lep, jets, jets_fwd, jets_nolepveto, this_weight*weight_cf, this_weight*weight_err_cf);
            if(isOffZ){
              FillDiLeptonPlot(this_suffix+"_OffZ_SS", shifted_lep, jets, jets_fwd, jets_nolepveto, this_weight*weight_cf, this_weight*weight_err_cf);
            }
            else{
              FillDiLeptonPlot(this_suffix+"_OnZ_SS", shifted_lep, jets, jets_fwd, jets_nolepveto, this_weight*weight_cf, this_weight*weight_err_cf);
            }
            
          }
        }

      } // passing this region
    } // Region loop


  } // Suffix (Channel) loop



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


void DiLeptonAnalyzer::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
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

  histname = "Cutflow_"+histname;

  if(GetHist(histname)){
    GetHist(histname)->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms(histname, 12,0.,12.);

    GetHist(histname)->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist(histname)->GetXaxis()->SetBinLabel(2,"MET_PV_Trig");
    GetHist(histname)->GetXaxis()->SetBinLabel(3,"TwoLeptons");
    GetHist(histname)->GetXaxis()->SetBinLabel(4,"NoExtraLepton");
    GetHist(histname)->GetXaxis()->SetBinLabel(5,"NoExtraFlavourLepton");
    GetHist(histname)->GetXaxis()->SetBinLabel(6,"InclusiveTwoJets");
    GetHist(histname)->GetXaxis()->SetBinLabel(7,"NoBJet");
    GetHist(histname)->GetXaxis()->SetBinLabel(8,"NoBJet_nolepveto");
    GetHist(histname)->GetXaxis()->SetBinLabel(9,"MET50");
    GetHist(histname)->GetXaxis()->SetBinLabel(10,"Mjj200");


  }

  if(!IsDATA){
    FillCutFlowByName("Unweighted_"+histname, cut, w, true);
  }

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
  FillHist("Nvtx_"+histsuffix, n_vtx, thisweight, 0., 100., 100);

  double HT(0.), ST(0.);
  for(unsigned int i=0; i<jets.size(); i++)     HT += jets.at(i).Pt();
  for(unsigned int i=0; i<jets_fwd.size(); i++) HT += jets_fwd.at(i).Pt();
  ST = HT;
  for(unsigned int i=0; i<leptons.size(); i++)  ST += leptons.at(i).Pt();

  FillHist("HT_"+histsuffix, HT, thisweight, 0., 2000., 2000);
  FillHist("ST_"+histsuffix, ST, thisweight, 0., 2000., 2000);

  FillHist("m_ll_"+histsuffix, (leptons.at(0)+leptons.at(1)).M(), thisweight, 0., 2000., 2000);

  for(int i=0; i<leptons.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"Lepton_Pt_"+histsuffix,  leptons.at(i).Pt(), thisweight, 0., 1000., 1000);
    FillHist(leporder[i]+"Lepton_Eta_"+histsuffix, leptons.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(int i=0; i<jets.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"Jet_Pt_"+histsuffix,  jets.at(i).Pt(), thisweight, 0., 1000., 1000);
    FillHist(leporder[i]+"Jet_Eta_"+histsuffix, jets.at(i).Eta(), thisweight, -3., 3., 60);
  }
  for(int i=0; i<jets_fwd.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"ForwardJet_Pt_"+histsuffix,  jets_fwd.at(i).Pt(), thisweight, 0., 1000., 1000);
    FillHist(leporder[i]+"ForwardJet_Eta_"+histsuffix, jets_fwd.at(i).Eta(), thisweight, -5., 5., 100);
  }
  for(int i=0; i<jets_nolepveto.size(); i++){
    if(i==4) break;
    FillHist(leporder[i]+"NoLepVetoJet_Pt_"+histsuffix, jets_nolepveto.at(i).Pt(), thisweight, 0., 1000., 1000);
    FillHist(leporder[i]+"NoLepVetoJet_Eta_"+histsuffix, jets_nolepveto.at(i).Eta(), thisweight, -3., 3., 60);
  }


  if(jets.size() >= 2){
    FillHist("m_jj_"+histsuffix, (jets.at(0)+jets.at(1)).M(), thisweight, 0., 1000., 1000);
    FillHist("m_Leadljj_"+histsuffix, (leptons.at(0)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_SubLeadljj_"+histsuffix, (leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
    FillHist("m_lljj_"+histsuffix, (leptons.at(0)+leptons.at(1)+jets.at(0)+jets.at(1)).M(), thisweight, 0., 2000., 2000);
  }

  if(thieweighterr!=0.){
    FillDiLeptonPlot(histsuffix+"_up",   leptons, jets, jets_fwd, jets_nolepveto, thisweight + thieweighterr, 0.);
    FillDiLeptonPlot(histsuffix+"_down", leptons, jets, jets_fwd, jets_nolepveto, thisweight - thieweighterr, 0.);
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

double DiLeptonAnalyzer::GetDijetMassClosest(std::vector<snu::KJet> js, double mass){

  if(js.size()<2) return -999;

  double m_close(9999999999999999.);
  for(unsigned int i=0; i<js.size()-1; i++){
    for(unsigned int j=i+1; j<js.size(); j++){

      double m_tmp = (js.at(i)+js.at(j)).M();
      if( fabs(m_tmp-mass) < fabs(m_close-mass) ){
        m_close = m_tmp;
      }

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
    if(!isT.at(i)) fr_onlyLoose.push_back( a.at(i) );
  }

  //==== Initialise weight
  float this_weight=-1.;

  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }

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




