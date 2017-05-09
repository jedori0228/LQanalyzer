// $Id: FakeRateCalculator_El_dxysig.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_El_dxysig Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_El_dxysig.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_El_dxysig);

FakeRateCalculator_El_dxysig::FakeRateCalculator_El_dxysig() :  AnalyzerCore(){
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_El_dxysig");
  
  Message("In FakeRateCalculator_El_dxysig constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_El_dxysig::InitialiseAnalysis() throw( LQError ) {
  
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

/*
  TDirectory* origDir = gDirectory;

  string lqdir = getenv("LQANALYZER_DIR");
  TFile *file_tmp = new TFile( (lqdir+"/data/Fake/80X/FR_sampleA.root").c_str() );

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

  FR_sampleA = (TH2D*)file_tmp->Get("FR_sampleA")->Clone();

  file_tmp->Close();
  delete file_tmp;

  origDir->cd();
*/

  return;
}


void FakeRateCalculator_El_dxysig::ExecuteEvents()throw( LQError ){

  double thisdXYCut = 0.01;

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

  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);

  int n_jets = jetColl_hn.size();
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
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

  //================================================
  //==== Apply corrections (PU reweight only here)
  //================================================

  if(!isData && !k_running_nonprompt){
    weight*=pileup_reweight;
    weight*=GetKFactor();
  }

  //============================
  //==== Call Loosest Electron
  //============================

  std::vector<snu::KElectron> electrontriNodXYCutVLooseColl_raw = GetElectrons(false, true, "ELECTRON_HN_TRI_NODXYCUT_VLOOSE");

  //==================
  //==== FR binnings
  //==================

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  //==== 2) back-to-back dijet topology

  //==== tag jet collections
  //==== pt > 40 GeV
  //==== |eta| > 2.4
  //==== LeptonVeto
  std::vector<snu::KJet> jetColl_tag = GetJets("JET_HN", 40., 2.4);
  std::vector<snu::KElectron> TEST_tight = GetElectrons("ELECTRON_HN_TRI_TIGHT");
  std::vector<snu::KElectron> TEST_loose = GetElectrons("ELECTRON_HN_TRI_LOOSE");

  if( PassTrigger("HLT_Mu17_v") || PassTrigger("HLT_Mu8_v") ){

    //==== dijet event seletion
    if(jetColl_tag.size() != 0){
      if(TEST_loose.size() == 1){

        Double_t this_weight_Loose = weight*GetPrescale(TEST_loose, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        KLepton electron = TEST_loose.at(0);

        float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;
        bool histfilled = false;
        for(unsigned int i=0; i<jetColl_tag.size(); i++){

          if(histfilled) break;
          dR = jetColl_tag.at(i).DeltaR(electron);
          dPhi = fabs(jetColl_tag.at(i).DeltaPhi(electron));
          ptmu_ptjet = electron.Pt() / jetColl_tag.at(i).Pt();
          FillHist("SingleElectronTrigger_Dijet_dR", dR, this_weight_Loose, 0., 10., 100);

          if( dR > 1.0 && ptmu_ptjet < 1. && dPhi > 2.5 ){
            FillHist("SingleElectronTrigger_Dijet_eta_F0", electron.Eta(), this_weight_Loose, -3, 3, 30);
            FillHist("SingleElectronTrigger_Dijet_pt_F0", electron.Pt(), this_weight_Loose, 0., 200., 200);
            FillHist("SingleElectronTrigger_Dijet_RelIso_F0", electron.RelIso(), this_weight_Loose, 0., 1., 100);
            FillHist("SingleElectronTrigger_Dijet_dXY_F0", fabs(electron.dXY()), this_weight_Loose, 0., 1., 1000);
            FillHist("SingleElectronTrigger_Dijet_dXYSig_F0", fabs(electron.dXYSig()), this_weight_Loose, 0., 15., 150);
            FillHist("SingleElectronTrigger_Dijet_dZ_F0", fabs(electron.dZ()), this_weight_Loose, 0., 0.5, 50);
            FillHist("SingleElectronTrigger_Dijet_events_F0", electron.Pt(), fabs(electron.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
            if(TEST_tight.size() == 1){
              histfilled=true;
              FillHist("SingleElectronTrigger_Dijet_eta_F", electron.Eta(), this_weight_Loose, -3, 3, 30);
              FillHist("SingleElectronTrigger_Dijet_pt_F", electron.Pt(), this_weight_Loose, 0., 200., 200);
              FillHist("SingleElectronTrigger_Dijet_RelIso_F", electron.RelIso(), this_weight_Loose, 0., 1., 100);
              FillHist("SingleElectronTrigger_Dijet_dXY_F", fabs(electron.dXY()), this_weight_Loose, 0., 1, 1000);
              FillHist("SingleElectronTrigger_Dijet_dXYSig_F", fabs(electron.dXYSig()), this_weight_Loose, 0., 15., 150);
              FillHist("SingleElectronTrigger_Dijet_dZ_F", fabs(electron.dZ()), this_weight_Loose, 0., 0.5, 50);
              FillHist("SingleElectronTrigger_Dijet_events_F", electron.Pt(), fabs(electron.Eta()), this_weight_Loose, ptarray, 9, etaarray, 4);
            }
          }

        } // tag jet for loop

      } // only one electron
    } // tag jet exist
  } // trigger

  //=============================================================
  //==== Large dXYSig Electron definitions for systematic study
  //=============================================================

  double dXYMin_central = 4.0;
  double RelIsoMax_central = 0.5;
  double TightIso = 0.05;

  const int n_dXYMins = 3;
  double dXYMins[n_dXYMins] = {3.0, 4.0, 5.0};
  const int n_RelIsoMaxs = 6;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.2, 0.4, 0.5, 0.6, 0.8, 1.0};

  //==== Central values only for test
  //const int n_dXYMins = 1;
  //double dXYMins[n_dXYMins] = {4.0};
  //const int n_RelIsoMaxs = 1;
  //double RelIsoMaxs[n_RelIsoMaxs] = {0.4};

  //=================
  //==== Start Loop
  //=================

  for(int aaa=0; aaa<n_dXYMins; aaa++){
    for(int bbb=0; bbb<n_RelIsoMaxs; bbb++){

      int dXY_Digit1 = int(dXYMins[aaa]);
      int dXY_Digit0p1 = 10*dXYMins[aaa]-10*dXY_Digit1;
      TString str_dXYCut = "dXYSigMin_"+TString::Itoa(dXY_Digit1,10)+"p"+TString::Itoa(dXY_Digit0p1,10);

      int iso_Digit1 = int(RelIsoMaxs[bbb]);
      int iso_Digit0p1 = 10*RelIsoMaxs[bbb]-10*iso_Digit1;
      TString str_iso = "LooseRelIsoMax_"+TString::Itoa(iso_Digit1,10)+"p"+TString::Itoa(iso_Digit0p1,10);

      str_dXYCut = str_dXYCut+"_"+str_iso;

      //===================================
      //==== Prepare Electron Collections
      //===================================

      //==== NodXYCut electrons
      std::vector<snu::KElectron> electrontriNodXYCutLooseColl_raw;
      std::vector<snu::KElectron> electrontriNodXYCutTightColl_raw;
      std::vector<snu::KElectron> electrontriNodXYCutLooseColl;
      std::vector<snu::KElectron> electrontriNodXYCutTightColl;
      //==== Small dXYSig electrons
      std::vector<snu::KElectron> electrontriLooseColl_raw;
      std::vector<snu::KElectron> electrontriTightColl_raw;
      std::vector<snu::KElectron> electrontriLooseColl;
      std::vector<snu::KElectron> electrontriTightColl;
      //==== Large dXYSig electrons
      std::vector<snu::KElectron> electrontriHighdXYLooseColl_raw;
      std::vector<snu::KElectron> electrontriHighdXYTightColl_raw;
      std::vector<snu::KElectron> electrontriHighdXYLooseColl;
      std::vector<snu::KElectron> electrontriHighdXYTightColl;

      for(unsigned int j=0; j<electrontriNodXYCutVLooseColl_raw.size(); j++){
        KLepton this_electron = electrontriNodXYCutVLooseColl_raw.at(j);
        snu::KElectron eltopush = electrontriNodXYCutVLooseColl_raw.at(j);
        //==== MaxRelIso for Loose
        if( this_electron.RelIso() >= RelIsoMaxs[bbb] ) continue;

        //==== Fill NodXYCut
        electrontriNodXYCutLooseColl_raw.push_back( eltopush ); 
        if( this_electron.RelIso() < TightIso ) electrontriNodXYCutTightColl_raw.push_back( eltopush );
        //==== Fill Small dXYSig
        if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
          electrontriLooseColl_raw.push_back( eltopush );
          if( this_electron.RelIso() < TightIso ) electrontriTightColl_raw.push_back( eltopush );
        }
        //==== Fill Large dXYSig
        if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
          electrontriHighdXYLooseColl_raw.push_back( eltopush );
          if( this_electron.RelIso() < TightIso ) electrontriHighdXYTightColl_raw.push_back( eltopush );
        }

        //==== for data, or QCD
        //==== do the same for prompt collection
        if( k_isdata || k_sample_name.Contains("QCD") ){
          //==== Fill NodXYCut
          electrontriNodXYCutLooseColl.push_back( eltopush );
          if( this_electron.RelIso() < TightIso ) electrontriNodXYCutTightColl.push_back( eltopush );
          //==== Fill Small dXYSig
          if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
            electrontriLooseColl.push_back( eltopush );
            if( this_electron.RelIso() < TightIso ) electrontriTightColl.push_back( eltopush );
          }
          //==== Fill Large dXYSig
          if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
            electrontriHighdXYLooseColl.push_back( eltopush );
            if( this_electron.RelIso() < TightIso ) electrontriHighdXYTightColl.push_back( eltopush );
          }
        }

        //==== for MC,
        //==== fill only MCIsPrompt() electrons
        else{
          if( this_electron.GetElectronPtr()->MCIsPrompt() ){
            //==== Fill NodXYCut
            electrontriNodXYCutLooseColl.push_back( eltopush );
            if( this_electron.RelIso() < TightIso ) electrontriNodXYCutTightColl.push_back( eltopush );
            //==== Fill Small dXYSig
            if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
              electrontriLooseColl.push_back( eltopush );
              if( this_electron.RelIso() < TightIso ) electrontriTightColl.push_back( eltopush );
            }
            //==== Fill Large dXYSig
            if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
              electrontriHighdXYLooseColl.push_back( eltopush );
              if( this_electron.RelIso() < TightIso ) electrontriHighdXYTightColl.push_back( eltopush );
            }
          }
        }


      } //==== Loosest electron loop

      //===================
      //==== dXYCut study
      //===================

      if( k_sample_name.Contains("DY") || k_sample_name.Contains("HN") ){

        for(unsigned int i=0; i<electrontriNodXYCutLooseColl.size(); i++){
          KLepton this_electron = electrontriNodXYCutLooseColl.at(i);
          if(fabs(this_electron.dXY()) < 1.){
            FillHist(str_dXYCut+"_prompt_Loose_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_prompt_Loose_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            if(this_electron.RelIso()<TightIso){
              FillHist(str_dXYCut+"_prompt_Tight_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_prompt_Tight_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }
      if( k_sample_name.Contains("TTJets_aMC") || k_sample_name.Contains("QCD") ){

        for(unsigned int i=0; i<electrontriNodXYCutLooseColl_raw.size(); i++){
          KLepton this_electron = electrontriNodXYCutLooseColl_raw.at(i);
          if( fabs(this_electron.dXY()) < 1. && !this_electron.GetElectronPtr()->MCIsPrompt() ){
            FillHist(str_dXYCut+"_fake_Loose_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_fake_Loose_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            if(this_electron.RelIso()<TightIso){
              FillHist(str_dXYCut+"_fake_Tight_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_fake_Tight_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }

      //=========================
      //==== SingleElectron Trigger
      //=========================

      bool highpt_trigger(false), lowpt_trigger(false);
      highpt_trigger = PassTrigger("HLT_Mu17_v");
      lowpt_trigger = PassTrigger("HLT_Mu8_v");

      if( highpt_trigger || lowpt_trigger ){

        Double_t this_weight_Loose = weight*GetPrescale(electrontriLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        Double_t this_weight_HighdXYLoose = weight*GetPrescale(electrontriHighdXYLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        Double_t this_weight_NodXYCutLoose = weight*GetPrescale(electrontriNodXYCutLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        //==== 1) LooseElectron study

        if( electrontriLooseColl.size() == 1 ){
          KLepton electron = electrontriLooseColl.at(0);
          double LeptonRelIso = electron.RelIso();
          FillHist(str_dXYCut+"_SingleElectronTrigger_LooseElectron_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleElectronTrigger_LooseElectron_dXY", fabs(electron.dXY()), this_weight_Loose, 0., 0.1, 100);
          FillHist(str_dXYCut+"_SingleElectronTrigger_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_Loose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleElectronTrigger_LooseElectron_dZ", fabs(electron.dZ()), this_weight_Loose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleElectronTrigger_LooseElectron_onebin", 0., this_weight_Loose, 0., 1., 1);
        }

        //==== 2) SingleElectron HighdXY LooseElectron study

        if( electrontriHighdXYLooseColl.size() == 1 ){
          KLepton electron = electrontriHighdXYLooseColl.at(0);
          double LeptonRelIso = electron.RelIso();
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_LooseElectron_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_LooseElectron_dXY", fabs(electron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_LooseElectron_dZ", fabs(electron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_LooseElectron_onebin", 0., this_weight_HighdXYLoose, 0., 1., 1);
        }

        //==== 2') SingleElectron NodXY Loose Electron study
        if( electrontriNodXYCutLooseColl.size() == 1 ){
          KLepton electron = electrontriNodXYCutLooseColl.at(0);
          double LeptonRelIso = electron.RelIso();
          FillHist(str_dXYCut+"_SingleElectronTrigger_NodXY_LooseElectron_RelIso", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleElectronTrigger_NodXY_LooseElectron_dXY", fabs(electron.dXY()), this_weight_NodXYCutLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleElectronTrigger_NodXY_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_NodXYCutLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleElectronTrigger_NodXY_LooseElectron_dZ", fabs(electron.dZ()), this_weight_NodXYCutLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleElectronTrigger_NodXY_LooseElectron_onebin", 0., this_weight_NodXYCutLoose, 0., 1., 1);
        }


        //==== 3) large dXY electron method

        if(electrontriHighdXYLooseColl.size()==1){
          KLepton HighdXYelectron = electrontriHighdXYLooseColl.at(0);
          double LeptonRelIso = HighdXYelectron.RelIso();

/*
          //==== half sample test
          if(dXYMins[aaa]== dXYMin_central && RelIsoMaxs[bbb]==RelIsoMax_central){

            int EventNumber = eventbase->GetEvent().EventNumber();

            double pt_for_FR = HighdXYelectron.Pt();
            double eta_for_FR = fabs(HighdXYelectron.Eta());
            if(pt_for_FR>=60.) pt_for_FR=59.;
            if(eta_for_FR>=2.5) eta_for_FR=2.2;

            int this_FR_bin = FR_sampleA->FindBin(pt_for_FR,eta_for_FR);
            double this_FR_sampleA = FR_sampleA->GetBinContent(this_FR_bin);

            snu::KEvent Evt = eventbase->GetEvent();
            double MET = Evt.MET();

            //==== sample A
            if(abs(EventNumber%2==0)){
              FillHist("TEST_HalfSample", 0., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleA";
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F0", MET, this_weight_HighdXYLoose, 0., 500., 500);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);

              if( LeptonRelIso < TightIso ){
                FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              }

              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_Predicted", MET, this_weight_HighdXYLoose*this_FR_sampleA, 0., 500., 500);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_Predicted", n_jets, this_weight_HighdXYLoose*this_FR_sampleA, 0., 10., 10);


            }
            //==== sample B
            else{
              FillHist("TEST_HalfSample", 1., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleB";
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F0", MET, this_weight_HighdXYLoose, 0., 500., 500);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              if( LeptonRelIso < TightIso ){
                FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
                FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F", MET, this_weight_HighdXYLoose, 0., 500., 500);
                FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              }
            }

          }
*/

          //==== all jet

          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_eta_F0", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_pt_F0", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dXY_F0", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dXYSig_F0", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dZ_F0", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
          FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
          if( LeptonRelIso < TightIso ){
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_eta_F", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_pt_F", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dXY_F", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dXYSig_F", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_dZ_F", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
          }

          //==== no jet

          if( n_jets == 0){
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_eta_F0", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_pt_F0", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dXY_F0", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dXYSig_F0", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dZ_F0", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < TightIso ){
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_eta_F", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_pt_F", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dXY_F", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dXYSig_F", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_dZ_F", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0jet_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
          }

          //==== with jet

          else{
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_eta_F0", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_pt_F0", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dXY_F0", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dXYSig_F0", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dZ_F0", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < TightIso ){
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_eta_F", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_pt_F", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dXY_F", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dXYSig_F", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_dZ_F", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withjet_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
          }

          //==== no b-jet

          if( n_bjets == 0){
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_eta_F0", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_pt_F0", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dXY_F0", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dXYSig_F0", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dZ_F0", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < TightIso ){
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_eta_F", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_pt_F", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dXY_F", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dXYSig_F", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_dZ_F", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_0bjet_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
          }

          //==== with b-jet

          else{
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_eta_F0", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_pt_F0", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_RelIso_F0", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100); 
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dXY_F0", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000); 
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dXYSig_F0", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dZ_F0", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_onebin_F0", 0., this_weight_HighdXYLoose, 0., 1., 1);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_n_jets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
            FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_events_F0", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            if( LeptonRelIso < TightIso ){
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_eta_F", HighdXYelectron.Eta(), this_weight_HighdXYLoose, -3., 3., 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_pt_F", HighdXYelectron.Pt(), this_weight_HighdXYLoose, 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_RelIso_F", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dXY_F", fabs(HighdXYelectron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dXYSig_F", fabs(HighdXYelectron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_dZ_F", fabs(HighdXYelectron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_onebin_F", 0., this_weight_HighdXYLoose, 0., 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_n_jets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              FillHist(str_dXYCut+"_SingleElectronTrigger_HighdXY_withbjet_events_F", HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 9, etaarray, 4);
            }
          }
        }

        //==== 5) MC Truth
        //==== here, we pick FAKE electrons using GEN info.
        if( !k_isdata ){

          for(unsigned int i=0; i<electrontriNodXYCutLooseColl_raw.size(); i++){

            KLepton electron = electrontriNodXYCutLooseColl_raw.at(i);

            //==== if prompt, skip
            if( electron.GetElectronPtr()->MCIsPrompt() ) continue;

            double LeptonRelIso = electron.RelIso();

            //=================================================
            //==== 1) |dXY| < 1 cm, |dXY/err| > dXYMins[aaa].
            //=================================================

            if( fabs( electron.dXY() ) < 1. &&
                fabs( electron.dXYSig() ) > dXYMins[aaa] ){

              //==== all jet

              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_eta_F0", electron.Eta(), 1., -3, 3, 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_pt_F0", electron.Pt(), 1., 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dXY_F0", fabs(electron.dXY()),1., 0., 1., 1000);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dXYSig_F0", fabs(electron.dXYSig()),1., 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dZ_F0", fabs(electron.dZ()), 1., 0., 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_onebin_F0", 0., 1., 0., 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
              if( LeptonRelIso < TightIso ){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_eta_F", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_pt_F", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dXY_F", fabs(electron.dXY()), 1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_dZ_F", fabs(electron.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_onebin_F", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
              }

              //==== no jet

              if(n_jets==0){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dXY_F0", fabs(electron.dXY()),1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dXYSig_F0", fabs(electron.dXYSig()),1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dZ_F0", fabs(electron.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_onebin_F0", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dXY_F", fabs(electron.dXY()), 1., 0., 1., 1000);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_dZ_F", fabs(electron.dZ()), 1., 0., 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_onebin_F", 0., 1., 0., 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0jet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== with jet

              else{
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dXY_F0", fabs(electron.dXY()),1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dXYSig_F0", fabs(electron.dXYSig()),1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dZ_F0", fabs(electron.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_onebin_F0", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dXY_F", fabs(electron.dXY()), 1., 0., 1., 1000);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_dZ_F", fabs(electron.dZ()), 1., 0., 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_onebin_F", 0., 1., 0., 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== no b-jet

              if(n_bjets==0){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dXY_F0", fabs(electron.dXY()),1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dXYSig_F0", fabs(electron.dXYSig()),1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dZ_F0", fabs(electron.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_onebin_F0", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dXY_F", fabs(electron.dXY()), 1., 0., 1., 1000);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_dZ_F", fabs(electron.dZ()), 1., 0., 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_onebin_F", 0., 1., 0., 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_0bjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== with b-jet

              else{
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dXY_F0", fabs(electron.dXY()),1., 0., 1., 1000);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dXYSig_F0", fabs(electron.dXYSig()),1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dZ_F0", fabs(electron.dZ()), 1., 0., 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_onebin_F0", 0., 1., 0., 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dXY_F", fabs(electron.dXY()), 1., 0., 1., 1000);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_dZ_F", fabs(electron.dZ()), 1., 0., 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_onebin_F", 0., 1., 0., 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_HighdXY_withbjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

            } // small dXYSig region


            //=========================================
            //==== 2) |dXY| < thisdXYCut cm, |dXY/err| < 3.
            //=========================================

            if( fabs( electron.dXY() ) < thisdXYCut &&
                fabs( electron.dXYSig() ) < 3. ){

              //==== all jet

              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_eta_F0", electron.Eta(), 1., -3, 3, 30);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_pt_F0", electron.Pt(), 1., 0., 200., 200);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dXY_F0", fabs(electron.dXY()), 1., 0., 0.1, 100);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dXYSig_F0", fabs(electron.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dZ_F0", fabs(electron.dZ()), 1., 0, 0.5, 50);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_onebin_F0", 0., 1., 0, 1., 1);
              FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
              if( LeptonRelIso < TightIso ){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_eta_F", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_pt_F", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dXY_F", fabs(electron.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_dZ_F", fabs(electron.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_onebin_F", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
              }

              //==== no jet

              if(n_jets==0){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dXY_F0", fabs(electron.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dXYSig_F0", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dZ_F0", fabs(electron.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_onebin_F0", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dXY_F", fabs(electron.dXY()), 1., 0., 0.1, 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_dZ_F", fabs(electron.dZ()), 1., 0, 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_onebin_F", 0., 1., 0, 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0jet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== with jet

              else{
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dXY_F0", fabs(electron.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dXYSig_F0", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dZ_F0", fabs(electron.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_onebin_F0", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dXY_F", fabs(electron.dXY()), 1., 0., 0.1, 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_dZ_F", fabs(electron.dZ()), 1., 0, 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_onebin_F", 0., 1., 0, 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== no b-jet

              if(n_bjets==0){
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dXY_F0", fabs(electron.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dXYSig_F0", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dZ_F0", fabs(electron.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_onebin_F0", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dXY_F", fabs(electron.dXY()), 1., 0., 0.1, 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_dZ_F", fabs(electron.dZ()), 1., 0, 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_onebin_F", 0., 1., 0, 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_0bjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

              //==== with b-jet

              else{
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_eta_F0", electron.Eta(), 1., -3, 3, 30);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_pt_F0", electron.Pt(), 1., 0., 200., 200);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_RelIso_F0", LeptonRelIso, 1., 0., 1., 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dXY_F0", fabs(electron.dXY()), 1., 0., 0.1, 100);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dXYSig_F0", fabs(electron.dXYSig()), 1., 0., 15., 150);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dZ_F0", fabs(electron.dZ()), 1., 0, 0.5, 50);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_onebin_F0", 0., 1., 0, 1., 1);
                FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_events_F0", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                if( LeptonRelIso < TightIso ){
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_eta_F", electron.Eta(), 1., -3, 3, 30);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_pt_F", electron.Pt(), 1., 0., 200., 200);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_RelIso_F", LeptonRelIso, 1., 0., 1., 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dXY_F", fabs(electron.dXY()), 1., 0., 0.1, 100);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dXYSig_F", fabs(electron.dXYSig()), 1., 0., 15., 150);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_dZ_F", fabs(electron.dZ()), 1., 0, 0.5, 50);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_onebin_F", 0., 1., 0, 1., 1);
                  FillHist(str_dXYCut+"_SingleElectronTrigger_MCTruth_withbjet_events_F", electron.Pt(), fabs(electron.Eta()), 1., ptarray, 9, etaarray, 4);
                }
              }

            } // Large dXYSig Region


          } // END Electron loop
        } // END k_isdata


      } // SingleElectron trigger fired


      //==================
      //==== TagZ method
      //==================

      std::vector<TString> triggerlist;
      triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

      if( PassTriggerOR(triggerlist) ){
        float trigger_ps_weight = WeightByTrigger(triggerlist, TargetLumi);
        double this_weight = weight*trigger_ps_weight;

        int n_nodxycut_loose = electrontriNodXYCutLooseColl.size();
        if( n_nodxycut_loose >= 3){

          //==== Find Z-pair
          double m_OS_Zclosest = (91.2+15.0);
          int z_pair_index[2];
          bool isZTagged = false;

          for(unsigned i=0; i<n_nodxycut_loose-1; i++){
            for(unsigned j=i+1; j<n_nodxycut_loose; j++){
              if( electrontriNodXYCutLooseColl.at(i).Charge() == electrontriNodXYCutLooseColl.at(j).Charge() ) continue;
          bool i_electron_tight =  fabs(electrontriNodXYCutLooseColl.at(i).dxy())<thisdXYCut && fabs(electrontriNodXYCutLooseColl.at(i).dxySig())<3.0 && electrontriNodXYCutLooseColl.at(i).PFRelIso(0.3)<TightIso;
          bool j_electron_tight =  fabs(electrontriNodXYCutLooseColl.at(j).dxy())<thisdXYCut && fabs(electrontriNodXYCutLooseColl.at(j).dxySig())<3.0 && electrontriNodXYCutLooseColl.at(j).PFRelIso(0.3)<TightIso;

              double m_thisOS = ( electrontriNodXYCutLooseColl.at(i)+electrontriNodXYCutLooseColl.at(j) ).M();
              if(  i_electron_tight &&  j_electron_tight && fabs(m_thisOS-91.2) < fabs(m_OS_Zclosest-91.2) ){
                isZTagged = true;
                m_OS_Zclosest = m_thisOS;
                z_pair_index[0] = i;
                z_pair_index[1] = j;
              }
            }
          }
          if(!isZTagged) continue;
          for(unsigned i=0; i<n_nodxycut_loose; i++){
            if(i==z_pair_index[0] || i==z_pair_index[1]) continue;

            KLepton this_electron =  electrontriNodXYCutLooseColl.at(i);

            //==== Large dXYSig
            if( fabs(this_electron.dXY()) < 1.0 && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_eta_F0", this_electron.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_F0", this_electron.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_RelIso_F0", this_electron.RelIso(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXY_F0", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig_F0", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dZ_F0", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_F0", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);


              if( this_electron.RelIso() < TightIso ){
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_eta_F", this_electron.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_F", this_electron.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_RelIso_F", this_electron.RelIso(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXY_F", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig_F", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dZ_F", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_F", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);

              }
            }

            //==== Small dXYSig
            if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_eta_F0", this_electron.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_F0", this_electron.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_RelIso_F0", this_electron.RelIso(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXY_F0", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig_F0", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dZ_F0", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_F0", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              if( this_electron.RelIso() < TightIso ){
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_eta_F", this_electron.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_F", this_electron.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_RelIso_F", this_electron.RelIso(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXY_F", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig_F", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dZ_F", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_F", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              }
            }

          }

        } //== END at least three electrons


      } //==== DiElectron Trigger



    } //==== RelIso loop
  } //==== dXYMin loop


/*
  if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ){
    float trigger_ps_weight = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);
    double this_weight = weight*trigger_ps_weight;
    std::vector<snu::KMuon> PtEtaMuons = GetMuons("MUON_PTETA", false);
    for(unsigned int i=0; i<PtEtaMuons.size(); i++){
      KLepton this_muon = PtEtaMuons.at(i);
      if( this_muon.GetElectronPtr()->MCIsPrompt() ){
        FillHist("DiElectronTrigger_PromptRate_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
        FillHist("DiElectronTrigger_PromptRate_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
        FillHist("DiElectronTrigger_PromptRate_RelIso_F0", this_muon.RelIso(), this_weight, 0., 1., 100);
        FillHist("DiElectronTrigger_PromptRate_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
        FillHist("DiElectronTrigger_PromptRate_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
        FillHist("DiElectronTrigger_PromptRate_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
        FillHist("DiElectronTrigger_PromptRate_onebin_F0", 0., this_weight, 0, 1., 1);
        FillHist("DiElectronTrigger_PromptRate_events_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
        if( this_muon.RelIso() < TightIso ){
          FillHist("DiElectronTrigger_PromptRate_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
          FillHist("DiElectronTrigger_PromptRate_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
          FillHist("DiElectronTrigger_PromptRate_RelIso_F", this_muon.RelIso(), this_weight, 0., 1., 100);
          FillHist("DiElectronTrigger_PromptRate_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
          FillHist("DiElectronTrigger_PromptRate_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
          FillHist("DiElectronTrigger_PromptRate_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
          FillHist("DiElectronTrigger_PromptRate_onebin_F", 0., this_weight, 0, 1., 1);
          FillHist("DiElectronTrigger_PromptRate_events_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray, 9, etaarray, 4);
        }
      }
    }
  }
*/

   return;
}// End of execute event loop
  


void FakeRateCalculator_El_dxysig::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_El_dxysig::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_El_dxysig::~FakeRateCalculator_El_dxysig() {
  
  Message("In FakeRateCalculator_El_dxysig Destructor" , INFO);
  
}


void FakeRateCalculator_El_dxysig::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 3,0.,3.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");

  }
}


void FakeRateCalculator_El_dxysig::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_El_dxysig::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_El_dxysigCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_El_dxysig::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_electrons.clear();
}

float FakeRateCalculator_El_dxysig::GetPrescale(std::vector<snu::KElectron> electron, bool passlow, bool passhigh){
  float prescale_trigger = 0.;

  if(electron.size() == 1){

    if(electron.at(0).Pt() >= 20.){
      if(passhigh){
        prescale_trigger = WeightByTrigger("HLT_Mu17_v", TargetLumi) ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = WeightByTrigger("HLT_Mu8_v", TargetLumi) ;
      }
      else prescale_trigger = 0.;
    }
  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}

