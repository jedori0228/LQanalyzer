// $Id: FakeRateCalculator_El_dxysig_DILEP.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_El_dxysig_DILEP Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_El_dxysig_DILEP.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_El_dxysig_DILEP);

FakeRateCalculator_El_dxysig_DILEP::FakeRateCalculator_El_dxysig_DILEP() :  AnalyzerCore() {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_El_dxysig_DILEP");
  
  Message("In FakeRateCalculator_El_dxysig_DILEP constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_El_dxysig_DILEP::InitialiseAnalysis() throw( LQError ) {
  
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
  TFile *file_tmp = new TFile( lqdir+"/JskimData/HalfSample/FR_sampleA.root" );

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

  return;
}


void FakeRateCalculator_El_dxysig_DILEP::ExecuteEvents()throw( LQError ){

  double thisdXYCut = 0.01;
  double TightISO = 0.08;

  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.MET();

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

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_NOLEPTONVETO", 30., 5.);

  int n_jets = jetColl_hn.size();
  int n_bjets=0;

  int For_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v=0;
/*
  //==== HLT_Mu3_PFJet40_v PFJet Pt Check
  if( PassTrigger("HLT_Mu3_PFJet40_v") ){
    if(jetColl_hn.size()!=0){
      FillHist("HLT_Mu3_PFJet40_v_LeadingJetPt", jetColl_hn.at(0).Pt(), 1., 0., 200., 200);
    }
  }
  return;
*/
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jetColl_hn.at(j).Pt() > 40. ) For_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v++;
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

  //==========================
  //==== Call Loosest Electron
  //==========================

  std::vector<snu::KElectron> electrontriNodXYCutVLooseColl_raw = GetElectrons(false, true, "Electron_HN_NODXYCUT_VLOOSE_lowestPtCut");

  //==================
  //==== FR binnings
  //==================

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {0., 5., 10., 15., 20., 25., 30., 40., 45., 50., 70., 100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  //===========
  //==== HLTs
  //===========

  // HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v 6.992
  // HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v 6.162
  // HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v 30.397
  // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v 16.43

  std::map< TString, std::vector<double> > HLT_ptrange;
  HLT_ptrange["HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"].push_back(10.);
  HLT_ptrange["HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"].push_back(15.);
  HLT_ptrange["HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v"].push_back(15.);
  HLT_ptrange["HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v"].push_back(20.);
  HLT_ptrange["HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v"].push_back(20.);
  HLT_ptrange["HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v"].push_back(25.);
  HLT_ptrange["HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v"].push_back(25.);
  HLT_ptrange["HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v"].push_back(9999.);

  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }

  //==== 2) back-to-back dijet topology

  //==== tag jet collections
  //==== pt > 40 GeV
  //==== |eta| > 2.4
  //==== LeptonVeto
  std::vector<snu::KJet> jetColl_tag = GetJets("JET_HN", 40., 2.4);
  std::vector<snu::KElectron> TEST_tight = GetElectrons(false, false, "ELECTRON_HN_TIGHTv4");
  std::vector<snu::KElectron> TEST_loose = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");

  if( PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") || PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v") ){

    //==== dijet event seletion
    if(jetColl_tag.size() != 0){
      if(TEST_loose.size() == 1){

        Double_t this_weight_Loose = weight*GetPrescale(TEST_loose, PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"), PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v"));

        snu::KElectron electron = TEST_loose.at(0);

        float dR=999.9, dPhi=999.9, ptmu_ptjet=999.9;
        bool histfilled = false;
        for(unsigned int i=0; i<jetColl_tag.size(); i++){

          if(histfilled) break;
          dR = jetColl_tag.at(i).DeltaR(electron);
          dPhi = fabs(jetColl_tag.at(i).DeltaPhi(electron));
          ptmu_ptjet = electron.Pt() / jetColl_tag.at(i).Pt();
          FillHist("SingleElectronTrigger_Dijet_dR", dR, this_weight_Loose, 0., 10., 100);

          if( dR > 1.0 && ptmu_ptjet < 1. && dPhi > 2.5 ){
            FillDenAndNum("SingleElectronTrigger_Dijet", electron, this_weight_Loose, (TEST_tight.size() == 1) );
            if(TEST_tight.size() == 1){
              histfilled=true;
            }
          }

        } // tag jet for loop

      } // only one electron
    } // tag jet exist
  } // trigger

  //=========================================================
  //==== Large dXYSig Electron definitions for systematic study
  //=========================================================

  //FIXME decide after first run
  double dXYMin_central = 4.0;
  double RelIsoMax_central = 0.4;

  const int n_dXYMins = 3;
  double dXYMins[n_dXYMins] = {4.0, 4.5, 5.0};
  const int n_RelIsoMaxs = 6;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.2, 0.3, 0.4, 0.6, 0.8, 1.0};

/*
  //==== Central values only for test
  const int n_dXYMins = 1;
  double dXYMins[n_dXYMins] = {4.0};
  const int n_RelIsoMaxs = 1;
  double RelIsoMaxs[n_RelIsoMaxs] = {0.4};
*/
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

      //===============================
      //==== Prepare Electron Collections
      //===============================

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
        snu::KElectron this_electron = electrontriNodXYCutVLooseColl_raw.at(j);
        //==== MaxRelIso for Loose
        if( this_electron.PFRelIso(0.3) >= RelIsoMaxs[bbb] ) continue;

        //==== Fill NodXYCut
        electrontriNodXYCutLooseColl_raw.push_back( this_electron ); 
        if( this_electron.PFRelIso(0.3) < TightISO ) electrontriNodXYCutTightColl_raw.push_back( this_electron );
        //==== Fill Small dXYSig
        if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
          electrontriLooseColl_raw.push_back( this_electron );
          if( this_electron.PFRelIso(0.3) < TightISO ) electrontriTightColl_raw.push_back( this_electron );
        }
        //==== Fill Large dXYSig
        if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
          electrontriHighdXYLooseColl_raw.push_back( this_electron );
          if( this_electron.PFRelIso(0.3) < TightISO ) electrontriHighdXYTightColl_raw.push_back( this_electron );
        }

        //==== for data, or QCD
        //==== do the same for prompt collection
        if( k_isdata || k_sample_name.Contains("QCD") ){
          //==== Fill NodXYCut
          electrontriNodXYCutLooseColl.push_back( this_electron );
          if( this_electron.PFRelIso(0.3) < TightISO ) electrontriNodXYCutTightColl.push_back( this_electron );
          //==== Fill Small dXYSig
          if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
            electrontriLooseColl.push_back( this_electron );
            if( this_electron.PFRelIso(0.3) < TightISO ) electrontriTightColl.push_back( this_electron );
          }
          //==== Fill Large dXYSig
          if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
            electrontriHighdXYLooseColl.push_back( this_electron );
            if( this_electron.PFRelIso(0.3) < TightISO ) electrontriHighdXYTightColl.push_back( this_electron );
          }
        }

        //==== for MC,
        //==== fill only MCMatched() electrons
        else{
          if( TruthMatched(this_electron,false) ){
            //==== Fill NodXYCut
            electrontriNodXYCutLooseColl.push_back( this_electron );
            if( this_electron.PFRelIso(0.3) < TightISO ) electrontriNodXYCutTightColl.push_back( this_electron );
            //==== Fill Small dXYSig
            if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
              electrontriLooseColl.push_back( this_electron );
              if( this_electron.PFRelIso(0.3) < TightISO ) electrontriTightColl.push_back( this_electron );
            }
            //==== Fill Large dXYSig
            if( fabs(this_electron.dXY()) < 1. && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
              electrontriHighdXYLooseColl.push_back( this_electron );
              if( this_electron.PFRelIso(0.3) < TightISO ) electrontriHighdXYTightColl.push_back( this_electron );
            }
          }
        }


      } //==== Loosest electron loop

      //===================
      //==== dXYCut study
      //===================

      if( k_sample_name.Contains("DY") || k_sample_name.Contains("HN") ){

        for(unsigned int i=0; i<electrontriNodXYCutLooseColl.size(); i++){
          snu::KElectron this_electron = electrontriNodXYCutLooseColl.at(i);
          if(fabs(this_electron.dXY()) < 1.){
            FillHist(str_dXYCut+"_prompt_Loose_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_prompt_Loose_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            if(this_electron.PFRelIso(0.3)<TightISO){
              FillHist(str_dXYCut+"_prompt_Tight_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_prompt_Tight_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }
      if( k_sample_name.Contains("TTJets_aMC") || k_sample_name.Contains("QCD") ){

        for(unsigned int i=0; i<electrontriNodXYCutLooseColl_raw.size(); i++){
          snu::KElectron this_electron = electrontriNodXYCutLooseColl_raw.at(i);
          if( fabs(this_electron.dXY()) < 1. && !TruthMatched(this_electron, false) ){
            FillHist(str_dXYCut+"_fake_Loose_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_fake_Loose_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            if(this_electron.PFRelIso(0.3)<TightISO){
              FillHist(str_dXYCut+"_fake_Tight_dXYSig", fabs(this_electron.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_fake_Tight_dXY", fabs(this_electron.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }

      //=========================
      //==== SingleElectron Trigger
      //=========================

      if( PassTriggerOR(AllHLTs) ){

        std::map<TString, double> this_weight_Loose, this_weight_HighdXYLoose, this_weight_NodXYCutLoose;
        for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
          this_weight_Loose[it->first]         = weight*GetTriggerWeightByPtRange(it->first, it->second, electrontriLooseColl, For_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
          this_weight_HighdXYLoose[it->first]  = weight*GetTriggerWeightByPtRange(it->first, it->second, electrontriHighdXYLooseColl, For_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
          this_weight_NodXYCutLoose[it->first] = weight*GetTriggerWeightByPtRange(it->first, it->second, electrontriNodXYCutLooseColl, For_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
        }

        //Double_t this_weight_Loose = weight*GetPrescale(electrontriLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_HighdXYLoose = weight*GetPrescale(electrontriHighdXYLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_NodXYCutLoose = weight*GetPrescale(electrontriNodXYCutLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        //==== 1) LooseElectron study

        if( electrontriLooseColl.size() == 1 ){
          snu::KElectron electron = electrontriLooseColl.at(0);
          double LeptonRelIso = electron.PFRelIso(0.3);
          double conept = ElectronConePt(electron,TightISO);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_eta", electron.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_pt", electron.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_dXY", fabs(electron.dXY()), this_weight_Loose, 0., 0.1, 100);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_Loose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_dZ", fabs(electron.dZ()), this_weight_Loose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_LooseElectron_onebin", 0., this_weight_Loose, 0., 1., 1);
        }

        //==== 2) SingleElectron HighdXY LooseElectron study

        if( electrontriHighdXYLooseColl.size() == 1 ){
          snu::KElectron electron = electrontriHighdXYLooseColl.at(0);
          double LeptonRelIso = electron.PFRelIso(0.3);
          double conept = ElectronConePt(electron,TightISO);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_eta", electron.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_pt", electron.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_dXY", fabs(electron.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_dZ", fabs(electron.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseElectron_onebin", 0., this_weight_HighdXYLoose, 0., 1., 1);
        }

        //==== 2') SingleElectron NodXY Loose Electron study
        if( electrontriNodXYCutLooseColl.size() == 1 ){
          snu::KElectron electron = electrontriNodXYCutLooseColl.at(0);
          double LeptonRelIso = electron.PFRelIso(0.3);
          double conept = ElectronConePt(electron,TightISO);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_eta", electron.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_pt", electron.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_RelIso", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_dXY", fabs(electron.dXY()), this_weight_NodXYCutLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_dXYSig", fabs(electron.dXYSig()), this_weight_NodXYCutLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_dZ", fabs(electron.dZ()), this_weight_NodXYCutLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseElectron_onebin", 0., this_weight_NodXYCutLoose, 0., 1., 1);
        }


        //==== 3) large dXY electron method

        if(electrontriHighdXYLooseColl.size()==1){
          snu::KElectron HighdXYelectron = electrontriHighdXYLooseColl.at(0);
          double LeptonRelIso = HighdXYelectron.PFRelIso(0.3);

          //==== half sample test
          if(dXYMins[aaa]== dXYMin_central && RelIsoMaxs[bbb]==RelIsoMax_central){

            int EventNumber = eventbase->GetEvent().EventNumber();

            double pt_for_FR = HighdXYelectron.Pt();
            double eta_for_FR = fabs(HighdXYelectron.Eta());
            if(pt_for_FR>=60.) pt_for_FR=59.; //FIXME should be 70
            if(eta_for_FR>=2.5) eta_for_FR=2.2;

            snu::KEvent Evt = eventbase->GetEvent();
            double MET = Evt.MET();

            //==== sample A
            //==== obtain FR_A(pt, eta)
            if(abs(EventNumber%2==0)){
              FillHist("TEST_HalfSample", 0., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleA";
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0", 
                                HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              }

            }
            //==== sample B
            //==== 1) obtain FR_measured(MET)
            //==== 2-1) obtain MET distribution := MET_F0
            //==== 2-2) obtain MET weighted by FR_A(pt_B, eta_B) := MET_F_predicted
            //==== 2-3) MET_F_predicted / MET_F0 := FR_predicted(MET)
            else{
              FillHist("TEST_HalfSample", 1., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleB";
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0",
                                HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F0", MET, this_weight_HighdXYLoose, 0., 500., 500);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYelectron.Pt(), fabs(HighdXYelectron.Eta()), this_weight_HighdXYLoose, ptarray, 11, etaarray, 4);
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F", MET, this_weight_HighdXYLoose, 0., 500., 500);
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              }

              int this_FR_bin = FR_sampleA->FindBin(pt_for_FR,eta_for_FR);
              double this_FR_sampleA = FR_sampleA->GetBinContent(this_FR_bin);

              std::map<TString, double> this_weight_HighdXYLoose_FRweighted;
              for(std::map< TString, double >::iterator it=this_weight_HighdXYLoose.begin(); it!=this_weight_HighdXYLoose.end(); it++){
                this_weight_HighdXYLoose_FRweighted[it->first] = (it->second)*this_FR_sampleA;
              }

              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_Predicted", MET, this_weight_HighdXYLoose_FRweighted, 0., 500., 500);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_Predicted", n_jets, this_weight_HighdXYLoose_FRweighted, 0., 10., 10);

            }

          }

          bool isThisTight = (LeptonRelIso < TightISO);

          //==== all jet
          FillDenAndNum(str_dXYCut+"_HighdXY_alljet", HighdXYelectron, this_weight_HighdXYLoose, isThisTight);

          //==== no jet
          if( n_jets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0jet", HighdXYelectron, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withjet", HighdXYelectron, this_weight_HighdXYLoose, isThisTight);
          }

          //==== no b-jet
          if( n_bjets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0bjet", HighdXYelectron, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with b-jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withbjet", HighdXYelectron, this_weight_HighdXYLoose, isThisTight);
          }


        } // END one HighdXYSig Electron

        //==== 5) MC Truth
        //==== here, we pick FAKE electrons using GEN info.
        if( !k_isdata ){

          for(unsigned int i=0; i<electrontriNodXYCutLooseColl_raw.size(); i++){

            snu::KElectron electron = electrontriNodXYCutLooseColl_raw.at(i);

            //==== if prompt, skip
            if( TruthMatched(electron, false) ) continue;

            double LeptonRelIso = electron.PFRelIso(0.3);
            bool isThisTight = (LeptonRelIso < TightISO);

            //=================================================
            //==== 1) |dXY| < 1 cm, |dXY/err| > dXYMins[aaa].
            //=================================================

            if( fabs( electron.dXY() ) < 1. &&
                fabs( electron.dXYSig() ) > dXYMins[aaa] ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_alljet", electron, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0jet", electron, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withjet", electron, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0bjet", electron, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withbjet", electron, 1., isThisTight);
              }

            } // small dXYSig region


            //=========================================
            //==== 2) |dXY| < thisdXYCut cm, |dXY/err| < 4.
            //=========================================

            if( fabs( electron.dXY() ) < thisdXYCut &&
                fabs( electron.dXYSig() ) < 4. ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_alljet", electron, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0jet", electron, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withjet", electron, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0bjet", electron, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withbjet", electron, 1., isThisTight);
              }

            } // Large dXYSig Region


          } // END Electron loop
        } // END k_isdata


      } // SingleElectron trigger fired


      //==================
      //==== TagZ method
      //==================

      std::vector<TString> triggerlist;
      triggerlist.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

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
          bool i_electron_tight =  fabs(electrontriNodXYCutLooseColl.at(i).dXY())<thisdXYCut && fabs(electrontriNodXYCutLooseColl.at(i).dXYSig())<3.0 && electrontriNodXYCutLooseColl.at(i).PFRelIso(0.3)<TightISO;
          bool j_electron_tight =  fabs(electrontriNodXYCutLooseColl.at(j).dXY())<thisdXYCut && fabs(electrontriNodXYCutLooseColl.at(j).dXYSig())<3.0 && electrontriNodXYCutLooseColl.at(j).PFRelIso(0.3)<TightISO;

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

            snu::KElectron this_electron =  electrontriNodXYCutLooseColl.at(i);
            double conept = ElectronConePt(this_electron,TightISO);

            //==== Large dXYSig
            if( fabs(this_electron.dXY()) < 1.0 && fabs(this_electron.dXYSig()) > dXYMins[aaa] ){
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_eta_F0", this_electron.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_F0", this_electron.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_RelIso_F0", this_electron.PFRelIso(0.3), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXY_F0", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig_F0", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dZ_F0", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_pt_vs_eta_F0", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_pt_cone_vs_eta_F0", conept, fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);


              if( this_electron.PFRelIso(0.3) < TightISO ){
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_eta_F", this_electron.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_F", this_electron.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_RelIso_F", this_electron.PFRelIso(0.3), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXY_F", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dXYSig_F", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_dZ_F", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_pt_vs_eta_F", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Large_events_pt_cone_vs_eta_F", conept, fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);

              }
            }

            //==== Small dXYSig
            if( fabs(this_electron.dXY()) < thisdXYCut && fabs(this_electron.dXYSig()) < 3.0 ){
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_eta_F0", this_electron.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_F0", this_electron.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_RelIso_F0", this_electron.PFRelIso(0.3), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXY_F0", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig_F0", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dZ_F0", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_pt_vs_eta_F0", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_pt_cone_vs_eta_F0", conept, fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              if( this_electron.PFRelIso(0.3) < TightISO ){
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_eta_F", this_electron.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_F", this_electron.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_RelIso_F", this_electron.PFRelIso(0.3), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXY_F", fabs(this_electron.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dXYSig_F", fabs(this_electron.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_dZ_F", fabs(this_electron.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_pt_vs_eta_F", this_electron.Pt(), fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiElectronTrigger_ZTag_Small_events_pt_cone_vs_eta_F", conept, fabs(this_electron.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              }
            }

          }

        } //== END at least three electrons


      } //==== DiElectron Trigger



    } //==== RelIso loop
  } //==== dXYMin loop


   return;
}// End of execute event loop
  


void FakeRateCalculator_El_dxysig_DILEP::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_El_dxysig_DILEP::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");

  
  return;
  
}

FakeRateCalculator_El_dxysig_DILEP::~FakeRateCalculator_El_dxysig_DILEP() {
  
  Message("In FakeRateCalculator_El_dxysig_DILEP Destructor" , INFO);
  
}


void FakeRateCalculator_El_dxysig_DILEP::FillCutFlow(TString cut, float weight){

  
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


void FakeRateCalculator_El_dxysig_DILEP::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_El_dxysig_DILEP::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_El_dxysig_DILEPCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_El_dxysig_DILEP::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
}

float FakeRateCalculator_El_dxysig_DILEP::GetPrescale(std::vector<snu::KElectron> electron, bool passlow, bool passhigh){
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

double FakeRateCalculator_El_dxysig_DILEP::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, std::vector<snu::KElectron> electrons, int npfjet50){

  double prescale_trigger = 0.;

  if(electrons.size()==1){
    snu::KElectron electron = electrons.at(0);
    double minpt = ptrange.at(0);
    double maxpt = ptrange.at(1);

    if(PassTrigger(hltname)){
      if(electron.Pt() >= minpt && electron.Pt() < maxpt){
        prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
      }

      if(hltname=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
        if(npfjet50==0) prescale_trigger = 0.;
      }
    }

  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}

void FakeRateCalculator_El_dxysig_DILEP::FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value, it->second, xmin, xmax, nbins);
  }

}
void FakeRateCalculator_El_dxysig_DILEP::FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value1, value2, it->second, x, nbinsx, y, nbinsy);
  }

}

void FakeRateCalculator_El_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KElectron electron, double thisweight, bool isTight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {0., 5., 10., 15., 20., 25., 30., 35., 45., 60., 80., 100.};

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  double TightISO = 0.07;
  double conept = ElectronConePt(electron,TightISO);

  FillHist(prefix+"_eta_F0", electron.Eta(), thisweight, -3., 3., 30);
  FillHist(prefix+"_pt_F0", electron.Pt(), thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_F0", conept, thisweight, 0., 200., 200);
  FillHist(prefix+"_RelIso_F0", electron.PFRelIso(0.3), thisweight, 0., 1., 100);
  FillHist(prefix+"_dXY_F0", fabs(electron.dXY()), thisweight, 0., 1., 1000);
  FillHist(prefix+"_dXYSig_F0", fabs(electron.dXYSig()), thisweight, 0., 15., 150);
  FillHist(prefix+"_dZ_F0", fabs(electron.dZ()), thisweight, 0., 0.5, 50);
  FillHist(prefix+"_onebin_F0", 0., thisweight, 0., 1., 1);
  FillHist(prefix+"_events_pt_vs_eta_F0", electron.Pt(), fabs(electron.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(electron.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_PFMET_F0", METauto, thisweight, 0., 1000., 1000);

  if( isTight ){
    FillHist(prefix+"_eta_F", electron.Eta(), thisweight, -3., 3., 30);
    FillHist(prefix+"_pt_F", electron.Pt(), thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_F", conept, thisweight, 0., 200., 200);
    FillHist(prefix+"_RelIso_F", electron.PFRelIso(0.3), thisweight, 0., 1., 100);
    FillHist(prefix+"_dXY_F", fabs(electron.dXY()), thisweight, 0., 1., 1000);
    FillHist(prefix+"_dXYSig_F", fabs(electron.dXYSig()), thisweight, 0., 15., 150);
    FillHist(prefix+"_dZ_F", fabs(electron.dZ()), thisweight, 0., 0.5, 50);
    FillHist(prefix+"_onebin_F", 0., thisweight, 0., 1., 1);
    FillHist(prefix+"_events_pt_vs_eta_F", electron.Pt(), fabs(electron.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(electron.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);

  }


}

void FakeRateCalculator_El_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KElectron electron, std::map<TString, double> hltweight, bool isTight){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){

    FillDenAndNum(it->first+"_"+prefix, electron, it->second, isTight);

  }

}


















