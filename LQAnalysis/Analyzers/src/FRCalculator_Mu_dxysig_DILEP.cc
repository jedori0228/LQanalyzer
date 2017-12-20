// $Id: FRCalculator_Mu_dxysig_DILEP.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFRCalculator_Mu_dxysig_DILEP Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FRCalculator_Mu_dxysig_DILEP.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FRCalculator_Mu_dxysig_DILEP);

FRCalculator_Mu_dxysig_DILEP::FRCalculator_Mu_dxysig_DILEP() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FRCalculator_Mu_dxysig_DILEP");
  
  Message("In FRCalculator_Mu_dxysig_DILEP constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FRCalculator_Mu_dxysig_DILEP::InitialiseAnalysis() throw( LQError ) {
  
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


void FRCalculator_Mu_dxysig_DILEP::ExecuteEvents()throw( LQError ){

  //if(!(eventbase->GetEvent().RunNumber()==273450 && eventbase->GetEvent().EventNumber()==725619133)) return;

  double thisdXYCut = 0.005;
  double TightISO = 0.07;

  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.MET();
  METphiauto = Evt.METPhi();

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

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_NOLEPTONVETO", 20., 2.5);

  int n_jets = jetColl_hn.size();
  int n_bjets=0;

  int For_HLT_Mu3_PFJet40_v=0;
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jetColl_hn.at(j).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
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
  //==== Call Loosest Muon
  //==========================

  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO");
  std::vector<snu::KElectron> electrons_veto = GetElectrons("ELECTRON_HN_VETO");

  if(electrons_veto.size()!=0) return;

  std::vector<snu::KMuon> muontriNodXYCutVLooseColl_raw = GetMuons("MUON_HN_NODXYCUT_VLOOSE_lowestPtCut", true);

  bool DoHighdXY = std::find(k_flags.begin(), k_flags.end(), "DoHighdXY") != k_flags.end();
  bool DoSmallHighdXY = std::find(k_flags.begin(), k_flags.end(), "DoSmallHighdXY") != k_flags.end();

  TString LooseID = "MUON_HN_LOOSEv7_SIP3";
  if(DoHighdXY){
    if(DoSmallHighdXY) LooseID = "MUON_HN_Loose_HighdXY_Small";
    else               LooseID = "MUON_HN_Loose_HighdXY";
  }
  //cout << "LooseID = " << LooseID << endl;
  //LooseID = "MUON_HN_LOOSE_8TeV";

  std::vector<snu::KMuon> hnloose_raw = GetMuons(LooseID, true);

  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSE", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_LoosenSIP", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_POGIP", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv5", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv6", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv7_SIP10", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv7_SIP3p5", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv7_SIP3", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_Loose_HighdXY", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv7_SIP9", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv8", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_POGIP_SIP6", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_POGIP_SIP8", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_POGIP_SIP5", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv2_POGIP_SIP4p5", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv3", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv4", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv9", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv10", true);
  //std::vector<snu::KMuon> hnloose_raw = GetMuons("MUON_HN_LOOSEv11", true);

  JSCorrectedMETRochester(muons_veto, METauto, METphiauto);

  //==================
  //==== FR binnings
  //==================

  const int n_eta = 3;
  float etaarray[n_eta+1] = {0.0, 0.8, 1.479, 2.5};

  //==== Data Binning
  const int n_pt = 8;
  float ptarray[n_pt+1] = {5., 10., 15., 20., 30., 40., 50., 60., 70.};

/*
  //==== MC Binning
  const int n_pt = 9;
  float ptarray[n_pt+1] = {5., 10., 15., 20., 25., 30., 40., 50., 60., 70.};
*/

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};

  //===========
  //==== HLTs
  //===========

  std::map< TString, std::vector<double> > HLT_ptconerange;
  std::map< TString, double > HLT_ptmin;
  HLT_ptconerange.clear();

  HLT_ptconerange["HLT_Mu3_PFJet40_v"].push_back(5.); // HLT_Mu3_PFJet40_v 7.408
  HLT_ptconerange["HLT_Mu3_PFJet40_v"].push_back(15.);
  HLT_ptmin["HLT_Mu3_PFJet40_v"] = 5.; // -> 6.65 GeV

  HLT_ptconerange["HLT_Mu8_TrkIsoVVL_v"].push_back(15.); // HLT_Mu8_TrkIsoVVL_v 7.832
  HLT_ptconerange["HLT_Mu8_TrkIsoVVL_v"].push_back(30.);
  HLT_ptmin["HLT_Mu8_TrkIsoVVL_v"] = 10.; // -> 13.3 GeV

  HLT_ptconerange["HLT_Mu17_TrkIsoVVL_v"].push_back(30.); // HLT_Mu17_TrkIsoVVL_v 217.553
  HLT_ptconerange["HLT_Mu17_TrkIsoVVL_v"].push_back(9999.);
  HLT_ptmin["HLT_Mu17_TrkIsoVVL_v"] = 20.; //-> 26.6 GeV

  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptconerange.begin(); it!=HLT_ptconerange.end(); it++){
    AllHLTs.push_back(it->first);
  }

  //==== 2) back-to-back dijet topology

  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end();
  bool DijetPrompt= std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end();

  //==== Norm Cehck for each trigger
  if(DijetPrompt && !DoHighdXY){

    double m_Z = 91.1876;

    std::vector< snu::KMuon > tightmuons = GetMuons("MUON_HN_TIGHT", true, 10., 2.4);

    double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", tightmuons, 0);
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(tightmuons);

    TLorentzVector metvec;
    metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto );

    for(unsigned int i=0; i<AllHLTs.size(); i++){
      TString ThisTrigger = AllHLTs.at(i);
      if( PassTrigger(ThisTrigger) && (tightmuons.size()!=0) ){

        if(HLT_ptmin[ThisTrigger] > tightmuons.at(0).Pt()) continue;

        double triggerweight = JSWeightByTrigger(ThisTrigger, TargetLumi);
        if(ThisTrigger.Contains("PFJet40")){
          if(For_HLT_Mu3_PFJet40_v==0) triggerweight = 0.;
        }
        double this_weight = weight*muon_id_iso_sf*MuTrkEffSF*triggerweight;
        //cout << ThisTrigger << "\t" << weight << "\t" << this_weight << endl; //FIXME
        //==== 1) Z-Peak
        if(tightmuons.size()==2){
          double mll = (tightmuons.at(0)+tightmuons.at(1)).M();
          if( (tightmuons.at(0).Pt() > 20.) && (tightmuons.at(1).Pt() > 10.) && (fabs(mll-m_Z) < 10.) ){
            FillHist(ThisTrigger+"_ZPeak_mll", mll, this_weight, 0., 200., 200);
            FillHist(ThisTrigger+"_ZPeak_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_ZPeak_subleadpt", tightmuons.at(1).Pt(), this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_ZPeak_PFMET", METauto, this_weight, 0., 500., 500);
          }
        }
        //==== 2) W
        if(tightmuons.size()==1){
          double MTval = AnalyzerCore::MT( tightmuons.at(0), metvec );
          //if( (METauto>50.) && (MTval>50.) && (tightmuons.at(0).Pt() > 20.) ){
          if( (METauto>40.) && (tightmuons.at(0).Pt() > 20.) ){
            FillHist(ThisTrigger+"_W_PFMET", METauto, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_MT", MTval, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
          }
          if( (METauto>40.) && (MTval>60.) && (MTval< 100.) && (tightmuons.at(0).Pt() > 20.) ){
            FillHist(ThisTrigger+"_W_John_PFMET", METauto, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_MT", MTval, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_leadeta", tightmuons.at(0).Eta(), this_weight, -3., 3., 60);
          }

        }


      }
    }


  }

  //==== tag jet collections
  //==== pt > 40 GeV
  //==== LeptonVeto
  std::vector<snu::KJet> jetColl_tag = GetJets("JET_HN");
  std::vector<snu::KJet> jetColl_nolepveto = GetJets("JET_NOLEPTONVETO");

  std::vector<snu::KMuon> hnloose;
  hnloose.clear();
  for(unsigned int i=0; i<hnloose_raw.size(); i++){

    if(k_isdata){
      hnloose.push_back( hnloose_raw.at(i) );
    }

    else{
      if(DijetFake){
        if( !TruthMatched( hnloose_raw.at(i)) ){
          hnloose.push_back( hnloose_raw.at(i) );
        }
      }
      else if(DijetPrompt){
        if( TruthMatched( hnloose_raw.at(i)) ){
          hnloose.push_back( hnloose_raw.at(i) );
        }
      }
      else{
        return;
      }
    }

  }

/*
  //==== 2-0) Loose Scan

  //std::vector<snu::KMuon> nocutmu_raw = GetMuons("MUON_HN_NOCUT", true);
  std::vector<snu::KMuon> nocutmu_raw = GetMuons("MUON_HN_NOCUT_IPsame", true);
  std::vector<snu::KMuon> nocutmu;
  nocutmu.clear();

  for(unsigned int i=0; i<nocutmu_raw.size(); i++){
    if( !TruthMatched( nocutmu_raw.at(i) ) ){
      nocutmu.push_back( nocutmu_raw.at(i) );
    }
  }

  for(unsigned int i=0; i<nocutmu.size(); i++){

    snu::KMuon muon = nocutmu.at(i);

    //==== find closeset jet
    double dr = 0.4;
    bool jetfound=false;
    snu::KJet cljet;
    for(unsigned int j=0; j<jetColl_nolepveto.size(); j++){

      snu::KJet jet = jetColl_nolepveto.at(j);
      if( muon.DeltaR( jet ) < dr ){
        dr = muon.DeltaR( jet );
        jetfound = true;
        cljet = jet;
      }

    }

    if(jetfound){

      int hf = cljet.HadronFlavour();
      TString flavour = "Heavy";

      if(hf<4) flavour = "Light";

      FillHist(flavour+"_Chi2", muon.GlobalChi2(), 1., 0., 50., 50);
      FillHist(flavour+"_dXY", fabs(muon.dXY()), 1., 0., 0.2, 200);
      FillHist(flavour+"_dZ", fabs(muon.dZ()), 1., 0., 0.1, 100);
      FillHist(flavour+"_dXYSig", fabs(muon.dXYSig()), 1., 0., 10., 100);
      FillHist(flavour+"_validHits", muon.validHits(), 1., 0., 50., 50);
      FillHist(flavour+"_validPixHits", muon.validPixHits(), 1., 0., 50., 50);
      FillHist(flavour+"_validStations", muon.validStations(), 1., 0., 50., 50);
      FillHist(flavour+"_ActiveLayer", muon.ActiveLayer(), 1., 0., 50., 50);

    }

  }

  //==== Make x:iso, y:mva 2D one-binned FR

  double MaxSIP = 3, dMaxSIP = 0.5;
  for(int a=0; a<60; a++){

    double MaxChi2 = 10., dMaxChi2 = 10.;
    for(int b=0; b<15; b++){

      int N_thismu(0);
      snu::KMuon muon;
      for(unsigned int i=0; i<nocutmu.size(); i++){
        if( (fabs(nocutmu.at(i).dXYSig()) < MaxSIP) && (nocutmu.at(i).GlobalChi2() < MaxChi2) ){
          N_thismu++;
          muon = nocutmu.at(i);
        }
      } 

      if(N_thismu==1){

        //==== find closeset jet
        double dr = 0.4;
        bool jetfound=false;
        snu::KJet cljet;
        for(unsigned int j=0; j<jetColl_nolepveto.size(); j++){

          snu::KJet jet = jetColl_nolepveto.at(j);
          if( muon.DeltaR( jet ) < dr ){
            dr = muon.DeltaR( jet );
            jetfound = true;
            cljet = jet;
          }

        }
        if(jetfound){

          //bool IsThisTight = PassID( muon, "MUON_HN_TIGHT" );
          bool IsThisTight = PassID( muon, "MUON_HN_TIGHTv2" );
          int hf = cljet.HadronFlavour();

          std::vector<snu::KMuon> onemu;
          onemu.push_back(muon);
          double this_weight = weight;

          double bin_MaxSIP = MaxSIP + dMaxSIP/2.;
          double bin_MaxChi2 = MaxChi2 + dMaxChi2/2.;

          double TightISO = 0.07;
          double conept = MuonConePt(muon,TightISO);

          TString SIP_Chi2 = TString::Itoa(MaxSIP*10,10)+"_"+TString::Itoa(MaxChi2,10);
          //cout << MaxSIP<<"\t"<<MaxChi2<<" ==> " << SIP_Chi2 << endl;

          if( hf >= 4 ){
            FillHist("HeavyFlavour_F0", bin_MaxSIP, bin_MaxChi2, 1., 0., 40., 80, 0., 200., 200);
            FillHist(SIP_Chi2+"_HeavyFlavour_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), 1., ptarray, n_pt, etaarray, n_eta);
            if(IsThisTight){
             FillHist("HeavyFlavour_F", bin_MaxSIP, bin_MaxChi2, 1., 0., 40., 80, 0., 200., 200);
             FillHist(SIP_Chi2+"_HeavyFlavour_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), 1., ptarray, n_pt, etaarray, n_eta);
            }
          }
          else{
            FillHist("LightFlavour_F0", bin_MaxSIP, bin_MaxChi2, 1., 0., 40., 80, 0., 200., 200);
            FillHist(SIP_Chi2+"_LightFlavour_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), 1., ptarray, n_pt, etaarray, n_eta);
            if(IsThisTight){
             FillHist("LightFlavour_F", bin_MaxSIP, bin_MaxChi2, 1., 0., 40., 80, 0., 200., 200);
             FillHist(SIP_Chi2+"_LightFlavour_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), 1, ptarray, n_pt, etaarray, n_eta);
            }
          }

        } // closest jet found

      } // only one muon

      MaxChi2 += dMaxChi2;

    } // chi2 loop

    MaxSIP += dMaxSIP;

  } // iso loop
  return;
*/

  double AwayjetPts[] = {20, 30, 40, 60};

  if( PassTriggerOR(AllHLTs) ){

    if( (hnloose_raw.size() == 1) && (hnloose.size() == 1) && (muons_veto.size() == 1) ){

      snu::KMuon muon = hnloose.at(0);

      double dr = 0.5;
      bool HasCloseBjet_Medium=false, HasCloseBjet_Loose=false;
      for(unsigned int j=0; j<jetColl_nolepveto.size(); j++){
        
        snu::KJet jet = jetColl_nolepveto.at(j);
        if( muon.DeltaR( jet ) < dr ){
          if(IsBTagged(jet, snu::KJet::CSVv2, snu::KJet::Medium)){
            HasCloseBjet_Medium = true;
          }
          if(IsBTagged(jet, snu::KJet::CSVv2, snu::KJet::Loose)){
            HasCloseBjet_Loose = true;
          }
        }
      }

      for(std::map< TString, std::vector<double> >::iterator it=HLT_ptconerange.begin(); it!=HLT_ptconerange.end(); it++){

        double weight_by_pt_cone = GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], hnloose, For_HLT_Mu3_PFJet40_v);

        TString TightID = "MUON_HN_TIGHT";
        if(DoHighdXY && !DoSmallHighdXY) TightID = "MUON_HN_Tight_HighdXY";

        //TightID = "MUON_HN_TIGHT_8TeV";

        bool IsThisTight = PassID( muon, TightID );

        if(DijetFake){

          double this_weight = weight_by_pt_cone*MCweight;

          double AwayjetPt = 40; // just using central value..

          TString HISTPREFIX = "SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10);
          if(DoHighdXY) HISTPREFIX = "SingleMuonTrigger_HighdXY";

          FillDenAndNum(HISTPREFIX, muon, this_weight, IsThisTight);
          if(HasCloseBjet_Medium) FillDenAndNum(HISTPREFIX+"_withbjet_Medium", muon, this_weight, IsThisTight);
          else                    FillDenAndNum(HISTPREFIX+"_withoutbjet_Medium", muon, this_weight, IsThisTight);
          if(HasCloseBjet_Loose) FillDenAndNum(HISTPREFIX+"_withbjet_Loose", muon, this_weight, IsThisTight);
          else                   FillDenAndNum(HISTPREFIX+"_withoutbjet_Loose", muon, this_weight, IsThisTight);

          double ptweight = JSWeightByTrigger(it->first, TargetLumi);
          ptweight *= MCweight;
          if( !PassTrigger(it->first) ) ptweight = 0.;
          if( muon.MiniAODPt() < HLT_ptmin[it->first] ) ptweight = 0.;
          if( (it->first).Contains("PFJet40_v") ){
            if( For_HLT_Mu3_PFJet40_v == 0 ) ptweight = 0.;
          }

          HISTPREFIX = (it->first)+"_SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10);
          if(DoHighdXY) HISTPREFIX = (it->first)+"_SingleMuonTrigger_HighdXY";

          double TightISO = 0.07;
          double conept = MuonConePt(muon,TightISO);

          FillHist(HISTPREFIX+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
          if(HasCloseBjet_Medium) FillHist(HISTPREFIX+"_withbjet_Medium_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
          else                    FillHist(HISTPREFIX+"_withoutbjet_Medium_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
          if(HasCloseBjet_Loose) FillHist(HISTPREFIX+"_withbjet_Loose_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
          else                    FillHist(HISTPREFIX+"_withoutbjet_Loose_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);

          if(IsThisTight){
            FillHist(HISTPREFIX+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
            if(HasCloseBjet_Medium) FillHist(HISTPREFIX+"_withbjet_Medium_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
            else                    FillHist(HISTPREFIX+"_withoutbjet_Medium_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
            if(HasCloseBjet_Loose) FillHist(HISTPREFIX+"_withbjet_Loose_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
            else                    FillHist(HISTPREFIX+"_withoutbjet_Loose_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
          } 

        }
        else{

          double this_weight = weight_by_pt_cone*weight;

          double pu_reweight = 1.;
          double trigger_additional_sf = 1.;
          double muon_scale_sf = 1.;
          if(!k_isdata){
            //pu_reweight = mcdata_correction->GetVtxReweight(it->first, eventbase->GetEvent().nVertices());

            double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(hnloose);
            muon_scale_sf *= MuTrkEffSF;

            if(IsThisTight){
/*
              if((it->first).Contains("HLT_Mu3_PFJet40_v")){
                trigger_additional_sf = 0.970059;
              }
              else if((it->first).Contains("HLT_Mu8_TrkIsoVVL_v")){
                trigger_additional_sf = 0.918045;
              }
              else if((it->first).Contains("HLT_Mu17_TrkIsoVVL_v")){
                trigger_additional_sf = 0.942808;
              }
*/
              muon_scale_sf *= mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, 0);
            }
          }

          this_weight *= pu_reweight*trigger_additional_sf*muon_scale_sf;


          TLorentzVector metvec;
          metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto );
          double MTval = AnalyzerCore::MT( muon, metvec );

          int n_awayjets = 4;
          if(DoHighdXY) n_awayjets = 1;

          for(int j=0; j<n_awayjets; j++){

            double AwayjetPt = AwayjetPts[j];

            int n_tagjets = jetColl_tag.size();
            if(DoHighdXY) n_tagjets = 1;

            bool histfilled = false; //Fill only one event at most
            for(unsigned int i=0; i<n_tagjets; i++){

              if(histfilled) break;

              bool UseEvent = false;

              if(!DoHighdXY){
                snu::KJet jet = jetColl_tag.at(i);
                if( jet.Pt() < AwayjetPt ) continue;
                double dPhi = muon.DeltaPhi( jet );
                //==== If QCD, don't have to require MET/MT
                if( DijetFake )        UseEvent = (dPhi > 2.5) && (jet.Pt()/muon.Pt() > 1.0);
                //==== If not, use it to remove W events
                else if( DijetPrompt ){
                  UseEvent = (dPhi > 2.5) && (jet.Pt()/muon.Pt() > 1.0) && (METauto < 80.) && (MTval < 25.);
/*
                  if(UseEvent && IsThisTight){
                  cout << "RunNumber = " << Evt.RunNumber() << endl;
                  cout << "EventNumber = " << Evt.EventNumber() << endl;
                  cout << "trigger = " << it->first << endl;
                  cout << "muon.Pt() = " << muon.Pt() << endl;
                  cout << "muon.RelIso04() = " << muon.RelIso04() << endl;
                  cout << "muon.dXYSig() = " << muon.dXYSig() << endl;
                  cout << "jet.Pt() = " << jet.Pt() << endl;
                  cout << "AwayjetPt = " << AwayjetPt << endl;
                  cout << "dPhi = " << dPhi << endl;
                  cout << "jet.Pt()/muon.Pt() = " << jet.Pt()/muon.Pt() << endl;
                  cout << "METauto = " << METauto << endl;
                  cout << "MTval = " << MTval << endl;
                  cout << "this_weight = " << this_weight << endl;
                  cout << "pu_reweight = " << pu_reweight << endl;
                  cout << "trigger_additional_sf = " << trigger_additional_sf << endl;
                  cout << "muon_scale_sf = " << muon_scale_sf << endl;
                  cout << "weight_by_pt_cone = " << weight_by_pt_cone << endl;
                  cout << "weight = " << weight << endl;
                  cout << "--> UseEvent = " << UseEvent << endl << endl;
                  }
*/
                }
              }
              else UseEvent = true;

              if( UseEvent ){

                TString HISTPREFIX = "SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10);
                if(DoHighdXY) HISTPREFIX = "SingleMuonTrigger_HighdXY";

                FillDenAndNum(HISTPREFIX, muon, this_weight, IsThisTight);
                if(HasCloseBjet_Medium) FillDenAndNum(HISTPREFIX+"_withbjet_Medium", muon, this_weight, IsThisTight);
                else                    FillDenAndNum(HISTPREFIX+"_withoutbjet_Medium", muon, this_weight, IsThisTight);
                if(HasCloseBjet_Loose) FillDenAndNum(HISTPREFIX+"_withbjet_Loose", muon, this_weight, IsThisTight);           
                else                   FillDenAndNum(HISTPREFIX+"_withoutbjet_Loose", muon, this_weight, IsThisTight);

                double ptweight = JSWeightByTrigger(it->first, TargetLumi);
                ptweight *= weight*pu_reweight*trigger_additional_sf*muon_scale_sf;
                if( !PassTrigger(it->first) ) ptweight = 0.;
                if( muon.MiniAODPt() < HLT_ptmin[it->first] ) ptweight = 0.;
                if( (it->first).Contains("PFJet40_v") ){
                  if( For_HLT_Mu3_PFJet40_v == 0 ) ptweight = 0.;
                }

                HISTPREFIX = (it->first)+"_SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10);
                if(DoHighdXY) HISTPREFIX = (it->first)+"_SingleMuonTrigger_HighdXY";

                double TightISO = 0.07;
                double conept = MuonConePt(muon,TightISO);

                FillHist(HISTPREFIX+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                if(HasCloseBjet_Medium) FillHist(HISTPREFIX+"_withbjet_Medium_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                else                    FillHist(HISTPREFIX+"_withoutbjet_Medium_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                if(HasCloseBjet_Loose) FillHist(HISTPREFIX+"_withbjet_Loose_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                else                    FillHist(HISTPREFIX+"_withoutbjet_Loose_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);

                if(IsThisTight){
                  FillHist(HISTPREFIX+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                  if(HasCloseBjet_Medium) FillHist(HISTPREFIX+"_withbjet_Medium_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                  else                    FillHist(HISTPREFIX+"_withoutbjet_Medium_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                  if(HasCloseBjet_Loose) FillHist(HISTPREFIX+"_withbjet_Loose_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                  else                    FillHist(HISTPREFIX+"_withoutbjet_Loose_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), ptweight, ptarray, n_pt, etaarray, n_eta);
                }

                histfilled = true;

              } // Use this event 

            } // END Tag jet loop

          } // Awayjet loop

        } // DijetPrompt

      } // Trigger loop

    } // Tag Jet and muon exist

  } // END PassTriggerOR

  if(DijetFake || DijetPrompt) return;

  //=========================================================
  //==== Large dXYSig Muon definitions for systematic study
  //=========================================================

  double dXYMin_central = 4.0;
  double RelIsoMax_central = 0.4;

  const int n_dXYMins = 3;
  double dXYMins[n_dXYMins] = {3.0, 4.0, 5.0};
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
      //==== Prepare Muon Collections
      //===============================

      //==== NodXYCut muons
      std::vector<snu::KMuon> muontriNodXYCutLooseColl_raw;
      std::vector<snu::KMuon> muontriNodXYCutTightColl_raw;
      std::vector<snu::KMuon> muontriNodXYCutLooseColl;
      std::vector<snu::KMuon> muontriNodXYCutTightColl;
      //==== Small dXYSig muons
      std::vector<snu::KMuon> muontriLooseColl_raw;
      std::vector<snu::KMuon> muontriTightColl_raw;
      std::vector<snu::KMuon> muontriLooseColl;
      std::vector<snu::KMuon> muontriTightColl;
      //==== Large dXYSig muons
      std::vector<snu::KMuon> muontriHighdXYLooseColl_raw;
      std::vector<snu::KMuon> muontriHighdXYTightColl_raw;
      std::vector<snu::KMuon> muontriHighdXYLooseColl;
      std::vector<snu::KMuon> muontriHighdXYTightColl;

      for(unsigned int j=0; j<muontriNodXYCutVLooseColl_raw.size(); j++){
        snu::KMuon this_muon = muontriNodXYCutVLooseColl_raw.at(j);
        //==== MaxRelIso for Loose
        if( this_muon.RelIso04() >= RelIsoMaxs[bbb] ) continue;

        //==== Fill NodXYCut
        muontriNodXYCutLooseColl_raw.push_back( this_muon ); 
        if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl_raw.push_back( this_muon );
        //==== Fill Small dXYSig
        if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
          muontriLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriTightColl_raw.push_back( this_muon );
        }
        //==== Fill Large dXYSig
        if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
          muontriHighdXYLooseColl_raw.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl_raw.push_back( this_muon );
        }

        //==== for data, or QCD
        //==== do the same for prompt collection
        if( k_isdata || k_sample_name.Contains("QCD") ){
          //==== Fill NodXYCut
          muontriNodXYCutLooseColl.push_back( this_muon );
          if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl.push_back( this_muon );
          //==== Fill Small dXYSig
          if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
            muontriLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriTightColl.push_back( this_muon );
          }
          //==== Fill Large dXYSig
          if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
            muontriHighdXYLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl.push_back( this_muon );
          }
        }

        //==== for MC,
        //==== fill only MCMatched() muons
        else{
          if( TruthMatched(this_muon) ){
            //==== Fill NodXYCut
            muontriNodXYCutLooseColl.push_back( this_muon );
            if( this_muon.RelIso04() < TightISO ) muontriNodXYCutTightColl.push_back( this_muon );
            //==== Fill Small dXYSig
            if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
              muontriLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < TightISO ) muontriTightColl.push_back( this_muon );
            }
            //==== Fill Large dXYSig
            if( fabs(this_muon.dXY()) < 1. && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              muontriHighdXYLooseColl.push_back( this_muon );
              if( this_muon.RelIso04() < TightISO ) muontriHighdXYTightColl.push_back( this_muon );
            }
          }
        }


      } //==== Loosest muon loop

      //===================
      //==== dXYCut study
      //===================

      if( k_sample_name.Contains("DY") || k_sample_name.Contains("HN") ){

        for(unsigned int i=0; i<muontriNodXYCutLooseColl.size(); i++){
          snu::KMuon this_muon = muontriNodXYCutLooseColl.at(i);
          if(fabs(this_muon.dXY()) < 1.){
            FillHist(str_dXYCut+"_prompt_Loose_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_prompt_Loose_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            if(this_muon.RelIso04()<TightISO){
              FillHist(str_dXYCut+"_prompt_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_prompt_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }
      if( k_sample_name.Contains("TTJets_aMC") || k_sample_name.Contains("QCD") ){

        for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){
          snu::KMuon this_muon = muontriNodXYCutLooseColl_raw.at(i);
          if( fabs(this_muon.dXY()) < 1. && !TruthMatched(this_muon) ){
            FillHist(str_dXYCut+"_fake_Loose_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
            FillHist(str_dXYCut+"_fake_Loose_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            if(this_muon.RelIso04()<TightISO){
              FillHist(str_dXYCut+"_fake_Tight_dXYSig", fabs(this_muon.dXYSig()), 1., 0., 15., 150);
              FillHist(str_dXYCut+"_fake_Tight_dXY", fabs(this_muon.dXY()), 1., 0., 0.1, 100);
            }
          }
        }

      }

      //=========================
      //==== SingleMuon Trigger
      //=========================

      if( PassTriggerOR(AllHLTs) ){

        std::map<TString, double> this_weight_Loose, this_weight_HighdXYLoose, this_weight_NodXYCutLoose;
        for(std::map< TString, std::vector<double> >::iterator it=HLT_ptconerange.begin(); it!=HLT_ptconerange.end(); it++){
          this_weight_Loose[it->first]         = weight*GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], muontriLooseColl, For_HLT_Mu3_PFJet40_v);
          this_weight_HighdXYLoose[it->first]  = weight*GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], muontriHighdXYLooseColl, For_HLT_Mu3_PFJet40_v);
          this_weight_NodXYCutLoose[it->first] = weight*GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], muontriNodXYCutLooseColl, For_HLT_Mu3_PFJet40_v);
        }

        //Double_t this_weight_Loose = weight*GetPrescale(muontriLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_HighdXYLoose = weight*GetPrescale(muontriHighdXYLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));
        //Double_t this_weight_NodXYCutLoose = weight*GetPrescale(muontriNodXYCutLooseColl, PassTrigger("HLT_Mu8_v"), PassTrigger("HLT_Mu17_v"));

        //==== 1) LooseMuon study

        if( muontriLooseColl.size() == 1 ){
          snu::KMuon muon = muontriLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_RelIso", LeptonRelIso, this_weight_Loose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_Loose, 0., 50., 50);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dXY", fabs(muon.dXY()), this_weight_Loose, 0., 0.1, 100);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_Loose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_dZ", fabs(muon.dZ()), this_weight_Loose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_LooseMuon_onebin", 0., this_weight_Loose, 0., 1., 1);
        }

        //==== 2) SingleMuon HighdXY LooseMuon study

        if( muontriHighdXYLooseColl.size() == 1 ){
          snu::KMuon muon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_RelIso", LeptonRelIso, this_weight_HighdXYLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_HighdXYLoose, 0, 50., 50);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_HighdXYLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_HighdXYLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_HighdXYLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_HighdXY_LooseMuon_onebin", 0., this_weight_HighdXYLoose, 0., 1., 1);
        }

        //==== 2') SingleMuon NodXY Loose Muon study
        if( muontriNodXYCutLooseColl.size() == 1 ){
          snu::KMuon muon = muontriNodXYCutLooseColl.at(0);
          double LeptonRelIso = muon.RelIso04();
          double conept = MuonConePt(muon,TightISO);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_eta", muon.Eta(), this_weight_Loose, -3., 3., 30);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_pt", muon.Pt(), this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_pt_cone", conept, this_weight_Loose, 0., 200., 200);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_RelIso", LeptonRelIso, this_weight_NodXYCutLoose, 0., 1., 100);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_Chi2", muon.GlobalChi2(), this_weight_NodXYCutLoose, 0, 50., 50);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dXY", fabs(muon.dXY()), this_weight_NodXYCutLoose, 0., 1., 1000);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dXYSig", fabs(muon.dXYSig()), this_weight_NodXYCutLoose, 0., 15., 150);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_dZ", fabs(muon.dZ()), this_weight_NodXYCutLoose, 0., 0.5, 50);
          FillHistByTrigger(str_dXYCut+"_NodXY_LooseMuon_onebin", 0., this_weight_NodXYCutLoose, 0., 1., 1);
        }


        //==== 3) large dXY muon method

        if(muontriHighdXYLooseColl.size()==1){
          snu::KMuon HighdXYmuon = muontriHighdXYLooseColl.at(0);
          double LeptonRelIso = HighdXYmuon.RelIso04();

          //==== half sample test
          if(dXYMins[aaa]== dXYMin_central && RelIsoMaxs[bbb]==RelIsoMax_central){

            int EventNumber = eventbase->GetEvent().EventNumber();

            double pt_for_FR = HighdXYmuon.Pt();
            double eta_for_FR = fabs(HighdXYmuon.Eta());
            if(pt_for_FR>=60.) pt_for_FR=59.;
            if(eta_for_FR>=2.5) eta_for_FR=2.2;

            snu::KEvent Evt = eventbase->GetEvent();
            double MET = Evt.MET();

            //==== sample A
            //==== obtain FR_A(pt, eta)
            if(abs(EventNumber%2==0)){
              FillHist("TEST_HalfSample", 0., 1., 0., 2., 2);
              TString HalfSampleIndex = "SampleA";
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F0", 
                                HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, n_pt, etaarray, n_eta);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, n_pt, etaarray, n_eta);
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
                                HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, n_pt, etaarray, n_eta);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_PFMET_F0", MET, this_weight_HighdXYLoose, 0., 500., 500);
              FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_njets_F0", n_jets, this_weight_HighdXYLoose, 0., 10., 10);
              if( LeptonRelIso < TightISO ){
                FillHistByTrigger(str_dXYCut+"_HighdXY_HalfSample_"+HalfSampleIndex+"_events_F",
                                  HighdXYmuon.Pt(), fabs(HighdXYmuon.Eta()), this_weight_HighdXYLoose, ptarray, n_pt, etaarray, n_eta);
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
          FillDenAndNum(str_dXYCut+"_HighdXY_alljet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);

          //==== no jet
          if( n_jets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0jet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }

          //==== no b-jet
          if( n_bjets == 0){
            FillDenAndNum(str_dXYCut+"_HighdXY_0bjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }
          //==== with b-jet
          else{
            FillDenAndNum(str_dXYCut+"_HighdXY_withbjet", HighdXYmuon, this_weight_HighdXYLoose, isThisTight);
          }


        } // END one HighdXYSig Muon 

        //==== 5) MC Truth
        //==== here, we pick FAKE muons using GEN info.
        if( !k_isdata ){

          for(unsigned int i=0; i<muontriNodXYCutLooseColl_raw.size(); i++){

            snu::KMuon muon = muontriNodXYCutLooseColl_raw.at(i);

            //==== if prompt, skip
            if( TruthMatched(muon) ) continue;

            double LeptonRelIso = muon.RelIso04();
            bool isThisTight = (LeptonRelIso < TightISO);

            //=================================================
            //==== 1) |dXY| < 1 cm, |dXY/err| > dXYMins[aaa].
            //=================================================

            if( fabs( muon.dXY() ) < 1. &&
                fabs( muon.dXYSig() ) > dXYMins[aaa] ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_alljet", muon, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0jet", muon, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withjet", muon, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_0bjet", muon, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_HighdXY_withbjet", muon, 1., isThisTight);
              }

            } // small dXYSig region


            //=========================================
            //==== 2) |dXY| < thisdXYCut cm, |dXY/err| < 3.
            //=========================================

            if( fabs( muon.dXY() ) < thisdXYCut &&
                fabs( muon.dXYSig() ) < 3. ){

              //==== all jet
              FillDenAndNum(str_dXYCut+"_MCTruth_alljet", muon, 1., isThisTight);

              //==== no jet
              if(n_jets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0jet", muon, 1., isThisTight);
              }
              //==== with jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withjet", muon, 1., isThisTight);
              }

              //==== no b-jet
              if(n_bjets==0){
                FillDenAndNum(str_dXYCut+"_MCTruth_0bjet", muon, 1., isThisTight);
              }
              //==== with b-jet
              else{
                FillDenAndNum(str_dXYCut+"_MCTruth_withbjet", muon, 1., isThisTight);
              }

            } // Large dXYSig Region


          } // END Muon loop
        } // END k_isdata


      } // SingleMuon trigger fired


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

        int n_nodxycut_loose = muontriNodXYCutLooseColl.size();
        if( n_nodxycut_loose >= 3){

          //==== Find Z-pair
          double m_OS_Zclosest = (91.2+15.0);
          int z_pair_index[2];
          bool isZTagged = false;

          for(unsigned i=0; i<n_nodxycut_loose-1; i++){
            for(unsigned j=i+1; j<n_nodxycut_loose; j++){
              if( muontriNodXYCutLooseColl.at(i).Charge() == muontriNodXYCutLooseColl.at(j).Charge() ) continue;
          bool i_muon_tight =  fabs(muontriNodXYCutLooseColl.at(i).dXY())<thisdXYCut && fabs(muontriNodXYCutLooseColl.at(i).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(i).RelIso04()<TightISO;
          bool j_muon_tight =  fabs(muontriNodXYCutLooseColl.at(j).dXY())<thisdXYCut && fabs(muontriNodXYCutLooseColl.at(j).dXYSig())<3.0 && muontriNodXYCutLooseColl.at(j).RelIso04()<TightISO;

              double m_thisOS = ( muontriNodXYCutLooseColl.at(i)+muontriNodXYCutLooseColl.at(j) ).M();
              if(  i_muon_tight &&  j_muon_tight && fabs(m_thisOS-91.2) < fabs(m_OS_Zclosest-91.2) ){
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

            snu::KMuon this_muon =  muontriNodXYCutLooseColl.at(i);
            double conept = MuonConePt(this_muon,TightISO);

            //==== Large dXYSig
            if( fabs(this_muon.dXY()) < 1.0 && fabs(this_muon.dXYSig()) > dXYMins[aaa] ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_vs_eta_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_cone_vs_eta_F0", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);


              if( this_muon.RelIso04() < TightISO ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_vs_eta_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Large_events_pt_cone_vs_eta_F", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);

              }
            }

            //==== Small dXYSig
            if( fabs(this_muon.dXY()) < thisdXYCut && fabs(this_muon.dXYSig()) < 3.0 ){
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 0, this_weight, 0., 2., 2);

              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F0", this_muon.Eta(), this_weight, -3, 3, 30);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F0", this_muon.Pt(), this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_cone_F0", conept, this_weight, 0., 200., 200);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F0", this_muon.RelIso04(), this_weight, 0., 1., 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F0", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F0", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F0", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F0", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F0", 0., this_weight, 0, 1., 1);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_vs_eta_F0", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_cone_vs_eta_F0", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              if( this_muon.RelIso04() < TightISO ){
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig", 1, this_weight, 0., 2., 2);

                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_eta_F", this_muon.Eta(), this_weight, -3, 3, 30);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_F", this_muon.Pt(), this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_pt_cone_F", conept, this_weight, 0., 200., 200);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_RelIso_F", this_muon.RelIso04(), this_weight, 0., 1., 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_Chi2_F", this_muon.GlobalChi2(), this_weight, 0, 50., 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXY_F", fabs(this_muon.dXY()), this_weight, 0., 0.1, 100);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dXYSig_F", fabs(this_muon.dXYSig()), this_weight, 0., 15., 150);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_dZ_F", fabs(this_muon.dZ()), this_weight, 0, 0.5, 50);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_onebin_F", 0., this_weight, 0, 1., 1);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_vs_eta_F", this_muon.Pt(), fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
                FillHist(str_dXYCut+"_DiMuonTrigger_ZTag_Small_events_pt_cone_vs_eta_F", conept, fabs(this_muon.Eta()), this_weight, ptarray_2, 3, etaarray_2, 2);
              }
            }

          }

        } //== END at least three muons


      } //==== DiMuon Trigger



    } //==== RelIso loop
  } //==== dXYMin loop


   return;
}// End of execute event loop
  


void FRCalculator_Mu_dxysig_DILEP::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FRCalculator_Mu_dxysig_DILEP::BeginCycle() throw( LQError ){
  
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

FRCalculator_Mu_dxysig_DILEP::~FRCalculator_Mu_dxysig_DILEP() {
  
  Message("In FRCalculator_Mu_dxysig_DILEP Destructor" , INFO);
  
}


void FRCalculator_Mu_dxysig_DILEP::FillCutFlow(TString cut, float weight){

  
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


void FRCalculator_Mu_dxysig_DILEP::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FRCalculator_Mu_dxysig_DILEP::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FRCalculator_Mu_dxysig_DILEPCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FRCalculator_Mu_dxysig_DILEP::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

float FRCalculator_Mu_dxysig_DILEP::GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh){
  float prescale_trigger = 0.;

  if(muon.size() == 1){

    if(muon.at(0).Pt() >= 20.){
      if(passhigh){
        prescale_trigger = JSWeightByTrigger("HLT_Mu17_v", TargetLumi) ; //// 20 + GeV bins
      }
      else prescale_trigger = 0.;
    }
    else{
      if(passlow){
        prescale_trigger = JSWeightByTrigger("HLT_Mu8_v", TargetLumi) ;
      }
      else prescale_trigger = 0.;
    }
  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;
}

double FRCalculator_Mu_dxysig_DILEP::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KMuon> muons, int npfjet50){

  double prescale_trigger = 0.;

  if(muons.size()==1){
    snu::KMuon muon = muons.at(0);
    double min_cone_pt = ptrange.at(0);
    double max_cone_pt = ptrange.at(1);

    bool SafePt = (muon.MiniAODPt() > trigger_safe_pt);

    if(SafePt && PassTrigger(hltname)){

      double TightISO = 0.07;
      double conept = MuonConePt(muon,TightISO);

      if(conept >= min_cone_pt && conept < max_cone_pt){
        prescale_trigger = JSWeightByTrigger(hltname, TargetLumi) ;
      }

      if(hltname.Contains("PFJet40")){
        if(npfjet50==0) prescale_trigger = 0.;
      }
    }

  }

  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}

void FRCalculator_Mu_dxysig_DILEP::FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value, it->second, xmin, xmax, nbins);
  }

}
void FRCalculator_Mu_dxysig_DILEP::FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){
    FillHist(it->first+"_"+histname, value1, value2, it->second, x, nbinsx, y, nbinsy);
  }

}

void FRCalculator_Mu_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight){

  const int n_eta = 3;
  float etaarray[n_eta+1] = {0.0, 0.8, 1.479, 2.5};

  //==== Data Binning
  const int n_pt = 8;
  float ptarray[n_pt+1] = {5., 10., 15., 20., 30., 40., 50., 60., 70.};

/*
  //==== MC Binning
  const int n_pt = 9;
  float ptarray[n_pt+1] = {5., 10., 15., 20., 25., 30., 40., 50., 60., 70.};
*/

  double TightISO = 0.07;
  double conept = MuonConePt(muon,TightISO);

  TLorentzVector METvec;
  METvec.SetPtEtaPhiE(METauto, 0, METphiauto, METauto);
  double this_mt = MT(muon, METvec);

  FillHist(prefix+"_eta_F0", muon.Eta(), thisweight, -3., 3., 30);
  FillHist(prefix+"_pt_F0", muon.Pt(), thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_F0", conept, thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_FRbinned_F0", conept, thisweight, ptarray, n_pt);
  FillHist(prefix+"_pt_cone_FRbinned_w1_F0", conept, 1., ptarray, n_pt);
  FillHist(prefix+"_RelIso_F0", muon.RelIso04(), thisweight, 0., 1., 100);
  FillHist(prefix+"_Chi2_F0", muon.GlobalChi2(), thisweight, 0., 200., 200);
  FillHist(prefix+"_dXY_F0", fabs(muon.dXY()), thisweight, 0., 1., 1000);
  FillHist(prefix+"_dXYSig_F0", fabs(muon.dXYSig()), thisweight, 0., 40., 400);
  FillHist(prefix+"_dZ_F0", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
  FillHist(prefix+"_Type_F0", muon.GetType(), thisweight, 0., 50., 50);
  FillHist(prefix+"_onebin_F0", 0., thisweight, 0., 1., 1);
  FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
  FillHist(prefix+"_PFMET_F0", METauto, thisweight, 0., 1000., 1000);
  FillHist(prefix+"_MT_F0", this_mt, thisweight, 0., 1000., 1000);

  if( isTight ){
    FillHist(prefix+"_eta_F", muon.Eta(), thisweight, -3., 3., 30);
    FillHist(prefix+"_pt_F", muon.Pt(), thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_F", conept, thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_FRbinned_F", conept, thisweight, ptarray, n_pt);
    FillHist(prefix+"_pt_cone_FRbinned_w1_F", conept, 1., ptarray, n_pt);
    FillHist(prefix+"_RelIso_F", muon.RelIso04(), thisweight, 0., 1., 100);
    FillHist(prefix+"_Chi2_F", muon.GlobalChi2(), thisweight, 0., 200., 200);
    FillHist(prefix+"_dXY_F", fabs(muon.dXY()), thisweight, 0., 1., 1000);
    FillHist(prefix+"_dXYSig_F", fabs(muon.dXYSig()), thisweight, 0., 40., 400);
    FillHist(prefix+"_dZ_F", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
    FillHist(prefix+"_Type_F", muon.GetType(), thisweight, 0., 50., 50);
    FillHist(prefix+"_onebin_F", 0., thisweight, 0., 1., 1);
    FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);
    FillHist(prefix+"_MT_F", this_mt, thisweight, 0., 1000., 1000);

  }


}

void FRCalculator_Mu_dxysig_DILEP::FillDenAndNum(TString prefix, snu::KMuon muon, std::map<TString, double> hltweight, bool isTight){

  for(std::map<TString, double>::iterator it=hltweight.begin(); it!=hltweight.end(); it++){

    FillDenAndNum(it->first+"_"+prefix, muon, it->second, isTight);

  }

}

float FRCalculator_Mu_dxysig_DILEP::JSWeightByTrigger(TString triggername, float tlumi){

  if(isData) return 1.;

  double NormSF(1.);
/*
  if(triggername=="HLT_Mu3_PFJet40_v"){
    NormSF = 0.72799;
  }
  if(triggername=="HLT_Mu8_TrkIsoVVL_v"){
    NormSF = 1.39876;
  }
*/
  return NormSF*WeightByTrigger(triggername, tlumi);


}
















