// $Id: trilepton_mumumu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu);

trilepton_mumumu::trilepton_mumumu() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu");
  
  Message("In trilepton_mumumu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  m_HNgenmatch->SetDrawHist(true);

  MakeCleverHistograms(hntrilephist, "cut0");
  MakeCleverHistograms(hntrilephist, "cutWlow");
  MakeCleverHistograms(hntrilephist, "cutWhigh");

  int signal_masses[] = {5, 10, 20, 30, 40, 50, 60, 70, 90, 100, 150, 200, 300, 400, 500, 700, 1000};
  for(int i=0; i<17; i++){
    TString thiscut = "cutHN"+TString::Itoa(signal_masses[i],10);
    MakeCleverHistograms(hntrilephist, thiscut);
  }

  MakeCleverHistograms(hntrilephist, "MuMuE");



}


void trilepton_mumumu::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_mumumu::ExecuteEvents()throw( LQError ){

/*
  //==== Mu3 test

  if(PassTrigger("HLT_Mu3_PFJet40_v")){
    cout << "HLT_Mu3_PFJet40_v fired" << endl;
    std::vector<snu::KMuon> testmuon = GetMuons("MUON_PTETA");
    std::vector<snu::KJet> testjet = GetJets("JET_HN");
    for(unsigned int i=0; i<testmuon.size(); i++){
      FillHist("Mu3Test_muon_pt", testmuon.at(i).Pt(), 1., 0., 50., 50);
    }
    for(unsigned int i=0; i<testjet.size(); i++){
      FillHist("Mu3Test_jet_pt", testjet.at(i).Pt(), 1., 0., 200., 200);
    }

    std::vector<snu::KJet> testjet40 = GetJets("JET_HN",40);
    FillHist("Mu3Test_jet40_size", testjet40.size(), 1, 0., 10., 10);
  }
  return;
*/

/*

  //==== Electron Charge flip

  float etaarray [] = {0.0, 0.9, 1.4442, 1.556, 2.5};
  float ptarray [] = {20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 450., 500.};
  std::vector<snu::KElectron> testel = GetElectrons(true, true, "ELECTRON_HN_TIGHT");
  for(unsigned int i=0; i<testel.size(); i++){
    snu::KElectron el = testel.at(i);
    FillHist("TEST_Electron_CF_F0", fabs(el.Eta()), el.Pt(), 1., etaarray, 4, ptarray, 25);
    if(MCIsCF(el)) FillHist("TEST_Electron_CF_F", fabs(el.Eta()), el.Pt(), 1., etaarray, 4, ptarray, 25);
  }
  return;
*/

/*

  //==== Muon from tau prompt check

  std::vector<snu::KMuon> testmuon = GetHNTriMuonsByLooseRelIso(0.4, true);
  for(unsigned int i=0;i<testmuon.size();i++){
    snu::KMuon muon = testmuon.at(i);
    if(!muon.MCMatched() && muon.MCFromTau()){
      FillHist("TEST_FromTauNotMCMatched", 0., 1., 0., 2., 2);
      if(muon.RelIso04()<0.1){
        FillHist("TEST_FromTauNotMCMatched", 1., 1., 0., 2., 2);
      }

      int mpdgid = muon.MotherPdgId();
      if(mpdgid<50){
        FillHist("TEST_FromTauNotMCMatched_mpdgidless50", 0., 1., 0., 2., 2);
        if(muon.RelIso04()<0.1){
          FillHist("TEST_FromTauNotMCMatched_mpdgidless50", 1., 1., 0., 2., 2);
        }
      }
      else{
        FillHist("TEST_FromTauNotMCMatched_mpdgidgreater50", 0., 1., 0., 2., 2);
        if(muon.RelIso04()<0.1){
          FillHist("TEST_FromTauNotMCMatched_mpdgidgreater50", 1., 1., 0., 2., 2);
        }
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
  FillHist("cutflow_MuMuE", 0., 1., 0., 10., 10);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  //====================================
  //==== For signal, run HNGenMatching
  //====================================

  if( k_sample_name.Contains("HN") ){

    std::vector<snu::KTruth> truthColl;
    eventbase->GetTruthSel()->Selection(truthColl);
    m_HNgenmatch->SetAllGenParticles(truthColl);
    m_HNgenmatch->SetSignalMass(GetSignalMass());
    m_HNgenmatch->SetHNpdgids(9900012);
    m_HNgenmatch->FindGenParticles();

    //==== check OS same-flavour mass
    bool isOSSF4GeV=false;
    if( m_HNgenmatch->gen_l_1.PdgId()+m_HNgenmatch->gen_l_2.PdgId() == 0 ){
      if( (m_HNgenmatch->gen_l_1+m_HNgenmatch->gen_l_2).M() < 4 ) isOSSF4GeV = true;
    }
    if( m_HNgenmatch->gen_l_1.PdgId()+m_HNgenmatch->gen_l_3.PdgId() == 0 ){
      if( (m_HNgenmatch->gen_l_1+m_HNgenmatch->gen_l_3).M() < 4 ) isOSSF4GeV = true;
    }
    if( m_HNgenmatch->gen_l_2.PdgId()+m_HNgenmatch->gen_l_3.PdgId() == 0 ){
      if( (m_HNgenmatch->gen_l_2+m_HNgenmatch->gen_l_3).M() < 4 ) isOSSF4GeV = true;
    }
    FillHist("GEN_isOSSF4GeV", isOSSF4GeV, 1., 0., 2., 2); 

    //==== snu::KTruth gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;
    //==== gen_l_1 and gen_l_2 are always same sign

    //==== mu mu mu
    if( fabs(m_HNgenmatch->gen_l_1.PdgId())==13 &&
        fabs(m_HNgenmatch->gen_l_2.PdgId())==13 &&
        fabs(m_HNgenmatch->gen_l_3.PdgId())==13 ){

      FillHist("SignalLeptonFlavour", 1., 1., 0., 40., 40);
      FillHist("MuMuMu_nu_pdgid", m_HNgenmatch->gen_nu.PdgId(), 1., -20., 20., 40);

    }
    //==== e mu mu
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==13 ){

      //==== SS
      if( m_HNgenmatch->gen_l_2.PdgId() == m_HNgenmatch->gen_l_3.PdgId() ){
        FillHist("SignalLeptonFlavour", 2., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 3., 1., 0., 40., 40);
      }

    }
    //==== mu e mu
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==13 ){

      //==== SS
      if( m_HNgenmatch->gen_l_1.PdgId() == m_HNgenmatch->gen_l_3.PdgId() ){
        FillHist("SignalLeptonFlavour", 4., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 5., 1., 0., 40., 40);
      }

    }
    //==== mu mu e
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==11 ){

      //==== SS
      if( m_HNgenmatch->gen_l_1.PdgId() == m_HNgenmatch->gen_l_2.PdgId() ){
        FillHist("SignalLeptonFlavour", 6., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 7., 1., 0., 40., 40);
      }

    }
    //==== e e mu
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==13 ){

      //==== SS
      if( m_HNgenmatch->gen_l_1.PdgId() == m_HNgenmatch->gen_l_2.PdgId() ){
        FillHist("SignalLeptonFlavour", 8., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 9., 1., 0., 40., 40);
      }

    }
    //==== e mu e
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==11 ){

      //==== SS
      if( m_HNgenmatch->gen_l_1.PdgId() == m_HNgenmatch->gen_l_3.PdgId() ){
        FillHist("SignalLeptonFlavour", 10., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 11., 1., 0., 40., 40);
      }

    }
    //==== mu e e
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==13 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==11 ){

      //==== SS
      if( m_HNgenmatch->gen_l_2.PdgId() == m_HNgenmatch->gen_l_3.PdgId() ){
        FillHist("SignalLeptonFlavour", 12., 1., 0., 40., 40);
      }
      //==== OS
      else{
        FillHist("SignalLeptonFlavour", 13., 1., 0., 40., 40);
      }

    }
    //==== e e e
    else if( fabs(m_HNgenmatch->gen_l_1.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_2.PdgId())==11 &&
             fabs(m_HNgenmatch->gen_l_3.PdgId())==11 ){

      FillHist("SignalLeptonFlavour", 14., 1., 0., 40., 40);

    }
    else{
      FillHist("SignalLeptonFlavour", 0., 1., 0., 40., 40);
        //cout << m_HNgenmatch->gen_l_1.PdgId() << "\t" << m_HNgenmatch->gen_l_1.Pt() << "\t" << m_HNgenmatch->gen_l_1.Eta() << endl;
        //cout << m_HNgenmatch->gen_l_2.PdgId() << "\t" << m_HNgenmatch->gen_l_2.Pt() << "\t" << m_HNgenmatch->gen_l_2.Eta() << endl;
        //cout << m_HNgenmatch->gen_l_3.PdgId() << "\t" << m_HNgenmatch->gen_l_3.Pt() << "\t" << m_HNgenmatch->gen_l_3.Eta() << endl;
    }

    if(m_HNgenmatch->allgenfound){
      FillHist("GenFound", 1., 1., 0., 2., 2);
    }
    else{
      FillHist("GenFound", 0., 1., 0., 2., 2);
    }

  }

  //======================
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  FillHist("cutflow_MuMuE", 1., 1., 0., 10., 10);

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
    FillHist("cutflow_MuMuE", 2., 1., 0., 10., 10);
    m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  }

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  //float weight_trigger_sf = TriggerScaleFactor(electrontriLooseColl, muonTightColl, "HLT_IsoMu20");

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
  FillHist("cutflow_MuMuE", 3., 1., 0., 10., 10);

  //======================
  //==== Prepare Leptons
  //======================

  double this_RelIso = 0.4;

  std::vector<snu::KMuon> muontriLooseColl;
  std::vector<snu::KElectron> electrontriLooseColl;

  bool diboson_had = std::find(k_flags.begin(), k_flags.end(), "diboson_had") != k_flags.end();
  //==== signal
  if( k_sample_name.Contains("HN") ){

    snu::KTruth gen_l_1 = m_HNgenmatch->gen_l_1;
    snu::KTruth gen_l_2 = m_HNgenmatch->gen_l_2;
    snu::KTruth gen_l_3 = m_HNgenmatch->gen_l_3;

    std::vector<snu::KMuon> muontriLooseColl_raw = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
    std::vector<snu::KElectron> electrontriLooseColl_raw = GetElectrons(false, true, "ELECTRON_HN_LOWDXY_FAKELOOSE");

    std::vector<KLepton> leptontriLooseColl_raw;
    for(unsigned int i=0; i<muontriLooseColl_raw.size(); i++) leptontriLooseColl_raw.push_back( muontriLooseColl_raw.at(i) );
    for(unsigned int i=0; i<electrontriLooseColl_raw.size(); i++) leptontriLooseColl_raw.push_back( electrontriLooseColl_raw.at(i) );

    //==== find gen_l_1
    //cout << "[gen_l_1] : pt = " << gen_l_1.Pt() << ", eta = " << gen_l_1.Eta() << endl;
    //cout << "[gen_l_2] : pt = " << gen_l_2.Pt() << ", eta = " << gen_l_2.Eta() << endl;
    //cout << "[gen_l_3] : pt = " << gen_l_3.Pt() << ", eta = " << gen_l_3.Eta() << endl;

    std::vector<int> loose_used;
    loose_used.clear();
    int loose_l_1_index = find_genmatching(gen_l_1, leptontriLooseColl_raw, loose_used);
    int loose_l_2_index = find_genmatching(gen_l_2, leptontriLooseColl_raw, loose_used);
    int loose_l_3_index = find_genmatching(gen_l_3, leptontriLooseColl_raw, loose_used);

    std::vector<KLepton> leptontriLooseColl_genorder;
    if(loose_l_1_index!=-1){
      leptontriLooseColl_genorder.push_back( leptontriLooseColl_raw.at(loose_l_1_index) );

      if(leptontriLooseColl_raw.at(loose_l_1_index).LeptonFlavour() == KLepton::MUON){
        FillHist("TEST_gen_l_1_dXY", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXY()), 1., 0., 0.1, 1000);
        FillHist("TEST_gen_l_1_dXYSig", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXYSig()), 1., 0., 10., 100);
        FillHist("TEST_gen_l_1_dZ", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dZ()), 1., 0., 0.5, 500);
        if(fabs(leptontriLooseColl_raw.at(loose_l_1_index).Eta()) < 1.4442){
          FillHist("TEST_gen_l_1_dXY_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_1_dXYSig_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_1_dZ_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dZ()), 1., 0., 0.5, 500);
        }
        else{
          FillHist("TEST_gen_l_1_dXY_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_1_dXYSig_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_1_dZ_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_1_index).dZ()), 1., 0., 0.5, 500);
        }
      }

    }
    if(loose_l_2_index!=-1){
      leptontriLooseColl_genorder.push_back( leptontriLooseColl_raw.at(loose_l_2_index) );

      if(leptontriLooseColl_raw.at(loose_l_2_index).LeptonFlavour() == KLepton::MUON){
        FillHist("TEST_gen_l_2_dXY", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXY()), 1., 0., 0.1, 1000);
        FillHist("TEST_gen_l_2_dXYSig", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXYSig()), 1., 0., 10., 100);
        FillHist("TEST_gen_l_2_dZ", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dZ()), 1., 0., 0.5, 500);
        if(fabs(leptontriLooseColl_raw.at(loose_l_2_index).Eta()) < 1.4442){
          FillHist("TEST_gen_l_2_dXY_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_2_dXYSig_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_2_dZ_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dZ()), 1., 0., 0.5, 500);
        }
        else{
          FillHist("TEST_gen_l_2_dXY_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_2_dXYSig_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_2_dZ_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_2_index).dZ()), 1., 0., 0.5, 500);
        }
      }

    }
    if(loose_l_3_index!=-1){
      leptontriLooseColl_genorder.push_back( leptontriLooseColl_raw.at(loose_l_3_index) );

      if(leptontriLooseColl_raw.at(loose_l_3_index).LeptonFlavour() == KLepton::MUON){
        FillHist("TEST_gen_l_3_dXY", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXY()), 1., 0., 0.1, 1000);
        FillHist("TEST_gen_l_3_dXYSig", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXYSig()), 1., 0., 10., 100);
        FillHist("TEST_gen_l_3_dZ", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dZ()), 1., 0., 0.5, 500);
        if(fabs(leptontriLooseColl_raw.at(loose_l_3_index).Eta()) < 1.4442){
          FillHist("TEST_gen_l_3_dXY_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_3_dXYSig_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_3_dZ_Barrel", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dZ()), 1., 0., 0.5, 500);
        }
        else{
          FillHist("TEST_gen_l_3_dXY_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXY()), 1., 0., 0.1, 1000);
          FillHist("TEST_gen_l_3_dXYSig_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dXYSig()), 1., 0., 10., 100);
          FillHist("TEST_gen_l_3_dZ_Endcap", fabs(leptontriLooseColl_raw.at(loose_l_3_index).dZ()), 1., 0., 0.5, 500);
        }
      }

    }

    //==== now sort leptontriLooseColl_genorder to ptorder, and replace
    leptontriLooseColl_genorder = sort_leptons_ptorder( leptontriLooseColl_genorder );
    for(unsigned int i=0; i<leptontriLooseColl_genorder.size(); i++){
      if(leptontriLooseColl_genorder.at(i).LeptonFlavour() == KLepton::MUON) muontriLooseColl.push_back( *leptontriLooseColl_genorder.at(i).GetMuonPtr() );
      else if(leptontriLooseColl_genorder.at(i).LeptonFlavour() == KLepton::ELECTRON) electrontriLooseColl.push_back( *leptontriLooseColl_genorder.at(i).GetElectronPtr() );
      else{
        cout << "LeptonFlavour is wrong" << endl;
      }
    }

  }
  //==== non-prompt : keep fake
  else if( k_sample_name.Contains("DY") || k_sample_name.Contains("WJets") || k_sample_name.Contains("TTJets") || k_sample_name.Contains("QCD") || k_sample_name.Contains("TTLL_powheg") ){

    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
    electrontriLooseColl = GetElectrons(false, true, "ELECTRON_HN_LOWDXY_FAKELOOSE");

/*
    //==== below is for MC fake test
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      snu::KMuon this_muon = muontriLooseColl.at(i);

      if(this_muon.MCIsFromConversion()){
        FillHist("TEST_Muon_FR_Conversion", 0., 1., 0., 2., 2);
        if(this_muon.RelIso04()<0.1){
          FillHist("TEST_Muon_FR_Conversion", 1., 1., 0., 2., 2);
        }
      }

      if(this_muon.MCFromTau()){
        FillHist("TEST_Muon_FR_Tau", 0., 1., 0., 2., 2);
        if(this_muon.RelIso04()<0.1){
          FillHist("TEST_Muon_FR_Tau", 1., 1., 0., 2., 2);
        }
      }

      if(!this_muon.MCMatched() && !this_muon.MCIsFromConversion() && !this_muon.MCFromTau()){

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
          if(jetColl_hn_nearby.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)){
            n_bjets_nearby++;
          }
        }

        for(unsigned int j=0; j<jetColl_hn_nearby.size(); j++){
          FillHist("dRNearByJetFakeMuon", jetColl_hn_nearby.at(j).DeltaR(this_muon), weight, 0., 6., 60);
          if(jetColl_hn_nearby.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)){
            FillHist("dRNearByBJetFakeMuon", jetColl_hn_nearby.at(j).DeltaR(this_muon), weight, 0., 6., 60);
          }
        }

        FillHist("TEST_Muon_FR_Fakable", 0., 1., 0., 2., 2);
        if(this_muon.RelIso04()<0.1){
          FillHist("TEST_Muon_FR_Fakable", 1., 1., 0., 2., 2);
        }
      }

    }

    std::vector<snu::KElectron> electrontriLooseColl = GetElectrons(true, true, "ELECTRON_HN_LOWDXY_FAKELOOSE");
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      snu::KElectron this_electron = electrontriLooseColl.at(i);

      if(this_electron.MCIsFromConversion()){
        FillHist("TEST_Electron_FR_Conversion", 0., 1., 0., 2., 2);
        if(eventbase->GetElectronSel()->ElectronPass(this_electron,"ELECTRON_HN_LOWDXY_TIGHT")){
          FillHist("TEST_Electron_FR_Conversion", 1., 1., 0., 2., 2);
        }
      }

      if(this_electron.MCFromTau()){
        FillHist("TEST_Electron_FR_Tau", 0., 1., 0., 2., 2);
        if(eventbase->GetElectronSel()->ElectronPass(this_electron,"ELECTRON_HN_LOWDXY_TIGHT")){
          FillHist("TEST_Electron_FR_Tau", 1., 1., 0., 2., 2);
        }
      }

      if(!this_electron.MCMatched() && !this_electron.MCIsFromConversion() && !this_electron.MCFromTau()){
        FillHist("TEST_Electron_FR_Fakable", 0., 1., 0., 2., 2);
        if(eventbase->GetElectronSel()->ElectronPass(this_electron,"ELECTRON_HN_LOWDXY_TIGHT")){
          FillHist("TEST_Electron_FR_Fakable", 1., 1., 0., 2., 2);
        }
      }
    }

    return;
*/


  }
  //==== diboson, but hadronic decays
  else if( (k_sample_name.Contains("WZ") || k_sample_name.Contains("ZZ") || k_sample_name.Contains("WW") ) && diboson_had ){
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
    electrontriLooseColl = GetElectrons(false, false, "ELECTRON_HN_LOWDXY_FAKELOOSE");
    if(muontriLooseColl.size()+electrontriLooseColl.size()==3) return;
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
    electrontriLooseColl = GetElectrons(false, true, "ELECTRON_HN_LOWDXY_FAKELOOSE");
  }
  //==== otherwise
  else{
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
    electrontriLooseColl = GetElectrons(false, false, "ELECTRON_HN_LOWDXY_FAKELOOSE");
  }

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);

  //=========================== 
  //==== Get Muon Corrections
  //===========================

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muontriLooseColl);

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

  //===========================
  //==== Trigger Scale Factor
  //===========================

  double trigger_sf = 1.;
  if(!k_isdata){
    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, muontriLooseColl, 0, 0, 0);
    double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrontriLooseColl, muontriLooseColl, 0, 1, 0);
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
    weight*=muon_id_iso_sf;
    weight*=trigger_sf;
    weight*=MuTrkEffSF;
    weight*=trigger_ps_weight;
    weight*=pileup_reweight;
    weight*=GetKFactor();
  }

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
    if(PassID(electrontriLooseColl.at(i), "ELECTRON_HN_LOWDXY_TIGHT")) n_triTight_electrons++;
  }

  int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
  int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

/*

  //==== Tight RelIso Study

  for(int i=0; i<100; i++){
    double test_reliso = 0.01*(i+1);
    std::vector<snu::KMuon> TestMuon = GetHNTriMuonsByLooseRelIso(test_reliso, true);
    if(TestMuon.size()==3){
      FillHist("TEST_TightRelIso", i, 1., 0., 100., 100);
    }
  }
  return;
*/

  //==== ppp to TTL/TLL/LLL ?
  if(!k_isdata){
    if(n_triLoose_muons==3){
      FillHist("PPP_nTight", n_triTight_muons, weight, 0., 4., 4);
    }
  }

  FillHist("n_loose_muon", n_triLoose_muons, 1., 0., 10., 10);
  FillHist("n_tight_muon", n_triTight_muons, 1., 0., 10., 10);

/*

  //==== ZG study

  if(k_sample_name.Contains("ZG")){

    std::vector<snu::KElectron> allel = GetElectrons(false, true, "ELECTRON_HN_LOWDXY_TIGHT");
    bool isthisit=false;
    for(unsigned int i=0; i<allel.size(); i++){
      if(!allel.at(i).MCMatched() && allel.at(i).MCIsFromConversion()){
        cout << "pt = " << allel.at(i).Pt() << "\t" << allel.at(i).Eta() << endl;
        isthisit = true;
      }
    }

    if(!isthisit) return;

    std::vector<snu::KTruth> truthColl;
    eventbase->GetTruthSel()->Selection(truthColl);
    cout << "=========================================================" << endl;
    cout << "RunNumber = " << eventbase->GetEvent().RunNumber() << endl;
    cout << "EventNumber = " << eventbase->GetEvent().EventNumber() << endl;
    cout << "truth size = " << truthColl.size() << endl;
    cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
    for(int i=2; i<truthColl.size(); i++){
      cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << "\t" << truthColl.at(i).Pt() << "\t" << truthColl.at(i).Eta() << endl;
    }

    return;
  }
*/

  //=============================
  //==== [CUT] Three Tight mons
  //=============================

  //==== Three Tight Muons, and no fourth loose lepton
  bool isThreeMuon     = (n_triLoose_leptons == 3)
                         && (n_triLoose_muons == 3 && n_triTight_muons == 3);
  bool isTwoMuonOneElectron = (n_triLoose_leptons == 3)
                              && (n_triLoose_muons == 2 && n_triTight_muons == 2)
                              && (n_triLoose_electrons == 1 && n_triTight_electrons == 1);

  if(!isThreeMuon && !isTwoMuonOneElectron) return;

  double MinLeadingMuonPt = 20;
  if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) return;
  FillCutFlow("3muon", 1.);
  FillHist("cutflow_MuMuE", 4., 1., 0., 10., 10);

  KLepton lep[3];
  snu::KParticle HN[4];
  for(unsigned int i=0;i<n_triLoose_muons;i++){
    lep[i] = muontriLooseColl.at(i);
  }
  for(unsigned int i=0;i<n_triLoose_electrons;i++){
    lep[i+n_triLoose_muons] = electrontriLooseColl.at(i);
  }

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

    FillHist("cutflow_MuMuE", 5., 1., 0., 10., 10);
    
  }
  else return;

  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowssllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);

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
  //MET = m_HNgenmatch->gen_nu.Pt();
  //METphi = m_HNgenmatch->gen_nu.Phi();
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

  if(isThreeMuon){
    FillHist("CutStudy_m_Z_candidate", z_candidate.M(), 1., 0., 1000., 1000);
    bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
    if(!VetoZResonance) return;
    FillCutFlow("ZVeto", 1.);

    FillHist("CutStudy_mlll", (lep[0] + lep[1] + lep[2]).M(), 1., 0., 1000., 1000);
    bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;
    if(!mllloffZ) return;
    FillCutFlow("mllloffZ", 1.);
  }

  FillHist("CutStudy_nbjet", n_bjets, 1., 0., 10., 10);

  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);
  FillHist("cutflow_MuMuE", 6., 1., 0., 10., 10);

  //==== preselection is done

/*

  //==== Trigger Check...

  FillHist("TEST_Triggers", 0., 1., 0., 10., 10);
  std::vector<TString> TESTtriggerlist1;
  TESTtriggerlist1.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  TESTtriggerlist1.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  if(PassTriggerOR(TESTtriggerlist1)){
    FillHist("TEST_Triggers", 1., 1., 0., 10., 10);
  }
  std::vector<TString> TESTtriggerlist2;
  TESTtriggerlist2.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  TESTtriggerlist2.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(PassTriggerOR(TESTtriggerlist2)){
    FillHist("TEST_Triggers", 2., 1., 0., 10., 10);
  }
  std::vector<TString> TESTtriggerlist3;
  TESTtriggerlist3.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  TESTtriggerlist3.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  TESTtriggerlist3.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  TESTtriggerlist3.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  if(PassTriggerOR(TESTtriggerlist3)){
    FillHist("TEST_Triggers", 3., 1., 0., 10., 10);
  }
  return;
*/

  if(DoCutOp){

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
      FillHist("TriggerStudy_weighted", 0., weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_1)){
      FillHist("TriggerStudy_unweighted", 1., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 1., weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_2)){
      FillHist("TriggerStudy_unweighted", 2., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 2., weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_3)){
      FillHist("TriggerStudy_unweighted", 3., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 3., weight, 0., 5., 5);
    }
    if(PassTriggerOR(triggerlist_123)){
      FillHist("TriggerStudy_unweighted", 4., 1., 0., 5., 5);
      FillHist("TriggerStudy_weighted", 4., weight, 0., 5., 5);
    }

  }

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);
  SetPlotHNTriLepMetInfo(MET, METphi);
  SetPlotHNTriLepParticleInfo(HN, W_pri_lowmass, W_pri_highmass, W_sec);
  SetPlotHNTriLepChargeSign(OppSign, SameSign[0], SameSign[1]);
  SetPlotHNTriBJet(n_bjets);

  bool isLowMass = (W_pri_lowmass.M() < 150.);
  bool isHighMass = (MET > 20.);

  if(isThreeMuon){

    FillCLHist(hntrilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);

    if( isLowMass ){
      FillCLHist(hntrilephist, "cutWlow", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
      FillCutFlow("LowMass", 1.);
    }

    if( isHighMass ){
      FillCLHist(hntrilephist, "cutWhigh", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
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

      //cout << "Tring PassOptimizedCut" << endl;
      bool pass_op = PassOptimizedCut(signal_masses[i],
        muontriLooseColl.at(0).Pt(), muontriLooseColl.at(1).Pt(), muontriLooseColl.at(2).Pt(),
        this_W_pri_mass, hnmass,
        deltaR_OS_min, gamma_star.M(),
        MET
      );
      //cout << "==>Done" << endl;

      if(pass_op){
        FillCLHist(hntrilephist, thiscut, eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
      }

    }

  }
  else if(isTwoMuonOneElectron){
    FillHist("TEST_MuMuE_nevent", 0., weight, 0., 1., 1);
    FillCLHist(hntrilephist, "MuMuE", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);

    if(DoCutOp){
      bool PassIsoMu24 = PassTrigger("HLT_IsoMu24_v") || PassTrigger("HLT_IsoTkMu24_v");
      if(PassIsoMu24){
        FillHist("TEST_MuMuE_IsoMu24", 0., weight, 0., 1., 1);
      }
      if(PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_DiMu", 0., weight, 0., 1., 1);
      }
      if(PassIsoMu24 || PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_IsoMu24_OR_DiMu", 0., weight, 0., 1., 1);
      }
      if(PassIsoMu24 && !PassTriggerOR(triggerlist)){
        FillHist("TEST_MuMuE_IsoMu24_AND_NOT_DiMu", 0., weight, 0., 1., 1);

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
  


void trilepton_mumumu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  cout
  << "n_gen_pass = " << m_HNgenmatch->n_gen_pass << endl
  << "best = " << m_HNgenmatch->sol_sel_chi2_best/m_HNgenmatch->n_gen_pass << endl
  << "plus = " << m_HNgenmatch->sol_sel_chi2_plus/m_HNgenmatch->n_gen_pass << endl
  << "minus = " << m_HNgenmatch->sol_sel_chi2_minus/m_HNgenmatch->n_gen_pass << endl
  << "smaller = " << m_HNgenmatch->sol_sel_chi2_smaller/m_HNgenmatch->n_gen_pass << endl
  << "larger = " << m_HNgenmatch->sol_sel_chi2_larger/m_HNgenmatch->n_gen_pass << endl;

  TH1D* GEN_solution_selection_chi2 = new TH1D("GEN_solution_selection_chi2", "", 6, 0, 6);
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(1, "n_gen_pass");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(2, "best");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(3, "plus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(4, "minus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(5, "smaller");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(6, "larger");
  GEN_solution_selection_chi2->SetBinContent(1, m_HNgenmatch->n_gen_pass);
  GEN_solution_selection_chi2->SetBinContent(2, m_HNgenmatch->sol_sel_chi2_best);
  GEN_solution_selection_chi2->SetBinContent(3, m_HNgenmatch->sol_sel_chi2_plus);
  GEN_solution_selection_chi2->SetBinContent(4, m_HNgenmatch->sol_sel_chi2_minus);
  GEN_solution_selection_chi2->SetBinContent(5, m_HNgenmatch->sol_sel_chi2_smaller);
  GEN_solution_selection_chi2->SetBinContent(6, m_HNgenmatch->sol_sel_chi2_larger);

  m_outputFile->cd();
  GEN_solution_selection_chi2->Write();


}


void trilepton_mumumu::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu::~trilepton_mumumu() {
  
  Message("In trilepton_mumumu Destructor" , INFO);
  
}


void trilepton_mumumu::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumuCore::MakeHistograms() to make new hists for your analysis
   **/

}


void trilepton_mumumu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

int trilepton_mumumu::GetSignalMass(){

  //==== NEW : HN_MuMuMu_5
  //==== OLD : HN40_mumumu_VmuN_0p1

  bool OldSig = k_sample_name.Contains("VmuN_0p1");  
  if(OldSig){
    if(k_sample_name.Contains("HN40_mumumu_VmuN_0p1_")) return 40;
    if(k_sample_name.Contains("HN60_mumumu_VmuN_0p1_")) return 60;
    if(k_sample_name.Contains("HN150_mumumu_VmuN_0p1_")) return 150;
    if(k_sample_name.Contains("HN700_mumumu_VmuN_0p1_")) return 700;
  }
  else{
    if(k_sample_name.Contains("HN_MuMuMu_5_")) return 5;
    if(k_sample_name.Contains("HN_MuMuMu_10_")) return 10;
    if(k_sample_name.Contains("HN_MuMuMu_20_")) return 20;
    if(k_sample_name.Contains("HN_MuMuMu_30_")) return 30;
    if(k_sample_name.Contains("HN_MuMuMu_40_")) return 40;
    if(k_sample_name.Contains("HN_MuMuMu_50_")) return 50;
    if(k_sample_name.Contains("HN_MuMuMu_60_")) return 60;
    if(k_sample_name.Contains("HN_MuMuMu_70_")) return 70;
    if(k_sample_name.Contains("HN_MuMuMu_90_")) return 90;
    if(k_sample_name.Contains("HN_MuMuMu_100_")) return 100;
    if(k_sample_name.Contains("HN_MuMuMu_150_")) return 150;
    if(k_sample_name.Contains("HN_MuMuMu_200_")) return 200;
    if(k_sample_name.Contains("HN_MuMuMu_300_")) return 300;
    if(k_sample_name.Contains("HN_MuMuMu_400_")) return 400;
    if(k_sample_name.Contains("HN_MuMuMu_500_")) return 500;
    if(k_sample_name.Contains("HN_MuMuMu_700_")) return 700;
    if(k_sample_name.Contains("HN_MuMuMu_1000_")) return 1000;
  }
  return 0;
  

}


