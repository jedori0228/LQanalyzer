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

  //====================
  //==== Prepare Muons
  //====================

  double this_RelIso = 0.4;
  std::vector<snu::KMuon> muontriLooseColl;
  bool diboson_had = std::find(k_flags.begin(), k_flags.end(), "diboson_had") != k_flags.end();
  //==== signal
  if( k_sample_name.Contains("HN") ){

    snu::KTruth gen_l_1 = m_HNgenmatch->gen_l_1;
    snu::KTruth gen_l_2 = m_HNgenmatch->gen_l_2;
    snu::KTruth gen_l_3 = m_HNgenmatch->gen_l_3;

    std::vector<snu::KMuon> muontriLooseColl_raw = GetHNTriMuonsByLooseRelIso(this_RelIso, true);

    //==== find gen_l_1
    //cout << "[gen_l_1] : pt = " << gen_l_1.Pt() << ", eta = " << gen_l_1.Eta() << endl;
    //cout << "[gen_l_2] : pt = " << gen_l_2.Pt() << ", eta = " << gen_l_2.Eta() << endl;
    //cout << "[gen_l_3] : pt = " << gen_l_3.Pt() << ", eta = " << gen_l_3.Eta() << endl;
    std::vector<int> loose_used;
    loose_used.clear();
    int loose_l_1_index = find_genmatching(gen_l_1, muontriLooseColl_raw, loose_used);
    int loose_l_2_index = find_genmatching(gen_l_2, muontriLooseColl_raw, loose_used);
    int loose_l_3_index = find_genmatching(gen_l_3, muontriLooseColl_raw, loose_used);

    std::vector<snu::KMuon> muontriLooseColl_genorder;
    if(loose_l_1_index!=-1){
      muontriLooseColl_genorder.push_back( muontriLooseColl_raw.at(loose_l_1_index) );
      FillHist("CH_F0", 1./muontriLooseColl_raw.at(loose_l_1_index).Pt(), 1., 0., 0.1, 100);
      if( muontriLooseColl_raw.at(loose_l_1_index).Charge() * gen_l_1.PdgId() > 0 ){
        FillHist("CH_F", 1./muontriLooseColl_raw.at(loose_l_1_index).Pt(), 1., 0., 0.1, 100);
      }
      FillHist("resPt", (gen_l_1.Pt()-muontriLooseColl_raw.at(loose_l_1_index).Pt())/muontriLooseColl_raw.at(loose_l_1_index).Pt(), 1., -1.0, 1.0, 200);
    }
    if(loose_l_2_index!=-1){
      muontriLooseColl_genorder.push_back( muontriLooseColl_raw.at(loose_l_2_index) );
      FillHist("CH_F0", 1./muontriLooseColl_raw.at(loose_l_2_index).Pt(), 1., 0., 0.1, 100);
      if( muontriLooseColl_raw.at(loose_l_2_index).Charge() * gen_l_2.PdgId() > 0 ){
        FillHist("CH_F", 1./muontriLooseColl_raw.at(loose_l_2_index).Pt(), 1., 0., 0.1, 100);
      }
      FillHist("resPt", (gen_l_2.Pt()-muontriLooseColl_raw.at(loose_l_2_index).Pt())/muontriLooseColl_raw.at(loose_l_2_index).Pt(), 1., -1.0, 1.0, 200);
    }
    if(loose_l_3_index!=-1){
      muontriLooseColl_genorder.push_back( muontriLooseColl_raw.at(loose_l_3_index) );
      FillHist("CH_F0", 1./muontriLooseColl_raw.at(loose_l_3_index).Pt(), 1., 0., 0.1, 100);
      if( muontriLooseColl_raw.at(loose_l_3_index).Charge() * gen_l_3.PdgId() > 0 ){
        FillHist("CH_F", 1./muontriLooseColl_raw.at(loose_l_3_index).Pt(), 1., 0., 0.1, 100);
      }
      FillHist("resPt", (gen_l_3.Pt()-muontriLooseColl_raw.at(loose_l_3_index).Pt())/muontriLooseColl_raw.at(loose_l_3_index).Pt(), 1., -1.0, 1.0, 200);
    }

    muontriLooseColl = sort_muons_ptorder( muontriLooseColl_genorder );

  }
  //==== non-prompt : keep fake
  else if( k_sample_name.Contains("DY") || k_sample_name.Contains("WJets") || k_sample_name.Contains("TTJets") || k_sample_name.Contains("QCD") ){
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);

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
    if(muontriLooseColl.size()==3) return;
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  }
  //==== otherwise
  else{
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
  }

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
      FillHist("bjet_pt", jetColl_hn.at(j).Pt(), 1., 0., 200., 200);
    }
  }

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrontriLooseColl = GetElectrons(false, false, "ELECTRON_HN_LOWDXY_FAKELOOSE");

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

  if(!isThreeMuon) return;

  double MinLeadingMuonPt = 20;
  if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) return;
  FillCutFlow("3muon", 1.);

  snu::KParticle lep[3], HN[4];
  for(int i=0;i<3;i++){
    lep[i] = muontriLooseColl.at(i);
  }

  int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
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

  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowosllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowssllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("CutStudy_lowllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);

  //==== MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
  //==== MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
  //==== POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
  if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
      (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) return;
  FillCutFlow("mllsf4", 1.);

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET();
  double METphi = Evt.METPhi();
  CorrectedMETRochester(muontriLooseColl, MET, METphi);
  m_HNgenmatch->SetMETInfo(MET, METphi);
  if(k_sample_name.Contains("HN") && m_HNgenmatch->allgenfound) m_HNgenmatch->solution_selection_study(muontriLooseColl);

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
  int l_3_index = find_mlmet_closest_to_W(lep, nu_highmass);
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

  if(l_3_index == OppSign){
     
    HN[2] = W_sec + lep[SameSign[1]]; // [class3]
    HN[3] = W_sec + lep[SameSign[0]]; // [class4]
 
  }
  else{
      HN[2] = W_sec + lep[OppSign]; // [class3]
      HN[3] = W_sec + lep[OppSign]; // [class4]
  }

  FillHist("CutStudy_m_Z_candidate", z_candidate.M(), 1., 0., 1000., 1000);

  bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
  if(!VetoZResonance) return;
  FillCutFlow("ZVeto", 1.);

  FillHist("CutStudy_mlll", (lep[0] + lep[1] + lep[2]).M(), 1., 0., 1000., 1000);

  bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;
  if(!mllloffZ) return;
  FillCutFlow("mllloffZ", 1.);

  FillHist("CutStudy_nbjet", n_bjets, 1., 0., 10., 10);

  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);

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

    return;

  }

  muontriLooseColl = sort_muons_ptorder(muontriLooseColl);
  SetPlotHNTriLepMetInfo(MET, METphi);

  bool isLowMass = (W_pri_lowmass.M() < 150.);
  bool isHighMass = (MET > 20.);

  FillHist("HN_mass_class1_cut0", HN[0].M(), weight, 0., 2000., 2000);
  FillHist("HN_mass_class2_cut0", HN[1].M(), weight, 0., 2000., 2000);
  FillHist("HN_mass_class3_cut0", HN[2].M(), weight, 0., 2000., 2000);
  FillHist("HN_mass_class4_cut0", HN[3].M(), weight, 0., 2000., 2000);
  FillHist("W_pri_lowmass_mass_cut0", W_pri_lowmass.M(), weight, 0., 2000., 2000);
  FillHist("W_pri_highmass_mass_cut0", W_pri_highmass.M(), weight, 0., 2000., 2000);
  FillHist("W_sec_highmass_mass_cut0", W_sec.M(), weight, 0., 2000., 2000);
  FillHist("deltaR_OS_min_cut0", deltaR_OS_min, weight, 0., 5., 50);
  FillHist("gamma_star_mass_cut0", gamma_star.M(), weight, 0., 200., 200);
  FillHist("z_candidate_mass_cut0", z_candidate.M(), weight, 0., 200., 200);
  FillCLHist(hntrilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
  FillHist("n_events_cut0", 0, weight, 0., 1., 1);

  if( isLowMass ){
    FillHist("HN_mass_class1_cutWlow", HN[0].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class2_cutWlow", HN[1].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class3_cutWlow", HN[2].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class4_cutWlow", HN[3].M(), weight, 0., 2000., 2000);
    FillHist("W_pri_lowmass_mass_cutWlow", W_pri_lowmass.M(), weight, 0., 2000., 2000);
    FillHist("W_pri_highmass_mass_cutWlow", W_pri_highmass.M(), weight, 0., 2000., 2000);
    FillHist("W_sec_highmass_mass_cutWlow", W_sec.M(), weight, 0., 2000., 2000);
    FillHist("deltaR_OS_min_cutWlow", deltaR_OS_min, weight, 0., 5., 50);
    FillHist("gamma_star_mass_cutWlow", gamma_star.M(), weight, 0., 200., 200);
    FillHist("z_candidate_mass_cutWlow", z_candidate.M(), weight, 0., 200., 200);
    FillCLHist(hntrilephist, "cutWlow", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
    FillHist("n_events_cutWlow", 0, weight, 0., 1., 1);

    FillCutFlow("LowMass", 1.);

  }

  if( isHighMass ){
    FillHist("HN_mass_class1_cutWhigh", HN[0].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class2_cutWhigh", HN[1].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class3_cutWhigh", HN[2].M(), weight, 0., 2000., 2000);
    FillHist("HN_mass_class4_cutWhigh", HN[3].M(), weight, 0., 2000., 2000);
    FillHist("W_pri_lowmass_mass_cutWhigh", W_pri_lowmass.M(), weight, 0., 2000., 2000);
    FillHist("W_pri_highmass_mass_cutWhigh", W_pri_highmass.M(), weight, 0., 2000., 2000);
    FillHist("W_sec_highmass_mass_cutWhigh", W_sec.M(), weight, 0., 2000., 2000);
    FillHist("deltaR_OS_min_cutWhigh", deltaR_OS_min, weight, 0., 5., 50);
    FillHist("gamma_star_mass_cutWhigh", gamma_star.M(), weight, 0., 200., 200);
    FillHist("z_candidate_mass_cutWhigh", z_candidate.M(), weight, 0., 200., 200);
    FillCLHist(hntrilephist, "cutWhigh", eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
    FillHist("n_events_cutWhigh", 0, weight, 0., 1., 1);

    FillCutFlow("HighMass", 1.);

  }

  int signal_masses[] = {5, 10, 20, 30, 40, 50, 60, 70, 90, 100, 150, 200, 300, 400, 500, 700, 1000};
  for(int i=0; i<17; i++){

    TString thiscut = "cutHN"+TString::Itoa(signal_masses[i],10);
    double this_W_pri_mass = W_pri_lowmass.M();
    if( signal_masses[i] > 80 ) this_W_pri_mass = W_pri_highmass.M();

    if( PassOptimizedCut(signal_masses[i], muontriLooseColl.at(0).Pt(), muontriLooseColl.at(1).Pt(), muontriLooseColl.at(2).Pt(), this_W_pri_mass, MET) ){
      FillHist("HN_mass_class1_"+thiscut, HN[0].M(), weight, 0., 2000., 2000);
      FillHist("HN_mass_class2_"+thiscut, HN[1].M(), weight, 0., 2000., 2000);
      FillHist("HN_mass_class3_"+thiscut, HN[2].M(), weight, 0., 2000., 2000);
      FillHist("HN_mass_class4_"+thiscut, HN[3].M(), weight, 0., 2000., 2000);
      FillHist("W_pri_lowmass_mass_"+thiscut, W_pri_lowmass.M(), weight, 0., 2000., 2000);
      FillHist("W_pri_highmass_mass_"+thiscut, W_pri_highmass.M(), weight, 0., 2000., 2000);
      FillHist("W_sec_highmass_mass_"+thiscut, W_sec.M(), weight, 0., 2000., 2000);
      FillHist("deltaR_OS_min_"+thiscut, deltaR_OS_min, weight, 0., 5., 50);
      FillHist("gamma_star_mass_"+thiscut, gamma_star.M(), weight, 0., 200., 200);
      FillHist("z_candidate_mass_"+thiscut, z_candidate.M(), weight, 0., 200., 200);
      FillCLHist(hntrilephist, thiscut, eventbase->GetEvent(), muontriLooseColl, electrontriLooseColl, jetColl_hn, weight);
      FillHist("n_events_"+thiscut, 0, weight, 0., 1., 1);
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

int trilepton_mumumu::find_genmatching(snu::KParticle gen, std::vector<snu::KMuon> recos, std::vector<int>& used_index){

  double mindr = 0.1;
  double maxPtDiff = 0.05;
  int found=-1;
  for(unsigned int i=0; i<recos.size(); i++){
    //cout << "["<<i<<"th reco] : pt = " << recos.at(i).Pt() << ", eta = " << recos.at(i).Eta() << endl;
    double dr = gen.DeltaR(recos.at(i));
    bool ptmatched(true);
    //if( fabs(gen.Pt() - recos.at(i).Pt()) / TMath::Min( gen.Pt(), recos.at(i).Pt() ) >= maxPtDiff ) ptmatched = false;
    bool isthisused = std::find(used_index.begin(), used_index.end(), i) != used_index.end();
    if(dr < mindr && ptmatched && !isthisused){
      mindr = dr;
      found = i;
    }
  }
  used_index.push_back(found);
  return found;
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

bool trilepton_mumumu::PassOptimizedCut(int sig_mass, double first_pt, double second_pt, double third_pt, double W_pri_mass, double PFMET){

  double cut_first_pt(0.), cut_second_pt(0.), cut_third_pt(0.), cut_W_pri_mass(0.), cut_PFMET(0.);

  if(sig_mass == 5){
    cut_first_pt = 60.;
    cut_second_pt = 45.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 125.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 10){
    cut_first_pt = 55.;
    cut_second_pt = 40.;
    cut_third_pt = 35.;
    cut_W_pri_mass = 130.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 20){
    cut_first_pt = 50.;
    cut_second_pt = 40.;
    cut_third_pt = 40.;
    cut_W_pri_mass = 130.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 30){
    cut_first_pt = 45.;
    cut_second_pt = 40.;
    cut_third_pt = 35.;
    cut_W_pri_mass = 130.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 40){
    cut_first_pt = 35.;
    cut_second_pt = 30.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 130.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 50){
    cut_first_pt = 30.;
    cut_second_pt = 30.;
    cut_third_pt = 30.;
    cut_W_pri_mass = 125.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 60){
    cut_first_pt = 30.;
    cut_second_pt = 25.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 130.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 70){
    cut_first_pt = 35.;
    cut_second_pt = 30.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 125.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 90){
    cut_first_pt = 45.;
    cut_second_pt = 40.;
    cut_third_pt = 15.;
    cut_W_pri_mass = 80.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 100){
    cut_first_pt = 30.;
    cut_second_pt = 15.;
    cut_third_pt = 15.;
    cut_W_pri_mass = 110.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 150){
    cut_first_pt = 45.;
    cut_second_pt = 40.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 160.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 200){
    cut_first_pt = 65.;
    cut_second_pt = 55.;
    cut_third_pt = 30.;
    cut_W_pri_mass = 250.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 300){
    cut_first_pt = 120.;
    cut_second_pt = 75.;
    cut_third_pt = 45.;
    cut_W_pri_mass = 350.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 400){
    cut_first_pt = 120.;
    cut_second_pt = 65.;
    cut_third_pt = 50.;
    cut_W_pri_mass = 480.;
    cut_PFMET = 40.;
  }
  else if(sig_mass == 500){
    cut_first_pt = 150.;
    cut_second_pt = 100.;
    cut_third_pt = 50.;
    cut_W_pri_mass = 530.;
    cut_PFMET = 50.;
  }
  else if(sig_mass == 700){
    cut_first_pt = 200.;
    cut_second_pt = 100.;
    cut_third_pt = 45.;
    cut_W_pri_mass = 760.;
    cut_PFMET = 50.;
  }
  else if(sig_mass == 1000){
    cut_first_pt = 290.;
    cut_second_pt = 180.;
    cut_third_pt = 50.;
    cut_W_pri_mass = 920.;
    cut_PFMET = 50.;
  }
  else{
    cout << "[trilepton_mumumu::PassOptimizedCut] Signal mass wrong" << endl;
  }

  bool pass = true;
  if(sig_mass < 80){
    if( !(first_pt < cut_first_pt) ) pass = false;
    if( !(second_pt < cut_second_pt) ) pass = false;
    if( !(third_pt < cut_third_pt) ) pass = false;
    if( !(W_pri_mass < cut_W_pri_mass) ) pass = false;
  }
  else{
    if( !(first_pt > cut_first_pt) ) pass = false;
    if( !(second_pt > cut_second_pt) ) pass = false;
    if( !(third_pt > cut_third_pt) ) pass = false;
    if( !(PFMET > cut_PFMET) ) pass = false;
    if( !(W_pri_mass > cut_W_pri_mass) ) pass = false;
  }

  return pass;

}


