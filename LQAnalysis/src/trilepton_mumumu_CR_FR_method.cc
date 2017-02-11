// $Id: trilepton_mumumu_CR_FR_method.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_CR_FR_method Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_CR_FR_method.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_CR_FR_method);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu_CR_FR_method::trilepton_mumumu_CR_FR_method() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_CR_FR_method");
  
  Message("In trilepton_mumumu_CR_FR_method constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_CR_FR_method::InitialiseAnalysis() throw( LQError ) {
  
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

  this_dXYSig = 4.0;
  this_RelIso = 0.4;

  return;

}


void trilepton_mumumu_CR_FR_method::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerlist;
  //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  //triggerlist.push_back("HLT_TripleMu_12_10_5_v");
  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  bool trigger_pass = false;
  for(unsigned int i=0; i<triggerlist.size(); i++){
    if( PassTrigger(triggerlist.at(i)) ){
      trigger_pass = true;
      break;
    }
  }
  if(!trigger_pass) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", weight);


  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  //std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  //std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  //std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  //std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed
  std::vector<snu::KMuon> muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);

  //CorrectMuonMomentum(muontriLooseColl); //FIXME do this for v8-0-4

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", true, 30., 2.4);
   
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  numberVertices = eventbase->GetEvent().nVertices();   

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    pileup_reweight = eventbase->GetEvent().PileUpWeight();
    //pileup_reweight = TempPileupWeight();
  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();

  //==== No normalization for MC Closure
  if(DoMCClosure){
    weight = 1.*MCweight;
    m_fakeobj->SetUseQCDFake(true); 
  }

  int n_triTight_muons(0);
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    if(muontriLooseColl.at(i).RelIso04() < 0.1) n_triTight_muons++;
  }
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_hn.size();

  FillHist("GenWeight_NJet" , n_jets*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  // CR related variables //
  FillHist("n_tight_muons_control", n_triTight_muons, weight*pileup_reweight, 0., 10., 10);
  FillHist("n_loose_muons_control", n_triLoose_muons, weight*pileup_reweight, 0., 10., 10);
  if( n_triTight_muons == 2 && n_triLoose_muons == 2){
    int isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ? 1 : 0;
    FillHist("2Muons_OS0_SS1_control", isSS, weight*pileup_reweight, 0., 2., 2);
  }
  FillHist("n_jets_control", n_jets, weight*pileup_reweight, 0., 10., 10);
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)) n_bjets++;
  }
  FillHist("n_bjets_control", n_bjets, weight*pileup_reweight, 0., 10., 10);

  //================
  //==== define CR
  //================

  int n_muons(0);

  bool isTwoMuon   = n_triLoose_muons == 2 && n_triTight_muons != 2; // no TT case
  bool isThreeMuon = n_triLoose_muons == 3 && n_triTight_muons != 3; // no TTT case
  if(n_triLoose_muons == 2 && n_triTight_muons == 0) FillHist("LL_TL_TT", 0., 1., 0., 3., 3);
  if(n_triLoose_muons == 2 && n_triTight_muons == 1) FillHist("LL_TL_TT", 1., 1., 0., 3., 3);
  if(n_triLoose_muons == 2 && n_triTight_muons == 2) FillHist("LL_TL_TT", 2., 1., 0., 3., 3);

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();

  //==== CR with Two Muons
  if(DoMCClosure && isTwoMuon){
    snu::KMuon lep[2];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);

    bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();

    double m_Z = 91.1876;
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();
    bool ZResonance = fabs(m_dimuon-m_Z) < 10.;

    std::map< TString, bool > map_whichCR_to_isCR;
    map_whichCR_to_isCR.clear();
    map_whichCR_to_isCR["DiMuon"] = isTwoMuon && leadPt20;
    map_whichCR_to_isCR["SSDiMuon"] = isTwoMuon && leadPt20 && isSS;
    map_whichCR_to_isCR["OSDiMuon"] = isTwoMuon && leadPt20 && !isSS;
    map_whichCR_to_isCR["OSDiMuon_Z_10GeV"] = isTwoMuon && leadPt20 && !isSS && ZResonance;

    //==== fake method weighting
    std::vector<snu::KElectron> empty_electron;
    empty_electron.clear();
    double FR_reweight = Get_DataDrivenWeight(false, muontriLooseColl, empty_electron, 2, 0);
    double weight_err = Get_DataDrivenWeight(true, muontriLooseColl, empty_electron, 2, 0);

    for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
      TString this_suffix = it->first;
      if(it->second){
        FillHist("weight_"+this_suffix, weight*FR_reweight, 1., -1., 1., 1000);
        FillUpDownHist("n_events_"+this_suffix+"", 0, weight*FR_reweight, weight_err, 0., 1., 1);
        FillUpDownHist("n_jets_"+this_suffix+"", n_jets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("n_bjets_"+this_suffix+"", n_bjets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("PFMET_"+this_suffix+"", MET, weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("PFMET_phi_"+this_suffix+"", METphi, weight*FR_reweight, weight_err, -3.2, 3.2, 64);
        FillUpDownHist("mll_"+this_suffix+"", m_dimuon , weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("n_vertices_"+this_suffix, numberVertices, weight, weight_err, 0., 50., 50);
        FillUpDownHist("leadingLepton_Pt_"+this_suffix+"", lep[0].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("leadingLepton_Eta_"+this_suffix+"", lep[0].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("leadingLepton_RelIso_"+this_suffix+"", lep[0].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("leadingLepton_Chi2_"+this_suffix+"", lep[0].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);
        FillUpDownHist("leadingLepton_dXY_"+this_suffix+"", fabs(lep[0].dXY()) , weight*FR_reweight, weight_err, 0., 0.1, 100);
        FillUpDownHist("leadingLepton_dXYSig_"+this_suffix+"", fabs(lep[0].dXYSig()) , weight*FR_reweight, weight_err, 0., 4., 40);
        FillUpDownHist("secondLepton_Pt_"+this_suffix+"", lep[1].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("secondLepton_Eta_"+this_suffix+"", lep[1].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("secondLepton_RelIso_"+this_suffix+"", lep[1].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("secondLepton_Chi2_"+this_suffix+"", lep[1].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);
        FillUpDownHist("secondLepton_dXY_"+this_suffix+"", fabs(lep[1].dXY()) , weight*FR_reweight, weight_err, 0., 0.1, 100);
        FillUpDownHist("secondLepton_dXYSig_"+this_suffix+"", fabs(lep[1].dXYSig()) , weight*FR_reweight, weight_err, 0., 4., 40);

      }
    } 

  } // is TwoMuon

  if(isThreeMuon){
    snu::KMuon lep[3];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);
    lep[2] = muontriLooseColl.at(2);

    //==== fake method weighting
    std::vector<snu::KElectron> empty_electron;
    empty_electron.clear();
    double FR_reweight =  Get_DataDrivenWeight(false, muontriLooseColl, empty_electron, 3, 0);
    double weight_err = Get_DataDrivenWeight(true, muontriLooseColl, empty_electron, 3, 0);

    //bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.; // This will be done for the Z-candidate muons
    bool AllSameCharge = ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) &&
                         ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(2).Charge() );

    if( !AllSameCharge ){
      snu::KMuon OS, SS[2];
      if     ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ){
        SS[0] = muontriLooseColl.at(0);
        SS[1] = muontriLooseColl.at(1);
        OS    = muontriLooseColl.at(2);
      }
      else if( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(2).Charge() ){
        SS[0] = muontriLooseColl.at(0);
        SS[1] = muontriLooseColl.at(2);
        OS    = muontriLooseColl.at(1);
      }
      else if( muontriLooseColl.at(1).Charge() == muontriLooseColl.at(2).Charge() ){
        SS[0] = muontriLooseColl.at(1);
        SS[1] = muontriLooseColl.at(2);
        OS    = muontriLooseColl.at(0);
      }
      else Message("?", INFO);

      double m_dimuon[2], m_Z = 91.1876;
      m_dimuon[0] = ( SS[0] + OS ).M();
      m_dimuon[1] = ( SS[1] + OS ).M();

      snu::KParticle Z_candidate;
      snu::KMuon ZMuon, WMuon;
      if( fabs(m_dimuon[0]-m_Z) < fabs(m_dimuon[1]-m_Z) ){
        Z_candidate = SS[0] + OS;
        ZMuon = SS[0];
        WMuon = SS[1];
      }
      else{
        Z_candidate = SS[1] + OS;
        ZMuon = SS[1];
        WMuon = SS[0];
      }

      snu::KMuon ZMuon_leading, ZMuon_subleading;
      if( ZMuon.Pt() > OS.Pt() ){
        ZMuon_leading = ZMuon;
        ZMuon_subleading = OS;
      }
      else{
        ZMuon_leading = OS;
        ZMuon_subleading = ZMuon;
      }

      bool ZMuonPtCut = (ZMuon.Pt() > 20.) || (OS.Pt() > 20.);
      bool isZresonance = (fabs(Z_candidate.M()-m_Z) < 10.);
      bool PtCutOnWMuon = (WMuon.Pt() > 20.);
      bool METCut = (MET > 30.);
      bool mlllCut = ((SS[0]+SS[1]+OS).M() > 100.);
      bool mll4 = (m_dimuon[0] < 4.) || (m_dimuon[1] < 4.);
      bool electronveto = (electronLooseColl.size() == 0);
      bool bjetveto = (n_bjets == 0);

      FillUpDownHist("m_Z_candidate_before_cut_WZ", Z_candidate.M(), weight*FR_reweight, weight_err, 0., 150., 150);
      FillUpDownHist("m_lll_before_cut_WZ", (SS[0]+SS[1]+OS).M(), weight*FR_reweight, weight_err, 0., 500., 500);
      FillUpDownHist("PFMET_before_cut_WZ", MET, weight*FR_reweight, weight_err, 0., 500., 500);
      FillUpDownHist("n_electrons_before_cut_WZ", electronLooseColl.size(), weight*FR_reweight, weight_err, 0., 10., 10);
      FillUpDownHist("n_bjets_before_cut_WZ", n_bjets, weight*FR_reweight, weight_err, 0., 10., 10);
      FillUpDownHist("n_vertices_before_cut_WZ", eventbase->GetEvent().nVertices(), weight*FR_reweight, weight_err, 0., 50., 50);

      //==== N-1 plots
      vector<bool> WZ_cuts;
      WZ_cuts.push_back( fabs(Z_candidate.M()-m_Z) < 15. );
      WZ_cuts.push_back( (SS[0]+SS[1]+OS).M() > 100. );
      WZ_cuts.push_back( n_bjets == 0 );
      WZ_cuts.push_back( MET > 30. );
      if( ZMuonPtCut && PtCutOnWMuon && !mll4 && electronveto ){
        FillHist("N1_preselection_WZ", 0, weight*pileup_reweight, 0., 1., 1);
        for(unsigned int i=0; i<WZ_cuts.size(); i++){
          bool this_bool = true;
          for(unsigned int j=0; j<WZ_cuts.size(); j++){
            if(j==i) continue;
            if(!WZ_cuts.at(j)) this_bool = false;
          }
          if(this_bool){
            if(i==0) FillUpDownHist("N1_Z_mass_WZ", fabs(Z_candidate.M()-m_Z), weight*FR_reweight, weight_err, 0., 60., 60);
            if(i==1) FillUpDownHist("N1_mlll_WZ",  (SS[0]+SS[1]+OS).M(), weight*FR_reweight, weight_err, 0., 500., 25);
            if(i==2) FillUpDownHist("N1_n_bjets_WZ", n_bjets, weight*FR_reweight, weight_err, 0., 4., 4);
            if(i==3) FillUpDownHist("N1_PFMET_WZ", MET, weight*FR_reweight, weight_err, 0., 200., 20);
          }
        }
      }

      if( ZMuonPtCut && isZresonance && PtCutOnWMuon && METCut && mlllCut && !mll4 && electronveto && bjetveto ){
        TString this_suffix = "WZ";
        snu::KParticle nu;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
        snu::KParticle W_candidate = nu+WMuon;

        FillUpDownHist("n_events_"+this_suffix+"", 0, weight*FR_reweight, weight_err, 0., 1., 1);
        FillUpDownHist("n_vertices_"+this_suffix+"", eventbase->GetEvent().nVertices(), weight*FR_reweight, weight_err, 0., 50., 50);
        FillUpDownHist("n_jets_"+this_suffix+"", n_jets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("n_bjets_"+this_suffix+"", n_bjets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("PFMET_"+this_suffix+"", MET, weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("PFMET_phi_"+this_suffix+"", METphi, weight*FR_reweight, weight_err, -3.2, 3.2, 64);
        FillUpDownHist("osllmass_"+this_suffix+"", m_dimuon[0], weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("osllmass_"+this_suffix+"", m_dimuon[1], weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("m_Z_candidate_"+this_suffix+"", Z_candidate.M(), weight*FR_reweight, weight_err, 0., 150., 150);
        FillUpDownHist("mt_W_candidate_"+this_suffix+"", MT(nu,WMuon), weight*FR_reweight, weight_err, 0., 300., 300);
        FillUpDownHist("m_lll_"+this_suffix+"", (SS[0]+SS[1]+OS).M(), weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("Z_candidate_Pt_"+this_suffix+"", Z_candidate.Pt(), weight*FR_reweight, weight_err, 0., 400., 400);
        FillUpDownHist("W_candidate_Pt_"+this_suffix+"", W_candidate.Pt(), weight*FR_reweight, weight_err, 0., 400., 400);
        FillUpDownHist("n_electron_"+this_suffix+"", electronColl.size(), weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("dRZMuonWMuon_"+this_suffix+"", ZMuon.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);
        FillUpDownHist("dRZMuonWMuon_"+this_suffix+"", OS.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);
        FillUpDownHist("dRMETWMuon_"+this_suffix+"", nu.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);

        FillUpDownHist("leadingLepton_Pt_"+this_suffix+"", lep[0].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("leadingLepton_Eta_"+this_suffix+"", lep[0].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("leadingLepton_RelIso_"+this_suffix+"", lep[0].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("leadingLepton_Chi2_"+this_suffix+"", lep[0].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);
        FillUpDownHist("secondLepton_Pt_"+this_suffix+"", lep[1].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("secondLepton_Eta_"+this_suffix+"", lep[1].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("secondLepton_RelIso_"+this_suffix+"", lep[1].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("secondLepton_Chi2_"+this_suffix+"", lep[1].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);
        FillUpDownHist("thirdLepton_Pt_"+this_suffix+"", lep[2].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("thirdLepton_Eta_"+this_suffix+"", lep[2].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("thirdLepton_RelIso_"+this_suffix+"", lep[2].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("thirdLepton_Chi2_"+this_suffix+"", lep[2].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);

        FillUpDownHist("ZMuon_leading_Pt_"+this_suffix+"", ZMuon_leading.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("ZMuon_leading_Eta_"+this_suffix+"", ZMuon_leading.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("ZMuon_leading_RelIso_"+this_suffix+"", ZMuon_leading.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("ZMuon_leading_Chi2_"+this_suffix+"", ZMuon_leading.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("ZMuon_leading_dXY_"+this_suffix+"", fabs(ZMuon_leading.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);
        FillUpDownHist("ZMuon_subleading_Pt_"+this_suffix+"", ZMuon_subleading.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("ZMuon_subleading_Eta_"+this_suffix+"", ZMuon_subleading.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("ZMuon_subleading_RelIso_"+this_suffix+"", ZMuon_subleading.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("ZMuon_subleading_Chi2_"+this_suffix+"", ZMuon_subleading.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);
        FillUpDownHist("ZMuon_subleading_dXY_"+this_suffix+"", fabs(ZMuon_subleading.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);
        FillUpDownHist("WMuon_Pt_"+this_suffix+"", WMuon.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("WMuon_Eta_"+this_suffix+"", WMuon.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("WMuon_RelIso_"+this_suffix+"", WMuon.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("WMuon_Chi2_"+this_suffix+"", WMuon.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("WMuon_dXY_"+this_suffix+"", fabs(WMuon.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);

      } // Z Resonance

      //==== Z+Jets selection
      if( ZMuonPtCut && isZresonance && (MET < 20.) && mlllCut && !mll4 && electronveto && bjetveto ){
        TString this_suffix = "ZJets";
        snu::KParticle nu;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
        snu::KParticle W_candidate = nu+WMuon;

        FillUpDownHist("n_events_"+this_suffix+"", 0, weight*FR_reweight, weight_err, 0., 1., 1);
        FillUpDownHist("n_vertices_"+this_suffix+"", eventbase->GetEvent().nVertices(), weight*FR_reweight, weight_err, 0., 50., 50);
        FillUpDownHist("n_jets_"+this_suffix+"", n_jets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("n_bjets_"+this_suffix+"", n_bjets, weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("PFMET_"+this_suffix+"", MET, weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("PFMET_phi_"+this_suffix+"", METphi, weight*FR_reweight, weight_err, -3.2, 3.2, 64);
        FillUpDownHist("osllmass_"+this_suffix+"", m_dimuon[0], weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("osllmass_"+this_suffix+"", m_dimuon[1], weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("m_Z_candidate_"+this_suffix+"", Z_candidate.M(), weight*FR_reweight, weight_err, 0., 150., 150);
        FillUpDownHist("mt_W_candidate_"+this_suffix+"", MT(nu, WMuon), weight*FR_reweight, weight_err, 0., 300., 300);
        FillUpDownHist("m_lll_"+this_suffix+"", (SS[0]+SS[1]+OS).M(), weight*FR_reweight, weight_err, 0., 500., 500);
        FillUpDownHist("Z_candidate_Pt_"+this_suffix+"", Z_candidate.Pt(), weight*FR_reweight, weight_err, 0., 400., 400);
        FillUpDownHist("W_candidate_Pt_"+this_suffix+"", W_candidate.Pt(), weight*FR_reweight, weight_err, 0., 400., 400);
        FillUpDownHist("n_electron_"+this_suffix+"", electronColl.size(), weight*FR_reweight, weight_err, 0., 10., 10);
        FillUpDownHist("dRZMuonWMuon_"+this_suffix+"", ZMuon.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);
        FillUpDownHist("dRZMuonWMuon_"+this_suffix+"", OS.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);
        FillUpDownHist("dRMETWMuon_"+this_suffix+"", nu.DeltaR(WMuon), weight*FR_reweight, weight_err, 0., 6., 60);

        FillUpDownHist("leadingLepton_Pt_"+this_suffix+"", lep[0].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("leadingLepton_Eta_"+this_suffix+"", lep[0].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("leadingLepton_RelIso_"+this_suffix+"", lep[0].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("leadingLepton_Chi2_"+this_suffix+"", lep[0].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("secondLepton_Pt_"+this_suffix+"", lep[1].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("secondLepton_Eta_"+this_suffix+"", lep[1].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("secondLepton_RelIso_"+this_suffix+"", lep[1].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("secondLepton_Chi2_"+this_suffix+"", lep[1].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("thirdLepton_Pt_"+this_suffix+"", lep[2].Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("thirdLepton_Eta_"+this_suffix+"", lep[2].Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("thirdLepton_RelIso_"+this_suffix+"", lep[2].RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("thirdLepton_Chi2_"+this_suffix+"", lep[2].GlobalChi2() , weight*FR_reweight, weight_err, 0., 10, 100);

        FillUpDownHist("ZMuon_leading_Pt_"+this_suffix+"", ZMuon_leading.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("ZMuon_leading_Eta_"+this_suffix+"", ZMuon_leading.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("ZMuon_leading_RelIso_"+this_suffix+"", ZMuon_leading.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("ZMuon_leading_Chi2_"+this_suffix+"", ZMuon_leading.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("ZMuon_leading_dXY_"+this_suffix+"", fabs(ZMuon_leading.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);
        FillUpDownHist("ZMuon_subleading_Pt_"+this_suffix+"", ZMuon_subleading.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("ZMuon_subleading_Eta_"+this_suffix+"", ZMuon_subleading.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("ZMuon_subleading_RelIso_"+this_suffix+"", ZMuon_subleading.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("ZMuon_subleading_Chi2_"+this_suffix+"", ZMuon_subleading.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("ZMuon_subleading_dXY_"+this_suffix+"", fabs(ZMuon_subleading.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);
        FillUpDownHist("WMuon_Pt_"+this_suffix+"", WMuon.Pt() , weight*FR_reweight, weight_err, 0., 200., 200);
        FillUpDownHist("WMuon_Eta_"+this_suffix+"", WMuon.Eta() , weight*FR_reweight, weight_err, -3., 3., 60);
        FillUpDownHist("WMuon_RelIso_"+this_suffix+"", WMuon.RelIso04() , weight*FR_reweight, weight_err, 0., 1.0, 100);
        FillUpDownHist("WMuon_Chi2_"+this_suffix+"", WMuon.GlobalChi2() , weight*FR_reweight, weight_err, 0., 10., 100);
        FillUpDownHist("WMuon_dXY_"+this_suffix+"", fabs(WMuon.dXY()) , weight*FR_reweight, weight_err, 0., 0.01, 100);

      }

    } // Not All Same Charge

  } // isThreeMuon


 
  return;

}// End of execute event loop
  


void trilepton_mumumu_CR_FR_method::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_CR_FR_method::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_CR_FR_method::~trilepton_mumumu_CR_FR_method() {
  
  Message("In trilepton_mumumu_CR_FR_method Destructor" , INFO);
  
}


void trilepton_mumumu_CR_FR_method::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_CR_FR_method::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_CR_FR_method::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_CR_FR_methodCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_mumumu_CR_FR_method::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



