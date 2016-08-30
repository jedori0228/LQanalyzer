// $Id: trilepton_mumumu_CR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_CR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_CR.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_CR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu_CR::trilepton_mumumu_CR() :  AnalyzerCore(), out_muons(0) {
  
  rmcor = new rochcor2015();
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_CR");
  
  Message("In trilepton_mumumu_CR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(trilephist,"cut0");
  MakeCleverHistograms(trilephist,"cut0_PU");
  MakeCleverHistograms(trilephist,"cutdR");
  MakeCleverHistograms(trilephist,"cutdR_PU");
  MakeCleverHistograms(trilephist,"cutdR_cutW");
  MakeCleverHistograms(trilephist,"cutdR_cutW_PU");

}


void trilepton_mumumu_CR::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_mumumu_CR::ExecuteEvents()throw( LQError ){

  
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  

  /// Acts on data to remove bad reconstructed event 
  if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;
  

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment

  if(!PassBasicEventCuts()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");

  //float trigger_ps_weight= ApplyPrescale("HLT_IsoMu20", TargetLumi,lumimask);

  if(!PassTrigger(triggerslist, prescale)) return;
  FillCutFlow("TriggerCut", weight);
  // Trigger matching is done using KMuon::TriggerMatched(TString) which returns a bool

  /* // #### CAT::: trigger matching information is stored for muons and electrons for:
  ///HLT_IsoMu24_eta2p1_v
  ///HLT_Mu17_Mu8_DZ_v
  ///HLT_Mu17_TkMu8_DZ_v
  ///HLT_IsoMu20
  ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
  ///HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele12_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
  ///HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_
  ///HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v
  ///HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v
  */

  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;


  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", weight);


  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed
  std::vector<snu::KMuon> muontriTightColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_TIGHT); 
  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons(BaseSelection::MUON_HN_TRI_LOOSE);

  std::vector<snu::KMuon> muontriTightColl, muontriLooseColl;
  for(unsigned int i=0; i<muontriTightColl_raw.size(); i++){
    snu::KMuon thismuon = muontriTightColl_raw.at(i);
    if(isData) muontriTightColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriTightColl.push_back(thismuon);
    }
  }
  for(unsigned int i=0; i<muontriLooseColl_raw.size(); i++){
    snu::KMuon thismuon = muontriLooseColl_raw.at(i);
    if(isData) muontriLooseColl.push_back(thismuon);
    else{
      if( thismuon.GetParticleType() == snu::KMuon::PROMPT ) muontriLooseColl.push_back(thismuon);
    }
  }
   
  CorrectMuonMomentum(muonTightColl);
  float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonTightColl,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl             = GetJets(BaseSelection::JET_NOLEPTONVETO); // All jets
  std::vector<snu::KJet> jetColl_loose       = GetJets(BaseSelection::JET_LOOSE); // pt > 10; eta < 5. ; PFlep veto
  std::vector<snu::KJet> jetColl_tight       = GetJets(BaseSelection::JET_TIGHT);// pt > 20 ; eta < 2.5; PFlep veto
  std::vector<snu::KJet> jetColl_hn          = GetJets(BaseSelection::JET_HN);// pt > 20 ; eta < 2.5; PFlep veto; pileup ID
   
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  numberVertices = eventbase->GetEvent().nVertices();   

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    pileup_reweight = eventbase->GetEvent().PileUpWeight(lumimask);

  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);


  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    weight*=pileup_reweight;
    //weight*=weight_trigger_sf;
    //weight*=trigger_ps_weight;
  }

  int n_triTight_muons = muontriTightColl.size();
  int n_triLoose_muons = muontriLooseColl.size();
  int n_jets = jetColl_loose.size();

  FillHist("GenWeight_NJet" , n_jets*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  // CR related variables //
  FillHist("n_tight_muons_control_PU", n_triTight_muons, weight*pileup_reweight, 0, 10, 10);
  FillHist("n_loose_muons_control_PU", n_triLoose_muons, weight*pileup_reweight, 0, 10, 10);
  if( n_triTight_muons == 2 && n_triLoose_muons == 2){
    int isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ? 1 : 0;
    FillHist("2Muons_OS0_SS1_control_PU", isSS, weight*pileup_reweight, 0, 2, 2);
  }
  FillHist("n_jets_control_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if(jetColl_loose.at(j).IsBTagged(snu::KJet::cMVAv2, snu::KJet::Medium)) n_bjets++;
  }
  FillHist("n_bjets_control_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);

  //================
  //==== define CR
  //================

  std::map< TString, bool > map_whichCR_to_isCR;

  int n_muons(0);

  bool isTwoMuon   = n_triLoose_muons == 2 && n_triTight_muons == 2;
  bool isThreeMuon = n_triLoose_muons == 3 && n_triTight_muons == 3;

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();

  //==== CR with Two Muons
  if(isTwoMuon){
    snu::KMuon lep[2];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);

    bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
    bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
    double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();

    map_whichCR_to_isCR["SS_0jet_vetoLowRes"] = leadPt20 && isSS && (n_jets == 0) && (m_dimuon > 15.);
    map_whichCR_to_isCR["SS_AL1bjet_vetoLowRes"] = leadPt20 && isSS && (n_bjets > 0) && (m_dimuon > 15.);
    map_whichCR_to_isCR["SS_AL2bjet_vetoLowRes"] = leadPt20 && isSS && (n_bjets > 1) && (m_dimuon > 15.);

    for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
      TString this_suffix = it->first;
      if(it->second){
				FillHist("n_events_"+this_suffix+"_PU", 0, weight*pileup_reweight, 0, 1, 1);
				FillHist("n_jets_"+this_suffix+"_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
				FillHist("n_bjets_"+this_suffix+"_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);
				FillHist("PFMET_"+this_suffix+"_PU", MET, weight*pileup_reweight, 0, 500, 500);
				FillHist("mll_"+this_suffix+"_PU", m_dimuon , weight*pileup_reweight, 0, 500, 500);
				FillHist("leadingLepton_Pt_"+this_suffix+"_PU", lep[0].Pt() , weight*pileup_reweight, 0, 200, 200);
				FillHist("leadingLepton_Eta_"+this_suffix+"_PU", lep[0].Eta() , weight*pileup_reweight, -3, 3, 60);
				FillHist("leadingLepton_RelIso_"+this_suffix+"_PU", lep[0].RelIso04() , weight*pileup_reweight, 0, 1.0, 100);
				FillHist("leadingLepton_Chi2_"+this_suffix+"_PU", lep[0].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
				FillHist("secondLepton_Pt_"+this_suffix+"_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
				FillHist("secondLepton_Eta_"+this_suffix+"_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
				FillHist("secondLepton_RelIso_"+this_suffix+"_PU", lep[1].RelIso04() , weight*pileup_reweight, 0, 1.0, 100);
				FillHist("secondLepton_Chi2_"+this_suffix+"_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
      }
    } 

  } // isTwoMuon

  //==== CR with Three Muons
  if(isThreeMuon){
    snu::KMuon lep[3];
    lep[0] = muontriLooseColl.at(0);
    lep[1] = muontriLooseColl.at(1);
    lep[2] = muontriLooseColl.at(2);

    bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
    bool AllSameCharge = ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) &&
                         ( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(2).Charge() );

    if( leadPt20 && !AllSameCharge ){
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
      snu::KMuon ExtraMuon;
      if( fabs(m_dimuon[0]-m_Z) < fabs(m_dimuon[1]-m_Z) ){
        Z_candidate = SS[0] + OS;
        ExtraMuon = SS[1];
      }
      else{
        Z_candidate = SS[1] + OS;
        ExtraMuon = SS[0];
      }

      bool isZresonance = fabs(Z_candidate.M()-m_Z) < 20.;
      bool PtCutOnExtraMuon = ExtraMuon.Pt() > 20.;
      bool METCut = MET > 20.;

      if( isZresonance && PtCutOnExtraMuon && METCut ){
        TString this_suffix = "WZ";
        snu::KParticle nu;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET); 
        snu::KParticle W_candidate = nu+ExtraMuon;

        FillHist("n_events_"+this_suffix+"_PU", 0, weight*pileup_reweight, 0, 1, 1);
        FillHist("n_jets_"+this_suffix+"_PU", n_jets, weight*pileup_reweight, 0, 10, 10);
        FillHist("n_bjets_"+this_suffix+"_PU", n_bjets, weight*pileup_reweight, 0, 10, 10);
        FillHist("PFMET_"+this_suffix+"_PU", MET, weight*pileup_reweight, 0, 500, 500);
        FillHist("m_Z_candidate_"+this_suffix+"_PU", Z_candidate.M(), weight*pileup_reweight, 0, 150, 150);
        FillHist("mt_W_candidate_"+this_suffix+"_PU", W_candidate.Mt(), weight*pileup_reweight, 0, 300, 300);

        FillHist("leadingLepton_Pt_"+this_suffix+"_PU", lep[0].Pt() , weight*pileup_reweight, 0, 200, 200);
        FillHist("leadingLepton_Eta_"+this_suffix+"_PU", lep[0].Eta() , weight*pileup_reweight, -3, 3, 60);
        FillHist("leadingLepton_RelIso_"+this_suffix+"_PU", lep[0].RelIso04() , weight*pileup_reweight, 0, 1.0, 100);
        FillHist("leadingLepton_Chi2_"+this_suffix+"_PU", lep[0].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
        FillHist("secondLepton_Pt_"+this_suffix+"_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
        FillHist("secondLepton_Eta_"+this_suffix+"_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
        FillHist("secondLepton_RelIso_"+this_suffix+"_PU", lep[1].RelIso04() , weight*pileup_reweight, 0, 1.0, 100);
        FillHist("secondLepton_Chi2_"+this_suffix+"_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);
        FillHist("thirdLepton_Pt_"+this_suffix+"_PU", lep[2].Pt() , weight*pileup_reweight, 0, 200, 200);
        FillHist("thirdLepton_Eta_"+this_suffix+"_PU", lep[2].Eta() , weight*pileup_reweight, -3, 3, 60);
        FillHist("thirdLepton_RelIso_"+this_suffix+"_PU", lep[2].RelIso04() , weight*pileup_reweight, 0, 1.0, 100);
        FillHist("thirdLepton_Chi2_"+this_suffix+"_PU", lep[2].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);

      } // Z Resonance
 

    } // Not All Same Charge

  } // isThreeMuon


 

/*
  if(n_muons==3){
    FillHist("TT_thirdLepton_Pt_control_PU", lep[1].Pt() , weight*pileup_reweight, 0, 200, 200);
    FillHist("TT_thirdLepton_Eta_control_PU", lep[1].Eta() , weight*pileup_reweight, -3, 3, 60);
    FillHist("TT_thirdLepton_RelIso_control_PU", LeptonRelIso[1] , weight*pileup_reweight, 0, 1.0, 100);
    FillHist("TT_thirdLepton_Chi2_control_PU", lep[1].GlobalChi2() , weight*pileup_reweight, 0, 10, 100);

    FillHist("TT_thirdLepton_Pt_1_control_PU", lep[2].Pt(), 1., 0, 200, 200);
    FillHist("TT_thirdLepton_Eta_1_control_PU", lep[2].Eta(), 1., -3, 3, 60);
    FillHist("TT_thirdLepton_RelIso_1_control_PU", LeptonRelIso[2], 1., 0, 1.0, 100);
    FillHist("TT_thirdLepton_Chi2_1_control_PU", lep[2].GlobalChi2(), 1., 0, 10, 100);
  }
*/




  return;

}// End of execute event loop
  


void trilepton_mumumu_CR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_CR::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

trilepton_mumumu_CR::~trilepton_mumumu_CR() {
  
  Message("In trilepton_mumumu_CR Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_mumumu_CR::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_CR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_CR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_CRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_mumumu_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

