// $Id: Validation_trilepton.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQValidation_trilepton Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "Validation_trilepton.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (Validation_trilepton);

Validation_trilepton::Validation_trilepton() :  AnalyzerCore(), out_muons(0) {
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("Validation_trilepton");
  
  Message("In Validation_trilepton constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void Validation_trilepton::InitialiseAnalysis() throw( LQError ) {
  
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


void Validation_trilepton::ExecuteEvents()throw( LQError ){

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

  if(!PassTriggerOR(triggerlist)) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

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

  //=====================
  //==== Get MET,METPhi
  //=====================

  snu::KEvent Evt = eventbase->GetEvent();

  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);
  int n_jets = jetColl_hn.size();
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)){
      n_bjets++;
      FillHist("bjet_pt", jetColl_hn.at(j).Pt(), 1., 0., 200., 200);
    }
  }

  //==========================
  //==== MuonIDs to validate
  //==========================

  vector<TString> muonids;
  muonids.push_back("MUON_HN_TRI_TIGHT");
  muonids.push_back("MUON_POG_TIGHT"); 

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight_CATTools = 1.0;
  float pileup_reweight_John = 1.0;
  
  if(!k_isdata){
    //==== CATTools reweight
    pileup_reweight_CATTools = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    pileup_reweight_John = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  numberVertices = eventbase->GetEvent().nVertices();

  for(unsigned int i_id=0; i_id<muonids.size(); i_id++){

    double this_weight = weight;

    TString muonid = muonids.at(i_id);

    std::vector<snu::KMuon> muoncoll = GetMuons(muonid);
    int n_triTight_muons = muoncoll.size();

    //==================
    //==== Correct MET
    //==================

    double MET = Evt.MET();
    double METphi = Evt.METPhi();
    CorrectedMETRochester(muoncoll, MET, METphi);

    //=========================== 
    //==== Get Muon Corrections
    //===========================

    double muon_id_sf = mcdata_correction->MuonScaleFactor(muonid, muoncoll, 0);
    double muon_iso_sf = 1.;
    if(muonid=="MUON_POG_TIGHT") mcdata_correction->MuonISOScaleFactor(muonid,muoncoll, 0);
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muoncoll);

    if(!isData && !k_running_nonprompt){
      this_weight*=muon_id_sf;
      this_weight*=muon_iso_sf;
      this_weight*=trigger_ps_weight;
      this_weight*=MuTrkEffSF;
    }

    if( muoncoll.size() == 2 ){

      std::map< TString, bool > map_whichCR_to_isCR;

      snu::KMuon lep[2];
      lep[0] = muoncoll.at(0);
      lep[1] = muoncoll.at(1);

      bool leadPt20 = muoncoll.at(0).Pt() > 20.;
      bool isSS = muoncoll.at(0).Charge() == muoncoll.at(1).Charge();

      double m_Z = 91.1876;
      double m_dimuon = ( muoncoll.at(0) + muoncoll.at(1) ).M();
      bool ZResonance = fabs(m_dimuon-m_Z) < 10.;

      map_whichCR_to_isCR[muonid+"_DiMuon"] = leadPt20;
      map_whichCR_to_isCR[muonid+"_SSDiMuon"] = leadPt20 && isSS;
      map_whichCR_to_isCR[muonid+"_OSDiMuon"] = leadPt20 && !isSS;
      map_whichCR_to_isCR[muonid+"_OSDiMuon_Z_10GeV"] = leadPt20 && !isSS && ZResonance;

      double trigger_sf = 1.;
      if(!k_isdata && muonid=="MUON_HN_TRI_TIGHT"){
        std::vector<snu::KElectron> empty_el;
        double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(empty_el, "ELECTRON_HN_TIGHTv4", muoncoll, "MUON_HN_TRI_TIGHT", 0, 0, 0);
        double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(empty_el, "ELECTRON_HN_TIGHTv4", muoncoll, "MUON_HN_TRI_TIGHT", 0, 1, 0);
        trigger_sf = trigger_eff_Data/trigger_eff_MC;
      }

      for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
        TString this_suffix = it->first;
        if(it->second){

          //==== CATTools PU reweight

          FillHist("weight_"+this_suffix+"_PUCATTools", this_weight*pileup_reweight_CATTools, 1., -1., 1., 1000);
          FillHist("n_events_"+this_suffix+"_PUCATTools", 0, this_weight*pileup_reweight_CATTools, 0., 1., 1);
          FillHist("n_jets_"+this_suffix+"_PUCATTools", n_jets, this_weight*pileup_reweight_CATTools, 0., 10., 10);
          FillHist("n_bjets_"+this_suffix+"_PUCATTools", n_bjets, this_weight*pileup_reweight_CATTools, 0., 10., 10);
          FillHist("PFMET_"+this_suffix+"_PUCATTools", MET, this_weight*pileup_reweight_CATTools, 0., 500., 500);
          FillHist("mll_"+this_suffix+"_PUCATTools", m_dimuon , this_weight*pileup_reweight_CATTools, 0., 200., 200);
          FillHist("n_vertices_"+this_suffix+"_PUCATTools", numberVertices, this_weight*pileup_reweight_CATTools, 0., 50., 50);
          FillHist("leadingLepton_Pt_"+this_suffix+"_PUCATTools", lep[0].Pt() , this_weight*pileup_reweight_CATTools, 0., 200., 200);
          FillHist("leadingLepton_Eta_"+this_suffix+"_PUCATTools", lep[0].Eta() , this_weight*pileup_reweight_CATTools, -3., 3., 60);
          FillHist("leadingLepton_RelIso_"+this_suffix+"_PUCATTools", lep[0].RelIso04() , this_weight*pileup_reweight_CATTools, 0., 1.0, 100);
          FillHist("leadingLepton_Chi2_"+this_suffix+"_PUCATTools", lep[0].GlobalChi2() , this_weight*pileup_reweight_CATTools, 0., 10., 100);
          FillHist("leadingLepton_dXY_"+this_suffix+"_PUCATTools"+"", fabs(lep[0].dXY()) , this_weight*pileup_reweight_CATTools, 0., 0.1, 100);
          FillHist("leadingLepton_dXYSig_"+this_suffix+"_PUCATTools"+"", fabs(lep[0].dXYSig()) , this_weight*pileup_reweight_CATTools, 0., 4., 40);
          FillHist("secondLepton_Pt_"+this_suffix+"_PUCATTools", lep[1].Pt() , this_weight*pileup_reweight_CATTools, 0., 200., 200);
          FillHist("secondLepton_Eta_"+this_suffix+"_PUCATTools", lep[1].Eta() , this_weight*pileup_reweight_CATTools, -3., 3., 60);
          FillHist("secondLepton_RelIso_"+this_suffix+"_PUCATTools", lep[1].RelIso04() , this_weight*pileup_reweight_CATTools, 0., 1.0, 100);
          FillHist("secondLepton_Chi2_"+this_suffix+"_PUCATTools", lep[1].GlobalChi2() , this_weight*pileup_reweight_CATTools, 0., 10., 100);
          FillHist("secondLepton_dXY_"+this_suffix+"_PUCATTools"+"", fabs(lep[1].dXY()) , this_weight*pileup_reweight_CATTools, 0., 0.1, 100);
          FillHist("secondLepton_dXYSig_"+this_suffix+"_PUCATTools"+"", fabs(lep[1].dXYSig()) , this_weight*pileup_reweight_CATTools, 0., 4., 40);

          //==== John PU reweight

          FillHist("weight_"+this_suffix+"_PUJohn", this_weight*pileup_reweight_John, 1., -1., 1., 1000);
          FillHist("n_events_"+this_suffix+"_PUJohn", 0, this_weight*pileup_reweight_John, 0., 1., 1);
          FillHist("n_jets_"+this_suffix+"_PUJohn", n_jets, this_weight*pileup_reweight_John, 0., 10., 10);
          FillHist("n_bjets_"+this_suffix+"_PUJohn", n_bjets, this_weight*pileup_reweight_John, 0., 10., 10);
          FillHist("PFMET_"+this_suffix+"_PUJohn", MET, this_weight*pileup_reweight_John, 0., 500., 500);
          FillHist("mll_"+this_suffix+"_PUJohn", m_dimuon , this_weight*pileup_reweight_John, 0., 200., 200);
          FillHist("n_vertices_"+this_suffix+"_PUJohn", numberVertices, this_weight*pileup_reweight_John, 0., 50., 50);
          FillHist("leadingLepton_Pt_"+this_suffix+"_PUJohn", lep[0].Pt() , this_weight*pileup_reweight_John, 0., 200., 200);
          FillHist("leadingLepton_Eta_"+this_suffix+"_PUJohn", lep[0].Eta() , this_weight*pileup_reweight_John, -3., 3., 60);
          FillHist("leadingLepton_RelIso_"+this_suffix+"_PUJohn", lep[0].RelIso04() , this_weight*pileup_reweight_John, 0., 1.0, 100);
          FillHist("leadingLepton_Chi2_"+this_suffix+"_PUJohn", lep[0].GlobalChi2() , this_weight*pileup_reweight_John, 0., 10., 100);
          FillHist("leadingLepton_dXY_"+this_suffix+"_PUJohn"+"", fabs(lep[0].dXY()) , this_weight*pileup_reweight_John, 0., 0.1, 100);
          FillHist("leadingLepton_dXYSig_"+this_suffix+"_PUJohn"+"", fabs(lep[0].dXYSig()) , this_weight*pileup_reweight_John, 0., 4., 40);
          FillHist("secondLepton_Pt_"+this_suffix+"_PUJohn", lep[1].Pt() , this_weight*pileup_reweight_John, 0., 200., 200);
          FillHist("secondLepton_Eta_"+this_suffix+"_PUJohn", lep[1].Eta() , this_weight*pileup_reweight_John, -3., 3., 60);
          FillHist("secondLepton_RelIso_"+this_suffix+"_PUJohn", lep[1].RelIso04() , this_weight*pileup_reweight_John, 0., 1.0, 100);
          FillHist("secondLepton_Chi2_"+this_suffix+"_PUJohn", lep[1].GlobalChi2() , this_weight*pileup_reweight_John, 0., 10., 100);
          FillHist("secondLepton_dXY_"+this_suffix+"_PUJohn"+"", fabs(lep[1].dXY()) , this_weight*pileup_reweight_John, 0., 0.1, 100);
          FillHist("secondLepton_dXYSig_"+this_suffix+"_PUJohn"+"", fabs(lep[1].dXYSig()) , this_weight*pileup_reweight_John, 0., 4., 40);

          //==== CATTools PU reweight Trigger Scale Factor

          FillHist("weight_"+this_suffix+"_PUCATTools_TriggerSF", this_weight*trigger_sf, 1., -1., 1., 1000);
          FillHist("n_events_"+this_suffix+"_PUCATTools_TriggerSF", 0, this_weight*trigger_sf, 0., 1., 1);
          FillHist("n_jets_"+this_suffix+"_PUCATTools_TriggerSF", n_jets, this_weight*trigger_sf, 0., 10., 10);
          FillHist("n_bjets_"+this_suffix+"_PUCATTools_TriggerSF", n_bjets, this_weight*trigger_sf, 0., 10., 10);
          FillHist("PFMET_"+this_suffix+"_PUCATTools_TriggerSF", MET, this_weight*trigger_sf, 0., 500., 500);
          FillHist("mll_"+this_suffix+"_PUCATTools_TriggerSF", m_dimuon , this_weight*trigger_sf, 0., 200., 200);
          FillHist("n_vertices_"+this_suffix+"_PUCATTools_TriggerSF", numberVertices, this_weight*trigger_sf, 0., 50., 50);
          FillHist("leadingLepton_Pt_"+this_suffix+"_PUCATTools_TriggerSF", lep[0].Pt() , this_weight*trigger_sf, 0., 200., 200);
          FillHist("leadingLepton_Eta_"+this_suffix+"_PUCATTools_TriggerSF", lep[0].Eta() , this_weight*trigger_sf, -3., 3., 60);
          FillHist("leadingLepton_RelIso_"+this_suffix+"_PUCATTools_TriggerSF", lep[0].RelIso04() , this_weight*trigger_sf, 0., 1.0, 100);
          FillHist("leadingLepton_Chi2_"+this_suffix+"_PUCATTools_TriggerSF", lep[0].GlobalChi2() , this_weight*trigger_sf, 0., 10., 100);
          FillHist("leadingLepton_dXY_"+this_suffix+"_PUCATTools_TriggerSF"+"", fabs(lep[0].dXY()) , this_weight*trigger_sf, 0., 0.1, 100);
          FillHist("leadingLepton_dXYSig_"+this_suffix+"_PUCATTools_TriggerSF"+"", fabs(lep[0].dXYSig()) , this_weight*trigger_sf, 0., 4., 40);
          FillHist("secondLepton_Pt_"+this_suffix+"_PUCATTools_TriggerSF", lep[1].Pt() , this_weight*trigger_sf, 0., 200., 200);
          FillHist("secondLepton_Eta_"+this_suffix+"_PUCATTools_TriggerSF", lep[1].Eta() , this_weight*trigger_sf, -3., 3., 60);
          FillHist("secondLepton_RelIso_"+this_suffix+"_PUCATTools_TriggerSF", lep[1].RelIso04() , this_weight*trigger_sf, 0., 1.0, 100);
          FillHist("secondLepton_Chi2_"+this_suffix+"_PUCATTools_TriggerSF", lep[1].GlobalChi2() , this_weight*trigger_sf, 0., 10., 100);
          FillHist("secondLepton_dXY_"+this_suffix+"_PUCATTools_TriggerSF"+"", fabs(lep[1].dXY()) , this_weight*trigger_sf, 0., 0.1, 100);
          FillHist("secondLepton_dXYSig_"+this_suffix+"_PUCATTools_TriggerSF"+"", fabs(lep[1].dXYSig()) , this_weight*trigger_sf, 0., 4., 40);

        }
      } 

    } // isTwoMuon

  }

  return;

}// End of execute event loop
  


void Validation_trilepton::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Validation_trilepton::BeginCycle() throw( LQError ){
  
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

Validation_trilepton::~Validation_trilepton() {
  
  Message("In Validation_trilepton Destructor" , INFO);
  
}


void Validation_trilepton::FillCutFlow(TString cut, float weight){

  
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


void Validation_trilepton::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void Validation_trilepton::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this Validation_trileptonCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void Validation_trilepton::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}




