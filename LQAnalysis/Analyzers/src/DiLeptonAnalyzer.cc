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


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DiLeptonAnalyzer);

DiLeptonAnalyzer::DiLeptonAnalyzer() :  AnalyzerCore(), MET(-999), METphi(-999), nbjets(-999), nbjets_fwd(-999), nbjets_nolepveto(-999), n_vtx(-999)
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
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);

  //====================
  //==== [CUT] Trigger //if(!PassTriggerOR(triggerlist_DiMuon)) return;
  //====================

  std::vector<TString> triggerlist_DiMuon;
  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  std::vector<TString> triggerlist_DiElectron;
  triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  //triggerlist_DiElectron.push_back("HLT_Ele27_WPTight_Gsf_v");

  std::vector<TString> triggerlist_EMu;
  triggerlist_EMu.push_back("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v");

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
  std::vector< snu::KMuon > muons = GetMuons("MUON_HN_TRI_LOOSE", false);
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muons, 0);
  double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muons);

  //==== Electrons
  std::vector< snu::KElectron > electrons = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrons, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);

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
  int n_triLoose_muons = muons.size();
  int n_triTight_muons(0);
  for(unsigned int i=0; i<muons.size(); i++){
    if(PassID(muons.at(i), "MUON_HN_TRI_TIGHT")) n_triTight_muons++;
  }

  int n_triLoose_electrons = electrons.size();
  int n_triTight_electrons(0);
  for(unsigned int i=0; i<electrons.size(); i++){
    if(PassID(electrons.at(i), "ELECTRON_HN_TIGHT")) n_triTight_electrons++;
  }

  int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
  int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

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

  if(n_triLoose_leptons != 2) return;
  //if(jets.size() <2) return;

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight=(1.0);
  if(!isData){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  if(!isData && !k_running_nonprompt){
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
  std::vector< bool > isTTs, isLOOSEs;

  bool RunningNonPromptData = k_running_nonprompt && isData;

  Suffixs.push_back("DiMuon");
  Triggers.push_back(triggerlist_DiMuon);
  isTTs.push_back( isTwoMuon_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isTwoMuon_Loose && RunningNonPromptData );

  Suffixs.push_back("DiElectron");
  Triggers.push_back(triggerlist_DiElectron);
  isTTs.push_back( isTwoElectron_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isTwoElectron_Loose && RunningNonPromptData );
/*
  Suffixs.push_back("EMu");
  Triggers.push_back(triggerlist_EMu);
  isTTs.push_back( isEMu_TT && !RunningNonPromptData );
  isLOOSEs.push_back( isEMu_Loose && RunningNonPromptData );
*/
  for(unsigned int i=0; i<Suffixs.size(); i++){

    TString Suffix = Suffixs.at(i);

    if(!PassTriggerOR( Triggers.at(i) )) continue;
    if(!isTTs.at(i) && !isLOOSEs.at(i)) continue;
    if(isLOOSEs.at(i) && !isData) continue;
    if(isData){
      if(Suffix=="DiMuon"){
        if(k_channel != "DoubleMuon") continue;
      }
      if(Suffix=="DiElectron"){
        if(k_channel != "DoubleEG") continue;
      }
      if(Suffix=="EMu"){
        if(k_channel != "MuonEG") continue;
      }
    }

    double trigger_ps_weight= WeightByTrigger(Triggers.at(i), TargetLumi);
    double this_weight = weight*trigger_ps_weight;

    double trigger_sf = 1.;
    if(!isData && Suffix=="DiMuon"){
      double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, muons, 0, 0, 0);
      double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, muons, 0, 1, 0);
      trigger_sf = trigger_eff_Data/trigger_eff_MC;
    }
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

    //==== Lepton Pt Cut
    if(Suffix=="DiMuon"){
      if(lep.at(0).Pt() < 20. || lep.at(1).Pt() < 10.) continue;
    }
    if(Suffix=="DiElectron"){
/*
      bool PassDiElTrig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      bool PassSingleElTrig = PassTrigger("HLT_Ele27_WPTight_Gsf_v");
      if(PassDiElTrig && !PassSingleElTrig){
        if(lep.at(0).Pt() < 25. || lep.at(1).Pt() < 15.) continue;
      }
      else if(!PassDiElTrig && PassSingleElTrig){
        if(lep.at(0).Pt() < 30.) continue;
      }
      else if(PassDiElTrig && PassSingleElTrig){
        if(lep.at(0).Pt() < 25.) continue;
      }
*/
      if(lep.at(0).Pt() < 25. || lep.at(1).Pt() < 15.) continue;
    }
    if(Suffix=="EMu"){
      //==== lep.at(0) : muon
      //==== lep.at(1) : eletron
      if(lep.at(0).Pt() < 10. || lep.at(1).Pt() < 15.) continue;
    }

    bool isSS = lep.at(0).Charge() == lep.at(1).Charge();
    double m_Z = 91.1876;
    bool isOffZ = fabs( (lep.at(0)+lep.at(1)).M() - m_Z ) > 15.;

    //==== mll Cut Study
    FillHist("CutStudy_m_ll_"+Suffix, ( lep.at(0)+lep.at(1) ).M(), 1., 0., 20., 200);
    bool mll4GeV = ( lep.at(0)+lep.at(1) ).M() < 12.;

    //if(!isSS && !isOffZ) continue;
    if(mll4GeV) continue;

    double this_weight_err(0.);
    if( isLOOSEs.at(i) ){
      //this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muons, "MUON_HN_TRI_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHT", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      //this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TRI_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHT", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
      m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TRI_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHT", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");

      this_weight = m_datadriven_bkg->GetFakeObj()->k_weight;
      this_weight_err = m_datadriven_bkg->GetFakeObj()->k_weight_err;
      
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

    for(std::map< TString, bool >::iterator it = map_Region_to_Bool.begin(); it != map_Region_to_Bool.end(); it++){
      TString this_suffix = it->first;
      if(it->second){

        FillDiLeptonPlot(this_suffix+"_AllCharge", lep, jets, jets_fwd, jets_nolepveto, this_weight);
        if( isLOOSEs.at(i) ){
          FillDiLeptonPlot(this_suffix+"_AllCharge_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
          FillDiLeptonPlot(this_suffix+"_AllCharge_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
        }

        if(isSS){
          FillDiLeptonPlot(this_suffix+"_SS", lep, jets, jets_fwd, jets_nolepveto, this_weight);
          if( isLOOSEs.at(i) ){
            FillDiLeptonPlot(this_suffix+"_SS_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
            FillDiLeptonPlot(this_suffix+"_SS_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
          }
        }
        else{
          FillDiLeptonPlot(this_suffix+"_OS", lep, jets, jets_fwd, jets_nolepveto, this_weight);
          if( isLOOSEs.at(i) ){
            FillDiLeptonPlot(this_suffix+"_OS_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
            FillDiLeptonPlot(this_suffix+"_OS_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
          }
        }

        if(isOffZ){

          FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge", lep, jets, jets_fwd, jets_nolepveto, this_weight);
          if( isLOOSEs.at(i) ){
            FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
            FillDiLeptonPlot(this_suffix+"_OffZ_AllCharge_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
          }

          if(isSS){
            FillDiLeptonPlot(this_suffix+"_OffZ_SS", lep, jets, jets_fwd, jets_nolepveto, this_weight);
            if( isLOOSEs.at(i) ){
              FillDiLeptonPlot(this_suffix+"_OffZ_SS_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
              FillDiLeptonPlot(this_suffix+"_OffZ_SS_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
            }
          }
          else{
            FillDiLeptonPlot(this_suffix+"_OffZ_OS", lep, jets, jets_fwd, jets_nolepveto, this_weight);
            if( isLOOSEs.at(i) ){
              FillDiLeptonPlot(this_suffix+"_OffZ_OS_up", lep, jets, jets_fwd, jets_nolepveto, this_weight+this_weight_err);
              FillDiLeptonPlot(this_suffix+"_OffZ_OS_down", lep, jets, jets_fwd, jets_nolepveto, this_weight-this_weight_err);
            }
          }

        } // Off-Z



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
  double thisweight
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

}














