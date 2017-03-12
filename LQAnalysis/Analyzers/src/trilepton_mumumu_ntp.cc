// $Id: trilepton_mumumu_ntp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_ntp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_ntp.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_ntp);

trilepton_mumumu_ntp::trilepton_mumumu_ntp() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_ntp");
  
  Message("In trilepton_mumumu_ntp constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_ntp::InitialiseAnalysis() throw( LQError ) {
  
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


void trilepton_mumumu_ntp::ExecuteEvents()throw( LQError ){

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

  //==================================
  //==== Prepare the lowest pt muons
  //==================================

  std::vector<snu::KMuon> muontriVLooseColl_lowestPtCut;
  double this_RelIso = 0.4;
  //==== signal
  if( k_sample_name.Contains("HN") ){

    std::vector<snu::KTruth> truthColl;
    eventbase->GetTruthSel()->Selection(truthColl);
    m_HNgenmatch->SetAllGenParticles(truthColl);
    m_HNgenmatch->SetSignalMass(GetSignalMass());
    m_HNgenmatch->SetHNpdgids(9900012);
    m_HNgenmatch->FindGenParticles();

    //==== save gen particles @ snu::KTruth gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;

    snu::KTruth gen_l_1 = m_HNgenmatch->gen_l_1;
    snu::KTruth gen_l_2 = m_HNgenmatch->gen_l_2;
    snu::KTruth gen_l_3 = m_HNgenmatch->gen_l_3;

    std::vector<snu::KMuon> muontriLooseColl_lowestPtCut_raw = GetMuons("MUON_HN_TRI_VLOOSE_lowestPtCut", true);

    //==== find gen_l_1
    //cout << "[gen_l_1] : pt = " << gen_l_1.Pt() << ", eta = " << gen_l_1.Eta() << endl;
    //cout << "[gen_l_2] : pt = " << gen_l_2.Pt() << ", eta = " << gen_l_2.Eta() << endl;
    //cout << "[gen_l_3] : pt = " << gen_l_3.Pt() << ", eta = " << gen_l_3.Eta() << endl;
    std::vector<int> loose_used;
    loose_used.clear();
    int loose_l_1_index = find_genmatching(gen_l_1, muontriLooseColl_lowestPtCut_raw, loose_used);
    int loose_l_2_index = find_genmatching(gen_l_2, muontriLooseColl_lowestPtCut_raw, loose_used);
    int loose_l_3_index = find_genmatching(gen_l_3, muontriLooseColl_lowestPtCut_raw, loose_used);

    std::vector<snu::KMuon> muontriLooseColl_genorder;
    if(loose_l_1_index!=-1) muontriLooseColl_genorder.push_back( muontriLooseColl_lowestPtCut_raw.at(loose_l_1_index) );
    if(loose_l_2_index!=-1) muontriLooseColl_genorder.push_back( muontriLooseColl_lowestPtCut_raw.at(loose_l_2_index) );
    if(loose_l_3_index!=-1) muontriLooseColl_genorder.push_back( muontriLooseColl_lowestPtCut_raw.at(loose_l_3_index) );

    muontriVLooseColl_lowestPtCut = sort_muons_ptorder( muontriLooseColl_genorder );

  }
  //==== non-prompt : keep fake
  else if( k_sample_name.Contains("DY") || k_sample_name.Contains("WJets") || k_sample_name.Contains("TTJets") || k_sample_name.Contains("QCD") ){
    muontriVLooseColl_lowestPtCut = GetMuons("MUON_HN_TRI_VLOOSE_lowestPtCut", true);
  }
  //==== otherwise
  else{
    muontriVLooseColl_lowestPtCut = GetMuons("MUON_HN_TRI_VLOOSE_lowestPtCut", false);
  }

  //=================================================
  //==== Get Muon Corrections (Rochestor only here)
  //=================================================


  //====================================
  //==== Get Jets (the lowest pt jets)
  //====================================

  std::vector<snu::KJet> jetColl_hn_lowestPtCut = GetJets("JET_HN", 20., 2.4);

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrontriLooseColl = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");

  //===============================
  //==== Get Electron Corrections
  //===============================

  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHT", electrontriLooseColl, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrontriLooseColl);

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight(1.0), pileup_reweight_down(1.0), pileup_reweight_up(1.0);
  if(!k_isdata){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    pileup_reweight_down = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);
    pileup_reweight_up = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),+1);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  //========================
  //==== Apply corrections
  //========================

  if(!k_isdata){
    weight*=trigger_ps_weight;
    //weight*=pileup_reweight;
    weight*=electron_sf;
    weight*=electron_RecoSF;
  }

  //====================================
  //==== Systematic source loop starts
  //==== 1) Muon Energy Scale
  //==== 2) Jet Energy Scale
  //==== 3) Jet Energy Resolution
  //==== 4) Unclustered Energy
  //==== 5) Muon ID Scale Factor
  //==== 6) Pileup
  //====================================

  double m_Z = 91.1876;

  int N_sys = 2*6+1;
  for(int it_sys=0; it_sys<N_sys; it_sys++){

    //==== MET
    //==== also set string for this systematic type
    snu::KEvent Evt = eventbase->GetEvent();
    double MET = Evt.MET();
    double METphi = Evt.METPhi();
    TString this_syst;
    if(it_sys==0){
      this_syst = "MuonEn_up";
      //MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "MuonEn_down";
      //MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::down);
    }
    else if(it_sys==2){
      this_syst = "JetEn_up";
      MET = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::up);
    }
    else if(it_sys==3){
      this_syst = "JetEn_down";
      MET = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::down);
    }
    else if(it_sys==4){
      this_syst = "JetRes_up";
      MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);
    }
    else if(it_sys==5){
      this_syst = "JetRes_down";
      MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);
    }
    else if(it_sys==6){
      this_syst = "Unclustered_up";
      MET = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
    }
    else if(it_sys==7){
      this_syst = "Unclustered_down";
      MET = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
    }
    else if(it_sys==8){
      this_syst = "Central";
    }
    else if(it_sys==9){
      this_syst = "MuonIDSF_up";
    }
    else if(it_sys==10){
      this_syst = "MuonIDSF_down";
    }
    else if(it_sys==11){
      this_syst = "PU_down";
    }
    else if(it_sys==12){
      this_syst = "PU_up";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    //==== Jet
    std::vector<snu::KJet> jetColl_hn;
    if(this_syst == "JetEn_up"){
      for(unsigned int i=0; i<jetColl_hn_lowestPtCut.size(); i++){
        snu::KJet this_jet = jetColl_hn_lowestPtCut.at(i);
        double this_E = this_jet.E()*this_jet.ScaledUpEnergy();
        double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
        double this_3p_sf = this_3p/this_jet.P();
        this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);
        if(this_jet.Pt() >= 30.) jetColl_hn.push_back(this_jet);
      }
    }
    else if(this_syst == "JetEn_down"){
      for(unsigned int i=0; i<jetColl_hn_lowestPtCut.size(); i++){
        snu::KJet this_jet = jetColl_hn_lowestPtCut.at(i);
        double this_E = this_jet.E()*this_jet.ScaledDownEnergy();
        double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
        double this_3p_sf = this_3p/this_jet.P();
        this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);
        if(this_jet.Pt() >= 30.) jetColl_hn.push_back(this_jet);
      }
    }
    else{
      for(unsigned int i=0; i<jetColl_hn_lowestPtCut.size(); i++){
        snu::KJet this_jet = jetColl_hn_lowestPtCut.at(i);
        if(this_jet.Pt() >= 30.) jetColl_hn.push_back(this_jet);
      }
    }
    int n_bjets=0;
    for(int j=0; j<jetColl_hn.size(); j++){
      if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) n_bjets++;
    }

    float BTagSF = BTagScaleFactor_1a_Weighted(jetColl_hn, snu::KJet::CSVv2, snu::KJet::Medium);

    //==== Muon
    std::vector<snu::KMuon> muontriLooseColl;
    if(this_syst == "MuonEn_up"){
      for(unsigned int i=0; i<muontriVLooseColl_lowestPtCut.size(); i++){
        snu::KMuon this_muon = muontriVLooseColl_lowestPtCut.at(i);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedUp();
        if( this_muon.Pt() >= 10. && new_RelIso < this_RelIso ) muontriLooseColl.push_back( this_muon );
      }
    }
    else if(this_syst == "MuonEn_down"){
      for(unsigned int i=0; i<muontriVLooseColl_lowestPtCut.size(); i++){
        snu::KMuon this_muon = muontriVLooseColl_lowestPtCut.at(i);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_RelIso = this_muon.RelIso04()/this_muon.PtShiftedDown();
        if( this_muon.Pt() >= 10. && new_RelIso < this_RelIso ) muontriLooseColl.push_back( this_muon );
      }
    }
    //==== normal muons
    else{
      for(unsigned int i=0; i<muontriVLooseColl_lowestPtCut.size(); i++){
        snu::KMuon this_muon = muontriVLooseColl_lowestPtCut.at(i);
        if( this_muon.Pt() >= 10. && this_muon.RelIso04() < this_RelIso ) muontriLooseColl.push_back( this_muon );
      }
    }

    //==== MET is calculated with No-Rochestor-Corrected Muons
    //==== In this step, muons are 
    //==== 1) Rochestor corrected & Up/Down
    //==== 2) Rochestor corrected
    //==== Both cases, we can correct MET (w.r.t. muon) using
    //==== AnalyzerCore::CorrectedMETRochester(std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi)
    CorrectedMETRochester(muontriLooseColl, MET, METphi);

    //==== Muon ID SF
    double muon_id_iso_sf;
    if(this_syst=="MuonIDSF_up"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor_Weighted("MUON_HN_TRI_TIGHT", muontriLooseColl, 1.); 
    }
    else if(this_syst=="MuonIDSF_down"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor_Weighted("MUON_HN_TRI_TIGHT", muontriLooseColl, -1.);
    }
    else{
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor_Weighted("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
    }

    //==== Muon TrackEff SF
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muontriLooseColl); //FIXME should add syst for this

    //==== this weight
    double this_weight = weight*muon_id_iso_sf*MuTrkEffSF*BTagSF;

    //==== now apply pileup reweight
    if(!k_isdata){
      if(this_syst == "PU_down"){
        this_weight *= pileup_reweight_down;
      }
      else if(this_syst == "PU_up"){
        this_weight *= pileup_reweight_up;
      }
      else{
        this_weight *= pileup_reweight;
      }
    }

    //==== number of leptons
    int n_triLoose_muons = muontriLooseColl.size();
    int n_triTight_muons(0);
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      if(eventbase->GetMuonSel()->MuonPass(muontriLooseColl.at(i), "MUON_HN_TRI_TIGHT")) n_triTight_muons++;
    }

    int n_triLoose_electrons = electrontriLooseColl.size();
    int n_triTight_electrons(0);
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      if(eventbase->GetElectronSel()->ElectronPass(electrontriLooseColl.at(i), "ELECTRON_HN_TIGHT")) n_triTight_electrons++;
    }
    int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
    int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

    bool isThreeMuon   = (n_triLoose_leptons == 3)
                         && (n_triLoose_muons == 3 && n_triTight_muons == 3);
    bool isThreeLepton = (n_triLoose_leptons == 3) && (n_triTight_leptons == 3);
    bool isFourLepton  = (n_triLoose_leptons == 4)
                         && (
                           (n_triLoose_muons == 4 && n_triTight_muons == 4) ||
                           (n_triLoose_electrons == 4 && n_triTight_electrons == 4) ||
                           (n_triLoose_muons == 2 && n_triTight_muons == 2 && n_triLoose_electrons == 2 && n_triTight_electrons == 2)
                         );

    if(n_triLoose_leptons < 3) continue;

    //==============================
    //==== Signal Region Variables
    //==============================

    bool isPreselection(false);
    snu::KParticle HN[4], W_pri_lowmass, W_pri_highmass, W_sec;

    if(isThreeMuon){

      snu::KParticle mu[3];
      for(int i=0;i<3;i++){
        mu[i] = muontriLooseColl.at(i);
      }

      bool AllSameCharge = false;
      int OppSign=0, SameSign[2]={1,2}; // SameSign[0].Pt() > SameSign[1].Pt()
      if(mu[0].Charge() * mu[1].Charge() > 0){ // Q(0) = Q(1)
        if(mu[1].Charge() * mu[2].Charge() < 0){ // Q(1) != Q(2)
          OppSign = 2;
          SameSign[0] = 0;
          SameSign[1] = 1;
        }
        else AllSameCharge = true; // veto Q(0) = Q(1) = Q(2)
      }
      else{ // Q(0) != Q(1)
        if(mu[0].Charge() * mu[2].Charge() > 0){ // Q(0) = Q(2)
          OppSign = 1;
          SameSign[0] = 0;
          SameSign[1] = 2;
        }
        else if(mu[1].Charge() * mu[2].Charge() > 0){ // Q(1) = Q(2)
          OppSign = 0;
          SameSign[0] = 1;
          SameSign[1] = 2;
        }
      } // Find l2 and assign l1&l3 in ptorder 

      bool mllsf4 = false;
      if( (mu[SameSign[0]]+mu[OppSign]).M() <= 4. ||
          (mu[SameSign[1]]+mu[OppSign]).M() <= 4.     ) mllsf4 = true;

      //==== m(HN) < 80 GeV region

      snu::KParticle nu_lowmass;
      nu_lowmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
      double pz_sol_lowmass[2];
      pz_sol_lowmass[0] = solveqdeq(80.4, mu[0]+mu[1]+mu[2], MET, METphi, "m"); // 0 = minus
      pz_sol_lowmass[1] = solveqdeq(80.4, mu[0]+mu[1]+mu[2], MET, METphi, "p"); // 1 = plus
      //PutNuPz(&selection_nu[0], solveqdeq(80.4, mu[0]+mu[1]+mu[2], MET, METphi, "m"));
      //PutNuPz(&selection_nu[1], solveqdeq(80.4, mu[0]+mu[1]+mu[2], MET, METphi, "p")); // 0 = minus, 1 = plus

      int solution_selection_lowmass = 0;
      if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
        //==== take the one with smaller magnitude
        if( fabs( pz_sol_lowmass[0] ) > fabs( pz_sol_lowmass[1] ) ){
          solution_selection_lowmass = 1;
        }
      }
      
      //==== reconstruct HN and W_real 4-vec with selected Pz solution
      PutNuPz(&nu_lowmass, pz_sol_lowmass[solution_selection_lowmass] );
      //==== SameSign[0] : leading among SS
      //==== SameSign[1] : subleading among SS
      //==== [class1]
      //==== HN40, HN50 - SS_leading is primary
      //==== [class2]
      //==== HN60       - SS_subleading is primary

      HN[0] = mu[OppSign] + mu[SameSign[1]] + nu_lowmass; // [class1]
      HN[1] = mu[OppSign] + mu[SameSign[0]] + nu_lowmass; // [class2]
      W_pri_lowmass = mu[0] + mu[1] + mu[2] + nu_lowmass;
      
      //==== m(HN) > 80 GeV region

      snu::KParticle nu_highmass;
      nu_highmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
      int l_3_index = find_mlmet_closest_to_W(mu, nu_highmass);
      double pz_sol_highmass[2]; 
      pz_sol_highmass[0] = solveqdeq(80.4, mu[l_3_index], MET, METphi, "m"); // 0 = minus
      pz_sol_highmass[1] = solveqdeq(80.4, mu[l_3_index], MET, METphi, "p"); // 1 = plus
      int solution_selection_highmass = 0;
      if( pz_sol_highmass[0] != pz_sol_highmass[1] ){ 
        //==== take the one with smaller magnitude
        if( fabs( pz_sol_highmass[0] ) > fabs( pz_sol_highmass[1] ) ){
          solution_selection_highmass = 1;
        }
      }
      PutNuPz( &nu_highmass, pz_sol_highmass[solution_selection_highmass] );

      W_pri_highmass = mu[0] + mu[1] + mu[2] + nu_highmass;

      //==== [class3]
      //==== m(HN) : 90 ~ 1000 GeV - primary lepton has larger pT
      //==== [class4]
      //==== m(HN) > 1000 GeV - primary lepton has smaller pT

      W_sec = mu[l_3_index] + nu_highmass;

      if(l_3_index == OppSign){
        HN[2] = W_sec + mu[SameSign[1]]; // [class3]
        HN[3] = W_sec + mu[SameSign[0]]; // [class4]
      }
      else{
          HN[2] = W_sec + mu[OppSign]; // [class3]
          HN[3] = W_sec + mu[OppSign]; // [class4]
      }

      //==== make mZ for Z veto

      double mz_SR = 0.;
      if( fabs( (mu[OppSign] + mu[SameSign[0]]).M() - 91.1876 ) <
          fabs( (mu[OppSign] + mu[SameSign[1]]).M() - 91.1876 )   ){
        mz_SR = (mu[OppSign] + mu[SameSign[0]]).M();
      }
      else{
        mz_SR = (mu[OppSign] + mu[SameSign[1]]).M();
      }

      bool LeadMuonPtCut = muontriLooseColl.at(0).Pt() > 20.;
      bool bjetveto = (n_bjets==0);
      //==== save nbjet in tree
      bjetveto = true;

      bool VetoZResonance = fabs(mz_SR-91.1876) > 15.;

      //==== Finally, boolean for Preselection
      isPreselection = isThreeMuon && LeadMuonPtCut && (!AllSameCharge) && (!mllsf4) && bjetveto && VetoZResonance;

    }

    //===============================
    //==== Control Region Variables
    //===============================

    bool isWZ(false), isZJets(false), isZLep(false), isZGamma(false);
    int ThreeLeptonConfig = -999;
    bool isZZ(false);
    int FourLeptonConfig = -999;

    if(isThreeLepton){

      std::vector<KLepton> lep;
      for(unsigned int i=0; i<muontriLooseColl.size(); i++){
        KLepton this_lep( muontriLooseColl.at(i) );
        lep.push_back( this_lep );
      }
      for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
        KLepton this_lep( electrontriLooseColl.at(i) );
        lep.push_back( this_lep );
      }

      bool AllSameCharge(false);
      //==== 3 Muon : ThreeLeptonConfig = 0;
      //==== 2 Muon + 1 Electron : ThreeLeptonConfig = 1;
      //==== 1 Muon + 2 Electron ; ThreeLeptonConfig = 2;
      //==== 3 Electron : ThreeLeptonConfig = 3;

      if(muontriLooseColl.size()==3 || electrontriLooseColl.size()==3){

        if(muontriLooseColl.size()==3)     ThreeLeptonConfig = 0;
        if(electrontriLooseColl.size()==3) ThreeLeptonConfig = 3;

        if( ( lep.at(0).Charge() == lep.at(1).Charge() ) &&
            ( lep.at(0).Charge() == lep.at(2).Charge() ) ) AllSameCharge = true;
      }
      else if(muontriLooseColl.size()==2){

        ThreeLeptonConfig = 1;

        if( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) AllSameCharge = true;
      }
      else if(electrontriLooseColl.size()==2){

        ThreeLeptonConfig = 2;

        if( electrontriLooseColl.at(0).Charge() == electrontriLooseColl.at(1).Charge() ) AllSameCharge = true;
      }
      else{
        Message("?", INFO);
      }

      if( !AllSameCharge ){

        snu::KParticle Z_candidate;
        KLepton ZLepton_leading, ZLepton_subleading, WLepton;
        double m_OSSF[2];

        if(ThreeLeptonConfig == 0 || ThreeLeptonConfig == 3){
          KLepton OS, SS[2];
          if     ( lep.at(0).Charge() == lep.at(1).Charge() ){
            SS[0] = lep.at(0);
            SS[1] = lep.at(1);
            OS    = lep.at(2);
          }
          else if( lep.at(0).Charge() == lep.at(2).Charge() ){
            SS[0] = lep.at(0);
            SS[1] = lep.at(2);
            OS    = lep.at(1);
          }
          else if( lep.at(1).Charge() == lep.at(2).Charge() ){
            SS[0] = lep.at(1);
            SS[1] = lep.at(2);
            OS    = lep.at(0);
          }
          else Message("?", INFO);

          m_OSSF[0] = ( SS[0] + OS ).M();
          m_OSSF[1] = ( SS[1] + OS ).M();

          KLepton SS_ZMuon;
          if( fabs(m_OSSF[0]-m_Z) < fabs(m_OSSF[1]-m_Z) ){
            Z_candidate = SS[0] + OS;
            SS_ZMuon = SS[0];
            WLepton = SS[1];
          }
          else{
            Z_candidate = SS[1] + OS;
            SS_ZMuon = SS[1];
            WLepton = SS[0];
          }

          if( SS_ZMuon.Pt() > OS.Pt() ){
            ZLepton_leading = SS_ZMuon;
            ZLepton_subleading = OS;
          }
          else{
            ZLepton_leading = OS;
            ZLepton_subleading = SS_ZMuon;
          }

        }
        else if(ThreeLeptonConfig == 1){
          Z_candidate = muontriLooseColl.at(0)+muontriLooseColl.at(1);
          ZLepton_leading = muontriLooseColl.at(0);
          ZLepton_subleading = muontriLooseColl.at(1);
          WLepton = electrontriLooseColl.at(0);

          m_OSSF[0] = Z_candidate.M();
          m_OSSF[1] = Z_candidate.M();

        }
        else if(ThreeLeptonConfig == 2){
          Z_candidate = electrontriLooseColl.at(0)+electrontriLooseColl.at(1);
          ZLepton_leading = electrontriLooseColl.at(0);
          ZLepton_subleading = electrontriLooseColl.at(1);
          WLepton = muontriLooseColl.at(0);

          m_OSSF[0] = Z_candidate.M();
          m_OSSF[1] = Z_candidate.M();
        }

        double mlll = (lep.at(0)+lep.at(1)+lep.at(2)).M();
        bool ZLeptonPtCut = (ZLepton_leading.Pt() > 20.);
        bool WLeptonPtCut = (WLepton.Pt() > 20.);
        bool isZresonance = (fabs(Z_candidate.M()-m_Z) < 10.);
        bool METCut = (MET > 30.);
        bool mlllCut = (mlll > 100.);
        bool mll4 = (m_OSSF[0] < 4.) || (m_OSSF[1] < 4.);
        bool bjetveto = (n_bjets == 0);
        //==== save nbjet in tree
        bjetveto = true;

        snu::KParticle nu;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
        snu::KParticle W_candidate = nu+WLepton;

        isWZ    = ZLeptonPtCut && isZresonance && WLeptonPtCut && METCut      && mlllCut   && !mll4 && bjetveto;
        isZJets = ZLeptonPtCut && isZresonance                 && (MET < 20.) && mlllCut   && !mll4 && bjetveto && MT(nu, WLepton) < 30.;
        isZLep  = ZLeptonPtCut && isZresonance                                && mlllCut   && !mll4 && bjetveto;
        isZGamma= ZLeptonPtCut && (fabs(Z_candidate.M()-m_Z) > 15.) && (MET < 50.) && (fabs(mlll-m_Z) < 10.) && !mll4 && bjetveto;

      } // Not All Same Charge

    } // isThreeLepton

    if(isFourLepton){

     std::vector<KLepton> lep;
      for(unsigned int i=0; i<muontriLooseColl.size(); i++){
        KLepton this_lep( muontriLooseColl.at(i) );
        lep.push_back( this_lep );
      }
      for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
        KLepton this_lep( electrontriLooseColl.at(i) );
        lep.push_back( this_lep );
      }

      //==== 4 Muon : FourLeptonConfig = 0;
      //==== 2 Muon + 2 Electron : FourLeptonConfig = 1;
      //==== 4 Electron : FourLeptonConfig = 2;

      if(muontriLooseColl.size()==4 || electrontriLooseColl.size()==4){

        if(muontriLooseColl.size()==4)     FourLeptonConfig = 0;
        if(electrontriLooseColl.size()==4) FourLeptonConfig = 2;

      }
      else if(muontriLooseColl.size()==2){

        FourLeptonConfig = 1;

      }
      else{
        Message("?", INFO);
      }

      std::vector<KLepton> LepPlus, LepMinus;
      for(unsigned int i=0; i<lep.size(); i++){
        if(lep.at(i).Charge() > 0) LepPlus.push_back(lep.at(i));
        else                       LepMinus.push_back(lep.at(i));
      }
      if( (LepPlus.size() == 2) && (LepMinus.size() == 2) ){

        snu::KParticle ll_case1_1 = LepPlus.at(0)+LepMinus.at(0);
        snu::KParticle ll_case1_2 = LepPlus.at(1)+LepMinus.at(1);
        bool TwoOnZ_case1 = ( fabs( ll_case1_1.M() - m_Z ) < 10. ) && (LepPlus.at(0).LeptonFlavour()==LepMinus.at(0).LeptonFlavour()) &&
                            ( fabs( ll_case1_2.M() - m_Z ) < 10. ) && (LepPlus.at(1).LeptonFlavour()==LepMinus.at(1).LeptonFlavour());

        snu::KParticle ll_case2_1 = LepPlus.at(0)+LepMinus.at(1);
        snu::KParticle ll_case2_2 = LepPlus.at(1)+LepMinus.at(0);
        bool TwoOnZ_case2 = ( fabs( ll_case2_1.M() - m_Z ) < 10. ) && (LepPlus.at(0).LeptonFlavour()==LepMinus.at(1).LeptonFlavour()) &&
                            ( fabs( ll_case2_2.M() - m_Z ) < 10. ) && (LepPlus.at(1).LeptonFlavour()==LepMinus.at(0).LeptonFlavour());

        if( (TwoOnZ_case1 && min(LepPlus.at(0).Pt(),LepMinus.at(0).Pt())>20. && min(LepPlus.at(1).Pt(),LepMinus.at(1).Pt())>20. ) ||
            (TwoOnZ_case2 && min(LepPlus.at(0).Pt(),LepMinus.at(1).Pt())>20. && min(LepPlus.at(1).Pt(),LepMinus.at(0).Pt())>20. )     ){

          isZZ = true;

        }

      } // 2OS

    } // isFourLepton


    double pt0(0.), pt1(0.), pt2(0.);
    if(isPreselection){
      pt0 = muontriLooseColl.at(0).Pt();
      pt1 = muontriLooseColl.at(1).Pt();
      pt2 = muontriLooseColl.at(2).Pt();
    }

    double cutop[100];
    cutop[0] = pt0;
    cutop[1] = pt1;
    cutop[2] = pt2;
    cutop[3] = HN[0].M();
    cutop[4] = HN[1].M();
    cutop[5] = HN[2].M();
    cutop[6] = HN[3].M();
    cutop[7] = W_pri_lowmass.M();
    cutop[8] = W_pri_highmass.M();
    cutop[9] = this_weight;
    cutop[10] = W_sec.M();
    cutop[11] = MET;
    cutop[12] = 0.; //weight_err
    cutop[13] = isPreselection;
    cutop[14] = isWZ;
    cutop[15] = isZJets;
    cutop[16] = isZLep;
    cutop[17] = isZGamma;
    cutop[18] = isZZ;
    cutop[19] = ThreeLeptonConfig;
    cutop[20] = FourLeptonConfig;
    cutop[21] = n_bjets;

    FillNtp("Ntp_"+this_syst,cutop);

  }

  return;

}// End of execute event loop
  


void trilepton_mumumu_ntp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_ntp::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_ntp::~trilepton_mumumu_ntp() {
  
  Message("In trilepton_mumumu_ntp Destructor" , INFO);
  
}


void trilepton_mumumu_ntp::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 4,0.,4.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    
  }
}


void trilepton_mumumu_ntp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_ntp::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_ntpCore::MakeHistograms() to make new hists for your analysis
   **/

  MakeNtp("Ntp_MuonEn_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_MuonEn_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_JetEn_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_JetEn_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_JetRes_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_JetRes_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_Unclustered_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_Unclustered_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_Central", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_MuonIDSF_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_MuonIDSF_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_PU_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");
  MakeNtp("Ntp_PU_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets");

}


void trilepton_mumumu_ntp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

int trilepton_mumumu_ntp::find_genmatching(snu::KParticle gen, std::vector<snu::KMuon> recos, std::vector<int>& used_index){

  double mindr = 0.1;
  int found=-1;
  for(unsigned int i=0; i<recos.size(); i++){
    //cout << "["<<i<<"th reco] : pt = " << recos.at(i).Pt() << ", eta = " << recos.at(i).Eta() << endl;
    double dr = gen.DeltaR(recos.at(i));
    bool isthisused = std::find(used_index.begin(), used_index.end(), i) != used_index.end();
    if(dr < mindr && !isthisused){
      mindr = dr;
      found = i;
    }
  }
  used_index.push_back(found);
  return found;
}

std::vector<snu::KMuon> trilepton_mumumu_ntp::sort_muons_ptorder(std::vector<snu::KMuon> muons){

  std::vector<snu::KMuon> outmuon;
  while(outmuon.size() != muons.size()){
    double this_maxpt = 0.;
    int index(0);
    for(unsigned int i=0; i<muons.size(); i++){
      bool isthisused = std::find( outmuon.begin(), outmuon.end(), muons.at(i) ) != outmuon.end();
      if(isthisused) continue;
      if( muons.at(i).Pt() > this_maxpt ){
        index = i;
        this_maxpt = muons.at(i).Pt();
      }
    }
    outmuon.push_back( muons.at(index) );
  }
  return outmuon;
 

}


int trilepton_mumumu_ntp::GetSignalMass(){

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


