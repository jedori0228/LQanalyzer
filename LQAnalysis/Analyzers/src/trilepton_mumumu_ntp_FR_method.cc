// $Id: trilepton_mumumu_ntp_FR_method.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_ntp_FR_method Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_ntp_FR_method.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_ntp_FR_method);

trilepton_mumumu_ntp_FR_method::trilepton_mumumu_ntp_FR_method() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_ntp_FR_method");
  
  Message("In trilepton_mumumu_ntp_FR_method constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_ntp_FR_method::InitialiseAnalysis() throw( LQError ) {
 
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


void trilepton_mumumu_ntp_FR_method::ExecuteEvents()throw( LQError ){

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

  std::vector<snu::KMuon>  muontriVLooseColl_lowestPtCut = GetMuons("MUON_HN_TRI_VLOOSE_lowestPtCut", true);

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

  std::vector<snu::KElectron> electrontriVLooseColl_lowestPtCut = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE");

  //====================================
  //==== Systematic source loop starts
  //==== 1) Muon Energy Scale
  //==== 2) Jet Energy Scale
  //==== 3) Jet Energy Resolution
  //==== 4) Unclustered Energy
  //==== 5) Muon ID Scale Factor
  //==== 6) Pileup
  //==== 7) Fake Rate HalfSampleTest Error
  //==== 8) Electron Energy Scale
  //====================================

  double m_Z = 91.1876;

  int N_sys = 2*8+1;
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
      //MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);
    }
    else if(it_sys==5){
      this_syst = "JetRes_down";
      //MET = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);
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
    else if(it_sys==13){
      this_syst = "FR_HalfSample_down";
    }
    else if(it_sys==14){
      this_syst = "FR_HalfSample_up";
    }
    else if(it_sys==15){
      this_syst = "ElectronEn_up";
    }
    else if(it_sys==16){
      this_syst = "ElectronEn_down";
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
    int n_jets = jetColl_hn.size();
    int n_bjets=0;
    for(int j=0; j<jetColl_hn.size(); j++){
      if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    }

    //m_datadriven_bkg->GetFakeObj()->SetNJet(n_jets);
    //m_datadriven_bkg->GetFakeObj()->SetNBJet(n_bjets);

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

    muontriLooseColl = sort_muons_ptorder(muontriLooseColl);

    //==== Electron
    std::vector<snu::KElectron> electrontriLooseColl;
    int ElEnDir = 0;
    if(this_syst == "ElectronEn_up"){
      ElEnDir = 1;
      for(unsigned int i=0; i<electrontriVLooseColl_lowestPtCut.size(); i++){
        snu::KElectron this_electron = electrontriVLooseColl_lowestPtCut.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
        this_electron.SetPFRelIso(0.3, new_RelIso);
        if( this_electron.Pt() >= 15. && PassID(this_electron, "ELECTRON_MVA_FAKELOOSE") ) electrontriLooseColl.push_back( this_electron );
      }
    }
    else if(this_syst == "ElectronEn_down"){
      ElEnDir = -1;
      for(unsigned int i=0; i<electrontriVLooseColl_lowestPtCut.size(); i++){
        snu::KElectron this_electron = electrontriVLooseColl_lowestPtCut.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_RelIso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
        this_electron.SetPFRelIso(0.3, new_RelIso);
        if( this_electron.Pt() >= 15. && PassID(this_electron, "ELECTRON_MVA_FAKELOOSE") ) electrontriLooseColl.push_back( this_electron );
      }
    }
    //==== normal electrons
    else{
      for(unsigned int i=0; i<electrontriVLooseColl_lowestPtCut.size(); i++){
        snu::KElectron this_electron = electrontriVLooseColl_lowestPtCut.at(i);
        if( this_electron.Pt() >= 15. && PassID(this_electron, "ELECTRON_MVA_FAKELOOSE") ) electrontriLooseColl.push_back( this_electron );
      }
    }

    JSCorrectedMETElectron(ElEnDir, electrontriLooseColl, MET, METphi);

    electrontriLooseColl = sort_electrons_ptorder(electrontriLooseColl);

    //==== MET is calculated with No-Rochestor-Corrected Muons
    //==== In this step, muons are 
    //==== 1) Rochestor corrected & Up/Down
    //==== 2) Rochestor corrected
    //==== Both cases, we can correct MET (w.r.t. muon) using
    //==== AnalyzerCore::JSCorrectedMETRochester(std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi)
    JSCorrectedMETRochester(muontriLooseColl, MET, METphi);

    //==== This is a data-driven run,
    //==== So IDSF and TrkEffSF are 1.
    //==== Just to see if they are 1..

    //==== Muon ID SF
    double muon_id_iso_sf;
    if(this_syst=="MuonIDSF_up"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 1.);
    }
    else if(this_syst=="MuonIDSF_down"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, -1.);
    }
    else{
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
    }

    //==== Muon TrackEff SF
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(muontriLooseColl); //FIXME should add syst for this

    //==== number of leptons
    int n_triLoose_muons = muontriLooseColl.size();
    int n_triTight_muons(0);
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      if(PassID(muontriLooseColl.at(i), "MUON_HN_TRI_TIGHT")) n_triTight_muons++;
    }

    int n_triLoose_electrons = electrontriLooseColl.size();
    int n_triTight_electrons(0);
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      if(PassID(electrontriLooseColl.at(i), "ELECTRON_MVA_TIGHT")) n_triTight_electrons++;
    }
    int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
    int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

    bool isThreeMuon   = (n_triLoose_leptons == 3)
                         && (n_triLoose_muons == 3 && n_triTight_muons != 3);
    bool isTwoMuonOneElectron = (n_triLoose_leptons == 3)
                                && (n_triLoose_muons == 2)
                                && (n_triLoose_electrons == 1)
                                && (n_triTight_leptons != 3);
    bool isThreeLepton = (n_triLoose_leptons == 3) && (n_triTight_leptons != 3);
    bool isFourLepton  = (n_triLoose_leptons == 4)
                         && (
                           (n_triLoose_muons == 4 && n_triTight_muons != 4) ||
                           (n_triLoose_electrons == 4 && n_triTight_electrons != 4) ||
                           (n_triLoose_muons == 2 && n_triTight_muons != 2 && n_triLoose_electrons == 2 && n_triTight_electrons != 2)
                         );

    if(n_triLoose_leptons < 3) continue;

    std::vector<KLepton> lep;
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      KLepton this_lep( muontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }
    for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
      KLepton this_lep( electrontriLooseColl.at(i) );
      lep.push_back( this_lep );
    }

    int ThreeLeptonConfig = -999, FourLeptonConfig = -999;

    if(isThreeLepton){

      //==== 3 Muon : ThreeLeptonConfig = 0;
      //==== 2 Muon + 1 Electron : ThreeLeptonConfig = 1;
      //==== 1 Muon + 2 Electron ; ThreeLeptonConfig = 2;
      //==== 3 Electron : ThreeLeptonConfig = 3;

      if(muontriLooseColl.size()==3 || electrontriLooseColl.size()==3){

        if(muontriLooseColl.size()==3)     ThreeLeptonConfig = 0;
        if(electrontriLooseColl.size()==3) ThreeLeptonConfig = 3;

      }
      else if(muontriLooseColl.size()==2){

        ThreeLeptonConfig = 1;

      }
      else if(electrontriLooseColl.size()==2){

        ThreeLeptonConfig = 2;

      }
      else{

        Message("?", INFO);

      }
    }

    if(isFourLepton){

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
    }

    //==============================
    //==== Signal Region Variables
    //==============================

    bool isPreselection(false);
    snu::KParticle HN[4], W_pri_lowmass, W_pri_highmass, W_sec, gamma_star;
    double deltaR_OS_min(-999);

    if(isThreeMuon || isTwoMuonOneElectron){

     bool AllSameCharge = false;
      int OppSign=0, SameSign[2]={1,2}; // SameSign[0].Pt() > SameSign[1].Pt()
      if(isThreeMuon){

        if(lep[0].Charge() * lep[1].Charge() > 0){ // Q(0) = Q(1)
          if(lep[1].Charge() * lep[2].Charge() < 0){ // Q(1) != Q(2)
            OppSign = 2;
            SameSign[0] = 0;
            SameSign[1] = 1;
          }
          else AllSameCharge = true; // veto Q(0) = Q(1) = Q(2)
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

      }
      else if(isTwoMuonOneElectron){

        //==== We want SS muon.
        //==== If OS, return
        if(muontriLooseColl.at(0).Charge() != muontriLooseColl.at(1).Charge()) AllSameCharge = true;
        //==== Remaining electron should have opposite sign
        //==== If SS, return
        if(muontriLooseColl.at(0).Charge() == electrontriLooseColl.at(0).Charge()) AllSameCharge = true;

        OppSign = 2;
        SameSign[0] = 0;
        SameSign[1] = 1;

      }
      else{}

      bool mllsf4 = false;
      if(isThreeMuon){
        if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
            (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) mllsf4 = true;
      }

      //==== m(HN) < 80 GeV region

      snu::KParticle nu_lowmass;
      nu_lowmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
      double pz_sol_lowmass[2];
      pz_sol_lowmass[0] = solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "m"); // 0 = minus
      pz_sol_lowmass[1] = solveqdeq(80.385, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 1 = plus

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

      HN[0] = lep[OppSign] + lep[SameSign[1]] + nu_lowmass; // [class1]
      HN[1] = lep[OppSign] + lep[SameSign[0]] + nu_lowmass; // [class2]
      W_pri_lowmass = lep[0] + lep[1] + lep[2] + nu_lowmass;

      //==== m(HN) > 80 GeV region

      snu::KParticle nu_highmass;
      nu_highmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
      int l_3_index(0);
      if(isThreeMuon) l_3_index = find_mlmet_closest_to_W(lep, nu_highmass);
      else if(isTwoMuonOneElectron) l_3_index = 2;
      double pz_sol_highmass[2];
      pz_sol_highmass[0] = solveqdeq(80.385, lep[l_3_index], MET, METphi, "m"); // 0 = minus
      pz_sol_highmass[1] = solveqdeq(80.385, lep[l_3_index], MET, METphi, "p"); // 1 = plus
      int solution_selection_highmass = 0;
      if( pz_sol_highmass[0] != pz_sol_highmass[1] ){
        //==== take the one with smaller magnitude
        if( fabs( pz_sol_highmass[0] ) > fabs( pz_sol_highmass[1] ) ){
          solution_selection_highmass = 1;
        }
      }
      PutNuPz( &nu_highmass, pz_sol_highmass[solution_selection_highmass] );

      W_pri_highmass = lep[0] + lep[1] + lep[2] + nu_highmass;

      //==== [class3]
      //==== m(HN) : 90 ~ 1000 GeV - primary lepton has larger pT
      //==== [class4]
      //==== m(HN) > 1000 GeV - primary lepton has smaller pT

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

      //==== make mZ for Z veto

      if( lep[OppSign].DeltaR(lep[SameSign[0]]) < lep[OppSign].DeltaR(lep[SameSign[1]]) ){
        deltaR_OS_min = lep[OppSign].DeltaR(lep[SameSign[0]]);
        gamma_star = lep[OppSign] + lep[SameSign[0]];
      }
      else{
        deltaR_OS_min = lep[OppSign].DeltaR(lep[SameSign[1]]);
        gamma_star = lep[OppSign] + lep[SameSign[1]];
      }

      double mz_SR = 0.;
      if( fabs( (lep[OppSign] + lep[SameSign[0]]).M() - 91.1876 ) <
          fabs( (lep[OppSign] + lep[SameSign[1]]).M() - 91.1876 )   ){
        mz_SR = (lep[OppSign] + lep[SameSign[0]]).M();
      }
      else{
        mz_SR = (lep[OppSign] + lep[SameSign[1]]).M();
      }

      bool LeadMuonPtCut = muontriLooseColl.at(0).Pt() > 20.;
      bool bjetveto = (n_bjets==0);
      //==== save nbjet in tree
      bjetveto = true;

      bool VetoZResonance = fabs(mz_SR-91.1876) > 15.;
      if(isTwoMuonOneElectron) VetoZResonance = true;

      bool mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;
      if(isTwoMuonOneElectron) mllloffZ = true;

      //==== Finally, boolean for Preselection
      isPreselection = LeadMuonPtCut && (!AllSameCharge) && (!mllsf4) && bjetveto && VetoZResonance && mllloffZ;

    }

    //===============================
    //==== Control Region Variables
    //===============================

    bool isWZ(false), isZJets(false), isZLep(false), isZGamma(false);
    bool isZZ(false);

    if(isThreeLepton){

      bool AllSameCharge(false);

      //==== 3 Muon : ThreeLeptonConfig = 0;
      //==== 2 Muon + 1 Electron : ThreeLeptonConfig = 1;
      //==== 1 Muon + 2 Electron ; ThreeLeptonConfig = 2;
      //==== 3 Electron : ThreeLeptonConfig = 3;

      if(ThreeLeptonConfig==0 || ThreeLeptonConfig==3){

        if( ( lep.at(0).Charge() == lep.at(1).Charge() ) &&
            ( lep.at(0).Charge() == lep.at(2).Charge() ) ) AllSameCharge = true;

      }
      else if(ThreeLeptonConfig==1){

        if( muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge() ) AllSameCharge = true;

      }
      else if(ThreeLeptonConfig==2){

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

        if( (TwoOnZ_case1 && max(LepPlus.at(0).Pt(),LepMinus.at(0).Pt())>20. && max(LepPlus.at(1).Pt(),LepMinus.at(1).Pt())>20. ) ||
            (TwoOnZ_case2 && max(LepPlus.at(0).Pt(),LepMinus.at(1).Pt())>20. && max(LepPlus.at(1).Pt(),LepMinus.at(0).Pt())>20. )     ){

          isZZ = true;

        }

      } // 2OS

    } // isFourLepton

    lep = sort_leptons_ptorder(lep);

    double pt0(0.), pt1(0.), pt2(0.);
    if(isPreselection){

      pt0 = lep.at(0).Pt();
      pt1 = lep.at(1).Pt();
      pt2 = lep.at(2).Pt();

    }

    //==== fake method weighting
    int Dir_HalfSample = 0;
    if(this_syst=="FR_HalfSample_down") Dir_HalfSample = -1;
    if(this_syst=="FR_HalfSample_up") Dir_HalfSample = +1;

    double this_weight =     m_datadriven_bkg->Get_DataDrivenWeight(false, muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40", Dir_HalfSample);
    double this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40", Dir_HalfSample);

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
    cutop[22] = deltaR_OS_min;
    cutop[23] = gamma_star.M();

    FillNtp("Ntp_"+this_syst,cutop);

  }

  return;

}// End of execute event loop
  


void trilepton_mumumu_ntp_FR_method::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_ntp_FR_method::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_ntp_FR_method::~trilepton_mumumu_ntp_FR_method() {
  
  Message("In trilepton_mumumu_ntp_FR_method Destructor" , INFO);
  
}


void trilepton_mumumu_ntp_FR_method::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_ntp_FR_method::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_ntp_FR_method::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_ntp_FR_methodCore::MakeHistograms() to make new hists for your analysis
   **/

MakeNtp("Ntp_MuonEn_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_MuonEn_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_JetEn_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_JetEn_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_JetRes_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_JetRes_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_Unclustered_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_Unclustered_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_Central", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_MuonIDSF_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_MuonIDSF_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_PU_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_PU_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_FR_HalfSample_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_FR_HalfSample_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_ElectronEn_up", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");
MakeNtp("Ntp_ElectronEn_down", "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets:isZLep:isZGamma:isZZ:ThreeLeptonConfig:FourLeptonConfig:nbjets:deltaR_OS_min:gamma_star_mass");


}


void trilepton_mumumu_ntp_FR_method::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}







