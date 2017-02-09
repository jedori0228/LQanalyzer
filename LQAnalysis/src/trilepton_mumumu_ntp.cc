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
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_ntp);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu_ntp::trilepton_mumumu_ntp() :  AnalyzerCore(),
out_muons(0)
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

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", 1.);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerlist;
  //triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  double MinLeadingMuonPt = 20;
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
  FillCutFlow("VertexCut", 1.);


  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  //std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT);  /// No cuts applied
  //std::vector<snu::KMuon> muonVetoColl = GetMuons(BaseSelection::MUON_HN_VETO);  // veto selection
  //std::vector<snu::KMuon> muonLooseColl = GetMuons(BaseSelection::MUON_HN_FAKELOOSE);  // loose selection
  //std::vector<snu::KMuon> muonTightColl = GetMuons(BaseSelection::MUON_HN_TIGHT,false); // tight selection : NonPrompt MC lep removed

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

  //CorrectMuonMomentum(muontriVLooseColl_lowestPtCut); //FIXME do this for v8-0-4

  std::vector<snu::KJet> jetColl_hn_lowestPtCut = GetJets("JET_HN", true, 20., 2.4);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  numberVertices = eventbase->GetEvent().nVertices();

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    //pileup_reweight = eventbase->GetEvent().PileUpWeight();
    //pileup_reweight = eventbase->GetEvent().AltPileUpWeight();
    pileup_reweight = TempPileupWeight();
  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);


  if(!isData && !k_running_nonprompt){
    weight*=trigger_ps_weight;
    weight*=pileup_reweight;
  }

  int N_sys = 2*5+1;
  for(int it_sys=0; it_sys<N_sys; it_sys++){

    //==== MET
    //==== also set string for this systematic type
    snu::KEvent Evt = eventbase->GetEvent();
    double MET;
    TString this_syst;
    if(it_sys==0){
      this_syst = "MuonEn_up";
      MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "MuonEn_down";
      MET = Evt.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::down);
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
      MET = Evt.MET();
    }
    else if(it_sys==9){
      this_syst = "MuonIDSF_up";
      MET = Evt.MET();
    }
    else if(it_sys==10){
      this_syst = "MuonIDSF_down";
      MET = Evt.MET();
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
      if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)) n_bjets++;
    }

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

    //==== Muon ID SF
    double muon_id_iso_sf;
    if(this_syst=="MuonIDSF_up"){
      muon_id_iso_sf = MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 1.); 
    }
    else if(this_syst=="MuonIDSF_down"){
      muon_id_iso_sf = MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, -1.);
    }
    else{
      muon_id_iso_sf = MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
    }

    //==== Muon TrackEff SF
    double MuTrkEffSF =  MuonTrackingEffScaleFactor(muontriLooseColl); //FIXME should add syst for this

    //==== this weight
    double this_weight = weight*muon_id_iso_sf*MuTrkEffSF;

    double METphi = Evt.METPhi();

    int n_triTight_muons(0);
    for(unsigned int i=0; i<muontriLooseColl.size(); i++){
      if(muontriLooseColl.at(i).RelIso04() < 0.1) n_triTight_muons++;
    }
    int n_triLoose_muons = muontriLooseColl.size();
    int n_jets = jetColl_hn.size();

    if( n_triLoose_muons != 3 ) continue;
    if( n_triTight_muons != 3 ) continue;

    if( muontriLooseColl.at(0).Pt() < MinLeadingMuonPt ) continue;
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
      else continue; // veto Q(0) = Q(1) = Q(2)
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

    // MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
    // MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
    // POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
    if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
        (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) continue;
    FillCutFlow("mllsf4", 1.);

    ///////////////////////////////////////////
    ////////// m(HN) < 80 GeV region //////////
    ///////////////////////////////////////////

    snu::KParticle W_pri_lowmass, nu_lowmass, gamma_star, z_candidate;
    nu_lowmass.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);
    double pz_sol_lowmass[2];
    pz_sol_lowmass[0] = solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"); // 0 = minus
    pz_sol_lowmass[1] = solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 1 = plus
    //PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
    //PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus

    int solution_selection_lowmass = 0;
    if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
      // take the one with smaller magnitude
      if( fabs( pz_sol_lowmass[0] ) > fabs( pz_sol_lowmass[1] ) ){
        solution_selection_lowmass = 1;
      }
    }
    
    // reconstruct HN and W_real 4-vec with selected Pz solution
    PutNuPz(&nu_lowmass, pz_sol_lowmass[solution_selection_lowmass] );
    // SameSign[0] : leading among SS
    // SameSign[1] : subleading among SS
    // [class1]
    // HN40, HN50 - SS_leading is primary
    // [class2]
    // HN60       - SS_subleading is primary

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
    pz_sol_highmass[0] = solveqdeq(80.4, lep[l_3_index], MET, METphi, "m"); // 0 = minus
    pz_sol_highmass[1] = solveqdeq(80.4, lep[l_3_index], MET, METphi, "p"); // 1 = plus
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
    // m(HN) : 90 ~ 500 GeV - primary lepton has larger pT
    // [class4]
    // m(HN) : 700 ~ 1000 GeV - primary lepton has smaller pT

    W_sec = lep[l_3_index] + nu_highmass;

    if(l_3_index == OppSign){
       
      HN[2] = W_sec + lep[SameSign[1]]; // [class3]
      HN[3] = W_sec + lep[SameSign[0]]; // [class4]
   
    }
    else{
        HN[2] = W_sec + lep[OppSign]; // [class3]
        HN[3] = W_sec + lep[OppSign]; // [class4]
    }

    if(n_bjets>0) continue;

    bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
    bool UseZResonance = fabs(z_candidate.M()-91.1876) < 10.;

    double m_dimuon[2], m_Z = 91.1876;
    m_dimuon[0] = ( lep[SameSign[0]] + lep[OppSign] ).M();
    m_dimuon[1] = ( lep[SameSign[1]] + lep[OppSign] ).M();
    
    snu::KParticle Z_candidate;
    snu::KParticle ZMuon, WMuon;
    if( fabs(m_dimuon[0]-m_Z) < fabs(m_dimuon[1]-m_Z) ){
      Z_candidate = lep[SameSign[0]] + lep[OppSign];
      ZMuon = lep[SameSign[0]];
      WMuon = lep[SameSign[1]];
    }
    else{
      Z_candidate = lep[SameSign[1]] + lep[OppSign];
      ZMuon = lep[SameSign[1]];
      WMuon = lep[SameSign[0]];
    }

    snu::KParticle ZMuon_leading, ZMuon_subleading;
    if( ZMuon.Pt() > lep[OppSign].Pt() ){
      ZMuon_leading = ZMuon;
      ZMuon_subleading = lep[OppSign];
    }
    else{
      ZMuon_leading = lep[OppSign];
      ZMuon_subleading = ZMuon;
    } 

    bool ZMuonPtCut = (ZMuon.Pt() > 20.) || (lep[OppSign].Pt() > 20.);
    bool PtCutOnWMuon = (WMuon.Pt() > 20.);
    bool METCut = (MET > 30.);
    bool mlllCut = ((lep[SameSign[0]]+lep[SameSign[1]]+lep[OppSign]).M() > 100.);
    bool electronveto = (electronLooseColl.size() == 0);

    bool isPreselection = VetoZResonance;
    bool isWZ = ZMuonPtCut && UseZResonance && PtCutOnWMuon && METCut && mlllCut && electronveto;
    bool isZJets = ZMuonPtCut && UseZResonance && (MET < 20.) && mlllCut && electronveto;

    double cutop[100];
    cutop[0] = muontriLooseColl.at(0).Pt();
    cutop[1] = muontriLooseColl.at(1).Pt();
    cutop[2] = muontriLooseColl.at(2).Pt();
    cutop[3] = deltaR_OS_min;
    cutop[4] = HN[0].M();
    cutop[5] = HN[1].M();
    cutop[6] = HN[2].M();
    cutop[7] = HN[3].M();
    cutop[8] = W_pri_lowmass.M();
    cutop[9] = W_pri_highmass.M();
    cutop[10] = this_weight;
    cutop[11] = W_sec.M();
    cutop[12] = MET;
    cutop[13] = 0.; //weight_err
    cutop[14] = isPreselection;
    cutop[15] = isWZ;
    cutop[16] = isZJets;
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
    AnalyzerCore::MakeHistograms("cutflow", 7,0.,7.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3muon");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS"); 
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"mllsf4");
    
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

  MakeNtp("Ntp_MuonEn_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_MuonEn_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_JetEn_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_JetEn_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_JetRes_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_JetRes_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_Unclustered_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_Unclustered_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_Central", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_MuonIDSF_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");
  MakeNtp("Ntp_MuonIDSF_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:W_sec_highmass_mass:PFMET:weight_err:isPreselection:isWZ:isZJets");

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


