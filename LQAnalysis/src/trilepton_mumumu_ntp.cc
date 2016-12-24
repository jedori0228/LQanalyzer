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
n_gen_pass(0), sol_sel_chi2_best(0), sol_sel_chi2_plus(0), sol_sel_chi2_minus(0), sol_sel_chi2_smaller(0), sol_sel_chi2_larger(0),
allgenfound(false),
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
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  //triggerlist.push_back("HLT_TripleMu_12_10_5_v");
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
    //==== save gen particles @ snu::KParticle gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;
    find_genparticles();

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

  //CorrectMuonMomentum(muonTightColl);
  //float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muontriTightColl, 0); ///MUON_POG_TIGHT == MUON_HN_TIGHT
  //muon_id_iso_sf *= MuonISOScaleFactor(BaseSelection::MUON_POG_TIGHT, muontriTightColl, 0);

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN");// pt > 20 ; eta < 2.5; PFlep veto; pileup ID

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
    //pileup_reweight = eventbase->GetEvent().PileUpWeight();
    //pileup_reweight = eventbase->GetEvent().AltPileUpWeight();
    pileup_reweight = TempPileupWeight();
  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);


  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    //weight*=weight_trigger_sf;
    weight*=trigger_ps_weight;
    weight*=pileup_reweight;
  }

  int N_sys = 2*4+1;
  for(int it_sys=0; it_sys<=N_sys; it_sys++){

    //==== MET
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
    else{
      Message("it_sys out of range!" , INFO);
      return;
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

    if(k_sample_name.Contains("HN") && allgenfound) solution_selection_stduy(muontriLooseColl);

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
    cutop[10] = weight;
    FillNtp("Ntp_"+this_syst,cutop);

  }


  return;

}// End of execute event loop
  


void trilepton_mumumu_ntp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  cout
  << "n_gen_pass = " << n_gen_pass << endl
  << "best = " << sol_sel_chi2_best/n_gen_pass << endl
  << "plus = " << sol_sel_chi2_plus/n_gen_pass << endl
  << "minus = " << sol_sel_chi2_minus/n_gen_pass << endl
  << "smaller = " << sol_sel_chi2_smaller/n_gen_pass << endl
  << "larger = " << sol_sel_chi2_larger/n_gen_pass << endl;

  TH1F* GEN_solution_selection_chi2 = new TH1F("GEN_solution_selection_chi2", "", 6, 0, 6);
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(1, "n_gen_pass");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(2, "best");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(3, "plus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(4, "minus");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(5, "smaller");
  GEN_solution_selection_chi2->GetXaxis()->SetBinLabel(6, "larger");
  GEN_solution_selection_chi2->SetBinContent(1, n_gen_pass);
  GEN_solution_selection_chi2->SetBinContent(2, sol_sel_chi2_best);
  GEN_solution_selection_chi2->SetBinContent(3, sol_sel_chi2_plus);
  GEN_solution_selection_chi2->SetBinContent(4, sol_sel_chi2_minus);
  GEN_solution_selection_chi2->SetBinContent(5, sol_sel_chi2_smaller);
  GEN_solution_selection_chi2->SetBinContent(6, sol_sel_chi2_larger);

  m_outputFile->cd();
  GEN_solution_selection_chi2->Write();

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

  MakeNtp("Ntp_MuonEn_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_MuonEn_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_JetEn_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_JetEn_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_JetRes_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_JetRes_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_Unclustered_up", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_Unclustered_down", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");
  MakeNtp("Ntp_Central", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight");

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

void trilepton_mumumu_ntp::find_genparticles(){

  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  //==== print truth info
  //cout << "=========================================================" << endl;
  //cout << "truth size = " << truthColl.size() << endl;
  //cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
  //for(int i=2; i<truthColl.size(); i++){
  //  cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << endl;
  //}

  //==== W_pri : first W
  //==== l_1 : lepton from first W
  //==== l_2 : lepton from HN
  //==== W_sec : W from HN
  //==== l_3 : lepton from second W

  int truthmax = truthColl.size();
  vector<int> gen_HN_indices, gen_W_pri_indices, gen_W_sec_indices, gen_l_1_indices, gen_nu_indices, gen_l_3_indices, gen_l_2_indices;
  bool W_sec_in_truth=false, isLowMass = false;

  //==== check if this is low/high mass region
  if(k_sample_name.Contains("HN40_") || k_sample_name.Contains("HN50_") || k_sample_name.Contains("HN60_")) isLowMass = true;
 
  //==== find HN index
  for(int i=2;i<truthmax;i++){
    if(truthColl.at(i).PdgId() == 80000002){
      gen_HN_indices.push_back(i);
      find_decay(truthColl, i, gen_HN_indices);
      break;
    }
  }
  //print_all_indices("gen_HN", gen_HN_indices);
  if(gen_HN_indices.size() == 0) return;
  int HN_motherindex = truthColl.at( gen_HN_indices.at(0) ).IndexMother();
  FillHist("TEST_HN_mother_pdgid", abs(truthColl.at(HN_motherindex).PdgId()), 1., 0., 30., 30);

  //======================
  //==== low mass region
  //======================
  if(isLowMass){
    //==== for low mass, W_pri is on-shell
    //==== so we can find them in gen particle collections
    //==== find W_pri index
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 24){
        gen_W_pri_indices.push_back(i);
        find_decay(truthColl, i, gen_W_pri_indices);
        break;
      }
    }
    //print_all_indices("gen_W_pri", gen_W_pri_indices);
    if(gen_W_pri_indices.size() == 0) return;

    //==== find l_1 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_pri_indices.size();j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = W_pri
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_pri_indices.at(j) ){
          gen_l_1_indices.push_back(i);
          find_decay(truthColl, i, gen_l_1_indices);
          break;
        }
      }
      if(gen_l_1_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_1", gen_l_1_indices);
    if(gen_l_1_indices.size() == 0) return;

    //==== As m(HN) goes closer to W mass, W_sec starting to appear in truth coll
    //==== if W_sec is virtual, it may not appear in the gen particle collections
    //==== check W_sec exists
    for(int i=2;i<truthmax;i++){
      //==== 1) |PID| = 24
      //==== 2) Mother's PID = 80000002
      if(fabs(truthColl.at(i).PdgId()) == 24 && truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 80000002){
        gen_W_sec_indices.push_back(i);
        find_decay(truthColl, i, gen_W_sec_indices);
        W_sec_in_truth = true;
        break;
      }
    }
    //print_all_indices("gen_W_sec", gen_W_sec_indices);

    //==== no W_sec in truthColl index case
    if(!W_sec_in_truth){
      //cout << "[No W_sec!]" << endl;

      //==== find nu at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = 14
          //==== 2) Mother = {HN}
          if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_nu_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_nu_indices.size() != 0) break;
      }
      //print_all_indices("gen_nu", gen_nu_indices);
      if(gen_nu_indices.size() == 0) return;

      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = PID of l_1 (SS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_l_3_indices);
            break;
          }
        }
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0) return;

      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = - (PID of l_1) (OS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
        if(gen_l_2_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0) return;
    }
    //==== W_sec in truthColl index case
    else{
      //cout << "[W_sec!]" << endl;
      //==== find nu at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
          //==== 1) |PID| = 14
          //==== 2) Mother = {W_sec}
          if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
            gen_nu_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_nu_indices.size() != 0) break;
      }
      //print_all_indices("gen_nu", gen_nu_indices);
      if(gen_nu_indices.size() == 0) return;

      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
          //==== 1) |PID| = 13
          //==== 2) Mother = {W_sec}
          if(fabs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_nu_indices);
            break;
          }
        }
        if(gen_l_3_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0) return;

      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) |PID| = 13
          //==== 2) Mother = {HN}
          if(fabs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
        if(gen_l_2_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0) return;
      
    }

    gen_l_1 = truthColl.at( gen_l_1_indices.back() );
    gen_l_2 = truthColl.at( gen_l_2_indices.back() );
    gen_l_3 = truthColl.at( gen_l_3_indices.back() );
    gen_nu = truthColl.at( gen_nu_indices.back() );
    allgenfound = true;

  }

  //=======================
  //==== high mass region
  //=======================
  else{
    //==== find l_1 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<1;j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = Mother of HN.at(0) : becase they are generated at the same time.
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == truthColl.at(gen_HN_indices.at(j)).IndexMother()){
          gen_l_1_indices.push_back(i);
          find_decay(truthColl, i, gen_l_1_indices);
          break;
        }
      }
      if(gen_l_1_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_1", gen_l_1_indices);
    if(gen_l_1_indices.size() == 0) return;

    //==== fine l_2 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_HN_indices.size();j++){
      //==== 1) |PID| = 13
      //==== 2) Mother = {HN}
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
          gen_l_2_indices.push_back(i);
          find_decay(truthColl, i, gen_l_2_indices);
          break;
        }
      }
      if(gen_l_2_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_2", gen_l_2_indices);    
    if(gen_l_2_indices.size() == 0) return;

    //==== find W_sec
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_HN_indices.size();j++){
        //==== 1) |PID| = 24
        //==== 2) Mother = {HN}
        if(abs(truthColl.at(i).PdgId()) == 24 && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
          gen_W_sec_indices.push_back(i);
          find_decay(truthColl, i, gen_W_sec_indices);
          break;
        }
      }
      if(gen_W_sec_indices.size() != 0) break;
    }
    //print_all_indices("gen_W_sec", gen_W_sec_indices);
    if(gen_W_sec_indices.size() == 0) return;

    //==== find nu at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
        //==== 1) |PID| = 14
        //==== 2) Mother = {W_sec}
        if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
          gen_nu_indices.push_back(i);
          find_decay(truthColl, i, gen_nu_indices);
          break;
        }
      }
      if(gen_nu_indices.size() != 0) break;
    }
    //print_all_indices("gen_nu", gen_nu_indices);
    if(gen_nu_indices.size() == 0) return;

    //==== find l_3 at gen. level
    for(int i=2;i<truthmax;i++){
      for(unsigned int j=0;j<gen_W_sec_indices.size();j++){
        //==== 1) |PID| = 13
        //==== 2) Mother = {W_sec}
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_sec_indices.at(j) ){
          gen_l_3_indices.push_back(i);
          find_decay(truthColl, i, gen_l_3_indices);
          break;
        }
      }
      if(gen_l_3_indices.size() != 0) break;
    }
    //print_all_indices("gen_l_3", gen_l_3_indices);
    if(gen_l_3_indices.size() == 0) return;

    gen_l_1 = truthColl.at( gen_l_1_indices.back() );
    gen_l_2 = truthColl.at( gen_l_2_indices.back() );
    gen_l_3 = truthColl.at( gen_l_3_indices.back() );
    gen_nu = truthColl.at( gen_nu_indices.back() );
    allgenfound = true;

  }

}

void trilepton_mumumu_ntp::solution_selection_stduy(std::vector<snu::KMuon> recomuons){

  snu::KParticle reco_lep[3];
  for(int i=0;i<3;i++){
    reco_lep[i] = recomuons.at(i);
  }

  int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
  if(reco_lep[0].Charge() * reco_lep[1].Charge() > 0){ // Q(0) = Q(1)
    if(reco_lep[1].Charge() * reco_lep[2].Charge() < 0){ // Q(1) != Q(2)
      OppSign = 2;
      SameSign[0] = 0;
      SameSign[1] = 1;
    }
    else return; // veto Q(0) = Q(1) = Q(2)
  }
  else{ // Q(0) != Q(1)
    if(reco_lep[0].Charge() * reco_lep[2].Charge() > 0){ // Q(0) = Q(2)
      OppSign = 1;
      SameSign[0] = 0;
      SameSign[1] = 2;
    }
    else if(reco_lep[1].Charge() * reco_lep[2].Charge() > 0){ // Q(1) = Q(2)
      OppSign = 0;
      SameSign[0] = 1;
      SameSign[1] = 2;
    }
  } // Find l2 and assign l1&l3 in ptorder

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
  snu::KParticle reco_MET;
  reco_MET.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, MET);

  //==== W_pri : first W
  //==== l_1 : lepton from first W
  //==== l_2 : lepton from HN
  //==== W_sec : W from HN
  //==== l_3 : lepton from second W
  
  bool isLowMass = false;

  //==== check if this is low/high mass region
  if(k_sample_name.Contains("HN40_") || k_sample_name.Contains("HN50_") || k_sample_name.Contains("HN60_")) isLowMass = true;
 
  //======================
  //==== low mass region
  //======================
  if(isLowMass){

    //==== solution selection
    double pz_sol_lowmass[2];
    pz_sol_lowmass[0] = solveqdeq(80.4, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "m"); // 0 = minus
    pz_sol_lowmass[1] = solveqdeq(80.4, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "p"); // 1 = plus
    if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
      FillHist("GEN_all_found", 0, 1., 0., 1., 1);
      n_gen_pass++;
      int best_sel = fabs(pz_sol_lowmass[0]-gen_nu.Pz()) < fabs(pz_sol_lowmass[1]-gen_nu.Pz()) ? 0 : 1;
      int smaller = fabs(pz_sol_lowmass[0]) < fabs(pz_sol_lowmass[1]) ? 0 : 1;
      int larger = smaller == 0 ? 1 : 0;
      //cout
      //<< "gen_nu.Pz() = " << gen_nu.Pz() << endl
      //<< "pz_sol_lowmass[0] = " << pz_sol_lowmass[0] << endl
      //<< "pz_sol_lowmass[1] = " << pz_sol_lowmass[1] << endl;
      FillHist("GEN_chi2_best", pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_chi2_plus", pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_chi2_minus", pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_chi2_smaller", pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_chi2_larger", pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);

      if(best_sel == 0) FillHist("GEN_solsel_minus_0_plus_1", 0, 1., 0., 2., 2);
      else              FillHist("GEN_solsel_minus_0_plus_1", 1, 1., 0., 2., 2);
      if(best_sel == smaller) FillHist("GEN_solsel_smaller_0_larger_1", 0, 1., 0., 2., 2);
      else                    FillHist("GEN_solsel_smaller_0_larger_1", 1, 1., 0., 2., 2);

      sol_sel_chi2_best += pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2);
      sol_sel_chi2_plus += pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2);
      sol_sel_chi2_minus += pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2);
      sol_sel_chi2_smaller += pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2);
      sol_sel_chi2_larger += pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2);
    }


  }

  //=======================
  //==== high mass region
  //=======================
  else{

    snu::KParticle reco_lep_tlv[3];
    for(int i=0; i<3; i++) reco_lep_tlv[i] = reco_lep[i];
    int l_3_cand = find_mlmet_closest_to_W(reco_lep_tlv, reco_MET);

    FillHist("GEN_reco_MET", reco_MET.Pt(), 1., 0., 120., 120);
    FillHist("GEN_reco_lep_1_MET", (reco_lep[0] + reco_MET).M() - 80.4, 1., -60., 60., 120);
    FillHist("GEN_reco_lep_2_MET", (reco_lep[1] + reco_MET).M() - 80.4, 1., -60., 60., 120);
    FillHist("GEN_reco_lep_3_MET", (reco_lep[2] + reco_MET).M() - 80.4, 1., -60., 60., 120);

    FillHist("GEN_l_3_cand", l_3_cand, 1., 0., 3., 3);
    //if( gen_l_3.DeltaR(reco_lep[l_3_cand]) < 0.1 ) FillHist("GEN_highmass_mlmet_Wmass_match_gen_l_3", 1, 1., 0., 2., 2);
    if( GenMatching(gen_l_3, reco_lep[l_3_cand], 0.1, 0.05)  ) FillHist("GEN_highmass_mlmet_Wmass_match_gen_l_3", 1, 1., 0., 2., 2);
    else FillHist("GEN_highmass_mlmet_Wmass_match_gen_l_3", 0, 1., 0., 2., 2);

    //==== 1) pt ordering firstly done = m1
    
    int l_1_cand_m1 = SameSign[0], l_SS_rem = SameSign[1], signal_class = 3;
    //==== pt ordering reversed for m(HN) >= 700 GeV
    if( k_sample_name.Contains("HN700_") || k_sample_name.Contains("HN1000_") ){
      signal_class = 4;
      l_1_cand_m1 = SameSign[1];
      l_SS_rem = SameSign[0];
    }
    //if( reco_lep[l_1_cand_m1].DeltaR(gen_l_1) < 0.1 ){
    if( GenMatching(reco_lep[l_1_cand_m1], gen_l_1, 0.1, 0.05) ){
      FillHist("GEN_pt_order_first", 1, 1., 0., 2., 2);

      int l_2_cand_m1, l_3_cand_m1;
      if( fabs( (reco_lep[OppSign]+reco_MET).M() - 80.4 ) < fabs( (reco_lep[l_SS_rem]+reco_MET).M() - 80.4 ) ){
        l_3_cand_m1 = OppSign;
        l_2_cand_m1 = l_SS_rem;
      }
      else{
        l_3_cand_m1 = l_SS_rem;
        l_2_cand_m1 = OppSign;
      }
      //if( gen_l_2.DeltaR( reco_lep[l_2_cand_m1] ) < 0.1 && gen_l_3.DeltaR( reco_lep[l_3_cand_m1] ) < 0.1 ) FillHist("GEN_pt_order_first_mlmet_next", 1, 1., 0., 2., 2);
      if( GenMatching(gen_l_2, reco_lep[l_2_cand_m1], 0.1, 0.05) && GenMatching(gen_l_3, reco_lep[l_3_cand_m1], 0.1, 0.05) ) FillHist("GEN_pt_order_first_mlmet_next", 1, 1., 0., 2., 2);
      else FillHist("GEN_pt_order_first_mlmet_next", 0, 1., 0., 2., 2);

    }
    else{
      FillHist("GEN_pt_order_first", 0, 1., 0., 2., 2);
    }

    //==== 2) mlmet first = m2
    
    //if( gen_l_3.DeltaR(reco_lep[l_3_cand]) < 0.1 ){
    if( GenMatching(gen_l_3, reco_lep[l_3_cand], 0.1, 0.05) ){
      FillHist("GEN_mlmet_first", 1, 1., 0., 2., 2);

      int l_1_cand_m2, l_2_cand_m2;
      if( l_3_cand == OppSign ){
        if( signal_class == 3){
          l_1_cand_m2 = SameSign[0];
          l_2_cand_m2 = SameSign[1];
        }
        else{
          l_1_cand_m2 = SameSign[1];
          l_2_cand_m2 = SameSign[0];
        }
      }
      else{
        l_2_cand_m2 = OppSign;
        if( l_3_cand == SameSign[0] ) l_1_cand_m2 = SameSign[1];
        else l_1_cand_m2 = SameSign[0];
      }

      //if( gen_l_1.DeltaR( reco_lep[l_1_cand_m2] ) < 0.1 && gen_l_2.DeltaR( reco_lep[l_2_cand_m2] ) < 0.1 ) FillHist("GEN_mlmet_first_pt_order_next", 1, 1., 0., 2., 2);
      if( GenMatching(gen_l_1, reco_lep[l_1_cand_m2], 0.1, 0.05) && GenMatching(gen_l_2, reco_lep[l_2_cand_m2], 0.1, 0.05) ) FillHist("GEN_mlmet_first_pt_order_next", 1, 1., 0., 2., 2);
      else FillHist("GEN_mlmet_first_pt_order_next", 0, 1., 0., 2., 2);

    }
    else{
      FillHist("GEN_mlmet_first", 0, 1., 0., 2., 2);
    }


  }

  //==== histograms
 
  FillHist("GEN_matching_validation_W_pri", (gen_l_1+gen_l_2+gen_l_3+gen_nu).M(), 1., 0., 1100., 110);
  FillHist("GEN_matching_validation_HN", (gen_l_2+gen_l_3+gen_nu).M(), 1., 0., 1100., 110);
  FillHist("GEN_matching_validation_W_sec", (gen_l_3+gen_nu).M(), 1., 0., 100., 100);

  FillHist("GEN_gen_l_1_Pt", gen_l_1.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_l_2_Pt", gen_l_2.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_l_3_Pt", gen_l_3.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_nu_Pt", gen_nu.Pt(), 1., 0., 1500., 1500);
 
  snu::KParticle gen_l_SS;
  if( gen_l_1.Charge() == gen_l_2.Charge() ) gen_l_SS = gen_l_2;
  else gen_l_SS = gen_l_3;

  FillHist("GEN_gen_l_SS_Pt", gen_l_SS.Pt(), 1., 0., 1500., 1500);
  
  //==== check in gen level
  if( gen_l_1.Pt() > gen_l_SS.Pt() ) FillHist("GEN_gen_pri_lep_pt_greater", 1, 1., 0., 2., 2);
  else FillHist("GEN_gen_pri_lep_pt_greater", 0, 1., 0., 2., 2);

  FillHist("TEST_leadingSS_pt", reco_lep[SameSign[0]].Pt(), 1., 0., 1500., 1500);
  FillHist("TEST_DeltaR_gen_l_1_AND_leadingSS", reco_lep[SameSign[0]].DeltaR(gen_l_1), 1., 0., 6., 60);
  FillHist("TEST_DeltaR_gen_l_1_AND_subleadingSS", reco_lep[SameSign[1]].DeltaR(gen_l_1), 1., 0., 6., 60);

  //if( reco_lep[SameSign[0]].DeltaR(gen_l_1) < 0.1 ){
  if( GenMatching(reco_lep[SameSign[0]], gen_l_1, 0.1, 0.05) ){
    FillHist("GEN_reco_leading_SS_match_gen_l_1", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_reco_leading_SS_match_gen_l_1", 0, 1., 0., 2., 2); 
  }

  //if( reco_lep[SameSign[1]].DeltaR(gen_l_1) < 0.1 ){
  if( GenMatching(reco_lep[SameSign[1]], gen_l_1, 0.1, 0.05) ){
    FillHist("GEN_reco_subleading_SS_match_gen_l_1", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_reco_subleading_SS_match_gen_l_1", 0, 1., 0., 2., 2);
  }

  if( gen_l_1.Charge() == gen_l_2.Charge() ){

    //if( reco_lep[SameSign[0]].DeltaR(gen_l_1) < 0.1 ){
    if( GenMatching(reco_lep[SameSign[0]], gen_l_1, 0.1, 0.05) ){
      FillHist("l1l2SS_GEN_reco_leading_SS_match_gen_l_1", 1, 1., 0., 2., 2);
    }
    else{
      FillHist("l1l2SS_GEN_reco_leading_SS_match_gen_l_1", 0, 1., 0., 2., 2);
    }

    if( gen_l_1.Pt() > gen_l_2.Pt() ) FillHist("l1l2SS_gen_l_1_leading", 1, 1., 0., 2., 2);
    else FillHist("l1l2SS_gen_l_1_leading", 0, 1., 0., 2., 2);

  }

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

void trilepton_mumumu_ntp::find_decay(std::vector<snu::KTruth> truthcoll, int target_index, std::vector<int>& indices){

  for(unsigned int i=target_index+1; i<truthcoll.size(); i++){ 
    if( truthcoll.at(i).IndexMother() == target_index && truthcoll.at(i).PdgId() == truthcoll.at(target_index).PdgId() ){
      indices.push_back(i);
      find_decay(truthcoll, i, indices);
    }
  }

}

void trilepton_mumumu_ntp::print_all_indices(TString particle, std::vector<int> vec){

  cout << particle+" indices" << endl;
  for(unsigned int i=0; i<vec.size(); i++) cout << " " << vec.at(i) << endl;

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
