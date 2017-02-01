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
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu::trilepton_mumumu() :  AnalyzerCore(),
n_gen_pass(0), sol_sel_chi2_best(0), sol_sel_chi2_plus(0), sol_sel_chi2_minus(0), sol_sel_chi2_smaller(0), sol_sel_chi2_larger(0),
allgenfound(false),
out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu");
  
  Message("In trilepton_mumumu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

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

  if( k_sample_name.Contains("HN") ){
    find_genparticles();
    if(allgenfound){
      FillHist("GenFound", 1., 1., 0., 2., 2);
    }
    else{
      FillHist("GenFound", 0., 1., 0., 2., 2);
    }
  }

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

  bool DoCutOp = std::find(k_flags.begin(), k_flags.end(), "cutop") != k_flags.end();

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

  double this_RelIso = 0.4;
  std::vector<snu::KMuon> muontriLooseColl;
  bool diboson_had = std::find(k_flags.begin(), k_flags.end(), "diboson_had") != k_flags.end();
  //==== signal
  if( k_sample_name.Contains("HN") ){
    //==== save gen particles @ snu::KParticle gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;

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
  }
  //==== diboson, but hadronica decays
  else if( (k_sample_name.Contains("WZ") || k_sample_name.Contains("ZZ") || k_sample_name.Contains("WW") ) && diboson_had ){
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
    if(muontriLooseColl.size()==3) return;
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, true);
  }
  //==== otherwise
  else{
    muontriLooseColl = GetHNTriMuonsByLooseRelIso(this_RelIso, false);
  }

  //CorrectMuonMomentum(muontriLooseColl); FIXME do this for v8-0-4
  double muon_id_iso_sf= MuonScaleFactor("MUON_HN_TRI_TIGHT", muontriLooseColl, 0);
  double MuTrkEffSF =  MuonTrackingEffScaleFactor(muontriLooseColl);

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", true, 30., 2.4);

  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  int n_bjets=0;
  for(int j=0; j<njet; j++){
    if(jetColl_hn.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)){
      n_bjets++;
      FillHist("TEST_bjet_pt", jetColl_hn.at(j).Pt(), 1., 0., 200., 200);
    }
  }

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
    weight*=muon_id_iso_sf;
    //weight*=weight_trigger_sf;
    weight*=MuTrkEffSF;
    weight*=trigger_ps_weight;
    weight*=pileup_reweight;
  }

  int n_triTight_muons(0);
  for(unsigned int i=0; i<muontriLooseColl.size(); i++){
    if(muontriLooseColl.at(i).RelIso04() < 0.1) n_triTight_muons++;
  }
  int n_triLoose_muons = muontriLooseColl.size();

  //==== ppp to TTL/TLL/LLL ?
  if(!k_isdata){
    if(n_triLoose_muons==3){
      FillHist("PPP_nTight", n_triTight_muons, weight, 0., 4., 4);
    }
  }

  FillHist("n_loose_muon", n_triLoose_muons, 1., 0., 10., 10);
  FillHist("n_tight_muon", n_triTight_muons, 1., 0., 10., 10);

  if( n_triLoose_muons != 3 ) return;
  if( n_triTight_muons != 3 ) return;

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

  FillHist("lowosllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("lowosllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowssllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[OppSign]+lep[SameSign[0]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[OppSign]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);
  FillHist("lowllmass", ( lep[SameSign[0]]+lep[SameSign[1]] ).M(), 1., 0., 20., 200);

  //==== MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
  //==== MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
  //==== POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
  if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
      (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) return;
  FillCutFlow("mllsf4", 1.);

  if(k_sample_name.Contains("HN") && allgenfound) solution_selection_stduy(muontriLooseColl);

  ///////////////////////////////////////////
  ////////// m(HN) < 80 GeV region //////////
  ///////////////////////////////////////////

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
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

  bool VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
  if(!VetoZResonance) return;
  FillCutFlow("ZVeto", 1.);

  if(n_bjets>0) return;
  FillCutFlow("bjetVeto", 1.);

  //==== preselection is done

  if(DoCutOp){
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
    cutop[11] = muontriLooseColl.at(0).dXY();
    cutop[12] = muontriLooseColl.at(1).dXY();
    cutop[13] = muontriLooseColl.at(2).dXY();
    cutop[14] = muontriLooseColl.at(0).dZ();
    cutop[15] = muontriLooseColl.at(1).dZ();
    cutop[16] = muontriLooseColl.at(2).dZ();
    cutop[17] = muontriLooseColl.at(0).RelIso04();
    cutop[18] = muontriLooseColl.at(1).RelIso04();
    cutop[19] = muontriLooseColl.at(2).RelIso04();
    cutop[20] = W_sec.M();
    cutop[21] = MET;
    cutop[22] = 0.; // weight_err

    FillNtp("cutop",cutop);
    return;
  }

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
  FillCLHist(trilephist, "cut0", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
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
    FillCLHist(trilephist, "cutWlow", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
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
    FillCLHist(trilephist, "cutWhigh", eventbase->GetEvent(), muontriLooseColl, electronColl, jetColl_hn, weight);
    FillHist("n_events_cutWhigh", 0, weight, 0., 1., 1);

    FillCutFlow("HighMass", 1.);

  }


   return;
}// End of execute event loop
  


void trilepton_mumumu::EndCycle()throw( LQError ){
  
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
    AnalyzerCore::MakeHistograms("cutflow", 11,0.,11.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3muon");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS"); 
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"mllsf4");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"ZVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"bjetVeto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"LowMass");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"HighMass");
    
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

  MakeNtp("cutop", "first_pt:second_pt:third_pt:deltaR_OS_min:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:weight:first_dXY:second_dXY:third_dXY:first_dZ:second_dZ:third_dZ:first_RelIso:second_RelIso:third_RelIso:W_sec_highmass_mass:PFMET:weight_err");

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

void trilepton_mumumu::find_genparticles(){

  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  bool OldSig = k_sample_name.Contains("VmuN_0p1");
  int PID_HN = 80000002;
  if(!OldSig) PID_HN = 9900012;

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
  vector<int> gen_HN_indices, gen_W_sec_indices, gen_l_1_indices, gen_nu_indices, gen_l_3_indices, gen_l_2_indices;
  gen_HN_indices.clear();
  gen_W_sec_indices.clear();
  gen_l_1_indices.clear();
  gen_nu_indices.clear();
  gen_l_3_indices.clear();
  gen_l_2_indices.clear();
  bool isLowMass(true);

  //==== check if this is low/high mass region
  if(GetSignalMass() < 80) isLowMass = true;
  else isLowMass = false;
 
  //==== find HN index
  for(int i=2;i<truthmax;i++){
    if(abs(truthColl.at(i).PdgId()) == PID_HN){
      gen_HN_indices.push_back(i);
      find_decay(truthColl, i, gen_HN_indices);
      break;
    }
  }
  //print_all_indices("gen_HN", gen_HN_indices);
  if(gen_HN_indices.size() == 0){
    Message("[GEN] HN not found", INFO);
    return;
  }
  int HN_motherindex = truthColl.at( gen_HN_indices.at(0) ).IndexMother();
  FillHist("GEN_HN_mother_pdgid", abs(truthColl.at(HN_motherindex).PdgId()), 1., 0., 30., 30);

  //======================
  //==== low mass region
  //======================
  if(isLowMass){
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
    if(gen_l_1_indices.size() == 0){
      Message("[GEN][LowMass] l_1 not found", INFO);
      return;
    }

    //==== As m(HN) goes closer to W mass, W_sec starting to appear in truth coll
    //==== if W_sec is virtual, it may not appear in the gen particle collections
    //==== check W_sec exists
    bool LowMass_Wsec(false), LowMass_Z(false);
    vector<int> gen_Z_indices;
    gen_Z_indices.clear();

    for(int i=2;i<truthmax;i++){
      //==== 1) |PID| = 24
      //==== 2) Mother's PID = PID_HN
      if( fabs(truthColl.at(i).PdgId()) == 24 && (abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId())) == PID_HN ){
        gen_W_sec_indices.push_back(i);
        find_decay(truthColl, i, gen_W_sec_indices);
        LowMass_Wsec = true;
        break;
      }
    }
    //print_all_indices("gen_W_sec", gen_W_sec_indices);

    //==== For low mass, HN can also decays via virtual Z
    for(int i=2;i<truthmax;i++){
      //==== 1) |PID| = 23
      //==== 2) Mother's PID = PID_HN
      if( fabs(truthColl.at(i).PdgId()) == 23 && (abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId())) == PID_HN ){
        gen_Z_indices.push_back(i);
        find_decay(truthColl, i, gen_Z_indices);
        LowMass_Z = true;
        break;
      }
    }
    //print_all_indices("gen_Z", gen_Z_indices);

    //============================
    //==== 1) HN to three particles
    //============================

    if( (!LowMass_Wsec) && (!LowMass_Z) ){
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
      if(gen_nu_indices.size() == 0){
        Message("[GEN][LowMass][1] nu not found", INFO);
        return;
      }

      //==== Let l_1 and l_2 are SS (considering dilep channel :D)
      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = PID of l_1 (SS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0){
        Message("[GEN][LowMass][1] l_2 not found", INFO);
        return;
      }

      //==== Let l_1 and l_3 are OS
      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) PID = - (PID of l_1) (OS)
          //==== 2) Mother = {HN}
          if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_HN_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_l_3_indices);
            break;
          }
        }
        if(gen_l_3_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0){
        Message("[GEN][LowMass][1] l_3 not found", INFO);
        return;
      }
    }

    //=========================
    //==== 2) HN decays via Z
    //=========================

    if( LowMass_Z ){
      //==== find nu at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_HN_indices.size();j++){
          //==== 1) |PID| = 14
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
      if(gen_nu_indices.size() == 0){
        Message("[GEN][LowMass][2] nu not found", INFO);
        return;
      }

      //==== Let l_1 and l_2 are SS (condiering dilep channel :D)
      //==== find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_Z_indices.size();j++){
          //==== 1) PID = PID of l_1 (SS)
          //==== 2) Mother = {Z}
          if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_Z_indices.at(j) ){
            gen_l_2_indices.push_back(i);
            find_decay(truthColl, i, gen_l_2_indices);
            break;
          }
        }
      }
      //print_all_indices("gen_l_2", gen_l_2_indices);
      if(gen_l_2_indices.size() == 0){
        Message("[GEN][LowMass][2] l_2 not found", INFO);
        return;
      }

      //==== Let l_1 and l_3 are OS
      //==== find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        for(unsigned int j=0;j<gen_Z_indices.size();j++){
          //==== 1) PID = - (PID of l_1) (OS)
          //==== 2) Mother = {Z}
          if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_indices.at(0)).PdgId() && truthColl.at(i).IndexMother() == gen_Z_indices.at(j) ){
            gen_l_3_indices.push_back(i);
            find_decay(truthColl, i, gen_l_3_indices);
            break;
          }
        }
        if(gen_l_3_indices.size() != 0) break;
      }
      //print_all_indices("gen_l_3", gen_l_3_indices);
      if(gen_l_3_indices.size() == 0){
        Message("[GEN][LowMass][2] l_3 not found", INFO);
        return;
      }

    }

    //====================================
    //==== 3) HN decays in to on-shell W
    //====================================

    if(LowMass_Wsec){
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
      if(gen_nu_indices.size() == 0){
        Message("[GEN][LowMass][3] nu not found", INFO);
        return;
      }

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
      if(gen_l_3_indices.size() == 0){
        Message("[GEN][LowMass][3] l_3 not found", INFO);
        return;
      }

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
      if(gen_l_2_indices.size() == 0){
        Message("[GEN][LowMass][3] l_3 not found", INFO);
        return;
      }
      
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
    gen_W_sec = truthColl.at( gen_W_sec_indices.back() );
    gen_l_3 = truthColl.at( gen_l_3_indices.back() );
    gen_nu = truthColl.at( gen_nu_indices.back() );
    allgenfound = true;

  }

}

void trilepton_mumumu::solution_selection_stduy(std::vector<snu::KMuon> recomuons){

  GENMatchingdR = 0.1;
  GENMatchingdPt = 0.05;

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
  if(GetSignalMass() < 80) isLowMass = true;

  //======================
  //==== low mass region
  //======================
  if(isLowMass){

    //==== solution selection
    double pz_sol_lowmass[2];
    pz_sol_lowmass[0] = solveqdeq(80.385, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "m"); // 0 = minus
    pz_sol_lowmass[1] = solveqdeq(80.385, reco_lep[0]+reco_lep[1]+reco_lep[2], MET, METphi, "p"); // 1 = plus
    if( pz_sol_lowmass[0] != pz_sol_lowmass[1] ){
      n_gen_pass++;
      int best_sel = fabs(pz_sol_lowmass[0]-gen_nu.Pz()) < fabs(pz_sol_lowmass[1]-gen_nu.Pz()) ? 0 : 1;
      int smaller = fabs(pz_sol_lowmass[0]) < fabs(pz_sol_lowmass[1]) ? 0 : 1;
      int larger = smaller == 0 ? 1 : 0;

      //=================
      //==== Chi2 study
      //=================

      FillHist("GEN_solsel_chi2_best", pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_solsel_chi2_plus", pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_solsel_chi2_minus", pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_solsel_chi2_smaller", pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);
      FillHist("GEN_solsel_chi2_larger", pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2), 1., 0., 10000., 10000);

      if(best_sel == 0) FillHist("GEN_solsel_chi2_minus_0_plus_1", 0, 1., 0., 2., 2);
      else              FillHist("GEN_solsel_chi2_minus_0_plus_1", 1, 1., 0., 2., 2);
      if(best_sel == smaller) FillHist("GEN_solsel_chi2_smaller_0_larger_1", 0, 1., 0., 2., 2);
      else                    FillHist("GEN_solsel_chi2_smaller_0_larger_1", 1, 1., 0., 2., 2);

      sol_sel_chi2_best += pow( (gen_nu.Pz() - pz_sol_lowmass[best_sel])/gen_nu.Pz() , 2);
      sol_sel_chi2_plus += pow( (gen_nu.Pz() - pz_sol_lowmass[1])/gen_nu.Pz() , 2);
      sol_sel_chi2_minus += pow( (gen_nu.Pz() - pz_sol_lowmass[0])/gen_nu.Pz() , 2);
      sol_sel_chi2_smaller += pow( (gen_nu.Pz() - pz_sol_lowmass[smaller])/gen_nu.Pz() , 2);
      sol_sel_chi2_larger += pow( (gen_nu.Pz() - pz_sol_lowmass[larger])/gen_nu.Pz() , 2);

      //===================
      //==== deltaR study
      //===================

      PutNuPz(&reco_MET, pz_sol_lowmass[best_sel]);
      FillHist("GEN_solsel_dR_best", gen_nu.DeltaR( reco_MET  ), 1., 0., 5., 50);
      PutNuPz(&reco_MET, pz_sol_lowmass[1]);
      FillHist("GEN_solsel_dR_plus", gen_nu.DeltaR( reco_MET  ), 1., 0., 5., 50);
      PutNuPz(&reco_MET, pz_sol_lowmass[0]);
      FillHist("GEN_solsel_dR_minus", gen_nu.DeltaR( reco_MET  ), 1., 0., 5., 50);
      PutNuPz(&reco_MET, pz_sol_lowmass[smaller]);
      FillHist("GEN_solsel_dR_smaller", gen_nu.DeltaR( reco_MET  ), 1., 0., 5., 50);
      PutNuPz(&reco_MET, pz_sol_lowmass[larger]);
      FillHist("GEN_solsel_dR_larger", gen_nu.DeltaR( reco_MET  ), 1., 0., 5., 50);

    }


  }

  //=======================
  //==== high mass region
  //=======================
  else{

    snu::KParticle reco_lep_tlv[3];
    for(int i=0; i<3; i++) reco_lep_tlv[i] = reco_lep[i];
    int l_3_cand = find_mlmet_closest_to_W(reco_lep_tlv, reco_MET);

    FillHist("GEN_highmass_reco_MET", reco_MET.Pt(), 1., 0., 120., 120);
    FillHist("GEN_highmass_MT_gen_l_1_MET", (gen_l_1 + reco_MET).M() - 80.385, 1., -60., 60., 120);
    FillHist("GEN_highmass_MT_gen_l_2_MET", (gen_l_2 + reco_MET).M() - 80.385, 1., -60., 60., 120);
    FillHist("GEN_highmass_MT_gen_l_3_MET", (gen_l_3 + reco_MET).M() - 80.385, 1., -60., 60., 120);
    FillHist("GEN_highmass_l_3_cand", l_3_cand, 1., 0., 3., 3);
    FillHist("GEN_highmass_gen_W_sec_pt", gen_W_sec.Pt(), 1., 0., 1000., 1000);
    FillHist("GEN_highmass_dR_gen_l_1_gen_nu", gen_l_1.DeltaR(gen_nu), 1., 0., 5., 50);
    FillHist("GEN_highmass_dR_gen_l_2_gen_nu", gen_l_2.DeltaR(gen_nu), 1., 0., 5., 50);
    FillHist("GEN_highmass_dR_gen_l_3_gen_nu", gen_l_3.DeltaR(gen_nu), 1., 0., 5., 50);

    //==== 1) l_1 first by pt ordering
    //==== signal_class = 3 : gen_l_1 is leading SS
    //==== signal_class = 4 : gen_l_1 is subleading SS (m > 1000 GeV)
    
    int l_1_cand = SameSign[0], l_SS_rem = SameSign[1];
    int signal_class = 3;
    //==== pt ordering reversed for m(HN) > 1000 GeV
    if( GetSignalMass() > 1000 ){ //FIXME study signal class
      signal_class = 4;
      l_1_cand = SameSign[1];
      l_SS_rem = SameSign[0];
    }
    //if( reco_lep[l_1_cand].DeltaR(gen_l_1) < 0.1 ){
    if( GenMatching(gen_l_1, reco_lep[l_1_cand], GENMatchingdR, GENMatchingdPt) ){
      FillHist("GEN_highmass_gen_l_1_first", 1, 1., 0., 2., 2);

      int l_2_cand_m1, l_3_cand_m1;
      //if( fabs( (reco_lep[OppSign]+reco_MET).M() - 80.385 ) < fabs( (reco_lep[l_SS_rem]+reco_MET).M() - 80.385 ) ){
      if( fabs( MT(reco_lep[OppSign], reco_MET) - 80.385 ) < fabs( MT(reco_lep[l_SS_rem], reco_MET) - 80.385 ) ){ 
        l_3_cand_m1 = OppSign;
        l_2_cand_m1 = l_SS_rem;
      }
      else{
        l_3_cand_m1 = l_SS_rem;
        l_2_cand_m1 = OppSign;
      }
      //if( gen_l_2.DeltaR( reco_lep[l_2_cand_m1] ) < 0.1 && gen_l_3.DeltaR( reco_lep[l_3_cand_m1] ) < 0.1 ) FillHist("GEN_pt_order_first_mlmet_next", 1, 1., 0., 2., 2);
      if( GenMatching(gen_l_2, reco_lep[l_2_cand_m1], GENMatchingdR, GENMatchingdPt) && GenMatching(gen_l_3, reco_lep[l_3_cand_m1], GENMatchingdR, GENMatchingdPt) ) FillHist("GEN_highmass_gen_l_1_first_mlmet_next", 1, 1., 0., 2., 2);
      else FillHist("GEN_highmass_gen_l_1_first_mlmet_next", 0, 1., 0., 2., 2);

    }
    else{
      FillHist("GEN_highmass_gen_l_1_first", 0, 1., 0., 2., 2);
    }

    //==== 2) gen_l_3 first by mlmet
    
    //if( gen_l_3.DeltaR(reco_lep[l_3_cand]) < 0.1 ){
    if( GenMatching(gen_l_3, reco_lep[l_3_cand], GENMatchingdR, GENMatchingdPt) ){
      FillHist("GEN_highmass_gen_l_3_first", 1, 1., 0., 2., 2);

      vector<snu::KParticle> lep_rem;
      for(unsigned int i=0; i<3; i++){
        if(i!=l_3_cand) lep_rem.push_back( reco_lep[i] );
      }

      //==== 2-1) Leading goes to l_1
      if( GenMatching(gen_l_1, lep_rem[0], GENMatchingdR, GENMatchingdPt) ){
        FillHist("GEN_highmass_gen_l_3_first_gen_l_1_leading_next", 1, 1., 0., 2., 2);
      }
      else{
        FillHist("GEN_highmass_gen_l_3_first_gen_l_1_leading_next", 0, 1., 0., 2., 2);
      }
      //==== 2-2) Leading goes to l_2
      if( GenMatching(gen_l_2, lep_rem[0], GENMatchingdR, GENMatchingdPt) ){
        FillHist("GEN_highmass_gen_l_3_first_gen_l_2_leading_next", 1, 1., 0., 2., 2);
      }
      else{
        FillHist("GEN_highmass_gen_l_3_first_gen_l_2_leading_next", 0, 1., 0., 2., 2);
      }

      //==== Combination
      int l_1_cand_m2, l_2_cand_m2;
      //==== if l_3 is OS : then, use pt ordering
      if( l_3_cand == OppSign ){
        //==== m(HN) : 90 ~ 200 GeV
        if( signal_class == 3){
          l_1_cand_m2 = SameSign[0];
          l_2_cand_m2 = SameSign[1];
        }
        //==== m(HN) : 300 ~ 1000 GeV
        else{
          l_1_cand_m2 = SameSign[1];
          l_2_cand_m2 = SameSign[0];
        }
      }
      //==== if l_3 is SS : then, remaing SS should be l_1
      else{
        l_2_cand_m2 = OppSign;
        if( l_3_cand == SameSign[0] ) l_1_cand_m2 = SameSign[1];
        else l_1_cand_m2 = SameSign[0];
      }

      //if( gen_l_1.DeltaR( reco_lep[l_1_cand_m2] ) < 0.1 && gen_l_2.DeltaR( reco_lep[l_2_cand_m2] ) < 0.1 ) FillHist("GEN_mlmet_first_pt_order_next", 1, 1., 0., 2., 2);
      if( GenMatching(gen_l_1, reco_lep[l_1_cand_m2], GENMatchingdR, GENMatchingdPt) && GenMatching(gen_l_2, reco_lep[l_2_cand_m2], GENMatchingdR, GENMatchingdPt) ) FillHist("GEN_highmass_gen_l_3_first_pt_order_next", 1, 1., 0., 2., 2);
      else FillHist("GEN_highmass_gen_l_3_first_pt_order_next", 0, 1., 0., 2., 2);

    }
    else{
      FillHist("GEN_highmass_gen_l_3_first", 0, 1., 0., 2., 2);
    }


  }

  //==== histograms
 
  FillHist("GEN_matching_validation_W_pri", (gen_l_1+gen_l_2+gen_l_3+gen_nu).M(), 1., 0., 1100., 1100);
  FillHist("GEN_matching_validation_HN", (gen_l_2+gen_l_3+gen_nu).M(), 1., 0., 1100., 1100);
  FillHist("GEN_matching_validation_W_sec", (gen_l_3+gen_nu).M(), 1., 0., 100., 1000);

  FillHist("GEN_gen_l_1_Pt", gen_l_1.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_l_2_Pt", gen_l_2.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_l_3_Pt", gen_l_3.Pt(), 1., 0., 1500., 1500);
  FillHist("GEN_gen_nu_Pt", gen_nu.Pt(), 1., 0., 1500., 1500);
 
  snu::KParticle gen_l_SS;
  if( gen_l_1.PdgId() == gen_l_2.PdgId() ) gen_l_SS = gen_l_2;
  else gen_l_SS = gen_l_3;
  FillHist("GEN_gen_SS_Pt", gen_l_SS.Pt(), 1., 0., 1500., 1500);

  //==== gen_l_1 : leadingSS
  if( gen_l_1.Pt() > gen_l_SS.Pt() ){
    FillHist("GEN_gen_l_1_leadingSS", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_gen_l_1_leadingSS", 0, 1., 0., 2., 2);
  }
  //==== gen_l_1 : leading
  if( gen_l_1.Pt() > gen_l_2.Pt() && gen_l_1.Pt() > gen_l_3.Pt() ){
    FillHist("GEN_gen_l_1_leading", 1., 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_gen_l_1_leading", 0., 1., 0., 2., 2);
  }
  //==== gen_l_2 : leading
  if( gen_l_2.Pt() > gen_l_1.Pt() && gen_l_2.Pt() > gen_l_3.Pt() ){
    FillHist("GEN_gen_l_2_leading", 1., 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_gen_l_2_leading", 0., 1., 0., 2., 2);
  }

  //if( reco_lep[SameSign[0]].DeltaR(gen_l_1) < 0.1 ){
  if( GenMatching(gen_l_1, reco_lep[SameSign[0]], GENMatchingdR, GENMatchingdPt) ){
    FillHist("GEN_reco_leading_SS_match_gen_l_1", 1, 1., 0., 2., 2);
    //cout << "=======================================" << endl;
    //cout << "== gen_l_1 is matched to leading SS  ==" << endl;
    //cout << gen_l_1.Pt() << '\t' << gen_l_1.PdgId()/-13 << endl;
    //cout << gen_l_2.Pt() << '\t' << gen_l_2.PdgId()/-13 << endl;
    //cout << gen_l_3.Pt() << '\t' << gen_l_3.PdgId()/-13 << endl;
    //cout << reco_lep[SameSign[0]].Pt() << '\t' << reco_lep[SameSign[0]].Charge() << endl;
    //cout << reco_lep[SameSign[1]].Pt() << '\t' << reco_lep[SameSign[0]].Charge() << endl;
  }
  else{
    FillHist("GEN_reco_leading_SS_match_gen_l_1", 0, 1., 0., 2., 2); 
  }

  //if( reco_lep[SameSign[1]].DeltaR(gen_l_1) < 0.1 ){
  if( GenMatching(gen_l_1, reco_lep[SameSign[1]], GENMatchingdR, GENMatchingdPt) ){
    FillHist("GEN_reco_subleading_SS_match_gen_l_1", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_reco_subleading_SS_match_gen_l_1", 0, 1., 0., 2., 2);
  }

  //if( reco_lep[0].DeltaR(gen_l_1) < 0.1 ){
  if( GenMatching(gen_l_1, reco_lep[0], GENMatchingdR, GENMatchingdPt) ){
    FillHist("GEN_reco_leading_match_gen_l_1", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_reco_leading_match_gen_l_1", 0, 1., 0., 2., 2);
  }

  //if( reco_lep[0].DeltaR(gen_l_2) < 0.1 ){
  if( GenMatching(gen_l_2, reco_lep[0], GENMatchingdR, GENMatchingdPt) ){
    FillHist("GEN_reco_leading_match_gen_l_2", 1, 1., 0., 2., 2);
  }
  else{
    FillHist("GEN_reco_leading_match_gen_l_2", 0, 1., 0., 2., 2);
  }


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

void trilepton_mumumu::find_decay(std::vector<snu::KTruth> truthcoll, int target_index, std::vector<int>& indices){

  for(unsigned int i=target_index+1; i<truthcoll.size(); i++){ 
    if( truthcoll.at(i).IndexMother() == target_index && truthcoll.at(i).PdgId() == truthcoll.at(target_index).PdgId() ){
      indices.push_back(i);
      find_decay(truthcoll, i, indices);
    }
  }

}

void trilepton_mumumu::print_all_indices(TString particle, std::vector<int> vec){

  cout << particle+" indices" << endl;
  for(unsigned int i=0; i<vec.size(); i++) cout << " " << vec.at(i) << endl;

}

std::vector<snu::KMuon> trilepton_mumumu::sort_muons_ptorder(std::vector<snu::KMuon> muons){

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


