// $Id: trilepton_mumumu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond     <jalmond@cern.ch>        - SNU
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
  *  This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu::trilepton_mumumu() :  AnalyzerCore(), out_muons(0)  {
  
  
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
  // You can out put messages simply with Message function. Message( "comment", output_level)  output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  //// default is silver, excpept in v7-6-2 ... set to gold if you use gold json in analysis
  /// only available in v7-6-X branch and newer
  //lumimask = snu::KEvent::gold;

  return;
}


void trilepton_mumumu::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  weight*=MCweight;
  
  /// Acts on data to remove bad reconstructed event 
  if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;
  

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment
  
  if(!PassBasicEventCuts()) return;    /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  //ListTriggersAvailable();

  if(!PassTrigger(triggerslist, prescale)) return;
  FillCutFlow("TriggerCut", weight);
  // Trigger matching is done using KMuon::TriggerMatched(TString) which returns a bool

  /* // #### CAT::: trigger matching information is stored for muons and electrons for:
  ///HLT_IsoMu24_eta2p1_v
  ///HLT_Mu17_Mu8_DZ_v
  ///HLT_Mu17_TkMu8_DZ_v
  ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
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
  //  ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", weight);
  

  /// List of preset muon collections : Can call also POGSoft/POGLoose/POGMedium/POGTight
  std::vector<snu::KMuon> muontriColl = GetMuons(BaseSelection::MUON_HN_TRI); 
  
  
  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets
  std::vector<snu::KJet> jetColl_hn       = GetJets(BaseSelection::JET_HN);// pt > 20 ; eta < 2.5; PFlep veto; pileup ID
  
  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 
  std::vector<snu::KElectron> electronColl         = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronLooseColl      = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  
  int n_jet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , n_jet*MCweight + MCweight*0.1, 1., -6. , 6., 12);

  numberVertices = eventbase->GetEvent().nVertices();  
  
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    if(eventbase->GetEvent().CatVersion().empty()) pileup_reweight = reweightPU->GetWeight(int(eventbase->GetEvent().nVertices()), k_mcperiod);
    else if(eventbase->GetEvent().CatVersion().find("v7-4") !=std::string::npos)  pileup_reweight = reweightPU->GetWeight(int(eventbase->GetEvent().nVertices()), k_mcperiod);
    else   pileup_reweight = eventbase->GetEvent().PileUpWeight(lumimask,snu::KEvent::central);
  }

  FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10);
  
  int n_loose_muon = muontriColl.size();
  if( n_loose_muon != 3 ) return;
  FillCutFlow("3loosemuon", weight);

  snu::KParticle lep[3], nu, selection_nu[2], HN, W_on_shell, gamma_star, z_candidate;
  for(int i=0;i<3;i++){
    lep[i] = muontriColl.at(i);
  }

  int OppSign, LepCand[2]; // LepCand[0].Pt() > LepCand[1].Pt()
  if(lep[0].Charge() * lep[1].Charge() > 0){ // Q(0) = Q(1)
    if(lep[1].Charge() * lep[2].Charge() < 0){ // Q(1) != Q(2)
      OppSign = 2;
      LepCand[0] = 0;
      LepCand[1] = 1;
    }
    else return; // veto Q(0) = Q(1) = Q(2)
  }
  else{ // Q(0) != Q(1)
    if(lep[0].Charge() * lep[2].Charge() > 0){ // Q(0) = Q(2)
      OppSign = 1;
      LepCand[0] = 0;
      LepCand[1] = 2;
    }
    else if(lep[1].Charge() * lep[2].Charge() > 0){ // Q(1) = Q(2)
      OppSign = 0;
      LepCand[0] = 1;
      LepCand[1] = 2;
    }
  } // Find l2 and assign l1&l3 in ptorder 
  FillCutFlow("2SS1OS", weight);

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.MET(), METphi = Evt.METPhi();
  nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
  PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
  PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus

  int solution_selection;
  if( selection_nu[0].Pz() == selection_nu[1].Pz() ){ // solution selection => smaller
    solution_selection = 0; // 0, 1 상관 없으므로
  }
  else{
    // take the one with smaller magnitude
    if(fabs(selection_nu[0].Pz()) > fabs(selection_nu[1].Pz())){
      solution_selection = 1;
    }
    else{
      solution_selection = 0;
    }
  }

  // reconstruct HN and W_real 4-vec with selected Pz solution
  PutNuPz(&nu, selection_nu[solution_selection].Pz());
  HN = lep[1] + lep[2] + nu;
  W_on_shell = HN + lep[0];

  double deltaR_OS_min;
  if( lep[0].DeltaR(lep[1]) < lep[2].DeltaR(lep[1]) ){
    deltaR_OS_min = lep[0].DeltaR(lep[1]);
    gamma_star = lep[1] + lep[0];
  }
  else{
    deltaR_OS_min = lep[2].DeltaR(lep[1]);
    gamma_star = lep[1] + lep[2];
  }

  if( fabs( (lep[0]+lep[1]).M() - 91.1876 ) <
      fabs( (lep[2]+lep[1]).M() - 91.1876 )   ){
    z_candidate = lep[1] + lep[0];
  }
  else{
    z_candidate = lep[1] + lep[2];
  }

  SetBinInfo(0);
  FillHist("cut0_HN_mass", HN.M(), weight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
  FillHist("cut0_W_on_shell_mass", W_on_shell.M(), weight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
  FillHist("cut0_deltaR_OS_min", deltaR_OS_min, weight, 0, 5, 5./0.1);
  FillHist("cut0_gamma_star_mass", gamma_star.M(), weight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
  FillHist("cut0_z_candidate_mass", z_candidate.M(), weight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
  FillHist("cut0_n_jet", n_jet, weight, 0, 10, 10);
  FillHist("cut0_HN_mass_PU", HN.M(), weight*pileup_reweight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
  FillHist("cut0_W_on_shell_mass_PU", W_on_shell.M(), weight*pileup_reweight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
  FillHist("cut0_deltaR_OS_min_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 5./0.1);
  FillHist("cut0_gamma_star_mass_PU", gamma_star.M(), weight*pileup_reweight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
  FillHist("cut0_z_candidate_mass_PU", z_candidate.M(), weight*pileup_reweight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
  FillHist("cut0_n_jet_PU", n_jet, weight*pileup_reweight, 0, 10, 10);
  if(deltaR_OS_min > 0.5){
    SetBinInfo(1);
    FillHist("cutdR_HN_mass", HN.M(), weight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
    FillHist("cutdR_W_on_shell_mass", W_on_shell.M(), weight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
    FillHist("cutdR_deltaR_OS_min", deltaR_OS_min, weight, 0, 5, 5./0.1);
    FillHist("cutdR_gamma_star_mass", gamma_star.M(), weight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
    FillHist("cutdR_z_candidate_mass", z_candidate.M(), weight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
    FillHist("cutdR_n_jet", n_jet, weight, 0, 10, 10);
    FillHist("cutdR_HN_mass_PU", HN.M(), weight*pileup_reweight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
    FillHist("cutdR_W_on_shell_mass_PU", W_on_shell.M(), weight*pileup_reweight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
    FillHist("cutdR_deltaR_OS_min_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 5./0.1);
    FillHist("cutdR_gamma_star_mass_PU", gamma_star.M(), weight*pileup_reweight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
    FillHist("cutdR_z_candidate_mass_PU", z_candidate.M(), weight*pileup_reweight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
    FillHist("cutdR_n_jet_PU", n_jet, weight*pileup_reweight, 0, 10, 10);
    if(W_on_shell.M() < 100){
      SetBinInfo(2);
      FillHist("cutdR_cutW_HN_mass", HN.M(), weight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
      FillHist("cutdR_cutW_W_on_shell_mass", W_on_shell.M(), weight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
      FillHist("cutdR_cutW_deltaR_OS_min", deltaR_OS_min, weight, 0, 5, 5./0.1);
      FillHist("cutdR_cutW_gamma_star_mass", gamma_star.M(), weight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
      FillHist("cutdR_cutW_z_candidate_mass", z_candidate.M(), weight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
      FillHist("cutdR_cutW_n_jet", n_jet, weight, 0, 10, 10);
      FillHist("cutdR_cutW_HN_mass_PU", HN.M(), weight*pileup_reweight, HN_x_min, HN_x_max, (HN_x_max-HN_x_min)/HN_dx);
      FillHist("cutdR_cutW_W_on_shell_mass_PU", W_on_shell.M(), weight*pileup_reweight, W_on_shell_x_min, W_on_shell_x_max, (W_on_shell_x_max-W_on_shell_x_min)/W_on_shell_dx);
      FillHist("cutdR_cutW_deltaR_OS_min_PU", deltaR_OS_min, weight*pileup_reweight, 0, 5, 5./0.1);
      FillHist("cutdR_cutW_gamma_star_mass_PU", gamma_star.M(), weight*pileup_reweight, gamma_star_x_min, gamma_star_x_max, (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx);
      FillHist("cutdR_cutW_z_candidate_mass_PU", z_candidate.M(), weight*pileup_reweight, z_candidate_x_min, z_candidate_x_max, (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx);
      FillHist("cutdR_cutW_n_jet_PU", n_jet, weight*pileup_reweight, 0, 10, 10);
    }
  }

  
  return;
}// End of execute event loop
  


void trilepton_mumumu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

trilepton_mumumu::~trilepton_mumumu() {
  
  Message("In trilepton_mumumu Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_mumumu::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
   GetHist("cutflow")->Fill(cut,weight);
  
  }
  else{
   AnalyzerCore::MakeHistograms("cutflow", 6,0.,6.);

   GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3loosemuon");
   GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"2SS1OS");  
   
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


void trilepton_mumumu::SetBinInfo(int cut){
  if(cut==0){ // no cut
    HN_x_min=0; HN_x_max=500; HN_dx=10;
    W_on_shell_x_min=70; W_on_shell_x_max=500; W_on_shell_dx=10;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1;
  }
  if(cut==1){ // "dR_OS_min > 0.5"
    HN_x_min=0; HN_x_max=500; HN_dx=10;
    W_on_shell_x_min=70; W_on_shell_x_max=500; W_on_shell_dx=10;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1;
  }
  if(cut==2){ // "dR_OS_min > 0.5", "W_on_shell < 100 GeV"
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_on_shell_x_min=70; W_on_shell_x_max=120; W_on_shell_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }
  if(cut==3){
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_on_shell_x_min=70; W_on_shell_x_max=120; W_on_shell_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }
  if(cut==4){
    HN_x_min=0; HN_x_max=100; HN_dx=10;
    W_on_shell_x_min=70; W_on_shell_x_max=120; W_on_shell_dx=5;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5;
  }

}
