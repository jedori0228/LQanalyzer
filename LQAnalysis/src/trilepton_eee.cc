// $Id: trilepton_eee.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_eee Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_eee.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_eee);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
trilepton_eee::trilepton_eee() :  AnalyzerCore() {
  
  
  // To have the correct name in the log:
  SetLogName("trilepton_eee");
  
  Message("In trilepton_eee constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
}


void trilepton_eee::InitialiseAnalysis() throw( LQError ) {
  
  
  /// Initialise histograms
  MakeHistograms();
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  return;
}


void trilepton_eee::ExecuteEvents()throw( LQError ){
  
  
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    pileup_reweight = reweightPU->GetWeight(eventbase->GetEvent().PileUpInteractionsTrue())* MCweight; // others
    //pileup_reweight = MCweight;   // HN40, HN50, HN60, w_llln
	 }
  
  FillCutFlow("NoCut", weight*pileup_reweight);
  
  if(!PassBasicEventCuts()) return;     /// Initial event cuts
  FillCutFlow("EventCut", weight*pileup_reweight);
  
  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  if(!PassTrigger(triggerslist, prescale)) return;
  FillCutFlow("TriggerCut", weight*pileup_reweight);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("VertexCut", weight*pileup_reweight);
  
	 numberVertices = eventbase->GetEvent().nVertices();
  

  
  std::vector<snu::KMuon> muonTightColl;
  eventbase->GetMuonSel()->HNTightMuonSelection(muonTightColl);
  std::vector<snu::KMuon> muonHighPtColl;
  eventbase->GetMuonSel()->HNTightHighPtMuonSelection(muonHighPtColl);
  /// Correct the muon momentum with rochester corrections
  CorrectMuonMomentum(muonTightColl);
  CorrectMuonMomentum(muonHighPtColl);
  std::vector<snu::KMuon> muonLooseColl;
  eventbase->GetMuonSel()->HNLooseMuonSelection(muonLooseColl);
  std::vector<snu::KMuon> muonVetoColl;
  eventbase->GetMuonSel()->HNVetoMuonSelection(muonVetoColl);
  
  std::vector<snu::KElectron> electronTightColl;
  eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);
  std::vector<snu::KElectron> electronVetoColl;
  eventbase->GetElectronSel()->HNVetoElectronSelection(electronVetoColl);
  std::vector<snu::KElectron> electronLooseColl;
  eventbase->GetElectronSel()->HNLooseElectronSelection(electronLooseColl);
  
  std::vector<snu::KJet> jetColl_lepveto;
  eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
  eventbase->GetJetSel()->SetPt(20.);
  eventbase->GetJetSel()->SetEta(2.5);
  eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonTightColl, electronTightColl);
  
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  double Event[100];
  double Electron[100];
  double Jet[100];
  double Gen[100];
  
  int n_loose_electron = electronLooseColl.size(), n_jet = jetColl_lepveto.size();
  
  if( n_loose_electron != 3 ) return;   // three loose elecron cut
  FillCutFlow("3Leptons", weight*pileup_reweight);
  
  snu::KParticle lep[3];
  for(int i=0;i<3;i++){
    lep[i] = electronLooseColl.at(i);
  }
  int OppSign, LepCand[2];
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
  FillCutFlow("OppSign", weight*pileup_reweight); // one opposite charge cut
  
  int order[3] = {LepCand[0], OppSign, LepCand[1]};
  Electron[0] = n_loose_electron;
  snu::KElectron el;
  double PHONH_03[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14};
  double rho = eventbase->GetEvent().JetRho();
  double ElectronIsoDR03, LeptonRelIsoDR03;
  int ifid, n_elvar = 13;
  float id_scalefactor;
  
  for(int i=0; i<3; i++){
				el = electronLooseColl.at(order[i]);
				if (fabs(el.SCEta()) < 1.0) ifid = 0;
        else if (fabs(el.SCEta()) < 1.479) ifid = 1;
        else if (fabs(el.SCEta()) < 2.0) ifid = 2;
        else if (fabs(el.SCEta()) < 2.2) ifid = 3;
        else if (fabs(el.SCEta()) < 2.3) ifid = 4;
        else if (fabs(el.SCEta()) < 2.4) ifid = 5;
        else ifid = 6;
				ElectronIsoDR03 =  el.PFChargedHadronIso03() + max( el.PFNeutralHadronIso03() + el.PFPhotonIso03() - rho * PHONH_03[ifid],  0.);
				LeptonRelIsoDR03 = ElectronIsoDR03/  el.Pt();
				Electron[1+i*n_elvar] = el.HasMatchedConvPhot();
				Electron[2+i*n_elvar] = el.MissingHits();
				Electron[3+i*n_elvar] = LeptonRelIsoDR03;
				Electron[4+i*n_elvar] = el.GsfCtfScPixChargeConsistency();
				Electron[5+i*n_elvar] = el.Eta();
				Electron[6+i*n_elvar] = el.Px();
				Electron[7+i*n_elvar] = el.Py();
				Electron[8+i*n_elvar] = el.Pz();
				Electron[9+i*n_elvar] = el.E();
				Electron[10+i*n_elvar] = el.dxy();
    
    id_scalefactor = 1.0;
    if(!isData) id_scalefactor *=  ElectronScaleFactor(el.Eta(), el.Pt());
    Electron[11+i*n_elvar] = id_scalefactor;
    Electron[12+i*n_elvar] = el.Charge();
    Electron[13+i*n_elvar] = el.Pt();
  }
  
  int n_bjet = 0;
  for(int i=0; i<n_jet; i++){
    if(jetColl_lepveto.at(i).CombinedSecVertexBtag() > 0.679){
      n_bjet++;
    }
  }
  
  snu::KEvent Evt = eventbase->GetEvent();
  double METv = Evt.PFMET();
  double METphi = Evt.PFMETphi();
  Event[0] = Evt.PFMET();
  Event[1] = Evt.PFMETphi();
  Event[2] = n_bjet;
  Event[3] = weight;
  Event[4] = pileup_reweight;
  Event[5] = MCweight;
  
  snu::KJet jet;
  int n_jetvar = 9;
  Jet[0] = n_jet; // <- only n_loose_electron <= 4 will be saved in the ntuple
  for(int i=0; i<n_jet; i++){
    jet = jetColl_lepveto.at(i);
    Jet[1+i*n_jetvar] = jet.Px();
    Jet[2+i*n_jetvar] = jet.Py();
    Jet[3+i*n_jetvar] = jet.Pz();
    Jet[4+i*n_jetvar] = jet.E();
    Jet[5+i*n_jetvar] = jet.Eta();
    Jet[6+i*n_jetvar] = jet.Pt();
    Jet[7+i*n_jetvar] = jet.PFJetTrackCountingHighPurBTag();
    Jet[8+i*n_jetvar] = jet.BtagProb();
    Jet[9+i*n_jetvar] = jet.CombinedSecVertexBtag();
  }
  
  FillNtp("Electron", Electron);
  FillNtp("Event", Event);
  FillNtp("Jet", Jet);
  
  /*
   ////////////////////////
   /// Generation Level ///
   ////////////////////////
   
   int truthmax = truthColl.size();
   int gen_l_1_index, gen_W_real_index, gen_HN_index, gen_W_virtual_index, gen_l_3_index;
   snu::KParticle gen_nu, gen_W_real, gen_HN, gen_W_virtual;
   snu::KParticle gen_l_1, gen_l_2, gen_l_3;
   bool isWV=false;
   
   
   
   for(int i=2;i<truthmax;i++){ // find W_real index
   if(abs(truthColl.at(i).PdgId()) == 24){
   gen_W_real_index = i;
   gen_W_real = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){ // find HN index
   if(truthColl.at(i).PdgId() == 80000002){
   gen_HN_index = i;
   gen_HN = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){ // check W_virtual exists
   if(fabs(truthColl.at(i).PdgId()) == 24 && truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 80000002){
   gen_W_virtual_index = i;
   gen_W_virtual = truthColl.at(i);
   isWV = true;
   break;
   }
   }
   
   if(!isWV){ // no W_virtual in truthColl index case
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){ // find l_1 at gen. level
   gen_l_1 = truthColl.at(i);
   gen_l_1_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_HN_index){ // find nu at gen. level
   gen_nu = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_index).PdgId() && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_3 at gen. level
   gen_l_3 = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_index).PdgId() && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_2 at gen. level
   gen_l_2 = truthColl.at(i);
   break;
   }
   }
   }
   else{ // W_virtual in truthColl index case
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){ // find l_1 at gen. level
   gen_l_1 = truthColl.at(i);
   gen_l_1_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_virtual_index){ // find nu at gen. level
   gen_nu = truthColl.at(i);
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(fabs(truthColl.at(i).PdgId()) == fabs(truthColl.at(gen_l_1_index).PdgId()) && truthColl.at(i).IndexMother() == gen_W_virtual_index){ // find l_3 at gen. level
   gen_l_3 = truthColl.at(i);
   gen_l_3_index = i;
   break;
   }
   }
   for(int i=2;i<truthmax;i++){
   if(fabs(truthColl.at(i).PdgId()) == fabs(truthColl.at(gen_l_1_index).PdgId()) && truthColl.at(i).IndexMother() == gen_HN_index){ // find l_2 at gen. level
   gen_l_2 = truthColl.at(i);
   }
   }
   if(truthColl.at(gen_l_3_index).PdgId() == -truthColl.at(gen_l_1_index).PdgId()){
   snu::KParticle gen_temp = gen_l_2;
   gen_l_2 = gen_l_3;
   gen_l_3 = gen_temp;
   }
   }
   
   Gen[0] = gen_l_1.Px();
   Gen[1] = gen_l_1.Py();
   Gen[2] = gen_l_1.Pz();
   Gen[3] = gen_l_1.E();
   Gen[4] = gen_l_2.Px();
   Gen[5] = gen_l_2.Py();
   Gen[6] = gen_l_2.Pz();
   Gen[7] = gen_l_2.E();
   Gen[8] = gen_l_3.Px();
   Gen[9] = gen_l_3.Py();
   Gen[10] = gen_l_3.Pz();
   Gen[11] = gen_l_3.E();
   Gen[12] = gen_nu.Px();
   Gen[13] = gen_nu.Py();
   Gen[14] = gen_nu.Pz();
   Gen[15] = gen_nu.E();
   Gen[16] = gen_W_real.Px();
   Gen[17] = gen_W_real.Py();
   Gen[18] = gen_W_real.Pz();
   Gen[19] = gen_W_real.E();
   Gen[20] = gen_HN.Px();
   Gen[21] = gen_HN.Py();
   Gen[22] = gen_HN.Pz();
   Gen[23] = gen_HN.E();
   
   
   FillNtp("Gen", Gen);
   */
  return;
  
}// End of execute event loop


void trilepton_eee::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  
}


void trilepton_eee::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "MyDataPileupHistogram.root").c_str());
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");
  
  
  return;
  
}

trilepton_eee::~trilepton_eee() {
  
  Message("In trilepton_eee Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void trilepton_eee::FillCutFlow(TString cut, float weight){
  
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
    
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 6,0.,6.);
    
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3Leptons");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OppSign");
  }
}


void trilepton_eee::BeginEvent( )throw( LQError ){
  
  Message("In BeginEvent() " , DEBUG);
  
  return;
}



void trilepton_eee::MakeHistograms(){
  //// Additional plots to make
  
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  
  MakeNtp("Electron", "n_LooseElectron:0_HasMatchedConvPhot:0_MissingHits:0_LeptonRelIsoDR03:0_GsfCtfScPixChargeConsistency:0_Eta:0_Px:0_Py:0_Pz:0_E:0_dxy:0_id_scalefactor:0_Charge:0_Pt:1_HasMatchedConvPhot:1_MissingHits:1_LeptonRelIsoDR03:1_GsfCtfScPixChargeConsistency:1_Eta:1_Px:1_Py:1_Pz:1_E:1_dxy:1_id_scalefactor:1_Charge:1_Pt:2_HasMatchedConvPhot:2_MissingHits:2_LeptonRelIsoDR03:2_GsfCtfScPixChargeConsistency:2_Eta:2_Px:2_Py:2_Pz:2_E:2_dxy:2_id_scalefactor:2_Charge:2_Pt");
  MakeNtp("Event", "MET:METphi:n_bjet:weight:PU_reweight:MCweight");
  //  MakeNtp("Gen", "mu_0_Px:mu_0_Py:mu_0_Pz:mu_0_E:mu_1_Px:mu_1_Py:mu_1_Pz:mu_1_E:mu_2_Px:mu_2_Py:mu_2_Pz:mu_2_E:nu_Px:nu_Py:nu_Pz:nu_E:W_Px:W_Py:W_Pz:W_E:HN_Px:HN_Py:HN_Pz:HN_E");
    MakeNtp("Jet","n_Jet:0_Px:0_Py:0_Pz:0_E:0_Eta:0_Pt:0_PFJetTrackCountingHighPurBTag:0_BtagProb:0_CombinedSecVertexBtag:1_Px:1_Py:1_Pz:1_E:1_Eta:1_Pt:1_PFJetTrackCountingHighPurBTag:1_BtagProb:1_CombinedSecVertexBtag:2_Px:2_Py:2_Pz:2_E:2_Eta:2_Pt:2_PFJetTrackCountingHighPurBTag:2_BtagProb:2_CombinedSecVertexBtag:3_Px:3_Py:3_Pz:3_E:3_Eta:3_Pt:3_PFJetTrackCountingHighPurBTag:3_BtagProb:3_CombinedSecVertexBtag:4_Px:4_Py:4_Pz:4_E:4_Eta:4_Pt:4_PFJetTrackCountingHighPurBTag:4_BtagProb:4_CombinedSecVertexBtag:5_Px:5_Py:5_Pz:5_E:5_Eta:5_Pt:5_PFJetTrackCountingHighPurBTag:5_BtagProb:5_CombinedSecVertexBtag:6_Px:6_Py:6_Pz:6_E:6_Eta:6_Pt:6_PFJetTrackCountingHighPurBTag:6_BtagProb:6_CombinedSecVertexBtag:7_Px:7_Py:7_Pz:7_E:7_Eta:7_Pt:7_PFJetTrackCountingHighPurBTag:7_BtagProb:7_CombinedSecVertexBtag:8_Px:8_Py:8_Pz:8_E:8_Eta:8_Pt:8_PFJetTrackCountingHighPurBTag:8_BtagProb:8_CombinedSecVertexBtag:9_Px:9_Py:9_Pz:9_E:9_Eta:9_Pt:9_PFJetTrackCountingHighPurBTag:9_BtagProb:9_CombinedSecVertexBtag");
  /**
   *  Remove//Overide this trilepton_eeeCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void trilepton_eee::ClearOutputVectors() throw(LQError) {
  
  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list.
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



