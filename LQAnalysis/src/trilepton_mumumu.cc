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
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  return;
}


void trilepton_mumumu::ExecuteEvents()throw( LQError ){
  
  
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
 
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    pileup_reweight = reweightPU->GetWeight(eventbase->GetEvent().PileUpInteractionsTrue())* MCweight; // others + HNXX_new
    //pileup_reweight = MCweight;   // HN40, HN50, HN60, Wtollln
   }
  
	FillCutFlow("NoCut", weight*pileup_reweight);

  if(!PassBasicEventCuts()) return;     /// Initial event cuts
  FillCutFlow("EventCut", weight*pileup_reweight);
  
  std::vector<TString> triggerslist;
  triggerslist.push_back("HLT_Mu17_TkMu8_v");
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
  std::vector<snu::KMuon> muonLooseColl;
  eventbase->GetMuonSel()->HNLooseMuonSelection(muonLooseColl);
  std::vector<snu::KMuon> muonVetoColl;
  eventbase->GetMuonSel()->HNVetoMuonSelection(muonVetoColl);
  /// Correct the muon momentum with rochester corrections
  CorrectMuonMomentum(muonTightColl);
  CorrectMuonMomentum(muonHighPtColl);
  //CorrectMuonMomentum(muonLooseColl);
	//CorrectMuonMomentum(muonVetoColl);  

  std::vector<snu::KElectron> electronTightColl;
	 eventbase->GetElectronSel()->HNTightElectronSelection(electronTightColl);
  
  std::vector<snu::KElectron> electronVetoColl;
	 eventbase->GetElectronSel()->HNVetoElectronSelection(electronVetoColl);
  
  std::vector<snu::KJet> jetColl_lepveto;
  eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
  eventbase->GetJetSel()->SetPt(20.);
  eventbase->GetJetSel()->SetEta(2.5);
  eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonTightColl, electronTightColl);
	//eventbase->GetJetSel()->JetHNSelection(jetColl_lepveto, muonLooseColl, electronLooseColl); 
 
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  double Event[100];
  double Muon[100];
	double Jet[100];
  double Gen[100]; 

	int n_loose_muon = muonLooseColl.size(), n_jet = jetColl_lepveto.size();
 
  if( muonLooseColl.size() != 3 ) return;   // three loose muon cut (with iso 0.1)
  FillCutFlow("3Leptons", weight*pileup_reweight);
 
/* 
  snu::KParticle lep[3];
  for(int i=0;i<3;i++){
    lep[i] = muonLooseColl.at(i);
  }
*/
/*
	//TEST//
	snu::KMuon test = muonLooseColl.at(0);
	FillHist("Pt()", test.Pt(), 1, 0., 300., 300);
	FillHist("CocktailPt()", test.MuonCocktailPt(), 1, 0., 300., 300);
	FillHist("MSPt()", test.MuonMSPt(), 1, 0., 300., 300);
	FillHist("IDPt()", test.MuonIDPt(), 1, 0., 300., 300);
*/
/*
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
*/
  Muon[0] = muonLooseColl.size();
  snu::KMuon mu;
  double LeptonRelIso;
	int n_muonvar = 20;
  for(int i=0; i<3; i++){
    //mu = muonLooseColl.at(order[i]);
    mu = muonLooseColl.at(i);
    LeptonRelIso = (mu.SumIsoCHDR03() + std::max(0.0, mu.SumIsoNHDR03() + mu.SumIsoPHDR03() - 0.5*mu.SumPUIsoR03()))/mu.Pt();
    Muon[1+n_muonvar*i] = mu.Px();
    Muon[2+n_muonvar*i] = mu.Py();
    Muon[3+n_muonvar*i] = mu.Pz();
    Muon[4+n_muonvar*i] = mu.E();
    Muon[5+n_muonvar*i] = LeptonRelIso;
    Muon[6+n_muonvar*i] = mu.Eta();
    Muon[7+n_muonvar*i] = mu.dXY();
    Muon[8+n_muonvar*i] = mu.dZ();
    Muon[9+n_muonvar*i] = mu.GlobalChi2();
    Muon[10+n_muonvar*i] = mu.IsoHcalVeto();
    Muon[11+n_muonvar*i] = mu.IsoEcalVeto();
    Muon[12+n_muonvar*i] = mu.IsGlobal();
    Muon[13+n_muonvar*i] = mu.validHits();
    Muon[14+n_muonvar*i] = mu.validPixHits();
    Muon[15+n_muonvar*i] = mu.validStations();
    Muon[16+n_muonvar*i] = mu.ActiveLayer();
    Muon[17+n_muonvar*i] = mu.IsPF();
    Muon[18+n_muonvar*i] = mu.IsTracker();
    Muon[19+n_muonvar*i] = mu.Charge();
    Muon[20+n_muonvar*i] = mu.Pt();
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
 
  FillNtp("Muon", Muon);
  FillNtp("Event", Event);
	FillNtp("Jet", Jet);


  ////////////////////////
  /// Generation Level ///
  ////////////////////////
  
  // l_1 : lepton from first W
  // l_2 : lepton from HN
  // l_3 : lepton from second W
  // W_real : on-shell W (First W for low mass, second W for high mass)
  
  int truthmax = truthColl.size();
  int gen_l_1_index, gen_l_2_index, gen_l_3_index, gen_W_real_index, gen_HN_index, gen_W_virtual_index, gen_nu_index;
  snu::KParticle gen_nu, gen_W_real, gen_HN, gen_W_virtual;
  snu::KParticle gen_l_1, gen_l_2, gen_l_3;
  bool isWV=false, isLowMass = true;
  
  // check if this is low/high mass region //
  // find HN index
  for(int i=2;i<truthmax;i++){
    if(truthColl.at(i).PdgId() == 80000002){
      if(truthColl.at(i).M() > 80) isLowMass = false;
      gen_HN_index = i;
      gen_HN = truthColl.at(i);
      break;
    }
  }
  
  // low mass region //
  if(isLowMass){
    // find W_real index
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 24){
        gen_W_real_index = i;
        gen_W_real = truthColl.at(i);
        break;
      }
    }
    // check W_virtual exists
    for(int i=2;i<truthmax;i++){
      if(fabs(truthColl.at(i).PdgId()) == 24 && truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 80000002){
        gen_W_virtual_index = i;
        gen_W_virtual = truthColl.at(i);
        isWV = true;
        break;
      }
    }
    // no W_virtual in truthColl index case
    if(!isWV){
      // find l_1 at gen. level
      for(int i=2;i<truthmax;i++){
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){
          gen_l_1 = truthColl.at(i);
          gen_l_1_index = i;
          break;
        }
      }
      // find nu at gen. level
      for(int i=2;i<truthmax;i++){
        if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_HN_index){
          gen_nu = truthColl.at(i);
          break;
        }
      }
      // find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        if(truthColl.at(i).PdgId() == truthColl.at(gen_l_1_index).PdgId() && truthColl.at(i).IndexMother() == gen_HN_index){
          gen_l_3 = truthColl.at(i);
          break;
        }
      }
      // find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        if(truthColl.at(i).PdgId() == -truthColl.at(gen_l_1_index).PdgId() && truthColl.at(i).IndexMother() == gen_HN_index){
          gen_l_2 = truthColl.at(i);
          break;
        }
      }
    }
    // W_virtual in truthColl index case
    else{
      // find l_1 at gen. level
      for(int i=2;i<truthmax;i++){
        if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_W_real_index){
          gen_l_1 = truthColl.at(i);
          gen_l_1_index = i;
          break;
        }
      }
      // find nu at gen. level
      for(int i=2;i<truthmax;i++){
        if(abs(truthColl.at(i).PdgId()) == 14 && truthColl.at(i).IndexMother() == gen_W_virtual_index){
          gen_nu = truthColl.at(i);
          break;
        }
      }
      // find l_3 at gen. level
      for(int i=2;i<truthmax;i++){
        if(fabs(truthColl.at(i).PdgId()) == fabs(truthColl.at(gen_l_1_index).PdgId()) && truthColl.at(i).IndexMother() == gen_W_virtual_index){
          gen_l_3 = truthColl.at(i);
          gen_l_3_index = i;
          break;
        }
      }
      // find l_2 at gen. level
      for(int i=2;i<truthmax;i++){
        if(fabs(truthColl.at(i).PdgId()) == fabs(truthColl.at(gen_l_1_index).PdgId()) && truthColl.at(i).IndexMother() == gen_HN_index){
          gen_l_2 = truthColl.at(i);
        }
      }
      if(truthColl.at(gen_l_3_index).PdgId() == -truthColl.at(gen_l_1_index).PdgId()){
        snu::KParticle gen_temp = gen_l_2;
        gen_l_2 = gen_l_3;
        gen_l_3 = gen_temp;
      }
    }
  }
  // high mass region //
  else{
    // find l_1 at gen. level
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == truthColl.at(gen_HN_index).IndexMother()){
        gen_l_1 = truthColl.at(i);
        gen_l_1_index = i;
        break;
      }
    }
    // fine l_2 at gen. level
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == gen_HN_index ){
        gen_l_2 = truthColl.at(i);
        gen_l_2_index = i;
      }
    }
    // find W_real
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 24 && truthColl.at(i).IndexMother() == truthColl.at(gen_l_2_index).IndexMother()){
        gen_W_real = truthColl.at(i);
        gen_W_real_index =i;
        break;
      }
    }
    // find nu at gen. level
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 14){
        gen_nu_index = i;
        gen_nu = truthColl.at(i) ;
        break;
      }
    }
    // find l_3 at gen. level
    for(int i=2;i<truthmax;i++){
      if(abs(truthColl.at(i).PdgId()) == 13 && truthColl.at(i).IndexMother() == truthColl.at(gen_nu_index).IndexMother()){
        gen_l_3 = truthColl.at(i);
        gen_l_3_index = i;
        break;
      }
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

	return;
  
}// End of execute event loop


void trilepton_mumumu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  
}


void trilepton_mumumu::BeginCycle() throw( LQError ){
  
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
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"3Leptons");
    //GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"OppSign");
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
  
  MakeNtp("Muon", "n_LooseMuon:0_Px:0_Py:0_Pz:0_E:0_LeptonRelIso:0_Eta:0_dXY:0_dZ:0_GlobalChi2:0_IsoHcalVeto:0_IsoEcalVeto:0_IsGlobal:0_validHits:0_validPixHits:0_validStations:0_ActiveLayer:0_IsPf:0_IsTracker:0_Charge:0_Pt:1_Px:1_Py:1_Pz:1_E:1_LeptonRelIso:1_Eta:1_dXY:1_dZ:1_GlobalChi2:1_IsoHcalVeto:1_IsoEcalVeto:1_IsGlobal:1_validHits:1_validPixHits:1_validStations:1_ActiveLayer:1_IsPf:1_IsTracker:1_Charge:1_Pt:2_Px:2_Py:2_Pz:2_E:2_LeptonRelIso:2_Eta:2_dXY:2_dZ:2_GlobalChi2:2_IsoHcalVeto:2_IsoEcalVeto:2_IsGlobal:2_validHits:2_validPixHits:2_validStations:2_ActiveLayer:2_IsPf:2_IsTracker:2_Charge:2_Pt");
	MakeNtp("Event", "MET:METphi:n_bjet:weight:PU_reweight:MCweight");
  MakeNtp("Gen", "mu_0_Px:mu_0_Py:mu_0_Pz:mu_0_E:mu_1_Px:mu_1_Py:mu_1_Pz:mu_1_E:mu_2_Px:mu_2_Py:mu_2_Pz:mu_2_E:nu_Px:nu_Py:nu_Pz:nu_E:W_Px:W_Py:W_Pz:W_E:HN_Px:HN_Py:HN_Pz:HN_E");
  MakeNtp("Jet","n_Jet:0_Px:0_Py:0_Pz:0_E:0_Eta:0_Pt:0_PFJetTrackCountingHighPurBTag:0_BtagProb:0_CombinedSecVertexBtag:1_Px:1_Py:1_Pz:1_E:1_Eta:1_Pt:1_PFJetTrackCountingHighPurBTag:1_BtagProb:1_CombinedSecVertexBtag:2_Px:2_Py:2_Pz:2_E:2_Eta:2_Pt:2_PFJetTrackCountingHighPurBTag:2_BtagProb:2_CombinedSecVertexBtag:3_Px:3_Py:3_Pz:3_E:3_Eta:3_Pt:3_PFJetTrackCountingHighPurBTag:3_BtagProb:3_CombinedSecVertexBtag:4_Px:4_Py:4_Pz:4_E:4_Eta:4_Pt:4_PFJetTrackCountingHighPurBTag:4_BtagProb:4_CombinedSecVertexBtag:5_Px:5_Py:5_Pz:5_E:5_Eta:5_Pt:5_PFJetTrackCountingHighPurBTag:5_BtagProb:5_CombinedSecVertexBtag:6_Px:6_Py:6_Pz:6_E:6_Eta:6_Pt:6_PFJetTrackCountingHighPurBTag:6_BtagProb:6_CombinedSecVertexBtag:7_Px:7_Py:7_Pz:7_E:7_Eta:7_Pt:7_PFJetTrackCountingHighPurBTag:7_BtagProb:7_CombinedSecVertexBtag:8_Px:8_Py:8_Pz:8_E:8_Eta:8_Pt:8_PFJetTrackCountingHighPurBTag:8_BtagProb:8_CombinedSecVertexBtag:9_Px:9_Py:9_Pz:9_E:9_Eta:9_Pt:9_PFJetTrackCountingHighPurBTag:9_BtagProb:9_CombinedSecVertexBtag"); 
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



