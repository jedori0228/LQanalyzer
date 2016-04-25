// $id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiMuon Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNDiMuon.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNDiMuon);

/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNDiMuon::HNDiMuon() :  AnalyzerCore(),  out_electrons(0) {

  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiMuon");

  Message("In HNDiMuon constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  
  k_met=0., k_emumass=0., k_emujjmass=0., k_l1jjmass=0., k_l2jjmass=0., k_njet=0;
  k_nbjet_m=-1;
  k_emujjmass_lowmass=0., k_l1jjmass_lowmass=0., k_l2jjmass_lowmass=0.;
  k_weight=0.;
  k_l1pt=0., k_l2pt=0., k_j1pt=0., k_jjmass=0., k_jjmass_lowmass=0.;
  k_l1eta=0., k_l2eta=0.;

  MakeCleverHistograms(sighist,"TChannel");
  /*
  MakeCleverHistograms(sighist,"TriLepEECR");
  MakeCleverHistograms(sighist,"TriLepMMCR");
  MakeCleverHistograms(sighist,"TriLepCR");
  MakeCleverHistograms(sighist,"OS_2Jet");

  MakeCleverHistograms(sighist,"SS_1Jet");
  MakeCleverHistograms(sighist,"SS_0bjet");
  MakeCleverHistograms(sighist,"SS_bjet");
  MakeCleverHistograms(sighist,"SS_DiJet");
  MakeCleverHistograms(sighist,"SS_DiJet_iso1");
  MakeCleverHistograms(sighist,"SS_DiJet_iso2");
  MakeCleverHistograms(sighist,"SS_DiJet_up");
  MakeCleverHistograms(sighist,"SS_DiJet_down");
  MakeCleverHistograms(sighist,"SS_lowmass");
  MakeCleverHistograms(sighist,"SS_lowmass2");
  MakeCleverHistograms(sighist,"SS_lowmass3");
  MakeCleverHistograms(sighist,"SS_lowmass4");
  MakeCleverHistograms(sighist,"SS_lowmass_40");
  MakeCleverHistograms(sighist,"SS_lowmass_50");
  MakeCleverHistograms(sighist,"SS_lowmass_60");
  MakeCleverHistograms(sighist,"SS_lowmass_70");
  MakeCleverHistograms(sighist,"SS_lowmass_80");
  MakeCleverHistograms(sighist,"SS_lowmass_80_2");
  MakeCleverHistograms(sighist,"SS_lowmassCR");

  MakeCleverHistograms(sighist,"SS_highmass");
  MakeCleverHistograms(sighist,"SS_highmass1");
  MakeCleverHistograms(sighist,"SS_highmass1_90");
  MakeCleverHistograms(sighist,"SS_highmass_90");
  MakeCleverHistograms(sighist,"SS_highmass_90b");
  MakeCleverHistograms(sighist,"SS_highmass_90c");
  MakeCleverHistograms(sighist,"SS_highmass1_100");
  MakeCleverHistograms(sighist,"SS_highmass_100");
  MakeCleverHistograms(sighist,"SS_highmass_100b");
  MakeCleverHistograms(sighist,"SS_highmass_100c");
  MakeCleverHistograms(sighist,"SS_highmass_125");
  MakeCleverHistograms(sighist,"SS_highmass_125b");
  MakeCleverHistograms(sighist,"SS_highmass_125c");
  MakeCleverHistograms(sighist,"SS_highmass_150");
  MakeCleverHistograms(sighist,"SS_highmass_150b");
  MakeCleverHistograms(sighist,"SS_highmass_150c");
  MakeCleverHistograms(sighist,"SS_highmass_175");
  MakeCleverHistograms(sighist,"SS_highmass_200");
  MakeCleverHistograms(sighist,"SS_highmass_250");
  MakeCleverHistograms(sighist,"SS_highmass_300");
  MakeCleverHistograms(sighist,"SS_highmass_exc");
  MakeCleverHistograms(sighist,"SS_highmass_350");
  MakeCleverHistograms(sighist,"SS_highmass_400");
  MakeCleverHistograms(sighist,"SS_highmass_500");


  MakeCleverHistograms(sighist,"SS_highmassCR");


  MakeCleverHistograms(sighist,"SSemu_1Jet");
  MakeCleverHistograms(sighist,"SSemu_DiJet");
  MakeCleverHistograms(sighist,"SSemu_DiJet_up");
  MakeCleverHistograms(sighist,"SSemu_DiJet_down");
  MakeCleverHistograms(sighist,"SSmue_1Jet");
  MakeCleverHistograms(sighist,"SSmue_DiJet");
  MakeCleverHistograms(sighist,"SSmue_DiJet_up");
  MakeCleverHistograms(sighist,"SSmue_DiJet_down");

  MakeCleverHistograms(sighist,"SSmue_2Jet");
  MakeCleverHistograms(sighist,"SSmue_3Jet");
  MakeCleverHistograms(sighist,"SSmue_4Jet");
  MakeCleverHistograms(sighist,"SSmue_BJet");
  MakeCleverHistograms(sighist,"SSmue_0BJet");

  MakeCleverHistograms(sighist,"SSemu_2Jet");
  MakeCleverHistograms(sighist,"SSemu_3Jet");
  MakeCleverHistograms(sighist,"SSemu_4Jet");
  MakeCleverHistograms(sighist,"SSemu_BJet");
  MakeCleverHistograms(sighist,"SSemu_0BJet");
  */
}


void HNDiMuon::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   return;
}


void HNDiMuon::ExecuteEvents()throw( LQError ){
    


  k_met=0., k_emumass=0., k_emujjmass=0., k_l1jjmass=0., k_l2jjmass=0., k_njet=0;
  k_emujjmass_lowmass=0., k_l1jjmass_lowmass=0., k_l2jjmass_lowmass=0.;

  k_nbjet_m=-1;
  k_weight=0.;
  k_l1pt=0., k_l2pt=0., k_j1pt=0., k_jjmass=0.,k_jjmass_lowmass=0.;
  k_l1eta=0., k_l2eta=0.;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  if(!PassBasicEventCuts())  throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  
  FillEventCutFlow("EventCut", "",weight);
  /// Trigger List 
  std::vector<TString> triggerslist;  
  //triggerslist.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  //triggerslist.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  
  triggerslist.push_back("HLT_Mu17_TkMu8_v");
  triggerslist.push_back("HLT_Mu17_Mu8_v");
  if(!PassTrigger(triggerslist, prescale)) throw LQError( "Fails basic cuts",  LQError::SkipEvent );   

  //// if the trigger that fired the event is prescaled you can reweight the event accordingly using the variable prescale
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  FillEventCutFlow("VertexCut","", weight);
    
  TString fake_loose_region = "looseregion2";
  TString fake_loose_label = "HNTight_loosereg2";

  //// Get the collection of electrons
  std::vector<snu::KElectron> electronAnalysisColl = GetElectrons(true,  true, fake_loose_label , weight);
  std::vector<snu::KMuon> muons = GetMuons("tight");
  
  if(!isData){
    for(std::vector<snu::KMuon>::iterator it = muons.begin(); it != muons.end(); it++){
      weight *= MuonScaleFactor(it->Eta(), it->Pt(), 0);
    }
 
  }
  
  

  vector<snu::KTruth> truth =  eventbase->GetTruth();
  
  std::vector<snu::KElectron> electronVetoColl       = GetElectrons(false, false, "veto"); 
  std::vector<snu::KElectron> electronLooseColl      = GetElectrons(false, false, "loose"); 

  std::vector<snu::KMuon> muonVetoColl  = GetMuons("veto");

  FillHist("Nveto_electrons", electronVetoColl.size() ,1., 0. , 4.,4);
  FillHist("Nveto_muons", muonVetoColl.size() ,1., 0. , 4., 4);
  FillHist("Nveto_leptons", electronVetoColl.size() + muonVetoColl.size()  ,1., 0. , 4., 4);


  //m_logger << INFO << "Number of veto el = " << electronVetoColl.size() << LQLogger::endmsg;
  //m_logger << INFO << "Number of veto muons = " << muonVetoColl.size() << LQLogger::endmsg;

  /// JETS

  std::vector<snu::KJet> jetColl             = GetJets("NoLeptonVeto");
  std::vector<snu::KJet> jetColl_lepveto     = GetJets("ApplyLeptonVeto");
  std::vector<snu::KJet> jetColl_lepveto_mva = GetJets("ApplyPileUpID");

  //  m_logger << INFO << "Number of jets1  = " << jetColl.size() << LQLogger::endmsg;
  //m_logger << INFO << "Number of jets2  = " << jetColl_lepveto.size() << LQLogger::endmsg;

  //RunMCCLosureTestEMU("loosereg2", jetColl_lepveto_mva,"",weight);


  ///// count number of bjets in the event (using cvs medium WP)
  int nbjet_m=0;

  for(unsigned int ij=0; ij <jetColl_lepveto_mva.size(); ij++){
    if(jetColl_lepveto_mva.at(ij).CombinedSecVertexBtag() > 0.679) nbjet_m++;
  }

 
  if(!(electronAnalysisColl.size() == 0 && muons.size() == 2)) throw LQError( "Fails basic cuts",  LQError::SkipEvent );

  if(k_running_nonprompt&&isData){
   
    weight      *= Get_DataDrivenWeight_MM(muons);

  }
  
 
 
  if(! ((muons.at(0).Pt() > 20.) &&  (muons.at(1).Pt() > 15.))) throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  

  snu::KParticle mumu =   muons.at(0)+ muons.at(1);                                                                                                                                                                                                            

  if(mumu.M()  < 10.) throw LQError( "Fails basic cuts",  LQError::SkipEvent );                                                                                                                                                                                            


  if(!SameCharge(muons)) throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  
  if ((electronVetoColl.size() + muonVetoColl.size()) >2) throw LQError( "Fails basic cuts",  LQError::SkipEvent );  
  
  cout << "SS MM " << endl;
 
  if(jetColl_lepveto_mva.size() > 3 ) {
    bool has_forward_jet(false), has_back_jet(false);
    for(unsigned int ij = 0 ; ij < jetColl_lepveto_mva.size(); ij++){
      if(jetColl_lepveto_mva.at(ij).Eta() > 1.5) has_forward_jet=true;
      if(jetColl_lepveto_mva.at(ij).Eta() < -1.5) has_back_jet=true;
      cout << "Passes selection ll jjjj " << endl;
      for(unsigned int ij1=0; ij1 < jetColl_lepveto_mva.size(); ij1++){
	cout << jetColl_lepveto_mva.at(ij1).Eta() << endl;
      }
    }
    if(has_forward_jet && has_back_jet)  FillCLHist(sighist, "TChannel", eventbase->GetEvent(), muons,electronAnalysisColl,jetColl_lepveto_mva, weight);
  }
  
  return;

    
  }// End of exeucte event loop



float HNDiMuon::WeightCFEvent(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, bool runchargeflip, bool useoldrates){

  if(electrons.size()!=1) return 0.;
  if(muons.size()!=1) return 0.;
  
  if(runchargeflip) {
    if(electrons.at(0).Charge() != muons.at(0).Charge()) {
      float cf1=  CFRate(electrons.at(0), useoldrates);

      return  cf1;
    }// OS requirement
    else return 0.;
  }// cf requirement
  else {
    if(electrons.at(0).Charge() != muons.at(0).Charge()) return 0.;
  }
  
  return 1.;
  
}


void HNDiMuon::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << INFO << "Number of os mc events = " << m_os_Z_nw  << LQLogger::endmsg; 
  m_logger << INFO << "Number of os mc events (weighted) = " << m_os_Z  << LQLogger::endmsg; 
  m_logger << INFO << "Number of ss mc events = " << m_ss_Z_nw  << LQLogger::endmsg; 
  m_logger << INFO << "Number of ss mc events (weighted)= " << m_ss_Z  << LQLogger::endmsg; 
}


void HNDiMuon::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  
  if(!k_isdata) reweightPU = new Reweight((analysisdir + "MyDataPileupHistogram_69400.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  
  DeclareVariable(k_met, "met", "MyTree");
  DeclareVariable(k_emumass, "emu_mass", "MyTree");
  DeclareVariable(k_emujjmass, "emujj_mass", "MyTree");
  DeclareVariable(k_l1jjmass, "l1jj_mass", "MyTree");
  DeclareVariable(k_l2jjmass, "l2jj_mass", "MyTree");
  DeclareVariable(k_emujjmass_lowmass, "emujj_mass_lowmass", "MyTree");
  DeclareVariable(k_l1jjmass_lowmass, "l1jj_mass_lowmass", "MyTree");
  DeclareVariable(k_l2jjmass_lowmass, "l2jj_mass_lowmass", "MyTree");
  DeclareVariable(k_njet, "njet", "MyTree");
  DeclareVariable(k_nbjet_m , "nbjet_m",  "MyTree");
  DeclareVariable(k_jjmass, "jj_mass", "MyTree");
  DeclareVariable(k_jjmass_lowmass, "jj_mass_lowmass", "MyTree");
  DeclareVariable(k_l1pt, "l1_pt", "MyTree");
  DeclareVariable(k_l1eta, "l1_eta", "MyTree");
  DeclareVariable(k_l2eta, "l2_eta", "MyTree");
  DeclareVariable(k_l2pt, "l2_pt", "MyTree");

  DeclareVariable(k_j1pt, "jet1_pt", "MyTree");
  DeclareVariable(k_weight, "weight", "MyTree");

  return;
  
}

HNDiMuon::~HNDiMuon() {
  
  Message("In HNDiMuon Destructor" , INFO);
  if(!k_isdata)delete reweightPU;

 }
     

void HNDiMuon::FillEventCutFlow(TString cut, TString label , float weight){

  if(GetHist(label + "_eventcutflow")) {
    GetHist(label + "_eventcutflow")->Fill(cut,weight);

  }
  else{
    AnalyzerCore::MakeHistograms(label + "_eventcutflow",19,0.,19.);

    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(2,"NoCut_w");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(3,"eventcut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(5,"VertexCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(6,"DiEl");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(7,"eedR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(8,"SSDiEl");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(9,"SS_lepveto");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(10,"DiJet");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(11,"Presel");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(12,"Presel_noZ");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(13,"Presel_nobjet");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(14,"lowmass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(15,"lowmassCR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(16,"mediummass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(17,"mediummassCR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(18,"highmass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(19,"highmassCR");
  }
  
}

     
void HNDiMuon::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow",16,0.,16.);
    
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"SS_NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"SS_Tight");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"SS_Tight_convveto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"SS_Tight_d0veto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"SS_Tight_reliso");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"SS_Medium_chargeconst");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"SS_Tight_chargeconst");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"SS_Tight_noclosejet");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"SS_anal_el");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"Signal_anal");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"Signal_Tightlooseiso_d0");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"Signal_Mediumlooseiso_d0");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(13,"Signal_drcut1");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(14,"Signal_drcut2");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(15,"Signal_anal_dr1");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(16,"Signal_anal_dr2");
  }
}
     
     
void HNDiMuon::FillIsoCutFlow(TString cut, float weight){
       
  
  if(GetHist("isocutflow")) {
    GetHist("isocutflow")->Fill(cut,weight);
    
  }
  else{
    AnalyzerCore::MakeHistograms("isocutflow",36,0.,36.);
    
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(1,"iso_d0_03_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(2,"iso_d0_03_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(3,"iso_d0_03_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(4,"iso_d0_03_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(5,"iso_d0_03_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(6,"iso_d0_03_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(7,"iso_d0_03_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(8,"iso_d0_03_iso3_075");
    
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(10,"iso_d0_02_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(11,"iso_d0_02_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(12,"iso_d0_02_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(13,"iso_d0_02_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(14,"iso_d0_02_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(15,"iso_d0_02_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(16,"iso_d0_02_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(17,"iso_d0_02_iso3_075");

    GetHist("isocutflow")->GetXaxis()->SetBinLabel(19,"iso_d0_01_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(20,"iso_d0_01_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(21,"iso_d0_01_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(22,"iso_d0_01_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(23,"iso_d0_01_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(24,"iso_d0_01_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(25,"iso_d0_01_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(26,"iso_d0_01_iso3_075");

    GetHist("isocutflow")->GetXaxis()->SetBinLabel(28,"iso_d0_005_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(29,"iso_d0_005_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(30,"iso_d0_005_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(31,"iso_d0_005_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(32,"iso_d0_005_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(33,"iso_d0_005_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(34,"iso_d0_005_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(35,"iso_d0_005_iso3_075");

    
    
  }
}


void HNDiMuon::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNDiMuon::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiMuonCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNDiMuon::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}




//  LocalWords:  masscuts jetResdown
