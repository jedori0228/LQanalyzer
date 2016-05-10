// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzerDiElectron Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzerDiElectron.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ExampleAnalyzerDiElectron);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
ExampleAnalyzerDiElectron::ExampleAnalyzerDiElectron() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("ExampleAnalyzerDiElectron");

  Message("In ExampleAnalyzerDiElectron constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void ExampleAnalyzerDiElectron::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist

   MakeCleverHistograms(sighist_ee, "DiElectron");
   MakeCleverHistograms(sighist_ee, "DiElectron_noTM");
   MakeCleverHistograms(sighist_ee, "DiElectronID");
   MakeCleverHistograms(sighist_ee, "DiElectronIDRECO");
   MakeCleverHistograms(sighist_ee, "SIGNAL");
   MakeCleverHistograms(sighist_ee, "SS_SIGNAL");
   MakeCleverHistograms(sighist_ee, "OS_SIGNAL");
   MakeCleverHistograms(sighist_ee, "DiElectron_HLT23");
   MakeCleverHistograms(sighist_ee, "DiElectronNoPRW");
   MakeCleverHistograms(sighist_ee, "DiElectronNoTrigger");
   MakeCleverHistograms(sighist_ee, "DiElectron_Zpeak");
   MakeCleverHistograms(sighist_ee, "DiElectron_Zpeak_chargeconsistency");
   MakeCleverHistograms(sighist_ee, "DiElectron_EE");
   MakeCleverHistograms(sighist_ee, "DiElectron_EB");
   MakeCleverHistograms(sighist_ee, "DiElectron_BB");
   MakeCleverHistograms(sighist_ee, "DiElectron_BJet");
   MakeCleverHistograms(sighist_ee, "SSElectron");
   MakeCleverHistograms(sighist_ee, "SSElectron_DiJet");
   MakeCleverHistograms(sighist_ee, "SSElectron_DiJet_ChargeConsistency");

   MakeCleverHistograms(trilephist,"TriElectron");
   MakeCleverHistograms(trilephist,"TriElectron_nomet");
   MakeCleverHistograms(trilephist,"TriElectron_noB");


   /// only available in v7-6-X branch and newer
   //// default lumimask is silver ////
   //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
   ///When METNoHF isfixed the default will be back to silver
   /// set to gold if you want to use gold json in analysis
   /// To set uncomment the line below:
   //ResetLumiMask(snu::KEvent::gold);
   
   return;
}


void ExampleAnalyzerDiElectron::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  std::vector<snu::KElectron> electronColl_nocut_truth             =  GetElectrons(BaseSelection::ELECTRON_NOCUT);
  
  FillHist("Nelectrons" , electronColl_nocut_truth.size(), 1.,  0. , 5., 5);
  
  int n_prompt(0);
  int n_cf(0);
  int n_conv(0);
  int n_tau(0);
  int n_fake(0);
  for(unsigned int iel=0; iel < electronColl_nocut_truth.size(); iel++){
    if(electronColl_nocut_truth.at(iel).MCIsPrompt()) n_prompt++;
    else n_fake++;

    if(electronColl_nocut_truth.at(iel).MCIsCF()) n_cf++;
    if(electronColl_nocut_truth.at(iel).MCIsFromConversion()) n_conv++;
    if(electronColl_nocut_truth.at(iel).MCFromTau()) n_tau++;
    FillHist("el_pdgid", fabs(electronColl_nocut_truth.at(iel).MCMatchedPdgId()),  fabs(electronColl_nocut_truth.at(iel).MotherPdgId()), 1., 0., 1000., 1000, 0., 1000., 1000); 
  }
  
  FillHist("electron_prompt_breakdown", n_prompt,n_fake, 1., 0., 4, 4, 0., 4., 4);
  FillHist("electron_cf_breakdown", n_cf, n_prompt,1., 0., 4, 4, 0., 5., 5);
  FillHist("electron_conv_breakdown", n_conv,n_prompt, 1., 0., 4, 4, 0., 5., 5);
  FillHist("electron_tau_breakdown", n_tau,n_prompt, 1., 0., 4, 4, 0., 5., 5);
  
  /// Can count number of bjets using IsBTagged function in KJet class 
  int nbjet_just_using_discriminant=0;
  for(unsigned int ij =0; ij < GetJets(BaseSelection::JET_HN).size(); ij++){
    if(GetJets(BaseSelection::JET_HN).at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight)) nbjet_just_using_discriminant++;
  }

  
  /// Updated way to cound bjets using NBJet function
  /// NBJet counts number of bjets, but varies the value of btag disciminant as expained in 
  /// 2a) on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

  /// Allowed input for taggers are:
  /// snu::KJet::CSVv2
  /// snu::KJet::cMVAv2
  /// Allowed values for WP are:
  /// snu::KJet::Loose
  /// snu::KJet::Medium
  /// snu::KJet::Tight
  int nbjet_using_btagsf = NBJet(GetJets(BaseSelection::JET_HN), snu::KJet::CSVv2, snu::KJet::Tight);
  
  /// Can also check invidual jets using IsBTagged function in AnalyzerCore
  int nbjet_just_using_sf(0);
  for(unsigned int ij =0; ij < GetJets(BaseSelection::JET_HN).size(); ij++){
    if(IsBTagged(GetJets(BaseSelection::JET_HN).at(ij), snu::KJet::CSVv2, snu::KJet::Tight)) nbjet_just_using_sf++;
  }
  

  std::vector<snu::KMuon> muons = GetMuons(BaseSelection::MUON_POG_TIGHT);
  CorrectMuonMomentum(muons);

  /*if(n_prompt > 3){
    for(unsigned int iel=0; iel < electronColl_nocut_truth.size(); iel++){
      cout << "eta " << electronColl_nocut_truth.at(iel).Eta() << endl;
      cout << "phi " << electronColl_nocut_truth.at(iel).Phi() << endl;
      cout << "pt " << electronColl_nocut_truth.at(iel).Pt() << endl;
      cout << "isprompt= " << electronColl_nocut_truth.at(iel).MCIsPrompt() << endl;
      cout << "mc pdgid = " << electronColl_nocut_truth.at(iel).MCMatchedPdgId() << endl;
      cout << "mc index = " <<  electronColl_nocut_truth.at(iel).MCTruthIndex() << endl;
      cout << "mother pdgid = " <<  electronColl_nocut_truth.at(iel).MotherPdgId()<< endl;
      cout << "motherindex= " <<electronColl_nocut_truth.at(iel).MotherTruthIndex()<< endl;

      cout << "Mother type = " << electronColl_nocut_truth.at(iel).GetMotherType() << endl;
      cout << "Particle type = " << electronColl_nocut_truth.at(iel).GetParticleType() << endl;

      for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
	if(eventbase->GetTruth().at(ig).Pt() != eventbase->GetTruth().at(ig).Pt()) continue;
	if(eventbase->GetTruth().at(ig).Pt() < 0.1) continue;
	if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
	if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
	cout << ig << " " << eventbase->GetTruth().at(ig).PdgId() << " : mother " <<  eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId() << " :  " << eventbase->GetTruth().at(ig).IndexMother() << endl;
	cout << eventbase->GetTruth().at(ig).Eta() << " " << eventbase->GetTruth().at(ig).Phi() << " " << eventbase->GetTruth().at(ig).Pt() <<endl;
      }
    }
    }*/

  std::vector<snu::KMuon> muonColl_truth = GetMuons(BaseSelection::MUON_POG_TIGHT);
  for(unsigned int im=0; im < muonColl_truth.size() ; im++){
  }
  /*cout << "mc pdgid = " <<muonColl_truth.at(im).MCMatchedPdgId() << endl;
    cout << "mc index = " << muonColl_truth.at(im).MCTruthIndex() << endl;
    cout << "mother pdgid = " << muonColl_truth.at(im).MotherPdgId()<< endl;
    cout << "motherindex= " <<muonColl_truth.at(im).MotherTruthIndex()<< endl;
    cout << "Is from conversion = " << muonColl_truth.at(im).MCIsFromConversion()<< endl;
  }
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).Pt() != eventbase->GetTruth().at(ig).Pt()) continue;
    if(eventbase->GetTruth().at(ig).Pt() < 0.1) continue;
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    cout << ig << " " << eventbase->GetTruth().at(ig).PdgId() << " : mother " <<  eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId() << " :  " << eventbase->GetTruth().at(ig).IndexMother() << endl;
    cout << eventbase->GetTruth().at(ig).Eta() << " " << eventbase->GetTruth().at(ig).Phi() << " " << eventbase->GetTruth().at(ig).Pt() <<endl;
  }
  */

  /// Apply MC weight for MCatnlo samples
  // MC weight = gen weight * lumimask weight
  // gen wieght = 1, -1
  // lumimask weight = 1 for silver json and lumi_gold/lumi_silver for gold

  if(!isData)weight*= MCweight;
  
  /// Apply json file if gold json is used. if lumimask == silver this does nothing  
  if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;
  
  /// FillCutFlow(cut, weight) fills a basic TH1 called cutflow. It is used to check number of events passing different cuts
  /// The string cut must match a bin label in FillCutFlow function
  FillCutFlow("NoCut", weight);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);

  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if(!PassBasicEventCuts()) return;     /// Initial event cuts  
  FillCutFlow("EventCut", weight);
  
  TString analysis_trigger="HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);
  
  std::vector<TString> triggerslist_23; // (prescaled)
  triggerslist_23.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");

  // This trigger will be  Unprescaled till L = 7E33 (may be used for this whole year) 
  // https://indico.cern.ch/event/370510/contribution/1/attachments/1161160/1671811/EleTriggers_Arun_28Sept_v1.pdf
  //if(!PassTrigger(triggerslist, prescale)) return;
  
  /// Target lumi = total lumi in json file. 
  /// ApplyPrescale reweights the MC to the luminosity of the trigger you are using


  if(PassTrigger(triggerslist, prescale)){
    FillCutFlow("TriggerCut", weight);
  }
  
  /// trigger_weight is for MC only: retruns 1 if data.
  /// Checks the luminosity of the trigger and returns weight that applied to 'weight' will correct for difference in luinosity of json file used in data
  float trigger_ps_weight= ApplyPrescale(analysis_trigger, TargetLumi,lumimask);

  FillHist("PSWeight" , trigger_ps_weight, 1., 0. , 2., 200);


  /// (0.998 is SF for trigger) https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgHLTScaleFactorMeasurements
  
  if(PassTrigger(triggerslist, prescale)){
    FillCutFlow("TriggerWeight", weight);
  }

  float weight_trigger_23=1.;
  if(!isData){
    weight_trigger_23 =  ApplyPrescale("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", TargetLumi,lumimask);
  }
  
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  
  //ListTriggersAvailable(); // uncomment this line to list off available triggers in the same

  // Trigger matching is done using KElectron::TriggerMatched(TString) which returns a bool

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("VertexCut", weight);


  /// Use the number of vertices in the event to check effect of pileup reweighting
  numberVertices = eventbase->GetEvent().nVertices();   
  
  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////

  
  // Get loose muons for veto: Can call  POGSoft/POGLoose/POGMedium/POGTight/HNVeto/HNLoose/HNMedium/HNTight
  std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_NOCUT); // loose selection
  
  /// Get tight jets : Can call NoLeptonVeto/Loose/Medium/Tight/HNJets
  std::vector<snu::KJet> jetColl_hn  = GetJets(BaseSelection::JET_HN);// pt > 20 ; eta < 2.5; PFlep veto; NO pileup ID
  std::vector<snu::KJet> jetColl_nlv  = GetJets(BaseSelection::JET_NOLEPTONVETO);
  std::vector<snu::KJet> jetColl_loose  = GetJets(BaseSelection::JET_LOOSE);

  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  
  cout << " " << endl;
  
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;

    if(eventbase->GetTruth().at(ig).PdgId() == 90 && eventbase->GetTruth().at(ig).GenStatus() == 22)  FillHist("N_mass", eventbase->GetTruth().at(ig).M() ,weight, 0. , 2500., 2500);
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 24 && eventbase->GetTruth().at(ig).GenStatus()== 22){
      cout << eventbase->GetTruth().at(ig).PdgId() << " " << eventbase->GetTruth().at(ig).GenStatus() << " " << eventbase->GetTruth().at(ig).M() << " " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< endl;
      if(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId() == 90)  FillHist("W2_mass", eventbase->GetTruth().at(ig).M() ,weight, 0. , 1500., 500);
      else if(eventbase->GetTruth().at(ig).M() == eventbase->GetTruth().at(ig).M())   FillHist("W1_mass", eventbase->GetTruth().at(ig).M() ,weight, 0. , 5000., 500);
    }
  }
  
      // Get POG electrons :  Can call POGVeto/POGLoose/POGMedium/POGTight/HNVeto/HNLoose/HNMedium/HNTight                                                                                              
  std::vector<snu::KElectron> electronLooseColl        = GetElectrons(BaseSelection::ELECTRON_POG_LOOSE);

  std::vector<snu::KElectron> electronColl_nocut             =  GetElectrons(BaseSelection::ELECTRON_NOCUT);

  //for(unsigned int iel = 0 ; iel < electronColl_nocut.size() ; iel++){
  // cout << "RECO " << electronColl_nocut.at(iel).Eta() << " " << electronColl_nocut.at(iel).Phi() << " " << electronColl_nocut.at(iel).Pt() << endl;
  // }

  
  std::vector<snu::KElectron> electronColl             = GetElectrons(false, false, BaseSelection::ELECTRON_POG_TIGHT);
  std::vector<snu::KElectron> electronColl_all             = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);
  
  FillHist("TruthMatchingAll", weight, electronColl_all.size(), 0., 6.,6);
  FillHist("TruthMatching", weight, electronColl.size(), 0., 6.,6);
 
  
  float weight_trigger_sf = TriggerScaleFactor(electronColl, muonColl, analysis_trigger);
  FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200);
  
  // Sets weight to weight if not running chargeflip bkg estimate or events are S
  //if(k_running_chargeflip) weight              *= WeightCFEvent(electronColl, k_running_chargeflip);
  

  std::vector<snu::KElectron> electronHNLooseColl  = GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE);
  std::vector<snu::KElectron> electronHNVetoColl   = GetElectrons(BaseSelection::ELECTRON_HN_VETO);
  std::vector<snu::KElectron> electronHNTightColl   = GetElectrons(BaseSelection::ELECTRON_HN_TIGHT);
  
  FillHist("NJets_nlv" , jetColl_nlv.size(), weight, 0., 5., 5);
  FillHist("NJets_loose" , jetColl_loose.size(), weight, 0., 5., 5);

  FillHist("NElectrons_nocut", electronColl_nocut.size(), weight, 0., 5., 5);
  FillHist("NElectrons_hnloose" ,  electronHNLooseColl.size(), weight, 0., 5., 5);
  FillHist("NElectrons_hnveto" ,   electronHNVetoColl.size(), weight, 0., 5., 5);
  FillHist("NElectrons_hntight" ,  electronHNTightColl.size(), weight, 0., 5., 5);
  FillHist("NElectrons_pogloose" ,  electronLooseColl.size(), weight, 0., 5., 5);
  FillHist("NElectrons_pogtight" ,  electronColl.size(), weight, 0., 5., 5);

  
  int njet = jetColl_hn.size();
  FillHist("GenWeight_NJet" , njet*MCweight + MCweight*0.1, 1., -6. , 6., 12);
 

  /// Correct MC for pileup   
  
  float pileup_reweight (1.);
  if (!k_isdata) {
    /// use silver or gold
    pileup_reweight = eventbase->GetEvent().PileUpWeight(lumimask);
  }
  FillHist("PileupWeight" , pileup_reweight, 1.,  0. , 2., 200);
  
  float id_weight=1.;
  float reco_weight=1.;
  
  if(!isData){
    for(unsigned int iel=0; iel < electronColl.size(); iel++){
      id_weight*= ElectronScaleFactor(BaseSelection::ELECTRON_POG_TIGHT, electronColl);
      reco_weight *= ElectronRecoScaleFactor(electronColl);
      
    }
  }
  FillHist("IDWeight" ,  id_weight,1.,  0. , 2., 200);
  FillHist("RecoWeight" ,  reco_weight, 1., 0. , 2., 200);

  FillCLHist(sighist_ee, "SIGNAL", eventbase->GetEvent(), muonColl,electronColl_nocut,jetColl_hn, weight*pileup_reweight);
  if(SameCharge(electronColl_nocut))    FillCLHist(sighist_ee, "SS_SIGNAL", eventbase->GetEvent(), muonColl,electronColl_nocut,jetColl_hn, weight*pileup_reweight);
  else     FillCLHist(sighist_ee, "OS_SIGNAL", eventbase->GetEvent(), muonColl,electronColl_nocut,jetColl_hn, weight*pileup_reweight);
  

  
  if(PassTrigger(triggerslist_23, prescale)){
    if(electronColl.size() ==2) {
      if(electronColl.at(0).Pt() > 25. && electronColl.at(1).Pt() > 15. ){
	FillCLHist(sighist_ee, "DiElectron_HLT23", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*weight_trigger_23);
      }
    }
  }
  

  if(electronColl_all.size() == 2 ) {

    if(electronColl_all.at(0).Pt() > 20. && electronColl_all.at(1).Pt() > 15. ){
      if(PassTrigger(triggerslist, prescale)){

        FillCLHist(sighist_ee, "DiElectron_noTM", eventbase->GetEvent(), muonColl,electronColl_all,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
      }
    }
  }

      
  if(electronColl.size() == 2 ) {
    
   if(electronColl.at(0).Pt() > 20. && electronColl.at(1).Pt() > 15. ){
      
      FillHist("Njets_dilepton", jetColl_hn.size() ,weight, 0. , 5., 5);
      FillCutFlow("DiEl_tight", weight);
      
      /// Method of plotting single histogram
      FillHist("zpeak_ee_noPUrw", GetDiLepMass(electronColl), weight, 0., 200.,400);
      FillHist("zpeak_ee", GetDiLepMass(electronColl), weight*pileup_reweight, 0., 200.,400);
      
      /// Standard set of histograms for muons/jets/electrons.. with no corrections
      
      FillCLHist(sighist_ee, "DiElectronNoPRW", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
      FillCLHist(sighist_ee, "DiElectronNoTrigger", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight);

      if(PassTrigger(triggerslist, prescale)){
	
	FillCLHist(sighist_ee, "DiElectron", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf);
	FillCLHist(sighist_ee, "DiElectronID" , eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight	);
	FillCLHist(sighist_ee, "DiElectronIDRECO" , eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);

	if(GetDiLepMass(electronColl) < 120. && GetDiLepMass(electronColl)  > 60. ){
	  if(!SameCharge(electronColl)){
	    FillCLHist(sighist_ee, "DiElectron_Zpeak", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	    if(electronColl.at(0).GsfCtfScPixChargeConsistency() && electronColl.at(1).GsfCtfScPixChargeConsistency()){
	      FillCLHist(sighist_ee, "DiElectron_Zpeak_chargeconsistency", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	      
	      
	      if(electronColl.at(0).IsEEFiducial() && electronColl.at(1).IsEEFiducial()) 
		FillCLHist(sighist_ee, "DiElectron_EE", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	      
	      if(electronColl.at(0).IsEBFiducial() && electronColl.at(1).IsEBFiducial())
		FillCLHist(sighist_ee, "DiElectron_BB", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	      
	      if(electronColl.at(0).IsEBFiducial() && electronColl.at(1).IsEEFiducial())
		FillCLHist(sighist_ee, "DiElectron_EB", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	      if(electronColl.at(1).IsEBFiducial() && electronColl.at(0).IsEEFiducial())
		FillCLHist(sighist_ee, "DiElectron_EB", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight*pileup_reweight*trigger_ps_weight*weight_trigger_sf*id_weight*reco_weight);
	    }
	  }
	}
      }
    }
  }
  
  
    
  
  if(electronLooseColl.size() == 3 && muonColl.size() == 0) {
    if(electronLooseColl.at(0).Pt() > 25. && electronLooseColl.at(2).Pt() > 25. ){
      FillHist("MET", eventbase->GetEvent().PFMET() ,weight, 0. , 100., 20);
      if( NBJet(jetColl_hn) == 0)FillCLHist(trilephist, "TriElectron_nomet", eventbase->GetEvent(), muonColl,electronLooseColl,jetColl_hn, weight*pileup_reweight);
      if(eventbase->GetEvent().PFMET() > 30){
	FillCLHist(trilephist, "TriElectron", eventbase->GetEvent(), muonColl,electronLooseColl,jetColl_hn, weight*pileup_reweight);
	if( NBJet(jetColl_hn) == 0)       FillCLHist(trilephist, "TriElectron_noB", eventbase->GetEvent(), muonColl,electronLooseColl,jetColl_hn, weight*pileup_reweight);
      }
    }
  }
    

  return;
}// End of execute event loop
  


void ExampleAnalyzerDiElectron::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

}


void ExampleAnalyzerDiElectron::BeginCycle() throw( LQError ){
  
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

ExampleAnalyzerDiElectron::~ExampleAnalyzerDiElectron() {
  
  Message("In ExampleAnalyzerDiElectron Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void ExampleAnalyzerDiElectron::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 10,0.,10.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerWeight");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"DiEl_tight");
   
    
  }
}


void ExampleAnalyzerDiElectron::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void ExampleAnalyzerDiElectron::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ExampleAnalyzerDiElectronCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ExampleAnalyzerDiElectron::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



