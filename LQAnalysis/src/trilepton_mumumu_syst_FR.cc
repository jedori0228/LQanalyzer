// $Id: trilepton_mumumu_syst_FR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQtrilepton_mumumu_syst_FR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "trilepton_mumumu_syst_FR.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_syst_FR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
trilepton_mumumu_syst_FR::trilepton_mumumu_syst_FR() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("trilepton_mumumu_syst_FR");
  
  Message("In trilepton_mumumu_syst_FR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void trilepton_mumumu_syst_FR::InitialiseAnalysis() throw( LQError ) {
 
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  TString lqdir = getenv("JSKIMROOTFILES");
  cout << lqdir+"/FRs.root" << endl;

  TFile *file_FRs = new TFile(lqdir+"/FRs.root");

  cout << "get scan values" << endl;
  TH1D *hist_dXYMins = (TH1D*)file_FRs->Get("hist_dXYMins");
  TH1D *hist_RelIsoMaxs = (TH1D*)file_FRs->Get("hist_RelIsoMaxs");
  for(int i=1; i<=hist_dXYMins->GetXaxis()->GetNbins(); i++) dXYMins.push_back( hist_dXYMins->GetBinContent(i) );
  for(int i=1; i<=hist_RelIsoMaxs->GetXaxis()->GetNbins(); i++) RelIsoMaxs.push_back( hist_RelIsoMaxs->GetBinContent(i) );

  cout << "get FRs" << endl;
  for(int aaa=0; aaa<dXYMins.size(); aaa++){
    for(int bbb=0; bbb<RelIsoMaxs.size(); bbb++){

      int dXY_Digit1 = int(dXYMins.at(aaa));
      int dXY_Digit0p1 = 10*dXYMins.at(aaa)-10*dXY_Digit1;
      TString str_dXYCut = "dXYSigMin_"+TString::Itoa(dXY_Digit1,10)+"p"+TString::Itoa(dXY_Digit0p1,10);

      int iso_Digit1 = int(RelIsoMaxs.at(bbb));
      int iso_Digit0p1 = 10*RelIsoMaxs.at(bbb)-10*iso_Digit1;
      TString str_iso = "LooseRelIsoMax_"+TString::Itoa(iso_Digit1,10)+"p"+TString::Itoa(iso_Digit0p1,10);

      str_dXYCut = str_dXYCut+"_"+str_iso;
      cout << "  "<<str_dXYCut<< endl;

      hist_trimuon_FR[str_dXYCut] = (TH2D*)file_FRs->Get(str_dXYCut+"_FR");
      hist_trimuon_FR_QCD[str_dXYCut] = (TH2D*)file_FRs->Get(str_dXYCut+"_FR_QCD");
      hist_trimuon_FRSF_QCD[str_dXYCut] = (TH2D*)file_FRs->Get(str_dXYCut+"_FRSF_QCD");
      //==== multiply SF
      hist_trimuon_FR_QCDSFed[str_dXYCut] = (TH2D*)hist_trimuon_FR[str_dXYCut]->Clone();
      hist_trimuon_FR_QCDSFed[str_dXYCut]->Multiply( hist_trimuon_FRSF_QCD[str_dXYCut] );
    }
  
  }


  cout << "get hist bins" << endl;
  TH1I* hist_bins = (TH1I*)file_FRs->Get("hist_bins");
  FR_n_pt_bin = hist_bins->GetBinContent(1);
  FR_n_eta_bin = hist_bins->GetBinContent(2);
  delete hist_bins;
  file_FRs->Close();
  delete file_FRs;

  return;

}


void trilepton_mumumu_syst_FR::ExecuteEvents()throw( LQError ){

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
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  //triggerlist.push_back("HLT_TripleMu_12_10_5_v");
  float trigger_ps_weight = WeightByTrigger(triggerlist, TargetLumi);
  bool trigger_pass = false;
  for(unsigned int i=0; i<triggerlist.size(); i++){
    if( PassTrigger(triggerlist.at(i)) ){
      trigger_pass = true;
      break;
    }
  }

  if(!DoCutOp){
    if(!trigger_pass) return;
    FillCutFlow("TriggerCut", 1.);
    m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;
  }


  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  /// Has Good Primary vertex:
  /// if ( vtx.ndof() > 4 &&
  //   ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= 24 ) &&
  //( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= 2 ) &&
  //!(vtx.isFake() ) ){
  FillCutFlow("VertexCut", 1.);


  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons("MUON_HN_TRI_VLOOSE");
 
 // CorrectMuonMomentum(muonTightColl);
  //float muon_id_iso_sf= MuonScaleFactor(BaseSelection::MUON_POG_TIGHT, muonTightColl,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT

  /// List of preset jet collections : NoLeptonVeto/Loose/Medium/Tight/TightLepVeto/HNJets

  /// can call POGVeto/POGLoose/POGMedium/POGTight/ HNVeto/HNLoose/HNTight/NoCut/NoCutPtEta 

  //float weight_trigger_sf = TriggerScaleFactor(electronColl, muonTightColl, "HLT_IsoMu20");

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    // check if catversion is empty. i.ie, v-7-4-X in which case use reweight class to get weight. In v-7-6-X+ pileupweight is stored in KEvent class, for silver/gold json
    //pileup_reweight = eventbase->GetEvent().PileUpWeight();
    pileup_reweight = TempPileupWeight();
  }

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();

  if(!isData && !k_running_nonprompt){
    //weight*=muon_id_iso_sf;
    //weight*=weight_trigger_sf;
    weight*=pileup_reweight;
    weight*=trigger_ps_weight;
    if(DoMCClosure){
      weight = 1.*MCweight;
    }
  }

  for(int bbb=0; bbb<RelIsoMaxs.size(); bbb++){

    std::vector<snu::KMuon> muontriLooseColl;
    muontriLooseColl.clear();
    for(unsigned int j=0; j<muontriLooseColl_raw.size(); j++){
      snu::KMuon this_muon = muontriLooseColl_raw.at(j);
      if( this_muon.RelIso04() < RelIsoMaxs.at(bbb) ) muontriLooseColl.push_back( this_muon );
    }

    for(int aaa=0; aaa<dXYMins.size(); aaa++){

      //==== prepare string

      int dXY_Digit1 = int(dXYMins.at(aaa));
      int dXY_Digit0p1 = 10*dXYMins.at(aaa)-10*dXY_Digit1;
      TString str_dXYCut = "dXYSigMin_"+TString::Itoa(dXY_Digit1,10)+"p"+TString::Itoa(dXY_Digit0p1,10);

      int iso_Digit1 = int(RelIsoMaxs.at(bbb));
      int iso_Digit0p1 = 10*RelIsoMaxs.at(bbb)-10*iso_Digit1;
      TString str_iso = "LooseRelIsoMax_"+TString::Itoa(iso_Digit1,10)+"p"+TString::Itoa(iso_Digit0p1,10);

      str_dXYCut = str_dXYCut+"_"+str_iso;

      int n_triTight_muons(0);
      for(unsigned int i=0; i<muontriLooseColl.size(); i++){
        if(muontriLooseColl.at(i).RelIso04() < 0.1) n_triTight_muons++;
      }
      int n_triLoose_muons = muontriLooseColl.size();

      bool isTwoMuon   = n_triLoose_muons == 2 && n_triTight_muons != 2; // no TT case
      bool isThreeMuon = n_triLoose_muons == 3 && n_triTight_muons != 3; // no TTT case

      //==== MCClosure

      if(isTwoMuon){

        snu::KMuon lep[2];
        lep[0] = muontriLooseColl.at(0);
        lep[1] = muontriLooseColl.at(1);

        double this_weight = weight;

        bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;
        bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();
        double m_dimuon = ( muontriLooseColl.at(0) + muontriLooseColl.at(1) ).M();

        std::map< TString, bool > map_whichCR_to_isCR;
        map_whichCR_to_isCR.clear();
        map_whichCR_to_isCR["DiMuon"] = isTwoMuon && leadPt20;
        map_whichCR_to_isCR["SSDiMuon"] = isTwoMuon && leadPt20 && isSS;

        vector<double> FR_muon, FR_error_muon;
        FR_muon.clear();
        FR_error_muon.clear();
        for(int i=0;i<2;i++){
          snu::KMuon this_muon = muontriLooseColl.at(i);
          //==== find loose but not tight muon ( 0.1 < RelIso (< 0.6) )
          if( this_muon.RelIso04() > 0.1 ){
            FR_muon.push_back( get_FR(this_muon, str_dXYCut, false) );
            FR_error_muon.push_back( get_FR(this_muon, str_dXYCut, true) );
          }
        }
        for(unsigned int i=0; i<FR_muon.size(); i++){
          this_weight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
        }
        if( FR_muon.size() == 2 ) this_weight *= -1.; // minus sign for LL
        //==== weight error
        double this_weight_err(0.);
        if( FR_muon.size() == 1 ){
          double fr1 = FR_muon.at(0);
          double fr1_err = FR_error_muon.at(0);
          this_weight_err = fr1_err/pow(fr1-1,2);
        }
        else if( FR_muon.size() == 2 ){
          double fr1 = FR_muon.at(0);
          double fr1_err = FR_error_muon.at(0);
          double fr2 = FR_muon.at(1);
          double fr2_err = FR_error_muon.at(1);
          this_weight_err = sqrt( pow( fr1_err*fr2*(1-fr2),2) +
                             pow( fr2_err*fr1*(1-fr1),2)   ) / pow( (1-fr1)*(1-fr2), 2 );
        }
        else{
          Message("?", INFO);
        }
        for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
          TString this_suffix = it->first;
          if(it->second){
            FillUpDownHist(str_dXYCut+"_n_events_"+this_suffix, 0., this_weight, this_weight_err, 0., 1., 1);
            FillHist(str_dXYCut+"_n_OnlyLoose_"+this_suffix, FR_muon.size(), 1., 0., 4., 4);
          }
        }

      } // is TwoMuon


      //==== Preselection

      if(isThreeMuon){
     
        bool leadPt20 = muontriLooseColl.at(0).Pt() > 20.;

        snu::KParticle lep[3];
        vector<double> FR_muon, FR_error_muon;
        FR_muon.clear();
        FR_error_muon.clear();
        //for(int i=0;i<k_flags.size();i++) cout << "k_flags = " << k_flags.at(i) << endl;
        for(int i=0;i<3;i++){
          lep[i] = muontriLooseColl.at(i);
          //==== find loose but not tight muon ( 0.1 < RelIso (< 0.6) )
          if( muontriLooseColl.at(i).RelIso04() > 0.1 ){
            FR_muon.push_back( get_FR(lep[i], str_dXYCut, false) );
            FR_error_muon.push_back( get_FR(lep[i], str_dXYCut, true) );
          }
        }

        double this_weight = weight;

        for(unsigned int i=0; i<FR_muon.size(); i++){
          this_weight *= FR_muon.at(i)/( 1.-FR_muon.at(i) );
        }
        if( FR_muon.size() == 2 ) this_weight *= -1.; // minus sign for TLL
        //==== weight error
        double this_weight_err(0.);
        if( FR_muon.size() == 1 ){
          double fr1 = FR_muon.at(0);
          double fr1_err = FR_error_muon.at(0);
          this_weight_err = fr1_err/pow(fr1-1,2);
        }
        else if( FR_muon.size() == 2 ){
          double fr1 = FR_muon.at(0);
          double fr1_err = FR_error_muon.at(0);
          double fr2 = FR_muon.at(1);
          double fr2_err = FR_error_muon.at(1);
          this_weight_err = sqrt( pow( fr1_err*fr2*(1-fr2),2) +
                             pow( fr2_err*fr1*(1-fr1),2)   ) / pow( (1-fr1)*(1-fr2), 2 );
        }
        else if( FR_muon.size() == 3 ){
          double fr1 = FR_muon.at(0);
          double fr1_err = FR_error_muon.at(0);
          double fr2 = FR_muon.at(1);
          double fr2_err = FR_error_muon.at(1);
          double fr3 = FR_muon.at(2);
          double fr3_err = FR_error_muon.at(2);
          this_weight_err = sqrt( pow( fr1_err*fr2*(1-fr2)*fr3*(1-fr3), 2) +
                             pow( fr2_err*fr3*(1-fr3)*fr1*(1-fr1), 2) +
                             pow( fr3_err*fr1*(1-fr1)*fr2*(1-fr2), 2)   ) / pow( (1-fr1)*(1-fr2)*(1-fr3) ,2 );
        }
        else{
          Message("?", INFO);
        }

        bool isOS = true;

        int OppSign, SameSign[2]; // SameSign[0].Pt() > SameSign[1].Pt()
        if(lep[0].Charge() * lep[1].Charge() > 0){ // Q(0) = Q(1)
          if(lep[1].Charge() * lep[2].Charge() < 0){ // Q(1) != Q(2)
            OppSign = 2;
            SameSign[0] = 0;
            SameSign[1] = 1;
          }
          else isOS = false; // veto Q(0) = Q(1) = Q(2)
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

        // MC samples has m(OS)_saveflavour > 4 GeV cut at gen level
        // MADGRAPH : https://github.com/cms-sw/genproductions/blob/master/bin/MadGraph5_aMCatNLO/cards/production/13TeV/WZTo3LNu01j_5f_NLO_FXFX/WZTo3LNu01j_5f_NLO_FXFX_run_card.dat#L130
        // POWHEG   : https://github.com/cms-sw/genproductions/blob/master/bin/Powheg/production/WZTo3lNu_NNPDF30_13TeV/WZ_lllnu_NNPDF30_13TeV.input#L2
        bool NOTmll4 = true;
        if(isOS){
          if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
              (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) NOTmll4 = false;
        }

        if( leadPt20 && isOS && NOTmll4 ){
          FillUpDownHist(str_dXYCut+"_n_events_"+"Preselection", 0., this_weight, this_weight_err, 0., 1., 1);
          FillHist(str_dXYCut+"_n_OnlyLoose_"+"Preselection", FR_muon.size(), 1., 0., 4., 4);
        }

      }





    }
  }



  return;

}// End of execute event loop
  


void trilepton_mumumu_syst_FR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void trilepton_mumumu_syst_FR::BeginCycle() throw( LQError ){
  
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

trilepton_mumumu_syst_FR::~trilepton_mumumu_syst_FR() {
  
  Message("In trilepton_mumumu_syst_FR Destructor" , INFO);
  
}


void trilepton_mumumu_syst_FR::FillCutFlow(TString cut, float weight){

  
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


void trilepton_mumumu_syst_FR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void trilepton_mumumu_syst_FR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this trilepton_mumumu_syst_FRCore::MakeHistograms() to make new hists for your analysis
   **/

}


void trilepton_mumumu_syst_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

double trilepton_mumumu_syst_FR::get_FR(snu::KParticle muon, TString whichFR, bool geterror){

  double this_pt = muon.Pt();
  double this_eta = fabs( muon.Eta() );

  // FR_n_pt_bin = 7
  // array index      0    1    2    3    4    5    6    7
  // bin numbe          1    2     3    4    5   6     7
  // ptarray[7+1] = {10., 15., 20., 25., 30., 35., 45., 60.}; 

  double ptarray[FR_n_pt_bin+1], etaarray[FR_n_eta_bin+1];
  //cout << "FR_n_pt_bin = " << FR_n_pt_bin << endl;
  for(int i=0; i<FR_n_pt_bin; i++){
    ptarray[i] = hist_trimuon_FR[whichFR]->GetXaxis()->GetBinLowEdge(i+1);
    //cout << " " << ptarray[i] << endl;
    if(i==FR_n_pt_bin-1){
      ptarray[FR_n_pt_bin] = hist_trimuon_FR[whichFR]->GetXaxis()->GetBinUpEdge(i+1);
      //cout << " " << ptarray[FR_n_pt_bin] << endl;
    }
  }
  //cout << "FR_n_eta_bin = " << FR_n_eta_bin << endl;
  for(int i=0; i<FR_n_eta_bin; i++){
    etaarray[i] = hist_trimuon_FR[whichFR]->GetYaxis()->GetBinLowEdge(i+1);
    //cout << " " << etaarray[i] << endl;
    if(i==FR_n_eta_bin-1){
      etaarray[FR_n_eta_bin] = hist_trimuon_FR[whichFR]->GetYaxis()->GetBinUpEdge(i+1);
      //cout << " " << etaarray[FR_n_eta_bin] << endl;
    }
  }

  int this_pt_bin;
  if( this_pt >= ptarray[FR_n_pt_bin] ) this_pt_bin = FR_n_pt_bin;
  else{
    for(int i=0; i<FR_n_pt_bin; i++){
      if( ptarray[i] <= this_pt && this_pt < ptarray[i+1] ){
        this_pt_bin = i+1;
        break;
      }
    }
  }
  //cout << "this pt bin = " << this_pt_bin << endl;
  int this_eta_bin;
  if( this_eta >= etaarray[FR_n_eta_bin] ) this_eta_bin = FR_n_eta_bin;
  else{
    for(int i=0; i<FR_n_eta_bin; i++){
      if( etaarray[i] <= this_eta && this_eta < etaarray[i+1] ){
        this_eta_bin = i+1;
        break;
      }
    }
  }
  //cout << "this eta bin = " << this_eta_bin << endl;

  //==== FR
  double this_FR = hist_trimuon_FR[whichFR]->GetBinContent(this_pt_bin, this_eta_bin);
  double FRSF_QCD = hist_trimuon_FRSF_QCD[whichFR]->GetBinContent(this_pt_bin, this_eta_bin);
  double this_FR_QCDSFed = hist_trimuon_FR_QCDSFed[whichFR]->GetBinContent(this_pt_bin, this_eta_bin);
  //==== error
  double this_FR_error = hist_trimuon_FR[whichFR]->GetBinError(this_pt_bin, this_eta_bin);
  double this_FR_QCDSFed_error = hist_trimuon_FR_QCDSFed[whichFR]->GetBinError(this_pt_bin, this_eta_bin);

  //==== FR QCD and error
  double this_FR_QCD = hist_trimuon_FR_QCD[whichFR]->GetBinContent(this_pt_bin, this_eta_bin);
  double this_FR_QCD_error = hist_trimuon_FR_QCD[whichFR]->GetBinError(this_pt_bin, this_eta_bin);

  //==== bool for MCClosure
  bool doMCClosure = std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();

  //cout << "==== FR ====" << endl;
  //cout << "this_FR = " << this_FR << endl;
  //cout << "FRSF_QCD = " << FRSF_QCD << endl;
  //cout << "==> multiply by hand gives = " << this_FR*FRSF_QCD << endl;
  //cout << "==> this_FR_QCDSFed = " << this_FR_QCDSFed << endl;
  //cout << "==== FR Error ====" << endl;
  //cout << "this_FR_error = " << this_FR_error << endl;
  //cout << "==> multiply by hand gives = " << this_FR_error*FRSF_QCD << endl;
  //cout << "==> this_FR_QCDSFed_error = " << this_FR_QCDSFed_error << endl << endl;

  if(geterror){
    if(doMCClosure) return this_FR_QCD_error;
    return this_FR_QCDSFed_error;
  }
  else{
    if(doMCClosure) return this_FR_QCD;
    return this_FR_QCDSFed;
  }

}



