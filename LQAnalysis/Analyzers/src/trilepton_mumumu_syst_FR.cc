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
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (trilepton_mumumu_syst_FR);

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

  return;

}


void trilepton_mumumu_syst_FR::ExecuteEvents()throw( LQError ){

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

  //=============================================
  //==== Prepare Muons and Rochestor Correction
  //=============================================

  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons("MUON_HN_TRI_VLOOSE");

  //===============
  //==== Get Jets
  //===============

  std::vector<snu::KJet> jetColl_hn = GetJets("JET_HN", 30., 2.4);
  int n_jets = jetColl_hn.size();
  int n_bjets=0;
  for(int j=0; j<n_jets; j++){
    if( IsBTagged(jetColl_hn.at(j), snu::KJet::CSVv2, snu::KJet::Medium) ){
      n_bjets++;
      FillHist("bjet_pt", jetColl_hn.at(j).Pt(), 1., 0., 200., 200);
    }
  }

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrontriLooseColl = GetElectrons(false, false, "ELECTRON_MVA_FAKELOOSE", 15., 2.5);

  int n_triLoose_electrons = electrontriLooseColl.size();
  int n_triTight_electrons(0);
  for(unsigned int i=0; i<electrontriLooseColl.size(); i++){
    if(PassID(electrontriLooseColl.at(i), "ELECTRON_MVA_TIGHT")) n_triTight_electrons++;
  }

  //=======================================================
  //==== For MC Closure test, let's not normalise to Lumi
  //=======================================================

  bool DoMCClosure = std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();
  //==== No normalization for MC Closure
  if(DoMCClosure){
    weight = 1.*MCweight;
    //m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true); //FIXME this functino is removedi n v8-0-7.27!!
  }
  dXYMins.clear();
  RelIsoMaxs.clear();
  //dXYMins = m_datadriven_bkg->GetFakeObj()->GetdXYMins(); //FIXME this functino is removedi n v8-0-7.27!!
  //RelIsoMaxs = m_datadriven_bkg->GetFakeObj()->GetRelIsoMaxs(); //FIXME this functino is removedi n v8-0-7.27!!

  bool UsePtCone = std::find(k_flags.begin(), k_flags.end(), "UsePtCone") != k_flags.end();
  m_datadriven_bkg->SetUsePtCone(UsePtCone);

  //===============================================
  //==== Large dXYSig muon definition loop starts
  //===============================================

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
        if(PassID(muontriLooseColl.at(i), "MUON_HN_TRI_TIGHT")) n_triTight_muons++;
      }
      int n_triLoose_muons = muontriLooseColl.size();

      int n_triLoose_leptons = n_triLoose_muons+n_triLoose_electrons;
      int n_triTight_leptons = n_triTight_muons+n_triTight_electrons;

      bool isTwoMuon = (n_triLoose_leptons == 2)
                       && (n_triLoose_muons == 2 && n_triTight_muons != 2);
      bool isThreeMuon = (n_triLoose_leptons == 3)
                         && (n_triLoose_muons == 3 && n_triTight_muons != 3);

      //==== MCClosure

      if(isTwoMuon && DoMCClosure){

        snu::KMuon lep[2];
        lep[0] = muontriLooseColl.at(0);
        lep[1] = muontriLooseColl.at(1);

        if( muontriLooseColl.at(0).Pt() < 20. ) continue;
        bool isSS = muontriLooseColl.at(0).Charge() == muontriLooseColl.at(1).Charge();

        std::map< TString, bool > map_whichCR_to_isCR;
        map_whichCR_to_isCR.clear();
        map_whichCR_to_isCR["DiMuon"] = true;
        map_whichCR_to_isCR["SSDiMuon"] = isSS;

        //==== fake method weighting
        //m_datadriven_bkg->GetFakeObj()->SetTrilepWP(dXYMins.at(aaa), RelIsoMaxs.at(bbb)); //FIXME this functino is removedi n v8-0-7.27!!
        std::vector<snu::KElectron> empty_electron;
        empty_electron.clear();
        double this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muontriLooseColl, "MUON_HN_TRI_TIGHT", 2, empty_electron, "ELECTRON_MVA_TIGHT", 0);
        double this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muontriLooseColl, "MUON_HN_TRI_TIGHT", 2, empty_electron, "ELECTRON_MVA_TIGHT", 0);

        for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
          TString this_suffix = it->first;
          if(it->second){
            FillUpDownHist(str_dXYCut+"_n_events_"+this_suffix, 0., this_weight, this_weight_err, 0., 1., 1);
            //FillHist(str_dXYCut+"_n_OnlyLoose_"+this_suffix, m_datadriven_bkg->GetFakeObj()->GetNLooseNotTight(), 1., 0., 4., 4);  //FIXME this functino is removedi n v8-0-7.27!!
          }
        }

      } // is TwoMuon


      //==== Preselection

      if(isThreeMuon){
     
        if( muontriLooseColl.at(0).Pt() < 20. ) continue;

        std::vector<KLepton> lep;
        for(unsigned int i=0; i<muontriLooseColl.size(); i++){
          KLepton this_lep( muontriLooseColl.at(i) );
          lep.push_back( this_lep );
        }

        snu::KParticle HN[4];;

        //==== fake method weighting
        //m_datadriven_bkg->GetFakeObj()->SetTrilepWP(dXYMins.at(aaa), RelIsoMaxs.at(bbb)); //FIXME this functino is removedi n v8-0-7.27!!
        std::vector<snu::KElectron> empty_electron;
        empty_electron.clear();

        //==== fake method weighting
        double this_weight     = m_datadriven_bkg->Get_DataDrivenWeight(false, muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
        double this_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muontriLooseColl, "MUON_HN_TRI_TIGHT", muontriLooseColl.size(), electrontriLooseColl, "ELECTRON_MVA_TIGHT", electrontriLooseColl.size(), "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");

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

        bool VetoZResonance;
        bool mllloffZ;
        bool isLowMass;
        bool isHighMass;

        //==== mll4
        if( (lep[SameSign[0]]+lep[OppSign]).M() <= 4. ||
            (lep[SameSign[1]]+lep[OppSign]).M() <= 4.     ) continue;

        ///////////////////////////////////////////
        ////////// m(HN) < 80 GeV region //////////
        ///////////////////////////////////////////

        snu::KEvent Evt = eventbase->GetEvent();
        double MET = Evt.MET();
        double METphi = Evt.METPhi();
        CorrectedMETRochester(muontriLooseColl, MET, METphi);

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

        VetoZResonance = fabs(z_candidate.M()-91.1876) > 15.;
        mllloffZ = fabs( (lep[0] + lep[1] + lep[2]).M() - 91.1876 ) > 15.;

        isLowMass = (W_pri_lowmass.M() < 150.);
        isHighMass = (MET > 20.);

        std::map< TString, bool > map_whichCR_to_isCR;
        map_whichCR_to_isCR.clear();
        map_whichCR_to_isCR["BasicSelection"] = true;
        map_whichCR_to_isCR["ZVeto"]          = VetoZResonance;
        map_whichCR_to_isCR["ZVeto_mllloffZ"] = VetoZResonance && mllloffZ;
        map_whichCR_to_isCR["Preselection"]   = VetoZResonance && mllloffZ && n_bjets==0;
        map_whichCR_to_isCR["LowMass"]        = map_whichCR_to_isCR["Preselection"] && isLowMass;
        map_whichCR_to_isCR["HighMass"]       = map_whichCR_to_isCR["Preselection"] && isHighMass;

        for(std::map< TString, bool >::iterator it = map_whichCR_to_isCR.begin(); it != map_whichCR_to_isCR.end(); it++){
          TString this_suffix = it->first;
          if(it->second){
            FillUpDownHist(str_dXYCut+"_n_events_"+this_suffix, 0., this_weight, this_weight_err, 0., 1., 1);
            //FillHist(str_dXYCut+"_n_OnlyLoose_"+this_suffix, m_datadriven_bkg->GetFakeObj()->GetNLooseNotTight(), 1., 0., 4., 4); //FIXME this functino is removedi n v8-0-7.27!!
          }
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




