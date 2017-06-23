// $Id: ChargeFlipRateCalculator.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQChargeFlipRateCalculator Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ChargeFlipRateCalculator.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ChargeFlipRateCalculator);

ChargeFlipRateCalculator::ChargeFlipRateCalculator() :  AnalyzerCore(), out_muons(0)
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ChargeFlipRateCalculator");
  
  Message("In ChargeFlipRateCalculator constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void ChargeFlipRateCalculator::InitialiseAnalysis() throw( LQError ) {
  
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


void ChargeFlipRateCalculator::ExecuteEvents()throw( LQError ){

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
  FillHist("cutflow_MuMuE", 0., 1., 0., 10., 10);
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  //======================
  //==== [CUT] METFilter
  //======================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", 1.);
  FillHist("cutflow_MuMuE", 1., 1., 0., 10., 10);

  //====================
  //==== [CUT] Trigger
  //====================
/*
  std::vector<TString> triggerlist;
  triggerlist.push_back("");

  if(!PassTriggerOR(triggerlist)) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
*/
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

  //======================
  //==== Pileup Reweight
  //======================

  float pileup_reweight=(1.0);
  if(!k_isdata){
    //==== CATTools reweight
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    //==== John reweight
    //pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
  }

  //====================
  //==== Get Electrons
  //====================

  std::vector<snu::KElectron> electrons_raw = GetElectrons(true, false, "ELECTRON_MVA_TIGHT");

  //===============================
  //==== Get Electron Corrections
  //===============================

  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electrons_raw, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons_raw);

  //========================
  //==== Apply corrections
  //========================

  if(!isData && !k_running_nonprompt){
    //weight*=trigger_ps_weight;
    weight*=pileup_reweight;
    weight*=electron_sf;
    weight*=electron_RecoSF;
  }

  //==== Charge Flip Rate Measurement with MC
  bool GetCFMC = std::find(k_flags.begin(), k_flags.end(), "GetCFMC") != k_flags.end();

  if(GetCFMC && !isData){

    float etaarray [] = {0.0, 0.9, 1.4442, 1.556, 2.5};
    float ptarray [] = {20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 450., 500.};

    for(unsigned int i=0; i<electrons_raw.size(); i++){
      snu::KElectron el = electrons_raw.at(i);

      bool isBarrel = fabs(el.SCEta())<1.479;
      FillDenAndNum("AllEtaRegion_", el, 1.);
      if(isBarrel){
        FillDenAndNum("Barrel_", el, 1.);
      }
      else{
        FillDenAndNum("Endcap_", el, 1.);
      }

    }

    //==== Z event
    double m_Z = 91.1876;
    if(electrons_raw.size()==2){
      if(electrons_raw.at(1).Pt() >= 20) {

        double MET = eventbase->GetEvent().MET();
        snu::KParticle Z_cand = electrons_raw.at(0)+electrons_raw.at(1);
        bool isZResonance = fabs(Z_cand.M() - m_Z) < 10.;

        if( MET < 30. && isZResonance ){
          for(int i=0; i<2; i++){
            snu::KElectron el = electrons_raw.at(i);

            bool isBarrel = fabs(el.SCEta())<1.479;
            FillDenAndNum("ZEvent_AllEtaRegion_", el, 1.);
            if(isBarrel){
              FillDenAndNum("ZEvent_Barrel_", el, 1.);
            }
            else{
              FillDenAndNum("ZEvent_Endcap_", el, 1.);
            }


          }
        }

      } // pt > 20 Gev
    } // two electrons

  }


  return;


}// End of execute event loop
  


void ChargeFlipRateCalculator::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ChargeFlipRateCalculator::BeginCycle() throw( LQError ){
  
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

ChargeFlipRateCalculator::~ChargeFlipRateCalculator() {
  
  Message("In ChargeFlipRateCalculator Destructor" , INFO);
  
}


void ChargeFlipRateCalculator::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12,0.,12.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    
  }
}


void ChargeFlipRateCalculator::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ChargeFlipRateCalculator::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ChargeFlipRateCalculatorCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ChargeFlipRateCalculator::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

void ChargeFlipRateCalculator::FillDenAndNum(TString prefix, snu::KElectron el, double thisweight){

  float etaarray [] = {0.0, 0.9, 1.4442, 1.556, 2.5};
  float ptarray [] = {20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 450., 500.};

  double invpt = 1./el.Pt();

  FillHist(prefix+"events_F0", fabs(el.Eta()), el.Pt(), thisweight, etaarray, 4, ptarray, 25);
  FillHist(prefix+"invpt_F0", invpt, thisweight, 0., 0.1, 100);
  FillHist(prefix+"pt_F0", el.Pt(), thisweight, 0., 500., 500);
  FillHist(prefix+"eta_F0", el.Eta(), thisweight, -3., 3., 60);
  FillHist(prefix+"mva_F0", el.MVA(), thisweight, -1., 1., 200);
  FillHist(prefix+"dXY_F0", fabs(el.dxy()), thisweight, 0., 0.1, 100);
  FillHist(prefix+"dZ_F0", fabs(el.dxy()), thisweight, 0., 0.2, 200);
  FillHist(prefix+"Type_F0", el.GetType(), thisweight, 0., 100., 100);

  if(MCIsCF(el)){
    FillHist(prefix+"events_F", fabs(el.Eta()), el.Pt(), thisweight, etaarray, 4, ptarray, 25);
    FillHist(prefix+"invpt_F", invpt, thisweight, 0., 0.1, 100);
    FillHist(prefix+"pt_F", el.Pt(), thisweight, 0., 500., 500);
    FillHist(prefix+"eta_F", el.Eta(), thisweight, -3., 3., 60);
    FillHist(prefix+"mva_F", el.MVA(), thisweight, -1., 1., 200);
    FillHist(prefix+"dXY_F", fabs(el.dxy()), thisweight, 0., 0.1, 100);
    FillHist(prefix+"dZ_F", fabs(el.dxy()), thisweight, 0., 0.2, 200);
    FillHist(prefix+"Type_F", el.GetType(), thisweight, 0., 100., 100);
  }


}









