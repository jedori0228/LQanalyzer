// $Id: AnalyzerSkeleton.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAnalyzerSkeleton Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "AnalyzerSkeleton.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (AnalyzerSkeleton);

AnalyzerSkeleton::AnalyzerSkeleton() :  AnalyzerCore()
{
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("AnalyzerSkeleton");
  
  Message("In AnalyzerSkeleton constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void AnalyzerSkeleton::InitialiseAnalysis() throw( LQError ) {
  
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


void AnalyzerSkeleton::ExecuteEvents()throw( LQError ){

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
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  if(!PassTriggerOR(triggerlist)) return;
  FillCutFlow("TriggerCut", 1.);
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);

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


  return;
} // End of execute event loop
  


void AnalyzerSkeleton::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void AnalyzerSkeleton::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  return;
  
}

AnalyzerSkeleton::~AnalyzerSkeleton() {
  
  Message("In AnalyzerSkeleton Destructor" , INFO);
  
}


void AnalyzerSkeleton::FillCutFlow(TString cut, float w){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,w);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 12,0.,12.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    
  }
}


void AnalyzerSkeleton::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void AnalyzerSkeleton::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this AnalyzerSkeletonCore::MakeHistograms() to make new hists for your analysis
   **/

}


void AnalyzerSkeleton::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //


}

