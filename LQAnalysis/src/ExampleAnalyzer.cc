// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzer.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ExampleAnalyzer);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
ExampleAnalyzer::ExampleAnalyzer() :  AnalyzerCore(), out_muons(0), out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("ExampleAnalyzer");

  Message("In ExampleAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
}


void ExampleAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist
   MakeCleverHistograms(sighist, "Zmuons_jlv");
   MakeCleverHistograms(sighist, "Zelectrons_jlv");
   MakeCleverHistograms(sighist, "Sigmuons");
   MakeCleverHistograms(sighist, "Sigelectrons");

   return;
 }


void ExampleAnalyzer::ExecuteEvents()throw( LQError ){

 //==== TEST 8 TeV DiMuon 
  std::vector<snu::KMuon> testmuons = GetMuons("tight");
  std::vector<snu::KMuon> testloosemuons = GetMuons("loose");
  std::vector<snu::KJet> testjets = GetJets("ApplyPileUpID");
  int n_bjets = NBJet(testjets);
  if(testmuons.size()>=1){
    FillHist("jets_size", testjets.size(), weight, 0., 10., 10);
    FillHist("jets_nbjet", n_bjets, weight, 0., 10., 10);
    int n_genpu=0;
    for(int i=0; i<testjets.size(); i++){
      FillHist("jets_PartonFlavour", testjets.at(i).PartonFlavour(), weight, -50, 50, 100);
      FillHist("jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
      FillHist("jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
      FillHist("jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
      if(testjets.at(i).PartonFlavour()==0){
        n_genpu++;
        FillHist("GENPU_jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
        FillHist("GENPU_jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
        FillHist("GENPU_jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
      }
    }
    FillHist("GENPU_jets_size", n_genpu, weight, 0., 10., 10);

    if(testloosemuons.size()==2){

      bool IsSS = (testloosemuons.at(0).Charge())==(testloosemuons.at(1).Charge());
      if(k_sample_name.Contains("Wjets")) IsSS = true;

      //==== SS
      if(IsSS && testloosemuons.at(0).Pt() >= 20. && testloosemuons.at(1).Pt() >= 15.){

        FillHist("TL_jets_size", testjets.size(), weight, 0., 10., 10);
        FillHist("TL_jets_nbjet", n_bjets, weight, 0., 10., 10);

        //==== SS+2jet
        if(testjets.size()>=2){
          FillHist("TL2jet_jets_size", testjets.size(), weight, 0., 10., 10);
          FillHist("TL2jet_jets_nbjet", n_bjets, weight, 0., 10., 10);
          for(int i=0; i<testjets.size(); i++){
            FillHist("TL2jet_jets_PartonFlavour", testjets.at(i).PartonFlavour(), weight, -50, 50, 100);
            FillHist("TL2jet_jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
            FillHist("TL2jet_jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
            FillHist("TL2jet_jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
          }

          //==== TT
          if(testmuons.size()==2){
            FillHist("TT2jet_jets_size", testjets.size(), weight, 0., 10., 10);
            FillHist("TT2jet_jets_nbjet", n_bjets, weight, 0., 10., 10);
            for(int i=0; i<testjets.size(); i++){
              FillHist("TT2jet_jets_PartonFlavour", testjets.at(i).PartonFlavour(), weight, -50, 50, 100);
              FillHist("TT2jet_jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
              FillHist("TT2jet_jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
              FillHist("TT2jet_jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
            }
          }

        }
        n_genpu=0;
        for(int i=0; i<testjets.size(); i++){
          FillHist("TL_jets_PartonFlavour", testjets.at(i).PartonFlavour(), weight, -50, 50, 100);
          FillHist("TL_jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
          FillHist("TL_jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
          FillHist("TL_jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
          if(testjets.at(i).PartonFlavour()==0){
            n_genpu++;
            FillHist("TL_GENPU_jets_pt", testjets.at(i).Pt(), weight, 0., 200., 200);
            FillHist("TL_GENPU_jets_eta", testjets.at(i).Eta(), weight, -3., 3., 60);
            FillHist("TL_GENPU_jets_PUMVA", testjets.at(i).PileupJetIDMVA(), weight, -1., 1., 200);
          }
        }
        FillHist("TL_GENPU_jets_size", n_genpu, weight, 0., 10., 10);

      }
    }

  }
  return;
  
}// End of execute event loop
  


void ExampleAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ExampleAnalyzer::BeginCycle() throw( LQError ){
  
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

ExampleAnalyzer::~ExampleAnalyzer() {
  
  Message("In ExampleAnalyzer Destructor" , INFO);
  if(!k_isdata)delete reweightPU;
  
}


void ExampleAnalyzer::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 5,0.,5.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"DiEl_tight");
   
    
  }
}


void ExampleAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void ExampleAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ExampleAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ExampleAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



