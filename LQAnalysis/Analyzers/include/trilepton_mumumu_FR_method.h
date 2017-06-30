#ifndef trilepton_mumumu_FR_method_h
#define trilepton_mumumu_FR_method_h

#include "AnalyzerCore.h"

class trilepton_mumumu_FR_method : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu_FR_method();
  ~trilepton_mumumu_FR_method();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  double this_dXYSig, this_RelIso;

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu_FR_method, 1);
};
#endif
