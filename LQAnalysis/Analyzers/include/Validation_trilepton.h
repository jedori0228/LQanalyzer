#ifndef Validation_trilepton_h
#define Validation_trilepton_h

#include "AnalyzerCore.h"

class Validation_trilepton : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Validation_trilepton();
  ~Validation_trilepton();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( Validation_trilepton, 1);
};
#endif
