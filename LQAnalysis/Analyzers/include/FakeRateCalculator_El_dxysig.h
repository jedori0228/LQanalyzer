#ifndef FakeRateCalculator_El_dxysig_h
#define FakeRateCalculator_El_dxysig_h

#include "AnalyzerCore.h"


class FakeRateCalculator_El_dxysig : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_El_dxysig();
  ~FakeRateCalculator_El_dxysig();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);


  float GetPrescale(std::vector<snu::KElectron> electron, bool passlow, bool passhigh);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KElectron> out_electrons;
  TH2D *FR_sampleA;


  ClassDef ( FakeRateCalculator_El_dxysig, 1);
};
#endif
