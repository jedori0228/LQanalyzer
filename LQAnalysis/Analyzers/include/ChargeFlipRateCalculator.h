#ifndef ChargeFlipRateCalculator_h
#define ChargeFlipRateCalculator_h

#include "AnalyzerCore.h"

class ChargeFlipRateCalculator : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ChargeFlipRateCalculator();
  ~ChargeFlipRateCalculator();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void FillDenAndNum(TString prefix, snu::KElectron el, double thisweight);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ChargeFlipRateCalculator, 1);
};
#endif
