#ifndef trilepton_mumumu_CR_h
#define trilepton_mumumu_CR_h

#include "AnalyzerCore.h"
#include "Trilepton.h"

class trilepton_mumumu_CR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu_CR();
  ~trilepton_mumumu_CR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void SetBinInfo(int cut);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu_CR, 1);
};
#endif
