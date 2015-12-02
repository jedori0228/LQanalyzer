#ifndef trilepton_mumue_h
#define trilepton_mumue_h

#include "AnalyzerCore.h"
#include <TNtupleD.h>

class trilepton_mumue : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumue();
  ~trilepton_mumue();
  
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
	TNtupleD *test;

  ClassDef (trilepton_mumue, 1);
};
#endif
