#ifndef trilepton_eemu_h
#define trilepton_eemu_h

#include "AnalyzerCore.h"
#include <TNtupleD.h>

class trilepton_eemu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_eemu();
  ~trilepton_eemu();
  
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

  ClassDef (trilepton_eemu, 1);
};
#endif
