#ifndef AnalyzerSkeleton_h
#define AnalyzerSkeleton_h

#include "AnalyzerCore.h"

class AnalyzerSkeleton : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  AnalyzerSkeleton();
  ~AnalyzerSkeleton();

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

  ClassDef ( AnalyzerSkeleton, 1);
};
#endif
