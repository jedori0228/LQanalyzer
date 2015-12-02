#ifndef trilepton_eee_h
#define trilepton_eee_h

#include "AnalyzerCore.h"
#include <TMinuit.h>
#include <TNtupleD.h>

class trilepton_eee : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_eee();
  ~trilepton_eee();

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

  void mysort(std::vector<snu::KJet> jetColl,int jidx[]);
  TNtupleD *test;

  ClassDef ( trilepton_eee, 1);
};
#endif
