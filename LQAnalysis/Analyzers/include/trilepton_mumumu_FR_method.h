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
  bool PassOptimizedCut(int sig_mass, double first_pt, double second_pt, double third_pt, double W_pri_mass, double PFMET);

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
