#ifndef trilepton_mumumu_FR_method_h
#define trilepton_mumumu_FR_method_h

#include "AnalyzerCore.h"
#include "Trilepton.h"

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

  TH2F* hist_trimuon_FR[2];
  int FR_n_pt_bin, FR_n_eta_bin;
  double get_FR(snu::KParticle muon, int n_jets);

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
