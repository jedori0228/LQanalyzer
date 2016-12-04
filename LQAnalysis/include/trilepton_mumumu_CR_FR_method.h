#ifndef trilepton_mumumu_CR_FR_method_h
#define trilepton_mumumu_CR_FR_method_h

#include "AnalyzerCore.h"

class trilepton_mumumu_CR_FR_method : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu_CR_FR_method();
  ~trilepton_mumumu_CR_FR_method();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  TH2F* hist_trimuon_FR[10];
  TH2F* hist_trimuon_FRSF;
  TH1F* hist_trimuon_FRSF_pt;
  int FR_n_pt_bin[10], FR_n_eta_bin[10];
  double get_FR(snu::KParticle muon, TString whichFR, int n_jets, bool geterror);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu_CR_FR_method, 1);
};
#endif
