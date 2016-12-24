#ifndef trilepton_mumumu_syst_FR_h
#define trilepton_mumumu_syst_FR_h

#include "AnalyzerCore.h"

class trilepton_mumumu_syst_FR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu_syst_FR();
  ~trilepton_mumumu_syst_FR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  vector<double> dXYMins, RelIsoMaxs;
  std::map<TString, TH2D*> hist_trimuon_FR;
  std::map<TString, TH2D*> hist_trimuon_FR_QCD;
  std::map<TString, TH2D*> hist_trimuon_FRSF_QCD;
  std::map<TString, TH2D*> hist_trimuon_FR_QCDSFed;
  int FR_n_pt_bin, FR_n_eta_bin;
  double get_FR(snu::KParticle muon, TString whichFR, bool geterror);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu_syst_FR, 1);
};
#endif
