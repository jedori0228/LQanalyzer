#ifndef trilepton_mumumu_CR_FR_method_h
#define trilepton_mumumu_CR_FR_method_h

#include "AnalyzerCore.h"
#include "Trilepton.h"

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
  void SetBinInfo(int cut);
  double HN_x_min, HN_x_max, HN_dx,
         W_pri_lowmass_x_min, W_pri_lowmass_x_max, W_pri_lowmass_dx,
         dR_x_min, dR_x_max, dR_dx,
         gamma_star_x_min, gamma_star_x_max, gamma_star_dx,
         z_candidate_x_min, z_candidate_x_max, z_candidate_dx;

  TH2F* hist_trimuon_FR;
  int FR_n_pt_bin, FR_n_eta_bin;
  double get_FR(snu::KParticle muon);

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
