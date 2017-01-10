#ifndef trilepton_mumumu_ntp_h
#define trilepton_mumumu_ntp_h

#include "AnalyzerCore.h"
class trilepton_mumumu_ntp : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu_ntp();
  ~trilepton_mumumu_ntp();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void find_genparticles();
  int find_genmatching(snu::KParticle gen, std::vector<snu::KMuon> recos, std::vector<int>& used_index);
  void solution_selection_stduy(std::vector<snu::KMuon> recomuons);
  void find_decay(std::vector<snu::KTruth> truthcoll, int target_index, std::vector<int>& indices);
  void print_all_indices(TString particle, std::vector<int> vec);
  int GetSignalMass();
  int n_gen_pass;
  double sol_sel_chi2_best, sol_sel_chi2_plus, sol_sel_chi2_minus, sol_sel_chi2_smaller, sol_sel_chi2_larger;
  bool allgenfound;
  snu::KParticle gen_nu, gen_W_pri, gen_HN, gen_W_sec, gen_l_1, gen_l_2, gen_l_3;
  std::vector<snu::KMuon> sort_muons_ptorder(std::vector<snu::KMuon> muons);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu_ntp, 1);
};
#endif
