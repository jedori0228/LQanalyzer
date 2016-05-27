#ifndef trilepton_mumumu_h
#define trilepton_mumumu_h

#include "AnalyzerCore.h"
#include "Trilepton.h"

class trilepton_mumumu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  trilepton_mumumu();
  ~trilepton_mumumu();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void gen_matching();
  void find_decay(std::vector<snu::KTruth> truthcoll, int target_index, std::vector<int>& indices);
  void print_all_indices(TString particle, std::vector<int> vec);
  int n_gen_pass;
  double sol_sel_chi2_best, sol_sel_chi2_plus, sol_sel_chi2_minus, sol_sel_chi2_smaller, sol_sel_chi2_larger;

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( trilepton_mumumu, 1);
};
#endif
