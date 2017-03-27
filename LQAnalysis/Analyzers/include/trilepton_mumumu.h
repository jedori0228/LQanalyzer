#ifndef trilepton_mumumu_h
#define trilepton_mumumu_h

#include "AnalyzerCore.h"

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
  int find_genmatching(snu::KParticle gen, std::vector<snu::KMuon> recos, std::vector<int>& used_index);
  int GetSignalMass();
  std::vector<snu::KMuon> sort_muons_ptorder(std::vector<snu::KMuon> muons);
  bool PassOptimizedCut(int sig_mass, double first_pt, double second_pt, double third_pt, double W_pri_mass, double PFMET);

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
