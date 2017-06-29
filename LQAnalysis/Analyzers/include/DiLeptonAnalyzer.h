#ifndef DiLeptonAnalyzer_h
#define DiLeptonAnalyzer_h

#include "AnalyzerCore.h"

class DiLeptonAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  DiLeptonAnalyzer();
  ~DiLeptonAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  void FillDiLeptonPlot(
    TString histprefix,
    std::vector< KLepton> leptons,
    std::vector< snu::KJet > jets,
    std::vector< snu::KJet > jets_fwd,
    std::vector< snu::KJet > jets_nolepveto,
    double thisweight
  );

  double MET, METphi;
  int nbjets, nbjets_fwd, nbjets_nolepveto;
  int n_vtx;

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs

  ClassDef ( DiLeptonAnalyzer, 1);
};
#endif
