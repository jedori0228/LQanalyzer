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
  void FillCutFlowByName(TString histname, TString cut, float w, bool IsDATA);

  void FillDiLeptonPlot(
    TString histprefix,
    std::vector< KLepton > leptons,
    std::vector< snu::KJet > jets,
    std::vector< snu::KJet > jets_fwd,
    std::vector< snu::KJet > jets_nolepveto,
    double thisweight,
    double thieweighterr
  );

  //==== FR
  double weight_fr, weight_err_fr;
  TH2D *hist_Muon_FR, *hist_Muon_PR, *hist_Electron_FR, *hist_Electron_PR;
  double GetMuonFR(bool geterr, float pt,  float eta);
  double GetMuonPR(bool geterr, float pt,  float eta);
  double GetElectronFR(bool geterr, float pt,  float eta);
  double GetElectronPR(bool geterr, float pt,  float eta);
  void get_eventweight(std::vector<TLorentzVector> muons, TString muid, std::vector<TLorentzVector> electrons, vector<TString> elcut, std::vector<bool> isT, int HalfSampleErrorDir);

  //==== CF
  TH2D *hist_CF;
  void GetCFWeight(KLepton lep1, KLepton lep2);
  double GetCF(KLepton lep, bool geterr);
  double weight_cf, weight_err_cf;

  double MET, METphi;
  int nbjets, nbjets_fwd, nbjets_nolepveto;
  int n_vtx;

  double GetDijetMassClosest(std::vector<snu::KJet> js, double mass);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs

  ClassDef ( DiLeptonAnalyzer, 1);
};
#endif
