#ifndef FRCalculator_Mu_dxysig_DILEP_h
#define FRCalculator_Mu_dxysig_DILEP_h

#include "AnalyzerCore.h"


class FRCalculator_Mu_dxysig_DILEP : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FRCalculator_Mu_dxysig_DILEP();
  ~FRCalculator_Mu_dxysig_DILEP();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  double GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KMuon> muons, int npfjet50);
  void FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins=0);
  void FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy);
  void FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight);
  void FillDenAndNum(TString prefix, snu::KMuon muon, std::map<TString, double> hltweight, bool isTight);

  float GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh);
  float JSWeightByTrigger(TString triggername, float tlumi);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
  TH2D *FR_sampleA;
  double METauto, METphiauto;


  ClassDef ( FRCalculator_Mu_dxysig_DILEP, 1);
};
#endif
