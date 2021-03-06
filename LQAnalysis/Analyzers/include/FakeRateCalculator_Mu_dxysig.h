#ifndef FakeRateCalculator_Mu_dxysig_h
#define FakeRateCalculator_Mu_dxysig_h

#include "AnalyzerCore.h"


class FakeRateCalculator_Mu_dxysig : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_Mu_dxysig();
  ~FakeRateCalculator_Mu_dxysig();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  double GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, std::vector<snu::KMuon> muons, int npfjet50);
  void FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins=0);
  void FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy);
  void FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight);
  void FillDenAndNum(TString prefix, snu::KMuon muon, std::map<TString, double> hltweight, bool isTight);

  float GetPrescale(std::vector<snu::KMuon> muon, bool passlow, bool passhigh);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
  TH2D *FR_sampleA;
  double METauto;


  ClassDef ( FakeRateCalculator_Mu_dxysig, 1);
};
#endif
