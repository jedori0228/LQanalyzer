#ifndef FakeRateCalculator_Mu_h
#define FakeRateCalculator_Mu_h

#include "AnalyzerCore.h"


class FakeRateCalculator_Mu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_Mu();
  ~FakeRateCalculator_Mu();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  double GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, std::vector<snu::KMuon> muons);
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


  ClassDef ( FakeRateCalculator_Mu, 1);
};
#endif
