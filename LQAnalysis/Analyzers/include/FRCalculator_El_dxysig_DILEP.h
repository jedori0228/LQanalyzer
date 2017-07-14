#ifndef FRCalculator_El_dxysig_DILEP_h
#define FRCalculator_El_dxysig_DILEP_h

#include "AnalyzerCore.h"


class FRCalculator_El_dxysig_DILEP : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FRCalculator_El_dxysig_DILEP();
  ~FRCalculator_El_dxysig_DILEP();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);

  double GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, std::vector<snu::KElectron> electrons, int npfjet50);
  void FillHistByTrigger(TString histname, float value, std::map<TString, double> hltweight, float xmin, float xmax, int nbins=0);
  void FillHistByTrigger(TString histname, float value1, float value2, std::map<TString, double> hltweight , float x[], int nbinsx, float y[], int nbinsy);
  void FillDenAndNum(TString prefix, snu::KElectron electron, double thisweight, bool isTight);
  void FillDenAndNum(TString prefix, snu::KElectron electron, std::map<TString, double> hltweight, bool isTight);

  float GetPrescale(std::vector<snu::KElectron> electron, bool passlow, bool passhigh);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  TH2D *FR_sampleA;
  double METauto;


  ClassDef ( FRCalculator_El_dxysig_DILEP, 1);
};
#endif
