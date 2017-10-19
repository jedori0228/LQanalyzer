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
    std::vector< snu::KFatJet > fatjets,
    double thisweight,
    double thieweighterr
  );

  TString MuonTightID, ElectronTightID;
  TString MuonLooseID_loosest, ElectronLooseID_loosest;
  TString MuonVetoID_loosest, ElectronVetoID_loosest;

  //==== FR
  double weight_fr, weight_err_fr;
  int NTightLeptons;
  //==== Prompt Rate
  TH2D *hist_Muon_PR, *hist_Electron_PR;
  //==== Fake Rate
  TH2D *hist_Muon_FR,                 *hist_Electron_FR;
  TH2D *hist_Muon_FR_withbjet,        *hist_Electron_FR_withbjet;
  TH2D *hist_Muon_FR_withoutbjet,     *hist_Electron_FR_withoutbjet;
  TH2D *hist_Muon_FR_QCD,             *hist_Electron_FR_QCD;
  TH2D *hist_Muon_FR_QCD_withbjet,    *hist_Electron_FR_QCD_withbjet;
  TH2D *hist_Muon_FR_QCD_withoutbjet, *hist_Electron_FR_QCD_withoutbjet;
  //==== awaypt change
  TString MuFR_key, ElFR_key;
  std::map< TString, TH2D* > hist_Muon_FR_syst;
  std::map< TString, TH2D* > hist_Electron_FR_syst;

  double CorrPt(TLorentzVector lep, double T_iso);
  double CorrPt(KLepton lep, double T_iso);
  double CorrPt(snu::KMuon lep, double T_iso);
  double CorrPt(snu::KElectron lep, double T_iso);
  double GetMuonFR(bool geterr, float pt,  float eta, int NearBjet);
  double GetMuonPR(bool geterr, float pt,  float eta);
  double GetElectronFR(bool geterr, float pt,  float eta, int NearBjet);
  double GetElectronPR(bool geterr, float pt,  float eta);
  void get_eventweight(std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<bool> isT, std::vector<int> NearBjet, int HalfSampleErrorDir);

  //==== CF
  void GetCFWeight(KLepton lep1, KLepton lep2);
  double GetCF(KLepton lep, bool geterr);
  double weight_cf, weight_err_cf;

  double MET, METphi, ST, HT, LT, contramass;
  int nbjets, nbjets_fwd, nbjets_nolepveto;
  int n_vtx;
  int index_jjW_j1, index_jjW_j2;
  int index_lljjW_j1, index_lljjW_j2;
  int index_fjW;

  bool JSFatJetID(snu::KFatJet fatjet);
  bool IsAwayFromFatjet(snu::KJet jet, vector<snu::KFatJet> fatjets);

  bool RunNtp;

  double GetDijetMassClosest(std::vector<snu::KJet> js, double mass, int& m, int& n);
  double GetDileptonDijetMassClosest(std::vector<KLepton> leps, std::vector<snu::KJet> js, double mass, int& m, int& n);
  double GetFatjetMassClosest(std::vector<snu::KFatJet> fjs, double mass, int& m);

  //==== PDFs
  vector<float> *ForTree_PdfWeights;
  vector<float> *ForTree_ScaleWeights;

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs

  ClassDef ( DiLeptonAnalyzer, 1);
};
#endif
