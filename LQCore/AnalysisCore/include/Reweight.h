#ifndef Reweight_h
#define Reweight_h

class TFile;
class TH1D;
class TString;
class TDirectory;

class Reweight {

 public:
  
  //// constructors
  Reweight(TString filename);
  ~Reweight();
  
  double GetWeight(int nvtx);
  TDirectory* getTemporaryDirectory(void) const;

 private:

  TFile* fileData_;
  TFile* fileMC_;

  TH1D* h_MCmod_;
  TH1D* h_Data_;

};

#endif

