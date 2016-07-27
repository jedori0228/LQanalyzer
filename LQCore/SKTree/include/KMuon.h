#ifndef _SKTree_KMuon_H__
#define _SKTree_KMuon_H__

/// Local includes
#include "KParticle.h"

// STD includes
#include <string>

namespace snu {
  
  class KMuon : public KParticle {
  public:
    
    KMuon();
    
    ///Copy constructor
    KMuon(const KMuon& muon);
    
    ///Destructor    
    virtual ~KMuon() ;
    
    KMuon& operator= (const KMuon& obj);

    ///Return the type of this object, i.e. KMuon.h              
    virtual std::string Type() const;
    
    void SetPtErr(double pterror);
    void SetEtaErr(double etaerror);
    void SetISOR03ChargedHad(double isor03ch );
    void SetISOR03NeutralHad(double isor03n );
    void SetISOR03Photon(double isor03ph );
    void SetISOR04ChargedHad(double isor04ch );
    void SetISOR04NeutralHad(double isor04n );
    void SetISOR04Photon(double isor04ph );
    void SetIsolationEcalVeto(double isoEcalveto );
    void SetIsolationHcalVeto(double isoHcalveto );
    void SetPileUp_R03(double pileupr03);
    void SetPileUp_R04(double pileupr04);
    void SetTrackVx(double vtx);
    void SetTrackVy(double vty);
    void SetTrackVz(double vtz);
    void Setdz(double dz);
    void Setdxy(double dxy);
    void Setdxy_pat(double dxypat);
    void Setdxyerr_pat(double dxyerrpat);
    void SetD0(double d0);
    void SetD0Error(double d0err);
    void SetGlobalchi2(double glob_chi2);
    void SetValidHits(int validhits);
    void SetPixelValidHits(int valid_pix_hits);
    void SetValidStations(int validstations);
    void SetLayersWithMeasurement(int layer_with_meas);
    void SetISPF(int ispf);
    void SetIsGlobal(int isglobal);
    void SetIsTracker(int istracker);
    void SetMuonMatchedGenParticlePt(double muonmatch_genpt);
    void SetMuonMatchedGenParticlePhi(double muonmatch_genphi);
    void SetMuonMatchedGenParticleEta(double muonmatch_geneta);
    void SetMuonVtxIndex(int ivertex);
    void SetVertexDistXY(double vdistxy);
    
    void SetHLTDoubleMuMatched(bool match);
    void SetHLTSingleMuMatched(bool match);
    void SetHLTSingleMuMatched5(bool match);
    void SetHLTSingleMuMatched8(bool match);
    void SetHLTSingleMuMatched12(bool match);
    void SetHLTSingleMuMatched17(bool match);
    void SetHLTSingleMuMatched24(bool match);
    void SetHLTEMuMatched8(bool match);
    void SetHLTEMuMatched17(bool match);
    void SetHLTSingleMuIsoMatched(bool match);
    
    
    void SetMuonIDPt(float pt);
    void SetMuonIDEta(float eta);
    void SetMuonIDPhi(float phi);
    void SetMuonIDCharge(int charge);
    
    void SetMuonMSPt(float pt);
    void SetMuonMSE(float e);
    void SetMuonMSEta(float eta);
    void SetMuonMSPhi(float phi);
    void SetMuonMSCharge(int charge);
    
    void SetMuonCocktailPt(double pt);
    void SetMuonCocktailEta(double eta);
    void SetMuonCocktailPhi(double phi);
    void SetMuonCocktailGlobalChi2(double chi2);
    void SetMuonCocktailTrkD0(double d0);
    void SetMuonCocktailTrkDz(double dz);
    void SetMuonCocktailCharge(double charge);
    

    void SetShiftedEUp(double pt_up);
    void SetShiftedEDown(double pt_down);

    inline Int_t IsPF() const {return k_muon_ispf;}
    inline Int_t IsGlobal() const {return k_muon_isglobal;}
    inline Int_t IsTracker() const {return k_muon_istracker;}
    inline Int_t validHits() const {return k_muon_valid_hits;}
    inline Int_t validPixHits() const {return k_muon_valid_pixhits;}
    inline Int_t validStations() const {return k_muon_valid_stations;}
    inline Int_t ActiveLayer() const {return k_muon_layer_with_meas;}
    
    //// Truth
    inline Double_t MuonMatchedGenParticlePt() const {return k_muongen_pt;}
    inline Double_t MuonMatchedGenParticleEta() const {return k_muongen_eta;}
    inline Double_t MuonMatchedGenParticlePhi() const {return k_muongen_phi;}
      
    inline Double_t PtError() const {return k_pterror;}
    inline Double_t EtaError() const {return k_etaerror;}

    inline Double_t SumIsoCHDR03() const {return k_isor03ch;} /// charged hadron
    inline Double_t SumIsoNHDR03() const {return k_isor03n;} /// neutral hadron
    inline Double_t SumIsoPHDR03() const {return k_isor03ph;} /// photon
    inline Double_t SumIsoCHDR04() const {return k_isor04ch;} /// charged hadron
    inline Double_t SumIsoNHDR04() const {return k_isor04n;} /// neutral hadron
    inline Double_t SumIsoPHDR04() const {return k_isor04ph;} /// photon
    inline Double_t IsoEcalVeto() const {return k_isoEcalveto;}
    inline Double_t IsoHcalVeto() const {return k_isoHcalveto;}
    inline Double_t SumPUIsoR03() const {return k_MuonPFIsoR03PU;}
    inline Double_t SumPUIsoR04() const {return k_MuonPFIsoR04PU;}
    inline Double_t muonVtx() const {return k_muonVtx;}
    inline Double_t muonVty() const {return k_muonVty;}
    inline Double_t muonVtz() const {return k_muonVtz;}
    inline Double_t dZ() const {return k_dz;}
    inline Double_t dXY() const {return k_dxy;}
    inline Double_t dXYPat() const {return k_dxy_pat;}
    inline Double_t dXYErrPat() const {return k_dxyerr_pat;}
    inline Double_t D0() const {return k_d0;}
    inline Double_t D0Err() const {return k_d0err;}
    inline Double_t GlobalChi2() const {return k_globmuon_chi2;}
    inline Double_t VertexDistXY() const {return k_vtxdistxy;}
    inline Int_t MuonVertexIndex() const {return i_muonVtx;}

    /// Trigger Matching
    inline Bool_t   MatchedDiMuonTrigger()  const {return k_dimuon_trig_match;}
    inline Bool_t   MatchedSingleMuonTrigger()  const {return k_muon_trig_match;}
    inline Bool_t   MatchedSingleMuonTrigger5()  const {return k_muon_trig_match5;}
    inline Bool_t   MatchedSingleMuonTrigger8()  const {return k_muon_trig_match8;}
    inline Bool_t   MatchedSingleMuonTrigger12()  const {return k_muon_trig_match12;}
    inline Bool_t   MatchedSingleMuonTrigger17()  const {return k_muon_trig_match17;}
    inline Bool_t   MatchedSingleMuonTrigger24()  const {return k_muon_trig_match24;}
    inline Bool_t   MatchedElMuonTrigger8()  const {return k_emu_trig_match8;}
    inline Bool_t   MatchedElMuonTrigger17()  const {return k_emu_trig_match17;}
    inline Bool_t   MatchedSingleIsoMuonTrigger()  const {return k_isomuon_trig_match;}
      
      
    inline Double_t MuonMSPt() const {return muon_ms_pt;}
    inline Double_t MuonMSEta() const {return muon_ms_eta;}
    inline Double_t MuonMSE() const {return muon_ms_e;}
    inline Double_t MuonMSPhi() const {return muon_ms_phi;}
    inline Double_t MuonMSCharge() const {return muon_ms_charge;}
    
    inline Double_t MuonIDPt() const {return muon_id_pt;}
    inline Double_t MuonIDEta() const {return muon_id_eta;}
    inline Double_t MuonIDPhi() const {return muon_id_phi;}
    inline Double_t MuonIDCharge() const {return muon_id_charge;}

    inline Double_t MuonCocktailPt() const {return muon_ck_pt;}
    inline Double_t MuonCocktailEta() const {return muon_ck_eta;}
    inline Double_t MuonCocktailPhi() const {return muon_ck_phi;}
    inline Double_t MuonCocktailGlobalChi2() const {return muon_ck_chi2;}
    inline Double_t MuonCocktailTrkD0() const {return muon_ck_trkd0;}
    inline Double_t MuonCocktailTrkDz() const {return muon_ck_trkdz;}
    inline Double_t MuonCocktailCharge() const {return muon_ck_charge;}
    
    inline Double_t IsoTerm() const {return (k_isor03ch + std::max(0.0, k_isor03n + k_isor03ph - k_MuonPFIsoR03PU));}
    
    inline Double_t PtShiftedUp() const {return muon_pt_up;}
    inline Double_t PtShiftedDown() const {return muon_pt_down;}

    Double_t LeptonRelIso();
    
  protected:
    /// Reset function.                                                                  
    virtual void Reset();    
    
  private:
    /// decalre private functions
  
    Double_t k_pterror, k_etaerror , k_isor03ch, k_isor03n, k_isor03ph, k_isor04ch, k_isor04n, k_isor04ph;
    Double_t k_isoEcalveto , k_isoHcalveto , k_MuonPFIsoR03PU,  k_MuonPFIsoR04PU, k_muonVtx, k_muonVty, k_muonVtz,k_muongen_pt,k_muongen_eta,k_muongen_phi;
    Double_t k_dz, k_dxy , k_d0, k_d0err, k_globmuon_chi2,  k_dxy_pat,  k_dxyerr_pat,k_vtxdistxy,k_reliso;
    Int_t k_muon_valid_hits, k_muon_valid_pixhits, k_muon_valid_stations, k_muon_layer_with_meas;
    Int_t k_muon_ispf, k_muon_isglobal, k_muon_istracker;
    Int_t i_muonVtx;
    Double_t muon_ck_pt, muon_ck_eta, muon_ck_phi, muon_ck_chi2, muon_ck_trkd0, muon_ck_trkdz, muon_ck_charge;

    Bool_t k_dimuon_trig_match,k_muon_trig_match,k_isomuon_trig_match, k_muon_trig_match5,k_muon_trig_match8,k_muon_trig_match12,k_muon_trig_match17,k_muon_trig_match24,k_emu_trig_match8, k_emu_trig_match17;
    Double_t muon_ms_pt, muon_ms_eta, muon_ms_phi, muon_ms_charge,  muon_ms_e;
    Double_t muon_id_pt, muon_id_eta, muon_id_phi, muon_id_charge;
    
    Double_t muon_pt_up, muon_pt_down;

    ClassDef(KMuon,9)
  };   
}//namespace snu

#endif
