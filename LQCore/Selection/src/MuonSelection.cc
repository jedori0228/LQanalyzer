#include "MuonSelection.h"

using namespace snu;

MuonSelection::MuonSelection(LQEvent ev) :
  BaseSelection()
{
  k_lqevent = ev;  
};

MuonSelection::~MuonSelection() {};



void MuonSelection::BasicSelection( std::vector<KMuon>& leptonColl, bool m_debug) {

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++)
    {

      if(muit->Pt() == 0.) continue;

      bool pass_selection = true;

      /// ONLY CUT ON PT/ETA/LOOSE ID
      if( muit->Pt() < pt_cut_min ) {
	pass_selection = false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Pt cut " << endl; 
      }
      if( fabs(muit->Eta()) > eta_cut){
	pass_selection =false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Eta cut " << endl; 
      }
      if(! (PassID(MUON_POG_LOOSE, *muit, m_debug))){
	pass_selection =false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Loose Selection" << endl;
      }
      
      if(pass_selection) leptonColl.push_back(*muit);
    }
  if(m_debug) cout << "BasicSelection::Number of muons = " << leptonColl.size() << endl;
  return;
}


void MuonSelection::SkimSelection( std::vector<KMuon>& leptonColl , bool m_debug) {

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++)
    {
      
      if(muit->Pt() == 0.) continue;

      bool pass_selection = true;
      
      /// ONLY CUT ON PT/ETA/LOOSE ID
      if( muit->Pt() < pt_cut_min ) {
        pass_selection = false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Pt cut " << endl;
      }
      if( fabs(muit->Eta()) > eta_cut){
        pass_selection =false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Eta cut " << endl;
      }
      if(! (PassID(MUON_POG_LOOSE, *muit, m_debug))){
        pass_selection =false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Loose Selection" << endl;
      }

      if (pass_selection) leptonColl.push_back(*muit);  
    }
  if(m_debug) cout << "SkimSelection::Number of muons = " << leptonColl.size() << endl;
  return;

}

void MuonSelection::Selection( std::vector<KMuon>& leptonColl, bool m_debug) {
  

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++)
    { 

      bool pass_selection(true);      
      if(muit->Pt() == 0.) continue;
      
      //// Calculate PF isolation
      /// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
      LeptonRelIso = (muit->RelIso03());
      
      if(apply_relisocut && !( LeptonRelIso < relIso_cut)) pass_selection = false;
      if(m_debug&&apply_relisocut && !( LeptonRelIso < relIso_cut))  cout << "Fails Selection::reliso cut " << endl;
      
      if(apply_ptcut && ! ( muit->Pt() > pt_cut_min )) pass_selection = false;
      if(m_debug&&apply_ptcut && ! (muit->Pt() >= pt_cut_min && muit->Pt() < pt_cut_max)) cout << "Fails Selection::pt cut " << endl;
      
      if(apply_etacut && !(fabs(muit->Eta()) < eta_cut)) pass_selection =false;
      if(m_debug&&apply_etacut && !(fabs(muit->Eta()) < eta_cut))  cout << "Fails Selection::eta cut " << endl;
  
      /// impact parameter cuts
      // Uses fabs(recoMu.muonBestTrack()->dxy(vertex->position())) as described in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon 
      // Also stores dB() on pat::Muon  as dxy_pat
      // Also stores D0 

      if(apply_dzcut && !(fabs(muit->dZ())<  dz_cut )) pass_selection = false;
      if(m_debug&&apply_dzcut && !(fabs(muit->dZ())<  dz_cut ))  cout << "Fails Selection::dz cut " << endl;
      if(apply_dxycut && !(fabs(muit->dXY())< dxy_cut )) pass_selection = false;
      if(m_debug&&apply_dxycut && !(fabs(muit->dXY())< dxy_cut ))cout << "Fails Selection::dxy cut " << endl;
      
      if(apply_ID && !PassID(k_id, *muit,m_debug)) pass_selection =false;
      if(m_debug&& apply_ID && !PassID(MUON_POG_TIGHT, *muit)) cout << "Fails Selection::ID cut " << endl;
      
      
      if(apply_chi2cut && !(muit->GlobalChi2() <chiNdof_cut) && !( muit->GlobalChi2()  >=chiNdofMIN_cut)) pass_selection = false;
      if(m_debug&&apply_chi2cut && !(muit->GlobalChi2() <chiNdof_cut) && !( muit->GlobalChi2()  >=chiNdofMIN_cut)) cout << "Fails chi2 cut " << endl;
      
      //// ADD EXTRA  cut on D0sig? or same vertex?     
      if(pass_selection) leptonColl.push_back(*muit);    

  }/// muon loop end

  BaseSelection::reset();
  return;
}  


////////// PREDEFINED MUON SELECTIONS
  
void MuonSelection::SelectMuons(std::vector<KMuon>& leptonColl, ID muid, float ptcut, float etacut){
  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++){

    bool pass_selection(true);
    if(muit->Pt() == 0.)   continue;

    MuonID = PassUserID(muid, *muit);
    if(!MuonID)  pass_selection = false;

    if(( muit->Pt() < ptcut )) pass_selection = false;
    if(!(fabs(muit->Eta()) < etacut)) pass_selection = false;
    if(pass_selection)  leptonColl.push_back(*muit);
  }
  return;
}


bool MuonSelection::PassUserID(ID id, snu::KMuon mu){

  if ( id == MUON_POG_LOOSE) return POGID(mu, MUON_POG_LOOSE);
  if ( id == MUON_POG_MEDIUM) return POGID(mu, MUON_POG_MEDIUM); 
  if ( id == MUON_POG_TIGHT)return POGID(mu, MUON_POG_TIGHT);
  if ( id == MUON_HN_VETO) return HNVetoMuonSelection(mu);
  if ( id == MUON_HN_FAKELOOSE) return HNLooseMuonSelection(mu);
  if ( id == MUON_HN_TIGHT) return HNTightMuonSelection(mu);
  if ( id == MUON_TOP_VETO) return TopVetoMuonSelection(mu);
  if ( id == MUON_TOP_LOOSE) return TopLooseMuonSelection(mu);
  if ( id == MUON_TOP_TIGHT) return TopTightMuonSelection(mu);
  if ( id == MUON_HN_TRI_NODXYCUT_TIGHT) return HNtriNodXYCutTightMuonSelection(mu);
  if ( id == MUON_HN_TRI_NODXYCUT_LOOSE) return HNtriNodXYCutLooseMuonSelection(mu);
  if ( id == MUON_HN_TRI_TIGHT) return HNtriTightMuonSelection(mu);
  if ( id == MUON_HN_TRI_LOOSE) return HNtriLooseMuonSelection(mu);
  if ( id == MUON_HN_TRI_HIGHDXY_TIGHT) return HNtriHighdXYTightMuonSelection(mu);
  if ( id == MUON_HN_TRI_HIGHDXY_LOOSE) return HNtriHighdXYLooseMuonSelection(mu);
  return false;
  
}


bool MuonSelection::HNVetoMuonSelection(KMuon mu) {

  bool pass_selection(true);
  LeptonRelIso = mu.RelIso03();
  
  if(!( LeptonRelIso < 0.6)) pass_selection = false;
  if(!(mu.GlobalChi2() < 500.)) pass_selection = false;
  if(!(fabs(mu.dZ())< 100.  ))  pass_selection = false;
  if(!(fabs(mu.dXY())< 10.0 )) pass_selection = false;
  if(!PassID(MUON_POG_LOOSE, mu))  pass_selection =false;
  
  //// Make Loose selection
  return pass_selection;
}

bool MuonSelection::HNLooseMuonSelection(KMuon mu) {
  
  //### THIS SELECTION IS USED FOR MUON FAKES STUDIES
  
  bool pass_selection(true);
  
  LeptonRelIso = (mu.RelIso04());
  if(!POGID(mu, MUON_POG_TIGHT)) pass_selection =false;
  if(!(LeptonRelIso < 0.6)) pass_selection = false;
  if(( mu.Pt() < 10. ))  pass_selection = false;
  if(( fabs(mu.Eta()) < 2.4 ))  pass_selection = false;

  return pass_selection;

}

bool MuonSelection::HNIsTight(KMuon muon, bool m_debug){
  
  bool pass_selection(true);
  
  if(muon.Pt() == 0.) return false;
  
  LeptonRelIso = (muon.RelIso04());
  
  /// TIGHT MUON SELECTION
  if(( muon.Pt() < 10. ))  pass_selection = false;
  if(!(fabs(muon.Eta()) < 2.4))   pass_selection =false;
  if(!( LeptonRelIso < 0.1))   pass_selection = false;
  if(!(fabs(muon.dXY())< 0.05 ))  pass_selection = false;

  /// TIGHT MUON from muon POG
  if(!POGID(muon, MUON_POG_TIGHT)) pass_selection =false;

  return pass_selection;
}


bool MuonSelection::POGID(KMuon muon, TString ID){

  bool pass_selection(true);

  if(muon.Pt() == 0.) return false;

  LeptonRelIso = (muon.RelIso04());

  if(ID == MUON_POG_LOOSE){
    if(! muon.IsLoose ())  pass_selection = false;
    if(!( LeptonRelIso < 0.25)) {
      pass_selection = false;
    }
  }
  if(ID== MUON_POG_MEDIUM){
    if(! muon.IsMedium ())  pass_selection = false;
    if(!( LeptonRelIso < 0.25)) {
      pass_selection = false;
    }
  }

  if(ID== MUON_POG_TIGHT){
    if(! muon.IsTight ())  pass_selection = false;
    if(!( LeptonRelIso < 0.15)) {
      pass_selection = false;
    }
  }

  return pass_selection;
}



bool MuonSelection::HNTightMuonSelection(KMuon mu) {
  
  if(HNIsTight(mu))  return true;
  else return false;
}

bool MuonSelection::HNtriNodXYCutTightMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(muon.Pt() == 0.) return false;

  LeptonRelIso = muon.RelIso04();

  /// TIGHT MUON SELECTION
  if(( muon.Pt() < 10. )) pass_selection = false;
  if(!(fabs(muon.Eta()) < 2.4)) pass_selection =false;
  if(!( LeptonRelIso < 0.05)) pass_selection = false;

  /// TIGHT MUON from muon POG
  if(!PassID(MUON_POG_TIGHT, muon)) pass_selection =false;

  return pass_selection;

}

bool MuonSelection::HNtriNodXYCutLooseMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(muon.Pt() == 0.) return false;

  LeptonRelIso = muon.RelIso04();

  /// TIGHT MUON SELECTION
  if(( muon.Pt() < 10. )) pass_selection = false;
  if(!(fabs(muon.Eta()) < 2.4)) pass_selection =false;
  if(!( LeptonRelIso < 0.6)) pass_selection = false;

  /// TIGHT MUON from muon POG
  if(!PassID(MUON_POG_TIGHT, muon)) pass_selection =false;

  return pass_selection;

}

bool MuonSelection::HNtriTightMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(!HNtriNodXYCutTightMuonSelection(muon)) pass_selection = false;
  if(!(fabs(muon.dXY()) < 0.005 )) pass_selection = false;
  if(!(fabs(muon.dXYSig()) < 3. )) pass_selection = false;
  if(!(fabs(muon.dZ()) < 0.10 )) pass_selection = false;

  return pass_selection;

}

bool MuonSelection::HNtriLooseMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(!HNtriNodXYCutLooseMuonSelection(muon)) pass_selection = false;
  if(!(fabs(muon.dXY()) < 0.005 )) pass_selection = false;
  if(!(fabs(muon.dXYSig()) < 3. )) pass_selection = false;
  if(!(fabs(muon.dZ()) < 0.10 )) pass_selection = false;

  return pass_selection;

}

bool MuonSelection::HNtriHighdXYTightMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(!HNtriNodXYCutTightMuonSelection(muon)) pass_selection = false;
  if(!(fabs(muon.dXY()) < 1. )) pass_selection = false;
  if(!(fabs(muon.dXYSig()) > 4. )) pass_selection = false;
  if(!(fabs(muon.dZ()) < 0.10 )) pass_selection = false;

  return pass_selection;

}

bool MuonSelection::HNtriHighdXYLooseMuonSelection(KMuon muon) {

  bool pass_selection(true);

  if(!HNtriNodXYCutLooseMuonSelection(muon)) pass_selection = false;
  if(!(fabs(muon.dXY()) < 1. )) pass_selection = false;
  if(!(fabs(muon.dXYSig()) > 4. )) pass_selection = false;
  if(!(fabs(muon.dZ()) < 0.10 )) pass_selection = false;

  return pass_selection;

}

////////// PREDEFINED MUON SELECTIONS FOR TOP ANALYSIS


bool MuonSelection::TopVetoMuonSelection(KMuon mu) {
  
  bool pass_selection(true);
  LeptonRelIso = (mu.RelIso03());
  
  if(!( LeptonRelIso < 0.2))     pass_selection = false;
  if(!PassID(MUON_POG_LOOSE, mu))      pass_selection =false;
  //// Make Loose selection
  
  return pass_selection;

}

// Intended use for di-lepton channel
bool MuonSelection::TopLooseMuonSelection(KMuon mu) {

  bool pass_selection(true);
  
  if(mu.Pt() == 0.) return false;
  LeptonRelIso = (mu.RelIso03());

  if(!PassID(MUON_POG_LOOSE, mu)) pass_selection =false;
  if(!( LeptonRelIso < 0.2)) pass_selection = false;

  
  return pass_selection;

}

bool MuonSelection::TopTightMuonSelection(KMuon mu) {

  bool pass_selection(true);
  
  if(mu.Pt() == 0.) return false;
  LeptonRelIso = (mu.RelIso03());
  
  if(!( LeptonRelIso < 0.1))    pass_selection = false;
  if(!(fabs(mu.dZ())< 0.10  ))  pass_selection = false;
  if(!(fabs(mu.dXY())< 0.005 ))   pass_selection = false;
  
  /// TIGHT MUON from muon POG
  if(!PassID(MUON_POG_TIGHT, mu)) pass_selection =false;
  if(mu.Pt() == 0.) return false;
  return pass_selection;

}

/// NO LONGER NEEDED
bool MuonSelection::PassID(ID id, snu::KMuon mu, bool m_debug){
  
  
  /// Taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIdRun2
  bool passID(true);
  if (id == MUON_POG_LOOSE) {
    if(!(mu.IsPF() == 1)) {
      passID = false;
      if(m_debug)cout << "PassID: Fail isPF" << endl;
   }
    if(!(mu.IsGlobal()==1 || mu.IsTracker() == 1 )){
      passID = false; 
      if(m_debug){
	cout << "PassID: Fail isGlobal||isTracker" << endl;
	cout << "PassID: mu.IsGlobal()=  " << mu.IsGlobal() << endl;
	cout << "PassID: mu.IsTracker()= " << mu.IsTracker() << endl;
      }
    }
  }


  else if (id == MUON_POG_TIGHT) {
    if(!(mu.IsPF() == 1        )){
      passID = false;
      if(m_debug)cout << "PassID: Fail isPF" << endl;
    }
    if(!(mu.IsGlobal()==1      )) {
      passID = false;
      if(m_debug)cout << "PassID: Fail isGlobal" << endl;
    }
    if( mu.validHits() == 0     ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validHit" << endl;
    }
    if( mu.validPixHits() == 0) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validPixelHit" << endl;
    }
    if( mu.validStations() <= 1 ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validStations" << endl;
    }
    if( mu.ActiveLayer() <= 5   ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail ActiveLayer " << endl;
    }
    //if( fabs(mu.dXY())    >= 0.2) {
    //  passID = false;
    //  if(m_debug)cout << "PassID: Fail dXY" << endl;
    //}
    //if( fabs(mu.dZ())    >= 0.5) {
    //  passID = false;
    //  if(m_debug)cout << "PassID: Fail dZ" << endl;
    //}
    if( mu.GlobalChi2() >=  10.){
      passID = false;
      if(m_debug) cout << "PassID: Fail  Chi2" << endl;
    }
  }

  else{
    cout << "Invalid ID set for muon selection" << endl;
  }
  return passID;
}



void MuonSelection::SetDeposits(Double_t ECalDeposit , Double_t HCalDeposit) {
    ECalDeposit ? ECalDeposit_max = ECalDeposit : ECalDeposit_max=4.0;
    HCalDeposit ? HCalDeposit_max = HCalDeposit : HCalDeposit_max=6.0;
    ECalDeposit_min = 0.0;
    HCalDeposit_min = 0.0;
}

void MuonSelection::SetDeposits(Double_t ECalDeposit1 , Double_t HCalDeposit1, Double_t ECalDeposit2 , Double_t HCalDeposit2) {
    ECalDeposit1 ? ECalDeposit_min = ECalDeposit1 : ECalDeposit_min=0.0;
    HCalDeposit1 ? HCalDeposit_min = HCalDeposit1 : HCalDeposit_min=0.0;
    ECalDeposit2 ? ECalDeposit_max = ECalDeposit2 : ECalDeposit_max=4.0;
    HCalDeposit2 ? HCalDeposit_max = HCalDeposit2 : HCalDeposit_max=6.0;
}


MuonSelection& MuonSelection::operator= (const MuonSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms);
    k_lqevent = ms.k_lqevent;  
  }
  return *this;
};

MuonSelection::MuonSelection(const MuonSelection& ms):
  BaseSelection(ms)
{
  k_lqevent = ms.k_lqevent;  
};

