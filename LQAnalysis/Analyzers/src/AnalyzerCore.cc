/***************************************************************************
 * @Project: LQAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes 
#include "AnalyzerCore.h"
#include "EventBase.h"

//Plotting                                                      
#include "MuonPlots.h"
#include "ElectronPlots.h"
#include "JetPlots.h"
#include "SignalPlotsEE.h"
#include "SignalPlotsMM.h"
#include "SignalPlotsEM.h"
#include "TriLeptonPlots.h"
#include "HNpairPlotsMM.h"
#include "HNTriLeptonPlots.h"

//ROOT includes
#include <TFile.h>
#include "TStyle.h"

AnalyzerCore::AnalyzerCore() : LQCycleBase(), n_cutflowcuts(0), MCweight(-999.),reset_lumi_mask(false),changed_target_lumi(false), k_reset_period(false), a_mcperiod(-1) {

  k_debugmode=false;
  IDSetup=false;  
  setupDDBkg=false;

  //// set rare collections to false. Can be set truth in analysis code
  k_usegenjet=true;
  k_usetruth=true;
  k_usephotons=true;
  k_usefatjet=true;

  TH1::SetDefaultSumw2(true);  
  /// clear list of triggers stored in KTrigger
  triggerlist.clear();
  cutflow_list.clear();
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("AnalyzerCore");

  Message("In AnalyzerCore constructor", INFO);
  
  /////////////////////////////////////////////////////////////////////// 
  //////// For HN analysis  /////////////////////////////////////////////  
  //////////////////////////////////////////////////////////////////////  
  //// MC Data corrections
  //////////////////////////////////////////////////////////////////////                                                                                                   
  
  //mcdata_correction = new MCDataCorrections();
  
  string lqdir =  getenv("LQANALYZER_DIR");

  string username = getenv("USER");

  if( TString(getenv("CATDEBUG")) == "True") k_debugmode=true;

  
  cout << "                                                  " << endl;
  cout << "   ########    ###       ###   ###  ###           " << endl;
  cout << "   ########    ####      ###   ###  ###           " << endl;
  cout << "   ###         #####     ###   ###  ###           " << endl;
  cout << "   ###         ### ##    ###   ###  ###           " << endl;
  cout << "   ########    ###  ##   ###   ###  ###           " << endl;
  cout << "   ########    ###   ##  ###   ###  ###           " << endl;
  cout << "        ###    ###    ## ###   ###  ###           " << endl;
  cout << "        ###    ###     #####   ###  ###           " << endl;
  cout << "   ########    ###      ####   ########           " << endl;
  cout << "   ########    ###       ###   ########           " << endl;
  cout << "                                                  " << endl;

  cout << "##################################################" << endl;
  compmap.clear();
  
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    //==== HN Gen Matching Class
    m_HNgenmatch = new HNGenMatching();
  }

  cout <<  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;


}


bool AnalyzerCore::CheckEventComparison(TString user, TString label){
  
  if(compmap.size() ==0){
    ifstream comp(( "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user+"/" + label + ".txt"));
    if(!comp) {
      exit(EXIT_FAILURE);
    }
    
    string lline;
    while(getline(comp,lline) ){
      std::istringstream is( lline );
      TString blank;
      int run;
      TString tmp;
      int ev;
      is >>blank;
      is >> run;
      is >> tmp;
      is >> ev;
      if(blank!=user) break;
      
      compmap[ev] = run;
      continue;
    }
  }
  
  for(map<int,int>::iterator mit = compmap.begin(); mit != compmap.end(); mit++){
    if(mit->first == eventbase->GetEvent().RunNumber() && mit->second == eventbase->GetEvent().EventNumber() ) return true;
  }
  
  return false;
}
void AnalyzerCore::FillEventComparisonFile(TString label){
  
  //// Make TEX file                                                                                                                                                        
  ofstream ofile_tex;
  string lqdir = getenv("LQANALYZER_DIR");

  string compfile = "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  string(getenv("USER")) + "/"+string(label)+ ".txt";     

  ofile_tex.open(compfile.c_str());
  ofile_tex.setf(ios::fixed,ios::floatfield);
  ofile_tex <<getenv("USER") << " " << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << endl;

}

vector<TString >  AnalyzerCore::GetHNDiLepElTriggers(){

  vector<TString> triglist;
  if(isData){
    if(k_channel.Contains("SingleElectron"))triglist.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    if(k_channel.Contains("DoubleEG")) triglist.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  else{
    triglist.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    triglist.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  return triglist;

}
bool AnalyzerCore::FailHNDataSetCheck(){

  bool _singleEG =(k_channel.Contains("SingleElectron"));
  bool _singleMuon =(k_channel.Contains("SingleMuon"));
  TString analysis_trigger_eg="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  TString analysis_trigger_muon="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  if(isData && _singleEG && PassTrigger(analysis_trigger_eg)) return false;
  if(isData && _singleMuon && PassTrigger(analysis_trigger_muon)) return false;

  return true;
}

void AnalyzerCore::setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.1);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.07, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  
  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

  

}

void AnalyzerCore::SetupLuminosityMap(bool initialsetup, TString forceperiod){
  if(isData) return ;

  TString lumitriggerpath="";
  TString singleperiod = getenv("CATAnalyzerPeriod");
  
  if(!initialsetup) {
    singleperiod=forceperiod;
    trigger_lumi_map_cat2016.clear();
  }
  string lqdir = getenv("LQANALYZER_DIR");
  
  cout << "CATAnalyzerPeriod = " << singleperiod << endl;
  if(singleperiod.Contains("None")){
    lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+".txt";
  }
  else{
    if(singleperiod== "B")   lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_272007_275376.txt";
    else if(singleperiod=="C")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_275657_276283.txt";
    else if(singleperiod=="D")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_277772_278808.txt";
    else if(singleperiod=="E")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276831_277420.txt";
    else if(singleperiod=="F")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276315_276811.txt";
    else if(singleperiod=="G")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_280919_284044.txt";
    else if(singleperiod.Contains("H"))lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_278820_280385.txt";
    else if(singleperiod=="GH")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_280919_280385.txt";
    else {  cerr << "Wrong period setting in SetupLuminosityMap"<< endl;  exit(EXIT_FAILURE);}
  }
  ifstream triglumi2016(lumitriggerpath.Data());
  if(!triglumi2016) {
    cerr << "Did not find period  "+ lumitriggerpath + " exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  
  string sline2016;
  

  if(k_debugmode)cout << "Trigname : Lumi pb-1" << endl;
  while(getline(triglumi2016,sline2016) ){
    std::istringstream is( sline2016 );
    
    string trigname;
    float trig_lumi;
    is >> trigname;
    if(trigname=="###" ) continue;
    
    
    if(trigname=="END") break;
    if(!TString(trigname).Contains("Lumi")){
      is >> trig_lumi;
      if(k_debugmode)cout << trigname << " " << trig_lumi << endl;
    }
    else{
      string tmp;
      is >> tmp;
      is >> trig_lumi;
      if(k_debugmode)cout << trigname << " " << trig_lumi << endl;
    }
    trigger_lumi_map_cat2016[TString(trigname)] = trig_lumi;
    continue;
  }
  triglumi2016.close();
}


bool  AnalyzerCore::Check(float val){

  if(abs(val) == 999.) return false;

  return true;
}



float AnalyzerCore::GetKFactor(){
  
  if(k_sample_name.Contains("ZZTo4L_powheg")) {
    // http://arxiv.org/abs/1405.2219
    //  1.16[3] brings pp->ZZ from NLO to NNLO
    return 1.16;
  }
  if(k_sample_name.Contains("ZZTo2L2Nu_Powheg")) {
    // http://arxiv.org/abs/1405.2219                                                                                     
    //  1.16[3] brings pp->ZZ from NLO to NNLO                                                                            
    return 1.16;
  }
  if(k_sample_name.Contains("ZZTo2L2Q_Powheg")) {
    // http://arxiv.org/abs/1405.2219                                                                                        //  1.16[3] brings pp->ZZ from NLO to NNLO 
    
    return 1.16;
  }

  if(k_sample_name.Contains("ggZZto")){
    //  1.67 brings gg->ZZ from LO to NLO (http://arxiv.org/abs/1509.06734)
    return 1.67;
  }
  if(k_sample_name.Contains("ggHtoZZ")){
    return 1.67;
    //AN2016_359
  }

  return 1.;
    
}

float AnalyzerCore::CorrectedMETRochester( std::vector<snu::KMuon> muall,bool update_met){

  /// function returns corrected met + can be used to set event met to corrected met

  float met_x =eventbase->GetEvent().PFMETx(); 
  float met_y =eventbase->GetEvent().PFMETy();
  
  float px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  for(unsigned int im=0; im < muall.size() ; im++){
      
      px_orig+= muall.at(im).MiniAODPt()*TMath::Cos(muall.at(im).Phi());
      py_orig+= muall.at(im).MiniAODPt()*TMath::Sin(muall.at(im).Phi());
      px_corrected += muall.at(im).Px();
      py_corrected += muall.at(im).Py();
      
  }
  met_x = met_x + px_orig - px_corrected;	
  met_y = met_y + py_orig - py_corrected;	
  
  if(update_met){
    if(!eventbase->GetEvent().PropagatedRochesterToMET()){
      snu::KEvent tempev = eventbase->GetEvent();
      tempev.SetMET(snu::KEvent::pfmet,  sqrt(met_x*met_x + met_y*met_y), eventbase->GetEvent().METPhi(), eventbase->GetEvent().SumET());
      tempev.SetPFMETx(met_x);
      tempev.SetPFMETy(met_y);
      tempev.SetPropagatedRochesterToMET(true);
      eventbase->SetEventBase(tempev);
    }
  }
  return sqrt(met_x*met_x + met_y*met_y);
}   


void AnalyzerCore::CorrectedMETRochester(std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi){

  float met_x = OrignialMET*TMath::Cos(OriginalMETPhi);
  float met_y = OrignialMET*TMath::Sin(OriginalMETPhi);

  float px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  for(unsigned int im=0; im < muall.size() ; im++){

      px_orig+= muall.at(im).MiniAODPt()*TMath::Cos(muall.at(im).Phi());
      py_orig+= muall.at(im).MiniAODPt()*TMath::Sin(muall.at(im).Phi());

      px_corrected += muall.at(im).Px();
      py_corrected += muall.at(im).Py();

  }
  met_x = met_x + px_orig - px_corrected;
  met_y = met_y + py_orig - py_corrected;

  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_x,met_y);

}



float AnalyzerCore::CorrectedMETElectron(std::vector<snu::KElectron> elall, int sys){

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int iel=0; iel < elall.size() ; iel++){


    px_orig+= elall.at(iel).Px();
    py_orig+= elall.at(iel).Py();
    if(sys==1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedDown();
    }


  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;


  return sqrt(met_x*met_x + met_y*met_y);

}

float AnalyzerCore::CorrectedMETMuon( std::vector<snu::KMuon> muall, int sys){
  
  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();
  
  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int imu=0; imu < muall.size() ; imu++){
    
    px_orig+= muall.at(imu).Px();
    py_orig+= muall.at(imu).Py();
    if(sys==1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedDown();
    }  
  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;
  
  
  return sqrt(met_x*met_x + met_y*met_y);
  
}

void AnalyzerCore::CorrectedMETMuon(int sys, std::vector<snu::KMuon> muall, double& OrignialMET, double& OriginalMETPhi){

  if(sys==0) return;

  float met_x = OrignialMET*TMath::Cos(OriginalMETPhi);
  float met_y = OrignialMET*TMath::Sin(OriginalMETPhi);

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int imu=0; imu < muall.size() ; imu++){

    px_orig+= muall.at(imu).Px();
    py_orig+= muall.at(imu).Py();
    if(sys==1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedDown();
    }
  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;

  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_x,met_y);

}

void AnalyzerCore::CorrectedMETElectron(int sys, std::vector<snu::KElectron> elall, double& OrignialMET, double& OriginalMETPhi){

  if(sys==0) return;

  float met_x = OrignialMET*TMath::Cos(OriginalMETPhi);
  float met_y = OrignialMET*TMath::Sin(OriginalMETPhi);

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int imu=0; imu < elall.size() ; imu++){

    px_orig+= elall.at(imu).Px();
    py_orig+= elall.at(imu).Py();
    if(sys==1){
      px_shifted += elall.at(imu).Px()*elall.at(imu).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += elall.at(imu).Px()*elall.at(imu).PtShiftedDown();
    }
  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;

  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_x,met_y);

}

snu::KJet AnalyzerCore::GetCorrectedJetCloseToLepton(snu::KElectron el, snu::KJet jet, bool usem){
  //jet_LepAwareJECv2 = (raw_jet * L1 - lepton) * L2L3Res + lepton
  
  float rawpt= jet.RawPt();  
  float rawe= jet.RawE();
  float L1corr = jet.L1JetCorr();
  float l2l3res = jet.L2L3ResJetCorr();
  float leppt = el.Pt();
  float lepe = el.E();
  float corr_pt = (rawpt*L1corr - leppt)*l2l3res + leppt;
  float corr_e = (rawe*L1corr - lepe)*l2l3res + lepe;
  
  snu::KJet jet_corr(jet);
  
  if(usem){
    TLorentzVector v;
    v.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M());
    v=v*(corr_pt/jet.Pt());
    jet_corr.SetPtEtaPhiM(v.Pt(), v.Eta(), v.Phi(), v.M());
    return jet_corr;
  }
  else jet_corr.SetPtEtaPhiE(corr_pt, jet.Eta(), jet.Phi(),corr_e);

  return jet_corr;
}

snu::KJet AnalyzerCore::GetCorrectedJetCloseToLepton(snu::KMuon mu, snu::KJet jet){

  //jet_LepAwareJECv2 = (raw_jet * L1 - lepton) * L2L3Res + lepton                                                                                                          
  
  float rawpt= jet.RawPt();
  float rawe= jet.RawE();
  float L1corr = jet.L1JetCorr();
  float l2l3res = jet.L2L3ResJetCorr();
  float leppt = mu.Pt();
  float lepe = mu.E();
  float corr_pt = (rawpt*L1corr - leppt)*l2l3res + leppt;
  float corr_e = (rawe*L1corr - lepe)*l2l3res + lepe;

  snu::KJet jet_corr(jet);
  TLorentzVector v;
  v.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M());
  v=v*(corr_pt/jet.Pt());
  jet_corr.SetPtEtaPhiM(v.Pt(), v.Eta(), v.Phi(), v.M());
  return jet_corr;
}


float AnalyzerCore::GetPtRelLepTJet(snu::KElectron electron, std::vector<snu::KJet> jets, bool usecorrectedpt){

  if(jets.size() == 0) return -999.;
  if(electron.Pt() < 10.) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(electron,jets.at(ijet));
      closejet=jets.at(ijet);
      mindR=electron.DeltaR(jets.at(ijet));
    }
  }

  if(mindR==0.7) return 0.;

  FillHist(("ptrel_dr"),mindR, weight, 0., 4., 100);

  TVector3 el3=  electron.Vect();
  TVector3 jet3= closejet.Vect();
  TVector3 lepjetrel = jet3-el3;
  FillHist(("ptrel_lepjetmag"),lepjetrel.Mag(), weight, 0., 100., 100);
  FillHist(("ptrel_crosslepjetmag"), (lepjetrel.Cross(el3)).Mag(), weight, 0., 100., 100);
  float ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();

  return ptrel;
}


float AnalyzerCore::GetPtRelLepTJet(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){

  if(jets.size() == 0) return -999.;
  if(muon.Pt() < 10.) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(muon,jets.at(ijet));
      else closejet = jets.at(ijet);
      mindR=muon.DeltaR(jets.at(ijet));
    }
  }

  if(mindR==0.7) return 0.;

  FillHist(("ptrel_dr"),mindR, weight, 0., 4., 100);

  TVector3 el3=  muon.Vect();
  TVector3 jet3= closejet.Vect();
  TVector3 lepjetrel = jet3-el3;
  FillHist(("ptrel_lepjetmag"),lepjetrel.Mag(), weight, 0., 100., 100);
  FillHist(("ptrel_crosslepjetmag"), (lepjetrel.Cross(el3)).Mag(), weight, 0., 100., 100);
  float ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();

  return ptrel;
}


float AnalyzerCore::GetJetsCloseToLeptonPt(snu::KElectron electron, std::vector<snu::KJet> jets,bool usecorrectedpt){

  float mindR=.4;
  float jetpT=-999.;

  if(electron.Pt() < 10.) return 0.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      mindR=electron.DeltaR(jets.at(ijet));
      if(usecorrectedpt)      jetpT=GetCorrectedJetCloseToLepton(electron,jets.at(ijet)).Pt();
      else jetpT=jets.at(ijet).Pt();
    }
  }

  return jetpT;
}



float AnalyzerCore::GetJetsCloseToLeptonPt(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){
  float mindR=.4;
  float jetpT=-999.;

  if(muon.Pt() < 10.) return 0.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      mindR=muon.DeltaR(jets.at(ijet));
      if(usecorrectedpt)jetpT=GetCorrectedJetCloseToLepton(muon,jets.at(ijet)).Pt();
      else jetpT=jets.at(ijet).Pt();
    }
  }
  return jetpT;
}





float AnalyzerCore::MassDrop(snu::KElectron electron, std::vector<snu::KJet> jets,bool usecorrectedpt){
  if(jets.size() == 0) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;
  if(electron.Pt() < 10.) return -999.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(electron,jets.at(ijet));
      else  closejet= jets.at(ijet);
      
      mindR=mindR;
    }
  }

  if(mindR >= 0.7)  return -999.;

  snu::KParticle lj = closejet+electron;

  return (lj.M() - closejet.M());


}

float AnalyzerCore::MassDrop(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){
  if(jets.size() == 0) return -999.;
  snu::KParticle closejet;
  float mindR=.7;

  if(muon.Pt() < 10.) return -999.;
  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(muon,jets.at(ijet));
      else  closejet= jets.at(ijet);
      mindR=mindR;
    }
  }
  if(mindR >= 0.7)  return -999.;

  snu::KParticle lj = closejet+muon;

  return (lj.M() - closejet.M());


}



// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
///// FUNCTION USED TO CREATE BTAG EFFICIENCIES USED BY BTAGSF.cxx CLass
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void AnalyzerCore::MakeBTagEfficiencyPlots(){
  GetJetTaggerEfficiences("CSVv2M",snu::KJet::CSVv2, snu::KJet::Medium); 
  GetJetTaggerEfficiences("CSVv2T",snu::KJet::CSVv2, snu::KJet::Tight); 
  GetJetTaggerEfficiences("cMVAv2M",snu::KJet::cMVAv2, snu::KJet::Medium); 
  GetJetTaggerEfficiences("cMVAv2L",snu::KJet::cMVAv2, snu::KJet::Loose); 
  GetJetTaggerEfficiences("cMVAv2T",snu::KJet::cMVAv2, snu::KJet::Tight);       
}
void AnalyzerCore::GetJetTaggerEfficiences(TString taggerWP, KJet::Tagger tag,  KJet::WORKING_POINT wp){
  // taken frmo https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#Example_code_AN5  (USE HADRON FLAVOUR)

  // Make tagger efficiency file for btagging 
  Float_t ptbins[7] = { 20., 40., 60., 80., 100., 120., 3000.};
  Float_t etabins[5] = { 0., 0.6, 1.2, 1.8, 2.4};


  for(unsigned int ij =0; ij < GetJets("JET_NOCUT").size(); ij++){
    
    snu::KJet j = GetJets("JET_NOCUT").at(ij);
    int hadronFlavor = j.HadronFlavour();
    if( abs(hadronFlavor)==5 ){
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_b" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp)) FillHist("h2_BTaggingEff_"+taggerWP+"_Num_b" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
    else if( abs(hadronFlavor)==4 ){
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_c" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp))          FillHist("h2_BTaggingEff_"+taggerWP+"_Num_c" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
    else {
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_udsg" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp) )         FillHist("h2_BTaggingEff_"+taggerWP+"_Num_udsg" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
  }
}


double AnalyzerCore::GetIsoCut(snu::KElectron el, TString cutlable){

  if(fabs(el.SCEta())<1.479){
    if(cutlable.Contains("b050")) return 0.05;
    if(cutlable.Contains("b0525")) return 0.0525;
    if(cutlable.Contains("b055")) return 0.055;
    if(cutlable.Contains("b060")) return 0.06;
    if(cutlable.Contains("b065")) return 0.065;
    if(cutlable.Contains("b075")) return 0.075;
    if(cutlable.Contains("b100")) return 0.1;
    if(cutlable.Contains("b125")) return 0.125;
  }
  else { 
    if(cutlable.Contains("e050")) return 0.05;
    if(cutlable.Contains("e0525")) return 0.0525;
    if(cutlable.Contains("e055")) return 0.055;
    if(cutlable.Contains("e060")) return 0.06;
    if(cutlable.Contains("e065")) return 0.065;
    if(cutlable.Contains("e075")) return 0.075;
    if(cutlable.Contains("e100")) return 0.1;
    if(cutlable.Contains("e125")) return 0.125;
  }   
  
  return 0.;
}
  
double AnalyzerCore::GetDXYCut(snu::KElectron el, TString cutlable){

  if(fabs(el.SCEta())<1.479){
    if(cutlable.Contains("b050")) return 0.05;
    if(cutlable.Contains("b040")) return 0.04;
    if(cutlable.Contains("b030")) return 0.03;
    if(cutlable.Contains("b025")) return 0.025;
    if(cutlable.Contains("b020")) return 0.02;
    if(cutlable.Contains("b019")) return 0.019;
    if(cutlable.Contains("b018")) return 0.018;      
    if(cutlable.Contains("b017")) return 0.017;
    if(cutlable.Contains("b016")) return 0.016;
    if(cutlable.Contains("b015")) return 0.015;
    if(cutlable.Contains("b014")) return 0.014;
    if(cutlable.Contains("b013")) return 0.013;
    if(cutlable.Contains("b012")) return 0.012;
    if(cutlable.Contains("b011")) return 0.011;
    if(cutlable.Contains("b010")) return 0.01;
  }
  else{
    if(cutlable.Contains("e100")) return 0.05;
    if(cutlable.Contains("e050")) return 0.05;
    if(cutlable.Contains("e040")) return 0.04;
    if(cutlable.Contains("e035")) return 0.035;
    if(cutlable.Contains("e025")) return 0.025;
    if(cutlable.Contains("e020")) return 0.02;
    if(cutlable.Contains("e019")) return 0.019;
    if(cutlable.Contains("e018")) return 0.018;
    if(cutlable.Contains("e017")) return 0.017;
    if(cutlable.Contains("e016")) return 0.016;
    if(cutlable.Contains("e015")) return 0.015;
    if(cutlable.Contains("e014")) return 0.014;
    if(cutlable.Contains("e013")) return 0.013;
    if(cutlable.Contains("e012")) return 0.012;
    if(cutlable.Contains("e011")) return 0.011;
    if(cutlable.Contains("e010")) return 0.01;
    
  }
  return 0.;
}

int AnalyzerCore::GetDataPeriod(){
  
  /// returns 1 for peiord B.... 7 for period H
  if(eventbase->GetEvent().RunNumber() < 272007) return -1;
  else  if(eventbase->GetEvent().RunNumber() <= 275376) return 1; 
  else  if(eventbase->GetEvent().RunNumber() <= 276283) return 2; 
  else  if(eventbase->GetEvent().RunNumber() <= 276811) return 3; 
  else  if(eventbase->GetEvent().RunNumber() <= 277420) return 4; 
  else  if(eventbase->GetEvent().RunNumber() <= 278808) return 5; 
  else  if(eventbase->GetEvent().RunNumber() <= 280385) return 6; 
  else  if(eventbase->GetEvent().RunNumber() <= 284044) return 7; 
  else return -1;
}

int AnalyzerCore::GetPeriod(){

  if(isData) return GetDataPeriod();
  else return GetMCPeriod();

}

int AnalyzerCore::GetMCPeriod(){
  /// This function returns a period B-H for MC events. 
  /// It uses a random number and retrunds a period based on the luminosity of each period
  /// It assumes the trigger used is unprescaled
  if(isData) return -1;
  if(!k_reset_period) return a_mcperiod;
  k_reset_period=false;
  


  a_mcperiod = k_mcperiod;
  return a_mcperiod;
  
}

int AnalyzerCore::GetMCPeriodRandom(){
  if(isData) return -1;
  double r =gRandom->Rndm(); /// random number between 0 and 1
  
  // https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=1689385956
  // using single muon luminosities (luminosities differ slightky for each dataset but difference is neglibable)
  double lumi_periodB = 5.929001722;
  double lumi_periodC = 2.645968083;
  double lumi_periodD = 4.35344881;
  double lumi_periodE = 4.049732039;
  double lumi_periodF = 3.157020934;
  double lumi_periodG = 7.549615806;
  double lumi_periodH = 8.545039549 + 0.216782873;
  double total_lumi = (lumi_periodB+lumi_periodC + lumi_periodD + lumi_periodE + lumi_periodF + lumi_periodG + lumi_periodH) ;
  
  vector<double> cum_lumi;
  cum_lumi.push_back(lumi_periodB/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH)/total_lumi); 
  
  /// returns an int

  /// r = 1       |  period B
  /// r = 2       |  period C
  /// r = 3       |  period D
  /// r = 4       |  period E
  /// r = 5       |  period F
  /// r = 6       |  period G
  /// r = 7       |  period H

  for(unsigned int i=0; i < cum_lumi.size(); i++){
    if ( r < cum_lumi.at(i)) {
      a_mcperiod =  (i+1);
      return a_mcperiod;
    }
  }

  /// return period H is for some reason r > cum_lumi.at(max) 'should not happen'
  return  cum_lumi.size();

}


void AnalyzerCore::SetupSelectionJet(std::string path_sel){
  Message("SetupSelectionJet", DEBUG);
  ifstream jetselconfig(path_sel.c_str());
  if(!jetselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "Found " + path_sel << endl;
   
  }
  string jetline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(jetselconfig,jetline) ){
    vector<pair<TString,TString> > string_jetsel;
    vector<pair<TString,float> > float_jetsel;
    TString idlabel;
    if (TString(jetline).Contains("webpage")) continue;
    if (TString(jetline).Contains("###")) continue;
    if (TString(jetline).Contains("ncut")) {
      std::istringstream is( jetline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(jetline).Contains("ptmin")) {
      std::istringstream is( jetline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( jetline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 6 && x < 18){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
        else if ( x ==1) {is >> idlabel; cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

        }
      }
    }
    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfJet.find(idlabel);

    if(idlabel!=""){
      
      if(fit != selectionIDMapfJet.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }

    selectionIDMapsJet[idlabel] = string_jetsel;
    selectionIDMapfJet[idlabel] = float_jetsel;
  }
}


void AnalyzerCore::SetupSelectionFatJet(std::string path_sel){
  Message("SetupSelectionJet", DEBUG);
  ifstream jetselconfig(path_sel.c_str());
  if(!jetselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string jetline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(jetselconfig,jetline) ){
    vector<pair<TString,TString> > string_jetsel;
    vector<pair<TString,float> > float_jetsel;
    TString idlabel;
    if (TString(jetline).Contains("webpage")) continue;
    if (TString(jetline).Contains("###")) continue;
    if (TString(jetline).Contains("ncut")) {
      std::istringstream is( jetline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(jetline).Contains("ptmin")) {
      std::istringstream is( jetline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( jetline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 8 && x < 20){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
        else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;
	  
        }
      }
    }

    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfFatJet.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfFatJet.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }

    selectionIDMapsFatJet[idlabel] = string_jetsel;
    selectionIDMapfFatJet[idlabel] = float_jetsel;
  }
}

void AnalyzerCore::SetupSelectionMuon(std::string path_sel){
  Message("SetupSelectionMuon", DEBUG);
  ifstream muonselconfig(path_sel.c_str());
  if(!muonselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string muonline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(muonselconfig,muonline) ){
    vector<pair<TString,TString> > string_muonsel;
    vector<pair<TString,float> > float_muonsel;
    TString idlabel;
    if (TString(muonline).Contains("webpage")) continue;
    if (TString(muonline).Contains("###")) continue;
    if (TString(muonline).Contains("ncut")) {
      std::istringstream is( muonline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(muonline).Contains("ptmin")) {
      std::istringstream is( muonline );
      string tmp;
      for (int x =0; x < ncuts; x++){
	is >> tmp;
	cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( muonline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
	if ( x%2 ==0) {
	  is >> tmp;
	  continue;
	}

	if (x > 12 && x < 18){
	  is >> tmp;
	  string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

	}
	else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
	else if (x > 30 && x < 32){
          is >> tmp;
          string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
	else {
	  is >> tmpf;
	  float_muonsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

	}
      }
    }
    
    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfMuon.end()){
	cerr << "Repeated ID " << idlabel << endl;
	exit(EXIT_FAILURE);
	
      }
    }
    selectionIDMapsMuon[idlabel] = string_muonsel;
    selectionIDMapfMuon[idlabel] = float_muonsel;
  }
}


void AnalyzerCore::SetupSelectionElectron(std::string path_sel){

  Message("SetupSelectionElectron", DEBUG);

  //// currently hard coded set of cuts

  ifstream elselconfig(path_sel.c_str());
  if(!elselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string elline;
  vector<TString> cutnames;
  int ncuts=0;
  while(getline(elselconfig,elline) ){
    vector<pair<TString,TString> > string_elsel;
    vector<pair<TString,float> > float_elsel;
    TString idlabel;
    if (TString(elline).Contains("webpage")) continue;
    if (TString(elline).Contains("###")) continue;
    if (TString(elline).Contains("ncut")) {
      std::istringstream is( elline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }

    if (TString(elline).Contains("ptmin")) {
      std::istringstream is( elline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( elline );
      string tmp;
      float tmpf;

      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 10 && x < 18){
          is >> tmp;
	  cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 26 && x < 30){
          is >> tmp;
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 34 && x < 46){
          is >> tmp;
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
        else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

          float_elsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
      }
    }

    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfElectron.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }
    selectionIDMapsElectron[idlabel] = string_elsel;
    selectionIDMapfElectron[idlabel] = float_elsel;
  }
}




std::map<TString,BTagSFUtil*> AnalyzerCore::SetupBTagger(std::vector<TString> taggers, std::vector<TString> wps){

  //// Btagging code for Moriond17 samples
  
  //// Current use is for HN analyses only.
  //// HN analysis uses method 2 a) from twiki:
  //// https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
  

  /// If users want to use another method please add additional function and modify BTag/ dir in LQLAnalysis/src/

  std::map<TString,BTagSFUtil*>  tmpmap;
  for(std::vector<TString>::const_iterator it = taggers.begin(); it != taggers.end(); it++){
    for(std::vector<TString>::const_iterator it2 = wps.begin(); it2 != wps.end(); it2++){
      if (it->Contains("CSVv2")){
	tmpmap[*it + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", (*it + "_BtoF").Data(), (*it + "_GtoH").Data(),  it2->Data());
	tmpmap[*it +  "_" + *it2 + "_hf"]= new BTagSFUtil("mujets", (*it + "_BtoF").Data(),  (*it + "_GtoH").Data(), it2->Data());
	// tmpmap[*it +  "_" + *it2 + "_hfcomb"]= new BTagSFUtil("comb", it->Data(), it2->Data());                /// SWITCH ON IF USER NEEDS THIS METHOD  
	// tmpmap[*it +  "_" + *it2 + "iterativefit"]= new BTagSFUtil("iterativefit", it->Data(), it2->Data());   /// SWITCH ON IF USER NEEDS THIS METHOD
      }
      if (it->Contains("cMVA")){
	tmpmap[*it  + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", (*it + "_BtoF").Data(),  (*it + "_GtoH").Data(), it2->Data());
	tmpmap[*it  +  "_" + *it2 + "_hf"]= new BTagSFUtil("ttbar", (*it + "_BtoF").Data(),  (*it + "_GtoH").Data(), it2->Data());

	// tmpmap[*it +  "_" + *it2 + "iterativefit"]= new BTagSFUtil("iterativefit", it->Data(), it2->Data());   /// SWITCH ON IF USER NEEDS THIS METHOD       
      }
      if (it->Contains("DeepCVS")){
	/// TO ADD If requested
      }
    }
  }

  return tmpmap;
}


//################################################################################################  
//@@@@@@@@@@@@@@@@@@@  ANALYSIS FUNCTIONALITY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

float AnalyzerCore::GetDiLepMass(std::vector<snu::KElectron> electrons){

  if(electrons.size() != 2) return 0.;
  snu::KParticle p = electrons.at(0) + electrons.at(1);
  return p.M();
}

float AnalyzerCore::GetDiLepMass(std::vector<snu::KMuon> muons){

  if(muons.size() != 2) return 0.;
  snu::KParticle p = muons.at(0) + muons.at(1);
  return p.M();
}


bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KElectron> electrons){
  if(electrons.size() != 2) return false;

  if(reg.Contains("EE")){
    if(fabs(electrons.at(0).Eta()) < 1.5) return false;
    if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(fabs(electrons.at(0).Eta()) > 1.5) return false;
    if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(fabs(electrons.at(0).Eta()) > 1.5){
      if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    }
    if(fabs(electrons.at(0).Eta()) < 1.5){
      if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    }
    return true;
  }
  return false;
  
}


bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KMuon> muons){
  
  if(muons.size() != 2) return false;
  if(reg.Contains("EE")){
    if(muons.at(0).Eta() < 1.5) return false;
    if(muons.at(1).Eta() < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(muons.at(0).Eta() > 1.5) return false;
    if(muons.at(1).Eta() > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(muons.at(0).Eta() > 1.5){
      if(muons.at(1).Eta() > 1.5) return false;
    }
    if(muons.at(0).Eta() < 1.5){
      if(muons.at(1).Eta() < 1.5) return false;
    }
    return true;
  }
  return false;
}



//################################################################################################ 
//@@@@@@@@@@@@@@@@@@@  GET SKTREE  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    


std::vector<snu::KFatJet> AnalyzerCore::GetFatJets(BaseSelection::ID jetid, float ptcut, float etacut){
  return GetFatJets(GetStringID(jetid), ptcut, etacut);
}


std::vector<snu::KJet> AnalyzerCore::GetJets(BaseSelection::ID jetid, float ptcut, float etacut){
  return GetJets(GetStringID(jetid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, bool keepfake, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), keepfake, ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(true, true, GetStringID(electronid), ptcut, etacut);
}




std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake,   BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(keepcf, keepfake, GetStringID(electronid), ptcut, etacut);
}



TString AnalyzerCore::GetStringID(BaseSelection::ID id){
  if(id == BaseSelection::ELECTRON_POG_VETO ) return "ELECTRON_POG_VETO";
  if(id == BaseSelection::ELECTRON_POG_LOOSE) return "ELECTRON_POG_LOOSE";
  if(id == BaseSelection::ELECTRON_POG_MEDIUM) return "ELECTRON_POG_MEDIUM";
  if(id == BaseSelection::ELECTRON_POG_TIGHT           ) return "ELECTRON_POG_TIGHT";
  if(id == BaseSelection::ELECTRON_POG_MVATrig) return "ELECTRON_POG_MVATrig";
  if(id == BaseSelection::ELECTRON_POG_MVANonTrig    ) return "ELECTRON_POG_MVANonTrig";
  if(id == BaseSelection::ELECTRON_ECAL_FIDUCIAL     ) return "ELECTRON_ECAL_FIDUCIAL";
  if(id == BaseSelection::ELECTRON_HN_VETO           ) return "ELECTRON_HN_VETO";
  if(id == BaseSelection::ELECTRON_HN_TIGHT          ) return "ELECTRON_HN_TIGHT";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE      ) return "ELECTRON_HN_FAKELOOSE";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0 ) return "ELECTRON_HN_FAKELOOSE_NOD0";
  if(id == BaseSelection::ELECTRON_TOP_VETO          ) return "ELECTRON_TOP_VETO";
  if(id == BaseSelection::ELECTRON_TOP_LOOSE         ) return "ELECTRON_TOP_LOOSE";
  if(id == BaseSelection::ELECTRON_TOP_TIGHT         ) return "ELECTRON_TOP_TIGHT";
  if(id == BaseSelection::ELECTRON_PTETA             ) return "ELECTRON_PTETA";
  if(id == BaseSelection::ELECTRON_NOCUT             ) return "ELECTRON_NOCUT";
  if(id == BaseSelection::MUON_POG_LOOSE             ) return "MUON_POG_LOOSE";
  if(id == BaseSelection::MUON_POG_MEDIUM            ) return "MUON_POG_MEDIUM";
  if(id == BaseSelection::MUON_POG_TIGHT             ) return "MUON_POG_TIGHT";
  if(id == BaseSelection::MUON_HN_VETO               ) return "MUON_HN_VETO";
 if(id == BaseSelection::MUON_HN_FAKELOOSE          ) return "MUON_HN_FAKELOOSE";
 if(id == BaseSelection::MUON_HN_TIGHT              ) return "MUON_HN_TIGHT";
 if(id == BaseSelection::MUON_FAKELOOSE             ) return "MUON_FAKELOOSE";
 if(id == BaseSelection::MUON_TOP_VETO              ) return "MUON_TOP_VETO";
 if(id == BaseSelection::MUON_TOP_LOOSE             ) return "MUON_TOP_LOOSE";
 if(id == BaseSelection::MUON_TOP_TIGHT             ) return "MUON_TOP_TIGHT";
 if(id == BaseSelection::MUON_PTETA                 ) return "MUON_PTETA";
 if(id == BaseSelection::MUON_NOCUT                 ) return "MUON_NOCUT";
 if(id == BaseSelection::PFJET_LOOSE                ) return "PFJET_LOOSE";
 if(id == BaseSelection::PFJET_MEDIUM               ) return "PFJET_MEDIUM";
 if(id == BaseSelection::PFJET_TIGHT                ) return "PFJET_TIGHT";
 if(id == BaseSelection::JET_HN                     ) return "JET_HN";
 if(id == BaseSelection::JET_HN_TChannel            ) return "JET_HN_TChannel";
 if(id == BaseSelection::JET_NOLEPTONVETO           ) return "JET_NOLEPTONVETO";
 if(id == BaseSelection::JET_LOOSE                  ) return "JET_LOOSE";
 if(id == BaseSelection::JET_TIGHT                  ) return "JET_TIGHT";
 if(id == BaseSelection::PHOTON_POG_LOOSE           ) return "PHOTON_POG_LOOSE";
 if(id == BaseSelection::PHOTON_POG_MEDIUM          ) return "PHOTON_POG_MEDIUM";
 if(id == BaseSelection::PHOTON_POG_TIGHT           ) return "PHOTON_POG_TIGHT";
 cerr << " ID [--] not found" << endl; exit(EXIT_FAILURE);

}

std::vector<snu::KJet> AnalyzerCore::GetJets(TString jetid,float ptcut, float etacut){
  
  std::vector<snu::KJet> jetColl;
  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsJet.find(jetid);
  std::map<TString, vector<pair<TString,float > > >::iterator fit = selectionIDMapfJet.find(jetid);
  if(it== selectionIDMapsJet.end()){
    cerr << "Jet ID ["+jetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }
    
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(jetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetJetSel()->SelectJets(jetColl, GetMuons(muontag), GetElectrons(eltag) , it->second, fit->second, ptcut,etacut);
  }
  return jetColl;
  
}

std::vector<snu::KJet> AnalyzerCore::GetJetsWFT(TString jetid,TString fatjetid,float ptcut, float etacut){

  std::vector<snu::KJet> jetColl;

  std::vector<snu::KFatJet> fatjetColl = GetFatJets(fatjetid);
  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsJet.find(jetid);
  std::map<TString, vector<pair<TString,float > > >::iterator fit = selectionIDMapfJet.find(jetid);
  if(it== selectionIDMapsJet.end()){
    cerr << "Jet ID ["+jetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }

    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(jetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetJetSel()->SelectJets(jetColl, fatjetColl,GetMuons(muontag), GetElectrons(eltag) , it->second, fit->second, ptcut,etacut);
  }
  return jetColl;

}




std::vector<snu::KFatJet> AnalyzerCore::GetFatJets(TString fatjetid, float ptcut, float etacut){

  std::vector<snu::KFatJet> fatjetColl;

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsFatJet.find(fatjetid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfFatJet.find(fatjetid);
  if(it== selectionIDMapsFatJet.end()){
    cerr << "FatJet ID ["+fatjetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetFatJetSel()->SelectFatJets(fatjetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetFatJetSel()->SelectFatJets(fatjetColl, GetMuons(muontag), GetElectrons(eltag) ,  it->second, fit->second, ptcut,etacut);
    
  }

  return fatjetColl;

}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, float ptcut, float etacut){
  return GetMuons(muid, true, ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, bool keepfakes, float ptcut, float etacut){

  std::vector<snu::KMuon> muonColl;
  
  if(muid.Contains("NONE")) return muonColl;

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsMuon.find(muid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(muid);
  
  if(1){
    //// This method was 10% faster in processing 40K MC events accessing 1000 vectors                                                                                          
    if(it== selectionIDMapsMuon.end()){
      cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    if(fit== selectionIDMapfMuon.end()){
      cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
    }

    if (ptcut == -999.)  eventbase->GetMuonSel()->SelectMuons(muonColl ,muid, it->second, fit->second);
    else eventbase->GetMuonSel()->SelectMuons(muonColl ,muid, it->second, fit->second, ptcut, etacut);
    return  GetTruePrompt(muonColl, keepfakes);
    
  }
  
  return  GetTruePrompt(muonColl, keepfakes);
  
}


std::vector<snu::KElectron> AnalyzerCore::GetElectrons(TString elid,float ptcut, float etacut){

  if(k_classname.Contains("HNDiElectronOpt")){
    if(k_running_chargeflip)  return GetElectrons( true, false,   elid, ptcut, etacut);
    return GetElectrons(false, false,elid, ptcut, etacut);
  }
  /// if cf flag set and MC keep CF and conversion electrons
  if(k_running_chargeflip)  return GetElectrons( true, false, elid, ptcut, etacut);
  
  
  //// by default keep all electrons
  return GetElectrons(true,  true ,elid, ptcut, etacut);

}



std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake,  TString elid,float ptcut, float etacut){

  std::vector<snu::KElectron> electronColl;
  
  if(elid.Contains("NONE")) return electronColl;

  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(elid);
  
  if(1){
    //// This method was 10% faster in processing 40K MC events accessing 1000 vectors
    if(it== selectionIDMapsElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    if(fit== selectionIDMapfElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    
    if(ptcut == -999.)eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, it->second, fit->second);
    else eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, it->second, fit->second,ptcut, etacut);
    
    return  GetTruePrompt(electronColl, keepcf, keepfake);
  }

  else{
    if(it== selectionIDMapsElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    else {
      
      bool check_cc(false);
      bool check_cv(false);
      TString el_id="";
      vector<pair<TString,TString> > v_string =  it->second;
      for(unsigned int iv=0; iv < v_string.size(); iv++){
	if(v_string[iv].second == "false") continue;
	if(v_string[iv].first.Contains("(POG)")) el_id=v_string[iv].first;
	if(v_string[iv].first.Contains("(MVA)") ) el_id=v_string[iv].first;
	if(v_string[iv].first == "GsfCtfScPix") check_cc=true;
	if(v_string[iv].first == "convveto") check_cv=true;
      }
      
      
      vector<pair<TString,float> > v_float =  fit->second;
      float isomax_b(-999.); 
      float isomax_e (-999.);
      float dxymax_b (-999.);
      float dxymax_e (-999.);
      float dzmax_b (-999.);
      float dzmax_e (-999.);
      
      float dxysigmax(-999.);
      float dxysigmin(-999.);
      
      for(unsigned int iv=0; iv < v_float.size(); iv++){
	if(!Check(v_float[iv].second)) continue;
	if(v_float[iv].first == "ptmin") {
	  if (ptcut == -999.) ptcut =v_float[iv].second;
	}
	if(v_float[iv].first == "|etamax|"){
	  if (etacut == -999.) etacut = v_float[iv].second;
	}
	if(v_float[iv].first == "isomax03_b")isomax_b =v_float[iv].second;
	if(v_float[iv].first == "isomax03_e") isomax_e=v_float[iv].second;
	if(v_float[iv].first == "|dxymax_b|") dxymax_b=v_float[iv].second;
	if(v_float[iv].first == "|dxymax_e|") dxymax_e=v_float[iv].second;
	if(v_float[iv].first == "|dzmax_b|") dzmax_b=v_float[iv].second;
	if(v_float[iv].first == "|dzmax_e|") dzmax_e=v_float[iv].second;
	if(v_float[iv].first == "|dxysigmax|") dxysigmax=v_float[iv].second;
	if(v_float[iv].first == "|dxysigmin|") dxysigmin=v_float[iv].second;
      }
      
      
      eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, el_id,check_cc,check_cv, isomax_b,isomax_e,dxymax_b,dxymax_e,dzmax_b,dzmax_e, dxysigmax,dxysigmin, ptcut, etacut);
    }
  }
  
  //if(elid == "ELECTRON_HN_TIGHT"){
    /// This is the vector of electrons with optimie cuts
    //std::vector<snu::KElectron> _electronColl;
    //if(k_running_nonprompt) eventbase->GetElectronSel()->SelectElectrons(_electronColl, BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0, 15., 2.5);
    //else eventbase->GetElectronSel()->SelectElectrons(_electronColl,BaseSelection::ELECTRON_HN_TIGHT, 15., 2.5);
    //electronColl =ShiftElectronEnergy(_electronColl, k_running_chargeflip);
  // }
   

  return  GetTruePrompt(electronColl, keepcf, keepfake); 

}




bool AnalyzerCore::HasCloseBJet(snu::KElectron el, KJet::Tagger tag, KJet::WORKING_POINT wp, int period){

  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO");

  if(period == 0) {
    Message("period not set in AnalyzerCore::HasCloseBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetPeriod();
  }

  bool cl = false;
  for(unsigned int ij =0; ij < alljets.size(); ij++){
    if(el.DeltaR(alljets.at(ij)) < 0.5){

      if(IsBTagged(alljets.at(ij), tag, wp, period))cl = true;

    }
  }

  return cl;
}



bool AnalyzerCore::TriggerMatch(TString trigname, vector<snu::KMuon> mu){
  
  if(mu.size() == 2){
    if(!mu.at(0).TriggerMatched(trigname)) return false;
    if(!mu.at(1).TriggerMatched(trigname)) return false;
  }
  return true;
}

bool AnalyzerCore::Is2015Analysis(){
  if(TString(getenv("running2015")).Contains("True")) return true;
  else return false;
}



float AnalyzerCore::WeightByTrigger(vector<TString> triggernames, float tlumi){
  
  if(isData){
    for(unsigned int i=0; i < triggernames.size() ; i++){
      //// code here sets weight to -99999. if user tries to use incorrect dataset in data 
      /// datasets for each trigger in 2015 can be seen in following googledoc:
      //// https://docs.google.com/spreadsheets/d/1BkgAHCC4UtP5sddTZ5G5iWY16BxleuK7rqT-Iz2LHiM/pubhtml?gid=0&single=true
       if(!k_channel.Contains("DoubleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_SameSign_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_Mu10_SameSign_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_TripleMu_12_10_5")) return -99999.;
      }
      if(!k_channel.Contains("SingleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu8_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu50_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu22_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu22_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu20_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu24_eta2p1_v" )) return -99999.;
      }
      if(!k_channel.Contains("DoubleEG")){
	if(triggernames.at(i).Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;

      }
      if(!k_channel.Contains("SingleElectron")){
	// Single Electon                                                                                                                                                       
	if(triggernames.at(i).Contains("HLT_Ele23_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_eta2p1_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele25_eta2p1_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele45_WPLoose_Gsf_v")) return -99999.;
      }
      if(!k_channel.Contains("MuonEG")){
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")) return -99999.;
      }
    }
    return 1.;
  }


  float trigps= -1.;
  for(unsigned int i=0; i < triggernames.size() ; i++){
    if(WeightByTrigger(triggernames.at(i), tlumi) > trigps) trigps = WeightByTrigger(triggernames.at(i),tlumi) ;
  }

  //if(trigps  > 1.) {m_logger << ERROR << "Error in getting weight for trigger prescale. It cannot be > 1, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  if(trigps  < 0.) {m_logger << ERROR << "Error in getting weight for trigger " << triggernames.at(0) << " prescale. It cannot be < 0, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  
  return trigps;
 
}

float AnalyzerCore::WeightByTrigger(TString triggername, float tlumi){

  /// Function applies weight to MC 
  /// Depends on trigger 
  /// 

  if(isData) return 1.;
  //  brilcalc lumi -u /pb 
  // --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json 
  // -i jsonfiles/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt 
  //--hltpath ""
  
  /// NUMBERS FROM GetLumi_Triggers.py SCRIPT
  
  /// In v766 path lumi is corrected for removal of bad beamspot LS
  // https://github.com/vallot/CATTools/commit/aae3e60b194b1bacf2595a33c8fa27f411dac16b
  for(map<TString, float>::iterator mit = trigger_lumi_map_cat2016.begin(); mit != trigger_lumi_map_cat2016.end(); mit++){
    if(triggername.Contains(mit->first)) return (mit->second / tlumi);
  }
  m_logger << ERROR << "Error in getting weight for trigger  " << triggername << "  prescale. Trigname is not correct or not in map"  << LQLogger::endmsg; exit(0);

  return 1.;
}



void AnalyzerCore::AddTriggerToList(TString triggername){
  
  triggerlist.push_back(triggername);
}

AnalyzerCore::~AnalyzerCore(){
  

  Message("In AnalyzerCore Destructor" , INFO);

  trigger_lumi_map_cat2016.clear();
  compmap.clear();
  for(map<TString, TH1*>::iterator it = maphist.begin(); it!= maphist.end(); it++){
    delete it->second;
  }
  maphist.clear();

  for(map<TString, TH2*>::iterator it = maphist2D.begin(); it!= maphist2D.end(); it++){
    delete it->second;
  }
  maphist2D.clear();
  
  for(map<TString, TH3*>::iterator it = maphist3D.begin(); it!= maphist3D.end(); it++){
    delete it->second;
  }
  maphist3D.clear();


  for(map<TString, MuonPlots*>::iterator it = mapCLhistMu.begin(); it != mapCLhistMu.end(); it++){
    delete it->second;
  }
  mapCLhistMu.clear();
  

  for(map<TString, JetPlots*>::iterator it = mapCLhistJet.begin(); it != mapCLhistJet.end(); it++){
    delete it->second;
  }
  mapCLhistJet.clear();

  for(map<TString, ElectronPlots*>::iterator it = mapCLhistEl.begin(); it != mapCLhistEl.end(); it++){
    delete it->second;
  }
  mapCLhistEl.clear();

  for(map<TString, SignalPlotsEE*>::iterator it = mapCLhistSigEE.begin(); it != mapCLhistSigEE.end(); it++){
    delete it->second;
  }
  mapCLhistSigEE.clear();

  for(map<TString, SignalPlotsMM*>::iterator it = mapCLhistSigMM.begin(); it != mapCLhistSigMM.end(); it++){
    delete it->second;
  }
  mapCLhistSigMM.clear();


  for(map<TString, SignalPlotsEM*>::iterator it = mapCLhistSigEM.begin(); it != mapCLhistSigEM.end(); it++){
    delete it->second;
  }
  mapCLhistSigEM.clear();

  
  for(map<TString, TriLeptonPlots*>::iterator it = mapCLhistTriLep.begin(); it != mapCLhistTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistTriLep.clear();
  
  for(map<TString, HNpairPlotsMM*>::iterator it = mapCLhistHNpairMM.begin(); it != mapCLhistHNpairMM.end(); it++){
    delete it->second;
  }
  mapCLhistHNpairMM.clear();

  for(map<TString,TNtupleD*>::iterator it = mapntp.begin(); it!= mapntp.end(); it++){ 
    delete it->second;
  }
  mapntp.clear();


  if(!k_classname.Contains("SKTreeMaker")){
    for(std::map<TString,BTagSFUtil*>::iterator it = MapBTagSF.begin(); it!= MapBTagSF.end(); it++){
      delete it->second;
    }
    MapBTagSF.clear();
  }

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistHNTriLep.clear();

  //// New class functions for databkg+corrections
  if(k_classname == "SKTreeMaker")   delete mcdata_correction;
  if(!k_classname.Contains("SKTreeMaker")){
    delete mcdata_correction;
  }

  if(k_running_nonprompt || k_running_chargeflip)  delete m_datadriven_bkg;

  //==== HN Gen Matching
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    delete m_HNgenmatch;
  }
}

//###
//###   IMPORTANT BASE FUNCTION: SETS UP EVENT FOR ALL CYCLES
//###

void AnalyzerCore::SetupID(){

  if(IDSetup) return;

  // IF STANDARD SKTREE CODE NO NEED FOR ID
  if(!k_classname.Contains("HN")){
    if(k_classname.Contains("SKTreeMaker")){
      IDSetup=true;
    return;
    }
  }

  string lqdir =  getenv("LQANALYZER_DIR");

  string username = getenv("USER");

  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/muons.sel");
  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/user_muons.sel");

  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/electrons.sel");
  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/user_electrons.sel");
  //if(k_classname.Contains("HNDiElectron"))SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  //if(k_classname.Contains("FakeRateCalculator_El")) SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  //if(k_classname.Contains("ElectronTypes")) SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/jets.sel");
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/user_jets.sel");

  SetupSelectionFatJet(lqdir + "/CATConfig/SelectionConfig/fatjets.sel");
  SetupSelectionFatJet(lqdir + "/CATConfig/SelectionConfig/user_fatjets.sel");

  IDSetup=true;
  if(k_debugmode){
    for( map<TString,vector<pair<TString,float> > >::iterator it=  selectionIDMapfMuon.begin() ; it !=  selectionIDMapfMuon.end(); it++){
      cout << it->first << endl;
      for (unsigned int i=0 ; i < it->second.size(); i++){
        cout << it->second.at(i).first << " " << it->second.at(i).second << endl;
      }
    }
  }
  Message("SetupSelection DONE", DEBUG);

}

void AnalyzerCore::SetupDDBkg(){
  if(k_running_nonprompt || k_running_chargeflip)m_datadriven_bkg = new DataDrivenBackgrounds();
  setupDDBkg=true;
  
  /// setup correction class at teh same time
  
  /// not needed for sktreemaker
  // save time as code does not need to setup and save files
  
  string lqdir =  getenv("LQANALYZER_DIR");

  // List of working points                                                                                                                                                                                                                                                  

  if(k_classname == "SKTreeMaker")  mcdata_correction = new MCDataCorrections();

  if(!k_classname.Contains("SKTreeMaker")){

    mcdata_correction = new MCDataCorrections();
    
    std::vector<TString> vtaggers;
    vtaggers.push_back("CSVv2Moriond17_2017_1_26");
    vtaggers.push_back("cMVAv2Moriond17_2017_1_26");

    std::vector<TString> v_wps;
    v_wps.push_back("Loose");
    v_wps.push_back("Medium");
    v_wps.push_back("Tight");
    MapBTagSF = SetupBTagger(vtaggers,v_wps);
    
    if(1){
      ifstream runlumi((lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt").c_str());
      if(!runlumi) {
	cerr << "Did not find "+lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt'), exiting ..." << endl;
	exit(EXIT_FAILURE);
      }
      string lline;
      int x=1;
      while(getline(runlumi,lline) ){
	std::istringstream is( lline );

	string trigname;
	float trig_lumi;
	int run;
	is >> trigname;
	if(trigname=="###" ) continue;
	if(trigname=="END") break;
	if(trigname=="run" ){
	  is >> run;
	  is >> trig_lumi;
	  if(k_debugmode)
	    cout << "Run number "<< run <<" ; Muon trigger (unprescaled) luminosity  " << trig_lumi << ";" << endl;

	  mapLumi2016[run] = trig_lumi;
	  continue;
	}
	if(trigname=="block" ){
	  is >> run;
	  is >> trig_lumi;
	  if(k_debugmode)cout << "mapLumi[" << run <<" ] = " << trig_lumi << ";" << endl;

	  mapLumiPerBlock2016[run] = trig_lumi;
	  ostringstream ss;
	  ss << x;
	  mapLumiNamePerBlock2016[run]="Lumi"+ss.str();
	  x++;
	  continue;
	}

      }
      runlumi.close();
    }

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    
    if(k_debugmode)cout << "Reading Luminosity File" << endl;

    SetupLuminosityMap(true);


  }

}

void AnalyzerCore::SetUpEvent(Long64_t entry, float ev_weight) throw( LQError ) {
  

  if(!IDSetup)   SetupID();
  if(!setupDDBkg)SetupDDBkg();
  

  Message("In SetUpEvent(Long64_t entry) " , DEBUG);
  m_logger << DEBUG << "This is entry " << entry << LQLogger::endmsg;
  if (!fChain) throw LQError( "Chain is not initialized",  LQError::SkipCycle );     
  
  if(LQinput){

    m_logger << DEBUG << "k_isdata = " << k_isdata << " and isData = " << isData << LQLogger::endmsg;
    if(k_isdata != isData) throw LQError( "!!! Event is confused. It does not know if it is data or MC", LQError::SkipCycle );
  }
  else isData = k_isdata;
  
  if (!(entry % output_interval)) {
    m_logger << INFO <<  "Processing entry " << entry <<  "/" << nentries << LQLogger::endmsg;

  }
  
  snu::KEvent eventinfo = GetEventInfo();
  
  if(k_isdata){
    if(ev_weight!=1.) Message("ERROR in setting weights. This is Data...", INFO);
    MCweight=1.;
    weight = 1.;
  }
  else {
    MCweight = eventinfo.MCWeight(); 
    weight= ev_weight; 
  }
 //
  // creates object that stores all SKTree classes	
  //                                                                                                        
  
  snu::KTrigger triggerinfo = GetTriggerInfo(triggerlist);
    
  std::vector<snu::KJet> skjets= GetAllJets();
  std::vector<snu::KFatJet> skfatjets;
  std::vector<snu::KGenJet> skgenjets;
  if(k_usegenjet) skgenjets=GetAllGenJets();
  if(k_usefatjet) skfatjets= GetAllFatJets();

  /// np == numberof particles you want to store at truth info. 30 is default unless running nocut sktree OR signal
  int np =  AssignnNumberOfTruth();
  
  //// Tmp vect to speed jobs when not used
  
  vector<snu::KPhoton> photons ;
  if(k_usephotons) photons =  GetAllPhotons();
  vector<snu::KTruth> gen;
  if(k_usetruth) gen=GetTruthParticles(np);
  
  vector<snu::KMuon>  muons = GetAllMuons();

  if(!k_classname.Contains("SKTreeMaker")){
    if(k_skim="FLATCAT"){
      SetCorrectedMomentum(muons, gen);
    }
  }
  LQEvent lqevent(muons, GetAllElectrons(), photons, skjets,skfatjets, skgenjets, gen, triggerinfo,eventinfo);
  

  //  eventbase is master class to use in analysis 
  //
  
  eventbase = new EventBase(lqevent);
  
  /*eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron);
    eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);*/
  eventbase->GetJetSel()->SetIDSMap(selectionIDMapsJet);
  eventbase->GetJetSel()->SetIDFMap(selectionIDMapfJet);
  eventbase->GetFatJetSel()->SetIDSMap(selectionIDMapsFatJet);
  eventbase->GetFatJetSel()->SetIDFMap(selectionIDMapfFatJet);

  bool setupFullIDinSelection=false;
  /// setting this to true makes the code run slower

  if(k_running_nonprompt || k_running_chargeflip || setupFullIDinSelection){
    //m_datadriven_bkg needs ID maps to get isTight for IDs 
    eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron); 
    eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);
    eventbase->GetMuonSel()->SetIDSMap(selectionIDMapsMuon);
    eventbase->GetMuonSel()->SetIDFMap(selectionIDMapfMuon);
    m_datadriven_bkg->SetEventBase(eventbase);
  }
  
  if(!k_isdata){
    if(!changed_target_lumi){
      changed_target_lumi=true;
    }
  }
  

  /// Setup correction class
  k_reset_period=true;
  

  if(!k_classname.Contains("SKTreeMaker")){
    mcdata_correction->SetPeriod(GetPeriod());
    mcdata_correction->SetIsData(isData);
  }
  if (k_classname == "SKTreeMaker"){
    mcdata_correction->SetPeriod(GetPeriod());
    mcdata_correction->SetIsData(isData);
    
  }
  
  
}


int AnalyzerCore::VersionStamp(TString cversion){
  
  if(cversion.Contains("v7-4-4")) return 1;
  else if(cversion.Contains("v7-4-5")) return 2;
  else if(cversion.Contains("v7-6-2") || cversion.Contains("v7-6-3") || cversion.Contains("v7-6-4")   ) return 3;
  else if((cversion.Contains("v7-6-5") || cversion.Contains("v7-6-6"))) return 4;
  else if((cversion.Contains("v8-0-1"))) return 5;
  else if((cversion.Contains("v8-0-2"))) return 6;
  else if((cversion.Contains("v8-0-3"))) return 7;
  else if((cversion.Contains("v8-0-4") || cversion.Contains("v8-0-5"))) return 8;
  
  return 5;
 
}


snu::KTruth AnalyzerCore::GetTruthMatchedParticle(snu::KElectron el){
  
  if(el.MCTruthIndex() >  eventbase->GetTruth().size()){
    snu::KTruth tr;
    return tr;
  }
  return   eventbase->GetTruth()[el.MCTruthIndex()];
}


int AnalyzerCore::AssignnNumberOfTruth(){
  int np = 1000;
  if(k_classname.Contains("SKTreeMaker")) np = 1000;
  if(k_classname.Contains("SKTreeMakerDiLep")) np = 0;
  if(k_classname.Contains("SKTreeMakerTriLep")) np = 1000;

  if(k_classname.Contains("SKTreeMaker")){
    if(k_sample_name.Contains("QCD") && !k_sample_name.Contains("mad")) np = 0;
  }

  /// List of signal samples
  /// G.Yu needs to add signal here
  
  if(IsSignal()) np = 1000;
  
  return np;
}



bool AnalyzerCore::IsSignal(){
  
  if(k_sample_name.Contains("Majornana")) return true;
  if(k_sample_name.Contains("HN")) return true;
  return false;
}


void AnalyzerCore::ClassInfo(){
  
  /*if(eventinfo.CatVersion().empty()){ 
    m_logger << INFO << "Catuple version is v7-4-X. Only basic infomation is available." << LQLogger::endmsg;
    
  }
    
  else if(TString(eventinfo.CatVersion()).Contains("v7-6-2")){
    m_logger << INFO <<  "Running on catuples version " << eventinfo.CatVersion() << LQLogger::endmsg;
    
    
    }*/

}

float AnalyzerCore::SumPt( std::vector<snu::KJet> particles){

  float sumpt=0.;
  for(std::vector<snu::KJet>::iterator it = particles.begin(); it != particles.end(); it++){
    sumpt += it->Pt();
  }
  return sumpt;
}


float AnalyzerCore::SumPt( std::vector<snu::KFatJet> particles){

  float sumpt=0.;
  for(std::vector<snu::KFatJet>::iterator it = particles.begin(); it != particles.end(); it++){
    sumpt += it->Pt();
  }
  return sumpt;
}

  
bool AnalyzerCore::IsDiEl(){
  if(isData) return false;
  int iel(0);
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11) iel++;
  }
  if(iel >1) return true;
  else return false;
}

bool AnalyzerCore::ISCF(snu::KElectron el){
  
  if(el.GetType() == 4) return true;
  if(el.GetType() == 5)return true;
  if(el.GetType() == 6)return true;
  if(el.GetType() == 13&&el.MCMatched())return true;
  if(el.GetType() == 19)return true;
  if(el.GetType() == 20)return true;
  if(el.GetType() == 21)return true;
  return false;
}

bool AnalyzerCore::TruthMatched(snu::KElectron el, bool keepCF){
  bool pass=false;
  if(!keepCF && ISCF(el)) return false;
  if(keepCF && ISCF(el)) return true;
  if((keepCF && !ISCF(el)) || !keepCF) {
    
    if(el.GetType() ==1)   pass=true; /// Z/W
    if(el.GetType() ==2)   pass=true; /// Z/W
    if(el.GetType() ==3)   pass=true; /// Z/W 
    if(el.GetType() == 11) pass=true; /// Tau
    if(el.GetType() == 14) pass=true; /// Z*
    if(el.GetType() == 15) pass=true; /// W*
    if(el.GetType() == 17) pass=true; /// * CF
    if(el.GetType() == 18) pass=true; /// * CF  
    //if(el.GetType() == 23) pass=true;
    if(el.GetType() == 35) pass=true;
    if(el.GetType() == 40) pass=true;
    if(el.GetType() == 16) pass=true;
  }

  return pass;
}

int AnalyzerCore::IsFakeEvent(vector<snu::KElectron> els ){
  
  if( ((k_sample_name.Contains("TT"))&&SameCharge(els))) return 2;
  if( ((k_sample_name.Contains("DY"))&&SameCharge(els))) return 2;
  if( ((k_sample_name.Contains("WJets"))&&els.size()==2)) return 2;
  if( ((k_sample_name.Contains("qcd"))&&els.size()>0)) return els.size();
  if( ((k_sample_name.Contains("QCD"))&&els.size()>1)) return els.size();
  
  return -1;

}

int  AnalyzerCore::IsFakeEvent(vector<snu::KMuon> mus ){

  if( ((k_sample_name.Contains("TT"))&&SameCharge(mus))) return 2;
  if( ((k_sample_name.Contains("DY"))&&SameCharge(mus))) return 2;
  if( ((k_sample_name.Contains("WJets"))&&mus.size()==2)) return 2;
  if( ((k_sample_name.Contains("qcd"))&&mus.size()>0)) return mus.size();
  if( ((k_sample_name.Contains("QCD"))&&mus.size()>1)) return mus.size();

  return -1;



}


bool AnalyzerCore::NonPrompt(snu::KElectron el){
  
  if(el.GetType() == 7) return true;
  return false;

}
bool AnalyzerCore::NonPrompt(snu::KMuon mu){

  if(mu.GetType() == 2) return true;
  return false;

}



bool AnalyzerCore::TruthMatched(std::vector<snu::KElectron> el, bool tightdxy, bool allowCF){
  
  bool pass=false;
  
  return pass;
}

bool AnalyzerCore::TruthMatched(snu::KMuon mu){

  bool pass=false;
  
  if(mu.GetType() ==1) pass=true;
  if(mu.GetType() ==6) pass=true;
  if(mu.GetType() ==8) pass=true;
  if(mu.GetType() ==9) pass=true;
  if(mu.GetType() ==12) pass=true;
  if(mu.GetType() ==25) pass=true;  //// HAS CHANGED AFTER SKTREE were remade
  //if(mu.GetType() ==28) pass=true;
  // 23>?? 28?? sshould these be included?
  return pass;
}

float AnalyzerCore::GetVirtualMass(int pdg, bool includenu, bool includeph){
  if(isData) return -999.;
  vector<KTruth> es1;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdg){
      if(eventbase->GetTruth().at(ig).GenStatus() ==1){
        es1.push_back(eventbase->GetTruth().at(ig));
      }
    }
    else   if(includenu){
      if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdg+1){
	if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	  es1.push_back(eventbase->GetTruth().at(ig));
	}
      }
    }

    else if(includeph){
      if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	es1.push_back(eventbase->GetTruth().at(ig));
      }
    }
  }

  if(!includeph){
    if(!includenu){
      if(es1.size()==2){
	snu::KParticle ll = es1[0]  + es1[1];
	return ll.M();
      }
    }
    else  if(es1.size()==2){
      snu::KParticle ll = es1[0]  + es1[1];
      return ll.M();
    }
  }
  else{
    snu::KParticle ll;
    for(unsigned int iel=0; iel < es1.size(); iel++){
      ll+= es1[iel];
    }
    return ll.M();
  }
  
  return -999.;

}




float AnalyzerCore::GetVirtualMassConv(int cmindex,int nconvindx){

  if(isData) return -999.;
  cout << "cmindex = " << cmindex << endl;
  vector<KTruth> es1;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    
    if(eventbase->GetTruth().at(ig).IndexMother() == eventbase->GetTruth().at(cmindex).IndexMother()){
      if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11){
	if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	  es1.push_back(eventbase->GetTruth().at(ig));
	}
      }
    }
  }
  
  if(es1.size()==3){
    if(nconvindx==0){
      snu::KParticle ll = es1[0]  + es1[1];
      return ll.M();
    }
    if(nconvindx==1){
      snu::KParticle ll = es1[0]  + es1[2];
      return ll.M();
    }
    if(nconvindx==2){
      snu::KParticle ll = es1[1]  + es1[2];
      return ll.M();
    }
  }
  return -999.;
}



void AnalyzerCore::TruthPrintOut(){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;



  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    //if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |       |         |" << endl;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size())){
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  ---  |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |  --- |" << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;

    }
    else{
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |   " << eventbase->GetTruth().at(ig).IndexMother()  << " " << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl; 
    }
  }

}


bool AnalyzerCore::isPrompt(long pdgid) {
  /// mother pdgid
  pdgid = abs(pdgid);
  if (pdgid == 24) return true; // Z
  else if (pdgid == 23) return true; // W
  else if (pdgid == 15) return true; // taus
  else if (pdgid == 90) return true; // N
  else return false;
}

void AnalyzerCore::EndEvent()throw( LQError ){
  
  delete eventbase;                                                                                                            

}
  

void AnalyzerCore::ListTriggersAvailable(){
  cout << "Set of triggers you can use are: " << endl;
  for(unsigned int i=0; i < eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().size(); i++){
    cout << eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().at(i)<< " has prescale " << eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales().at(i)<< endl;
  }
  return;
}


bool AnalyzerCore::PassJets(std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets ){

  if(jets.size() < 2 && fatjets.size() == 0) return false;
  
  return true;
}

//################################################################################################
//@@@@@@@@@@@@@@@@@@@  TRIGGER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     

bool AnalyzerCore::PassTrigger(vector<pair<TString,TString> > list){
  
  vector<TString> triglist;
  for(vector<pair<TString,TString>  >::iterator it = list.begin(); it != list.end(); it++){
    if(k_channel.Contains(it->second)){
      for(unsigned int i=0; i <  triglist.size(); i++){
	if(PassTrigger(it->first)) return false;
      }
      if(PassTrigger(it->first)) return true;
      triglist.push_back(it->first);
    }
  }
  return false;
}



bool AnalyzerCore::PassTrigger(TString trig){
  vector<TString> list;
  list.push_back(trig);
  int  prescaler=1.;
  return TriggerSelector(list, eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerDecisions(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales(), prescaler);
 
 
}


bool AnalyzerCore::PassTriggerOR(vector<TString> list){

  int  prescaler=1.;
  return TriggerSelector(list, eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerDecisions(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales(), prescaler);
  
}


////###############################################################################################
/// @@@@@@@@@@@@@@@@@@@@@@@@@ MISC   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TDirectory* AnalyzerCore::GetTemporaryDirectory(void) const
{
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:                                              
    std::stringstream dirname;
    dirname << "CATAnalzer_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                                                       
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}


std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KElectron> el){return mcdata_correction->MakeTLorentz(el);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KMuon> mu){return mcdata_correction->MakeTLorentz(mu);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KJet> jet){return mcdata_correction->MakeTLorentz(jet);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KFatJet> jet){return mcdata_correction->MakeTLorentz(jet);}

/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2  FUNCTIONS USED BY TRILEPTON ANALYSIS: jskim

//==== HN Trilepton stuffs 

void AnalyzerCore::PutNuPz(TLorentzVector *nu, double Pz){
  double Px, Py;
  Px = nu->Px();
  Py = nu->Py();
  nu->SetPxPyPzE(Px, Py, Pz, TMath::Sqrt(Px*Px+Py*Py+Pz*Pz));
}

void AnalyzerCore::PutNuPz(snu::KParticle *nu, double Pz){
  double Px, Py;
  Px = nu->Px();
  Py = nu->Py();
  nu->SetPxPyPzE(Px, Py, Pz, TMath::Sqrt(Px*Px+Py*Py+Pz*Pz));
}

double AnalyzerCore::solveqdeq(double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm){
  TLorentzVector met;
  met.SetPxPyPzE(MET*cos(METphi),
		 MET*sin(METphi),
		 0,
		 MET);

  Double_t d = (W_mass*W_mass)-(l1l2l3.M())*(l1l2l3.M())+2.0*l1l2l3.Px()*met.Px()+2.0*l1l2l3.Py()*met.Py();
  Double_t a = l1l2l3.E()*l1l2l3.E() - l1l2l3.Pz()*l1l2l3.Pz();
  Double_t b = d*l1l2l3.Pz();
  Double_t c = l1l2l3.E()*l1l2l3.E()*met.E()*met.E()-d*d/4.0;
  if(b*b-4*a*c<0){
    return b/(2*a);
  }
  else{
    if(pm=="p") return (b+TMath::Sqrt(b*b-4*a*c))/(2*a);
    else if(pm=="m")  return (b-TMath::Sqrt(b*b-4*a*c))/(2*a);
    else return 0;
  }
}

int AnalyzerCore::find_mlmet_closest_to_W(snu::KParticle lep[], snu::KParticle MET, int n_lep){
  double m_diff[n_lep];
  double m_diff_min = 999999999.;
  int outindex = 0;
  for(int i=0; i<n_lep; i++){
    double dphi = lep[i].DeltaPhi(MET);
    double mt2 = 2.*lep[i].Pt()*MET.Pt()*(1.-TMath::Cos(dphi));
    m_diff[i] = fabs( sqrt(mt2) - 80.385 );
    if( m_diff[i] < m_diff_min ){
      m_diff_min = m_diff[i];
      outindex = i;
    }
  }
  return outindex;
  
}

int AnalyzerCore::find_mlmet_closest_to_W(std::vector<KLepton> lep, snu::KParticle  MET){

  int n_lep = lep.size();

  double m_diff[n_lep];
  double m_diff_min = 999999999.;
  int outindex = 0;
  for(int i=0; i<n_lep; i++){
    double dphi = lep[i].DeltaPhi(MET);
    double mt2 = 2.*lep[i].Pt()*MET.Pt()*(1.-TMath::Cos(dphi));
    m_diff[i] = fabs( sqrt(mt2) - 80.385 );
    if( m_diff[i] < m_diff_min ){
      m_diff_min = m_diff[i];
      outindex = i;
    }
  }
  return outindex;

}

double AnalyzerCore::MT(TLorentzVector a, TLorentzVector b){
  double dphi = a.DeltaPhi(b);
  return TMath::Sqrt( 2.*a.Pt()*b.Pt()*(1.- TMath::Cos(dphi) ) );
  
}

bool AnalyzerCore::GenMatching(snu::KParticle gen, snu::KParticle reco, double maxDeltaR, double maxPtDiff){
  bool matched = true;
  if(reco.Pt() > 200.) maxPtDiff = 0.1;
  
  if( gen.DeltaR(reco) >= maxDeltaR ) matched = false;
  if( fabs(gen.Pt() - reco.Pt()) / reco.Pt() >= maxPtDiff ) matched = false;
  
  return matched;
  
}

std::vector<snu::KMuon> AnalyzerCore::GetHNTriMuonsByLooseRelIso(double LooseRelIsoMax, bool keepfake, float ptcut, float etacut){
 
  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons("MUON_HN_TRI_VLOOSE", keepfake, ptcut, etacut);
  std::vector<snu::KMuon> muontriLooseColl;
  muontriLooseColl.clear();
  for(unsigned int j=0; j<muontriLooseColl_raw.size(); j++){
    snu::KMuon this_muon = muontriLooseColl_raw.at(j);
    if( this_muon.RelIso04() < LooseRelIsoMax ) muontriLooseColl.push_back( this_muon );
  }
  return muontriLooseColl;
  
}

void AnalyzerCore::PrintTruth(){
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  cout << "=========================================================" << endl;
  cout << "truth size = " << truthColl.size() << endl;
  cout << "index" << '\t' << "pdgid" << '\t' << "status" << '\t' << "mother" << '\t' << "mother pid" << "\t" << "pt" << "\t" << "eta" << "\t" << "mass" << endl;
  for(int i=2; i<truthColl.size(); i++){
    cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).GenStatus() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << "\t" << truthColl.at(i).Pt() << "\t" << truthColl.at(i).Eta() << "\t" << truthColl.at(i).M() << endl;
  }

}


void AnalyzerCore::Message(TString message, LQMsgType type){
  m_logger <<  type << message << LQLogger::endmsg;
}



////###############################################################################################                                                                           
 /// @@@@@@@@@@@@@@@@@@@@@@@@@ HIST   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                                                             

void AnalyzerCore::MakeCleverHistograms(histtype type, TString clhistname ){
  
  //// ELECTRON PLOTs                                                                                          
  if(type==elhist) mapCLhistEl[clhistname] = new ElectronPlots(clhistname);
  //// BaseSelection::MUON PLOTs                                                                                              
  if(type==muhist) mapCLhistMu[clhistname] = new MuonPlots(clhistname);
  /// JET PLOTs                                                                                                
  if(type==jethist) mapCLhistJet[clhistname] = new JetPlots(clhistname);
  /// Signal plots                                                                                             
  if(type==sighist_e)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 1);
  if(type==sighist_ee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 2);
  if(type==sighist_eee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 3);
  if(type==sighist_eeee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname,-1);
  if(type==sighist_m)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,1);
  if(type==sighist_mm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,2);
  if(type==sighist_mmm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,3);
  if(type==sighist_mmmm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,-1);
  if(type==sighist_em)  mapCLhistSigEM[clhistname] = new SignalPlotsEM(clhistname);

  if(type==trilephist)  mapCLhistTriLep[clhistname] = new TriLeptonPlots(clhistname);
  if(type==hnpairmm) mapCLhistHNpairMM[clhistname] = new HNpairPlotsMM(clhistname);
  if(type==hntrilephist)  mapCLhistHNTriLep[clhistname] = new HNTriLeptonPlots(clhistname);
    
  return;
}

void AnalyzerCore::MakeHistograms(){
  //// Additional plots to make                                                                                
  maphist.clear();
  maphist2D.clear();
  maphist3D.clear();

    
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xbins[], TString label){
  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xbins);
  maphist[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xmin, float xmax, TString label){

  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xmin,xmax);
  //maphist[hname]->GetXaxis()->SetTitle("TEST");
  maphist[hname]->GetXaxis()->SetTitle(label);
}


void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
}


void AnalyzerCore::MakeHistograms3D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, int nbinsz, float zmin, float zmax, TString label) {
  
  maphist3D[hname] =  new TH3D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax, nbinsz, zmin, zmax);
  maphist3D[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx,  float xbins[], int nbinsy,  float ybins[], TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx , xbins, nbinsy,ybins);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
}


void AnalyzerCore::FillHist(TString histname, float value, float w, float xbins[], int nbins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;

  if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xbins, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }

}


void AnalyzerCore::FillHistPerLumi(TString histname, float value, float w, float xmin, float xmax,int nbins, int nlumibins){
  
  if(!Is2015Analysis()){
    if(nlumibins==10){
      
      if(!GetHist(histname+"_perlumi")) {
	MakeHistograms(histname+"_perlumi", 9, 0., 9.);
	int nbin=0;

	for(std::map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	  nbin++;
	  GetHist(histname+"_perlumi")->GetXaxis()->SetBinLabel(nbin,it->second);
	}
      }
      
      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	
        if(!GetHist(histname+"_"+it->second)) {
          MakeHistograms(histname+"_"+it->second, nbins, xmin, xmax);
	}
      } 


      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	if(eventbase->GetEvent().RunNumber()  < it->first) {
	  map<int,float>::iterator it2 = mapLumiPerBlock.find(it->first);
	  
	  if(isData){
	    float neww= w /it2->second;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	  else{
	    //	    float neww = w * (it2->second/TargetLumi);
	    float neww = w /TargetLumi;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	} 
      }
    }// nbins
  }

}

void AnalyzerCore::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins , TString label){
  
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist(histname)) GetHist(histname)->Fill(value, w);  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xmin, xmax, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }
  
}

void AnalyzerCore::FillHist(TString histname, float value1, float value2, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , TString label){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  }

}

void AnalyzerCore::FillHist(TString histname, float valuex, float valuey, float w, float xbins[], int nxbins, float ybins[], int nybins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex,valuey, w);

  else{
    if (nxbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nxbins, xbins, nybins, ybins , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex, valuey, w);
  }

}


void AnalyzerCore::FillHist(TString histname, float value1, float value2,float value3, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , float zmin, float zmax, int nbinxz, TString label){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist3D(histname)) GetHist3D(histname)->Fill(value1,value2,value3, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms3D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , nbinxz, zmin, zmax,label);
    if(GetHist3D(histname)) GetHist3D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist3D(histname)) GetHist3D(histname)->Fill(value1,value2,value3, w);
  }

}


void AnalyzerCore::FillHist(TString histname, float value, float w , TString label){

  if(GetHist(histname)){
    GetHist(histname)->Fill(value, w);  /// Plots Z peak        
    GetHist(histname)->GetXaxis()->SetTitle(label);
  }
  
  else m_logger << INFO << histname << " was NOT found. Will add the histogram to the hist map on first event." << LQLogger::endmsg;
  
  
  return;
}

void AnalyzerCore::FillUpDownHist(TString histname, float value, float w, float w_err, float xmin, float xmax, int nbins , TString label){

  FillHist(histname+"_up", value, w+w_err, xmin, xmax, nbins, label);
  FillHist(histname+"_down", value, w-w_err, xmin, xmax, nbins, label);
  FillHist(histname, value, w, xmin, xmax, nbins, label);

}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons, double w){

  if(type==muhist){
    map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.find(hist);
    if(mupit != mapCLhistMu.end()) mupit->second->Fill(w,muons);
    else m_logger << INFO  << hist << " not found in mapCLhistMu" << LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to muhist, is this a mistake?" << LQLogger::endmsg;

}


void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons,double w){

  if(type==elhist){
    map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.find(hist);
    if(elpit !=mapCLhistEl.end()) elpit->second->Fill(w,electrons);
    else m_logger << INFO  << hist << " not found in mapCLhistEl" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to elhist, is this a mistake?" << LQLogger::endmsg;
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KJet> jets, double w){
  
  if(type==jethist){
    map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.find(hist);
    if(jetpit !=mapCLhistJet.end()) jetpit->second->Fill(w,jets);
    else m_logger << INFO  << hist << " not found in mapCLhistJet" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  <<"Type not set to jethist, is this a mistake?" << LQLogger::endmsg;

}
void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double w, int nbjet){
  if(type==hnpairmm){
    map<TString, HNpairPlotsMM*>::iterator HNpairmmit = mapCLhistHNpairMM.find(hist);
    if(HNpairmmit !=mapCLhistHNpairMM.end()) HNpairmmit->second->Fill(ev, muons, electrons, jets, w, nbjet);
    else {
      mapCLhistHNpairMM[hist] = new HNpairPlotsMM(hist);
      HNpairmmit = mapCLhistHNpairMM.find(hist);
      HNpairmmit->second->Fill(ev, muons, electrons, jets, w, nbjet);
    }
  }
}


void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double w){
  
  vector<snu::KFatJet> fatjets;
  vector<snu::KJet> alljets;
  FillCLHist(type, hist, ev, muons,electrons ,jets, alljets,fatjets,w);
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,vector<snu::KJet> alljets,double w){

  vector<snu::KFatJet> fatjets;
  FillCLHist(type, hist, ev, muons,electrons ,jets, alljets, fatjets,w);
}


void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets, vector<snu::KJet> alljets, vector<snu::KFatJet> fatjets,double w){

  if(type==trilephist){

    map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.find(hist);
    if(trilepit !=mapCLhistTriLep.end()) trilepit->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistTriLep[hist] = new TriLeptonPlots(hist);
      trilepit = mapCLhistTriLep.find(hist);
      trilepit->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==hntrilephist){

    map<TString, HNTriLeptonPlots*>::iterator hntrilepit = mapCLhistHNTriLep.find(hist);
    if(hntrilepit !=mapCLhistHNTriLep.end()) hntrilepit->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistHNTriLep[hist] = new HNTriLeptonPlots(hist);
      hntrilepit = mapCLhistHNTriLep.find(hist);
      hntrilepit->second->Fill(ev, muons, electrons, jets,w);
    }
  }
 
  else if(type==sighist_ee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,  alljets,fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,2);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }
  else if(type==sighist_eee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist, 3);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets,fatjets,w);
    }
  }
  else if(type==sighist_eeee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,-1);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }
  else if(type==sighist_e){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,1);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }


  else if(type==sighist_m){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,1);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);

    }
  }
  else if(type==sighist_mm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,2);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
      
    }
  }
  else if(type==sighist_mmm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,3);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
      
    }
  }
  else if(type==sighist_mmmm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,-1);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
      
    }
  }

  else if(type==sighist_em){

    map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.find(hist);
    if(sigpit_em !=mapCLhistSigEM.end()) sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigEM[hist] = new SignalPlotsEM(hist);
      sigpit_em = mapCLhistSigEM.find(hist);
      sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    }
  }
 else  m_logger << INFO  <<"Type not set to sighist, is this a mistake?" << LQLogger::endmsg;


}


void AnalyzerCore::WriteHistograms() throw (LQError){
  // This function is called after the cycle is ran. It wrues all histograms to the output file. This function is not used by user. But by the contrioller code.
  WriteHists();
  WriteCLHists();
  WriteNtp();
}

  
void AnalyzerCore::WriteCLHists(){

  for(map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.begin(); mupit != mapCLhistMu.end(); mupit++){

    Dir = m_outputFile->mkdir(mupit->first);
    m_outputFile->cd( Dir->GetName() );
    mupit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.begin(); elpit != mapCLhistEl.end(); elpit++)\
    {

      Dir = m_outputFile->mkdir(elpit->first);
      m_outputFile->cd( Dir->GetName() );
      elpit->second->Write();
      m_outputFile->cd();
    }

  for(map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.begin(); jetpit != mapCLhistJet.end(); jetpit++)\
    {
      
      Dir = m_outputFile->mkdir(jetpit->first);
      m_outputFile->cd( Dir->GetName() );
      jetpit->second->Write();
      m_outputFile->cd();
    }
  for(map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.begin(); sigpit_ee != mapCLhistSigEE.end(); sigpit_ee++){
    
    Dir = m_outputFile->mkdir(sigpit_ee->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_ee->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.begin(); sigpit_mm != mapCLhistSigMM.end(); sigpit_mm++){

    Dir = m_outputFile->mkdir(sigpit_mm->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_mm->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.begin(); sigpit_em != mapCLhistSigEM.end(); sigpit_em++){

    Dir = m_outputFile->mkdir(sigpit_em->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_em->second->Write();
    m_outputFile->cd();
  }
  
  for(map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.begin(); trilepit != mapCLhistTriLep.end(); trilepit++){

    Dir = m_outputFile->mkdir(trilepit->first);
    m_outputFile->cd( Dir->GetName() );
    trilepit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, HNpairPlotsMM*>::iterator HNpairmmit = mapCLhistHNpairMM.begin(); HNpairmmit != mapCLhistHNpairMM.end(); HNpairmmit++){
    
    Dir = m_outputFile->mkdir(HNpairmmit->first);
    m_outputFile->cd( Dir->GetName() );
    HNpairmmit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, HNTriLeptonPlots*>::iterator hntrilepit = mapCLhistHNTriLep.begin(); hntrilepit != mapCLhistHNTriLep.end(); hntrilepit++){

    //==== (jskim)I don't need director!
    //Dir = m_outputFile->mkdir(hntrilepit->first);
    //m_outputFile->cd( Dir->GetName() );
    hntrilepit->second->Write();
    m_outputFile->cd();
  }



  return;
}

void AnalyzerCore::WriteHists(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    
    
    
    if(mapit->first.Contains("closejet")){
      if(!m_outputFile->GetDirectory( "closejet" )){
	Dir = m_outputFile->mkdir("closejet");
	m_outputFile->cd( Dir->GetName() );
      }
      else  m_outputFile->cd("closejet");
      mapit->second->Write();
      m_outputFile->cd();
    }

    
    
    else {
      mapit->second->Write();
    }
  }
  
  for(map<TString, TH2*>::iterator mapit = maphist2D.begin(); mapit != maphist2D.end(); mapit++){
    mapit->second->Write();
  }
  for(map<TString, TH3*>::iterator mapit = maphist3D.begin(); mapit != maphist3D.end(); mapit++){
    mapit->second->Write();
  }

  //==== HN Gen Matching
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    m_HNgenmatch->WriteHNGenHists();
  }
  return;
}

TH1* AnalyzerCore::GetHist(TString hname){

  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}



TH2* AnalyzerCore::GetHist2D(TString hname){

  TH2* h = NULL;
  std::map<TString, TH2*>::iterator mapit = maphist2D.find(hname);
  if(mapit != maphist2D.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}


TH3* AnalyzerCore::GetHist3D(TString hname){

  TH3* h = NULL;
  std::map<TString, TH3*>::iterator mapit = maphist3D.find(hname);
  if(mapit != maphist3D.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}




void AnalyzerCore::FillCutFlow(TString cut, float weight){

  bool cut_exists=false;
  m_logger << DEBUG  << "FillCutFlow" <<  LQLogger::endmsg;
  for(unsigned int i=0; i < cutflow_list.size();  i++){
    if (cut == cutflow_list.at(i)) cut_exists=true;
  }

  if (!cut_exists){
    m_logger << DEBUG  << "FillCutFlow: " << cut <<  LQLogger::endmsg;

    n_cutflowcuts=n_cutflowcuts+1;
    if(GetHist("cutflow")) {
      m_logger << DEBUG  << "FillCutFlow: mod " << cut <<  LQLogger::endmsg;

      map<TString, TH1*>::iterator it = maphist.find("cutflow");
      vector<float> counters;
      for(int j = 0; j  < n_cutflowcuts ; j++){
	counters.push_back(it->second->GetBinContent(1+j));
      }
      delete it->second;
      
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      for(unsigned int i=0;i < cutflow_list.size();  i++){
	GetHist("cutflow")->GetXaxis()->SetBinLabel(i+1,cutflow_list.at(i));
	GetHist("cutflow")->Fill(cutflow_list.at(i), weight);
      }
      GetHist("cutflow")->GetXaxis()->SetBinLabel(n_cutflowcuts,cut);
      cutflow_list.push_back(cut);
      
    }
    else{
      m_logger << DEBUG  << "FillCutFlow:Fill " << cut <<  LQLogger::endmsg;
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      GetHist("cutflow")->GetXaxis()->SetBinLabel(1,cut);
      GetHist("cutflow")->Fill(cut,weight);
      cutflow_list.push_back(cut);
    }
  }
  else {
    if(GetHist("cutflow")) {
      GetHist("cutflow")->Fill(cut,weight);
    }
  }
}

////###############################################################################################
/// @@@@@@@@@@@@@@@@@@@@@@@@@ JETPILEUP   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



////############################################################################################### 
/// @@@@@@@@@@@@@@@@@@@@@@@@@ CUTS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

bool AnalyzerCore::PassMETFilter(){

  bool pass (true);
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_8011_ICHEP_dataset

  if (!eventbase->GetEvent().PassTightHalo2016Filter()) {
    pass = false;    m_logger << DEBUG << "Event Fails PassTightHalo2016Filter " << LQLogger::endmsg;
  }
  if (!eventbase->GetEvent().PassHBHENoiseFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassHBHENoiseFilter " << LQLogger::endmsg;
  }
  if (!eventbase->GetEvent().PassHBHENoiseIsoFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassHBHENoiseIsoFilter " << LQLogger::endmsg;
  }
  if(!eventbase->GetEvent().PassEcalDeadCellTriggerPrimitiveFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassEcalDeadCellTriggerPrimitiveFilter" << LQLogger::endmsg;
  }
  
  if(!eventbase->GetEvent().PassBadChargedCandidateFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassBadChargedCandidateFilterr" << LQLogger::endmsg;
  }

  if(!eventbase->GetEvent().PassBadPFMuonFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails  PassBadPFMuonFilter" << LQLogger::endmsg;
  }

  if (isData){
    if(!eventbase->GetEvent().PassBadEESupercrystalFilter()) {
      pass = false;
      m_logger << DEBUG << "Event FailsPassBadEESupercrystalFilter" << LQLogger::endmsg;
    }
  }
  

  return pass;
}



bool AnalyzerCore::Zcandidate(std::vector<snu::KMuon> muons, float interval, bool require_os){
  
  if(muons.size()!=2) return false;
  if(require_os&&SameCharge(muons)) return false;
  
  KParticle Z = muons.at(0) + muons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;
  
}
  
bool AnalyzerCore::SameCharge(std::vector<snu::KMuon> muons){
  
  if(muons.size()!=2) return false;
  if(muons.at(0).Charge() == muons.at(1).Charge()) return true;
  return false;
}


bool AnalyzerCore::Zcandidate(std::vector<snu::KElectron> electrons, float interval, bool require_os){

  if(electrons.size()!=2) return false;
  if(require_os&&SameCharge(electrons)) return false;

  KParticle Z = electrons.at(0) + electrons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;

}

bool AnalyzerCore::SameCharge(std::vector<snu::KElectron> electrons, bool runningcf){

  if(electrons.size()!=2) return false;


  if(!runningcf){
    if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;

  return false;
}


bool AnalyzerCore::OppositeCharge(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons){
  
  if(electrons.size() != 1) return false;
  if(muons.size() != 1) return false;
  
  if(electrons[0].Charge() == muons[0].Charge()) return false;
  return true;
}

bool AnalyzerCore::OppositeCharge(std::vector<snu::KElectron> electrons, bool runningcf){

  if(electrons.size()!=2) return false;

  if(!runningcf){
    if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;

  return false;
}


int AnalyzerCore::NBJet(std::vector<snu::KJet> jets,  KJet::Tagger tag, KJet::WORKING_POINT wp, int period){

  int nbjet=0;

  if(period == 0) {
    Message("period not set in AnalyzerCore::NBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";

  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26";
  
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";
   
  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);
  
  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis 
  if ( tag == snu::KJet::JETPROB) return -999;
  for(unsigned int ij=0; ij <jets.size(); ij++){

    if(IsBTagged(jets.at(ij), tag, wp)) nbjet++;
    continue;

    bool isBtag=false;
    if (isData) {

      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  -999999, jets.at(ij).Pt(), jets.at(ij).Eta(),period))
	isBtag=true;
    }
    else if (jets.at(ij).HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta(),period))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta(),period))
	isBtag=true;
    }
    
    if(isBtag )nbjet++;
  }
  return nbjet;
}


bool AnalyzerCore::IsBTagged(snu::KJet jet,  KJet::Tagger tag, KJet::WORKING_POINT wp, int mcperiod){

  if(mcperiod == 0) {
    Message("mcperiod not set in AnalyzerCore::IsBTagged. Will assign mcperiod for you but this may not give correct behaviour", WARNING);      
    mcperiod=GetPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26";
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";



  /// Data applied no correction. So only mcperiod is set

  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis
  if ( tag == snu::KJet::JETPROB) return -999;
  
  bool isBtag=false;
  if (isData) {
    
    if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  -999999, jet.Pt(), jet.Eta(), mcperiod))
      isBtag=true;
  }
    else if (jet.HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta(),mcperiod))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta(),mcperiod))
        isBtag=true;
    }
  
  return isBtag;
}


float AnalyzerCore::BTagScaleFactor_1a(std::vector<snu::KJet> jetColl, KJet::Tagger tag, KJet::WORKING_POINT wp, int mcperiod){

  //BTag SF from 1a method.
  //This is coded for H+->WA analysis. I'm fine with anybody else using this function, but be aware that HN analyses decided to use 2a method.
  //And I currently have no plan to use multiple WP. so I just coded to work only for single WP regime.

  if(isData) return 1.;

  if(mcperiod == 0) {
    Message("FYI : mcperiod not set in AnalyzerCore::BTagScaleFactor_1a: meaning auto-set", DEBUG);
    mcperiod=GetPeriod();
  }

  TString btag_key_lf(""), btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)  wp_string = "Loose";
  if(wp == snu::KJet::Medium) wp_string = "Medium";
  if(wp == snu::KJet::Tight)  wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2)  tag_string ="CSVv2Moriond17_2017_1_26";
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";

  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  if ( tag == snu::KJet::JETPROB) return -999;

  float BTagSF=1.;
  for(unsigned int i=0; i<jetColl.size(); i++){
    if(jetColl.at(i).IsBTagged(tag, wp)){
      if(jetColl.at(i).HadronFlavour()==5 || jetColl.at(i).HadronFlavour()==4){
        BTagSF *= it_hf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
      }
      else{
        BTagSF *= it_lf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
      }
    }
    else{
      float SFj=1., Effj=1.;
      if(jetColl.at(i).HadronFlavour()==5 || jetColl.at(i).HadronFlavour()==4){
        SFj  = it_hf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
        Effj = it_hf->second->JetTagEfficiency(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta());
      }
      else{
        SFj  = it_lf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
        Effj = it_lf->second->JetTagEfficiency(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta());
      }

      if( (1.-Effj)==0. ) return 0.;
      BTagSF *= (1.-SFj*Effj)/(1.-Effj);
    }
  }

  return BTagSF;

}



double AnalyzerCore::MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w){
  
  if(mu.size()< 2) return 0.;
  snu::KParticle Z = mu.at(0) + mu.at(1);
  
  double factor (1.);
  if(Z.M() > 90.){
    factor = 8.37401e-01 + 1.61277e-03*Z.M();
  }
  return w*factor;
}


bool AnalyzerCore::IsTight(snu::KMuon muon){
  /// ADD TIGHT BaseSelection::MUON REQUIREMENT
  float reliso= muon.RelIso04();

  if(( reliso >= 0.1)) return false;
  if(fabs(muon.dXY()) >= 0.05) return false; 
  return true;
}

bool AnalyzerCore::PassID(snu::KMuon mu, TString muid){


  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsMuon.find(muid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(muid);
  if(it== selectionIDMapsMuon.end()){
    cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {

    return eventbase->GetMuonSel()->PassUserID(muid,mu, it->second,fit->second);
  }
  return true;
}


bool AnalyzerCore::PassID(snu::KElectron el, TString elid){
  

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(elid);
  if(it== selectionIDMapsElectron.end()){
    cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {

    return eventbase->GetElectronSel()->PassUserID(elid,el, it->second, fit->second);
  }
  return true;
}

bool AnalyzerCore::IsCF(snu::KElectron el){
  vector<snu::KTruth> truth =  eventbase->GetTruth();
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11){
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 23 ||
	 fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 24){
	if(eventbase->GetTruth().at(ig).PdgId() * el.Charge() > 0 ) return true;
	else return false;
      }
    }
  }
  return false;
}

vector<snu::KElectron> AnalyzerCore::GetTruePrompt(vector<snu::KElectron> electrons, bool keep_chargeflip, bool keepfake){
  
  if(electrons.size() == 0)
    return electrons;
  

  vector<snu::KElectron> prompt_electrons;
  for(unsigned int i = 0; i < electrons.size(); i++){

    if(!k_isdata){
      if(k_running_taudecays){
	
	if(electrons.at(i).MCFromTau())  prompt_electrons.push_back(electrons.at(i));
      }
      else{
	
	bool ismatched = TruthMatched(electrons.at(i),  keep_chargeflip);                                                                                                                            
	//electrons.at(i).MCMatched();                                                                                                                                                                         
	//TruthMatched(electrons.at(i),  keep_chargeflip);
	//electrons.at(i).MCMatched();
	if(electrons.at(i).MCFromTau()) ismatched=false;
	if(keepfake&&keep_chargeflip) prompt_electrons.push_back(electrons.at(i));
	else if(keep_chargeflip&& ismatched) prompt_electrons.push_back(electrons.at(i));
	else if(keepfake&&! MCIsCF(electrons.at(i))) prompt_electrons.push_back(electrons.at(i)); 
	else if(ismatched && !MCIsCF(electrons.at(i))) prompt_electrons.push_back(electrons.at(i));
      }
    }// Data
    else prompt_electrons.push_back(electrons.at(i));
  }/// loop

  return prompt_electrons;


}
bool  AnalyzerCore::MCIsCF(snu::KElectron el){
  
  if(el.GetType() == 4) return true;
  if(el.GetType() == 5) return true;
  if(el.GetType() == 6) return true;
  if(el.GetType() == 13) return true;
  if(el.GetType() == 19) return true;
  if(el.GetType() == 20) return true;
  if(el.GetType() == 21) return true;
  
  return false;
}
vector<snu::KMuon> AnalyzerCore::GetTruePrompt(vector<snu::KMuon> muons, bool keepfake){
  if(muons.size()==0)return muons;

  vector<snu::KMuon> prompt_muons;

  for(unsigned int i = 0; i < muons.size(); i++){
    if(!k_isdata){
      
      if(keepfake) prompt_muons.push_back(muons.at(i));
      else if(TruthMatched(muons.at(i))) prompt_muons.push_back(muons.at(i));
    }
    // Data
    else prompt_muons.push_back(muons.at(i));
  }/// loop
  return prompt_muons;
  
}


void AnalyzerCore::CorrectMuonMomentum(vector<snu::KMuon>& k_muons){
  
  mcdata_correction->CorrectMuonMomentum(k_muons,eventbase->GetTruth());
  Message("END CorrectMuonMomentum",DEBUG);    
  return;
    
}

void AnalyzerCore::SetCorrectedMomentum(vector<snu::KMuon>& k_muons, vector<snu::KTruth> truth){

  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){

    if(it->RochPt() < 0.){
      it->SetRochPt(mcdata_correction->GetCorrectedMuonMomentum(*it, truth));
    }
  }

}


void AnalyzerCore::SetCorrectedMomentum(vector<snu::KMuon>& k_muons){
  
  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){
    //if(k_classname=="SKTreeMaker" && (it->RochPt() >0.)) exit(EXIT_FAILURE);
    if(k_classname!="SKTreeMaker" && k_classname.Contains("SKTreeMaker")){
      
    if(it->RochPt() < 0.) {
	cerr << "Roch Pt wrongly set in dilep skim" << endl;
	exit(EXIT_FAILURE);
      }
    }
    
    if(it->RochPt() < 0.){
      it->SetRochPt(mcdata_correction->GetCorrectedMuonMomentum(*it, eventbase->GetTruth()));
    }    
  }
  
}

void AnalyzerCore::MakeNtp(TString hname, TString myvar){

  mapntp[hname] =  new TNtupleD(hname.Data(),hname.Data(),myvar.Data());
}


void AnalyzerCore::FillNtp(TString hname, Double_t myinput[]){

  if (GetNtp(hname)) GetNtp(hname)->Fill(myinput);
  else m_logger << INFO << hname << " was NOT found. Check you ntp. " << LQLogger::endmsg;

  return;
}

// //
void AnalyzerCore::WriteNtp(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TNtupleD*>::iterator mapit = mapntp.begin(); mapit != mapntp.end(); mapit++){
    mapit->second->Write();
  }

  return;
}
// //



// //
TNtupleD* AnalyzerCore::GetNtp(TString hname){

  TNtupleD* n = NULL;
  std::map<TString, TNtupleD*>::iterator mapit = mapntp.find(hname);
  if (mapit != mapntp.end()) return mapit->second;
  else m_logger << INFO << hname << " was not found in map" << LQLogger::endmsg;

  return n;
}

void AnalyzerCore::FillLeptonKinematicPlot(std::vector<KLepton> lep, TString suffix, double w){

  //==== prefix : "ZLepton"
  //==== suffix : "WZ"

  TString OrderStr[5] = {"leading", "second", "third", "fourth", "fifth"};

  int counter_muon(0), counter_electron(0);

  for(unsigned int i=0; i<lep.size(); i++){

    //==== return more than 5 leptons
    if(i>=5) break;

    //==== if lepton flavour is not set, return
    if(lep.at(i).LeptonFlavour()==KLepton::NOTSET){
      m_logger << INFO << "lepton flavour not set" << LQLogger::endmsg;
      continue;
    }


    //==== Lepton plot

    TString this_order = OrderStr[i];

    FillHist(this_order+"Lepton_Pt_"+suffix, lep.at(i).Pt(), w, 0., 200., 200);
    FillHist(this_order+"Lepton_Eta_"+suffix, lep.at(i).Eta(), w, -3., 3., 60);
    FillHist(this_order+"Lepton_RelIso_"+suffix, lep.at(i).RelIso(), w, 0., 1.0, 100);
    FillHist(this_order+"Lepton_dXY_"+suffix, fabs(lep.at(i).dXY()), w, 0., 0.01, 100);
    FillHist(this_order+"Lepton_dXYSig_"+suffix, fabs(lep.at(i).dXYSig()), w, 0., 8., 80);
    FillHist(this_order+"Lepton_dZ_"+suffix, fabs(lep.at(i).dZ()), w, 0., 0.5, 500);

    //==== Muon plot

    if(lep.at(i).LeptonFlavour()==KLepton::MUON){

      const snu::KMuon* this_muon = lep.at(i).GetMuonPtr();
      TString this_muon_order = OrderStr[counter_muon];

      FillHist(this_muon_order+"Muon_Pt_"+suffix, this_muon->Pt(), w, 0., 200., 200);
      FillHist(this_muon_order+"Muon_Eta_"+suffix, this_muon->Eta(), w, -3., 3., 60);
      FillHist(this_muon_order+"Muon_RelIso_"+suffix, this_muon->RelIso04(), w, 0., 1.0, 100);
      FillHist(this_muon_order+"Muon_dXY_"+suffix, fabs(this_muon->dXY()), w, 0., 0.01, 100);
      FillHist(this_muon_order+"Muon_dXYSig_"+suffix, fabs(this_muon->dXYSig()), w, 0., 8., 80);
      FillHist(this_muon_order+"Muon_dZ_"+suffix, fabs(this_muon->dZ()), w, 0., 0.5, 500);
      FillHist(this_muon_order+"Muon_GlobalChi2_"+suffix, this_muon->GlobalChi2(), w, 0., 100., 100);

      counter_muon++;
    }

    //==== Electron plot

    else if(lep.at(i).LeptonFlavour()==KLepton::ELECTRON){

      const snu::KElectron* this_electron = lep.at(i).GetElectronPtr();
      TString this_electron_order = OrderStr[counter_electron];

      FillHist(this_electron_order+"Electron_Pt_"+suffix, this_electron->Pt(), w, 0., 200., 200);
      FillHist(this_electron_order+"Electron_Eta_"+suffix, this_electron->Eta(), w, -3., 3., 60);
      FillHist(this_electron_order+"Electron_RelIso_"+suffix, this_electron->PFRelIso(0.3), w, 0., 1.0, 100);
      FillHist(this_electron_order+"Electron_dXY_"+suffix, fabs(this_electron->dxy()), w, 0., 0.01, 100);
      FillHist(this_electron_order+"Electron_dXYSig_"+suffix, fabs(this_electron->dxySig()), w, 0., 8., 80);
      FillHist(this_electron_order+"Electron_dZ_"+suffix, fabs(this_electron->dz()), w, 0., 0.5, 500);

      counter_electron++;

    }

    else{

    }

  }

}

void AnalyzerCore::FillUpDownLeptonKinematicPlot(std::vector<KLepton> lep, TString suffix, double w, double w_err){

  FillLeptonKinematicPlot(lep, suffix+"_up", w+w_err);
  FillLeptonKinematicPlot(lep, suffix+"_down", w-w_err);
  FillLeptonKinematicPlot(lep, suffix, w);

}


void AnalyzerCore::SetPlotHNTriLepMetInfo(double met, double metphi){

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    it->second->SetMetInfo(met,metphi);
  }

}

void AnalyzerCore::SetPlotHNTriLepParticleInfo(snu::KParticle *hn, snu::KParticle w_pri_lowmass, snu::KParticle w_pri_highmass, snu::KParticle w_sec){

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    it->second->SetParticleInfo(hn, w_pri_lowmass, w_pri_highmass, w_sec);
  }

}

void AnalyzerCore::SetPlotHNTriLepChargeSign(double os, double ss0, double ss1){

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    it->second->SetChargeSign(os, ss0, ss1);
  }

}

void AnalyzerCore::SetPlotHNTriBJet(int nbjet){

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    it->second->SetBJet(nbjet);
  }

}

std::vector<snu::KMuon> AnalyzerCore::sort_muons_ptorder(std::vector<snu::KMuon> muons){

  std::vector<snu::KMuon> outmuon;
  while(outmuon.size() != muons.size()){
    double this_maxpt = 0.;
    int index(0);
    for(unsigned int i=0; i<muons.size(); i++){
      bool isthisused = std::find( outmuon.begin(), outmuon.end(), muons.at(i) ) != outmuon.end();
      if(isthisused) continue;
      if( muons.at(i).Pt() > this_maxpt ){
        index = i;
        this_maxpt = muons.at(i).Pt();
      }
    }
    outmuon.push_back( muons.at(index) );
  }
  return outmuon;
 

}

std::vector<snu::KElectron> AnalyzerCore::sort_electrons_ptorder(std::vector<snu::KElectron> electrons){

  std::vector<snu::KElectron> outelectron;
  while(outelectron.size() != electrons.size()){
    double this_maxpt = 0.;
    int index(0);
    for(unsigned int i=0; i<electrons.size(); i++){
      bool isthisused = std::find( outelectron.begin(), outelectron.end(), electrons.at(i) ) != outelectron.end();
      if(isthisused) continue;
      if( electrons.at(i).Pt() > this_maxpt ){
        index = i;
        this_maxpt = electrons.at(i).Pt();
      }
    }
    outelectron.push_back( electrons.at(index) );
  }
  return outelectron;



}

std::vector<KLepton> AnalyzerCore::sort_leptons_ptorder(std::vector<KLepton> leptons){

  std::vector<KLepton> outlepton;
  while(outlepton.size() != leptons.size()){
    double this_maxpt = 0.;
    int index(0);
    for(unsigned int i=0; i<leptons.size(); i++){
      bool isthisused = std::find( outlepton.begin(), outlepton.end(), leptons.at(i) ) != outlepton.end();
      if(isthisused) continue;
      if( leptons.at(i).Pt() > this_maxpt ){
        index = i;
        this_maxpt = leptons.at(i).Pt();
      }
    }
    outlepton.push_back( leptons.at(index) );
  }
  return outlepton;


}

bool AnalyzerCore::PassOptimizedCut(
  int sig_mass, double first_pt, double second_pt, double third_pt,
  double W_pri_mass, double HN_mass, 
  double deltaR_OS_min, double gamma_star_mass,
  double PFMET){

  double cut_first_pt(0.), cut_second_pt(0.), cut_third_pt(0.),
         cut_W_pri_mass(0.), cut_HN_mass(0.), 
         cut_deltaR_OS_min(0.), cut_gamma_star_mass(0.),
         cut_PFMET(0.);

  if(sig_mass == 5){
    cut_first_pt = 40.;
    cut_second_pt = 35.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 110.;
    cut_HN_mass = 35.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 10){
    cut_first_pt = 40.;
    cut_second_pt = 30.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 135.;
    cut_HN_mass = 40.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 20){
    cut_first_pt = 40.;
    cut_second_pt = 35.;
    cut_third_pt = 20.;
    cut_W_pri_mass = 130.;
    cut_HN_mass = 40.;
    cut_deltaR_OS_min = 0.4;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 30){
    cut_first_pt = 40.;
    cut_second_pt = 40.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 115.;
    cut_HN_mass = 40.;
    cut_deltaR_OS_min = 0.5;
    cut_gamma_star_mass = 10.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 40){
    cut_first_pt = 40.;
    cut_second_pt = 35.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 100.;
    cut_HN_mass = 50.;
    cut_deltaR_OS_min = 0.7;
    cut_gamma_star_mass = 10.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 50){
    cut_first_pt = 30.;
    cut_second_pt = 30.;
    cut_third_pt = 30.;
    cut_W_pri_mass = 125.;
    cut_HN_mass = 80.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 20.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 60){
    cut_first_pt = 30.;
    cut_second_pt = 30.;
    cut_third_pt = 30.;
    cut_W_pri_mass = 100.;
    cut_HN_mass = 80.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 20.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 70){
    cut_first_pt = 35.;
    cut_second_pt = 35.;
    cut_third_pt = 30.;
    cut_W_pri_mass = 100.;
    cut_HN_mass = 80.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 15.;
    cut_PFMET = 0.;
  }
  else if(sig_mass == 90){
    cut_first_pt = 45.;
    cut_second_pt = 40.;
    cut_third_pt = 10.;
    cut_W_pri_mass = 80.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 100){
    cut_first_pt = 25.;
    cut_second_pt = 15.;
    cut_third_pt = 15.;
    cut_W_pri_mass = 110.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 150){
    cut_first_pt = 45.;
    cut_second_pt = 40.;
    cut_third_pt = 20.;
    cut_W_pri_mass = 160.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 0;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 200){
    cut_first_pt = 65.;
    cut_second_pt = 50.;
    cut_third_pt = 25.;
    cut_W_pri_mass = 250.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 0.7;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 300){
    cut_first_pt = 100.;
    cut_second_pt = 55.;
    cut_third_pt = 35.;
    cut_W_pri_mass = 330.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 1.4;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 400){
    cut_first_pt = 120.;
    cut_second_pt = 65.;
    cut_third_pt = 50.;
    cut_W_pri_mass = 480.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 1.4;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 500){
    cut_first_pt = 130.;
    cut_second_pt = 80.;
    cut_third_pt = 40.;
    cut_W_pri_mass = 550.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 1.4;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 40.;
  }
  else if(sig_mass == 700){
    cut_first_pt = 200.;
    cut_second_pt = 80.;
    cut_third_pt = 45.;
    cut_W_pri_mass = 760.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 1.4;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 20.;
  }
  else if(sig_mass == 1000){
    cut_first_pt = 290.;
    cut_second_pt = 180.;
    cut_third_pt = 50.;
    cut_W_pri_mass = 900.;
    cut_HN_mass = 9999999.;
    cut_deltaR_OS_min = 0.7;
    cut_gamma_star_mass = 0.;
    cut_PFMET = 50.;
  }
  else{
    cout << "[trilepton_mumumu::PassOptimizedCut] Signal mass wrong" << endl;
  }

  bool pass = true;
  if(sig_mass < 80){
    if( ! (first_pt < cut_first_pt) ) pass = false;
    if( ! (second_pt < cut_second_pt) ) pass = false;
    if( ! (third_pt < cut_third_pt) ) pass = false;
    if( ! (W_pri_mass < cut_W_pri_mass) ) pass = false;
    if( ! (HN_mass < cut_HN_mass ) ) pass = false;
    if( ! (deltaR_OS_min > cut_deltaR_OS_min) ) pass = false;
    if( ! (gamma_star_mass > cut_gamma_star_mass) ) pass = false;
  }
  else{
    if( ! (first_pt > cut_first_pt) ) pass = false;
    if( ! (second_pt > cut_second_pt) ) pass = false;
    if( ! (third_pt > cut_third_pt) ) pass = false;
    if( ! (PFMET > cut_PFMET) ) pass = false;
    if( ! (W_pri_mass > cut_W_pri_mass) ) pass = false;
    if( ! (deltaR_OS_min > cut_deltaR_OS_min) ) pass = false;
  }

  return pass;

}

int AnalyzerCore::find_genmatching(snu::KTruth gen, std::vector<KLepton> recos, std::vector<int>& used_index){

  double mindr = 0.1;
  double maxPtDiff = 9999.;
  int found=-1;
  for(unsigned int i=0; i<recos.size(); i++){
    if(abs(gen.PdgId())==13){
      if(recos.at(i).LeptonFlavour() != KLepton::MUON) continue;
    }
    else if(abs(gen.PdgId())==11){
      if(recos.at(i).LeptonFlavour() != KLepton::ELECTRON) continue;
    }
    else{
      cout << "HNGenMatched lepton is not lepton?!" << endl;
      cout << "  gen.PdgId() = " << gen.PdgId() << endl;
    }
    
    //cout << "["<<i<<"th reco] : pt = " << recos.at(i).Pt() << ", eta = " << recos.at(i).Eta() << endl;
    double dr = gen.DeltaR(recos.at(i));
    bool ptmatched(true);
    //if( fabs(gen.Pt() - recos.at(i).Pt()) / TMath::Min( gen.Pt(), recos.at(i).Pt() ) >= maxPtDiff ) ptmatched = false;
    bool isthisused = std::find(used_index.begin(), used_index.end(), i) != used_index.end();
    if(dr < mindr && ptmatched && !isthisused){
      mindr = dr;
      found = i;
    }
  }
  used_index.push_back(found);
  return found;
}


double AnalyzerCore::MuonConePt(snu::KMuon muon, double tightiso){

  double mu_pt_corr = muon.Pt()*(1+max(0.,(muon.RelIso04()-tightiso))) ;

  return mu_pt_corr;

}














void AnalyzerCore::SetHNTriCutOp(TString filepath){

  string cutline;
  ifstream in(filepath);
  while(getline(in,cutline)){
    std::istringstream is( cutline );
    TString channel;
    TString mass;
    TString var;
    double value;
    is >> channel;
    is >> mass;
    is >> var;
    is >> value;

    TString thiskey = channel+"_"+mass+"_"+var;

    map_HNTriChannl_cutop[thiskey] = value;
  }

}

bool AnalyzerCore::PassOptimizedCut( 
  TString channel, int sig_mass,
  double first_pt, double second_pt, double third_pt,
  double W_pri_mass, double HN_mass, 
  double deltaR_OS_min, double gamma_star_mass,
  double PFMET){

  double cut_first_pt(0.), cut_second_pt(0.), cut_third_pt(0.),
         cut_W_pri_mass(0.), cut_HN_mass(0.), 
         cut_deltaR_OS_min(0.), cut_gamma_star_mass(0.),
         cut_PFMET(0.);

  TString thiskey_prefix = channel+"_"+TString::Itoa(sig_mass,10);

  std::vector<TString> variables;

  variables.push_back("first_pt");
  variables.push_back("second_pt");
  variables.push_back("third_pt");
  variables.push_back("W_pri_mass");
  variables.push_back("HN_mass");
  variables.push_back("deltaR_OS_min");
  variables.push_back("gamma_star_mass");
  variables.push_back("PFMET");

  std::map< TString, double > cut_variables;

  for(unsigned int i=0; i<variables.size(); i++){
    TString key = thiskey_prefix+"_"+variables.at(i);
    //cout << "[filling] "<<key<<endl;

    std::map< TString, double >::const_iterator mapit;
    mapit = map_HNTriChannl_cutop.find(key);
    bool valuefound = (mapit!=map_HNTriChannl_cutop.end());

    if(valuefound){
      cut_variables[variables.at(i)] = map_HNTriChannl_cutop[key];
    }
    else{
      cerr << "Did not find optimized cut value for "<<key << endl;
      exit(EXIT_FAILURE);
    }
  }

  bool pass = true;

//  for(int i=0; i<n_var; i++){
//    cout << channel << "\t" << sig_mass << "\t" << variables[i] << "\t" << cut_variables[i] << endl;
//  }

  cut_first_pt =  cut_variables["first_pt"];
  cut_second_pt = cut_variables["second_pt"];
  cut_third_pt = cut_variables["third_pt"];
  cut_W_pri_mass = cut_variables["W_pri_mass"];
  cut_HN_mass = cut_variables["HN_mass"];
  cut_deltaR_OS_min = cut_variables["deltaR_OS_min"];
  cut_gamma_star_mass = cut_variables["gamma_star_mass"];
  cut_PFMET = cut_variables["PFMET"];

  if(channel == "MuMuMu"){
    if(sig_mass < 80){
      if( ! (first_pt < cut_first_pt) ) pass = false;
      if( ! (second_pt < cut_second_pt) ) pass = false;
      if( ! (third_pt < cut_third_pt) ) pass = false;
      if( ! (W_pri_mass < cut_W_pri_mass) ) pass = false;
      if( ! (HN_mass < cut_HN_mass ) ) pass = false;
      if( ! (deltaR_OS_min > cut_deltaR_OS_min) ) pass = false;
      if( ! (gamma_star_mass > cut_gamma_star_mass) ) pass = false;
    }
    else{
      if( ! (first_pt > cut_first_pt) ) pass = false;
      if( ! (second_pt > cut_second_pt) ) pass = false;
      if( ! (third_pt > cut_third_pt) ) pass = false;
      if( ! (W_pri_mass > cut_W_pri_mass) ) pass = false;
      if( ! (PFMET > cut_PFMET) ) pass = false;
      if( ! (deltaR_OS_min > cut_deltaR_OS_min) ) pass = false;
    }
  }

  if(channel == "SSSF_MuMuE"){
    if(sig_mass < 80){
      if( ! (first_pt < cut_first_pt) ) pass = false;
      if( ! (second_pt < cut_second_pt) ) pass = false;
      if( ! (third_pt < cut_third_pt) ) pass = false;
      if( ! (W_pri_mass < cut_W_pri_mass) ) pass = false;
      //if( ! (HN_mass < cut_HN_mass ) ) pass = false;
      if( ! (PFMET < cut_PFMET) ) pass = false;
    }
    else{
      if( ! (first_pt > cut_first_pt) ) pass = false;
      if( ! (second_pt > cut_second_pt) ) pass = false;
      if( ! (third_pt > cut_third_pt) ) pass = false;
      if( ! (W_pri_mass > cut_W_pri_mass) ) pass = false;
      if( ! (PFMET > cut_PFMET) ) pass = false;
    }
  }

  return pass;

}










