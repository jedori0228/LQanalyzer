// STD includes
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <sstream>

/// Local includes
#include "HistUtils.hpp"

///ROOT includes
#include "TH1.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include <TStyle.h>
#include "TLatex.h"


using namespace std;

/// LIST OF FUNCIONS FOR PLOTTER CODE
int MakePlots(std::string hist);
void MakeCutFlow(std::string hist);
int MakeCutFlow_Plots(string configfile);
void PrintCanvas(TCanvas* c1, std::string folder, std::string title);
bool repeat(string hname);
TLegend* MakeLegend(map<TString, TH1*> legmap,TH1* h_legdata, bool rundata, bool log);
TH1* MakeDataHist(string name, double xmin, double xmax, TH1* h_up,bool ylog , int rebin);
void CheckHist(TH1* h);
void CheckSamples(int nsamples);
vector<pair<TString,float> >  InitSample (TString sample);
THStack* MakeStack(vector<pair<pair<vector<pair<TString,float> >, int >, TString > > sample, TString type, string name, float xmin, float xmax,map<TString, TH1*>& legmap, int rebin);
void SetErrors(TH1* hist, float normerr);
TH1* MakeStackUp(map<TString, TH1*> map_of_stacks, TString clonename);
TH1* MakeStackDown(map<TString, TH1*> map_of_stacks, TString clonename);
TH1* MakeSumHist(THStack* thestack);
float  GetMaximum(TH1* h_data, TH1* h_up, bool ylog, string name);
void SetTitles(TH1* hist, string name);
bool HistInGev(string name);
void fixOverlay();
void setTDRStyle();
TCanvas* CompDataMC(TH1* hdata, vector<THStack*> mcstack,TH1* hup, TH1* hdown,TLegend* legend, const string hname, const int rebin, double xmin, double xmax,double ymin, double ymax,string path , string folder, bool logy, bool usedata, TString channel);
TH1* MakeSumHist2(THStack* thestack);
TH1* MakeErrorBand(TH1* hnom, TH1* hup, TH1* hdown);
void SetNomBinError(TH1* hnom, TH1* hup, TH1* hdown);
void MakeLabel( float rhcol_x, float rhcol_y);
std::map<std::string,std::string> _htmls;

float GetSyst(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples );
float GetSystPercent(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples );
void setZ(bool useAlpgen);
float Calculate(TString cut, TString variance,  pair<vector<pair<TString,float> >, TString >  samples);
void  SetUpConfig(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples , vector<string>& cut_label);
void  SetUpMasterConfig(std::string filename);


////// For cutflow
float Error(TH1* h);
float GetTotal(TString cut, vector<pair<TString,float> > sample);
float GetStatError(TString cut, vector<pair<TString,float> > sample);
float GetStatError2(TString cut, vector<pair<TString,float> > sample);
float GetIntegral(TString cut, TString isample, TString type);
float GetNormErr(TString cut,  vector<pair<TString,float> > samples);
float GetNormErr2(TString cut,  vector<pair<TString,float> > samples);
float GetErr(TString cut,  vector<pair<TString,float> > samples, TString err_type,TString var);
float GetErr2(TString cut,  vector<pair<TString,float> > samples, TString err_type,TString var);
float GetError(TString cut, TString isample, TString type);


//// GLOBAL VARIABLES
int isig=0;
map<string,int> norepeatplot;
TString columnname="";
TString caption="";
  
std::string hist;
bool showdata=true;
std::string cutfile;
std::string histfile;
bool ylog;
bool usenp(false);

TString channel;

std::string path;
std::string message;
std::string fileprefix="";
std::string filepostfix = "";

std::ofstream page;
std::ofstream histpage;

vector<string> cuts; 
vector<string> allcuts; 
vector<string> listofsamples; 

//// Standard bkg folders 
string  mcloc="";
/// Data folder
string dataloc = "";
/// data driven
string datadrivenloc= "";
string plotloc ="";
string cutloc ="";

string histdir="";






//##################################################################
//###### MAIN CALL
//##################################################################

int main(int argc, char *argv[]) {
  
  /////////////////////////////////////////////////////
  //
  //  3 stages to setup : 
  //  1: McatNloWZ specifies with WZ MC we plot
  //  2: reg: specifies which fake region is plotted
  //  3: BM specifies if we plot 1.2 or 4.7 fb-1
  //
  //////////////////////////////////////////////////////
  
  TH1::SetDefaultSumw2(true);
  
  if(argc == 1) {
    cout << "No config file set" << endl;
    return 0;  
  }
  
  /// Set Plotting style
  setTDRStyle();
  gStyle->SetPalette(1);
  
  //read in config
  
  for (int i = 1; i < argc; ++i) {
    string configfile = argv[i];   
    
    SetUpMasterConfig(configfile);
    
    int a =MakeCutFlow_Plots(configfile);
  }
  
  return 0;
}


void setTDRStyle() {
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
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
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
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
  
  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  //  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);

  
  // tdrStyle->SetDateY(Float_t y = 0.01);
  
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
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);
  
  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

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


  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}



// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}


int MakeCutFlow_Plots(string configfile){
  
  std::string pname = "/home/jalmond/WebPlots/"+ path + "/indexCMS.html";
  std::string phistname = "/home/jalmond/WebPlots/"+ path + "/histograms/" + histdir  + "/indexCMS.html";
  system(("mkdir /home/jalmond/WebPlots/" + path).c_str());
  system(("mkdir /home/jalmond/WebPlots/" + path+ "/histograms/").c_str());
  system(("mkdir /home/jalmond/WebPlots/" + path+"/histograms/" + histdir + "/").c_str());
  cout << "HIST page is set to " << phistname.c_str() << endl;

  histpage.open(phistname.c_str());
  page.open(pname.c_str());
  
  page << "<html><font face=\"Helvetica\"><head><title> HvyN Analysis </title></head>" << endl;
  page << "<body>" << endl;
  page << "<h1> HvyN Analysis Plots </h1>" << endl;
  page << "<br> <font size=\"4\"><b> " << message <<  " </b></font> <br><br>" << endl;
  
  page << "<a href=\"histograms/" +histdir + "/indexCMS.html\">"+ histdir + "</a><br>"; 


  int M=MakePlots(histdir);  
  //MakeCutFlow(histdir);

  return M;

}



int MakePlots(string hist) {

  
  cout << "\n ---------------------------------------- " << endl;
  cout << "MakeDataMCCompPlots::MakePlots(string hist) " << endl;

  ////////////////////// ////////////////
  ////  MAIN PART OF CODE for user/
  ///////////////////////////////////////
  //// What samples to use in histogram
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > samples;  
  vector<string> cut_label;
  //// Sets flags for using CF/NP/logY axis/plot data/ and which mc samples to use
  
  SetUpConfig( samples, cut_label);  

  cuts.clear();

  // ----------Get list of cuts to plot  ----------------------
  ifstream cut_name_file(cutfile.c_str());
  if(!cut_name_file) {
    cerr << "Did not find " + cutfile + ", exiting ..." << endl;
    return 1;
  }
  while(!cut_name_file.eof()) {
    string cutname;
    cut_name_file >> cutname;
    if(cutname=="END") break;
    allcuts.push_back(cutname);
    cout << "Added " << cutname << endl;
  }
  

  ifstream histo_name_file(histfile.c_str());
  if(!histo_name_file) {
    cerr << "Did not find " << histfile << ", exiting ..." << endl;
    return 1;
  }
  
  histpage << "<table border = 1><tr>"
	   << "<th> <a name=\"Cut : Plot\">Cut : Plot</a> </th>"
	   << "<th> Plots </th>"
	   << "</tr>" << endl;
  
  while(!histo_name_file.eof()) {
    string h_name;
    int rebin;
    double xmin,xmax;
    histo_name_file >> h_name;
    if(repeat(h_name))continue;
    if(h_name=="END") break;
    histo_name_file >> rebin;
    histo_name_file >> xmin;
    histo_name_file >> xmax;
    
    if(h_name.find("#")!=string::npos) continue;
    
    for(unsigned int ncut=0; ncut<allcuts.size();  ncut++){
      string name = h_name+allcuts.at(ncut);
       
	cout << "\n------------------------------------------------------- \n" << endl;
	cout << "Making histogram " << name << endl;
	
	
	/// Make nominal histogram stack
	map<TString, TH1*> legmap;
	
	cout << "Making Nominal histogram " << endl;
	THStack* mstack=  MakeStack(samples , "Nominal",name, xmin, xmax, legmap, rebin );
	
	//// mhist sets error config
	map<TString,TH1*> mhist;
	mhist["Nominal"] = MakeSumHist(mstack);
	
	TH1* hup = MakeStackUp(mhist, name+"UP");
	TH1* hdown = MakeStackDown(mhist, name+"DOWN");
		
	cout << "Final Background Integral = " <<  MakeSumHist(mstack)->Integral() << " : Up = " << hup->Integral() << " : Down= " << hdown->Integral() << endl;
	
	/// Make data histogram
	TH1* hdata = MakeDataHist(name, xmin, xmax, hup, ylog, rebin);
	CheckHist(hdata);	
	float ymin (0.), ymax( 1000000.);
	ymax = GetMaximum(hdata, hup, ylog, name);
  
	if(showdata)cout << "Total data = " <<  hdata->Integral() << endl;
	
	/// Make legend
	TLegend* legend = MakeLegend(legmap, hdata, showdata, ylog);       		
	
        vector<THStack*> vstack;		
	vstack.push_back(mstack);   	
	

	cout << " Making canvas" << endl;
	
	TCanvas* c = CompDataMC(hdata, vstack,hup,hdown, legend,name,rebin,xmin,xmax, ymin,ymax, path, histdir,ylog, showdata, channel);      	
	cout << " Made canvas" << endl;
	PrintCanvas(c, histdir, c->GetName());
    }
  }            
  page.close();
  
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///   CODE FOR CUTFLOW
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void MakeCutFlow(string type){
  
  return ;
  vector<string> cut_label;  
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > cfsamples;  
  SetUpConfig( cfsamples, cut_label);

  cuts.clear();    
  // ----------Get list of cuts to plot  ----------------------
  ifstream cut_name_file(cutfile.c_str());
  while(!cut_name_file.eof()) {
    string cutname;
    cut_name_file >> cutname;
    if(cutname=="END") break;
    cuts.push_back((hist + cutname).c_str());    
    cout << "Making cutflow for MuonPlots/mu1_eta"<< cutname << endl;
  }

 
  vector<float> totalnumbers;
  vector<float> totalnumbersup;
  vector<float> totalnumbersdown;  
  
  int i_cut(0);
  for(vector<string>::iterator it = cuts.begin(); it!=cuts.end(); it++, i_cut++){
    
    vector<pair<vector<pair<TString,float> >, TString> > samples;   
    for(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it2 = cfsamples.begin(); it2!=cfsamples.end(); it2++){
      samples.push_back(make_pair(it2->first.first,it2->second));      
    }
    
      
    /// Vectors for cutflow	table
    map<TString,float> samples_numbers;
    map<TString,float> samples_numbers_staterr;
    map<TString,float> samples_numbers_up;
    map<TString,float> samples_numbers_down;
    
    /// Vector for systematic table
    map<TString,float> syst_stat;
    map<TString,float> syst_JESup;
    map<TString,float> syst_JESdown;
    map<TString,float> syst_MUONISOup;
    map<TString,float> syst_MUONISOdown;
    map<TString,float> syst_MUONRECOup;
    map<TString,float> syst_MUONRECOdown;
    map<TString,float> syst_MUONSCALEup;
    map<TString,float> syst_MUONSCALEdown;
    map<TString,float> syst_ELECTRONISOup;
    map<TString,float> syst_ELECTRONISOdown;
    map<TString,float> syst_ELECTRONIDup;
    map<TString,float> syst_ELECTRONIDdown;
    map<TString,float> syst_ELECTRONRECOup;
    map<TString,float> syst_ELECTRONRECOdown;
    map<TString,float> syst_ELECTRONSCALEup;
    map<TString,float> syst_ELECTRONSCALEdown;
    map<TString,float> syst_ELECTRONSMEARup;
    map<TString,float> syst_ELECTRONSMEARdown;
    map<TString,float> syst_CFup;
    map<TString,float> syst_CFdown;
    map<TString,float> syst_JVFup;
    map<TString,float> syst_JVFdown;
    map<TString,float> syst_JER;
    map<TString,float> syst_norm;
    map<TString,float> syst_total;

    
    float totalbkg(0.),totalbkgdown(0.), totalbkgup(0.); 
    for(vector<pair<vector<pair<TString,float> >, TString> >::iterator it2 = samples.begin() ; it2!= samples.end(); it2++){
      TString cutname = *it;
      totalbkg+= Calculate(*it,"Normal",*it2);
      totalbkgup+= Calculate(*it,"Up",*it2);
      totalbkgdown+= Calculate(*it,"Down",*it2);
      
      samples_numbers[it2->second] = Calculate(cutname,"Normal",*it2);
      samples_numbers_up[it2->second] = Calculate(cutname,"Down",*it2);
      samples_numbers_down[it2->second] = Calculate(cutname,"Up",*it2);
      samples_numbers_staterr[it2->second] = Calculate(cutname,"StatErr",*it2);
    }	
    
    
    float totaldata(0.);
    if(showdata) totaldata = GetIntegral(*it,"data","data");
    float errdata(0.);
    if(showdata) errdata= GetError(*it,"data","data");
    
    
    float totalerr_up(0.),totalerr_down(0.),totalerrup(0.),totalerrdown(0.),total_staterr(0.);
    for(map<TString,float>::iterator mapit = samples_numbers.begin(); mapit!= samples_numbers.end(); mapit++){
      
      map<TString,float>::iterator mapit_up;
      mapit_up = samples_numbers_up.find(mapit->first);
      map<TString,float>::iterator mapit_down;
      mapit_down = samples_numbers_down.find(mapit->first);
      map<TString,float>::iterator mapit_stat;
      mapit_stat = samples_numbers_staterr.find(mapit->first);
      totalerr_up += ( mapit_up->second*mapit_up->second+ mapit_stat->second*mapit_stat->second); 
      totalerr_down += (mapit_down->second*mapit_down->second + mapit_stat->second*mapit_stat->second); 
      totalerrup += (mapit_up->second*mapit_up->second); 
      totalerrdown += (mapit_down->second*mapit_down->second); 
      total_staterr += mapit_stat->second*mapit_stat->second;
      cout << mapit->first << " background = " << mapit->second << " +- " << mapit_stat->second << " + " << mapit_up->second << " - " << mapit_down->second <<  endl;      
   
      
    }
  
    
    totalerr_up = sqrt(totalerr_up);
    totalerr_down = sqrt(totalerr_down);
    total_staterr = sqrt(total_staterr);
    totalerrup = sqrt(totalerrup);
    totalerrdown = sqrt(totalerrdown);
    
    cout << "Total Bkg   = " << totalbkg << "+- " << total_staterr << " + " << totalerrup << " - " << totalerrdown << endl;
    if(showdata)cout << "Total Data  = " << totaldata << endl;
    cout << "-------------" << endl;
    if(totaldata > totalbkg)cout <<"Significance = " << (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_up*totalerr_up))) << endl;
    else cout <<"Significance = " << (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_down*totalerr_down))) << endl;
    float significance = (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_up*totalerr_up))) ;
    
    if(significance < 0.) significance = (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_down*totalerr_down))) ;
    
    ofstream ofile;
    
    cout << cut_label.size() << endl;
    string latex =  "Tables/" + cut_label.at(i_cut) + "Table.txt";
    
    ofile.open(latex.c_str());
    ofile.setf(ios::fixed,ios::floatfield); 
    
    ofile.precision(1);
    ofile << "\\begin{table}[h]" << endl;
    ofile << "\\begin{center}" << endl;
    ofile << "\\begin{tabular}{lr@{\\hspace{0.5mm}}c@{\\hspace{0.5mm}}c@{\\hspace{0.5mm}}l}" << endl;
    ofile << "\\hline" << endl;
    ofile << "\\hline" << endl;   
    //ofile << "Source & \\multicolumn{4}{c}{$\\mu\\mu\\mu$} \\"<<"\\" << endl;
    ofile << "Source & \\multicolumn{4}{c}{" << columnname << "} \\"<<"\\" << endl;
    ofile << "\\hline" << endl;   
	 
    
    for(map<TString,float>::iterator mapit = samples_numbers.begin(); mapit!= samples_numbers.end(); mapit++){
      
      map<TString,float>::iterator mapit_up;
      mapit_up = samples_numbers_up.find(mapit->first);
      map<TString,float>::iterator mapit_down;
      mapit_down = samples_numbers_down.find(mapit->first);
      map<TString,float>::iterator mapit_stat;
      mapit_stat = samples_numbers_staterr.find(mapit->first);
      
      if(mapit->second!=0.0){
	ofile << mapit->first + "&" <<  mapit->second << "& $\\pm$& "  << mapit_stat->second <<  "&$^{+" <<  mapit_up->second << "}_{-" <<  mapit_down->second  << "}$" ; 
	ofile  <<  "\\"  << "\\" << endl;	   
      }
    }
    ofile << "\\hline" << endl;
    ofile << "Total&" << totalbkg << "& $\\pm$&"  << total_staterr << "&$^{+" << totalerrup  << "}_{-" << totalerrdown << "}$" ; 
    ofile  <<  "\\"  << "\\" << endl;
    ofile << "\\hline" << endl;
    
    ofile << "Data&  \\multicolumn{4}{c}{$" << totaldata << "$}\\" << "\\" <<endl;
    ofile << "\\hline" << endl;
    if(significance < 0) ofile << "Signficance&  \\multicolumn{4}{c}{$" << significance << "\\sigma$}\\" << "\\" <<endl;
    if(significance > 0) ofile << "Signficance&  \\multicolumn{4}{c}{$+" << significance << "\\sigma$}\\" << "\\" <<endl;
    ofile << "\\hline" << endl;   
    ofile << "\\hline" << endl;   
    ofile << "\\end{tabular}" << endl;
    ofile << "\\caption{" << caption << "}" << endl;
    ofile << "\\end{center}" << endl;
    ofile << "\\end{table}" << endl;    
  }
  
  
  string latex_command = "latex Tables/" + cut_label.at(0) +".tex";
  string dvi_command = "dvipdf " + cut_label.at(0) +".dvi";
  string mv_command = "mv " + cut_label.at(0) +".pdf /home/jalmond/WebPlots/" + path +"/histograms/"+ histdir ;
  
  system((latex_command.c_str()));
  system((dvi_command.c_str()));
  system((mv_command.c_str()));
  system(("rm *aux"));
  system(("rm *log"));
  system(("rm *dvi"));
  
  string cftitle = cut_label.at(0);
  
  histpage << "<tr><td>"<< "cutflow " + cut_label.at(0)  <<"</td>"<<endl;
  histpage <<"<td>"<<endl;
  histpage << "<a href=\"" << cut_label.at(0)  << ".pdf\">";
  histpage << "<img src=\"" << cut_label.at(0)  << ".pdf\" width=\"100%\"/>";
  histpage << "</td>" << endl;

  
  return;
   }

bool repeat (string hname){
  map<string,int>::iterator mit = norepeatplot.find(hname);
  if(mit!=norepeatplot.end())return true;
  else{
    norepeatplot[hname]=1;
    return false;
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintCanvas(TCanvas* c1, string folder, string title){

  std::string tpdf = "/home/jalmond/WebPlots/"+ path +  "/histograms/"+folder+"/"+title;
  string plot_description;

  plot_description = title;
  
  
  if(plot_description.empty())plot_description=title;
  histpage << "<tr><td>"<< plot_description <<"</td>"<<endl;
  histpage <<"<td>"<<endl;
  histpage << "<a href=\"" << title.c_str() << ".png\">";
  histpage << "<img src=\"" << title.c_str() << ".png\" width=\"100%\"/>";
  histpage << "</td>" << endl;
  
  return;
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


TLegend* MakeLegend(map<TString, TH1*> map_legend,TH1* hlegdata,  bool rundata , bool logy){
  
  double x1 = 0.75;
  double y1 = 0.6;
  double x2 = 0.9;
  double y2 = 0.85;

  if(logy){
    
    x1 = 0.75;
    y1 = 0.7;
    x2 = 0.92;
    y2 = 0.9;

  }

  if(!showdata){
    x1= 0.6;
    x2= 0.85;
    y1= 0.6;
    y2= 0.9;
  }

  TLegend* legendH = new TLegend(x1,y1,x2,y2);
  legendH->SetFillColor(10);
  legendH->SetBorderSize(0);
  legendH->SetTextSize(0.035);
  
  
  if(rundata) 	legendH->AddEntry(hlegdata,"Data","pl");
  
  for(map<TString, TH1*>::iterator it = map_legend.begin(); it!= map_legend.end(); it++){
    legendH->AddEntry(it->second,it->first.Data(),"f");    
  }

  return legendH;
  
}




TH1* MakeDataHist(string name, double xmin, double xmax, TH1* hup, bool ylog, int rebin){

  /// Make data histogram
  TFile* file_data =  TFile::Open((dataloc).c_str());
  TH1* hdata = dynamic_cast<TH1*> ((file_data->Get(name.c_str()))->Clone());
  
  hdata->Rebin(rebin);

  float ymin (0.), ymax( 1000000.);
  if(ylog) ymin= 1.;
  ymax = GetMaximum(hdata, hup, ylog, name);
  
  /// Set Ranges / overflows
  FixOverUnderFlows(hdata, xmax);  
  
  cout << "Ymax = " << ymax << endl;
  hdata->GetXaxis()->SetRangeUser(xmin,xmax);
  hdata->GetYaxis()->SetRangeUser(ymin, ymax);

  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.2);
	
 
  //// X title  
  SetTitles(hdata, name);
  
  return hdata;

}

void CheckHist(TH1* h){
  if(!h) {
    cout << "Not able to find data histogram" << endl;
    exit(1);
  }
}


vector<pair<TString,float> >  InitSample (TString sample){
  
  vector<pair<TString,float> > list;  
  
  if(sample.Contains("DY")){
    list.push_back(make_pair("DY10to50",0.2));    
    list.push_back(make_pair("DY50plus",0.2));    
  }
  
  ///// Top samples //////////////////    
  if(sample.Contains("top")){
    list.push_back(make_pair("ttbar",0.08));
  }
  //////// Diboson ////////
  if(sample.Contains("wz")){    
    list.push_back(make_pair("WZ",0.09));
  }
  
  if(sample.Contains("zz")){
    list.push_back(make_pair("ZZ",0.09));
  }
  
  if(sample.Contains("ww")){      
    list.push_back(make_pair("WW",0.09));
  }
  
  if(sample.Contains("wjet")){
    list.push_back(make_pair("Wjets",0.09));
  }
  
  
  

  //////// SS WW /////////  
  if(sample.Contains("ssww")){         
    list.push_back(make_pair("W-W-",0.22));              
    list.push_back(make_pair("W+W+",0.22));              
  }


  if(list.size()==0) cout << "Error in making lists" << endl;
  
  return list;
}

void CheckSamples(int nsamples){
  if(nsamples==0) {
    cout << "No sample in this vector" << endl;
    exit(1);
  }
  return;
}




THStack* MakeStack(vector<pair<pair<vector<pair<TString,float> >, int >, TString > > sample, TString type, string name, float xmin, float xmax,map<TString, TH1*>& legmap , int rebin){
  
  string clonename = name;	
  THStack* stack = new THStack();
  
  bool debug(false);
  if(type.Contains("Nominal")) debug=true;

  TString fileloc = "";
 
  TDirectory* origDir = gDirectory;

  cout << "File prefix = " << fileprefix << " postfix= "  << filepostfix << endl;

  float sum_integral=0.;
  for(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it = sample.begin() ; it!= sample.end(); it++){
	
    if(type.Contains("Nominal")) fileloc = mcloc;
        
    if(!type.Contains("Nominal")) {
      if(it->first.first.at(0).first.Contains("datadriven"))fileloc=mcloc;
    }    
    
    
    CheckSamples( it->first.first.size() );

    TFile* file =  TFile::Open((fileloc+ fileprefix + it->first.first.at(0).first + filepostfix).Data());

   
    gROOT->cd();
    TDirectory* tempDir = 0;
    int counter = 0;
    while (not tempDir) {
      std::stringstream dirname;
      dirname << "WRHNCommonLeptonFakes_%i" << counter;
      if (gROOT->GetDirectory((dirname.str()).c_str())) {
	++counter;
	continue;
      }      
      tempDir = gROOT->mkdir((dirname.str()).c_str());      
    }
            

    tempDir->cd();

    TH1* h_tmp = dynamic_cast<TH1*> ((file->Get(name.c_str()))->Clone(clonename.c_str()));

    CheckHist(h_tmp);

    if(debug)cout <<  it->second <<  "  contribution " << 1 << "/" << it->first.first.size()  << " is from " << filepostfix + it->first.first.at(0).first + filepostfix <<" : Integral = " <<h_tmp->Integral() << " " << fileloc << endl;
    
    
    for(unsigned int i=1; i < it->first.first.size(); i++){	    
      clonename+="A";
             
      origDir->cd();
      
      TFile* file_loop =  TFile::Open((fileloc+ fileprefix + it->first.first.at(i).first + filepostfix).Data());	    

      tempDir->cd();
      TH1* h_loop = dynamic_cast<TH1*> ((file_loop->Get(name.c_str()))->Clone(clonename.c_str()));	    	    	    
      CheckHist(h_loop);
      h_tmp->Add(h_loop);	  	    	    
      
      if(debug)cout <<  it->second <<  "  contribution " <<i+1 <<"/" << it->first.first.size()  << " is from ExampleAnalyzer_SK" << it->first.first.at(i).first << ".NTUP_SMWZ.Reco.root : Integral = " <<h_loop->Integral() << " sum integral = " << h_tmp->Integral()    << endl;
      file_loop->Close();
    }	  	  
    
    /// TH1* is now made. Now make pretty
    FixOverUnderFlows(h_tmp, xmax);	  
    ///Set colors
    h_tmp->SetFillColor(it->first.second);
    h_tmp->SetLineColor(it->first.second);	  
    
    
    if(!stack->GetHists()) {
      stack->SetName( (string("s_") + name).c_str() );
      stack->SetTitle( (string("s_") + name).c_str() );
      SetTitles(h_tmp, name);

  
    }//stack empt   
    

    h_tmp->Rebin(rebin);
    SetErrors(h_tmp, it->first.first.at(0).second);
           
    
    stack->Add(h_tmp);
    sum_integral+=h_tmp->Integral();
    
    if(type.Contains("Nominal")) {
      legmap[it->second] = h_tmp;
    }
    
    file->Close();	  
    origDir->cd();
  }
  
  cout << type << " has integral = " << sum_integral << endl;
  
  return stack;
}


TH1* MakeStackUp(map<TString, TH1*> map_of_stacks, TString clonename){
  
  map<TString, TH1*>::iterator it =  map_of_stacks.find("Nominal");
  
  float norm=1;
  TH1* h_up = dynamic_cast<TH1*>(it->second->Clone(clonename.Data())); // copy of nominal
  
  for(int binx=1 ; binx < h_up->GetNbinsX()+1; binx++){
    float nom_content = h_up->GetBinContent(binx);
    float nom_error = h_up->GetBinError(binx);
    
    float errup2 =  nom_error*nom_error ;
    
    /// add rest of systs

    float new_bin = nom_content * 1.075; //sqrt(errup2);
        
    h_up->SetBinContent(binx,new_bin);
    
  }
  
  return  h_up;
  
}


TH1* MakeStackDown(map<TString, TH1*> map_of_stacks, TString clonename){
  
  map<TString, TH1*>::iterator it =  map_of_stacks.find("Nominal");
  
  float norm=1;
  TH1* h_down = dynamic_cast<TH1*>(it->second->Clone(clonename.Data())); // copy of nominal
  
  for(int binx=1 ; binx < h_down->GetNbinsX()+1; binx++){
    float nom_content = h_down->GetBinContent(binx);
    float nom_error = h_down->GetBinError(binx);
    
    
    //// nom_error = stat err + normalisation error, set previously on nom hist
    float errdown2 =  nom_error*nom_error;
    
    float new_bin = nom_content  * (1./1.075) ;//- sqrt(errdown2);
    
        
    h_down->SetBinContent(binx,new_bin);
    
  }
  
  return  h_down;


}


TH1* MakeSumHist(THStack* thestack){
  
  TH1* hsum=0;  
  TList* list = thestack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    
    if(!hsum) hsum = (TH1*)h->Clone( (string(h->GetName()) + "_sum").c_str() );
    else {
      hsum->Add(h, 1.0);
    }
  }//hist loop
  
  return hsum;
}


void SetErrors(TH1* hist, float normerr){

  
  for(int binx =1; binx < hist->GetNbinsX()+1; binx++){
    float newbinerr = hist->GetBinError(binx)*hist->GetBinError(binx) + hist->GetBinContent(binx)*hist->GetBinContent(binx)*normerr*normerr;
    hist->SetBinError(binx, sqrt(newbinerr));
  }
  
  return;

}



void SetTitles(TH1* hist, string name){
  
  string xtitle ="";
  string ytitle ="Entries";

  float binedge_up = hist->GetBinLowEdge(2);
  float binedge_down = hist->GetBinLowEdge(1);
  
  float width = binedge_up - binedge_down;
  
  std::ostringstream str_width;
  str_width<< width;

  if(HistInGev(name)) ytitle = "Entries / " +str_width.str() + " GeV";
  
  if(name.find("MET")!=string::npos)xtitle="E^{miss}_{T} [GeV]"; 

  if(name.find("mu_eta")!=string::npos)xtitle="Muon #eta";
  if(name.find("mu_pt")!=string::npos)xtitle="Muon p_{T} [GeV]";
  if(name.find("mu1_pt")!=string::npos)xtitle="Lead p_{T} [GeV]";
  if(name.find("mu2_pt")!=string::npos)xtitle="Second p_{T} [GeV]";
  if(name.find("mu3_pt")!=string::npos)xtitle="Third p_{T} [GeV]";


  if(name.find("el_eta")!=string::npos)xtitle="Electron #eta";
  if(name.find("el_pt")!=string::npos)xtitle="Electron p_{T} [GeV]";
  if(name.find("el1_pt")!=string::npos)xtitle="Lead p_{T} [GeV]";
  if(name.find("el2_pt")!=string::npos)xtitle="Second p_{T} [GeV]";
  if(name.find("el3_pt")!=string::npos)xtitle="Third p_{T} [GeV]";
  
  if(name.find("charge")!=string::npos)xtitle="sum of lepton charge";

  if(name.find("leaddimuma")!=string::npos)xtitle="m(#mu#mu) [GeV]";
  if(name.find("leaddielma")!=string::npos)xtitle="m(ee) [GeV]";
  if(name.find("leademuma")!=string::npos)xtitle="m(e#mu) [GeV]";
  
  if(name.find("jet_eta")!=string::npos)xtitle="jet #eta";
  if(name.find("1jet_eta")!=string::npos)xtitle="Leading jet #eta";
  if(name.find("2jet_eta")!=string::npos)xtitle="2^{nd} Leading jet #eta";
  if(name.find("njet")!=string::npos)xtitle="Number of jets";
  if(name.find("nbjet")!=string::npos)xtitle="Number of bjets";

  if(name.find("muall")!=string::npos)xtitle="m(#mu#mu#mu(#mu)) [GeV]";

  if(name.find("el1jet_mindr")!=string::npos)xtitle="min#Delta R(e_{1}j)";
  if(name.find("el2jet_mindr")!=string::npos)xtitle="min#Delta R(e_{2}j)";

  if(name.find("mu1jet_mindr")!=string::npos)xtitle="min#Delta R(#mu_{1}j)";
  if(name.find("mu2jet_mindr")!=string::npos)xtitle="min#Delta R(#mu_{2}j)";
  if(name.find("mujj_massle")!=string::npos)xtitle="m(#mu_{1}jj) [GeV]";
  if(name.find("mujj_masstr")!=string::npos)xtitle="m(#mu_{2}jj) [GeV]";
  if(name.find("mumujj_mass")!=string::npos)xtitle="m(#mu#mujj) [GeV]";

  if(name.find("ejj_massle")!=string::npos)xtitle="m(e_{1}jj) [GeV]";
  if(name.find("ejj_masstr")!=string::npos)xtitle="m(e_{2}jj) [GeV]";
  if(name.find("eejj_mass")!=string::npos)xtitle="m(eejj) [GeV]";

  if(name.find("muon_deta_")!=string::npos)xtitle="#Delta #eta (#mu,#mu)";
  if(name.find("el_deta_")!=string::npos)xtitle="#Delta #eta (e,e)";
  if(name.find("leaddimudeltaR_")!=string::npos)xtitle="#Delta R (#mu,#mu)";
  if(name.find("leaddieldeltaR_")!=string::npos)xtitle="#Delta R (e,e)";

  if(name.find("leaddijetma")!=string::npos)xtitle="m(j_{1}j_{2}) [GeV]";
  if(name.find("leaddijetdr")!=string::npos)xtitle="#Delta R(j_{1}j_{2})";
  if(name.find("jet_pt")!=string::npos)xtitle="jet p_{T} [GeV]";



  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());

  return;
}


bool HistInGev(string name){
  
  bool ingev=false;
  if(name.find("_pt_")!=string::npos)ingev=true;
  if(name.find("mass_")!=string::npos)ingev=true;
  
  return ingev;

}


float  GetMaximum(TH1* h_data, TH1* h_up, bool ylog, string name){

  float yscale= 1.4;
  if(!showdata) yscale = 1.2;
  if(ylog) yscale = 100000.;
  
  cout << name << endl;
  if(name.find("eta")!=string::npos) yscale*=1.5;
  if(name.find("MET")!=string::npos) yscale*=1.;
  if(name.find("charge")!=string::npos) yscale*=2.5;
  if(name.find("deltaR")!=string::npos) yscale*=2.;
  
  float max_data = h_data->GetMaximum()*yscale;
  float max_bkg = h_up->GetMaximum()*yscale;

  
  if(max_data > max_bkg) return max_data;
  else return max_bkg;
  
  return -1000.;
}



float GetTotal(TString cut, vector<pair<TString,float> > samples){
  
  float total(0.);
  for(vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++){
    total += GetIntegral(cut,(*it).first,"MC");
  }
  
  return total;    
}


float GetStatError2(TString cut, vector<pair<TString,float> > samples){  
  float err = GetStatError(cut,samples);
  err = err*err;
  return err;
}

float GetStatError(TString cut, vector<pair<TString,float> > samples){  

  TString path  = "/home/jalmond/LQanalyzer/data/output/Electron/";
  
  TFile* f0 =  TFile::Open((path+ fileprefix + samples.at(0).first +  filepostfix).Data());  
  TH1* h_tmp = dynamic_cast<TH1*> ((f0->Get(cut.Data()))->Clone());
  float stat_error(-99999.);

  for(unsigned int i=1; i < samples.size(); i++){
    
    
    TFile* f =  TFile::Open((path+ fileprefix + samples.at(i).first +  filepostfix).Data());
    TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data()))->Clone());
    h_tmp->Add(h);
    f->Close();
    if(i == (samples.size()- 1)){
      stat_error = Error(h_tmp);
    }
  }

  if(samples.size()==1) stat_error = Error(h_tmp);
  
  f0->Close();
  return stat_error;
  
}




float GetIntegral(TString cut, TString isample, TString type){

  TFile* f =  TFile::Open(( mcloc  + fileprefix + isample + filepostfix).Data());  
  
  
  if(!((f->Get(cut.Data())))){
    cout << "Histogram " << cut << " in "  << mcloc+  fileprefix + isample + filepostfix << " not found" << endl;
    exit(0);
  }

  TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data())->Clone()));
  
  
  if(!h) {
    cout << "Histogram " << cut << " in "  <<  (mcloc + fileprefix + isample + filepostfix) << " not found" << endl;
    exit(1);
  }
  
  float integral = h->Integral();
 
  
  if(!h) {
    f->Close();
    cout << "Systematic file does not exist. Setting this error to zero" << endl;
    type = "MC";
    path  = "/home/jalmond/LQanalyzer/data/output" + type + "/";

    TFile* f_tmp =  TFile::Open((mcloc+ fileprefix + isample +  filepostfix).Data());
    TH1* h_tmp = dynamic_cast<TH1*> ((f_tmp->Get(cut.Data())->Clone()));
    integral = h_tmp->Integral();
    f_tmp->Close();
    return integral;
  }
  

  f->Close();

  return integral;
  
}

float GetError(TString cut, TString isample, TString type){
  
  TString path  = mcloc + type + "/";
    
  TFile* f =  TFile::Open((path+ fileprefix + isample + filepostfix).Data());
  TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data())->Clone()));
  cout << h << endl;
  
  if(!h) {
    cout << "Histogram " << cut << " in "  << (path+ fileprefix + isample + filepostfix) << " not found" << endl;
    exit(1);
  }
  
  float err = Error(h);
  
  f->Close();

  return err;
  
}


float GetNormErr(TString cut, vector<pair<TString,float> > samples){
  
  float norm_err(0.);
  int i=0;
  for( vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++, i++){
    norm_err += GetIntegral(cut,it->first,"MC")* it->second*GetIntegral(cut,it->first,"MC")* it->second;
  }
  
  return sqrt(norm_err);
}

float GetNormErr2(TString cut,  vector<pair<TString,float> > samples){

  float err = GetNormErr(cut,samples);
  err = err*err;
  return err;
  
}



float GetErr2(TString cut, vector<pair<TString,float> > samples, TString err_type,TString var){
  
  float err = GetErr(cut,samples,err_type,var);
  err = err*err;
  return err;

}

float GetErr(TString cut, vector<pair<TString,float> > samples, TString err_type,TString var){

  float total(0.);
  float total_witherr(0.);
  
  int i=0;
  for(vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++, i++){
    total += GetIntegral(cut,it->first,("MC/"));
    total_witherr += GetIntegral(cut,it->first,("MC_" + err_type));
  }
  
  
  float err=  total - total_witherr;      
  

  if(err_type.Contains("MM")){    
    if((var.Contains("UP")||var.Contains("Up"))) {
      if(( ( total - total_witherr)< 0.)) return fabs(err);
      else return total*0.1;
    }
    if((var.Contains("DOWN")||var.Contains("Down")) && ( ( total - total_witherr > 0.))) return err;
    else return total*0.1;  
  }
  
  if((var.Contains("UP")||var.Contains("Up"))) {
    if(( ( total - total_witherr)< 0.)) return fabs(err);
    else return 0.;
  }
  if((var.Contains("DOWN")||var.Contains("Down"))){
    if( ( total - total_witherr > 0.)) return err;
    else return 0.;  
  }

    return err;

}



float GetSystPercent(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples ){

  return ( 100.*(GetSyst(cut, syst,samples )/ Calculate(cut,"Normal",samples)));
}



float GetSyst(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples ){
  
  if(syst.Contains("Stat"))return GetStatError(cut,samples.first) ;    
  if(syst.Contains("Normalisation")) return GetNormErr(cut,samples.first);

  if(syst.Contains("TOTAL"))return sqrt(  GetStatError2(cut,samples.first)+
					  GetNormErr2(cut,samples.first));  
  return -9999.;
}


float Calculate(TString cut, TString variance, pair<vector<pair<TString,float> >,TString > samples ){
  
  cout << cut << " " << variance << " " << samples.second << endl;
  
  if(samples.second.Contains("NonPrompt")){
    if(variance.Contains("Normal"))  return GetTotal(cut,samples.first) ;  
    if(variance.Contains("StatErr")) return GetStatError(cut,samples.first) ;  
    //if(variance.Contains("Up")) return sqrt(GetErr2(cut,samples.first,"MMFDown",variance)+GetErr2(cut,samples.first,"MMFUp",variance)+GetErr2(cut,samples.first,"MMRUp",variance)+GetErr2(cut,samples.first,"MMRDown",variance));
    //if(variance.Contains("Down")) return sqrt(GetErr2(cut,samples.first,"MMFDown",variance)+GetErr2(cut,samples.first,"MMFUp",variance)+GetErr2(cut,samples.first,"MMRUp",variance)+GetErr2(cut,samples.first,"MMRDown",variance) );     

  }
  
  
  if(variance.Contains("Up") || variance.Contains("UP")||variance.Contains("Down") || variance.Contains("DOWN")  ){
    bool debug =false; 
    if(debug){
      cout << "Variance = " << variance << endl;
      cout << "Normalisation =" <<  GetNormErr(cut,samples.first) << endl;
      cout << "JESUp = " <<  GetErr(cut,samples.first,"JESup",variance) << endl;
      cout << "JESDown = " <<  GetErr(cut,samples.first,"JESdown",variance) << endl;
      cout << "MUONISOUp = " <<  GetErr(cut,samples.first,"MUONISOup",variance) << endl;
      cout << "MUONISODown = " <<  GetErr(cut,samples.first,"MUONISOdown",variance) << endl;
      cout << "JVFup = " <<  GetErr(cut,samples.first,"JVFup",variance) << endl;
      cout << "JVFdown = " <<  GetErr(cut,samples.first,"JVFdown",variance) << endl;
      cout << "JER = " <<  GetErr(cut,samples.first,"JER",variance) << endl;
      cout << "CFup = " <<  GetErr(cut,samples.first,"CFup",variance) << endl;
      cout << "CFdown = " <<  GetErr(cut,samples.first,"CFdown",variance) << endl;
    }

    return fabs(sqrt( (GetNormErr2(cut,samples.first))));
    
  }
  
  if(variance.Contains("Normal")) return GetTotal(cut,samples.first) ;    
  
  if(variance.Contains("StatErr")) return GetStatError(cut,samples.first) ;    
    
  return -999999.;

}




float Error(TH1* h){
  double err ;
  double integral = h->IntegralAndError(0,h->GetNbinsX(),err,"");
  
  return err;
}



void SetUpMasterConfig(string name){
  
  // Get list of cuts to plot
  ifstream master_config_name_file(name.c_str());
  if(!master_config_name_file) {
    cerr << "Did not find " + name + ", exiting ..." << endl;
    return;
  }
  while(!master_config_name_file.eof()) {
    string tmp;
    string tmppath;
    master_config_name_file >> tmp;
    master_config_name_file >> tmppath;
    
    if(tmp=="END") break;
    if(tmp.find("#")!=string::npos) continue;

    if(tmp=="mcpath") mcloc = tmppath;
    if(tmp=="datapath") dataloc = tmppath;
    if(tmp=="datadrivenpath") datadrivenloc = tmppath;
    
    if(tmp=="prefix") fileprefix = tmppath;
    if(tmp=="postfix") filepostfix = tmppath;
    
    if(tmp=="plottingpath") plotloc = tmppath;
    if(tmp=="cutpath")  cutloc = tmppath;
    if(tmp=="cutpath")  cerr << "tmppath = " << tmppath << std::endl;

    if(tmp=="outputdir")    path = tmppath;

    
    if(tmp=="showdata")    {
      if (tmppath == "true") showdata=true;
      else showdata=false;
    }
    if(tmp=="ylog") {
      if (tmppath == "true")ylog=true;
      else ylog = false;
    }

    if(tmp=="usenp"){
      if (tmppath == "true")usenp = true;
      else usenp=false;
    }
    
    if(tmp=="samples"){
      listofsamples.push_back(tmppath);
    }
    
    if(tmp=="histdir") histdir = tmppath;
    
    cutfile = cutloc;
    histfile =  plotloc;

  }

}

void  SetUpConfig(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples, vector<string>& cut_label){
    
  cout << " /// MakeDataMCComplots::SetUpConfig " << endl;
  /// colours of histograms
  int tcol(0), zzcol(0), fcol(0), zcol(0), wzcol(0), sscol(0),  wwcol(0), wcol(0);
  
  // Get list of cuts to plot  
  ifstream colour_name_file("Config/colour.txt");
  if(!colour_name_file) {
    cerr << "Did not find Config/colour.txt, exiting ..." << endl;
    return;
  }
  while(!colour_name_file.eof()) {        
    string histname;
    int col;
    colour_name_file >> histname;    
    if(histname=="END") break;
    colour_name_file >> col;

    if(histname=="tcol") tcol =col;
    if(histname=="zzcol") zzcol =col;
    if(histname=="fcol") fcol =col;
    if(histname=="zcol") zcol =col;
    if(histname=="wzcol")wzcol =col;
    if(histname=="sscol") sscol =col;
    if(histname=="wwcol") wwcol =col;
    if(histname=="wcol") wcol =col;
    cout << "Set sample " << histname << " with colour " << col  << endl;
  }
  
  /// Setup list of samples: grouped into different processes 
  //// MC (truth only)
  vector<pair<TString,float> > top = InitSample("top");
  vector<pair<TString,float> > wz = InitSample("wz");
  vector<pair<TString,float> > zz = InitSample("zz");
  vector<pair<TString,float> > ww = InitSample("ww");
  vector<pair<TString,float> > ssww = InitSample("ssww");
  vector<pair<TString,float> > z = InitSample("DY");
  vector<pair<TString,float> > w = InitSample("wjet");
  
  /// NP is datadriven
  vector<pair<TString,float> > np;
  np.push_back(make_pair("datadriven",0.));
  
  for( unsigned int i = 0; i < listofsamples.size(); i++){
    if(listofsamples.at(i) =="WW")samples.push_back(make_pair(make_pair(ww,wwcol),"WW")); 
    if(listofsamples.at(i) =="ZZ")samples.push_back(make_pair(make_pair(zz,zzcol),"ZZ"));
    if(listofsamples.at(i) =="WZ")samples.push_back(make_pair(make_pair(wz,wzcol),"WZ"));
    if(listofsamples.at(i) =="DY")samples.push_back(make_pair(make_pair(z,zcol),"DY"));
    if(listofsamples.at(i) =="Top")samples.push_back(make_pair(make_pair(top,tcol),"Top"));
    if(listofsamples.at(i) =="Wjet")samples.push_back(make_pair(make_pair(w,wcol),"Wjet"));
    if(listofsamples.at(i) =="NonPrompt")samples.push_back(make_pair(make_pair(np,fcol),"NonPrompt"));   
  }

  ///// Fix cut flow code
  caption="";
  cut_label.push_back("DY");
  hist = "/mu1_eta";
  columnname="";

  
  cout << "Configured as :: " << endl;
  cout << "- channel = " << channel  << endl;
  cout << "- usenp = " << usenp << endl;
  cout << "- showdata = " << showdata << endl;
  cout << "- Y axis set log = " << ylog << endl;
  cout << "- cutfile = " << cutfile << endl;
  cout << "- histogram name = " << hist << endl;
  cout << "- Samples include: " << endl;

  for( vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it = samples.begin(); it!=samples.end();++it){
    cout << it->second << endl;
  }
  

  return;

}




TCanvas* CompDataMC(TH1* hdata, vector<THStack*> mcstack,TH1* hup, TH1* hdown,TLegend* legend, const string hname, const  int rebin, double xmin, double xmax,double ymin, double ymax,string path , string folder, bool logy, bool usedata, TString channel) {
  
  std::cout << "start " << std::endl;
  

  string cname;
  if(hdata) cname= string("c_") + hdata->GetName();
  else cname = string("c_") + ((TNamed*)mcstack.at(0)->GetHists()->First())->GetName();

  //Create Canvases
  TCanvas* canvas = new TCanvas((cname+"significance").c_str(), (cname+"significance").c_str(), 800, 600);

  
  std::string title=canvas->GetName();
  std::string tpdf = "/home/jalmond/WebPlots/"+ path + "/histograms/"+folder+"/"+title+".png";
  
  ///####################   Standard plot
  if(logy)canvas->SetLogy();
  canvas->cd();

  
  //// %%%%%%%%%% TOP HALF OF PLOT %%%%%%%%%%%%%%%%%%
  TH1* h_nominal = MakeSumHist2(mcstack.at(0));

  MakeLabel(0.2,0.8);

  TH1* errorband = MakeErrorBand(h_nominal,hup, hdown) ;

  SetNomBinError(h_nominal, hup, hdown);

  if(usedata){
    hdata->Draw("p");
    mcstack.at(0)->Draw("HIST same");
    hdata->Draw("p same");
    hdata->Draw("axis same");
    errorband->Draw("E2same");
  }
  else{
    errorband->GetXaxis()->SetRangeUser(xmin,xmax);
    errorband->GetYaxis()->SetRangeUser(ymin,ymax);
    errorband->Draw("E2");
    mcstack.at(0)->Draw("same HIST");
    errorband->Draw("E2same");

    
  }
  legend->Draw("same");

  if(usedata){
    //// %%%%%%%%%% BOTTOM (SIGNIFICANCE) HALF OF PLOT %%%%%%%%%%%%%%%%%%

    /// Make significance hist

    TH1* h_significance=(TH1F*)hdata->Clone();
    TH1* h_divup=(TH1F*)hup->Clone();
    TH1* h_divdown=(TH1F*)hdown->Clone();

    TH1* errorbandratio = (TH1*)h_nominal->Clone("AA");

    hdata->GetXaxis()->SetLabelSize(0.); ///
    hdata->GetXaxis()->SetTitle("");

    h_divup->Divide(h_nominal);
    h_divdown->Divide(h_nominal);

    for(int i=1; i < errorbandratio->GetNbinsX()+1; i++){

      float bc = ((h_divup->GetBinContent(i)+h_divdown->GetBinContent(i))/2.);
      float bd = ((h_divup->GetBinContent(i)-h_divdown->GetBinContent(i))/2.);

      errorbandratio->SetBinContent(i,bc);
      errorbandratio->SetBinError(i,bd);
    }

    errorbandratio->SetFillStyle(3354);
    errorbandratio->SetFillColor(kBlue-8);
    errorbandratio->SetMarkerStyle(0);

    for(int i=1; i < h_significance->GetNbinsX()+1; i++){
      float num = h_significance->GetBinContent(i) - h_nominal->GetBinContent(i);
      float denom = sqrt( (h_significance->GetBinError(i)*h_significance->GetBinError(i) + h_nominal->GetBinError(i)*h_nominal->GetBinError(i)));

      float sig = 0.;
      if(denom!=0.) sig = num / denom;

      /// For  now plot data/mc ...
      if(h_nominal->GetBinContent(i)!=0. ) sig =   h_significance->GetBinContent(i)/ h_nominal->GetBinContent(i);
      h_significance->SetBinContent(i,sig);
    }


    // How large fraction that will be taken up by the data/MC ratio part
    double FIGURE2_RATIO = 0.35;
    double SUBFIGURE_MARGIN = 0.15;
    canvas->SetBottomMargin(FIGURE2_RATIO);
    TPad *p = new TPad( "p_test", "", 0, 0, 1, 1.0 - SUBFIGURE_MARGIN, 0, 0, 0);  // create new pad, fullsize to have equal font-sizes in both plots
    p->SetTopMargin(1-FIGURE2_RATIO);   // top-boundary (should be 1 - thePad->GetBottomMargin() )
    p->SetFillStyle(0);     // needs to be transparent
    p->Draw();
    p->cd();
    
    //h_significance->SetFillColor(kGray+1);
    //h_significance->SetLineColor(kGray+1);

    h_significance->GetYaxis()->SetNdivisions(10204);
    //h_significance->GetYaxis()->SetTitle("Significance");
    h_significance->GetYaxis()->SetTitle("Data/MC");
    //h_significance->GetYaxis()->SetRangeUser(-4., 4.);
    h_significance->GetYaxis()->SetRangeUser(0., 2.);
    h_significance->GetXaxis()->SetRangeUser(xmin, xmax);
    h_significance->Draw("hist");
    TLine *line = new TLine(h_significance->GetBinLowEdge(h_significance->GetXaxis()->GetFirst()),1.0,h_significance->GetBinLowEdge(h_significance->GetXaxis()->GetLast()+1),1.0);
    
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
    h_significance->Draw("HISTsame");
    
  }
  
  canvas->Print(tpdf.c_str(), ".png");


  return canvas;



}



TH1* MakeSumHist2(THStack* thestack){

  TH1* hsum=0;
  TList* list = thestack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    if(!hsum) hsum = (TH1*)h->Clone( (string(h->GetName()) + "_sum").c_str() );
    else {
      hsum->Add(h, 1.0);
    }
  }//hist loop

  return hsum;
}



void SetNomBinError(TH1* hnom, TH1* hup, TH1* hdown){

  for(int i=1; i < hnom->GetNbinsX()+1; i++){

    float err1 = fabs(hnom->GetBinContent(i)- hup->GetBinContent(i));
    float err2 = fabs(hnom->GetBinContent(i)- hdown->GetBinContent(i));

    if(err1 > err2 ) hnom->SetBinError(i, err1);
    if(err2 > err1 ) hnom->SetBinError(i, err2);
  }
  return;
}



TH1* MakeErrorBand(TH1* hnom, TH1* hup, TH1* hdown){

  TH1* errorband = (TH1*)hnom->Clone("aa");

  for(int i=1; i < errorband->GetNbinsX()+1; i++){

    float bin_content = (hup->GetBinContent(i)+ hdown->GetBinContent(i))/2.;
    float bin_error = (hup->GetBinContent(i)- hdown->GetBinContent(i))/2.;

    errorband->SetBinContent(i,bin_content);
    errorband->SetBinError(i,bin_error);
  }

  errorband->SetFillStyle(3354);
  errorband->SetFillColor(kBlue-8);
  errorband->SetMarkerSize(0);
  errorband->SetMarkerStyle(0);
  errorband->SetLineColor(kWhite);
  errorband->Draw("E2Same");

  return errorband;

}


void MakeLabel(float rhcol_x, float rhcol_y){
  TLatex label;
  label.SetTextSize(0.04);
  label.SetTextColor(2);
  label.SetTextFont(42);
  label.SetNDC();
  label.SetTextColor(1);
  label.DrawLatex(rhcol_x,rhcol_y,"#int L dt = 20.4 fb^{-1}");
  label.DrawLatex(rhcol_x + 0.2,rhcol_y ,"#sqrt{s}= 8 TeV");
  label.SetTextSize(0.045);

  label.DrawLatex(rhcol_x+0.115, rhcol_y + 0.09,"Work In Progress");
  label.SetTextFont(72);
  label.DrawLatex(rhcol_x, rhcol_y + 0.09,"CMS");

  return;
}


