#ifndef Trilepton_h
#define Trilepton_h

#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <TMath.h>

///////////////////////////////////////
/// Put Pz to Neutirno's 4-vector ///
///////////////////////////////////////
void PutNuPz(TLorentzVector *nu, double Pz);


////////////////////////////////
/// Solve quadratic equation ///
////////////////////////////////
double solveqdeq(bool *isNegative, double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm);

////////////////////////////////
/// Solve quadratic equation ///
////////////////////////////////
double solveqdeq(double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm);

/////////////////////
/// SetNeutrinoPz ///
/////////////////////
void SetNeutrinoPz(TLorentzVector *nu, double Pz);

/////////////////////////
/// Angle in XY plane ///
/////////////////////////
double XYAngle(TLorentzVector a, TLorentzVector b);

int find_mlmet_closest_to_W(TLorentzVector* lep, TLorentzVector MET);

double Mt3(TLorentzVector a, TLorentzVector b);

#endif
