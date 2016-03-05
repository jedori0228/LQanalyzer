#include "Trilepton.h"

using namespace std;

///////////////////////////////////////
/// Put Pz to Neutirno's 4-vector ///
///////////////////////////////////////
void PutNuPz(TLorentzVector *nu, double Pz){
  double Px, Py;
  Px = nu->Px();
  Py = nu->Py();
  nu->SetPxPyPzE(Px, Py, Pz, TMath::Sqrt(Px*Px+Py*Py+Pz*Pz));
}


////////////////////////////////
/// Solve quadratic equation ///
////////////////////////////////
double solveqdeq(bool *isNegative, double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm){
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
    *isNegative = true;
    return b/(2*a);
  }
  else{
    *isNegative = false;
    if(pm=="p") return (b+TMath::Sqrt(b*b-4*a*c))/(2*a);
    else if(pm=="m")  return (b-TMath::Sqrt(b*b-4*a*c))/(2*a);
    else{
      cout << "wrong plus/minus option" << endl;
      return 0;
    }
  }
}

////////////////////////////////
/// Solve quadratic equation ///
////////////////////////////////
double solveqdeq(double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm){
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

/////////////////////
/// SetNeutrinoPz ///
/////////////////////
void SetNeutrinoPz(TLorentzVector *nu, double Pz){
  double MET = TMath::Sqrt((*nu).Px()*(*nu).Px() + (*nu).Py()*(*nu).Py());
  (*nu).SetPxPyPzE((*nu).Px(), (*nu).Py(), Pz, TMath::Sqrt(MET*MET+Pz*Pz));

}

/////////////////////////
/// Angle in XY plane ///
/////////////////////////
double XYAngle(TLorentzVector a, TLorentzVector b){
  TLorentzVector a_temp, b_temp;
  a_temp.SetPxPyPzE(a.Px(), a.Py(), 0, 0);
  b_temp.SetPxPyPzE(b.Px(), b.Py(), 0, 0);
  return a_temp.Angle(b_temp.Vect());
}
