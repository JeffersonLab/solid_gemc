#ifndef __ECAL_HH
#define __ECAL_HH
#include <functional>
//reference - http://stackoverflow.com/questions/3973659/c-unordered-map-compiling-issue-with-g
#include <tr1/unordered_map>

#include <map>
#include <cassert>

#include <TVector2.h>
#include <TMath.h>

#define __IO_MAXHIT 10000

using namespace std;
using namespace std::tr1;

//ecal geometry (m) 
Double_t blockside = 6.26/100;//6.25
Double_t rmax = 265.0/100;
Double_t rmin = 110.0/100;

Double_t xspace = blockside*3;//18.75 cm
Double_t yspace = blockside*sqrt(3.0);// 10.825 cm



struct hashFunc{
  size_t operator()(const TVector2 &k) const{
    size_t h1 = hash<double>()(k.Px());
    size_t h2 = hash<double>()(k.Py());

    //Cantor pairing function:
    size_t A = h1 >= 0 ? 2 * h1 : -2 * h1 - 1;
    size_t B = h2 >= 0 ? 2 * h2 : -2 * h2 - 1;

    return ((A + B) * (A + B + 1) / 2 + A);
  }
};


struct equalsFunc{
  bool operator()( const TVector2& lhs, const TVector2& rhs ) const{
    return ((lhs.Px() == rhs.Px()) && (lhs.Py() == rhs.Py()));
    //return ( (abs(lhs.Px() - rhs.Px()) < 0.02) && (abs(lhs.Py() - rhs.Py()) < 0.02) );
  }
};


typedef unordered_map<TVector2, double, hashFunc, equalsFunc> TPoint2DMap;

struct singledethit{
  Int_t fGenDetHit_det;
  Int_t fGenDetHit_trid;
  Int_t fGenDetHit_mtrid;
  Int_t fGenDetHit_pid;
  Double_t fGenDetHit_P;
  Double_t fGenDetHit_PZ;
  Double_t fGenDetHit_X;
  Double_t fGenDetHit_Y;
  Double_t fGenDetHit_Z;
  Double_t fGenDetHit_VZ;
  Double_t fGenDetHit_T;
};

struct singledetsum{
  Int_t fCalDetSum_det;
  Double_t fCalDetSum_edep;
  Double_t fCalDetSum_x;
  Double_t fCalDetSum_y;
  Double_t fCalDetPos_X;
  Double_t fCalDetPos_Y;
  Double_t fCalDetPos_Z;
  Double_t fCalDet_ST;
  Double_t fCalDet_ET;  
};


struct singleevent{
  Int_t fEv_num;
  Int_t fNGenDetHit;
  Int_t fNCalDetSum;
  vector< singledethit > hitlist;
  vector< singledetsum > sumlist;
};

#endif//__ECAL_HH

/*
Hash function for two input parameters,

Cantor pairing function: (reference http://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way)

A = a >= 0 ? 2 * a : -2 * a - 1;
B = b >= 0 ? 2 * b : -2 * b - 1;
(A + B) * (A + B + 1) / 2 + A;




 */
