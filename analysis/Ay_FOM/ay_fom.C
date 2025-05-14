#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void ay_fom()
{

gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
//   gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);  
  gStyle->SetPadRightMargin(0.15);  
  
  gStyle->SetPadColor(0);

  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.04,"xyz");   
  gStyle->SetTitleSize(0.07,"t");    
//   gStyle->SetPaintTextFormat("4.1f");       
    
//counts estimation based the following
//15uA e- beam on 40cm 10amg He3 target, approved 48 PAC days
//this study uses e- crosssection from Eric Christy and Peter Boosted 2021 fit, even though rate on He3 still need improvement (eAll at https://github.com/JeffersonLab/evgen_inclusive_e). proposal use e- crosssection from pdf fit (eDIS at https://github.com/JeffersonLab/evgen_inclusive) which is a few percent lower than eAll
//this study used acceptance from SoLID large angle detector with window collimator version 201701 at https://github.com/JeffersonLab/solid_gemc/tree/master/analysis/acceptance/result_SIDIS_He3, proposal used an older version 201402 which is a few percent higher
//cut Q2>1 GeV2, W>2 GeV and P cut

//count in each Q2 bin
const int n=8;    
char *Q2_center[n]={"2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5"};
// double count[n]={1.1761e+06,1.32249e+09,2.96848e+09,2.43187e+09,1.15184e+09,3.77492e+08,6.77e+07,  1.27412e+06};   //11GeV beam, P>3.5GeV, rate=total/3600/24/48=2khz
// double count[n]={1.1761e+06+5.51285e+08,1.32249e+09+2.35171e+09,2.96848e+09+ 1.00804e+09,2.43187e+09+2.1516e+08,1.15184e+09+5.88678e+06,3.77492e+08,6.77e+07,1.27412e+06};   //11GeV+8.8GeV beam, P>3.5GeV
// double count[n]={2.97546e+08,3.15867e+09,4.01386e+09,2.98397e+09,1.50698e+09,5.18632e+08,1.38184e+08,2.57073e+07};  //11GeV beam, P>3.0GeV, rate=total/3600/24/48=2.9khz
double count[n]={2.97546e+08+1.87335e+09,3.15867e+09+3.39065e+09,4.01386e+09+1.43267e+09,2.98397e+09+2.79171e+08,1.50698e+09,5.18632e+08+6.00997e+06,1.38184e+08,2.57073e+07};  //11GeV+8.8GeV beam, P>3.0GeV, rate=total/3600/24/48=2.9khz

//Assumption
double Ay=-0.03; // close to PRL data
// double Ay=-0.01; // same assumption as in proposal
double Pt=0.6; //He3 polarization
double Pn=0.86; //effective neutron polarization
double dhe3=0.85; //He3 dilution
double dn=0.2; //neutron dilution in He3

TCanvas *c_Aut = new TCanvas("c_Aut","c_Aut",1900,800);
c_Aut->Divide(n/2,2);

TH1F* h_Ay=new TH1F("Ay","SoLID He3 projection;Q^{2} (GeV^{2});A_{y} (neutron)",10,0,10);

TH1F* h_Aut[n];
for(int i=0;i<n;i++){
   char hstname[100];   
   sprintf(hstname,"Aut_%i",i);
   h_Aut[i]=new TH1F(hstname,Form("Q^{2}=%sGeV^{2};#phi (rad);A_{UT}",Q2_center[i]),12,0,3.1416*2);
   
   for(int j=0;j<12;j++){
    double Aut=Ay*sin(3.1416*2/12*(j+0.5)); // with simple sin form
    double Aut_error_stat=1/sqrt(count[i]/12.)/Pt/Pn/dhe3/dn; //abs stat error, average in 12 phi_s bins, note SoLID acceptance is roughly symmetric in phi_s with some theta dependence.
    double Aut_error_sys=0.07*Aut+1e-4+1e-3; //abs sys error according to proposal
    double Aut_error_total=Aut_error_stat; // not using sys error yet
    double Aut_smeared=gRandom->Gaus(Aut,Aut_error_total); // smearing to prevent overfitting
   
//     h_Aut[i]->SetBinContent(j+1,Aut); 
    h_Aut[i]->SetBinContent(j+1,Aut_smeared);
    h_Aut[i]->SetBinError(j+1,Aut_error_total); 
   }
   
   c_Aut->cd(i+1);
   h_Aut[i]->SetMaximum(Ay+abs(Ay)*3);
   h_Aut[i]->SetMinimum(Ay-abs(Ay)*1);   
   h_Aut[i]->Draw("E1 X0");
   
   TF1 *f = new TF1("f", "[0]*sin(x)", 0, 10);
   f->SetParameter(0,Ay); 
   h_Aut[i]->Fit(f);
   double Ay_error=f->GetParError(0);
   double Ay_fit=f->GetParameter(0);       
   
   TPaveText *pt = new TPaveText(.5, .8, .95, .9,"NDC");
   pt->AddText(Form("Ay %f #pm %f",Ay_fit,Ay_error));
   pt->Draw();
   
   if (i!=7){
    h_Ay->SetBinContent(i+3,Ay);
    h_Ay->SetBinError(i+3,Ay_error);
    cout << Ay_error/Ay << endl;
   }
}

TCanvas *c_Ay = new TCanvas("c_Ay","c_Ay",1000,800);
   
h_Ay->SetMaximum(Ay+0.035);
h_Ay->SetMinimum(Ay-0.035);
// h_Ay->SetMaximum(Ay+0.004);
// h_Ay->SetMinimum(Ay-0.004);
h_Ay->SetMarkerColor(kRed);
h_Ay->SetLineColor(kRed);
h_Ay->SetMarkerStyle(21);
h_Ay->SetMarkerSize(0.5);
h_Ay->Draw("E1 X0");

// He3 11GeV, Fig 33 and Fig 37 of proposal
// double Ay_error_proposal[9]={1.47e-5,2.58e-5,4.88e-5,8.57e-5,1.27e-4,2.04e-4,4.11e-4,1.07e-3,8.68e-3};
double Ay_error_proposal[8]={1.47e-5,2.58e-5,4.88e-5,8.57e-5,1.27e-4,2.04e-4,4.11e-4,1.07e-3};  
TH1F* h_Ay_proposal=new TH1F("Ay_proposal","",10,0.1,10.1);
for(int i=0;i<8;i++){
h_Ay_proposal->SetBinContent(i+2,Ay);
h_Ay_proposal->SetBinError(i+2,Ay_error_proposal[i]);
}
h_Ay_proposal->SetMarkerStyle(21);
h_Ay_proposal->SetMarkerSize(0.5);
h_Ay_proposal->SetMarkerColor(kBlack);
h_Ay_proposal->SetLineColor(kBlack);
// h_Ay_proposal->Draw("same E1 X0");

const int m=5;
Double_t x[m]  = {3.24,2.56,2.08,1.58,1.05};
Double_t y[m]  = {-3.87e-2,-3.89e-2,-1.08e-2,-3.84e-2,-0.64e-2};  
Double_t ex[m] = {0,0,0,0};
Double_t ey[m] = {1.55e-2,0.96e-2,1.18e-2,2e-2,0.41e-2};
TGraphErrors *gr = new TGraphErrors(m,x,y,ex,ey);
gr->SetMarkerColor(kBlue);
gr->SetLineColor(kBlue);
gr->SetMarkerStyle(21);
gr->SetMarkerSize(0.5);
gr->Draw("P");

TLegend* leg = new TLegend(0.6, 0.8, 0.9, 0.9);
leg->AddEntry(h_Ay,Form("SoLID LD"),"l");
// leg->AddEntry(h_Ay_proposal,Form("FD+LD in proposal"),"l");
leg->AddEntry(gr,Form("HallA data (PRL 2014)"),"l");
leg->Draw();

}
