#include <string>
#include <iterator>
#include <algorithm>
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <stdio.h>
#include <vector>
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void plot_parameter()
{
gROOT->Reset();
gStyle->SetPalette(1);
// gStyle->SetOptStat(0);
gStyle->SetOptStat(111111);

const double DEG=180./3.1415926;

const int n=41;

double PhotonEnergy[n]={
2.04358, 2.0664, 2.09046, 2.14023, 
2.16601, 2.20587, 2.23327, 2.26137, 
2.31972, 2.35005, 2.38116, 2.41313, 
2.44598, 2.47968, 2.53081, 2.58354, 
2.6194, 2.69589, 2.73515, 2.79685, 
2.86139, 2.95271, 3.04884, 3.12665, 
3.2393, 3.39218, 3.52508, 3.66893,
3.82396, 3.99949, 4.13281, 4.27679, 
4.48244, 4.65057, 4.89476, 5.02774, 
5.16816, 5.31437, 5.63821, 5.90401, 
6.19921
};  // in ev
double QE_H8500_03[n] = {
0.008, 0.0124, 0.0157, 0.02125, 
0.0275, 0.034, 0.04, 0.048, 
0.062, 0.0753, 0.09, 0.1071, 
0.12144, 0.1428, 0.15, 0.16429, 
0.17857, 0.1928, 0.2, 0.2125,
0.225, 0.2375, 0.25, 0.2625, 
0.275, 0.275, 0.275, 0.275, 
0.275, 0.275, 0.2625, 0.25, 
0.2375, 0.2125, 0.192859, 0.185716, 
0.178573, 0.15714, 0.13572, 0.1143,
0.09  
}; 
double QE_H12700_03[n] = {
0.016,0.02,0.025,0.033,
0.042,0.048,0.056,0.06,
0.075,0.085,0.096,0.121,
0.147,0.166,0.182,0.194,
0.203,0.22,0.238,0.253,
0.269,0.287,0.3,0.31,
0.32,0.33,0.335,0.335,
0.335,0.33,0.325,0.31,
0.296,0.282,0.257,0.237,
0.22,0.197,0.165,0.139,
0.114
};
double QE_H12700_03_WLS_meas[n] = {
0.016,0.02,0.0243455,0.0349796,0.0400769,0.0495496,0.054666,0.0612895,0.0758019,0.0853365,0.100662,0.121331,0.144678,0.162644,0.180719,0.194414,0.202599,0.224051,0.235051,0.253334,0.268143,0.285398,0.30002,0.309013,0.319247,0.328839,0.333333,0.335,0.33337,0.327161,0.321697,0.328776,0.333637,0.318123,0.313051,0.326953,0.331335,0.331335,0.331335,0.331335,0.331335
};

double Wavelength[n];
for(int i=0;i<n;i++){
Wavelength[i]=1.24/PhotonEnergy[i]*1e3; // in nm
QE_H12700_03[i]=QE_H12700_03[i]*100;
QE_H12700_03_WLS_meas[i]=QE_H12700_03_WLS_meas[i]*100;
QE_H8500_03[i]=QE_H8500_03[i]*100;
}

TGraph *g_QE_H12700_03_E=new TGraph(n,PhotonEnergy,QE_H12700_03);
g_QE_H12700_03_E->SetTitle(";E(eV);QE(%)");
g_QE_H12700_03_E->SetName("H12700_03");
TGraph *g_QE_H12700_03_WL=new TGraph(n,Wavelength,QE_H12700_03);
g_QE_H12700_03_WL->SetTitle(";wavelength(nm);QE(%)");
g_QE_H12700_03_WL->SetName("H12700_03");

TGraph *g_QE_H8500_03_E=new TGraph(n,PhotonEnergy,QE_H8500_03);
g_QE_H8500_03_E->SetTitle(";E(eV);QE(%)");
g_QE_H8500_03_E->SetName("QE_H8500_03");
TGraph *g_QE_H8500_03_WL=new TGraph(n,Wavelength,QE_H8500_03);
g_QE_H8500_03_WL->SetTitle(";wavelength(nm);QE(%)");
g_QE_H8500_03_WL->SetName("QE_H8500_03");

TGraph *g_QE_H12700_03_WLS_meas_E=new TGraph(n,PhotonEnergy,QE_H12700_03_WLS_meas);
g_QE_H12700_03_WLS_meas_E->SetTitle(";E(eV);QE(%)");
g_QE_H12700_03_WLS_meas_E->SetName("H12700_03_WLS_meas");
TGraph *g_QE_H12700_03_WLS_meas_WL=new TGraph(n,Wavelength,QE_H12700_03_WLS_meas);
g_QE_H12700_03_WLS_meas_WL->SetTitle(";wavelength(nm);QE(%)");
g_QE_H12700_03_WLS_meas_WL->SetName("H12700_03_WLS_meas");

double Wavelength_Hamamatsu[51] = {
200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700
};

double QE_H12700_03_Hamamatsu[51] = {
11.50,14.12,16.83,18.91,22.45,25.18,27.22,28.86,30.15,31.35,32.87,33.37,33.63,33.73,33.46,33.12,32.74,32.42,32.45,31.54,30.76,29.76,28.77,27.59,26.16,24.31,22.48,20.93,19.67,18.49,16.50,13.05,10.11,8.36,7.15,6.12,5.19,4.29,3.46,2.69,2.00,1.41,0.95,0.60,0.35,0.19,0.10,0.05,0.02,0.01,0.00  
};

TGraph *g_QE_H12700_03_Hamamatsu_WL=new TGraph(51,Wavelength_Hamamatsu,QE_H12700_03_Hamamatsu);
g_QE_H12700_03_Hamamatsu_WL->SetTitle(";wavelength(nm);QE(%)");
g_QE_H12700_03_Hamamatsu_WL->SetName("H12700_03_Hamamatsu");

// LGC optical parameters from ../../geometry/lgc_moved/solid_lgc_OptMats.pl and solid_lgc_mirror.pl
// LGC materials use the HGC photon energy table plus two more points
const int n_LGC=43;
double PhotonEnergy_LGC[n_LGC];
for(int i=0;i<n;i++) PhotonEnergy_LGC[i]=PhotonEnergy[i];
PhotonEnergy_LGC[41]=6.5254848;
PhotonEnergy_LGC[42]=6.8927794367;

double QE_LGC_H12700[n_LGC] = {  // SL_H12700 "efficiency"
0.016, 0.02, 0.024, 0.035,
0.040, 0.050, 0.055, 0.061,
0.077, 0.085, 0.101, 0.121,
0.145, 0.163, 0.181, 0.194,
0.203, 0.224, 0.235, 0.253,
0.268, 0.285, 0.300, 0.309,
0.319, 0.329, 0.333, 0.335,
0.333, 0.327, 0.322, 0.329,
0.334, 0.318, 0.313, 0.327,
0.331, 0.331, 0.331, 0.331,
0.331, 0.331, 0.331
};

double Wavelength_LGC[n_LGC];
for(int i=0;i<n_LGC;i++){
Wavelength_LGC[i]=1.24/PhotonEnergy_LGC[i]*1e3; // in nm
QE_LGC_H12700[i]=QE_LGC_H12700[i]*100;
}

TGraph *g_QE_LGC_H12700_E=new TGraph(n_LGC,PhotonEnergy_LGC,QE_LGC_H12700);
g_QE_LGC_H12700_E->SetTitle(";E(eV);QE(%)");
g_QE_LGC_H12700_E->SetName("LGC_H12700");
TGraph *g_QE_LGC_H12700_WL=new TGraph(n_LGC,Wavelength_LGC,QE_LGC_H12700);
g_QE_LGC_H12700_WL->SetTitle(";wavelength(nm);QE(%)");
g_QE_LGC_H12700_WL->SetName("LGC_H12700");


TLegend* leg_QE_E = new TLegend(0.8, 0.8, .95, .95);
leg_QE_E->AddEntry(g_QE_H12700_03_E,g_QE_H12700_03_E->GetName());
leg_QE_E->AddEntry(g_QE_H12700_03_WLS_meas_E,g_QE_H12700_03_WLS_meas_E->GetName());
leg_QE_E->AddEntry(g_QE_H8500_03_E,g_QE_H8500_03_E->GetName());
leg_QE_E->AddEntry(g_QE_LGC_H12700_E,g_QE_LGC_H12700_E->GetName());

TLegend* leg_QE_WL = new TLegend(0.8, 0.8, .95, .95);
leg_QE_WL->AddEntry(g_QE_H12700_03_WL,g_QE_H12700_03_WL->GetName());
leg_QE_WL->AddEntry(g_QE_H12700_03_WLS_meas_WL,g_QE_H12700_03_WLS_meas_WL->GetName());
leg_QE_WL->AddEntry(g_QE_H8500_03_WL,g_QE_H8500_03_WL->GetName());
leg_QE_WL->AddEntry(g_QE_H12700_03_Hamamatsu_WL,g_QE_H12700_03_Hamamatsu_WL->GetName());
leg_QE_WL->AddEntry(g_QE_LGC_H12700_WL,g_QE_LGC_H12700_WL->GetName());


TCanvas *c_QE = new TCanvas("QE","QE",1600,800);
c_QE->Divide(2,1);
c_QE->cd(1);
// gPad->SetLogy(1);
gPad->SetGrid(1);
g_QE_H12700_03_E->SetMaximum(35);
g_QE_H12700_03_E->SetMinimum(0);
g_QE_H12700_03_E->GetXaxis()->SetRangeUser(2,7);
g_QE_H12700_03_E->SetMarkerColor(kBlue);
g_QE_H12700_03_E->Draw("AC*");
g_QE_H12700_03_WLS_meas_E->SetMarkerColor(kGreen);
g_QE_H12700_03_WLS_meas_E->Draw("C*");
g_QE_H8500_03_E->SetMarkerColor(kRed);
g_QE_H8500_03_E->Draw("C*");
g_QE_LGC_H12700_E->SetMarkerColor(kMagenta);
g_QE_LGC_H12700_E->SetMarkerStyle(24);
g_QE_LGC_H12700_E->Draw("P");
leg_QE_E->Draw();
c_QE->cd(2);
// gPad->SetLogy(1);
gPad->SetGrid(1);
g_QE_H12700_03_WL->SetMaximum(35);
g_QE_H12700_03_WL->SetMinimum(0);
g_QE_H12700_03_WL->GetXaxis()->SetRangeUser(200,800);
g_QE_H12700_03_WL->SetMarkerColor(kBlue);
g_QE_H12700_03_WL->Draw("AC*");
g_QE_H12700_03_WLS_meas_WL->SetMarkerColor(kGreen);
g_QE_H12700_03_WLS_meas_WL->Draw("C*");
g_QE_H8500_03_WL->SetMarkerColor(kRed);
g_QE_H8500_03_WL->Draw("C*");
g_QE_H12700_03_Hamamatsu_WL->SetMarkerColor(kBlack);
g_QE_H12700_03_Hamamatsu_WL->Draw("C*");
g_QE_LGC_H12700_WL->SetMarkerColor(kMagenta);
g_QE_LGC_H12700_WL->SetMarkerStyle(24);
g_QE_LGC_H12700_WL->Draw("P");
leg_QE_WL->Draw();


double reflectivity[n] = {
0.8678125, 0.8651562, 0.8639063, 0.8637500,
0.8640625, 0.8645313, 0.8643750, 0.8656250,
0.8653125, 0.8650000, 0.8648437, 0.8638281, 
0.8635156, 0.8631250, 0.8621875, 0.8617188,
0.8613281, 0.8610156, 0.8610938, 0.8616016,
0.8623047, 0.8637500, 0.8655859, 0.8673828,
0.8700586, 0.8741992, 0.8781055, 0.8825195,
0.8876172, 0.8937207, 0.8981836, 0.9027441,
0.9078369, 0.9102002, 0.9093164, 0.9061743,
0.9004223, 0.8915210, 0.8599536, 0.8208313,
0.7625024
};

TGraph *g_reflectivity_E=new TGraph(n,PhotonEnergy,reflectivity);
g_reflectivity_E->SetTitle(";E(eV);reflectivity");
TGraph *g_reflectivity_WL=new TGraph(n,Wavelength,reflectivity);
g_reflectivity_WL->SetTitle(";wavelength(nm);reflectivity");
g_reflectivity_E->SetName("reflectivity");
g_reflectivity_WL->SetName("reflectivity");

double reflectivity_LGC[n] = {  // LGC_Mirror and LGC_WinstonCone "reflectivity"
0.8678125, 0.8651562, 0.8639063, 0.8637500,
0.8640625, 0.8645313, 0.8643750, 0.8656250,
0.8653125, 0.8650000, 0.8648437, 0.8638281,
0.8635156, 0.8631250, 0.8621875, 0.8617188,
0.8613281, 0.8610156, 0.8610938, 0.8616016,
0.8623047, 0.8637500, 0.8655859, 0.8673828,
0.8700586, 0.8741992, 0.8781055, 0.8825195,
0.8876172, 0.8937207, 0.8981836, 0.9027441,
0.9078369, 0.9102002, 0.9093164, 0.9061743,
0.9004223, 0.8915210, 0.8599536, 0.8208313,
0.7625024
};

TGraph *g_reflectivity_LGC_E=new TGraph(n,PhotonEnergy,reflectivity_LGC);
g_reflectivity_LGC_E->SetName("LGC mirror, Winston cone");
TGraph *g_reflectivity_LGC_WL=new TGraph(n,Wavelength,reflectivity_LGC);
g_reflectivity_LGC_WL->SetName("LGC mirror, Winston cone");

TLegend* leg_reflectivity_E = new TLegend(0.6, 0.8, .95, .95);
leg_reflectivity_E->AddEntry(g_reflectivity_E,g_reflectivity_E->GetName());
leg_reflectivity_E->AddEntry(g_reflectivity_LGC_E,g_reflectivity_LGC_E->GetName());

TLegend* leg_reflectivity_WL = new TLegend(0.6, 0.8, .95, .95);
leg_reflectivity_WL->AddEntry(g_reflectivity_WL,g_reflectivity_WL->GetName());
leg_reflectivity_WL->AddEntry(g_reflectivity_LGC_WL,g_reflectivity_LGC_WL->GetName());

TCanvas *c_reflectivity = new TCanvas("reflectivity","reflectivity",1600,800);
c_reflectivity->Divide(2,1);
c_reflectivity->cd(1);
g_reflectivity_E->Draw("AC*");
g_reflectivity_E->SetMaximum(1);
g_reflectivity_E->SetMinimum(0.5);
g_reflectivity_E->GetXaxis()->SetRangeUser(2,7);
g_reflectivity_LGC_E->SetMarkerColor(kMagenta);
g_reflectivity_LGC_E->SetMarkerStyle(24);
g_reflectivity_LGC_E->Draw("P");
leg_reflectivity_E->Draw();
c_reflectivity->cd(2);
g_reflectivity_WL->Draw("AC*");
g_reflectivity_WL->SetMaximum(1);
g_reflectivity_WL->SetMinimum(0.5);
g_reflectivity_WL->GetXaxis()->SetRangeUser(200,800);
g_reflectivity_LGC_WL->SetMarkerColor(kMagenta);
g_reflectivity_LGC_WL->SetMarkerStyle(24);
g_reflectivity_LGC_WL->Draw("P");
leg_reflectivity_WL->Draw();

double RefractiveIndex_C4F8O_oneandhalfatm[n] = {
1.00205, 1.00205, 1.00205, 1.00206,
1.00206, 1.00206, 1.00206, 1.00206,
1.00206, 1.00206, 1.00206, 1.00207,
1.00207, 1.00207, 1.00207, 1.00207,
1.00207, 1.00208, 1.00208, 1.00208,
1.00209, 1.00209, 1.00209, 1.0021,
1.00211, 1.00211, 1.00212, 1.00213,
1.00214, 1.00215, 1.00216, 1.00217,
1.00219, 1.0022, 1.00222, 1.00223,
1.00224, 1.00225, 1.00228, 1.00231,
1.00234
}; //at oneandhalfatm

double RefractiveIndex_C4F8O_oneatm[n] = {
(1.00205-1)/1.5+1, (1.00205-1)/1.5+1, (1.00205-1)/1.5+1, (1.00206-1)/1.5+1,
(1.00206-1)/1.5+1, (1.00206-1)/1.5+1, (1.00206-1)/1.5+1, (1.00206-1)/1.5+1,
(1.00206-1)/1.5+1, (1.00206-1)/1.5+1, (1.00206-1)/1.5+1, (1.00207-1)/1.5+1,
(1.00207-1)/1.5+1, (1.00207-1)/1.5+1, (1.00207-1)/1.5+1, (1.00207-1)/1.5+1,
(1.00207-1)/1.5+1, (1.00208-1)/1.5+1, (1.00208-1)/1.5+1, (1.00208-1)/1.5+1,
(1.00209-1)/1.5+1, (1.00209-1)/1.5+1, (1.00209-1)/1.5+1, (1.0021-1)/1.5+1,
(1.00211-1)/1.5+1, (1.00211-1)/1.5+1, (1.00212-1)/1.5+1, (1.00213-1)/1.5+1,
(1.00214-1)/1.5+1, (1.00215-1)/1.5+1, (1.00216-1)/1.5+1, (1.00217-1)/1.5+1,
(1.00219-1)/1.5+1, (1.0022-1)/1.5+1, (1.00222-1)/1.5+1, (1.00223-1)/1.5+1,
(1.00224-1)/1.5+1, (1.00225-1)/1.5+1, (1.00228-1)/1.5+1, (1.00231-1)/1.5+1,
(1.00234-1)/1.5+1
}; //at 1atm

double RefractiveIndex_C4F10_oneandhalfatm[n] = { 
// 1.0012258504, 1.001225763, 1.0012265232, 1.0012264533, 
(1.0012272134-1)*1.5+1, (1.0012288123-1)*1.5+1, (1.001232115-1)*1.5+1, (1.0012294938-1)*1.5+1, 
(1.001230254-1)*1.5+1, (1.0012310141-1)*1.5+1, (1.0012309355-1)*1.5+1, (1.0012316956-1)*1.5+1, 
(1.0012324557-1)*1.5+1, (1.0012349196-1)*1.5+1, (1.0012356798-1)*1.5+1, (1.0012364486-1)*1.5+1, 
(1.0012363525-1)*1.5+1, (1.0012379689-1)*1.5+1, (1.0012395853-1)*1.5+1, (1.0012411842-1)*1.5+1,
(1.0012436481-1)*1.5+1, (1.0012444169-1)*1.5+1, (1.0012477021-1)*1.5+1, (1.001248471-1)*1.5+1, 
(1.0012492399-1)*1.5+1, (1.001251695-1)*1.5+1, (1.0012533027-1)*1.5+1, (1.001256614-1)*1.5+1, 
(1.0012599254-1)*1.5+1, (1.0012623719-1)*1.5+1, (1.0012656745-1)*1.5+1, (1.0012698334-1)*1.5+1, 
(1.0012731361-1)*1.5+1, (1.0012781425-1)*1.5+1, (1.0012848264-1)*1.5+1, (1.0012889853-1)*1.5+1, 
(1.001295678-1)*1.5+1, (1.0013032182-1)*1.5+1, (1.0013099196-1)*1.5+1, (1.0013174685-1)*1.5+1,
(1.0013292375-1)*1.5+1, (1.0013401678-1)*1.5+1, (1.0013553268-1)*1.5+1, (1.0013713595-1)*1.5+1,
(1.0013941548-1)*1.5+1, 
// 1.0014118738, 1.0014465255, 1.0014761009, 1.0015183889  
};

double RefractiveIndex_C4F10_oneatm[n] = { 
// 1.0012258504, 1.001225763, 1.0012265232, 1.0012264533, 
1.0012272134, 1.0012288123, 1.001232115, 1.0012294938, 
1.001230254, 1.0012310141, 1.0012309355, 1.0012316956, 
1.0012324557, 1.0012349196, 1.0012356798,1.0012364486, 
1.0012363525, 1.0012379689, 1.0012395853, 1.0012411842,
1.0012436481, 1.0012444169, 1.0012477021, 1.001248471, 
1.0012492399, 1.001251695, 1.0012533027, 1.001256614, 
1.0012599254, 1.0012623719, 1.0012656745, 1.0012698334, 
1.0012731361, 1.0012781425, 1.0012848264, 1.0012889853, 
1.001295678, 1.0013032182, 1.0013099196, 1.0013174685,
1.0013292375, 1.0013401678, 1.0013553268, 1.0013713595,
1.0013941548, 
// 1.0014118738, 1.0014465255, 1.0014761009, 1.0015183889  
};

TGraph *g_RefractiveIndex_C4F8O_oneandhalfatm_E=new TGraph(n,PhotonEnergy,RefractiveIndex_C4F8O_oneandhalfatm);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->SetTitle(";E(eV);RefractiveIndex");
g_RefractiveIndex_C4F8O_oneandhalfatm_E->SetName("n C4F8O 1.5atm");
TGraph *g_RefractiveIndex_C4F8O_oneandhalfatm_WL=new TGraph(n,Wavelength,RefractiveIndex_C4F8O_oneandhalfatm);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->SetTitle(";wavelength(nm);RefractiveIndex");
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->SetName("n C4F8O 1.5atm");

TGraph *g_RefractiveIndex_C4F8O_oneatm_E=new TGraph(n,PhotonEnergy,RefractiveIndex_C4F8O_oneatm);
g_RefractiveIndex_C4F8O_oneatm_E->SetTitle(";E(eV);RefractiveIndex");
g_RefractiveIndex_C4F8O_oneatm_E->SetName("n C4F8O 1.0atm");
TGraph *g_RefractiveIndex_C4F8O_oneatm_WL=new TGraph(n,Wavelength,RefractiveIndex_C4F8O_oneatm);
g_RefractiveIndex_C4F8O_oneatm_WL->SetTitle(";wavelength(nm);RefractiveIndex");
g_RefractiveIndex_C4F8O_oneatm_WL->SetName("n C4F8O 1.0atm");

TGraph *g_RefractiveIndex_C4F10_oneandhalfatm_E=new TGraph(n,PhotonEnergy,RefractiveIndex_C4F10_oneandhalfatm);
g_RefractiveIndex_C4F10_oneandhalfatm_E->SetTitle(";E(eV);RefractiveIndex");
g_RefractiveIndex_C4F10_oneandhalfatm_E->SetName("n C4F10 1.5atm");
TGraph *g_RefractiveIndex_C4F10_oneandhalfatm_WL=new TGraph(n,Wavelength,RefractiveIndex_C4F10_oneandhalfatm);
g_RefractiveIndex_C4F10_oneandhalfatm_WL->SetTitle(";wavelength(nm);RefractiveIndex");
g_RefractiveIndex_C4F10_oneandhalfatm_WL->SetName("n C4F10 1.5atm");

TGraph *g_RefractiveIndex_C4F10_oneatm_E=new TGraph(n,PhotonEnergy,RefractiveIndex_C4F10_oneatm);
g_RefractiveIndex_C4F10_oneatm_E->SetTitle(";E(eV);RefractiveIndex");
g_RefractiveIndex_C4F10_oneatm_E->SetName("n C4F10 1.0atm");
TGraph *g_RefractiveIndex_C4F10_oneatm_WL=new TGraph(n,Wavelength,RefractiveIndex_C4F10_oneatm);
g_RefractiveIndex_C4F10_oneatm_WL->SetTitle(";wavelength(nm);RefractiveIndex");
g_RefractiveIndex_C4F10_oneatm_WL->SetName("n C4F10 1.0atm");

double RefractiveIndex_LGC_CO2_oneatm[n_LGC] = {  // SL_LGCCgas_SIDIS, used by LGC in both SIDIS and PVDIS
1.000418, 1.000418, 1.000418, 1.000419,
1.000419, 1.000419, 1.000419, 1.000419,
1.000420, 1.000420, 1.000421, 1.000421,
1.000421, 1.000421, 1.000422, 1.000422,
1.000423, 1.000423, 1.000424, 1.000425,
1.000425, 1.000426, 1.000427, 1.000428,
1.000429, 1.000431, 1.000433, 1.000435,
1.000437, 1.000440, 1.000442, 1.000445,
1.000448, 1.000451, 1.000456, 1.000459,
1.000462, 1.000466, 1.000473, 1.000481,
1.000489, 1.000489, 1.000489
};

TGraph *g_RefractiveIndex_LGC_CO2_oneatm_E=new TGraph(n_LGC,PhotonEnergy_LGC,RefractiveIndex_LGC_CO2_oneatm);
g_RefractiveIndex_LGC_CO2_oneatm_E->SetTitle(";E(eV);RefractiveIndex");
g_RefractiveIndex_LGC_CO2_oneatm_E->SetName("n LGC CO2 1.0atm");
TGraph *g_RefractiveIndex_LGC_CO2_oneatm_WL=new TGraph(n_LGC,Wavelength_LGC,RefractiveIndex_LGC_CO2_oneatm);
g_RefractiveIndex_LGC_CO2_oneatm_WL->SetTitle(";wavelength(nm);RefractiveIndex");
g_RefractiveIndex_LGC_CO2_oneatm_WL->SetName("n LGC CO2 1.0atm");

TLegend* leg_Re_E = new TLegend(0.8, 0.8, .95, .95);
leg_Re_E->AddEntry(g_RefractiveIndex_C4F8O_oneandhalfatm_E,g_RefractiveIndex_C4F8O_oneandhalfatm_E->GetName());
leg_Re_E->AddEntry(g_RefractiveIndex_C4F8O_oneatm_E,g_RefractiveIndex_C4F8O_oneatm_E->GetName());
leg_Re_E->AddEntry(g_RefractiveIndex_C4F10_oneandhalfatm_E,g_RefractiveIndex_C4F10_oneandhalfatm_E->GetName());
leg_Re_E->AddEntry(g_RefractiveIndex_C4F10_oneatm_E,g_RefractiveIndex_C4F10_oneatm_E->GetName());
leg_Re_E->AddEntry(g_RefractiveIndex_LGC_CO2_oneatm_E,g_RefractiveIndex_LGC_CO2_oneatm_E->GetName());

TLegend* leg_Re_WL = new TLegend(0.8, 0.8, .95, .95);
leg_Re_WL->AddEntry(g_RefractiveIndex_C4F8O_oneandhalfatm_WL,g_RefractiveIndex_C4F8O_oneandhalfatm_WL->GetName());
leg_Re_WL->AddEntry(g_RefractiveIndex_C4F8O_oneatm_WL,g_RefractiveIndex_C4F8O_oneatm_WL->GetName());
leg_Re_WL->AddEntry(g_RefractiveIndex_C4F10_oneandhalfatm_WL,g_RefractiveIndex_C4F10_oneandhalfatm_WL->GetName());
leg_Re_WL->AddEntry(g_RefractiveIndex_C4F10_oneatm_WL,g_RefractiveIndex_C4F10_oneatm_WL->GetName());
leg_Re_WL->AddEntry(g_RefractiveIndex_LGC_CO2_oneatm_WL,g_RefractiveIndex_LGC_CO2_oneatm_WL->GetName());

TCanvas *c_RefractiveIndex = new TCanvas("RefractiveIndex","RefractiveIndex",1600,800);
c_RefractiveIndex->Divide(2,1);
c_RefractiveIndex->cd(1);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->SetMaximum(1.0025);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->SetMinimum(1.0003);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->GetXaxis()->SetRangeUser(2,7);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->SetMarkerColor(kBlack);
g_RefractiveIndex_C4F8O_oneandhalfatm_E->Draw("AC*");
g_RefractiveIndex_C4F8O_oneatm_E->SetMarkerColor(kBlack);
g_RefractiveIndex_C4F8O_oneatm_E->Draw("C* same");
g_RefractiveIndex_C4F10_oneandhalfatm_E->SetMarkerColor(kRed);
g_RefractiveIndex_C4F10_oneandhalfatm_E->Draw("C* same");
g_RefractiveIndex_C4F10_oneatm_E->SetMarkerColor(kRed);
g_RefractiveIndex_C4F10_oneatm_E->Draw("C* same");
g_RefractiveIndex_LGC_CO2_oneatm_E->SetMarkerColor(kMagenta);
g_RefractiveIndex_LGC_CO2_oneatm_E->Draw("C* same");
leg_Re_E->Draw();
c_RefractiveIndex->cd(2);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->SetMaximum(1.0025);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->SetMinimum(1.0003);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->GetXaxis()->SetRangeUser(200,800);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->SetMarkerColor(kBlack);
g_RefractiveIndex_C4F8O_oneandhalfatm_WL->Draw("AC*");
g_RefractiveIndex_C4F8O_oneatm_WL->SetMarkerColor(kBlack);
g_RefractiveIndex_C4F8O_oneatm_WL->Draw("C* same");
g_RefractiveIndex_C4F10_oneandhalfatm_WL->SetMarkerColor(kRed);
g_RefractiveIndex_C4F10_oneandhalfatm_WL->Draw("C* same");
g_RefractiveIndex_C4F10_oneatm_WL->SetMarkerColor(kRed);
g_RefractiveIndex_C4F10_oneatm_WL->Draw("C* same");
g_RefractiveIndex_LGC_CO2_oneatm_WL->SetMarkerColor(kMagenta);
g_RefractiveIndex_LGC_CO2_oneatm_WL->Draw("C* same");
leg_Re_WL->Draw();

double RefractiveIndex_LGC_UVwin[n_LGC] = {  // SL_H12700UVwin, LGC PMT UV window
1.51568, 1.513, 1.51024, 1.50472,
1.50196, 1.49783, 1.49507, 1.49231,
1.48679, 1.48403, 1.48128, 1.47852,
1.47576, 1.473, 1.46901, 1.46995,
1.46667, 1.46799, 1.46733, 1.46767,
1.46934, 1.47201, 1.47434, 1.47838,
1.4868, 1.49445, 1.5008, 1.50321,
1.504, 1.504, 1.504, 1.504,
1.504, 1.504, 1.504, 1.504,
1.504, 1.504, 1.504, 1.504,
1.504, 1.504, 1.504
};

TGraph *g_RefractiveIndex_LGC_UVwin_E=new TGraph(n_LGC,PhotonEnergy_LGC,RefractiveIndex_LGC_UVwin);
g_RefractiveIndex_LGC_UVwin_E->SetTitle("LGC PMT UV window;E(eV);RefractiveIndex");
TGraph *g_RefractiveIndex_LGC_UVwin_WL=new TGraph(n_LGC,Wavelength_LGC,RefractiveIndex_LGC_UVwin);
g_RefractiveIndex_LGC_UVwin_WL->SetTitle("LGC PMT UV window;wavelength(nm);RefractiveIndex");

TCanvas *c_RefractiveIndex_LGC_UVwin = new TCanvas("RefractiveIndex_LGC_UVwin","RefractiveIndex_LGC_UVwin",1600,800);
c_RefractiveIndex_LGC_UVwin->Divide(2,1);
c_RefractiveIndex_LGC_UVwin->cd(1);
g_RefractiveIndex_LGC_UVwin_E->SetMarkerColor(kMagenta);
g_RefractiveIndex_LGC_UVwin_E->Draw("AC*");
g_RefractiveIndex_LGC_UVwin_E->GetXaxis()->SetRangeUser(2,7);
c_RefractiveIndex_LGC_UVwin->cd(2);
g_RefractiveIndex_LGC_UVwin_WL->SetMarkerColor(kMagenta);
g_RefractiveIndex_LGC_UVwin_WL->Draw("AC*");
g_RefractiveIndex_LGC_UVwin_WL->GetXaxis()->SetRangeUser(200,800);

// double abs_length_C4F8O[n] = {
// 87.7642, 87.7642, 87.7642, 87.7642, 
// 87.7642, 87.7642, 87.7642, 87.7642, 
// 87.7642, 87.7642, 87.7642, 87.7642, 
// 87.7642, 87.7642, 87.616, 90.6032, 
// 93.7711, 102.171, 106.728, 113.497, 
// 119.095, 123.247, 122.872, 119.913, 
// 113.597, 104.84, 98.9037, 94.5673, 
// 91.9657, 90.7634, 90.408, 89.8311, 
// 87.4725, 83.3764, 73.8737, 67.4882, 
// 60.427, 53.2042, 39.2524, 30.5066, 
// 23.3447
// }; //in m

// double CO2_1atm_AbsLen_alt[n] = {
double abs_length_C4F8O[n] = {  
         70316.5, 66796.2, 63314.0, 56785.7,
         53726.5, 49381.2, 46640.7, 44020.0,
         39127.2, 36845.7, 34671.4, 32597.4,
         30621.3, 28743.4, 26154.3, 23775.1,
         22306.7, 19526.3, 18263.4, 16473.0,
         14823.5, 12818.8, 11053.4, 9837.32,
         8351.83, 6747.67, 5648.87, 4694.87,
         3876.99, 3150.27, 2706.97, 2310.46,
         1859.36, 1568.2, 1237.69, 1093.38,
         962.586, 846.065, 643.562, 80.0,
         4.0
};

TGraph *g_abs_length_C4F8O_E=new TGraph(n,PhotonEnergy,abs_length_C4F8O);
// g_abs_length_C4F8O_E->SetTitle(";E(eV);abs_length_C4F8O (m)");
g_abs_length_C4F8O_E->SetTitle(";E(eV);abs_length_CO2 (m)");
TGraph *g_abs_length_C4F8O_WL=new TGraph(n,Wavelength,abs_length_C4F8O);
// g_abs_length_C4F8O_WL->SetTitle(";wavelength(nm);abs_length_C4F8O (m)");
g_abs_length_C4F8O_WL->SetTitle(";wavelength(nm);abs_length_CO2 (m)");
g_abs_length_C4F8O_E->SetName("abs_length");
g_abs_length_C4F8O_WL->SetName("abs_length");

double abs_length_LGC_CO2[n_LGC] = {  // SL_LGCCgas_SIDIS, from CO2_1atm_AbsLen_alt
         70316.5, 66796.2, 63314.0, 56785.7,
         53726.5, 49381.2, 46640.7, 44020.0,
         39127.2, 36845.7, 34671.4, 32597.4,
         30621.3, 28743.4, 26154.3, 23775.1,
         22306.7, 19526.3, 18263.4, 16473.0,
         14823.5, 12818.8, 11053.4, 9837.32,
         8351.83, 6747.67, 5648.87, 4694.87,
         3876.99, 3150.27, 2706.97, 2310.46,
         1859.36, 1568.2, 1237.69, 1093.38,
         962.586, 846.065, 643.562, 80.0,
         4.0, 0.7465970, 0
};

TGraph *g_abs_length_LGC_CO2_E=new TGraph(n_LGC,PhotonEnergy_LGC,abs_length_LGC_CO2);
g_abs_length_LGC_CO2_E->SetName("LGC CO2 1.0atm");
TGraph *g_abs_length_LGC_CO2_WL=new TGraph(n_LGC,Wavelength_LGC,abs_length_LGC_CO2);
g_abs_length_LGC_CO2_WL->SetName("LGC CO2 1.0atm");

TLegend* leg_abs_length_E = new TLegend(0.6, 0.8, .95, .95);
leg_abs_length_E->AddEntry(g_abs_length_C4F8O_E,g_abs_length_C4F8O_E->GetName());
leg_abs_length_E->AddEntry(g_abs_length_LGC_CO2_E,g_abs_length_LGC_CO2_E->GetName());

TLegend* leg_abs_length_WL = new TLegend(0.6, 0.8, .95, .95);
leg_abs_length_WL->AddEntry(g_abs_length_C4F8O_WL,g_abs_length_C4F8O_WL->GetName());
leg_abs_length_WL->AddEntry(g_abs_length_LGC_CO2_WL,g_abs_length_LGC_CO2_WL->GetName());
 


TCanvas *c_abs_length = new TCanvas("abs_length","abs_length",1600,800);
c_abs_length->Divide(2,1);
c_abs_length->cd(1);
g_abs_length_C4F8O_E->Draw("AC*");
// g_abs_length_C4F8O_E->SetMaximum(150);
g_abs_length_C4F8O_E->SetMinimum(1);
g_abs_length_C4F8O_E->GetXaxis()->SetRangeUser(2,7);
g_abs_length_LGC_CO2_E->SetMarkerColor(kMagenta);
g_abs_length_LGC_CO2_E->SetMarkerStyle(24);
g_abs_length_LGC_CO2_E->Draw("P");
leg_abs_length_E->Draw();
c_abs_length->cd(2);
g_abs_length_C4F8O_WL->Draw("AC*");
// g_abs_length_C4F8O_WL->SetMaximum(150);
g_abs_length_C4F8O_WL->SetMinimum(1);
g_abs_length_C4F8O_WL->GetXaxis()->SetRangeUser(200,800);
g_abs_length_LGC_CO2_WL->SetMarkerColor(kMagenta);
g_abs_length_LGC_CO2_WL->SetMarkerStyle(24);
g_abs_length_LGC_CO2_WL->Draw("P");
leg_abs_length_WL->Draw();

// save all canvases in the current dir, pdf works also inside the container
c_QE->SaveAs("plot_parameter_QE.pdf");
c_reflectivity->SaveAs("plot_parameter_reflectivity.pdf");
c_RefractiveIndex->SaveAs("plot_parameter_RefractiveIndex.pdf");
c_RefractiveIndex_LGC_UVwin->SaveAs("plot_parameter_RefractiveIndex_LGC_UVwin.pdf");
c_abs_length->SaveAs("plot_parameter_abs_length.pdf");

}