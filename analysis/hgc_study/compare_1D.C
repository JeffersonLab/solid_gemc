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

void compare_1D()
{
  gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// const int m=4;
// char* input_filename[m]={
// "output_z0_theta8_blockoff_fieldoff_1e4_output.root",
// "output_z0_theta8_blockon_fieldoff_1e4_output.root",
// "output_z0_theta8_blockoff_fieldon_1e4_output.root",
// "output_z0_theta8_blockon_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "blockoff and fieldoff",
// "blockon and fieldoff",
// "blockoff and fieldon",
// "blockon and fieldon",
// };
  
//     const int n=6;
//     const int m=13;
// char* input_filename[n][m]={
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta8_phi0_blockoff_fieldoff_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta8_phi0_blockon_fieldon_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta14.8_phi0_blockon_fieldon_1e3_output.root",
// 
// };
// char* label[n]={
// "z0_theta8_phi0_blockoff_fieldoff",  
// "z0_theta14.8_phi0_blockoff_fieldoff",
//   
// "z0_theta8_phi0_blockoff_fieldon",  
// "z0_theta14.8_phi0_blockoff_fieldon",
//   
// "z0_theta8_phi0_blockon_fieldon",  
// "z0_theta14.8_phi0_blockon_fieldon",
// };


//     const int n=4;
//     const int m=13;
// char* input_filename[n][m]={
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p2.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p2.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p3.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p3.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p4.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p4.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p5.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p5.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p6.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p6.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p7.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p7.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p8.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// 
// 
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p2.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p2.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p3.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p3.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p4.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p4.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p5.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p5.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p6.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p6.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p7.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p7.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_alone_pim/output_p8.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// 
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p2.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p3.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p4.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p5.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p6.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p7.5_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// "data/gemc2.1_geant4.9.6.p03_full_pim/output_p8.0_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",
// 
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p2.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p2.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p3.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p3.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p4.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p4.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p5.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p5.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p6.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p6.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p7.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p7.5_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// // "data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p8.0_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",
// 
// 
// };
// 
// char* label[n]={
// "alone,#pi-,z=0cm,#theta=8deg",    
// "full,#pi-,z=0cm,#theta=8deg",
// "alone,#pi-,z=0cm,#theta=14.8deg",    
// "full,#pi-,z=0cm,#theta=14.8deg",
// };


// char* label[m]={
// };
// int color[m]={4,1};

//   const int m=2;
// char* input_filename[m]={
// "output_z0_theta8_blockoff_fieldoff_1e4_output.root",
// "output_z0_theta8_blockoff_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "theta8_blockoff_fieldoff",
// "theta8_blockoff_fieldon",
// };
// int color[m]={1,2};
  
//     const int m=2;
// char* input_filename[m]={
// "output_z0_theta8_blockoff_fieldon_1e4_output.root",
// "output_z0_theta8_blockon_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "theta8_blockoff_fieldon",
// "theta8_blockon_fieldon",
// };
// int color[m]={2,3};


// const int n=4;
// const int m=13;  
// char* var[m]={"2.0","2.5","3.0","3.5","4.0","4.5","5.0","5.5","6.0","6.5","7.0","7.5","8.0"};    
// string input_filename[n][m];    
// for(int i=0;i<m;i++) {
//   char name[200];
//   sprintf(name,Form("data/gemc2.1_geant4.9.6.p03_alone_pim/output_p%s_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",var[i]));   
//   input_filename[0][i]=name;
//   sprintf(name,Form("data/gemc2.1_geant4.9.6.p03_alone_pim/output_p%s_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",var[i]));   
//   input_filename[1][i]=name;
//   sprintf(name,Form("data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p%s_z0_theta8_phi0_blockoff_fieldon_1e3_output.root",var[i]));   
//   input_filename[2][i]=name;
//   sprintf(name,Form("data/gemc2.1_geant4.9.6.p03_alone_pim_old/output_p%s_z0_theta14.8_phi0_blockoff_fieldon_1e3_output.root",var[i]));   
//   input_filename[3][i]=name;
// }
// 
// char* label[n]={
// "alone,#pi-,z=0cm,#theta=8deg",
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// "alone,#pi-,z=0cm,#theta=8deg",    
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// };

// const int m=11;  
// char* var[m]={"2.5","3.0","3.5","4.0","4.5","5.0","5.5","6.0","6.5","7.0","7.5"};    

// const int n=6;
// string input_filename[n][m];

// for(int i=0;i<m;i++) {
//   char name[200]; 
// 
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[0][i]=name;  
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[1][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[2][i]=name;  
//   
// 
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmoveNO_coneYES/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[3][i]=name;     
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmoveNO_coneYES/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[4][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmoveNO_coneYES/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[5][i]=name;    
//   
// }
// char* label[n]={
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",  
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// "(old) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",
// "(old) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(old) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// };
// int MarkerStyle[n]={4,4,4,8,8,8};

// const int n=6;
// string input_filename[n][m];
// for(int i=0;i<m;i++) {
//   char name[200]; 
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[0][i]=name;  
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[1][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[2][i]=name;  
//   
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03_WLS_meas/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[3][i]=name;  
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03_WLS_meas/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[4][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.2/pmtmove10cm_cone10cmtilt65deg/QE_H12700_03_WLS_meas/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[5][i]=name;  
// }
// char* label[n]={
// "(uncoated) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",  
// "(uncoated) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(uncoated) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// "(coated) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",
// "(coated) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(coated) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// };
// int MarkerStyle[n]={4,4,4,26,26,26};
// int color[n]={1,2,3,1,2,3};

// const int n=6;
// const int m=13;  
// char* var[m]={"2.0","2.5","3.0","3.5","4.0","4.5","5.0","5.5","6.0","6.5","7.0","7.5","8.0"};    
// string input_filename[n][m];
// for(int i=0;i<m;i++) {
//   char name[200]; 
// // pmtmove10cm_cone10cmtilt65deg
//   
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_shift20cmdown_new/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[0][i]=name;  
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_shift20cmdown_new/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[1][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_shift20cmdown_new/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[2][i]=name; 
//   
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_what/output_pim_z0_p%s_theta7.5_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[3][i]=name;  
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_what/output_pim_z0_p%s_theta8.0_phi0_blockoff_fieldoff_1e3_output.root",var[i]));  
//   input_filename[4][i]=name;
//   sprintf(name,Form("data/JLAB_VERSION_1.3/He3_what/output_pim_z0_p%s_theta14.8_phi0_blockoff_fieldoff_1e3_output.root",var[i]));   
//   input_filename[5][i]=name;    
// }
// char* label[n]={
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",  
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(new) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// "(NH3_shift20cmdown) alone,#pi-,field off,block off,z=0cm,#theta=7.5deg,#phi=0deg",
// "(NH3_shift20cmdown) alone,#pi-,field off,block off,z=0cm,#theta=8.0deg,#phi=0deg",
// "(NH3_shift20cmdown) alone,#pi-,field off,block off,z=0cm,#theta=14.8deg,#phi=0deg",  
// };
// int MarkerStyle[n]={4,4,4,26,26,26};
// int color[n]={1,2,3,1,2,3};

// const int n=4;
// const int m=2;
// char* var[m]={"2.5","7.5"};
// string input_filename[n][m];
// for(int i=0;i<m;i++) {
//   char name[1000];
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p%s_theta8.0_phi45-145_noblock_1e3_output.root",var[i]));  
//   input_filename[0][i]=name;
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p%s_theta14.8_phi45-145_noblock_1e3_output.root",var[i]));  
//   input_filename[1][i]=name;
//   
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p%s_theta8.0_phi45-145_yesblock_1e3_output.root",var[i]));  
//   input_filename[2][i]=name;
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p%s_theta14.8_phi45-145_yesblock_1e3_output.root",var[i]));  
//   input_filename[3][i]=name;  
// }
// char* label[n]={
// "(no) alone,#pi-,field on,z=0cm,#theta=8.0deg,#phi=0deg",  
// "(nu) alone,#pi-,field on,z=0cm,#theta=14.8deg,#phi=0deg",
// "(yes) alone,#pi-,field on,z=0cm,#theta=8.0deg,#phi=0deg",  
// "(yes) alone,#pi-,field on,z=0cm,#theta=14.8deg,#phi=0deg",
// };
// int MarkerStyle[n]={4,4,26,26};
// int color[n]={1,2,1,2};

const int n=3;
const int m=2;
string input_filename[n][m];

  char name[1000];
  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi45-145_noblock_1e3_output.root"));  
  input_filename[0][0]=name;
  
  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta14.8_phi45-145_noblock_1e3_output.root"));  
  input_filename[0][1]=name;
  
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi55-65_yesblock_1e3_output.root"));  
//   input_filename[1][0]=name;
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta14.8_phi65-75_yesblock_1e3_output.root"));  
//   input_filename[1][1]=name;  
// 
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi55-65_yesblockmirror1_1e3_output.root"));  
//   input_filename[2][0]=name;
//   sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta14.8_phi65-75_yesblockmirror1_1e3_output.root"));  
//   input_filename[2][1]=name;  

  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi57-63_yesblock_1e4_output.root"));  
  input_filename[1][0]=name;
  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta14.8_phi67-73_yesblock_1e4_output.root"));  
  input_filename[1][1]=name;  

  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi57-63_yesblockmirror1_1e4_output.root"));  
  input_filename[2][0]=name;
  sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testblock/hgc_SIDIS_He3_pim_z350_p2.5_theta14.8_phi67-73_yesblockmirror1_1e4_output.root"));  
  input_filename[2][1]=name;  

char* label[n]={
"(no) alone,#pi-,field on,z=0cm,#theta=8.0deg,#phi=0deg",  
"(yes) alone,#pi-,field on,z=0cm,#theta=14.8deg,#phi=0deg",
"(mirror1) alone,#pi-,field on,z=0cm,#theta=14.8deg,#phi=0deg",
};
int MarkerStyle[n]={4,26,27};
int color[n]={1,2,3};

TH1F *hcount_total_p[n];
for(int j=0;j<n;j++){
hcount_total_p[j]=new TH1F(Form("hcount_total_p_%i",j),"photoelectron;P (GeV);count",m,1.75,8.25);
// hcount_total_p[j]=new TH1F(Form("hcount_total_p_%i",j),";P (GeV);photoelectron count (sim*0.5)",m,2.25,7.75);  
// hcount_total_p->SetAxisRange(2,8,"X");
hcount_total_p[j]->SetAxisRange(0,150,"Y");   
hcount_total_p[j]->SetMarkerStyle(MarkerStyle[j]);
hcount_total_p[j]->SetMarkerSize(2);
hcount_total_p[j]->SetMarkerColor(color[j]);
hcount_total_p[j]->SetLineColor(color[j]);
}


TCanvas *c = new TCanvas("compare","compare",1600,800);
c->Divide(n,1);
TFile *input[n][m];
TH1F *h[n][m];  
TLegend* leg[n];
for(int j=0;j<n;j++){
  c->cd(j+1);  
  leg[j]= new TLegend(0.7, 0.95-0.05*m, 0.95, 0.95);  
for(int i=0;i<m;i++){
//   cout << j << " " << i << endl;  
//     if(i==3 || i>6) continue;
  input[j][i]=new TFile(input_filename[j][i].c_str());
//   cout << " " << input_filename[j][i] << endl;
  if (input[j][i]->IsZombie()) {
    cout << "Error opening ratefile " << input_filename[j][i] << endl;
    exit(-1);
  }
  else cout << "ok open file " << input_filename[j][i] << endl;
  
  
//   char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
//   cout << hstname << endl;
//   h[j][i]=(TH1F*) input[j][i]->Get("hcount");
  h[j][i]=(TH1F*) input[j][i]->Get("hnpe_no0");  
//   h[i]->SetAxisRange(ymin,ymax,"Y");    
//   h[i]->SetAxisRange(xmin,xmax,"X");  
//   h[i]->SetTitle(Form("%s %s",title,label[i]));
//   h[j][i]->Rebin(2);    
  h[j][i]->SetLineColor(i+1);  
  if (i==0) h[j][i]->Draw();  
  else h[j][i]->Draw("same");
//   h[i]->SetMarkerStyle(8);
//   g[i]->SetMarkerSize(0.15*(m-i));    
//   g[i]->SetMarkerColor(color[i]);
//   g[i]->SetLineColor(color[i]);
//   g[i]->Draw("same P");
  hcount_total_p[j]->SetBinContent(i+1,h[j][i]->GetMean());
  hcount_total_p[j]->SetBinError(i+1,h[j][i]->GetRMS());  
  
//   cout << h[i]->Integral() << endl;
//     input.Close();
//   leg->AddEntry(h[i], Form("%s   %02f",input_filename[i],h[i]->GetMean()),"l");  
  leg[j]->AddEntry(h[j][i], Form("%0.1f %0.1f",h[j][i]->GetMean(),h[j][i]->GetRMS()),"l");    
//   leg->AddEntry(g[i], label[i],"l");    
  
}
leg[j]->Draw();
}


TLegend* legend;
  legend= new TLegend(0.5, 0.98-0.05*n, 0.95, 0.98);  
TCanvas *c1 = new TCanvas("c","c",1000,800);
for(int j=0;j<n;j++){
if (j==0) hcount_total_p[j]->Draw("E1");
else hcount_total_p[j]->Draw("E same");
  legend->AddEntry(hcount_total_p[j], Form("%s",label[j]),"lep");    
}
legend->Draw();


}