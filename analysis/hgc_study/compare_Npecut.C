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
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void compare_Npecut()
{
gROOT->Reset();
// gStyle->SetPalette(57);
gStyle->SetOptStat(0);
gStyle->SetPadColor(0);
  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);  
  
  gStyle->SetLabelSize(0.05,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetTitleOffset(0.7,"y");
  gStyle->SetTitleOffset(1,"x");    
  gStyle->SetTitleSize(0.07,"t");   

// char* label[4]={
// "alone,#pi-,z=0cm,#theta=8deg",
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// "alone,#pi-,z=0cm,#theta=8deg",    
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// };

// const int n=4;

// char* input_filename[n]={
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_pim_z350_p2.5_theta8.0_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_km_z350_p2.5_theta8.0_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_pim_z350_p7.5_theta14.5_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_km_z350_p7.5_theta14.5_1e5_output.root",
// };
// char* label[2]={
// "#pi^{-}",
// "K^{-}",    
// };
// char* dir1="solid_SIDIS_He3_full_lgc/";
//   char* dir1="solid_SIDIS_He3_full_lgcec/";
char* dir1="solid_PVDIS_LD2_full_lgcec/";

// char* dir2="backgroundno_nodecaytrack/";
// char* dir2="background3sector_nodecaytrack/";
char* dir2="background3sectordouble_nodecaytrack/";

const int n=2;
char* input_filename[n]={
// "solid_SIDIS_He3_moved_full_ele_z350_p2.5_theta8.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p2.5_theta8.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_ele_z350_p4.5_theta14.5_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p4.5_theta14.5_1e6_output.root",

// "solid_SIDIS_He3_moved_full_ele_z350_p3.0_theta12.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p3.0_theta12.0_1e6_output.root",

"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output.root",
"solid_PVDIS_LD2_moved_full_pim_z10_p2.5_theta22_1e6_output.root",    
};
char* label[2]={
"e^{-}",
"#pi^{-}",    
};

// char* dir1="solid_PVDIS_LD2_full_lgc/";
// 
// char* dir2="backgroundno_nodecaytrack/";
// // char* dir2="background3sector_nodecaytrack/";
// // char* dir2="background3sectordouble_nodecaytrack/";
// 
// const int n=4;
// char* input_filename[n]={
// "solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output.root",
// "solid_PVDIS_LD2_moved_full_pim_z10_p2.5_theta22_1e6_output.root",
// "solid_PVDIS_LD2_moved_full_ele_z10_p4.5_theta34_1e6_output.root",
// "solid_PVDIS_LD2_moved_full_pim_z10_p4.5_theta34_1e6_output.root",
// };
// char* label[2]={
// "e^{-}",
// "#pi^{-}",    
// };


TCanvas *c = new TCanvas("compare","compare",n/2*800,800);
c->Divide(n/2,2);
TFile *input[n];
TH1F *h[n];
TH1F *h_cut[n];
TLegend* leg[n];
for(int j=0;j<n/2;j++){
  c->cd(j+1); 
  gPad->SetLogy(1);
  leg[j]= new TLegend(0.75, 0.75, 0.85, 0.9);  
for(int i=0;i<2;i++){
//     if(i==3 || i>6) continue;
  input[2*j+i]=new TFile(Form("%s%s%s",dir1,dir2,input_filename[2*j+i]));
  if (input[2*j+i]->IsZombie()) {
    cout << "Error opening file " << input_filename[2*j+i] << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename[2*j+i] << endl;
  
//   char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
//   cout << hstname << endl;
  h[2*j+i]=(TH1F*) input[2*j+i]->Get("hnpe");

  h[2*j+i]->SetAxisRange(1,1e6,"Y");    
//   h[2*j+i]->SetAxisRange(0,100,"X");  
  h[2*j+i]->SetTitle(";Npe;counts");  
  
  if (i==0){
    h[2*j+i]->SetLineColor(kRed);  
    h[2*j+i]->Draw();
  }
  else{
    h[2*j+i]->SetLineColor(kGreen);      
    h[2*j+i]->Scale(0.1); //manually make the bad particles 10 times smaller
    h[2*j+i]->Draw("same");
  }
  
//   h[i]->SetMarkerStyle(8);
//   g[i]->SetMarkerSize(0.15*(m-i));    
//   g[i]->SetMarkerColor(color[i]);
//   g[i]->SetLineColor(color[i]);
//   g[i]->Draw("same P");
//   hcount_total_p[j]->SetBinContent(i+1,h[2*j+i]->GetMean());
//   hcount_total_p[j]->SetBinError(i+1,h[2*j+i]->GetRMS());  
  
//   cout << h[i]->Integral() << endl;
//     input.Close();
//   leg->AddEntry(h[i], Form("%s   %02f",input_filename[i],h[i]->GetMean()),"l");  
//   leg[j]->AddEntry(h[2*j+i], Form("%02f %02f",h[2*j+i]->GetMean(),h[2*j+i]->GetRMS()),"l");    
  leg[j]->AddEntry(h[2*j+i], label[i],"l");    

  h_cut[2*j+i]=(TH1F*) h[2*j+i]->Clone();
  h_cut[2*j+i]->Reset("ICESM");  
  int Nbin=h[2*j+i]->GetNbinsX();
//   int Ntotal=h[2*j+i]->GetEntries();
  int Ntotal=h[2*j+i]->GetEntries();  
  cout << Nbin << " " << Ntotal << endl;
  for(int bin=1;bin<Nbin+1;bin++){
    if (i==0) h_cut[2*j+i]->SetBinContent(bin,h[2*j+i]->Integral(bin,Nbin)/Ntotal);
    else h_cut[2*j+i]->SetBinContent(bin,h[2*j+i]->Integral(1,bin)/Ntotal*10);  //10 because manually make the bad particles 10 times smaller    
//     cout << h[2*j+i]->Integral(1,bin) << endl;
  }
}
leg[j]->Draw();

c->cd(j+1+n/2); 
for(int i=0;i<2;i++){
//   h_cut[2*j+i]->SetAxisRange(0,1,"Y");
//   h_cut[2*j+i]->SetAxisRange(0.8,1,"Y");        
  h_cut[2*j+i]->SetAxisRange(0,1,"Y");          
  h_cut[2*j+i]->SetAxisRange(0,100,"X");  
  h_cut[2*j+i]->SetTitle(";Npe cut;FOM");
  
  if (i==0) {
    h_cut[2*j+i]->SetLineColor(kRed);  
    h_cut[2*j+i]->Draw("");
  }
  else{
    h_cut[2*j+i]->SetLineColor(kGreen);      
    h_cut[2*j+i]->Draw("same");
  }
}
leg[j]->Draw();

}

// TLegend* legend;
//   legend= new TLegend(0.5, 0.9-0.05*n, 0.95, 0.95);  
// TCanvas *c1 = new TCanvas("c","c",1000,800);
// for(int j=0;j<n;j++){
// if (j==0) hcount_total_p[j]->Draw("E1");
// else hcount_total_p[j]->Draw("E1 same");
//   legend->AddEntry(hcount_total_p[j], Form("%s",label[j]),"l");    
// }
// legend->Draw();


c->SaveAs("compare.png");

}