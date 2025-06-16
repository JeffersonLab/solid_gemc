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

void compare_testmirror()
{
  gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

const int n=3;
const int m=3;
string input_filename[n][m];

char name[1000];
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row+5cm_1e3_output.root"));  
input_filename[0][0]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row+10cm_1e3_output.root"));  
input_filename[1][0]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row+15cm_1e3_output.root"));  
input_filename[2][0]=name;

sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row-5cm_1e3_output.root"));  
input_filename[0][1]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row-10cm_1e3_output.root"));  
input_filename[1][1]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_row-15cm_1e3_output.root"));  
input_filename[2][1]=name;

sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_nochange_1e3_output.root"));  
input_filename[0][2]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_nochange_1e3_output.root"));  
input_filename[1][2]=name;
sprintf(name,Form("/work/halla/solid/zwzhao/solid/hgc_sim/JLAB_VERSION_1.3/testmirror/hgc_SIDIS_He3_pim_z350_p2.5_theta8.0_phi360_nochange_1e3_output.root"));  
input_filename[2][2]=name;

char* label[n][m]={
"roll+5cm",
"roll-5cm",
"nochange",
"roll+10cm",
"roll-10cm",
"nochange",
"roll+15cm",
"roll-15cm",
"nochange",
};
int color[m]={1,2,3};
int MarkerStyle[n]={4,26,27};
int linestyle[n]={1,2,9};

TH1F *hcount_total_p[n];
for(int j=0;j<n;j++){
hcount_total_p[j]=new TH1F(Form("hcount_total_p_%i",j),"photoelectron;P (GeV);count",m,1.75,8.25);
// hcount_total_p[j]=new TH1F(Form("hcount_total_p_%i",j),";P (GeV);photoelectron count (sim*0.5)",m,2.25,7.75);  
// hcount_total_p->SetAxisRange(2,8,"X");
hcount_total_p[j]->SetAxisRange(0,100,"Y");   
hcount_total_p[j]->SetMarkerStyle(MarkerStyle[j]);
hcount_total_p[j]->SetMarkerSize(2);
hcount_total_p[j]->SetMarkerColor(color[j]);
hcount_total_p[j]->SetLineColor(color[j]);
}


TCanvas *c = new TCanvas("compare","compare",1600,800);
// c->Divide(n,1);
TFile *input[n][m];
TH1F *h[n][m];  
TLegend* leg= new TLegend(0.5, 0.95-0.05*m*n, 0.95, 0.95);  

for(int i=0;i<m;i++){
  for(int j=0;j<n;j++){
  //   c->cd(j+1);  
//   cout << j << " " << i << endl;  
    if(i==2 && j!=0) continue;
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
  h[j][i]->SetTitle("Testing HGC mirror;Npe;normalized event count");
  h[j][i]->Rebin(2);    
  h[j][i]->Scale(1./h[j][i]->Integral());   
  h[j][i]->SetAxisRange(0,100,"X");  
  h[j][i]->SetAxisRange(0,0.2,"Y");     
  h[j][i]->SetLineStyle(linestyle[j]);  
  h[j][i]->SetLineColor(color[i]);  
  if (i==0 && j==0) h[j][i]->Draw("HIST");  
  else h[j][i]->Draw("HIST same");
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
  leg->AddEntry(h[j][i], Form("%s, Mean %0.1f, RMS %0.1f",label[j][i],h[j][i]->GetMean(),h[j][i]->GetRMS()),"l");    
//   leg->AddEntry(g[i], label[i],"l");    
  
}
leg->Draw();
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