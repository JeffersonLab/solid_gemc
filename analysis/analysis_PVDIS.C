  // TF1* fq22 = new TF1("fq22",Form("%f",q2min),0,1);
  // fq22->SetLineColor(2);
  // fq22->Draw("same");
  // TLatex* tq22 = new TLatex (0.2,7,Form ("Q^{2}>%3.1f GeV^{2}", q2min));
  // tq22->SetTextColor(2);
  // tq22->SetTextSize(0.035);
  // tq22->Draw();
  // TF1* fw22 = new TF1("fw22",Form ("x*(%f-.938*.938)/(1-x)", wmin*wmin),0,1);
  // fw22->SetLineColor(6);
  // fw22->Draw("same");
  // TLatex* tw22 = new TLatex (0.25,2.,Form ("W>%3.1f GeV", wmin));
  // tw22->SetTextColor(6);
  // tw22->SetTextSize(0.035);
  // tw22->Draw();
  // TLine* fxbj2 = new TLine(xbjmin,0,xbjmin,14);
  // fxbj2->SetLineWidth(2);
  // fxbj2->SetLineColor(1);
  // fxbj2->Draw("same");
  // TLatex* txbj2 = new TLatex (.57,2.,Form("x_{bj}>%4.2f", xbjmin));
  // txbj2->SetTextColor(1);
  // txbj2->SetTextSize(0.035);
  // txbj2->Draw();
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 

#include <TArc.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>

#include "analysis_tree_solid_ec.C"
#include "analysis_tree_solid_lgc.C"
//#include "PVDIS_LD2_FAEC_electron_trigger_HallDRakitha.C"
#include "Get_PVDIS_trigger_efficiency_EC-1.C"

using namespace std;


void ContoursPT (Double_t ebeam, Bool_t writec = false)
{
  vector <Double_t> Q2v;
  vector <Double_t> Wv;
  vector <Double_t> xv;
  Q2v.push_back (3.0);
  Q2v.push_back (4.5);
  Q2v.push_back (6.0);
  Q2v.push_back (7.5);
  Q2v.push_back (9.0);
  Wv.push_back (1.0);
  Wv.push_back (1.5);
  Wv.push_back (2.0);
  Wv.push_back (3.0);
  xv.push_back (0.35);
  xv.push_back (0.55);
  xv.push_back (0.75);
  xv.push_back (0.95);

  TString Q2formula (TString ("%f") + Form("/(2*%f*(1-cos(x*3.14159/180)))", ebeam));
  TString Wformula (TString (Form("(.938*.938+2*.938*%f-", ebeam)) + "%f" + Form (")/(2*(.938+%f*(1-cos(x*3.14159/180))))", ebeam));
  TString xformula (TString("(.938*%f") + Form ("*%f)/(%f*(1-cos(x*3.14159/180))+.938*", ebeam, ebeam) + "%f)");

  TLegend* tl = new TLegend (.65, .6, .85, .9);

  for (UInt_t i = 0; i < Q2v.size(); ++i)
    {
      TF1* fq2 = new TF1("fq2", Form(Q2formula, Q2v[i]), 10, 50);
      fq2->SetLineColor (kBlue);
      fq2->SetLineStyle (i+1);
      fq2->SetLineWidth (3);
      fq2->Draw("same");

      if (writec)
	fq2->Write (Form ("fq2_%d", i));
      tl->AddEntry (fq2, Form ("Q^{2}=%f", Q2v[i]), "l");
    }
  for (UInt_t i = 0; i < Wv.size(); ++i)
    {
      TF1* fw = new TF1("fw", Form(Wformula, pow (Wv[i],2)), 10, 50);
      fw->SetLineColor (kMagenta);
      fw->SetLineStyle (i+1);
      fw->SetLineWidth (3);
      fw->Draw("same");      
      if (writec)
	fw->Write (Form ("fw_%d", i));
      tl->AddEntry (fw, Form ("W=%f", Wv[i]), "l");
    }
  for (UInt_t i = 0; i < xv.size(); ++i)
    {
      TF1* fx = new TF1("fx", Form(xformula, xv[i], xv[i]), 10, 50);
      fx->SetLineColor (kBlack);
      fx->SetLineStyle (i+1);
      fx->SetLineWidth (3);
      fx->Draw("same");
      if (writec)
	fx->Write (Form ("fx_%d", i));
      tl->AddEntry (fx, Form ("x=%f", xv[i]), "l");
    }
  tl->SetTextSize (0.015);
  tl->Draw();
}


void ContoursQ2x (Double_t ebeam, Bool_t writec = false)
{
  vector <Double_t> Q2v;
  vector <Double_t> Wv;
  vector <Double_t> xv;
  Q2v.push_back (3.0);
  Q2v.push_back (4.5);
  Q2v.push_back (6.0);
  Q2v.push_back (7.5);
  Q2v.push_back (9.0);
  Wv.push_back (1.0);
  Wv.push_back (1.5);
  Wv.push_back (2.0);
  Wv.push_back (3.0);
  xv.push_back (0.35);
  xv.push_back (0.55);
  xv.push_back (0.75);
  xv.push_back (0.95);

  TString Wformula ("x*(%f-.938*.938)/(1-x)");

  TLegend* tl = new TLegend (.15, .6, .35, .9);

  for (UInt_t i = 0; i < Q2v.size(); ++i)
    {
      TLine* lq2 = new TLine (0.0, Q2v[i], 1.0, Q2v[i]);
      lq2->SetLineColor (kBlue);
      lq2->SetLineStyle (i+1);
      lq2->SetLineWidth (3);
      lq2->Draw("same");

      tl->AddEntry (lq2, Form ("Q^{2}=%f", Q2v[i]), "l");
    }
  for (UInt_t i = 0; i < Wv.size(); ++i)
    {
      TF1* fw = new TF1("fw", Form(Wformula, pow (Wv[i],2)), 0.0, 1.0);
      fw->SetLineColor (kMagenta);
      fw->SetLineStyle (i+1);
      fw->SetLineWidth (3);
      fw->Draw("same");      
      if (writec)
	fw->Write (Form ("fw2_%d", i));
      tl->AddEntry (fw, Form ("W=%f", Wv[i]), "l");
    }
  for (UInt_t i = 0; i < xv.size(); ++i)
    {
      TLine* lx = new TLine (xv[i], 0., xv[i], 14.);
      lx->SetLineColor (kBlack);
      lx->SetLineStyle (i+1);
      lx->SetLineWidth (3);
      lx->Draw("same");
      tl->AddEntry (lx, Form ("x=%f", xv[i]), "l");
    }
  tl->SetTextSize (0.015);
  tl->Draw();
}

void analysis_PVDIS (string input_filename,
		     UInt_t nev = 1e9, 
		     UInt_t sev = 0, 
		     bool kinecut=false, 
		     bool lgcpecut=true, 
		     bool sixpointsix=false)
{
  // input_filename is normally DIS (no kill) file
  // debug = true will analyze 1% of file
  // kinecut = true will impose kinematics cuts, useful for LGC efficiency, not for acceptance
  // lgcpecut = true will impose PE >= 3 requirement for acceptance, otherwise imposes
  //    lgc fluxg counter requirement

  Double_t ebeam = sixpointsix ? 6.6 : 11.0;
  Double_t q2min = sixpointsix ? 3.0 : 6.0;
  Double_t wmin = sixpointsix ? 1.5 : 2.0;
  Double_t xbjmin = sixpointsix ? 0.50 : 0.55;

  gROOT->Reset();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1111111);

  gStyle->SetPadColor(0);
  gStyle->SetPalette(1);

  gStyle->SetLabelSize(0.05,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.08,"xyz");  
  gStyle->SetTitleOffset(0.7,"y");
  gStyle->SetTitleOffset(1,"x");    
  gStyle->SetTitleSize(0.07,"t");  
  
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);  
  
  TRandom3 rand;
  rand.SetSeed(0);
  
  const double DEG=180./3.1415926;
  Int_t nphecut = 3;
  
  char the_filename[200];
  sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());
  
  char output_filename[200];
  sprintf(output_filename, "%s_output.root",the_filename);
  TFile *outputfile=new TFile(output_filename, "recreate");
  
  const int m=12;
  const char *detector_name[m]={"GEM 1","GEM 2","GEM 3","GEM 4","GEM 5","GEM 6","LGC","HGC","FASPD","LASPD","FAEC","LAEC"};
  TH2F *hflux_hitxy[m];
  for (Int_t i=0;i<m;i++) {
    hflux_hitxy[i]=new TH2F(Form("hflux_hitxy_%i",i),Form("flux_hitxy %s;x(cm);y(cm)",detector_name[i]),600,-300,300,600,-300,300);
  }
  
  TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","gen_ThetaP",40,10,50,110,0,11);     
  hgen_ThetaP->SetStats(0);
  TH2F *hacceptance_ThetaP;
  hacceptance_ThetaP=new TH2F("acceptance_ThetaP_FA","PVDIS acceptance",40,10,50,110,0,11);     
  hacceptance_ThetaP->SetStats(0);
  TH2F *hacceptanceeff_ThetaP;
  hacceptanceeff_ThetaP=new TH2F("acceptanceeff_ThetaP_FA","PVDIS electron acceptance & efficiency",40,10,50,110,0,11);     
  hacceptanceeff_ThetaP->SetStats(0);
  TH2F *hacceptance_RhitPhit;
  hacceptance_RhitPhit=new TH2F("acceptance_RhitPhit_FA","PVDIS acceptance",60,0.,300.,110,0,11);     
  hacceptance_RhitPhit->SetStats(0);
  TH2F *hacceptanceeff_RhitPhit;
  hacceptanceeff_RhitPhit=new TH2F("acceptanceeff_RhitPhit_FA","PVDIS electron acceptance & efficiency",60,0.,300.,110,0,11);     
  hacceptanceeff_RhitPhit->SetStats(0);
  TH1F *hacceptance_P;
  hacceptance_P=new TH1F("acceptance_P_FA","PVDIS acceptance",110,0,11);     
  hacceptance_P->SetStats(0);
  TH1F *hacceptanceeff_P;
  hacceptanceeff_P=new TH1F("acceptanceeff_P_FA","PVDIS electron acceptance & efficiency",110,0,11);     
  hacceptanceeff_P->SetStats(0);
  
  TH1F *hdist_x=new TH1F("dist_x","dist_x", 50, 0. , 1.);
  TH1F *hdist_Q2=new TH1F("dist_Q2","dist_Q2", 70, 0. , 14.);
  
  TH2F *hgen_Q2x=new TH2F("gen_Q2x","gen_Q2x", 50, 0. , 1., 70, 0., 14.);
  hgen_Q2x->SetStats(0);
  TH2F *hacceptance_Q2x;
  hacceptance_Q2x=new TH2F("acceptance_Q2x_FA","PVDIS acceptance", 50, 0. , 1., 70, 0., 14.);   
  hacceptance_Q2x->SetStats(0);
  TH2F *hacceptanceeff_Q2x;
  hacceptanceeff_Q2x=new TH2F("acceptanceeff_Q2x_FA","PVDIS electron acceptance & efficiency", 50, 0. , 1., 70, 0., 14.);   
  hacceptanceeff_Q2x->SetStats(0);
  
  
  TH1F *hnphe_lgc=new TH1F("hnphe_lgc","Photoelectrons",60,-0.5,59.5);
  hnphe_lgc->SetStats(0);
  TH1F *hsectoring_ec_lgc=new TH1F("hsectoring_ec_lgc","PE rate per segment",30,-14.5,15.5);
  hsectoring_ec_lgc->SetStats(0);
  
  // TH1F *htotEdep_ec=new TH1F("htotEdep_ec","htotEdep_ec",100,0,2000);
  
  // TH2F *htotEdep_ec_gen=new TH2F("htotEdep_ec_gen","htotEdep_ec_gen",100,0,2000,110,0,11000);
  
  TFile *file=new TFile(input_filename.c_str());
  if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);  }
  else cout << "open file " << input_filename << endl;    
  
  outputfile->cd();

  TTree *tree_header = (TTree*) file->Get("header");
  vector <string> *header_time=0;
  vector <int> *header_evn=0,*header_evn_type=0;
  vector <double> *header_beamPol=0;
  vector <double> *header_var1=0,*header_var2=0,*header_var3=0,*header_var4=0,*header_var5=0,*header_var6=0,*header_var7=0,*header_var8=0;
  tree_header->SetBranchAddress("time",&header_time);
  tree_header->SetBranchAddress("evn",&header_evn);
  tree_header->SetBranchAddress("evn_type",&header_evn_type);
  tree_header->SetBranchAddress("beamPol",&header_beamPol);
  tree_header->SetBranchAddress("var1",&header_var1);
  tree_header->SetBranchAddress("var2",&header_var2);
  tree_header->SetBranchAddress("var3",&header_var3);
  tree_header->SetBranchAddress("var4",&header_var4);
  tree_header->SetBranchAddress("var5",&header_var5);
  tree_header->SetBranchAddress("var6",&header_var6);
  tree_header->SetBranchAddress("var7",&header_var7);
  tree_header->SetBranchAddress("var8",&header_var8);
  // if(debug){
  // char *branchname_header[12]={"time","evn","evn_type","beamPol","var1","var2","var3","var4","var5","var6","var7","var8"};
  // cout << endl << "tree_header" << endl;
  // for (Int_t i=0;i<12;i++) { 
  // cout << branchname_header[i] << " " <<  tree_header->GetBranch(branchname_header[i])->GetLeaf(branchname_header[i])->GetTypeName() << ",";
  // }
  // }
  
  TTree *tree_generated = (TTree*) file->Get("generated");
  vector <int> *gen_pid=0;
  vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
  tree_generated->SetBranchAddress("pid",&gen_pid);
  tree_generated->SetBranchAddress("px",&gen_px);
  tree_generated->SetBranchAddress("py",&gen_py);
  tree_generated->SetBranchAddress("pz",&gen_pz);
  tree_generated->SetBranchAddress("vx",&gen_vx);
  tree_generated->SetBranchAddress("vy",&gen_vy);
  tree_generated->SetBranchAddress("vz",&gen_vz);
  // if(debug){
  // char *branchname_generated[7]={"pid","px","py","pz","vx","vy","vz"};
  // cout << endl << "tree_generated" << endl;
  // for (Int_t i=0;i<7;i++) { 
  // cout << branchname_generated[i] << " " <<  tree_generated->GetBranch(branchname_generated[i])->GetLeaf(branchname_generated[i])->GetTypeName() << ",";
  // }
  // }
  
  TTree *tree_flux = (TTree*) file->Get("flux");
  vector<int> *flux_id=0,*flux_hitn=0;
  vector<int> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
  vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
  tree_flux->SetBranchAddress("hitn",&flux_hitn);
  tree_flux->SetBranchAddress("id",&flux_id);
  tree_flux->SetBranchAddress("pid",&flux_pid);
  tree_flux->SetBranchAddress("mpid",&flux_mpid);
  tree_flux->SetBranchAddress("tid",&flux_tid);
  tree_flux->SetBranchAddress("mtid",&flux_mtid);
  tree_flux->SetBranchAddress("otid",&flux_otid);
  tree_flux->SetBranchAddress("trackE",&flux_trackE);
  tree_flux->SetBranchAddress("totEdep",&flux_totEdep);
  tree_flux->SetBranchAddress("avg_x",&flux_avg_x);
  tree_flux->SetBranchAddress("avg_y",&flux_avg_y);
  tree_flux->SetBranchAddress("avg_z",&flux_avg_z);
  tree_flux->SetBranchAddress("avg_lx",&flux_avg_lx);
  tree_flux->SetBranchAddress("avg_ly",&flux_avg_ly);
  tree_flux->SetBranchAddress("avg_lz",&flux_avg_lz);
  tree_flux->SetBranchAddress("px",&flux_px);
  tree_flux->SetBranchAddress("py",&flux_py);
  tree_flux->SetBranchAddress("pz",&flux_pz);
  tree_flux->SetBranchAddress("vx",&flux_vx);
  tree_flux->SetBranchAddress("vy",&flux_vy);
  tree_flux->SetBranchAddress("vz",&flux_vz);
  tree_flux->SetBranchAddress("mvx",&flux_mvx);
  tree_flux->SetBranchAddress("mvy",&flux_mvy);
  tree_flux->SetBranchAddress("mvz",&flux_mvz);
  tree_flux->SetBranchAddress("avg_t",&flux_avg_t);
  // if(debug){
  // char *branchname_flux[26]={"hitn","id","pid","mpid","tid","mtid","otid","trackE","totEdep","trackE","avg_x","avg_y","avg_z","avg_lx","avg_ly","avg_lz","px","py","pz","vx","vy","vz","mvx","mvy","mvz","avg_t"};
  // cout << endl << "tree_flux" << endl;
  // for (Int_t i=0;i<26;i++) { 
  // cout << branchname_flux[i] << " " <<  tree_flux->GetBranch(branchname_flux[i])->GetLeaf(branchname_flux[i])->GetTypeName() << ",";
  // }
  // }
  
  TTree *tree_solid_ec = (TTree*) file->Get("solid_ec");
  setup_tree_solid_ec(tree_solid_ec);
  
  TTree *tree_solid_lgc = (TTree*) file->Get("solid_lgc");
  setup_tree_solid_lgc(tree_solid_lgc);

  double filenum = 1;
  if (input_filename.find("_filenum",0) != string::npos) 
    {
      filenum = atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
      cout << "filenum " << filenum << " for additional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
    }
  else 
    cout << "This file has no filenum, please check if you need filenum for additional normalization" << endl;

  

  // Event loop =====================================

  int nselected = 0;

  Double_t n_acc = 0; // accepted events
  Double_t n_trig[3] = {0, 0, 0};  // events triggering LGC, EC, both

  UInt_t nevmax = tree_generated->GetEntries();
  if (sev > nevmax)
    return;
  if (nev + sev > nevmax)
    nev = nevmax - sev;
  
  for (UInt_t i = sev; i < sev + nev; ++i) {
    
    cout << i << "\r";
    //   cout << i << "\n";
    
    tree_header->GetEntry(i);
    double rate=header_var8->at(0) / filenum;
    
    tree_generated->GetEntry(i);  
    int pid_gen=0;
    double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;      
    Double_t Q2 = 0.0;
    Double_t W2 = 0.0;
    Double_t xbj = 0.0;
    
    //       cout << "gen_pid->size() " << gen_pid->size() << endl;        
    for (UInt_t j=0;j<gen_pid->size();j++) {
      //       cout << gen_pid->at(j) << " " << gen_px->at(j) << endl;//<< " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
      pid_gen=gen_pid->at(j);
      px_gen=gen_px->at(j);
      py_gen=gen_py->at(j);
      pz_gen=gen_pz->at(j);
      vx_gen=gen_vx->at(j);
      vy_gen=gen_vy->at(j);
      vz_gen=gen_vz->at(j);
    }
    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    theta_gen=acos(pz_gen/p_gen);
    phi_gen=atan2(py_gen,px_gen);

    Q2 = 2 * ebeam*1e3 * p_gen * (1 - cos (theta_gen));
    W2 = (.938e3*.938e3+2*.938e3*(ebeam*1e3 - p_gen)-Q2);
    xbj = Q2 / 2 / 0.938e3 / (ebeam*1e3 - p_gen);
    Bool_t is_kine = (Q2 > q2min*1e6 && W2 > wmin*wmin*1e6 && xbj > xbjmin);
    if (kinecut && !is_kine)
      continue;
    
    hgen_ThetaP->Fill(theta_gen*DEG,p_gen/1e3, rate);
    hgen_Q2x->Fill (xbj, Q2/1e6, rate);


    // Process flux hits =============================
    
    bool Is_acc=false,Is_ec=false,Is_gem[5]={false,false,false,false,false},Is_lgc=false;
    
    tree_flux->GetEntry(i);  
    
    int sec_ec=0;
    Double_t p_ec = 0.;
    Double_t r_ec = 0.;
    Double_t phi_ec = 0.;

    for (UInt_t j=0;j<flux_hitn->size();j++) {
      //       cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
      
      int detector_ID=flux_id->at(j)/1000000;
      int subdetector_ID=(flux_id->at(j)%1000000)/100000;
      int subsubdetector_ID=((flux_id->at(j)%1000000)%100000)/10000;
      int component_ID=flux_id->at(j)%10000;      
      
      //       if (detector_ID==5 && subdetector_ID == 1 && subsubdetector_ID == 1)   cout << "particle mom entering SPD " << flux_trackE->at(j) << endl;   
      
      //       if (detector_ID==4 && subdetector_ID == 1 && subsubdetector_ID == 1)   cout << "particle mom entering MRPC " << flux_trackE->at(j) << endl;   
      
      //       if (detector_ID==3 && subdetector_ID == 1 && subsubdetector_ID == 1)   cout << "particle mom entering EC " << flux_trackE->at(j) << endl;         
      double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2));
      double_t hit_momen = sqrt (flux_px->at(j) *flux_px->at(j) +
				 flux_py->at(j) *flux_py->at(j) +
				 flux_pz->at(j) *flux_pz->at(j));
      double hit_y=flux_avg_y->at(j),hit_x=flux_avg_x->at(j),hit_z=flux_avg_z->at(j);          
      double hit_phi=atan2(hit_y,hit_x)*DEG;
      
      int hit_id=-1;
      if (detector_ID==1 && subdetector_ID == 1 && subsubdetector_ID == 1) hit_id=0;
      if (detector_ID==1 && subdetector_ID == 2 && subsubdetector_ID == 1) hit_id=1;	  
      if (detector_ID==1 && subdetector_ID == 3 && subsubdetector_ID == 1) hit_id=2;	  
      if (detector_ID==1 && subdetector_ID == 4 && subsubdetector_ID == 1) hit_id=3;	  
      if (detector_ID==1 && subdetector_ID == 5 && subsubdetector_ID == 1) hit_id=4;	  
      if (detector_ID==1 && subdetector_ID == 6 && subsubdetector_ID == 1) hit_id=5;	        
      if (detector_ID==2 && subdetector_ID == 1 && subsubdetector_ID == 1) hit_id=6;
      if (detector_ID==2 && subdetector_ID == 2 && subsubdetector_ID == 1) hit_id=7;	              
      if (detector_ID==5 && subdetector_ID == 1 && subsubdetector_ID == 1) hit_id=8;
      if (detector_ID==5 && subdetector_ID == 2 && subsubdetector_ID == 1) hit_id=9;	                          
      if (detector_ID==3 && subdetector_ID == 1 && subsubdetector_ID == 1) hit_id=10;
      if (detector_ID==3 && subdetector_ID == 2 && subsubdetector_ID == 1) hit_id=11;	                    
      
      if (0<=hit_id && hit_id<=11) hflux_hitxy[hit_id]->Fill(flux_avg_x->at(j)/10.,flux_avg_y->at(j)/10., rate);
      //       else cout << "flux_id->at(j) " << flux_id->at(j) << endl;
      
      //check hit on EC and find sec_ec
      if(hit_id==10 && flux_tid->at(j)==1){
	if (110<=hit_r/1e1 && hit_r/1e1<=250) { //cut on EC hit
	  Is_ec=true;
	  int sec_shift=1.7;  // shift to match electron turning in field
	  if (hit_phi > 90+sec_shift) sec_ec=int((hit_phi-90-sec_shift)/12+1);
	  else sec_ec=int((hit_phi+360-90-sec_shift)/12+1);
	  	// cout << " hit_phi " << hit_phi << " sec_ec " << sec_ec << endl;	
	  p_ec = sqrt (pow (flux_px->at(j),2) +
		       pow (flux_py->at(j),2) +
		       pow (flux_pz->at(j),2)) / 1000.0;  //covert MeV to GeV
	  r_ec = hit_r /  10.0;  //convert mm to cm
	  phi_ec = fabs(atan(flux_avg_y->at(j)/flux_avg_x->at(j))/3.1416*180);	
	}
      }
      
      //check hit on GEM
      if (detector_ID==1 && flux_tid->at(j)==1) {
	// some low mom tracks spiral and travel back in field and go through one plane twice or more,   flux bank average these steps to get hit position which is wrong and can be outside of the virtual plane.
	//       for PVDIS
	//  my @Rin = (48,59,65,105,109);
	//  my @Rout = (122,143,143,230,237);	      
	//       for SIDIS
	//        my @Rin = (36,21,25,32,42,55);
	//        my @Rout = (87,98,112,135,100,123);
	double Rin[6]={0,0,0,0,0,0},Rout[6]={0,0,0,0,0,0};
	//       if (Is_PVDIS){
	Rin[0]=48;Rin[1]=59;Rin[2]=65;Rin[3]=105;Rin[4]=109;Rin[5]=0;     
        Rout[0]=122;Rout[1]=143;Rout[2]=143;Rout[3]=230;Rout[4]=237;Rout[5]=300;
	//       }
	//       else {
	// 	Rin[0]=36;Rin[1]=21;Rin[2]=25;Rin[3]=32;Rin[4]=42;Rin[5]=55;     
	//         Rout[0]=87;Rout[1]=98;Rout[2]=112;Rout[3]=135;Rout[4]=100;Rout[5]=123;
	//       }
	
	if (Rin[subdetector_ID-1]<=hit_r/1e1 && hit_r/1e1<Rout[subdetector_ID-1]) {
	  Is_gem[subdetector_ID-1]=true;
	  // 	cout << flux_id->at(j) << endl; 	
	  continue;	
	}  
	
      }
      
      //check hit on LGC (flux counter)
      if (detector_ID==2 && flux_tid->at(j)==1) {
	// some low mom tracks spiral and travel back in field and go through one plane twice or more,   flux bank average these steps to get hit position which is wrong and can be outside of the virtual plane.
	double RinLGC = 65.0;
	double RoutLGC = 144.0;
	
	if (RinLGC<=hit_r/1e1 && hit_r/1e1<RoutLGC) 
	  {
	    Is_lgc = true;
	    continue;	
	  }  
	
      }
      
    } // end loop over flux hits

    // Process LGC ==================================
    
    //   tree_solid_ec->GetEntry(i);    
    //   double totEdep_ec=process_tree_solid_ec(tree_solid_ec);
    //   cout << "totEdep_ec " << totEdep_ec << endl;
    
    tree_solid_lgc->GetEntry(i);
    
    Int_t nphe_lgc[30]={0};
    process_tree_solid_lgc(tree_solid_lgc,nphe_lgc);
    
    Int_t ntrigsecs = 0;
    Bool_t is_lgc_trigger = false;
    Int_t trigger_lgc[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // Check for LGC triggers
    process_tree_solid_lgc_trigger (NULL, trigger_lgc, ntrigsecs);
    

    // Now handle event if it had an EC hit  ====================
    
    Bool_t is_ec_trigger = false;

    double EC_efficiency=0;
    if (Is_ec){
      //check how number of p.e. for all LGC sector and plot it at 0 when it's same sector like sec_ec
      for (int j=0;j<30;j++){    
	int lgcsec = j+1;
	int index;
	if(-14<=lgcsec-sec_ec && lgcsec-sec_ec<16) index=lgcsec-sec_ec;
	else if (lgcsec-sec_ec>=16) index=lgcsec-sec_ec-30;
	else if (lgcsec-sec_ec<-14) index=lgcsec-sec_ec+30;    
	//     if (nphe_lgc[index]<0) cout << sec_ec << " " << lgcsec << " " << index << " " << nphe_lgc[index] << endl;    
	hsectoring_ec_lgc->Fill(index,nphe_lgc[j]*rate);
      }	
      
      // Check triggering
      is_lgc_trigger = trigger_lgc[sec_ec-1] == 1;
	
      //using Rakitha's curve, no need to have phi dependent
      EC_efficiency = GetElectronTriggerEffi(GetRadiusIndex(r_ec), GetMomentumIndex(p_ec));

      // // using Jin's curves
      // TString region;		
      // if (phi_ec - int (phi_ec / 12.)*12<6) region="high";
      // else region="low";
      //       EC_efficiency =Get_PVDIS_trigger_efficiency_EC (region, "electron", r_ec, p_ec);

      double test_num = rand.Uniform(0,1);
      if (test_num <= EC_efficiency)
	{
	  is_ec_trigger = true;
	}    
      
      // cout << r_ec << " " << p_ec << " " << p_gen << " " << EC_efficiency
      // 	   <<" EC " << (is_ec_trigger?"T":"F") << endl;
	
    }
    
    if (Is_ec)
      {
	if(Is_gem[0] && Is_gem[1] && Is_gem[2] && Is_gem[3] && Is_gem[4])
	  {
	    //sum number of p.e. from LGC sector sec_ec-sec_width_sum to sec_ec+sec_width_sum
	    int sec_width_sum=0; 
	    int nphe_lgc_total=0;  
	    for (int j=sec_ec-sec_width_sum;j<=sec_ec+sec_width_sum;j++)
	      {    
		int index;
		if (0<j && j<=30) index=j-1;
		else if (j>30) index=j-30-1;
		else if (j<=0) index=j+30-1;
		else cout << "something wrong with sec" << endl;      
		nphe_lgc_total += nphe_lgc[index];
	      }      
	    hnphe_lgc->Fill(nphe_lgc_total,rate);  
	    
	    if (lgcpecut)
	      {
		// Require minimal LGC signal for acceptance
		if (nphe_lgc_total>nphecut){  // cut on lgc
		  Is_acc=true;
		}
	      }
	    else
	      {
		// Require LGC flux counter for acceptance
		if (Is_lgc)
		  Is_acc = true;
	      }
	  }      
      }

    if (Is_acc) 
      {
	hacceptance_ThetaP->Fill(theta_gen*DEG,p_gen/1e3,rate);  
	hacceptance_P->Fill(p_gen/1e3,rate);  
	hacceptance_RhitPhit->Fill(r_ec,p_ec,rate);  
	hdist_x->Fill (xbj, rate);
	hdist_Q2->Fill (Q2/1e6, rate);
	hacceptance_Q2x->Fill (xbj, Q2/1e6,rate);

	if (is_kine)
	  {
	    n_acc += rate;
	    if (is_lgc_trigger)
	      n_trig[0] += rate;
	    if (is_ec_trigger)
	      n_trig[1] += rate;
	  }

	if (is_lgc_trigger && is_ec_trigger)
	  {
	    if (is_kine)
	      n_trig[2] += rate;
	    hacceptanceeff_ThetaP->Fill(theta_gen*DEG,p_gen/1e3,rate);  
	    hacceptanceeff_P->Fill(p_gen/1e3,rate);  
	    hacceptanceeff_Q2x->Fill (xbj, Q2/1e6,rate);
	    hacceptanceeff_RhitPhit->Fill(r_ec,p_ec,rate);  
	  }
      }
    
  } // end event loop
  file->Close();

  if (kinecut) cout << "(with kine cuts:)" << endl;
  cout << "accepted    " << n_acc << endl;
  cout << "LGC trig    " << n_trig[0] << " " << float (n_trig[0]) / n_acc << endl;
  cout << "EC trig     " << n_trig[1] << " " << float (n_trig[1]) / n_acc << endl;
  cout << "LGC+EC trig " << n_trig[2] << " " << float (n_trig[2]) / n_acc << endl;
  
  TCanvas *c_flux_hitxy = new TCanvas("flux_hitxy","flux_hitxy",1800,900);
  c_flux_hitxy->Divide(5,2);
  for (Int_t i=0;i<m;i++) {
    c_flux_hitxy->cd(i+1);
    gPad->SetLogz();
    hflux_hitxy[i]->Draw("colz");
  }
  
  TCanvas *c_npe_lgc = new TCanvas("npe_lgc","npe_lgc",1600,900);
  c_npe_lgc->Divide(2,1);
  c_npe_lgc->cd(1);
  
  hnphe_lgc->SetFillColor(1);
  hnphe_lgc->SetBarWidth(0.5);
  hnphe_lgc->SetBarOffset(0.25);
  hnphe_lgc->Draw("b");
  
  Double_t intall = hnphe_lgc->Integral();
  Double_t intcut = hnphe_lgc->Integral(nphecut+2, hnphe_lgc->GetNbinsX());
  cout << intcut/intall << endl;
  
  TPaveText* tpt = new TPaveText (0.6, 0.5, 0.85, 0.6, "NDC");
  tpt->AddText (Form ("Eff = %4.1f%%", intcut/intall*100));
  tpt->Draw();
  
  TLine* tln = new TLine (hnphe_lgc->GetBinLowEdge(nphecut+2), 0, 
			 hnphe_lgc->GetBinLowEdge(nphecut+2), hnphe_lgc->GetMaximum()/10);
  tln->SetLineColor(2);
  tln->SetLineWidth(3);
  tln->SetLineStyle(2);
  tln->Draw();
  
  c_npe_lgc->cd(2);
  
  hsectoring_ec_lgc->SetFillColor(1);
  hsectoring_ec_lgc->SetBarWidth(0.5);
  hsectoring_ec_lgc->SetBarOffset(0.25);
  hsectoring_ec_lgc->Draw("b");
  gPad->SetLogy(1);
  
  TCanvas *c_acc = new TCanvas("acc","acc",1800,750);
  c_acc->Divide(2,1);
  
  c_acc->cd(1);
  hacceptance_ThetaP->Divide(hacceptance_ThetaP,hgen_ThetaP);  
  hacceptance_ThetaP->SetMinimum(0);  
  hacceptance_ThetaP->SetMaximum(1);    
  hacceptance_ThetaP->GetXaxis()->SetTitle ("#theta [deg]");
  hacceptance_ThetaP->GetXaxis()->SetTitleSize (0.06);
  hacceptance_ThetaP->GetYaxis()->SetTitle ("p [GeV/c]");
  hacceptance_ThetaP->GetYaxis()->SetTitleSize (0.06);
  hacceptance_ThetaP->Draw("colz");
  
  ContoursPT (ebeam, true);

  // TF1* fq2 = new TF1("fq2",Form("%f/(2*%f*(1-cos(x*3.14159/180)))", q2min, ebeam),10,50);
  // fq2->SetLineColor(2);
  // fq2->Write ("fq2");
  // fq2->Draw("same");
  // TLatex* tq2 = new TLatex (15,9, Form ("Q^{2}>%3.1f GeV^{2}", q2min));
  // tq2->SetTextColor(2);
  // tq2->SetTextSize(0.035);
  // tq2->Write ("tq2");
  // tq2->Draw();
  // TF1* fw2 = new TF1("fw2",Form("(.938*.938+2*.938*%f-%f)/(2*(.938+%f*(1-cos(x*3.14159/180))))", ebeam, wmin*wmin, ebeam),10,50);
  // fw2->SetLineColor(6);
  // fw2->Write ("fw2");
  // fw2->Draw("same");
  // TLatex* tw2 = new TLatex (11,5.5, Form("W>%3.1f GeV", wmin));
  // tw2->SetTextColor(6);
  // tw2->SetTextSize(0.035);
  // tw2->Write ("tw2");
  // tw2->Draw();
  // TF1* fxbj = new TF1("fxbj",Form("(.938*%f*%f)/(%f*(1-cos(x*3.14159/180))+.938*%f)", xbjmin, ebeam, ebeam, xbjmin),10,50);
  // fxbj->SetLineColor(1);
  // fxbj->Write ("fxbj");
  // fxbj->Draw("same");
  // TLatex* txbj = new TLatex (41,2.25,Form("x_{bj}>%4.2f", xbjmin));
  // txbj->SetTextColor(1);
  // txbj->SetTextSize(0.035);
  // txbj->Write ("txbj");
  // txbj->Draw();
  
  c_acc->cd(2);
  hacceptance_Q2x->Divide(hacceptance_Q2x,hgen_Q2x);  
  hacceptance_Q2x->SetMinimum(0);  
  hacceptance_Q2x->SetMaximum(1);    
  hacceptance_Q2x->GetXaxis()->SetTitle ("x_{bj}");
  hacceptance_Q2x->GetXaxis()->SetTitleSize (0.06);
  hacceptance_Q2x->GetYaxis()->SetTitle ("Q^{2} [(GeV/c)^{2}]");
  hacceptance_Q2x->GetYaxis()->SetTitleSize (0.06);
  hacceptance_Q2x->Draw("colz");

  ContoursQ2x (ebeam, true);
  
  TCanvas *c_acceff = new TCanvas("acceff","acceff",1800,750);
  c_acceff->Divide(2,1);
  
  c_acceff->cd(1);
  hacceptanceeff_ThetaP->Divide(hacceptanceeff_ThetaP,hgen_ThetaP);  
  hacceptanceeff_ThetaP->SetMinimum(0);  
  hacceptanceeff_ThetaP->SetMaximum(1);    
  hacceptanceeff_ThetaP->GetXaxis()->SetTitle ("#theta [deg]");
  hacceptanceeff_ThetaP->GetXaxis()->SetTitleSize (0.06);
  hacceptanceeff_ThetaP->GetYaxis()->SetTitle ("p [GeV/c]");
  hacceptanceeff_ThetaP->GetYaxis()->SetTitleSize (0.06);
  hacceptanceeff_ThetaP->Draw("colz");
  
  ContoursPT (ebeam);
  
  c_acceff->cd(2);
  hacceptanceeff_Q2x->Divide(hacceptanceeff_Q2x,hgen_Q2x);  
  hacceptanceeff_Q2x->SetMinimum(0);  
  hacceptanceeff_Q2x->SetMaximum(1);    
  hacceptanceeff_Q2x->GetXaxis()->SetTitle ("x_{bj}");
  hacceptanceeff_Q2x->GetXaxis()->SetTitleSize (0.06);
  hacceptanceeff_Q2x->GetYaxis()->SetTitle ("Q^{2} [(GeV/c)^{2}]");
  hacceptanceeff_Q2x->GetYaxis()->SetTitleSize (0.06);
  hacceptanceeff_Q2x->Draw("colz");

  ContoursQ2x (ebeam);

  TCanvas *c_acc2 = new TCanvas("acc_q2x","acc",900,750);
  hacceptance_Q2x->Draw("colz");

  ContoursQ2x (ebeam);

  TCanvas *c_acceff2 = new TCanvas("acceff_q2x","acceff",900,750);
  hacceptanceeff_Q2x->Draw("colz");

  ContoursQ2x (ebeam);

  TCanvas *c_eff = new TCanvas("eff","eff",900,750);
  
  TH2F* heff_ThetaP = new TH2F (*hacceptanceeff_ThetaP);
  heff_ThetaP->SetStats(0);
  heff_ThetaP->SetName ("eff_ThetaP_FA");
  heff_ThetaP->SetTitle ("Efficiency");
  heff_ThetaP->Divide(hacceptanceeff_ThetaP,hacceptance_ThetaP);
  heff_ThetaP->GetXaxis()->SetTitle ("#theta [deg]");
  heff_ThetaP->GetXaxis()->SetTitleSize (0.06);
  heff_ThetaP->GetYaxis()->SetTitle ("p [GeV/c]");
  heff_ThetaP->GetYaxis()->SetTitleSize (0.06);
  heff_ThetaP->Draw("colz");
  
  ContoursPT (ebeam);
  
  TCanvas *c_eff1 = new TCanvas("eff1","eff",900,750);
  
  TH1F* heff_P = new TH1F (*hacceptanceeff_P);
  heff_P->SetStats(0);
  heff_P->SetName ("eff_P_FA");
  heff_P->SetTitle ("Efficiency");
  heff_P->Divide(hacceptanceeff_P,hacceptance_P);
  heff_P->GetXaxis()->SetTitle ("p [GeV/c]");
  heff_P->GetXaxis()->SetTitleSize (0.06);
  heff_P->Draw();

  TCanvas *c_effrhph = new TCanvas("effrhph","effrhph",900,750);
  
  TH2F* heff_RhitPhit = new TH2F (*hacceptanceeff_RhitPhit);
  heff_RhitPhit->SetStats(0);
  heff_RhitPhit->SetName ("eff_RhitPhit_FA");
  heff_RhitPhit->SetTitle ("Efficiency");
  heff_RhitPhit->Divide(hacceptanceeff_RhitPhit,hacceptance_RhitPhit);
  heff_RhitPhit->GetXaxis()->SetTitle ("r_{hit} [cm]");
  heff_RhitPhit->GetXaxis()->SetTitleSize (0.06);
  heff_RhitPhit->GetYaxis()->SetTitle ("p_{hit} [GeV/c]");
  heff_RhitPhit->GetYaxis()->SetTitleSize (0.06);
  heff_RhitPhit->Draw("colz");


  outputfile->Write();
  outputfile->Flush();
  
}
