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
#include <TRandom3.h>
#include <TLegend.h>

#include "analysis_tree_solid_hgc.C"
#include "analysis_tree_solid_lgc.C"
#include "analysis_tree_solid_ec.C"
#include "../ec_layout/ec_layout.h"

using namespace std;

void analysis_aiml(string input_filename,string det,bool Is_simsafe=false,bool Is_textfile=false, bool Is_decaycut=false, int factor_background=0,string the_dir=".",int Nreduce=1)
{
gROOT->Reset();
// gStyle->SetPalette(57);
gStyle->SetOptStat(0);
// gStyle->SetOptStat(1);
// gStyle->SetOptStat("ioue");
// gStyle->SetOptStat(11111111);
gStyle->SetOptFit(111111);

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

TRandom3 rand;
rand.SetSeed(0);		  

const double DEG=180./3.1415926;

LoadEC_map("../ec_layout");
const int ecmodule=1748;

int vplaneid;
int ch_pmt;
if (det=="lgc"){
  vplaneid=2110000;
  ch_pmt=270;
}
else if (det=="lgcec") {
  vplaneid=3110000;  //FA only
  ch_pmt=270;
}
else if (det=="hgc") {
  vplaneid=2210000;  
  ch_pmt=480;
}
else {cout << "bad detector name" << endl; return; }

//   ch_quad=1920; ch_pixel=30720; for hgc
//   ch_quad=1080; ch_pixel=17280; for lgc
int ch_quad=ch_pmt*4,ch_pixel=ch_pmt*64;
int sen_pmt=ch_pmt/30,sen_quad=ch_quad/30,sen_pixel=ch_pixel/30;
int sentrans_pmt = sqrt(sen_pmt),sentrans_quad = sqrt(sen_quad),sentrans_pixel = sqrt(sen_pixel);

// Float_t npe[ch],npe[ch];
Float_t npe_total=0,npe_total_trigged=0;
	
//occupancy threshold
double occ_threshold=0; //in N_p.e.

double PEthresh=2; // pe threshold for each pmt
double PMTthresh=2; // pmt threshold, at least 2pmts are fired in each sector

// cout << eff_PMT[0] << " " << eff_PMT[1] << " " << eff_PMT[2] << " " << eff_PMT[3] << " " << endl;

// std::size_t found = input_filename.rfind("cache");
// string input_filename_original=input_filename;
// if (found!=std::string::npos) input_filename.replace(found,5,"work");

char the_filename[400];
string input_filename_short=input_filename.substr(input_filename.rfind("/")+1);
sprintf(the_filename, "%s",input_filename_short.substr(0,input_filename_short.find(".root")).c_str());
sprintf(the_filename,Form("%s/%s",the_dir.c_str(),the_filename));
// cout << " the_filename " << the_filename << endl;

string solid_setup;
if(input_filename.rfind("SIDIS_He3") !=std::string::npos) solid_setup="SIDIS_He3";
else if(input_filename.rfind("PVDIS_LD2") !=std::string::npos) solid_setup="PVDIS_LD2";
else{cout <<"unknown solid_setup "<< solid_setup<<endl; return;}

ofstream textfile_row_npe,textfile_row_pmt,textfile_row_quad,textfile_row_pixel;
char textfile_name[100];
sprintf(textfile_name,"%s_row_npe.csv",the_filename);
if (Is_textfile) textfile_row_npe.open(textfile_name,std::ofstream::trunc);
sprintf(textfile_name,"%s_row_pmt.csv",the_filename);
if (Is_textfile) textfile_row_pmt.open(textfile_name,std::ofstream::trunc);
sprintf(textfile_name,"%s_row_quad.csv",the_filename);
if (Is_textfile) textfile_row_quad.open(textfile_name,std::ofstream::trunc);
sprintf(textfile_name,"%s_row_pixel.csv",the_filename);
if (Is_textfile) textfile_row_pixel.open(textfile_name,std::ofstream::trunc);

int Ngood_cc=0,Nrand_cc=0;
int Ngood_ec=0,Nrand_ec=0;
TFile *file_cc_bg,*file_ec_bg;
TH2F *hevent_hit_npe_bg, *hevent_hit_pmt_bg, *hevent_hit_quad_bg, *hevent_hit_pixel_bg;
TH2F *hevent_ec_shower_bg,*hevent_ec_preshower_bg;
if(factor_background){
// here is an example what I think how to merge signal and background
// 1.assume beam current 15ua and 50ns time window. we have 15e-6/1.6e-19*50e-9=4.7e6 beam e- within 50ns time window
// 2.assume beamontarget file is from 1e10 e- with 100x actual target density and it has 1e5 number of events with hits in at least one detector.
// 3.beamontarget file has 1e10/4.7e6 = 2127 number of time window.
// for each signal particle, we should add 1e5/100/2127=0.5 background event from beamontarget file randomly for every 4.detector channel, then output in csv.
// 5.for GEM, we need to use 300ns time window. all other detector use 50ns. so GEM needs to add 6 time more background events for each signal event
  
  cout << "add background" << endl;
  if (det=="lgc" || det=="lgcec"){
	if(solid_setup=="SIDIS_He3"){
//  /cache/halla/solid/sim/solid_gemc/SIDIS_He3_JLAB_VERSION_2.5/pass3/solid_SIDIS_He3_moved_full_BeamOnTarget_5e10_filter.root    
    //     total filtered event 34324461
    file_cc_bg=new TFile("/group/solid/www/solid/html/files/AIML/background/solid_SIDIS_He3_moved_full_BeamOnTarget_5e10_filter_output_lgc.root");
    Ngood_cc=9992; //good events of lgc background
    Nrand_cc=int(1*factor_background); //for lgc background within 50ns, 9992/(5e10/(15e-6/1.6e-19*50e-9))=0.94
    file_ec_bg=new TFile("/group/solid/www/solid/html/files/AIML/background/solid_SIDIS_He3_moved_full_BeamOnTarget_5e10_filter_output_ec_merge100.root");
    Ngood_ec=48576; //4857574 good event of ec background, after merging 100
    Nrand_ec=int(5*factor_background); //for ec background within 50ns, 48576/(5e10/(15e-6/1.6e-19*50e-9))=4.6
	}
	else if(solid_setup=="PVDIS_LD2"){
// /cache/halla/solid/sim/solid_gemc/PVDIS_JLAB_VERSION_2.5/pass3/solid_PVDIS_LD2_moved_full_BeamOnTarget_1e10_filter.root        
    //     total filtered event  520781336
    file_cc_bg=new TFile("/group/solid/www/solid/html/files/AIML/background/solid_PVDIS_LD2_moved_full_BeamOnTarget_1e10_filter_output_lgc.root");
    Ngood_cc=3557;
    Nrand_cc=int(6*factor_background); //for lgc background within 50ns, 3557/(1e10/(50e-6/1.6e-19*50e-9))=5.6
    file_ec_bg=new TFile("/group/solid/www/solid/html/files/AIML/background/solid_PVDIS_LD2_moved_full_BeamOnTarget_1e10_filter_output_ec_merge100.root");
    Ngood_ec=13925; //1392446 good event of ec background, after merging 100
    Nrand_ec=int(22*factor_background); //for ec background within 50ns, 13925/(1e10/(50e-6/1.6e-19*50e-9))=22
	}
	else{cout <<"unknown solid_setup "<< solid_setup<<endl; return;}
  }  
  else if (det=="hgc") {
//     file_cc_bg=new TFile("solid_SIDIS_He3_moved_full_BeamOnTarget_5e9_filter_output_hgc.root");    
// /cache/halla/solid/sim/solid_gemc/SIDIS_He3_JLAB_VERSION_1.3/pass8/farm_solid_SIDIS_He3_moved_BeamOnTarget_0.561e10_skim_HGCwinCF1.root
//     total filtered event 3990487    
    file_cc_bg=new TFile("/group/solid/www/solid/html/files/AIML/background/farm_solid_SIDIS_He3_moved_BeamOnTarget_0.561e10_skim_HGCwinCF1_output.root");
    Ngood_cc=9000; //good events of hgc background
    Nrand_cc=int(7.5*factor_background); //for hgc background within 50ns, 9000/(0.561e10/(15e-6/1.6e-19*50e-9))=7.5    
  }
  else {cout << "bad detector name" << endl; return;}
 
  hevent_hit_npe_bg=(TH2F*) file_cc_bg->Get("hevent_hit_npe");
  hevent_hit_pmt_bg=(TH2F*) file_cc_bg->Get("hevent_hit_pmt");	
  hevent_hit_quad_bg=(TH2F*) file_cc_bg->Get("hevent_hit_quad");	
  hevent_hit_pixel_bg=(TH2F*) file_cc_bg->Get("hevent_hit_pixel");	  
  
  if(det=="lgcec"){
    hevent_ec_shower_bg=(TH2F*) file_ec_bg->Get("hevent_ec_shower");
    hevent_ec_preshower_bg=(TH2F*) file_ec_bg->Get("hevent_ec_preshower");	      
  }
  
}

char output_filename[200];
sprintf(output_filename, "%s_output.root",the_filename);
TFile *outputfile=new TFile(output_filename, "recreate");

TH1F *hnpe_no0=new TH1F("hnpe_no0",";Npe;event count",200,0,200);

TH1F *hnpe=new TH1F("hnpe",";Npe;event count",200,0,200);

TH1F *hnpe_QE=new TH1F("hnpe_QE",";Npe/Nphoton;event count",100,0,1);

TH1F *hnpe_3secratio=new TH1F("hnpe_3secratio",";ratio;percent",120,0,1.2);
TH1F *hnpe_2secratioleft=new TH1F("hnpe_2secratioleft",";ratio;percent",120,0,1.2);
TH1F *hnpe_2secratioright=new TH1F("hnpe_2secratioright",";ratio;percent",120,0,1.2);
TH1F *hnpe_1secratio=new TH1F("hnpe_1secratio",";ratio;percent",120,0,1.2);

TH2F *hnpe_pmt_2D=new TH2F("hnpe_pmt_2D","Npe; sensor_x; sensor_y",3*sentrans_pmt,0,3*sentrans_pmt,sentrans_pmt,0,sentrans_pmt);
TH2F *hnpe_quad_2D=new TH2F("hnpe_quad_2D","Npe; sensor_x; sensor_y",3*sentrans_quad,0,3*sentrans_quad,sentrans_quad,0,sentrans_quad);
TH2F *hnpe_pixel_2D=new TH2F("hnpe_pixel_2D","Npe; sensor_x; sensor_y",3*sentrans_pixel,0,3*sentrans_pixel,sentrans_pixel,0,sentrans_pixel);

TH2F *hnpe_pmt_2D_phi[12];
TH2F *hnpe_quad_2D_phi[12];
TH2F *hnpe_pixel_2D_phi[12];

for(int i=0;i<12;i++){
   char hstname[100];   
   sprintf(hstname,"hnpe_pmt_2D_phi%i",i);
   hnpe_pmt_2D_phi[i]=new TH2F(hstname,"Npe; sensor_x; sensor_y",3*sentrans_pmt,0,3*sentrans_pmt,sentrans_pmt,0,sentrans_pmt);
   sprintf(hstname,"hnpe_quad_2D_phi%i",i);   
   hnpe_quad_2D_phi[i]=new TH2F(hstname,"Npe; sensor_x; sensor_y",3*sentrans_quad,0,3*sentrans_quad,sentrans_quad,0,sentrans_quad);
   sprintf(hstname,"hnpe_pixel_2D_phi%i",i);
   hnpe_pixel_2D_phi[i]=new TH2F(hstname,"Npe; sensor_x; sensor_y",3*sentrans_pixel,0,3*sentrans_pixel,sentrans_pixel,0,sentrans_pixel);
}


// TH1F *hnpe[2];
// for(int i=0;i<2;i++){
//    char hstname[100];   
//    sprintf(hstname,"npe_%i",i);
//   hnpe[i]=new TH1F(hstname,";Npe;event count",100,0,100);
// }

TH1F *hpe=new TH1F("pe","N_{p.e.};sensor;",sen_pmt,0,sen_pmt);
TH2F *hpe_2D=new TH2F("pe_2D","N_{p.e.};sensor_lx;sensor_ly",sentrans_pmt,0,sentrans_pmt,sentrans_pmt,0,sentrans_pmt);	

	TH1F *htime_photon=new TH1F("time_photon","time_photon;t (ns)",100,0,100);
	TH1F *hwl_photon=new TH1F("wl_photon","wavelength_photon;wl (nm)",1000,0,1000);	
	TH1F *hr_photon=new TH1F("r_photon","radius_photon;radius (cm);",100,0,10);

	TH1F *hmotherP=new TH1F("motherP","motherP;P (GeV);",11000,0,11);

	TH1F *hmother_procid=new TH1F("mother_procid","mother_procid",150,-0.5,149.5);

// const int Nbin_Theta=10,Nbin_Phi=18;
const int Nbin_Theta=100,Nbin_Phi=180;
double Max_Theta=50;
TH1F *hnpe_ThetaPhi[Nbin_Theta][Nbin_Phi];
for(int i=0;i<Nbin_Theta;i++){
  for(int j=0;j<Nbin_Phi;j++){
   char hstname[100];   
   sprintf(hstname,"hnpe_ThetaPhi_%i_%i",i,j);
   hnpe_ThetaPhi[i][j]=new TH1F(hstname,";number of photoelectron;count",100,0,50);
  }
}
TH2F *havg_pe=new TH2F("havg_pe","avg number of photoelectron;#theta(deg);#phi(deg)",Nbin_Theta,0,Max_Theta,Nbin_Phi,-180,180);

TH2F *hhitxy=new TH2F("hhitxy","p.e. pattern; r (mm); #phi (mm)",32,-102,102,32,-102,102);

TH2F *hhitxy_mirror=new TH2F("hhitxy_mirror","photon pattern on mirror; r (mm); #phi (mm)",140,800,2200,60,-300,300);
TH2F *hhitxy_cone=new TH2F("hhitxy_cone","photon pattern on cone; r (mm); #phi (mm)",60,-300,300,60,-300,300);

TH1F *hedep_ec_shower_sum=new TH1F("hedep_ec_shower_sum","hedep_ec_shower_sum;Edep(MeV);count",1000,0,1000);
TH1F *hedep_ec_preshower_sum=new TH1F("hedep_ec_preshower_sum","hedep_ec_preshower_sum;Edep(MeV);count",200,0,200);
TH2F *hedep_ec_sum=new TH2F("hedep_ec_sum","hedep_ec_sum;shower+preshower Edep(MeV);preshower Edep(MeV)",1000,0,1000,200,0,200);
// TH2F *hedep_ec_sum=new TH2F("hedep_ec_sum","hedep_ec_sum;shower+preshower Edep(MeV);preshower Edep(MeV)",100,0,1,200,0,200);
TH3F *hnpe_lgc_edep_ec_sum=new TH3F("hnpe_lgc_edep_ec_sum","hnpe_lgc_edep_ec_sum;shower+preshower Edep(MeV);preshower Edep(MeV);Npe LGC",1000,0,1000,200,0,200,200,0,200);
TH1F *hedep_ec_shower_cluster=new TH1F("hedep_ec_shower_cluster","hedep_ec_shower_cluster;Edep(MeV);count",1000,0,1000);
TH1F *hedep_ec_preshower_cluster=new TH1F("hedep_ec_preshower_cluster","hedep_ec_preshower_cluster;Edep(MeV);count",200,0,200);
TH2F *hedep_ec_cluster=new TH2F("hedep_ec_cluster","hedep_ec_cluster;shower+preshower Edep(MeV);preshower Edep(MeV)",1000,0,1000,200,0,200);
// TH2F *hedep_ec_cluster=new TH2F("hedep_ec_cluster","hedep_ec_cluster;shower+preshower Edep(MeV);preshower Edep(MeV)",100,0,1,200,0,200);
TH3F *hnpe_lgc_edep_ec_cluster=new TH3F("hnpe_lgc_edep_ec_cluster","hnpe_lgc_edep_ec_cluster;shower+preshower Edep(MeV);preshower Edep(MeV);Npe LGC",1000,0,1000,200,0,200,200,0,200);
TH1F *hedep_ec_clusterquality=new TH1F("hedep_ec_clusterquality","hedep_ec_clusterquality;Edep cluster/sum;count",100,0,1);


//save background histogram
// int cceventmax=10000; //save large histograms will crash, only enable this when using it
int cceventmax=1;
TH2F *hevent_hit_npe=new TH2F("hevent_hit_npe","Npe of all events; npe; event",1,0,1,cceventmax,0,cceventmax);
TH2F *hevent_hit_pmt=new TH2F("hevent_hit_pmt","Npe of all events; pmt; event",ch_pmt,0,ch_pmt,cceventmax,0,cceventmax);
TH2F *hevent_hit_quad=new TH2F("hevent_hit_quad","Npe of all events; quad; event",ch_quad,0,ch_quad,cceventmax,0,cceventmax);
TH2F *hevent_hit_pixel=new TH2F("hevent_hit_pixel","Npe of all events; pixel; event",ch_pixel,0,ch_pixel,cceventmax,0,cceventmax);

// int eceventmax=50000; //save large histograms will crash, only enable this when using it, sidis_he3
// int eceventmax=20000; //save large histograms will crash, only enable this when using it, pvdis_ld2
int eceventmax=1;
TH2F *hevent_ec_shower=new TH2F("hevent_ec_shower","EC shower Edep of all events; id; event",ecmodule,0,ecmodule,eceventmax,0,eceventmax);
TH2F *hevent_ec_preshower=new TH2F("hevent_ec_preshower","EC preshower Edep of all events; id; event",ecmodule,0,ecmodule,eceventmax,0,eceventmax);

TFile *file=new TFile(input_filename.c_str());
if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);
}
else cout << "open file " << input_filename << endl;    

TTree *tree_header = (TTree*) file->Get("header");
vector <int> *evn=0,*evn_type=0;
vector <double> *beamPol=0;
vector <int> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
tree_header->SetBranchAddress("evn",&evn);
tree_header->SetBranchAddress("evn_type",&evn_type);
tree_header->SetBranchAddress("beamPol",&beamPol);
// tree_header->SetBranchAddress("var1",&var1);
// tree_header->SetBranchAddress("var2",&var2);
// tree_header->SetBranchAddress("var3",&var3);
// tree_header->SetBranchAddress("var4",&var4);
// tree_header->SetBranchAddress("var5",&var5);
// tree_header->SetBranchAddress("var6",&var6);
// tree_header->SetBranchAddress("var7",&var7);
// tree_header->SetBranchAddress("var8",&var8);

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

TTree *tree_solid_ec, *tree_solid_ec_ps;
tree_solid_ec = (TTree*) file->Get("solid_ec");
setup_tree_solid_ec(tree_solid_ec);
tree_solid_ec_ps = (TTree*) file->Get("solid_ec_ps");
setup_tree_solid_ec_ps(tree_solid_ec_ps);

TTree *tree_solid_cc;
if (det=="lgc" || det=="lgcec"){
tree_solid_cc = (TTree*) file->Get("solid_lgc");
setup_tree_solid_lgc(tree_solid_cc);
}
else if (det=="hgc") {
tree_solid_cc = (TTree*) file->Get("solid_hgc");
setup_tree_solid_hgc(tree_solid_cc);
}
else {cout << "bad detector name" << endl; }

int nevent = (int)tree_generated->GetEntries();
int nselected = 0;
cout << "nevent " << nevent << endl;

int sensor_good=0;
int event_cc_good=0,event_ec_good=0,event_ec_good_merge=0;
int decay_no=0,decay_first=0,decay_lost=0;
double edep_ec_preshower_FA_tmp[ecmodule]={0},edep_ec_shower_FA_tmp[ecmodule]={0};


for (Int_t i=0;i<nevent/Nreduce;i++) { 
// for (Int_t i=0;i<nevent/100;i++) {
// for (Int_t i=0;i<1e4;i++) { 
  cout << "event " << i << "\r";
//   cout << "event " << i << "\n";

  tree_header->GetEntry(i);
  
  tree_generated->GetEntry(i);  
  int pid_gen=0;
  double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;      
  for (int j=0;j<gen_pid->size();j++) {
//       cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
      pid_gen=gen_pid->at(j);
      px_gen=gen_px->at(j)/1e3;
      py_gen=gen_py->at(j)/1e3;
      pz_gen=gen_pz->at(j)/1e3;
      vx_gen=gen_vx->at(j)/1e1;
      vy_gen=gen_vy->at(j)/1e1;
      vz_gen=gen_vz->at(j)/1e1;
      p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
      theta_gen=acos(pz_gen/p_gen)*DEG;
      phi_gen=atan2(py_gen,px_gen)*DEG;
      
//       cout << "gen " << p_gen << " "  << theta_gen <<  " " << phi_gen << endl;
  }

    tree_flux->GetEntry(i);  
    bool Is_decay=true,Is_decayfirst=false,Is_decaylost=false;
    double hit_px=0,hit_py=0,hit_pz=0,hit_p=0,hit_x=0,hit_y=0,hit_z=0,hit_r=0,hit_phi=0,hit_phip=0;
    for (Int_t j=0;j<flux_hitn->size();j++) {
//       cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  

      int detector_ID=flux_id->at(j)/1000000;
      int subdetector_ID=(flux_id->at(j)%1000000)/100000;
      int subsubdetector_ID=((flux_id->at(j)%1000000)%100000)/10000;
      int component_ID=flux_id->at(j)%10000;      
           
//       if (detector_ID==3 && subdetector_ID == 1 && subsubdetector_ID == 1)   cout << "particle mom entering EC " << flux_trackE->at(j) << endl;         
      
//       if (flux_pid->at(j)==0){       
// 	if (flux_id->at(j)==2211023) hhitxy_mirror->Fill(flux_avg_lx->at(j),flux_avg_ly->at(j));	
// 	if (flux_id->at(j)==2212023) hhitxy_cone->Fill(flux_avg_lx->at(j),flux_avg_ly->at(j));	
//       }
      
      
      if (flux_id->at(j)==vplaneid){  //for solid lgc or hgc front virtual plane
	if (flux_tid->at(j)==1) {  //  original track
	  decay_no++;
	  
	  hit_px=flux_px->at(j)/1e3;
	  hit_py=flux_py->at(j)/1e3;
	  hit_pz=flux_pz->at(j)/1e3;	  
	  hit_p=sqrt(hit_px*hit_px+hit_py*hit_py+hit_pz*hit_pz);
	  hit_x=flux_avg_x->at(j)/10.;
	  hit_y=flux_avg_y->at(j)/10.;	  
	  hit_z=flux_avg_z->at(j)/10.;
	  hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2));
	  hit_phi=atan2(hit_y,hit_x)*DEG;       //rad to deg		  	  	  
	  hit_phip=atan2(hit_py,hit_px)*DEG;       //rad to deg			  
// 	   cout << "hit_p " << hit_p << " "  << hit_phi <<  " " << hit_r << endl;

// 	  Is_decay=false;
	  if (hit_p/p_gen>0.8) Is_decay=false; // this mom cut reduce particles lost most mom
	}	
	else if (flux_mtid->at(j)==1 && (abs(flux_pid->at(j))==211 || abs(flux_pid->at(j))==13 || abs(flux_pid->at(j))==11) ) {
	  decay_first++;
// 	  if (-3700<flux_vz->at(j) && flux_vz->at(j)<-3300 && fabs(flux_vx->at(j))<50 && fabs(flux_vy->at(j))<50) Is_decayfirst=true;
	  Is_decayfirst=true;	  
	}
	else{}
// 	else  cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << endl;

// 	  cout << "hit Mom " <<  hit_px << ","<< hit_py<< ","<<hit_pz<<","<<hit_phip<< endl;
// 	  cout << "hit pos " << hit_x << ","<<hit_y<< ","<<hit_z<<","<<hit_phi<< endl;
      }
      
    }// end of flux	    
    
	if (Is_decay && !Is_decayfirst){
// 	if (Is_decay && Is_decayfirst){
	  decay_lost++;
	  Is_decaylost=true;
	  
// 	  for (Int_t j=0;j<flux_hitn->size();j++) {
// 	    cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
// 	  }
	  
	}
    
    if (Is_decaycut && Is_decay) continue;	 // skip decay particle

		//-----------------------	
		//--- lgc or hgc
		//-----------------------
		tree_solid_cc->GetEntry(i);

		//use hgc value for array size, as lgc values are smaller
    		double hit_pmt[480]={0};		
		double hit_quad[1920]={0};
		double hit_pixel[30720]={0};						
		int trigger[30]={0};
		int ntrigsecs=0;
		int photon_mtid=0;
		int nphoton=0;

		if (det=="lgc" || det=="lgcec") process_tree_solid_lgc(tree_solid_cc,Is_simsafe,hit_pmt,hit_quad,hit_pixel,trigger,ntrigsecs,nphoton,PMTthresh,PEthresh);
		else if (det=="hgc") process_tree_solid_hgc(tree_solid_cc,Is_simsafe,hit_pmt,hit_quad,hit_pixel,trigger,ntrigsecs,PMTthresh,PEthresh);
		else {cout << "bad detector name" << endl; return;}

		npe_total=0;
		for(int index=0;index<ch_pmt;index++) {
		  npe_total += hit_pmt[index];
// 		  if (hit_pmt[index]>0) cout << "pmt index " << index << "  " << hit_pmt[index] << endl;		  
		}
// 		cout << " npe_total " << npe_total << endl;
		
		//find sector by max Npe		
// 		int sec=0;		
// 		double hit_sec[30]={0};		  		  
// 		for(int index_sec=0;index_sec<30;index_sec++){
// 		  for(int index_pmt=0;index_pmt<16;index_pmt++){
// 		    hit_sec[index_sec] += hit_pmt[16*index_sec+index_pmt];
// 		  }
// 		}				
// 		for(int index_sec=0;index_sec<30;index_sec++){
// 		  if (hit_sec[index_sec]>hit_sec[sec]) sec=index_sec;
// 		}
		
		//find sector by track	
		int sec=0;
		if (hit_phi>=90) sec=int((hit_phi-90)/12);
		else sec=int((hit_phi+360-90)/12);
// 		cout << hit_phi << " " << sec << endl;
		if (sec<0 || sec > 29) {cout<< "bad sec " << sec << " " << hit_phi << endl; continue;}
		
		int sec_left=sec-1,sec_right=sec+1;
		if (sec_left==-1) sec_left=29;
		if (sec_right==30) sec_right=0;
// 		cout << " npe_total " << npe_total << endl;
// 		<<  hit_sec[sec_left]+hit_sec[sec]+hit_sec[sec_right]  << " " <<  hit_sec[sec_left] << " " <<  hit_sec[sec]  << " " <<  hit_sec[sec_right] << endl;
			
		//add background
		if(factor_background){
		  //add npe_total
// 		  cout << "npe_total before " << npe_total << endl;
// 		  if (rand.Uniform(0,1)<0.25){  // the chance to have a good lgc event in 1 sector within 50ns
// 		  if (rand.Uniform(0,1)<0.5){  // the chance to have a good lgc event in 2 sector within 50ns
// 		  if (rand.Uniform(0,1)<0.75){  // the chance to have a good lgc event in 3 sector within 50ns	    	    
// 		    int ibgevent=int(rand.Uniform(0,1)*9000);
// // 		    cout<< "ibgevent " << ibgevent << " npe_bg " << hevent_hit_npe_bg->GetBinContent(1,ibgevent) << endl;	
// 		    npe_total += hevent_hit_npe_bg->GetBinContent(1,ibgevent);
// 		  }
		  		  
		  //add pmt,quad,pixel	  
		  for(int ibgcount=0;ibgcount<Nrand_cc;ibgcount++){
		    int ibgevent=int(rand.Uniform(0,1)*Ngood_cc);
		    for(int index=0;index<ch_pmt;index++) hit_pmt[index] += hevent_hit_pmt_bg->GetBinContent(index+1,ibgevent);
		    for(int index=0;index<ch_quad;index++) hit_quad[index] += hevent_hit_quad_bg->GetBinContent(index+1,ibgevent);
		    for(int index=0;index<ch_pixel;index++) hit_pixel[index] += hevent_hit_pixel_bg->GetBinContent(index+1,ibgevent);
		  }		  
		  
		  //add npe for 3 sectors
		  int Nsensor=sen_pmt;
		  int npe_tmp=0;
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec_left+index];
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec+index];
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec_right+index]; 		  
		  npe_total=npe_tmp;
		}

		int Nsensor=sen_pmt;
		int npe_tmp_left=0,npe_tmp=0,npe_tmp_right=0;
		for(int index=0;index<Nsensor;index++) npe_tmp_left += hit_pmt[Nsensor*sec_left+index];
		for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec+index];
		for(int index=0;index<Nsensor;index++) npe_tmp_right += hit_pmt[Nsensor*sec_right+index]; 		  
		
		hnpe_3secratio->Fill(double(npe_tmp_left+npe_tmp+npe_tmp_right)/double(npe_total));
		hnpe_2secratioleft->Fill(double(npe_tmp_left+npe_tmp)/double(npe_total));		  hnpe_2secratioright->Fill(double(npe_tmp+npe_tmp_right)/double(npe_total));
		hnpe_1secratio->Fill(double(npe_tmp)/double(npe_total));		
		
		if (npe_total>0) hnpe_no0->Fill(npe_total);
				
		hnpe->Fill(npe_total);	
		
		hnpe_QE->Fill(double(npe_total)/double(nphoton));
		
		//tdc for pixel
		for(int index=0;index<ch_pixel;index++) if(hit_pixel[index]>0) hit_pixel[index]=1;  
		
		//save background
		if (npe_total>0){
		  event_cc_good++;
// 		  cout<< " ok " << endl;
		  hevent_hit_npe->Fill(0.5,event_cc_good,npe_total);		  
		  for(int index=0;index<ch_pmt;index++)	hevent_hit_pmt->Fill(index,event_cc_good,hit_pmt[index]);		
		  for(int index=0;index<ch_quad;index++) hevent_hit_quad->Fill(index,event_cc_good,hit_quad[index]);	
		  for(int index=0;index<ch_pixel;index++) hevent_hit_pixel->Fill(index,event_cc_good,hit_pixel[index]);  
		}

		//fill 2D plot
		{
// 		  int bin_phi=int(phi_gen+180)%12;
		  int bin_phi=int(hit_phi+180)%12;
// 		  cout<< " hit_phi " << hit_phi<< " bin_phi " << bin_phi<< endl;		  
		  if (bin_phi<0 || bin_phi > 11) {cout<< "bad bin_phi " << bin_phi<< endl; continue;}
		  
		  int Nsensor;
		  Nsensor=ch_pmt/30;
		  for(int index=0;index<Nsensor;index++) {
		    int thetax=index/int(sqrt(Nsensor));  //along solid theta angle direction
		    hnpe_pmt_2D->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec_left+index]);
		    hnpe_pmt_2D->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec+index]);
		    hnpe_pmt_2D->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec_right+index]);
		    
		    hnpe_pmt_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec_left+index]);
		    hnpe_pmt_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec+index]);
		    hnpe_pmt_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pmt[Nsensor*sec_right+index]);		    
		  }
		  Nsensor=ch_quad/30;
		  for(int index=0;index<Nsensor;index++) {
		    hnpe_quad_2D->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec_left+index]);
		    hnpe_quad_2D->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec+index]);
		    hnpe_quad_2D->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec_right+index]);	
		    
		    hnpe_quad_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec_left+index]);
		    hnpe_quad_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec+index]);
		    hnpe_quad_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_quad[Nsensor*sec_right+index]);	
		  }
		  Nsensor=ch_pixel/30;
		  for(int index=0;index<Nsensor;index++) {
		    hnpe_pixel_2D->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec_left+index]);
		    hnpe_pixel_2D->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec+index]);
		    hnpe_pixel_2D->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec_right+index]);	
		    
		    hnpe_pixel_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec_left+index]);
		    hnpe_pixel_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec+index]);
		    hnpe_pixel_2D_phi[bin_phi]->Fill(index/int(sqrt(Nsensor))+2*int(sqrt(Nsensor)),index%int(sqrt(Nsensor)),hit_pixel[Nsensor*sec_right+index]);			    
		  }
		}
		
		//-----------------------	
		//--- ec
		//-----------------------
		tree_solid_ec->GetEntry(i);
		tree_solid_ec_ps->GetEntry(i);

		double edep_ec_preshower_FA[ecmodule]={0},edep_ec_shower_FA[ecmodule]={0};
		double edep_ec_preshower_LA[ecmodule]={0},edep_ec_shower_LA[ecmodule]={0};				
		process_tree_solid_ec(tree_solid_ec,tree_solid_ec_ps,edep_ec_preshower_FA,edep_ec_preshower_LA,edep_ec_shower_FA,edep_ec_shower_LA);
		
		//add background
		if(factor_background){
		  if(det=="lgcec"){
		    for(int ibgcount=0;ibgcount<Nrand_ec;ibgcount++){
		      int ibgevent=int(rand.Uniform(0,1)*Ngood_ec);
		      for(int index=0;index<ecmodule;index++) edep_ec_preshower_FA[index] += hevent_ec_preshower_bg->GetBinContent(index+1,ibgevent);
		      for(int index=0;index<ecmodule;index++) edep_ec_shower_FA[index] += hevent_ec_shower_bg->GetBinContent(index+1,ibgevent);		      
		    }		  	    
		  }		
		}
		
		int ec_shower_E_tmp=0;
		int ecal_maxEdepID=-1;
		double edep_ec_shower_sum=0,edep_ec_preshower_sum=0;
		for(int index=0;index<ecmodule;index++){		  
		  edep_ec_shower_sum += edep_ec_shower_FA[index];
		  edep_ec_preshower_sum += edep_ec_preshower_FA[index];
		
		  //find EC module from largest shower energy deposit		
		  if(edep_ec_shower_FA[index]>ec_shower_E_tmp){
		    ec_shower_E_tmp=edep_ec_shower_FA[index];
		    ecal_maxEdepID=index;
		  }
// 		  if(edep_ec_shower_LA[index]>ec_shower_E_tmp){
// 		    ec_shower_E_tmp=edep_ec_shower_FA[index];
// 		    ecal_hitblockID=index;
// 		  }
		}
		
		int ecal_hitblockID=GetECALBlock_id(hit_x,hit_y); //get ec id from hit on virtual plane
		
// 		cout << "ecal_maxEdepID " <<  ecal_maxEdepID << " ecal_hitblockID " << ecal_hitblockID  << " " << hit_x << " " << hit_y << " " << hit_z << endl;
// 		cout << "ecal_maxEdepID " <<  ecal_maxEdepID << endl;
// 		
		//save background
		if(edep_ec_shower_sum>0 || edep_ec_preshower_sum>0){
		  event_ec_good++;
		  if(event_ec_good%100==0){ //merge every 100 events together to reduce output size
		    event_ec_good_merge++;
		    for(int index=0;index<ecmodule;index++) {
// 		      if(index==0) cout <<" event_ec_good " << event_ec_good << " " << event_ec_good_merge << endl;
		      hevent_ec_shower->Fill(index,event_ec_good_merge,edep_ec_shower_FA_tmp[index]);
		      hevent_ec_preshower->Fill(index,event_ec_good_merge,edep_ec_preshower_FA_tmp[index]);
		      edep_ec_shower_FA_tmp[index] = 0;
		      edep_ec_preshower_FA_tmp[index] = 0;
		    }
		  }else{
		    for(int index=0;index<ecmodule;index++) {
		      edep_ec_shower_FA_tmp[index] += edep_ec_shower_FA[index];
		      edep_ec_preshower_FA_tmp[index] += edep_ec_preshower_FA[index];
		    }
		  }
		}
		
		hedep_ec_shower_sum->Fill(edep_ec_shower_sum);
		hedep_ec_preshower_sum->Fill(edep_ec_preshower_sum);		
		hedep_ec_sum->Fill(edep_ec_shower_sum+edep_ec_preshower_sum,edep_ec_preshower_sum);
// 		hedep_ec_sum->Fill((edep_ec_shower_sum+edep_ec_preshower_sum)/hit_p/1e3,edep_ec_preshower_sum);		
		hnpe_lgc_edep_ec_sum->Fill(edep_ec_shower_sum+edep_ec_preshower_sum,edep_ec_preshower_sum,npe_total);
		
		const int Nneighbor=19; //hexagon has 6+1 or 18+1 layout
		Int_t cluster_edep_blockid[Nneighbor]={0};	
		int block_count = GetECALCluserAll(Nneighbor,ecal_hitblockID,cluster_edep_blockid); 		
// 		cout << "ecal_hitblockID " << ecal_hitblockID << " block_count " << block_count <<  endl;		
		
		double edep_ec_shower_cluster=0, edep_ec_preshower_cluster=0;
		for(Int_t j=0;j<Nneighbor;j++){
// 		  cout << " cluster_edep_blockid " << cluster_edep_blockid[j] << endl; 
		  if(cluster_edep_blockid[j] > 0){
		    edep_ec_shower_cluster += edep_ec_shower_FA[cluster_edep_blockid[j]];
		    edep_ec_preshower_cluster += edep_ec_preshower_FA[cluster_edep_blockid[j]];		  
		  }
		}

		hedep_ec_shower_cluster->Fill(edep_ec_shower_cluster);
		hedep_ec_preshower_cluster->Fill(edep_ec_preshower_cluster);		
		hedep_ec_cluster->Fill(edep_ec_shower_cluster+edep_ec_preshower_cluster,edep_ec_preshower_cluster);
// 		hedep_ec_cluster->Fill((edep_ec_shower_cluster+edep_ec_preshower_cluster)/hit_p/1e3,edep_ec_preshower_cluster);		
		hnpe_lgc_edep_ec_cluster->Fill(edep_ec_shower_cluster+edep_ec_preshower_cluster,edep_ec_preshower_cluster,npe_total);
		
		hedep_ec_clusterquality->Fill((edep_ec_shower_cluster+edep_ec_preshower_cluster)/(edep_ec_shower_sum+edep_ec_preshower_sum));
		
		//output textfile
		if (Is_textfile){		  
		  //output 3 sectors	
		  int Nsensor=ch_pmt/30;
		  int npe_tmp=0;
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec_left+index];
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec+index];
		  for(int index=0;index<Nsensor;index++) npe_tmp += hit_pmt[Nsensor*sec_right+index]; 
		  textfile_row_npe << npe_tmp << ",";
		  textfile_row_npe << hit_px << ","<< hit_py<< ","<<hit_pz<< ","<<hit_x<< ","<<hit_y<< ","<<hit_z<< ",";	
		  if (det=="lgcec") {
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_npe << edep_ec_preshower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_npe << 0 <<",";
		    }
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_npe << edep_ec_shower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_npe << 0 <<",";
		    }
		  }
		  if (Is_simsafe) textfile_row_npe << 1 << endl;
		  else 	textfile_row_npe << 0 << endl;
// 		  if (Is_simsafe) {
// 		    if (Is_decay)	textfile_row_npe << 2 << endl;
// 		    else 	textfile_row_npe << 1 << endl;
// 		  }
// 		  else {
// 		    if (Is_decay)	textfile_row_npe << -1 << endl;
// 		    else 	textfile_row_npe << 0 << endl;
// 		  }
		  
		  Nsensor=ch_pmt/30;
		  for(int index=0;index<Nsensor;index++) textfile_row_pmt << hit_pmt[Nsensor*sec_left+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_pmt << hit_pmt[Nsensor*sec+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_pmt << hit_pmt[Nsensor*sec_right+index] << ","; 
		  textfile_row_pmt << hit_px << ","<< hit_py<< ","<<hit_pz<< ","<<hit_x<< ","<<hit_y<< ","<<hit_z<< ",";
 		  if (det=="lgcec") {
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_pmt << edep_ec_preshower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_pmt << 0 <<",";
		    }
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_pmt << edep_ec_shower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_pmt << 0 <<",";
		    }
		  }
		  if (Is_simsafe) textfile_row_pmt << 1 << endl;
		  else 	textfile_row_pmt << 0 << endl;
// 		  if (Is_simsafe) {
// 		    if (Is_decay)	textfile_row_pmt << 2 << endl;
// 		    else 	textfile_row_pmt << 1 << endl;
// 		  }
// 		  else {
// 		    if (Is_decay)	textfile_row_pmt << -1 << endl;
// 		    else 	textfile_row_pmt << 0 << endl;
// 		  }

		  Nsensor=ch_quad/30;
		  for(int index=0;index<Nsensor;index++) textfile_row_quad << hit_quad[Nsensor*sec_left+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_quad << hit_quad[Nsensor*sec+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_quad << hit_quad[Nsensor*sec_right+index] << ","; 
		  textfile_row_quad << hit_px << ","<< hit_py<< ","<<hit_pz<< ","<<hit_x<< ","<<hit_y<< ","<<hit_z<< ",";	
		  if (det=="lgcec") {
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_quad << edep_ec_preshower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_quad << 0 <<",";
		    }
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_quad << edep_ec_shower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_quad << 0 <<",";
		    }
		  }
		  if (Is_simsafe) textfile_row_quad << 1 << endl;
		  else 	textfile_row_quad << 0 << endl;		  
// 		  if (Is_simsafe) {
// 		    if (Is_decay)	textfile_row_quad << 2 << endl;
// 		    else 	textfile_row_quad << 1 << endl;
// 		  }
// 		  else {
// 		    if (Is_decay)	textfile_row_quad << -1 << endl;
// 		    else 	textfile_row_quad << 0 << endl;
// 		  }

		  Nsensor=ch_pixel/30;
		  for(int index=0;index<Nsensor;index++) textfile_row_pixel << hit_pixel[Nsensor*sec_left+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_pixel << hit_pixel[Nsensor*sec+index] << ",";
		  for(int index=0;index<Nsensor;index++) textfile_row_pixel << hit_pixel[Nsensor*sec_right+index] << ","; 
		  textfile_row_pixel << hit_px << ","<< hit_py<< ","<<hit_pz<< ","<<hit_x<< ","<<hit_y<< ","<<hit_z<< ",";	
		  if (det=="lgcec") {
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_pixel << edep_ec_preshower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_pixel << 0 <<",";
		    }
		    for(Int_t j=0;j<Nneighbor;j++){
		      if(cluster_edep_blockid[j] > 0) textfile_row_pixel << edep_ec_shower_FA[cluster_edep_blockid[j+1]]<<",";	
		      else textfile_row_pixel << 0 <<",";
		    }
		  }
		  if (Is_simsafe) textfile_row_pixel << 1 << endl;
		  else 	textfile_row_pixel << 0 << endl;
// 		  if (Is_simsafe) {
// 		    if (Is_decay)	textfile_row_pixel << 2 << endl;
// 		    else 	textfile_row_pixel << 1 << endl;
// 		  }
// 		  else {
// 		    if (Is_decay)	textfile_row_pixel << -1 << endl;
// 		    else 	textfile_row_pixel << 0 << endl;
// 		  }
		  
		}
		
// 		if (2.5<p_gen && p_gen <3.0 && 8<theta_gen && theta_gen<9) hnpe[0]->Fill(npe_total);
// 		if (7.0<p_gen && p_gen <7.5 && 14<theta_gen && theta_gen<15) hnpe[1]->Fill(npe_total);
 
		    
// 		npe_total=0;		
// // 		if(ntrigsecs){
// 		if(true){
// 
// 		for(int sensor_id=0;sensor_id<ch;sensor_id++){
// 		      int sensor_sec=sensor_id/sensor;		  
// 		      int sensor_num=sensor_id%sensor;
// // 		      int sensor_x=sensor_num%sensor_trans,sensor_y=sensor_trans-1-sensor_num/sensor_trans;
// 		      int sensor_x=sensor_num%sensor_trans,sensor_y=sensor_num/sensor_trans;		      
// 		  
// 		      if(hit[sensor_id]>0){
// 			sensor_good++;			
// // 			cout << sensor_num << " " << sensor_x << " " << sensor_y << endl;
// 			
// 			hpe->Fill(sensor_num,hit[sensor_id]);	    
// 			hpe_2D->Fill(sensor_x,sensor_y,hit[sensor_id]);	
// 			
// 			npe_total += hit[sensor_id];
// 		      }
// 		      
// 		}		
// 		      
// 		if (npe_total>0){
// 		    event_cc_good++;
// 		    hnpe->Fill(npe_total);
// 	
// 		      if (theta_gen<=Max_Theta) {
// 			int bin_Theta=int(theta_gen/(Max_Theta/Nbin_Theta));
// 			int bin_Phi=int((phi_gen-(-180))/(360./Nbin_Phi));
// 		// 	 cout << "bin_Theta " << bin_Theta << "theta_gen " << theta_gen << endl;	
// 			hnpe_ThetaPhi[bin_Theta][bin_Phi]->Fill(npe_total);       
// 		      }
// 		      else cout << "theta_gen too large " << theta_gen << endl;
// 		}
// 			
// 		}
 
}
file->Close();

outputfile->Write();
outputfile->Flush();

cout <<" sensor_good " << sensor_good << endl;
cout <<" event_cc_good " << event_cc_good << endl;
cout <<" event_ec_good " << event_ec_good << endl;
cout <<" decay_no " << decay_no <<" decay_first " << decay_first << " decay_lost " << decay_lost << endl;

// TCanvas *c = new TCanvas("c","c",2000,1000);
// c->Divide(Nbin_Theta,Nbin_Phi);
// for(int j=0;j<Nbin_Phi;j++){
//   for(int i=0;i<Nbin_Theta;i++){    
//     c->cd(j*Nbin_Theta+i+1);
//     gPad->SetLogy();
//     hnpe_ThetaPhi[i][j]->Draw();
// //     cout << hnpe_ThetaPhi[i][j]->GetMean() << " "; 
//     havg_pe->SetBinContent(i+1,j+1,hnpe_ThetaPhi[i][j]->GetMean());
//     havg_pe->SetBinError(i+1,j+1,hnpe_ThetaPhi[i][j]->GetRMS());
//   }
// //   cout << endl;
// }
// 
// TCanvas *c_havg_pe = new TCanvas("havg_pe","havg_pe",1000,1000);
// havg_pe->SetMaximum(50);
// havg_pe->Draw("colz");
// c_havg_pe->SaveAs(Form("%s_avg_pe.png",the_filename));


if (Is_textfile) textfile_row_npe.close();
if (Is_textfile) textfile_row_pmt.close();
if (Is_textfile) textfile_row_quad.close();
if (Is_textfile) textfile_row_pixel.close();

TCanvas *c_check = new TCanvas("check","check",1600,1000);
c_check->Divide(3,2);
c_check->cd(1);
hnpe_no0->Draw();
c_check->cd(2);
hpe->Draw();  
c_check->cd(3);
hpe_2D->Draw("colz");
c_check->cd(4);
htime_photon->Draw();
c_check->cd(5);
hr_photon->Draw();
c_check->cd(6);
hwl_photon->Draw();
// c_check->SaveAs(Form("%s_count.png",the_filename));

TCanvas *c_npe_QE = new TCanvas("c_npe_QE","c_npe_QE",1600,1000);
hnpe_QE->Draw();

TLegend* leg= new TLegend(0.85, 0.75, 0.95, 0.9);
TCanvas *c_npe_secratio = new TCanvas("c_npe_secratio","c_npe_secratio",1600,1000);
hnpe_3secratio->SetLineColor(1);
// hnpe_3secratio->SetMarkerStyle(24);
hnpe_3secratio->Scale(1./hnpe_3secratio->GetEntries());
hnpe_3secratio->Draw("E");
leg->AddEntry(hnpe_3secratio,"hnpe_3secratio","l");  
hnpe_2secratioleft->SetLineColor(2);
// hnpe_2secratioleft->SetMarkerStyle(25);
hnpe_2secratioleft->Scale(1./hnpe_2secratioleft->GetEntries());
hnpe_2secratioleft->Draw("E same");
leg->AddEntry(hnpe_2secratioleft,"hnpe_2secratioleft","l");  
hnpe_2secratioright->SetLineColor(3);
// hnpe_2secratioright->SetMarkerStyle(26);
hnpe_2secratioright->Scale(1./hnpe_2secratioright->GetEntries());
hnpe_2secratioright->Draw("E same");
leg->AddEntry(hnpe_2secratioright,"hnpe_2secratioright","l");  
hnpe_1secratio->SetLineColor(4);
// hnpe_1secratio->SetMarkerStyle(27);
hnpe_1secratio->Scale(1./hnpe_1secratio->GetEntries());
hnpe_1secratio->Draw("E same");
leg->AddEntry(hnpe_1secratio,"hnpe_1secratio","l");  
leg->Draw();
hnpe_1secratio->SetMaximum(1);


// TCanvas *c_npe = new TCanvas("npe","npe",1600,1000);
// c_npe->Divide(2,2);
// c_npe->cd(1);
// hnpe[0]->Draw();
// c_npe->cd(2);
// hnpe[1]->Draw();

TCanvas *c_npe = new TCanvas("npe","npe",1600,1000);
c_npe->cd(1);
hnpe->Draw();

TCanvas *c_npe_2D = new TCanvas("npe_2D","npe_2D",1000,1000);
c_npe_2D->Divide(1,3);
c_npe_2D->cd(1);
hnpe_pmt_2D->Draw("colz");
c_npe_2D->cd(2);
hnpe_quad_2D->Draw("colz");
c_npe_2D->cd(3);
hnpe_pixel_2D->Draw("colz");

TCanvas *c_npe_2D_phi = new TCanvas("npe_2D_phi","npe_2D_phi",1900,1000);
c_npe_2D_phi->Divide(12,3);
for(int i=0;i<12;i++){
c_npe_2D_phi->cd(i+1);
hnpe_pmt_2D_phi[i]->Draw("colz");
c_npe_2D_phi->cd(12+i+1);
hnpe_quad_2D_phi[i]->Draw("colz");
c_npe_2D_phi->cd(24+i+1);
hnpe_pixel_2D_phi[i]->Draw("colz");  
}

TCanvas *c_ec = new TCanvas("ec","ec",1900,1000);
c_ec->Divide(3,3);
c_ec->cd(1);
hedep_ec_shower_sum->Draw();
c_ec->cd(2);
hedep_ec_preshower_sum->Draw();
c_ec->cd(3);
hedep_ec_sum->Draw("colz");
c_ec->cd(4);
hedep_ec_shower_cluster->Draw();
c_ec->cd(5);
hedep_ec_preshower_cluster->Draw();
c_ec->cd(6);
hedep_ec_cluster->Draw("colz");
c_ec->cd(7);
hedep_ec_clusterquality->Draw();

// TCanvas *c_hitxy = new TCanvas("hhitxy","hhitxy",1000,1000);
// hhitxy->Draw("colz");
// c_hitxy->SaveAs(Form("%s_hitpmt.png",the_filename));
// 
// TCanvas *c_hitxy_mirror = new TCanvas("hitxy_mirror","hitxy_mirror",1000,1000);
// hhitxy_mirror->Draw("colz");
// c_hitxy_mirror->SaveAs(Form("%s_hitmirror.png",the_filename));
// 
// TCanvas *c_hitxy_cone = new TCanvas("hitxy_cone","hitxy_cone",1000,1000);
// hhitxy_cone->Draw("colz");
// c_hitxy_cone->SaveAs(Form("%s_hitcone.png",the_filename));
 
// cout << hhitxy->GetSum() << endl;

}

