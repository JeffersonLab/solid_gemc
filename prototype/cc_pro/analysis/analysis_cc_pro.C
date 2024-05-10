#include <iostream> 
#include <fstream>
#include <cmath> 
#include "math.h" 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"

using namespace std;

#include "analysis_tree_solid_hgc.C"

// some numbers to be hard coded 
// make sure they are correct while using this script
//################################################################################################################################################## 
// const double filenum=50; //file numbers while running GEMC in order to be correct for normalization
const int loop_time=1;   //electron to be 1, pion to be many times to take advantage of statistics, pion has low efficiency on EC
const int add_norm=1; // additional normalization factor

//distance between two ec cluster to have coincidance trigger
// const double threshold_distance=0;
// const double threshold_distance=0.1;
const double threshold_distance=32.5; 

const int with_background_on_lgc=0;     //0: no background on lgc, 1: yes background on lgc

//trigger threshold
// lgc 
const double PEthresh_lgc=2; //lgc pe threshold for each pmt
const double PMTthresh_lgc=2; //lgc pmt threshold, at least 2pmts are fired in each sector

// hgc 
// const double PEthresh_hgc=2; //hgc pe threshold for each pmt
// const double PMTthresh_hgc=2; //hgc pmt threshold, at least 2pmts are fired in each sector
const double PEthresh_hgc=1; //hgc pe threshold for each pmt
const double PMTthresh_hgc=2; //hgc pmt threshold, at least 2pmts are fired in each sector

//spd 
// const double trigger_threshold_spd_FA=0.35;         //in MeV
const double trigger_threshold_spd_FA=0.5;         //in MeV
const double trigger_threshold_spd_LA=1.5;         //in MeV
//mrpc threshold
const double mrpc_block_threshold_FA=5;  //how many layers are required to be fired

//occupancy threshold
double occ_threshold_lgc=0,occ_threshold_hgc=0; //in N_p.e.
double occ_threshold_spd_FA=trigger_threshold_spd_FA/5.,occ_threshold_spd_LA=trigger_threshold_spd_LA/5.; //in MeV
double occ_threshold_ec_preshower=0.4,occ_threshold_ec_shower=6; //in MeV

// double occ_threshold_lgc=1,occ_threshold_hgc=1; //in N_p.e.
// // double occ_threshold_spd_FA=trigger_threshold_spd_FA/2.,occ_threshold_spd_LA=trigger_threshold_spd_LA/2.; //in MeV
// double occ_threshold_spd_FA=0.5,occ_threshold_spd_LA=3.; 
// double occ_threshold_ec_preshower=0.8,occ_threshold_ec_shower=12; //in MeV

//EC radius cut for physics result
double rout_cut_FA_phys=220;
double rin_cut_FA_phys=105;
double rout_cut_LA_phys=127;
double rin_cut_LA_phys=83; 
//EC radius cut for trigger
double rout_cut_FA_trigger=235;
double rin_cut_FA_trigger=105;
double rout_cut_LA_trigger=140; 
double rin_cut_LA_trigger=80;

bool Is_debug=false;

const double DEG=180./3.1415926;   //rad to degree

//#####################################################################################################################################################

int analysis_cc_pro(string inputfile_name,string runmode="trigger", bool Is_tellorig=false,string filetype="",bool Is_new=true){

// gStyle->SetOptStat(11111111);
// gStyle->SetOptStat("ioue");
gStyle->SetOptStat(0);

gStyle->SetPalette(57);

double rout_cut_FA=0,rin_cut_FA=0,rout_cut_LA=0,rin_cut_LA=0;
if (runmode=="phys"){
 cout << "runmode: phys" << endl;  
 rout_cut_FA=rout_cut_FA_phys;
 rin_cut_FA=rin_cut_FA_phys;
 rout_cut_LA=rout_cut_LA_phys;
 rin_cut_LA=rin_cut_LA_phys;   
}else if(runmode=="trigger"){
 cout << "runmode: trigger" << endl;
 rout_cut_FA=rout_cut_FA_trigger;
 rin_cut_FA=rin_cut_FA_trigger;
 rout_cut_LA=rout_cut_LA_trigger;
 rin_cut_LA=rin_cut_LA_trigger;   
}
else {cout << "need to know runmode: phys or trigger" << endl; return 0;}
  
bool Is_singlefile=false;
bool Is_pi0=false;
if(Is_tellorig){
if(filetype.find("single",0) != string::npos) {
  Is_singlefile=true;
  cout << "this is a single file" << endl;  
}
else if(filetype.find("sidis",0) != string::npos) {
  Is_singlefile=false;
  cout << "this is a sidis file" << endl;      
}
else {cout << "unknown file type, choose either single or sidis" << endl;return 0;}

if (inputfile_name.find("pi0",0) != string::npos) {
  Is_pi0=true;
  cout << "this is a pi0 file" << endl;  
}
else {cout << "this is NOT a pi0 file" << endl;}
}

string filemode;
double event_actual=1;
if (inputfile_name.find("BeamOnTargetEM",0) != string::npos) {
  filemode="BeamOnTargetEM";
  cout << "this is a BeamOnTargetEM file" << endl;  
  
  event_actual=atof(inputfile_name.substr(inputfile_name.find("BeamOnTargetEM",0)+15,inputfile_name.find("_")).c_str());
  cout << "event_actual " << event_actual <<  endl;  
}
else if (inputfile_name.find("BeamOnTarget",0) != string::npos) {
  filemode="BeamOnTarget";
  cout << "this is a BeamOnTarget file" << endl;  
  
  event_actual=atof(inputfile_name.substr(inputfile_name.find("BeamOnTarget",0)+13,inputfile_name.find("_")).c_str());
  cout << "event_actual " << event_actual <<  endl;  
}
else if (inputfile_name.find("even",0) != string::npos) {
  filemode="even";
  cout << "this is a evenly distributed file" << endl;  
}
else {
  filemode="rate";  
  cout << "this is rate dependent file" << endl;  
}

if (event_actual<1) {cout << "wrong event_actual" << endl; return 0;}

double filenum=1;
if (inputfile_name.find("_filenum",0) != string::npos) {
  filenum=atof(inputfile_name.substr(inputfile_name.find("_filenum")+8,inputfile_name.find("_")).c_str());
    cout << "filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
}
else {
  if (filemode=="rate"){
    cout << "this file is rate dependent, but has no filenum, something is wrong" << endl;  
    return 0;
  }
  else{
    cout << "this file has no filenum, please check if you need filenum for addtional normalization" << endl;      
  }
}

// bool Is_He3=false,Is_C=false,Is_NOtarget=false;
// if(inputfile_name.find("_He3_",0) != string::npos) {
//   Is_He3=true;
//   cout << "He3 setup" << endl;  
// }
// else if(inputfile_name.find("_C_",0) != string::npos) {
//   Is_C=true;
//   cout << "C setup" << endl;  
// }
// else if(inputfile_name.find("_NOtarget",0) != string::npos) {
//   Is_NOtarget=true;
//   cout << "NOtarget setup" << endl;  
// }
// else {
//     cout << "Not He3 or C or NOtarget" << endl;    
//     return 0;
// }

//Cherenkov sensor for 30 sectors
const int ch_lgc=270;
// const int ch_hgc=480;    	//pmt readout
const int ch_hgc=1920;	//quad readout
// const int ch_hgc=30720;		//pixel readout

const int sensor_hgc = ch_hgc/30;
const int sensor_trans_hgc = sqrt(sensor_hgc);		

Float_t npe_lgc[ch_lgc],npe_hgc[ch_hgc];
Float_t npe_lgc_total=0,npe_lgc_total_trigged=0;
Float_t npe_hgc_total=0,npe_hgc_total_trigged=0;
		  
// char mlpfile_lgc_name[200];
// sprintf(mlpfile_lgc_name, "%s_mlp_lgc.root",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());
// TFile *mlpfile_lgc=new TFile(mlpfile_lgc_name, "recreate");
// TTree *mlptree_lgc = new TTree("lgc", "npe");
// for(int i=0;i<ch_lgc;i++) mlptree_lgc->Branch(Form("npe_lgc%i",i), &npe_lgc[i], Form("npe_lgc%i/F",i));
// mlptree_lgc->Branch("npe_lgc_total", &npe_lgc_total, "npe_lgc_total/F");
// 
// char mlpfile_hgc_name[200];
// sprintf(mlpfile_hgc_name, "%s_mlp_hgc.root",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());
// TFile *mlpfile_hgc=new TFile(mlpfile_hgc_name, "recreate");
// TTree *mlptree_hgc = new TTree("hgc", "npe");
// for(int i=0;i<ch_hgc;i++) mlptree_hgc->Branch(Form("npe_hgc%i",i), &npe_hgc[i], Form("npe_hgc%i/F",i));
// mlptree_hgc->Branch("npe_hgc_total", &npe_hgc_total, "npe_hgc_total/F");

TFile *file=new TFile(inputfile_name.c_str());

// 	TString background_inputfile_name="parametrized_lgc.root";      //h_pe is here	
// 	TFile *background_file=new TFile(background_inputfile_name);
// 	TH1F *h_pe=(TH1F*)background_file->Get("h_pe");

std::size_t found = inputfile_name.rfind("cache");
if (found!=std::string::npos)  inputfile_name.replace(found,5,"work");

char the_filename[500];
sprintf(the_filename, "%s",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());

char outputfile_name[200];
sprintf(outputfile_name, "%s_output.root",the_filename);
TFile *outputfile=new TFile(outputfile_name, "recreate");

// prepare for outputs
// define histograms, output txt files etc...

// 	TH1F *hpe_hgc=new TH1F("pe_hgc","pe;sensor;",sensor_hgc,0,sensor_hgc);
// 	TH1F *hhit_hgc=new TH1F("hit_hgc","hit;sensor;N_{p.e.} rate(kHz)",sensor_hgc,0,sensor_hgc);
// 	TH1F *hocc_hgc=new TH1F("occ_hgc","occupancy;sensor;rate(kHz)",sensor_hgc,0,sensor_hgc);
// 	
// 	TH2F *hpe_hgc_2D=new TH2F("pe_hgc_2D","pe;sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);	
// 	TH2F *hhit_hgc_2D=new TH2F("hit_hgc_2D","hit,N_{p.e.} rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
// 	TH2F *hocc_hgc_2D=new TH2F("occ_hgc_2D","occupancy,rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);

	TH1F *hpe_hgc[4],*hhit_hgc[4],*hocc_hgc[4];
	TH2F *hpe_hgc_2D[4],*hhit_hgc_2D[4],*hocc_hgc_2D[4];	
	for(int i=0;i<4;i++){
	  char hstname[200]; 	  
	  sprintf(hstname,"pe_hgc_%i",i);	  
	  hpe_hgc[i]=new TH1F(hstname,"N_{p.e.};sensor;",sensor_hgc,0,sensor_hgc);
	  sprintf(hstname,"hit_hgc_%i",i);	  
	  hhit_hgc[i]=new TH1F(hstname,"N_{p.e.}*rate(kHz);sensor;N_{p.e.} rate(kHz)",sensor_hgc,0,sensor_hgc);
	  sprintf(hstname,"occ_hgc_%i",i);	  
	  hocc_hgc[i]=new TH1F(hstname,"rate(kHz);sensor; rate(kHz)",sensor_hgc,0,sensor_hgc);
	
	  sprintf(hstname,"pe_hgc_2D_%i",i);	  	  
	  hpe_hgc_2D[i]=new TH2F(hstname,"N_{p.e.};sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);	
	  sprintf(hstname,"hit_hgc_2D_%i",i);	  	  
	  hhit_hgc_2D[i]=new TH2F(hstname,"N_{p.e.}*rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
	  sprintf(hstname,"occ_hgc_2D_%i",i);	  	  
	  hocc_hgc_2D[i]=new TH2F(hstname,"rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
	}
	
// 	TH1F *htime_photon[2];
// 	htime_photon[0]=new TH1F("time_photon_0","time",400,0,1000);
// 	htime_photon[1]=new TH1F("time_photon_1","time",400,0,1000);	
	TH1F *htime_photon=new TH1F("time_photon","time_photon;t (ns)",100,0,100);

	TH1F *hmotherP=new TH1F("motherP","mother particle rate (kHz);P (GeV);",5500,0,11);

	TH1F *hmother_procid=new TH1F("mother_procid","mother_procid",150,-0.5,149.5);

	TH1F *hmother_tid=new TH1F("mother_tid","mother_tid",100,0,100);	
	TH1F *hmother_pid=new TH1F("mother_pid","mother_pid",200,-100,100);
	TH1F *hmother_vz=new TH1F("mother_vz","mother particle rate (kHz);vz (cm);",1300,-100,1200);	
	TH2F *hmother_vzvx=new TH2F("mother_vzvx","mother particle rate (kHz);vz (cm); vx (cm)",1300,-100,1200,1000,-20,20);		

	TH1F *hgrandma_tid=new TH1F("grandma_tid","grandma_tid",100,0,100);	
	TH1F *hgrandma_pid=new TH1F("grandma_pid","grandma_pid",200,-100,100);
	TH1F *hgrandma_vz=new TH1F("grandma_vz","grandma particle rate (kHz);vz (cm);",1300,-100,1200);	
	TH2F *hgrandma_vzvx=new TH2F("grandma_vzvx","grandma particle rate (kHz);vz (cm); vx (cm)",1300,-100,1200,1000,-20,20);	
	
	TH1F *hnpe_hgc_count[2],*hnpe_hgc[2];
	TH1F *hnpe_sensor_hgc_count[2][sensor_hgc],*hnpe_sensor_hgc[2][sensor_hgc];		
	TH2F *hnpe_nsensor_hgc_count[2],*hnpe_nsensor_hgc[2];	
	for(int j=0;j<2;j++){	
	  hnpe_hgc_count[j]=new TH1F(Form("npe_hgc_count_%i",j),";Npe;count",100,0,100);		
	  hnpe_hgc[j]=new TH1F(Form("npe_hgc_%i",j),";Npe;rate(kHz)",100,0,100);
	  hnpe_nsensor_hgc_count[j]=new TH2F(Form("npe_nsensor_hgc_count_%i",j),";Npe;Nsensor",100,0,100,sensor_hgc+1,-0.5,sensor_hgc+0.5);
	  hnpe_nsensor_hgc[j]=new TH2F(Form("npe_nsensor_hgc_%i",j),";Npe;Nsensor",100,0,100,sensor_hgc,0.5,sensor_hgc+0.5);
	  for(int i=0;i<sensor_hgc;i++){	
	    hnpe_sensor_hgc_count[j][i]=new TH1F(Form("npe_sensor_hgc_count_%i_%i",j,i),";Npe;count",30,0,30);		
	    hnpe_sensor_hgc[j][i]=new TH1F(Form("npe_sensor_hgc_%i_%i",j,i),";Npe;rate(kHz)",30,0,30);
	  }	  
	}
	
	TH1F *hnpe_hgc_cut[sensor_hgc][2];
	for(int i=0;i<sensor_hgc;i++){
	  for(int j=0;j<2;j++){	
	  char hstname[200]; 
	  sprintf(hstname,"npe_hgc_cut_%i_%i",i,j);
	  hnpe_hgc_cut[i][j]=new TH1F(hstname,";Npe;rate(kHz)",200,0,200);
	  }	  
	}	

	TH1F *hnpe_hgc_mom[10],*hnpe_hgc_mom_theta[10][10];
	for(int i=0;i<10;i++){
	  char hstname[200]; 
	  sprintf(hstname,"npe_hgc_mom_%i",i);
	  hnpe_hgc_mom[i]=new TH1F(hstname,";Npe;count",200,0,200);
	  for(int j=0;j<10;j++){
	    sprintf(hstname,"npe_hgc_mom_theta_%i_%i",i,j);
	    hnpe_hgc_mom_theta[i][j]=new TH1F(hstname,";Npe;count",200,0,200);
	  }	  
	}
	
	TH2F *hvertex_rz=new TH2F("hvertex_rz","hvertex_rz",1800,-400,500,600,0,300);
	
	TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","generated events;vertex #theta (deg);vertex P (GeV)",100,0,50,110,0,11);  
	TH2F *hgen_ThetaPhi=new TH2F("gen_ThetaPhi","generated events;vertex #theta (deg);vertex #phi (deg)",100,0,50,360,-180,180);     
	TH2F *hgen_PhiP=new TH2F("gen_PhiP","generated events;vertex #phi (deg);vertex P (GeV)",360,-180,180,110,0,11);	
	
	TH3F *hgen_ThetaPhiP=new TH3F("gen_ThetaPhiP","gen_ThetaPhiP",50,0,50,180,-180,180,55,0,11);   
	
	const int n=4;
	char *detname[n]={"sc_front","sc_back","ec"};	
	TH1F *hhit_vz[n];	
	TH2F *hhit_xy[n],*hhit_PhiR[n];
	TH1F *hhit_E[n],*hhit_E_mip[n],*hhit_E_photonele[n],*hhit_E_ele[n],*hhit_Edep[n];		  	
	for(int i=0;i<n;i++){
	  char hstname[100];
	  sprintf(hstname,"hit_vz_%i",i);
	  hhit_vz[i]=new TH1F(hstname,detname[i],1600,-100,1500);  
	  sprintf(hstname,"hit_xy_%i",i);
	  hhit_xy[i]=new TH2F(hstname,detname[i],400,-20,20,400,-20,20);        
	  sprintf(hstname,"hit_PhiR_%i",i);
	  hhit_PhiR[i]=new TH2F(hstname,detname[i],360,-180,180,300,0,300);
	  sprintf(hstname,"hit_E_%i",i);
	  hhit_E[i]=new TH1F(hstname,detname[i],1100,0,11);
	  sprintf(hstname,"hit_E_mip_%i",i);
	  hhit_E_mip[i]=new TH1F(hstname,detname[i],1100,0,11);
	  sprintf(hstname,"hit_E_photonele_%i",i);
	  hhit_E_photonele[i]=new TH1F(hstname,detname[i],1100,0,11);
	  sprintf(hstname,"hit_E_ele_%i",i);
	  hhit_E_ele[i]=new TH1F(hstname,detname[i],1100,0,11);	  
	  sprintf(hstname,"hit_Edep_%i",i);
	  hhit_Edep[i]=new TH1F(hstname,detname[i],1000,0,0.1);	  
	}
	
// 	TH1F *hhit_Edepsc1=new TH1F("hit_Edepsc1","SC1;Edep (GeV);rate(kHz)",1000,0,0.1);
// 	TH1F *hhit_Edepsc2=new TH1F("hit_Edepsc2","SC2;Edep (GeV);rate(kHz)",1000,0,0.1);
// 	TH1F *hhit_Eec=new TH1F("hit_Eec","EC;(GeV);rate(kHz)",1100,0,11);
// 	TH1F *hhit_Eec_photonele=new TH1F("hit_Eec_photonele","EC (photon,e-/e+);E(GeV);rate(kHz)",1100,0,11);
// 	TH1F *hhit_Eec_ele=new TH1F("hit_Eec_ele","EC (e-/e+);E(GeV);rate(kHz)",1100,0,11);			
// 	TH1F *hhit_Pec_pion=new TH1F("hit_Pec_pion","EC;P(GeV);rate(kHz)",1100,0,11);
// 	TH1F *hhit_Pec_proton=new TH1F("hit_Pec_proton","EC;P(GeV);rate(kHz)",1100,0,11);	

	TH1F *hprocid_sc1=new TH1F("procid_sc1","SC1;ProcID;rate(kHz)",150,-0.5,149.5);
	TH1F *hprocid_sc2=new TH1F("procid_sc2","SC2;ProcID;rate(kHz)",150,-0.5,149.5);
	TH1F *hprocid_ec=new TH1F("procid_ec","EC;ProcID;rate(kHz)",150,-0.5,149.5);

	TH1F *hhit_Edepsc1=new TH1F("hit_Edepsc1","SC1;log10(Edep)(GeV);rate(kHz)",100,-6,-1);
	TH1F *hhit_Edepsc2=new TH1F("hit_Edepsc2","SC2;log10(Edep)(GeV);rate(kHz)",100,-6,-1);
	TH1F *hhit_Eec=new TH1F("hit_Eec","EC (all);log10(E)(GeV);rate(kHz)",70,-6,1);
	TH1F *hhit_Eec_photonele=new TH1F("hit_Eec_photonele","EC (photon,e-/e+);log10(E)(GeV);rate(kHz)",70,-6,1);
	TH1F *hhit_Eec_ele=new TH1F("hit_Eec_ele","EC (e-/e+);log10(E)(GeV);rate(kHz)",70,-6,1);			
	TH1F *hhit_Pec_pion=new TH1F("hit_Pec_pion","EC;log10(P)(GeV);rate(kHz)",70,-6,1);
	TH1F *hhit_Pec_proton=new TH1F("hit_Pec_proton","EC;log10(P)(GeV);rate(kHz)",70,-6,1);	
	
	double hit_Edepsc1_count[3],hit_Edepsc2_count[3];
	double hit_Eec_count[3],hit_Eec_photonele_count[3],hit_Eec_ele_count[3];
	
	double hit_Pec_pion_count;
	double hit_Pec_proton_count;
	
	TH1F *hhit_mom_hgc_e=new TH1F("hit_mom_hgc_e","hit_mom_hgc_e",100,0,0.05);		
	//-------------------------
	//   get trees in the real data file
	//-------------------------
	
	//---header tree
	TTree *tree_header = (TTree*) file->Get("header");
	vector <int> *evn=0,*evn_type=0,*runNo=0;
	vector <double> *beamPol=0;
	vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
	tree_header->SetBranchAddress("evn",&evn);      // number 
	tree_header->SetBranchAddress("evn_type",&evn_type);  // evn_type==-1 for simulated events
	tree_header->SetBranchAddress("beamPol",&beamPol);   //beam polarization
	tree_header->SetBranchAddress("runNo",&runNo);  // run number	  	
	if (!Is_new){	
	tree_header->SetBranchAddress("var1",&var1);     // W+ rate
	tree_header->SetBranchAddress("var2",&var2);     // W- rate
	tree_header->SetBranchAddress("var3",&var3);     // target pol
	tree_header->SetBranchAddress("var4",&var4);     //x
	tree_header->SetBranchAddress("var5",&var5);     //y
	tree_header->SetBranchAddress("var6",&var6);     //w
	tree_header->SetBranchAddress("var7",&var7);     //Q2
	tree_header->SetBranchAddress("var8",&var8);     //rate, Hz, should check the input file of the simulation
	}
	
	TTree *tree_userHeader;
	vector <double> *userVar001=0,*userVar002=0,*userVar003=0,*userVar004=0,*userVar005=0,*userVar006=0,*userVar007=0,*userVar008=0,*userVar009=0,*userVar010=0;
	if (Is_new){
	tree_userHeader = (TTree*) file->Get("userHeader");
	tree_userHeader->SetBranchAddress("userVar001",&userVar001);
	tree_userHeader->SetBranchAddress("userVar002",&userVar002);
	tree_userHeader->SetBranchAddress("userVar003",&userVar003);
	tree_userHeader->SetBranchAddress("userVar004",&userVar004);
	tree_userHeader->SetBranchAddress("userVar005",&userVar005);
	tree_userHeader->SetBranchAddress("userVar006",&userVar006);
	tree_userHeader->SetBranchAddress("userVar007",&userVar007);
	tree_userHeader->SetBranchAddress("userVar008",&userVar008);
	tree_userHeader->SetBranchAddress("userVar009",&userVar009);
	tree_userHeader->SetBranchAddress("userVar010",&userVar010);
	}
	
	//---generated tree
	//particle generated with certain momentum at certain vertex
	TTree *tree_generated = (TTree*) file->Get("generated");
	vector <int> *gen_pid=0;
	vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	tree_generated->SetBranchAddress("pid",&gen_pid);   //particle ID 
	tree_generated->SetBranchAddress("px",&gen_px);     //momentum of the generated particle at target
	tree_generated->SetBranchAddress("py",&gen_py);
	tree_generated->SetBranchAddress("pz",&gen_pz);
	tree_generated->SetBranchAddress("vx",&gen_vx);    //vertex of the generated particle at target
	tree_generated->SetBranchAddress("vy",&gen_vy);
	tree_generated->SetBranchAddress("vz",&gen_vz);

	//--- flux
	//the real deal output from the GEMC simulation
	TTree *tree_flux = (TTree*) file->Get("flux");
	vector<int> *flux_id=0,*flux_hitn=0;
	vector<int> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0,*flux_procID=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
	tree_flux->SetBranchAddress("hitn",&flux_hitn);     // hit number
	tree_flux->SetBranchAddress("id",&flux_id);         //hitting detector ID
	tree_flux->SetBranchAddress("pid",&flux_pid);       //pid
	tree_flux->SetBranchAddress("mpid",&flux_mpid);     // mother pid
	tree_flux->SetBranchAddress("tid",&flux_tid);       // track id
	tree_flux->SetBranchAddress("mtid",&flux_mtid);     // mother track id
	tree_flux->SetBranchAddress("otid",&flux_otid);     // original track id
	tree_flux->SetBranchAddress("trackE",&flux_trackE);  // track energy of 1st step,  track here is G4 track
	tree_flux->SetBranchAddress("totEdep",&flux_totEdep); // totEdep in all steps, track here is G4 track
	tree_flux->SetBranchAddress("avg_x",&flux_avg_x);     //average x, weighted by energy deposition in each step
	tree_flux->SetBranchAddress("avg_y",&flux_avg_y);     //average y
	tree_flux->SetBranchAddress("avg_z",&flux_avg_z);     //average z
	tree_flux->SetBranchAddress("avg_lx",&flux_avg_lx);   // local average x 
	tree_flux->SetBranchAddress("avg_ly",&flux_avg_ly);   // local average y
	tree_flux->SetBranchAddress("avg_lz",&flux_avg_lz);   // local average z
	tree_flux->SetBranchAddress("px",&flux_px);          // px of 1st step
	tree_flux->SetBranchAddress("py",&flux_py);          // py of 1st step
	tree_flux->SetBranchAddress("pz",&flux_pz);          // pz of 1st step
	tree_flux->SetBranchAddress("vx",&flux_vx);          // x coordinate of 1st step
	tree_flux->SetBranchAddress("vy",&flux_vy);          // y coordinate of 1st step
	tree_flux->SetBranchAddress("vz",&flux_vz);          // z coordinate of 1st step
	tree_flux->SetBranchAddress("mvx",&flux_mvx);        // mother
	tree_flux->SetBranchAddress("mvy",&flux_mvy);
	tree_flux->SetBranchAddress("mvz",&flux_mvz);
	tree_flux->SetBranchAddress("avg_t",&flux_avg_t);     //average time stamp
	tree_flux->SetBranchAddress("procID",&flux_procID);     // process id

	//---hgc
	//information recorded by hgc
	TTree* tree_solid_hgc= (TTree*) file->Get("solid_hgc");
	setup_tree_solid_hgc(tree_solid_hgc);	
	
	TRandom3 rand;
	rand.SetSeed(0);

	int sensor_good=0;
	int event_good=0,event_trig_good=0;
	
// 	long int N_events = (long int)tree_header->GetEntries();
	long int N_events = (long int)tree_generated->GetEntries();	

	cout << "total number of events : " << N_events << endl;	

	//----------------------------
	//      loop trees
	//---------------------------
	for(int loop_id=1;loop_id<=loop_time;loop_id++){
		cout<<"loop.....  "<<loop_id<<endl;
	
	for(long int i=0;i<N_events;i++){	  		
// 	for(long int i=1;i<N_events-1;i++){	  				
// 	for(long int i=0;i<N_events/100;i++){	  
// 	for(long int i=N_events/2;i<N_events;i++){	  		
// 	for(long int i=520;i<521;i++){  //pip event
// 	for(long int i=5289;i<5290;i++){	  // background event			  
// 			cout<<"" << i<<endl;
		if (i%1000==0) cout<<i<<"\r";
// 		if (i%1000==0) cout<<i<<"\n";
		
		//---
		//---header tree
		//---
		double rate=0;
		if (Is_new){
			tree_userHeader->GetEntry(i);
			rate=userVar006->at(0); // new eDIS and eAll generator
	// 		cout << "rate  " << rate   << endl;
		}
		else {
			tree_header->GetEntry(i);
			rate=var8->at(0); // old eDIS and eAll generator
		}

		if (filemode=="BeamOnTargetEM" || filemode=="BeamOnTarget") {
		  rate=1e-6/1.6e-19/event_actual/loop_time*add_norm;  //1ua		  
// 		  if(Is_He3) rate=1e-6/1.6e-19/event_actual/loop_time*add_norm;
// 		  else if(Is_C) rate=1e-6/1.6e-19/event_actual/loop_time*add_norm;		  
// 		  else if(Is_NOtarget) rate=1e-6/1.6e-19/event_actual/loop_time*add_norm;		  		  
// 		  else {
// 		    cout << "Not He3 or C or NOtarget setup" << endl;    
// 		    return 0;
// 		  }
		}
		else if (filemode=="even") rate=1;
		else if (filemode=="rate") rate=rate/filenum/loop_time*add_norm;
		else {
		    cout << "Not right filemode" << endl;    
		    return 0; 
		}
		
/*		double x=var4->at(0);	
		double y=var5->at(0);
		double W=var6->at(0);		
		double Q2=var7->at(0);	*/	
		//cout<<"header tree: "<<rate<<endl;	

		//---
		//---generated tree
		//---
		tree_generated->GetEntry(i);
		int n_gen=gen_pid->size();
		//cout<<"generated : "<<n_gen<<endl;
		int pid_gen=0;
		double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;      
	      //       cout << "gen_pid->size() " << gen_pid->size() << endl;        
		for (int j=0;j<gen_pid->size();j++) {
// 	            cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
		    pid_gen=gen_pid->at(j);
		    px_gen=gen_px->at(j);
		    py_gen=gen_py->at(j);
		    pz_gen=gen_pz->at(j);
		    vx_gen=gen_vx->at(j);
		    vy_gen=gen_vy->at(j);
		    vz_gen=gen_vz->at(j);
		    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
		    theta_gen=acos(pz_gen/p_gen)*DEG;
		    phi_gen=atan2(py_gen,px_gen)*DEG;
// 	            cout << "p_gen " << p_gen << endl; 
		    
		}		
// 		TVector3 *vec_p_gen(px_gen,py_gen,pz_gen),*vec_v_gen(vx_gen,vy_gen,vz_gen);
		
		// 		if (phi_gen<-1 || phi_gen>1) continue;
// 			cout<<"" << i<<endl;
		
// 		if (vz_gen<-3550 || -3450<vz_gen) continue;

		hgen_ThetaP->Fill(theta_gen,p_gen/1e3,rate);
		hgen_ThetaPhi->Fill(theta_gen,phi_gen,rate);                  		
		hgen_PhiP->Fill(phi_gen,p_gen/1e3,rate);                  				
		hgen_ThetaPhiP->Fill(theta_gen,phi_gen,p_gen/1e3,rate);                  			
		///////////////////////////////////////////////////////////////////////////////////////
		//       do trigger
		////////////////////////////////////////////////////////////////////////////////////////
		//---	
		//---flux tree
		//---
		tree_flux->GetEntry(i);		  
		  
		  int sec_hgc=0;		

		bool Is_trig=false;					  
		
		//check on EC and other by flux
// 		 cout << "flux_hitn  " << flux_hitn->size() << endl;
		double Eec=0,Eec_photonele=0,Eec_ele=0,Edepsc1=0,Edepsc2=0;
		
		for (Int_t j=0;j<flux_hitn->size();j++) {
// 	          cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  

		double hit_vr=sqrt(pow(flux_vx->at(j),2)+pow(flux_vy->at(j),2))/1e1; //mm to cm
		double hit_vy=flux_vy->at(j)/1e1,hit_vx=flux_vx->at(j)/1e1,hit_vz=flux_vz->at(j)/1e1;           //mm to cm		  
		double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))/1e1; //mm to cm
		double hit_y=flux_avg_y->at(j)/1e1,hit_x=flux_avg_x->at(j)/1e1,hit_z=flux_avg_z->at(j)/1e1;           //mm to cm		
		double hit_phi=atan2(hit_y,hit_x)*DEG;       //rad to  deg
		double hit_p=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j))/1e3;  //MeV to GeV
		
// 		TVector3 *vec_hit(hit_x,hit_y,hit_z), *vec_v(hit_vx,hit_vy,hit_vz),*vec_p(flux_px->at(j),flux_py->at(j),flux_pz->at(j));
// 		TVector3 *vec_path=vec_hit-vec_v_gen;
		  
		  int hit_id=-1;
		  if(flux_id->at(j)==1) hit_id=0; // SC front,SC1
		  else if(flux_id->at(j)==2) hit_id=1; // SC back, SC2
		  else if(flux_id->at(j)==3) hit_id=2; // EC
		  else if(flux_id->at(j)==10) hit_id=3; // Cherenkov front window
		  else cout << "wrong flux_id" << flux_id->at(j) << endl;
// 		  if (hit_id==-1) {/cout << flux_id->at(j) << " " << flux_avg_z->at(j) << endl;
		  
		  hhit_xy[hit_id]->Fill(flux_avg_lx->at(j)/1e1,flux_avg_ly->at(j)/1e1,rate/1e3);
		  hhit_vz[hit_id]->Fill(hit_vz,rate/1e3);

		  double E=flux_trackE->at(j)/1e3,Edep=flux_totEdep->at(j)/1e3;		  
		  hhit_Edep[hit_id]->Fill(Edep,rate/1e3);		  
		  hhit_E[hit_id]->Fill(E,rate/1e3);
		  if (abs(flux_pid->at(j)) == 211 || flux_pid->at(j)==2212|| flux_pid->at(j)==13) hhit_E_mip[hit_id]->Fill(E,rate/1e3);
		  if (abs(flux_pid->at(j)) == 11 || flux_pid->at(j)==22) hhit_E_photonele[hit_id]->Fill(E,rate/1e3);
		  if (abs(flux_pid->at(j)) == 11) hhit_E_ele[hit_id]->Fill(E,rate/1e3);  
		  
// 		  if (hit_id==1 && 1e-3<Edep && Edep <2e-3) {		  
// 		  if (hit_id==2 && 0.93<E && E <1) {
// 		  cout << "flux " << " !!! " << " " << flux_pid->at(j) << " " <<  flux_totEdep->at(j) << endl;
// 		  cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
// 		  }
		  
		  if (hit_id==0) {
		    Edepsc1 += Edep;
		    hprocid_sc1->Fill(flux_procID->at(j),rate/1e3);
		  }
		  if (hit_id==1) {
		    Edepsc2 += Edep;
		    hprocid_sc2->Fill(flux_procID->at(j),rate/1e3);
		  }		  
		  if (hit_id==2) {
// 		    cout <<  "x,y " << fabs(flux_avg_lx->at(j)) << " " << fabs(flux_avg_ly->at(j)) << endl;
		    //EC central module (made of 4 small module)
		    if (fabs(flux_avg_lx->at(j))/10.<4.25*2.54/2. && fabs(flux_avg_ly->at(j))/10.<4.25*2.54/2.){
		    //EC central c9_4 small module https://logbooks.jlab.org/entry/3853160		    
// 		     if (fabs(flux_avg_lx->at(j)/1e1-(-4.25*2.54/2/2))<4.25*2.54/2/2 && fabs(flux_avg_ly->at(j)/1e1-(-4.25*2.54/2/2))<4.25*2.54/2/2){
		    
// 		      if (flux_pid->at(j) == 22) 
			Eec += E;
		    }		    
		    if (abs(flux_pid->at(j)) == 11 || flux_pid->at(j)==22) Eec_photonele +=E;
		    if (abs(flux_pid->at(j)) == 11) Eec_ele +=E;	
		    
		    if (abs(flux_pid->at(j)) == 211 ) {hhit_Pec_pion->Fill(log10(hit_p),rate/1e3); hit_Pec_pion_count += rate/1e3;}	    
		    if (abs(flux_pid->at(j)) == 2212) {hhit_Pec_proton->Fill(log10(hit_p),rate/1e3); hit_Pec_proton_count += rate/1e3;}	    
		    
		    hprocid_ec->Fill(flux_procID->at(j),rate/1e3);
		  }		  
	  
		  //fill every hits
// 		  if (hit_id==0) {
// 		    Edepsc1 += Edep;
// 		    hprocid_sc1->Fill(flux_procID->at(j),rate/1e3);
// 		    
// 		    hhit_Edepsc1->Fill(log10(Edep),rate/1e3);
// 		  }
// 		  if (hit_id==1) {
// 		    Edepsc2 += Edep;
// 		    hprocid_sc2->Fill(flux_procID->at(j),rate/1e3);
// 
// 		    hhit_Edepsc2->Fill(log10(Edep),rate/1e3);		    
// 		  }		  
// 		  if (hit_id==2) {  //new change 2021/02/21
// 		    Eec += E;
// 		    hhit_Eec->Fill(log10(E),rate/1e3);		    
// 		    
// 		    if (abs(flux_pid->at(j)) == 11 || flux_pid->at(j)==22) {
// 		      Eec_photonele +=E;
// 		      hhit_Eec_photonele->Fill(log10(E),rate/1e3);
// 		    }
// 		    if (abs(flux_pid->at(j)) == 11) {
// 		      Eec_ele +=E;	
// 		      hhit_Eec_ele->Fill(log10(E),rate/1e3);		      
// 		    }
// 		    
// 		    if (abs(flux_pid->at(j)) == 211 ) {hhit_Pec_pion->Fill(log10(hit_p),rate/1e3); hit_Pec_pion_count += rate/1e3;}	    
// 		    if (abs(flux_pid->at(j)) == 2212) {hhit_Pec_proton->Fill(log10(hit_p),rate/1e3); hit_Pec_proton_count += rate/1e3;}	    
// 		    
// 		    hprocid_ec->Fill(flux_procID->at(j),rate/1e3);
// 		  }		  
		    
		}	// end of flux		

		//fill every events
// 		  hhit_Edepsc1->Fill(Edepsc1,rate/1e3);
// 		  hhit_Edepsc2->Fill(Edepsc2,rate/1e3);
// 		  hhit_Eec->Fill(Eec,rate/1e3);
// 		  hhit_Eec_photonele->Fill(Eec_photonele,rate/1e3);		  
// 		  hhit_Eec_ele->Fill(Eec_ele,rate/1e3);	

		  hhit_Edepsc1->Fill(log10(Edepsc1),rate/1e3);
		  hhit_Edepsc2->Fill(log10(Edepsc2),rate/1e3);
		  hhit_Eec->Fill(log10(Eec),rate/1e3);
		  hhit_Eec_photonele->Fill(log10(Eec_photonele),rate/1e3);		  
		  hhit_Eec_ele->Fill(log10(Eec_ele),rate/1e3);	
  
// 		  if (Edepsc2>0.7e-3 && Eec>0.1) Is_trig=true;
		  if (Eec>0.01) Is_trig=true;	  
		  
		//-----------------------	
		//--- hgc 
		//-----------------------
		tree_solid_hgc->GetEntry(i);

		double hit_hgc[ch_hgc]={0};
		int trigger_hgc[30]={0};
		int ntrigsecs_hgc=0;
		int photon_mtid=0;
		process_tree_solid_hgc(tree_solid_hgc,hit_hgc,trigger_hgc,ntrigsecs_hgc,PMTthresh_hgc,PEthresh_hgc,ch_hgc,htime_photon,photon_mtid);

		for (Int_t j=0;j<flux_hitn->size();j++) {
		  if(flux_id->at(j)==1) {
// 		      cout << "photon_mtid " << photon_mtid << endl;		      		    
		    if(flux_tid->at(j)==photon_mtid) {
		      hmotherP->Fill(sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j))/1e3,rate/1e3);
		      hmother_tid->Fill(flux_tid->at(j),rate/1e3);		      		      
		      hmother_pid->Fill(flux_pid->at(j),rate/1e3);
		      hmother_vz->Fill(flux_vz->at(j)/10.,rate/1e3);
		      hmother_vzvx->Fill(flux_vz->at(j)/10.,flux_vx->at(j)/10.,rate/1e3);		      		      

		      hmother_procid->Fill(flux_procID->at(j),rate/1e3);	  		      
		      
		      hgrandma_tid->Fill(flux_mtid->at(j),rate/1e3);		      		      
		      hgrandma_pid->Fill(flux_mpid->at(j),rate/1e3);	
		      hgrandma_vz->Fill(flux_mvz->at(j)/10.,rate/1e3);
		      hgrandma_vzvx->Fill(flux_mvz->at(j)/10.,flux_mvx->at(j)/10.,rate/1e3);		      		      
		    }
		    
// 		    else cout << "photon_mtid " << photon_mtid << endl;
		  }
		}
				
// 		if(ntrigsecs_hgc){
		if(true){
  
		for(int sensor_id=0;sensor_id<ch_hgc;sensor_id++){
		      int sensor_sec=sensor_id/sensor_hgc;		  
		      int sensor_num=sensor_id%sensor_hgc;
// 		      int sensor_x=sensor_num%sensor_trans_hgc,sensor_y=sensor_trans_hgc-1-sensor_num/sensor_trans_hgc;
		      int sensor_x=sensor_num%sensor_trans_hgc,sensor_y=sensor_num/sensor_trans_hgc;		      
		  
		      if(hit_hgc[sensor_id]>0){
			sensor_good++;			
// 			cout << sensor_num << " " << sensor_x << " " << sensor_y << endl;
			
			hpe_hgc[0]->Fill(sensor_num,hit_hgc[sensor_id]);	    
			hpe_hgc_2D[0]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]);			
			hhit_hgc[0]->Fill(sensor_num,hit_hgc[sensor_id]*rate/1e3);	    
			hhit_hgc_2D[0]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]*rate/1e3);
			hocc_hgc[0]->Fill(sensor_num,rate/1e3);
			hocc_hgc_2D[0]->Fill(sensor_x,sensor_y,rate/1e3);
			if (Edepsc2>0.7e-3){
			  hpe_hgc[1]->Fill(sensor_num,hit_hgc[sensor_id]);	    
			  hpe_hgc_2D[1]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]);			
			  hhit_hgc[1]->Fill(sensor_num,hit_hgc[sensor_id]*rate/1e3);	    
			  hhit_hgc_2D[1]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]*rate/1e3);
			  hocc_hgc[1]->Fill(sensor_num,rate/1e3);
			  hocc_hgc_2D[1]->Fill(sensor_x,sensor_y,rate/1e3);			  
			}
// 			if (Eec>0.3){
			if (Eec>0.1){		  
			  hpe_hgc[2]->Fill(sensor_num,hit_hgc[sensor_id]);	    
			  hpe_hgc_2D[2]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]);			
			  hhit_hgc[2]->Fill(sensor_num,hit_hgc[sensor_id]*rate/1e3);	    
			  hhit_hgc_2D[2]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]*rate/1e3);
			  hocc_hgc[2]->Fill(sensor_num,rate/1e3);
			  hocc_hgc_2D[2]->Fill(sensor_x,sensor_y,rate/1e3);			  
			}	
// 			if (Edepsc2>0.7e-3 && Eec>0.3){
			if (Edepsc2>0.7e-3 && Eec>0.1){
			  hpe_hgc[3]->Fill(sensor_num,hit_hgc[sensor_id]);	    
			  hpe_hgc_2D[3]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]);			
			  hhit_hgc[3]->Fill(sensor_num,hit_hgc[sensor_id]*rate/1e3);	    
			  hhit_hgc_2D[3]->Fill(sensor_x,sensor_y,hit_hgc[sensor_id]*rate/1e3);
			  hocc_hgc[3]->Fill(sensor_num,rate/1e3);
			  hocc_hgc_2D[3]->Fill(sensor_x,sensor_y,rate/1e3);			  
			}			    
		      }	

		}		
		
		} //pass hgc trigger in offline		
		
// 		if (1<=sec_hgc && sec_hgc<=30){
// 		    int this_sec=sec_hgc;
// 		    int last_sec=-1,next_sec=-1;
// 		    if (this_sec==1) {last_sec=30;next_sec=2;}
// 		    else if (this_sec==30) {last_sec=29;next_sec=1;}
// 		    else {last_sec=this_sec-1;next_sec=this_sec+1;}		    		  
// 		  double npe_hgc=0;
// 		  for(int sensor_id=0;sensor_id<16;sensor_id++) npe_hgc += hit_hgc[16*(sec_hgc-1)+sensor_id];
// 		  for(int sensor_id=0;sensor_id<16;sensor_id++) npe_hgc += hit_hgc[16*(last_sec-1)+sensor_id];
// 		  for(int sensor_id=0;sensor_id<16;sensor_id++) npe_hgc += hit_hgc[16*(next_sec-1)+sensor_id];

		  npe_hgc_total=0;
		  npe_hgc_total_trigged=0;		  
		  for(int sensor_id=0;sensor_id<ch_hgc;sensor_id++) {	
		    if (hit_hgc[sensor_id]>0){
// 		      cout << "" << i << " hgc sensor_id " << sensor_id << " npe " << hit_hgc[sensor_id] << endl;
		      npe_hgc_total += hit_hgc[sensor_id];  
		      
		      if (Is_trig) npe_hgc_total_trigged += hit_hgc[sensor_id];  		      
		      
// 		      if(sensor_id>=16) cout << "why " << sensor_id << endl;
// 		      npe_hgc[sensor_id] += hit_hgc[sensor_id];
// 		      npe_hgc[sensor_id]=1;
		    }
		  }
		  
		  if (npe_hgc_total>0) {		  
// 		  if (npe_hgc_total>0 && Is_trig) {
		    event_good++;
		    
		    hnpe_hgc_count[0]->Fill(npe_hgc_total);		    
		    hnpe_hgc[0]->Fill(npe_hgc_total,rate/1e3);
		    for(int sensor_id=0;sensor_id<sensor_hgc;sensor_id++) {	
		      if (hit_hgc[sensor_id]>0) {
			hnpe_sensor_hgc_count[0][sensor_id]->Fill(hit_hgc[sensor_id]);		    
			hnpe_sensor_hgc[0][sensor_id]->Fill(hit_hgc[sensor_id],rate/1e3);		    
		      }
		    }
		    
		    hnpe_nsensor_hgc_count[0]->Fill(npe_hgc_total,trigger_hgc[0]);		    
		    hnpe_nsensor_hgc[0]->Fill(npe_hgc_total,trigger_hgc[0],rate/1e3);		    
		    if (1<=trigger_hgc[0] && trigger_hgc[0]<=sensor_hgc){
		      hnpe_hgc_cut[trigger_hgc[0]-1][0]->Fill(npe_hgc_total,rate/1e3);
		    }
		    else cout << "wrong trigger_hgc[0] " << trigger_hgc[0] << endl;
		  }
		  if (npe_hgc_total_trigged>0) {
		    event_trig_good++;
		    
		    hnpe_hgc_count[1]->Fill(npe_hgc_total_trigged);		    
		    hnpe_hgc[1]->Fill(npe_hgc_total_trigged,rate/1e3);	
		    for(int sensor_id=0;sensor_id<sensor_hgc;sensor_id++) {	
		      if (hit_hgc[sensor_id]>0) {		      
			hnpe_sensor_hgc_count[1][sensor_id]->Fill(hit_hgc[sensor_id]);		    
			hnpe_sensor_hgc[1][sensor_id]->Fill(hit_hgc[sensor_id],rate/1e3);		    
		      }
		    }		    
		    
		    hnpe_nsensor_hgc_count[1]->Fill(npe_hgc_total_trigged,trigger_hgc[0]);		    
		    hnpe_nsensor_hgc[1]->Fill(npe_hgc_total_trigged,trigger_hgc[0],rate/1e3);		    		    
		    if (1<=trigger_hgc[0] && trigger_hgc[0]<=sensor_hgc){
		      hnpe_hgc_cut[trigger_hgc[0]-1][1]->Fill(npe_hgc_total_trigged,rate/1e3);
		    }
		    else cout << "wrong trigger_hgc[0] " << trigger_hgc[0] << endl;
		  }
		  
// 		  if (npe_hgc_total>0) cout<<"" << i<< " npe_hgc_total " << npe_hgc_total << endl;	
		  
// 		  for(int sensor_id=0;sensor_id<ch_hgc;sensor_id++) cout << npe_hgc[sensor_id] << "\t";
// 		  if (npe_hgc_total>0) mlptree_hgc->Fill();
		  
// 		if (-365< vz_gen/10. && vz_gen/10. < -335){
		if (7<=theta_gen && theta_gen <=15){ //has to cut to print(-1)=0 !!!!!
		if (2500<=p_gen && theta_gen <=7500){	//has to cut to print(-1)=0	  
		  
		  int index_mom=int((p_gen/1e3-2.5)/0.5);  //int(-1)=0 Warning !!!!!!!!!!!!!!!!!!
		  if (0<=index_mom && index_mom<10){
// 		    hnpe_hgc->Fill(npe_hgc_total*0.5);
		    hnpe_hgc_mom[index_mom]->Fill(npe_hgc_total*0.5);
		    int index_theta=int((theta_gen-7));    //int(-1)=0 Warning !!!!!!!!!!!!!!!!!!
		    if (0<=index_theta && index_theta<10){
		      hnpe_hgc_mom_theta[index_mom][index_theta]->Fill(npe_hgc_total*0.5);	
		    }
		  }
		}
		}
// 		}
		
		
// 		  if (npe_hgc==0){
// 		    cout<<"" << i<<endl;
// 		    cout << "npe_hgc " << npe_hgc << endl;		    
// 		    cout << "sec_hgc " << sec_hgc << endl;
// 		for (Int_t j=0;j<flux_hitn->size();j++) {
// 	          if (flux_id->at(j)==2210000 && (flux_tid->at(j)==1 || abs(flux_pid->at(j))==221 || abs(flux_pid->at(j))==13))   cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
// 		}		  
// 		  }						
		 		
		

	} //end loop
	
	} //end loop time

cout <<" sensor_good " << sensor_good << endl;
cout <<" event_good " << event_good << endl;
cout <<" event_trig_good " << event_trig_good << endl;

//do outputs

outputfile->Write();	
outputfile->Flush();

cout << "hocc_hgc[0]/sensor_hgc " << hocc_hgc[0]->Integral()/sensor_hgc << endl;

TCanvas *c0 = new TCanvas("c0", "c0",1900,900);
c0->Divide(3,2);
c0->cd(1);
hhit_hgc_2D[0]->Draw("colz");
c0->cd(2);
hocc_hgc_2D[0]->Draw("colz");
c0->cd(3);
hpe_hgc_2D[0]->Draw("colz");
c0->cd(4);
hhit_hgc[0]->Draw("HIST");
c0->cd(5);
hocc_hgc[0]->Draw("HIST");
c0->cd(6);
hpe_hgc[0]->Draw("HIST");
c0->SaveAs("c0.png");

TCanvas *c1 = new TCanvas("c1", "c1",1000,900);
hnpe_hgc_count[0]->SetAxisRange(0,100);  
hnpe_hgc_count[0]->Draw("HIST");
c1->SaveAs("c1.png");


TCanvas *c2 = new TCanvas("c2", "c2",1000,900);
// hpe_hgc_2D[0]->Scale(1/1e3);
hpe_hgc_2D[0]->Draw("colz");
// c2->SaveAs("c2.png");

TCanvas *c3 = new TCanvas("c3", "c3",1000,900);
hocc_hgc[0]->Draw();
c3->SaveAs("c3.png");

TCanvas *c_npe_nsensor_hgc = new TCanvas("npe_nsensor_hgc", "npe_nsensor_hgc",1800,1000);
c_npe_nsensor_hgc->Divide(2,2);
c_npe_nsensor_hgc->cd(1);
hnpe_nsensor_hgc_count[0]->Draw("colz");
c_npe_nsensor_hgc->cd(2);
hnpe_nsensor_hgc_count[1]->Draw("colz");
c_npe_nsensor_hgc->cd(3);
hnpe_nsensor_hgc[0]->Draw("colz");
c_npe_nsensor_hgc->cd(4);
hnpe_nsensor_hgc[1]->Draw("colz");

TCanvas *c_npe_nsensor_hgc_before = new TCanvas("npe_nsensor_hgc_before", "npe_nsensor_hgc_before",1800,1000);
hnpe_nsensor_hgc[1]->SetTitle(";;");
hnpe_nsensor_hgc[1]->Draw("colz");

TCanvas *c_npe_nsensor_hgc_after = new TCanvas("npe_nsensor_hgc_after", "npe_nsensor_hgc_after",1800,1000);
TH2F *hnpe_nsensor_hgc_after=(TH2F*) hnpe_nsensor_hgc[1]->Clone();
for(int x=0;x<100;x++) {	
hnpe_nsensor_hgc_after->SetBinContent(x,1,0);
hnpe_nsensor_hgc_after->SetBinContent(x,2,0);
// hnpe_nsensor_hgc_after->SetBinContent(x,3,0);  //for quad
}
hnpe_nsensor_hgc_after->SetTitle(";;");
hnpe_nsensor_hgc_after->Draw("colz");

TCanvas *c_npe = new TCanvas("npe", "npe",1000,900);
// hpe_hgc_2D[0]->Scale(1/1e3);
hnpe_hgc[0]->Draw();
c_npe->SaveAs("c_npe.png");
// cout << "npe_hgc_0 " << " " << hnpe_hgc[0]->Integral() << endl;

TCanvas *c_npe_sensor = new TCanvas("npe_sensor", "npe_sensor",1000,900);
c_npe_sensor->Divide(sensor_trans_hgc,sensor_trans_hgc);
for(int sensor_id=0;sensor_id<sensor_hgc;sensor_id++) {	
  c_npe_sensor->cd(sensor_id+1);
  hnpe_sensor_hgc[0][sensor_id]->Draw();
//   cout << "npe_sensor_hgc_0 " << sensor_id  << " "  << hnpe_sensor_hgc[0][sensor_id]->Integral() << endl;  
}
c_npe_sensor->SaveAs("c_npe_sensor.png");


TCanvas *c_6 = new TCanvas("c6", "c6",1000,600);
c_6->Divide(3,1);
c_6->cd(1);
htime_photon->Draw("HIST");
c_6->cd(2);
gPad->SetLogx();
hmotherP->Draw("HIST");
c_6->cd(3);
// gPad->SetLogx();
hmother_procid->Draw("HIST");

TCanvas *c_mother = new TCanvas("mother", "mother",1800,1000);
c_mother->Divide(4,2);
c_mother->cd(1);
gPad->SetLogy();
hmother_tid->Draw("HIST");
c_mother->cd(2);
gPad->SetLogy();
hmother_pid->Draw("HIST");
c_mother->cd(3);
gPad->SetLogy();
hmother_vz->Draw("HIST");
c_mother->cd(4);
gPad->SetLogz();
hmother_vzvx->Draw("colz");
c_mother->cd(5);
gPad->SetLogy();
hgrandma_tid->Draw("HIST");
c_mother->cd(6);
gPad->SetLogy();
hgrandma_pid->Draw("HIST");
c_mother->cd(7);
gPad->SetLogy();
hgrandma_vz->Draw("HIST");
c_mother->cd(8);
gPad->SetLogz();
hgrandma_vzvx->Draw("colz");

TCanvas *c_mother_final= new TCanvas("mother_final", "mother_final",1900,1000);
c_mother_final->Divide(3,1);
c_mother_final->cd(1);
gPad->SetLogy();
hmother_vz->Draw("HIST");
c_mother_final->cd(2);
gPad->SetLogz();
hmother_vzvx->Draw("colz");
c_mother_final->cd(3);
gPad->SetLogx();
hmotherP->Draw("HIST");
c_mother_final->SaveAs("mother_final.png");
c_mother_final->SaveAs("mother_final.pdf");


// TCanvas *c_sc_count = new TCanvas("sc_count", "sc_count",1800,1000);
// c_sc_count->Divide(2,1);
// c_sc_count->cd(1);
// hhit_Edepsc1_count->Draw();
// cout << hhit_Edepsc1_count->GetBinContent(1) << " " << hhit_Edepsc1_count->GetBinContent(2) << " " << hhit_Edepsc1_count->GetBinContent(3) << endl;
// c_sc_count->cd(2);
// hhit_Edepsc2_count->Draw();
// cout << hhit_Edepsc2_count->GetBinContent(1) << " " << hhit_Edepsc2_count->GetBinContent(2) << " " << hhit_Edepsc2_count->GetBinContent(3) << endl;
// 
// TCanvas *c_ec_count = new TCanvas("ec_count", "ec_count",1800,1000);
// c_ec_count->Divide(3,1);
// c_ec_count->cd(1);
// hhit_Eec_count->Draw();
// cout << hhit_Eec_count->GetBinContent(1) << " " << hhit_Eec_count->GetBinContent(2) << " " << hhit_Eec_count->GetBinContent(3) << endl;
// c_ec_count->cd(2);
// hhit_Eec_photonele_count->Draw();
// cout << hhit_Eec_photonele_count->GetBinContent(1) << " " << hhit_Eec_photonele_count->GetBinContent(2) << " " << hhit_Eec_photonele_count->GetBinContent(3) << endl;
// c_ec_count->cd(3);
// hhit_Eec_ele_count->Draw();
// cout << hhit_Eec_ele_count->GetBinContent(1) << " " << hhit_Eec_ele_count->GetBinContent(2) << " " << hhit_Eec_ele_count->GetBinContent(3) << endl;

TCanvas *c_ecsc = new TCanvas("ecsc", "ecsc",1800,1000);
c_ecsc->Divide(3,2);
c_ecsc->cd(1);
gPad->SetLogy();
hhit_Edepsc1->SetLineColor(kBlack);
hhit_Edepsc1->Draw();
hhit_Edepsc2->SetLineColor(kRed);
hhit_Edepsc2->Draw("same");
c_ecsc->cd(2);
gPad->SetLogy();
hhit_Eec->SetLineColor(kBlack);
hhit_Eec->Draw();
hhit_Eec_photonele->SetLineColor(kRed);
hhit_Eec_photonele->Draw("same");
hhit_Eec_ele->SetLineColor(kBlue);
hhit_Eec_ele->Draw("same");
c_ecsc->cd(3);
gPad->SetLogy();
hhit_Pec_pion->SetLineColor(kBlack);
hhit_Pec_pion->Draw();
hhit_Pec_proton->SetLineColor(kRed);
hhit_Pec_proton->Draw("same");
c_ecsc->cd(4);
hprocid_sc1->Draw("HIST");
c_ecsc->cd(5);
hprocid_sc1->Draw("HIST");
c_ecsc->cd(6);
hprocid_ec->Draw("HIST");

cout << hhit_Edepsc1->Integral(hhit_Edepsc1->FindBin(-5),hhit_Edepsc1->GetNbinsX())/11 << " " << hhit_Edepsc1->Integral(hhit_Edepsc1->FindBin(-4),hhit_Edepsc1->GetNbinsX())/11 << " " <<hhit_Edepsc1->Integral(hhit_Edepsc1->FindBin(-3),hhit_Edepsc1->GetNbinsX())/11 << endl;

cout << hhit_Edepsc2->Integral(hhit_Edepsc2->FindBin(-5),hhit_Edepsc2->GetNbinsX())/11 << " " << hhit_Edepsc2->Integral(hhit_Edepsc2->FindBin(-4),hhit_Edepsc2->GetNbinsX())/11 << " " <<hhit_Edepsc2->Integral(hhit_Edepsc2->FindBin(-3),hhit_Edepsc2->GetNbinsX())/11 << endl;

cout << hhit_Eec->Integral(hhit_Eec->FindBin(-2),hhit_Eec->GetNbinsX())/9 << " " << hhit_Eec->Integral(hhit_Eec->FindBin(-1),hhit_Eec->GetNbinsX())/9 << " " <<hhit_Eec->Integral(hhit_Eec->FindBin(log10(0.2)),hhit_Eec->GetNbinsX())/9 << endl;

cout << hhit_Eec_photonele->Integral(hhit_Eec_photonele->FindBin(-2),hhit_Eec_photonele->GetNbinsX())/9 << " " << hhit_Eec_photonele->Integral(hhit_Eec_photonele->FindBin(-1),hhit_Eec_photonele->GetNbinsX())/9 << " " <<hhit_Eec_photonele->Integral(hhit_Eec_photonele->FindBin(log10(0.2)),hhit_Eec_photonele->GetNbinsX())/9 << endl;

cout << hhit_Eec_ele->Integral(hhit_Eec_ele->FindBin(-2),hhit_Eec_ele->GetNbinsX())/9 << " " << hhit_Eec_ele->Integral(hhit_Eec_ele->FindBin(-1),hhit_Eec_ele->GetNbinsX())/9 << " " <<hhit_Eec_ele->Integral(hhit_Eec_ele->FindBin(log10(0.2)),hhit_Eec_ele->GetNbinsX())/9 << endl;

cout << hhit_Pec_pion->Integral()/9 << endl;

cout << hhit_Pec_proton->Integral()/9 << endl;

// cout << hit_Edepsc1_count[0]/11 << " " << hit_Edepsc1_count[1]/11 << " " << hit_Edepsc1_count[2]/11 << endl;
// 
// cout << hit_Edepsc2_count[0]/11 << " " << hit_Edepsc2_count[1]/11 << " " << hit_Edepsc2_count[2]/11 << endl;
// 
// cout << hit_Eec_count[0]/9 << " " << hit_Eec_count[1]/9 << " " << hit_Eec_count[2]/9 << endl;
// 
// cout << hit_Eec_photonele_count[0]/9 << " " << hit_Eec_photonele_count[1]/9 << " " << hit_Eec_photonele_count[2]/9 << endl;
// 
// cout << hit_Eec_ele_count[0]/9 << " " << hit_Eec_ele_count[1]/9 << " " << hit_Eec_ele_count[2]/9 << endl;
// 
// cout << hit_Pec_pion_count/9 << endl;
// 
// cout << hit_Pec_proton_count/9 << endl;



// mlpfile_lgc->Write();	
// mlpfile_lgc->Flush();
// mlpfile_hgc->Write();	
// mlpfile_hgc->Flush();

}
