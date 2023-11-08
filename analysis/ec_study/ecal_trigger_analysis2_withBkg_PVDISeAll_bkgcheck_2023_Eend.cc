/*
Author : Rakitha Beminiwattha, rakithab@jlab.org
Wed Jun 29 12:58:25 EDT 2016
file name : ecal_energy_res.cc

This is the first script to access GEMC based ECAL. Ultimate goal is to get energy resolution for ECAL

to run,
./ecal_energy_res 

It is problematic to access Trees with std::vector data in it using TChain. The issues fixed by following the thread discussion at, https://root.cern.ch/phpBB3/viewtopic.php?t=4939

*/

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <new>
#include <cstdlib>
#include <math.h>

#include <TRandom.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TSystem.h>

#include <TH2F.h>
#include <TTree.h>
#include <TF1.h>
#include <TProfile.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h> 
#include <TDatime.h>
#include <TStopwatch.h>
#include <stdexcept>
#include <time.h>
#include <cstdio>
#include <map>
#include <cassert>

#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>

#include <TCanvas.h>
#include <TLine.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFrame.h>
#include <TObjArray.h>
#include <TVector2.h>
#include <TVirtualFitter.h>

#include "ecal.h"

#define __IO_MAXHIT 10000
#define __OCTANTS 30
Bool_t trigger_state_1[__OCTANTS]={kFALSE}; //octant level 6+1 blocks trigger state for the event
Bool_t trigger_state_2[__OCTANTS]={kFALSE}; //octant level 2+1 blocks trigger state for the event
Bool_t trigger_state_PS[__OCTANTS]={kFALSE};//if single event is above MIP then trigger state for the window will be true
Bool_t trigger_state_3={kFALSE};
Bool_t window_trigger_state_1=kFALSE;//trigger window level 6+1 blocks trigger state for the event
Bool_t window_trigger_state_2=kFALSE;//trigger window level 2+1 blocks trigger state for the event
Bool_t bDisableBkg=kFALSE;//disable background signal loading
Bool_t bEarlyBeakBkg = kFALSE;//early break background loop
Bool_t kSaveRootFile=kTRUE;
const Double_t DEG=180./3.1415926;
Double_t fTotalRate[3][5];
const Double_t TriggerWindow =30;// ns
const Double_t DelayedHitTimeLimit = 50;// ns//or 30 ns was set before Tue Jun 23 11:02:34 EDT 2015
int event_bkg_N=0;
TFile * rootfile;
//header bank
vector <int> *evn=0,*evn_type=0;
vector <double> *beamPol=0;
vector <int> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
//Generated input bank
vector<Int_t> *fGen_pid=0;
vector<Double_t> *fGen_Px=0;
vector<Double_t> *fGen_Py=0;
vector<Double_t> *fGen_Pz=0;
vector<Double_t> *fGen_vx=0;
vector<Double_t> *fGen_vy=0;
vector<Double_t> *fGen_vz=0;
Int_t fTimeWindow;
//FLUX bank
vector<Int_t> *fFluxHit_n=0;
vector<Int_t> *fFluxHit_id=0;
vector<Int_t> *fFluxHit_tid=0;
vector<Int_t> *fFluxHit_mtid=0;
vector<Int_t> *fFluxHit_otid=0;
vector<Int_t> *fFluxHit_pid=0;
vector<Int_t> *fFluxHit_mpid=0;
vector<Double_t> *fFluxHit_trackE=0;
vector<Double_t> *fFluxHit_totEdep=0;
vector<Double_t> *fFluxHit_Px=0;
vector<Double_t> *fFluxHit_Py=0;
vector<Double_t> *fFluxHit_Pz=0;
vector<Double_t> *fFluxHit_T=0;
vector<Double_t> *fFluxHit_Avg_x=0;
vector<Double_t> *fFluxHit_Avg_y=0;
vector<Double_t> *fFluxHit_Avg_z=0;
vector<Double_t> *fFluxHit_Avg_lx=0;
vector<Double_t> *fFluxHit_Avg_ly=0;
vector<Double_t> *fFluxHit_Avg_lz=0;
vector<Double_t> *fFluxHit_vx=0;
vector<Double_t> *fFluxHit_vy=0;
vector<Double_t> *fFluxHit_vz=0;
vector<Double_t> *fFluxHit_mvx=0;
vector<Double_t> *fFluxHit_mvy=0;
vector<Double_t> *fFluxHit_mvz=0;
vector<Double_t> *fFluxHit_Avg_t=0;
Double_t edep_PS_sum[__OCTANTS]={0};
TBranch *bFluxHit_id = 0;
TBranch *bFluxHit_n = 0;
TBranch *bFluxHit_pid = 0;
TBranch *bFluxHit_mpid = 0;
TBranch *bFluxHit_tid = 0;
TBranch *bFluxHit_mtid = 0;
TBranch *bFluxHit_otid = 0;
TBranch *bFluxHit_trackE = 0;
TBranch *bFluxHit_totEdep = 0;
TBranch *bFluxHit_Avg_x = 0;
TBranch *bFluxHit_Avg_y = 0;
TBranch *bFluxHit_Avg_z = 0;
TBranch *bFluxHit_Avg_lx = 0;
TBranch *bFluxHit_Avg_ly = 0;
TBranch *bFluxHit_Avg_lz = 0;
TBranch *bFluxHit_Px = 0;
TBranch *bFluxHit_Py = 0;
TBranch *bFluxHit_Pz = 0;
TBranch *bFluxHit_vx = 0;
TBranch *bFluxHit_vy = 0;
TBranch *bFluxHit_vz = 0;
TBranch *bFluxHit_mvx = 0;
TBranch *bFluxHit_mvy = 0;
TBranch *bFluxHit_mvz = 0;
TBranch *bFluxHit_t = 0;

Int_t rndTimeWindow;//to save the random number of the trigger window
//EC bank
vector<int> *fECHit_n=0;
vector<int> *fECHit_id=0;
vector<Int_t> *fECHit_tid=0;
vector<Int_t> *fECHit_mtid=0;
vector<Int_t> *fECHit_otid=0;
vector<Int_t> *fECHit_pid=0;
vector<Int_t> *fECHit_mpid=0;
vector<Double_t> *fECHit_trackE=0;
vector<Double_t> *fECHit_totEdep=0;
vector<Double_t> *fECHit_totEdepB=0;
vector<Double_t> *fECHit_totEend=0;
vector<Double_t> *fECPSHit_totEdepB=0;
vector<Double_t> *fECPSHit_totEend=0;
vector<Double_t> *fECHit_Px=0;
vector<Double_t> *fECHit_Py=0;
vector<Double_t> *fECHit_Pz=0;
vector<Double_t> *fECHit_Avg_x=0;
vector<Double_t> *fECHit_Avg_y=0;
vector<Double_t> *fECHit_Avg_z=0;
vector<Double_t> *fECHit_Avg_lx=0;
vector<Double_t> *fECHit_Avg_ly=0;
vector<Double_t> *fECHit_Avg_lz=0;
vector<Double_t> *fECHit_vx=0;
vector<Double_t> *fECHit_vy=0;
vector<Double_t> *fECHit_vz=0;
vector<Double_t> *fECHit_mvx=0;
vector<Double_t> *fECHit_mvy=0;
vector<Double_t> *fECHit_mvz=0;
vector<Double_t> *fECHit_Avg_t=0;
vector<int> *fECPSHit_n=0;
vector<int> *fECPSHit_id=0;
vector<Int_t> *fECPSHit_tid=0;
vector<Int_t> *fECPSHit_mtid=0;
vector<Int_t> *fECPSHit_otid=0;
vector<Int_t> *fECPSHit_pid=0;
vector<Int_t> *fECPSHit_mpid=0;
vector<Double_t> *fECPSHit_trackE=0;
vector<Double_t> *fECPSHit_totEdep=0;
vector<Double_t> *fECPSHit_Px=0;
vector<Double_t> *fECPSHit_Py=0;
vector<Double_t> *fECPSHit_Pz=0;
vector<Double_t> *fECPSHit_Avg_x=0;
vector<Double_t> *fECPSHit_Avg_y=0;
vector<Double_t> *fECPSHit_Avg_z=0;
vector<Double_t> *fECPSHit_Avg_lx=0;
vector<Double_t> *fECPSHit_Avg_ly=0;
vector<Double_t> *fECPSHit_Avg_lz=0;
vector<Double_t> *fECPSHit_vx=0;
vector<Double_t> *fECPSHit_vy=0;
vector<Double_t> *fECPSHit_vz=0;
vector<Double_t> *fECPSHit_mvx=0;
vector<Double_t> *fECPSHit_mvy=0;
vector<Double_t> *fECPSHit_mvz=0;
vector<Double_t> *fECPSHit_Avg_t=0;
vector<double>  *userVar001;
vector<double>  *userVar002;
vector<double>  *userVar003;
vector<double>  *userVar004;
vector<double>  *userVar005;
vector<double>  *userVar006;
vector<double>  *userVar007;
vector<double>  *userVar008;
vector<double>  *userVar009;
vector<double>  *userVar010;
TBranch        *b_userVar001;   //!
TBranch        *b_userVar002;   //!
TBranch        *b_userVar003;   //!
TBranch        *b_userVar004;   //!
TBranch        *b_userVar005;   //!
TBranch        *b_userVar006;   //!
TBranch        *b_userVar007;   //!
TBranch        *b_userVar008;   //!
TBranch        *b_userVar009;   //!
TBranch        *b_userVar010;   //!

TBranch *bECHit_id = 0;
TBranch *bECHit_n = 0;
TBranch *bECHit_pid = 0;
TBranch *bECHit_mpid = 0;
TBranch *bECHit_tid = 0;
TBranch *bECHit_mtid = 0;
TBranch *bECHit_otid = 0;
TBranch *bECHit_trackE = 0;
TBranch *bECHit_totEdep = 0;
TBranch *bECHit_Avg_x = 0;
TBranch *bECHit_Avg_y = 0;
TBranch *bECHit_Avg_z = 0;
TBranch *bECHit_Avg_lx = 0;
TBranch *bECHit_Avg_ly = 0;
TBranch *bECHit_Avg_lz = 0;
TBranch *bECHit_Px = 0;
TBranch *bECHit_Py = 0;
TBranch *bECHit_Pz = 0;
TBranch *bECHit_vx = 0;
TBranch *bECHit_vy = 0;
TBranch *bECHit_vz = 0;
TBranch *bECHit_mvx = 0;
TBranch *bECHit_mvy = 0;
TBranch *bECHit_mvz = 0;
TBranch *bECHit_t = 0;
Int_t fEv_simpid;
Double_t fEvTimeStamp;
Int_t fEvTimeWindow;
TChain * TGEMC_EC;
Bool_t kSaveCanvas=kTRUE;

Int_t pid_gen=0;
Double_t theta_gen=0,phi_gen=0,p_gen=0,E_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;
Int_t FluxHit_detector_ID,FluxHit_subdetector_ID,FluxHit_subsubdetector_ID,FluxHit_component_ID;
Int_t ECHit_detector_ID,ECHit_subdetector_ID,ECHit_subsubdetector_ID,ECHit_component_ID;
Int_t ECPSHit_detector_ID,ECPSHit_subdetector_ID,ECPSHit_subsubdetector_ID,ECPSHit_component_ID;
Double_t trig_low_R[8]={110.0 ,130.0 ,150.0 ,170.0 ,190.0 ,210.0 ,230.0 ,250.0};
Double_t trig_high_R[8]={130.0 ,150.0 ,170.0 ,190.0 ,210.0 ,230.0 ,270.0 ,270.0};
Double_t trig_thresh_6p1[8] = {375.156,310.137,260.115,244.673,222.564,178.257,151.943,151.943};
Double_t trig_thresh_2p1[6] = {501.5 ,471.9 ,412.8 ,340.5 };//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0  };//not in use
Double_t trig_thresh_PS[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{20.9 ,28.2 ,28.3 ,27.7 ,27.5 ,29.0 ,31.7 ,17.7};//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0  };//not in use
Double_t trig_thresh_PS_6p1[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{20.9 ,28.2 ,28.3 ,27.7 ,27.5 ,29.0 ,31.7 ,17.7 };//{126.6 ,128.7 ,123.6 ,124.4 ,124.1 ,125.0 ,126.4 ,119.7};////Using 6p1 cluster sum as total PS sum threshold
Double_t trig_thresh_PS_2p1[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{121.9 ,120.8 ,115.0 ,115.3 ,116.4 ,116.3 ,117.5 ,131.4};

//routine to load ecal block id,coordinates, and sector information
void LoadEC_map(TString map_file);
//routine to return ecal block x,y
TVector2 GetECALBlock_coord(Int_t block_id);
Int_t GetECALBlock_sector(Int_t block_id);
//return an pointer to energy deposit for 6 neighboring ecal blocks
Int_t GetECALCluser(Int_t block_id,Int_t *cluster_edep_blockid);
//Int_t GetECALCluser(Double_t hit_x, Double_t hit_y, Int_t *cluster_edep_blockid);
//save x,y,id and sectors and status for which moudles are active only for SIDIS case
Int_t sector[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t id[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t num_module_in_row[54]={0};
Double_t y_bak[54]={100000};
Double_t x[54][50]={{100000}};           
Double_t y[54][50]={{100000}}; 
Int_t status[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
std::pair<Int_t,Int_t> block_map[2000];
Int_t GetRadiusIndex(Double_t radius);
Double_t GetThreshold6p1(Double_t radius);
Double_t Calfactor(Double_t edep);
//routines for properly access vector based TTree using TChain
void SetFluxBranchAddresses();
void GetFluxEntryByBranch(Long64_t local);

void SetECBranchAddresses();
void GetECEntryByBranch(Long64_t local);
Double_t edep_6p1_max1[__OCTANTS]={0};
Double_t radius_6p1_max1[__OCTANTS]={0};
Double_t edep_6p1_preshower_max1[__OCTANTS]={0};
Double_t radius_6p1_preshower_max1[__OCTANTS]={0};
Double_t edep_6p1_max1_bkg=0;  
Double_t edep_6p1_premax1_bkg=0;  
Int_t block_count_bkg=0;
Int_t cluster_edep_blockid_bkg[7]={0};
Double_t cluster_edep_bkg[7]={0};
Double_t cluster_PS_edepbkg[7]={0};
Double_t edep_6p1_x[__OCTANTS]={0};
Double_t edep_6p1_y[__OCTANTS]={0};
Double_t edep_6p1_x_preshower[__OCTANTS]={0};
Double_t edep_6p1_y_preshower[__OCTANTS]={0};
//get energy resolution
TGraphErrors * GetRMSReolution(TH2F * Difference);
struct event {
	Int_t pid;
	Double_t x;
	Double_t y;
	Double_t pf; 
	//  Double_t hit_time; 
	Int_t blockID;
	Double_t Q2;
};

struct trigger_event {
	TString pid;
	Double_t pf;  
	Double_t r;
	Double_t x;
	Double_t y;
};

struct ecal_event {
	Double_t edep;  
	Double_t r;
	Double_t x;
	Double_t y;
};

Bool_t bEarlyBeak = kFALSE;//kTRUE;
vector< vector<trigger_event> > event_summary;
vector< vector<ecal_event> > ecal_summary;
trigger_event trigger_event1;
ecal_event ecal_event1;
using namespace std;
using namespace std::tr1;
void LoadBackgrounds();
void getBkgEdepPS(Int_t rndnum);//this routine will fill ecalPSMap with background signals
void getBkgEdepSh(Int_t rndnum);//this routine will fill ecalShMap with background signals
vector<std::map<Int_t,Double_t> >total_ecalShMap;
vector<std::map<Int_t,Double_t> > total_ecalPSMap;
vector<std::map<Int_t,Double_t> >total_ecalShMapbkg;
vector<std::map<Int_t,Double_t> > total_ecalPSMapbkg;
std::map<Int_t,Double_t> ecalPSMap;//map for pre-shower
std::map<Int_t,Double_t> ecalShMap;//map for shower
std::map<Int_t,Double_t> ecalPSMapbkg;//map for pre-shower
std::map<Int_t,Double_t> ecalShMapbkg;//map for shower
void set_plot_style();

std::map<int,int> pidmap;
vector<event> event_list;//use this to store events in a single time window
int main(Int_t argc,Char_t* argv[]) { 
	//set the timer
	TStopwatch watch;
	watch.Start();


	TVirtualFitter::SetDefaultFitter("Minuit2");

	TApplication theApp("App",&argc,argv);
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat("eMr");
	gStyle->SetOptFit();
	set_plot_style();

	Int_t pid_ec;

	//bin resolution for histograms
	Int_t mom_bins = 100;
	//pid map to access histograms pid=0 is returned when invalid pid is passed
	pidmap[11]=1;
	pidmap[-211]=2;
	pidmap[22]=3;
	pidmap[-11]=4;
	pidmap[211]=5; 
	TString rootfilename="ecal_pimbkg_eAll_trigger_Eend_pL1800_11GeV_2023_Q2L2X30p1600.root";

	if (kSaveRootFile){
		TString rootfilestatus="RECREATE";
		rootfile = new TFile(rootfilename, rootfilestatus);
		rootfile->cd();
	} 

	TH2F *hPlog_R=new TH2F("hPlog_R","hPlog_R",300, 0, 300, 200,-6,1.3);    
	hPlog_R->SetTitle(";R (cm);log(P) (GeV)");          


	TH2F *hElog_R=new TH2F("hElog_R","hElog_R",300, 0, 300, 200,-6,1.3);    
	hElog_R->SetTitle(";R (cm);log(E) (GeV)");
	TH2F *hEklog_R=new TH2F("hEklog_R","hEklog_R",300, 0, 300, 200,-6,1.3);
	hEklog_R->SetTitle(";R (cm);log(Ek) (GeV)");  
	TH1F *hfluxR=new TH1F("hfluxR","hfluxR",60,0,300);
	hfluxR->SetTitle(";R (cm);flux (kHz/mm2)");
	TH1F *hEfluxR=new TH1F("hEfluxR","hEfluxR",30,0,300);
	hEfluxR->SetTitle(";R (cm);Eflux (1e3*GeV/100cm2/s)");

	TH1F *Gen_mom_distr = new TH1F("Gen_mom_distr","Primary Track Momentum;Momentum (MeV)",mom_bins ,1000,8000);
	TH1F *Gen_theta_distr = new TH1F("Gentheta_distr","Primary Track Theta;#theta (deg)",mom_bins ,10,50);  
	TH2F *Gen_theta_vs_p = new TH2F("Gen_theta_vs_p","#pi^{-} Theta vs Mom;theta (deg);p (GeV/c)",100,22,35,100,0,7);  
	TH2F *Gen_theta_vs_p_EC = new TH2F("Gen_theta_vs_p_EC","#pi^{-} Theta vs Mom;theta (deg);p (GeV/c)",100,22,35,100,0,7);  
	TH2F *Gen_theta_vs_p_EC_accep = new TH2F("Gen_theta_vs_p_EC_accep","#pi^{-} Theta vs Mom; r_EC (cm);p (GeV/c)",200,100,300,100,0,7);  
	TH2F *Gen_theta_vs_p_Xcut = new TH2F("Gen_theta_vs_p_Xcut","#pi^{-} Theta vs Mom;theta (deg);p (GeV/c)",100,22,35,100,0,7);  
	TH2F *Gen_theta_vs_p_pcut = new TH2F("Gen_theta_vs_p_pcut","#pi^{-} Theta vs Mom;theta (deg);p (GeV/c)",100,22,35,100,0,7);  
	TH2F *Gen_theta_vs_p_Q2Xpcut = new TH2F("Gen_theta_vs_p_Q2Xpcut","#pi^{-} Theta vs Mom;theta (deg);p (GeV/c)",100,22,35,100,0,7);  
	TH1F *ECALDet_mom_distr = new TH1F("ECALDet_mom_distr","ECAL Front All e^{-} Track Momentum;Momentum (MeV)",mom_bins ,-8000,1000);
	TH1F *htotEdep_ec_shower=new TH1F("htotEdep_ec_shower","ec shower;totEdep(MeV);",5000,0,5000);
	TH1F *htotEdep_ec_preshower=new TH1F("htotEdep_ec_preshower","ec preshower;totEdep(MeV);",100,0,500);
	TH1F *BlockX_hist=new TH1F("BlockX_hist",";#Delta X(cm);",1000,-500,500);
	TH1F *BlockY_hist=new TH1F("BlockY_hist",";#Delta Y(cm);",1000,-500,500);
	TString strig[3] = {"No. Trig","6+1 Trig.","2+1 Trig."};
	TString sradii[9] = {"1.1 - 1.3 m","1.3 - 1.5 m","1.5 - 1.7 m","1.7 - 1.9 m","1.9 - 2.1 m","2.1 - 2.3 m","2.3 - 2.7 m","2.5 - 2.7 m","1.1 - 2.7 m"};
	//Only index 0 - e+/-, 1 - pions. 2 - gamma are used in count_pid
	Int_t count_pid[30][3][5]={{{0}}};//eventually 5 pid type (e-, pi-,gamma,e+, and pi+) to count total, 6+1 trigger and 2+1 trigger would be used
	TH1F *histo_edep_Sh_radius[3][9];//total,6+1 trigger, 2+1 trigger for 8 radius bins and all radiifluxHit_totEdep
	TH1F *histo_edep_Sh_cluster_radius[3][9];//total,6+1 trigger, 2+1 trigger for 8 radius bins and all radii
	TH1F *histo_edep_PS_radius[3][9];
	TH1F *histo_edep_PS_cluster_radius[3][9];
	TH1F *histo_edep_Sh_radius_threshold[9];
	TH1F *histo_edep_Sh_radius_nothreshold[9];
	TH1F *histo_edep_Sh_radius_nothreshold_nobkg[9];
	TH1F *histo_E_Sh_radius_threshold_nobkg[9];
	TH1F *histo_trig_window_residue_r[9]; //
	TH1F *histo_trig_window_residue_x[9]; //
	TH1F *histo_trig_window_residue_y[9]; //
	TH1F *histo_edep_Sh_radius_threshold_cut[9];
	TH1F *histo_edep_Sh_radius_nothreshold_cut[9];
	TH1F *histo_edep_Sh_radius_threshold_nobkg[9];
	TH2F *Shower_vs_timewindows_hist;
	TFile *file[5];
	TTree *tree_generated[5];
	TTree *chain_Flux[5];
	TTree *TGEMC_EC[5];
	TTree *tree_header[5];
	TTree *TGEMC_ECPS[5]; 
	Shower_vs_timewindows_hist  = new TH2F("timewindow_vs_sum6p1Ebkg_hist",";E_sum (GeV); Timewondow number",120,-0.5,11.5,108,-0.5,107.5);

	for(Int_t k=0;k<8;k++){
		histo_edep_Sh_radius_threshold[k]  = new TH1F(Form("histo_edep_Sh_radius_threshold_%d",k),Form("shower_6p1E_R%d",k),8000,0,8000);//for Pf<1 change bin limit to 1500
		histo_edep_Sh_radius_threshold_nobkg[k]  = new TH1F(Form("histo_edep_Sh_radius_threshold_nobkg_%d",k),Form("shower_6p1E_R%d",k),8000,0,8000);//for Pf<1 change bin limit to 1500
		histo_edep_Sh_radius_nothreshold[k]  = new TH1F(Form("histo_edep_Sh_radius_nothreshold_%d",k),Form("PShower_6p1E_R_%d",k),300,0,300);
		histo_edep_Sh_radius_nothreshold_nobkg[k]  = new TH1F(Form("histo_edep_Sh_radius_nothreshold_nobkg_%d",k),Form("PShower_6p1E_R_nobkg%d",k),300,0,300);
		histo_E_Sh_radius_threshold_nobkg[k]  = new TH1F(Form("histo_E_Sh_radius_threshold_nobkg_%d",k),Form("Eflux_shower_R%d",k),11000,0,11000);//for Pf<1 change bin limit to 1500
		histo_edep_Sh_radius_threshold[k]->GetXaxis()->SetTitle("E_6p1 [MeV]");
		histo_edep_Sh_radius_threshold[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_edep_Sh_radius_nothreshold[k]->GetXaxis()->SetTitle("E_pre [MeV]");
		histo_edep_Sh_radius_nothreshold[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_edep_Sh_radius_nothreshold_nobkg[k]->GetXaxis()->SetTitle("E_pre [MeV]");
		histo_edep_Sh_radius_nothreshold_nobkg[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_edep_Sh_radius_threshold_cut[k]  = new TH1F(Form("histo_edep_Sh_radius_threshold_cut_%d",k),Form("shower_6p1E_R%d",k),8000,0,8000);//for Pf<1 change bin limit to 1500
		histo_edep_Sh_radius_nothreshold_cut[k]  = new TH1F(Form("histo_edep_Sh_radius_nothreshold_cut_%d",k),Form("PShower_6p1E_R_%d",k),300,0,300);
		histo_edep_Sh_radius_threshold_cut[k]->GetXaxis()->SetTitle("E_6p1 [MeV]");
		histo_edep_Sh_radius_threshold_cut[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_edep_Sh_radius_nothreshold_cut[k]->GetXaxis()->SetTitle("E_6p1 [MeV]");
		histo_edep_Sh_radius_nothreshold_cut[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_edep_Sh_radius_nothreshold[k]->GetXaxis()->SetTitle("E_pre [MeV]");
		histo_edep_Sh_radius_nothreshold[k]->GetYaxis()->SetTitle("rate [Hz]");
		histo_E_Sh_radius_threshold_nobkg[k]->GetXaxis()->SetTitle("E_flux [MeV/c]");
		histo_E_Sh_radius_threshold_nobkg[k]->GetYaxis()->SetTitle("rate [Hz]");
	}
	Int_t ecal_bin_limit[4]={2000,2000,400,400};// for blocker shower total //{2500,1000,500,400};//for no photon blocker //{2000,600,400,250};// for blocker shower total,shower 6+1, PS total, PS 6+1 for P>1 1000,750,200,200 and for p<1 or total {2500,1000,500,400}
for(Int_t i=0;i<3;i++){
	for(Int_t j=0;j<9;j++){
		//shower total energy all clusters in the octant. The total sum variation with R is difficult to define since R is extracted based on R of the max cluster
		histo_edep_Sh_radius[i][j] = new TH1F(Form("histo_edep_Sh_radius_%d_%d",i,j),Form("Sh. Total energy deposit in 30 ns (R : %s, %s);edep (MeV)",sradii[j].Data(),strig[i].Data()),100,0,ecal_bin_limit[0]);
		histo_edep_Sh_cluster_radius[i][j] = new TH1F(Form("histo_edep_Sh_cluster_radius_%d_%d",i,j),Form("Sh. 6+1 energy deposit in 30 ns (R : %s, %s);edep (MeV)",sradii[j].Data(),strig[i].Data()),100,0,ecal_bin_limit[1]);
		//preshower total energy all clusters in the octant. The total sum variation with R is difficult to define since R is extracted based on R of the max cluster
		histo_edep_PS_radius[i][j] = new TH1F(Form("histo_edep_PS_radius_%d_%d",i,j),Form("PS. Total energy deposit in 30 ns (R : %s, %s);edep (MeV)",sradii[j].Data(),strig[i].Data()),100,0,ecal_bin_limit[2]);
		histo_edep_PS_cluster_radius[i][j] = new TH1F(Form("histo_edep_PS_cluster_radius_%d_%d",i,j),Form("PS. 6+1 energy deposit in 30 ns (R : %s, %s);edep (MeV)",sradii[j].Data(),strig[i].Data()),100,0,ecal_bin_limit[3]);
		if (i==0){
			histo_trig_window_residue_r[j]= new TH1F(Form("histo_trig_window_residue_r_%d",j),Form("Track vs. Cluster R Variance for %s; R_{track} - R_{6p1} (cm)",sradii[j].Data()),100,-50,50);
			histo_trig_window_residue_x[j]= new TH1F(Form("histo_trig_window_residue_x_%d",j),Form("Track vs. Cluster X Variance for %s; X_{track} - X_{6p1} (cm)",sradii[j].Data()),100,-50,50);
			histo_trig_window_residue_y[j]= new TH1F(Form("histo_trig_window_residue_y_%d",j),Form("Track vs. Cluster Y Variance for %s; Y_{track} - Y_{6p1} (cm)",sradii[j].Data()),100,-50,50);
		}
	}
}

TH1F *histo_pid[3][5];//momentum distr on incident particles on ecal sens detector for background 
TH1F *histo_pid2[3][5][9];//momentum distr on incident particles on ecal sens detector for flat elec. + background merged
TH1F *histo_pid2_cut[3][5][9];//momentum distr on incident particles on ecal sens detector for flat elec. + background merged
TH1F *histo_pid_time[3][5];//hit time distr on incident particles on ecal sens detector 

TH1F *histo_delayed_pid[3][5];//momentum distr on incident particles on ecal sens detector 
TH1F *histo_delayed_pid_time[3][5];//hit time distr on incident particles on ecal sens detector 
TH1F *histo_pid_time_energyW[5];//hit time distr on incident particles on ecal sens detector 


TH2F *histo_rphi[3][5];//plot xy hit distr on ecal sens detector for e+-, any pion, gamma
TString spid[5] = {"Electron","Pion","Pi0 Gamma","Gamma","empty"};//e+-,Pi+-,gamma
Int_t ecal_bin_limit2[5]={8000,8000,1200,1200,8000};//for  blocker;  //for no blocker //{1000,5000,1200,1200,8000};//for blocker //{1000,1000,1000,1000,1000};//{10,10,10,10,10};//for no blocker {1200,5000,2000,1200,8000}
//for use with histograms for flat elec. + background merged simulations
Int_t ecal_bin_limit3[5]={1000,8000,2000,8000,8000};//for  blocker;
Int_t ecal_bin_limit3_low[5]={0,1000,0,1000,0};
TString spid3[5] = {"Bkg. e^{#pm}","Pion","Gamma","e^{-}","empty"};//e+-,Pi+-,gamma
//cout<<"FAL======"<<spid3[1].Data()<<endl;
for(Int_t i=0;i<3;i++){
	for(Int_t j=0;j<5;j++){
		histo_pid[i][j]  = new TH1F(Form("Histo_pid_%d_%d",i,j),Form("%s Momentum (%s);Momentum (MeV)",spid[j].Data(),strig[i].Data()),mom_bins ,0,ecal_bin_limit2[j]);//for Pf<1 change bin limit to 1500
		for(Int_t k=0;k<9;k++){
			histo_pid2[i][j][k]  = new TH1F(Form("Histo_pid2_%d_%d_%d",i,j,k),Form("%s Momentum (%s);Momentum (MeV) (R %s)",spid3[j].Data(),strig[i].Data(),sradii[k].Data()),35,ecal_bin_limit3_low[j],ecal_bin_limit3[j]);//for Pf<1 change bin limit to 1500
			histo_pid2_cut[i][j][k]  = new TH1F(Form("Histo_pid2_cut_%d_%d_%d",i,j,k),Form("%s Momentum (%s);Momentum (MeV) (R %s)",spid3[j].Data(),strig[i].Data(),sradii[k].Data()),35,ecal_bin_limit3_low[j],ecal_bin_limit3[j]);//for Pf<1 change bin limit to 1500
		}
		histo_pid_time[i][j]  = new TH1F(Form("Histo_pid_time_%d_%d",i,j),Form("%s Hit Time  (%s);Time (ns)",spid[j].Data(),strig[i].Data()),60 ,0,30);
		histo_rphi[i][j]  = new TH2F(Form("Histo_rphi_%d_%d",i,j),Form("%s Hit Distribution (%s); #phi (deg); R (cm)",spid[j].Data(),strig[i].Data()),48,0,12,200 ,100,300);//phi range -6,6
		fTotalRate[i][j] = 0;

		histo_delayed_pid[i][j]  = new TH1F(Form("Histo_delayed_pid_%d_%d",i,j),Form("%s Momentum (%s) Delayed Hit;Momentum (MeV)",spid[j].Data(),strig[i].Data()),mom_bins ,0,100);
		histo_delayed_pid_time[i][j]  = new TH1F(Form("Histo_delayed_pid_time_%d_%d",i,j),Form("%s Delayed Hit Time  (%s);Time (ns)",spid[j].Data(),strig[i].Data()),1000 ,0,10000);
	}    
}
for(Int_t j=0;j<5;j++){
	histo_pid_time_energyW[j] = new TH1F(Form("Histo_pid_time_energyW_%d",j),Form("%s Hit Time Energy Weighted ;Time (ns)",spid[j].Data()),200 ,0,500);
}




TH1F *histo_PS_edep = new TH1F("histo_PS_edep","PS energy deposit over 30 ns;edep (MeV)",100,0,500);


TH1F *histo_octant_trig[2];//for 6+1 and 2+1
histo_octant_trig[0] = new TH1F("histo_octant_trig_6p1","Octant Level 6+1 Trigger Multiplicity;No. of Octants Trig.",30,0,30);
histo_octant_trig[1] = new TH1F("histo_octant_trig_2p1","Octant Level 2+1 Trigger Multiplicity;No. of Octants Trig.",30,0,30);
Int_t ocantcount_6p1;//count no.of octants 6+1 fired
Int_t ocantcount_2p1;//count no.of octants 2+1 fired

TH1F *histo_pid_count[3][5];
Int_t pid_bin_limit[5]={20,10,250,100,100};//for pf>1 change gamma range to 10 for e- and gamma //for total and pf < 1 set e- 20 gamma 25
TString spid2[5] = {"e^{#pm}","#pi^{#pm}","#pi^{0} #gamma","#gamma","#pi^{+}"};
for(Int_t i=0;i<3;i++){
	for(Int_t j=0;j<5;j++){
		histo_pid_count[i][j]=new TH1F(Form("Histo_pid_count_%d_%d",i,j),Form("%s Count (%s);Count per Window",spid2[j].Data(),strig[i].Data()),10 ,0,pid_bin_limit[j]);
	}
}

//Energy resolution plots
Double_t EoverPElec_low_limit = 0;//0. for electron shower, -0.1 for pre shower//-2. for pions
Double_t EoverPElec_up_limit = 0.25;//0.5 for electronshower, 0.1 for pre shower//2. for pions
Double_t E_low_limit = 0;//0. for electron shower, -0.1 for pre shower//-2. for pions
Double_t E_up_limit = 0;//0.5 for electronshower, 0.1 for pre shower//2. for pions
EoverPElec_low_limit = 1.2;//0.90; defalut from remoll
EoverPElec_up_limit = 1.5;//1.05; defalut from remoll
E_low_limit = 1000;//0.90; defalut from remoll
E_up_limit = 12000;//1.05; defalut from remoll
//for Pre-Shower+Shower Cluster
TH2F * EoverPElec_total_ecal_fullsum = new TH2F("EoverPElec_total_ecal_fullsum","Total PS+Sh Edep over Pf Ratio;Momentum (GeV); Edep/Pf",mom_bins,2,7,100,EoverPElec_low_limit,EoverPElec_up_limit);
//EoverPElec_low_limit = 0.98;//for sampling fraction related use  0.85 to 1.1 limit//for all simulation based estimation use 0.985 to 1.0
//EoverPElec_up_limit = 0.995;  //with holes use
//EoverPElec_low_limit = 0.85;

TH2F * EoverPElec_total_ecal_6p1sum = new TH2F("EoverPElec_total_ecal_6p1sum","6+1 PS+Sh Edep over Pf Ratio;Momentum (GeV); Edep/Pf",mom_bins,2,7,100,EoverPElec_low_limit,EoverPElec_up_limit);

//sampling frac for lead
Double_t f_samp = 0.23;//0.199776;
Double_t f_samp_PS = 0.290808;
event single_event;//for GEM hits in a single time window
event_list.clear();//empty the event list vector
TH2F * E_total_ecal_fullsum = new TH2F("E_total_ecal_fullsum","Total PS+Sh Edep;Momentum (GeV); Edep",mom_bins,2,7,1000,E_low_limit,E_up_limit);
//EoverPElec_low_limit = 0.98;//for sampling fraction related use  0.85 to 1.1 limit//for all simulation based estimation use 0.985 to 1.0
//EoverPElec_up_limit = 0.995;  //with holes use
//EoverPElec_low_limit = 0.85;

TH2F * E_total_ecal_6p1sum = new TH2F("E_total_ecal_6p1sum","6+1 PS+Sh Edep;Momentum (GeV); Edep",mom_bins,2,7,1000,E_low_limit,E_up_limit);

if (!bDisableBkg)
	LoadBackgrounds();
	else
	printf("Backgrounds are disabled! \n");
	// rod 6mm eAll
	file[0]=new TFile("/cache/halla/solid/user/tianye/container/eAll_commited4fe_rod_6mm/eAll_commited4fe_rod_6mm_95files.root");
	file[1]=new TFile("/cache/halla/solid/user/tianye/container/eAll_commited4fe_rod_6mm_2/eAll_commited4fe_rod_6mm_2_100files_1e6.root");
	int totalN[2]={95,100};
double event_time=0.0;
double time_windows=0;
int N_loss=0;
int N_total=0;
Int_t nentries[5];
double areaR=0;
double weightR=0;
double r_flux=0;
double flux_E=0;
double dE=0;
Double_t pf;
double total_rate_eDIS=0;
for(int n=0;n<2;n++){

	tree_header[n] = (TTree*) file[n]->Get("userHeader");
	tree_header[n]->SetBranchAddress("userVar001", &userVar001, &b_userVar001);
	tree_header[n]->SetBranchAddress("userVar002", &userVar002, &b_userVar002);
	tree_header[n]->SetBranchAddress("userVar003", &userVar003, &b_userVar003);
	tree_header[n]->SetBranchAddress("userVar004", &userVar004, &b_userVar004);
	tree_header[n]->SetBranchAddress("userVar005", &userVar005, &b_userVar005);
	tree_header[n]->SetBranchAddress("userVar006", &userVar006, &b_userVar006);
	tree_header[n]->SetBranchAddress("userVar007", &userVar007, &b_userVar007);
	tree_header[n]->SetBranchAddress("userVar008", &userVar008, &b_userVar008);
	tree_header[n]->SetBranchAddress("userVar009", &userVar009, &b_userVar009);
	tree_header[n]->SetBranchAddress("userVar010", &userVar010, &b_userVar010);

	tree_generated[n] = (TTree*) file[n]->Get("generated");

	tree_generated[n]->SetBranchAddress("pid",&fGen_pid);
	tree_generated[n]->SetBranchAddress("px",&fGen_Px);
	tree_generated[n]->SetBranchAddress("py",&fGen_Py);
	tree_generated[n]->SetBranchAddress("pz",&fGen_Pz);
	tree_generated[n]->SetBranchAddress("vx",&fGen_vx);
	tree_generated[n]->SetBranchAddress("vy",&fGen_vy);
	tree_generated[n]->SetBranchAddress("vz",&fGen_vz);

	chain_Flux[n] = (TTree*) file[n]->Get("flux");
	chain_Flux[n]->SetBranchAddress("hitn",&fFluxHit_n);
	chain_Flux[n]->SetBranchAddress("id",&fFluxHit_id);
	chain_Flux[n]->SetBranchAddress("pid",&fFluxHit_pid);
	chain_Flux[n]->SetBranchAddress("mpid",&fFluxHit_mpid);
	chain_Flux[n]->SetBranchAddress("mtid",&fFluxHit_mtid);
	chain_Flux[n]->SetBranchAddress("tid",&fFluxHit_tid);
	chain_Flux[n]->SetBranchAddress("trackE",&fFluxHit_trackE);
	chain_Flux[n]->SetBranchAddress("totEdep",&fFluxHit_totEdep);
	chain_Flux[n]->SetBranchAddress("avg_x",&fFluxHit_Avg_x);
	chain_Flux[n]->SetBranchAddress("avg_y",&fFluxHit_Avg_y);
	chain_Flux[n]->SetBranchAddress("avg_z",&fFluxHit_Avg_z);
	chain_Flux[n]->SetBranchAddress("px",&fFluxHit_Px);
	chain_Flux[n]->SetBranchAddress("py",&fFluxHit_Py);
	chain_Flux[n]->SetBranchAddress("pz",&fFluxHit_Pz);
	chain_Flux[n]->SetBranchAddress("avg_t",&fFluxHit_T);
	chain_Flux[n]->SetBranchAddress("vz",&fFluxHit_vz);

	TGEMC_EC[n] = (TTree*) file[n]->Get("solid_ec");

	TGEMC_EC[n]->SetBranchAddress("id",&fECHit_id);
	TGEMC_EC[n]->SetBranchAddress("totEdep",&fECHit_totEdep);
	TGEMC_EC[n]->SetBranchAddress("totEdepB",&fECHit_totEdepB);
	TGEMC_EC[n]->SetBranchAddress("totEend",&fECHit_totEend);
	TGEMC_EC[n]->SetBranchAddress("avg_z",&fECHit_Avg_z);
	TGEMC_EC[n]->SetBranchAddress("mpid",&fECHit_mpid);
	TGEMC_EC[n]->SetBranchAddress("tid",&fECHit_tid);
	TGEMC_EC[n]->SetBranchAddress("mtid",&fECHit_mtid); 
	TGEMC_EC[n]->SetBranchAddress("avg_lx",&fECHit_Avg_lx);
	TGEMC_EC[n]->SetBranchAddress("avg_ly",&fECHit_Avg_ly);
	TGEMC_EC[n]->SetBranchAddress("avg_lz",&fECHit_Avg_lz);
	TGEMC_EC[n]->SetBranchAddress("avg_x",&fECHit_Avg_x);
	TGEMC_EC[n]->SetBranchAddress("avg_y",&fECHit_Avg_y);
	TGEMC_ECPS[n] = (TTree*) file[n]->Get("solid_ec_ps");

	TGEMC_ECPS[n]->SetBranchAddress("id",&fECPSHit_id);
	TGEMC_ECPS[n]->SetBranchAddress("totEdep",&fECPSHit_totEdep);
	TGEMC_ECPS[n]->SetBranchAddress("totEdepB",&fECPSHit_totEdepB);
	TGEMC_ECPS[n]->SetBranchAddress("totEend",&fECPSHit_totEend);
	TGEMC_ECPS[n]->SetBranchAddress("avg_z",&fECPSHit_Avg_z);
	TGEMC_ECPS[n]->SetBranchAddress("mpid",&fECPSHit_mpid);
	TGEMC_ECPS[n]->SetBranchAddress("tid",&fECPSHit_tid);
	TGEMC_ECPS[n]->SetBranchAddress("mtid",&fECPSHit_mtid); 
	TGEMC_ECPS[n]->SetBranchAddress("avg_lx",&fECPSHit_Avg_lx);
	TGEMC_ECPS[n]->SetBranchAddress("avg_ly",&fECPSHit_Avg_ly);
	TGEMC_ECPS[n]->SetBranchAddress("avg_lz",&fECPSHit_Avg_lz);
	TGEMC_ECPS[n]->SetBranchAddress("avg_x",&fECPSHit_Avg_x);
	TGEMC_ECPS[n]->SetBranchAddress("avg_y",&fECPSHit_Avg_y);

	nentries[n] = (Int_t)tree_generated[n]->GetEntries();
	printf("Entries = %i \n",nentries[n]);
	Int_t treenumber=-1;
	Int_t ECtreenumber=-1;

	Bool_t Fillecal=kFALSE;  

	//flux hits
	Double_t pf;//for ECAL front det
	Double_t pf_flux;
	Double_t E_flux;

	Double_t th;//for ECAL front det
	Double_t r;//for ECAL front det
	Double_t phi;
	Int_t i_pf=0;//for GEM plane hit index  
	TRandom *r3 = new TRandom3();
	ecalPSMap.clear();
	ecalPSMap.clear();
	ecalShMapbkg.clear();
	ecalShMapbkg.clear();

	//load EC map file
	LoadEC_map("../layout/map_FAEC_ANL_20130628.txt");

	//vector to store ecal hit coordinates
	TVector2 vec_ecalBlock;
	Double_t ecalBlock_x;
	Double_t ecalBlock_y;

	TVector2 vec_ecalBlock_preshower;
	Double_t ecalBlock_x_preshower;
	Double_t ecalBlock_y_preshower;

	TVector2 vec_ecalBlock_bkg;
	Double_t ecalBlock_x_bkg;
	Double_t ecalBlock_y_bkg;
	TVector2 vec_ecalCluster;
	Double_t ecalCluster_x;
	Double_t ecalCluster_y;
	Int_t ecal_clusterID;
	Int_t ecal_blockID;
	Int_t ecal_blockID_bkg;

	Int_t ecal_blockID_psh;  
	Int_t ecal_hitblockID_psh;  
	Int_t ecal_hitblockID;
	Double_t ecalBlock_x0;
	Double_t ecalBlock_y0;
	Int_t ecal_blockID0;
	Int_t ecal_blockID0_preshower;
	Double_t ecalBlock_z0;
	Double_t ecalPos_x0;//store main(center) block location of the 6+1 cluster
	Double_t ecalPos_y0;
	Double_t ecalBlock_x0_2p1;//store main(center) block location of the 2+1 cluster
	Double_t ecalBlock_y0_2p1;
	Double_t ecal_offset=3.5;//global center of the ecal wrt remoll setup in m
	Double_t targ_offset=0.1;//in m
	Int_t octantno;
	Double_t ecalClustBlock_x;
	Double_t ecalClustBlock_y;
	Double_t cluster_edep[7]={0};
	Int_t cluster_edep_blockid[7]={0};
	Double_t cluster_PS_edep[7]={0};
	Int_t cluster_edep_PS_blockid[7]={0};
	Int_t block_count=0;
	Double_t edep_6p1_sum_bkg;  //sum of 6+1 blocks 
	Double_t edep_6p1_presum_bkg;  //sum of 6+1 blocks 
	Int_t block_count_PS=0; 
	Int_t block_sector=0;
	Double_t edep_6p1_sum;  //sum of 6+1 blocks
	Double_t edep_6p1_PS_sum;  //sum of 6+1 blocks  
	Double_t edep_total_sum;  //sum ofall the ecal blocks 
	Double_t edep_1block_sum;//photon sum of the main ecal block
	Double_t edep_1block_sum_preshower;//photon sum of the main ecal block
	Double_t edep_1block_sum_bkg; 

	// Double_t edep_6p1_max=0;
	Double_t radius_6p1_max=0;
	Double_t radius_6p1=0;
	Double_t radius_6p1_bkg=0;
	Double_t radius_6p1_preshower=0;
	Double_t edep_sum0;  //photon sum of the  ecalBlock_x0,ecalBlock_y0 ecal block
	Double_t DeltaRrDiff;
	Double_t edep_total_PS_sum;  //sum ofall the ecal blocks 
	//Double_t edep_2p1_PS_sum;  //sum of 6+1 blocks 
	Double_t X;
	Double_t Y;
	int eventid=0;
	Bool_t processWindow=kFALSE;//when all the events in a single window are filled to the ecal, start cluster processing
	Fillecal=kFALSE;

	//start dissecting this ecal with a primary track
	fTimeWindow=0;//use to keep track of current time window
	edep_1block_sum = 0;
	edep_1block_sum_bkg= 0;
	edep_1block_sum_preshower = 0; 
	edep_sum0 = 0;
	edep_total_sum = 0;
	Double_t totEdep_shower=0;
	Double_t totEdep_preshower=0;
	//PS ecal cluster information
	edep_6p1_PS_sum=0;
	edep_total_PS_sum=0;	
	double pidmass=0.511/1000;//MeV electron
	double EC_averagex=0.0;
	double EC_averagey=0.0;
	double EC_localx=0.0;
	double EC_localy=0.0;
	int flux_size=0;
	double phi_flux=0;
	double E_flux_max= 0;


	//double edep_total_PS_sum=0;
	double edep_total_Sh_sum=0;
	// double event_rate=50e-6/1.6e-19/1e9;//EM rate
	double event_rate=0;//pi0 rate
	double Q2=0;//pi0 rate
	double Xjb=0;//pi0 rate
	int event_sector=0;
	/// in mm2
	for (Int_t i=0; i<nentries[n]; i++) {
		double pf_satisfy=0;
		Long64_t local;

		tree_header[n]->GetEntry(i);
		//event_rate= var8->at(0)/500;
		event_rate= userVar006->at(0)/(195.);
		Q2= userVar005->at(0);
		Xjb= userVar002->at(0);
		//event_rate= userVar010->at(0)/1000;
		pf=0;
		total_rate_eDIS += event_rate;
		tree_generated[n]->GetEntry(i);
		for (int j=0;j<fGen_pid->size();j++) {
			pid_gen=fGen_pid->at(j);
			px_gen=fGen_Px->at(j);
			py_gen=fGen_Py->at(j);
			pz_gen=fGen_Pz->at(j);
			vx_gen=fGen_vx->at(j);
			vy_gen=fGen_vy->at(j);
			vz_gen=fGen_vz->at(j);
			p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
			E_gen=sqrt(pow(p_gen,2)-pow(pidmass,2))*0.001;
			theta_gen=TMath::ACos(pz_gen/p_gen)*DEG;
			phi_gen=TMath::ATan2(py_gen,px_gen)*DEG;
			Gen_mom_distr->Fill(p_gen);
			Gen_theta_distr->Fill(theta_gen);
			Gen_theta_vs_p->Fill(theta_gen,p_gen*0.001,event_rate);
			if(p_gen>1800){
				Gen_theta_vs_p_pcut->Fill(theta_gen,p_gen*0.001,event_rate);
				if(Xjb>0.3){
					Gen_theta_vs_p_Xcut->Fill(theta_gen,p_gen*0.001,event_rate);
					if(Q2>2){
						Gen_theta_vs_p_Q2Xpcut->Fill(theta_gen,p_gen*0.001,event_rate);
					}
				}
			}
		}
		chain_Flux[n]->GetEntry(i);
		E_flux_max=0;
		for (int j = 0; j<fFluxHit_id->size(); j++){
			// 
			FluxHit_detector_ID= fFluxHit_id->at(j)/1000000;
			FluxHit_subdetector_ID=(fFluxHit_id->at(j)%1000000)/100000;
			FluxHit_subsubdetector_ID=((fFluxHit_id->at(j)%1000000)%100000)/10000;
			FluxHit_component_ID=fFluxHit_id->at(j)%10000;        
			pf_flux=TMath::Sqrt(TMath::Power(fFluxHit_Px->at(j),2)+TMath::Power(fFluxHit_Py->at(j),2)+TMath::Power(fFluxHit_Pz->at(j),2));
			r_flux = TMath::Sqrt(TMath::Power(fFluxHit_Avg_x->at(j),2)+TMath::Power(fFluxHit_Avg_y->at(j),2));
			E_flux=fFluxHit_trackE->at(j);

			if (FluxHit_detector_ID==3 && FluxHit_subdetector_ID == 1 && FluxHit_subsubdetector_ID == 1 && Xjb>0.3 && Q2>2 /*&& theta_gen>22 && theta_gen<35 */ &&  pf_flux>=1800 &&  pf_flux<=12000 ){
				th=TMath::ATan((r_flux*0.1)/(fFluxHit_Avg_z->at(j)*0.1 - targ_offset))*DEG;
				ECALDet_mom_distr->Fill(fFluxHit_vz->at(j));
				pf_satisfy+=pf_flux;
				flux_E= fFluxHit_totEdep->at(j);
				flux_size =fFluxHit_id->size(); 
				//cout<<"event="<<i<<"pf="<<pf_flux*0.001<<"theta="<<th<<endl; 
				Gen_theta_vs_p_EC->Fill(theta_gen,p_gen*0.001,event_rate);
				Gen_theta_vs_p_EC_accep->Fill(r_flux,pf_flux*0.001,event_rate);
				//save all the hits per event
				single_event.pid=fFluxHit_pid->at(j);
				single_event.x=fFluxHit_Avg_x->at(j)*0.1;
				single_event.y=fFluxHit_Avg_y->at(j)*0.1;
				X=fFluxHit_Avg_x->at(j)*0.1;
				Y=fFluxHit_Avg_y->at(j)*0.1;
				single_event.pf=pf_flux;
				single_event.Q2=Q2;
				single_event.blockID =FluxHit_component_ID;    
				event_list.push_back(single_event);          
				if(fFluxHit_Pz->at(j)>0){  
					if (E_flux >= E_flux_max){
						E_flux_max=E_flux;
					}
				}
				//event_list.push_back(single_event);
				Fillecal=kTRUE;
				processWindow=kTRUE;
				break;
			}
		}

		if (Fillecal){
			TGEMC_EC[n]->GetEntry(i);
			TGEMC_ECPS[n]->GetEntry(i);
			totEdep_shower=0;
			totEdep_preshower=0;
			EC_averagex=0.0;
			EC_averagey=0.0;
			EC_localx=0.0;
			EC_localy=0.0;

			for (int j = 0; j<fECHit_id->size(); j++){
				ECHit_detector_ID=fECHit_id->at(j)/1000000;
				ECHit_subdetector_ID=(fECHit_id->at(j)%1000000)/100000;
				ECHit_subsubdetector_ID=((fECHit_id->at(j)%1000000)%100000)/10000;
				ECHit_component_ID=fECHit_id->at(j)%10000;
				if (ECHit_detector_ID==3 && ECHit_subdetector_ID == 1 && ECHit_subsubdetector_ID == 0){//shower 
					EC_averagex=fECHit_Avg_x->at(j)*0.1;
					EC_averagey=fECHit_Avg_y->at(j)*0.1;
					EC_localx=fECHit_Avg_lx->at(j)*0.1;
					EC_localy=fECHit_Avg_ly->at(j)*0.1;
					totEdep_shower+=fECHit_totEend->at(j);
					if (ecalShMap.count(ECHit_component_ID)){
						ecalShMap[ECHit_component_ID]+=fECHit_totEend->at(j);// MeV
					} else {
						ecalShMap[ECHit_component_ID]=fECHit_totEend->at(j);// MeV
					}
				}
			}// EC loop
			for (int k = 0; k<fECPSHit_id->size(); k++){
				ECPSHit_detector_ID=fECPSHit_id->at(k)/1000000;
				ECPSHit_subdetector_ID=(fECPSHit_id->at(k)%1000000)/100000;
				ECPSHit_subsubdetector_ID=((fECPSHit_id->at(k)%1000000)%100000)/10000;
				ECPSHit_component_ID=fECPSHit_id->at(k)%10000;

				if (ECPSHit_detector_ID==3 && ECPSHit_subdetector_ID == 1 && ECPSHit_subsubdetector_ID ==1 ){//shower 
					totEdep_preshower+=fECPSHit_totEend->at(k);
					if (ecalPSMap.count(ECPSHit_component_ID)){
						ecalPSMap[ECPSHit_component_ID]+=fECPSHit_totEend->at(k);// MeV
					} else {
						ecalPSMap[ECPSHit_component_ID]=fECPSHit_totEend->at(k);// MeV
					}
				}
			}//end PSEChit loop

			htotEdep_ec_shower->Fill(totEdep_shower);
			htotEdep_ec_preshower->Fill(totEdep_preshower);

		} //fillcall
		//process hits in this time window
		if (processWindow){  

			edep_6p1_sum=0;
			for (Int_t c=0;c<7;c++){//init cluster array
				cluster_edep[c]=0;
				cluster_edep_blockid[c]=0;
				cluster_PS_edepbkg[c]=0;
				cluster_edep_bkg[c]=0;
			}
			for (Int_t c=0;c<30;c++){
				edep_6p1_max1[c] = 0;
				radius_6p1_max1[c]={0};
				edep_6p1_preshower_max1[c] = 0;
				radius_6p1_preshower_max1[c]={0};	
				edep_6p1_x[c]={0};
				edep_6p1_y[c]={0};
				edep_6p1_x_preshower[c]={0};
				edep_6p1_y_preshower[c]={0};	
			}   
			//edep_6p1_max = 0;
			octantno = -1;

			for (std::map<Int_t,Double_t>::iterator it = ecalShMap.begin(); it != ecalShMap.end(); ++it ){
				ecal_blockID = it->first;	
				edep_1block_sum = it->second;
				vec_ecalBlock = GetECALBlock_coord(ecal_blockID);
				ecalBlock_x = vec_ecalBlock.Px();
				ecalBlock_y = vec_ecalBlock.Py();
				radius_6p1=TMath::Sqrt(TMath::Power(ecalBlock_x,2)+TMath::Power(ecalBlock_y,2));
				octantno = 0;//Int_t(phi/12);//octant number 0 to 29 for 30 octants
				if (octantno>29)
					octantno=29;//for phi is 360       
				ecal_hitblockID=ecal_blockID;
				//find 6+1 blocks based on hit position x and y
				block_count = GetECALCluser(ecal_hitblockID,cluster_edep_blockid);
				cluster_edep[0]=ecalShMap[cluster_edep_blockid[0]];
				edep_6p1_sum = cluster_edep[0];
				for(Int_t j=0;j<block_count;j++){
					cluster_edep[j+1]=ecalShMap[cluster_edep_blockid[j+1]];
					edep_6p1_sum +=cluster_edep[j+1];
				}
				if (edep_6p1_sum >= edep_6p1_max1[octantno]){
					ecalBlock_x0=ecalBlock_x;
					ecalBlock_y0=ecalBlock_y;
					ecal_blockID0 = ecal_blockID;
					edep_6p1_x[octantno]=ecalBlock_x;
					edep_6p1_y[octantno]=ecalBlock_y;
					radius_6p1_max1[octantno]=TMath::Sqrt(TMath::Power(ecalBlock_x,2)+TMath::Power(ecalBlock_y,2));
					edep_6p1_max1[octantno]=edep_6p1_sum;
				}
			}// endl shower

			vec_ecalBlock_preshower = GetECALBlock_coord(ecal_blockID);
			ecalBlock_x_preshower = vec_ecalBlock_preshower.Px();
			ecalBlock_y_preshower = vec_ecalBlock_preshower.Py();
			radius_6p1_preshower=TMath::Sqrt(TMath::Power(ecalBlock_x_preshower,2)+TMath::Power(ecalBlock_y_preshower,2));
			octantno = 0;//Int_t(phi/12);//octant number 0 to 29 for 30 octants
			if (octantno>29)
				octantno=29;//for phi is 360           

			ecal_hitblockID_psh = ecal_blockID0;
			block_count_PS = GetECALCluser(ecal_hitblockID_psh,cluster_edep_PS_blockid);
			cluster_PS_edep[0]=ecalPSMap[cluster_edep_PS_blockid[0]];
			edep_6p1_PS_sum = cluster_PS_edep[0];
			for(Int_t j=0;j<block_count_PS;j++){
				cluster_PS_edep[j+1]=ecalPSMap[cluster_edep_PS_blockid[j+1]];
				edep_6p1_PS_sum +=cluster_PS_edep[j+1];
			}

			//	  if (edep_6p1_PS_sum >= edep_6p1_preshower_max1[octantno]){
			ecal_blockID0_preshower = ecal_blockID_psh;
			edep_6p1_x_preshower[octantno]=ecalBlock_x_preshower;
			edep_6p1_y_preshower[octantno]=ecalBlock_y_preshower;
			edep_6p1_preshower_max1[octantno]=edep_6p1_PS_sum;
			for (Int_t oc=0;oc<__OCTANTS;oc++){//access the octant level cluster sums to check trigger condition
				//shower clusters
				if (edep_6p1_max1[oc]>0){//check only if max is greater than or equal zero
					//Ri=GetRadiusIndex(radius_6p1_max1[oc]);
					if(edep_6p1_max1[oc] >= GetThreshold6p1(radius_6p1_max1[oc])){//use r to get radius at the GEM && trigger_state_PS[oc]
						eventid=i;
						trigger_state_1[oc]=kTRUE;//octant level trigger
						break;
					} else{
						trigger_state_1[oc]=kFALSE;	 
					} 
				} 
			}//end of oc loop 
			if (!bDisableBkg){
				edep_total_sum=0;
				rndTimeWindow = (Int_t)(r3->Uniform()*(total_ecalShMapbkg.size()-1));
				edep_total_sum=0;
				edep_6p1_sum_bkg=0;  
				edep_6p1_presum_bkg=0;  
				edep_6p1_max1_bkg=0;
				edep_6p1_premax1_bkg=0;
				getBkgEdepPS(rndTimeWindow);
				getBkgEdepSh(rndTimeWindow);
				int showermap_N=0;
				ecal_blockID_bkg = ecal_blockID0; 
				block_count_bkg = GetECALCluser(ecal_blockID_bkg,cluster_edep_blockid_bkg);
				cluster_edep_bkg[0]=ecalShMapbkg[cluster_edep_blockid_bkg[0]];
				edep_6p1_sum_bkg = cluster_edep_bkg[0];
				for(Int_t j=0;j<block_count_bkg;j++){
					cluster_edep_bkg[j+1]=ecalShMapbkg[cluster_edep_blockid_bkg[j+1]];
					edep_6p1_sum_bkg +=cluster_edep_bkg[j+1];
				}
				edep_6p1_max1_bkg=edep_6p1_sum_bkg+edep_6p1_max1[octantno];
				cluster_PS_edepbkg[0]=ecalPSMapbkg[cluster_edep_blockid_bkg[0]];
				edep_6p1_presum_bkg = cluster_PS_edepbkg[0];
				for(Int_t j=0;j<block_count_bkg;j++){
					cluster_PS_edepbkg[j+1]=ecalPSMapbkg[cluster_edep_blockid_bkg[j+1]];
					edep_6p1_presum_bkg +=cluster_PS_edepbkg[j+1];
				}

				edep_6p1_premax1_bkg=edep_6p1_presum_bkg+edep_6p1_preshower_max1[octantno];
			}

			// cout<<"event_pass="<<eventid<<endl;    
			Int_t hit_pid_index;
			Int_t hit_pid_index2;
			Double_t hit_phi;
			Double_t event_p;
			Double_t event_E;
			Double_t event_Ek;
			Double_t DeltaX; 
			Double_t DeltaY;
			Double_t DeltaR;
			Int_t hit_octant;
			Int_t Ri;
			for(Int_t e=0;e<event_list.size();e++){

				event_p = event_list[e].pf/1000;
				Ri=GetRadiusIndex(TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2)));
				dE=1.0; 
				DeltaRrDiff = TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2));
				hit_pid_index =  pidmap[event_list[e].pid];//pid index starts from 1 but array index starts from 0
				hit_octant = 0;//since there will only be one event per window I disabled sector id//Int_t(phi/12);
				if (hit_octant>29)
					hit_octant=29;     
				if (hit_pid_index>0 && hit_pid_index<=5 && (event_p <= 12)){
					if (hit_pid_index==1 || hit_pid_index==4)//e+-
						hit_pid_index2 = 0;
					if (hit_pid_index==2 || hit_pid_index==5)//any pion
						hit_pid_index2 = 1;
					if (hit_pid_index==3 ||hit_pid_index==1 ||hit_pid_index==4){//gamma
						if (event_list[e].pid==111)//no simpid in the electron+bkg merged root file
							hit_pid_index2 = 2;//photons from pi0 simulation
						else
							hit_pid_index2 = 3;//photons from other sources
					}


					histo_pid[0][hit_pid_index2]->Fill(event_list[e].pf);//momentum distribution

					if (hit_pid_index2 == 0 && event_list[e].pf< 1000)//bkg e+/- less than 1 GeV
						histo_pid2[0][0][8]->Fill(event_list[e].pf);//momentum distribution
					if (hit_pid_index == 1 && event_list[e].pf >= 1000){//high energy electrons
						histo_pid2[0][3][8]->Fill(event_list[e].pf);//momentum distribution        
						histo_pid2[0][3][Ri]->Fill(event_list[e].pf);
					}
					Gen_theta_vs_p_EC_accep->Fill(radius_6p1_max1[octantno],event_list[e].pf*0.001,event_rate);
					if (event_list[e].pf >= 1800){//any pion
						histo_pid2[0][1][8]->Fill(event_list[e].pf);//momentum distribution
						histo_pid2[0][1][Ri]->Fill(event_list[e].pf);
						histo_E_Sh_radius_threshold_nobkg[7]->Fill(E_flux_max,event_rate); 
						if(110<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=130  ){
							histo_edep_Sh_radius_threshold[0]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[0]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[0]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[0]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[0]->Fill(E_flux_max,event_rate); 

						}else if(130<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=150  ){
							histo_edep_Sh_radius_threshold[1]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[1]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[1]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[1]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[1]->Fill(E_flux_max,event_rate); 

						}else if(150<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=170  ){
							histo_edep_Sh_radius_threshold[2]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[2]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[2]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[2]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[2]->Fill(E_flux_max,event_rate); 

						}else if(170<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=190){
							histo_edep_Sh_radius_threshold[3]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[3]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[3]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[3]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[3]->Fill(E_flux_max,event_rate); 

						}else if(190<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=210 ){
							histo_edep_Sh_radius_threshold[4]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[4]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[4]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[4]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[4]->Fill(E_flux_max,event_rate); 

						}else if(210<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=230 ){
							histo_edep_Sh_radius_threshold[5]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[5]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[5]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[5]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[5]->Fill(E_flux_max,event_rate); 

						}else if(230<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=270 ){
							histo_edep_Sh_radius_threshold[6]->Fill(edep_6p1_max1_bkg,event_rate); 
							histo_edep_Sh_radius_nothreshold_nobkg[6]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 
							histo_edep_Sh_radius_threshold_nobkg[6]->Fill(edep_6p1_max1[octantno],event_rate); 
							histo_edep_Sh_radius_nothreshold[6]->Fill(edep_6p1_premax1_bkg,event_rate); 
							histo_E_Sh_radius_threshold_nobkg[6]->Fill(E_flux_max,event_rate); 

						}

					}
					if (hit_pid_index == 3)//any gamma
						histo_pid2[0][2][8]->Fill(event_list[e].pf);//momentum distribution


					if(edep_6p1_max1_bkg>GetThreshold6p1(radius_6p1_max1[octantno])){
						DeltaR = (TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2)) - TMath::Sqrt(TMath::Power(edep_6p1_x[hit_octant],2)+TMath::Power(edep_6p1_y[hit_octant],2)));
						//cout<<"DeltaR="<<DeltaR<<endl;
						histo_pid[1][hit_pid_index2]->Fill(event_list[e].pf);//update the total hits
						if (hit_pid_index2 == 0 && event_list[e].pf < 1000)//bkg e+/- less than 1 GeV
							histo_pid2[1][0][8]->Fill(event_list[e].pf);//momentum distribution
						if (hit_pid_index == 1 && event_list[e].pf >= 1000){//high energy electrons
							histo_pid2[1][3][8]->Fill(event_list[e].pf);//momentum distribution
							histo_pid2[1][3][Ri]->Fill(event_list[e].pf);
						}
						if (event_list[e].pf >= 1000){//any pion
							histo_pid2_cut[0][1][8]->Fill(event_list[e].pf);//momentum distribution
							histo_pid2_cut[0][1][Ri]->Fill(event_list[e].pf);
							if(110<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=130  ){
								histo_edep_Sh_radius_threshold_cut[0]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[0]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(130<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=150  ){
								histo_edep_Sh_radius_threshold_cut[1]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[1]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(150<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=170  ){
								histo_edep_Sh_radius_threshold_cut[2]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[2]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(170<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=190){
								histo_edep_Sh_radius_threshold_cut[3]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[3]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(190<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=210 ){
								histo_edep_Sh_radius_threshold_cut[4]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[4]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(210<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=230 ){
								histo_edep_Sh_radius_threshold_cut[5]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[5]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}else if(230<radius_6p1_max1[octantno] && radius_6p1_max1[octantno]<=270 ){
								histo_edep_Sh_radius_threshold_cut[6]->Fill(edep_6p1_max1[hit_octant],event_rate); 
								histo_edep_Sh_radius_nothreshold_cut[6]->Fill(edep_6p1_preshower_max1[hit_octant],event_rate); 

							}
						}
						if (hit_pid_index2 == 1){//any pion
							histo_pid2[1][1][8]->Fill(event_list[e].pf);//momentum distribution
							histo_pid2[1][1][Ri]->Fill(event_list[e].pf);
						}
						if (hit_pid_index == 3)//any gamma
							histo_pid2[1][2][8]->Fill(event_list[e].pf);//momentum distribution	
					}    
				} 	     


			}//end eventlist
			processWindow=kFALSE;//reset till next window is filled
			ecalPSMap.clear();
			ecalShMap.clear();
			event_list.clear();
			Fillecal=kFALSE;

		}//endProcesswindow
		}//end entries
	}
	cout<<"N_toal="<<N_total<<"  "<<"N_loss="<<N_loss<<"total_rate="<<total_rate_eDIS<<endl;
	//efficiency plots
	Double_t p[9][100] = { {0} },p_2[9][100] = { {0} };
	Double_t eff[9][100] = { {0} },eff_2[9][100] = { {0} };
	Double_t efferr[9][100] = { {0} },efferr_2[9][100] = { {0} };
	int cnt[9] = {0};
	int cnt_2[9] = {0};
	TGraphErrors * geffi_total_pion[9];
	//  printf("for Pions \n Momentum \t Efficiency \t Error \n");
	for (Int_t Ri=0;Ri<9;Ri++){
		for ( int i = 1; i <= histo_pid2[1][1][Ri]->GetNbinsX(); i++ ){
			if ( histo_pid2[1][1][Ri]->GetBinContent(i) <= histo_pid2[0][1][Ri]->GetBinContent(i) && histo_pid2[0][1][Ri]->GetBinContent(i) > 20 ){
				p[Ri][cnt[Ri]] = histo_pid2[1][1][Ri]->GetBinCenter(i)/1000;
				eff[Ri][cnt[Ri]] = histo_pid2[1][1][Ri]->Integral(i,i) / histo_pid2[0][1][Ri]->Integral(i,i);
				efferr[Ri][cnt[Ri]] = TMath::Sqrt(eff[Ri][cnt[Ri]] * (1 - eff[Ri][cnt[Ri]]) / histo_pid2[0][1][Ri]->GetBinContent(i));
				printf("%5i %13.3f %13.3f %13.3f \n",Ri,p[Ri][cnt[Ri]],eff[Ri][cnt[Ri]],efferr[Ri][cnt[Ri]]);
				cnt[Ri]++;
			}
		}
		geffi_total_pion[Ri] = new TGraphErrors(cnt[Ri], p[Ri], eff[Ri], 0, efferr[Ri]);
		geffi_total_pion[Ri]->SetName(Form("geffi_total_pion_%d",Ri));
		geffi_total_pion[Ri]->SetTitle(Form("Pion Efficiency at %s; Momentum (GeV);Efficiency",sradii[Ri].Data()));

		//cnt = 0;
	}

	//save trigger efficiencies

	printf("Double_t R_low[8] = {");
	for(Int_t i=0;i<8;i++)
		printf("%4.1f,",trig_low_R[i]);
	printf("}; \n");
	printf("Double_t R_up[8] = {");
	for(Int_t i=0;i<8;i++)
		printf("%4.1f,",trig_high_R[i]);
	printf("}; \n");

	Double_t mom_bin_width = 0.28;
	printf("Double_t P_low[%d] = {",cnt[8]);
	for(Int_t i=0;i<cnt[8];i++)
		printf("%4.1f,",p[8][i]-mom_bin_width/2);
	printf("}; \n");
	printf("Double_t P_high[%d] = {",cnt[8]);
	for(Int_t i=0;i<cnt[8];i++)
		printf("%4.1f,",p[8][i]+mom_bin_width/2);
	printf("}; \n");


	printf("Double_t trig_pi_eff[9][%d] = { ",cnt[8]);
	for(Int_t Ri=0;Ri<9;Ri++){
		printf("{ ");
		if (cnt[8] != cnt[Ri]) {
			printf("momentum array length mistmatch, [8] and [%d]",cnt[Ri]);
			break;
		}
		for(Int_t Pi=0;Pi<cnt[Ri];Pi++)
			printf("%4.3f,", eff[Ri][Pi]);
		printf("},");
	}
	printf("}; \n");

	printf("Double_t trig_pi_eff_err[9][%d] = { ",cnt[8]);
	for(Int_t Ri=0;Ri<9;Ri++){
		printf("{ ");
		if (cnt[8] != cnt[Ri]) {
			printf("momentum array length mistmatch, [8] and [%d]",cnt[Ri]);
			break;
		}
		for(Int_t Pi=0;Pi<cnt[Ri];Pi++)
			printf("%4.3f,", efferr[Ri][Pi]);
		printf("},");
	}
	printf("}; \n");

	//cnt = 0;
	printf("for Electrons \n Momentum \t Efficiency \t Error \n");
	TGraphErrors * geffi_total_electron[9];
	for (Int_t Ri=0;Ri<9;Ri++){
		for ( int i = 1; i <= histo_pid2[1][3][Ri]->GetNbinsX(); i++ ){
			if ( histo_pid2[1][3][Ri]->GetBinContent(i) <= histo_pid2[0][3][Ri]->GetBinContent(i) && histo_pid2[0][3][Ri]->GetBinContent(i) > 20){
				p_2[Ri][cnt_2[Ri]] = histo_pid2[1][3][Ri]->GetBinCenter(i)/1000;
				eff_2[Ri][cnt_2[Ri]] = histo_pid2[1][3][Ri]->GetBinContent(i)/ histo_pid2[0][3][Ri]->GetBinContent(i);
				efferr_2[Ri][cnt_2[Ri]] = TMath::Sqrt(eff_2[Ri][cnt_2[Ri]] * (1 - eff_2[Ri][cnt_2[Ri]]) / histo_pid2[0][3][Ri]->GetBinContent(i));
				printf("%5i %13.3f %13.3f %13.3f \n",Ri,p_2[Ri][cnt_2[Ri]],eff_2[Ri][cnt_2[Ri]],efferr_2[Ri][cnt_2[Ri]]);
				cnt_2[Ri]++;
			}
		}
		geffi_total_electron[Ri] = new TGraphErrors(cnt_2[Ri], p_2[Ri], eff_2[Ri], 0, efferr_2[Ri]);
		geffi_total_electron[Ri]->SetName(Form("geffi_total_electron_%d",Ri));
		geffi_total_electron[Ri]->SetTitle(Form("Electron Efficiency at %s; Momentum (GeV);Efficiency",sradii[Ri].Data()));

		//cnt = 0;
	}

	printf("Double_t trig_ele_eff[9][%d] = { ",cnt_2[8]);
	for(Int_t Ri=0;Ri<9;Ri++){
		printf("{ ");
		if (cnt_2[8] != cnt_2[Ri]) {
			printf("momentum array length mistmatch, [8] and [%d]",cnt_2[Ri]);
			break;
		}
		for(Int_t Pi=0;Pi<cnt_2[Ri];Pi++)
			printf("%4.3f,", eff_2[Ri][Pi]);
		printf("},");
	}
	printf("}; \n");

	printf("Double_t trig_ele_eff_err[9][%d] = { ",cnt_2[8]);
	for(Int_t Ri=0;Ri<9;Ri++){
		printf("{ ");
		if (cnt_2[8] != cnt_2[Ri]) {
			printf("momentum array length mistmatch, [8] and [%d]",cnt_2[Ri]);
			break;
		}
		for(Int_t Pi=0;Pi<cnt_2[Ri];Pi++)
			printf("%4.3f,", efferr_2[Ri][Pi]);
		printf("},");
	}
	printf("}; \n");


	TCanvas *canvas_eff_pion = new TCanvas("canvas_eff_pion","canvas_eff_pion", 1000, 700);
	geffi_total_pion[8]->Draw("AP*");
	geffi_total_pion[8]->SetTitle(Form("Pion Efficiency ;Momentum (GeV);Efficiency"));
	gPad->Update();
	TCanvas *canvas_eff_electron = new TCanvas("canvas_eff_electron","canvas_eff_electron", 1000, 700);
	geffi_total_electron[8]->Draw("AP*");
	geffi_total_electron[8]->SetTitle(Form("Electron Efficiency ;Momentum (GeV);Efficiency"));
	gPad->Update();
	TMultiGraph *mg_eff_electron = new TMultiGraph();

	TLegend *leg1 = new TLegend(0.71,0.31,0.89,0.59,NULL,"brNDC");//0.6,0.78,0.9,0.90//numbers obtained by first moving the legend to a prefered location and then looking at the inspect pop-up window
	leg1->SetTextFont(62);
	//leg1->SetTextSize(0.03);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(1001);

	for(Int_t i=0;i<7;i++){
		geffi_total_electron[i]->SetLineStyle(i+1);
		geffi_total_electron[i]->SetLineColor(1+i);//high contrast colors    
		geffi_total_electron[i]->SetLineWidth(2);//high contrast colors
		mg_eff_electron->Add(geffi_total_electron[i],"lp");    
		leg1->AddEntry(geffi_total_electron[i],Form("%s",sradii[i].Data()),"l");
	}
	mg_eff_electron->SetName("mg_eff_electron");


	TCanvas *canvas_mg_eff_electron = new TCanvas("canvas_mg_eff_electron","canvas_mg_eff_electron", 1000, 700);
	canvas_mg_eff_electron->SetGrid();
	mg_eff_electron->SetTitle("Electron Efficiency ;Momentum (GeV);Efficiency");
	mg_eff_electron->SetMaximum(1.05);
	mg_eff_electron->SetMinimum(0.0);
	//mg_eff_electron->GetXaxis()->SetRangeUser(0.98,8.2);
	mg_eff_electron->Draw("A");
	leg1->Draw();
	canvas_mg_eff_electron->RedrawAxis();
	gPad->Update();
	TMultiGraph *mg_eff_pion = new TMultiGraph();

	TLegend *leg2 = new TLegend(0.11,0.71,0.30,0.89,NULL,"brNDC");//0.6,0.78,0.9,0.90
	leg2->SetTextFont(62);
	//leg2->SetTextSize(0.03);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(1001);

	for(Int_t i=0;i<8;i++){
		geffi_total_pion[i]->SetLineStyle(i+1);
		geffi_total_pion[i]->SetLineColor(1+i);//high contrast colors
		geffi_total_pion[i]->SetLineWidth(2);//high contrast colors
		mg_eff_pion->Add(geffi_total_pion[i],"lp");
		leg2->AddEntry(geffi_total_pion[i],Form("%s",sradii[i].Data()),"l");
	}
	mg_eff_pion->SetName("mg_eff_pion");

	TCanvas *canvas_mg_eff_pion = new TCanvas("canvas_mg_eff_pion","canvas_mg_eff_pion", 1000, 700);
	mg_eff_pion->SetTitle("Pion Efficiency ;Momentum (GeV);Efficiency");
	mg_eff_pion->Draw("A");
	leg2->Draw();
	gPad->Update();

	if (kSaveCanvas){
		canvas_eff_electron->SaveAs("/home/tianye/solid/plots/canvas_eff_electron_R_Track_PVDIS.png");
		canvas_mg_eff_electron->SaveAs("/home/tianye/solid/plots/canvas_eff_radii_electron_R_Track_PVDIS.png");
	}



	/*****************************************************************************************/
	TCanvas *shower_6p1_time = new TCanvas("shower_6p1_time","shower_6p1_time",1000,1000);
	Shower_vs_timewindows_hist->Draw("COLZ");
	TLegend *leg4 = new TLegend(0.71,0.31,0.89,0.59,NULL,"brNDC");//0.6,0.78,0.9,0.90//numbers obtained by first moving the legend to a prefered location and then looking at the inspect pop-up window
	leg4->SetTextFont(62);
	//leg1->SetTextSize(0.03);
	leg4->SetFillColor(0);
	leg4->SetFillStyle(1001);
	TCanvas * canvas_6plus1_threshold = new TCanvas("canvas_6plus1_threshold","canvas_6plus1_threshold",1400,400);
	//canvas_6plus1_threshold->Divide(2,3);
	TLine *line[6];
	histo_edep_Sh_radius_threshold[0]->SetLineWidth(2);
	histo_edep_Sh_radius_threshold[0]->SetLineColor(1); 
	leg4->AddEntry(histo_edep_Sh_radius_threshold[0],Form("%s",sradii[0].Data()),"l");
	histo_edep_Sh_radius_threshold[0]->Draw();
	for(int i=1;i<4;i++){
		//canvas_6plus1_threshold->cd(i+1);
		line[i] = new TLine(trig_thresh_6p1[i],0,trig_thresh_6p1[i],30000);
		//histo_pid2[1][3][i]->Draw();
		histo_edep_Sh_radius_threshold[i]->SetLineWidth(2);
		histo_edep_Sh_radius_threshold[i]->SetLineColor(i+1); 
		leg4->AddEntry(histo_edep_Sh_radius_threshold[i],Form("%s",sradii[i].Data()),"l");
		histo_edep_Sh_radius_threshold[i]->Draw("same");
		line[i]->SetLineColor(kRed);
		// line[i]->Draw("same");
	} 
	TLegend *leg5 = new TLegend(0.71,0.31,0.89,0.59,NULL,"brNDC");//0.6,0.78,0.9,0.90//numbers obtained by first moving the legend to a prefered location and then looking at the inspect pop-up window
	leg5->SetTextFont(62);
	//leg1->SetTextSize(0.03);
	leg5->SetFillColor(0);
	leg5->SetFillStyle(1001);
	TCanvas * canvas_6plus1_nothreshold = new TCanvas("canvas_6plus1_nothreshold","canvas_6plus1_bothreshold",1400,400);
	//canvas_6plus1_threshold->Divide(2,3);
	histo_edep_Sh_radius_nothreshold[0]->SetLineWidth(2);
	histo_edep_Sh_radius_nothreshold[0]->SetLineColor(1); 
	leg5->AddEntry(histo_edep_Sh_radius_nothreshold[0],Form("%s",sradii[0].Data()),"l");
	histo_edep_Sh_radius_nothreshold[0]->Draw();
	for(int i=1;i<4;i++){
		//canvas_6plus1_threshold->cd(i+1);
		line[i] = new TLine(trig_thresh_6p1[i],0,trig_thresh_6p1[i],30000);
		//histo_pid2[1][3][i]->Draw();
		histo_edep_Sh_radius_nothreshold[i]->SetLineWidth(2);
		histo_edep_Sh_radius_nothreshold[i]->SetLineColor(i+1); 
		leg5->AddEntry(histo_edep_Sh_radius_nothreshold[i],Form("%s",sradii[i].Data()),"l");
		histo_edep_Sh_radius_nothreshold[i]->Draw("same");
		line[i]->SetLineColor(kRed);
		// line[i]->Draw("same");
	}         
	TCanvas * canvas_kinematics_ecal_front = new TCanvas("canvas_kinematics_ecal_front","canvas_kinematics_ecal_front",1400,400);
	canvas_kinematics_ecal_front->Divide(3,1);
	canvas_kinematics_ecal_front->cd(1);
	Gen_mom_distr->Draw();
	canvas_kinematics_ecal_front->cd(2);
	Gen_theta_distr->Draw();
	canvas_kinematics_ecal_front->cd(3);
	ECALDet_mom_distr->Draw();

	TCanvas * canvas_DeltaX = new TCanvas("canvas_DeltaXY","canvas_DeltaXY",1000,400);
	canvas_DeltaX->Divide(2,1);   
	canvas_DeltaX->cd(1);
	BlockX_hist->Draw();
	canvas_DeltaX->cd(2);
	BlockY_hist->Draw();      

	//For PS+Sh 
	TCanvas * canvas_run_level_e_res_PS_shower = new TCanvas("canvas_run_level_e_res_PS_shower","canvas_run_level_e_res_PS_shower",1000,1400);
	canvas_run_level_e_res_PS_shower->Divide(2,2);
	//canvas_run_level_e_res_PS_shower->Divide(2,2);
	canvas_run_level_e_res_PS_shower->cd(1);
	EoverPElec_total_ecal_fullsum->Draw("colz");
	canvas_run_level_e_res_PS_shower->cd(2);
	EoverPElec_total_ecal_6p1sum->Draw("colz");
	canvas_run_level_e_res_PS_shower->cd(3);
	TProfile * EoverPElec_total_ecal_fullsum_x = EoverPElec_total_ecal_fullsum->ProfileX();
	EoverPElec_total_ecal_fullsum_x->Draw();
	EoverPElec_total_ecal_fullsum_x->Fit("pol0");
	canvas_run_level_e_res_PS_shower->cd(4);
	TProfile * EoverPElec_total_ecal_6p1sum_x = EoverPElec_total_ecal_6p1sum->ProfileX();
	EoverPElec_total_ecal_6p1sum_x->Draw();
	EoverPElec_total_ecal_6p1sum_x->Fit("pol0");

	TCanvas * canvas_run_level_e_res_PS_shower_final = new TCanvas("canvas_run_level_e_res_PS_shower_final","canvas_run_level_e_res_PS_shower_final",1000,700);
	canvas_run_level_e_res_PS_shower_final->Divide(2,1);

	canvas_run_level_e_res_PS_shower_final->cd(1);
	TGraph * EoverpRMS_total_fullsum = GetRMSReolution(EoverPElec_total_ecal_fullsum);//0.94,1.02 for limited angle simulation //0.96,1.02 for all angle range
	EoverpRMS_total_fullsum->SetTitle("ECAL PS+Sh Total Energy Resolution VS p; Momentum (GeV/c); #DeltaE/E");
	TF1 * total_f1 = new TF1("total_Eres1", "sqrt(pow([0]/sqrt(x),2)+pow([1],2)+pow([2]/x,2))", 2, 7);//1,10 is this in GeV?
	total_f1->SetLineColor(kRed);
	EoverpRMS_total_fullsum->Fit(total_f1, "M");
	canvas_run_level_e_res_PS_shower_final->cd(2);
	TGraph * EoverpRMS_total_6p1sum = GetRMSReolution(EoverPElec_total_ecal_fullsum);//combined PS+sh 0.84,0.98//if actual edep on material used 0.93,0.995
	EoverpRMS_total_6p1sum->SetTitle("ECAL PS+Sh 6+1 Energy Resolution VS p; Momentum (GeV/c); #DeltaE/E");
	TF1 * total_f2 = new TF1("total_Eres2", "sqrt(pow([0]/sqrt(x),2)+pow([1],2)+pow([2]/x,2))", 2, 7);
	total_f2->SetLineColor(kRed);
	EoverpRMS_total_6p1sum->Fit(total_f2, "M");
	if (kSaveRootFile){
		rootfile->Write();
	}      
	watch.Stop();
	//printf("Total Time %3.4f min \n Completed, Exit now :) \n",watch.RealTime()/60);
	theApp.Run();

	if (kSaveRootFile){
		rootfile->Close();
	}
	return(1);
}

void GetFluxEntryByBranch(Long64_t local){
	bFluxHit_id->GetEntry(local);  
	bFluxHit_n->GetEntry(local);
	bFluxHit_pid->GetEntry(local);
	bFluxHit_tid->GetEntry(local);  
	bFluxHit_trackE->GetEntry(local); 
	bFluxHit_Avg_x->GetEntry(local);  
	bFluxHit_Avg_y->GetEntry(local);  
	bFluxHit_Avg_z->GetEntry(local);
	bFluxHit_Px->GetEntry(local);  
	bFluxHit_Py->GetEntry(local);  
	bFluxHit_Pz->GetEntry(local);  
};

void SetECBranchAddresses(){
	TGEMC_EC->SetBranchAddress("id",&fECHit_id,&bECHit_id);
	TGEMC_EC->SetBranchAddress("totEdep",&fECHit_totEdep,&bECHit_totEdep);
	TGEMC_EC->SetBranchAddress("avg_x",&fECHit_Avg_x,&bECHit_Avg_x);
	TGEMC_EC->SetBranchAddress("avg_y",&fECHit_Avg_y,&bECHit_Avg_y);
	TGEMC_EC->SetBranchAddress("avg_z",&fECHit_Avg_z,&bECHit_Avg_z);
};

void GetECEntryByBranch(Long64_t local){
	bECHit_id->GetEntry(local);  
	bECHit_totEdep->GetEntry(local);  
	bECHit_Avg_x->GetEntry(local);  
	bECHit_Avg_y->GetEntry(local);  
	bECHit_Avg_z->GetEntry(local);  
};

void LoadEC_map(TString map_file){
	//std. map file in the svn repo  ../layout/map_FAEC_ANL_20130628.txt
	ifstream INPUT_file;
	INPUT_file.open(map_file.Data());
	if(!INPUT_file){ 
		printf("ERROR!!! Can't open %s \n",map_file.Data());	
		exit(1);
	}



	//load x,y locations and block id
	Double_t total_module=0;
	Int_t counter_id=0;
	for(int i=0;i<54;i++){ //54 rows
		INPUT_file>>num_module_in_row[i];
		//std::cout<<num_module_in_row[i]<<" ";
		num_module_in_row[i]=num_module_in_row[i]-1;  //first one is y coordinate
		total_module+=num_module_in_row[i];
		Double_t tmp_y;
		INPUT_file>>tmp_y;
		y_bak[i]=tmp_y;       //make a backup in order to judge which row a certain particle hits the EC
		for(Int_t j=0;j<num_module_in_row[i];j++){
			INPUT_file>>x[i][j];
			y[i][j]=tmp_y;     // in each row , y coordinate is the same
			counter_id++;
			id[i][j]=counter_id;
			//update the pair array
			block_map[counter_id].first = i;
			block_map[counter_id].second = j;
			//sector
			TVector2 vec(x[i][j],y[i][j]);
			Double_t phi_module=vec.Phi();
			for(int k=0;k<30;k++){   //30 sectors
				if(phi_module>=k*12.0/180.0*3.141592 && phi_module<(k+1.0)*12.0/180.0*3.141592){  //sector k
					sector[i][j]=k+1;
				}
			}

		}
	}
	printf("Total ECAL Modules loaded %i \n",counter_id);

	Int_t total_module_active=0;
	for(Int_t i=0;i<54;i++){
		for(Int_t j=0;j<num_module_in_row[i];j++){
			status[i][j]=1; //all active for PVDIS FAEC
			if (status[i][j]==1) total_module_active++;
		}
	}




};

TVector2 GetECALBlock_coord(Int_t block_id){
	Int_t idx=-1,idy=-1;
	TVector2 vec_ecalBlock;
	idx = block_map[block_id].second;
	idy = block_map[block_id].first;
	vec_ecalBlock.Set(floor(x[idy][idx]*100+0.5)/100,floor(y[idy][idx]*100+0.5)/100); 
	//printf("DEBUG : ecal id %i [%i,%i],[%f,%f] \n",block_id,idx,idy,x[idy][idx],y[idy][idx]);//vec_ecalBlock.Px(),vec_ecalBlock.Py()
	return vec_ecalBlock;
};

Int_t GetECALBlock_sector(Int_t block_id){
	Int_t idx = block_map[block_id].second;
	Int_t idy = block_map[block_id].first;

	return sector[idy][idx];
};

Int_t GetECALCluser(Int_t block_id,Int_t *cluster_edep_blockid){
	Int_t hit_idx = block_map[block_id].second;
	Int_t hit_idy = block_map[block_id].first;
	Int_t hit_around_idx[6]={100000};      //around 6 modules, a variable label indicates how many surrounded modules are around the hitted module
	Int_t hit_around_idy[6]={100000};
	//__________________________________find the surrounded other 6 modules________________________
	Int_t tmp_idx[15]={hit_idx-7, hit_idx-6, hit_idx-5, hit_idx-4, hit_idx-3, hit_idx-2, hit_idx-1, hit_idx, hit_idx+1, hit_idx+2, hit_idx+3, hit_idx+4, hit_idx+5, hit_idx+6, hit_idx+7};
	Int_t tmp_idy[3]={hit_idy-1, hit_idy, hit_idy+1};
	Int_t label=0;
	for(Int_t i=0;i<15;i++){// x scan
		for(Int_t j=0;j<3;j++){  //y scan
			if(tmp_idy[j]>=0 && tmp_idy[j]<54 && tmp_idx[i]>=0 && tmp_idx[i]<num_module_in_row[tmp_idy[j]] && (tmp_idx[i]!=hit_idx || tmp_idy[j]!=hit_idy) ){ // in range
				if(sqrt( pow( (x[tmp_idy[j]][tmp_idx[i]]-x[hit_idy][hit_idx]),2 )+ pow((y[tmp_idy[j]][tmp_idx[i]]-y[hit_idy][hit_idx]),2 ))<15.0){
					hit_around_idx[label]=tmp_idx[i];
					hit_around_idy[label]=tmp_idy[j];
					label++;
				}
			}
		}
	}
	cluster_edep_blockid[0] = id[hit_idy][hit_idx];
	for(int l=0;l<label;l++){
		cluster_edep_blockid[l+1]=id[hit_around_idy[l]][hit_around_idx[l]];
	}
	return label;
};


TGraphErrors * GetRMSReolution(TH2F * Difference){
	assert(Difference);

	TH1D * hpt = Difference->ProjectionX();

	int cnt = 0;
	Double_t pt[1000] = { 0 };
	Double_t resolution[2][1000]={0};
	Double_t resolution_error[2][1000]={0};
	Double_t fit_range_up[1000] = { 0 };
	Double_t fit_range_low[1000] = { 0 };
	for ( int i = 0; i < Difference->GetNbinsX() / 4; i++ )
	{

		TH1D * proj = Difference->ProjectionY("_py", i * 10 + 1, (i + 1) * 10);

		if ( proj->GetSum() > 100 )//check for entries
		{
			pt[cnt] = Difference->GetXaxis()->GetBinCenter(i * 10 + 10 / 2);//in GeV
			Double_t fit_range_up=1.2;
			Double_t fit_range_low=1.5;   

			TF1 * gaus = new TF1("gausEoverP", "gaus", fit_range_low, fit_range_up);
			proj -> Fit(gaus, "R0");
			resolution[0][cnt] = gaus->GetParameter(2)/gaus->GetParameter(1);//RMS[cnt];
			resolution[1][cnt] = proj->GetRMS()/proj->GetMean();
			resolution_error[0][cnt] = TMath::Sqrt(TMath::Power(gaus->GetParError(1)/gaus->GetParameter(1),2)+TMath::Power(gaus->GetParError(2)/gaus->GetParameter(2),2))*resolution[0][cnt];//RMSErr[cnt];
			resolution_error[1][cnt] = TMath::Sqrt(TMath::Power(proj->GetRMSError(),2)+TMath::Power(proj->GetMeanError(),2));//proj->GetRMSError()/proj->GetMean();
			cnt++;

		}

		delete proj;
	}

	//    delete Difference;
	TGraphErrors * gptRMS = new TGraphErrors(cnt, pt, resolution[0], 0, resolution_error[0]);
	gptRMS -> SetLineWidth(3);
	gptRMS -> Draw("*A");
	gptRMS -> GetYaxis()->CenterTitle();

	// degub print
	for(Int_t i=0;i<cnt;i++){
		// printf(" %13.2f  %13.3f %13.3f %13.3f %13.5f \n",pt[i],resolution[0][i],resolution_error[0][i],resolution[1][i],resolution_error[1][i]);
	}

	return gptRMS;  
};

Int_t GetRadiusIndex(Double_t radius){
	Int_t Ri=-1;
	for (Int_t i =0; i<8;i++){
		if (radius>trig_low_R[i] && radius<=trig_high_R[i]){
			Ri=i;
			break;
		}
	}
	return Ri;
};

Double_t GetThreshold6p1(Double_t radius){
	Double_t thresh=0;
	if (radius<trig_low_R[0] || radius>trig_high_R[7]){
		thresh=11000;//set a high number so that R below the ecal limit will be rejected
		//thresh=0;
		return thresh;
	} else {
		for(Int_t i=0;i<8;i++){
			if(radius>trig_low_R[i] && radius<=trig_high_R[i]){
				thresh=trig_thresh_6p1[i];
				break;
			}       
		}
	}
	return thresh;
};
Double_t Calfactor(Double_t edep){
	Double_t factor=0;
	Double_t convertE=edep/1000.0;
	if (edep<0){
		factor=0;

		return factor;
	} else {    
		factor=(0.190176+0.00815213*convertE-0.00044322*convertE*convertE);	
	}       

	return factor;
};

void set_plot_style()
{
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}
void LoadBackgrounds(){
	//background tree for now this works with merged bkg.
	ecalPSMap.clear();
	ecalShMap.clear();

	TVector2 vec_ecalBlock;
	Int_t nentries[5];
	Bool_t Fillecal=kFALSE;  

	Double_t pf;//for GEM
	Int_t i_pf=0;//for GEM plane hit index
	Double_t th;//for GEM
	Double_t r;//for GEM
	Double_t X[2];//for GEM hit coordinates


	Int_t evnum=-1;
	TFile *filebkg[5];
	TTree *TGEMC_Generate_bkg[5];
	TTree *TGEMC_Flux_bkg[5];
	TTree *TGEMC_EM_bkg[5];
	TTree *TGEMC_ECPS_bkg[5]; 
	const Double_t DEG=180./3.1415926;
	Int_t counttracks=0;
	//double detlaT = 6.13711543300417939e-05;//EM only
	double detlaT = 0.035183270; //HallD pim
	// double detlaT = 0.032784305; //hallD pi0
	// double detlaT = 0.035166937; //HallD pip
	double event_time=0.0;
	double time_windows=0;

	//    double dE=0;
	Bool_t processWindow=kFALSE;//when all the events in a single window are filled to the ecal, start cluster processing
	filebkg[0]=new TFile("/cache/halla/solid/user/tianye/container/pimBggen_all_rod_6mm/pimBggen_all_rod_6mm_995files.root");
	filebkg[1]=new TFile("/cache/halla/solid/user/tianye/container/pimBggen_all_rod_6mm_2/pimBggen_all_rod_6mm_2_991file.root");
	filebkg[2]=new TFile("/cache/halla/solid/user/tianye/container/pimBggen_all_rod_6mm_3/pimBggen_all_rod_6mm_993files.root");
	//event info
	total_ecalShMap.clear();
	total_ecalPSMap.clear();
	total_ecalShMapbkg.clear();
	total_ecalPSMapbkg.clear();
	double pid_gen_bkg=0;
	double bkg_number=0;
	for(int n=0;n<3;n++){
		TGEMC_Generate_bkg[n] = (TTree*) filebkg[n]->Get("generated");

		TGEMC_Generate_bkg[n]->SetBranchAddress("pid",&fGen_pid);
		TGEMC_Generate_bkg[n]->SetBranchAddress("px",&fGen_Px);
		TGEMC_Generate_bkg[n]->SetBranchAddress("py",&fGen_Py);
		TGEMC_Generate_bkg[n]->SetBranchAddress("pz",&fGen_Pz);
		TGEMC_Generate_bkg[n]->SetBranchAddress("vx",&fGen_vx);
		TGEMC_Generate_bkg[n]->SetBranchAddress("vy",&fGen_vy);
		TGEMC_Generate_bkg[n]->SetBranchAddress("vz",&fGen_vz);

		TGEMC_Flux_bkg[n]  = (TTree*) filebkg[n]->Get("flux");

		TGEMC_Flux_bkg[n]->SetBranchAddress("hitn",&fFluxHit_n);
		TGEMC_Flux_bkg[n]->SetBranchAddress("id",&fFluxHit_id);
		TGEMC_Flux_bkg[n]->SetBranchAddress("pid",&fFluxHit_pid);
		TGEMC_Flux_bkg[n]->SetBranchAddress("mpid",&fFluxHit_mpid);
		TGEMC_Flux_bkg[n]->SetBranchAddress("mtid",&fFluxHit_mtid);
		TGEMC_Flux_bkg[n]->SetBranchAddress("tid",&fFluxHit_tid);
		TGEMC_Flux_bkg[n]->SetBranchAddress("trackE",&fFluxHit_trackE);
		TGEMC_Flux_bkg[n]->SetBranchAddress("avg_x",&fFluxHit_Avg_x);
		TGEMC_Flux_bkg[n]->SetBranchAddress("avg_y",&fFluxHit_Avg_y);
		TGEMC_Flux_bkg[n]->SetBranchAddress("avg_z",&fFluxHit_Avg_z);
		TGEMC_Flux_bkg[n]->SetBranchAddress("px",&fFluxHit_Px);
		TGEMC_Flux_bkg[n]->SetBranchAddress("py",&fFluxHit_Py);
		TGEMC_Flux_bkg[n]->SetBranchAddress("pz",&fFluxHit_Pz);
		TGEMC_Flux_bkg[n]->SetBranchAddress("avg_t",&fFluxHit_T);
		TGEMC_Flux_bkg[n]->SetBranchAddress("vz",&fFluxHit_vz);

		TGEMC_EM_bkg[n] = (TTree*) filebkg[n]->Get("solid_ec");

		TGEMC_EM_bkg[n]->SetBranchAddress("id",&fECHit_id);
		TGEMC_EM_bkg[n]->SetBranchAddress("totEdep",&fECHit_totEdep);
		TGEMC_EM_bkg[n]->SetBranchAddress("totEdepB",&fECHit_totEdepB);
		TGEMC_EM_bkg[n]->SetBranchAddress("totEend",&fECHit_totEend);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_z",&fECHit_Avg_z);
		TGEMC_EM_bkg[n]->SetBranchAddress("mpid",&fECHit_mpid);
		TGEMC_EM_bkg[n]->SetBranchAddress("tid",&fECHit_tid);
		TGEMC_EM_bkg[n]->SetBranchAddress("mtid",&fECHit_mtid); 
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_lx",&fECHit_Avg_lx);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_ly",&fECHit_Avg_ly);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_lz",&fECHit_Avg_lz);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_x",&fECHit_Avg_x);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_y",&fECHit_Avg_y);
		TGEMC_ECPS_bkg[n] = (TTree*) filebkg[n]->Get("solid_ec_ps");

		TGEMC_ECPS_bkg[n]->SetBranchAddress("id",&fECPSHit_id);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("totEdep",&fECPSHit_totEdep);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("totEdepB",&fECPSHit_totEdepB);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("totEend",&fECPSHit_totEend);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_z",&fECPSHit_Avg_z);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("mpid",&fECPSHit_mpid);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("tid",&fECPSHit_tid);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("mtid",&fECPSHit_mtid); 
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_lx",&fECPSHit_Avg_lx);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_ly",&fECPSHit_Avg_ly);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_lz",&fECPSHit_Avg_lz);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_x",&fECPSHit_Avg_x);
		TGEMC_ECPS_bkg[n]->SetBranchAddress("avg_y",&fECPSHit_Avg_y);
		//total_ecalShMap.clear();
		//total_ecalPSMap.clear();
		nentries[n] = (Int_t)TGEMC_Generate_bkg[n]->GetEntries();
		printf("Merged Background Entries = %i \n",nentries[n]);
		Double_t totEdep_shower=0;
		Double_t totEdep_preshower=0;
		Fillecal=kFALSE;  

		processWindow=kFALSE;//when all the events in a single window are filled to the ecal, start cluster processing
		Double_t pf_flux_bkg;
		Double_t r_flux_bkg;
		for (int i=0; i<nentries[n]; i++) {
			TGEMC_Flux_bkg[n]->GetEntry(i); 
			TGEMC_EM_bkg[n]->GetEntry(i); 
			TGEMC_ECPS_bkg[n]->GetEntry(i); 
			TGEMC_Generate_bkg[n]->GetEntry(i);  
			event_time +=detlaT;
			time_windows = int(event_time/TriggerWindow);
			if(fTimeWindow<time_windows){

				fTimeWindow=time_windows;
				i--;
				processWindow=kTRUE;
				//  cout<<"event="<<i<<"timewindows="<<fTimeWindow<<"Evtime_windows="<<time_windows<<endl;
			}else{
				Fillecal=kFALSE;
				for (int j = 0; j<fFluxHit_id->size(); j++){
					FluxHit_detector_ID=fFluxHit_id->at(j)/1000000;
					FluxHit_subdetector_ID=(fFluxHit_id->at(j)%1000000)/100000;
					FluxHit_subsubdetector_ID=((fFluxHit_id->at(j)%1000000)%100000)/10000;
					FluxHit_component_ID=fFluxHit_id->at(j)%10000;  
					pf_flux_bkg=TMath::Sqrt(TMath::Power(fFluxHit_Px->at(j),2)+TMath::Power(fFluxHit_Py->at(j),2)+TMath::Power(fFluxHit_Pz->at(j),2));
					r_flux_bkg = TMath::Sqrt(TMath::Power(fFluxHit_Avg_x->at(j),2)+TMath::Power(fFluxHit_Avg_y->at(j),2));

					if (FluxHit_detector_ID==3 && FluxHit_subdetector_ID == 1 && FluxHit_subsubdetector_ID == 1 /*&& r_flux_bkg>r_min && r_flux_bkg<=r_max && fFluxHit_pid->at(j)==-211 && fFluxHit_tid->at(j)==1 && pf_flux_bkg>0 && pf_flux_bkg<=12000 && fFluxHit_Px->at(j)>0*/){
						Fillecal=kTRUE;

						break;
					}
				}

				if (Fillecal){
					bkg_number+=1;
					//processWindow=kTRUE;
					totEdep_shower=0;
					totEdep_preshower=0; 
					for (int j = 0; j<fECHit_id->size(); j++){
						ECHit_detector_ID=fECHit_id->at(j)/1000000;
						ECHit_subdetector_ID=(fECHit_id->at(j)%1000000)/100000;
						ECHit_subsubdetector_ID=((fECHit_id->at(j)%1000000)%100000)/10000;
						ECHit_component_ID=fECHit_id->at(j)%10000;
						if (ECHit_detector_ID==3 && ECHit_subdetector_ID == 1 && ECHit_subsubdetector_ID == 0){//shower 
							totEdep_shower+=fECHit_totEend->at(j);
							//cout<<"ECHit_component_ID="<<ECHit_component_ID<<"totEdep_shower="<<totEdep_shower<<endl;
							if (ecalShMapbkg.count(ECHit_component_ID)){
								ecalShMapbkg[ECHit_component_ID]+=fECHit_totEend->at(j);//5794;//106;// MeV
							} else {
								ecalShMapbkg[ECHit_component_ID]=fECHit_totEend->at(j);//5794;//106;// MeV
							}
							//cout<<"totEdep_shower="<<totEdep_shower<<endl;
						}
					}
					for (int k = 0; k<fECPSHit_id->size(); k++){
						ECPSHit_detector_ID=fECPSHit_id->at(k)/1000000;
						ECPSHit_subdetector_ID=(fECPSHit_id->at(k)%1000000)/100000;
						ECPSHit_subsubdetector_ID=((fECPSHit_id->at(k)%1000000)%100000)/10000;
						ECPSHit_component_ID=fECPSHit_id->at(k)%10000;
						//cout<<"event="<<i<<"E_preshower="<<totEdep_shower<<endl;
						if (ECPSHit_detector_ID==3 && ECPSHit_subdetector_ID == 1 && ECPSHit_subsubdetector_ID == 1){//pershower
							//if (solid_ec_avg_z[j]>2000){
							totEdep_preshower+=fECPSHit_totEend->at(k);

							if (ecalPSMapbkg.count(ECPSHit_component_ID)){
								ecalPSMapbkg[ECPSHit_component_ID]+=fECPSHit_totEend->at(k);// MeV
							} else {
								ecalPSMapbkg[ECPSHit_component_ID]=fECPSHit_totEend->at(k);// MeV
							}
							//}//end preshower
						}

					}//end EChit loop
				}//end of fillecal
				//after filling ecal object with ecal hits for 30 ns window now add these objects to the vector
				if (processWindow){ 
					//add ecal objects 
					total_ecalShMapbkg.push_back(ecalShMapbkg);
					total_ecalPSMapbkg.push_back(ecalPSMapbkg);
					processWindow=kFALSE;//reset till next window is filled

					//clear objects that reset after 30 ns window
					ecalShMapbkg.clear();//empty the ecalMap
					ecalPSMapbkg.clear();//empty the ecalMap
				}//end of processWindow

			}//end of if (fTimeWindow<fEvTimeWindow)

			if (i>100000 && bEarlyBeakBkg)
				break;

			if (i%100000==0)
				printf("Background Event %d \n",i);

		}//end of nentries
	}

	// cout<<"total_ecalPSMap.size()="<<total_ecalPSMap.size()<<endl;
	printf("Total Background Trigger windows (@ %3.0f ns) available %d \n",TriggerWindow, total_ecalPSMapbkg.size());

	//end of background. It's in the memory
};

void getBkgEdepPS(Int_t rndnum){//this routine will fill ecalPSMap with background signals
	ecalPSMapbkg = total_ecalPSMapbkg[rndnum];
};
void getBkgEdepSh(Int_t rndnum){//this routine will fill ecalShMap with background signals
	ecalShMapbkg = total_ecalShMapbkg[rndnum];
};	   
