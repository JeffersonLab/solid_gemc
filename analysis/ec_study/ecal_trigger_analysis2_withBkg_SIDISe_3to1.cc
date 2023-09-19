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
#include <stdio.h>
#include <limits.h>
#include "ecal.h"

#define __IO_MAXHIT 10000
#define __OCTANTS 30
Bool_t trigger_state_1[__OCTANTS]={kFALSE}; //octant level 6+1 blocks trigger state for the event
Bool_t trigger_state_2[__OCTANTS]={kFALSE}; //octant level 2+1 blocks trigger state for the event
Bool_t trigger_state_PS[__OCTANTS]={kFALSE};//if single event is above MIP then trigger state for the window will be true
Bool_t trigger_state_3={kFALSE};
Bool_t window_trigger_state_1=kFALSE;//trigger window level 6+1 blocks trigger state for the event
Bool_t window_trigger_state_2=kFALSE;//trigger window level 2+1 blocks trigger state for the event
Bool_t bDisableBkg=kFALSE;//kTRUE;//disable background signal loading
Bool_t bEarlyBeakBkg = kFALSE;//early break background loop
const Double_t DEG=180./3.1415926;
Bool_t processmodules= kFALSE;
Double_t fTotalRate[3][5];
const Double_t module_factor=1.0;
const Double_t TriggerWindow = 30;// ns
const Double_t DelayedHitTimeLimit = 50;// ns//or 30 ns was set before Tue Jun 23 11:02:34 EDT 2015
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
vector<int> *fFluxHit_n=0;
vector<int> *fFluxHit_id=0;
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

TChain * TGEMC_Flux;
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
TChain * TGEMC_EC;
Bool_t kSaveCanvas=kTRUE;

Int_t pid_gen=0;
Double_t theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;
Int_t FluxHit_detector_ID,FluxHit_subdetector_ID,FluxHit_subsubdetector_ID,FluxHit_component_ID;
Int_t ECHit_detector_ID,ECHit_subdetector_ID,ECHit_subsubdetector_ID,ECHit_component_ID;
Double_t trig_low_R[5]={90.0 ,105.0 ,115.0 ,130.0,150.0};
Double_t trig_high_R[5]={105.0 ,115.0 ,130.0 ,150.0,200.0};
Double_t trig_thresh_6p1[5] = {990.095,762.6,557.97,355.25,185.87};
Double_t trig_thresh_2p1[6] = {501.5 ,471.9 ,412.8 ,340.5 };//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0  };//not in use
Double_t trig_thresh_PS[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{20.9 ,28.2 ,28.3 ,27.7 ,27.5 ,29.0 ,31.7 ,17.7};//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0  };//not in use
Double_t trig_thresh_PS_6p1[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{20.9 ,28.2 ,28.3 ,27.7 ,27.5 ,29.0 ,31.7 ,17.7 };//{126.6 ,128.7 ,123.6 ,124.4 ,124.1 ,125.0 ,126.4 ,119.7};////Using 6p1 cluster sum as total PS sum threshold
Double_t trig_thresh_PS_2p1[6] = {0.0, 0.0, 0.0, 0.0,0.0,0.0};//{121.9 ,120.8 ,115.0 ,115.3 ,116.4 ,116.3 ,117.5 ,131.4};

//routine to load ecal block id,coordinates, and sector information
void LoadEC_map(TString map_file);
//routine to return ecal block x,y
TVector2 GetECALBlock_coord(Int_t block_id);
Int_t GetECALBlock_sector(Int_t block_id);
Int_t GetECALCluser(Int_t block_id,Int_t *cluster_edep_blockid);
Int_t GetECALlargeCluser(Int_t block_id,Int_t *cluster_edep_largeblockid,Int_t *cluster_edep_largegroupid);
Int_t sector[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t id[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t num_module_in_row[54]={0};
Double_t y_bak[54]={100000};
Double_t x[54][50]={{100000}};           
Double_t y[54][50]={{100000}}; 
Int_t status[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
int group_x[54][50]={100000};
int num_module_in_row_group[54]={0};
std::pair<Int_t,Int_t> block_map[2000];
Int_t GetRadiusIndex(Double_t radius);
Double_t GetThreshold6p1(Double_t radius);
Double_t Calfactor(Double_t edep);
//routines for properly access vector based TTree using TChain
void SetFluxBranchAddresses();
void GetFluxEntryByBranch(Long64_t local);
//void print2largest();
void SetECBranchAddresses();
void GetECEntryByBranch(Long64_t local);
Double_t edep_6p1_max1[__OCTANTS]={0};
Double_t edep_2p1_max1[__OCTANTS]={0};
Double_t radius_6p1_max1[__OCTANTS]={0};
int module_2p1_group1[9]={0};
Double_t edep_2p1_around[__OCTANTS]={0}; 
Double_t edep_6p1_x[__OCTANTS]={0};
Double_t edep_6p1_y[__OCTANTS]={0};
//get energy resolution
TGraphErrors * GetRMSReolution(TH2F * Difference);
struct event {
	Int_t pid;
	Double_t x;
	Double_t y;
	Double_t pf; 
	//  Double_t hit_time; 
	Int_t blockID;
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

void top3Repeated(int arr[], int n,vector<int> &group_number)
{
	if (n < 3) {
		cout << "Invalid Input";
		return;
	}

	std::map<int, int> fre;
	fre.clear();
	for (int i = 0; i < n; i++){
		fre[arr[i]]++;
	}
	for (std::map<int,int>::iterator it=fre.begin(); it!=fre.end(); ++it){
		if (it->second ==3) {
			group_number.push_back(it->first);
		}
	}
}



void  print2largest(double arr[], int arr_size, double &first, double &second,double &third, int &first_index,int &second_index, int &third_index)
{
	int i;

	if (arr_size < 3)
	{
		printf(" Invalid Input ");
		return;
	}

	third = first = second = 0;
	for (i = 0; i < arr_size ; i ++)
	{

		if (arr[i] > first)
		{
			third = second;
			third_index= second_index;
			second = first;
			second_index = first_index;
			first = arr[i];
			first_index =i;
		}
		else if (arr[i] > second)
		{
			third = second;
			third_index = second_index;
			second = arr[i];
			second_index = i;
		}
		else if (arr[i] > third){
			third = arr[i];
			third_index=i;
		}

	}
}


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
std::map<Int_t,Double_t> ecalPSMap;//map for pre-shower
std::map<Int_t,Double_t> ecalShMap;//map for shower
vector<std::map<Int_t,Double_t> >total_ecalShMapbkg;
vector<std::map<Int_t,Double_t> > total_ecalPSMapbkg;
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

	TH1F *Gen_mom_distr = new TH1F("Gen_mom_distr","Primary Track Momentum;Momentum (MeV)",mom_bins ,1000,8000);
	TH1F *Gen_theta_distr = new TH1F("Gentheta_distr","Primary Track Theta;#theta (deg)",mom_bins ,10,50);  
	TH1F *ECALDet_mom_distr = new TH1F("ECALDet_mom_distr","ECAL Front All e^{-} Track Momentum;Momentum (MeV)",mom_bins ,-8000,1000);
	TH1F *htotEdep_ec_shower=new TH1F("htotEdep_ec_shower","ec shower;totEdep(MeV);",100,0,4000);
	TH1F *htotEdep_ec_preshower=new TH1F("htotEdep_ec_preshower","ec preshower;totEdep(MeV);",100,0,500);
	TH1F *BlockX_hist=new TH1F("BlockX_hist",";#Delta X(cm);",1000,-500,500);
	TH1F *BlockY_hist=new TH1F("BlockY_hist",";#Delta Y(cm);",1000,-500,500);
	TString strig[3] = {"No. Trig","6+1 Trig.","2+1 Trig."};
	TString sradii[7] = {"0.9 - 1.05 m","1.05 - 1.15 m","1.15 - 1.30 m","1.30 - 1.50 m","1.50 - 2.0 m","2.0 - 2.3 m","0.9 - 2.3 m"};
	//Only index 0 - e+/-, 1 - pions. 2 - gamma are used in count_pid
	Int_t count_pid[30][3][5]={{{0}}};//eventually 5 pid type (e-, pi-,gamma,e+, and pi+) to count total, 6+1 trigger and 2+1 trigger would be used
	TH1F *histo_edep_Sh_radius[3][5];//total,6+1 trigger, 2+1 trigger for 8 radius bins and all radii
	TH1F *histo_edep_Sh_cluster_radius[3][5];//total,6+1 trigger, 2+1 trigger for 8 radius bins and all radii
	TH1F *histo_edep_PS_radius[3][5];
	TH1F *histo_edep_PS_cluster_radius[3][5];
	TH1F *histo_edep_Sh_radius_threshold[6];
	TH2F *histo_edep_Sh_radius_nothreshold[6];
	TH1F *histo_trig_window_residue_r[5]; //
	TH1F *histo_trig_window_residue_x[5]; //
	TH1F *histo_trig_window_residue_y[5]; //
	for(Int_t k=0;k<6;k++){
		histo_edep_Sh_radius_threshold[k]  = new TH1F(Form("histo_edep_Sh_radius_threshold_%d",k),Form("histo_edep_Sh_radius_threshold_%d",k),500,0,5000);//for Pf<1 change bin limit to 1500
		histo_edep_Sh_radius_nothreshold[k]  = new TH2F(Form("histo_edep_Sh_radius_nothreshold_%d",k),Form("histo_edep_Sh_radius_nothreshold_%d",k),500,0,5000,12000,0,12000);
		histo_edep_Sh_radius_threshold[k]->GetXaxis()->SetTitle("E_6p1 [MeV]");
		histo_edep_Sh_radius_threshold[k]->GetYaxis()->SetTitle("flux_p [MeV]");
		histo_edep_Sh_radius_nothreshold[k]->GetXaxis()->SetTitle("E_6p1 [MeV]");
		histo_edep_Sh_radius_nothreshold[k]->GetYaxis()->SetTitle("flux_p [MeV]");
	}
	Int_t ecal_bin_limit[4]={2000,2000,400,400};// for blocker shower total //{2500,1000,500,400};//for no photon blocker //{2000,600,400,250};// for blocker shower total,shower 6+1, PS total, PS 6+1 for P>1 1000,750,200,200 and for p<1 or total {2500,1000,500,400}
for(Int_t i=0;i<3;i++){
	for(Int_t j=0;j<5;j++){
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
TH1F *histo_pid2[3][5][7];//momentum distr on incident particles on ecal sens detector for flat elec. + background merged
TH1F *histo_pid_time[3][5];//hit time distr on incident particles on ecal sens detector 

TH1F *histo_delayed_pid[3][5];//momentum distr on incident particles on ecal sens detector 
TH1F *histo_delayed_pid_time[3][5];//hit time distr on incident particles on ecal sens detector 
TH1F *histo_pid_time_energyW[5];//hit time distr on incident particles on ecal sens detector 


TH2F *histo_rphi[3][5];//plot xy hit distr on ecal sens detector for e+-, any pion, gamma
TString spid[5] = {"Electron","Pion","Pi0 Gamma","Gamma","empty"};//e+-,Pi+-,gamma
Int_t ecal_bin_limit2[5]={8000,8000,1200,1200,8000};//for  blocker;  //for no blocker //{1000,5000,1200,1200,8000};//for blocker //{1000,1000,1000,1000,1000};//{10,10,10,10,10};//for no blocker {1200,5000,2000,1200,8000}
Int_t ecal_bin_limit3[5]={1000,7833.33,2000,7833.33,8000};//for  blocker;
Int_t ecal_bin_limit3_low[5]={0,1166.67,0,1166.67,0};
//
TString spid3[5] = {"Bkg. e^{#pm}","Pion","Gamma","e^{-}","empty"};//e+-,Pi+-,gamma
//cout<<"FAL======"<<spid3[1].Data()<<endl;
for(Int_t i=0;i<3;i++){
	for(Int_t j=0;j<5;j++){
		histo_pid[i][j]  = new TH1F(Form("Histo_pid_%d_%d",i,j),Form("%s Momentum (%s);Momentum (MeV)",spid[j].Data(),strig[i].Data()),mom_bins ,0,ecal_bin_limit2[j]);//for Pf<1 change bin limit to 1500
		for(Int_t k=0;k<7;k++){
			histo_pid2[i][j][k]  = new TH1F(Form("Histo_pid2_%d_%d_%d",i,j,k),Form("%s Momentum (%s);Momentum (MeV) (R %s)",spid3[j].Data(),strig[i].Data(),sradii[k].Data()),20,ecal_bin_limit3_low[j],ecal_bin_limit3[j]);//for Pf<1 change bin limit to 1500
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

TFile *output_file=new TFile("pim_mip_threshod_study.root","RECREATE");  
if (!bDisableBkg)
	LoadBackgrounds();
	else
	printf("Backgrounds are disabled! \n");
	//Flux bank
	TFile *file=new TFile("/cache/halla/solid/sim/solid_gemc/SIDIS_He3_JLAB_VERSION_1.3/pass7/background_solid_SIDIS_He3_dirty_even_e_FAEClarge_vacuum_field_1e6.root");
	TTree *tree_generated = (TTree*) file->Get("generated");

	tree_generated->SetBranchAddress("pid",&fGen_pid);
	tree_generated->SetBranchAddress("px",&fGen_Px);
	tree_generated->SetBranchAddress("py",&fGen_Py);
	tree_generated->SetBranchAddress("pz",&fGen_Pz);
	tree_generated->SetBranchAddress("vx",&fGen_vx);
	tree_generated->SetBranchAddress("vy",&fGen_vy);
	tree_generated->SetBranchAddress("vz",&fGen_vz);

	TTree *TGEMC_Flux = (TTree*) file->Get("flux");

	TGEMC_Flux->SetBranchAddress("hitn",&fFluxHit_n);
	TGEMC_Flux->SetBranchAddress("id",&fFluxHit_id);
	TGEMC_Flux->SetBranchAddress("pid",&fFluxHit_pid);
	TGEMC_Flux->SetBranchAddress("mpid",&fFluxHit_mpid);
	TGEMC_Flux->SetBranchAddress("mtid",&fFluxHit_mtid);
	TGEMC_Flux->SetBranchAddress("tid",&fFluxHit_tid);
	TGEMC_Flux->SetBranchAddress("trackE",&fFluxHit_trackE);
	TGEMC_Flux->SetBranchAddress("avg_x",&fFluxHit_Avg_x);
	TGEMC_Flux->SetBranchAddress("avg_y",&fFluxHit_Avg_y);
	TGEMC_Flux->SetBranchAddress("avg_z",&fFluxHit_Avg_z);
	TGEMC_Flux->SetBranchAddress("px",&fFluxHit_Px);
	TGEMC_Flux->SetBranchAddress("py",&fFluxHit_Py);
	TGEMC_Flux->SetBranchAddress("pz",&fFluxHit_Pz);
	TGEMC_Flux->SetBranchAddress("avg_t",&fFluxHit_T);
	TGEMC_Flux->SetBranchAddress("vz",&fFluxHit_vz);

	TTree *TGEMC_EC = (TTree*) file->Get("solid_ec");

	TGEMC_EC->SetBranchAddress("id",&fECHit_id);
	TGEMC_EC->SetBranchAddress("totEdep",&fECHit_totEdep);
	TGEMC_EC->SetBranchAddress("avg_z",&fECHit_Avg_z);
	TGEMC_EC->SetBranchAddress("mpid",&fECHit_mpid);
	TGEMC_EC->SetBranchAddress("tid",&fECHit_tid);
	TGEMC_EC->SetBranchAddress("mtid",&fECHit_mtid); 
	TGEMC_EC->SetBranchAddress("avg_lx",&fECHit_Avg_lx);
	TGEMC_EC->SetBranchAddress("avg_ly",&fECHit_Avg_ly);
	TGEMC_EC->SetBranchAddress("avg_lz",&fECHit_Avg_lz);
	TGEMC_EC->SetBranchAddress("avg_x",&fECHit_Avg_x);
	TGEMC_EC->SetBranchAddress("avg_y",&fECHit_Avg_y);
	Int_t nentries = (Int_t)tree_generated->GetEntries();
	printf("Entries = %i \n",nentries);
	Int_t treenumber=-1;
	Int_t ECtreenumber=-1;

	Bool_t Fillecal=kFALSE;  

	//flux hits
	Double_t pf;//for ECAL front det
	Double_t pf_flux;

	Double_t th;//for ECAL front det
	Double_t r;//for ECAL front det
	Double_t phi;
	Int_t i_pf=0;//for GEM plane hit index  
	//2D hash map to store hits in time window interval and other ecal related parameters
	//TPoint2DMap ecalPSMap;//map for pre-shower
	//TPoint2DMap ecalShMap;//map for shower
	TRandom *r3 = new TRandom3();
	ecalPSMap.clear();
	ecalShMap.clear();
	ecalPSMapbkg.clear();
	ecalShMapbkg.clear();   
	//load EC map file
	LoadEC_map("../layout/map_FAEC_ANL_20130628.txt");

	//vector to store ecal hit coordinates
	TVector2 vec_ecalBlock;
	Double_t ecalBlock_x;
	Double_t ecalBlock_y;
	TVector2 vec_ecalBlock_bkg;
	Double_t ecalBlock_x_bkg;
	Double_t ecalBlock_y_bkg;
	TVector2 vec_ecalCluster;
	Double_t ecalCluster_x;
	Double_t ecalCluster_y;
	Int_t ecal_clusterID;
	Int_t ecal_blockID;
	Int_t ecal_hitblockID;
	Int_t ecal_hitblockIDbkg;
	Double_t ecalBlock_x0;
	Double_t ecalBlock_y0;
	Int_t ecal_blockID0;
	Double_t ecal_blockE0=0.0;
	Double_t ecal_blockE1=0.0;
	Double_t ecal_blockE2=0.0;
	Double_t ecal_blockE3=0.0;
	Double_t ecal_blockE4=0.0;
	Double_t ecal_blockE5=0.0;
	Double_t ecal_blockE6=0.0;

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
Double_t cluster_edepbkg[7]={0};
Int_t cluster_edep_blockidbkg[7]={0};


Double_t cluster_largeedep[40]={0};
Double_t cluster_largeedepbkg[40]={0};
Int_t cluster_edep_largeblockid[40]={0};
Int_t cluster_edep_largegroupid[40]={0};
Int_t block_count=0;
Int_t block_countbkg=0;
double Edep_maximum[10];
Int_t moreblock_count=0;
Int_t moreblock_countbkg=0;

Int_t block_sector=0;
Double_t edep_6p1_sum;  //sum of 6+1 blocks 
Double_t edep_total_sum;  //sum ofall the ecal blocks 
Double_t edep_1block_sum;//photon sum of the main ecal block
// Double_t edep_6p1_max=0;
Double_t radius_6p1_max=0;
Double_t radius_6p1=0;
Double_t edep_6p1_sumbkg=0;
Double_t edep_2p1_sumbkg=0;
Double_t radius_6p1_bkg=0;
Double_t edep_sum0;  //photon sum of the  ecalBlock_x0,ecalBlock_y0 ecal block
Double_t DeltaRrDiff;
Double_t edep_6p1_PS_sum;  //sum of 6+1 blocks 
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
edep_sum0 = 0;
edep_total_sum = 0;
Double_t totEdep_shower=0;
Double_t totEdep_preshower=0;
//PS ecal cluster information
edep_6p1_PS_sum=0;
edep_total_PS_sum=0;	
double EC_averagex=0.0;
double EC_averagey=0.0;
double EC_localx=0.0;
double EC_localy=0.0;
int maxmodule_index[6];
int flux_size=0;
int N_loss=0;
int N_total=0;
int diff_moduleID1=0;
int diff_moduleID2=0;
int diff_moduleID3=0;
int diff_moduleID4=0;
int diff_moduleID5=0;
int diff_moduleID6=0;
vector <int> moduleID_group;

for (Int_t i=0; i<nentries; i++) {
	double pf_satisfy=0;
	Long64_t local;
	tree_generated->GetEntry(i);
	for (int j=0;j<fGen_pid->size();j++) {
		pid_gen=fGen_pid->at(j);
		px_gen=fGen_Px->at(j);
		py_gen=fGen_Py->at(j);
		pz_gen=fGen_Pz->at(j);
		vx_gen=fGen_vx->at(j);
		vy_gen=fGen_vy->at(j);
		vz_gen=fGen_vz->at(j);
		p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
		theta_gen=TMath::ACos(pz_gen/p_gen)*DEG;
		phi_gen=TMath::ATan2(py_gen,px_gen)*DEG;
		Gen_mom_distr->Fill(p_gen);
		Gen_theta_distr->Fill(theta_gen);
		if(i==50649){
			// cout<<"event="<<i<<"  "<<"p_gen="<<"  "<<p_gen<<endl;
		}
	}

	//FLUX bank

	TGEMC_Flux->GetEntry(i);
	//if(i==334){

	//}
	for (int j = 0; j<fFluxHit_id->size(); j++){
		FluxHit_detector_ID=fFluxHit_id->at(j)/1000000;
		FluxHit_subdetector_ID=(fFluxHit_id->at(j)%1000000)/100000;
		FluxHit_subsubdetector_ID=((fFluxHit_id->at(j)%1000000)%100000)/10000;
		FluxHit_component_ID=fFluxHit_id->at(j)%10000;        
		r = 0.1*TMath::Sqrt(TMath::Power(fFluxHit_Avg_x->at(j),2)+TMath::Power(fFluxHit_Avg_y->at(j),2));
		pf_flux=TMath::Sqrt(TMath::Power(fFluxHit_Px->at(j),2)+TMath::Power(fFluxHit_Py->at(j),2)+TMath::Power(fFluxHit_Pz->at(j),2));
		if (FluxHit_detector_ID==3 && FluxHit_subdetector_ID == 1 && FluxHit_subsubdetector_ID == 1 && theta_gen>7.5 && theta_gen<14.85 /*&& fFluxHit_pid->at(j)==11 &&  pf_flux>1000 && pf_flux<=12000*/ && fFluxHit_tid->at(j)==1 ){
			th=TMath::ATan(r/(fFluxHit_Avg_z->at(j) - targ_offset))*DEG;
			ECALDet_mom_distr->Fill(fFluxHit_vz->at(j));
			pf_satisfy+=pf_flux;
			flux_size =fFluxHit_id->size();
			processWindow=kTRUE;//process this event
			single_event.pid=fFluxHit_pid->at(j);
			single_event.x=fFluxHit_Avg_x->at(j)*0.1;
			single_event.y=fFluxHit_Avg_y->at(j)*0.1;
			X=fFluxHit_Avg_x->at(j)*0.1;
			Y=fFluxHit_Avg_y->at(j)*0.1;
			//radius_6p1_max = TMath::Sqrt(TMath::Power(X,2)+TMath::Power(Y,2));
			single_event.pf=pf_flux;
			single_event.blockID =FluxHit_component_ID;    
			event_list.push_back(single_event);          
			//event_list.push_back(single_event);
			Fillecal=kTRUE;

			//cout<<"satisfyevent="<<i<<endl;
			break;
		} else{
			Fillecal=kFALSE;
			// processWindow=kFALSE;
		}
	} 
	if (Fillecal){

		TGEMC_EC->GetEntry(i);
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
				//cout<<"avg_z="<<fECHit_Avg_z->at(j)<<endl;
				if (fECHit_Avg_z->at(j)>2000){
					EC_averagex=fECHit_Avg_x->at(j)*0.1;
					EC_averagey=fECHit_Avg_y->at(j)*0.1;
					EC_localx=fECHit_Avg_lx->at(j)*0.1;
					EC_localy=fECHit_Avg_ly->at(j)*0.1;
					totEdep_shower+=fECHit_totEdep->at(j);
					// if(i==2504){

					// }
					if (ecalShMap.count(ECHit_component_ID)){
						ecalShMap[ECHit_component_ID]+=fECHit_totEdep->at(j);// MeV
					} else {
						ecalShMap[ECHit_component_ID]=fECHit_totEdep->at(j);// MeV
					}
				} //end shower 
			}
			//cout<<"event="<<i<<"E_preshower="<<totEdep_shower<<endl;
			if (ECHit_detector_ID==3 && ECHit_subdetector_ID == 1 && ECHit_subsubdetector_ID == 1){//pershower
				if (fECHit_Avg_z->at(j)>2000){
					totEdep_preshower+=fECHit_totEdep->at(j);
					if (ecalPSMap.count(ECHit_component_ID)){
						ecalPSMap[ECHit_component_ID]+=fECHit_totEdep->at(j);// MeV
					} else {
						ecalPSMap[ECHit_component_ID]=fECHit_totEdep->at(j);// MeV
					}
				}//end preshower
			}

		}//end EChit loop

	}
	//process hits in this time window
	if (processWindow){  
		edep_6p1_sum=0;
		for (Int_t c=0;c<7;c++){//init cluster array
			cluster_edep[c]=0;
			cluster_edep_blockid[c]=0;
			cluster_edepbkg[c]=0;
			cluster_edep_blockidbkg[c]=0;
			maxmodule_index[c]=0;
		}

		for (Int_t c=0;c<30;c++){
			cluster_largeedep[c]=0;
			edep_6p1_max1[c] = 0;
			edep_2p1_max1[c] = 0;
			radius_6p1_max1[c]={0};
			edep_6p1_x[c]={0};
			edep_6p1_y[c]={0};
		}   


		for (Int_t c=0;c<40;c++){
			cluster_largeedep[c]=0;  
			cluster_edep_largegroupid[c]=0;
			cluster_edep_largeblockid[c]=0;
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
			/********************method1**************************/
			cluster_edep[0]=ecalShMap[cluster_edep_blockid[0]];
			edep_6p1_sum = cluster_edep[0];
			for(Int_t j=0;j<block_count;j++){
				cluster_edep[j+1]=ecalShMap[cluster_edep_blockid[j+1]];
				edep_6p1_sum +=cluster_edep[j+1];
			}
			/**************************************************/



			if (edep_6p1_sum >= edep_6p1_max1[octantno]){
				ecalBlock_x0=ecalBlock_x;
				ecalBlock_y0=ecalBlock_y;
				edep_6p1_x[octantno]=ecalBlock_x;
				edep_6p1_y[octantno]=ecalBlock_y;
				radius_6p1_max1[octantno]=TMath::Sqrt(TMath::Power(ecalBlock_x,2)+TMath::Power(ecalBlock_y,2));
				edep_6p1_max1[octantno]=edep_6p1_sum;
				ecal_blockID0 = ecal_blockID;
				maxmodule_index[0]=cluster_edep_blockid[0];
				maxmodule_index[1]=cluster_edep_blockid[1];
				maxmodule_index[2]=cluster_edep_blockid[2];
				maxmodule_index[3]=cluster_edep_blockid[3];
				maxmodule_index[4]=cluster_edep_blockid[4];
				maxmodule_index[5]=cluster_edep_blockid[5];
				maxmodule_index[6]=cluster_edep_blockid[6];
				ecal_blockE0 = cluster_edep[0];
				ecal_blockE1 = cluster_edep[1];
				ecal_blockE2 = cluster_edep[2];
				ecal_blockE3 = cluster_edep[3];
				ecal_blockE4 = cluster_edep[4];
				ecal_blockE5 = cluster_edep[5];
				ecal_blockE6 = cluster_edep[6];
			}

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
				}//end of 6+1 shower cluster
			} 
			//     }//end of processing shower


			moreblock_count = GetECALlargeCluser(ecal_blockID0,cluster_edep_largeblockid, cluster_edep_largegroupid);
			//cout<<"hit_ID="<<ecal_blockID0<<"modules="<<moreblock_count<<endl;
			std::map<int,int> ecalMap;
			ecalMap.clear();
			vector <int> group_number;
			group_number.clear();
			int arr[40]={0};
			// vector <int> moduleID_group;
			moduleID_group.clear();
			for(int m=0;m<moreblock_count;m++){
				// if(cluster_edep_largegroupid[m]!=0){
				// cout<<"surround ID="<<cluster_edep_largeblockid[m]<<"   "<<"groupID="<<cluster_edep_largegroupid[m]<<endl;
				arr[m]=cluster_edep_largegroupid[m];
				ecalMap[cluster_edep_largeblockid[m]]=cluster_edep_largegroupid[m];
				// }
			}
			if(ecalMap[ecal_blockID0]!=0){
				processmodules= kTRUE;
			}else processmodules= kFALSE;
			int n = sizeof(arr) / sizeof(arr[0]);
			top3Repeated(arr, n,group_number); 
			for(int i = 0; i < group_number.size(); i++){
				for (std::map<int,int>::iterator it = ecalMap.begin(); it != ecalMap.end(); ++it ){
					if( it->second == group_number[i] ){
						moduleID_group.push_back(it->first);
						//   cout<<"group_number="<<group_number[i]<<"  "<<"cluster module ID="<<it->first<<endl;
					}
				}
			}
			//double Edep_maximum[7];
			for(int i = 0; i < moduleID_group.size(); i++){
				cluster_largeedep[i]=ecalShMap[moduleID_group[i]];
				// cout<<"moduleID="<<moduleID_group[i]<<"  "<<"Edep="<<cluster_largeedep[i]<<endl;
			}

	}//end of processing shower
	if(processmodules){
		Edep_maximum[0]=cluster_largeedep[0]+cluster_largeedep[1]+cluster_largeedep[2];
		Edep_maximum[1]=cluster_largeedep[3]+cluster_largeedep[4]+cluster_largeedep[5];
		Edep_maximum[2]=cluster_largeedep[6]+cluster_largeedep[7]+cluster_largeedep[8];
		Edep_maximum[3]=cluster_largeedep[9]+cluster_largeedep[10]+cluster_largeedep[11];
		Edep_maximum[4]=cluster_largeedep[12]+cluster_largeedep[13]+cluster_largeedep[14];
		Edep_maximum[5]=cluster_largeedep[15]+cluster_largeedep[16]+cluster_largeedep[17];
		Edep_maximum[6]=cluster_largeedep[18]+cluster_largeedep[19]+cluster_largeedep[20];
		Edep_maximum[7]=cluster_largeedep[21]+cluster_largeedep[22]+cluster_largeedep[23];
		Edep_maximum[8]=cluster_largeedep[24]+cluster_largeedep[25]+cluster_largeedep[26];
		Edep_maximum[9]=cluster_largeedep[27]+cluster_largeedep[28]+cluster_largeedep[29];


		// find out first three maximum numbers
		//double arr[] = {ecal_blockID0,ecal_blockE1, ecal_blockE2, ecal_blockE3, ecal_blockE4, ecal_blockE5, ecal_blockE6};
		int E_size = sizeof(Edep_maximum)/sizeof(Edep_maximum[0]);
		double first, second,third;
		int first_index, second_index,third_index;
		print2largest(Edep_maximum, E_size,first,second,third,first_index, second_index,third_index);
		edep_2p1_max1[octantno] = first+second+third;
		int index1=first_index*3;
		module_2p1_group1[0]=moduleID_group[index1];
		module_2p1_group1[1]=moduleID_group[index1+1]; 
		module_2p1_group1[2]=moduleID_group[index1+2];

		module_2p1_group1[3]=moduleID_group[second_index*3];
		module_2p1_group1[4]=moduleID_group[second_index*3+1];
		module_2p1_group1[5]=moduleID_group[second_index*3+2];


		module_2p1_group1[6]=moduleID_group[third_index*3];
		module_2p1_group1[7]=moduleID_group[third_index*3+1];
		module_2p1_group1[8]=moduleID_group[third_index*3+2];
	}else edep_2p1_max1[octantno]=0;
	/*********************************************load bkg map*************************************************************************************************/
	if (!bDisableBkg){
		edep_total_sum=0;
		rndTimeWindow = (Int_t)(r3->Uniform()*(total_ecalShMapbkg.size()-1));
		getBkgEdepPS(rndTimeWindow);
		getBkgEdepSh(rndTimeWindow);
		edep_2p1_sumbkg=0;
		for(Int_t j=0;j<9;j++){
			edep_2p1_sumbkg += ecalShMapbkg[module_2p1_group1[j]];
		}
	}

	/******************************************************************************************************************************************************/  
	Int_t hit_pid_index;
	Int_t hit_pid_index2;
	Double_t hit_phi;
	Double_t event_p;
	Double_t DeltaX; 
	Double_t DeltaY;
	Double_t DeltaR;
	Int_t hit_octant;
	Int_t Ri;
	if(edep_2p1_max1[octantno]>0){

		for(Int_t e=0;e<event_list.size();e++){

			event_p = event_list[e].pf/1000;
			//cout<<"eventlist=="<<i<<"pf="<<event_list[e].pf<<endl;
			hit_pid_index =  pidmap[event_list[e].pid];//pid index starts from 1 but array index starts from 0
			hit_octant = 0;//since there will only be one event per window I disabled sector id//Int_t(phi/12);
			if (hit_octant>29)
				hit_octant=29;     
			if (hit_pid_index>0 && hit_pid_index<=5 && (/*event_p > 1.0 &&*/ event_p <= 12)){
				if (hit_pid_index==1 || hit_pid_index==4)//e+-
					hit_pid_index2 = 0;
				if (hit_pid_index==2 || hit_pid_index==5)//any pion
					hit_pid_index2 = 1;
				if (hit_pid_index==3){//gamma
					if (event_list[e].pid==111)//no simpid in the electron+bkg merged root file
						hit_pid_index2 = 2;//photons from pi0 simulation
					else
						hit_pid_index2 = 3;//photons from other sources
				}


				histo_pid[0][hit_pid_index2]->Fill(event_list[e].pf);//momentum distribution

				//cout<<"Event="<<i<<"E61="<<edep_6p1_max<<"R="<<radius_6p1_max<<"ThreshE="<<GetThreshold6p1(radius_6p1_max)<<endl;
				DeltaRrDiff = TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2));
				if(DeltaRrDiff>230){
					//cout<<"R="<<DeltaRrDiff<<endl;

				}
				Ri=GetRadiusIndex(TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2)));
				double High_pf= 4.6*(edep_6p1_max1[hit_octant]+30)+350;
				double Low_pf= 4.25*(edep_6p1_max1[hit_octant]+50)-225;
				N_total+=1;
				N_loss+=1; 
				BlockX_hist->Fill(edep_6p1_x[octantno]-event_list[e].x);
				BlockY_hist->Fill(edep_6p1_y[octantno]-event_list[e].y);
				if (hit_pid_index2 == 0 && event_list[e].pf< 1000)//bkg e+/- less than 1 GeV
					histo_pid2[0][0][5]->Fill(event_list[e].pf);//momentum distribution
				if (hit_pid_index == 1 && event_list[e].pf >= 1000){//high energy electrons
					histo_pid2[0][3][5]->Fill(event_list[e].pf);//momentum distribution        
					histo_pid2[0][3][Ri]->Fill(event_list[e].pf);
					if(90<DeltaRrDiff && DeltaRrDiff<=105 ){

						histo_edep_Sh_radius_nothreshold[0]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}else if(105<DeltaRrDiff && DeltaRrDiff<=115  ){
						histo_edep_Sh_radius_nothreshold[1]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}else if(115<DeltaRrDiff && DeltaRrDiff<=130 ){
						histo_edep_Sh_radius_nothreshold[2]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}else if(130<DeltaRrDiff&& DeltaRrDiff<=150 ){
						histo_edep_Sh_radius_nothreshold[3]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}else if(150<DeltaRrDiff && DeltaRrDiff<=200 ){
						histo_edep_Sh_radius_nothreshold[4]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}else if(200<DeltaRrDiff && DeltaRrDiff<=230 ){
						histo_edep_Sh_radius_nothreshold[5]->Fill(edep_2p1_max1[octantno],event_list[e].pf);

					}
					if(90<DeltaRrDiff && DeltaRrDiff<=105 && event_list[e].pf<=5500 && event_list[e].pf>=4500 ){
						histo_edep_Sh_radius_threshold[0]->Fill(edep_2p1_max1[hit_octant]);

					}else if(105<DeltaRrDiff && DeltaRrDiff<=115 && event_list[e].pf<=4400 && event_list[e].pf>=3600 ){
						histo_edep_Sh_radius_threshold[1]->Fill(edep_2p1_max1[hit_octant]);

					}else if(115<DeltaRrDiff && DeltaRrDiff<=130 && event_list[e].pf<=3300 && event_list[e].pf>=2700 ){
						histo_edep_Sh_radius_threshold[2]->Fill(edep_2p1_max1[hit_octant]);

					}else if(130<DeltaRrDiff && DeltaRrDiff<=150 && event_list[e].pf<=2200 && event_list[e].pf>=1800){
						histo_edep_Sh_radius_threshold[3]->Fill(edep_2p1_max1[hit_octant]);

					}else if(150<DeltaRrDiff && DeltaRrDiff<=200 && event_list[e].pf<=1100 && event_list[e].pf>=900){
						histo_edep_Sh_radius_threshold[4]->Fill(edep_2p1_max1[hit_octant]);

					}else if(200<DeltaRrDiff && DeltaRrDiff<=240 && event_list[e].pf<=2200 && event_list[e].pf>=1800){
						histo_edep_Sh_radius_threshold[5]->Fill(edep_2p1_max1[hit_octant]);

					}

				}
				if (hit_pid_index2 == 1){//any pion
					histo_pid2[0][1][5]->Fill(event_list[e].pf);//momentum distribution
					//cout<<"Ri================="<<Ri<<"Radius=="<<DeltaR<<endl;
					histo_pid2[0][1][Ri]->Fill(event_list[e].pf);
				}
				if (hit_pid_index == 3)//any gamma
					histo_pid2[0][2][5]->Fill(event_list[e].pf);//momentum distribution


				if((edep_2p1_max1[octantno]+edep_2p1_sumbkg*module_factor) >GetThreshold6p1(DeltaRrDiff)){

					DeltaR = (TMath::Sqrt(TMath::Power(event_list[e].x,2)+TMath::Power(event_list[e].y,2)) - TMath::Sqrt(TMath::Power(edep_6p1_x[hit_octant],2)+TMath::Power(edep_6p1_y[hit_octant],2)));
					histo_pid[1][hit_pid_index2]->Fill(event_list[e].pf);//update the total hits
					if (hit_pid_index2 == 0 && event_list[e].pf < 1000)//bkg e+/- less than 1 GeV
						histo_pid2[1][0][5]->Fill(event_list[e].pf);//momentum distribution
					if (hit_pid_index == 1 && event_list[e].pf >= 1000){//high energy electrons
						histo_pid2[1][3][5]->Fill(event_list[e].pf);//momentum distribution
						histo_pid2[1][3][Ri]->Fill(event_list[e].pf);
					}
					if (hit_pid_index2 == 1){//any pion
						histo_pid2[1][1][5]->Fill(event_list[e].pf);//momentum distribution
						histo_pid2[1][1][Ri]->Fill(event_list[e].pf);
					}
					if (hit_pid_index == 3)//any gamma
						histo_pid2[1][2][5]->Fill(event_list[e].pf);//momentum distribution	
				}    
			}	     

		}//end eventlist
	}//6p1E>0
	processWindow=kFALSE;//reset till next window is filled
	ecalPSMap.clear();
	ecalShMap.clear();
	event_list.clear();
	Fillecal=kFALSE;

}//endProcesswindow
}//end entries
cout<<"N_toal="<<N_total<<"  "<<"N_loss="<<N_loss<<endl;
output_file->Close();

//efficiency plots
Double_t p[6][100] = { {0} },p_2[6][100] = { {0} };
Double_t eff[6][100] = { {0} },eff_2[6][100] = { {0} };
Double_t efferr[6][100] = { {0} },efferr_2[6][100] = { {0} };
int cnt[6] = {0};
int cnt_2[6] = {0};
TGraphErrors * geffi_total_pion[6];
//  printf("for Pions \n Momentum \t Efficiency \t Error \n");
for (Int_t Ri=0;Ri<6;Ri++){
	for ( int i = 1; i <= histo_pid2[1][1][Ri]->GetNbinsX(); i++ ){
		if ( histo_pid2[1][1][Ri]->GetBinContent(i) <= histo_pid2[0][1][Ri]->GetBinContent(i) && histo_pid2[0][1][Ri]->GetBinContent(i) > 20 ){
			p[Ri][cnt[Ri]] = histo_pid2[1][1][Ri]->GetBinLowEdge(i)/1000;
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

printf("Double_t R_low[5] = {");
for(Int_t i=0;i<5;i++)
printf("%4.1f,",trig_low_R[i]);
printf("}; \n");
printf("Double_t R_up[5] = {");
for(Int_t i=0;i<5;i++)
printf("%4.1f,",trig_high_R[i]);
printf("}; \n");

Double_t mom_bin_width = 0.28;
printf("Double_t P_low[%d] = {",cnt[5]);
for(Int_t i=0;i<cnt[5];i++)
printf("%4.1f,",p[5][i]-mom_bin_width/2);
printf("}; \n");
printf("Double_t P_high[%d] = {",cnt[5]);
for(Int_t i=0;i<cnt[5];i++)
printf("%4.1f,",p[5][i]+mom_bin_width/2);
printf("}; \n");


printf("Double_t trig_pi_eff[6][%d] = { ",cnt[5]);
for(Int_t Ri=0;Ri<6;Ri++){
	printf("{ ");
	if (cnt[5] != cnt[Ri]) {
		printf("momentum array length mistmatch, [5] and [%d]",cnt[Ri]);
		break;
	}
	for(Int_t Pi=0;Pi<cnt[Ri];Pi++)
		printf("%4.3f,", eff[Ri][Pi]);
	printf("},");
}
printf("}; \n");

printf("Double_t trig_pi_eff_err[6][%d] = { ",cnt[5]);
for(Int_t Ri=0;Ri<6;Ri++){
	printf("{ ");
	if (cnt[5] != cnt[Ri]) {
		printf("momentum array length mistmatch, [5] and [%d]",cnt[Ri]);
		break;
	}
	for(Int_t Pi=0;Pi<cnt[Ri];Pi++)
		printf("%4.3f,", efferr[Ri][Pi]);
	printf("},");
}
printf("}; \n");

//cnt = 0;
printf("for Electrons \n Momentum \t Efficiency \t Error \n");
TGraphErrors * geffi_total_electron[6];
for (Int_t Ri=0;Ri<6;Ri++){
	for ( int i = 1; i <= histo_pid2[1][3][Ri]->GetNbinsX(); i++ ){
		if ( histo_pid2[1][3][Ri]->GetBinContent(i) <= histo_pid2[0][3][Ri]->GetBinContent(i) && histo_pid2[0][3][Ri]->GetBinContent(i) > 20){
			p_2[Ri][cnt_2[Ri]] = histo_pid2[1][3][Ri]->GetBinLowEdge(i)/1000;
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

printf("Double_t trig_ele_eff[6][%d] = { ",cnt_2[5]);
for(Int_t Ri=0;Ri<6;Ri++){
	printf("{ ");
	if (cnt_2[5] != cnt_2[Ri]) {
		printf("momentum array length mistmatch, [5] and [%d]",cnt_2[Ri]);
		break;
	}
	for(Int_t Pi=0;Pi<cnt_2[Ri];Pi++)
		printf("%4.3f,", eff_2[Ri][Pi]);
	printf("},");
}
printf("}; \n");

printf("Double_t trig_ele_eff_err[6][%d] = { ",cnt_2[5]);
for(Int_t Ri=0;Ri<6;Ri++){
	printf("{ ");
	if (cnt_2[5] != cnt_2[Ri]) {
		printf("momentum array length mistmatch, [5] and [%d]",cnt_2[Ri]);
		break;
	}
	for(Int_t Pi=0;Pi<cnt_2[Ri];Pi++)
		printf("%4.3f,", efferr_2[Ri][Pi]);
	printf("},");
}
printf("}; \n");


TCanvas *canvas_eff_pion = new TCanvas("canvas_eff_pion","canvas_eff_pion", 1000, 700);
geffi_total_pion[5]->Draw("AP*");
geffi_total_pion[5]->SetTitle(Form("Pion Efficiency ;Momentum (GeV);Efficiency"));
gPad->Update();
TCanvas *canvas_eff_electron = new TCanvas("canvas_eff_electron","canvas_eff_electron", 1000, 700);
geffi_total_electron[5]->Draw("AP*");
geffi_total_electron[5]->SetTitle(Form("Electron Efficiency ;Momentum (GeV);Efficiency"));
gPad->Update();
TMultiGraph *mg_eff_electron = new TMultiGraph();

TLegend *leg1 = new TLegend(0.71,0.31,0.89,0.59,NULL,"brNDC");//0.6,0.78,0.9,0.90//numbers obtained by first moving the legend to a prefered location and then looking at the inspect pop-up window
leg1->SetTextFont(62);
//leg1->SetTextSize(0.03);
leg1->SetFillColor(0);
leg1->SetFillStyle(1001);

for(Int_t i=0;i<5;i++){
	geffi_total_electron[i]->SetLineStyle(i+1);
	geffi_total_electron[i]->SetLineColor(1+i);//high contrast colors    
	geffi_total_electron[i]->SetLineWidth(2);//high contrast colors
	mg_eff_electron->Add(geffi_total_electron[i],"lp");    
	leg1->AddEntry(geffi_total_electron[i],Form("%s",sradii[i].Data()),"l");
}
mg_eff_electron->SetName("mg_eff_electron");


TCanvas *canvas_mg_eff_electron = new TCanvas("canvas_mg_eff_electron","canvas_mg_eff_electron", 1000, 700);
mg_eff_electron->SetTitle("Electron Efficiency ;Momentum (GeV);Efficiency");
mg_eff_electron->SetMaximum(1.05);
mg_eff_electron->SetMinimum(0.0);
mg_eff_electron->Draw("A");
leg1->Draw();
gPad->Update();
TMultiGraph *mg_eff_pion = new TMultiGraph();

TLegend *leg2 = new TLegend(0.11,0.71,0.30,0.89,NULL,"brNDC");//0.6,0.78,0.9,0.90
leg2->SetTextFont(62);
leg2->SetFillColor(0);
leg2->SetFillStyle(1001);

for(Int_t i=0;i<5;i++){
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
	canvas_eff_electron->SaveAs("/home/tianye/solid/plots/canvas_eff_electron_R_Track_FAL.png");
	canvas_mg_eff_electron->SaveAs("/home/tianye/solid/plots/canvas_eff_radii_electron_R_Track_FAL.png");
}



/*****************************************************************************************/

TCanvas * canvas_6plus1_threshold = new TCanvas("canvas_6plus1_threshold","canvas_6plus1_threshold",1400,400);
canvas_6plus1_threshold->Divide(2,3);
TLine *line[5];

for(int i=0;i<5;i++){
	canvas_6plus1_threshold->cd(i+1);
	line[i] = new TLine(trig_thresh_6p1[i],0,trig_thresh_6p1[i],30000);
	//histo_pid2[1][3][i]->Draw();
	histo_edep_Sh_radius_threshold[i]->Draw();
	line[i]->SetLineColor(kRed);
	line[i]->Draw("same");
}    
TLine *line_cut[5];
TF1 *function1=new TF1("function1","4.25*(x+50)-225",0,3000);
TCanvas * canvas_6plus1_nothreshold = new TCanvas("canvas_6plus1_nothreshold","canvas_6plus1_nothreshold",1400,400);
canvas_6plus1_nothreshold->Divide(2,3);
for(int j=0;j<5;j++){
	canvas_6plus1_nothreshold->cd(j+1);
	line_cut[j] = new TLine(trig_thresh_6p1[j],0,trig_thresh_6p1[j],10000);
	histo_edep_Sh_radius_nothreshold[j]->Draw("COLZ");
	line_cut[j]->SetLineColor(kRed);
	line_cut[j]->Draw("same");
	function1->SetLineColor(1);
	function1->Draw("same");

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

watch.Stop();
//printf("Total Time %3.4f min \n Completed, Exit now :) \n",watch.RealTime()/60);
theApp.Run();
return(1);
}

void SetFluxBranchAddresses(){
	TGEMC_Flux->SetBranchAddress("hitn",&fFluxHit_n,&bFluxHit_n);
	TGEMC_Flux->SetBranchAddress("id",&fFluxHit_id,&bFluxHit_id);
	TGEMC_Flux->SetBranchAddress("pid",&fFluxHit_pid,&bFluxHit_pid);
	TGEMC_Flux->SetBranchAddress("tid",&fFluxHit_tid,&bFluxHit_tid);
	TGEMC_Flux->SetBranchAddress("trackE",&fFluxHit_trackE,&bFluxHit_trackE);
	TGEMC_Flux->SetBranchAddress("avg_x",&fFluxHit_Avg_x,&bFluxHit_Avg_x);
	TGEMC_Flux->SetBranchAddress("avg_y",&fFluxHit_Avg_y,&bFluxHit_Avg_y);
	TGEMC_Flux->SetBranchAddress("avg_z",&fFluxHit_Avg_z,&bFluxHit_Avg_z);
	TGEMC_Flux->SetBranchAddress("px",&fFluxHit_Px,&bFluxHit_Px);
	TGEMC_Flux->SetBranchAddress("py",&fFluxHit_Py,&bFluxHit_Py);
	TGEMC_Flux->SetBranchAddress("pz",&fFluxHit_Pz,&bFluxHit_Pz);

	TGEMC_Flux->SetBranchAddress("vx",&fFluxHit_vx,&bFluxHit_vx);
	TGEMC_Flux->SetBranchAddress("vy",&fFluxHit_vy,&bFluxHit_vy);
	TGEMC_Flux->SetBranchAddress("vz",&fFluxHit_vz,&bFluxHit_vz);
	TGEMC_Flux->SetBranchAddress("mvx",&fFluxHit_mvx,&bFluxHit_mvx);
	TGEMC_Flux->SetBranchAddress("mvy",&fFluxHit_mvy,&bFluxHit_mvy);
	TGEMC_Flux->SetBranchAddress("mvz",&fFluxHit_mvz,&bFluxHit_mvz);
	TGEMC_Flux->SetBranchAddress("avg_t",&fFluxHit_Avg_t,&bFluxHit_t);  

};

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
	ifstream infile;	
	infile.open("map_FAEC_ANL_20130628_group.txt");	

	//load x,y locations and block id
	Double_t total_module=0;
	Int_t counter_id=0;
	for(int i=0;i<54;i++){ //54 rows
		INPUT_file>>num_module_in_row[i];
		infile>>num_module_in_row_group[i];
		//std::cout<<num_module_in_row[i]<<" ";
		num_module_in_row[i]=num_module_in_row[i]-1;  //first one is y coordinate
		total_module+=num_module_in_row[i];
		Double_t tmp_y;
		INPUT_file>>tmp_y;
		double tmp_y_group;
		infile>>tmp_y_group;
		y_bak[i]=tmp_y;       //make a backup in order to judge which row a certain particle hits the EC
		for(Int_t j=0;j<num_module_in_row[i];j++){
			INPUT_file>>x[i][j];
			infile>>group_x[i][j];
			y[i][j]=tmp_y;     // in each row , y coordinate is the same
			counter_id++;
			id[i][j]=counter_id;
			//update the pair array
			block_map[counter_id].first = i;
			block_map[counter_id].second = j;
			//printf("DEBUG EC map update : ecal id %i [%i,%i],[%f,%f] \n",counter_id,j,i,x[i][j],y[i][j]);//vec_ecalBlock.Px(),vec_ecalBlock.Py()
			//sector
			TVector2 vec(x[i][j],y[i][j]);
			Double_t phi_module=vec.Phi();
			for(int k=0;k<30;k++){   //30 sectors
				if(phi_module>=k*12.0/180.0*3.141592 && phi_module<(k+1.0)*12.0/180.0*3.141592){  //sector k
					sector[i][j]=k+1;
					// 				    cout << sector[i][j] << endl;				
				}
			}

		}
		//cout<<num_module_in_row[i]<<"	"<<y[i][0]<<endl;
		//std::cout<<" "<<std::endl;
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
	Int_t hit_around_idx[6]={100000};      
	Int_t hit_around_idy[6]={100000};
	Int_t tmp_idx[15]={hit_idx-7, hit_idx-6, hit_idx-5, hit_idx-4, hit_idx-3, hit_idx-2, hit_idx-1, hit_idx, hit_idx+1, hit_idx+2, hit_idx+3, hit_idx+4, hit_idx+5, hit_idx+6, hit_idx+7};
	Int_t tmp_idy[3]={hit_idy-1, hit_idy, hit_idy+1};
	Int_t label=0;
	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<3;j++){  
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



Int_t GetECALlargeCluser(Int_t block_id,Int_t *cluster_edep_largeblockid, Int_t *cluster_edep_largegroupid){
	Int_t hit_idx = block_map[block_id].second;
	Int_t hit_idy = block_map[block_id].first;
	Int_t hit_around_idx[40]={100000};      
	Int_t hit_around_idy[40]={100000};
	Int_t tmp_idx[21]={hit_idx-10,hit_idx-9,hit_idx-8,hit_idx-7, hit_idx-6, hit_idx-5, hit_idx-4, hit_idx-3, hit_idx-2, hit_idx-1, hit_idx, hit_idx+1, hit_idx+2, hit_idx+3, hit_idx+4, hit_idx+5, hit_idx+6, hit_idx+7,hit_idx+8,hit_idx+9,hit_idx+10};
	int tmp_idy[7]={hit_idy-3,hit_idy-2,hit_idy-1, hit_idy, hit_idy+1,hit_idy+2,hit_idy+3}; 
	Int_t label=0;
	for(Int_t i=0;i<21;i++){
		for(Int_t j=0;j<7;j++){  
			if(tmp_idy[j]>=0 && tmp_idy[j]<54 && tmp_idx[i]>=0 && tmp_idx[i]<num_module_in_row[tmp_idy[j]] /*&& (tmp_idx[i]!=hit_idx || tmp_idy[j]!=hit_idy)*/ ){ // in range
				if(sqrt( pow( (x[tmp_idy[j]][tmp_idx[i]]-x[hit_idy][hit_idx]),2 )+ pow((y[tmp_idy[j]][tmp_idx[i]]-y[hit_idy][hit_idx]),2 ))<30.0){
					hit_around_idx[label]=tmp_idx[i];
					hit_around_idy[label]=tmp_idy[j];
					label++;
				}
			}
		}
	}
	for(int l=0;l<label;l++){
		cluster_edep_largeblockid[l]=id[hit_around_idy[l]][hit_around_idx[l]];
		cluster_edep_largegroupid[l] = group_x[hit_around_idy[l]][hit_around_idx[l]];

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
	//TCanvas *c2 = new TCanvas("GetEnergyRes2", "GetEnergyRes2", 800 * 1.5, 600 * 1);//DEBUG
	for ( int i = 0; i < Difference->GetNbinsX() / 4; i++ )
	{
		//pt[cnt] = Difference->GetXaxis()->GetBinCenter(i * 10 + 10 / 2);//in GeV

		TH1D * proj = Difference->ProjectionY("_py", i * 10 + 1, (i + 1) * 10);

		if ( proj->GetSum() > 100 )//check for entries
		{
			pt[cnt] = Difference->GetXaxis()->GetBinCenter(i * 10 + 10 / 2);//in GeV
			Double_t fit_range_up=1.2;
			Double_t fit_range_low=1.5;   

			//with no fits
			TF1 * gaus = new TF1("gausEoverP", "gaus", fit_range_low, fit_range_up);
			proj -> Fit(gaus, "R0");
			resolution[0][cnt] = gaus->GetParameter(2)/gaus->GetParameter(1);//RMS[cnt];
			resolution[1][cnt] = proj->GetRMS()/proj->GetMean();
			resolution_error[0][cnt] = TMath::Sqrt(TMath::Power(gaus->GetParError(1)/gaus->GetParameter(1),2)+TMath::Power(gaus->GetParError(2)/gaus->GetParameter(2),2))*resolution[0][cnt];//RMSErr[cnt];
			resolution_error[1][cnt] = TMath::Sqrt(TMath::Power(proj->GetRMSError(),2)+TMath::Power(proj->GetMeanError(),2));//proj->GetRMSError()/proj->GetMean();
			cnt++;
			//}


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
	for (Int_t i =0; i<5;i++){
		if (radius>trig_low_R[i] && radius<=trig_high_R[i]){
			Ri=i;
			break;
		}
	}
	return Ri;
};

Double_t GetThreshold6p1(Double_t radius){
	Double_t thresh=0;
	if (radius<trig_low_R[0] || radius>trig_high_R[4]){
		thresh=11000;//set a high number so that R below the ecal limit will be rejected
		//thresh=0;
		return thresh;
	} else {
		for(Int_t i=0;i<5;i++){
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
	//Load the background in to memory
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
	const Double_t DEG=180./3.1415926;
	Int_t counttracks=0;
	//double detlaT = 6.13711543300417939e-05;// PVDIS
	double detlaT = 1.51920503039169666e-02; //SIDIS
	double event_time=0.0;
	double time_windows=0;
	Bool_t processWindow=kFALSE;//when all the events in a single window are filled to the ecal, start cluster processing
	filebkg[0]=new TFile("/cache/halla/solid/sim/solid_gemc/SIDIS_He3_JLAB_VERSION_1.3/pass7/background_solid_SIDIS_He3_BeamOnTarget_0.999e10_skim.root");
	filebkg[1]=new TFile("/cache/halla/solid/sim/solid_gemc/PVDIS_LD2_JLAB_VERSION_1.3/pass5/solid_PVDIS_LD2_BeamOnTarget_2e9_skim_2.root");
	filebkg[2]=new TFile("/cache/halla/solid/sim/solid_gemc/PVDIS_LD2_JLAB_VERSION_1.3/pass5/solid_PVDIS_LD2_BeamOnTarget_2e9_skim_3.root");
	filebkg[3]=new TFile("/cache/halla/solid/sim/solid_gemc/PVDIS_LD2_JLAB_VERSION_1.3/pass5/solid_PVDIS_LD2_BeamOnTarget_2e9_skim_4.root");
	filebkg[4]=new TFile("/cache/halla/solid/sim/solid_gemc/PVDIS_LD2_JLAB_VERSION_1.3/pass5/solid_PVDIS_LD2_BeamOnTarget_2e9_skim_5.root");

	//event info
	total_ecalShMap.clear();
	total_ecalPSMap.clear();
	total_ecalShMapbkg.clear();
	total_ecalPSMapbkg.clear();
	double pid_gen_bkg=0;
	double bkg_number=0;
	for(int n=0;n<1;n++){
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
		TGEMC_EM_bkg[n]->SetBranchAddress("pid",&fECHit_pid);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_z",&fECHit_Avg_z);
		TGEMC_EM_bkg[n]->SetBranchAddress("mpid",&fECHit_mpid);
		TGEMC_EM_bkg[n]->SetBranchAddress("tid",&fECHit_tid);
		TGEMC_EM_bkg[n]->SetBranchAddress("mtid",&fECHit_mtid); 
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_lx",&fECHit_Avg_lx);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_ly",&fECHit_Avg_ly);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_lz",&fECHit_Avg_lz);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_x",&fECHit_Avg_x);
		TGEMC_EM_bkg[n]->SetBranchAddress("avg_y",&fECHit_Avg_y);

		TGEMC_EM_bkg[n]->SetBranchAddress("px",&fECHit_Px);
		TGEMC_EM_bkg[n]->SetBranchAddress("py",&fECHit_Py);
		TGEMC_EM_bkg[n]->SetBranchAddress("pz",&fECHit_Pz);


		nentries[n] = (Int_t)TGEMC_Generate_bkg[n]->GetEntries();
		printf("Merged Background Entries = %i \n",nentries[n]);
		Double_t totEdep_shower=0;
		Double_t totEdep_preshower=0;
		Fillecal=kFALSE;  

		processWindow=kFALSE;//when all the events in a single window are filled to the ecal, start cluster processing
		Double_t pf_flux_bkg;
		Double_t pf_EC_bkg;
		for (int i=0; i<nentries[n]; i++) {
			TGEMC_Flux_bkg[n]->GetEntry(i); 
			TGEMC_EM_bkg[n]->GetEntry(i); 
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
					if (FluxHit_detector_ID==3 && FluxHit_subdetector_ID == 1 && FluxHit_subsubdetector_ID == 1 && pf_flux_bkg>0 && pf_flux_bkg<=12000 && fFluxHit_Px->at(j)>0){
						if(fFluxHit_pid->at(j)==11 && pf_flux_bkg>1000){
							Fillecal=kFALSE;
						}else 
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
						//        cout<<"EC_pid="<<fECHit_pid->at(j)<<endl;
						pf_EC_bkg=TMath::Sqrt(TMath::Power(fECHit_Px->at(j),2)+TMath::Power(fECHit_Py->at(j),2)+TMath::Power(fECHit_Pz->at(j),2));
						if(fECHit_pid->at(j)= 11 && pf_EC_bkg>1000){
							ecalShMapbkg[ECHit_component_ID]=0;
							ecalPSMapbkg[ECHit_component_ID]=0;
						}else{
							if (ECHit_detector_ID==3 && ECHit_subdetector_ID == 1 && ECHit_subsubdetector_ID == 0){//shower 
								if (fECHit_Avg_z->at(j)>2000){
									totEdep_shower+=fECHit_totEdep->at(j);
									//cout<<"ECHit_component_ID="<<ECHit_component_ID<<"totEdep_shower="<<totEdep_shower<<endl;
									if (ecalShMapbkg.count(ECHit_component_ID)){
										ecalShMapbkg[ECHit_component_ID]+=fECHit_totEdep->at(j);//5794;//106;// MeV
									} else {
										ecalShMapbkg[ECHit_component_ID]=fECHit_totEdep->at(j);//5794;//106;// MeV
									}
								}
								//cout<<"totEdep_shower="<<totEdep_shower<<endl;
							}
							//cout<<"event="<<i<<"E_preshower="<<totEdep_shower<<endl;
							if (ECHit_detector_ID==3 && ECHit_subdetector_ID == 1 && ECHit_subsubdetector_ID == 1){//pershower
								if (fECHit_Avg_z->at(j)>2000){
									totEdep_preshower+=fECHit_totEdep->at(j);

									if (ecalPSMapbkg.count(ECHit_component_ID)){
										ecalPSMapbkg[ECHit_component_ID]+=fECHit_totEdep->at(j);// MeV
									} else {
										ecalPSMapbkg[ECHit_component_ID]=fECHit_totEdep->at(j);// MeV
									}
								}//end preshower
							}
						}//end EChit loop
					}


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

			if (i%1000000==0)
				printf("Background Event %d \n",i);

		}//end of nentries
	}

	// cout<<"total_ecalPSMap.size()="<<total_ecalPSMap.size()<<endl;
	printf("Total no 6p1 max Background Trigger windows (@ %3.0f ns) available %d \n",TriggerWindow, total_ecalPSMapbkg.size());

	//end of background. It's in the memory
};

void getBkgEdepPS(Int_t rndnum){//this routine will fill ecalPSMap with background signals
	ecalPSMapbkg = total_ecalPSMapbkg[rndnum];
};
void getBkgEdepSh(Int_t rndnum){//this routine will fill ecalShMap with background signals
	ecalShMapbkg = total_ecalShMapbkg[rndnum];
};   
