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

int analysis(string inputfile_name,string runmode="trigger", bool Is_tellorig=false,string filetype="",bool Is_new=true){

// gStyle->SetOptStat(11111111);
gStyle->SetOptStat("ioue");
// gStyle->SetOptStat(0);

// gStyle->SetPalette(57);

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
// if (inputfile_name.find("BeamOnTargetEM",0) != string::npos) {
//   filemode="BeamOnTargetEM";
//   cout << "this is a BeamOnTargetEM file" << endl;  
//   
//   event_actual=atof(inputfile_name.substr(inputfile_name.find("BeamOnTargetEM",0)+15,inputfile_name.find("_")).c_str());
//   cout << "event_actual " << event_actual <<  endl;  
// }
// else if (inputfile_name.find("BeamOnTarget",0) != string::npos) {
//   filemode="BeamOnTarget";
//   cout << "this is a BeamOnTarget file" << endl;  
//   
//   event_actual=atof(inputfile_name.substr(inputfile_name.find("BeamOnTarget",0)+13,inputfile_name.find("_")).c_str());
//   cout << "event_actual " << event_actual <<  endl;  
// }
// else if (inputfile_name.find("even",0) != string::npos) {
//   filemode="even";
//   cout << "this is a evenly distributed file" << endl;  
// }
// else {
//   filemode="rate";  
//   cout << "this is rate dependent file" << endl;  
// }

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
	
	TH2F *hvertex_rz=new TH2F("hvertex_rz","hvertex_rz",1800,-400,500,600,0,300);
	
	TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","generated events;vertex #theta (deg);vertex P (GeV)",100,0,50,110,0,11);  
	TH2F *hgen_ThetaPhi=new TH2F("gen_ThetaPhi","generated events;vertex #theta (deg);vertex #phi (deg)",100,0,50,360,-180,180);     
	TH2F *hgen_PhiP=new TH2F("gen_PhiP","generated events;vertex #phi (deg);vertex P (GeV)",360,-180,180,110,0,11);	
	
	TH3F *hgen_ThetaPhiP=new TH3F("gen_ThetaPhiP","gen_ThetaPhiP",50,0,50,180,-180,180,55,0,11);   
	
	const int n=10;
	string detname[n]={
	  "muon_FA_1_virt","muon_FA_2_virt","muon_FA_3_virt","muon_LA_1_virt","muon_barrel_1_virt",
	  "muon_FA_1_scint","muon_FA_2_scint","muon_FA_3_scint","muon_LA_1_scint","muon_barrel_1_scint",
	};
	TH2F *hhit_rz[n],*hhit_rz_orig[n];
	TH2F *hhit_xy[n],*hhit_xy_orig[n],*hhit_PhiR[n];
	TH1F *hhit_E[n],*hhit_E_mip[n],*hhit_E_photonele[n],*hhit_E_ele[n],*hhit_Edep[n];		  	
	for(int i=0;i<n;i++){
	  char hstname[100];
	  sprintf(hstname,"hit_rz_%i",i);
	  hhit_rz[i]=new TH2F(hstname,detname[i].c_str(),1200,-400,800,300,0,300);
	  sprintf(hstname,"hit_rz_orig_%i",i);
	  hhit_rz_orig[i]=new TH2F(hstname,detname[i].c_str(),1200,-400,800,300,0,300);
	  sprintf(hstname,"hit_xy_%i",i);
	  hhit_xy[i]=new TH2F(hstname,detname[i].c_str(),600,-300,300,600,-300,300);        
	  sprintf(hstname,"hit_xy_orig_%i",i);
	  hhit_xy_orig[i]=new TH2F(hstname,detname[i].c_str(),600,-300,300,600,-300,300);         
	  sprintf(hstname,"hit_PhiR_%i",i);
	  hhit_PhiR[i]=new TH2F(hstname,detname[i].c_str(),360,-180,180,300,0,300);
	  sprintf(hstname,"hit_E_%i",i);
	  hhit_E[i]=new TH1F(hstname,detname[i].c_str(),1100,0,11);
	  sprintf(hstname,"hit_E_mip_%i",i);
	  hhit_E_mip[i]=new TH1F(hstname,detname[i].c_str(),1100,0,11);
	  sprintf(hstname,"hit_E_photonele_%i",i);
	  hhit_E_photonele[i]=new TH1F(hstname,detname[i].c_str(),1100,0,11);
	  sprintf(hstname,"hit_E_ele_%i",i);
	  hhit_E_ele[i]=new TH1F(hstname,detname[i].c_str(),1100,0,11);	  
	  sprintf(hstname,"hit_Edep_%i",i);
	  hhit_Edep[i]=new TH1F(hstname,detname[i].c_str(),1000,0,0.1);	  
	}
	
        TH2F *hhit_EdepP_P_FAMD=new TH2F("hhit_EdepP_P_FAMD","hhit_EdepP_P_FAMD;P(GeV);Edep_FAMD/P(GeV)",110,0,22,100,0,0.1);
        TH2F *hhit_Edep_P_FAMD=new TH2F("hhit_Edep_P_FAMD","hhit_Edep_P_FAMD;P(GeV);Edep_FAMD(GeV)",110,0,22,100,0,0.1);
        TH1F *hhit_Edep_FAMD=new TH1F("hhit_Edep_FAMD","hhit_Edep_FAMD;Edep_FAMD(GeV);count",100,0,0.1);
        TH2F *hhit_Edep_hitr_FAMD=new TH2F("hhit_Edep_hitr_FAMD","hhit_Edep_hitr_FAMD;hitr_FAMD(cm);Edep_FAMD(GeV)",100,0,1e5,100,0,0.1);
        
	//-------------------------
	//   get trees in the data file
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
// 	for(long int i=0;i<N_events;i++){
	for(long int i=0;i<N_events/10;i++){	  
// 	for(long int i=0;i<1e4;i++){	  		
// 	for(long int i=520;i<521;i++){  //pip event
// 	for(long int i=5289;i<5290;i++){	  // background event			  
// 			cout<<"" << i<<endl;
		cout<<"event " << i << "\r";
// 		if (i%1000==0) cout<<i<<"\r";
// 		if (i%1000==0) cout<<i<<"\n";
		
		//---
		//---header tree
		//---
		double rate=1;
// 		cout << "rate " << rate << endl;
		
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
		
		//---	
		//---flux tree
		//---
		tree_flux->GetEntry(i);		  
		  	
// 		 cout << "flux_hitn  " << flux_hitn->size() << endl;
                
		bool Is_reachlast=false;
		for (Int_t j=0;j<flux_hitn->size();j++) {
//                     if(abs(flux_pid->at(j))==211 && flux_tid->at(j)!=1 && flux_id->at(j)==6103000) Is_reachlast=true;                                        
//                     if(abs(flux_pid->at(j))==211 && flux_id->at(j)==6103000) Is_reachlast=true;                                        
                    if(flux_tid->at(j)==1 && flux_id->at(j)==6103000) Is_reachlast=true;                                        
//                     if(flux_pid->at(j)==13 && flux_mtid->at(j)==1 && flux_id->at(j)==6103000)  Is_reachlast=true;
                }

		int count_layer=0;
		double Edep_FAMD=0,hit_r_FAMD=0;
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
		  if(flux_id->at(j)==6110000) hit_id=0; 
		  else if(flux_id->at(j)==6120000) hit_id=1;
		  else if(flux_id->at(j)==6130000) hit_id=2;
		  else if(flux_id->at(j)==6210000) hit_id=3;
		  else if(flux_id->at(j)==6310000) hit_id=4;
		  else if(flux_id->at(j)==6101000) hit_id=5; 
		  else if(flux_id->at(j)==6102000) hit_id=6;
		  else if(flux_id->at(j)==6103000) hit_id=7;
		  else if(flux_id->at(j)==6201000) hit_id=8;
		  else if(flux_id->at(j)==6301000) hit_id=9;
		  else {
//                       cout << "wrong flux_id " << flux_id->at(j) << endl;
                      continue;
                  }
// 		  if (hit_id==-1) {/cout << flux_id->at(j) << " " << flux_avg_z->at(j) << endl;
		  
		  hhit_xy[hit_id]->Fill(hit_x,hit_y,rate);
		  hhit_rz[hit_id]->Fill(hit_z,hit_r,rate);
		  
		  if (flux_tid->at(j) == 1){ // hit is from original particle
//                       cout << hit_x << " " << hit_y << endl;
		    hhit_xy_orig[hit_id]->Fill(hit_x,hit_y,rate);
		    hhit_rz_orig[hit_id]->Fill(hit_z,hit_r,rate);
		  }

		  double E=flux_trackE->at(j)/1e3,Edep=flux_totEdep->at(j)/1e3;		  
		  hhit_Edep[hit_id]->Fill(Edep,rate);
		  hhit_E[hit_id]->Fill(E,rate);
		  if (abs(flux_pid->at(j)) == 211 || flux_pid->at(j)==2212|| flux_pid->at(j)==13) hhit_E_mip[hit_id]->Fill(E,rate);
		  if (abs(flux_pid->at(j)) == 11 || flux_pid->at(j)==22) hhit_E_photonele[hit_id]->Fill(E,rate);
		  if (abs(flux_pid->at(j)) == 11) hhit_E_ele[hit_id]->Fill(E,rate);  
		  
// 		  if(flux_tid->at(j)==1 && (hit_id==5||hit_id==6||hit_id==7)) {count_layer++;Edep_FAMD += Edep;}
// 		  if(Is_reachlast && (hit_id==5||hit_id==6||hit_id==7)) {hit_r_FAMD += hit_r;Edep_FAMD += Edep;}
                  if(hit_id==5||hit_id==6||hit_id==7) {hit_r_FAMD += hit_r;Edep_FAMD += Edep;}
		}	// end of flux		
		
// 		if (Edep_FAMD>0){
		if (Is_reachlast && Edep_FAMD>0){
                    hhit_EdepP_P_FAMD->Fill(p_gen/1e3,Edep_FAMD/(p_gen/1e3));
                    hhit_Edep_P_FAMD->Fill(p_gen/1e3,Edep_FAMD);
                    hhit_Edep_FAMD->Fill(Edep_FAMD);
                    hhit_Edep_hitr_FAMD->Fill(hit_r_FAMD,Edep_FAMD);
		}
		

} //end loop
	

cout <<" sensor_good " << sensor_good << endl;
cout <<" event_good " << event_good << endl;
cout <<" event_trig_good " << event_trig_good << endl;

//do outputs

outputfile->Write();	
outputfile->Flush();

TCanvas *c_hit_xy_orig = new TCanvas("hit_xy_orig", "hit_xy_orig",1900,900);
c_hit_xy_orig->Divide(5,2);
for(int i=0;i<n;i++){
c_hit_xy_orig->cd(i+1);
hhit_xy_orig[i]->Draw("colz box");
}

TCanvas *c_hit_rz_orig = new TCanvas("hit_rz_orig", "hit_rz_orig",1900,900);
c_hit_rz_orig->Divide(5,2);
for(int i=0;i<n;i++){
c_hit_rz_orig->cd(i+1);
hhit_rz_orig[i]->Draw("colz");
}

TCanvas *c_Edep = new TCanvas("Edep", "Edep",1900,900);
c_Edep->Divide(5,2);
for(int i=0;i<5;i++){
c_Edep->cd(i+1);
gPad->SetLogy(1);
hhit_E[i]->Draw("colz");
}
for(int i=5;i<n;i++){
c_Edep->cd(i+1);
gPad->SetLogy(1);
hhit_Edep[i]->Draw("colz");
}

TCanvas *c_Edep_FAMD = new TCanvas("Edep_FAMD", "Edep_FAMD",1900,1000);
c_Edep_FAMD->Divide(2,2);
c_Edep_FAMD->cd(1);
hhit_EdepP_P_FAMD->Draw("colz");
c_Edep_FAMD->cd(2);
hhit_Edep_P_FAMD->Draw("colz");
c_Edep_FAMD->cd(3);
hhit_Edep_FAMD->Draw();
c_Edep_FAMD->cd(4);
hhit_Edep_hitr_FAMD->Draw("colz");
// c_Edep_FAMD->SaveAs("Edep_FAMD.png");

// exit(0);
}
