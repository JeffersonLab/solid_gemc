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
const double PEthresh_hgc=2; //hgc pe threshold for each pmt
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

int analysis_hgc(string inputfile_name,string runmode="trigger", bool Is_tellorig=false,string filetype=""){

// gStyle->SetOptStat(11111111);
gStyle->SetOptStat(0);

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
//     return 0;
  }
  else{
    cout << "this file has no filenum, please check if you need filenum for addtional normalization" << endl;      
  }
}

bool Is_SIDIS_He3=false,Is_SIDIS_NH3=false,Is_JPsi_LH2=false,Is_C=false;
if(inputfile_name.find("SIDIS_He3",0) != string::npos) {
  Is_SIDIS_He3=true;
  cout << "SIDIS_He3 setup" << endl;  
}
else if(inputfile_name.find("SIDIS_NH3",0) != string::npos) {
  Is_SIDIS_NH3=true;
  cout << "SIDIS_NH3 setup" << endl;  
}
else if(inputfile_name.find("JPsi_LH2",0) != string::npos) {
  Is_JPsi_LH2=true;
  cout << "JPsi_LH2 setup" << endl;  
}
else if(inputfile_name.find("_C",0) != string::npos) {
  Is_C=true;
  cout << "Carbon setup" << endl;  
}
else {
    cout << "Not SIDIS_He3 or SIDIS_NH3 or JPsi_LH2 setup" << endl;    
//     return 0;
}

//Cherenkov sensor for 30 sectors
const int ch_lgc=270;
// const int ch_hgc=480;    //use pmt readout
// const int ch_hgc=1920;		//use quad readout
const int ch_hgc=30720;		//use pixel readout

int sensor_hgc = ch_hgc/30;
int sensor_trans_hgc = sqrt(sensor_hgc);

Float_t npe_lgc[ch_lgc],npe_hgc[ch_hgc];
Float_t npe_lgc_total=0;
Float_t npe_hgc_total=0;
		  
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

char textfile_name[1000];
sprintf(textfile_name,"%s_pixel.ml",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());
ofstream textfile(textfile_name,std::ofstream::trunc);

char matrixfile_name[1000];
sprintf(matrixfile_name,"%s_pixel_matrix.ml",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());
ofstream matrixfile(matrixfile_name,std::ofstream::trunc);

TFile *file=new TFile(inputfile_name.c_str());

// 	TString background_inputfile_name="parametrized_lgc.root";      //h_pe is here	
// 	TFile *background_file=new TFile(background_inputfile_name);
// 	TH1F *h_pe=(TH1F*)background_file->Get("h_pe");

std::size_t found = inputfile_name.rfind("cache");
if (found!=std::string::npos)  inputfile_name.replace(found,5,"work");

char the_filename[500];
sprintf(the_filename, "%s",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());

// char outputfile_name[200];
// sprintf(outputfile_name, "%s_output.root",the_filename);
// TFile *outputfile=new TFile(outputfile_name, "recreate");

// prepare for outputs
// define histograms, output txt files etc...
	
	TH2F *hpattern_hgc=new TH2F("pattern_hgc","N_{p.e.};x;y",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
	
	TH1F *hhit_hgc=new TH1F("hit_hgc","hit;sensor;N_{p.e.} rate(kHz)",sensor_hgc,0,sensor_hgc);
	TH1F *hocc_hgc=new TH1F("occ_hgc","occupancy;sensor;event rate(kHz)",sensor_hgc,0,sensor_hgc);
	TH2F *hhit_hgc_2D=new TH2F("hit_hgc_2D","hit,N_{p.e.} rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
	TH2F *hocc_hgc_2D=new TH2F("occ_hgc_2D","occupancy,event rate(kHz);sensor_lx;sensor_ly",sensor_trans_hgc,0,sensor_trans_hgc,sensor_trans_hgc,0,sensor_trans_hgc);
	
	TH1F *hnpe_hgc=new TH1F("npe_hgc",";Npe;",400,0,200);
	TH1F *hnpe_hgc_mom[10],*hnpe_hgc_mom_theta[10][10];
	for(int i=0;i<10;i++){
	  char hstname[200]; 
	  sprintf(hstname,"npe_hgc_mom_%i",i);
	  hnpe_hgc_mom[i]=new TH1F(hstname,";Npe;",400,0,200);
	  for(int j=0;j<10;j++){
	    sprintf(hstname,"npe_hgc_mom_theta_%i_%i",i,j);
	    hnpe_hgc_mom_theta[i][j]=new TH1F(hstname,";Npe;",400,0,200);
	  }	  
	}	
	
	TH2F *hvertex_rz=new TH2F("hvertex_rz","hvertex_rz",1800,-400,500,600,0,300);
	
	TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","generated events;vertex #theta (deg);vertex P (GeV)",100,0,50,110,0,11);  
	TH2F *hgen_ThetaPhi=new TH2F("gen_ThetaPhi","generated events;vertex #theta (deg);vertex #phi (deg)",100,0,50,360,-180,180);     
	TH2F *hgen_PhiP=new TH2F("gen_PhiP","generated events;vertex #phi (deg);vertex P (GeV)",360,-180,180,110,0,11);	
	
	TH3F *hgen_ThetaPhiP=new TH3F("gen_ThetaPhiP","gen_ThetaPhiP",50,0,50,180,-180,180,55,0,11);   
	
	const int n=16;
	char *detname[n]={"GEM 1","GEM 2","GEM 3","GEM 4","GEM 5","GEM 6","LGC","HGC","FASPD","LASPD","FAEC","LAEC","FAMRPC","FAMUON","LAMUON","new"};	
// 	const int n=12;
// 	char *detname[n]={"GEM 1","GEM 2","GEM 3","GEM 4","GEM 5","GEM 6","LGC","HGC","FASPD","LASPD","FAEC","LAEC"};		
	TH2F *hhit_xy[n],*hhit_PhiR[n];
	TH1F *hhit_loss[n];
	TH2F *hhit_momloss[n];
	for(int i=0;i<n;i++){
	  char hstname[100];
	  sprintf(hstname,"hit_xy_%i",i);
	  hhit_xy[i]=new TH2F(hstname,detname[i],600,-300,300,600,-300,300);        
	  sprintf(hstname,"hit_PhiR_%i",i);
	  hhit_PhiR[i]=new TH2F(hstname,detname[i],360,-180,180,300,0,300);
	  sprintf(hstname,"hit_loss_%i",i);
	  hhit_loss[i]=new TH1F(hstname,detname[i],55,0,11);
	  sprintf(hstname,"hit_momloss_%i",i);
	  hhit_momloss[i]=new TH2F(hstname,detname[i],55,0,11,200,0,1);	  
	}
	
	TH1F *hhit_mom_hgc_e=new TH1F("hit_mom_hgc_e","hit_mom_hgc_e",100,0,0.05);		
	//-------------------------
	//   get trees in the real data file
	//-------------------------
	
	//---header tree
	TTree *tree_header = (TTree*) file->Get("header");
	vector <int> *evn=0,*evn_type=0;
	vector <double> *beamPol=0;
	vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
	tree_header->SetBranchAddress("evn",&evn);      // event number 
	tree_header->SetBranchAddress("evn_type",&evn_type);  // evn_type==-1 for simulated events
	tree_header->SetBranchAddress("beamPol",&beamPol);   //beam polarization
	tree_header->SetBranchAddress("var1",&var1);     // W+ rate
	tree_header->SetBranchAddress("var2",&var2);     // W- rate
	tree_header->SetBranchAddress("var3",&var3);     // target pol
	tree_header->SetBranchAddress("var4",&var4);     //x
	tree_header->SetBranchAddress("var5",&var5);     //y
	tree_header->SetBranchAddress("var6",&var6);     //w
	tree_header->SetBranchAddress("var7",&var7);     //Q2
	tree_header->SetBranchAddress("var8",&var8);     //rate, Hz, should check the input file of the simulation

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
	vector<int> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
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

	//---hgc
	//information recorded by hgc
	TTree* tree_solid_hgc= (TTree*) file->Get("solid_hgc");
	setup_tree_solid_hgc(tree_solid_hgc);	
	
	TRandom3 rand;
	rand.SetSeed(0);
	
// 	long int N_events = (long int)tree_header->GetEntries();
	long int N_events = (long int)tree_generated->GetEntries();	

	cout << "total number of events : " << N_events << endl;	

	//----------------------------
	//      loop trees
	//---------------------------
	for(int loop_id=1;loop_id<=loop_time;loop_id++){
		cout<<"loop.....  "<<loop_id<<endl;
	
	for(long int i=0;i<N_events;i++){	  		
// 	for(long int i=0;i<N_events/10;i++){	  
// 	for(long int i=N_events/2;i<N_events;i++){	  		
// 	for(long int i=520;i<521;i++){  //pip event
// 	for(long int i=5289;i<5290;i++){	  // background event			  
// 			cout<<"event " << i<<endl;
// 		if (i%1000==0) cout<<i<<"\r";
// 		if (i%1000==0) cout<<i<<"\n";
		
		//---
		//---header tree
		//---
		tree_header->GetEntry(i);
		double rate=var8->at(0);
		if (filemode=="BeamOnTargetEM" || filemode=="BeamOnTarget") {
		  if(Is_SIDIS_He3) rate=30e-6/1.6e-19/event_actual/loop_time*add_norm;
		  else if(Is_SIDIS_NH3) rate=100e-9/1.6e-19/event_actual/loop_time*add_norm;
		  else if(Is_C) rate=30e-6/1.6e-19/event_actual/loop_time*add_norm;		  		  
		  else if(Is_JPsi_LH2) rate=3e-6/1.6e-19/event_actual/loop_time*add_norm; 	  
		  else {
		    cout << "Not SIDIS_He3 or SIDIS_NH3  or JPsi_LH2 setup" << endl;    
		    return 0;
		  }
		}
		else if (filemode=="even") rate=1;
		else if (filemode=="rate") rate=rate/filenum/loop_time*add_norm;
		else {
		    cout << "Not right filemode" << endl;    
// 		    return 0; 
		}
		
		double x=var4->at(0);	
		double y=var5->at(0);
		double W=var6->at(0);		
		double Q2=var7->at(0);		
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
// 			cout<<"event " << i<<endl;
		
// 		if (vz_gen<-3550 || -3450<vz_gen) continue;

		hgen_ThetaP->Fill(theta_gen,p_gen/1e3,rate);
		hgen_ThetaPhi->Fill(theta_gen,phi_gen,rate);                  		
		hgen_PhiP->Fill(phi_gen,p_gen/1e3,rate);                  				
		hgen_ThetaPhiP->Fill(theta_gen,phi_gen,p_gen/1e3,rate);                  			

		
		
		bool Is_ok=false;
		
		//-----------------------	
		//--- hgc 
		//-----------------------
		tree_solid_hgc->GetEntry(i);

		int hit_hgc[ch_hgc]={0};
		int trigger_hgc[30]={0};
		int ntrigsecs_hgc=0;
		
		
		if(solid_hgc_hitn->size()>0) {
		  Is_ok=true; 
		  textfile << "event" << "\t" <<  i << endl;
		  textfile << "cher" << "\t" <<  solid_hgc_hitn->size() << endl;		  
		}
		process_tree_solid_hgc(tree_solid_hgc,hit_hgc,trigger_hgc,ntrigsecs_hgc,PMTthresh_hgc,PEthresh_hgc,ch_hgc,textfile);
		
// 		for(int index=0;index<ch_hgc;index++){
// 		  if (hit_hgc[index]>0) {
// 		    cout << i << "\t" << index << "\t" << hit_hgc[index] << endl;		    
// 		    textfile << i << "\t" << index << "\t" << hit_hgc[index] << endl;
// 		  }
// 		}

 		
// 		if(ntrigsecs_hgc){
		if(true){
  
		for(int index=0;index<ch_hgc;index++){    
// 		  cout << hit_hgc[index] << endl;
		      matrixfile << hit_hgc[index] << "\t";
		      
		      int pmt_sec=index/sensor_hgc;		  
		      int pmt_hgc=index%sensor_hgc;
// 		      int pmt_x=pmt_hgc%sensor_trans_hgc,pmt_y=sensor_trans_hgc-1-pmt_hgc/sensor_trans_hgc;
		      int pmt_x=pmt_hgc%sensor_trans_hgc,pmt_y=pmt_hgc/sensor_trans_hgc;		      
		  
		      if(hit_hgc[index]>occ_threshold_hgc){
// 			cout << pmt_hgc << " " << pmt_x << " " << pmt_y << endl;
			
			hpattern_hgc->Fill(pmt_x,pmt_y,hit_hgc[index]);		
			
			hhit_hgc->Fill(pmt_hgc,hit_hgc[index]*rate/1e3);	    
			hhit_hgc_2D->Fill(pmt_x,pmt_y,hit_hgc[index]*rate/1e3);

			hocc_hgc->Fill(pmt_hgc,rate/1e3);
			hocc_hgc_2D->Fill(pmt_x,pmt_y,rate/1e3);
		      }	

		}
		matrixfile << endl;		

		} //pass hgc trigger in offline		
		
// 		if (1<=sec_hgc && sec_hgc<=30){
// 		    int this_sec=sec_hgc;
// 		    int last_sec=-1,next_sec=-1;
// 		    if (this_sec==1) {last_sec=30;next_sec=2;}
// 		    else if (this_sec==30) {last_sec=29;next_sec=1;}
// 		    else {last_sec=this_sec-1;next_sec=this_sec+1;}		    		  
// 		  double npe_hgc=0;
// 		  for(int pmt_id=0;pmt_id<16;pmt_id++) npe_hgc += hit_hgc[16*(sec_hgc-1)+pmt_id];
// 		  for(int pmt_id=0;pmt_id<16;pmt_id++) npe_hgc += hit_hgc[16*(last_sec-1)+pmt_id];
// 		  for(int pmt_id=0;pmt_id<16;pmt_id++) npe_hgc += hit_hgc[16*(next_sec-1)+pmt_id];

		  npe_hgc_total=0;
		  for(int pmt_id=0;pmt_id<ch_hgc;pmt_id++) {
		    npe_hgc[pmt_id]=0;		    
		    if (hit_hgc[pmt_id]>0){
// 		      cout << "event " << i << " hgc pmt_id " << pmt_id << " npe " << hit_hgc[pmt_id] << endl;
		      npe_hgc_total += hit_hgc[pmt_id];  
		      
		      npe_hgc[pmt_id]=hit_hgc[pmt_id];
// 		      npe_hgc[pmt_id]=1;		      
		    }
		  }
		  hnpe_hgc->Fill(npe_hgc_total);

// 		  if (npe_hgc_total>0) cout<<"event " << i<< " npe_hgc_total " << npe_hgc_total << endl;		  
		  
// 		  for(int pmt_id=0;pmt_id<ch_hgc;pmt_id++) cout << npe_hgc[pmt_id] << "\t";
// 		  if (npe_hgc_total>0) mlptree_hgc->Fill();
		  
// 		if (-365< vz_gen/10. && vz_gen/10. < -335){
		if (7<=theta_gen && theta_gen <=15){ //has to cut to prevent int(-1)=0 !!!!!
		if (2500<=p_gen && theta_gen <=7500){	//has to cut to prevent int(-1)=0	  
		  
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
// 		    cout<<"event " << i<<endl;
// 		    cout << "npe_hgc " << npe_hgc << endl;		    
// 		    cout << "sec_hgc " << sec_hgc << endl;
// 		for (Int_t j=0;j<flux_hitn->size();j++) {
// 	          if (flux_id->at(j)==2210000 && (flux_tid->at(j)==1 || abs(flux_pid->at(j))==221 || abs(flux_pid->at(j))==13))   cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
// 		}		  
// 		  }						
		 		
		//---	
		//---flux tree
		//---
		tree_flux->GetEntry(i);
		
		if(Is_ok && flux_hitn->size()>0) {
		  textfile << "flux" << "\t" << flux_hitn->size() << endl;		  
		}		
		
		//check on EC and other by flux
		for (Int_t j=0;j<flux_hitn->size();j++) {
// 	          cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  
		  

// 		  if (Is_ok && flux_id->at(j)==10) {  // for cc_pro front window
		  if (Is_ok && flux_id->at(j)==2210000) {  // for solid hgc front virtual plane		    
		    textfile << flux_hitn->at(j) << "\t" <<  flux_pid->at(j) << "\t" <<  flux_mpid->at(j) << "\t" <<  flux_tid->at(j) << "\t" <<  flux_mtid->at(j) << "\t" <<  flux_trackE->at(j) << "\t" <<  flux_totEdep->at(j) << "\t" <<  flux_avg_lx->at(j) << "\t" <<  flux_avg_ly->at(j) << "\t" <<  flux_avg_lz->at(j) << "\t" <<  flux_px->at(j) << "\t" <<  flux_py->at(j) << "\t" <<  flux_pz->at(j) << "\t" <<   flux_avg_t->at(j) << endl;  

		  }
		  
		}			
		

	} //end event loop
	
	} //end loop time


//do outputs
textfile.close();
matrixfile.close();

TCanvas *c = new TCanvas("c", "c",1900,900);
c->Divide(2,1);
c->cd(1);
hhit_hgc_2D->Draw("colz");
c->cd(2);
hocc_hgc_2D->Draw("colz");
// c->cd(3);
// hnpe_hgc->Draw();

// outputfile->Write();	
// outputfile->Flush();

// mlpfile_lgc->Write();	
// mlpfile_lgc->Flush();
// mlpfile_hgc->Write();	
// mlpfile_hgc->Flush();

}
