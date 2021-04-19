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
#include "TMath.h"
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

// //the wiser curve
// #include "ec_triggerfile/SIDIS_He3_JinHuang/SIDIS_He3_FAEC_electron_trigger_WiserJinHuang.C"
// #include "ec_triggerfile/SIDIS_He3_JinHuang/SIDIS_He3_LAEC_electron_trigger_WiserJinHuang.C"
// #include "ec_triggerfile/SIDIS_He3_JinHuang/SIDIS_He3_FAEC_par_trigger_WiserJinHuang.C"
// // EC efficiency table inside
// 
// #include "analysis_tree_solid_lgc.C"
// #include "analysis_tree_solid_hgc.C"
// #include "analysis_tree_solid_spd.C"
// #include "analysis_tree_solid_mrpc.C"
// #include "analysis_tree_solid_ec.C"


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
const double PEthresh=2; //lgc pe shreshold for each pmt
const double PMTthresh=2; //lgc pmt shreshold, at least 2pmts are fired in each sector
//spd 
// const double trigger_threshold_spd_FA=0.35;         //in MeV
const double trigger_threshold_spd_FA=0.5;         //in MeV
const double trigger_threshold_spd_LA=1.5;         //in MeV
//mrpc threshold
const double mrpc_block_threshold_FA=5;  //how many layers are required to be fired

//occupancy threshold
// double occ_threshold_lgc=0,occ_threshold_hgc=0; //in N_p.e.
// double occ_threshold_spd_FA=trigger_threshold_spd_FA/5.,occ_threshold_spd_LA=trigger_threshold_spd_LA/5.; //in MeV
// double occ_threshold_ec_preshower=0.4,occ_threshold_ec_shower=6; //in MeV

double occ_threshold_lgc=1,occ_threshold_hgc=1; //in N_p.e.
// double occ_threshold_spd_FA=trigger_threshold_spd_FA/2.,occ_threshold_spd_LA=trigger_threshold_spd_LA/2.; //in MeV
double occ_threshold_spd_FA=0.5,occ_threshold_spd_LA=3.; 
double occ_threshold_ec_preshower=0.8,occ_threshold_ec_shower=12; //in MeV

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

int analysis_rate(string inputfile_name,string par="e",string runmode="trigger", bool Is_tellorig=false,string filetype="",string outputdir="",bool Is_new=true){

// gStyle->SetOptStat(11111111);
// gStyle->SetOptStat(111);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.04,"xyz");   
  gStyle->SetTitleSize(0.07,"t");    
  gStyle->SetPaintTextFormat("4.1f");   
  
  
  string par_title,e_title,par_filename;
  if (par=="pi0") {par_title="#pi^{0} e^{-}"; e_title="e^{-}"; par_filename="pion_0";}  
  else if (par=="pim") {par_title="#pi^{-}"; e_title="e^{-}"; par_filename="pion_m";}
  else if (par=="pip") {par_title="#pi^{+}"; e_title="e^{+}"; par_filename="pion_p";}
  else if (par=="p") {    par_title="p";  e_title="e^{+}"; par_filename="proton";  }
  else if (par=="eany") {    par_title="any e^{-}";  e_title="e^{-}"; par_filename="eany";  }
  else if (par=="e") {    par_title="e^{-}";  e_title="e^{-}"; par_filename="e";  }  
  else {
     cout << "unknown par" << endl;
     return 0;
  }  

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

bool Is_SIDIS_He3=false,Is_SIDIS_NH3=false,Is_JPsi_LH2=false,Is_PVDIS_LD2=false;
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
else if(inputfile_name.find("PVDIS_LD2",0) != string::npos) {
  Is_PVDIS_LD2=true;
  cout << "PVDIS_LD2 setup" << endl;  
}
else {
    cout << "Not SIDIS_He3 or SIDIS_NH3 or JPsi_LH2 setup" << endl;    
    return 0;
}

TFile *file=new TFile(inputfile_name.c_str());

// 	TString background_inputfile_name="parametrized_lgc.root";      //h_pe is here	
// 	TFile *background_file=new TFile(background_inputfile_name);
// 	TH1F *h_pe=(TH1F*)background_file->Get("h_pe");

std::size_t found = inputfile_name.rfind("cache");
char the_filename[400];
if (outputdir=="original") {
  if (found!=std::string::npos) inputfile_name.replace(found,5,"work");
}
else {
  inputfile_name=inputfile_name.substr(inputfile_name.rfind("/")+1);
}
sprintf(the_filename, "%s",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());  

char outputfile_name[200];
sprintf(outputfile_name, "%s_output.root",the_filename);
TFile *outputfile=new TFile(outputfile_name, "recreate");
cout << "output file " << outputfile_name << endl;

   TH2F *hPTheta=new TH2F("hPTheta","hPTheta",20,5,45,20,0,10);
   TH2F *hQ2x=new TH2F("hQ2x","hQ2x",10,0,1,14,0,14);   

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
	if (!Is_new){	
	tree_header->SetBranchAddress("runNo",&runNo);  // run number	  
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
	
	
	int event_good=0,event_trig_good=0;
	
// 	long int N_events = (long int)tree_header->GetEntries();
	long int N_events = (long int)tree_generated->GetEntries();	

	cout << "total number of events : " << N_events << endl;	
	for(int loop_id=1;loop_id<=loop_time;loop_id++){
		cout<<"loop.....  "<<loop_id<<endl;

	for(long int i=0;i<N_events;i++){	  		
		if (i%1000==0) cout<<i<<"\r";
// 		if (i%1000==0) cout<<i<<"\n";
		
		//---
		//---header tree
		//---
		double rate=0;
		if (Is_new){
			tree_userHeader->GetEntry(i);
			if (par=="e")	rate=userVar006->at(0); // new eDIS and eAll generator
			else rate=userVar010->at(0); // bggen
// 			cout << "rate  " << rate   << endl;
		}
		else {
			tree_header->GetEntry(i);
			rate=var8->at(0); // old eDIS and eAll generator
		}

		if (filemode=="BeamOnTargetEM" || filemode=="BeamOnTarget") {
		  if(Is_SIDIS_He3) rate=15e-6/1.6e-19/event_actual/loop_time*add_norm;
		  else if(Is_SIDIS_NH3) rate=100e-9/1.6e-19/event_actual/loop_time*add_norm;
		  else if(Is_JPsi_LH2) rate=3e-6/1.6e-19/event_actual/loop_time*add_norm; 	  
		  else if(Is_PVDIS_LD2)  rate=50e6/1.6e-19/event_actual/loop_time*add_norm;
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
		
		double Q2_gen=0,x_gen=0;
		int tid_theone=0;
		int pid_gen=0;
		double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;      
	      //       cout << "gen_pid->size() " << gen_pid->size() << endl;        
		for (int j=0;j<gen_pid->size();j++) {
// 	            cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
		    pid_gen=TMath::Nint(gen_pid->at(j));
		    px_gen=gen_px->at(j)/1e3;
		    py_gen=gen_py->at(j)/1e3;
		    pz_gen=gen_pz->at(j)/1e3;
		    vx_gen=gen_vx->at(j)/1e1;
		    vy_gen=gen_vy->at(j)/1e1;
		    vz_gen=gen_vz->at(j)/1e1;
		    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
		    theta_gen=acos(pz_gen/p_gen)*DEG;
		    phi_gen=atan2(py_gen,px_gen)*DEG;
// 	            cout << "p_gen " << p_gen << endl; 

		    if ((par=="e" && pid_gen==11) || (par=="pi0" && pid_gen==111) || (par=="pim" && pid_gen==-211) || 
				(par=="pip" && pid_gen==211) || (par=="p" && pid_gen==2212)) {
		      Q2_gen=4*11*p_gen*sin(theta_gen/2/180*3.1416)*sin(theta_gen/2/180*3.1416);
		      x_gen=Q2_gen/2/0.938/(11-p_gen);
		      tid_theone=j+1;
// 		      cout << "pid_gen " << pid_gen<< " tid_theone " <<  tid_theone << endl;
		    }
		    
		}		
		  
		
// 		TVector3 *vec_p_gen(px_gen,py_gen,pz_gen),*vec_v_gen(vx_gen,vy_gen,vz_gen);
		
		// 		if (phi_gen<-1 || phi_gen>1) continue;
// 			cout<<"" << i<<endl;
		
// 		if (vz_gen<-3550 || -3450<vz_gen) continue;
                			
		//---	
		//---flux tree
		//---
		tree_flux->GetEntry(i);		  
		  		
		//check on EC and other by flux
// 		 cout << "flux_hitn  " << flux_hitn->size() << endl;
		double Eec=0,Eec_photonele=0,Eec_ele=0,Edepsc1=0,Edepsc2=0;
		
		for (Int_t j=0;j<flux_hitn->size();j++) {
// 	          cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  

		double hit_vr=sqrt(pow(flux_vx->at(j),2)+pow(flux_vy->at(j),2))/1e1; //mm to cm
		double hit_vy=flux_vy->at(j)/1e1,hit_vx=flux_vx->at(j)/1e1,hit_vz=flux_vz->at(j)/1e1;           //mm to cm		  
		double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))/1e1; //mm to cm
		double hit_y=flux_avg_y->at(j)/1e1,hit_x=flux_avg_x->at(j)/1e1,hit_z=flux_avg_z->at(j)/1e1;     //mm to cm		
		double hit_phi=atan2(hit_y,hit_x)*DEG;       //rad to  deg
		double hit_p=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j))/1e3;  //MeV to GeV
			
// 		TVector3 *vec_hit(hit_x,hit_y,hit_z), *vec_v(hit_vx,hit_vy,hit_vz),*vec_p(flux_px->at(j),flux_py->at(j),flux_pz->at(j));
// 		TVector3 *vec_path=vec_hit-vec_v_gen;
		  
		  if (TMath::Nint(flux_id->at(j))==3110000) {
			//apply cuts
            if(hit_p<0.3) continue;      //particle less than 0.3GeV will not pass pre-lead (1.12cm), do not event consider this particle
			if(hit_r > rout_cut_FA || hit_r < rin_cut_FA)  continue;
			
// 		    cout << " ec " << flux_tid->at(j) << " " <<  flux_mtid->at(j) << endl;
		    int this_pid = TMath::Nint(flux_pid->at(j));
		    int this_tid = TMath::Nint(flux_tid->at(j));
		    int this_otid = TMath::Nint(flux_otid->at(j));
		    //find the original particle and plot its vertex values assuming tracking can reconstruct those
		    if ((par=="e" && this_pid==11) || (par=="pim" && this_tid==tid_theone && this_pid==-211) || 
				(par=="pip" && this_tid==tid_theone && this_pid==211) || (par=="p" && this_tid==tid_theone && this_pid==2212) ) {
// 		    if ((par=="e" && flux_tid->at(j)==tid_theone && flux_pid->at(j)==11) || (par=="pim" && flux_tid->at(j)==tid_theone && flux_pid->at(j)==-211) || (par=="pip" && flux_tid->at(j)==tid_theone && flux_pid->at(j)==211) || (par=="p" && flux_tid->at(j)==tid_theone && flux_pid->at(j)==2212)){
		      hPTheta->Fill(theta_gen,p_gen,rate/1e3);
		      hQ2x->Fill(x_gen,Q2_gen,rate/1e3);
		    }
		    
// 		    for pi0, find the e- from pi0 originally
		    if (par=="pi0" && this_otid==tid_theone && this_pid==11) {
// 		    for pi0, find the e- from direct decay only pi0 ->e-e+
// 		    if (par=="pi0" && flux_mtid->at(j)==tid_theone && flux_pid->at(j)==11){
// 		      cout << " ec " << tid_theone << " " <<  flux_tid->at(j) << " " <<  flux_mtid->at(j) << " " << flux_otid->at(j) << endl;		      
		      
		      // assume this e- is scattered e- beam and use flux value instead of vertex values
		      double hit_theta=atan2(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j),flux_pz->at(j)*flux_pz->at(j))/3.1416*180;
		      double hit_Q2=4*11*hit_p*sin(hit_theta/2/180*3.1416)*sin(hit_theta/2/180*3.1416);
		      double hit_x=hit_Q2/2/0.938/(11-hit_p);
		      hPTheta->Fill(hit_theta,hit_p,rate/1e3);
		      hQ2x->Fill(hit_x,hit_Q2,rate/1e3);
		    }
		    
		    if (par=="eany" && this_pid==11) {
		      // assume this e- is scattered e- beam and use flux value instead of vertex values
		      double hit_theta=atan2(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j),flux_pz->at(j)*flux_pz->at(j))/3.1416*180;
		      double hit_Q2=4*11*hit_p*sin(hit_theta/2/180*3.1416)*sin(hit_theta/2/180*3.1416);
		      double hit_x=hit_Q2/2/0.938/(11-hit_p);
		      hPTheta->Fill(hit_theta,hit_p,rate/1e3);
		      hQ2x->Fill(hit_x,hit_Q2,rate/1e3);
		    }
		    
		  }
		  
		}

	}	
	}
	
outputfile->Write();	
outputfile->Flush();

   string type="sim_solid_PVDIS_LD2";
   
   TCanvas *cPTheta = new TCanvas("cPTheta","cPTheta",1800,1000);
   hPTheta->SetTitle(Form("%s rate (kHz);#theta (deg);P (GeV/c)",par_title.c_str()));       
//    hPTheta->GetXaxis()->SetRangeUser(min_theta,max_theta);   
//    hPTheta->GetYaxis()->SetRangeUser(min_p,max_p);     
   hPTheta->Draw("text colz");
//    hPTheta->SetMaximum(1e13);
//    hPTheta->SetMinimum(1e-3);
   cPTheta->SaveAs(Form("%s_rate_PTheta_%s.pdf",type.c_str(),par.c_str()));  
   cPTheta->SaveAs(Form("%s_rate_PTheta_%s.png",type.c_str(),par.c_str()));         
	
   TCanvas *cQ2x = new TCanvas("cQ2x","cQ2x",1800,1000);
   hQ2x->SetTitle(Form("%s rate (kHz);x;Q^{2} (GeV/c)^{2}",par_title.c_str()));       
//    hQ2x->GetXaxis()->SetRangeUser(min_theta,max_theta);   
//    hQ2x->GetYaxis()->SetRangeUser(min_p,max_p);     
   hQ2x->Draw("text colz");
//    hQ2x->SetMaximum(1e13);
//    hQ2x->SetMinimum(1e-3);
   cQ2x->SaveAs(Form("%s_rate_Q2x_%s.pdf",type.c_str(),par.c_str()));  
   cQ2x->SaveAs(Form("%s_rate_Q2x_%s.png",type.c_str(),par.c_str()));         
	
	
}
