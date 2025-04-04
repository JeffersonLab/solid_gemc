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

int analysis_physics(string inputfile_name,string evgen="twopeg",string runmode="trigger", bool Is_tellorig=false,string filetype="",bool Is_new=true){

// gStyle->SetOptStat(11111111);
// gStyle->SetOptStat("ioue");
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

// string filemode;
// double event_actual=1;
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
// 
// if (event_actual<1) {cout << "wrong event_actual" << endl; return 0;}
// 
// double filenum=1;
// if (inputfile_name.find("_filenum",0) != string::npos) {
//   filenum=atof(inputfile_name.substr(inputfile_name.find("_filenum")+8,inputfile_name.find("_")).c_str());
//     cout << "filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
// }
// else {
//   if (filemode=="rate"){
//     cout << "this file is rate dependent, but has no filenum, something is wrong" << endl;  
//     return 0;
//   }
//   else{
//     cout << "this file has no filenum, please check if you need filenum for addtional normalization" << endl;      
//   }
// }

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

//      TString background_inputfile_name="parametrized_lgc.root";      //h_pe is here  
//      TFile *background_file=new TFile(background_inputfile_name);
//      TH1F *h_pe=(TH1F*)background_file->Get("h_pe");

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
        
        TH2F *hhit_mum_coin=new TH2F("hhit_mum_coin","mum coin rate(khz);x(cm);y(cm)",600,-300,300,600,-300,300);
        TH2F *hhit_mup_coin=new TH2F("hhit_mup_coin","mup coin rate(khz);x(cm);y(cm)",600,-300,300,600,-300,300);
        TH2F *hhit_mu_single=new TH2F("hhit_mu_single","mu single rate(khz);x(cm);y(cm)",600,-300,300,600,-300,300);
                                
        TH1F *hcount_InvM_raw=new TH1F("hcount_InvM_raw","e,mu+,mu- raw count;InvM(GeV);count/50MeV",80,0,4);
        TH1F *hcount_InvM=new TH1F("hcount_InvM","e,mu+,mu- event count;InvM(GeV);count/50MeV",80,0,4);        
        
        
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
        
//      long int N_events = (long int)tree_header->GetEntries();
        long int N_events = (long int)tree_generated->GetEntries();     

        cout << "total number of events : " << N_events << endl;        

        //----------------------------
        //      loop trees
        //---------------------------
        for(long int i=0;i<N_events;i++){                       
//      for(long int i=1;i<N_events-1;i++){                                     
//      for(long int i=0;i<N_events/100;i++){     
//      for(long int i=N_events/2;i<N_events;i++){                      
//      for(long int i=520;i<521;i++){  //pip event
//      for(long int i=5289;i<5290;i++){          // background event                     
//                      cout<<"" << i<<endl;
                cout<<"event " << i << "\r";
//              if (i%1000==0) cout<<i<<"\r";
//              if (i%1000==0) cout<<i<<"\n";
                
                //---
                //---header tree
                //---
                
                tree_userHeader->GetEntry(i);
                
                double W_tmp=userVar008->at(0);
                double Q2_tmp=userVar009->at(0); 
                double weight_tmp=userVar010->at(0);
                if (i<10) cout << i << " " << W_tmp << " " << Q2_tmp << " " << weight_tmp << endl;
                

      double rate_convert=0;
      double effxsec=0;
  if(evgen=="grape"){
    effxsec=weight_tmp/double(N_events);
//     effxsec=weight_tmp;      //before fixing old weight_tmp=crosssection/nevent
    
//  grape generator output unit pb = 1e-36 cm2, lumi 1.2e37/cm2/s, 50 days, 0.85 eff
//       rate_convert = 1e-36*1.2e37*0.85;  
      rate_convert = 1e-36*1.2e37*0.7;          
  }
  else if(evgen=="twopeg"){
    double Ebeam=11.;
// from twopeg generator, this is max allowed W and Q2 range for Ebeam=11 and no limit on scattered e-
// Minimum W  has been changed to 1.2375
// Minimum Q2 has been changed to 5e-05
// Maximum W  has been changed to 4.63921
// W^2=M^2+2M(E-E')-Q2
    double W_min=1.2375,W_max=sqrt(0.938*0.938+2*0.938*Ebeam-Q2_tmp);
    double Q2_min=5e-5,Q2_max=0.938*0.938+2*0.938*Ebeam-W_tmp*W_tmp;
//     effxsec=(W_max-W_min)*(Q2_max-Q2_min)*weight_tmp/double(N_events);
    effxsec=(W_max-W_min)*(Q2_max-Q2_min)*weight_tmp/double(1.1e8);    // have to use actual event number
          
// twopeg generator output unit ub = 1e-30 cm2, lumi 1.2e37/cm2/s, 50 days, 0.85 eff
//       rate_convert = 1e-30*1.2e37*0.85;
      rate_convert = 1e-30*1.2e37*0.7;    
  }
      
      double count_convert = rate_convert*3600*24*50;
      double rate=effxsec*rate_convert;
      double count=effxsec*count_convert;

                //---
                //---generated tree
                //---
                tree_generated->GetEntry(i);
                int n_gen=gen_pid->size();
//                 cout<<"generated : "<<n_gen<<endl;
                
                TLorentzVector mum_gen,mup_gen;
                int pid_gen=0;
                double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;      
              //       cout << "gen_pid->size() " << gen_pid->size() << endl;        
                for (int j=0;j<gen_pid->size();j++) {
//                  cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
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
//                  cout << "p_gen " << p_gen << endl; 
                    
                    if (pid_gen==13 || pid_gen==-211) mum_gen.SetXYZM(px_gen/1e3,py_gen/1e3,pz_gen/1e3,0.105658);
                    if (pid_gen==-13 || pid_gen==211) mup_gen.SetXYZM(px_gen/1e3,py_gen/1e3,pz_gen/1e3,0.105658);
                    
                }               
//              TVector3 *vec_p_gen(px_gen,py_gen,pz_gen),*vec_v_gen(vx_gen,vy_gen,vz_gen);
                
                //              if (phi_gen<-1 || phi_gen>1) continue;
//                      cout<<"" << i<<endl;
                
//              if (vz_gen<-3550 || -3450<vz_gen) continue;

                hgen_ThetaP->Fill(theta_gen,p_gen/1e3,rate);
                hgen_ThetaPhi->Fill(theta_gen,phi_gen,rate);                            
                hgen_PhiP->Fill(phi_gen,p_gen/1e3,rate);                                                
                hgen_ThetaPhiP->Fill(theta_gen,phi_gen,p_gen/1e3,rate);                                         
                
                //---   
                //---flux tree
                //---
                tree_flux->GetEntry(i);           
//                 cout<<"flux: "<<flux_hitn->size()<<endl;
                        
                bool Is_mum_gem1=false,Is_mum_gem2=false,Is_mum_gem3=false,Is_mum_gem4=false,Is_mum_gem5=false,Is_mum_gem6=false;
                bool Is_mup_gem1=false,Is_mup_gem2=false,Is_mup_gem3=false,Is_mup_gem4=false,Is_mup_gem5=false,Is_mup_gem6=false;
                bool Is_e_gem1=false,Is_e_gem2=false,Is_e_gem3=false,Is_e_gem4=false,Is_e_gem5=false,Is_e_gem6=false;
                bool Is_mum_FAMD=false,Is_mup_FAMD=false,Is_e_FAEC=false,Is_e_LAEC=false;
                bool Is_mum=false,Is_mup=false,Is_e=false;
                int flux_index_mum=0,flux_index_mup=0;
//               cout << "flux_hitn  " << flux_hitn->size() << endl;
                double Eec=0,Eec_photonele=0,Eec_ele=0,Edepsc1=0,Edepsc2=0;
                for (Int_t j=0;j<flux_hitn->size();j++) {
//                 if(flux_id->at(j)==6140000) cout << "flux " << " !!! " << flux_hitn->at(j) << " " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;  

                double hit_vr=sqrt(pow(flux_vx->at(j),2)+pow(flux_vy->at(j),2))/1e1; //mm to cm
                double hit_vy=flux_vy->at(j)/1e1,hit_vx=flux_vx->at(j)/1e1,hit_vz=flux_vz->at(j)/1e1;           //mm to cm                
                double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))/1e1; //mm to cm
                double hit_y=flux_avg_y->at(j)/1e1,hit_x=flux_avg_x->at(j)/1e1,hit_z=flux_avg_z->at(j)/1e1;           //mm to cm                
                double hit_phi=atan2(hit_y,hit_x)*DEG;       //rad to  deg
                double hit_p=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j))/1e3;  //MeV to GeV
                
//              TVector3 *vec_hit(hit_x,hit_y,hit_z), *vec_v(hit_vx,hit_vy,hit_vz),*vec_p(flux_px->at(j),flux_py->at(j),flux_pz->at(j));
//              TVector3 *vec_path=vec_hit-vec_v_gen;


                   //find two muon and one electron
//                    if (flux_pid->at(j) == 13 && flux_mpid->at(j) == 0) {                   
//                    if (flux_pid->at(j) == 13) {                   
                   if (flux_pid->at(j) == 13 || flux_pid->at(j) == -211) {                   
                    if (flux_id->at(j)==1110000) Is_mum_gem1=true;
                    else if (flux_id->at(j)==1210000) Is_mum_gem2=true;
                    else if (flux_id->at(j)==1310000) Is_mum_gem3=true;
                    else if (flux_id->at(j)==1410000) Is_mum_gem4=true;
                    else if (flux_id->at(j)==1510000) Is_mum_gem5=true;
                    else if (flux_id->at(j)==1610000) Is_mum_gem6=true;            
                    else if (flux_id->at(j)==6140000) {Is_mum_FAMD=true;flux_index_mum=j;}
                    else{};
                   }
//                    else if (flux_pid->at(j) == -13 && flux_mpid->at(j) == 0) {                                      
//                    else if (flux_pid->at(j) == -13) {
                   if (flux_pid->at(j) == -13 || flux_pid->at(j) == 211) {                   
                    if (flux_id->at(j)==1110000) Is_mup_gem1=true;
                    else if (flux_id->at(j)==1210000) Is_mup_gem2=true;
                    else if (flux_id->at(j)==1310000) Is_mup_gem3=true;
                    else if (flux_id->at(j)==1410000) Is_mup_gem4=true;
                    else if (flux_id->at(j)==1510000) Is_mup_gem5=true;
                    else if (flux_id->at(j)==1610000) Is_mup_gem6=true;            
                    else if (flux_id->at(j)==6140000) {Is_mup_FAMD=true;flux_index_mup=j;}
                    else{};
                   }
                   else if (flux_pid->at(j) == 11 && flux_mpid->at(j) == 0) {                                      
                    if (flux_id->at(j)==1110000) Is_mup_gem1=true;
                    else if (flux_id->at(j)==1210000) Is_e_gem2=true;
                    else if (flux_id->at(j)==1310000) Is_e_gem3=true;
                    else if (flux_id->at(j)==1410000) Is_e_gem4=true;
                    else if (flux_id->at(j)==1510000) Is_e_gem5=true;
                    else if (flux_id->at(j)==1610000) Is_e_gem6=true;            
                    else if (flux_id->at(j)==3110000) Is_e_FAEC=true;
                    else if (flux_id->at(j)==3210000) Is_e_LAEC=true;                    
                    else{};
                   }
                   else{};
                  
                }       // end of flux   
                
                if (Is_mum_gem2 && Is_mum_gem3 && Is_mum_gem4 && Is_mum_gem5 && Is_mum_gem6 && Is_mum_FAMD) Is_mum=true;
                if (Is_mup_gem2 && Is_mup_gem3 && Is_mup_gem4 && Is_mup_gem5 && Is_mup_gem6 && Is_mup_FAMD) Is_mup=true;
                if (Is_e_gem2 && Is_e_gem3 && Is_e_gem4 && Is_e_gem5 && Is_e_gem6 && Is_e_FAEC) Is_e=true;
                if (Is_e_gem1 && Is_e_gem2 && Is_e_gem3 && Is_e_gem4 && Is_e_LAEC) Is_e=true;                
                
                if (Is_mum && Is_mup && Is_e){                
//                 if (Is_mum && Is_mup){
                    hhit_mum_coin->Fill(flux_avg_x->at(flux_index_mum)/1e1,flux_avg_y->at(flux_index_mum)/1e1,rate/1e3);     
                    hhit_mup_coin->Fill(flux_avg_x->at(flux_index_mup)/1e1,flux_avg_y->at(flux_index_mup)/1e1,rate/1e3);                    
                    
                    hcount_InvM->Fill((mum_gen+mup_gen).M(),count);
                    hcount_InvM_raw->Fill((mum_gen+mup_gen).M());
                }
                else if(Is_mum) hhit_mu_single->Fill(flux_avg_x->at(flux_index_mum)/1e1,flux_avg_y->at(flux_index_mum)/1e1,rate/1e3); 
                else if(Is_mup) hhit_mu_single->Fill(flux_avg_x->at(flux_index_mup)/1e1,flux_avg_y->at(flux_index_mup)/1e1,rate/1e3);
                else{}
                
              
} //end loop
        

cout <<" sensor_good " << sensor_good << endl;
cout <<" event_good " << event_good << endl;
cout <<" event_trig_good " << event_trig_good << endl;

//do outputs

outputfile->Write();    
outputfile->Flush();

TCanvas *c_muon = new TCanvas("c_muon", "c_muon",1900,900);
c_muon->Divide(3,1);
c_muon->cd(1);
hhit_mum_coin->Draw("colz");
c_muon->cd(2);
hhit_mup_coin->Draw("colz");
c_muon->cd(3);
hhit_mu_single->Draw("colz");

TCanvas *c_mass = new TCanvas("c_mass", "c_mass",1000,900);
c_mass->Divide(1,1);
c_mass->cd(1);
hcount_InvM->Draw();

TCanvas *c_mass_raw = new TCanvas("c_mass_raw", "c_mass_raw",1000,900);
c_mass_raw->Divide(1,1);
c_mass_raw->cd(1);
hcount_InvM_raw->Draw();

}
