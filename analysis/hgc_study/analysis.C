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

using namespace std;

void analysis(string input_filename,string cherenkov,bool Is_decaycut=false,bool Is_simsafe=false,int pid_opticalphoton=-22)
{
gROOT->Reset();
// gStyle->SetPalette(57);
// gStyle->SetOptStat(0);
// gStyle->SetOptStat(1);
// gStyle->SetOptStat("ioue");
gStyle->SetOptStat(11111111);
gStyle->SetOptFit(111111);
gStyle->SetPadRightMargin(0.15); 

TRandom3 rand;
rand.SetSeed(0);		  

//safety factor for pion, manual reduce 2
// double factor=0.5;

const double DEG=180./3.1415926;

int vplaneid;
int ch_pmt;
if (cherenkov=="lgc"){
  vplaneid=2110000;
  ch_pmt=270;
}
else if (cherenkov=="hgc") {
  vplaneid=2210000;  
  ch_pmt=480;
}
else {cout << "bad Cherenkov name" << endl; return; }


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
sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());
// sprintf(the_filename, "%s",input_filename.substr(input_filename.rfind("/")+1,input_filename.rfind(".")).c_str());
// sprintf(the_filename,Form("%s/%s",the_dir.c_str(),the_filename));

cout << " the_filename " << the_filename << endl;

int Ngood=0,Nrand=0;

char output_filename[200];
sprintf(output_filename, "%s_output.root",the_filename);
TFile *outputfile=new TFile(output_filename, "recreate");

TH1F *hnpe_no0=new TH1F("hnpe_no0",";Npe;event count",200,0,200);

TH1F *hnpe=new TH1F("hnpe",";Npe;event count",200,0,200);

TH1F *hnpe_QE=new TH1F("hnpe_QE",";Npe/Nphoton;event count",100,0,1);

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

TTree *tree_solid_cc;
if (cherenkov=="lgc"){
tree_solid_cc = (TTree*) file->Get("solid_lgc");
setup_tree_solid_lgc(tree_solid_cc);
}
else if (cherenkov=="hgc") {
tree_solid_cc = (TTree*) file->Get("solid_hgc");
setup_tree_solid_hgc(tree_solid_cc);
}
else {cout << "bad Cherenkov name" << endl; }

int nevent = (int)tree_generated->GetEntries();
int nselected = 0;
cout << "nevent " << nevent << endl;

int sensor_good=0;
int event_good=0;
int decay_no=0,decay_first=0,decay_lost=0;

for (Int_t i=0;i<nevent;i++) { 
// for (Int_t i=0;i<nevent/100;i++) {
// for (Int_t i=0;i<10;i++) { 
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
    double hit_px=0,hit_py=0,hit_pz=0,hit_x=0,hit_y=0,hit_z=0,hit_r=0,hit_phi=0,hit_phip=0;
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
	  Is_decay=false;
	  
	  hit_px=flux_px->at(j);
	  hit_py=flux_py->at(j);
	  hit_pz=flux_pz->at(j);	  
	  hit_x=flux_avg_x->at(j);
	  hit_y=flux_avg_y->at(j);	  
	  hit_z=flux_avg_z->at(j);
	  hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2));
	  hit_phi=atan2(hit_y,hit_x)*DEG;       //rad to deg		  	  
	  
	  hit_phip=atan2(hit_py,hit_px)*DEG;       //rad to deg			  
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
		
		if (cherenkov=="lgc") process_tree_solid_lgc(tree_solid_cc,Is_simsafe,hit_pmt,hit_quad,hit_pixel,trigger,ntrigsecs,nphoton,PMTthresh,PEthresh);		
		else if (cherenkov=="hgc") process_tree_solid_hgc(tree_solid_cc,Is_simsafe,hit_pmt,hit_quad,hit_pixel,trigger,ntrigsecs,PMTthresh,PEthresh,pid_opticalphoton);
		else {cout << "bad Cherenkov name" << endl; return;}
		
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
			
		if (npe_total>0) hnpe_no0->Fill(npe_total);
				
		hnpe->Fill(npe_total);	
		
		hnpe_QE->Fill(double(npe_total)/double(nphoton));
		
		//tdc for pixel
		for(int index=0;index<ch_pixel;index++) if(hit_pixel[index]>0) hit_pixel[index]=1;  

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
// 		    event_good++;
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
cout <<" event_good " << event_good << endl;
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
