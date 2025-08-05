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
#include <TMath.h>

using namespace std;



vector<int> *solid_lgc_hitn=0;
vector<int> *solid_lgc_sector=0,*solid_lgc_pmt=0,*solid_lgc_pixel=0,*solid_lgc_nphe=0;
vector<double> *solid_lgc_avg_t=0;

void setup_tree_solid_lgc(TTree *tree_solid_lgc)
{  
tree_solid_lgc->SetBranchAddress("hitn",&solid_lgc_hitn);
tree_solid_lgc->SetBranchAddress("sector",&solid_lgc_sector);
tree_solid_lgc->SetBranchAddress("pmt",&solid_lgc_pmt);
tree_solid_lgc->SetBranchAddress("pixel",&solid_lgc_pixel);
tree_solid_lgc->SetBranchAddress("nphe",&solid_lgc_nphe);
tree_solid_lgc->SetBranchAddress("avg_t",&solid_lgc_avg_t);

return;
}

//Simple trigger, no timing information is used.  If at least 1 sector meets the criteria for trigger, then the trigger fires.
//Must imput the lgc_tree, the event number, and the PMT and PEperPMT thresholds (default is a 2x2 trigger).

// Bool_t process_tree_solid_lgc(TTree *tree_solid_lgc,double *hit_lgc_pmt,double *hit_lgc_pixel, Int_t *trigger_lgc, Int_t &ntrigsecs, Int_t PMTthresh = 2, Int_t PEthresh = 2)
Bool_t process_tree_solid_lgc(TTree *tree_solid_lgc,bool Is_simsafe,double *hit_lgc_pmt,double *hit_lgc_quad,double *hit_lgc_pixel, Int_t *trigger_lgc, Int_t &ntrigsecs,Int_t &nphoton, Int_t PMTthresh = 2, Int_t PEthresh = 2)
{
  double factor_packing=0.8; //from PMT packing ratio
  double factor=factor_packing;
  if (Is_simsafe) factor=factor*0.5; //and addition sim safety factor  

  if(!solid_lgc_hitn->size()) return 0;
   //if using root6, uncomment line below, and comment out following line
  //std::vector<std::vector<int>> sectorhits (30, std::vector<int>(9,0));  //initialize a 30x9 vector array
  double sectorhits_pmt[30][9] = {0};  //need to intialize to zero or bad stuff
  double sectorhits_quad[30][36] = {0};  //need to intialize to zero or bad stuff  
  double sectorhits_pixel[30][576] = {0};  //need to intialize to zero or bad stuff
  
  Int_t ntrigpmts =0;
  nphoton=0;
  for(UInt_t i = 0; i < solid_lgc_hitn->size(); i++){
    if (solid_lgc_nphe->at(i)==0) nphoton += 1;  // at least 1 photon
    else nphoton += solid_lgc_nphe->at(i);  // Nphoton >= Npe
    
    if(solid_lgc_nphe->at(i)>0){
//       cout << "solid_lgc " << " !!! " << solid_lgc_hitn->size() << ", " << solid_lgc_hitn->at(i) << " " << solid_lgc_sector->at(i) << " " << solid_lgc_pmt->at(i) << " " << solid_lgc_pixel->at(i) << " " << solid_lgc_nphe->at(i) << " " << solid_lgc_avg_t->at(i) << endl;      
      sectorhits_pmt[solid_lgc_sector->at(i)-1][solid_lgc_pmt->at(i)-1] += solid_lgc_nphe->at(i)*factor;
      sectorhits_pixel[solid_lgc_sector->at(i)-1][(solid_lgc_pmt->at(i)-1)*64+solid_lgc_pixel->at(i)-1] += solid_lgc_nphe->at(i)*factor;
      
      int quadid=1;
      if ((1<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=4) || (9<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=12) || (17<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=20) || (25<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=28)) quadid=1;
      else if ((5<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=8) || (13<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=16) || (21<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=24) || (19<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=32)) quadid=2;
      else if ((33<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=36) || (41<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=44) || (49<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=52) || (57<=solid_lgc_pixel->at(i) && solid_lgc_pixel->at(i)<=60)) quadid=3;
      else quadid=4;
      
      sectorhits_quad[solid_lgc_sector->at(i)-1][(solid_lgc_pmt->at(i)-1)*4+quadid-1] += solid_lgc_nphe->at(i)*factor;      

    }
  }
   
  for(UInt_t i = 0; i < 30; i++){
    ntrigpmts = 0;
    for(Int_t j = 0; j < 9; j++){      
      hit_lgc_pmt[i*9+j]=sectorhits_pmt[i][j];            

      if(sectorhits_pmt[i][j] >= PEthresh) ntrigpmts++;
    }   
    if(ntrigpmts >= PMTthresh) {
      ntrigsecs++;
      trigger_lgc[i]=1;
    }
    
    for(Int_t j = 0; j < 36; j++){      
      hit_lgc_quad[i*36+j]=sectorhits_quad[i][j];            
    }        
    for(Int_t j = 0; j < 576; j++){      
      hit_lgc_pixel[i*576+j]=sectorhits_pixel[i][j];            
    }    
  }
  if(ntrigsecs){
    return 1;
  }else{
    return 0;
  }
}


