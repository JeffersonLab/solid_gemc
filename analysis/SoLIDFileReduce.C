#include <vector>
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
#include <TColor.h>

using namespace std;

void SoLIDFileReduce(string setup, string input_filename, string format="1")
{
	if(setup=="PVDIS" || setup=="SIDIS" || setup=="HGC" || setup=="flux"){cout << "setup " << setup << endl;}
	else {cout << "can't recognize setup " << setup << endl; return;}
  
  	if(format=="1" || format=="2"){cout << "format " << format << endl;}
	else {cout << "can't recognize format " << format << endl; return;}
	
	TFile *f1 = new TFile(input_filename.c_str());
	TFile *f2 = new TFile(Form("%s_skim",input_filename.c_str()), "RECREATE");
	
	TTree *header;
	TTree *generated;
	TTree *userHeader;		
	TTree *ancestors;	
	TTree *chargeTime;	
	TTree *counter;	
	TTree *mirror;	
	TTree *rf;	
	TTree *allraws;  
	TTree *flux;	
	TTree *solid_gem;
	TTree *solid_ec;
	TTree *solid_ec_ps;	
	TTree *solid_lgc;
	TTree *solid_mrpc;
	TTree *solid_spd;
	TTree *solid_hgc;
	
	header = (TTree*)f1->Get("header");
	generated = (TTree*)f1->Get("generated");
	if(format=="2"){
	userHeader = (TTree*)f1->Get("userHeader");	
	ancestors = (TTree*)f1->Get("ancestors");
	chargeTime = (TTree*)f1->Get("chargeTime");	
	counter = (TTree*)f1->Get("counter");
	mirror = (TTree*)f1->Get("mirror");	
	rf = (TTree*)f1->Get("rf");
	}	
	if(setup=="flux"){	
	allraws = (TTree*)f1->Get("allraws");
	flux = (TTree*)f1->Get("flux");
	}
	if(setup=="PVDIS"){	
	allraws = (TTree*)f1->Get("allraws");
	flux = (TTree*)f1->Get("flux");
	solid_gem = (TTree*)f1->Get("solid_gem");
	solid_ec = (TTree*)f1->Get("solid_ec");
	solid_ec_ps = (TTree*)f1->Get("solid_ec_ps");	
	solid_lgc = (TTree*)f1->Get("solid_lgc");
	}	
	if(setup=="SIDIS"){
	allraws = (TTree*)f1->Get("allraws");
	flux = (TTree*)f1->Get("flux");
	solid_gem = (TTree*)f1->Get("solid_gem");
	solid_ec = (TTree*)f1->Get("solid_ec");
	solid_ec_ps = (TTree*)f1->Get("solid_ec_ps");		
	solid_lgc = (TTree*)f1->Get("solid_lgc");	  
	solid_mrpc = (TTree*)f1->Get("solid_mrpc");
	solid_spd = (TTree*)f1->Get("solid_spd");
	solid_hgc = (TTree*)f1->Get("solid_hgc");
	}
	if(setup=="HGC"){
	allraws = (TTree*)f1->Get("allraws");
	flux = (TTree*)f1->Get("flux");
	solid_hgc = (TTree*)f1->Get("solid_hgc");
	}
		
	
// 	cout <<   flux->GetBranch("hitn")->GetLeaf("hitn")->GetTypeName() << endl;
// 	cout <<   flux->FindLeaf("hitn")->GetValue(1) << endl;
	
	TTree *output_header;
	TTree *output_generated;
	TTree *output_userHeader;		
	TTree *output_ancestors;	
	TTree *output_chargeTime;	
	TTree *output_counter;	
	TTree *output_mirror;	
	TTree *output_rf;	
	TTree *output_allraws;  
	TTree *output_flux;	
	TTree *output_solid_gem;
	TTree *output_solid_ec;
	TTree *output_solid_ec_ps;	
	TTree *output_solid_lgc;
	TTree *output_solid_mrpc;
	TTree *output_solid_spd;
	TTree *output_solid_hgc;
	
	output_header = header->CloneTree(0);
	output_generated = generated->CloneTree(0);
	if(format=="2"){	
	output_userHeader = userHeader->CloneTree(0);
	output_ancestors = ancestors->CloneTree(0);
	output_chargeTime = chargeTime->CloneTree(0);
	output_counter = counter->CloneTree(0);
	output_mirror = mirror->CloneTree(0);
	output_rf = rf->CloneTree(0);	
	}
	if(setup=="flux"){	
	output_allraws = allraws->CloneTree(0);
	output_flux = flux->CloneTree(0);	
	}
	if(setup=="PVDIS"){	
	output_allraws = allraws->CloneTree(0);
	output_flux = flux->CloneTree(0);	
	output_solid_gem = solid_gem->CloneTree(0);
	output_solid_ec = solid_ec->CloneTree(0);
	output_solid_ec_ps = solid_ec_ps->CloneTree(0);		
	output_solid_lgc = solid_lgc->CloneTree(0);
	}	
	if(setup=="SIDIS"){
	output_allraws = allraws->CloneTree(0);
	output_flux = flux->CloneTree(0);	
	output_solid_gem = solid_gem->CloneTree(0);
	output_solid_ec = solid_ec->CloneTree(0);
	output_solid_ec_ps = solid_ec_ps->CloneTree(0);	
	output_solid_lgc = solid_lgc->CloneTree(0);	  
	output_solid_mrpc = solid_mrpc->CloneTree(0);
	output_solid_spd = solid_spd->CloneTree(0);
	output_solid_hgc = solid_hgc->CloneTree(0);
	}
	if(setup=="HGC"){
	output_allraws = allraws->CloneTree(0);
	output_flux = flux->CloneTree(0);	
	output_solid_hgc = solid_hgc->CloneTree(0);
	}	

//!!! have to use double instead of int, otherwise output tree will have weird vector out_of_range error	
vector<double> *flux_hitn=0;
vector<double> *allraws_hitn=0,*solid_gem_hitn=0,*solid_ec_hitn=0,*solid_ec_ps_hitn=0,*solid_lgc_hitn=0,*solid_mrpc_hitn=0,*solid_spd_hitn=0,*solid_hgc_hitn=0;

	if(setup=="flux"){
	allraws->SetBranchAddress("hitn", &allraws_hitn);
	flux->SetBranchAddress("hitn", &flux_hitn);
	}
	if(setup=="PVDIS"){	
	allraws->SetBranchAddress("hitn", &allraws_hitn);
	flux->SetBranchAddress("hitn", &flux_hitn);
	solid_gem->SetBranchAddress("hitn", &solid_gem_hitn);
	solid_ec->SetBranchAddress("hitn", &solid_ec_hitn);
	solid_ec_ps->SetBranchAddress("hitn", &solid_ec_ps_hitn);	
	solid_lgc->SetBranchAddress("hitn", &solid_lgc_hitn);
	}	
	if(setup=="SIDIS"){
	allraws->SetBranchAddress("hitn", &allraws_hitn);
	flux->SetBranchAddress("hitn", &flux_hitn);
	solid_gem->SetBranchAddress("hitn", &solid_gem_hitn);
	solid_ec->SetBranchAddress("hitn", &solid_ec_hitn);
	solid_ec_ps->SetBranchAddress("hitn", &solid_ec_ps_hitn);		
	solid_lgc->SetBranchAddress("hitn", &solid_lgc_hitn);
	solid_mrpc->SetBranchAddress("hitn", &solid_mrpc_hitn);
	solid_spd->SetBranchAddress("hitn", &solid_spd_hitn);
	solid_hgc->SetBranchAddress("hitn", &solid_hgc_hitn);
	}
	if(setup=="HGC"){
	allraws->SetBranchAddress("hitn", &allraws_hitn);
	flux->SetBranchAddress("hitn", &flux_hitn);
	solid_hgc->SetBranchAddress("hitn", &solid_hgc_hitn);
	}

	long long int fNentries = generated->GetEntries();
	
	cout << "total number of events " << fNentries << endl;
	
	int evncounter=0;
	for (int i=0; i<fNentries; i++){
// 		if (i % 1000000 == 0) std::cout<<i<<endl;
			header->GetEntry(i);
			generated->GetEntry(i);
			if(format=="2"){			
			userHeader->GetEntry(i);
			ancestors->GetEntry(i);
			chargeTime->GetEntry(i);
			counter->GetEntry(i);
			mirror->GetEntry(i);
			rf->GetEntry(i);				
			}
			if(setup=="flux"){				
			allraws->GetEntry(i);
			flux->GetEntry(i);
			}
			if(setup=="PVDIS"){
			allraws->GetEntry(i);
			flux->GetEntry(i);			  
			solid_gem->GetEntry(i);
			solid_ec->GetEntry(i);
			solid_ec_ps->GetEntry(i);			
			solid_lgc->GetEntry(i);
			}	
			if(setup=="SIDIS"){
			allraws->GetEntry(i);
			flux->GetEntry(i);			  
			solid_gem->GetEntry(i);
			solid_ec->GetEntry(i);
			solid_ec_ps->GetEntry(i);						
			solid_lgc->GetEntry(i);	  
			solid_mrpc->GetEntry(i);
			solid_spd->GetEntry(i);
			solid_hgc->GetEntry(i);
			}
			if(setup=="HGC"){
			allraws->GetEntry(i);
			flux->GetEntry(i);			  
			solid_hgc->GetEntry(i);
			}		
			
			if(setup=="flux"){				
			  if (flux_hitn->size() > 0){
			      evncounter++;			      
			      output_header->Fill();
			      output_generated->Fill();
			      if(format=="2"){			      
			      output_userHeader->Fill();
			      output_ancestors->Fill();
			      output_chargeTime->Fill();
			      output_counter->Fill();
			      output_mirror->Fill();
			      output_rf->Fill();			      
			      }
			      output_allraws->Fill();
			      output_flux->Fill();
			  }
			}			
			if(setup=="PVDIS"){	
			  if (flux_hitn->size() > 0 || 
			      solid_gem_hitn->size() > 0 ||			      
			      solid_ec_hitn->size() > 0 ||
			      solid_ec_ps_hitn->size() > 0 ||			      
			      solid_lgc_hitn->size() > 0){
			    
			      evncounter++;			      
			      output_header->Fill();
			      output_generated->Fill();
			      if(format=="2"){			  
			      output_userHeader->Fill();
			      output_ancestors->Fill();
			      output_chargeTime->Fill();
			      output_counter->Fill();
			      output_mirror->Fill();
			      output_rf->Fill();
			      }
			      output_allraws->Fill();
			      output_flux->Fill();		
			      output_solid_gem->Fill();
			      output_solid_ec->Fill();
			      output_solid_ec_ps->Fill();			      
			      output_solid_lgc->Fill();
			    }	
			}	
			if(setup=="SIDIS"){
			  if (flux_hitn->size() > 0 ||			      
			      solid_gem_hitn->size() > 0 ||			    
			      solid_ec_hitn->size() > 0 ||
			      solid_ec_ps_hitn->size() > 0 ||			      			      
			      solid_lgc_hitn->size() > 0 ||		  
			      solid_mrpc_hitn->size() > 0 ||
		              solid_spd_hitn->size() > 0 ||			      
		              solid_hgc_hitn->size() > 0){
			    
			      evncounter++;			      
			      output_header->Fill();
			      output_generated->Fill();
			      if(format=="2"){			  
			      output_userHeader->Fill();
			      output_ancestors->Fill();
			      output_chargeTime->Fill();
			      output_counter->Fill();
			      output_mirror->Fill();
			      output_rf->Fill();
			      }
			      output_allraws->Fill();
			      output_flux->Fill();		
			      output_solid_gem->Fill();
			      output_solid_ec->Fill();
			      output_solid_ec_ps->Fill();			      			      
			      output_solid_lgc->Fill();		
			      output_solid_mrpc->Fill();
			      output_solid_spd->Fill();
			      output_solid_hgc->Fill();
			    }	
			}
			if(setup=="HGC"){
			  if (flux_hitn->size() > 0 ||				      
		              solid_hgc_hitn->size() > 0){
			    
			      evncounter++;			      
			      output_header->Fill();
			      output_generated->Fill();
			      if(format=="2"){			  
			      output_userHeader->Fill();
			      output_ancestors->Fill();
			      output_chargeTime->Fill();
			      output_counter->Fill();
			      output_mirror->Fill();
			      output_rf->Fill();
			      }
			      output_allraws->Fill();
			      output_flux->Fill();	
			      output_solid_hgc->Fill();
			    }
			}
			
	}
	cout << "skim event counter " <<  evncounter << endl;
	
	f2->cd();
	
	output_header->Write();
	output_generated->Write();
	if(format=="2"){	
	output_userHeader->Write();
	output_ancestors->Write();
	output_chargeTime->Write();
	output_counter->Write();
	output_mirror->Write();
	output_rf->Write();	
	}
	if(setup=="flux"){				
	output_allraws->Write();
	output_flux->Write();
	}
	if(setup=="PVDIS"){	
	output_allraws->Write();
	output_flux->Write();	  
	output_solid_gem->Write();
	output_solid_ec->Write();
	output_solid_ec_ps->Write();	
	output_solid_lgc->Write();
	}	
	if(setup=="SIDIS"){
	output_allraws->Write();
	output_flux->Write();	  
	output_solid_gem->Write();
	output_solid_ec->Write();
	output_solid_ec_ps->Write();		
	output_solid_lgc->Write();  
	output_solid_mrpc->Write();
	output_solid_spd->Write();
	output_solid_hgc->Write();
	}
	if(setup=="HGC"){
	output_allraws->Write();
	output_flux->Write();	  
	output_solid_hgc->Write();
	}	
	
	f2->Close();
	f1->Close();
	
};
