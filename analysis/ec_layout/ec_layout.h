#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"

//routine to load ecal block id,coordinates, and sector information
// void LoadEC_map(TString map_file);
void LoadEC_map();
//routine to return ecal block x,y
TVector2 GetECALBlock_coord(Int_t block_id);
Int_t GetECALBlock_sector(Int_t block_id);
Int_t GetECALBlock_id(double hit_x,double hit_y);
//return an pointer to energy deposit for 6 neighboring ecal blocks
Int_t GetECALCluser(Int_t Nneighbor,Int_t block_id,Int_t *cluster_edep_blockid);
//Int_t GetECALCluser(Double_t hit_x, Double_t hit_y, Int_t *cluster_edep_blockid);
//save x,y,id and sectors and status for which moudles are active only for SIDIS case
Int_t GetECALCluserAll(Int_t Nneighbor,Int_t block_id,Int_t *cluster_edep_blockid);
Int_t sector[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t id[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
Int_t num_module_in_row[54]={0};
Double_t y_bak[54]={100000};
Double_t x[54][50]={{100000}};           
Double_t y[54][50]={{100000}}; 
Int_t status[54][50]={{100000}};    //the structure is the same as x[54][50], y[54][50]   54 is the number of y rows
std::pair<Int_t,Int_t> block_map[2000];

// void LoadEC_map(TString map_file){
// 	if (map_file == "") map_file="map_FAEC_ANL_20130628.txt";
void LoadEC_map(string path){
// 	cout << " this " << gSystem->Getenv("solid_gemc") << endl;
// 	TString map_file= string(gSystem->Getenv("solid_gemc"))+"/analysis/ec_layout/map_FAEC_ANL_20130628.txt";
        TString map_file= path+"/map_FAEC_ANL_20130628.txt";
	ifstream INPUT_file;
	INPUT_file.open(map_file.Data());
	if(!INPUT_file){ 
		printf("ERROR!!! Can't open %s \n",map_file.Data());	
		exit(1);
	}
	else printf("open ec_layout map %s \n",map_file.Data());



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
	int counter_id_new=0;
	
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

Int_t GetECALBlock_id(double hit_x=0,double hit_y=0){
	int hit_idx;                   //once we get hit id, we can easily locate the hitted module coordinate
	int hit_idy;
	
	int first_y=0;
	int last_y=54-1;
	while( last_y!=first_y+1 ){
		int mid_y=(int)(first_y+last_y)/2;
		if(hit_y>y_bak[mid_y]){
			first_y=mid_y;
		}else if(hit_y<y_bak[mid_y]){
			last_y=mid_y;
		}else{
			first_y=mid_y;
			last_y=first_y+1;
		}
	}
	if( fabs(hit_y-y_bak[first_y])<=fabs(hit_y-y_bak[last_y]) ){
		hit_idy=first_y;                                 //found y id
	}else{
		hit_idy=last_y;
	}

	//then go and find x id
	int first_x=0;
	int last_x=num_module_in_row[hit_idy]-1;   //last depends on which row the hit it is on
	while( last_x!=first_x+1 ){
		int mid_x=(int)(first_x+last_x)/2;
		if(hit_x>x[hit_idy][mid_x]){
			first_x=mid_x;
		}else if(hit_x<x[hit_idy][mid_x]){
			last_x=mid_x;
		}else{
			first_x=mid_x;
			last_x=first_x+1;
		}
	}
	if( fabs(hit_x-x[hit_idy][first_x])<=fabs(hit_x-x[hit_idy][last_x]) ){
		hit_idx=first_x;            //found x id
	}else{
		hit_idx=last_x;
	}
	
// 	cout << hit_x << " " << hit_y << " " << hit_idx+1 << " " << hit_idy+1 << " " << id[hit_idy][hit_idx] << endl;
	return id[hit_idy][hit_idx];
};

Int_t GetECALCluser(Int_t Nneighbor,Int_t block_id,Int_t *cluster_edep_blockid){
// 	cluster_edep_blockid has [0] as the central block and varying number of other actual neighbor blocks
//	it return number of actual blocks including the central one
  
	Int_t hit_idx = block_map[block_id].second;
	Int_t hit_idy = block_map[block_id].first;
	Int_t hit_around_idx[100]={100000};   //a variable label indicates how many surrounded modules are around the hitted module
	Int_t hit_around_idy[100]={100000};
	//__________________________________find the surrounded other 18 modules________________________
	Int_t tmp_idx[21]={hit_idx-10,hit_idx-9,hit_idx-8,hit_idx-7, hit_idx-6, hit_idx-5, hit_idx-4, hit_idx-3, hit_idx-2, hit_idx-1, hit_idx, hit_idx+1, hit_idx+2, hit_idx+3, hit_idx+4, hit_idx+5, hit_idx+6, hit_idx+7, hit_idx+8, hit_idx+9, hit_idx+10};
	Int_t tmp_idy[5]={hit_idy-2,hit_idy-1, hit_idy, hit_idy+1,hit_idy+2};
	Int_t label=0;
	
	//6.423cm is the side length of layout hexagon,not module hexagon which is smaller. 
// 	so distance between center and 6 neighbor layout hexagon is 11.176cm, which is smaller than 15cm
// 	so max distance between center and 18 neighbor layout hexagon is 22.343cm, which is smaller than 25cm	
	int dmax=0;	  	
	if(Nneighbor==7){
	  dmax=15;	  
	}else if (Nneighbor==19){
	  dmax=25;	  	  
	}else {cout << "wrong Nneighbor value" << endl;}
	  
	for(Int_t i=0;i<21;i++){// x scan
		for(Int_t j=0;j<5;j++){  //y scan
			if(tmp_idy[j]>=0 && tmp_idy[j]<54 && tmp_idx[i]>=0 && tmp_idx[i]<num_module_in_row[tmp_idy[j]] && (tmp_idx[i]!=hit_idx || tmp_idy[j]!=hit_idy) ){ // in range
				if(sqrt( pow( (x[tmp_idy[j]][tmp_idx[i]]-x[hit_idy][hit_idx]),2 )+ pow((y[tmp_idy[j]][tmp_idx[i]]-y[hit_idy][hit_idx]),2 ))<dmax){
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
	return label+1;
};


Int_t GetECALCluserAll(Int_t Nneighbor,Int_t block_id,Int_t *cluster_edep_blockid){
// 	cluster_edep_blockid have blocks in the order of from lower to higher in y and right to left in x. for none-existing block, it has 0 value
//	it return number of actual blocks including the central one
  
	Int_t hit_idx = block_map[block_id].second;
	Int_t hit_idy = block_map[block_id].first;

	//6.423cm is the side length of layout hexagon,not module hexagon which is smaller. 
// 	so distance between center and 6 neighbor layout hexagon is 11.176cm, which is smaller than 15cm
// 	so max distance between center and 18 neighbor layout hexagon is 22.343cm, which is smaller than 25cm		
	if(Nneighbor==7){
	for(Int_t i=hit_idy-1;i<=hit_idy+1;i++){
	  for(Int_t j=0;j<num_module_in_row[i];j++){
	    if(i==hit_idy-1 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[0]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[1]=id[i][j];
	    }
	    else if(i==hit_idy && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[2]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-0)<1) cluster_edep_blockid[3]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[4]=id[i][j];	      
	    }
	    else if(i==hit_idy+1 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[5]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[6]=id[i][j];
	    }    
	    
	  }
	}
	}
	else if(Nneighbor==19){
	for(Int_t i=hit_idy-2;i<=hit_idy+2;i++){
	  for(Int_t j=0;j<num_module_in_row[i];j++){
	    if(i==hit_idy-2 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[0]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1) cluster_edep_blockid[1]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[2]=id[i][j];
	    }
	    else if(i==hit_idy-1 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[3]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[4]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[5]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[6]=id[i][j];
	    }
	    else if(i==hit_idy && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[7]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[8]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-0)<1) cluster_edep_blockid[9]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[10]=id[i][j];	      
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[11]=id[i][j];	      
	    }	    
	    else if(i==hit_idy+1 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[12]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[13]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-11.176)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[14]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[15]=id[i][j];
	    }
	    else if(i==hit_idy+2 && (0<=i && i<=54)){
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]<x[hit_idy][hit_idx]) cluster_edep_blockid[16]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-19.269)<1) cluster_edep_blockid[17]=id[i][j];
	      if(fabs(sqrt(pow((x[i][j]-x[hit_idy][hit_idx]),2)+pow((y[i][j]-y[hit_idy][hit_idx]),2))-22.343)<1 && x[i][j]>x[hit_idy][hit_idx]) cluster_edep_blockid[18]=id[i][j];
	    }
	  }
	}
	}
		
	int label=0;
	for(int l=0;l<Nneighbor;l++){
	  if(cluster_edep_blockid[l] !=0) label++;
	}	
	
	return label;
};
