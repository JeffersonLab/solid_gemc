{
  gROOT->Reset();
// gStyle->SetPalette(1);
gStyle->SetOptStat(0);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);  
  
  gStyle->SetPadColor(0);

  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.04,"xyz");   
  gStyle->SetTitleSize(0.07,"t");    
  gStyle->SetPaintTextFormat("4.1f");   

char hstname[100],title[500],outname[100]; 

const int n=3;
char *filename[n]={ 
"solid_JPsi_DDVCS_LH2_moved_muon_even_mum_theta12_phi0_flux_1e6_output_muon.root",
"solid_JPsi_DDVCS_LH2_moved_muon_even_pim_theta12_phi0_flux_1e6_output_pionprimary.root",
"solid_JPsi_DDVCS_LH2_moved_muon_even_pim_theta12_phi0_flux_1e6_output_pionsecondary.root",
};
char *label[n]={"muon","primary pion","secondary pion"};
sprintf(hstname,"hhit_Edep_FAMD");
// sprintf(title,"count/50MeV;l^{+}l^{-} InvM (GeV);");
sprintf(title,";Edep (GeV);");
sprintf(outname,"compare");
int color[n]={1,2,4};
int style[n]={1,1,1};

// const int n=2;
// char *filename[n]={ 
// // "/volatile/halla/solid/sim/solid_gemc/JPsi_DDVCS_JLAB_VERSION_2.5/pass1/solid_JPsi_DDVCS_LH2_moved_full_grape_flux_2.5e6_output.root_2muFAMD",
// // "/volatile/halla/solid/sim/solid_gemc/JPsi_DDVCS_JLAB_VERSION_2.5/pass1/solid_JPsi_DDVCS_LH2_moved_full_grape_flux_2.5e6_output.root_2muFAMDGEM",
// 
// "/volatile/halla/solid/sim/solid_gemc/JPsi_DDVCS_JLAB_VERSION_2.5/pass1/solid_JPsi_DDVCS_LH2_moved_full_grape_flux_2.5e6_output.root",
// 
// // "/volatile/halla/solid/sim/solid_gemc/JPsi_DDVCS_JLAB_VERSION_2.5/pass1/solid_JPsi_DDVCS_LH2_moved_full_twopeg_flux_1e6_output.root",
// "/volatile/halla/solid/sim/solid_gemc/JPsi_DDVCS_JLAB_VERSION_2.5/pass1/solid_JPsi_DDVCS_LH2_moved_full_twopeg_flux_1.1e8_FAMD_output.root",
// 
// };
// char *label[n]={"grape di-muon","twopeg di-muon"};
// sprintf(hstname,"hcoin");
// // sprintf(title,"count/50MeV;l^{+}l^{-} InvM (GeV);");
// sprintf(title,";l^{+}l^{-} InvM (GeV);count/50MeV");
// sprintf(outname,"compare_solid_ddvcs_11GeV_22GeV_InvM_FALA_nosmear_mumup");
// int color[n]={1,2};
// int style[n]={1,1};

TFile *file[n];
TH1F *h[n];
TH1F *hadd_ISRon,*hadd_ISRoff;
TCanvas *c = new TCanvas("c","c",1200,800);
TLegend* leg = new TLegend(0.5, 0.5, 0.9, 0.9);
// TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
for (Int_t i=0;i<n;i++) {
//   if (i>2) continue;
//     if (i==3) sprintf(hstname,"lepIM1_2");
  cout << i << endl;
  
//   cout << Form("%s/grp_out.root",name[i]) << endl;
//   file[i]=new TFile(Form("%s/acceptance_forward/grp_out.root",name[i]));  
//   file[i]=new TFile(Form("%s/acceptance_forwardandlarge/grp_out.root",name[i]));    
//   file[i]=new TFile(Form("%s/grp_out.root",name[i]));
  file[i]=new TFile(Form("%s",filename[i]));        

  h[i]=(TH1F*) file[i]->Get(hstname);
  
      
   h[i]->SetLineColor(color[i]);
   h[i]->SetLineStyle(style[i]);  
//  h[i]->SetMaximum(5e5); 
//  h[i]->SetMaximum(1.5e-5); 
//  h[i]->SetMaximum(0.05);  
//  h[i]->SetMaximum(0.12);   
//  h[i]->SetMaximum(0.4);   
 h[i]->SetTitle(title);
 

  if (i==0) h[i]->Scale(0.02);
  
 if (i==0) h[i]->Draw("H");
 else h[i]->Draw("H same");
//  if (i==0) h[i]->Draw();
//  else h[i]->Draw("same"); 
  cout << "count " << h[i]->Integral() << endl;
 
//  if (i==0) hadd_ISRon=(TH1F*)h[i]->Clone();
//  if (i==1) hadd_ISRon->Add(h[i]); 
//  if (i==4) hadd_ISRoff=(TH1F*)h[i]->Clone();
//  if (i==5) hadd_ISRoff->Add(h[i]); 
 
// label = new TText(0.5,8e-6+i*2e-6,name[i]);
// label = new TText(0.5,8e-6+i*2e-6,name[i]); 
// label->SetTextColor(i+1);
// label->SetTextSize(0.05);
// label->Draw();
  leg->AddEntry(h[i],label[i],"l");  
//   leg->AddEntry(h[i],Form("%s, integral %1.2f",label[i],h[i]->Integral()),"l");   
}

//   hadd_ISRon->SetLineColor(6);
//   hadd_ISRon->SetLineStyle(1);  
//   hadd_ISRon->Draw("H same");
//   hadd_ISRoff->SetLineColor(6);
//   hadd_ISRoff->SetLineStyle(2);  
//   hadd_ISRoff->Draw("H same");
//   leg->AddEntry(hadd_ISRon,Form("ISRon BH+Compton, integral %1.2f",hadd_ISRon->Integral()),"l");   
//   leg->AddEntry(hadd_ISRoff,Form("ISRoff BH+Compton, integral %1.2f",hadd_ISRoff->Integral()),"l");   
  
  leg->Draw();
  
// c->SaveAs(Form("%s.png",outname));
// c->SaveAs(Form("%s.pdf",outname));

}