{
gROOT->Reset();
// gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// gStyle->SetOptStat(1111111);
// gStyle->SetPadLeftMargin(0.1);
// gStyle->SetPadRightMargin(0.32);
//   gStyle->SetLabelSize(0.1,"xyz"); // size of axis values
  
// char input_dir[200]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/old/withoutEC";
// char input_dir[200]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass1";

// char input_dir[200]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2";
// const int m=2;
// char* input_filename[m]={
// "cc_pro_He3_CO2LA_even_e_central_1e3_output.root",
// "cc_pro_He3_CO2LA_even_e_broad_1e5_output.root",
// };
// char *label[m]={"even e central","even e broad"};
// int MarkerStyle[m]={24,26};
// int color[m]={1,1};
// int style[m]={1,1};
// char *hst[m]={"pe_hgc_3","pe_hgc_3"};
// // char *hst[m]={"npe_hgc_1","npe_hgc_1"};
// double ymin=1,ymax=4e4;
// double xmin=1,xmax=100;
// char *title="TCD with CO2 at Large Angle";

// char input_dir[300]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2";
// const int m=6;
// char* input_filename[m]={
// "cc_pro_He3_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",  "cc_pro_Beup_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_Bedown_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_He3_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_Beup_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_Bedown_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// };
// char *title="photon time with CO2 at Large Angle";
// char *label[m]={
// "eDIS (He3)","eDIS (Be up)","eDIS (Be down)",
// "pi0Wiser (He3)","pi0Wiser (Be up)","pi0Wiser (Be down)",
// };
// int MarkerStyle[m]={
// 24,24,24,25,25,25,
// };
// int color[m]={
// 1,2,3,1,2,3,
// };
// int style[m]={
// 1,1,1,2,2,2,
// };
// char *hst[m]={
// "time_photon","time_photon","time_photon","time_photon","time_photon","time_photon",
// };
// double ymin=1e-6,ymax=1e6;

// // pass2
// // char input_dir[300]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2/output_SC0.7_EC0.3";
// // char input_dir[300]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2/output_SC0.7_EC0.1";
// // char input_dir[300]="/work/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2/output_trigSC0.7EC0.1/";
// // char input_dir[300]="/work/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2/output_trigno/";
// char input_dir[300]="/work/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass2/";
// 
// const int m=28;
// 
// char* input_filename[m]={
// "cc_pro_Beup_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_Alup_CO2LA_dirty_weighted_eDIS_filenum98_0.98e9_skim_output.root",
// "cc_pro_Aldown_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_Bedown_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_He3_N2up_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_He3_winup_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_He3_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// 
// "cc_pro_He3_windown_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_He3_N2down_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// "cc_pro_C1_CO2LA_dirty_weighted_eDIS_filenum99_0.99e9_skim_output.root",
// "cc_pro_C2_CO2LA_dirty_weighted_eDIS_filenum97_0.97e9_skim_output.root",
// "cc_pro_C3_CO2LA_dirty_weighted_eDIS_filenum94_0.94e9_skim_output.root",
// "cc_pro_C_N2_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// 
// "cc_pro_Beup_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_Alup_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_Aldown_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_Bedown_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_He3_N2up_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_He3_winup_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_He3_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_He3_windown_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_He3_N2down_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_C1_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// "cc_pro_C2_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // "cc_pro_C3_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // "cc_pro_C_N2_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// 
// // "cc_pro_He3_CO2LA_BeamOnTarget_0.966e11_skim_output.root",
// // "cc_pro_C_CO2LA_BeamOnTarget_0.99e10_skim_output.root",
// // "cc_pro_LH2_CO2LA_BeamOnTarget_1e10_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_0.878e10_skim_output.root",
// // "cc_pro_LH2x10_CO2LA_BeamOnTarget_0.999e10_skim_output.root",
// 
// // "cc_pro_Bedown_C4F8LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",  //test
// // "cc_pro_He3_C4F8FA_BeamOnTarget_0.98e10_skim_output.root",
// //"cc_pro_LH2x10_C4F8LA_BeamOnTarget_1e10_skim_output.root",
// //"cc_pro_LH2x100_C4F8LA_BeamOnTarget_0.99e9_new_skim_output.root",
// // "cc_pro_LH2x100_C4F8LA_BeamOnTarget_0.952e10_skim_output.root",
// // "cc_pro_LH2x100_C4F8LA_BeamOnTarget_1.947e10_skim_output.root",
// 
// //"cc_pro_Bedown_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root", //test
// // "cc_pro_He3_CO2FA_BeamOnTarget_0.97e10_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_0.988e9_new_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_0.878e10_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_1.951e10_skim_output.root",
// 
// // "output_old/cc_pro_LH2x100_CO2LA_BeamOnTarget_0.878e10_skim_output.root",
// // "output_old/cc_pro_LH2x10_CO2LA_BeamOnTarget_0.999e10_skim_output.root",
// 
// // "cc_pro_He3_CO2LA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// 
// // "cc_pro_LH2x100_C4F8LA_BeamOnTarget_1.947e10_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_1.951e10_skim_output.root",
// // "new/cc_pro_LH2x100_C4F8LA_BeamOnTarget_0.981e10_HIGHMOM5MeV_skim_output.root",
// // "new/cc_pro_LH2x100_CO2LA_BeamOnTarget_0.97e10_HIGHMOM5MeV_skim_output.root",
// 
// // "cc_pro_He3_overall_C4F8LA_dirty_normalized_pi0Wiser_filenum100_1e6_skim_output.root",
// // "cc_pro_He3_overall_CO2LA_dirty_normalized_pi0Wiser_filenum98_0.98e6_skim_output.root",
// // "cc_pro_He3_overall_C4F8LA_dirty_weighted_eDIS_filenum79_0.79e6_skim_output.root",
// // "cc_pro_He3_overall_CO2LA_dirty_weighted_eDIS_filenum73_0.73e6_skim_output.root",
// 
// "cc_pro_He3_CO2FA_BeamOnTarget_0.99e10_skim_output.root",
// "cc_pro_He3_CO2FA_dirty_normalized_pi0Wiser_filenum98_0.98e9_skim_output.root",
// "cc_pro_He3_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// 
// // "cc_pro_He3_CO2LA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_0.988e9_new_skim_output.root",
// // "cc_pro_LH2x100_CO2LA_BeamOnTarget_1.951e10_skim_output.root",
// 
// // "cc_pro_LH2x100_C4F8LA_BeamOnTarget_0.99e9_new_skim_output.root",
// // "cc_pro_LH2x100_C4F8LA_BeamOnTarget_1.947e10_skim_output.root",
// 
// // "cc_pro_He3_CO2FA_dirty_weighted_eDIS_filenum98_0.98e9_skim_output.root",
// 
// // "cc_pro_Alup_CO2FA_dirty_weighted_eDIS_filenum100_1e9_skim_output.root",
// // "cc_pro_Alup_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// 
// // "cc_pro_He3_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_He3_CO2FA_dirty_normalized_pi0Wiser_filenum98_0.98e9_skim_output.root",
// // "cc_pro_He3_CO2FA_BeamOnTarget_0.99e10_skim_output.root",
// 
// };
// 
// // char* input_filename[m]={
// // "cc_pro_Beup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_Alup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_Aldown_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_Bedown_CO2FA_dirty_weighted_eDIS_filenum97_0.97e8_skim_output.root",
// // "cc_pro_He3_N2up_CO2FA_dirty_weighted_eDIS_filenum10_1e8_skim_output.root",
// // "cc_pro_He3_winup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_He3_CO2FA_dirty_weighted_eDIS_filenum98_0.98e8_skim_output.root",
// // "cc_pro_He3_N2down_CO2FA_dirty_weighted_eDIS_filenum10_1e8_skim_output.root",
// // "cc_pro_He3_windown_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_C1_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_C2_CO2FA_dirty_weighted_eDIS_filenum99_0.99e8_skim_output.root",
// // "cc_pro_C3_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_C_N2_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // 
// // "cc_pro_Beup_CO2FA_dirty_normalized_pi0Wiser_filenum98_0.98e9_skim_output.root",
// // "cc_pro_Alup_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // "cc_pro_Aldown_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_Bedown_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_He3_N2up_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // "cc_pro_He3_winup_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_He3_CO2FA_dirty_normalized_pi0Wiser_filenum98_0.98e9_skim_output.root",
// // "cc_pro_He3_windown_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_He3_N2down_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // "cc_pro_C1_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_C2_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // "cc_pro_C3_CO2FA_dirty_normalized_pi0Wiser_filenum99_0.99e9_skim_output.root",
// // // "cc_pro_C_N2_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e9_skim_output.root",
// // 
// // // "cc_pro_He3_CO2FA_BeamOnTarget_0.97e10_skim_output.root",
// // // "cc_pro_C_CO2FA_BeamOnTarget_1e9_skim_output.root",
// // 
// // "cc_pro_He3_CO2FA_dirty_weighted_eDIS_filenum98_0.98e9_skim_output.root",
// // "cc_pro_He3_CO2FA_dirty_normalized_pi0Wiser_filenum98_0.98e9_skim_output.root",
// // "cc_pro_He3_CO2FA_BeamOnTarget_0.99e10_skim_output.root",
// // };
// 
// // char *title="Cherenkov";
// // char *title="Signal Sum Distribution";
// 
// char *label[m]={
// "eDIS (Be up)","eDIS (Al up)","eDIS (Al down)","eDIS (Be down)",
// "eDIS (N2 up)","eDIS (win up)",
// "pi0Wiser (C4F8)",
// "eDIS (win down)","eDIS (N2 down)",
// "eDIS (C1)","eDIS (C2)","eDIS (C3)","eDIS (C_N2)",
// 
// "pi0Wiser (Be up)","pi0Wiser (Al up)","pi0Wiser (Al down)","pi0Wiser (Be down)",
// "pi0Wiser (N2 up)","pi0Wiser (win up)",
// "pi0Wiser (CO2)",
// "pi0Wiser (win down)","pi0Wiser (N2 down)",
// "pi0Wiser (C1)","pi0Wiser (C2)",//"pi0Wiser (C3)",//"pi0Wiser (C_N2)",
// 
// // "BeamOnTarget (He3)","BeamOnTarget (C)",
// // "BeamOnTarget (LH2x10/10)","BeamOnTarget (LH2x100/100)",
// // "BeamOnTarget","BeamOnTarget",
// // "e^{-}","#pi^{0}","BeamOnTarget"
// // "C4F8 BeamOntarget","CO2 BeamOntarget","C4F8 pi0Wiser","CO2 pi0Wiser",
// "CO2 BeamOntarget","CO2 pi0Wiser","CO2 eDIS",""
// };
// 
// int MarkerStyle[m]={
// 24,24,24,24,24,24,24,24,24,24,24,24,24,
// 26,26,26,26,26,26,26,26,26,26,26,26,26,
// 25,25,
// };
// int color[m]={
// // 1,2,3,4,5,6,7,8,9,10,11,12,13,
// // 1,2,3,4,5,6,7,8,9,10,11,12,13,
// // 1,7,3,4,5,6,2,8,9,10,11,12,13,
// // 1,7,3,4,5,6,2,8,9,10,11,12,13,
// // 2,2,
// 1,7,3,4,5,6,3,8,9,10,11,12,13,
// 1,7,3,4,5,6,2,8,9,10,11,//12,//13,
// // 1,2,
// 2,1,2,1,
// };
// int style[m]={
// 1,1,1,1,1,1,1,1,1,1,1,1,1,
// 2,2,2,2,2,2,2,2,2,2,2,//2,2,
// 1,1,2,2,
// };
// // char *hst[m]={"occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0", "occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0",};
// // char *hst[m]={"hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_3","hit_hgc_2",};
// // double ymin=1e-5,ymax=1e3;
// // double ymin=1e-2,ymax=1e6;
// // char *hst[m]={"pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3","pe_hgc_3"};
// char *hst[m]={"npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1",};
// double xmin=0,xmax=100;
// double ymin=1e-8,ymax=2e1;
// // double ymin=1e-8,ymax=1e6;
// // double ymin=1e0,ymax=1e7;
// // double ymin=1e-8,ymax=1e6;
// char *title="";
// 
// // char *hst[m]={"motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP","motherP",};
// // double xmin=0,xmax=0.1;
// // double ymin=1e-4,ymax=1e2;
// // char *title="electron P;GeV;rate(kHz)";
// 
// // char *hst[m]={"hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele","hit_Eec_ele"};
// // double ymin=1e-6,ymax=1e6;
// // double ymin=1e-3,ymax=1e9;

// pass5
// char input_dir[300]="/work/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_1.3/pass5/";
// // char input_dir[300]="/volatile/halla/solid/sim/solid_gemc/cc_pro_JLAB_VERSION_devel/pass5/";
// 
// // const int m=11;
// const int m=10;
// 
// char* input_filename[m]={
// 
// // "cc_pro_NOtarget_Beup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Alup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_N2_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Bedown_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Aldown_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// // 
// // "cc_pro_NOtarget_Beup_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Alup_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_N2_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Aldown_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// // "cc_pro_NOtarget_Bedown_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// // 
// // "cc_pro_NOtarget_CO2FA_BeamOnTarget_1e10_skim_output.root",
// 
// // "cc_pro_NOtarget_CO2FA_BeamOnTarget_1e10_skim_output.root",
// // "cc_pro_NOtarget_CO2FA_BeamOnTarget_1e10_filter_output.root",
// // "cc_pro_NOtarget_CO2FA_HBfield_BeamOnTarget_1e10_filter_output.root",
//   
// "cc_pro_NOtarget_CO2FA_BeamOnTarget_1e10_skim_output.root",
// 
// "cc_pro_NOtarget_Beup_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// "cc_pro_NOtarget_N2_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// "cc_pro_NOtarget_Bedown_CO2FA_dirty_weighted_eDIS_filenum100_1e8_skim_output.root",
// 
// "cc_pro_NOtarget_Beup_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// "cc_pro_NOtarget_N2_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// "cc_pro_NOtarget_Bedown_CO2FA_dirty_normalized_pi0Wiser_filenum100_1e8_skim_output.root",
// 
// "cc_pro_NOtarget_Beup_CO2FA_dirty_normalized_allBggen_filenum300_3e7_skim_output.root",
// "cc_pro_NOtarget_N2_CO2FA_dirty_normalized_allBggen_filenum300_3e7_skim_output.root",
// "cc_pro_NOtarget_Bedown_CO2FA_dirty_normalized_allBggen_filenum300_3e7_skim_output.root",
// 
// 
// };
// // char *title="Cherenkov";
// // char *title="Signal Sum Distribution";
// 
// char *label[m]={
//   "BeamOnTarget",
//   
// // "eDIS (Be up)","eDIS (Al up)","eDIS (N2)","eDIS (Al down)","eDIS (Be down)",
// 
// // "pi0Wiser (Be up)","pi0Wiser (Al up)","pi0Wiser (N2)","pi0Wiser (Al down)","pi0Wiser (Be down)",
// 
// "eDIS (Be up)","eDIS (N2)","eDIS (Be down)",
// 
// "pi0 Wiser (Be up)","pi0 Wiser (N2)","pi0 Wiser (Be down)",
// 
// "all Bggen (Be up)","all Bggen (N2)","all Bggen (Be down)",
// 
// };
// 
// int MarkerStyle[m]={
// 24,24,24,24,24,
// 24,24,24,24,24,
// // 24,
// };
// int color[m]={
// 1,7,3,4,5,6,3,8,9,
// 2,
// // 2
// };
// int style[m]={
// 1,1,1,1,1,1,1,1,1,1
// // 2,2,2,2,2,
// // 1,
// };
// // char *hst[m]={
// //   "occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0",
// //   "occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0","occ_hgc_0",
// //   "occ_hgc_0",
// // };
// char *hst[m]={
//   "hit_Eec","hit_Eec","hit_Eec",
//   "hit_Eec","hit_Eec","hit_Eec",
//   "hit_Eec","hit_Eec","hit_Eec","hit_Eec",
// //   "hit_Eec",
// };
// // char *hst[m]={"npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1","npe_hgc_1",};
// // double xmin=0,xmax=100;
// // double ymin=1e-8,ymax=2e1;
// // double ymin=1e-8,ymax=1e6;
// // double ymin=1e0,ymax=1e7;
// // double ymin=1e-8,ymax=1e6;
// // char *title="SC back;log10(Edep)(GeV);rate (kHz)";
// // char *title="SC front;log10(Edep)(GeV);rate (kHz)";
// char *title="EC;log10(Edep)(GeV);rate (kHz)";

char input_dir[300]="farm_cc_pro_even";

const int m=3;

char* input_filename[m]={
// "cc_pro_CO2LA_tcdonly_even_ecentral_1e3_output_MAPMT.root",
// "cc_pro_CO2LA_tcdonly_even_ecentral_1e3_output_LAPPD.root",
// "cc_pro_C4F8LA_tcdonly_even_ecentral_1e3_output_LAPPD.root"
"cc_pro_CO2LA_tcdonly_even_ecentral_1e3_output_MAPMT_new.root",
"cc_pro_CO2LA_tcdonly_even_ecentral_1e3_output_LAPPD_new.root",
"cc_pro_C4F8LA_tcdonly_even_ecentral_1e3_output_LAPPD_new.root"
};
// char *title="Cherenkov";
// char *title="Signal Sum Distribution";

char *label[m]={
"CO2 MAPMT","CO2 LAPPD","C4F8 LAPPD",
};

int MarkerStyle[m]={
24,24,24,
};
int color[m]={
1,2,3
};
int style[m]={
1,1,1
};
char *hst[m]={
  "npe_hgc_count_0","npe_hgc_count_0","npe_hgc_count_0",
};
char *title=";Npe;";

double occ_threshold_spd_FA=0.5,occ_threshold_spd_LA=3.; 
// double occ_threshold_ec_preshower=0.8,occ_threshold_ec_shower=12; //in MeV
double occ_threshold_ec_preshower=0.8,occ_threshold_ec_shower=12; //in MeV

///option
bool Is_R=true,Is_Phi=false;
// bool Is_logy=false;
bool Is_logy=true;
bool Is_logx=false;
// bool Is_logx=true;
bool Is_cut=false;

double xmax_npe[8]={10,20,30,50,50,50,80,80};

TCanvas *c = new TCanvas("c","c",1600,800);
c->Divide(1,1);
// c->Divide(4,2);
// TPaveLabel* ctitle = new TPaveLabel(0.1,0.96,0.9,0.99,"number of photo-electrons");
// ctitle->Draw();

// TLegend* leg = new TLegend(0.70, 0.98-0.03*m, 0.98, 0.98);
TLegend* leg = new TLegend(0.8, 0.95-0.07*m, 0.98, 0.95);
TFile *input[m];
TH1F *h[m];
TH1F *h_sum;

double total_He3_BeamOnTarget=0,total_C_BeamOnTarget=0;
double total_He3=0,total_C=0;
for(int i=0;i<m;i++){  
//    if (i!=26 && i!=6 && i!=19) continue;
//    if (i!=6 && i!=19) continue;  
//    if (i !=25 && i!=26) continue;    
//    if (i < 24 || i>26) continue;      
//    if (i<10) continue;      
  
  input[i]=new TFile(Form("%s/%s",input_dir,input_filename[i]));
  if (input[i]->IsZombie()) {
    cout << "Error opening file " << input_filename[i] << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename[i] << endl;
 
// for(int j=0;j<8;j++){    
//   c->cd(j+1);
  
  char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
  sprintf(hstname,"%s",hst[i]);      
//   sprintf(hstname,"%s_cut_%i",hst[i],j);
//   sprintf(hstname,"%s_count",hst[i],j); 
  cout << hstname << endl;
  h[i]=(TH1F*) input[i]->Get(hstname);

//   if (Is_cut) {
//     for(int j=0;j<h[i]->GetXaxis()->FindBin(cut[i]);j++)  h[i]->SetBinContent(j+1,0);
//   }
  
  if(Is_logy) gPad->SetLogy();
  if(Is_logx) gPad->SetLogx();  
//   h[i]->SetLineColor(color[i]);
//   h[i]->SetLineWidth(3);  
//   h[i]->SetLineStyle(style[i]);
//    h[i]->Scale(scale[i]);    
//   h[i]->SetMinimum(ymin);
//   h[i]->SetMaximum(ymax);
//   h[i]->SetAxisRange(xmin,xmax);
//   h[i]->SetAxisRange(xmin,xmax_npe[j]);  
//     h[i]->SetAxisRange(0,50);
//     h[i]->SetAxisRange(0,xmax);     
  h[i]->SetTitle(title);  
//   if (i==0) h[i]->Draw("C");
//   else h[i]->Draw("C same");

// h[i]->SetMarkerStyle(MarkerStyle[i]);
// h[i]->SetMarkerSize(2);
h[i]->SetMarkerColor(color[i]);
h[i]->SetLineColor(color[i]);
// h[i]->GetXaxis()->SetLabelSize(0.06);
// h[i]->GetYaxis()->SetLabelSize(0.06);

// string input_filename_s=input_filename[i];
// if(input_filename_s.find("LH2x10_",0) != string::npos) {h[i]->Scale(1./150.); cout << "scale0.007" << endl;}
// else if(input_filename_s.find("LH2x100_",0) != string::npos) {h[i]->Scale(1./1500.); cout << "scale0.0007" << endl;}
// else {}
  
//scale from 30uA to 1uA
// if(input_filename_s.find("FA",0) != string::npos) {h[i]->Scale(1./30.); cout << "scale0.033" << endl;}

  //remove spike CO2LA pi0 output_trigSC0.7EC0.1
//   if (i==25) h[i]->SetBinContent(3+1,h[i]->GetBinContent(3));
//   if (i==25) h[i]->SetBinContent(9+1,h[i]->GetBinContent(9));  

  //remove spike CO2LA pi0 output_trigno 
//   if (i==25) h[i]->SetBinContent(6+1,h[i]->GetBinContent(5));
//   if (i==25) h[i]->SetBinContent(7+1,h[i]->GetBinContent(5));  
//   if (i==25) h[i]->SetBinContent(17+1,h[i]->GetBinContent(17));  

//   if (i==6) h[i]->Draw("P");
// //   if (i==0) h[i]->Draw("P");
//   else h[i]->Draw("P same");

  if (i==0) h[i]->Draw("HIST");
  else h[i]->Draw("same HIST");

//      h[i]->Draw();
     
//   if (i==6) h[i]->Draw();
// //   if (i==0) h[i]->Draw();
//   else h[i]->Draw("same");
// 
//   if (i==6) h_sum=(TH1F) h[i]->Clone();
//   if (i==19) h_sum->Add(h[i]);
// 
// //   if (i==6 || i==19) h[i]->Scale(6);
//   if (i==19) h[i]->Scale(6);  
// 
//   double a=h[i]->Integral();
//   
//   if ((0<=i && i<=3) || (0<=13 && i<=16))   {total_He3 += a;total_C += a;}
//   if ((4<=i && i<=8) || (0<=17 && i<=21))  {total_He3 += a;}
//   if ((9<=i && i<=12) || (0<=22 && i<=25))  {total_C += a;}
//   if (i==26) total_He3_BeamOnTarget =a;
//   if (i==27) total_C_BeamOnTarget =a;
  
  cout << "GetEntries " << h[i]->GetEntries() << " Integral " << h[i]->Integral()/16. << endl;  
  //   input.Close(); 
  
//   leg->AddEntry(h[i], label[i],"p");    
  leg->AddEntry(h[i], label[i],"l");  
  
// }
		
}
leg->Draw();

// h_sum->SetLineColor(kRed);
// h_sum->Draw("same");

cout << "======================================== " <<  endl;
cout << "He3_BeamOnTarget avg " << total_He3_BeamOnTarget/16 <<  endl;
cout << "H3 avg " << total_He3/16 <<  endl;
cout << "C_BeamOnTarget avg " << total_C_BeamOnTarget/16 <<  endl;
cout << "C avg " << total_C/16 <<  endl;

c->SaveAs("c.png");
c->SaveAs("c.C");

}
