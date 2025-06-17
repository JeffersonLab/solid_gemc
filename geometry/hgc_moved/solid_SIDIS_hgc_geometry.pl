use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;
use Getopt::Long;
use Math::Trig;
# use Math::VectorReal;
use VectorReal_mine;  

my $DetectorName = 'solid_SIDIS_hgc';
my $DetectorMother="root";

sub solid_SIDIS_hgc_geometry
{
make_chamber();
make_gas();
# make_block();
make_blockmirror1();
make_window_front();
# make_window_back();
make_cone();
make_shieldback();
make_pmt();
make_mirror();
}

# N.B.
## 1. All the lengths are in centimeters

# Constants
my $N=30; # number of sectors
my $ang_width=360/$N;  # =12, phi angle each sector covers in degree
my $sec_start = 90+$ang_width/2;  # =96, degree where the 1st sector starts
my $DEG=180/3.1415926;  # conversion factor between degrees and radians
my $Z_target=350;  # distance between the target center and the origin

# Parameters
## Chamber
my $Rmin0_chamber=190;  
my $Rmax0_chamber=240;  
my $Rmin1_chamber=77;  # inner radius of the chamber at the upstream side
my $Rmax1_chamber=240;  # outer radius of the chamber at the upstream side
my $Rmin2_chamber=87; 
my $Rmax2_chamber=240;
my $Rmin3_chamber=225;  # inner radius of the chamber at the downstream side
my $Rmax3_chamber=240;  # outer radius of the chamber at the downstream side

my $Z0_chamber=316;  # z position of the chamber at the upstream side
my $Zmin_chamber=326;  # z position of the chamber at the upstream side
my $Zmid_chamber=416;  # z position of the chamber at the upstream side
my $Zmax_chamber=426;  # z position of the chamber at the downstream side

## Gas
my $Rmin0_gas=$Rmin0_chamber+2.54*1.25;
my $Rmax0_gas=$Rmax0_chamber-2.54;  
my $Rmin1_gas=$Rmin1_chamber+2.54*1.25;  # inner radius of the gas at the upstream side
my $Rmax1_gas=$Rmax1_chamber-2.54;  # outer radius of the gas at the upstream side
my $Rmin2_gas=$Rmin2_chamber+2.54*1.25; 
my $Rmax2_gas=$Rmax2_chamber-2.54; 
my $Rmin3_gas=$Rmin3_chamber+2.54*1.25;  # inner radius of the gas at the downstream side
my $Rmax3_gas=$Rmax3_chamber-2.54;  # outer radius of the gas at the downstream side

## Windows
# my $halfthickness_window_front_1=0.0215;  # half thickness of the 1st front window
# my $halfthickness_window_front_2=0.0065;  # half thickness of the 2nd front window
my $halfthickness_window_front_1=2.54*0.04/2;  # half thickness of the 1st front window
my $halfthickness_window_back=2.54*0.25/2;  # half thickness of the back window
my $Z_window_front_1=$Zmin_chamber+$halfthickness_window_front_1;  # z position of the 1st front window
# my $Z_window_front_2=$Zmin_chamber+$halfthickness_window_front_1*2+$halfthickness_window_front_2;   # z position of the 2nd front window
# my $Z_window_back=$Zmax_chamber-$halfthickness_window_back;   # z position of the back window
my $Z_window_back=$Zmax_chamber-$halfthickness_window_back;   # z position of the back window

my $Rmin_window_front=$Rmin1_gas;  # inner radius of the front windows
my $Rmax_window_front=187;  # outer radius of the front windows
my $Rmin_window_back=$Rmin2_gas;  # inner radius of the back windows
my $Rmax_window_back=$Rmax2_gas;  # outer radius of the back windows

## Gas Cont.
# my $Zmin_gas=$Zmin_chamber+($halfthickness_window_front_1*2+$halfthickness_window_front_2*2);  # z position of the gas at the upstream side
my $Z0_gas=$Z0_chamber+1;
my $Zmin_gas=$Zmin_chamber+$halfthickness_window_front_1*2;  # z position of the gas at the upstream side
my $Zmid_gas=$Zmid_chamber;  
my $Zmax_gas=$Zmax_chamber-$halfthickness_window_back*2;  # z position of the gas at the downstream side

## Cone
# my $ang_upper = 126;  # angle between the PMT and the upper side of the cone in degrees
# my $ang_lower = 126;  # angle between the PMT and the lower side of the cone in degrees
# my $halflength_upper = 12;  # half length of the upper side of the cone
# my $halflength_lower_extend = $halflength_upper*cos(($ang_upper-90)/$DEG)/cos(($ang_lower-90)/$DEG);  # extended half length of the lower side of the cone, used in the definition of the larger cone to be cut
# my $halflength_lower = 

## Mirror
### Radius of the mirror
my $sphere_r_inner = 210;  # inner radius of the mirror
my $thickness_sphere = 0.3;  # thickness of the mirror
my $sphere_r_outer = $sphere_r_inner+$thickness_sphere;  # outer radius of the mirror

### Cone used to cut the mirror
my $Z_front = $Z_target+$Zmin_chamber;  # z position of the upstream side of the cone used to cut the mirror
my $Z_end = $Z_target+$Zmax_chamber;  # z position of the downstream side of the cone used to cut the mirror
my $Z_half_w = ($Z_end-$Z_front)/2.;  # half width of the cone used to cut the mirror
my $cons_z=($Zmin_chamber+$Zmax_chamber)/2;  # z position of the center of the cone used to cut the mirror
my $ang_min = 6.8;  # inner angle of the cone used to cut the mirror
my $ang_max = 16;  # outer anlge of the cone used to cut the mirror
my $R_front_min = $Z_front*tan($ang_min/$DEG);  # inner radius of the cone at the upstream side
my $R_front_max = $Z_front*tan($ang_max/$DEG);  # outer radius of the cone at the upstream side
my $R_end_min = $Z_end*tan($ang_min/$DEG);  # inner radius of the cone at the downstream side
my $R_end_max = $Z_end*tan($ang_max/$DEG);  # outer radius of the cone at the downstream side

## PMT
my $ang_tilt = 39;  # tilt angle between the PMT surface and the radial direction in degrees
my $half_width = 10.65;  # half width of the PMT
my $windowhalf_z = 0.001;  # half thickness of the front side of PMT
my $backendhalf_z = 0.001;   # half thickness of the back side of PMT
my $half_z = $windowhalf_z + $backendhalf_z;  # total half length of the PMT

# Hittype and materials
my $hitype="solid_hgc";  # alternative: "flux"
my $material_chamber="G4_Al";
my $material_gas="SL_HGC_C4F8O";  # alternative: "G4_Al"
my $material_block="G4_Al";  # alternative: "G4_Al"
# my $material_window_front_1 = "SL_HGC_kevlar";  # alternative: "G4_Al"
# my $material_window_front_2 = "SL_HGC_mylar";  # alternative: "G4_Al"
my $material_window_front_1 = "G4_Al";
my $material_window_back = "G4_Al";
my $material_cone= "G4_Fe";
my $material_pmt_surface = "SL_HGC_C4F8O";  # alternative: "G4_GLASS_PLATE" or "SL_HGC_pmt_surface"
my $material_pmt_backend= "Kryptonite";
my $material_mirror= "SL_HGC_CFRP";

#for round cone
my $rmin_w_front_cone = 10.65; 
my $rmax_w_front_cone = $rmin_w_front_cone+0.57;
my $rmin_w_end = 16; 
my $rmax_w_end = $rmin_w_end+0.57;

#for trapezoid cone
my $trd_dx1_inner = $half_width + 0.5;
my $trd_dy1_inner = $half_width + 0.5;
my $trd_dx2_inner = 22.41;
my $trd_dy2_inner = 16;
my $trd_dx1_outer = $trd_dx1_inner + 0.5;
my $trd_dy1_outer = $trd_dy1_inner + 0.5;
my $trd_dx2_outer = $trd_dx2_inner + 0.5;
my $trd_dy2_outer = $trd_dy2_inner + 0.5;
my $trd_dz = 8.09;
my $shieldback_dz = 10;

my $dis_P_PMT = 135;
my $z_angle_min = 390;
my $z_angle_max = 420;

my $Angle_in = 7;
my $Angle_out = 15;
my $Angle_opt = 6.8;

my $y_angle_min = ($z_angle_min+$Z_target)*tan($Angle_in/$DEG);
my $y_angle_max = ($z_angle_max+$Z_target)*tan($Angle_out/$DEG);
my $pos_A = vector(0, $y_angle_max, $z_angle_max);
my $pos_B = vector(0, $y_angle_min, $z_angle_min);
my $vec_BA = $pos_A - $pos_B;
my $dir_BA = $vec_BA->norm();
my $dir_BA_ortho = vector(0, $dir_BA->z(),-$dir_BA->y());
my $dis_BA = sqrt($vec_BA->x()*$vec_BA->x()+$vec_BA->y()*$vec_BA->y()+$vec_BA->z()*$vec_BA->z());
my $dis_C_BA = sqrt($sphere_r_inner*$sphere_r_inner-$dis_BA*$dis_BA/4);
my $pos_sphere_center = $pos_B+$vec_BA/2+$dir_BA_ortho*$dis_C_BA;

my $dis_C_PT = ($pos_sphere_center->y()-tan($Angle_opt/$DEG)*$pos_sphere_center->z()-tan($Angle_opt/$DEG)*350)/sqrt(1+tan($Angle_opt/$DEG)*tan($Angle_opt/$DEG));
my $ang_CPT = asin($dis_C_PT/$sphere_r_inner);
my $dir_P_PMT = vector(0,sin(2*$ang_CPT-$Angle_opt/$DEG),-cos(2*$ang_CPT-$Angle_opt/$DEG));
my $pos_P = vector(0,$pos_sphere_center->y()-$sphere_r_inner*sin($ang_CPT-$Angle_opt/$DEG),$pos_sphere_center->z()+$sphere_r_inner*cos($ang_CPT-$Angle_opt/$DEG));             
my $pos_PMT = $pos_P+$dis_P_PMT*$dir_P_PMT;
my $image_y = $pos_PMT->y();
my $image_z = $pos_PMT->z();
my $image_x = 0;
print "PMT position in xyz of $image_x $image_y $image_z\n";
my $pos_x_sphere=sprintf('%f',$pos_sphere_center->y()); # somehow this needs to be swtiched
my $pos_y_sphere=sprintf('%f',$pos_sphere_center->x()); # somehow this needs to be swtiched
my $pos_z_sphere=sprintf('%f',$pos_sphere_center->z());
print "Position of the center of sphere in xyz of $pos_x_sphere $pos_y_sphere $pos_z_sphere\n";

sub make_chamber
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_chamber";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Polycone";
#  $detector{"dimensions"}  = "0*deg 360*deg 3*counts $Rmin1_chamber*cm $Rmin2_chamber*cm $Rmin3_chamber*cm $Rmax1_chamber*cm $Rmax2_chamber*cm $Rmax3_chamber*cm $Zmin_chamber*cm $Zmid_chamber*cm $Zmax_chamber*cm";
 $detector{"dimensions"}  = "0*deg 360*deg 5*counts $Rmin0_chamber*cm $Rmin0_chamber*cm $Rmin1_chamber*cm $Rmin2_chamber*cm $Rmin3_chamber*cm $Rmax0_chamber*cm $Rmax1_chamber*cm $Rmax1_chamber*cm $Rmax2_chamber*cm $Rmax3_chamber*cm $Z0_chamber*cm $Zmin_chamber*cm $Zmin_chamber*cm $Zmid_chamber*cm $Zmax_chamber*cm";
 $detector{"material"}    = "$material_chamber";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}

sub make_gas
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_gas";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Polycone";
#  $detector{"dimensions"}  = "0*deg 360*deg 3*counts $Rmin1_gas*cm $Rmin2_gas*cm $Rmin3_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Rmax3_gas*cm $Zmin_gas*cm $Zmid_gas*cm $Zmax_gas*cm";
 $detector{"dimensions"}  = "0*deg 360*deg 5*counts $Rmin0_gas*cm $Rmin0_gas*cm $Rmin1_gas*cm $Rmin2_gas*cm $Rmin3_gas*cm $Rmax0_gas*cm $Rmax1_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Rmax3_gas*cm $Z0_gas*cm $Zmin_gas*cm $Zmin_gas*cm $Zmid_gas*cm $Zmax_gas*cm";
 $detector{"material"}    = "$material_gas";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
#       $detector{"sensitivity"} = "flux";
#       $detector{"hit_type"}    = "flux";
#       $detector{"identifiers"} = "id manual 1";
 print_det(\%configuration, \%detector);
}


sub make_window_front
{
    my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_window_front_1";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $Z_window_front_1*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_window_front*cm $Rmax_window_front*cm $halfthickness_window_front_1*cm 0*deg 360*deg";
 $detector{"material"}    = $material_window_front_1;
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
 
#  $detector{"name"}        = "$DetectorName\_window_front_2";
#  $detector{"mother"}      = "$DetectorName\_chamber";
#  $detector{"description"} = $detector{"name"};
#  $detector{"pos"}         = "0*cm 0*cm $Z_window_front_2*cm";
#  $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#  $detector{"color"}       = "CC33CC";
#  $detector{"type"}        = "Tube";
#  $detector{"dimensions"}  = "$Rmin_window_front*cm $Rmax_window_front*cm $halfthickness_window_front_2*cm 0*deg 360*deg";
#  $detector{"material"}    = $material_window_front_2;
#  $detector{"mfield"}      = "no";
#  $detector{"ncopy"}       = 1;
#  $detector{"pMany"}       = 1;
#  $detector{"exist"}       = 1;
#  $detector{"visible"}     = 1;
#  $detector{"style"}       = 1;
#  $detector{"sensitivity"} = "no";
#  $detector{"hit_type"}    = "no";
#  $detector{"identifiers"} = "no";
#  print_det(\%configuration, \%detector);  
}

sub make_window_back
{
    my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_window_back";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $Z_window_back*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_window_back*cm $Rmax_window_back*cm $halfthickness_window_back*cm 0*deg 360*deg";
 $detector{"material"}    = $material_window_back;
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}

sub make_block
{
     for(my $i=1; $i<=$N; $i++){
      my $sector_start=$sec_start+0.5*$ang_width+$ang_width*($i-1);
      
      if ($i!=30 && $i!=15) {next;}
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_block_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm 0*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "22CC33";  
      $detector{"type"}        = "Polycone";
      $detector{"dimensions"}  = "$sector_start*deg 0.5*deg 2*counts $Rmin1_gas*cm $Rmin2_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Zmin_gas*cm $Zmax_gas*cm";
      $detector{"material"}    = $material_block;     
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;          
      print_det(\%configuration, \%detector);
    }
}

sub make_blockmirror1
{
     for(my $i=1; $i<=$N; $i++){
      my $sector_start=$sec_start+0.5*$ang_width+$ang_width*($i-1);
      
      if ($i!=30 && $i!=15) {next;}
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_block_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm 0*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "22CC33";  
      $detector{"type"}        = "Polycone";
      $detector{"dimensions"}  = "$sector_start*deg 0.5*deg 2*counts $Rmin1_gas*cm $Rmin2_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Zmin_gas*cm $Zmax_gas*cm";
      $detector{"material"}    = "$material_mirror";      
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;
      $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
      $detector{"hit_type"}    = "mirror";
      my $id=2200000+$i*1000+3;
      $detector{"identifiers"} = "id manual $id";            
      print_det(\%configuration, \%detector);
    }
}

sub make_cone
{   
  my $pos_r = $image_y-sin($ang_tilt/$DEG)*$trd_dz;  
  my $pos_z = $image_z+cos($ang_tilt/$DEG)*$trd_dz;
     for(my $i=1; $i<=$N; $i++){
      my $pos_x = $pos_r*cos(($i-1)*$ang_width/$DEG+$sec_start/$DEG);
      my $pos_y = $pos_r*sin(($i-1)*$ang_width/$DEG+$sec_start/$DEG);    
      my $ang_zrot=-(($i+-1)*$ang_width+$sec_start);
      my $ang_xrot=0;
      my $ang_yrot=$ang_tilt; 
      
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_coneout_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Trd";
      $detector{"dimensions"}  = "$trd_dx1_outer*cm $trd_dx2_outer*cm $trd_dy1_outer*cm $trd_dy2_outer*cm $trd_dz*cm";   
      $detector{"material"}    = "Component";
      print_det(\%configuration, \%detector);  

      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_conein_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Trd";
      $detector{"dimensions"}  = "$trd_dx1_inner*cm $trd_dx2_inner*cm $trd_dy1_inner*cm $trd_dy2_inner*cm $trd_dz*cm";   
      $detector{"material"}    = "Component";
      print_det(\%configuration, \%detector);
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_cone_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Operation:@ $DetectorName\_coneout_$i - $DetectorName\_conein_$i";
      $detector{"dimensions"}  = "0";   
      $detector{"material"}    = $material_cone;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;
      $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
      $detector{"hit_type"}    = "mirror";
      my $id=2200000+$i*1000+2;
      $detector{"identifiers"} = "id manual $id";      
      print_det(\%configuration, \%detector);

               
    }
}

sub make_shieldback
{   
  my $pos_r = $image_y-sin($ang_tilt/$DEG)*(-$shieldback_dz-0.1);  
  my $pos_z = $image_z+cos($ang_tilt/$DEG)*(-$shieldback_dz-0.1);
     for(my $i=1; $i<=$N; $i++){
      my $pos_x = $pos_r*cos(($i-1)*$ang_width/$DEG+$sec_start/$DEG);
      my $pos_y = $pos_r*sin(($i-1)*$ang_width/$DEG+$sec_start/$DEG);    
      my $ang_zrot=-(($i+-1)*$ang_width+$sec_start);
      my $ang_xrot=0;
      my $ang_yrot=$ang_tilt; 
      
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_shieldbackout_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Box";
      $detector{"dimensions"}  = "$trd_dx1_outer*cm $trd_dy1_outer*cm $shieldback_dz*cm";        
      $detector{"material"}    = "Component";
      print_det(\%configuration, \%detector);  

      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_shieldbackin_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Box";
      $detector{"dimensions"}  = "$trd_dx1_inner*cm $trd_dy1_inner*cm $shieldback_dz*cm";        
      $detector{"material"}    = "Component";
      print_det(\%configuration, \%detector);
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_shieldback_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Operation:@ $DetectorName\_shieldbackout_$i - $DetectorName\_shieldbackin_$i";
      $detector{"dimensions"}  = "0";   
      $detector{"material"}    = $material_cone;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 0;
      $detector{"sensitivity"} = "no";
      $detector{"hit_type"}    = "no";
      $detector{"identifiers"} = "no";      
      print_det(\%configuration, \%detector);
    }
}

sub make_pmt
{ 
  my $pos_r = $image_y+sin($ang_tilt/$DEG)*$half_z;  
  my $pos_z = $image_z-cos($ang_tilt/$DEG)*$half_z;
     for(my $i=1; $i<=$N; $i++){
      my $pos_x = $pos_r*cos(($i-1)*$ang_width/$DEG+$sec_start/$DEG);
      my $pos_y = $pos_r*sin(($i-1)*$ang_width/$DEG+$sec_start/$DEG);    
      my $ang_zrot=-(($i-1)*$ang_width+$sec_start);
      my $ang_xrot=0;
      my $ang_yrot=$ang_tilt; 

      my %detector;      
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "000000";  #cyan
      $detector{"type"}        = "Box";
      $detector{"dimensions"}  = "$half_width*cm $half_width*cm $half_z*cm";
      $detector{"material"}    = $material_pmt_surface;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 0; 
      print_det(\%configuration, \%detector); 
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_pmt_surface_$i";     
      $detector{"mother"}      = "$DetectorName\_pmt_$i";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm $windowhalf_z*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "009999";  #cyan
      $detector{"type"}        = "Box";
      $detector{"dimensions"}  = "$half_width*cm $half_width*cm $windowhalf_z*cm";
      $detector{"material"}    = $material_pmt_surface;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;
      $detector{"sensitivity"} = $hitype;
      $detector{"hit_type"}    = $hitype;
#       my $id=2200000+$i*10000;      old and bad 
      my $id=2200000+$i*1000;      
      $detector{"identifiers"} = "id manual $id";      
      print_det(\%configuration, \%detector);
   
   
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt_backend_$i";
      $detector{"mother"}      = "$DetectorName\_pmt_$i";      
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm -$backendhalf_z*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "000000";  #cyan
      $detector{"type"}        = "Box";
      $detector{"dimensions"}  = "$half_width*cm $half_width*cm $backendhalf_z*cm";
      $detector{"material"}    = $material_pmt_backend;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;
      print_det(\%configuration, \%detector);
    }    
}

sub make_mirror
{
  my $ang_start=-0.5*$ang_width;  
  
#       // make a cone to intersect a sphere      
     for(my $i=1; $i<=$N; $i++){      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_mirror_cons_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm $cons_z*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "808080";  #gray
      $detector{"type"}        = "Cons";
      $detector{"dimensions"}  = "$R_front_min*cm $R_front_max*cm $R_end_min*cm $R_end_max*cm $Z_half_w*cm $ang_start*deg $ang_width*deg";
      $detector{"material"}    = "Component";     
      print_det(\%configuration, \%detector);
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_mirror_sphere_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x_sphere*cm $pos_y_sphere*cm $pos_z_sphere*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "808080";  #gray
      $detector{"type"}        = "Sphere";
      $detector{"dimensions"}  = "$sphere_r_inner*cm $sphere_r_outer*cm 0*deg 360*deg 0*deg 90*deg";
      $detector{"material"}    = "Component";            
      print_det(\%configuration, \%detector);      
      
      my $ang_zrot = -($sec_start+($i-1)*$ang_width);      
      
# Make the subtraction of the inner ellipsoid from the outer barrel:
# the "Operation:@" indicates that gemc should assume the coordinates
# and rotations of the mirror ellipsoid are given in its mother coordinate system,
# not relative to the outer barrel coordinate system:      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_mirror_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm $cons_z*cm";      
#       $detector{"rotation"}    = "0*deg 0*deg $ang_zrot*deg";
      $detector{"rotation"}    = "0*deg 0*deg $ang_zrot*deg";
      $detector{"color"}       = "808080";  #gray				
      $detector{"type"}        = "Operation:@ $DetectorName\_mirror_cons_$i * $DetectorName\_mirror_sphere_$i";
      $detector{"dimensions"}  = "0";      
      $detector{"material"}    = "$material_mirror";      
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;        
      $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
      $detector{"hit_type"}    = "mirror";
      my $id=2200000+$i*1000+1;      
      $detector{"identifiers"} = "id manual $id";      
      print_det(\%configuration, \%detector);  
      }
}