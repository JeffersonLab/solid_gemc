use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;
use Math::VectorReal;
  
my $DetectorName = 'solid_SIDIS_hgc';

my $DetectorMother="root";

sub solid_SIDIS_hgc_geometry
{
make_chamber();
make_gas();
make_block();
make_window_front();
make_window_back();
make_cone();
make_pmt();
make_mirror();
}

#BaBar
# Z(306,396)
# Rin(96,104)
# Rout(265,265)

# CLEOv8
# Z(306,396)
# Rin(86,98)
# Rout(265,265)

my $DEG=180/3.1415916;

my $hitype="solid_hgc";
# my $hitype="flux";

my $material_chamber="G4_Al";
# my $material_gas="G4_Al";
my $material_gas="SL_HGC_C4F8O";
my $material_block="SL_HGC_C4F8O";
# my $material_block="G4_Al";
my $material_window_front_1 = "SL_HGC_kevlar";
my $material_window_front_2 = "SL_HGC_mylar";
my $material_window_back = "G4_Al";
# my $material_window_front_1 = "G4_Al";
# my $material_window_front_2 = "G4_Al";
my $material_cone= "G4_GLASS_PLATE";
# my $material_pmt_surface = "SL_HGC_pmt_surface";
# my $material_pmt_surface = "G4_GLASS_PLATE";
my $material_pmt_surface = "SL_HGC_C4F8O";
my $material_pmt_backend= "Kryptonite";
my $material_mirror= "SL_HGC_CFRP";

my $N=30; # number of sectors
my $ang_width=360/$N;  #in degree
# my $ang_tilt = 44.; #degree  #  this is tilt angle of cone and PMTs around the focus axis
my $ang_tilt = 65; #degree  #  this is tilt angle of cone and PMTs around the focus axis

my $sec_start = 96; #degree where the 1st sector starts

#  all in cm
my $image_x = 0.;
my $image_y = 239.7;
my $image_z = 325;

my $halfthickness_window_front_1=0.0215;
my $halfthickness_window_front_2=0.0065;
my $halfthickness_window_back=0.05;

# my $Rmin1_chamber=86;
my $Rmin1_chamber=80;
my $Rmax1_chamber=265;
# my $Rmin2_chamber=98;
my $Rmin2_chamber=94;
my $Rmax2_chamber=265;
my $Rmin1_gas=$Rmin1_chamber+0.5;
my $Rmax1_gas=$Rmax1_chamber-0.5;
my $Rmin2_gas=$Rmin2_chamber+0.5;
my $Rmax2_gas=$Rmax2_chamber-0.5;

my $Zmin_chamber=306;
my $Zmax_chamber=406;
my $Zmin_gas=$Zmin_chamber+($halfthickness_window_front_1*2+$halfthickness_window_front_2*2);
my $Zmax_gas=$Zmax_chamber-$halfthickness_window_back*2;

my $Z_window_front_1=$Zmin_chamber+$halfthickness_window_front_1;
my $Z_window_front_2=$Zmin_chamber+$halfthickness_window_front_1*2+$halfthickness_window_front_2;
my $Rmin_window_front=$Rmin1_gas;
my $Rmax_window_front=187;

my $Z_window_back=$Zmax_chamber-$halfthickness_window_back;
my $Rmin_window_back=$Rmin2_gas;
my $Rmax_window_back=$Rmax2_gas;

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
#  $detector{"dimensions"}  = "0*deg 360*deg 2*counts 83*cm 96*cm 265*cm 265*cm 306*cm 406*cm";
 $detector{"dimensions"}  = "0*deg 360*deg 2*counts $Rmin1_chamber*cm $Rmin2_chamber*cm $Rmax1_chamber*cm $Rmax2_chamber*cm $Zmin_chamber*cm $Zmax_chamber*cm";
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
#  $detector{"dimensions"}  = "0*deg 360*deg 2*counts 83*cm 96*cm 265*cm 265*cm 306*cm 406*cm";
 $detector{"dimensions"}  = "0*deg 360*deg 2*counts $Rmin1_gas*cm $Rmin2_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Zmin_gas*cm $Zmax_gas*cm";
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
 
 $detector{"name"}        = "$DetectorName\_window_front_2";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $Z_window_front_2*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC33CC";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_window_front*cm $Rmax_window_front*cm $halfthickness_window_front_2*cm 0*deg 360*deg";
 $detector{"material"}    = $material_window_front_2;
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
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_block_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm 0*cm 0*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "22CC33";  
      $detector{"type"}        = "Polycone";
      $detector{"dimensions"}  = "$sector_start*deg 0.000005*deg 2*counts $Rmin1_gas*cm $Rmin2_gas*cm $Rmax1_gas*cm $Rmax2_gas*cm $Zmin_gas*cm $Zmax_gas*cm";      
      $detector{"material"}    = $material_block;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 0;
      $detector{"style"}       = 0;
      print_det(\%configuration, \%detector);
    }
}

sub make_cone
{
  my $rmin_w_front_cone = 14.16; 
  my $rmax_w_front_cone = 14.73; 
  my $rmin_w_end = 22.04; 
  my $rmax_w_end = 22.61;
#   my $z_w_half = 19. + 0.1;
  my $z_w_half = 18.5;  

#   my $Dz = 18.;
  
  my $pos_r = $image_y-sin($ang_tilt/$DEG)*$z_w_half+10;  
  my $pos_z = $image_z+cos($ang_tilt/$DEG)*$z_w_half+10;
print "cone $pos_r $pos_z\n";
     for(my $i=1; $i<=$N; $i++){
      my $pos_x = $pos_r*cos(($i-1)*$ang_width/$DEG+$sec_start/$DEG);
      my $pos_y = $pos_r*sin(($i-1)*$ang_width/$DEG+$sec_start/$DEG);    
#       my $ang_xrot=-($ang_tilt*sin(($i-1)*$ang_width/$DEG));
#       my $ang_yrot=($ang_tilt*cos(($i-1)*$ang_width/$DEG));
      my $ang_zrot=-(($i+-1)*$ang_width+$sec_start);
      my $ang_xrot=0;
      my $ang_yrot=$ang_tilt; 

#       my $z1=$pos_z-$z_w_half-0.002;
#       my $z2=$pos_z-$z_w_half-0.001;
#       my $z3=$pos_z-$z_w_half;
#       my $z4=$pos_z+$z_w_half;      
      
      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_cone_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";  
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
      $detector{"color"}       = "3dff84";  
      $detector{"type"}        = "Cons";
      $detector{"dimensions"}  = "$rmin_w_front_cone*cm $rmax_w_front_cone*cm $rmin_w_end*cm $rmax_w_end*cm $z_w_half*cm 0*deg 360*deg";   
#       $detector{"type"}        = "Polycone";
#       $detector{"dimensions"}  = "0*deg 360*deg 4*counts 0*cm 0*cm $rmin_w_front_cone*cm $rmin_w_end*cm $rmin_w_front_cone*cm $rmin_w_front_cone*cm $rmax_w_front_cone*cm $rmin_w_end*cm $z1*cm $z2*cm $z3*cm $z4*cm";
      $detector{"material"}    = $material_cone;
      $detector{"mfield"}      = "no";
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;
      $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
      $detector{"hit_type"}    = "mirror";
      $detector{"identifiers"} = "no";      
      print_det(\%configuration, \%detector);      
    }
}

sub make_pmt
{ 
#  my $z_w_half = 18.5;
#  my $pos_z = -($z_w_half-0.5);
 
 my $half_width = 10.16; 
#  my $windowhalf_z = 0.254;
 my $windowhalf_z = 0.001;
 my $backendhalf_z = 0.001; 
 my $half_z = $windowhalf_z + $backendhalf_z; 
 
  my $pos_r = $image_y+sin($ang_tilt/$DEG)*$half_z+10;  
  my $pos_z = $image_z-cos($ang_tilt/$DEG)*$half_z+10;
  print "pmt $pos_r $pos_z\n";
     for(my $i=1; $i<=$N; $i++){
      my $pos_x = $pos_r*cos(($i-1)*$ang_width/$DEG+$sec_start/$DEG);
      my $pos_y = $pos_r*sin(($i-1)*$ang_width/$DEG+$sec_start/$DEG);    
#       my $ang_xrot=-($ang_tilt*sin(($i-1)*$ang_width/$DEG));
#       my $ang_yrot=($ang_tilt*cos(($i-1)*$ang_width/$DEG));
#       my $ang_zrot=-($i-1)*$ang_width;
      my $ang_zrot=-(($i-1)*$ang_width+$sec_start);
      my $ang_xrot=0;
      my $ang_yrot=$ang_tilt; 

#       my $posz=-(18.5-$half_z);

      my %detector;      
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt_$i";
#       $detector{"mother"}      = "$DetectorName\_cone_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
      $detector{"rotation"}    =  "ordered: zxy $ang_zrot*deg $ang_xrot*deg $ang_yrot*deg";
#       $detector{"pos"}         = "0*cm 0*cm $posz*cm";
#       $detector{"rotation"}    = "0*deg 0*deg 0*deg";
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
#       $detector{"sensitivity"} = "mirror: SL_HGC_mirror_pmt";
#       $detector{"hit_type"}    = "mirror";
#       $detector{"identifiers"} = "no";      
      print_det(\%configuration, \%detector); 
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_pmt_surface_$i";
#       $detector{"mother"}      = "$DetectorName\_gas";      
      $detector{"mother"}      = "$DetectorName\_pmt_$i";
      $detector{"description"} = $detector{"name"};
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
#       $detector{"rotation"}    = "$ang_xrot*deg $ang_yrot*deg 0*deg";
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
#       $detector{"sensitivity"} = "mirror: SL_HGC_mirror_pmt";
#       $detector{"hit_type"}    = "mirror";
#       $detector{"identifiers"} = "no";
      $detector{"sensitivity"} = $hitype;
      $detector{"hit_type"}    = $hitype;
      my $id=2200000+$i*10000;      
      $detector{"identifiers"} = "id manual $id";      
      print_det(\%configuration, \%detector);
   
#       %detector=init_det();      
#       $detector{"name"}        = "$DetectorName\_pmt_backend_$i";
# #       $detector{"mother"}      = "$DetectorName\_gas";
#       $detector{"mother"}      = "$DetectorName\_pmt_$i";      
#       $detector{"description"} = $detector{"name"};
# #       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
# #       $detector{"rotation"}    = "$ang_xrot*deg $ang_yrot*deg 0*deg";
#       $detector{"pos"}         = "0*cm 0*cm -$windowhalf_z*cm";
#       $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#       $detector{"color"}       = "000000";  #cyan
#       $detector{"type"}        = "Box";
#       $detector{"dimensions"}  = "$half_width*cm $half_width*cm $backendhalf_z*cm";
#       $detector{"material"}    = "Air_Opt";
#       $detector{"mfield"}      = "no";
#       $detector{"ncopy"}       = 1;
#       $detector{"pMany"}       = 1;
#       $detector{"exist"}       = 1;
#       $detector{"visible"}     = 1;
#       $detector{"style"}       = 1;
#       $detector{"sensitivity"} = "flux";
#       $detector{"hit_type"}    = "flux";
#       my $id=2200000+$i*10000;      
#       $detector{"identifiers"} = "id manual $id";
#       print_det(\%configuration, \%detector);
   
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt_backend_$i";
#       $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"mother"}      = "$DetectorName\_pmt_$i";      
      $detector{"description"} = $detector{"name"};
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
#       $detector{"rotation"}    = "$ang_xrot*deg $ang_yrot*deg 0*deg";
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
#=== cons of mirror ========================================  
#   // max and min polar angles in SIDIS in degree
#   my $Angle_in = 7.5;
  my $Angle_in = 7.3;
  my $Angle_out = 15.5; 

#   // z dist. between the target and front/back walls of the tank: 350 + 306 & 350 + 396
#   // add 15*cm to R_in because the mirror is sticking out
  my $Z_front = 656;
  my $Z_end = 756;
  my $Z_half_w = ($Z_end-$Z_front)/2.;  

  my $R_front_in = $Z_front*tan($Angle_in/$DEG);
  my $R_front_out = $Z_front*tan($Angle_out/$DEG);
  my $R_end_in = $Z_end*tan($Angle_in/$DEG);
  my $R_end_out = $Z_end*tan($Angle_out/$DEG);
print "$R_front_in $R_front_out $R_end_in $R_end_out\n";  
#   86.3640091200709 181.92483680982 99.5292544127646 209.657281445463

  my $ang_start=-0.5*$ang_width;  
#   my $ang_start=0;  

#   // this is where the center of the cone should stay w.r.t. (0,0,0)
#   // i.e. in the middle of the tank: (306+406)/2=356*cm
  my $cons_z=356;
#===============================================================  
  
#=== sphere of mirror ========================================  
  my $Z_mirror = 388;  
  my $T_M1 = 0.3;
  my $P1 = vector(0,0,-350);    
  my $P2 = vector($image_x, $image_y, $image_z);
  my $V0 = vector(0., sin(11.5/$DEG), cos(11.5/$DEG));
  
my $V;
my $V_theta;
#         // Calculate the incident vector w.r.t. the source (i.e. target)
        my $Vi = (($Z_mirror - $P1->z())/$V0->z()) * $V0;
$V=$Vi; 
$V_theta=atan(sqrt(($V->x()*$V->x()+$V->y()*$V->y())/($V->z()*$V->z())))*$DEG;
print "$V_theta\n";  #11.5
#         // Crossing point on mirror plane w.r.t. the origin of coordinates
        my $Pm = $P1 + $Vi;
$V=$Pm;
$V_theta=atan(sqrt(($V->x()*$V->x()+$V->y()*$V->y())/($V->z()*$V->z())))*$DEG;
print "$V_theta\n"; #21.15
#         // Reflected vector: P2 w.r.t. the origin of coordinates; Vr doesn't matter
        my $Vr = $P2 - $Pm;  
$V=$Vr; 	#54.8
$V_theta=atan(sqrt(($V->x()*$V->x()+$V->y()*$V->y())/($V->z()*$V->z())))*$DEG;
print "$V_theta\n";                   
#         // Calculate the unitory normal vector
        my $Vn = $Vr->norm() - $Vi->norm();
        $Vn = $Vn->norm();
$V=$Vn; 
$V_theta=atan(sqrt(($V->x()*$V->x()+$V->y()*$V->y())/($V->z()*$V->z())))*$DEG;
print "$V_theta\n";    #21.69                       
#         // Calculate Angle
        my $ang_cos = -($Vi->norm().$Vn);  
my $ang=acos($ang_cos)*$DEG;                 
print "$ang\n";   #33.19
#         // Radius
        my $R = 2./$ang_cos/(1./$Vr->length() + 1./$Vi->length());        
# print "$R\n";        
#         // Spherical Center w.r.t. the origin of coordinates
        my $Position = $Pm + $R * $Vn;

# 	my $pos_x_sphere=sprintf('%f',$Position->x());
# 	my $pos_y_sphere=sprintf('%f',$Position->y());
# 	my $pos_z_sphere=sprintf('%f',$Position->z());	

	my $pos_x_sphere=sprintf('%f',$Position->y());
	my $pos_y_sphere=sprintf('%f',$Position->x());
	my $pos_z_sphere=sprintf('%f',$Position->z());
# print "$pos_x_sphere $pos_y_sphere $pos_z_sphere \n";
# # print "$Position->x()";
        
#         G4cout <<"Radius of Spherical Mirror: "<< R/cm << " " << "and pos mirror 1 is: " 
#                << Position/cm  << " " << "vi_mag is: " << Vi.mag() << " " << "Z_mirror: " 
#                << Z_mirror << " " << "P1z is: " << P1.z() << " " << "v0z is: " 
#                << V0.z() << " ang_cos is: " << ang_cos << G4endl;        
#
# Radius of Spherical Mirror: 228.451 and pos mirror 1 is: (0,234.568,175.719) vi_mag is: 7531.19 Z_mirror: 3880 P1z is: -3500 v0z is: 0.979925 ang_cos is: 0.83689
# rfin is: 86.364 rfout is: 181.925 rein is: 99.5293 reout is: 209.657 z_half is: 50

#   sphere radius
  my $R_in = $R;
  my $R_out = $R + $T_M1;  
#===============================================================
  
#       // make a cone to intersect a sphere      
     for(my $i=1; $i<=$N; $i++){      
      my %detector=init_det();
      $detector{"name"}        = "$DetectorName\_mirror_cons_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $cons_z*cm";
#       $detector{"rotation"}    = "$ang_xrot*deg $ang_yrot*deg 0*deg";
      $detector{"pos"}         = "0*cm 0*cm $cons_z*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "808080";  #gray
      $detector{"type"}        = "Cons";
      $detector{"dimensions"}  = "$R_front_in*cm $R_front_out*cm $R_end_in*cm $R_end_out*cm $Z_half_w*cm $ang_start*deg $ang_width*deg";
      $detector{"material"}    = "Component";
#       $detector{"material"}    = "$material_pmt_backend";
#       $detector{"mfield"}      = "no";
#       $detector{"ncopy"}       = 1;
#       $detector{"pMany"}       = 1;
#       $detector{"exist"}       = 1;
#       $detector{"visible"}     = 1;
#       $detector{"style"}       = 1;      
      print_det(\%configuration, \%detector);
      
      %detector=init_det();
      $detector{"name"}        = "$DetectorName\_mirror_sphere_$i";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $pos_z*cm";
#       $detector{"pos"}         = "0*cm 0*cm $Z_mirror*cm";
#       $detector{"pos"}         = "$pos_x*cm $pos_y*cm $cons_z*cm";
#       $detector{"rotation"}    = "$ang_xrot*deg $ang_yrot*deg 0*deg";
      $detector{"pos"}         = "$pos_x_sphere*cm $pos_y_sphere*cm $pos_z_sphere*cm";
      $detector{"rotation"}    = "0*deg 0*deg 0*deg";
      $detector{"color"}       = "808080";  #gray
      $detector{"type"}        = "Sphere";
      $detector{"dimensions"}  = "$R_in*cm $R_out*cm 0*deg 360*deg 0*deg 90*deg";
      $detector{"material"}    = "Component";
#       $detector{"material"}    = "$material_pmt_backend";      
#       $detector{"ncopy"}       = 1;
#       $detector{"pMany"}       = 1;
#       $detector{"exist"}       = 1;
#       $detector{"visible"}     = 1;
#       $detector{"style"}       = 1;            
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
      $detector{"rotation"}    = "0*deg 0*deg $ang_zrot*deg";
      $detector{"color"}       = "808080";  #gray				
      $detector{"type"}        = "Operation:@ $DetectorName\_mirror_cons_$i * $DetectorName\_mirror_sphere_$i";
      $detector{"dimensions"}  = "0";
#       $detector{"material"}    = "Component";      
      $detector{"material"}    = "$material_mirror";      
      $detector{"ncopy"}       = 1;
      $detector{"pMany"}       = 1;
      $detector{"exist"}       = 1;
      $detector{"visible"}     = 1;
      $detector{"style"}       = 1;        
#       $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
#       $detector{"hit_type"}    = "mirror";
      $detector{"sensitivity"} = "mirror: SL_HGC_mirror";
      $detector{"hit_type"}    = "mirror";
      $detector{"identifiers"} = "id manual 888888";      
      print_det(\%configuration, \%detector);           
    }
}
