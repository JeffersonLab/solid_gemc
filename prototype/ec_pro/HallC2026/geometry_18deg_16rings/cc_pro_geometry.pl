use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;
use Getopt::Long;
# use Math::Trig;
# use Math::VectorReal;  

my $DetectorName = 'cc_pro';
my $DetectorMother="root";

# N.B.
## 1. All the lengths are in centimeters

# Constants
my $DEG=180/3.1415926;  # conversion factor between degrees and radians

# Parameters
## Chamber
my $AngX_tcd=0.0;
# my $AngY_tcd=0.0;
# my $AngY_tcd=-74.55;  #2020 hallc test, large angle 74.55deg at left side of beam direction, angle from the upstream beamline to the Cherenkov is about 105.45 with a hard to estimate error,  
# my $AngY_tcd=3.5; #2020 hallc test, small angle 3.5deg at right side of beam direction
#my $AngY_tcd=-82.2; #2022 hallc test, large angle at left side of beam direction, according to survey and alignment to target center is very good within 0.3deg
my $AngY_tcd=18; #2022 hallc test, large angle at left side of beam direction, according to survey and alignment to target center is very good within 0.3deg
#my $AngY_tcd=7; #2022 hallc test, large angle at left side of beam direction, according to survey and alignment to target center is very good within 0.3deg

# at 2020 hallc test
# measurement done at front scintilator plane ($sc1_r) which is 5" before chamber front window
# 17 feet 10 inch (543.6cm) to pivot at large angle, 39 feet (1189cm) to pivot at small angle, by Jack Seagal
# my $rmin_chamber=543.6+5*2.54;  # z position of the chamber entrance at large angle,
# my $rmin_chamber=1189+5*2.54;  # z position of the chamber entrance at small angle

# at 2022 hallc test
#my $rmin_chamber=1086.33;  # r position of the chamber entrance at large angle according to survey, sqrt(8.4369^2+1.09^2)=850.7
#my $rmin_chamber=1981.2;  # r position of the chamber entrance at large angle according to survey, sqrt(8.4369^2+1.09^2)=850.7
my $rmin_chamber=2021.67;  # r position of the chamber entrance at large angle according to survey,19.95-30.48+45./2*2.54 

# my $halflength_chamber_l = 60*2.54/2;
my $halflength_chamber_l = 56*2.54/2;
my $rmax_chamber=$rmin_chamber+$halflength_chamber_l*2;  # z position of the chamber exit
my $r_chamber=$rmin_chamber+$halflength_chamber_l; # z position of the chamber center and is tcd center

my $x_chamber = sin(-$AngY_tcd/$DEG)*$r_chamber;
my $y_chamber = 0.0;
my $z_chamber = cos($AngY_tcd/$DEG)*$r_chamber;

my $Rmin_chamber=0;  # inner radius of the chamber
# my $Rmax_chamber=14/2*2.54;  # outer radius of the chamber, 14 inch diameter, 
my $Rmax_chamber=(7.625+1.2)/2*2.54;  # outer radius of the chamber

## Gas
my $Rmin_gas=0;  # inner radius of the gas
# my $Rmax_gas=12.25/2*2.54;  # outer radius of the gas, 12.25 inch diameter
my $Rmax_gas=7.625/2*2.54;  # outer radius of the gas

## Windows
my $halfthickness_window_front=0.003/2*2.54; # half thickness of the front window
my $halfthickness_window_back=0.003/2*2.54;  # half thickness of the back window

my $Z_window_front=-$halflength_chamber_l+$halfthickness_window_front;
my $Z_window_back=$halflength_chamber_l-$halfthickness_window_back;
my $Rmin_window_front=$Rmin_gas;  # inner radius of the front windows
my $Rmax_window_front=$Rmax_gas;  # outer radius of the front windows
my $Rmin_window_back=$Rmin_gas;  # inner radius of the back windows
my $Rmax_window_back=$Rmax_gas;  # outer radius of the back windows

# my $r_window_front=$rmin_chamber+$halfthickness_window_front;  # r position of the front window
# my $r_window_back=$rmax_chamber-$halfthickness_window_back;   # r position of the back window

## Gas Cont.
my $rmin_gas=$rmin_chamber+$halfthickness_window_front*2;  # z position of the gas at the upstream side
my $rmax_gas=$rmax_chamber-$halfthickness_window_back*2;  # z position of the gas at the downstream side

## Mirror
my $mirror_Rmin=0;
my $mirror_Rmax=$Rmax_gas;
my $mirror_halfthickness=0.005;
my $mirror_d=44*2.54-$halflength_chamber_l;  # mirror at 44" front the chamber front window
# my $mirror_R=$Rmin_chamber+$mirror_d;

## PMT
my $PMT_y = 13.5; #  12*2.54-17
# my $half_width = 10.65;  # half width of the PMT
my $half_width = 5.30;  # half width of the PMT
my $windowhalf_z = 0.1;  # half thickness of the front side of PMT
my $backendhalf_z = 0.1;   # half thickness of the back side of PMT
my $half_z = $windowhalf_z + $backendhalf_z;  # total half length of the PMT

# my $halflength_chamber_t = ($PMT_y-14/2)*2.54;
# my $halflength_chamber_t = $PMT_y-$Rmax1_chamber;
my $halflength_chamber_t = $PMT_y;

# Hittype and materials
my $hitype="solid_hgc";  # alternative: "flux"
my $material_chamber="G4_POLYVINYL_CHLORIDE";
my $material_window_front = "G4_Al";  
my $material_window_back = "G4_Al";  
# my $material_window_back = "G4_POLYVINYL_CHLORIDE";  #actually, the back will be PVC
# my $material_gas="SL_HGC_C4F8_oneatm";  
# my $material_pmt_surface = "SL_HGC_C4F8_oneatm";
my $material_gas="SL_LGCCgas";
my $material_pmt_surface = "SL_LGCCgas";
# my $material_gas="SL_N2gas";
# my $material_pmt_surface = "SL_N2gas";
my $material_pmt_backend= "Kryptonite";
my $material_mirror= "SL_HGC_CFRP";
my $hy1 = 6*2.54/2;
my $hx1	= 30.48/2;
my $hdz1= 10.16/2;
my $hdz2= 30.48/2;
my $hdz3= 40.64/2;
my $hx2	= 10.16/2;
my $hy2	= 6*2.54/2;
my $hx3	= 10.16/2;
my $hy3	= 6*2.54/2;
my $hz1	= 23;
my $hz2	= 22.+10./2.+30.48/2.;
my $hz3	= 56*2.54/2+4.5+8.5+6.12+5.62+6.89+5.12+40.64/2;

sub cc_pro_geometry
{
make_tcd();
make_chamber();
make_gas();
make_window_front();
make_window_back();
make_mirror();
make_pmt();
make_coli();
make_coli2();
make_coli3();
}

sub make_tcd
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_tcd";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x_chamber*cm $y_chamber*cm $z_chamber*cm";
 $detector{"rotation"}    = "$AngX_tcd*deg $AngY_tcd*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Box";
 #$detector{"dimensions"}  = "35*cm 60*cm 250*cm";
 $detector{"dimensions"}  = "50*cm 60*cm 250*cm";
# $detector{"dimensions"}  = "70*cm 60*cm 300*cm";
 $detector{"material"}    = "G4_AIR";
 #$detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 0;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}

sub make_chamber
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_chamber_l";
 $detector{"mother"}      = "$DetectorName\_tcd";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC00";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_chamber*cm $Rmax_chamber*cm $halflength_chamber_l*cm 0*deg 360*deg"; 
 $detector{"material"}    = "Component";
 print_det(\%configuration, \%detector);
 
  %detector=init_det(); 
 $detector{"name"}        = "$DetectorName\_chamber_t";
 $detector{"mother"}      = "$DetectorName\_tcd";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm $PMT_y*cm $mirror_d*cm";
 $detector{"rotation"}    = "90*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC00";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_chamber*cm $Rmax_chamber*cm $halflength_chamber_t*cm 0*deg 360*deg"; 
 $detector{"material"}    = "Component";
 print_det(\%configuration, \%detector); 

  %detector=init_det();
 $detector{"name"}        = "$DetectorName\_chamber";
 $detector{"mother"}      = "$DetectorName\_tcd";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg -90*deg";
 $detector{"color"}       = "ff0000";
 $detector{"type"}        = "Operation:@ $DetectorName\_chamber_l + $DetectorName\_chamber_t";
 $detector{"dimensions"}  = "0";      
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
 $detector{"name"}        = "$DetectorName\_gas_l";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "AACC33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_gas*cm $Rmax_gas*cm $halflength_chamber_l*cm 0*deg 360*deg"; 
 $detector{"material"}    = "Component";
 print_det(\%configuration, \%detector);
 
 %detector=init_det(); 
 $detector{"name"}        = "$DetectorName\_gas_t";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm $PMT_y*cm $mirror_d*cm";
 $detector{"rotation"}    = "90*deg 0*deg 0*deg";
 $detector{"color"}       = "AACC33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_gas*cm $Rmax_gas*cm $halflength_chamber_t*cm 0*deg 360*deg"; 
 $detector{"material"}    = "Component";
 print_det(\%configuration, \%detector); 

  %detector=init_det();
 $detector{"name"}        = "$DetectorName\_gas";
 $detector{"mother"}      = "$DetectorName\_chamber";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "AACC33";
 $detector{"type"}        = "Operation:@ $DetectorName\_gas_l + $DetectorName\_gas_t";
 $detector{"dimensions"}  = "0";      
 $detector{"material"}    = "$material_gas";
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


sub make_window_front
{
    my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_window_front";
 $detector{"mother"}      = "$DetectorName\_gas";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $Z_window_front*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCAA33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_window_front*cm $Rmax_window_front*cm $halfthickness_window_front*cm 0*deg 360*deg";
 $detector{"material"}    = $material_window_front;
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
  $detector{"sensitivity"} = "flux";
  $detector{"hit_type"}    = "flux";
  $detector{"identifiers"} = "id manual 10";
  print_det(\%configuration, \%detector);
}

sub make_window_back
{
    my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_window_back";
 $detector{"mother"}      = "$DetectorName\_gas";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $Z_window_back*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCAA33";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin_window_back*cm $Rmax_window_back*cm $halfthickness_window_back*cm 0*deg 360*deg";
 $detector{"material"}    = $material_window_back;
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

sub make_mirror
{
    my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_mirror";
 $detector{"mother"}      = "$DetectorName\_gas";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $mirror_d*cm";
 $detector{"color"}       = "aaaaff";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"type"}        = "CTube";
 $detector{"dimensions"}  = "$mirror_Rmin*cm $mirror_Rmax*cm $mirror_halfthickness*cm 0*deg 360*deg 0*cm 0.7*cm -0.7*cm 0*cm -0.7*cm 0.7*cm"; 
 $detector{"material"}    = $material_window_back;
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

sub make_pmt
{ 
      my %detector;      
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt";
      $detector{"mother"}      = "$DetectorName\_gas";
      $detector{"description"} = $detector{"name"};
      $detector{"pos"}         = "0*cm $PMT_y*cm $mirror_d*cm";
      $detector{"rotation"}    = "-90*deg 0*deg 0*deg";
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
      $detector{"name"}        = "$DetectorName\_pmt_surface";     
      $detector{"mother"}      = "$DetectorName\_pmt";
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
      my $id=2200000+10000;      
      $detector{"identifiers"} = "id manual $id";      
      print_det(\%configuration, \%detector);
   
   
      %detector=init_det();      
      $detector{"name"}        = "$DetectorName\_pmt_backend";
      $detector{"mother"}      = "$DetectorName\_pmt";      
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
sub make_coli
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli";
#	$detector{"mother"}      = "$DetectorMother";
        $detector{"mother"}      = "$DetectorName\_tcd";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "-28.5*cm 0*cm $hz1*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "1a4fff";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx1*cm $hy1*cm $hdz1*cm";    
	$detector{"material"}    = "G4_Pb";
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
sub make_coli2
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli2";
#	$detector{"mother"}      = "$DetectorMother";
        $detector{"mother"}      = "$DetectorName\_tcd";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "18.5*cm 0*cm $hz2*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "1a4fff";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx2*cm $hy2*cm $hdz2*cm";    
	$detector{"material"}    = "G4_Pb";
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
sub make_coli3
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli3";
        $detector{"mother"}      = "$DetectorName\_tcd";
#        $detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "22.5*cm 0*cm $hz3*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "1a4fff";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx3*cm $hy3*cm $hdz3*cm";    
	$detector{"material"}    = "G4_Pb";
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

cc_pro_geometry();
1;
