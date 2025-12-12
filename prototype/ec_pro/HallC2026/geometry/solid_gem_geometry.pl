# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_gem';

my $DetectorMother="cc_pro_tcd";

sub solid_gem_geometry
{
make_gem();
make_coli();
}

#my $z1	= -56*2.54/2-30.48-8.26-7.03-9.3-1.5/2;
my $z1	= -56*2.54/2-8.26-7.03-9.-1.5/2;
my $z2	=  56*2.54/2+4.5+1.5/2;
#my $z3	= -56*2.54/2-30.48-8.26-7.03-1.5/2;
my $z3	= -56*2.54/2-8.26-7.03-1.5/2;
my $z4	=  56*2.54/2+4.5+8.5+1.5/2;
my $hx	= 5.12;
my $hy	= 5.12;
my $hx1	= 12*2.54/2.;
my $hy1	= 12*2.54/2.;
#my $cz1	= -56*2.54/2-8.26-7.03-9.-5.5-2.54/2;
my $cz1	= -56*2.54/2-8.26-7.03-9.-4;

my @PlateZ = ($z1,$z3,$z2,$z4);

my @PlateZ_virt = ($z1-1,$z3-1,$z2-1,$z4-1);

sub make_gem
{
# == older version ====================================
#  my $Dz   = 0.48;
#  my $material="DCgas";
#  my $color="44ee11";
#  my $Dz   = 15.955/2;

# ==== old copy from SBS code ================================
#  * Describe the single GEM Chamber module (similar to COMPASS)
#  * see: "Construction Of GEM Detectors for the COMPASS experiment", CERN Tech Note TA1/00-03
#  *
#  * Consist of 15 layers of different size, material and position
#  *
#  *
#  * HoneyComb
#  *  0   NEMA G10 120 um
#  *  1   NOMEX    3 mm  #should be 3um?
#  *  2   NEMA G10 120 um
#  * Drift Cathode
#  *  3   Copper 5 um    #should exchange with 4?
#  *  4   Kapton 50 um   #should exchange with 3?
#  *  5   Air 3 mm
#  * GEM0
#  *  6   Copper 5 um
#  *  7   Kapton 50 um
#  *  8   Copper 5 um
#  *  9   Air 2 mm
#  * GEM1
#  * 10   Copper 5 um
#  * 11   Kapton 50 um
#  * 12   Copper 5 um
#  * 13   Air 2 mm
#  * GEM2
#  * 14   Copper 5 um
#  * 15   Kapton 50 um
#  * 16   Copper 5 um
#  * 17   Air 2 mm 
#  * Readout Board
#  * 18   Copper 10 um
#  * 19   Kapton 50 um
#  * 20   G10 120 um + 60 um (assume 60 um glue as G10)    # not implmented yet
#  * Honeycomb
#  * 21   NEMA G10 120 um
#  * 22   NOMEX    3 mm       #should be 3um?
#  * 23   NEMA G10 120 um
# 
# my $Dz   = 9.781/2; # unit in mm 
#  my $Nlayer = 23;
#  my @layer_thickness = (0.12,0.003,0.12,0.05,0.005,3,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.01,0.05,0.12,0.003,0.12);   # unit in mm
#  my @material = ("SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_NEMAG10","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_NEMAG10");
#  my @sens = ("no","no","no","no","$hittype","$hittype","$hittype","no","no","no","no","no","no","no","no","no","no","no","$hittype","no","no","no","no");
#  my @hitt = ("no","no","no","no","$hittype","$hittype","$hittype","no","no","no","no","no","no","no","no","no","no","no","$hittype","no","no","no","no");
#  my $color_NEMAG10 = "00ff00";
#  my $color_NOMEX = "ffse14";
#  my $color_Copper = "ffe731";
#  my $color_Kapton = "1a4fff";
#  my $color_Air = "ff33fc";
#  my @color = ($color_NEMAG10,$color_NOMEX,$color_NEMAG10,$color_Kapton,$color_Copper,$color_Air,$color_Copper,$color_Kapton,$color_Copper,$color_Air,$color_Copper,$color_Kapton,$color_Copper,$color_Air,$color_Copper,$color_Kapton,$color_Copper,$color_Air,$color_Copper,$color_Kapton,$color_NEMAG10,$color_NOMEX,$color_NEMAG10);
# ======================================

# =====SoLID GEM =================================
# refer  to https://hallaweb.jlab.org/wiki/index.php/Solid_Tracking#GEM_module.27s_geometry_and_material

my $Nlayer = 24;
my @layer_thickness =
(0.005,0.02,3,0.05,0.005,3,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.01,0.05,0.18,3,0.02,0.005);  # unit in mm
my $Dz   = 0;
$Dz += $_ for @layer_thickness;
$Dz = int($Dz*5000)/10000.0;  # half thickness in mm

my @material =
("G4_Al","SL_gem_mylar","SL_gem_GEMgas","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","SL_gem_NEMAG10","SL_gem_GEMgas","SL_gem_mylar","G4_Al");
my $color_NEMAG10 = "00ffdd";
my $color_Al = "00ff00";
my $color_mylar = "ffse14";
my $color_Copper = "ffe731";
my $color_Kapton = "1a4fff";
my $color_gas = "ff33fc"; 
my @color = ($color_Al,$color_mylar,$color_gas,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_NEMAG10,$color_gas,$color_mylar,$color_Al);
my @hittype = ("no","no","no","no","solid_gem","solid_gem","solid_gem","no","solid_gem","solid_gem","solid_gem","no","no","no","no","no","no","no","solid_gem","no","no","no","no","no");

 for(my $n=1; $n<=4; $n++)
 {
     my %detector=init_det();

    $detector{"name"}        = "$DetectorName\_virt_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ_virt[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "111111";
    $detector{"type"}       = "Box";
    $detector{"dimensions"} = "$hx*cm $hy*cm 0.01*mm";
    $detector{"material"}    = "G4_Galactic";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    my $ID = 1+$n;
    $detector{"identifiers"} = "id manual $ID";
    print_det(\%configuration, \%detector);

    $detector{"name"}        = "$DetectorName\_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "111111";
    $detector{"type"}       = "Box";
    $detector{"dimensions"} = "$hx*cm $hy*cm $Dz*mm";    
    $detector{"material"}   = "G4_AIR";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = $n;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

    for( my $sec = 1; $sec <= 1; $sec++ ){

# was
#      my $thisrot = -($sec-1)*12.0 + $offrot[$n-1];

# If offset is 0,
#   sector_number=int((phi-90)/$angle_width+1)
# So 
# sector 1, phi(90 to 90+$angle_width)deg
# sector 2, phi(90+$angle_width to 90+2*$angle_width)deg
# .............
# sector $Nsect, phi(450-$angle_width to 450)deg

# $thisrot is central angle of sector (with a minus sign)

#       my $thisrot = -(90.0 - $angle_start + ($sec-1) * $angle_width + $offrot[$n-1]);
      my $thisrot = 0;
      
      for(my $i=1; $i<=$Nlayer; $i++)
      {
	my $layerZ = -$Dz;
	for(my $k=1; $k<=$i-1; $k++)
	{	
	   $layerZ = $layerZ+$layer_thickness[$k-1];
	}
	$layerZ = $layerZ+$layer_thickness[$i-1]/2;
	
	my $DlayerZ=$layer_thickness[$i-1]/2;

	    my %detector=init_det();
	    $detector{"name"}        = "$DetectorName\_$n\_$sec\_$i";
	    $detector{"mother"}      = "$DetectorName\_$n";
	    $detector{"description"} = $detector{"name"};
	    $detector{"pos"}        = "0*cm 0*cm $layerZ*mm";
	    $detector{"rotation"}   = "0*deg 0*deg $thisrot*deg";
	    $detector{"color"}      = "$color[$i-1]";
	    $detector{"type"}       = "Box";
	    $detector{"dimensions"} = "$hx*cm $hy*cm $DlayerZ*mm";	    
	    $detector{"material"}   = "$material[$i-1]";
	    $detector{"mfield"}     = "no";
	    $detector{"ncopy"}      = 1;
	    $detector{"pMany"}       = 1;
	    $detector{"exist"}       = 1;
	    $detector{"style"}       = 1;
	    $detector{"sensitivity"} = "$hittype[$i-1]";
	    $detector{"hit_type"}    = "$hittype[$i-1]";
	    my $id=1000000+$n*100000+$sec*1000+$i;
	    $detector{"identifiers"} = "id manual $id";
	    print_det(\%configuration, \%detector);
	}
    }


 }
}
sub make_coli
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $cz1*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "ff0000";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx1*cm $hy1*cm 1.27*cm";    
	#$detector{"material"}    = "G4_Pb";
	$detector{"material"}    = "G4_POLYSTYRENE";
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

solid_gem_geometry();
1;
