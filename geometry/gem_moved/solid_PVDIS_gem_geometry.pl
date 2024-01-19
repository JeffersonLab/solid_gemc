use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_gem';

my $DetectorMother="root";

sub solid_PVDIS_gem_geometry
{
make_gem();
}

my $Nplate	= $parameters{"Nplate"};
my $Nsect       = $parameters{"Nsect"};
# my $angle_start=-5;
# my $angle_width=10;
# my $angle_start=-6;
# my $angle_width=12;
my $angle_width = 360.0 / $Nsect;
my $angle_start = -$angle_width * 0.5;
my $PlateZ1	= $parameters{"PlateZ1"};
my $PlateZ2	= $parameters{"PlateZ2"};
my $PlateZ3	= $parameters{"PlateZ3"};
my $PlateZ4	= $parameters{"PlateZ4"};
my $PlateZ5	= $parameters{"PlateZ5"};
my $Rin1	= $parameters{"Rin1"};
my $Rin2	= $parameters{"Rin2"};
my $Rin3	= $parameters{"Rin3"};
my $Rin4	= $parameters{"Rin4"};
my $Rin5	= $parameters{"Rin5"};
my $Rout1	= $parameters{"Rout1"};
my $Rout2	= $parameters{"Rout2"};
my $Rout3	= $parameters{"Rout3"};
my $Rout4	= $parameters{"Rout4"};
my $Rout5	= $parameters{"Rout5"};
my $offrot1	= $parameters{"offrot1"}; # angular offsets in degrees (w.r.t. start at 90+$angle_width*$n for $n-th sector) 
my $offrot2	= $parameters{"offrot2"};
my $offrot3	= $parameters{"offrot3"};
my $offrot4	= $parameters{"offrot4"};
my $offrot5	= $parameters{"offrot5"};

my @PlateZ = ($PlateZ1,$PlateZ2,$PlateZ3,$PlateZ4,$PlateZ5);
my @Rin    = ($Rin1,$Rin2,$Rin3,$Rin4,$Rin5);
my @Rout   = ($Rout1,$Rout2,$Rout3,$Rout4,$Rout5);
my @offrot = ($offrot1,$offrot2,$offrot3,$offrot4,$offrot5);

#my @offrot = (3.5,3.0,3.0,2.5,2.5);  
 
sub make_gem
{

# == older version ====================================
#  my $Dz   = 0.48;
#  my $material="DCgas";
#  my $color="44ee11";
#  my $Dz   = 15.955/2;

# ==== old copy from SBS code ================================

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

my $color_NOMEX = "00ff11";
my $color_NEMAG10 = "00ffdd";
my $color_Al = "00ff00";
my $color_mylar = "ffse14";
my $color_Copper = "ffe731";
my $color_Kapton = "1a4fff";
my $color_gas = "ff33fc"; 

# add a missing G10 layer at top according to Xinzhan Bai
my $Nlayer = 25;
my @layer_thickness = (0.05,0.125,3,0.18,0.05,0.005,3,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.01,0.05,0.18,3,0.125,0.05);  # unit in mm
my @material = ("G4_Al","SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_NEMAG10","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_NEMAG10","G4_Al");
my @color = ($color_Al,$color_NEMAG10,$color_NOMEX,$color_NEMAG10,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_NEMAG10,$color_NOMEX,$color_NEMAG10,$color_Al);
my @hittype = ("no","no","no","no","no","solid_gem","solid_gem","solid_gem","no","solid_gem","solid_gem","solid_gem","no","no","no","no","no","no","no","solid_gem","no","no","no","no","no");

# my $Nlayer = 24;

# change top and bottom 3 layers to 50um Al + 125um G10 + 3mm honey comb according to SBS GEM experience
# my @layer_thickness = (0.05,0.125,3,0.05,0.005,3,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.01,0.05,0.18,3,0.125,0.05);  # unit in mm
# my @material = ("G4_Al","SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","SL_gem_NEMAG10","SL_gem_NOMEX","SL_gem_NEMAG10","G4_Al");
# my @color = ($color_Al,$color_NEMAG10,$color_NOMEX,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_NEMAG10,$color_NOMEX,$color_NEMAG10,$color_Al);

# my @layer_thickness = (0.005,0.02,3,0.05,0.005,3,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.005,0.05,0.005,2,0.01,0.05,0.18,3,0.02,0.005);  # unit in mm
# my @material = ("G4_Al","SL_gem_mylar","SL_gem_GEMgas","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","G4_Cu","SL_gem_GEMgas","G4_Cu","SL_gem_Kapton","SL_gem_NEMAG10","SL_gem_GEMgas","SL_gem_mylar","G4_Al");
# my @color = ($color_Al,$color_mylar,$color_gas,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_Copper,$color_gas,$color_Copper,$color_Kapton,$color_NEMAG10,$color_gas,$color_mylar,$color_Al);

# my @hittype = ("no","no","no","no","solid_gem","solid_gem","solid_gem","no","solid_gem","solid_gem","solid_gem","no","no","no","no","no","no","no","solid_gem","no","no","no","no","no");

my $Dz   = 0;
$Dz += $_ for @layer_thickness;
$Dz = int($Dz*5000)/10000.0;  # half thickness in mm

 for(my $n=1; $n<=$Nplate; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "111111";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz*mm 0*deg 360*deg";
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

    for( my $sec = 1; $sec <= $Nsect; $sec++ ){

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

      my $thisrot = -(90.0 - $angle_start + ($sec-1) * $angle_width + $offrot[$n-1]);
    
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
	    $detector{"type"}       = "Tube";
	    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $DlayerZ*mm $angle_start*deg $angle_width*deg";
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
