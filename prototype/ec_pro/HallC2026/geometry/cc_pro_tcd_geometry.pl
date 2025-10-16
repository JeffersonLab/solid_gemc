use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;
use Getopt::Long;
# use Math::Trig;
# use Math::VectorReal;

my $DetectorName = 'cc_pro_tcd';
my $DetectorMother="root";

# N.B.
## 1. All the lengths are in centimeters

# Constants
my $Pi=3.1415926;
my $DEG=180/3.1415926;  # conversion factor between degrees and radians

# Parameters
## Chamber
# at 2022 hallc test
#my $rmin_chamber=1086.33;  # r position of the chamber entrance at large angle according to survey, sqrt(8.4369^2+1.09^2)=850.7
#my $rmin_chamber=1981.2;  # r position of the chamber entrance at large angle according to survey, sqrt(8.4369^2+1.09^2)=850.7
my $rmin_chamber=2021.67;  # r position of the chamber entrance at large angle according to survey,19.95-30.48+45./2*2.54

# my $halflength_chamber_l = 60*2.54/2;
my $halflength_chamber_l = 56*2.54/2;
my $rmax_chamber=$rmin_chamber+$halflength_chamber_l*2;  # z position of the chamber exit
my $r_chamber=$rmin_chamber+$halflength_chamber_l; # z position of the chamber center and is tcd center

# my $rmin_tcd=$rmin_chamber-500;
# my $rmax_tcd=$rmax_chamber+500;

my $Ang_tcd=18;

my $Ntcd=16; # make number of copy
# my $Ntcd=1; # just make 1

sub make_tcd
{
 my %detector=init_det();
 for(my $id=1; $id<=$Ntcd; $id++){
    my $x_chamber = -sin($Ang_tcd/$DEG)*cos(($id-1)*360/$Ntcd/$DEG)*$r_chamber;
    my $y_chamber = -sin($Ang_tcd/$DEG)*sin(($id-1)*360/$Ntcd/$DEG)*$r_chamber;
    my $z_chamber = cos($Ang_tcd/$DEG)*$r_chamber;

    my $AngX_tcd=$Ang_tcd*sin(-($id-1)*(360/$Ntcd/$DEG));
    my $AngY_tcd=$Ang_tcd*cos(($id-1)*(360/$Ntcd/$DEG));

#     my $AngX_tcd=0;
#     my $AngY_tcd=0;
    my $AngZ_tcd=0;

    if($id == 1){$detector{"name"}        = "$DetectorName";}
    else {$detector{"name"}        = "$DetectorName\_$id";}

    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
     $detector{"pos"}         = "$x_chamber*cm $y_chamber*cm $z_chamber*cm";
#     $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "$AngX_tcd*deg $AngY_tcd*deg $AngZ_tcd*deg";
#     $detector{"rotation"}    = "ordered: zyx $AngX_tcd*deg $AngY_tcd*deg $AngZ_tcd*deg";
    $detector{"color"}       = "CCCC33";

    if($id == 1){
    $detector{"type"}        = "Box";
    $detector{"dimensions"}  = "60*cm 60*cm 500*cm";
#     $detector{"type"}        = "Polycone";
#     $detector{"dimensions"}  = "0*deg 360*deg 2*counts 0*cm 0*cm 100*cm 100*cm $rmin_tcd*cm $rmax_tcd*cm";
    }
    else {$detector{"type"}       = "CopyOf $DetectorName";}

    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "G4_Galactic";
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

make_tcd();

1;
