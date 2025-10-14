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

sub make_tcd
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x_chamber*cm $y_chamber*cm $z_chamber*cm";
 $detector{"rotation"}    = "$AngX_tcd*deg $AngY_tcd*deg 0*deg";
 $detector{"color"}       = "CCCC33";
 $detector{"type"}        = "Box";
 #$detector{"dimensions"}  = "35*cm 60*cm 250*cm";
 $detector{"dimensions"}  = "100*cm 100*cm 500*cm";
# $detector{"dimensions"}  = "70*cm 60*cm 300*cm";
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

make_tcd();

1;
