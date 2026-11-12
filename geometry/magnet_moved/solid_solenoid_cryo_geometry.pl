#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_solenoid_cryo';

my $DetectorMother="root";
my $CoilMother="$DetectorName\_CryostatVacuum";

sub solid_solenoid_cryo
{
make_cryostat();
make_coilCLEO();
make_coilShell();
make_coilSecondary();
}

# refer to  https://userweb.jlab.org/~gen/jlab12gev/cleo_mag/
# from CLEO-II NIM paper "the stainless cryostat" "consists of a 12 mm thick outer cylinder, a 10 mm inner cylinder and two 20 mm thick end flanges which are bolted and sealed with O-rings."
# stainless steel cryostat 3780 mm long, IR 1440 mm, OR 1750 mm; the shells run between the end flanges

sub make_cryostat
{
 my $Length = 378.0;  # full length, cm
 my $IR     = 144.0;
 my $OR     = 175.0;
 my $tIn    = 1.0;    # inner shell thickness
 my $tOut   = 1.2;    # outer shell thickness
 my $tSide  = 2.0;    # end flange thickness

 my $Nplate  = 4;
 my @PlateZ  = (0, 0, -($Length/2-$tSide/2), $Length/2-$tSide/2);
 my @Rin  = ($IR, $OR-$tOut, $IR, $IR);
 my @Rout = ($IR+$tIn, $OR, $OR, $OR);
 my @Dz   = ($Length/2-$tSide, $Length/2-$tSide, $tSide/2, $tSide/2);
 my @name = ("CryostatInner","CryostatOuter","CryostatFlangeUpstream","CryostatFlangeDownstream");
 my $material="G4_STAINLESS-STEEL";
 my $color="ffffff";

 for(my $n=1; $n<=$Nplate; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
    $detector{"material"}   = "$material";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }

 # the space inside is insulating vacuum (CLEO-II NIM: "SS VACUUM VESSEL", indirect cooling, LHe only in pipes on the coil shell);
 # it is the mother of the coil, its secondary and its support cylinder
 my %detector=init_det();
 $detector{"name"}        = "$CoilMother";
 $detector{"mother"}      = "$DetectorMother" ;
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}        = "0*cm 0*cm 0*cm";
 $detector{"rotation"}   = "0*deg 0*deg 0*deg";
 $detector{"color"}      = "$color";
 $detector{"type"}       = "Tube";
 $detector{"dimensions"} = ($IR+$tIn)."*cm ".($OR-$tOut)."*cm ".($Length/2-$tSide)."*cm 0*deg 360*deg";
 $detector{"material"}   = "G4_Galactic";
 $detector{"mfield"}     = "no";
 $detector{"ncopy"}      = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 0;
 $detector{"style"}       = 1;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}

# coil and its support cylinder (cold dimensions):
# "2 layers of conductor wound inside of a support cylinder", conductor 16x5 mm2 Al with 4x1.2 mm2 superconductor,
# "the center of the two layers is at R=1533mm", "the coil length is 3475mm", so 3.2 cm thick, z +-173.75 (all Al here);
# support cylinder "3535mm long, ID=3100mm (cold), 40mm thick", Al, assumed centered like the coil
# from CLEO-II NIM paper sec 2.3, "a single turn secondary in the form of a 1.5 mm thick, high purity (RRR > 2000) aluminum sheet
# is first fastened to the inside of the shell": r 154.85-155.0, length assumed equal to the coil;
# the coil sits just inside it, r 151.65-154.85 (center 153.25)

sub make_coilCLEO
{
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_coilCLEO";
    $detector{"mother"}      = "$CoilMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "ff8000";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 2*counts 151.65*cm 151.65*cm 154.85*cm 154.85*cm -173.75*cm 173.75*cm";
    $detector{"material"}   = "G4_Al";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_coilShell
{
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_coilShell";
    $detector{"mother"}      = "$CoilMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "aaaaaa";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "155.0*cm 159.0*cm 176.75*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Al";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_coilSecondary
{
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_coilSecondary";
    $detector{"mother"}      = "$CoilMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "cccccc";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "154.85*cm 155.0*cm 173.75*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Al";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}
