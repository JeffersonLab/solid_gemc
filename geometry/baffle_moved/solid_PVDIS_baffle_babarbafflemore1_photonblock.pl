#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_baffle_babarbafflemore1_photonblock';

my $DetectorMother="root";

sub solid_PVDIS_baffle_babarbafflemore1_photonblock
{
make_ec_forwardangle_block();
}

my $color_baffle="00C0C0";

my $material="G4_Pb";#"Lead";

sub make_ec_forwardangle_block
{
 for(my $n=1; $n<=30; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 330*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_baffle";
    $detector{"type"}       = "Tube";
    my $phi_s=($n-1)*12+2.2;
    $detector{"dimensions"} = "110*cm 200*cm 2.5*cm $phi_s*deg 2.5*deg";
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
}
