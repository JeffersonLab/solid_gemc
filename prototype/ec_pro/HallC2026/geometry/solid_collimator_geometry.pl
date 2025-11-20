#!/usr/bin/perl -w
# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_collimator';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";

my $z	= -140; # before 11ring

sub solid_collimator
{
make();
}

sub make
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CCCC00";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "7.5*cm 15*cm 5*cm 0*deg 360*deg";
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

solid_collimator();
1;
