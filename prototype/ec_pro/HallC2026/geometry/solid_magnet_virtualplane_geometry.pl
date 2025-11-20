#!/usr/bin/perl -w
# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_magnet_virtualplane';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";

my $x1	= -35;
my $x2  = -35;
my $z1	= -100;
my $z2  = -200;
my $hx	= 50;
my $hy	= 50;

sub solid_magnet_virtualplane
{
make_1();
make_2();
}

sub make_1
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_1";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x1*cm 0*cm $z1*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 51;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}

sub make_2
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_2";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x2*cm 0*cm $z2*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 52;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}

solid_magnet_virtualplane();
1;
