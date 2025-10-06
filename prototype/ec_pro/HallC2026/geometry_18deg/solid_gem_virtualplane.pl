#!/usr/bin/perl -w
# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_gem';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";

#my $z1	= -56*2.54/2-30.48-8.26-7.03-9.3-1.5/2-1;
my $z1	= -56*2.54/2-8.26-7.03-9.-1.5/2-1;
my $z2	=  56*2.54/2+4.5+1.5/2-1;
#my $z3	= -56*2.54/2-30.48-8.26-7.03-1.5/2-1;
my $z3	= -56*2.54/2-8.26-7.03-1.5/2-1;
my $z4	=  56*2.54/2+4.5+8.5+1.5/2-1;
my $hx	= 5.12;
my $hy	= 5.12;

sub solid_gem_virtualplane
{
make_1();
make_2();
make_3();
make_4();
}

sub make_1
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_virtualplane_1";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z1*cm";
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
 my $ID = 5;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}

sub make_2
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_virtualplane_2";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z2*cm";
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
 my $ID = 6;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_3
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_virtualplane_3";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z3*cm";
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
 my $ID = 20;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_4
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_virtualplane_4";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z4*cm";
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
 my $ID = 21;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}

solid_gem_virtualplane();
1;
