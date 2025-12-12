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

my $x1	= 0;
my $x2  = 0;
my $x3  = 0;
my $x4  = -9.1;
my $x5  = 9.1;
my $y4  = 0;
my $y5  = 0;
my $z1	= -98;
my $z3  = -167;
my $z13  = -167-2.2;
my $z2  = -131;
my $z4  = -116.5;
my $z5  = -116.5;
my $z6  = -267;
my $z16  = -267-2.2;
my $z7  = -114.5;
my $hx	= 50;
my $hy	= 50;
my $hx4	= 15.5;
my $hy4	= 7.15;
my $hx5	= 15.5;
my $hy5	= 7.15;
my $Rmin= 7.1;
my $Rmax= 7.1+0.0001;
my $hx13= 30;
my $hy13= 30;
my $hx16= 30;
my $hy16= 30;
my $hx24= 15.5;
my $hy24= 9.15;
my $hx25= 15.5;
my $hy25= 9.15;

sub solid_magnet_virtualplane
{
make_1();
#make_2();
make_3();
make_colli1();
#make_4();
#make_5();
#make_6();
#make_colli2();
#make_7();
make_Tunnel_left();
make_Tunnel_right();
make_Tunnel_top();
make_Tunnel_bottom();
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
sub make_3
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_3";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x3*cm 0*cm $z3*cm";
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
 my $ID = 53;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_4
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_4";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x4*cm $y4*cm $z4*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx4*cm $hy4*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 54;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_5
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_5";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x5*cm $y5*cm $z5*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx5*cm $hy5*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 55;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_6
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_6";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x3*cm 0*cm $z6*cm";
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
 my $ID = 56;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_7
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_7";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z7*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Tube";
 $detector{"dimensions"}  = "$Rmin*cm $Rmax*cm 16.4*cm 0*deg 360*deg"; 
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 54;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_colli1{
    # Outer box: 30x30x4 cm³ Pb 
    my $outer_name = "${DetectorName}_collimator_box";
    my %outer = init_det();
    $outer{"name"}        = $outer_name;
    $outer{"mother"}      = $DetectorMother;
    $outer{"description"} = "60x60 cm pb collimator outer box";
    $outer{"pos"}         = "0*cm 0*cm $z13*cm";
    $outer{"rotation"}    = "0*deg 0*deg 0*deg";
    $outer{"color"}       = "cc6633";
    $outer{"type"}        = "Box";
    $outer{"dimensions"}  = "$hx13*cm $hy13*cm 2*cm";
    $outer{"material"}    = "G4_Pb";
    $outer{"ncopy"}       = 1;
    $outer{"pMany"}       = 1;
    $outer{"exist"}       = 1;
    $outer{"visible"}     = 0;
    print_det(\%configuration, \%outer);

    # Inner box (hole): 10x10x0.1 cm³ Air
    my $inner_name = "${DetectorName}_collimator_hole";
    my %inner = init_det();
    $inner{"name"}        = $inner_name;
    $inner{"mother"}      = $DetectorMother;
    $inner{"description"} = "10x10 cm air hole inside collimator";
    $inner{"pos"}         = "0*cm 0*cm $z13*cm";
    $inner{"rotation"}    = "0*deg 0*deg 0*deg";
    $inner{"color"}       = "00ccff";  # light blue for visibility
    $inner{"type"}        = "Box";
    $inner{"dimensions"}  = "5*cm 5*cm 2*cm";
    $inner{"material"}    = "G4_AIR";
    $inner{"ncopy"}       = 1;
    $inner{"pMany"}       = 1;
    $inner{"exist"}       = 1;
    $inner{"style"}       = 1;
    $inner{"visible"}     = 1;
    $inner{"mfield"}      = "no";
    $inner{"sensitivity"} = "no";
    $inner{"hit_type"}    = "no";
    $inner{"identifiers"} = "no";
    print_det(\%configuration, \%inner);

    # Boolean operation: subtract hole from the box
    my %collimator = init_det();
    $collimator{"name"}        = "${DetectorName}_collimator_final";
    $collimator{"mother"}      = $DetectorMother;
    $collimator{"description"} = "Final collimator with 10x10 cm air hole";
    $collimator{"pos"}         = "0*cm 0*cm $z13*cm";
    $collimator{"rotation"}    = "0*deg 0*deg 0*deg";
    $collimator{"color"}       = "cc6633";
    $collimator{"type"}        = "Operation: $outer_name - $inner_name";
    $collimator{"material"}    = "G4_Pb";
    $collimator{"ncopy"}       = 1;
    $collimator{"pMany"}       = 1;
    $collimator{"exist"}       = 1;
    $collimator{"style"}       = 1;
    $collimator{"visible"}     = 1;
    $collimator{"mfield"}      = "no";
    $collimator{"sensitivity"} = "no";
    $collimator{"hit_type"}    = "no";
    $collimator{"identifiers"} = "no";
    print_det(\%configuration, \%collimator);
}
sub make_colli2 {
    # Outer box: 30x30x4 cm³ Pb
    my $outer_name = "${DetectorName}_collimator_box_2";
    my %outer = init_det();
    $outer{"name"}        = $outer_name;
    $outer{"mother"}      = $DetectorMother;
    $outer{"description"} = "60x60 cm Pb collimator outer box 2";
    $outer{"pos"}         = "0*cm 0*cm $z16*cm";
    $outer{"rotation"}    = "0*deg 0*deg 0*deg";
    $outer{"color"}       = "cc6633";
    $outer{"type"}        = "Box";
    $outer{"dimensions"}  = "$hx16*cm $hx16*cm 2*cm";
    $outer{"material"}    = "G4_Pb";
    $outer{"ncopy"}       = 1;
    $outer{"pMany"}       = 1;
    $outer{"exist"}       = 1;
    $outer{"visible"}     = 0;
    print_det(\%configuration, \%outer);

    # Inner box (hole): 10x10x0.1 cm³ Air
    my $inner_name = "${DetectorName}_collimator_hole_2";
    my %inner = init_det();
    $inner{"name"}        = $inner_name;
    $inner{"mother"}      = $DetectorMother;
    $inner{"description"} = "10x10 cm air hole inside collimator 2";
    $inner{"pos"}         = "0*cm 0*cm $z16*cm";
    $inner{"rotation"}    = "0*deg 0*deg 0*deg";
    $inner{"color"}       = "00ccff";  # light blue for visibility
    $inner{"type"}        = "Box";
    $inner{"dimensions"}  = "5*cm 5*cm 2*cm";
    $inner{"material"}    = "G4_AIR";
    $inner{"ncopy"}       = 1;
    $inner{"pMany"}       = 1;
    $inner{"exist"}       = 1;
    $inner{"style"}       = 1;
    $inner{"visible"}     = 1;
    $inner{"mfield"}      = "no";
    $inner{"sensitivity"} = "no";
    $inner{"hit_type"}    = "no";
    $inner{"identifiers"} = "no";
    print_det(\%configuration, \%inner);

    # Boolean operation: subtract hole from the box
    my %collimator = init_det();
    $collimator{"name"}        = "${DetectorName}_collimator_final_2";
    $collimator{"mother"}      = $DetectorMother;
    $collimator{"description"} = "Final collimator with 10x10 cm air hole";
    $collimator{"pos"}         = "0*cm 0*cm $z16*cm";
    $collimator{"rotation"}    = "0*deg 0*deg 0*deg";
    $collimator{"color"}       = "cc6633";
    $collimator{"type"}        = "Operation: $outer_name - $inner_name";
    $collimator{"material"}    = "G4_Pb";
    $collimator{"ncopy"}       = 1;
    $collimator{"pMany"}       = 1;
    $collimator{"exist"}       = 1;
    $collimator{"style"}       = 1;
    $collimator{"visible"}     = 1;
    $collimator{"mfield"}      = "no";
    $collimator{"sensitivity"} = "no";
    $collimator{"hit_type"}    = "no";
    $collimator{"identifiers"} = "no";
   print_det(\%configuration, \%collimator);
}
sub make_Tunnel_right
{
 my %detector=init_det();
 $detector{"name"}        = "PbTunnel_right";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x4*cm $y4*cm $z4*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx4*cm $hy4*cm 2*cm";	    
 $detector{"material"}    = "G4_Pb";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
 $detector{"mfield"}      = "no";
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}
sub make_Tunnel_left 
{
 my %detector=init_det();
 $detector{"name"}        = "PbTunnel_left";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x5*cm $y5*cm $z5*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx5*cm $hy5*cm 2*cm";	    
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
sub make_Tunnel_top
{
 my %detector=init_det();
 $detector{"name"}        = "PbTunnel_top";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 9.1*cm $z4*cm";
 $detector{"rotation"}    = "90*deg 0*deg 0*deg";
 $detector{"color"}       = "999999";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hy24*cm $hx24*cm 2*cm";	    
 $detector{"material"}    = "G4_Pb";
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"mfield"}      = "no";
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}
sub make_Tunnel_bottom 
{
 my %detector=init_det();
 $detector{"name"}        = "PbTunnel_bottom";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm -9.1*cm $z5*cm";
 $detector{"rotation"}    = "90*deg 0*deg 0*deg";
 $detector{"color"}       = "999999";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hy25*cm $hx25*cm 2*cm";	    
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
solid_magnet_virtualplane();
1;
