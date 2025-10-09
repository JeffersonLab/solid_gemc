# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_spd';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";
my $z_shower = 56*2.54/2+4.5+1.5+14.7+5+2+0.6+1.5+2+15.1+1.1+2.4+2+45/2.; 

my $sc1_x_hlfln = 5./2;  
my $sc1_y_hlfln = 7.5/2;  
my $sc2_x_hlfln = 5.5;
my $sc2_y_hlfln = 6.35;
my $sc3_x_hlfln = 5./2;  
my $sc3_y_hlfln = 10/2;  
my $sc4_y1_hlfln = 3.5/2;  
my $sc4_y2_hlfln = 4.5/2;  
my $sc4_x_hlfln = 18./2;  
my $sc_z_hlfln = 1./2;
my $sc2_z_hlfln = 2./2;
my $Z_sc1 = -56*2.54/2-8.26-1./2;
my $Z_sc2 =  56*2.54/2+4.5+8.5+6.12+5.62+6.89+2./2;  

my $spd1_x_hlfln = 10./2;  
my $spd1_y_hlfln = 25./2;  
my $spd2_x_hlfln = 14./2;
my $spd2_y_hlfln = 40./2;
my $spd1_z_hlfln = 0.6/2;
my $spd2_z_hlfln = 2./2;

my $Z_spd1 = 56*2.54/2+4.5+1.5+14.7+5.+2+0.6/2;
my $Z_spd2 = 56*2.54/2+4.5+8.5+6.12+5.62+2./2;  

my $material_scint = "SL_spd_ScintillatorB";

sub spd_geometry
{
make_sc1();
make_sc2();
make_sc3();
make_sc4();
#make_spd1();
make_spd2();
}

sub make_sc1
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_sc1";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $Z_sc1*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$sc1_x_hlfln*cm $sc1_y_hlfln*cm $sc_z_hlfln*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 1";
	print_det(\%configuration, \%detector);
 
}

sub make_sc2
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_sc2";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $Z_sc2*cm"; 
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$sc2_x_hlfln*cm $sc2_y_hlfln*cm $sc2_z_hlfln*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 2";
	print_det(\%configuration, \%detector);
}
sub make_sc3
{
 my $z_sc3 = -56*2.54/2+30.48+243.84+22.86+1/2.;
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_sc3";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $z_sc3*cm"; 
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$sc3_x_hlfln*cm $sc3_y_hlfln*cm $sc_z_hlfln*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 11";
	print_det(\%configuration, \%detector);
}
sub make_sc4
{
 my $z_sc4 = 56*2.54/2+4.5+8.5+6.12+2./2;
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_sc4";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $z_sc4*cm"; 
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$sc4_x_hlfln*cm $sc4_y2_hlfln*cm $sc2_z_hlfln*cm";
	#$detector{"dimensions"}  = "1.8*cm 1.5*cm 32*cm 32*cm 10*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 12";
	print_det(\%configuration, \%detector);
}


sub make_spd1
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_spd1";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $Z_spd1*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$spd1_x_hlfln*cm $spd1_y_hlfln*cm $spd1_z_hlfln*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 3";
	print_det(\%configuration, \%detector);
 
}

sub make_spd2
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_spd2";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm $Z_spd2*cm"; 
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "CC6633";
	$detector{"type"}        = "Box";
	$detector{"dimensions"}  = "$spd2_x_hlfln*cm $spd2_y_hlfln*cm $spd2_z_hlfln*cm";
	$detector{"material"}    = $material_scint;
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "solid_spd";
	$detector{"hit_type"}    = "solid_spd";
	$detector{"identifiers"} = "id manual 4";
	print_det(\%configuration, \%detector);
}

spd_geometry();

1;
