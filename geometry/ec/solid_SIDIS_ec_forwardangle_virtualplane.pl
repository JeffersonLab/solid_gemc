#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_SIDIS_ec_forwardangle_virtualplane';

my $DetectorMother="root";

sub solid_SIDIS_ec_forwardangle_virtualplane
{
make_solid_SIDIS_ec_forwardangle_virtualplane_front();
# make_solid_SIDIS_ec_forwardangle_virtualplane_middle();
make_solid_SIDIS_ec_forwardangle_virtualplane_inner();
make_solid_SIDIS_ec_forwardangle_virtualplane_rear();
}

sub make_solid_SIDIS_ec_forwardangle_virtualplane_front
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_front";
 $detector{"mother"}      = $DetectorMother;
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 413*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 90;
  my $Rmax1 = 265;
  my $Rmin2 = 90;
  my $Rmax2 = 265;
  my $Dz    = 0.001/2;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 3110000";
 print_det(\%configuration, \%detector);
}

sub make_solid_SIDIS_ec_forwardangle_virtualplane_middle
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_middle";
 $detector{"mother"}      = $DetectorMother;
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 414.6*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 90;
  my $Rmax1 = 265;
  my $Rmin2 = 90;
  my $Rmax2 = 265;
  my $Dz    = 0.001/2;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 3120000";
 print_det(\%configuration, \%detector);
}

sub make_solid_SIDIS_ec_forwardangle_virtualplane_inner
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_inner";
 $detector{"mother"}      = $DetectorMother;
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 440*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 96;
  my $Rmax1 = 96.1;
  my $Rmin2 = 96;
  my $Rmax2 = 96.1;
  my $Dz    = 25;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 3130000";
 print_det(\%configuration, \%detector);
}

sub make_solid_SIDIS_ec_forwardangle_virtualplane_rear
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_rear";
 $detector{"mother"}      = $DetectorMother;
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 466*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 90;
  my $Rmax1 = 265;
  my $Rmin2 = 90;
  my $Rmax2 = 265;
  my $Dz    = 0.001/2;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 3140000";
 print_det(\%configuration, \%detector);
}

