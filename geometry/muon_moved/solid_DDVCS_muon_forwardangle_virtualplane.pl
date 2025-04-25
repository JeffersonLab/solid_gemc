#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_forwardangle_virtualplane';

my $DetectorMother="root";

sub solid_DDVCS_muon_forwardangle_virtualplane
{
make_solid_DDVCS_muon_forwardangle_virtualplane_1();
make_solid_DDVCS_muon_forwardangle_virtualplane_2();
make_solid_DDVCS_muon_forwardangle_virtualplane_3();
make_solid_DDVCS_muon_forwardangle_virtualplane_4();
make_solid_DDVCS_muon_forwardangle_virtualplane_5();
}

sub make_solid_DDVCS_muon_forwardangle_virtualplane_1
{
 my $z=570+50-5;    
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_1";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 100;
  my $Rmax1 = 285;
  my $Rmin2 = 100;
  my $Rmax2 = 285;
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
 $detector{"identifiers"} = "id manual 6110000";
 print_det(\%configuration, \%detector);
}

sub make_solid_DDVCS_muon_forwardangle_virtualplane_2
{
 my $z=570+50+18+18+5;    
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_2";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 100;
  my $Rmax1 = 285;
  my $Rmin2 = 100;
  my $Rmax2 = 285;
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
 $detector{"identifiers"} = "id manual 6120000";
 print_det(\%configuration, \%detector);
}

sub make_solid_DDVCS_muon_forwardangle_virtualplane_3
{
 my $z=570+50+18+18+10+18+18+5;    
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_3";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 100;
  my $Rmax1 = 285;
  my $Rmin2 = 100;
  my $Rmax2 = 285;
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
 $detector{"identifiers"} = "id manual 6130000";
 print_det(\%configuration, \%detector);
}

sub make_solid_DDVCS_muon_forwardangle_virtualplane_4
{
 my $z=570+50+18+18+10+18+18+10+18+18+5;    
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_4";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 100;
  my $Rmax1 = 285;
  my $Rmin2 = 100;
  my $Rmax2 = 285;
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
 $detector{"identifiers"} = "id manual 6140000";
 print_det(\%configuration, \%detector);
}

sub make_solid_DDVCS_muon_forwardangle_virtualplane_5
{
 my $z=570+50+18+18+10+18+18+10+18+18+5+10;    
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_5";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 100;
  my $Rmax1 = 285;
  my $Rmin2 = 100;
  my $Rmax2 = 285;
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
 $detector{"identifiers"} = "id manual 6150000";
 print_det(\%configuration, \%detector);
}
