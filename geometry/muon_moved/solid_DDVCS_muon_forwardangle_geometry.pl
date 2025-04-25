#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_forwardangle';

my $DetectorMother="root";

# 1 of 8 side of 3rd layer of CLEO is roughly 98.55"x14.21"x210" (250x36x534cm)
# echa side is made of 2 pieces of 62" and 37"
# we have 13 pieces: 7 at 62" x 14" x 210" and 6 at 37" x 14" x 210", while Cornell 3 for their use
# see drawings at
# https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/6052-303%20-Sh%2003.pdf
# https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/6052-303%20-Sh%2006.pdf
#  https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/
# 
#  photos at
#  https://solid.jlab.org/files/cleo_manual/cleo_scan1.pdf
#  https://solid.jlab.org/files/cleo_manual/cleo_scan2.pdf

 my $color="00ffff";
 my $material="G4_Fe";
 my $color_scint="ff0000"; 
 my $material_scint="G4_PLASTIC_SC_VINYLTOLUENE";
 
sub solid_DDVCS_muon_forwardangle
{
# make_solid_DDVCS_muon_forwardangle_w();
# make_solid_DDVCS_muon_forwardangle_0();
make_solid_DDVCS_muon_forwardangle_1();
make_solid_DDVCS_muon_forwardangle_2();
make_solid_DDVCS_muon_forwardangle_3();
}

sub make_solid_DDVCS_muon_forwardangle_w
{
    my $z=-200;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_w";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "bbffff";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "15*cm 55*cm 25*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_W";
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

sub make_solid_DDVCS_muon_forwardangle_0
{
    my $z=570+50/2.;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_0";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "bbffff";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "35*cm 285*cm 18*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Pb";
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

sub make_solid_DDVCS_muon_forwardangle_1
{
    my $z=570+50+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_11";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_12";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "-160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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

    my $z_virt_0=$z-18-1;    
    $detector{"name"}        = "$DetectorName\_0_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_virt_0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 0.01*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6100000";
    print_det(\%configuration, \%detector);
    
    my $z_virt=$z+18+1;    
    $detector{"name"}        = "$DetectorName\_1_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_virt*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 0.01*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6110000";
    print_det(\%configuration, \%detector);
    
    my $z_scint=$z+18+7;    
    $detector{"name"}        = "$DetectorName\_1_scint";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_scint*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color_scint;
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 2.5*cm 0*deg 360*deg";
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6101000"; 
    print_det(\%configuration, \%detector);   
}


sub make_solid_DDVCS_muon_forwardangle_2
{
    my $z=570+50+18+18+10+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_21";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 160*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "267*cm 125*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_22";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm -160*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "267*cm 125*cm 18*cm";
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
    
    my $z_virt=$z+18+1;    
    $detector{"name"}        = "$DetectorName\_2_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_virt*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 0.01*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6120000";
    print_det(\%configuration, \%detector);
    
    my $z_scint=$z+18+7;    
    $detector{"name"}        = "$DetectorName\_2_scint";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_scint*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color_scint;
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 2.5*cm 0*deg 360*deg";
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6102000";
    print_det(\%configuration, \%detector);      
}


sub make_solid_DDVCS_muon_forwardangle_3
{
    my $z=570+50+18+18+10+18+18+10+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_31";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_32";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "-160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
    
    my $z_virt=$z+18+1;    
    $detector{"name"}        = "$DetectorName\_3_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_virt*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 0.01*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6130000";
    print_det(\%configuration, \%detector);
    
    my $z_scint=$z+18+7;    
    $detector{"name"}        = "$DetectorName\_3_scint";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_scint*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color_scint;
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 2.5*cm 0*deg 360*deg";
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6103000";
    print_det(\%configuration, \%detector);     
    
    my $z_virt_4=$z+18+10;    
    $detector{"name"}        = "$DetectorName\_4_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_virt_4*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "100*cm 285*cm 0.01*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6140000";
    print_det(\%configuration, \%detector);    
}
