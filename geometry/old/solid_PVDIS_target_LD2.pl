#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_target_LD2';

my $DetectorMother="root";

sub solid_PVDIS_target_LD2
{
make_target_PVDIS_target();
}

# target offset in cm
my $targetoff=10;

# C --     Target  ==================================
# C
# GPARVOL21  'TACH'  209  'TARC'    0.    0.    0.    0  'TUBE'  3   0.    5.00   25.00    
# GPARVOL22  'TACV'  221  'TACH'    0.    0.    0.    0  'TUBE'  3   0.    4.95   24.95    
# GPARVOL23  'TAW1'  221  'TACH'    0.    0.  -24.975 0  'TUBE'  3   0.    1.22    0.025   
# GPARVOL24  'TAW2'  221  'TACH'    0.    0.   24.975 0  'TUBE'  3   0.    1.50    0.025   
# GPARVOL25  'TALU'  209  'TACV'    0.    0.    0.    0  'TUBE'  3   0.    1.918  20.0     
# GPARVOL26  'TLH2'  201  'TALU'    0.    0.    0.    0  'TUBE'  3   0.    1.900  19.982   

#40cm LH2/LD2 target

sub make_target_PVDIS_target
{
# 120um windows, 15mil wall, 4cm radius?
 my $NUM  = 6;
 my @z    = (0.+$targetoff,0.,-39.98,39.98,0.,0.);
 my @Rin  = (0.,0.,0.,0.,0.,0.);
 my @Rout = (3,2.9,2.9,2.9,2.538,2.5);
 my @Dz   = (40,39.96,0.02,0.02,20.012,20.);
 my @name = ("$DetectorName\_TACH","$DetectorName\_TACV","$DetectorName\_TAW1","$DetectorName\_TAW2","$DetectorName\_TALU","$DetectorName\_TAH2"); 
 my @mother=("$DetectorMother","$DetectorName\_TACH","$DetectorName\_TACH","$DetectorName\_TACH","$DetectorName\_TACV","$DetectorName\_TALU");
 #SL_Vacuum is the vacuum with certain air, or use G4_Galactic for pure vacuum 
 my @mat  = ("SL_Vacuum","SL_Vacuum","SL_Vacuum","SL_Vacuum","G4_Al","LD2");
 my @color = ("0000ff","808080","0000ff","0000ff","0000ff","ff0000");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
    $detector{"material"}   = $mat[$n-1];
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
