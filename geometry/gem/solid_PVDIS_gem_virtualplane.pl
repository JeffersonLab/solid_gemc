use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_gem_virtualplane';

my $DetectorMother="root";

sub solid_PVDIS_gem_virtualplane
{
make_gem_virtualplane();
make_gem_shield();
}

my $Nplate	= $parameters{"Nplate"};
my $PlateZ1	= $parameters{"PlateZ1"};
my $PlateZ2	= $parameters{"PlateZ2"};
my $PlateZ3	= $parameters{"PlateZ3"};
my $PlateZ4	= $parameters{"PlateZ4"};
my $PlateZ5	= $parameters{"PlateZ5"};
my $Rin1	= $parameters{"Rin1"};
my $Rin2	= $parameters{"Rin2"};
my $Rin3	= $parameters{"Rin3"};
my $Rin4	= $parameters{"Rin4"};
my $Rin5	= $parameters{"Rin5"};
my $Rout1	= $parameters{"Rout1"};
my $Rout2	= $parameters{"Rout2"};
my $Rout3	= $parameters{"Rout3"};
my $Rout4	= $parameters{"Rout4"};
my $Rout5	= $parameters{"Rout5"};

 my @PlateZ = ($PlateZ1-0.9,$PlateZ2-0.9,$PlateZ3-0.9,$PlateZ4-0.9,$PlateZ5-0.9);
 my @Rin    = ($Rin1,$Rin2,$Rin3,$Rin4,$Rin5);
 my @Rout   = ($Rout1,$Rout2,$Rout3,$Rout4,$Rout5);
 
sub make_gem_virtualplane
{
#  my $Nplate  = 4;
# my @PlateZ  = (155,185,295,310);  # as on p56 of pac34 proposal
# my @PlateZ  = (157.5,185.5,297,306); # as in Eugen's code
#  my @PlateZ  = (157.5,185.5,306,321);  # change for last two planes further back as Cherenkov needs 10cm more
#  my @Rin  = (55,65,105,115);
#  my @Rout = (115,140,200,215);
#  my @PlateZ = (157.5,185.5,306,315);
#  my $Nplate  = 5;
#  my @PlateZ = (156.8,184.8,188.8,305,314);

#cover target center at 10cm from 21 to 36 degree
#  my @Rin = (56,67,113,117);
#  my @Rout = (108,129,215,222);
#  my $Dz   = 0.5;
#cover 40cm long full target with center at 10cm from 21 to 36 degree
#  my @Rin = (48,59,105,109);
#  my @Rout = (122,143,230,237);
#  my @Rin = (45,55,63,90,90);
#  my @Rout = (144,144,144,270,270);
 my $Dz   = 0.001/2; 
#  my $color="44ee11";

 for(my $n=1; $n<=$Nplate; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    my $id=1000000+$n*100000+10000;
    $detector{"identifiers"} = "id manual $id";
    print_det(\%configuration, \%detector);
 }
}

sub make_gem_shield
{
 my $Dz   = 0.01436/2; 
 my @PlateZ_shield = ($PlateZ1-0.95);

 for(my $n=1; $n<=2; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$n";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $PlateZ_shield[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "44ee11";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Gu";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}
