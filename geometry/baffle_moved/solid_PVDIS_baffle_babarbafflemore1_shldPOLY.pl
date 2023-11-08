#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_baffle_babarbafflemore1_shldPOLY';

my $DetectorMother="root";

sub solid_PVDIS_baffle_babarbafflemore1_shldPOLY
{
make_baffle_plate_inner_absorber_shldPOLY();
make_baffle_plate_outer_absorber_shldPOLY();
make_baffle_plate_inner_shldPOLY();
make_baffle_plate_outer_shldPOLY();
make_baffle_plate_shldPOLY();
make_baffle_blocks_shldPOLY();
}

#6 plates, 9cm thick,every 30cm, 30 slits per plate, 20 blocks per slit

my $color_baffle="00C0C0";
my $color_baffleout="CCFFFF";
#my $color_abs="00C0C0";
my $color_abs="9999FF"; # in case of poly

my $material_baffle_within="G4_AIR";
my $material_baffle="SL_PVDIS_baffle_BorPoly";
my $sensitivity_baffle="no";
my $hit_baffle="no";

# In PVDIS PAC34 proposal, page 35 (http://hallaweb.jlab.org/collab/PAC/PAC34/PR-09-012-pvdis.pdf). It's mentioned that the baffle position is 30,60,90,120,150,180cm relative to BaBar solenoid coil center at 0. But it's actually relative to the target center which is at 10cm as seen on the Fig 3.3 on page 31.
# Later on, it seems Eugene has further tweaked them to 30,58,86,114,142,170. refer to http://hallaweb.jlab.org/12GeV/SoLID/download/sim/geant3/solid_comgeant_pvdis/pvdis_02_01_p_14_01/fort.22

my $targetoff=10; # target offset in cm

 my $Nplate  = 11;  # default 6;
 my @PlateZ  = (30+$targetoff,44+$targetoff,58+$targetoff,72+$targetoff,86+$targetoff,100+$targetoff,114+$targetoff,128+$targetoff,142+$targetoff,156+$targetoff,170+$targetoff); # (30+$targetoff,58+$targetoff,86+$targetoff,114+$targetoff,142+$targetoff,170+$targetoff);
 my $Dz   = 9.0/2.;
my $Dz_shld = 5.0; #full length (Polycone)

# gem
#  my @PlateZ = (157.5,185.5,190,306,315);
#  my @Rin = (48,59,65,105,109);
#  my @Rout = (122,143,143,230,237); 

sub make_baffle_plate_inner_absorber_shldPOLY
{
 #my @Rin  = (3.89,14.,19.,23.9,28.9,33.8); #for 6 baffles geometry
 #my @Rout = (3.9,15.3,26.6,37.9,49.2,60.4); #for 6 baffles geometry
 #my @Rin  = (3.89,8.95,14.00,16.50,19.00,21.45,23.90,26.40,28.90,31.35,33.80);
#  my @Rin  = (3.89,8.95,14.00,18.50,24.00,30.45,35.90,41.40,47.90,52.35,58.80);
#  my @Rout = (3.90,9.62,15.30,20.92,26.60,32.22,37.90,43.52,49.20,54.82,60.40); 
 my @Rout  = (8.95-0.1,8.95-0.1,14.00-0.1,14.00-0.1,18.50-0.1,18.50-0.1,24.00-0.1,24.00-0.1,30.45-0.1,30.45-0.1,35.90-0.1,35.90-0.1,41.40-0.1,41.40-0.1,47.90-0.1,47.90-0.1,52.35-0.1,52.35-0.1,58.80-0.1,58.80-0.1);
 my @Rin = (5.56,7.56,8.67,10.67,11.80,13.80,14.90,16.90,18.03,20.03,20.03,20.03,20.03,20.03,20.03,20.03,20.03,20.03,20.03,20.03);
 my @Zabs= ($PlateZ[1]-4.5,$PlateZ[1]+4.5,$PlateZ[2]-4.5,$PlateZ[2]+4.5,$PlateZ[3]-4.5,$PlateZ[3]+4.5,$PlateZ[4]-4.5,$PlateZ[4]+4.5,$PlateZ[5]-4.5,$PlateZ[5]+4.5,$PlateZ[6]-4.5,$PlateZ[6]+4.5,$PlateZ[7]-4.5,$PlateZ[7]+4.5,$PlateZ[8]-4.5,158.00,$PlateZ[9]-4.5,172.00,$PlateZ[10]-4.5,188.99);
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_plateinner_absorber";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 20*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rin[15]*cm $Rin[16]*cm $Rin[17]*cm $Rin[18]*cm $Rin[19]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $Rout[15]*cm $Rout[16]*cm $Rout[17]*cm $Rout[18]*cm $Rout[19]*cm $Zabs[0]*cm $Zabs[1]*cm $Zabs[2]*cm $Zabs[3]*cm $Zabs[4]*cm $Zabs[5]*cm $Zabs[6]*cm $Zabs[7]*cm $Zabs[8]*cm $Zabs[9]*cm $Zabs[10]*cm $Zabs[11]*cm $Zabs[12]*cm $Zabs[13]*cm $Zabs[14]*cm $Zabs[15]*cm $Zabs[16]*cm $Zabs[17]*cm $Zabs[18]*cm $Zabs[19]*cm";
    $detector{"material"}   = "$material_baffle";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_baffle";
    $detector{"hit_type"}    = "$hit_baffle";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_baffle_plate_inner_shldPOLY
{
 #my @Rin  = (3.89,14.,19.,23.9,28.9,33.8); #for 6 baffles geometry
 #my @Rout = (3.9,15.3,26.6,37.9,49.2,60.4); #for 6 baffles geometry
#  my @Rin  = (3.89,8.95,14.00,16.50,19.00,21.45,23.90,26.40,28.90,31.35,33.80);
 my @Rin  = (3.89,8.95,14.00,18.50,24.00,30.45,35.90,41.40,47.90,52.35,58.80);
 my @Rout = (3.90,9.62,15.30,20.92,26.60,32.22,37.90,43.52,49.20,54.82,60.40);

 for(my $n=2; $n<=($Nplate-3); $n++)
 {
    my $n_c     = cnumber($n-1, 1);
    my $z_c = $PlateZ[$n-1]+4.5;
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_plateinner$n_c";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_c*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 2*counts $Rin[$n-1]*cm $Rin[$n]*cm $Rout[$n-1]*cm $Rout[$n]*cm 0.0*cm $Dz_shld*cm";
    $detector{"material"}   = "$material_baffle";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_baffle";
    $detector{"hit_type"}    = "$hit_baffle";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}

sub make_baffle_plate_outer_absorber_shldPOLY
{
 #my @Rin  = (3.89,14.,19.,23.9,28.9,33.8); #for 6 baffles geometry
 #my @Rout = (3.9,15.3,26.6,37.9,49.2,60.4); #for 6 baffles geometry
 #my @Rin  = (3.89,8.95,14.00,16.50,19.00,21.45,23.90,26.40,28.90,31.35,33.80)
#  my @Rin = (52,52,62,62,72,72,82,82,92,92,102,102,112,112,120,120,132,132,140,140);
 my @Rin = (52,52,62,62,72,72,82,82,92,92,102,102,112,112,122,122,132,132,140,140);
 my @Rout = (140,140,140,140,140,140,140,140,140,140,140,140,140,140,140,140,140,140,140,140);
 my @Zabs= ($PlateZ[1]-4.5,$PlateZ[1]+4.5,$PlateZ[2]-4.5,$PlateZ[2]+4.5,$PlateZ[3]-4.5,$PlateZ[3]+4.5,$PlateZ[4]-4.5,$PlateZ[4]+4.5,$PlateZ[5]-4.5,$PlateZ[5]+4.5,$PlateZ[6]-4.5,$PlateZ[6]+4.5,$PlateZ[7]-4.5,$PlateZ[7]+4.5,$PlateZ[8]-4.5,158.00,$PlateZ[9]-4.5,172.00,$PlateZ[10]-4.5,188.99);
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_plateouter_absorber";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 20*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rin[15]*cm $Rin[16]*cm $Rin[17]*cm $Rin[18]*cm $Rin[19]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $Rout[15]*cm $Rout[16]*cm $Rout[17]*cm $Rout[18]*cm $Rout[19]*cm $Zabs[0]*cm $Zabs[1]*cm $Zabs[2]*cm $Zabs[3]*cm $Zabs[4]*cm $Zabs[5]*cm $Zabs[6]*cm $Zabs[7]*cm $Zabs[8]*cm $Zabs[9]*cm $Zabs[10]*cm $Zabs[11]*cm $Zabs[12]*cm $Zabs[13]*cm $Zabs[14]*cm $Zabs[15]*cm $Zabs[16]*cm $Zabs[17]*cm $Zabs[18]*cm $Zabs[19]*cm";
    $detector{"material"}   = "$material_baffle";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_baffle";
    $detector{"hit_type"}    = "$hit_baffle";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_baffle_plate_outer_shldPOLY
{
 #my @Rin  = (34.8,54.3,73.9,93.4,112.8,132.1); #for 6 baffles geometry
 #my @Rout = (140,140,140,140,140,140); #coil edge is at 142cm #for 6 baffles geometry

 my @Rin  = (34.70,44.54,54.30,64.01,73.80,83.48,93.30,102.95,112.80,122.43,132.00);
# my @Rout = (140,140,140,140,140,140,140,140,140,140,140); #coil edge is at 142cm
my @Rout = (42,52,62,72,82,92,102,112,120,132,140);
#  my @Rout = (42,143.9,143.9,143.9,143.9,143.9,143.9,143.9,143.9,143.9,143.9);

 for(my $n=2; $n<=($Nplate-3); $n++)
 {
    my $n_c     = cnumber($n-1, 1);
    my $z_c = $PlateZ[$n-1]+4.5;
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_plateouter$n_c";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_c*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 2*counts $Rin[$n-1]*cm $Rin[$n]*cm $Rout[$n-1]*cm $Rout[$n]*cm 0.0*cm $Dz_shld*cm";
    $detector{"material"}   = "$material_baffle";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_baffle";
    $detector{"hit_type"}    = "$hit_baffle";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}

sub make_baffle_plate_shldPOLY   #viutral container for baffle plate
{
 #my @Rin = (3.9,15.3,26.6,37.9,49.2,60.4); #for 6 baffles geometry
 #my @Rout  = (34.7,54.3,73.8,93.3,112.8,132.0); #for 6 baffles geometry

  my @Rin = (3.90,9.62,15.30,20.92,26.60,32.22,37.90,43.52,49.20,54.82,60.40);
  my @Rout  = (34.70,44.54,54.30,64.01,73.80,83.48,93.30,102.95,112.80,122.43,132.00);


 for(my $n=2; $n<=($Nplate-3); $n++)
 {
    my $n_c     = cnumber($n-1, 1);
    my $z_c = $PlateZ[$n-1]+4.5;
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_plate$n_c";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_c*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "000000";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 2*counts $Rin[$n-1]*cm $Rin[$n]*cm $Rout[$n-1]*cm $Rout[$n]*cm 0.0*cm $Dz_shld*cm";
    $detector{"material"}   = "$material_baffle_within";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_baffle";
    $detector{"hit_type"}    = "$hit_baffle";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}

#making buffle blocks, palte by plate, slit by slit,block by block
sub make_baffle_blocks_shldPOLY
{
my @x =( #Rin     Rout   SPhi    DPhi
#1 plate
[	 3.90 ,   5.44 ,  176.26 , 183.00 ,
	 5.44 ,   6.98 ,  176.44 , 182.13 ,
	 6.98 ,   8.52 ,  176.61 , 181.67 ,
	 8.52 ,  10.06 ,  176.75 , 181.78 ,
	10.06 ,  11.60 ,  176.89 , 181.90 ,
	11.60 ,  13.14 ,  177.03 , 182.01 ,
	13.14 ,  14.68 ,  177.18 , 182.12 ,
	14.68 ,  16.22 ,  177.32 , 182.24 ,
	16.22 ,  17.76 ,  177.46 , 182.35 ,
	17.76 ,  19.30 ,  177.60 , 182.46 ,
	19.30 ,  20.84 ,  177.74 , 182.58 ,
	20.84 ,  22.38 ,  177.88 , 182.69 ,
	22.38 ,  23.92 ,  178.02 , 182.80 ,
	23.92 ,  25.46 ,  178.16 , 182.92 ,
	25.46 ,  27.00 ,  178.30 , 183.03 ,
	27.00 ,  28.54 ,  178.44 , 183.15 ,
	28.54 ,  30.08 ,  178.58 , 183.26 ,
	30.08 ,  31.62 ,  178.72 , 183.37 ,
	31.62 ,  33.16 ,  178.86 , 183.49 ,
	33.16 ,  34.70 ,  179.00 , 183.60 	],
#plate 2 : I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
    [ 9.62,	11.36,	176.09,	182.91,
      11.36,	13.11,	176.29,	182.38,
      13.11,	14.86,	176.48,	182.11,
      14.86,	16.60,	176.67,	182.08,
      16.60,	18.35,	176.86,	182.06,
      18.35,	20.09,	177.05,	182.06,
      20.09,	21.84,	177.25,	182.11,
      21.84,	23.59,	177.43,	182.21,
      23.59,	25.33,	177.61,	182.31,
      25.33,	27.08,	177.77,	182.44,
      27.08,	28.82,	177.92,	182.57,
      28.82,	30.57,	178.06,	182.71,
      30.57,	32.32,	178.20,	182.83,
      32.32,	34.06,	178.36,	182.96,
      34.06,	35.81,	178.52,	183.09,
      35.81,	37.56,	178.68,	183.23,
      37.56,	39.30,	178.84,	183.35,
      39.30,	41.05,	179.00,	183.48,
      41.05,	42.79,	179.16,	183.61,
      42.79,	44.54,	179.32,	183.74 ],

#3 plate
[	15.30 ,  17.25 ,  175.93 , 182.78 ,
	17.25 ,  19.20 ,  176.14 , 182.56 ,
	19.20 ,  21.15 ,  176.35 , 182.36 ,
	21.15 ,  23.10 ,  176.58 , 182.19 ,
	23.10 ,  25.05 ,  176.81 , 182.05 ,
	25.05 ,  27.00 ,  177.06 , 181.93 ,
	27.00 ,  28.95 ,  177.28 , 181.98 ,
	28.95 ,  30.90 ,  177.46 , 182.12 ,
	30.90 ,  32.85 ,  177.64 , 182.26 ,
	32.85 ,  34.80 ,  177.82 , 182.41 ,
	34.80 ,  36.75 ,  178.00 , 182.55 ,
	36.75 ,  38.70 ,  178.18 , 182.70 ,
	38.70 ,  40.65 ,  178.35 , 182.84 ,
	40.65 ,  42.60 ,  178.53 , 182.98 ,
	42.60 ,  44.55 ,  178.71 , 183.13 ,
	44.55 ,  46.50 ,  178.89 , 183.27 ,
	46.50 ,  48.45 ,  179.07 , 183.41 ,
	48.45 ,  50.40 ,  179.25 , 183.56 ,
	50.40 ,  52.35 ,  179.42 , 183.70 ,
	52.35 ,  54.30 ,  179.60 , 183.84 	],
#plate 4: I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
    [ 20.92, 23.07,	175.82,	182.56,
      23.07,	25.23,	176.05,	182.20,
      25.23,	27.38,	176.28,	182.02,
      27.38,	29.54,	176.51,	182.01,
      29.54,	31.69,	176.74,	182.02,
      31.69,	33.85,	176.99,	182.04,
      33.85,	36.00,	177.23,	182.11,
      36.00,	38.15,	177.47,	182.20,
      38.15,	40.31,	177.71,	182.31,
      40.31,	42.46,	177.94,	182.45,
      42.46,	44.62,	178.16,	182.60,
      44.62,	46.77,	178.38,	182.76,
      46.77,	48.93,	178.57,	182.91,
      48.93,	51.08,	178.77,	183.06,
      51.08,	53.24,	178.97,	183.22,
      53.24,	55.39,	179.16,	183.38,
      55.39,	57.55,	179.36,	183.54,
      57.55,	59.70,	179.56,	183.70,
      59.70,	61.86,	179.77,	183.87,
      61.86,	64.01,	179.97,	184.03  ],


#5 plate    
      [ 26.60 ,  28.96 ,  175.71 , 182.49 ,                    
	28.96 ,  31.32 ,  175.95 , 182.40 ,                   
	31.32 ,  33.68 ,  176.20 , 182.32 ,                   
	33.68 ,  36.04 ,  176.46 , 182.27 ,                   
	36.04 ,  38.40 ,  176.73 , 182.23 ,                   
	38.40 ,  40.76 ,  177.02 , 182.22 ,                   
	40.76 ,  43.12 ,  177.31 , 182.22 ,                   
	43.12 ,  45.48 ,  177.62 , 182.23 ,                   
	45.48 ,  47.84 ,  177.93 , 182.28 ,                   
	47.84 ,  50.20 ,  178.15 , 182.46 ,                   
	50.20 ,  52.56 ,  178.36 , 182.63 ,                   
	52.56 ,  54.92 ,  178.58 , 182.81 ,                   
	54.92 ,  57.28 ,  178.80 , 182.98 ,                   
	57.28 ,  59.64 ,  179.01 , 183.16 ,                   
	59.64 ,  62.00 ,  179.23 , 183.33 ,                   
	62.00 ,  64.36 ,  179.45 , 183.51 ,                   
	64.36 ,  66.72 ,  179.67 , 183.68 ,                   
	66.72 ,  69.08 ,  179.88 , 183.85 ,                   
	69.08 ,  71.44 ,  180.10 , 184.03 ,                   
	71.44 ,  73.80 ,  180.32 , 184.21	],  

#plate 6: I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
    [  32.22, 34.78, 175.54, 182.21,
       34.78,	37.34,	175.81,	182.01,
       37.34,	39.91,	176.07,	181.94,
       39.91,	42.47,	176.35,	181.95,
       42.47,	45.03,	176.63,	181.98,
       45.03,	47.60,	176.92,	182.02,
       47.60,	50.16,	177.22,	182.10,
       50.16,	52.72,	177.52,	182.20,
       52.72,	55.29,	177.82,	182.32,
       55.29,	57.85,	178.11,	182.47,
       57.85,	60.41,	178.41,	182.63,
       60.41,	62.98,	178.69,	182.81,
       62.98,	65.54,	178.95,	182.98,
       65.54,	68.10,	179.18,	183.17,
       68.10,	70.67,	179.41,	183.35,
       70.67,	73.23,	179.65,	183.54,
       73.23,	75.79,	179.88,	183.73,
       75.79,	78.36,	180.13,	183.93,
       78.36,	80.92,	180.37,	184.13,
       80.92,	83.48,	180.62,	184.33  ],
                
#7 plate
      [ 37.90  , 40.67  , 175.38 , 182.05 ,
	40.67  , 43.44  , 175.66 , 182.04 ,
	43.44  , 46.21  , 175.95 , 182.04 ,
	46.21  , 48.98  , 176.25 , 182.05 ,
	48.98  , 51.75  , 176.57 , 182.08 ,
	51.75  , 54.52  , 176.89 , 182.13 ,
	54.52  , 57.29  , 177.22 , 182.19 ,
	57.29  , 60.06  , 177.57 , 182.27 ,
	60.06  , 62.83  , 177.92 , 182.36 ,
	62.83  , 65.60  , 178.29 , 182.47 ,
	65.60  , 68.37  , 178.66 , 182.66 ,
	68.37  , 71.14  , 178.90 , 182.84 ,
	71.14  , 73.91  , 179.15 , 183.02 ,
	73.91  , 76.68  , 179.39 , 183.22 ,
	76.68  , 79.45  , 179.65 , 183.42 ,
	79.45  , 82.22  , 179.90 , 183.62 ,
	82.22  , 84.99  , 180.16 , 183.84 ,
	84.99  , 87.76  , 180.43 , 184.06 ,
	87.76  , 90.53  , 180.71 , 184.28 ,
	90.53  , 93.30  , 180.98 , 184.51	],

#plate 8: I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
    [ 43.52, 46.49,	175.27,	181.86,
      46.49,	49.46,	175.56,	181.83,
      49.46,	52.43,	175.87,	181.85,
      52.43,	55.40,	176.19,	181.88,
      55.40,	58.38,	176.51,	181.93,
      58.38,	61.35,	176.86,	182.00,
      61.35,	64.32,	177.20,	182.09,
      64.32,	67.29,	177.56,	182.20,
      67.29,	70.26,	177.92,	182.32,
      70.26,	73.24,	178.28,	182.48,
      73.24,	76.21,	178.65,	182.66,
      76.21,	79.18,	179.00,	182.85,
      79.18,	82.15,	179.32,	183.06,
      82.15,	85.12,	179.58,	183.27,
      85.12,	88.10,	179.85,	183.48,
      88.10,	91.07,	180.13,	183.70,
      91.07,	94.04,	180.40,	183.92,
      94.04,	97.01,	180.69,	184.15,
      97.01,	99.98,	180.97,	184.38,
      99.98,	102.95,	181.26,	184.62 ],

#9 plate
      [ 49.20  , 52.38  , 175.06 , 181.61 ,
	52.38  , 55.56  , 175.37 , 181.64 ,
	55.56  , 58.74  , 175.70 , 181.70 ,
	58.74  , 61.92  , 176.04 , 181.76 ,
	61.92  , 65.10  , 176.39 , 181.84 ,
	65.10  , 68.28  , 176.76 , 181.94 ,
	68.28  , 71.46  , 177.13 , 182.05 ,
	71.46  , 74.64  , 177.52 , 182.17 ,
	74.64  , 77.82  , 177.91 , 182.31 ,
	77.82  , 81.00  , 178.32 , 182.46 ,
	81.00  , 84.18  , 178.75 , 182.63 ,
	84.18  , 87.36  , 179.20 , 182.85 ,
	87.36  , 90.54  , 179.48 , 183.06 ,
	90.54  , 93.72  , 179.76 , 183.28 ,
	93.72  , 96.90  , 180.05 , 183.51 ,
	96.90  ,100.08  , 180.35 , 183.75 ,
	100.08 ,103.26  , 180.65 , 183.99 ,
	103.26 ,106.44  , 180.96 , 184.24 ,
	106.44 ,109.62  , 181.27 , 184.49 ,
	109.62 ,112.80  , 181.59 , 184.75	],

#plate 10: I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
    [ 54.82, 58.20,	174.99,	181.52,
      58.20,	61.58,	175.32,	181.64,
      61.58,	64.96,	175.67, 181.76,
      64.96,	68.34,	176.03,	181.82,
      68.34,	71.72,	176.40,	181.89,
      71.72,	75.10,	176.79,	181.98,
      75.10,	78.48,	177.19,	182.08,
      78.48,	81.86,	177.60,	182.20,
      81.86,	85.24,	178.02,	182.33,
      85.24,	88.62,	178.45,	182.49,
      88.62,	92.00,	178.89,	182.69,
      92.00,	95.38,	179.31,	182.90,
      95.38,	98.76,	179.69,	183.13,
      98.76,	102.14,	179.99,	183.37,
      102.14,	105.52,	180.30,	183.61,
      105.52, 108.90,	180.61,	183.85,
      108.90,	112.28,	180.92,	184.11,
      112.28,	115.67,	181.25,	184.37,
      115.67,	119.05,	181.58,	184.64,
      119.05,	122.43,	181.91,	184.91  ],

#11 plate
      [ 60.40  , 63.98  , 174.92 , 181.35 ,
	63.98  , 67.56  , 175.27 , 181.43 ,
	67.56  , 71.14  , 175.63 , 181.52 ,
	71.14  , 74.72  , 176.00 , 181.63 ,
	74.72  , 78.30  , 176.38 , 181.75 ,
	78.30  , 81.88  , 176.78 , 181.88 ,
	81.88  , 85.46  , 177.18 , 182.02 ,
	85.46  , 89.04  , 177.60 , 182.18 ,
	89.04  , 92.62  , 178.03 , 182.35 ,
	92.62  , 96.20  , 178.47 , 182.53 ,
	96.20  , 99.78  , 178.92 , 182.73 ,
	99.78  ,103.36  , 179.39 , 182.94 ,
	103.36 , 106.94 ,  179.88,  183.19,
	106.94 , 110.52 ,  180.20,  183.43,
	110.52 , 114.10 ,  180.52,  183.69,
	114.10 , 117.68 ,  180.84,  183.94,
	117.68 , 121.26 ,  181.17,  184.21,
	121.26 , 124.84 ,  181.51,  184.48,
	124.84 , 128.42 ,  181.86,  184.76,
	128.42 , 132.00 ,  182.21,  185.04	]
);

 my $Nslit  = 30;
 my $Nblock  = 20;

# my @offset =(-5.6, -4.4, -3.3, -2.1, -0.9, 0.1);  #according to Seamus, it's substraction 
 my @offset =(-5.6, -5.00,-4.4, -3.85,-3.3, -2.70,-2.1, -1.55,-0.9, -0.4,0.1);  #according above values, I extrapolated for new z locations. Rakitha Thu Aug 29 15:05:57 EDT 2013
 

 for(my $n=2; $n<=($Nplate-3); $n++){
	my $n_c     = cnumber($n-1, 1);
  for(my $i=1; $i<=$Nslit; $i++){  # making all slits
	my $slit_rotation = ($i-1)*12-$offset[$n-1]; #note the minus sign here
	my $i_c     = cnumber($i-1, 10);
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_plate$n_c.slit$i_c";
	$detector{"mother"}      = "$DetectorName\_plate$n_c";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm 0*cm";
	$detector{"rotation"}    = "0*deg 0*deg $slit_rotation*deg";
	$detector{"color"}       = "000000";
	$detector{"type"}        = "Polycone";
	my $Rin  = $x[$n-1][0];
	my $Rout = $x[$n-1][77];
	my $Rin2  = $x[$n][0];
	my $Rout2 = $x[$n][77];
	my $Sphi = 0.5*($x[$n-1][2]+$x[$n-1][2]+$offset[$n-1]-$offset[$n]);
	my $Dphi = 12;
	$detector{"dimensions"}  = "$Sphi*deg $Dphi*deg 2*counts $Rin*cm $Rin2*cm $Rout*cm $Rout2*cm 0.0*cm $Dz_shld*cm";
	$detector{"material"}    = "$material_baffle_within";
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 0;
	$detector{"style"}       = 0;
	$detector{"sensitivity"} = "$sensitivity_baffle";
	$detector{"hit_type"}    = "$hit_baffle";
 	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);
    
      for(my $j=1; $j<=$Nblock; $j++){ # making blocks within slits
	my $j_c     = cnumber($j-1, 10);
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_plate$n_c.slit$i_c.block$j_c";
	$detector{"mother"}      = "$DetectorName\_plate$n_c.slit$i_c";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "0*cm 0*cm 0*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "$color_abs";
	$detector{"type"}        = "Polycone";
	my $Rin  = $x[$n-1][($j-1)*4+0];
	my $Rout = $x[$n-1][($j-1)*4+1];
	my $Rin2  = $x[$n][($j-1)*4+0];
	my $Rout2 = $x[$n][($j-1)*4+1];
	my $Sphi = 0.5*($x[$n-1][($j-1)*4+2]+$x[$n][($j-1)*4+2]);
 	my $Dphi = 0.5*($x[$n-1][($j-1)*4+3]-$x[$n-1][($j-1)*4+2]+1+$x[$n][($j-1)*4+3]-$x[$n][($j-1)*4+2]+1);
	$detector{"dimensions"}  = "$Sphi*deg $Dphi*deg 2*counts $Rin*cm $Rin2*cm $Rout*cm $Rout2*cm 0.0*cm $Dz_shld*cm";
	$detector{"material"}    = "$material_baffle";
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "$sensitivity_baffle";
	$detector{"hit_type"}    = "$hit_baffle";
 	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);
      }

    }
  }

#  for(my $n=1; $n<=$Nplate; $n++){
#   for(my $i=1; $i<=$Nslit; $i++){
#     for(my $j=1; $j<=$Nblock; $j++){# 
# 	my $n_c     = cnumber($n-1, 1);
# 	my $i_c     = cnumber($i-1, 10);
# 	my $j_c     = cnumber($j-1, 10);
# 
# 	$detector{"name"}        = "$DetectorName\_plate$n_c.slit$i_c.block$j_c";
# 	$detector{"mother"}      = "$DetectorName";
# 	$detector{"description"} = $detector{"name"};
# 	$detector{"pos"}         = "0*cm 0*cm $PlateZ[$n-1]*cm";
# 	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
# 	$detector{"color"}       = "$color_abs";
# 	$detector{"type"}        = "Polycone";
# 	my $Rin  = $x[$n-1][($j-1)*4+0];
# 	my $Rout = $x[$n-1][($j-1)*4+1];
# 	my $Sphi = $x[$n-1][($j-1)*4+2]+($i-1)*12;
#  	my $Dphi = $x[$n-1][($j-1)*4+3]-$x[$n-1][($j-1)*4+2];
# 	$detector{"dimensions"}  = "$Rin*cm $Rout*cm $PlateDz*cm $Sphi*deg $Dphi*deg";
# 	$detector{"material"}    = "$material_baffle";
# 	$detector{"mfield"}      = "no";
# 	$detector{"ncopy"}       = ($n-1)*$Nplate+($i-1)*$Nslit+$j;
# 	$detector{"pMany"}       = 1;
# 	$detector{"exist"}       = 1;
# 	$detector{"visible"}     = 1;
# 	$detector{"style"}       = 1;
# 	$detector{"sensitivity"} = "$hit_baffle";
# 	$detector{"hit_type"}    = "$hit_baffle";
#  	$detector{"identifiers"} = $detector{"name"};
# 
# 	print_det(\%configuration, \%detector);
#      }
#    }
#  }
}
