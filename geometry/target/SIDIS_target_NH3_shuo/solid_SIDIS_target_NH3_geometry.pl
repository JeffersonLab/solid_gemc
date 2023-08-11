use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_SIDIS_target_NH3';

my $DetectorMother="root";

sub solid_SIDIS_target_NH3
{
make_target_field();
make_scattering_chamber();
make_scattering_windows();
make_target();
make_target_endcaps();
make_target_LHe();
make_target_tail_nose();
make_target_4Kshield();
make_target_LN2shield();
#make_target_steel();
#make_target_coil();
#make_target_coil_box();
#make_target_coil_box_2();
}


sub make_target_field
{
 my $NUM  = 1;
 my @z    = (-350);
#  my @Rin  = (0);
#  my @Rout = (65);
#  my @Dz   = (65);
  my @Dx = (65);
  my @Dy = (65);
  my @Dz = (65);
 my @name = ("field"); 
 my @mother = ("$DetectorMother");
 my @mat  = ("G4_AIR");
 
 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "ff0000";
    $detector{"type"}       = "Box";
    $detector{"dimensions"} = "$Dx[$n-1]*cm $Dy[$n-1]*cm $Dz[$n-1]*cm";    
#     $detector{"type"}       = "Tube";
#     $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
#     $detector{"type"}       = "Polycone";
#     $detector{"dimensions"} = "0*deg 360*deg 4*counts 0*cm 0*cm 1.1*cm 3*cm 65*cm 65*cm 65*cm 65*cm -65*cm 25*cm 25*cm 65*cm";
    $detector{"material"}   = $mat[$n-1];
#     $detector{"mfield"}     = "solenoid_ptarget";
#     $detector{"mfield"}     = "g2p_ptarget";
    $detector{"mfield"}     = "oxford_ptarget";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 0;  # style 0 shows only borders
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}

# (according to Kalyan in 2014 and 2016)  Scattering chamber : cylinder of  74cm(height) with inner diameter of 45cm and thickness of 2.5cm made of Aluminum.
# (according to oxford 2012 design "CRABB DESIGN STUDY - October 2012-1.pdf"), in 2020, Zhiwen change to 76cm height and 68cm diameter

my $chamber_height=76;
my $chamber_diameter=68;
my $chamber_thk=2.5;
my $chamber_inner=$chamber_diameter/2-$chamber_thk;

sub make_scattering_chamber
{
 my $NUM  = 2;
 my @Rin  = (0,0);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = ($chamber_diameter/2,$chamber_inner);
 my @Dz   = ($chamber_height/2,$chamber_height/2); 
 my @name = ("SC_out","SC_in");
 my @mother = ("$DetectorName\_field","$DetectorName\_SC_out"); 
 my @mat  = ("G4_Al","G4_Galactic");
 my @rot  = (90,0);
 my @color = ("FF6600","FFFFFF");
 
 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "$rot[$n-1]*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];    
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
    $detector{"material"}   = $mat[$n-1];
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

# sub make_scattering_chamber
# {
#  my $NUM  = 4;
#  my @y    = (26.0,-26.0,0.0,0.0);
# #  my @z    = (0.0-350,0.0-350,0.0-350,0.0-350);
#  my @z    = (0,0,0,0); 
#  my @Rin  = (22.5,22.5,22.5,22.5);
#  my @Rout = (25.,25.,25.,25.);
#  my @Dz   = (11.0,11.0,15.0,15.0);
#  my @SPhi = (0.0,0.0,95.0,300.0);
#  my @DPhi = (360.0,360.0,145.0,145.0);
#  my @name = ("SC1","SC2","SC3","SC4");
#  my @mother = ("$DetectorName\_field","$DetectorName\_field","$DetectorName\_field","$DetectorName\_field"); 
#  my @mat  = ("G4_Al","G4_Al","G4_Al","G4_Al");
# 
#  for(my $n=1; $n<=$NUM; $n++)
#  {
#      my %detector=init_det(); 
#     $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#     $detector{"mother"}      = "$mother[$n-1]" ;
#     $detector{"description"} = "$DetectorName\_$name[$n-1]";
#     $detector{"pos"}        = "0*cm $y[$n-1]*cm $z[$n-1]*cm";
#     $detector{"rotation"}   = "0*deg 0*deg 0*deg";
#     $detector{"color"}      = "FF6600";
#     $detector{"type"}       = "Tube";
#     $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm $SPhi[$n-1]*deg $DPhi[$n-1]*deg";
#     $detector{"material"}   = $mat[$n-1];
#     $detector{"mfield"}     = "no";
#     $detector{"ncopy"}      = 1;
#     $detector{"pMany"}       = 1;
#     $detector{"exist"}       = 1;
#     $detector{"visible"}     = 1;
#     $detector{"style"}       = 1;
#     $detector{"sensitivity"} = "no";
#     $detector{"hit_type"}    = "no";
#     $detector{"identifiers"} = "no";
#     print_det(\%configuration, \%detector);
#  }
# }

# some info in 2014 below is not good 
#2. The scattering chamber has entrance and exit windows on them, made of thin Al ( 0.04cm)
#  a) beam entrance window dimensions can be small:  20cm wide  and 30cm  height
#  b) beam exit  window dimensions should be large enough to cover +/- 28deg in the scattering particles. You can keep the same height (30cm) but need at least 100cm width to cover angular range.

# (according to "g2p_target4.xls"), scattering chamber has entrance windows is Al (0.02cm), exit window is Al (0.05 cm) 

sub make_scattering_windows
{
 my $NUM  = 4;
#  my @z    = (0.0-350,0.0-350);
 my @z    = (0,0,0,0); 
#  my @Rin  = (22.5,24.96,22.5,24.96);
#  my @Rout = (24.96,25,24.96,25);
#  my @Dz   = (5,5,10,10);
#  my @SPhi = (85,85,245,245);
#  my @DPhi = (10,10,50,50);
 
 my @Rin  = ($chamber_inner,$chamber_inner+($chamber_thk-0.02),$chamber_inner,$chamber_inner+($chamber_thk-0.04));
 my @Rout = ($chamber_inner+($chamber_thk-0.02),$chamber_diameter/2,$chamber_inner+($chamber_thk-0.04),$chamber_diameter/2);
 my @Dz   = (6,6,16.6,16.6);
 my @SPhi = (80,80,242,242);
 my @DPhi = (20,20,56,56); 
 my @name = ("entrance_cut","entrance_win","exit_cut","exit_win");
 my @mother = ("$DetectorName\_SC_out","$DetectorName\_SC_out","$DetectorName\_SC_out","$DetectorName\_SC_out");
 my @mat  = ("G4_Galactic","G4_Al","G4_Galactic","G4_Al");
 my @color = ("FFFFFF","000000","FFFFFF","000000");

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm $SPhi[$n-1]*deg $DPhi[$n-1]*deg";
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

sub make_target
{
 my $NUM  = 1;
#  my @z    = (0.0-350);
 my @z    = (0);
 my @Rin  = (0.0);
 my @Rout = (1.413);
 my @Dz   = (1.423);
 my @name = ("target"); 
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("SL_target_NH3_NH3He4");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "ff0000";
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

# For rest of the items below I had to use  position in y-direction (as suppose to z) to place the material. Maybe due to mother volume rotated already..but I couldn't figure our how to use z direction to place the items ..yet. 
sub make_target_endcaps
{
 my $NUM  = 2;
 my @y    = (-1.426,1.426); # y instead of z !! ?
 my @Rin  = (0.0,0.0);
 my @Rout = (1.613,1.613);
 my @Dz   = (0.000889,0.000889);
 my @name = ("up_endcap","down_endcap"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";  # y instead of z !! ?
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "808080";
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

sub make_target_LHe
{
 my $NUM  = 2;
 my @y    = (-1.65,1.65);
 my @Rin  = (0.0,0.0);
 my @Rout = (1.613,1.613);
 my @Dz   = (0.218,0.218);
 my @name = ("up_LHe","down_LHe"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("SL_target_NH3_He4_liquid","SL_target_NH3_He4_liquid");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "00BFFF";
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

sub make_target_tail_nose
{
 my $NUM  = 2;
 my @y    = (-2.65,2.65);
 my @Rin  = (0.0,0.0);
 my @Rout = (1.613,1.613);
 my @Dz   = (0.006345,0.006345);
 my @name = ("up_tail_nose","down_tail_nose"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "8A2BE2";
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


sub make_target_4Kshield
{
 my $NUM  = 2;
 my @y    = (-3.65,3.65);
 my @Rin  = (0.0,0.0);
 my @Rout = (1.613,1.613);
 my @Dz   = (0.0006345,0.0006345);
 my @name = ("up_4Kshield","down_4Kshield"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "2F4F4F";
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


sub make_target_LN2shield
{
 my $NUM  = 2;
 my @y    = (-4.65,4.65);
 my @Rin  = (0.0,0.0);
 my @Rout = (1.613,1.613);
 my @Dz   = (0.0019,0.0019);
 my @name = ("up_LN2shield","down_LN2shield"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";
    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
    $detector{"color"}      = "FF8C00";
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

sub make_target_steel
{
 my $NUM  = 2;
 my @x    = (-5.55,5.55);
 my @Rin  = (5.9,5.9);
 my @Rout = (7.9,7.9);
 my @Dz   = (1.95,1.95);
 my @name = ("left_steel","right_steel"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_STAINLESS-STEEL","G4_STAINLESS-STEEL");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
    $detector{"color"}      = "00BFFF";
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

sub make_target_coil
{
 my $NUM  = 8;
 my @x    = (-13.5,-11.15,-9.1,-6.95,6.95,9.1,11.15,13.5);
 my @Rin  = (19.9,16.1,13,9.9,9.9,13,16.1,19.9);
 my @Rout = (22.9,18.6,15,11.4,11.4,15,18.6,22.9);
 my @Dz   = (2.3,1.85,1.5,1.15,1.15,1.5,1.85,2.3);
 my @name = ("left_coil_xl","left_coil_l","left_coil_m","left_coil_s","right_coil_s","right_coil_m","right_coil_l","right_coild_xl"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
    $detector{"color"}      = "FF8C00";
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

sub make_target_coil_box
{
 my $NUM  = 2;
 my @x    = (-3.6,3.6,-3.6,3.6);
# my @numZPlane = 15;
 my @Rot  = (90,-90,90,-90);
 my @Rin  = (7.9,    7.9,   5.88,  5.88,   5.88,   6.58125,11.685, 11.685,15.212,15.212, 15.212, 15.212, 15.212, 18.856, 18.856);
 my @Rout = (8.57803,17.015,17.015,17.853, 8.75,   8.75,   17.853, 23.114,23.114,23.8,   23.8,   25,     25,     25,     25);
 my @zPln = (0,      4,     4,     4.39751,4.39751,5.51756,4.39751,6.8915,6.8915,7.21668,8.27625,8.27625,9.28451,9.28451,12.14951);
 my @name = ("left_coil_box","right_coil_box"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"rotation"}   = "0*deg $Rot[$n-1]*deg 0*deg";
    $detector{"color"}      = "2F4F4F";
    $detector{"type"}       = "Polycone";
    if($n<=2){
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"dimensions"} = "0*deg 360*deg 6*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm";
  } 
    if($n>=3){
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"dimensions"} = "0*deg 360*deg 9*counts $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm";
    }
    
#"0*deg 360*deg 15*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm ";
#"0*deg 360*deg 3*counts 7.9*cm 7.9*cm 5.88*cm 8.58*cm 17.015*cm 17.015*cm 0*cm 4*cm 4*cm"; 
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

sub make_target_coil_box_2
{
 my $NUM  = 2;
 my @x    = (-3.6,3.6,-3.6,3.6);
# my @numZPlane = 15;
 my @Rot  = (90,-90,90,-90);
 my @Rin  = (7.9,    7.9,   5.88,  5.88,   5.88,   6.58125,11.685, 11.685,15.212,15.212, 15.212, 15.212, 15.212, 18.856, 18.856);
 my @Rout = (8.57803,17.015,17.015,17.853, 8.75,   8.75,   17.853, 23.114,23.114,23.8,   23.8,   25,     25,     25,     25);
 my @zPln = (0,      4,     4,     4.39751,4.39751,5.51756,4.39751,6.8915,6.8915,7.21668,8.27625,8.27625,9.28451,9.28451,12.14951);
 my @name = ("left_coil_box_2","right_coil_box_2"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"rotation"}   = "0*deg $Rot[$n-1]*deg 0*deg";
    $detector{"color"}      = "2F4F4F";
    $detector{"type"}       = "Polycone";
    # if($n<=2){
    # $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    # $detector{"dimensions"} = "0*deg 360*deg 6*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm";
    # } 
    # if($n>=3){
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"dimensions"} = "0*deg 360*deg 9*counts $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm";
    #}
    
#"0*deg 360*deg 15*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm ";
#"0*deg 360*deg 3*counts 7.9*cm 7.9*cm 5.88*cm 8.58*cm 17.015*cm 17.015*cm 0*cm 4*cm 4*cm"; 
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

sub make_target_coil_lid
{
 my $NUM  = 2;
 my @x    = (-8,8);
# my @numZPlane = 15;
 my @Rot  = (90,-90);
 my @Rin  = (8.75,  8.75,   8.825,  9.954, 9.954, 12.736,12.736,16.068,16.068,17.57, 17.75); 
 my @Rout = (11.685,11.685, 11.685, 11.685,15.212,15.212,18.851,18.851,25,    25,    24.8);
 my @zPln = (8,     9.384,  9.384,  10.494,10.494,12.887,12.887,15.752,15.752,16.852,17.052);
 my @name = ("left_coil_lid","right_coil_lid"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"rotation"}   = "0*deg $Rot[$n-1]*deg 0*deg";
    $detector{"color"}      = "2F4F4F";
    $detector{"type"}       = "Polycone";
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"dimensions"} = "0*deg 360*deg 11*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm"; 
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

