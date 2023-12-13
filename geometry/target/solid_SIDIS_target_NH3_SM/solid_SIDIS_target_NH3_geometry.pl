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
make_30K_shield();
make_30K_window();
make_target_LHe();
make_target_LHe_shield();
make_target();
make_target_endcaps();
make_target_cell();
make_target_steel();
make_target_coil_box();
make_target_coil_box_2();
#make_target_coil_1();#somehow the solid_slice.vis would change the rotation 
make_target_coil();
make_target_coil_lid();
make_magnet_support();
#Not used
#make_5T_JLab();#just for map the magnetic field purpose
#make_vacuum_4K();#dimensions are from the g2p target from ChaoGu
#make_vacuum_LN2();#from ChaoGu
}

my $target_l=2.82702;#cm/1.113"
my $target_r=1.36144;#cm/0.536"
my $cap_thick=0.001778;#cm/0.7mil
my $cell_thick=0.0889;#cm
my $LHe_r=2.10058;#nose_r in chaogu's file cm
my $LHe_l=22.86;#cm
my $LHe_shield_thick=0.0101;#cm

#Nov3,2023 update
#my $target_l=3;
#my $LHe_r=3;


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
#     $detector{"mfield"}     = "oxford_ptarget";
    $detector{"mfield"}     = "no";
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
#my $chamber_height=76;
#my $chamber_diameter=68;
#my $chamber_thk=2.5;
#my $chamber_inner=$chamber_diameter/2-$chamber_thk;
my $chamber_height=86.36;
my $chamber_out=95.885/2;
my $chamber_thk=2.54;
my $chamber_inner=$chamber_out-$chamber_thk;
my $window_beam_thick=0.02032;#from Chao Gu
my $window_beam_inner= $chamber_thk-$window_beam_thick;#for the vacuum inside the chamber 
my $window_exit_thick=0.05;#from Chao Gu
my $window_exit_inner=$chamber_thk-$window_exit_thick;#for the vacuum inside the chamber

sub make_scattering_chamber
{
 my $NUM  = 2;
 my @Rin  = (0,0);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = ($chamber_out,$chamber_inner);
 my @Dz   = ($chamber_height/2,$chamber_height/2); 
 my @name = ("SC_out","SC_in");
 my @mother = ("$DetectorName\_field","$DetectorName\_SC_out"); 
 my @mat  = ("G4_Al","G4_Galactic");
 my @rot  = (90,0);
 #my @color = ("ff0000","ff0000");
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
 
 my @Rin = ($chamber_inner,$chamber_inner+$window_beam_inner,$chamber_inner,$chamber_inner+$window_exit_inner);
 my @Rout = ($chamber_inner+$window_beam_inner,$chamber_out,$chamber_inner+$window_exit_inner,$chamber_out);
 #my @Rin  = ($chamber_inner,$chamber_inner+($window_beam_inner),$chamber_inner+$window_beam_inner+$window_beam_thick,$chamber_inner,$chamber_inner+($window_exit_inner),$chamber_inner+$window_exit_inner+$window_exit_thick);
 #my @Rout = ($chamber_inner+($window_beam_inner),$chamber_inner+$window_beam_inner+$window_beam_thick,$chamber_out,$chamber_inner+($window_exit_inner),$chamber_inner+$window_exit_inner+$window_exit_thick,$chamber_out);
 my @Dz   = (19,19,19,19);
 my @SPhi = (80,80,242,242);
 my @DPhi = (20,20,56,56); 
 my @name = ("entrance_cut","entrance_win","exit_cut","exit_win");
 my @mother = ("$DetectorName\_SC_out","$DetectorName\_SC_out","$DetectorName\_SC_out","$DetectorName\_SC_out");
 my @mat  = ("G4_Galactic","G4_Al","G4_Galactic","G4_Al");
 my @color = ("FFFFFF","000000","FFFFFF","000000");
 my @style = (0,1,0,1);

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
    $detector{"style"}       = $style[$n-1];
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}

sub make_target_LHe
{
 my $NUM  = 1;
 my @Rin  = (0);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = ($LHe_r);
 my @Dz   = ($LHe_l/2); 
 my @name = ("target_LHe");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("SL_target_NH3_He4_liquid");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("FF6600");

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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

sub make_target_LHe_shield
{
 my $NUM  = 1;
 my @Rin  = ($LHe_r);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = ($LHe_r+$LHe_shield_thick);
 my @Dz   = ($LHe_l/2); 
 my @name = ("target_LHe_shield");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("FFFFFF");

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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

#check:what is this volume in CAD
sub make_5T_JLab
{
 my $NUM  = 1;
 my @Rin  = (3.49);#Looks like a whole volume
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = (3.5);
 my @Dz   = ((116.94-56)/2); 
 my @name = ("target_5T");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("FF6600");
 my @z = (5);

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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

#check: from ChaoGu, this size exceed the opening length of the two iron ring
sub make_vacuum_4K
{
 my $NUM  = 1;
 my @Rin  = (3.81);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = (3.8138);
 my @Dz   = (135.89/2); 
 my @name = ("target_vacuum_4K");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("ff0000");

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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

#from ChaoGu
sub make_vacuum_LN2
{
 my $NUM  = 1;
 my @Rin  = (41.91);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = (41.9138);
 my @Dz   = (135.89/2); 
 my @name = ("target_vacuum_LN2");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("ff0000");

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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

my $shield_30K_Rin=38;#cm
my $shield_30K_Rout=39.5;#cm
my $shield_30K_height=75/2;#cm
my $shield_30K_window_thick=0.02;#cm 200um
sub make_30K_shield
{
 my $NUM  = 1;
 my @Rin  = ($shield_30K_Rin);
#  my @Rout = (25,22.5);
#  my @Dz   = (37,37); 
 my @Rout = ($shield_30K_Rout);
 my @Dz   = ($shield_30K_height); 
 #my @Dz   = (105/2); 
 my @name = ("30K");
 my @mother = ("$DetectorName\_SC_in"); 
 my @mat  = ("G4_Al");
 my @rot  = (90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("FF6600");
 my @z = (0);#move it to center
 #my @z = (5);#estimate

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
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
sub make_30K_window
{
 my $NUM  = 4;
 my @Rin  = ($shield_30K_Rin,$shield_30K_Rout-$shield_30K_window_thick,$shield_30K_Rin,$shield_30K_Rout-$shield_30K_window_thick);
 #my @Rout = ($shield_30K_Rout,$shield_30K_Rout);
 my @Rout = ($shield_30K_Rout-$shield_30K_window_thick,$shield_30K_Rout,$shield_30K_Rout-$shield_30K_window_thick,$shield_30K_Rout);
 my @Dz   = (10,10,20,20);#cm estimate. first is for beam entrance 
 #my @Dz   = (105/2); 
 my @name = ("30K_inwindow_vacuum","30K_inwindow_Al","30K_outwindow_vacuum","30K_outwindow_Al");
 #my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mother = ("$DetectorName\_30K","$DetectorName\_30K","$DetectorName\_30K","$DetectorName\_30K"); 
 my @mat  = ("G4_Galactic","G4_Al","G4_Galactic","G4_Al");
 my @rot  = (-90,-90,-90,-90);
 #my @color = ("ff0000","ff0000");
 #my @color = ("FF6600","FFFFFF");
 my @color = ("FFFFFF","000000","FFFFFF","000000");
 my @z = (0,0,0,0);#move it back to center
 #my @z = (-5,-5,-5,-5);#estimate
 my @SPhi = (80,80,242,242);
 my @DPhi = (20,20,56,56);
 my @style = (0,1,0,1);

 for(my $n=1; $n<=$NUM; $n++)
 {
     my %detector=init_det(); 
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
    $detector{"color"}      = $color[$n-1];    
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm $SPhi[$n-1]*deg $DPhi[$n-1]*deg";
    $detector{"material"}   = $mat[$n-1];
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = $style[$n-1];
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
 }
}
#This shield used same logic as the scattering chamber, that the inner volume is inside the out volume, then everything inside should be a daughter volume of the inner volume. 
#sub make_30K_shield
#{
# my $NUM  = 2;
# my @Rin  = (0,0);
##  my @Rout = (25,22.5);
##  my @Dz   = (37,37); 
# my @Rout = ($shield_30K_Rout,$shield_30K_Rin);
# my @Dz   = ($shield_30K_height,$shield_30K_height); 
# #my @Dz   = (105/2); 
# my @name = ("30K","30K_in");
# my @mother = ("$DetectorName\_SC_in","$DetectorName\_30K"); 
# my @mat  = ("G4_Al","G4_Galactic");
# my @rot  = (90,0);
# #my @color = ("ff0000","ff0000");
# #my @color = ("FF6600","FFFFFF");
# my @color = ("FF6600","FFFFFF");
# my @z = (0,0);#move it back to center
# #my @z = (5,5);#estimate
#
# for(my $n=1; $n<=$NUM; $n++)
# {
#     my %detector=init_det(); 
#    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#    $detector{"mother"}      = "$mother[$n-1]" ;
#    $detector{"description"} = "$DetectorName\_$name[$n-1]";
#    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
#    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
#    $detector{"color"}      = $color[$n-1];    
#    $detector{"type"}       = "Tube";
#    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
#    $detector{"material"}   = $mat[$n-1];
#    $detector{"mfield"}     = "no";
#    $detector{"ncopy"}      = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 0;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);
# }
#}
#sub make_30K_window#Not actually the window itself, but the vacuum inside the Al shield
#{
# my $NUM  = 2;
# my @Rin  = ($shield_30K_Rin,$shield_30K_Rin);
# #my @Rout = ($shield_30K_Rout,$shield_30K_Rout);
# my @Rout = ($shield_30K_Rout-$shield_30K_window_thick,$shield_30K_Rout-$shield_30K_window_thick);
# my @Dz   = (10,20);#cm estimate. first is for beam entrance 
# #my @Dz   = (105/2); 
# my @name = ("30K_inwindow","30K_outwindow");
# #my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
# my @mother = ("$DetectorName\_30K","$DetectorName\_30K"); 
# my @mat  = ("G4_Galactic","G4_Galactic");
# my @rot  = (-90,-90);
# #my @color = ("ff0000","ff0000");
# #my @color = ("FF6600","FFFFFF");
# my @color = ("FFFFFF","FFFFFF");
# my @z = (0,0);#move it back to center
# #my @z = (-5,-5);#estimate
# my @SPhi = (80,242);
# my @DPhi = (20,56);
#
# for(my $n=1; $n<=$NUM; $n++)
# {
#     my %detector=init_det(); 
#    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#    $detector{"mother"}      = "$mother[$n-1]" ;
#    $detector{"description"} = "$DetectorName\_$name[$n-1]";
#    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
#    $detector{"rotation"}   = "0*deg 0*deg $rot[$n-1]*deg";
#    $detector{"color"}      = $color[$n-1];    
#    $detector{"type"}       = "Tube";
#    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm $SPhi[$n-1]*deg $DPhi[$n-1]*deg";
#    $detector{"material"}   = $mat[$n-1];
#    $detector{"mfield"}     = "no";
#    $detector{"ncopy"}      = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);
# }
#}

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


sub make_target
{
 my $NUM  = 1;
#  my @z    = (0.0-350);
 my @z    = (0);
 my @Rin  = (0.0);
 my @Rout = ($target_r);
 my @Dz   = (($target_l-2*$cap_thick)/2);#(target_l-2*cap_thick)/2=(2.82702-2*0.001778)/2=1.4117cm
 my @name = ("target"); 
 my @mother = ("$DetectorName\_target_LHe"); 
 my @mat  = ("SL_target_NH3_NH3He4");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
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

#cap thick = 0.7mil from ChaoGu
sub make_target_endcaps
{
 my $NUM  = 2;
#  my @z    = (0.0-350);
 my @z    = ($target_l/2,-$target_l/2);#2.82702/2
 my @Rin  = (0.0,0.0);
 my @Rout = ($target_r,$target_r);#1.36144
 my @Dz   = ($cap_thick/2,$cap_thick/2);#0.7mil/2
 my @name = ("target_endcaps_u","target_endcaps_d"); 
 my @mother = ("$DetectorName\_target_LHe","$DetectorName\_target_LHe"); 
 my @mat  = ("G4_Al","G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$z[$n-1]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
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

sub make_target_cell
{
 my $NUM  = 1;
#  my @z    = (0.0-350);
 my @z    = (0);
 my @Rin  = ($target_r);
 my @Rout = ($target_r+$cell_thick);
 my @Dz   = (($target_l)/2);#(target_l-2*cap_thick)/2=(2.82702-2*0.001778)/2=1.4117cm
 my @name = ("target_cell"); 
 my @mother = ("$DetectorName\_target_LHe"); 
 my @mat  = ("G4_Al");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $z[$n-1]*cm 0*cm";
    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
    $detector{"color"}      = "FFFFFF";
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
#sub make_target_endcaps
#{
# my $NUM  = 2;
# my @y    = (-1.426,1.426); # y instead of z !! ?
# my @Rin  = (0.0,0.0);
# my @Rout = (1.613,1.613);
# my @Dz   = (0.000889,0.000889);
# my @name = ("up_endcap","down_endcap"); 
# my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
# my @mat  = ("G4_Al","G4_Al");
#
# for(my $n=1; $n<=$NUM; $n++)
# {
##     my $pnumber     = cnumber($n-1, 10);
#     my %detector=init_det();
#    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#    $detector{"mother"}      = "$mother[$n-1]" ;
#    $detector{"description"} = "$DetectorName\_$name[$n-1]";
#    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";  # y instead of z !! ?
#    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
#    $detector{"color"}      = "808080";
#    $detector{"type"}       = "Tube";
#    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
#    $detector{"material"}   = $mat[$n-1];
#    $detector{"mfield"}     = "no";
#    $detector{"ncopy"}      = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);
# }
#}

#sub make_target_LHe
#{
# my $NUM  = 2;
# my @y    = (-1.65,1.65);
# my @Rin  = (0.0,0.0);
# my @Rout = (1.613,1.613);
# my @Dz   = (0.218,0.218);
# my @name = ("up_LHe","down_LHe"); 
# my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
# my @mat  = ("SL_target_NH3_He4_liquid","SL_target_NH3_He4_liquid");
#
# for(my $n=1; $n<=$NUM; $n++)
# {
##     my $pnumber     = cnumber($n-1, 10);
#     my %detector=init_det();
#    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#    $detector{"mother"}      = "$mother[$n-1]" ;
#    $detector{"description"} = "$DetectorName\_$name[$n-1]";
#    $detector{"pos"}        = "0*cm $y[$n-1]*cm 0*cm";
#    $detector{"rotation"}   = "-90*deg 0*deg 0*deg";
#    $detector{"color"}      = "00BFFF";
#    $detector{"type"}       = "Tube";
#    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
#    $detector{"material"}   = $mat[$n-1];
#    $detector{"mfield"}     = "no";
#    $detector{"ncopy"}      = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);
# }
#}






sub make_target_steel
{
 my $NUM  = 2;
 my @x    = (-5.55,5.55);
 my @Rot  = (-90,90);
 my @Rin  = (5.9,5.9,5.9);
 my @Rout = (7.9,7.9,7.82);#25deg chamfer
 my @zPln = (0,1.91,1.95);#25deg chamfer 4mm*8mm
 my @name = ("left_steel","right_steel"); 
 my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Fe","G4_Fe");
 #my @mat  = ("G4_STAINLESS-STEEL","G4_STAINLESS-STEEL");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg $Rot[$n-1]*deg 0*deg";
    $detector{"color"}      = "00BFFF";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0*deg 360*deg 3*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm";
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
#sub make_target_steel
#{
# my $NUM  = 2;
# my @x    = (-5.55,5.55);
# my @Rin  = (5.9,5.9);
# my @Rout = (7.9,7.9);
# my @Dz   = (1.95,1.95);
# my @name = ("left_steel","right_steel"); 
# my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in"); 
# my @mat  = ("G4_Fe","G4_Fe");
# #my @mat  = ("G4_STAINLESS-STEEL","G4_STAINLESS-STEEL");
#
# for(my $n=1; $n<=$NUM; $n++)
# {
##     my $pnumber     = cnumber($n-1, 10);
#     my %detector=init_det();
#    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
#    $detector{"mother"}      = "$mother[$n-1]" ;
#    $detector{"description"} = "$DetectorName\_$name[$n-1]";
#    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
#    $detector{"rotation"}   = "0*deg -90*deg 0*deg";
#    $detector{"color"}      = "00BFFF";
#    $detector{"type"}       = "Tube";
#    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
#    $detector{"material"}   = $mat[$n-1];
#    $detector{"mfield"}     = "no";
#    $detector{"ncopy"}      = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);
# }
#}

sub make_target_coil
{
 my $NUM  = 8;
 my @x    = (9.5,7.15,5.1,2.95,2.95,5.1,7.15,9.5);
 #my @x    = (-13.5,-11.15,-9.1,-6.95,6.95,9.1,11.15,13.5);
 my @rot  = (0,0,0,0,0,0,0,0);
 #my @rot  = (90,90,90,90,-90,-90,-90,-90);
 my @Rin  = (19.9,16.1,13,9.9,9.9,13,16.1,19.9);
 my @Rout = (22.9,18.6,15,11.4,11.4,15,18.6,22.9);
 my @Dz   = (2.3,1.85,1.5,1.15,1.15,1.5,1.85,2.3);
 my @name = ("left_coil_xl","left_coil_l","left_coil_m","left_coil_s","right_coil_s","right_coil_m","right_coil_l","right_coild_xl"); 
 #my @mother = ("$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box"); 
 my @mother = ("$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2"); 
 #my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    #$detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"pos"}        = "0*cm 0*cm $x[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
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

#somehow the addcutwayplane in solid_slice.vis would change the position and rotation of the coils. This is the origin parameters without addcutwayplane option.
sub make_target_coil_1
{
 my $NUM  = 8;
 my @x    = (-13.5,-11.15,-9.1,-6.95,6.95,9.1,11.15,13.5);
 #my @rot  = (0,0,0,0,-180,-180,-180,-180);
 my @rot  = (90,90,90,90,-90,-90,-90,-90);
 my @Rin  = (19.9,16.1,13,9.9,9.9,13,16.1,19.9);
 my @Rout = (22.9,18.6,15,11.4,11.4,15,18.6,22.9);
 my @Dz   = (2.3,1.85,1.5,1.15,1.15,1.5,1.85,2.3);
 my @name = ("left_coil_xl","left_coil_l","left_coil_m","left_coil_s","right_coil_s","right_coil_m","right_coil_l","right_coild_xl"); 
 #my @mother = ("$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_right_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box","$DetectorName\_left_coil_box"); 
 my @mother = ("$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_left_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2","$DetectorName\_right_coil_box_2"); 
 #my @mother = ("$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in","$DetectorName\_SC_in"); 
 my @mat  = ("G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu","G4_Cu");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
     my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg $rot[$n-1]*deg 0*deg";
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
 my @x    = (-4,4,-4,4);
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
    $detector{"pos"}        = "$x[$n-1]*cm 0*cm 0*cm";
    $detector{"dimensions"} = "0*deg 360*deg 2*counts $Rin[4]*cm $Rin[5]*cm $Rout[4]*cm $Rout[5]*cm $zPln[4]*cm $zPln[5]*cm";
    
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
 my @x    = (-4,4,-4,4);
 #moved 30K shield back to center
 #my @z    = (-5,-5,-5,-5);#need to move it back to center due to the 30K shield position
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
    $detector{"dimensions"} = "0*deg 360*deg 13*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm";
    #$detector{"dimensions"} = "0*deg 360*deg 9*counts $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rin[10]*cm $Rin[11]*cm $Rin[12]*cm $Rin[13]*cm $Rin[14]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $Rout[10]*cm $Rout[11]*cm $Rout[12]*cm $Rout[13]*cm $Rout[14]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm $zPln[10]*cm $zPln[11]*cm $zPln[12]*cm $zPln[13]*cm $zPln[14]*cm";
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
 my @x    = (-8.39751,8.39751);
# my @numZPlane = 15;
 my @Rot  = (90,-90);
 my @Rin  = (8.75,  8.75,   8.825,  9.954, 9.954, 12.736,12.736,16.068,16.068,17.57, 17.75); 
 my @Rout = (11.685,11.685, 11.685, 11.685,15.212,15.212,18.851,18.851,25,    25,    24.8);
 my @zPln = (0,     1.384,  1.384,  2.494,2.494,4.887,4.887,7.752,7.752,8.852,9.052);
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

#The out surface is tube,need to be cutted more, but the square cut is too complicated
sub make_magnet_support_1
{
  #for the tube 
  my $NUM  = 1;
  my @Rin  = (47.6/2,47.6/2,23.22839568,23.22839568,39/2,39/2,23.22839568,23.22839568,47.6/2,47.6/2);
  my @Rout = (25,25,25,24,24,24,24,25,25,25);
  my @zPln = (0,1.05957,2.5,2.5,(24.31549-18.18608)/2,(24.31549-18.18608)/2+18.18608,24.31549-2.5,24.31549-2.5,24.31549-1.05957,24.31549);
  my @x    = (24.31549/2);
  my @name = ("magnet_support");
  my @mother = ("$DetectorName\_SC_in"); 
  my @mat  = ("G4_Al");
  my @rot  = (90);
  #my @color = ("ff0000","ff0000");
  my @color = ("FF6600");

  #for the cone holes
  my $N      = 4;
  my $R_out1 = 18.3592/2;
  my $R_out2 = 25*0.52159;
  my $pDz    = 7.71051/2;
  my @cone_x = (0,22.25,0,-22.25);
  my @cone_y = (22.25,0,-22.25,0);
  my @cone_rot = (90,0,90,0);

  my %detector=init_det(); 
  $detector{"name"}        = "$DetectorName\_polycone";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The tube";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Polycone";
  $detector{"dimensions"} = "0*deg 360*deg 10*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm"; 
  $detector{"material"}   = $mat[0];
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);

  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cone_1";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The holes cones";
  $detector{"pos"}        = "0*cm $cone_y[0]*cm $cone_x[0]*cm";
  $detector{"rotation"}   = "$cone_rot[0]*deg 0*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Cons";
  $detector{"dimensions"} = "0*cm $R_out1*cm 0*cm $R_out2*cm $pDz*cm 0*deg 360*deg"; 
  $detector{"material"}   = $mat[0];
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);

  #Operation@ indicates that the postion and rotation depends on the first volume( volume before the operation)
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cuts_1";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The holes";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"pos"}        = "0*cm 22.169365*cm 0*cm";
  #$detector{"rotation"}   = "90*deg 0*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Operation:@ $DetectorName\_polycone * $DetectorName\_cone_1";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = $mat[0];
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);
  
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cone_2";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The holes cones";
  $detector{"pos"}        = "0*cm $cone_y[1]*cm $cone_x[1]*cm";
  $detector{"rotation"}   = "$cone_rot[1]*deg 0*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Cons";
  $detector{"dimensions"} = "0*cm $R_out1*cm 0*cm $R_out2*cm $pDz*cm 0*deg 360*deg"; 
  $detector{"material"}   = $mat[0];
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cuts_2";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The holes";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"pos"}        = "0*cm 22.169365*cm 0*cm";
  #$detector{"rotation"}   = "90*deg 0*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Operation:@ $DetectorName\_polycone * $DetectorName\_cone_2";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = $mat[0];
  #$detector{"mfield"}     = "no";
  #$detector{"ncopy"}      = 1;
  #$detector{"pMany"}       = 1;
  $detector{"exist"}       = 0;
  #$detector{"visible"}     = 1;
  #$detector{"style"}       = 1;
  print_det(\%configuration, \%detector);

   %detector=init_det();
   $detector{"name"}        = "$DetectorName\_magnet_support";
   $detector{"mother"}      = "$mother[0]" ;
   $detector{"description"} = "The magnet support";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
   #$detector{"pos"}        = "0*cm 22.169365*cm 0*cm";
   #$detector{"rotation"}   = "90*deg 0*deg 0*deg";
   $detector{"color"}      = "FF6600";    
   $detector{"type"}       = "Operation:@ $DetectorName\_polycone - $DetectorName\_cuts_1";
   $detector{"dimensions"} = "0"; 
   $detector{"material"}   = $mat[0];
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

sub make_magnet_support
{
  #for the tube 
  my $NUM  = 1;
  my @Rin  = (47.6/2,47.6/2,23.22839568,23.22839568,39/2,39/2,23.22839568,23.22839568,47.6/2,47.6/2);
  my @Rout = (25,25,25,24,24,24,24,25,25,25);
  my @zPln = (0,1.05957,2.5,2.5,(24.31549-18.18608)/2,(24.31549-18.18608)/2+18.18608,24.31549-2.5,24.31549-2.5,24.31549-1.05957,24.31549);
  my @x    = (24.31549/2);
  my @name = ("magnet_support");
  my @mother = ("$DetectorName\_SC_in"); 
  my @mat  = ("G4_Al");
  my @rot  = (90);
  #my @color = ("ff0000","ff0000");
  my @color = ("FF6600");

  ##for the cone holes, from measuring
  #my $N      = 4;
  #my $R_out1 = 18.3592/2;
  #my $R_out2 = 25*0.52159;
  #my $pDz    = 7.71051/2;
  #my @cone_x = (0,21.1448,0,-21.1448);
  #my @cone_y = (21.1448,0,-21.1448,0);
  #my @cone_rot = (90,0,-90,180);
  
  #For the cone holes, assumming 25deg
  my $DEG    = 3.1415926/180;
  my $N      = 4;
  my $R_out1 = 19.5*sin(25*$DEG);#195*sin25deg
  my $R_out2 = 25*tan(25*$DEG);#25*tan25deg
  my $pDz    = (25-19.5*cos(25*$DEG))/2;
  my @cone_x = (0,19.5*cos(25*$DEG)+$pDz,0,-(19.5*cos(25*$DEG)+$pDz));
  my @cone_y = (19.5*cos(25*$DEG)+$pDz,0,-(19.5*cos(25*$DEG)+$pDz),0);
  my @cone_rot = (90,0,-90,180);

  #For the additional hole cuts
  my @rec_rot    = (0,90,0,-90);
  my @rec_z      = (0,20.58,0,-20.58);
  my @rec_y      = (20.58,0,-20.58,0);
  my @addtub_rot = (90,90,90,90);
  my @addtub_rot2 = (0,90,180,270);
  my @addtub_z   = (0,21.66,0,-21.66);
  my @addtub_y   = (21.66,0,-21.66,0);

  my %detector=init_det(); 
  $detector{"name"}        = "$DetectorName\_polycone";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The tube";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "808080";    
  $detector{"type"}       = "Polycone";
  $detector{"dimensions"} = "0*deg 360*deg 10*counts $Rin[0]*cm $Rin[1]*cm $Rin[2]*cm $Rin[3]*cm $Rin[4]*cm $Rin[5]*cm $Rin[6]*cm $Rin[7]*cm $Rin[8]*cm $Rin[9]*cm $Rout[0]*cm $Rout[1]*cm $Rout[2]*cm $Rout[3]*cm $Rout[4]*cm $Rout[5]*cm $Rout[6]*cm $Rout[7]*cm $Rout[8]*cm $Rout[9]*cm $zPln[0]*cm $zPln[1]*cm $zPln[2]*cm $zPln[3]*cm $zPln[4]*cm $zPln[5]*cm $zPln[6]*cm $zPln[7]*cm $zPln[8]*cm $zPln[9]*cm"; 
  $detector{"material"}   = "Component";
  print_det(\%configuration, \%detector);


  for(my $i=1; $i<=$N; $i++){
    #The holes
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_cone_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes cones";
    $detector{"pos"}        = "0*cm $cone_y[$i-1]*cm $cone_x[$i-1]*cm";
    $detector{"rotation"}   = "$cone_rot[$i-1]*deg 0*deg 0*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Cons";
    $detector{"dimensions"} = "0*cm $R_out1*cm 0*cm $R_out2*cm $pDz*cm 0*deg 360*deg"; 
    $detector{"material"}   = "Component";
    $detector{"exist"}       = 0;
    print_det(\%configuration, \%detector);
    
    #The cuts addition to the holes
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_rec_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes addition cuts";
    $detector{"pos"}        = "0*cm $rec_y[$i-1]*cm $rec_z[$i-1]*cm";
    $detector{"rotation"}   = "$rec_rot[$i-1]*deg 0*deg 0*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Box";
    $detector{"dimensions"} = "12.157*cm 1.08*cm 2.5*cm"; 
    $detector{"material"}   = "Component";
    $detector{"mfield"}     = "no";
    print_det(\%configuration, \%detector);
    #The cuts addition to the holes
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_addtub_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes addition cuts";
    $detector{"pos"}        = "0*cm $addtub_y[$i-1]*cm $addtub_z[$i-1]*cm";
    $detector{"rotation"}   = "0*deg $addtub_rot[$i-1]*deg  $addtub_rot2[$i-1]*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "0*cm 2.5*cm 12.157*cm 0*deg 180*deg"; 
    $detector{"material"}   = "Component";
    print_det(\%configuration, \%detector);
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_rectub_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes addition cuts";
    $detector{"pos"}        = "0*cm $rec_y[$i-1]*cm $rec_z[$i-1]*cm";
    $detector{"rotation"}   = "$rec_rot[$i-1]*deg 0*deg  0*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Operation:@ $DetectorName\_rec_$i + $DetectorName\_addtub_$i";
    $detector{"dimensions"} = "0"; 
    $detector{"material"}   = "Component";
    print_det(\%configuration, \%detector);
    
    #combine the rectangular tube and cone together, they corresponds to each other, lucky!
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_holes_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes all";
    $detector{"pos"}        = "0*cm $rec_y[$i-1]*cm $rec_z[$i-1]*cm";
    $detector{"rotation"}   = "$rec_rot[$i-1]*deg 0*deg  0*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Operation:@ $DetectorName\_rectub_$i + $DetectorName\_cone_$i";
    $detector{"dimensions"} = "0"; 
    $detector{"material"}   = "Component";
    #$detector{"material"}   = $mat[0];
    #$detector{"mfield"}     = "no";
    #$detector{"ncopy"}      = 1;
    #$detector{"pMany"}       = 1;
    #$detector{"exist"}       = 1;
    #$detector{"visible"}     = 1;
    #$detector{"style"}       = 1;
    #$detector{"sensitivity"} = "no";
    #$detector{"hit_type"}    = "no";
    #$detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

    ##Operation@ indicates that the postion and rotation depends on the first volume( volume before the operation)
    #%detector=init_det();
    #$detector{"name"}        = "$DetectorName\_cuts_$i";
    #$detector{"mother"}      = "$mother[0]" ;
    #$detector{"description"} = "The holes";
    #$detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
    #$detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
    ##$detector{"pos"}        = "0*cm 22.169365*cm 0*cm";
    ##$detector{"rotation"}   = "90*deg 0*deg 0*deg";
    #$detector{"color"}      = "808080";    
    #$detector{"type"}       = "Operation:@ $DetectorName\_polycone * $DetectorName\_holes_$i";
    #$detector{"dimensions"} = "0"; 
    #$detector{"material"}   = $mat[0];
    #$detector{"mfield"}     = "no";
    #$detector{"ncopy"}      = 1;
    #$detector{"pMany"}       = 1;
    #$detector{"exist"}       = 1;
    #$detector{"visible"}     = 1;
    #$detector{"style"}       = 1;
    #$detector{"sensitivity"} = "no";
    #$detector{"hit_type"}    = "no";
    #$detector{"identifiers"} = "no";
    #print_det(\%configuration, \%detector);
    #Operation@ indicates that the postion and rotation depends on the first volume( volume before the operation)
    %detector=init_det();
    $detector{"name"}        = "$DetectorName\_cuts_$i";
    $detector{"mother"}      = "$mother[0]" ;
    $detector{"description"} = "The holes";
    $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
    #$detector{"pos"}        = "0*cm 22.169365*cm 0*cm";
    #$detector{"rotation"}   = "90*deg 0*deg 0*deg";
    $detector{"color"}      = "808080";    
    $detector{"type"}       = "Operation:@ $DetectorName\_polycone * $DetectorName\_holes_$i";
    $detector{"dimensions"} = "0"; 
    $detector{"material"}   = "Component";
    $detector{"mfield"}     = "no";
    print_det(\%configuration, \%detector);
    #$detector{"material"}   = $mat[0];
    #$detector{"mfield"}     = "no";
    #$detector{"ncopy"}      = 1;
    #$detector{"pMany"}       = 1;
    #$detector{"exist"}       = 1;
    #$detector{"visible"}     = 1;
    #$detector{"style"}       = 1;
    #$detector{"sensitivity"} = "no";
    #$detector{"hit_type"}    = "no";
    #$detector{"identifiers"} = "no";
    #print_det(\%configuration, \%detector);

  }
  #%detector=init_det();
  #$detector{"name"}        = "$DetectorName\_magnet_support_1";
  #$detector{"mother"}      = "$mother[0]" ;
  #$detector{"description"} = "The magnet support";
  #$detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  #$detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"color"}      = "FF6600";    
  #$detector{"type"}       = "Operation:@ $DetectorName\_polycone - $DetectorName\_cuts_1";
  #$detector{"dimensions"} = "0"; 
  #$detector{"material"}   = "Component";
  #$detector{"mfield"}     = "no";
  #$detector{"exist"}       = 0;
  #print_det(\%configuration, \%detector);
  #%detector=init_det();
  #$detector{"name"}        = "$DetectorName\_magnet_support_2";
  #$detector{"mother"}      = "$mother[0]" ;
  #$detector{"description"} = "The magnet support";
  #$detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  #$detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"color"}      = "FF6600";    
  #$detector{"type"}       = "Operation:@ $DetectorName\_magnet_support_1 - $DetectorName\_cuts_2";
  #$detector{"dimensions"} = "0"; 
  #$detector{"material"}   = "Component";
  #$detector{"mfield"}     = "no";
  #$detector{"exist"}       = 0;
  #print_det(\%configuration, \%detector);
  #%detector=init_det();
  #$detector{"name"}        = "$DetectorName\_magnet_support_3";
  #$detector{"mother"}      = "$mother[0]" ;
  #$detector{"description"} = "The magnet support";
  #$detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  #$detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"color"}      = "FF6600";    
  #$detector{"type"}       = "Operation:@ $DetectorName\_magnet_support_2 - $DetectorName\_cuts_3";
  #$detector{"dimensions"} = "0"; 
  #$detector{"material"}   = "Component";
  #$detector{"mfield"}     = "no";
  #$detector{"exist"}       = 0;
  #print_det(\%configuration, \%detector);
  #%detector=init_det();
  #$detector{"name"}        = "$DetectorName\_magnet_support";
  #$detector{"mother"}      = "$mother[0]" ;
  #$detector{"description"} = "The magnet support";
  #$detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  #$detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  #$detector{"color"}      = "FF6600";    
  #$detector{"type"}       = "Operation:@ $DetectorName\_magnet_support_3 - $DetectorName\_cuts_4";
  #$detector{"dimensions"} = "0"; 
  #$detector{"material"}   = $mat[0];
  #$detector{"mfield"}     = "no";
  #$detector{"ncopy"}      = 1;
  #$detector{"pMany"}       = 1;
  #$detector{"exist"}       = 1;
  #$detector{"visible"}     = 1;
  #$detector{"style"}       = 1;
  #$detector{"sensitivity"} = "no";
  #$detector{"hit_type"}    = "no";
  #$detector{"identifiers"} = "no";
  #print_det(\%configuration, \%detector);
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cuts_12";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The magnet support";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "FF6600";    
  $detector{"type"}       = "Operation:@ $DetectorName\_cuts_1+ $DetectorName\_cuts_2";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = "Component";
  $detector{"mfield"}     = "no";
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cuts_123";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The magnet support";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "FF6600";    
  $detector{"type"}       = "Operation:@ $DetectorName\_cuts_12 + $DetectorName\_cuts_3";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = "Component";
  $detector{"mfield"}     = "no";
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_cuts_1234";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The magnet support";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "FF6600";    
  $detector{"type"}       = "Operation:@ $DetectorName\_cuts_123 + $DetectorName\_cuts_4";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = "Component";
  $detector{"mfield"}     = "no";
  $detector{"exist"}       = 0;
  print_det(\%configuration, \%detector);
  %detector=init_det();
  $detector{"name"}        = "$DetectorName\_magnet_support";
  $detector{"mother"}      = "$mother[0]" ;
  $detector{"description"} = "The magnet support";
  $detector{"pos"}        = "$x[0]*cm 0*cm 0*cm";
  $detector{"rotation"}   = "0*deg $rot[0]*deg 0*deg";
  $detector{"color"}      = "FF6600";    
  $detector{"type"}       = "Operation:@ $DetectorName\_polycone - $DetectorName\_cuts_1234";
  $detector{"dimensions"} = "0"; 
  $detector{"material"}   = $mat[0];
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
