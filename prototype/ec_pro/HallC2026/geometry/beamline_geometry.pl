use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'beamline';

my $DetectorMother="root";

# my $z_target			= $parameters{"z_target"};

#He3 target
my $target_length = 40;

#Carbon target 2.21 g/cm3 density
my $target_length_C = 0.04525;

# carbon optics target run
# my @y    = (0,0,(1.125/2-0.688/2)*2.54,(1.125/2-0.688/2)*2.54);  #beam pass through pos1 foil center, and pos4 and pos6

# carbon hole target run
my @y    = (0,-(0.906-0.688/2)*2.54,-(0.906-1.125/2)*2.54,-(0.906-1.125/2)*2.54);  #beam pass through pos4 and pos6 holes

# carbon foil target run
# my @y    = ???  #beam pass through pos4 and pos6 top
my $material_window = "Al_7075";
#my $hx1	= 2*2.54/2;
my $hx1	= 6*2.54/2;
my $hy1	= 4*2.54/2;
#my $hz1	= 130-30*2.54+41*2.54+15*2.54+36*2.54; //to match 80 in long
my $hz1	= 2180+10;
my $hz2	= $hz1+45*2.54+35*2.54 ;
#30in first try
#my $hdz= 30*2.54;
#my $hdz1= 80*2.54;
my $hdz1= 18*2.54/2.;
my $hdz2= 35*2.54;
# my $z14 = 6;
my $z14 = 20;
sub beamline
{
#make_beam_entrance_cryo();
#make_beam_exit_cryo();
#make_beam_coli();
#make_beam_coli2();
#make_beam_dummy();
# make_beam_LH2_CH();
make_target();
# make_beam_LH2();
# make_beam_LD2();
#make_beam_carbon();
# make_beam_LH2_CH_con();
# make_beam_entrance();
# make_beam_exit();
make_beam_exit_long();
# make_virutal();

}
sub make_beam_entrance_cryo
{
 my @z    = (-88.52,0.0);
 my @Rin  = (0.,0.);
 my @Rout = (2.75/2*2.54,2.067/2*2.54);
 my @Dz   = (72.9/2,72.9/2);
 
 my $NUM  = 2;
 my @name = ("BMP1","BMV1"); 
 my @mother = ("$DetectorMother","$DetectorName\_BMP1"); 
 my @mat  = ("G4_Al","G4_Galactic");
 my @color = ("ff00ff","808080");

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
sub make_beam_coli2
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli2";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "-10.5*cm 0*cm $hz2*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"color"}       = "ff00ff";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx1*cm $hy1*cm $hdz2*cm";    
#	$detector{"material"}    = "G4_Al";
	$detector{"material"}    = "G4_Pb";
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "flux";
	$detector{"hit_type"}    = "flux";
        my $ID = 21;
        $detector{"identifiers"} = "id manual $ID";
        print_det(\%configuration, \%detector);
}
sub make_beam_coli
{
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_coli";
	$detector{"mother"}      = "$DetectorMother";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}         = "-240.5*cm 0*cm $hz1*cm";
	$detector{"rotation"}    = "0*deg 7*deg 0*deg";
	$detector{"color"}       = "1a4fff";
        $detector{"type"}       = "Box";
        $detector{"dimensions"} = "$hx1*cm $hy1*cm $hdz1*cm";    
#	$detector{"material"}    = "G4_Al";
	$detector{"material"}    = "G4_Pb";
	$detector{"mfield"}      = "no";
	$detector{"ncopy"}       = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "flux";
	$detector{"hit_type"}    = "flux";
        my $ID = 20;
        $detector{"identifiers"} = "id manual $ID";
        print_det(\%configuration, \%detector);
}
sub make_beam_exit_cryo
{
 my @z    = (65.025,0.0,161.2,0.,359.4,0.,632.35,0.,838.05,0.,980.9,0.);
 my @x    = (-0.55,0.,-1.17,0.,-1.95,0.,-3.24,0.,-3.6,0.,-2.81,0.);
 my @y    = (-0.28,0.,-0.13,0.,-0.38,0.,-0.3,0.,0.02,0.,0.12,0.);
 my @Rin  = (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
 my @Rout = (7/2,2.067/2*2.54,5.1/2,(5.1-1.75)/2,12.7/2,(12.7-1.75)/2,12.7/2,(12.7-1.75)/2,27.3/2,(27.3-1.75)/2,27.3/2,(27.3-1.75)/2);
 my @Dz   = (15.75/2,15.75/2,(249.5-72.9)/2,(249.5-72.9)/2,(469.3-249.5)/2,(469.3-249.5)/2,(795.4-469.3)/2,(795.4-469.3)/2,(880.7-795.4)/2,(880.7-795.4)/2,(1081.1-880.7)/2,(1081.1-880.7)/2); 
 #my @Dz   = (15.75/2,15.75/2,88.3,88.3,109.9,109.9,163.05,163.05,42.65,42.65); 
 my $NUM  = 12;
 my @name = ("BMP2","BMV2","BMD2","BMW2","BMC2","BMF2","BMC3","BMF3","BMC4","BMF4","BMC5","BMF5"); 
 my @mother = ("$DetectorMother","$DetectorName\_BMP2","$DetectorMother","$DetectorName\_BMD2","$DetectorMother","$DetectorName\_BMC2","$DetectorMother","$DetectorName\_BMC3","$DetectorMother","$DetectorName\_BMC4","$DetectorMother","$DetectorName\_BMC5"); 
 my @mat  = ("G4_Al","G4_Galactic","G4_Al","G4_Galactic","G4_Al","G4_Galactic","G4_Al","G4_Galactic","G4_Al","G4_Galactic","G4_Al","G4_Galactic");
 my @color = ("ff00ff","808080","ff00ff","808080","00FFFF","808080","ff00ff","808080","00FFFF","808080","ff00ff","808080");

 for(my $n=1; $n<=$NUM; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "$x[$n-1]*cm $y[$n-1]*cm $z[$n-1]*cm";
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

sub make_beam_entrance
{
# upstream Be window .01thk, 12.80 inch from target center at pivot,and 0.0002" thk Al cover, inner diameter 2.067" outter diameter 2.375", length 200cm
 my @z    = (-12.80*2.54-100,0,-100+0.01*2.54/2,100-0.0002*2.54-0.01*2.54/2,100-0.0002*2.54/2);
 my @Rin  = (0.,0.,0.,0.,0);
 my @Rout = (2.375/2*2.54,2.067/2*2.54,2.067/2*2.54,2.067/2*2.54,2.067/2*2.54);
 my @Dz   = (100,100,0.01*2.54/2,0.01*2.54/2,0.0002*2.54/2);
 
 my $NUM  = 5;
 my @name = ("BMP1","BMV1","BMD1","BMW1","BMC1"); 
 my @mother = ("$DetectorMother","$DetectorName\_BMP1","$DetectorName\_BMV1","$DetectorName\_BMV1","$DetectorName\_BMV1"); 
 my @mat  = ("G4_Al","G4_Galactic","G4_Galactic","G4_Be","G4_Al");
 my @color = ("0000ff","808080","808080","00FFFF","FF00FF");

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

sub make_beam_exit
{
# downstream Be window .02thk, 17.25 inch from target center at pivot,and 0.001" thk Al cover, inner diameter 2.067" outter diameter 2.375"
 my @z    = (17.25*2.54+100,0,-100+0.001*2.54/2,-100+0.001*2.54+0.02*2.54/2,100-0.02*2.54/2,);
 my @Rin  = (0.,0.,0.,0.,0);
 my @Rout = (2.375/2*2.54,2.067/2*2.54,2.067/2*2.54,2.067/2*2.54,2.067/2*2.54);
 my @Dz   = (100,100,0.001*2.54/2,0.02*2.54/2,0.02*2.54/2,);

 my $NUM  = 5; 
 my @name = ("B3PP","B3PV","B3C1","B3W1","B3DM"); 
 my @mother=("$DetectorMother","$DetectorName\_B3PP","$DetectorName\_B3PV","$DetectorName\_B3PV","$DetectorName\_B3PV");
 my @mat  = ("G4_Al","G4_Galactic","G4_Al","G4_Be","G4_Galactic");
 my @color = ("0000ff","808080","FF00FF","00FFFF","808080");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
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

sub make_beam_exit_long
{
 my $NUM  = 4;
 my @z    = (10/2+1500,0,-1500+0.02*2.54/2,1500-0.01);
 my @Rin = (0,0,0,0);
 my @Rout1 = (2.375/2*2.54,2.067/2*2.54,2.067/2*2.54,10.067/2*2.54);
 my @Rout2 = (10.375/2*2.54,10.067/2*2.54,2.067/2*2.54,10.067/2*2.54);
 my @Dz   = (1500,1500,0.02*2.54/2,0.01);
 my @name = ("exit_long_outer","exit_long_inner","exit_long_front","exit_long_end");
 my @mother=("$DetectorMother","$DetectorName\_exit_long_outer","$DetectorName\_exit_long_inner","$DetectorName\_exit_long_inner");
 my @mat  = ("G4_Al","G4_Galactic","G4_Be","Kryptonite");
 my @color = ("ff00ff","808080","ff0000","ffffff");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Cons";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout1[$n-1]*cm $Rin[$n-1]*cm $Rout2[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg";
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
sub make_beam_LH2_CH
{ # Target windows 2024-T3 aluminum 0.02 inches thick. The Chamber has an inner and outer radisus of 41 and 45 inches respectively.        
 my $NUM  = 5;
 my @z    = (0,0,0,0,0);
# my @z    = (0,0,0,0);
 my @Rin  = (41/2*2.54,0,(45./2-0.02)*2.54,41/2*2.54,(45./2-0.02)*2.54);
 my @phi  = (-90,-90,0,0,0);
 #my @phi  = (-90,0,-90,0);
 my @thetaT  = (360.,360.,113.5,113.5,92.1);
 my @thetaS  = (0,0,66.5,66.5,87.89);
 my @Rout = (45./2*2.54,41/2*2.54,45./2*2.54,(45./2-0.02)*2.54,45./2*2.54);
 my @Dz   = (44.75/2*2.54,44.75/2*2.54,17.0/2*2.54,17.0/2*2.54,1.66/2*2.54);
 my @name = ("$DetectorName\_TACHR","$DetectorName\_TACHI","$DetectorName\_TACWI","$DetectorName\_TACWI2","$DetectorName\_TACWIH");
 my @mother=("$DetectorMother","$DetectorMother","$DetectorName\_TACHR","$DetectorName\_TACHR","$DetectorName\_TACWI");
 #my @mother=("$DetectorMother","$DetectorName\_TACHR","$DetectorMother","$DetectorName\_TACWI");
 my @mat  = ("G4_Al","G4_Galactic","Al_2024","G4_Galactic","G4_Galactic");
 my @color  = ("808080","ff00ff","808080","ff00ff","ff00ff");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "$phi[$n-1]*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm $thetaS[$n-1]*deg $thetaT[$n-1]*deg";
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
sub make_beam_LH2
{ #cro target 10cm loop, upstream (0.05'') and downstream (0.032'') target foils with Al_7075.For loop2--LH2 target, the entrance Al_7075 window thickness is 0.15mm and the exit window thickness is 0.19 mm. For loop 3 -LD2 target, the entrance window thickness is 0.116mm and the exit window thickness is 0.184mm.       
 my $NUM  = 4;
 my @z    = (0,-5-0.015,5+0.019,0.);
 my @phi  = (-90,0,0,0);
 my @Rin  = (0.,0.,0.,0.);
 my @Rout = (3.81/2,3.81/2-0.0278,3.81/2-0.0278,3.81/2-0.0278);
 my @Dz   = (5.0+0.15,0.015/2,0.019/2,10.0/2);
 my @name = ("$DetectorName\_TACH","$DetectorName\_TAEN","$DetectorName\_TAEX","$DetectorName\_TALH");
 my @mother=("$DetectorName\_TACHI","$DetectorName\_TACH","$DetectorName\_TACH","$DetectorName\_TACH");
 my @mat  = ("G4_Galactic","Al_7075","Al_7075","G4_lH2");
 my @color  = ("101011","00FFFF","00FFFF","808080");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "$phi[$n-1]*deg 0*deg 0*deg";
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
sub make_beam_carbon
{
# Carbon foil 100mg/cm2
 my $coolgas_z=(-12.80*2.54+17.25*2.54)/2;
 my $NUM  = 1;
 my @z    = (0);
 my @Dx = (0.625/2*2.54);
 my @Dy = (0.688/2*2.54);
 my @Dz   = ($target_length_C/2);
 my @name = ("target_carbon1"); 
 #my @mother = ("$DetectorMother");
 my @mother = ("$DetectorName\_TACHI");
 my @mat  = ("G4_GRAPHITE");
 my @color = ("101011");
 
 for(my $n=1; $n<=$NUM; $n++)
 {
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = "$DetectorName\_$name[$n-1]";
    $detector{"pos"}        = "0*cm $y[$n-1]*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "90*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Box";
    $detector{"dimensions"} = "$Dx[$n-1]*cm $Dy[$n-1]*cm $Dz[$n-1]*cm";
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
sub make_beam_LD2
{ #cro target 10cm loop, upstream (0.05'') and downstream (0.032'') target foils with Al_7075.For loop2--LH2 target, the entrance Al_7075 window thickness is 0.15mm and the exit window thickness is 0.19 mm. For loop 3 -LD2 target, the entrance window thickness is 0.116mm and the exit window thickness is 0.184mm.       
 my $NUM  = 4;
 my @z    = (0,-5-0.0116,5+0.0184,0.);
 my @phi  = (-90,0,0,0);
 my @Rin  = (0.,0.,0.,0.);
 my @Rout = (3.81/2,3.81/2-0.0278,3.81/2-0.0278,3.81/2-0.0278);
 my @Dz   = (5.0+0.15,0.0116/2,0.0184/2,10.0/2);
 my @name = ("$DetectorName\_TACH","$DetectorName\_TAEN","$DetectorName\_TAEX","$DetectorName\_TALH");
 my @mother=("$DetectorName\_TACHI","$DetectorName\_TACH","$DetectorName\_TACH","$DetectorName\_TACH");
 my @mat  = ("G4_Galactic","Al_7075","Al_7075","LD2");
 my @color  = ("101011","00FFFF","00FFFF","808080");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "$phi[$n-1]*deg 0*deg 0*deg";
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

sub make_beam_LH2_CH_con
{ # connector between target chamber and the front end of the beam pipe.       
 my $NUM  = 2;
 my @z    = (45./2*2.54+20.5/2,0);
 my @y    = (0,0);
 my @phi  = (0,0);
 my @Rin  = (0,0);
 my @Rout = (1.66/2*2.54,1.38/2*2.54);
 my @Dz   = (20.5/2,20.5/2);
 my @name = ("$DetectorName\_TACCN","$DetectorName\_TACCNV");
 my @mother=("$DetectorMother","$DetectorName\_TACCN");
 my @mat  = ("G4_Al","G4_Galactic");
 my @color  = ("101011","ff00ff");

 for(my $n=1; $n<=$NUM; $n++)
 {
#     my $pnumber     = cnumber($n-1, 10);
    my %detector=init_det();
    $detector{"name"}        = "$name[$n-1]";
    $detector{"mother"}      = "$mother[$n-1]" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm $y[$n-1]*cm $z[$n-1]*cm";
    $detector{"rotation"}   = "$phi[$n-1]*deg 0*deg 0*deg";
    $detector{"color"}      = $color[$n-1];
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "$Rin[$n-1]*cm $Rout[$n-1]*cm $Dz[$n-1]*cm 0*deg 360*deg ";
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
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_target";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "ff0000";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"} = "0*cm 1.88*cm 5*cm 0*deg 360*deg";
#     $detector{"material"}   = "G4_lH2";
     $detector{"material"}   = "LD2";
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
sub make_beam_dummy
{ #cro target 10cm loop, upstream (0.05'') and downstream (0.032'') target foils with Al_7075.For loop2--LH2 target, the entrance Al_7075 window thickness is 0.15mm and the exit window thickness is 0.19 mm. For loop 3 -LD2 target, the entrance window thickness is 0.116mm and the exit window thickness is 0.184mm.       
 my $NUM  = 4;
 my @z    = (0,-5-0.0861,5+0.109,0.);
 my @Rin  = (0.,0.,0.,0.);
 my @Rout = (3.81/2,3.81/2,3.81/2,3.81/2-0.0278);
 my @Dz   = (5.0+0.15,0.05/2*2.54,0.032/2*2.54,10.0/2);
 my @name = ("$DetectorName\_TACH","$DetectorName\_TAW1","$DetectorName\_TAW2","$DetectorName\_TALH");
 my @mother=("$DetectorMother","$DetectorName\_TACH","$DetectorName\_TACH","$DetectorName\_TACH");
 my @mat  = ("G4_Galactic","Al_7075","Al_7075","G4_lH2");
 my @color  = ("101011","FF00FF","FF0000","808080");

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

sub make_virutal 
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_virtualplane";
 $detector{"mother"}      = "$DetectorName\_TACHI";
 $detector{"description"} = $detector{"name"};
#  $detector{"pos"}         = "0*cm 0*cm $z14*cm";
 $detector{"pos"}         = "0*cm $z14*cm 0*cm";
 $detector{"rotation"}    = "-90*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Tube";
 $detector{"dimensions"} = "0.1*cm 30*cm 0.0001*cm 0*deg 360*deg";
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

beamline();
1;
