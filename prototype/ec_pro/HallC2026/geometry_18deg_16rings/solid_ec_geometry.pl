# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_ec';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";

my $ec_sen="solid_ec";
my $ec_ps_sen="solid_ec_ps";

my $color_abs="0000ff";
my $color_scint="f1ffe1";
my $color_gap="80ff77";
my $color_paint="77ff80";
my $color_support="123456";
my $color_wrap_sectorodd="00ffff";
my $color_wrap_sectoreven="ffff00";
my $color_rod="654321";

my $material_abs="G4_Pb";
my $material_scint="G4_POLYSTYRENE";
# my $material_scint="G4_PLASTIC_SC_VINYLTOLUENE";
my $material_gap="G4_MYLAR";
my $material_paint="SL_EC_paint";
my $material_support="G4_Al";
my $material_wrap="G4_MYLAR";
my $material_rod="G4_STAINLESS-STEEL";

# config mostly from CaloSim/include/Config.h
# Sampling ratio = 0.253017(EM) 0.285811(MIP) 0.34+-0.04(H)
# my $Nlayer=194;
# my $Thickness_lead = 0.05; #cm
# my $Thickness_scint = 0.15; #cm
# my $Thickness_gap = 0.024; #cm
# my $Thickness_shield = 0.5137*2; #cm, it should be 0.5612cm for X0
# my $Thickness_prescint = 2; #cm 
# my $Thickness_support = 2; #cm
# my $Thickness_preshower=$Thickness_shield+$Thickness_prescint; #3.0274mm
# my $z_shower=350;

#   my $Rmin = 110;
#   my $Rmax = 265;
#   my $Sphi  = 0;
#   my $Dphi  = 360;

#length in cm, angle in deg
my $Nlayer		= $parameters{"Nlayer"};
my $Thickness_lead	= $parameters{"Thickness_lead"};
my $Thickness_scint 	= $parameters{"Thickness_scint"};
my $Thickness_gap 	= $parameters{"Thickness_gap"};
my $Thickness_paint 	= $parameters{"Thickness_paint"};
my $Thickness_shield 	= $parameters{"Thickness_shield"};
my $Thickness_prescint 	= $parameters{"Thickness_prescint"};
my $Thickness_support 	= $parameters{"Thickness_support"};
my $Angle_module        = $parameters{"Angle_module"};
my $Radius_shower 	= $parameters{"Radius_shower"};
my $Radius_rod 		= $parameters{"Radius_rod"};
my $R_rod 		= $parameters{"R_rod"};
my $z_shower 		= $parameters{"z_shower"};
my $Rmin 		= $parameters{"Rmin"};
my $Rmax 		= $parameters{"Rmax"};
my $Sphi 		= $parameters{"Sphi"};
my $Dphi 		= $parameters{"Dphi"};
my $total_module 	= $parameters{"total_module"};
my $total_module_active	= $parameters{"total_module_active"};
my $thismodule		= $parameters{"thismodule"};


my @status;
my @idy;
my @idx;
my @y;
my @x;
my @sector;
for(my $i=1; $i<=$total_module; $i++)
{
  $status[$i-1] = $parameters{"id$i\_status"};
  $idy[$i-1] = $parameters{"id$i\_idy"};
  $idx[$i-1] = $parameters{"id$i\_idx"};
  $y[$i-1] = $parameters{"id$i\_y"};
  $x[$i-1] = $parameters{"id$i\_x"};  
  $sector[$i-1] = $parameters{"id$i\_sector"};
#   print "$idx[$i-1] \n";
}

my $vis_inner=0;

my $Thickness_layer=$Thickness_lead+$Thickness_scint+$Thickness_gap+$Thickness_paint; #0.234cm
my $Thickness_shower=$Thickness_layer*$Nlayer;  #45.395cm

my $z_support = $z_shower-$Thickness_shower/2-$Thickness_support/2;  #cm
my $z_prescint = $z_support-$Thickness_support/2-$Thickness_prescint/2; #cm 
my $z_shield = $z_prescint-$Thickness_prescint/2-$Thickness_shield/2;  #cm
print "z_shower z_support z_prescint z_shield $z_shower $z_support $z_prescint $z_shield\n";

# preshower lead is build as a while piece
# each EC shower module layer is made of sublayers "0.005cm paint + 0.05cm lead + 0.005cm paint + 0.012cm mylar + 0.15cm scint + 0.012cm mylar"
# we build the whole module with scint and be sensitive, then put all others sublayer volumes as its insensitive daughters 
# to minimize mem use, we use copy in geant4. we build 1 module and copy it. and within that module, we build 1st layer and copy it. we do that to preshower and support later also
# 6 stainless steel rods going through the whole EC module. we build them as daughter volume in each sublayer of 1st layer, then copy it. in each sublayer, we build 1st rod and copy it. rods in scint needs to build/copy one by one because scint is a whole volume. (we turn off rods in paint and mylar to save mem and speed up sim running for now, which is 15% rod mass)

sub solid_ec_geometry()
{
make_ec_shield(); # preshower lead
make_ec_module();
}

sub make_ec_shield
{
 my $Dz    = $Thickness_shield/2;

 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_shield";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z_shield*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "$color_abs";
 $detector{"type"}       = "Tube";
 $detector{"dimensions"} = "$Rmin*cm $Rmax*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = $material_abs;
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 1;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "no";
 $detector{"identifiers"} = "no";
 print_det(\%configuration, \%detector);
}

sub make_ec_module
{   
# build the module 424 x=-39.116cm y=-120.984cm first, at radius near 125cm, it exist in all three calorimeter, then make copy of it
# the original has to be made before copies, otherwise it can have segmentation fault. because GEMC make volume in alphabetic order of names, one way to make sure the original built first is to make the original name shorter

	  make_ec_module_shower($thismodule);
	  make_ec_module_prescint($thismodule);
	  make_ec_module_support($thismodule);
	  
#    for(my $i=$thismodule; $i<=$thismodule; $i++)   
   for(my $i=1; $i<=$total_module; $i++)   
    {
      if($i== $thismodule) {next;}
      if($status[$i-1] eq 1) {      
	my $R=sqrt($x[$i-1]*$x[$i-1]+$y[$i-1]*$y[$i-1]);
# 	if (130<$R && $R<140 ){
# 	  print "build this $i $x[$i-1] $y[$i-1]\n"
# 	}
	if (($Rmin)<$R && $R<($Rmax)){
# 	if (($Rmin)<$R){
	  make_ec_module_shower($i);
	  make_ec_module_prescint($i);
	  make_ec_module_support($i);
	}
# 	else {print "remove module $i $x[$i-1] $y[$i-1]\n";}
      }
    }
}

sub make_ec_module_shower()
{   
    my $id = $_[0];

    my %detector=init_det();
 
    
    my $Dz_shower = $Thickness_shower/2;

    if($id == $thismodule){ $detector{"name"}        = "$DetectorName\_showe$id";} #use shorter name so it can be built first
    else {$detector{"name"}        = "$DetectorName\_shower$id";}
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "$x[$id-1]*cm $y[$id-1]*cm $z_shower*cm";
    $detector{"rotation"}   = "0*deg 0*deg $Angle_module*deg";
    if ($sector[$id-1]%2 eq 0) {$detector{"color"}      = $color_wrap_sectoreven;}
    else {$detector{"color"}      = $color_wrap_sectorodd;}
    if($id == $thismodule){ 
    $detector{"type"}       = "Pgon";
    $detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_shower*cm $Dz_shower*cm";
    }    
    else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$thismodule";}         
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;   
    $detector{"sensitivity"} = $ec_sen;
    $detector{"hit_type"}    = $ec_sen;
    my $ID = $detID+$id;    
    $detector{"identifiers"} = "id manual $ID";
    print_det(\%configuration, \%detector); 

   if($id == $thismodule){
   for(my $j=1; $j<=$Nlayer; $j++){
	my $Dz_layer = $Thickness_layer/2;   
	my $layerZ = -$Dz_shower+($j-0.5)*$Thickness_layer;
	
	my $Dz_paint = $Thickness_paint/4;	
	my $Dz_lead = $Thickness_lead/2;	
	my $Dz_gap = $Thickness_gap/4;		
	my $Dz_scint = $Thickness_scint/2;		

	my $paintfrontZ = $layerZ-$Dz_layer+$Dz_paint;	
	my $leadZ = $layerZ-$Dz_layer+$Dz_paint*2+$Thickness_lead/2;	
	my $paintbackZ = $layerZ-$Dz_layer+$Dz_paint*2+$Thickness_lead+$Dz_paint;	
	my $gapfrontZ = $layerZ-$Dz_layer+$Dz_paint*2+$Thickness_lead+$Dz_paint*2+$Dz_gap;	
	my $scintZ = $layerZ-$Dz_layer+$Dz_paint*2+$Thickness_lead+$Dz_paint*2+$Dz_gap*2+$Thickness_scint/2;	
	my $gapbackZ = $layerZ-$Dz_layer+$Dz_paint*2+$Thickness_lead+$Dz_paint*2+$Dz_gap*2+$Thickness_scint+$Dz_gap;
	
        $detector{"name"}        = "$DetectorName\_showe$id\_lead$j";
        $detector{"mother"}      = "$DetectorName\_showe$id";
        $detector{"description"} = $detector{"name"};
        $detector{"pos"}        = "0*cm 0*cm $leadZ*cm";
        $detector{"rotation"}   = "0*deg 0*deg 0*deg";
        $detector{"color"}      = "$color_abs";
	if($j == 1){
	$detector{"type"}       = "Pgon";
	$detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_lead*cm $Dz_lead*cm";
	}
	else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_lead1";}	
        $detector{"material"}   = $material_abs;
        $detector{"mfield"}     = "no";
        $detector{"ncopy"}      = 1;
        $detector{"pMany"}       = 1;
        $detector{"exist"}       = 1;
	$detector{"visible"}     = $vis_inner;
	$detector{"style"}       = 0;     	
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";	
	print_det(\%configuration, \%detector);
	
	if($j == 1){
	  for(my $l=1; $l<=6; $l++){
	    $detector{"name"}        = "$DetectorName\_showe$id\_lead1_rod$l";
	    $detector{"mother"}      = "$DetectorName\_showe$id\_lead1";
	    $detector{"description"} = $detector{"name"};
	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm 0*cm";
	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
	    $detector{"color"}      = "$color_rod";
	    if($l == 1){
	    $detector{"type"}       = "Tube";
	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_lead*cm 0*deg 360*deg";    
	    }
	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_lead1_rod1";}	
	    $detector{"material"}   = $material_rod;
	    $detector{"mfield"}     = "no";
	    $detector{"ncopy"}      = 1;
	    $detector{"pMany"}       = 1;
	    $detector{"exist"}       = 1;
	    $detector{"visible"}     = $vis_inner;
	    $detector{"style"}       = 1;     	
	    $detector{"sensitivity"} = "no";
	    $detector{"hit_type"}    = "no";
	    $detector{"identifiers"} = "no";	
	    print_det(\%configuration, \%detector);
	  }
	}

	#build paint
        $detector{"name"}        = "$DetectorName\_showe$id\_paintfront$j";
        $detector{"mother"}      = "$DetectorName\_showe$id";
        $detector{"description"} = $detector{"name"};
        $detector{"pos"}        = "0*cm 0*cm $paintfrontZ*cm";
        $detector{"rotation"}   = "0*deg 0*deg 0*deg";
        $detector{"color"}      = "$color_paint";
	if($j == 1){
	$detector{"type"}       = "Pgon";
	$detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_paint*cm $Dz_paint*cm";
	}
	else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_paintfront1";}	
        $detector{"material"}   = "$material_paint";
        $detector{"mfield"}     = "no";
        $detector{"ncopy"}      = 1;
        $detector{"pMany"}       = 1;
        $detector{"exist"}       = 1;
	$detector{"visible"}     = $vis_inner;
	$detector{"style"}       = 0;
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);
	
# 	if($j == 1){
# 	  for(my $l=1; $l<=6; $l++){
# 	    $detector{"name"}        = "$DetectorName\_showe$id\_paintfront1_rod$l";
# 	    $detector{"mother"}      = "$DetectorName\_showe$id\_paintfront1";
# 	    $detector{"description"} = $detector{"name"};
# 	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
# 	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
# 	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm 0*cm";
# 	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
# 	    $detector{"color"}      = "$color_rod";
# 	    if($l == 1){
# 	    $detector{"type"}       = "Tube";
# 	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_paint*cm 0*deg 360*deg";    
# 	    }
# 	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_paintfront1_rod1";}	
# 	    $detector{"material"}   = $material_rod;
# 	    $detector{"mfield"}     = "no";
# 	    $detector{"ncopy"}      = 1;
# 	    $detector{"pMany"}       = 1;
# 	    $detector{"exist"}       = 1;
# 	    $detector{"visible"}     = $vis_inner;
# 	    $detector{"style"}       = 1;     	
# 	    $detector{"sensitivity"} = "no";
# 	    $detector{"hit_type"}    = "no";
# 	    $detector{"identifiers"} = "no";	
# 	    print_det(\%configuration, \%detector);
# 	  }
# 	}		
	
        $detector{"name"}        = "$DetectorName\_showe$id\_paintback$j";
        $detector{"mother"}      = "$DetectorName\_showe$id";
        $detector{"description"} = $detector{"name"};
        $detector{"pos"}        = "0*cm 0*cm $paintbackZ*cm";
        $detector{"rotation"}   = "0*deg 0*deg 0*deg";
        $detector{"color"}      = "$color_paint";
	if($j == 1){
	$detector{"type"}       = "Pgon";
	$detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_paint*cm $Dz_paint*cm";
	}
	else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_paintback1";}	
        $detector{"material"}   = "$material_paint";
        $detector{"mfield"}     = "no";
        $detector{"ncopy"}      = 1;
        $detector{"pMany"}       = 1;
        $detector{"exist"}       = 1;
	$detector{"visible"}     = $vis_inner;
	$detector{"style"}       = 0;
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);	
	
# 	if($j == 1){
# 	  for(my $l=1; $l<=6; $l++){
# 	    $detector{"name"}        = "$DetectorName\_showe$id\_paintback1_rod$l";
# 	    $detector{"mother"}      = "$DetectorName\_showe$id\_paintback1";
# 	    $detector{"description"} = $detector{"name"};
# 	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
# 	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
# 	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm 0*cm";
# 	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
# 	    $detector{"color"}      = "$color_rod";
# 	    if($l == 1){
# 	    $detector{"type"}       = "Tube";
# 	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_paint*cm 0*deg 360*deg";    
# 	    }
# 	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_paintback1_rod1";}	
# 	    $detector{"material"}   = $material_rod;
# 	    $detector{"mfield"}     = "no";
# 	    $detector{"ncopy"}      = 1;
# 	    $detector{"pMany"}       = 1;
# 	    $detector{"exist"}       = 1;
# 	    $detector{"visible"}     = $vis_inner;
# 	    $detector{"style"}       = 1;     	
# 	    $detector{"sensitivity"} = "no";
# 	    $detector{"hit_type"}    = "no";
# 	    $detector{"identifiers"} = "no";	
# 	    print_det(\%configuration, \%detector);
# 	  }
# 	}
	
	#build gap
        $detector{"name"}        = "$DetectorName\_showe$id\_gapfront$j";
        $detector{"mother"}      = "$DetectorName\_showe$id";
        $detector{"description"} = $detector{"name"};
        $detector{"pos"}        = "0*cm 0*cm $gapfrontZ*cm";
        $detector{"rotation"}   = "0*deg 0*deg 0*deg";
        $detector{"color"}      = "$color_gap";
	if($j == 1){
	$detector{"type"}       = "Pgon";
	$detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_gap*cm $Dz_gap*cm";
	}
	else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_gapfront1";}	
        $detector{"material"}   = "$material_gap";
        $detector{"mfield"}     = "no";
        $detector{"ncopy"}      = 1;
        $detector{"pMany"}       = 1;
        $detector{"exist"}       = 1;
	$detector{"visible"}     = $vis_inner;
	$detector{"style"}       = 0;
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);
	
# 	if($j == 1){
# 	  for(my $l=1; $l<=6; $l++){
# 	    $detector{"name"}        = "$DetectorName\_showe$id\_gapfront1_rod$l";
# 	    $detector{"mother"}      = "$DetectorName\_showe$id\_gapfront1";
# 	    $detector{"description"} = $detector{"name"};
# 	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
# 	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
# 	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm 0*cm";
# 	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
# 	    $detector{"color"}      = "$color_rod";
# 	    if($l == 1){
# 	    $detector{"type"}       = "Tube";
# 	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_gap*cm 0*deg 360*deg";    
# 	    }
# 	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_gapfront1_rod1";}	
# 	    $detector{"material"}   = $material_rod;
# 	    $detector{"mfield"}     = "no";
# 	    $detector{"ncopy"}      = 1;
# 	    $detector{"pMany"}       = 1;
# 	    $detector{"exist"}       = 1;
# 	    $detector{"visible"}     = $vis_inner;
# 	    $detector{"style"}       = 1;     	
# 	    $detector{"sensitivity"} = "no";
# 	    $detector{"hit_type"}    = "no";
# 	    $detector{"identifiers"} = "no";	
# 	    print_det(\%configuration, \%detector);
# 	  }
# 	}		

        $detector{"name"}        = "$DetectorName\_showe$id\_gapback$j";
        $detector{"mother"}      = "$DetectorName\_showe$id";
        $detector{"description"} = $detector{"name"};
        $detector{"pos"}        = "0*cm 0*cm $gapbackZ*cm";
        $detector{"rotation"}   = "0*deg 0*deg 0*deg";
        $detector{"color"}      = "$color_gap";
	if($j == 1){
	$detector{"type"}       = "Pgon";
	$detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_gap*cm $Dz_gap*cm";
	}
	else {  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_gapback1";}	
        $detector{"material"}   = "$material_gap";
        $detector{"mfield"}     = "no";
        $detector{"ncopy"}      = 1;
        $detector{"pMany"}       = 1;
        $detector{"exist"}       = 1;
	$detector{"visible"}     = $vis_inner;
	$detector{"style"}       = 0;
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);	
	
# 	if($j == 1){
# 	  for(my $l=1; $l<=6; $l++){
# 	    $detector{"name"}        = "$DetectorName\_showe$id\_gapback1_rod$l";
# 	    $detector{"mother"}      = "$DetectorName\_showe$id\_gapback1";
# 	    $detector{"description"} = $detector{"name"};
# 	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
# 	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
# 	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm 0*cm";
# 	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
# 	    $detector{"color"}      = "$color_rod";
# 	    if($l == 1){
# 	    $detector{"type"}       = "Tube";
# 	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_gap*cm 0*deg 360*deg";    
# 	    }
# 	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_gapback1_rod1";}	
# 	    $detector{"material"}   = $material_rod;
# 	    $detector{"mfield"}     = "no";
# 	    $detector{"ncopy"}      = 1;
# 	    $detector{"pMany"}       = 1;
# 	    $detector{"exist"}       = 1;
# 	    $detector{"visible"}     = $vis_inner;
# 	    $detector{"style"}       = 1;     	
# 	    $detector{"sensitivity"} = "no";
# 	    $detector{"hit_type"}    = "no";
# 	    $detector{"identifiers"} = "no";	
# 	    print_det(\%configuration, \%detector);
# 	  }
# 	}		
	
	#build scint rod
	  for(my $l=1; $l<=6; $l++){
	    if($j == 1){ $detector{"name"}        = "$DetectorName\_showe$id\_scint$j\_rod$l";}
	    else {$detector{"name"}        = "$DetectorName\_showe$id\_scint$j\_rod$l";}
	    $detector{"mother"}      = "$DetectorName\_showe$id";
	    $detector{"description"} = $detector{"name"};
	    my $rod_x=$Radius_rod*cos(($l-1)*60/180*3.1416);
	    my $rod_y=$Radius_rod*sin(($l-1)*60/180*3.1416);
	    $detector{"pos"}        = "$rod_x*cm $rod_y*cm $scintZ*cm";
	    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
	    $detector{"color"}      = "$color_rod";
	    if(($j == 1) && ($l == 1)){
	    $detector{"type"}       = "Tube";
	    $detector{"dimensions"} = "0*cm $R_rod*cm $Dz_scint*cm 0*deg 360*deg";    
	    }
	    else{  $detector{"type"}       = "CopyOf $DetectorName\_showe$id\_scint1_rod1";}	
	    $detector{"material"}   = $material_rod;
	    $detector{"mfield"}     = "no";
	    $detector{"ncopy"}      = 1;
	    $detector{"pMany"}       = 1;
	    $detector{"exist"}       = 1;
	    $detector{"visible"}     = $vis_inner;
	    $detector{"style"}       = 1;     	
	    $detector{"sensitivity"} = "no";
	    $detector{"hit_type"}    = "no";
	    $detector{"identifiers"} = "no";	
	    print_det(\%configuration, \%detector);
	  }
	
    }    
    }

}

sub make_ec_module_prescint()
{   
    my $id = $_[0];
 
    
    my $Dz_prescint = $Thickness_prescint/2;
    
    my %detector=init_det();

    if($id == $thismodule){ $detector{"name"}        = "$DetectorName\_prescin$id";} #use shorter name so it can be built first
    else {$detector{"name"}        = "$DetectorName\_prescint$id";}    
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "$x[$id-1]*cm $y[$id-1]*cm $z_prescint*cm";
    $detector{"rotation"}   = "0*deg 0*deg $Angle_module*deg";
    $detector{"color"}      = $color_scint;
    if($id == $thismodule){ 
    $detector{"type"}       = "Pgon";
    $detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_prescint*cm $Dz_prescint*cm";
    }    
    else {  $detector{"type"}       = "CopyOf $DetectorName\_prescin$thismodule";}      
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = $ec_ps_sen;
    $detector{"hit_type"}    = $ec_ps_sen;
    my $ID = $detID+10000+$id;
    $detector{"identifiers"} = "id manual $ID";
    print_det(\%configuration, \%detector); 
}

sub make_ec_module_support()
{   
    my $id = $_[0];
 
    
    my $Dz_support = $Thickness_support/2;
    
    my %detector=init_det();

    if($id == $thismodule){ $detector{"name"}        = "$DetectorName\_suppor$id";} #use shorter name so it can be built first
    else {$detector{"name"}        = "$DetectorName\_support$id";}    
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "$x[$id-1]*cm $y[$id-1]*cm $z_support*cm";
    $detector{"rotation"}   = "0*deg 0*deg $Angle_module*deg";
    $detector{"color"}      = $color_support;
    if($id == $thismodule){ 
    $detector{"type"}       = "Pgon";
    $detector{"dimensions"} = "$Sphi*deg $Dphi*deg 6*counts 2*counts 0*cm 0*cm $Radius_shower*cm $Radius_shower*cm -$Dz_support*cm $Dz_support*cm";
    }    
    else {  $detector{"type"}       = "CopyOf $DetectorName\_suppor$thismodule";}         
    $detector{"material"}   = "$material_support";
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

solid_ec_geometry();
1;