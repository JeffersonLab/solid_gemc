#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_PVDIS_baffle_CLEO2_shldPOLY';

my $DetectorMother="root";

# Global variables governing some of the variations
my @color_abs;
my @sensitivity_abs;
my @hit_abs;
my @mother;
my @material_abs;
my @material_abs_within;

# System parameters
my $Nplate; # # of baffles
my $Dz;  # thickness
my $Nslit;  # of slits
my $Nblock; # # of radial layers
my $zc0; # center of baffle system
my $Dzc;  # baffle spacing

# Baffle parameters
my @rinin;  # inner radius of inner ring
my @rinout; # outer radius of inner ring
my @routin;  # inner radius of outer ring
my @routout; # outer radius of outer ring
my @offset;  # phi offset

# Slit parameterization for CLEO2 baffles:
#
# Slits in a given plate are hexagons in r-phi space described by 9
# numbers:
#
# 3 radii r0, r1, r2
# 3 low angles phi00 phi01 phi02
# 3 high angles phi10 phi11 phi12
#
# (in most plates actually the high angle sides are colinear and the
# slit is a pentagon but we still use this parameterization.)
#
# Layer centered at r = ri has blocks starting at sphi0i and ending
# at sphi1i where
#
#   sphi0i = phi01 + (ri-r1) * (phi01-phi00)/(r1-r0)  for ri < r1
#   sphi0i = phi01 + (ri-r1) * (phi02-phi01)/(r2-r1)  for ri > r1
#
# and similarly for sphi1i.
#
# These 9 parameters, except for r0 and r2, are computed from a linear model:
#
# r1 = c0_r1 + (z-zc) * c1_r1  etc.
#
# where zc is the z coordinate of the center of the baffles and z
# is the z coordinate of the center of the plate.
#
# Exception: for first two plates slit is a quadrilateral, not a
# hexagon, in r-phi space. Phi for inner vertices (at r = r0) are
# obtained by extrapolation or interpolation of the line segments
# (r1, phi01) to (r2, phi02) and (r1, phi11) to (r2,
# phi12). Parameters phi00 and phi10 are ignored.

# Variation can be any mix of the strings:

# Enclosure     Puts 1st baffle inside solid_PVDIS_target_LD2_enclosure_TACV
#               with vacuum in the slits
# Kill          Kryptonite baffles with vacuum in the slits
# Copper        Copper baffles
# NoInner       Inner rings on all plates are omitted
# Zigzag        Even # plates (starting with second one) have slits opened up

# c0, c1, c2, are constant, linear, quadratic coefficients
# Slit parameters
my $c0_r1; # middle radius
my $c1_r1;
my $c2_r1;
my $c0_phi00; # phi0 at inner radius
my $c1_phi00;
my $c2_phi00;
my $c0_phi01; # phi0 at middle radius
my $c1_phi01;
my $c2_phi01;
my $c0_phi02; # phi0 at outer radius
my $c1_phi02;
my $c2_phi02;
my $c0_phi10; # phi1 at inner radius
my $c1_phi10;
my $c2_phi10;
my $c0_phi11; # phi1 at middle radius
my $c1_phi11;
my $c2_phi11;
my $c0_phi12; # phi1 at outer radius
my $c1_phi12;
my $c2_phi12;

my @PlateZ;
my @Dphi_adj; # adjustments to Dphi

# Absorber parameters
my $Dz_gap;  # gap between absorber and baffle
my $Dz_shld;

sub solid_PVDIS_baffle_CLEO2_shldPOLY_geometry
{
# the first argument to this function becomes the variation
    $configuration{"variation"} = shift;

    $Nplate  = $parameters{"Nplate"}; 
    for (my $i = 0; $i < $Nplate; $i++)
    {
	$Dphi_adj[$i] = 0;
	$mother[$i] = $DetectorMother;
	# adjustments -- 
	if ($configuration{"variation"} =~ /Zigzag/)
	{
	    $Dphi_adj[$i] = ($i % 2 == 0 ? 0 : -1); # Open up even # baffles
	}
	if ($configuration{"variation"} =~ /Kill/)
	{
	    # Inert Kryptonite
	    $material_abs[$i] = "Kryptonite";
	    $material_abs_within[$i] = "G4_Galactic";
	    $sensitivity_abs[$i] = "no";
	    $hit_abs[$i] = "no";
	}
	else
	{
	    # Inert borated poly
	    $material_abs[$i] = "SL_PVDIS_baffle_BorPoly";
	    $material_abs_within[$i] = "G4_AIR";
	    $sensitivity_abs[$i] = "no";
	    $hit_abs[$i] = "no";
	}
	if ($i == 0 && $configuration{"variation"} =~ /Enclosure/)
	{
	    $mother[$i] = "solid_PVDIS_target_enclosure_TACV";
	    $material_abs_within[$i] = "G4_Galactic";
	}
	$color_abs[$i] = "9999FF";
    }


    $Dz   = $parameters{"Dz"};  
    $Nslit  = $parameters{"Nslit"};
    $Nblock  = $parameters{"Nblock"};
    $zc0 = $parameters{"zc0"}; 
    $Dzc = $parameters{"Dzc"};
    $Dz_gap = 0.1; # space for baffle observers (full length)
    $Dz_shld =  $Dzc - 2 * $Dz - 2 * $Dz_gap ; #full length (Polycone)
    @rinin = ($parameters{"rinin1"}, $parameters{"rinin2"}, 
	      $parameters{"rinin3"}, $parameters{"rinin4"}, 
	      $parameters{"rinin5"}, $parameters{"rinin6"}, 
	      $parameters{"rinin7"}, $parameters{"rinin8"}, 
	      $parameters{"rinin9"}, $parameters{"rinin10"}, 
	      $parameters{"rinin11"});
    @rinout = ($parameters{"rinout1"}, $parameters{"rinout2"}, 
	       $parameters{"rinout3"}, $parameters{"rinout4"}, 
	       $parameters{"rinout5"}, $parameters{"rinout6"}, 
	       $parameters{"rinout7"}, $parameters{"rinout8"}, 
	       $parameters{"rinout9"}, $parameters{"rinout10"}, 
	       $parameters{"rinout11"});
    @routin = ($parameters{"routin1"}, $parameters{"routin2"}, 
	       $parameters{"routin3"}, $parameters{"routin4"}, 
	       $parameters{"routin5"}, $parameters{"routin6"}, 
	       $parameters{"routin7"}, $parameters{"routin8"}, 
	       $parameters{"routin9"}, $parameters{"routin10"}, 
	       $parameters{"routin11"});
    @routout = ($parameters{"routout1"}, $parameters{"routout2"}, 
		$parameters{"routout3"}, $parameters{"routout4"}, 
		$parameters{"routout5"}, $parameters{"routout6"}, 
		$parameters{"routout7"}, $parameters{"routout8"}, 
		$parameters{"routout9"}, $parameters{"routout10"}, 
		$parameters{"routout11"});
    @offset = ($parameters{"offset1"}, $parameters{"offset2"}, 
	       $parameters{"offset3"}, $parameters{"offset4"}, 
	       $parameters{"offset5"}, $parameters{"offset6"}, 
	       $parameters{"offset7"}, $parameters{"offset8"}, 
	       $parameters{"offset9"}, $parameters{"offset10"}, 
	       $parameters{"offset11"});

    $c0_r1 = $parameters{"c0_r1"};
    $c1_r1 = $parameters{"c1_r1"};
    $c2_r1 = $parameters{"c2_r1"};
    $c0_phi00 = $parameters{"c0_phi00"};
    $c1_phi00 = $parameters{"c1_phi00"};
    $c2_phi00 = $parameters{"c2_phi00"};
    $c0_phi01 = $parameters{"c0_phi01"};
    $c1_phi01 = $parameters{"c1_phi01"};
    $c2_phi01 = $parameters{"c2_phi01"};
    $c0_phi02 = $parameters{"c0_phi02"};
    $c1_phi02 = $parameters{"c1_phi02"};
    $c2_phi02 = $parameters{"c2_phi02"};
    $c0_phi10 = $parameters{"c0_phi10"};
    $c1_phi10 = $parameters{"c1_phi10"};
    $c2_phi10 = $parameters{"c2_phi10"};
    $c0_phi11 = $parameters{"c0_phi11"};
    $c1_phi11 = $parameters{"c1_phi11"};
    $c2_phi11 = $parameters{"c2_phi11"};
    $c0_phi12 = $parameters{"c0_phi12"};
    $c1_phi12 = $parameters{"c1_phi12"};
    $c2_phi12 = $parameters{"c2_phi12"};
    
    @PlateZ = ();
    for (my $i = 0; $i < $Nplate; ++$i)
    {
	my $z = $zc0 + $Dzc * ($i - ($Nplate-1) * 0.5);
	push @PlateZ, $z;
    }

    if ($configuration{"variation"} !~ /NoInner/)
    {
	make_CLEO2_baffle_shldPOLY_plate_inner_absorber();
	make_CLEO2_baffle_shldPOLY_plate_inner();
    }
    make_CLEO2_baffle_shldPOLY_plate_outer_absorber();
    make_CLEO2_baffle_shldPOLY_plate_outer();
    make_CLEO2_baffle_shldPOLY_plate();
    make_CLEO2_baffle_shldPOLY_blocks();
}

1; # return true


sub zvals1_shldPOLY
{
# Given $c0, $c1, $zc, @z, return array of $c0+$c1*($z[i]-$zc)
    my ($c0, $c1, $zc, @z) = @_;
    return map ($c0 + $c1 * ($_ - $zc), @z)
}

sub zvals2_shldPOLY
{
# Given $c0, $c1, $c2, $zc, @z, return array of $c0+$c1*($z[i]-$zc)+$c2*()^2
    my ($c0, $c1, $c2, $zc, @z) = @_;
    return map ($c0 + ($c1 + $c2 * ($_ - $zc)) * ($_ - $zc), @z)
}

sub make_CLEO2_baffle_shldPOLY_plate_inner_absorber
{
    my @Rin = ();
    my @Rout = ();
    my @Zabs = ();
 
# For now the inner radii are are calculated for a hardcoded beamline shape
# It's from 5.56 at 49.5 to 20.03 at 114.5, then 20.03 after that.

    for (my $n=2; $n<=$Nplate; $n++)
    {
	my $z0 = $PlateZ[$n-1]-$Dz-$Dz_gap;
	my $z1 = $n < 9 ? $PlateZ[$n-1]+$Dz+$Dz_gap
	    : $n < 11 ? $PlateZ[$n-1] + 6.0
	    : 188.99;
	my $rad0 = $z0 >= 114.5 ? 20.03 :
	    5.56 + (20.03-5.56) / (114.5-49.5) * ($z0-49.5);
	my $rad1 = $z1 >= 114.5 ? 20.03 :
	    5.56 + (20.03-5.56) / (114.5-49.5) * ($z1-49.5);
	push @Rin, $rad0, $rad1;
	push @Rout, $rinin[$n-1]-0.1, $rinin[$n-1]-0.1;
	push @Zabs, $z0, $z1;
    }

    my $dims0 = "0*deg 360*deg 20*counts ";
    my $dims1 = "";
    my $dims2 = "";
    my $dims3 = "";
    for (my $i = 0; $i < 20; ++$i)
    {
	$dims1 .= "$Rin[$i]*cm ";
	$dims2 .= "$Rout[$i]*cm ";
	$dims3 .= "$Zabs[$i]*cm ";
    }
    my %detector = init_det();
    $detector{"name"}        = "$DetectorName\_plateinner_absorber";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs[0]";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = $dims0 . $dims1 . $dims2 . $dims3;
    $detector{"material"}   = "$material_abs[0]";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_abs[0]";
    $detector{"hit_type"}    = "$hit_abs[0]";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_CLEO2_baffle_shldPOLY_plate_inner
{
    for (my $n=2; $n<=$Nplate-3; $n++)
    {
	next if ($rinout[$n-1] <= $rinin[$n-1]); # make no ring if zero (or neg) thickness
	my $n_c     = cnumber($n-1, 1);
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_plateinner$n_c";
	$detector{"mother"}      = "$mother[$n-1]" ;
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}        = "0*cm 0*cm " . ($PlateZ[$n-1] + $Dz + $Dz_gap) . "*cm";
	$detector{"rotation"}   = "0*deg 0*deg 0*deg";
	$detector{"color"}      = "$color_abs[$n-1]";
	$detector{"type"}       = "Polycone";
	$detector{"dimensions"} = "0*deg 360*deg 2*counts $rinin[$n-1]*cm $rinin[$n]*cm $rinout[$n-1]*cm $rinout[$n]*cm 0.0*cm $Dz_shld*cm";
	$detector{"material"}   = "$material_abs[$n-1]";
	$detector{"mfield"}     = "no";
	$detector{"ncopy"}      = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "$sensitivity_abs[$n-1]";
	$detector{"hit_type"}    = "$hit_abs[$n-1]";
	my $id = 1000000 + $n;
	$detector{"identifiers"} = "id manual $id";
	print_det(\%configuration, \%detector);
    }
}

sub make_CLEO2_baffle_shldPOLY_plate_outer_absorber
{
    my @Rin = ();
    my @Rout = ();
    my @Zabs = ();
 
    for (my $n=2; $n<=$Nplate; $n++)
    {
	my $z0 = $PlateZ[$n-1]-$Dz-$Dz_gap;
	my $z1 = $n < 9 ? $PlateZ[$n-1]+$Dz+$Dz_gap
	    : $n < 11 ? $PlateZ[$n-1] + 6.0
	    : 188.99;
	push @Rin, $routout[$n-1], $routout[$n-1];
	push @Rout, 140, 140;
	push @Zabs, $z0, $z1;
    }

    my $dims0 = "0*deg 360*deg " . ($#Rin+1) . "*counts ";
    my $dims1 = "";
    my $dims2 = "";
    my $dims3 = "";
    for (my $i = 0; $i <= $#Rin; ++$i)
    {
	$dims1 .= "$Rin[$i]*cm ";
	$dims2 .= "$Rout[$i]*cm ";
	$dims3 .= "$Zabs[$i]*cm ";
    }
    my %detector = init_det();
    $detector{"name"}        = "$DetectorName\_plateouter_absorber";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "$color_abs[0]";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = $dims0 . $dims1 . $dims2 . $dims3;
    $detector{"material"}   = "$material_abs[0]";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "$sensitivity_abs[0]";
    $detector{"hit_type"}    = "$hit_abs[0]";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}

sub make_CLEO2_baffle_shldPOLY_plate_outer
{
    for (my $n=2; $n<=$Nplate-3; $n++)
    {
	my $n_c     = cnumber($n-1, 1);
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_plateouter$n_c";
	$detector{"mother"}      = "$mother[$n-1]" ;
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}        = "0*cm 0*cm " . ($PlateZ[$n-1] + $Dz + $Dz_gap) . "*cm";
	$detector{"rotation"}   = "0*deg 0*deg 0*deg";
	$detector{"color"}      = "$color_abs[$n-1]";
	$detector{"type"}       = "Polycone";
	$detector{"dimensions"} = "0*deg 360*deg 2*counts $routin[$n-1]*cm $routin[$n]*cm $routout[$n-1]*cm $routout[$n]*cm 0.0*cm $Dz_shld*cm";
	$detector{"material"}   = "$material_abs[$n-1]";
	$detector{"mfield"}     = "no";
	$detector{"ncopy"}      = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 1;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "$sensitivity_abs[$n-1]";
	$detector{"hit_type"}    = "$hit_abs[$n-1]";
	my $id = 3000000 + $n;
	$detector{"identifiers"} = "id manual $id";
	print_det(\%configuration, \%detector);
    }
}

sub make_CLEO2_baffle_shldPOLY_plate   #virtual container for baffle plate
{
    for (my $n=2; $n<=$Nplate-3; $n++)
    {
	my $n_c     = cnumber($n-1, 1);
	my %detector=init_det();
	$detector{"name"}        = "$DetectorName\_plate$n_c";
	$detector{"mother"}      = "$mother[$n-1]";
	$detector{"description"} = $detector{"name"};
	$detector{"pos"}        = "0*cm 0*cm " . ($PlateZ[$n-1] + $Dz + $Dz_gap) . "*cm";
	$detector{"rotation"}   = "0*deg 0*deg 0*deg";
	$detector{"color"}      = "000000";
	$detector{"type"}       = "Polycone";
	$detector{"dimensions"} = "0*deg 360*deg 2*counts $rinout[$n-1]*cm $rinout[$n]*cm $routin[$n-1]*cm $routin[$n]*cm 0.0*cm $Dz_shld*cm";
	$detector{"material"}   = "$material_abs_within[$n-1]";
	$detector{"mfield"}     = "no";
	$detector{"ncopy"}      = 1;
	$detector{"pMany"}       = 1;
	$detector{"exist"}       = 1;
	$detector{"visible"}     = 0;
	$detector{"style"}       = 1;
	$detector{"sensitivity"} = "no";
	$detector{"hit_type"}    = "no";
	$detector{"identifiers"} = "no";
	print_det(\%configuration, \%detector);
    }
}

#making baffle blocks, plate by plate, slit by slit, block by block
sub make_CLEO2_baffle_shldPOLY_blocks
{
    my @R1 = &zvals2_shldPOLY ($c0_r1, $c1_r1, $c2_r1, $zc0, @PlateZ);
    my @Phi00 = &zvals2_shldPOLY ($c0_phi00, $c1_phi00, $c2_phi00, $zc0, @PlateZ);
    my @Phi01 = &zvals2_shldPOLY ($c0_phi01, $c1_phi01, $c2_phi01, $zc0, @PlateZ);
    my @Phi02 = &zvals2_shldPOLY ($c0_phi02, $c1_phi02, $c2_phi02, $zc0, @PlateZ);
    my @Phi10 = &zvals2_shldPOLY ($c0_phi10, $c1_phi10, $c2_phi10, $zc0, @PlateZ);
    my @Phi11 = &zvals2_shldPOLY ($c0_phi11, $c1_phi11, $c2_phi11, $zc0, @PlateZ);
    my @Phi12 = &zvals2_shldPOLY ($c0_phi12, $c1_phi12, $c2_phi12, $zc0, @PlateZ);

    my @x = ([], [], [], [], [], [], [], [], [], [], []);
    
    for (my $ip = 0; $ip < $Nplate; ++$ip)
    {
	my $dr = ($routin[$ip] - $rinout[$ip]) / $Nblock;

	foreach my $r (0..$Nblock-1)
	{
	    my $r0b = $rinout[$ip] + $r * $dr;
	    my $rc = $r0b + 0.5 * $dr;
	    my $r1b = $r0b + $dr;
	    my $phi00a = $Phi00[$ip];
	    my $phi01a = $Phi01[$ip];
	    my $phi02a = $Phi02[$ip];
	    my $phi10a = $Phi10[$ip];
	    my $phi11a = $Phi11[$ip];
	    my $phi12a = $Phi12[$ip];
	    if ($ip < 2)
	    {
		my $phi = $phi01a + ($rc - $R1[$ip]) *
		    ($phi02a-$phi01a) / ($routin[$ip]-$R1[$ip]);
		my $Dphi = $phi11a + ($rc - $R1[$ip]) *
		    ($phi12a-$phi11a) / ($routin[$ip]-$R1[$ip]) - $phi;
		$Dphi += 360.0/$Nslit - 12.0; # adjust for non 30 sectors
		push @{$x[$ip]}, ($r0b, $r1b, $phi, $phi+$Dphi);
	    }
	    else
	    {
		my $phi = $phi01a + ($rc - $R1[$ip]) *
		    ($rc < $R1[$ip] ? ($phi00a-$phi01a) / ($rinout[$ip]-$R1[$ip])
		     : ($phi02a-$phi01a) / ($routin[$ip]-$R1[$ip]));
		my $Dphi = $phi11a + ($rc - $R1[$ip]) *
		    ($rc < $R1[$ip] ? ($phi10a-$phi11a) / ($rinout[$ip]-$R1[$ip])
		     : ($phi12a-$phi11a) / ($routin[$ip]-$R1[$ip])) - $phi;
		$Dphi += 360.0/$Nslit - 12.0; # adjust for non 30 sectors
		push @{$x[$ip]}, ($r0b, $r1b, $phi, $phi+$Dphi);
	    }
	}
    }

    for (my $n=2; $n<=$Nplate-3; $n++)
    {
	my $n_c     = cnumber($n-1, 1);
	for(my $i=1; $i<=$Nslit; $i++){  # making all slits
# was
#	my $slit_rotation = ($i-1)*12-$offset[$n-1]; #note the minus sign here

	    my $slit_rotation = -(96.0 + ($i-1) * 360.0/$Nslit + $offset[$n-1]);
# 	    $slit_rotation -= 360.0 if $slit_rotation > 360.0;

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
	    my $Rin2  = $x[$n][0];
	    my $Rout = $x[$n-1][$Nblock*4-3];
	    my $Rout2 = $x[$n][$Nblock*4-3];
	    my $Sphi = 0.5*($x[$n-1][2]+$x[$n][2]+$offset[$n-1]-$offset[$n]);
	    my $Dphi = 360./$Nslit;
	    $detector{"dimensions"}  = "$Sphi*deg $Dphi*deg 2*counts $Rin*cm $Rin2*cm $Rout*cm $Rout2*cm 0.0*cm $Dz_shld*cm";
	    $detector{"material"}    = "$material_abs_within[$n-1]";
	    $detector{"mfield"}      = "no";
	    $detector{"ncopy"}       = 1;
	    $detector{"pMany"}       = 1;
	    $detector{"exist"}       = 1;
	    $detector{"visible"}     = 0;
	    $detector{"style"}       = 0;
	    $detector{"sensitivity"} = "no";
	    $detector{"hit_type"}    = "no";
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
		$detector{"color"}       = "$color_abs[$n-1]";
		$detector{"type"}        = "Polycone";
		my $Rin  = $x[$n-1][($j-1)*4+0];
		my $Rin2  = $x[$n][($j-1)*4+0];
		my $Rout = $x[$n-1][($j-1)*4+1];
		my $Rout2 = $x[$n][($j-1)*4+1];
		my $Sphi =  0.5*($x[$n-1][($j-1)*4+2]+$x[$n][($j-1)*4+2]);
		my $Dphi =  0.5*($x[$n-1][($j-1)*4+3]-$x[$n-1][($j-1)*4+2]+$Dphi_adj[$n-1]
				 +$x[$n][($j-1)*4+3]-$x[$n][($j-1)*4+2]+$Dphi_adj[$n]);
		$detector{"dimensions"}  = "$Sphi*deg $Dphi*deg 2*counts $Rin*cm $Rin2*cm $Rout*cm $Rout2*cm 0.0*cm $Dz_shld*cm";
		$detector{"material"}    = "$material_abs[$n-1]";
		$detector{"mfield"}      = "no";
		$detector{"ncopy"}       = 1;
		$detector{"pMany"}       = 1;
		$detector{"exist"}       = 1;
		$detector{"visible"}     = 1;
		$detector{"style"}       = 1;
		$detector{"sensitivity"} = "$sensitivity_abs[$n-1]";
		$detector{"hit_type"}    = "$hit_abs[$n-1]";
		my $id = 2000000 + $j*10000 + $i*100 + $n;
		$detector{"identifiers"} = "id manual $id";
		print_det(\%configuration, \%detector);
	    }
	}
    }
}
