#!/usr/bin/perl -w

use strict;
use lib ("$ENV{GEMC}/io");
use lib ("$ENV{GEMC}/api/perl");
use parameters;
use utils;

use geometry;
use hit;
use bank;
use math;

use Math::Trig;
# use Math::MatrixReal;
# use Math::VectorReal;

# system("rm meic_det1_simple__*txt");

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   detector.pl <configuration filename>\n";
 	print "   Will create the detector\n";
 	print "   Note: The passport and .visa files must be present to connect to MYSQL. \n\n";
	exit;
}

# Make sure the argument list is correct
if( scalar @ARGV == 0 ) 
{
	help();
	exit;
}


# Loading configuration file and paramters
my $config_file   = $ARGV[0];
our %configuration = load_configuration($config_file);

# my $type = $ARGV[1];
# echo $type;
# my $material_gas="";
# if (char $type=="C4F8") $material_gas="SL_HGC_C4F8_oneatm";  
# else if (char $type=="CO2") $material_gas="SL_LGCCgas"; 
# else exit;

#if($#ARGV < 0){
#    do 'configs/pvdis_CLEO_nominal.pl';
#}else{
#    do 'configs/'.$ARGV[0];
#}

# One can change the "variation" here if one is desired different from the config.dat
$configuration{"detector_name"} = "cc_pro";
$configuration{"variation"} = "CO2LA";

# To get the parameters proper authentication is needed.
# our %parameters    = get_parameters(%configuration);

#Geometry definition
require "./cc_pro_geometry.pl";

#materials
require "./cc_pro_materials.pl";

#mirror
require "./cc_pro_mirror.pl";

#hit definition
require "./solid_hgc_hit.pl";

# bank definition
require "./solid_hgc_bank.pl";
