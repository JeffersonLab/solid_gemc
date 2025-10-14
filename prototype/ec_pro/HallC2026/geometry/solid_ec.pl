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

# use Math::Trig;
# use Math::MatrixReal;
# use Math::VectorReal;

# system("rm meic_det1_simple__*txt");

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   detector.pl <detector name>\n";
 	print "   Will create the detector\n";
 	print "   Note: The passport and .visa files must be present to connect to MYSQL. \n\n";
	exit;
}

# Make sure the argument list is correct
if( scalar @ARGV != 1) 
{
	help();
	exit;
}

# Loading configuration file and paramters
my $config_file   = $ARGV[0];
our %configuration = load_configuration($config_file);

# One can change the "variation" here if one is desired different from the config.dat
$configuration{"detector_name"} = "solid_ec";
$configuration{"variation"} = "Original";

# To get the parameters proper authentication is needed.
our %parameters    = get_parameters(%configuration);

# our $DetectorName=$ARGV[0];
# print "DetectorName $DetectorName \n";

# our $detID=0;
# if($DetectorName eq "solid_ec"){
#     $detID=3100000;
# }
# else{}
# print "detID $detID \n";

our $DetectorName= "solid_ec";
our $detID=3100000;

#Geometry definition
require "./solid_ec_geometry.pl";
require "./solid_ec_virtualplane.pl";

#materials
require "./solid_ec_materials.pl";

#hit and bank definition Execute only when there are changes
#hit
require "./solid_ec_hit.pl";
require "./solid_ec_ps_hit.pl";

# banks
require "./solid_ec_bank.pl";
require "./solid_ec_ps_bank.pl";

