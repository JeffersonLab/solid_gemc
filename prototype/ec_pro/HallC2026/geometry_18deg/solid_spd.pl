#!/usr/bin/perl -w

# use strict;
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
$configuration{"detector_name"} = "solid_spd";
# $configuration{"variation"} = "Original";

# To get the parameters proper authentication is needed.
# our %parameters    = get_parameters(%configuration);

#Geometry definition
require "./solid_spd_geometry.pl";
require "./solid_spd_virtualplane.pl";

#materials
require "./solid_spd_materials.pl";

#hit and bank definition Execute only when there are changes
#hit
require "./solid_spd_hit.pl";

# banks
require "./solid_spd_bank.pl";

