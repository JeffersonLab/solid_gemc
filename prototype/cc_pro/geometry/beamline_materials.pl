#!/usr/bin/perl -w
###################################################################
# This script defines all materials used in the SoLID simulation. #
# If a new material is added, please specify your name, date and  #
# the purpose of adding this. Please make sure all materials      #
# have their names started with SL_mrpc_*                              #
# 
# Note: G4_* are for materials, so use percentage
#       For components,like C9H10, use the element "C" but not "G4_C"                                
#                                                                 #
#  -- Zhihong Ye, yez@jlab.org, 06/12/2014                        # 
###################################################################

############
#Note:
# (0) Define your new material in the section it belongs to,
#     DO-NOT just simply add it to the end of the file. 
#     Put your name and date near where you define your new items.
# (1) Pay attention to the density unit, which should be g/cm3
# (2) For elements, they should be like "He","C", ...
#     For materials, they should be like "G4_He", "G4_C"
# (3) "SL_mrpc_NewMaterial" is the newly defined material for SoLID    
# (4) If the new material is composed of elements, 
#     use "integer" to define the number of elements,
#     e.g "H 2 O 1"
# (5) If the new material is mixers of other materials,
#     use "mass fraction" to define the components, 
#     e.g. "G4_Si 0.70 G4_O 0.30"       
# (6) When you define a new material mixed by other materials, 
#     pay attention to how they are mixed,
#     e.g. by "mass fraction" or by "mole fraction" or volumn etc.
#
 
use strict;
use lib ("$ENV{GEMC}/io");
use utils;
use materials;

# Make sure the argument list is correct
# If not pring the help
if( scalar @ARGV != 1)
{
	help();
	exit;
}

# Loading configuration file from argument
our %configuration = load_configuration($ARGV[0]);

sub define_material
{

	my %mat;

	%mat = init_mat();
	$mat{"name"}          = "SL_beamline_N2_1atm";
	$mat{"density"}       = "0.00125053";  #in g/cm3
	$mat{"ncomponents"}   = "1";
	$mat{"components"}    = "G4_N 1";
	print_mat(\%configuration, \%mat);
	
}
define_material();


