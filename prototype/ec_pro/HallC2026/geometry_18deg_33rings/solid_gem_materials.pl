#!/usr/bin/perl -w
###################################################################
# This script defines all materials used in the SoLID simulation. #
# If a new material is added, please specify your name, date and  #
# the purpose of adding this. Please make sure all materials      #
# have their names started with SL_gem_mrpc_*                              #
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
# (3) "SL_gem_mrpc_NewMaterial" is the newly defined material for SoLID    
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
 
# use strict;
use lib ("$ENV{GEMC}/io");
use lib ("$ENV{GEMC}/api/perl");
use utils;
use materials;

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   this_script.pl <configuration filename>\n";
 	print "   Will create materials used in SoLID\n";
	exit;
}

# Make sure the argument list is correct
# If not pring the help
if( scalar @ARGV != 1)
{
	help();
	exit;
}

# Loading configuration file from argument
# our %configuration = load_configuration($ARGV[0]);

sub define_material
{

	my %mat;


##########################
# GEM Section 
##########################
	##############
	#GEM gas--basically ArCO2
	#
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_CO2";
	$mat{"description"}   = "CO2 mixed in GEM gas";
	$mat{"density"}       = "0.001977";  #in g/cm3
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "C 1 O 2";
	print_mat(\%configuration, \%mat);

	%mat = init_mat();
	$mat{"name"}          = "SL_gem_Argon";
	$mat{"description"}   = "Argon mixed in GEM gas";
	$mat{"density"}       = "0.0017823";  #in g/cm3
	$mat{"ncomponents"}   = "1";
	$mat{"components"}    = "G4_Ar 1";
	print_mat(\%configuration, \%mat);
		
#note		
# 	G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
# 	G4Material *ArCO2 = new G4Material("GEMgas", density_ArCO2, nel=2);
# 	ArCO2->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
# 	ArCO2->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;		
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_GEMgas";
	$mat{"description"}   = "GEM gas";
	$mat{"density"}       = "0.0018407";  #in g/cm3
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "SL_gem_CO2 0.32 SL_gem_Argon 0.68";
	print_mat(\%configuration, \%mat);

 
	##############
	#NEMAG10
	#
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_NEMAG10";
	$mat{"description"}   = "NEMA G10 in GEM";
	$mat{"density"}       = "1.700";  #in g/cm3
	$mat{"ncomponents"}   = "4";
	$mat{"components"}    = "Si 1 O 2 C 3 H 3";
	print_mat(\%configuration, \%mat);

	##############
	#NOMEX
	#
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_NOMEX_pure";
	$mat{"description"}   = "pure NOMEX";
	$mat{"density"}       = "1.380";  #in g/cm3
	$mat{"ncomponents"}   = "5";
	$mat{"components"}    = "G4_H 0.04 G4_C 0.54 G4_N 0.09 G4_O 0.10 G4_Cl 0.23";
	print_mat(\%configuration, \%mat);

	%mat = init_mat();
	$mat{"name"}          = "SL_gem_NOMEX";
	$mat{"description"}   = "NOMEX mixed with air in GEM";
	$mat{"density"}       = "0.040";  #in g/cm3
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "G4_AIR 0.55 SL_gem_NOMEX_pure 0.45";
	print_mat(\%configuration, \%mat);
	
	#Kapton
	#
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_Kapton";
	$mat{"description"}   = "Kapton";
	$mat{"density"}       = "1.420";  #in g/cm3
	$mat{"ncomponents"}   = "4";
	$mat{"components"}    = "G4_H 0.026362 G4_C 0.691133 G4_N 0.07327 G4_O 0.209235";
	print_mat(\%configuration, \%mat);	
	
	%mat = init_mat();
	$mat{"name"}          = "SL_gem_mylar";
	$mat{"description"}   = "SL_gem_mylar";
	$mat{"density"}       = "1.397";  # in g/cm3
	$mat{"ncomponents"}   = "3";
	$mat{"components"}    = "G4_H 0.042 G4_C 0.625 G4_O 0.333";
	print_mat(\%configuration, \%mat);		
}
define_material();


