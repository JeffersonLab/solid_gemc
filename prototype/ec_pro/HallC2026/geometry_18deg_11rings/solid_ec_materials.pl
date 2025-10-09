# #!/usr/bin/perl -w
###################################################################
# This script defines all materials used in the SoLID simulation. #
# If a new material is added, please specify your name, date and  #
# the purpose of adding this. Please make sure all materials      #
# have their names started with SL_*                              #
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
# (3) "SL_NewMaterial" is the newly defined material for SoLID    
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
# use lib ("$ENV{GEMC}/io");
# use lib ("$ENV{GEMC}/api/perl");
use utils;
use materials;

# # Help Message
# sub help()
# {
# 	print "\n Usage: \n";
# 	print "   materials.pl <configuration filename>\n";
#  	print "   Will create a simple scintillator material\n";
# 	exit;
# }
# 
# # Make sure the argument list is correct
# # If not pring the help
# if( scalar @ARGV != 1)
# {
# 	help();
# 	exit;
# }

# Loading configuration file from argument
# our %configuration = load_configuration($ARGV[0]);
# my $config_file   = "config.dat";
# our %configuration = load_configuration($config_file);

# One can change the "variation" here if one is desired different from the config.dat
# $configuration{"detector_name"} = "solid_PVDIS_ec_forwardangle";
# $configuration{"detector_name"} = "$ARGV[0]";
# $configuration{"variation"} = "Original";

sub define_materials
{
	my %mat;
	    
	%mat = init_mat();
	
# https://pdg.lbl.gov/2019/AtomicNuclearProperties/HTML/titanium_dioxide.html	
	$mat{"name"}          = "SL_EC_paint_TiO2";
	$mat{"description"}   = "SL_EC_paint_TiO2";
	$mat{"density"}       = "4.26";  # in g/cm3
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "Ti 1 O 2";
	print_mat(\%configuration, \%mat);

# https://pdg.lbl.gov/2019/AtomicNuclearProperties/HTML/barium_sulfate.html	
	$mat{"name"}          = "SL_EC_paint_BaSO4";
	$mat{"description"}   = "SL_EC_paint_BaSO4";
	$mat{"density"}       = "4.5";  # in g/cm3
	$mat{"ncomponents"}   = "3";
	$mat{"components"}    = "Ba 1 S 1 O 4";
	print_mat(\%configuration, \%mat);

# from Dong li
# 序号  组分	型号	重量	备注		Materials 	type		weight
# 1	聚酯	3#A	115			POLYESTER-3		3#A		115
# 2	环氧树脂	604	105			epoxy resin		604		105
# 3	环氧硅油	604硅油	10			Epoxy silicone oil	604silicone oil	10
# 4	硫酸钡	BaS04	95			barium sulfate		BaS04		95
# 5	钛白粉	R216	168			titanium pigment	R216		168
# 6	流平剂	922	5			flatting agent		922		5
# 7	光亮剂	701	5			brightener		701		5
# 8	消泡剂	安息香	2			antifoamer		benzoin		2
# 9	增电剂	SA2483	1.5			Electrostatic agent	SA2483		1.5

    #Chendi Shen  Here is the datasheet for the TiO2 powder paint that will be used for module SDU4: power painting.pdf. THU2 used a similar paint but that factory has been closed.
    # https://hallaweb.jlab.org/experiment/PVDIS/SoLID/EC/meetings/20180412/powder%20painting.pdf	
	
#     POLYESTER and epoxy (close to G4_POLYSTYRENE) fractionmass 0.45, TiO2 0.35, BaSO4 0.20
	$mat{"name"}          = "SL_EC_paint";
	$mat{"description"}   = "SL_EC_paint";
	$mat{"density"}       = "1.5";  # in g/cm3
	$mat{"ncomponents"}   = "3";
	$mat{"components"}    = "G4_POLYSTYRENE 0.45 SL_EC_paint_TiO2 0.35 SL_EC_paint_BaSO4 0.20";
	print_mat(\%configuration, \%mat);
    
}

define_materials();
1;
