use strict;
use warnings;

our %configuration;

# Variable Type is two chars.
# The first char:
#   R for raw integrated variables
#   D for dgt integrated variables
#   S for raw step by step variables
#   M for digitized multi-hit variables
#   V for voltage(time) variables
#
# The second char:
#   i for integers
#   d for doubles

my $bankID   = 2500;
my $bankname = "lmu_gem";

sub solid_lmu_gem_bank
{
    # LMU-style bank definition.
    #
    # In your adapted SoLID GEM/uRWell readout:
    #   sector    = GEM/uRWell module number
    #   layer     = readout layer, e.g. 1 = U, 2 = V
    #   component = strip number

    insert_bank_variable(\%configuration, $bankname, "bankid",    $bankID, "Di", "$bankname bank ID");

    insert_bank_variable(\%configuration, $bankname, "sector",    1, "Di", "module/sector number");
    insert_bank_variable(\%configuration, $bankname, "layer",     2, "Di", "readout layer: 1 = U, 2 = V");
    insert_bank_variable(\%configuration, $bankname, "component", 3, "Di", "strip number");

    insert_bank_variable(\%configuration, $bankname, "ADC_order", 4, "Di", "always 0");
    insert_bank_variable(\%configuration, $bankname, "ADC_ADC",   5, "Di", "ADC");
    insert_bank_variable(\%configuration, $bankname, "ADC_time",  6, "Dd", "time");
    insert_bank_variable(\%configuration, $bankname, "ADC_ped",   7, "Di", "pedestal");

    insert_bank_variable(\%configuration, $bankname, "hitn",     99, "Di", "hit number");
}

# Optional alias if your main driver expects define_bank()
sub define_bank
{
    solid_lmu_gem_bank();
}

1;
