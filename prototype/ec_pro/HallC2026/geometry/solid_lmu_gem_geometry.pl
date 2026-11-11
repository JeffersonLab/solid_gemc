# use strict;
# use warnings;
#
# LMU-size four-GEM geometry adapted for SoLID HallC2026 beam-test placement.
#
# This version keeps the SoLID GEM detector positions and flip settings from
# solid_lmu_gem_geometry.pl, but uses a compact SoLID mother/region Z size.
# The internal layer thicknesses match the LMU GEM region1 values from
# lmu__volumes_default.txt, and the z positions are calculated by thickness
# accumulation centered around z = 0.
#
# LMU GEM reference:
#   transverse half-size: 5.1 cm x 5.1 cm
#   compact SoLID half-Z: mother 0.60 cm, region 0.55 cm
#
# The GEM layer stack is:
#   window   = kapton, Al, gas
#   cathode  = kapton, Cu, gas   <-- sensitive urwell gas
#   gem1     = Cu, kapton, Cu
#   gem2     = Cu, kapton, Cu
#   gem3     = Cu, kapton, Cu, glue
#   readout1 = Cu, g10, Cu
#   readout2 = glue, g10, Cu
#
# This file is intended to replace only the geometry file. Keep the same
# materials/hit/bank files you already use for solid_lmu_gem.

our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorBaseName = "solid_lmu_gem";
my $DetectorMother   = "cc_pro_tcd";

# Keep cathode gas sensitive, as in the LMU generated geometry.
# Set to 0 for passive geometry tests.
my $USE_URWELL_SENSITIVE_GAS = 1;

# --------------------------------------------------------------------------
# Keep SoLID beam-test GEM positions and flip settings
# --------------------------------------------------------------------------

my @module_z = (
    -71.12-69.49,
    -71.12-5.28,
    71.12+5.84,
    71.12+38.85,
);

# GEMs 1 and 3 are flipped; GEMs 2 and 4 are normal.
my @module_flipped = (
    1,
    0,
    1,
    0,
);

# --------------------------------------------------------------------------
# Exact LMU GEM region1 Box dimensions from lmu__volumes_default.txt
# --------------------------------------------------------------------------

my $lmu_hx        = 5.1;
my $lmu_hy        = 5.1;

# Compact SoLID mother/region Z sizes.
# These remove the large empty LMU container space but keep the real LMU GEM
# layer thicknesses and transverse size.
my $lmu_mother_hz = 0.60;
my $lmu_region_hz = 0.55;

# The fourth value is the full layer thickness in cm.
# This version does not use the original LMU absolute z positions; instead,
# it centers the full stack by accumulating thickness from -total_thickness/2.
my @lmu_gem_layers = (
    ["window_kapton",  "window",   "kapton", 0.00250],
    ["window_Al",      "window",   "Al",     0.00050],
    ["window_gas",     "window",   "gas",    0.30000],

    ["cathode_kapton", "cathode",  "kapton", 0.00500],
    ["cathode_Cu",     "cathode",  "Cu",     0.00050],
    ["cathode_gas",    "cathode",  "gas",    0.30000],

    ["gem1_Cu_1",      "gem1",     "Cu",     0.00500],
    ["gem1_kapton",    "gem1",     "kapton", 0.05000],
    ["gem1_Cu_2",      "gem1",     "Cu",     0.00500],

    ["gem2_Cu_1",      "gem2",     "Cu",     0.00500],
    ["gem2_kapton",    "gem2",     "kapton", 0.05000],
    ["gem2_Cu_2",      "gem2",     "Cu",     0.00500],

    ["gem3_Cu_1",      "gem3",     "Cu",     0.00500],
    ["gem3_kapton",    "gem3",     "kapton", 0.05000],
    ["gem3_Cu_2",      "gem3",     "Cu",     0.00500],
    ["gem3_glue",      "gem3",     "glue",   0.00280],

    ["readout1_Cu_1",  "readout1", "Cu",     0.00370],
    ["readout1_g10",   "readout1", "g10",    0.10000],
    ["readout1_Cu_2",  "readout1", "Cu",     0.00370],

    ["readout2_glue",  "readout2", "glue",   0.01250],
    ["readout2_g10",   "readout2", "g10",    0.10000],
    ["readout2_Cu",    "readout2", "Cu",     0.00370],
);

my %colorGem = (
    "kapton" => "bf0000",
    "Al"     => "2a3158",
    "gas"    => "afb0ba",
    "Cu"     => "fd7f00",
    "glue"   => "14ce3d",
    "g10"    => "888808",
    "air"    => "ff66cc",
);

sub solid_lmu_gem_geometry
{
    make_lmu_size_gem_modules();
}

sub make_lmu_size_gem_modules
{
    for(my $imod = 0; $imod < 4; $imod++)
    {
        my $module_id   = $imod + 1;
        my $mother_name = "${DetectorBaseName}_${module_id}";
        my $region_name = "${mother_name}_region";
        my $is_flipped  = $module_flipped[$imod];

        # ------------------------------------------------------------
        # Outer mother volume, corresponding to LMU_region1
        # ------------------------------------------------------------
        my %detector = init_det();
        $detector{"name"}        = $mother_name;
        $detector{"mother"}      = $DetectorMother;
        $detector{"description"} = "LMU-size compact GEM module $module_id outer mother";
        $detector{"pos"}         = "0*cm 0*cm $module_z[$imod]*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"color"}       = "ff66cc";
        $detector{"type"}        = "Box";
        $detector{"dimensions"}  = "$lmu_hx*cm $lmu_hy*cm $lmu_mother_hz*cm";
        $detector{"material"}    = "G4_AIR";
        set_common(\%detector);
        $detector{"visible"}     = 1;
        $detector{"style"}       = 0;
        $detector{"sensitivity"} = "no";
        $detector{"hit_type"}    = "no";
        $detector{"identifiers"} = "no";
        print_det(\%configuration, \%detector);

        # ------------------------------------------------------------
        # Inner region volume, corresponding to region1
        # ------------------------------------------------------------
        %detector = init_det();
        $detector{"name"}        = $region_name;
        $detector{"mother"}      = $mother_name;
        $detector{"description"} = "LMU-size compact GEM module $module_id inner region";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"color"}       = "4f84f7";
        $detector{"type"}        = "Box";
        $detector{"dimensions"}  = "$lmu_hx*cm $lmu_hy*cm $lmu_region_hz*cm";
        $detector{"material"}    = "G4_AIR";
        set_common(\%detector);
        $detector{"visible"}     = 1;
        $detector{"style"}       = 0;
        $detector{"sensitivity"} = "no";
        $detector{"hit_type"}    = "no";
        $detector{"identifiers"} = "no";
        print_det(\%configuration, \%detector);

        # ------------------------------------------------------------
        # LMU GEM daughter layers.
        # Calculate z positions by thickness accumulation, centered at z = 0.
        # For flipped modules, reverse the layer order.
        # ------------------------------------------------------------
        my @layers_this_module = $is_flipped ? reverse(@lmu_gem_layers) : @lmu_gem_layers;

        my $total_thickness = 0.0;
        foreach my $L (@layers_this_module) {
            my ($local_name, $group, $mat, $thick) = @$L;
            $total_thickness += $thick;
        }

        my $z = -$total_thickness/2.0;

        for(my $ilay = 0; $ilay <= $#layers_this_module; $ilay++)
        {
            my ($local_name, $group, $mat, $thick) = @{$layers_this_module[$ilay]};
            my $zhalf = $thick/2.0;
            my $zpos  = $z + $zhalf;
            $z += $thick;

            my $lname = "${mother_name}_${local_name}";

            %detector = init_det();
            $detector{"name"}        = $lname;
            $detector{"mother"}      = $region_name;
            $detector{"description"} = "LMU-size compact GEM $module_id $group $mat layer";
            $detector{"pos"}         = "0*cm 0*cm $zpos*cm";
            $detector{"rotation"}    = "0*deg 0*deg 0*deg";
            $detector{"color"}       = exists $colorGem{$mat} ? $colorGem{$mat} : "aaaaaa";
            $detector{"type"}        = "Box";
            $detector{"dimensions"}  = "$lmu_hx*cm $lmu_hy*cm $zhalf*cm";
            $detector{"material"}    = $mat;
            set_common(\%detector);
            $detector{"visible"}     = 1;
            $detector{"style"}       = 1;

            if($USE_URWELL_SENSITIVE_GAS && $group eq "cathode" && $mat eq "gas")
            {
                $detector{"sensitivity"} = "lmu_gem";
                $detector{"hit_type"}    = "lmu_gem";
		my $sector_id = $module_id;
                $detector{"identifiers"} =
                    "region manual 1 sector manual $sector_id chamber manual 1 layer manual 1 component manual 1";
            }
            else
            {
                $detector{"sensitivity"} = "no";
                $detector{"hit_type"}    = "no";
                $detector{"identifiers"} = "no";
            }

            print_det(\%configuration, \%detector);
        }
    }
}

sub set_common
{
    my $d = shift;

    $d->{"mfield"} = "no";
    $d->{"ncopy"}  = 1;
    $d->{"pMany"}  = 0;
    $d->{"exist"}  = 1;
}

1;
