use strict;
use warnings;

use lib ("$ENV{GEMC}/io");
use lib ("$ENV{GEMC}/api/perl");

use geometry;

our %configuration;

my $DetectorBaseName = "solid_urgroove";
my $DetectorMother   = "cc_pro_tcd";

# Keep manual SoLID-style geometry now.
# Later, if a uRGroove/uRWell volume table is available, set this to 1
# and fill get_volume_info() near the bottom.
my $USE_VOLUME_TABLE = 0;

# If you have a dedicated C++ hit process named "urgroove", keep this as urgroove.
# If you want to temporarily reuse the already-compiled uRWell hit process,
# you may set it to "urwell" for testing.
my $HIT_TYPE = "urgroove";

my $NMODULES = 2;

my $x0 = 0;
my $y0 = 0;

# Positions of the two uRGroove modules.
# These use the first two positions from the current uRWell geometry style.
# Adjust as needed.
my @module_z = (
    -71.12 -38.66,
    -71.12 -22.15,
);

# No flip by default.
my @module_flipped = (
    0,
    0,
);

# GEM-like transverse size.
my $manual_hx = 5.12;  # cm
my $manual_hy = 5.12;  # cm

# uRGroove detail configuration.
# Layer stack and thicknesses follow the information from Jaydeep:
#   Kapton, uRGrooveGas, Copper, EM528, Carbon, Kapton,
#   Copper, uRGrooveGas, Copper, Kapton
# with thicknesses in cm:
#   0.0035, 0.3, 0.001, 0.0028, 0.00001,
#   0.005, 0.001, 0.3, 0.0005, 0.0035
#
# Note: the geometry uses valid GEMC material names:
#   uRGrooveGas -> urgroove_gas
#   Copper      -> Cu
#   EM528       -> EM528
#   Carbon      -> dlc
my $readout_pitch_cm    = 0.020;  # 200 um
my $groove_top_width_cm = 0.007;  # 70 um

# For 10.24 cm full width and 0.020 cm pitch, this is 512 strips.
my $approx_nstrips = int((2.0*$manual_hx)/$readout_pitch_cm + 0.5);

my %color = (
    "kapton"       => "bf0000",
    "Al"           => "2a3158",
    "urgroove_gas" => "afb0ba",
    "gas"          => "afb0ba",
    "Cu"           => "fd7f00",
    "dlc"          => "14b6ce",
    "EM528"        => "14ce3d",
    "glue"         => "14ce3d",
    "g10"          => "aa44d8"
);

sub set_common
{
    my ($d) = @_;

    $d->{"rotation"} = "0*deg 0*deg 0*deg" unless defined $d->{"rotation"};
    $d->{"mfield"}   = "no";
    $d->{"ncopy"}    = 1;
    $d->{"pMany"}    = 0;
    $d->{"exist"}    = 1;
    $d->{"visible"}  = 1 unless defined $d->{"visible"};
    $d->{"style"}    = 1 unless defined $d->{"style"};

    $d->{"sensitivity"} = "no" unless defined $d->{"sensitivity"};
    $d->{"hit_type"}    = "no" unless defined $d->{"hit_type"};
    $d->{"identifiers"} = "no" unless defined $d->{"identifiers"};
}

sub print_checked_det
{
    my ($d) = @_;

    foreach my $k (qw(name mother description pos rotation color type dimensions material mfield ncopy pMany exist visible style sensitivity hit_type identifiers)) {
        if (!defined $d->{$k}) {
            print "MISSING: $k in $d->{name}\n";
        }
        elsif ($d->{$k} eq "") {
            print "EMPTY: $k in $d->{name}\n";
        }
    }

    print_det(\%configuration, $d);
}

sub make_urgroove
{
    my %detector;

    # uRGroove layer stack from Jaydeep.
    # Format: [group, GEMC material name, thickness in cm, display label].
    # The display label is kept in the volume description so the generated geometry
    # file is easy to compare with the requested layer list.
    my @layers = (
        ["window",     "kapton",       0.0035,  "Kapton"],
        ["gas1",       "urgroove_gas", 0.3000,  "uRGrooveGas"],
        ["copper1",    "Cu",           0.0010,  "Copper"],
        ["adhesive",   "EM528",         0.0028,  "EM528"],
        ["carbon",     "dlc",          0.00001, "Carbon"],
        ["kapton_mid", "kapton",       0.0050,  "Kapton"],
        ["copper2",    "Cu",           0.0010,  "Copper"],
        ["gas2",       "urgroove_gas", 0.3000,  "uRGrooveGas"],
        ["copper3",    "Cu",           0.0005,  "Copper"],
        ["kapton_out", "kapton",       0.0035,  "Kapton"],
    );

    my $total_thickness = 0.0;
    foreach my $L (@layers) {
        $total_thickness += $L->[2];
    }

    my $mother_z_half = $total_thickness/2.0 + 0.20;

    for (my $imod = 0; $imod < $NMODULES; $imod++)
    {
        my $DetectorName = "${DetectorBaseName}_${imod}";
        my $z0 = $module_z[$imod];

        %detector = init_det();

        $detector{"name"}        = $DetectorName;
        $detector{"mother"}      = $DetectorMother;
        $detector{"description"} =
            "SoLID uRGroove module $imod mother; groove top ${groove_top_width_cm} cm; pitch ${readout_pitch_cm} cm; approx strips $approx_nstrips";
        $detector{"pos"}         = "$x0*cm $y0*cm $z0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"color"}       = "aaaaaa";
        $detector{"type"}        = "Box";
        $detector{"dimensions"}  = "$manual_hx*cm $manual_hy*cm $mother_z_half*cm";

        if ($USE_VOLUME_TABLE) {
            my ($type_from_table, $dim_from_table) =
                get_volume_info("mother", $imod, "mother", "air", 0, $module_flipped[$imod]);

            if (defined $type_from_table && defined $dim_from_table) {
                $detector{"type"}       = $type_from_table;
                $detector{"dimensions"} = $dim_from_table;
            }
        }

        $detector{"material"}    = "G4_AIR";
        $detector{"visible"}     = 1;
        $detector{"style"}       = 0;

        set_common(\%detector);
        print_checked_det(\%detector);

        my $is_flipped = $module_flipped[$imod];
        my @layers_this_module = $is_flipped ? reverse(@layers) : @layers;

        my $z = -$total_thickness/2.0;
        my $i = 0;

        foreach my $L (@layers_this_module)
        {
            my ($group, $mat, $thick, $label) = @$L;
            my $zhalf = $thick/2.0;
            my $zpos = $z + $thick/2.0;
            $z += $thick;

            %detector = init_det();

            $detector{"name"}        = "${DetectorName}_${group}_${mat}_${i}";
            $detector{"mother"}      = $DetectorName;
            $detector{"description"} = "uRGroove module $imod $group layer; requested material label $label; thickness $thick cm";
            $detector{"pos"}         = "0*cm 0*cm $zpos*cm";
            $detector{"rotation"}    = "0*deg 0*deg 0*deg";
            $detector{"color"}       = exists $color{$mat} ? $color{$mat} : "aaaaaa";

            $detector{"type"}        = "G4Trap";
            $detector{"dimensions"} =
                "$zhalf*cm 0*deg 0*deg " .
                "$manual_hy*cm $manual_hx*cm $manual_hx*cm 0*deg " .
                "$manual_hy*cm $manual_hx*cm $manual_hx*cm 0*deg";

            if ($USE_VOLUME_TABLE) {
                my ($type_from_table, $dim_from_table) =
                    get_volume_info("layer", $imod, $group, $mat, $i, $is_flipped);

                if (defined $type_from_table && defined $dim_from_table) {
                    $detector{"type"}       = $type_from_table;
                    $detector{"dimensions"} = $dim_from_table;
                }
            }

            $detector{"material"}    = $mat;
            $detector{"visible"}     = 1;
            $detector{"style"}       = 1;

            # Sensitive ionization volume:
            # only gas1 is sensitive, based on Jaydeep's confirmation.
            # uRgroove_HitProcess::processID() maps each sensitive gas hit to:
            #   layer 1 -> U strips
            #   layer 2 -> V strips
            if ($group eq "gas1" && $mat eq "urgroove_gas") {
                $detector{"sensitivity"} = $HIT_TYPE;
                $detector{"hit_type"}    = $HIT_TYPE;

                $detector{"identifiers"} =
                    "region manual 1 sector manual " . ($imod + 1) .
                    " chamber manual 1 layer manual 1 component manual 1";
            }

            set_common(\%detector);
            print_checked_det(\%detector);

            $i++;
        }
    }
}

sub get_volume_info
{
    my ($kind, $imod, $group, $mat, $layer_index, $is_flipped) = @_;

    # Future hook for a volume table.
    return (undef, undef);
}

make_urgroove();

1;
